#!/usr/bin/env pvbatch
# Headless, frame-parallel rendering of xy-plane (constant-z) slices of NGA2 AMReX field plotfiles (ParaView >= 6).
# Run:  pvbatch render_slices.py            (uses SETTINGS and FIELDS below)
#       pvbatch render_slices.py --help     (every SETTINGS entry is also a flag; a flag overrides its setting)

# =============================================== USER SETTINGS ===============================================
# Relative paths are taken from the folder you launch pvbatch from (here: the case folder).
# None means "not set". Each key matches the flag of the same name (underscores <-> dashes).
SETTINGS = dict(
    # --- data
    data="amrviz/drop",                     # folder of plt.nga2.cell.* plotfiles, or a single plotfile folder
    out="slices",                           # output folder; images go to <out>/<FIELD>_slice/
    fields=None,                            # subset of FIELDS to render, e.g. ["RHOG"]; None = all of FIELDS
    slice_z=None,                           # slice plane z [d0]; None = domain centre. A z on a cell face is
                                            # moved up half a finest cell so every point samples one cell
    frames=None,                            # 1-based 'A', 'A:B' or 'A:B:STEP'; None = all plotfiles
    overwrite=False,                        # True re-renders frames whose PNGs already exist
    # --- execution
    np=4,                                   # parallel pvbatch workers (each needs ~4 GB for a 17 GB plotfile)
    gl="egl",                               # "egl" = GPU; "mesa" = CPU (Kitware binary); "none" = OSMesa cluster builds
    pvbatch=None,                           # pvbatch for the workers; None = the one running this script
    # --- appearance
    res=[1920, 1080],                       # image width, height [px]
    background=[0.0, 0.0, 0.0],             # background RGB (0-1), seen outside the domain
    # --- camera: a .pvcc is applied first; any cam_* / parallel_scale / projection set here or as a flag
    #     overrides the matching .pvcc value (set them to None to use the .pvcc as saved)
    camera_file=None,                       # .pvcc saved from the ParaView GUI
    cam_position=[3.0, 0.0, 10.0],          # looking down -z onto the xy plane
    cam_focal_point=[3.0, 0.0, 0.0],
    cam_view_up=[0.0, 1.0, 0.0],
    cam_view_angle=30.0,                    # perspective view angle [deg]; unused in parallel projection
    parallel_scale=4.5,                     # parallel-projection half-height [d0] -> x in [-5, 11] at 16:9
    projection="parallel",                  # None, "perspective" or "parallel"
)

# Colour map of each rendered variable (name as in the plotfile Header). Fixed ranges keep colours comparable
# across frames. log=True uses a logarithmic colour scale (range must be > 0).
FIELDS = {
    "RHOG": dict(range=[0.15, 1.85],  log=False, colormap="Cool to Warm", title="Gas density"),
    "PG":   dict(range=[0.003, 0.5],  log=True,  colormap="Cool to Warm", title="Gas pressure"),
    "TG":   dict(range=[0.13, 1.5],   log=False, colormap="Cool to Warm", title="Gas temperature"),
}
# =============================================================================================================

DESCRIPTION = """Headless, frame-parallel rendering of constant-z slices of NGA2 AMReX field plotfiles (ParaView >= 6).
Defaults come from SETTINGS/FIELDS at the top of this script; any flag below overrides its setting for that run.

Data
  --data is a folder of plt.nga2.cell.* plotfiles (particle plotfiles in the same folder are ignored) or one
  plotfile folder. Frames are the plotfiles in time order; times come from each plotfile's Header. All AMR
  levels are read (only the rendered variables), and the plane is cut with ParaView's AMR-aware
  AxisAlignedSlice, which skips coarse cells covered by finer levels.

Output
  <out>/<FIELD>_slice/frame_NNNNNN_t<time>.png for every field, <out>/frames.txt (frame -> time, plotfile),
  <out>/logs/worker_NN.log. Each worker reads a plotfile once and writes all fields from it. Existing PNGs
  are skipped (rerun to resume), images are written atomically, and the exit status is nonzero if any
  requested image is missing at the end.

Examples
  pvbatch render_slices.py                                  # everything from SETTINGS/FIELDS
  pvbatch render_slices.py --frames 1 --fields RHOG         # one image, to check camera and range
  pvbatch render_slices.py --data ~/monitor_hang_drop6 --out slices_test --np 1
"""
__doc__ = DESCRIPTION

import argparse
import glob
import os
import shutil
import signal
import subprocess
import sys
import time as wallclock
import xml.etree.ElementTree as ET

MIN_PARAVIEW_MAJOR = 6


# ----------------------------------------------------------------------------------------------- inputs
def die(msg):
    print(f"render_slices: ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(2)


def check_fields_table():
    """Validate FIELDS once; return it with ranges as floats."""
    if not isinstance(FIELDS, dict) or not FIELDS:
        die("FIELDS must be a non-empty dict of {variable: dict(range=..., log=..., colormap=..., title=...)}")
    out = {}
    for name, spec in FIELDS.items():
        if not isinstance(spec, dict) or set(spec) - {"range", "log", "colormap", "title"}:
            die(f"FIELDS['{name}'] must be a dict with keys range, log, colormap, title")
        try:
            lo, hi = (float(v) for v in spec["range"])
        except (KeyError, TypeError, ValueError):
            die(f"FIELDS['{name}']['range'] must be [min, max]")
        log = spec.get("log", False)
        if not isinstance(log, bool):
            die(f"FIELDS['{name}']['log'] must be True or False")
        if not lo < hi:
            die(f"FIELDS['{name}']['range'] = {spec['range']}: min must be < max")
        if log and lo <= 0:
            die(f"FIELDS['{name}']: a log colour scale needs a range > 0, got {spec['range']}")
        out[name] = dict(range=[lo, hi], log=log, colormap=str(spec.get("colormap", "Cool to Warm")),
                         title=str(spec.get("title", name)))
    return out


def parse_args(argv):
    ap = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    io = ap.add_argument_group("data")
    io.add_argument("--data", help="folder of plt.nga2.cell.* plotfiles, or one plotfile (%(default)s)")
    io.add_argument("--out", help="output folder (%(default)s)")
    io.add_argument("--fields", nargs="+", help="variables to render, keys of FIELDS (%(default)s = all)")
    io.add_argument("--slice-z", type=float, help="slice plane z (%(default)s = domain centre)")
    io.add_argument("--frames", help="1-based inclusive selection 'A', 'A:B' or 'A:B:STEP' (%(default)s = all)")
    io.add_argument("--overwrite", action=argparse.BooleanOptionalAction,
                    help="re-render frames whose PNGs exist (%(default)s)")
    run = ap.add_argument_group("execution")
    run.add_argument("--np", type=int, help="number of parallel pvbatch workers (%(default)s)")
    run.add_argument("--gl", choices=("egl", "mesa", "none"),
                     help="headless OpenGL: egl = GPU; mesa = CPU via the Mesa bundled with Kitware binaries; "
                          "none = no backend flag, the build's default (use for OSMesa cluster builds) (%(default)s)")
    run.add_argument("--pvbatch", help="pvbatch executable for workers (%(default)s = the one running this)")
    look = ap.add_argument_group("appearance")
    look.add_argument("--res", type=int, nargs=2, metavar=("W", "H"), help="(%(default)s)")
    look.add_argument("--background", type=float, nargs=3, metavar=("R", "G", "B"),
                      help="background colour, 0-1 (%(default)s)")
    cam = ap.add_argument_group("camera (these override --camera-file values)")
    cam.add_argument("--camera-file", help=".pvcc camera file saved from the ParaView GUI (%(default)s)")
    cam.add_argument("--cam-position", type=float, nargs=3, metavar=("X", "Y", "Z"), help="(%(default)s)")
    cam.add_argument("--cam-focal-point", type=float, nargs=3, metavar=("X", "Y", "Z"), help="(%(default)s)")
    cam.add_argument("--cam-view-up", type=float, nargs=3, metavar=("X", "Y", "Z"), help="(%(default)s)")
    cam.add_argument("--cam-view-angle", type=float, help="perspective view angle [deg] (%(default)s)")
    cam.add_argument("--parallel-scale", type=float,
                     help="parallel projection half-height, implies parallel (%(default)s)")
    cam.add_argument("--projection", choices=("perspective", "parallel"),
                     help="force the projection type (%(default)s)")
    ap.add_argument("--worker", type=int, help=argparse.SUPPRESS)        # internal: worker id
    ap.add_argument("--worker-frames", help=argparse.SUPPRESS)          # internal: comma list of frame indices

    # SETTINGS become the defaults; argparse does not type-check defaults, so check/convert them here
    actions = {a.dest: a for a in ap._actions if a.dest not in ("help", "worker", "worker_frames")}
    unknown = sorted(set(SETTINGS) - set(actions))
    if unknown:
        die(f"unknown SETTINGS key(s) {unknown}; valid keys: {sorted(actions)}")
    defaults = {}
    for key, val in SETTINGS.items():
        act = actions[key]
        try:
            if val is None:
                pass
            elif key == "overwrite":
                if not isinstance(val, bool):
                    raise ValueError("must be True or False")
            elif key == "fields":
                if not isinstance(val, (list, tuple)) or not all(isinstance(v, str) for v in val):
                    raise ValueError("must be a list of variable names (or None)")
                val = list(val)
            elif act.type is None and act.choices is None:
                if not isinstance(val, str):
                    raise ValueError("must be a string (or None)")
            elif isinstance(act.nargs, int):
                if not isinstance(val, (list, tuple)) or len(val) != act.nargs:
                    raise ValueError(f"needs a list of {act.nargs} values")
                val = [act.type(v) for v in val]
            elif act.type is not None:
                if isinstance(val, bool):
                    raise ValueError("is a bool")
                val = act.type(val)
            if val is not None and act.choices is not None and val not in act.choices:
                raise ValueError(f"must be one of {list(act.choices)}")
        except (TypeError, ValueError) as err:
            die(f"SETTINGS['{key}'] = {SETTINGS[key]!r}: {err}")
        defaults[key] = val
    ap.set_defaults(**defaults)
    args = ap.parse_args(argv)

    for key in ("data", "out", "overwrite", "np", "gl", "res", "background"):
        if getattr(args, key) is None:
            die(f"'{key}' is not set: give it in SETTINGS or as --{key.replace('_', '-')}")
    if args.np < 1:
        die("--np must be >= 1")
    if min(args.res) < 1:
        die("--res must be positive")
    if any(not 0.0 <= c <= 1.0 for c in args.background):
        die("--background components must be in [0, 1]")
    table = check_fields_table()
    names = args.fields or list(table)
    bad = [n for n in names if n not in table]
    if bad:
        die(f"field(s) {bad} not in FIELDS {list(table)}; add them to FIELDS with a range")
    if len(set(names)) != len(names):
        die(f"--fields lists a variable twice: {names}")
    args.fields = names
    args.field_spec = {n: table[n] for n in names}
    return args


def read_plotfile_header(pf):
    """Return dict(time, names, finest, lo, hi, dx) from an AMReX plotfile Header, or None if not a field file."""
    try:
        with open(os.path.join(pf, "Header")) as f:
            lines = f.read().splitlines()
        if not lines[0].startswith("HyperCLaw"):
            return None
        nvar = int(lines[1])
        if nvar == 0:                                              # particle plotfile Header: no mesh data
            return None
        names = [s.strip() for s in lines[2:2 + nvar]]
        i = 2 + nvar
        dim, t, finest = int(lines[i]), float(lines[i + 1]), int(lines[i + 2])
        lo = [float(v) for v in lines[i + 3].split()]
        hi = [float(v) for v in lines[i + 4].split()]
        # i+5 refinement ratios, i+6 domain boxes, i+7 level steps, then one cell-size line per level
        dx = [[float(v) for v in lines[i + 8 + L].split()] for L in range(finest + 1)]
    except (OSError, IndexError, ValueError) as err:
        die(f"{pf}/Header: cannot parse ({err})")
    if dim != 3 or len(lo) != 3 or len(hi) != 3 or any(len(d) != 3 for d in dx):
        die(f"{pf}: expected a 3D plotfile (dimension {dim})")
    return dict(time=t, names=names, finest=finest, lo=lo, hi=hi, dx=dx)


def slice_plane_z(hdr, requested, pf):
    """Plane z: requested or domain centre; nudged up half a finest cell if it lies on a finest-level cell face."""
    z = 0.5 * (hdr["lo"][2] + hdr["hi"][2]) if requested is None else requested
    if not hdr["lo"][2] <= z <= hdr["hi"][2]:
        die(f"{pf}: slice z={z} is outside the domain z in [{hdr['lo'][2]}, {hdr['hi'][2]}]")
    dz = hdr["dx"][-1][2]
    n = (z - hdr["lo"][2]) / dz
    if abs(n - round(n)) < 1e-6:                                   # on a face of every level's cells
        z += 0.5 * dz if z + 0.5 * dz <= hdr["hi"][2] else -0.5 * dz
    return z


def build_frames(args):
    """Return the time-ordered frame list: dicts with index (1-based), time, plotfile, header, z."""
    path = args.data
    if not os.path.isdir(path):
        die(f"--data {path}: not a folder")
    if os.path.isfile(os.path.join(path, "Header")):
        candidates = [path]
    else:
        candidates = sorted(p for p in glob.glob(os.path.join(path, "plt*")) if os.path.isdir(p))
    items = []
    for pf in candidates:
        hdr = read_plotfile_header(pf)
        if hdr is not None:
            items.append((hdr["time"], pf, hdr))
    if not items:
        die(f"--data {path}: no AMReX field plotfiles found")
    items.sort(key=lambda it: it[0])
    for t, pf, hdr in items:
        missing = [f for f in args.fields if f not in hdr["names"]]
        if missing:
            die(f"{pf}: variable(s) {missing} not in this plotfile; it has {hdr['names']}")
    tags = [time_tag(t) for t, _, _ in items]
    if len(set(tags)) != len(tags):
        die("two plotfiles have the same time to the file-name precision; cannot name images uniquely")
    return [dict(index=i + 1, time=t, plotfile=pf, header=hdr, z=slice_plane_z(hdr, args.slice_z, pf))
            for i, (t, pf, hdr) in enumerate(items)]


def time_tag(t):
    return f"{t:.6f}"


def png_path(out, field, fr):
    return os.path.join(out, f"{field}_slice", f"frame_{fr['index']:06d}_t{time_tag(fr['time'])}.png")


def select_frames(frames, spec):
    if not spec:
        return frames
    try:
        parts = [int(x) for x in spec.split(":")]
        if len(parts) > 3 or any(x < 1 for x in parts):
            raise ValueError
    except ValueError:
        die(f"--frames '{spec}': expected positive integers 'A', 'A:B' or 'A:B:STEP'")
    a = parts[0]
    b = parts[1] if len(parts) > 1 else a
    step = parts[2] if len(parts) > 2 else 1
    chosen = [f for f in frames if a <= f["index"] <= b and (f["index"] - a) % step == 0]
    if not chosen:
        die(f"--frames '{spec}' selects nothing (there are {len(frames)} frames)")
    return chosen


def read_pvcc(path):
    """Return {property: [values]} from a ParaView .pvcc camera file."""
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as err:
        die(f"cannot parse camera file {path}: {err}")
    props = {}
    for prop in root.iter("Property"):
        vals = [float(e.attrib["value"]) for e in prop.iter("Element") if "value" in e.attrib]
        if vals:
            props[prop.attrib.get("name")] = vals
    if "CameraPosition" not in props:
        die(f"{path}: no CameraPosition found; is this a .pvcc camera file?")
    return props


def resolve_camera(args, verbose):
    """Merge .pvcc values and SETTINGS/flags into one validated camera dict."""
    cam = {}
    if args.camera_file:
        p = read_pvcc(args.camera_file)
        for key, name, n in (("position", "CameraPosition", 3), ("focal_point", "CameraFocalPoint", 3),
                             ("view_up", "CameraViewUp", 3), ("view_angle", "CameraViewAngle", 1),
                             ("parallel_scale", "CameraParallelScale", 1)):
            if name in p:
                cam[key] = p[name] if n == 3 else p[name][0]
        if "CameraParallelProjection" in p:
            cam["parallel"] = bool(p["CameraParallelProjection"][0])
    from_file = set(cam)
    overrides = {"cam_" + k: getattr(args, "cam_" + k) for k in ("position", "focal_point", "view_up", "view_angle")}
    overrides.update(parallel_scale=args.parallel_scale, projection=args.projection)
    replaced = [k for k, v in overrides.items() if v is not None and
                (k.removeprefix("cam_") in from_file or (k == "projection" and "parallel" in from_file))]
    if replaced and verbose:
        print(f"render_slices: NOTE: {replaced} (from SETTINGS or flags) override the values in "
              f"{args.camera_file}; set them to None to use the .pvcc as saved", flush=True)
    for key in ("position", "focal_point", "view_up", "view_angle"):
        val = getattr(args, "cam_" + key)
        if val is not None:
            cam[key] = val
    if args.parallel_scale is not None:
        cam["parallel_scale"] = args.parallel_scale
        cam["parallel"] = True
    if args.projection:
        cam["parallel"] = args.projection == "parallel"
    cam.setdefault("parallel", False)
    cam.setdefault("view_up", [0.0, 1.0, 0.0])
    cam.setdefault("view_angle", 30.0)

    if "position" not in cam or "focal_point" not in cam:
        die("camera needs a position and a focal point: give --camera-file or --cam-position and "
            "--cam-focal-point")
    d = [f - p for p, f in zip(cam["position"], cam["focal_point"])]
    dn = sum(x * x for x in d) ** 0.5
    un = sum(x * x for x in cam["view_up"]) ** 0.5
    if dn == 0.0:
        die("camera position equals the focal point")
    if un == 0.0:
        die("camera view-up is the zero vector")
    if abs(sum(a * b for a, b in zip(d, cam["view_up"]))) > (1 - 1e-6) * dn * un:
        die("camera view-up is parallel to the viewing direction")
    if cam["parallel"] and not cam.get("parallel_scale", 0) > 0:
        die("parallel projection needs a positive --parallel-scale (or one in the .pvcc)")
    if not cam["parallel"] and not 0.0 < cam["view_angle"] < 180.0:
        die("--cam-view-angle must be in (0, 180) degrees")
    return cam


def check_paraview_version():
    from paraview import servermanager as sm
    major, minor = sm.vtkSMProxyManager.GetVersionMajor(), sm.vtkSMProxyManager.GetVersionMinor()
    if major < MIN_PARAVIEW_MAJOR:
        die(f"ParaView {major}.{minor} found; this script requires ParaView >= {MIN_PARAVIEW_MAJOR}")
    return f"{major}.{minor}"


def check_colormaps(spec):
    from paraview import simple as pv
    presets = pv.servermanager.vtkSMTransferFunctionPresets.GetInstance()
    for name, s in spec.items():
        if not presets.HasPreset(s["colormap"]):
            die(f"FIELDS['{name}']: unknown ParaView colour map '{s['colormap']}'")


# ----------------------------------------------------------------------------------------------- master
def master(args, argv):
    version = check_paraview_version()
    check_colormaps(args.field_spec)
    cam = resolve_camera(args, verbose=True)
    frames = build_frames(args)
    chosen = select_frames(frames, args.frames)
    for f in args.fields:
        os.makedirs(os.path.join(args.out, f"{f}_slice"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "logs"), exist_ok=True)

    table = os.path.join(args.out, "frames.txt")                 # frame -> data map, rewritten each run
    with open(table + ".tmp", "w") as f:
        f.write(f"# {'frame':>6} {'time':>24} {'slice_z':>14}  plotfile\n")
        for fr in frames:
            f.write(f"  {fr['index']:6d} {fr['time']:24.16e} {fr['z']:14.8f}  {fr['plotfile']}\n")
    os.replace(table + ".tmp", table)

    def done(fr):
        return all(os.path.isfile(png_path(args.out, fld, fr)) for fld in args.fields)
    todo = [fr["index"] for fr in chosen if args.overwrite or not done(fr)]
    print(f"render_slices: script {os.path.abspath(__file__)}", flush=True)
    print(f"render_slices: data {os.path.abspath(args.data)} -> out {os.path.abspath(args.out)}, "
          f"overwrite {args.overwrite}, fields {args.fields}", flush=True)
    print(f"render_slices: ParaView {version}, {len(frames)} plotfiles found, {len(chosen)} selected, "
          f"{len(chosen) - len(todo)} already rendered, {len(todo)} to render", flush=True)
    print(f"render_slices: slice z {sorted({fr['z'] for fr in chosen})}, camera {cam}", flush=True)
    if not todo:
        return 0

    pvbatch = args.pvbatch                                    # else next to the running one, else on PATH
    if not pvbatch:
        here = os.path.join(os.path.dirname(sys.executable), "pvbatch")
        pvbatch = here if os.access(here, os.X_OK) else shutil.which("pvbatch")
    if not pvbatch or not os.access(pvbatch, os.X_OK):
        die(f"pvbatch executable not found ({pvbatch}); pass --pvbatch")
    backend = {"egl": ["--opengl-window-backend", "EGL"], "mesa": ["--mesa"], "none": []}[args.gl]
    nw = min(args.np, len(todo))
    env = dict(os.environ)
    if args.gl != "egl":                                      # CPU rendering: share the allocated cores
        ncpu = len(os.sched_getaffinity(0))                   # respects the batch job's CPU set
        env.setdefault("LP_NUM_THREADS", str(max(1, ncpu // nw)))
    script = os.path.abspath(__file__)

    procs = []
    def stop(signum, _frame):                                     # do not orphan workers
        for p, _ in procs:
            if p.poll() is None:
                p.terminate()
        sys.exit(128 + signum)
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    for w in range(nw):
        mine = todo[w::nw]
        log = open(os.path.join(args.out, "logs", f"worker_{w:02d}.log"), "w")
        cmd = [pvbatch, *backend, "--force-offscreen-rendering", script, *argv,
               "--worker", str(w), "--worker-frames", ",".join(map(str, mine))]
        procs.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env), log))
    print(f"render_slices: {nw} workers ({args.gl}), logs in {os.path.join(args.out, 'logs')}", flush=True)

    byindex = {fr["index"]: fr for fr in frames}
    t0, last = wallclock.time(), -1
    while any(p.poll() is None for p, _ in procs):
        wallclock.sleep(10)
        if args.overwrite:                                        # old PNGs would inflate the count
            continue
        n = sum(done(byindex[i]) for i in todo)
        if n != last:
            print(f"render_slices: {n}/{len(todo)} frames rendered ({wallclock.time() - t0:.0f} s)", flush=True)
            last = n
    for _, log in procs:
        log.close()

    bad_workers = [w for w, (p, _) in enumerate(procs) if p.returncode != 0]
    missing = [i for i in todo if not done(byindex[i])]
    print(f"render_slices: finished in {wallclock.time() - t0:.0f} s: {len(todo) - len(missing)}/{len(todo)} "
          f"frames rendered ({len(args.fields)} images each)", flush=True)
    if bad_workers:
        print(f"render_slices: workers with errors: {bad_workers} (see their logs)", flush=True)
    if missing:
        print(f"render_slices: frames with MISSING images: {missing[:20]}{' ...' if len(missing) > 20 else ''}",
              flush=True)
    return 1 if (missing or bad_workers) else 0


# ----------------------------------------------------------------------------------------------- worker
def worker(args):
    from paraview import simple as pv
    check_paraview_version()
    cam = resolve_camera(args, verbose=False)
    frames = {fr["index"]: fr for fr in build_frames(args)}
    try:
        mine = [frames[int(i)] for i in args.worker_frames.split(",")]
    except (KeyError, ValueError):
        die(f"worker {args.worker}: bad frame assignment {args.worker_frames}")
    tag = f"worker {args.worker}:"

    view = pv.CreateRenderView()
    view.ViewSize = list(args.res)
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.BackgroundColorMode = "Single Color"
    view.Background = list(args.background)
    if not view.GetRenderWindow().GetOffScreenRendering():
        die(f"{tag} render window is not offscreen; refusing to open a display window")

    reader = pv.AMReXBoxLibGridReader(FileNames=[mine[0]["plotfile"]])     # ParaView 6: FileNames list
    reader.UpdatePipelineInformation()
    reader.Level = mine[0]["header"]["finest"]                     # default reads only levels 0-1
    reader.CellArrayStatus = list(args.fields)                     # read only the rendered variables
    cut = pv.AxisAlignedSlice(Input=reader)                        # AMR aware: skips covered coarse cells
    cut.CutFunction.Normal = [0.0, 0.0, 1.0]
    cut.CutFunction.Origin = [0.0, 0.0, mine[0]["z"]]
    cut.Level = mine[0]["header"]["finest"]                        # default slices level 0 only
    rep = pv.Show(cut, view)
    rep.Representation = "Surface"
    rep.Ambient, rep.Diffuse, rep.Specular = 1.0, 0.0, 0.0         # unlit: pixels show the colour map exactly

    luts = {}
    for name, s in args.field_spec.items():
        lut = pv.GetColorTransferFunction(name)
        bar = pv.GetScalarBar(lut, view)
        bar.Title, bar.ComponentTitle = s["title"], ""
        bar.RangeLabelFormat = "%-#.3g"                            # end labels show the FIELDS range exactly
        bar.Visibility = 0
        luts[name] = (lut, bar)

    def color_by(name):
        """Colour by one field with its fixed FIELDS map. Must run after ColorBy, which re-initialises the
        colour map's range in ParaView 6.1 (setting it beforehand is silently lost)."""
        s = args.field_spec[name]
        lut, bar = luts[name]
        pv.ColorBy(rep, ("CELLS", name))
        rep.LookupTable = lut
        lut.AutomaticRescaleRangeMode = "Never"
        lut.UseLogScale = 0
        lut.ApplyPreset(s["colormap"], True)
        lut.RescaleTransferFunction(*s["range"])
        if s["log"]:
            lut.MapControlPointsToLogSpace()
            lut.UseLogScale = 1
        for other, (_, b) in luts.items():
            b.Visibility = int(other == name)
        return lut

    pv.Render(view)                                                  # first render may reset the camera,
    view.CameraPosition = list(cam["position"])                     # so the fixed camera is set after it
    view.CameraFocalPoint = list(cam["focal_point"])
    view.CameraViewUp = list(cam["view_up"])
    view.CameraViewAngle = cam["view_angle"]
    view.CameraParallelProjection = int(cam["parallel"])
    if cam["parallel"]:
        view.CameraParallelScale = cam["parallel_scale"]

    nfail = 0
    for fr in mine:
        t0 = wallclock.time()
        tmps = []
        try:
            reader.FileNames = [fr["plotfile"]]
            reader.UpdatePipelineInformation()
            reader.Level = fr["header"]["finest"]
            reader.CellArrayStatus = list(args.fields)
            cut.CutFunction.Origin = [0.0, 0.0, fr["z"]]
            cut.Level = fr["header"]["finest"]
            cut.UpdatePipeline()
            if cut.GetDataInformation().GetNumberOfCells() == 0:
                raise RuntimeError(f"slice at z={fr['z']} is empty")
            for name in args.fields:
                lut = color_by(name)
                lo, hi = args.field_spec[name]["range"]
                got = (lut.RGBPoints[0], lut.RGBPoints[-4])
                if abs(got[0] - lo) > 1e-9 * abs(lo) + 1e-300 or abs(got[1] - hi) > 1e-9 * abs(hi) + 1e-300:
                    raise RuntimeError(f"{name} colour range is {got}, expected {(lo, hi)}")
                out = png_path(args.out, name, fr)
                tmp = out[:-4] + f".tmp{args.worker}.png"
                tmps.append(tmp)
                pv.SaveScreenshot(tmp, view, ImageResolution=list(args.res), TransparentBackground=0)
                if not os.path.isfile(tmp) or os.path.getsize(tmp) == 0:
                    raise RuntimeError(f"SaveScreenshot wrote no image for {name}")
                os.replace(tmp, out)
            print(f"{tag} frame {fr['index']:6d}  t={fr['time']}  {fr['plotfile']}  "
                  f"{wallclock.time() - t0:.1f} s", flush=True)
        except Exception as err:                                     # report, keep going
            nfail += 1
            print(f"{tag} FAILED frame {fr['index']} ({fr['plotfile']}): {err!r}", flush=True)
            for tmp in tmps:
                if os.path.exists(tmp):
                    os.remove(tmp)
    print(f"{tag} done, {len(mine) - nfail}/{len(mine)} frames rendered", flush=True)
    pv.Delete(view)          # ParaView 6.1 aborts at exit (double free) if a render view is still alive
    del view
    return 1 if nfail else 0


if __name__ == "__main__":
    argv = sys.argv[1:]
    args = parse_args(argv)
    sys.exit(worker(args) if args.worker is not None else master(args, argv))
