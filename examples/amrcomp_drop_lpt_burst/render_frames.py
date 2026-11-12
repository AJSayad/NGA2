#!/usr/bin/env pvbatch
# Headless, frame-parallel rendering of NGA2 PLIC surfaces and spray particles to PNG (ParaView >= 6).
# Run:  pvbatch render_frames.py            (uses SETTINGS below)
#       pvbatch render_frames.py --help     (every setting is also a flag; a flag overrides its setting)

# =============================================== USER SETTINGS ===============================================
# Relative paths are taken from the folder you launch pvbatch from (here: the case folder).
# None means "not set". Each key matches the flag of the same name (underscores <-> dashes).
SETTINGS = dict(
    # --- data
    plic="amrviz/burst_tau_2.50000E-01",    # PLIC folder (plic.pvd + plic_*.vtp) or .pvd file; None = no PLIC
    spray="amrviz/burst_tau_2.50000E-01",   # folder of plt.part.* plotfiles; None = no particles
    out="frames/tau0.25",                   # output folder for the PNGs
    frames=None,                            # 1-based 'A', 'A:B' or 'A:B:STEP'; None = all frames
    time_tol=1e-6,                          # PLIC/spray time match tolerance, relative to max(1,|t|)
    overwrite=False,                        # True re-renders frames whose PNG already exists
    # --- execution
    np=16,                                  # number of parallel pvbatch workers
    gl="egl",                               # "egl" = GPU; "mesa" = CPU (Kitware binary); "none" = OSMesa cluster builds
    pvbatch=None,                           # pvbatch for the workers; None = the one running this script
    # --- appearance
    res=[1920, 1080],                       # image width, height [px]
    color=[1.0, 1.0, 1.0],                  # solid RGB (0-1) of PLIC and particles
    background=[0.0, 0.0, 0.0],             # background RGB (0-1)
    diameter_array="d",                     # particle diameter array (glyphs drawn at true diameter)
    glyph_res=16,                           # sphere glyph theta/phi resolution
    # --- camera: a .pvcc is applied first; any cam_* / parallel_scale / projection set here or as a flag
    #     overrides the matching .pvcc value (set them to None to use the .pvcc as saved)
    camera_file=None,                       # .pvcc saved from the ParaView GUI
    #     Camera from the earlier imaging script (AMRmpcomp/imaging/PLIC_imaging_pvbatch_projectile.py):
    #     lengths in d0; view direction atan(0.79/4.5) = 9.96 deg upstream (-x) of the perpendicular to the
    #     shot line (x axis), offset along +z. Parallel projection, so the 4.57 d0 camera distance does not
    #     affect the image; parallel_scale = 1.5 -> visible height 3.0 d0 (width 3.0*W/H).
    cam_position=[-0.79, 0.0, 4.5],
    cam_focal_point=[0.0, 0.0, 0.0],        # initial drop centre
    cam_view_up=[0.0, 1.0, 0.0],
    cam_view_angle=30.0,                    # perspective view angle [deg]; unused in parallel projection
    parallel_scale=1.5,                     # parallel-projection half-height [d0]; setting it implies parallel
    projection="parallel",                  # None, "perspective" or "parallel"
)
# =============================================================================================================

DESCRIPTION = """Headless, frame-parallel rendering of NGA2 PLIC surfaces and spray particles to PNG (ParaView >= 6).
Defaults come from SETTINGS at the top of this script; any flag below overrides its setting for that run.

Inputs
  --plic   a folder holding plic.pvd (times) and plic_*.vtp, a .pvd file, or a folder of plic_*.vtp only
  --spray  a folder holding plt.part.* (each needs the top-level 'Header': run make_particle_headers.py)
Either input may be omitted. When both are given, each PLIC frame is paired with the spray plotfile of
nearest time (within --time-tol); unmatched frames are rendered PLIC-only with a warning.

Rendering
  PLIC is drawn as a solid surface; particles as sphere glyphs of true diameter (array --diameter-array)
  in the same solid colour. The camera is fixed for all frames and fully user controlled: load a .pvcc
  saved from the ParaView GUI (File > Save Camera / Adjust Camera > Save) and/or give the individual
  --cam-* flags, which override the corresponding .pvcc values.

Parallelism
  This process only plans the work; it launches --np independent headless pvbatch workers, each rendering
  its own subset of frames (data per frame is small, so frame parallelism beats domain decomposition).
  Existing PNGs are skipped (rerun to resume), images are written atomically, failures are reported per
  frame, and the exit status is nonzero if any requested frame is missing at the end.

Examples
  pvbatch render_frames.py                     # everything from SETTINGS
  pvbatch render_frames.py --frames 1          # SETTINGS, but only frame 1 (camera check)
  pvbatch render_frames.py --plic amrviz/burst_tau_2.50000E-01 --spray amrviz/burst_tau_2.50000E-01 \\
      --out frames/tau0.25 --np 16 --camera-file view.pvcc
  pvbatch render_frames.py --plic amrviz/plic --spray amrviz/spray --out frames/main --np 8 --gl mesa \\
      --cam-position 3 2 4 --cam-focal-point 0.5 0 0 --cam-view-up 0 1 0 --cam-view-angle 30 --frames 1:1
Outputs: <out>/frame_NNNNNN.png (1-based, time ordered), <out>/frames.txt (frame -> time and files),
<out>/logs/worker_NN.log.
"""
__doc__ = DESCRIPTION

import argparse
import bisect
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
    print(f"render_frames: ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(2)


def parse_args(argv):
    ap = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    io = ap.add_argument_group("data")
    io.add_argument("--plic", help="PLIC folder (with plic.pvd or plic_*.vtp) or .pvd file (%(default)s)")
    io.add_argument("--spray", help="folder holding plt.part.* plotfiles (%(default)s)")
    io.add_argument("--out", help="output folder for PNGs (%(default)s)")
    io.add_argument("--frames", help="1-based inclusive selection 'A', 'A:B' or 'A:B:STEP' (%(default)s = all)")
    io.add_argument("--time-tol", type=float,
                    help="PLIC/spray time match tolerance, relative to max(1,|t|) (%(default)s)")
    io.add_argument("--overwrite", action=argparse.BooleanOptionalAction,
                    help="re-render frames whose PNG exists (%(default)s)")
    run = ap.add_argument_group("execution")
    run.add_argument("--np", type=int, help="number of parallel pvbatch workers (%(default)s)")
    run.add_argument("--gl", choices=("egl", "mesa", "none"),
                     help="headless OpenGL: egl = GPU; mesa = CPU via the Mesa bundled with Kitware binaries; "
                          "none = no backend flag, the build's default (use for OSMesa cluster builds) (%(default)s)")
    run.add_argument("--pvbatch", help="pvbatch executable for workers (%(default)s = the one running this)")
    look = ap.add_argument_group("appearance")
    look.add_argument("--res", type=int, nargs=2, metavar=("W", "H"), help="(%(default)s)")
    look.add_argument("--color", type=float, nargs=3, metavar=("R", "G", "B"),
                      help="solid colour of PLIC and particles, 0-1 (%(default)s)")
    look.add_argument("--background", type=float, nargs=3, metavar=("R", "G", "B"),
                      help="background colour, 0-1 (%(default)s)")
    look.add_argument("--diameter-array", help="particle diameter array (%(default)s)")
    look.add_argument("--glyph-res", type=int, help="sphere glyph theta/phi resolution (%(default)s)")
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

    for key in ("overwrite", "np", "gl", "time_tol", "res", "color", "background", "diameter_array", "glyph_res"):
        if getattr(args, key) is None:
            die(f"'{key}' is not set: give it in SETTINGS or as --{key.replace('_', '-')}")
    if not args.out:
        die("no output folder: set SETTINGS['out'] or give --out")
    if not args.plic and not args.spray:
        die("give --plic and/or --spray")
    if args.np < 1:
        die("--np must be >= 1")
    if min(args.res) < 1 or args.glyph_res < 3 or args.time_tol < 0:
        die("--res must be positive, --glyph-res >= 3, --time-tol >= 0")
    for name in ("color", "background"):
        if any(not 0.0 <= c <= 1.0 for c in getattr(args, name)):
            die(f"--{name} components must be in [0, 1]")
    return args


def read_pvd(pvd):
    """Return [(time, file)] from a .pvd, sorted by time; a repeated time keeps its last entry."""
    try:
        root = ET.parse(pvd).getroot()
    except (OSError, ET.ParseError) as err:
        die(f"cannot parse {pvd}: {err}")
    base = os.path.dirname(os.path.abspath(pvd))
    bytime = {}
    for ds in root.iter("DataSet"):
        try:
            t = float(ds.attrib["timestep"])
            f = ds.attrib["file"]
        except (KeyError, ValueError):
            die(f"{pvd}: malformed DataSet entry {ds.attrib}")
        bytime[t] = f if os.path.isabs(f) else os.path.join(base, f)
    if not bytime:
        die(f"{pvd}: no DataSet entries")
    return sorted(bytime.items())


def plic_frames(path):
    """Return [(time or None, vtp path)] in time order."""
    if os.path.isfile(path) and path.endswith(".pvd"):
        entries = read_pvd(path)
    elif os.path.isdir(path):
        pvd = os.path.join(path, "plic.pvd")
        if os.path.isfile(pvd):
            entries = read_pvd(pvd)
        else:
            entries = [(None, f) for f in sorted(glob.glob(os.path.join(path, "plic_*.vtp")))]
            if entries:
                print(f"render_frames: WARNING: {path} has no plic.pvd; PLIC times unknown, file order used")
    else:
        die(f"--plic {path}: not a folder or .pvd file")
    missing = [f for _, f in entries if not os.path.isfile(f)]
    if missing:
        print(f"render_frames: WARNING: {len(missing)} PLIC files listed but missing (skipped), "
              f"first: {missing[0]}")
    entries = [e for e in entries if os.path.isfile(e[1])]
    if not entries:
        die(f"--plic {path}: no PLIC files found")
    return entries


def spray_frames(path):
    """Return [(time, plt.part folder)] in time order; every folder must carry the top-level Header."""
    if not os.path.isdir(path):
        die(f"--spray {path}: not a folder")
    out, noheader = [], []
    for pf in sorted(p for p in glob.glob(os.path.join(path, "plt.part.*")) if os.path.isdir(p)):
        if not os.path.isfile(os.path.join(pf, "Header")):
            noheader.append(pf)
            continue
        try:
            with open(os.path.join(pf, "time")) as f:
                out.append((float(f.read().split()[0]), pf))
        except (OSError, IndexError, ValueError) as err:
            die(f"{pf}: cannot read 'time' ({err})")
    if noheader:
        die(f"{len(noheader)} spray plotfiles lack the top-level 'Header' (first: {noheader[0]}). Run:\n"
            f"    python3 make_particle_headers.py {path}")
    if not out:
        die(f"--spray {path}: no plt.part.* plotfiles found")
    out.sort()
    return out


def build_frames(args):
    """Return the full time-ordered frame list: dicts with index (1-based), time, plic, spray."""
    plic = plic_frames(args.plic) if args.plic else None
    spray = spray_frames(args.spray) if args.spray else None
    if plic is None:
        pairs = [(t, None, s) for t, s in spray]
    elif spray is None:
        pairs = [(t, p, None) for t, p in plic]
    elif plic[0][0] is None:                                   # no PLIC times: pair by order only
        if len(plic) != len(spray):
            die(f"PLIC has no plic.pvd and {len(plic)} files vs {len(spray)} spray plotfiles: cannot pair")
        pairs = [(st, p, s) for (_, p), (st, s) in zip(plic, spray)]
    else:
        stimes = [t for t, _ in spray]
        pairs, unmatched = [], 0
        for t, p in plic:
            k = bisect.bisect_left(stimes, t)
            best = min((j for j in (k - 1, k) if 0 <= j < len(stimes)), key=lambda j: abs(stimes[j] - t))
            ok = abs(stimes[best] - t) <= args.time_tol * max(1.0, abs(t))
            unmatched += not ok
            pairs.append((t, p, spray[best][1] if ok else None))
        if unmatched:
            print(f"render_frames: WARNING: {unmatched}/{len(plic)} PLIC frames have no spray plotfile within "
                  f"--time-tol; they are rendered without particles (see frames.txt)")
    return [dict(index=i + 1, time=t, plic=p, spray=s) for i, (t, p, s) in enumerate(pairs)]


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


def resolve_camera(args):
    """Merge .pvcc values and flags into one validated camera dict."""
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
    if replaced:
        print(f"render_frames: NOTE: {replaced} (from SETTINGS or flags) override the values in "
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


def png_path(out, index):
    return os.path.join(out, f"frame_{index:06d}.png")


def check_paraview_version():
    from paraview import servermanager as sm
    major, minor = sm.vtkSMProxyManager.GetVersionMajor(), sm.vtkSMProxyManager.GetVersionMinor()
    if major < MIN_PARAVIEW_MAJOR:
        die(f"ParaView {major}.{minor} found; this script requires ParaView >= {MIN_PARAVIEW_MAJOR}")
    return f"{major}.{minor}"


# ----------------------------------------------------------------------------------------------- master
def master(args, argv):
    version = check_paraview_version()
    cam = resolve_camera(args)
    frames = build_frames(args)
    chosen = select_frames(frames, args.frames)
    os.makedirs(os.path.join(args.out, "logs"), exist_ok=True)

    table = os.path.join(args.out, "frames.txt")                 # frame -> data map, rewritten each run
    with open(table + ".tmp", "w") as f:
        f.write(f"# {'frame':>6} {'time':>24}  plic  spray\n")
        for fr in frames:
            f.write(f"  {fr['index']:6d} {fr['time'] if fr['time'] is not None else float('nan'):24.16e}  "
                    f"{fr['plic'] or '-'}  {fr['spray'] or '-'}\n")
    os.replace(table + ".tmp", table)

    todo = [fr["index"] for fr in chosen if args.overwrite or not os.path.isfile(png_path(args.out, fr["index"]))]
    print(f"render_frames: ParaView {version}, {len(frames)} frames found, {len(chosen)} selected, "
          f"{len(chosen) - len(todo)} already rendered, {len(todo)} to render", flush=True)
    print(f"render_frames: camera {cam}", flush=True)
    if not todo:
        return 0

    pvbatch = args.pvbatch                                    # else next to the running one, else on PATH
    if not pvbatch:
        here = os.path.join(os.path.dirname(sys.executable), "pvbatch")
        pvbatch = here if os.access(here, os.X_OK) else shutil.which("pvbatch")
    if not pvbatch or not os.access(pvbatch, os.X_OK):
        die(f"pvbatch executable not found ({pvbatch}); pass --pvbatch")
    # egl: GPU; mesa: Kitware-binary bundled Mesa; none: the build's own default (e.g. OSMesa cluster builds)
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
    print(f"render_frames: {nw} workers ({args.gl}), logs in {os.path.join(args.out, 'logs')}", flush=True)

    t0, last = wallclock.time(), -1
    while any(p.poll() is None for p, _ in procs):
        wallclock.sleep(5)
        if args.overwrite:                                        # old PNGs would inflate the count
            continue
        done = sum(os.path.isfile(png_path(args.out, i)) for i in todo)
        if done != last:
            print(f"render_frames: {done}/{len(todo)} rendered ({wallclock.time() - t0:.0f} s)", flush=True)
            last = done
    for _, log in procs:
        log.close()

    bad_workers = [w for w, (p, _) in enumerate(procs) if p.returncode != 0]
    missing = [i for i in todo if not os.path.isfile(png_path(args.out, i))]
    print(f"render_frames: finished in {wallclock.time() - t0:.0f} s: {len(todo) - len(missing)}/{len(todo)} "
          f"frames rendered", flush=True)
    if bad_workers:
        print(f"render_frames: workers with errors: {bad_workers} (see their logs)", flush=True)
    if missing:
        print(f"render_frames: MISSING frames: {missing[:20]}{' ...' if len(missing) > 20 else ''}", flush=True)
    return 1 if (missing or bad_workers) else 0


# ----------------------------------------------------------------------------------------------- worker
def worker(args):
    from paraview import simple as pv
    check_paraview_version()
    cam = resolve_camera(args)
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
    view.UseFXAA = 1
    if not view.GetRenderWindow().GetOffScreenRendering():
        die(f"{tag} render window is not offscreen; refusing to open a display window")

    plic_rep = spray_rep = plic_rd = spray_rd = None
    first_plic = next((f["plic"] for f in mine if f["plic"]), None)
    first_spray = next((f["spray"] for f in mine if f["spray"]), None)
    if first_plic:
        plic_rd = pv.XMLPolyDataReader(FileName=[first_plic])
        plic_rd.CellArrayStatus = []                                # geometry only: less I/O and memory
        plic_rd.PointArrayStatus = []
        plic_rep = pv.Show(plic_rd, view)
        plic_rep.Representation = "Surface"
        plic_rep.ColorArrayName = ["POINTS", ""]                  # solid colour (ColorBy(None) fails in 6.1)
        plic_rep.AmbientColor = plic_rep.DiffuseColor = list(args.color)
    if first_spray:
        spray_rd = pv.AMReXBoxLibParticlesReader(FileNames=[first_spray])   # ParaView 6: FileNames list
        spray_rd.PointArrayStatus = [args.diameter_array]
        glyph = pv.Glyph(Input=spray_rd, GlyphType="Sphere")
        glyph.GlyphType.Radius = 0.5                                 # unit-diameter sphere scaled by d
        glyph.GlyphType.ThetaResolution = glyph.GlyphType.PhiResolution = args.glyph_res
        glyph.ScaleArray = ["POINTS", args.diameter_array]
        glyph.ScaleFactor = 1.0
        glyph.GlyphMode = "All Points"
        spray_rep = pv.Show(glyph, view)
        spray_rep.Representation = "Surface"
        spray_rep.ColorArrayName = ["POINTS", ""]                  # solid colour (ColorBy(None) fails in 6.1)
        spray_rep.AmbientColor = spray_rep.DiffuseColor = list(args.color)

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
        out = png_path(args.out, fr["index"])
        tmp = out[:-4] + f".tmp{args.worker}.png"
        t0 = wallclock.time()
        try:
            if plic_rep is not None:
                plic_rep.Visibility = int(fr["plic"] is not None)
                if fr["plic"]:
                    plic_rd.FileName = [fr["plic"]]
            npart = 0
            if spray_rep is not None:
                if fr["spray"]:
                    spray_rd.FileNames = [fr["spray"]]
                    spray_rd.UpdatePipelineInformation()
                    spray_rd.PointArrayStatus = [args.diameter_array]
                    spray_rd.UpdatePipeline()
                    info = spray_rd.GetDataInformation()
                    npart = info.GetNumberOfPoints()
                    if npart and spray_rd.PointData.GetArray(args.diameter_array) is None:
                        raise RuntimeError(f"{fr['spray']} has no particle array '{args.diameter_array}'")
                spray_rep.Visibility = int(npart > 0)
            pv.SaveScreenshot(tmp, view, ImageResolution=list(args.res), TransparentBackground=0)
            if not os.path.isfile(tmp) or os.path.getsize(tmp) == 0:
                raise RuntimeError("SaveScreenshot wrote no image")
            os.replace(tmp, out)
            print(f"{tag} frame {fr['index']:6d}  t={fr['time']}  particles={npart}  "
                  f"{wallclock.time() - t0:.2f} s", flush=True)
        except Exception as err:                                     # report, keep going
            nfail += 1
            print(f"{tag} FAILED frame {fr['index']} ({fr['plic']}, {fr['spray']}): {err!r}", flush=True)
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
