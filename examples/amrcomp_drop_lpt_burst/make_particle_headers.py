#!/usr/bin/env python3
"""Add the top-level AMReX 'Header' that ParaView needs to NGA2 particle plotfiles.

Temporary fix: amrlptviz writes plt.part.NNNNNN/{particles/, time} but no top-level
'Header'. Without it, ParaView's AMReX/BoxLib particles reader rejects the folder and
reports time 0.0 for every file, so a stack of plt.part.* has no time steps.
This script writes a minimal Header (format version, 0 mesh variables, dimension, time)
into each plt.part.* folder, taking the time from that folder's 'time' file.

Usage:
    python3 make_particle_headers.py amrviz/spray amrviz/burst_tau_*
    python3 make_particle_headers.py --force amrviz/spray     # overwrite existing Headers
Arguments may be series folders (containing plt.part.*) or single plt.part.* folders.
"""
import argparse
import glob
import os
import sys


def plotfiles(path):
    """Return the plt.part.* folders given a series folder or a single plotfile folder."""
    if os.path.basename(os.path.normpath(path)).startswith("plt.part."):
        return [path]
    return sorted(p for p in glob.glob(os.path.join(path, "plt.part.*")) if os.path.isdir(p))


def write_header(pf, force):
    """Write pf/Header; return 'written', 'skipped' or an error message."""
    header = os.path.join(pf, "Header")
    if os.path.exists(header) and not force:
        return "skipped"
    try:
        with open(os.path.join(pf, "time")) as f:
            time = f.read().split()[0]                 # keep the full-precision text as written
        with open(os.path.join(pf, "particles", "Header")) as f:
            dim = f.read().split()[1]                  # line 2 of the particle Header is the dimension
        float(time); int(dim)                          # fail loudly on malformed input
    except (OSError, IndexError, ValueError) as err:
        return f"ERROR ({err})"
    with open(header, "w") as f:
        f.write(f"HyperCLaw-V1.1\n0\n{dim}\n{time}\n")
    return "written"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="series folders or plt.part.* folders")
    ap.add_argument("--force", action="store_true", help="overwrite an existing Header")
    args = ap.parse_args()
    nerr = 0
    for path in args.paths:
        pfs = plotfiles(path)
        if not pfs:
            print(f"{path}: no plt.part.* folders found")
            nerr += 1
            continue
        counts = {"written": 0, "skipped": 0}
        for pf in pfs:
            status = write_header(pf, args.force)
            if status in counts:
                counts[status] += 1
            else:
                print(f"{pf}: {status}")
                nerr += 1
        print(f"{path}: {counts['written']} written, {counts['skipped']} skipped (Header exists), "
              f"{len(pfs) - sum(counts.values())} errors")
    sys.exit(1 if nerr else 0)


if __name__ == "__main__":
    main()
