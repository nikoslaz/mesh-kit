#!/usr/bin/env python3
"""Run meshkit_cli over every .msh file in a folder and report mesh quality.

For each mesh it writes a config file, runs the C++ tool (read -> smooth ->
write), then prints the minimum triangle angle before/after smoothing and a
histogram of the per-triangle minimum angles.

    python3 scripts/batch_quality.py
    python3 scripts/batch_quality.py --data tests/data --smooth 20 --writer nastran
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

BINS = [0, 10, 20, 30, 40, 50, 60]          # degrees; an equilateral triangle is 60
BAR_WIDTH = 30


def run_one(cli, msh, out_dir, writer, smooth):
    """Convert one mesh. Returns (stats dict, before angles, after angles)."""
    ext = {"vtk": "vtk", "nastran": "bdf"}[writer]
    stem = out_dir / msh.stem
    cfg = stem.with_suffix(".cfg")
    quality = stem.with_name(msh.stem + "_quality.csv")
    cfg.write_text(
        "reader = gmsh\n"
        f"writer = {writer}\n"
        f"input = {msh}\n"
        f"output = {stem.with_suffix('.' + ext)}\n"
        f"smooth = {smooth}\n"
        f"quality = {quality}\n"
    )

    result = subprocess.run([str(cli), str(cfg)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"exit code {result.returncode}")

    stats = {}
    for line in result.stdout.splitlines():           # lines look like "key value"
        key, _, value = line.partition(" ")
        stats[key] = float(value)

    before, after = [], []
    with open(quality, newline="") as f:
        for row in csv.DictReader(f):
            before.append(float(row["min_angle_before"]))
            after.append(float(row["min_angle_after"]))
    return stats, before, after


def histogram(values):
    counts = [0] * (len(BINS) - 1)
    for v in values:
        for i in range(len(counts)):
            if v < BINS[i + 1] or i == len(counts) - 1:
                counts[i] += 1
                break
    return counts


def print_histograms(before, after):
    hb, ha = histogram(before), histogram(after)
    biggest = max(hb + ha) or 1
    print(f"    {'min angle':>9}  {'before':<{BAR_WIDTH + 6}} after")
    for i in range(len(hb)):
        label = f"{BINS[i]:>2}-{BINS[i + 1]:<2}°"
        bar_b = "#" * round(hb[i] / biggest * BAR_WIDTH)
        bar_a = "#" * round(ha[i] / biggest * BAR_WIDTH)
        print(f"    {label:>9}  {bar_b:<{BAR_WIDTH}} {hb[i]:>5}  {bar_a:<{BAR_WIDTH}} {ha[i]:>5}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cli", default="build/meshkit_cli", help="path to meshkit_cli")
    ap.add_argument("--data", default="tests/data", help="folder with .msh files")
    ap.add_argument("--out", default="output/batch", help="folder for outputs")
    ap.add_argument("--writer", default="vtk", choices=["vtk", "nastran"])
    ap.add_argument("--smooth", type=int, default=10, help="Laplacian iterations")
    args = ap.parse_args()

    cli, data, out = Path(args.cli), Path(args.data), Path(args.out)
    if not cli.exists():
        sys.exit(f"meshkit_cli not found at {cli} (build first: cmake --build build)")
    meshes = sorted(data.glob("*.msh"))
    if not meshes:
        sys.exit(f"no .msh files in {data}")
    out.mkdir(parents=True, exist_ok=True)

    rows, failed = [], 0
    for msh in meshes:
        print(f"\n== {msh.name}")
        try:
            stats, before, after = run_one(cli, msh, out, args.writer, args.smooth)
        except Exception as e:                       # keep going with the other files
            print(f"    FAILED: {e}")
            failed += 1
            continue
        print(f"    {int(stats['points'])} points, {int(stats['triangles'])} triangles")
        print_histograms(before, after)
        rows.append((msh.name, stats))

    print(f"\nSummary (smooth = {args.smooth} iterations, outputs in {out}/)")
    print(f"{'mesh':<22}{'points':>8}{'tris':>8}{'min before':>12}{'min after':>11}{'mean before':>13}{'mean after':>12}")
    for name, s in rows:
        print(f"{name:<22}{int(s['points']):>8}{int(s['triangles']):>8}"
              f"{s['min_angle_before_min']:>12.2f}{s['min_angle_after_min']:>11.2f}"
              f"{s['min_angle_before_mean']:>13.2f}{s['min_angle_after_mean']:>12.2f}")
    if failed:
        print(f"\n{failed} file(s) failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
