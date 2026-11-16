#!/usr/bin/env python3
"""Generate sample Gmsh 2.2 meshes for the batch script.

Each sample is a unit-square grid split into triangles, with the interior
nodes moved randomly ("jitter") so the mesh has bad triangles to fix.

    python3 scripts/make_samples.py            # writes tests/data/grid_*.msh
"""
import argparse
import random
from pathlib import Path


def grid_mesh(n, jitter, seed):
    """n x n cells on [0,1]^2; interior nodes moved by up to jitter * cell size."""
    rng = random.Random(seed)
    h = 1.0 / n
    nodes = []
    for j in range(n + 1):
        for i in range(n + 1):
            x, y = i * h, j * h
            if 0 < i < n and 0 < j < n:
                x += rng.uniform(-jitter, jitter) * h
                y += rng.uniform(-jitter, jitter) * h
            nodes.append((x, y))

    def nid(i, j):                       # 1-based Gmsh node id
        return j * (n + 1) + i + 1

    tris = []
    for j in range(n):
        for i in range(n):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))
    return nodes, tris


def write_msh(path, nodes, tris):
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n{len(nodes)}\n")
        for k, (x, y) in enumerate(nodes, start=1):
            f.write(f"{k} {x:.10g} {y:.10g} 0\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n{len(tris)}\n")
        for k, (a, b, c) in enumerate(tris, start=1):
            f.write(f"{k} 2 2 0 1 {a} {b} {c}\n")    # type 2 = triangle, 2 tags
        f.write("$EndElements\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="tests/data", help="output folder")
    args = ap.parse_args()

    samples = [  # (name, cells per side, jitter)
        ("grid_10_regular", 10, 0.0),
        ("grid_10_jitter", 10, 0.3),
        ("grid_25_jitter", 25, 0.35),
        ("grid_50_jitter", 50, 0.35),
    ]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for seed, (name, n, jitter) in enumerate(samples, start=1):
        nodes, tris = grid_mesh(n, jitter, seed)
        path = out / f"{name}.msh"
        write_msh(path, nodes, tris)
        print(f"wrote {path}  ({len(nodes)} nodes, {len(tris)} triangles)")


if __name__ == "__main__":
    main()
