# Benchmarks

`meshkit_bench`, built with `-O2`. Times are wall-clock seconds. Absolute numbers depend on the
machine; compare before and after on the same one.

## Baseline (before optimisation)

Each insertion scans every triangle, so cost grows quadratically: doubling the points
multiplies the time by about 4.

```
Delaunay triangulation of random points in the unit square
  points   seconds   triangles
     500    0.0071         983
    1000    0.0267        1983
    2000    0.1066        3978
    4000    0.4344        7976

meshDomain (max radius-edge ratio 1.2) on a 10x10 plate with a circular hole
  hole sides   seconds    points   triangles
          16    0.0008        52          72
          32    0.0033       111         172
          64    0.0107       204         324
         128    0.0351       364         579
         256    0.1253       665        1055
```

## After: neighbour-based insertion

Every triangle now stores its three neighbours. `insertPoint` finds the triangle containing the
new point by walking across neighbours, grows the cavity by a breadth-first search over
neighbours, and replaces the cavity in place. It no longer scans every triangle.

The table below was measured with old and new code on the *same* machine (a Linux sandbox,
`-O2`, random points in the unit square), so the two columns can be compared with each other.
The baseline above was measured on a different machine, so do not compare it with these numbers.

```
Delaunay triangulation of random points in the unit square
  points   scan all triangles   neighbours   speed-up
    1000             0.0125        0.0022        5.6x
    2000             0.0517        0.0047         11x
    4000             0.2130        0.0104         21x
    8000             0.8932        0.0255         35x
   16000             3.7112        0.0638         58x
   32000            22.424         0.1661        135x
   64000                 -         0.4561           -
```

The old time grows by about 4x when the points double. The new time grows by about 2.5x, and
the output is the same triangulation (checked on random sets and exact grids).

## After: spatial insertion order and incremental refinement

Two more changes:

- The walk now starts from the last inserted triangle, and `delaunay()` inserts the points in
  a "snake" order over rows of cells (then renumbers them, so the output keeps the input
  numbering). Each walk is now only a few steps long.
- `refineMesh` no longer rescans the mesh on every step. New triangles go on a work list of
  possibly bad triangles, inherit the inside/outside flag of the triangle they replaced, and
  segments are checked for encroachment only against the vertices opposite them. Each of
  these shortcuts was checked against the full recomputation at every step while developing.

Same machine as the table above (Linux sandbox, `-O2`):

```
Delaunay, random points           before (neighbours only)   now
    16000                                   0.061 s        0.032 s
    64000                                   0.419 s        0.148 s
   256000                                   3.602 s        0.683 s

meshDomain, ratio 1.2, 10x10 plate with an n-sided hole
  hole sides         before            now     points now
          64       0.0063 s       0.0008 s           196
         128       0.0229 s       0.0013 s           339
         256       0.0853 s       0.0028 s           684
         512       0.3699 s       0.0060 s          1331
        1024       1.6932 s       0.0115 s          2549
        2048       7.9142 s       0.0262 s          5053
```

At 2048 sides refinement is about 300 times faster, and doubling the hole now roughly doubles
the time instead of multiplying it by 4.5. The meshes are not identical to before: bad
triangles are fixed in a different order, so the point count differs by a few percent. They
meet the same guarantees (every triangle within the ratio limit, the same area and topology).

## 3D Delaunay tetrahedralization

`delaunay3d`, random points in the unit cube, `-O2`.

```
  points   tets     seconds (random order, start at tet 0)   seconds (3D snake order, start at last tet)
   10000   66380    0.14                                     0.16
   50000   335372   0.96                                     0.94
  100000   671683   2.29                                     1.97
```

The walk start and the ordering help only a little: at 100,000 points about 1.4 s is insertion and
0.4 s is rebuilding the neighbours after the super tetrahedron is removed. Removing the
per-insertion allocations changed nothing, so the cost is spread over the predicates and the
cavity bookkeeping.

## 3D sliver removal

`meshkit_cli tets` on the cube corners plus random points (seed 42), threshold 0.1, `-O2`, the
same Linux sandbox as above. `improveSlivers` runs up to 20 rounds of flips and smoothing.

```
  points   tets before   slivers before   tets after   slivers after   flips   seconds
     300          1915              323         1853             240      76     0.03
    1000          6556              803         6372             544     204     0.13
    3000         19875             1721        19419            1139     506     0.55
   10000         66821             4181        65518            2494    1433     3.56
```

Worst radius ratio does not change (1.6e-6 at 300 points): the worst tets are flat tets against
the convex hull, and the hull is never changed. Slivers that touch a hull vertex barely move
(491 -> 473 at 10,000 points); the other slivers drop by about 45-50% (3690 -> 2021).

### Speed: local neighbour updates

The first version rebuilt the neighbours of the whole mesh after every flip (40 ms at 66,000 tets)
and made one flip per pass, so the time grew about as the square of the number of flips. Now
`apply23` / `apply32` re-link only the faces of the new tets and their outside neighbours, and
`apply32` removes a tet by moving the last tet into its slot instead of shifting the array. A pass
tries every sliver once, skipping tets that an earlier flip already replaced. The local update was
checked against a full `buildNeighbours` after each of 666 random flips (identical).

```
  points   one flip per pass + full rebuild   now     speed-up
     300                          0.098 s     0.029 s     3x
    2000                          3.18 s      0.31 s     10x
   10000                         75.1 s       3.58 s     21x
```

What is left is rebuilding the vertex-to-tets lists after each flip (about 1.7 ms at 66,000 tets).

### Delaunay property

Before removal no tet has a point inside its circumsphere. After removal they do: 693 of 1853 tets
at 300 points, 5739 of 19419 at 3000. The flips choose tets by quality, not by the
empty-circumsphere rule, so the result is no longer a Delaunay triangulation. That is a
trade-off, not a bug, but it must be kept in mind when using these meshes.
