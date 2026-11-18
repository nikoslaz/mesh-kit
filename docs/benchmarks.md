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
