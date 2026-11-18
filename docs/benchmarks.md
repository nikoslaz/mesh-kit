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

`meshDomain` is not much faster yet: its time is now spent in the refinement loop, which
recomputes the interior and scans all triangles and points on every iteration.
