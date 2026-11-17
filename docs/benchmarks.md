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
