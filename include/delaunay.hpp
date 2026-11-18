#ifndef DELAUNAY
#define DELAUNAY

#include "mesh.hpp"
#include "geometry.hpp"
#include "quality.hpp"

Mesh delaunay(const std::vector<Point>& points);
// What one call to insertPoint changed. Used by refinement to update its own data
// instead of rescanning the whole mesh.
struct InsertChange {
    std::vector<Triangle> removed;   // copies of the cavity triangles, as they were before
    std::vector<int> added;          // indices of the new triangles
    std::vector<int> addedFrom;      // for each new triangle, the old cavity triangle it replaced
};

// Inserts p (Bowyer-Watson). Needs valid neighbours and keeps them valid. If `change` is
// given, it is filled with what changed. A point equal to an existing vertex is added to
// meshPoint but changes no triangle.
int insertPoint(Mesh& mesh, Point p, InsertChange* change = nullptr);

// The triangles whose circumcircle contains p (the cavity p would create), without inserting
// it. Empty if p is outside the mesh.
std::vector<int> findCavityOf(const Mesh& mesh, const Point& p);
Mesh makeSuperTriangle(const std::vector<Point>& points);
void removeSuperTriangle(Mesh& mesh);

#endif