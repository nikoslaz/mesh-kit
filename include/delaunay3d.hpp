#ifndef DELAUNAY3D_HPP
#define DELAUNAY3D_HPP

#include "mesh.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh delaunay3d(const std::vector<Point3>& points);
// What one call to insertPoint changed. Used by refinement to update its own data
// instead of rescanning the whole mesh.
struct InsertChangeTet {
    std::vector<Tet> removed;       // copies of the cavity triangles, as they were before
    std::vector<int> added;          // indices of the new triangles
    std::vector<int> addedFrom;      // for each new triangle, the old cavity triangle it replaced
};

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

// Inserts p (Bowyer-Watson). Needs valid neighbours and keeps them valid. If `change` is
// given, it is filled with what changed. A point equal to an existing vertex is added to
// meshPoint but changes no triangle.
int insertPoint3d(TetMesh& mesh, Point3 p, InsertChangeTet* change = nullptr);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

// The triangles whose circumcircle contains p (the cavity p would create), without inserting
// it. Empty if p is outside the mesh.
std::vector<int> findCavityOf(const TetMesh& mesh, const Point3& p);
TetMesh makeSuperTet(const std::vector<Point3>& points);
void removeSuperTet(TetMesh& mesh);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif