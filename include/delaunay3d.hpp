#ifndef DELAUNAY3D_HPP
#define DELAUNAY3D_HPP

#include "mesh.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh delaunay3d(const std::vector<Point3>& points);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

// Inserts p (Bowyer-Watson 3D). Needs valid neighbours and keeps them valid.
// Returns the index of p, or -1 if p is outside the mesh.
int insertPoint3d(TetMesh& mesh, Point3 p);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh makeSuperTet(const std::vector<Point3>& points);
void removeSuperTet(TetMesh& mesh);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif