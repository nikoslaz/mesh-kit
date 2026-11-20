#ifndef DELAUNAY3D_HPP
#define DELAUNAY3D_HPP

#include "mesh.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh delaunay3d(const std::vector<Point3>& points);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

// Inserts p (Bowyer-Watson 3D) into a closed mesh with ghost tets (point 0 = the vertex at infinity),
// as built inside delaunay3d. Needs valid neighbours and keeps them valid. Returns the index of p.
int insertPoint3d(TetMesh& mesh, Point3 p);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif