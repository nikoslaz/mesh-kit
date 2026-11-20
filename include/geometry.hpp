#ifndef GEOMETRY
#define GEOMETRY

#include "mesh.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

double orientation(const Point& a, const Point& b, const Point& c);                         // 2D
double orientation(const Point3& a, const Point3& b, const Point3& c, const Point3& d);     // 3D

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

bool inCircumcircle(const Point& a, const Point& b, const Point& c, const Point& d); // 2D
bool inCircumcircle(const Point3& a, const Point3& b, const Point3& c, const Point3& d, const Point3& e); // 3D

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

bool encroaches(const Point& a, const Point& b, const Point& p);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

Point circumcenter(const Point& a, const Point& b, const Point& c); 

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

bool flip23Tets(const TetMesh& mesh, int t1, int t2, Tet out[3]);
bool flip32Tets(const TetMesh& mesh, int t1, int t2, int t3, Tet out[2]);
// Flips applied to the mesh. Need valid neighbours and keep them valid (updated locally).
// apply23 adds one tet at the end. apply32 removes t3 by moving the last tet into slot t3.
bool apply23(TetMesh& mesh, int t1, int t2);
bool apply32(TetMesh& mesh, int t1, int t2, int t3);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif