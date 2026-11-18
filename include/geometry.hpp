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

#endif