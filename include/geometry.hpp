#ifndef GEOMETRY
#define GEOMETRY

#include "mesh.hpp"

double orientation(const Point& a, const Point& b, const Point& c);
bool inCircumcircle(const Point& a, const Point& b, const Point& c, const Point& d);
bool encroaches(const Point& a, const Point& b, const Point& p);
Point circumcenter(const Point& a, const Point& b, const Point& c);

#endif