#ifndef GEOMETRY
#define GEOMETRY

#include "mesh.hpp"

double orientation(const Point& a, const Point& b, const Point& c);
bool inCircumcircle(const Point& a, const Point& b, const Point& c, const Point& d);

#endif