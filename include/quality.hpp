#ifndef QUALITY
#define QUALITY

#include "mesh.hpp"
#include "geometry.hpp"

double angleAt(const Point& corner, const Point& p, const Point& q);
double minAngle(const Point& a, const Point& b, const Point& c);
std::vector<double> minAngles(const Mesh& mesh);
double radiusEdgeRatio(const Point& a, const Point& b, const Point& c);

#endif