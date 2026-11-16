#ifndef DELAUNAY
#define DELAUNAY

#include "mesh.hpp"
#include "geometry.hpp"
#include "quality.hpp"

Mesh delaunay(const std::vector<Point>& points);

#endif