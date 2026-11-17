#ifndef DELAUNAY
#define DELAUNAY

#include "mesh.hpp"
#include "geometry.hpp"
#include "quality.hpp"

Mesh delaunay(const std::vector<Point>& points);
int insertPoint(Mesh& mesh, Point p);
Mesh makeSuperTriangle(const std::vector<Point>& points);
void removeSuperTriangle(Mesh& mesh);

#endif