#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include <map>
#include <utility>
#include <vector>
#include "mesh.hpp"
#include "domain.hpp"

std::map<std::pair<int,int>, std::vector<int>> buildEdgeMap(const Mesh& mesh);
int findTriangleContaining(const Mesh& mesh, const Point& p);
std::vector<bool> insideTriangles(const Mesh& mesh, const std::vector<Segment>& segments, const std::vector<Point>& holes);
Mesh extractInside(const Mesh& mesh, const std::vector<bool>& inside);

#endif