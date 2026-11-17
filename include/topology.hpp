#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include <map>
#include <utility>
#include <vector>
#include "mesh.hpp"

std::map<std::pair<int,int>, std::vector<int>> buildEdgeMap(const Mesh& mesh);
int findTriangleContaining(const Mesh& mesh, const Point& p);

#endif