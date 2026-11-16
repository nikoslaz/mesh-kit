#ifndef SMOOTH
#define SMOOTH

#include "mesh.hpp"
#include <algorithm>
#include <map>
#include <utility>
#include <set>

std::vector<bool> boundaryVertices(const Mesh& mesh);
std::vector<std::set<int>> buildNeighbours(const Mesh& mesh);
void laplacianSmooth(Mesh& mesh, int iterations);

#endif