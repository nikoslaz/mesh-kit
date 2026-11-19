#ifndef SMOOTH
#define SMOOTH

#include "mesh.hpp"
#include <algorithm>
#include <map>
#include <utility>
#include <set>
#include <vector>

std::vector<bool> boundaryVertices(const Mesh& mesh);
std::vector<std::set<int>> buildNeighbours(const Mesh& mesh);
void laplacianSmooth(Mesh& mesh, int iterations);
bool tryMoveVertex(TetMesh& mesh, int v, const Point3& target, const std::vector<int>& tetsAround);

// 3D guarded Laplacian smoothing. Hull vertices stay fixed. Returns the number of moves that were kept.
int laplacianSmooth3d(TetMesh& mesh, int iterations);
std::vector<std::vector<int>> vertexTets(const TetMesh& mesh);
std::vector<std::set<int>> vertexNeighbours(const TetMesh& mesh, const std::vector<std::vector<int>>& vt);
std::vector<bool> hullVertices(const TetMesh& mesh); // true if v lies on the boundary
int laplacianSmooth3d(TetMesh& mesh, int iterations);

#endif