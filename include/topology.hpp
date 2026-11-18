#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include <map>
#include <utility>
#include <vector>
#include "mesh.hpp"
#include "domain.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

std::map<std::pair<int,int>, std::vector<int>> buildEdgeMap(const Mesh& mesh); // 2D
std::map<std::array<int,3>, std::vector<int>> buildFaceMap(const TetMesh& mesh); // 3D

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

int findTriangleContaining(const Mesh& mesh, const Point& p);
std::vector<bool> insideTriangles(const Mesh& mesh, const std::vector<Segment>& segments, const std::vector<Point>& holes);
Mesh extractInside(const Mesh& mesh, const std::vector<bool>& inside);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

void buildNeighbours(Mesh& mesh); // 2D
void buildNeighbours(TetMesh& tet); // 3D

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

int walkToTriangle(const Mesh& mesh, int start, const Point& p); //2D
int walkToTet(const TetMesh& mesh, int start, const Point3& p); // 3D

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

// Set the neighbour slot of tri that faces edge (a, b) to `other`.
void setNeighbour(Triangle& tri, int a, int b, int other); // 2D
void setNeighbour(Tet& tet, int a, int b, int c,int other); // 3D

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif