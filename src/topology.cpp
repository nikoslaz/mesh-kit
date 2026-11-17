#include "topology.hpp"
#include "geometry.hpp"
#include <algorithm>

std::map<std::pair<int,int>, std::vector<int>> buildEdgeMap(const Mesh& mesh) {
    std::map<std::pair<int,int>, std::vector<int>> edgeMap;

    for (std::size_t i = 0; i < mesh.meshTriangle.size(); i++) {
        const Triangle& tri = mesh.meshTriangle[i];
        std::pair<int,int> key;
        // edge 1: between v0 and v1
        key = { std::min(tri.v0, tri.v1), std::max(tri.v0, tri.v1) };   
        edgeMap[key].push_back(static_cast<int>(i));

        // edge 2: between v1 and v2
        key = { std::min(tri.v1, tri.v2), std::max(tri.v1, tri.v2) };   
        edgeMap[key].push_back(static_cast<int>(i));

        // edge 3: between v2 and v0
        key = { std::min(tri.v2, tri.v0), std::max(tri.v2, tri.v0) };
        edgeMap[key].push_back(static_cast<int>(i));
    }

    return edgeMap;
}

int findTriangleContaining(const Mesh& mesh, const Point& p) {
    for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
        const Triangle& tri = mesh.meshTriangle[ti];
        const Point& a = mesh.meshPoint[tri.v0];
        const Point& b = mesh.meshPoint[tri.v1];
        const Point& c = mesh.meshPoint[tri.v2];

        double o1 = orientation(a, b, p);
        double o2 = orientation(b, c, p);
        double o3 = orientation(c, a, p);

        if (o1 >= 0 && o2 >= 0 && o3 >=0) {
            return static_cast<int>(ti);
        }
        else if (o1 <= 0 && o2 <= 0 && o3 <=0) {
            return static_cast<int>(ti);
        }
    }
    return -1;
}