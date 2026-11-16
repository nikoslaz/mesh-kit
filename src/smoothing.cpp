#include "smoothing.hpp"

std::vector<bool> boundaryVertices(const Mesh& mesh) {
    std::map<std::pair<int,int>, int> edgeCount;
    std::vector<bool> result(mesh.meshPoint.size(), false);

    for (const auto& tri : mesh.meshTriangle) {
        std::pair<int,int> key;
        // edge 1: between v0 and v1
        key = { std::min(tri.v0, tri.v1), std::max(tri.v0, tri.v1) };   
        edgeCount[key]++;

        // edge 2: between v1 and v2
        key = { std::min(tri.v1, tri.v2), std::max(tri.v1, tri.v2) };   
        edgeCount[key]++;

        // edge 3: between v2 and v0
        key = { std::min(tri.v2, tri.v0), std::max(tri.v2, tri.v0) };
        edgeCount[key]++;
    }

    for (const auto& entry : edgeCount) {
        if (entry.second == 1) {
            result[entry.first.first]  = true;
            result[entry.first.second] = true;
        }
    }

    return result;
}

std::vector<std::set<int>> buildNeighbours(const Mesh& mesh) { 
    std::vector<std::set<int>> neighbours(mesh.meshPoint.size());
    for (const auto& tri : mesh.meshTriangle) {
        neighbours[tri.v0].insert(tri.v1);
        neighbours[tri.v0].insert(tri.v2);
        neighbours[tri.v1].insert(tri.v0);
        neighbours[tri.v1].insert(tri.v2);
        neighbours[tri.v2].insert(tri.v0);
        neighbours[tri.v2].insert(tri.v1);
    }

    return neighbours;
}

void laplacianSmooth(Mesh& mesh, int iterations) {
    auto boundary   = boundaryVertices(mesh);
    auto neighbours = buildNeighbours(mesh);

    for (int iter = 0; iter < iterations; iter++) {
        std::vector<Point> next = mesh.meshPoint;
        
        for (std::size_t i = 0; i < next.size(); i++) {
            if (neighbours[i].empty()) continue;
            if (!boundary[i]) {
                double sumX = 0;
                double sumY = 0;
                for (const auto& neigh : neighbours[i]) {
                    sumX += mesh.meshPoint[neigh].x;
                    sumY += mesh.meshPoint[neigh].y;
                }
                next[i].x = sumX / neighbours[i].size();
                next[i].y = sumY / neighbours[i].size();
            }
        }

        mesh.meshPoint = next;
    }
}

