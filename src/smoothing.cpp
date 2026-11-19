#include "smoothing.hpp"
#include "geometry.hpp"
#include "quality.hpp"

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

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

std::vector<std::vector<int>> vertexTets(const TetMesh& mesh) { 
    std::vector<std::vector<int>> result(mesh.points.size());

    for (std::size_t i = 0; i < mesh.tets.size(); i++) {
        result[mesh.tets[i].v0].push_back(static_cast<int>(i));
        result[mesh.tets[i].v1].push_back(static_cast<int>(i));
        result[mesh.tets[i].v2].push_back(static_cast<int>(i));
        result[mesh.tets[i].v3].push_back(static_cast<int>(i));
    }

    return result;
}

std::vector<std::set<int>> vertexNeighbours(const TetMesh& mesh, const std::vector<std::vector<int>>& vt) {
    std::vector<std::set<int>> result(vt.size());

    for (std::size_t v = 0; v < vt.size(); ++v) {
        for (const auto& idx : vt[v]) {
            const Tet& tet = mesh.tets[idx];
            int c[4] = {tet.v0, tet.v1, tet.v2, tet.v3};
            for (std::size_t i = 0; i < 4; i ++) {
                if (c[i] != static_cast<int>(v)) {
                    result[v].insert(c[i]);
                }
            }
        }
    }

    return result;
}

std::vector<bool> hullVertices(const TetMesh& mesh) {
    std::vector<bool> result(mesh.points.size(), false);

    for (std::size_t i = 0; i < mesh.tets.size(); i++) {
        const Tet& tet = mesh.tets[i];

        int v[4] = {tet.v0, tet.v1, tet.v2, tet.v3};
        int n[4] = {tet.n0, tet.n1, tet.n2, tet.n3};

        for(std::size_t j = 0; j < 4; j++) {
            if (n[j] == -1) {
                for (std::size_t k = 0; k < 4; k++) {
                    if (k != j) result[v[k]] = true;
                }
            }
        }
    }

    return result;
}

bool tryMoveVertex(TetMesh& mesh, int v, const Point3& target, const std::vector<int>& tetsAround) {
    auto quality = [&mesh](int id) {
        const Tet& t = mesh.tets[id];
        return tetRadiusRatio(mesh.points[t.v0], mesh.points[t.v1], mesh.points[t.v2], mesh.points[t.v3]);
    };

    double worstBefore = 1.0;
    for (int id : tetsAround) worstBefore = std::min(worstBefore, quality(id));

    const Point3 oldPos = mesh.points[v];
    mesh.points[v] = target;

    bool ok = true;
    double worstAfter = 1.0;
    for (int id : tetsAround) {
        const Tet& t = mesh.tets[id];
        if (orientation(mesh.points[t.v0], mesh.points[t.v1], mesh.points[t.v2], mesh.points[t.v3]) <= 0) {
            ok = false;
            break;
        }
        worstAfter = std::min(worstAfter, quality(id));
    }

    if (!ok || worstAfter < worstBefore) {
        mesh.points[v] = oldPos;
        return false;
    }
    return true;
}
