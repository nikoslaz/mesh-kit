#include "delaunay.hpp"
#include "topology.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <utility>
#include <cstdlib>
#include <iostream>

struct CaviteFace {
    int a;
    int b;
    int c;
    int outside;
    int owner;     // the cavity tetra this face belongs to
};

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh makeSuperTet(const std::vector<Point3>& points) {
    TetMesh mesh;
    if (points.empty()) return mesh;

    double maxX = points[0].x;
    double minX = points[0].x;
    double maxY = points[0].y;
    double minY = points[0].y;
    double maxZ = points[0].z;
    double minZ = points[0].z;
    double width, height, m, depth;
    double cx, cy, cz;

    for (const auto& point: points) {
        maxX = std::max(maxX, point.x);
        minX = std::min(minX, point.x);
        maxY = std::max(maxY, point.y);
        minY = std::min(minY, point.y);
        maxZ = std::max(maxZ, point.z);
        minZ = std::min(minZ, point.z);
    }

    width = maxX - minX;
    height = maxY - minY;
    depth = maxZ - minZ;

    m = std::max({width, height, depth});
    if (m == 0) m = 1;
    cx = (minX + maxX) / 2;
    cy = (minY + maxY) / 2;
    cz = (minZ + maxZ) / 2;

    const double K = 1000;  
    double S  = 3 * (K + 1) * m;                       // size of the tetrahedron
    double ox = cx - K * m;                            // the corner "o"
    double oy = cy - K * m;
    double oz = cz - K * m;

    Point3 a{ox, oy, oz};
    Point3 b{ox + S, oy, oz};
    Point3 c{ox, oy + S, oz };
    Point3 d{ox, oy, oz + S};                  

    Tet t;
    t.v0 =  static_cast<int>(mesh.addPoint(a));
    t.v1 =  static_cast<int>(mesh.addPoint(c)); // changed c with b 
    t.v2 =  static_cast<int>(mesh.addPoint(b)); // for positive orientation
    t.v3 =  static_cast<int>(mesh.addPoint(d));

    mesh.addTet(t);

    return mesh;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //


std::vector<int> findCavityOf(const TetMesh& mesh, const Point3& p) {

}

static std::vector<int> findCavity(const TetMesh& mesh, int start, const Point3& p) {
    std::vector<int> cavity;
    std::vector<bool> visited(mesh.tets.size(), false);
    std::queue<int> q;

    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        int t = q.front();
        q.pop();
        cavity.push_back(t);

        const Tet& tet = mesh.tets[t];
        int nbrs[4] = {tet.n0, tet.n1, tet.n2, tet.n3};
        for (int k : nbrs) {
            if (k == -1 || visited[k]) continue;
            const Tet& nb   = mesh.tets[k];
            const Point3& a = mesh.points[nb.v0];
            const Point3& b = mesh.points[nb.v2];
            const Point3& c = mesh.points[nb.v1];
            const Point3& d = mesh.points[nb.v3];
            if (inCircumcircle(a, b, c, d, p)) {
                q.push(k);
                visited[k] = true;
            }
        }
    }
    return cavity;
}

static std::vector<CaviteFace> cavityBoundary(const TetMesh& mesh, const std::vector<int>& cavity) {
    std::vector<CaviteFace> caviteFace;
    std::vector<bool> inCavity(mesh.tets.size(), false);

    for (const auto& c : cavity) {
        inCavity[c] = true;
    }
    
    for (const auto& c : cavity) {
        const Tet& t = mesh.tets[c];

        if (t.n0 == -1 || !inCavity[t.n0]) {
            CaviteFace ce = {t.v1, t.v2, t.v3, t.n0, c};
            caviteFace.push_back(ce);
        }
        if (t.n1 == -1 || !inCavity[t.n1]) {
            CaviteFace ce = {t.v0, t.v2, t.v3, t.n1, c};
            caviteFace.push_back(ce);
        }
        if (t.n2 == -1 || !inCavity[t.n2]) {
            CaviteFace ce = {t.v0, t.v1, t.v3, t.n2, c};
            caviteFace.push_back(ce);
        }
        if (t.n3 == -1 || !inCavity[t.n3]) {
            CaviteFace ce = {t.v0, t.v1, t.v2, t.n3, c};
            caviteFace.push_back(ce);
        }
    }

    return caviteFace;
}