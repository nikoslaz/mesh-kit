#include "delaunay.hpp"
#include "topology.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <utility>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <iostream>
#include <cassert>
#include <functional>

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

void removeSuperTet(TetMesh& mesh) {
    std::vector<Tet> result;

    for (const auto& tet : mesh.tets) {
        if (tet.v0 >= 4 && tet.v1 >= 4 && tet.v2 >= 4 && tet.v3 >= 4) {
            result.push_back(tet);
        }
    }

    for (auto& tet : result) {
        tet.v0 -= 4;
        tet.v1 -= 4;
        tet.v2 -= 4;
        tet.v3 -= 4;
    }

    mesh.points.erase(mesh.points.begin(), mesh.points.begin() + 4);
    mesh.tets = result;
    mesh.lastTet = 0;   // the tet were renumbered
    buildNeighbours(mesh);   // tet indices changed, so the old neighbours are invalid
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

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

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

int insertPoint3d(TetMesh& mesh, Point3 p) {
    int start = walkToTet(mesh, 0, p);
    if (start == -1) return -1;
    int idx = mesh.addPoint(p);

    std::vector<int> cavity = findCavity(mesh, start, p);
    std::vector<CaviteFace> caviteFace = cavityBoundary(mesh,cavity);
    std::vector<int> newTets;

    for (size_t i = 0; i < caviteFace.size(); i++) {
        int a = caviteFace[i].a;
        int b = caviteFace[i].b;
        int c = caviteFace[i].c;
        Tet tet = {a, b, c, idx};

        Point3 pa = mesh.points[a];
        Point3 pb = mesh.points[b];
        Point3 pc = mesh.points[c];

        if (orientation(pa, pb, pc, p) < 0) {
            tet.v1 = c;
            tet.v2 = b;
        } 

        int slot;
        if (i < cavity.size()) { 
            slot = cavity[i]; 
            mesh.tets[slot] = tet; 
        }
        else { 
            mesh.addTet(tet); 
            slot = (int)mesh.tets.size() - 1; 
        }
        newTets.push_back(slot);
    }

    for (size_t i = 0; i < newTets.size(); i++) {
        if (caviteFace[i].outside != -1) {
            Tet& outsideTet = mesh.tets[caviteFace[i].outside];
            int a = caviteFace[i].a;
            int b = caviteFace[i].b;
            int c = caviteFace[i].c;
            setNeighbour(outsideTet, a, b, c, newTets[i]);
            setNeighbour(mesh.tets[newTets[i]], a, b, c, caviteFace[i].outside);
        }
    }

    std::map<std::pair<int,int>, std::vector<int>> edgeTets;
    for (size_t i = 0; i < newTets.size(); i++) {
        std::pair<int,int> key;
        int a = caviteFace[i].a;
        int b = caviteFace[i].b;
        int c = caviteFace[i].c;

        // edge 1: between a and b
        key = { std::min(a, b), std::max(a, b) };   
        edgeTets[key].push_back((newTets[i]));

        // edge 2: between b and c
        key = { std::min(b, c), std::max(b, c) };   
        edgeTets[key].push_back((newTets[i]));

        // edge 3: between c and a
        key = { std::min(c, a), std::max(c, a) };
        edgeTets[key].push_back((newTets[i]));
    }

    for (const auto& e : edgeTets) {
        assert(e.second.size() == 2);
        int x = e.first.first, y = e.first.second;
        int t1 = e.second[0], t2 = e.second[1];
        setNeighbour(mesh.tets[t1], x, y, idx, t2);
        setNeighbour(mesh.tets[t2], x, y, idx, t1);
    }

    // 3D can create fewer tets than the cavity held so we drop the unused cavity slots
    if (cavity.size() > newTets.size()) {
        std::vector<int> spare(cavity.begin() + newTets.size(), cavity.end());
        std::sort(spare.begin(), spare.end(), std::greater<int>());
        for (int s : spare) {
            int last = static_cast<int>(mesh.tets.size()) - 1;
            if (s != last) {
                mesh.tets[s] = mesh.tets[last];
                Tet& moved = mesh.tets[s];
                int* ns[4] = {&moved.n0, &moved.n1, &moved.n2, &moved.n3};
                for (int* n : ns) {
                    if (*n == -1) continue;
                    Tet& o = mesh.tets[*n];
                    int* os[4] = {&o.n0, &o.n1, &o.n2, &o.n3};
                    for (int* q : os) if (*q == last) *q = s;
                }
            }
            mesh.tets.pop_back();
        }
    }

    return idx;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh delaunay3d(const std::vector<Point3>& points) {
    if (points.empty()) return TetMesh{};

    TetMesh mesh = makeSuperTet(points);
    for (const auto& p : points) {
        insertPoint3d(mesh, p);
    }

    removeSuperTet(mesh);

    return mesh;
}

// if (points.empty()) return Mesh{};
//     std::vector<int> order = spatialOrder(points);

//     Mesh mesh = makeSuperTriangle(points);
//     for (int i : order) {
//         insertPoint(mesh, points[i]);
//     }
//     removeSuperTriangle(mesh);

//     // point k of the mesh is points[order[k]]: renumber so the output keeps the input numbering
//     for (auto& t : mesh.meshTriangle) {
//         t.v0 = order[t.v0];
//         t.v1 = order[t.v1];
//         t.v2 = order[t.v2];
//     }
//     mesh.meshPoint = points;
//     return mesh;