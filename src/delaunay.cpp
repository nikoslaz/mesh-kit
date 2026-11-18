#include "delaunay.hpp"
#include "topology.hpp"
#include <algorithm>
#include <map>
#include <queue>
#include <utility>
#include <cstdlib>
#include <iostream>

// One edge of the cavity outline, going a -> b with the cavity on its left,
// and the triangle on the other side (-1 if none).
struct CavityEdge {
    int a;
    int b;
    int outside;
};

Mesh makeSuperTriangle(const std::vector<Point>& points) {
    Mesh mesh;
    if (points.empty()) return mesh;

    double maxX = points[0].x;
    double minX = points[0].x;
    double maxY = points[0].y;
    double minY = points[0].y;
    double width, height, m;
    double cx, cy;

    for (const auto& point: points) {
        maxX = std::max(maxX, point.x);
        minX = std::min(minX, point.x);
        maxY = std::max(maxY, point.y);
        minY = std::min(minY, point.y);
    }

    width = maxX - minX;
    height = maxY - minY;
    m = std::max(width, height);
    if (m == 0) m = 1;
    cx = (minX + maxX) / 2;
    cy = (minY + maxY) / 2;

    Point a{cx - 1000*m, cy - m};
    Point b{cx + 1000*m, cy - m};
    Point c{cx, cy + 1000*m};

    Triangle t;
    t.v0 =  static_cast<int>(mesh.addPoint(a));
    t.v1 =  static_cast<int>(mesh.addPoint(b));
    t.v2 =  static_cast<int>(mesh.addPoint(c));

    mesh.addTriangle(t);

    return mesh;
}

static std::vector<int> findCavity(const Mesh& mesh, int start, const Point& p) {
    std::vector<int> cavity;
    std::vector<bool> visited(mesh.meshTriangle.size(), false);
    std::queue<int> q;

    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        int t = q.front();
        q.pop();
        cavity.push_back(t);

        const Triangle& tri = mesh.meshTriangle[t];
        int nbrs[3] = {tri.n0, tri.n1, tri.n2};
        for (int k : nbrs) {
            if (k == -1 || visited[k]) continue;
            const Triangle& nb = mesh.meshTriangle[k];
            const Point& a = mesh.meshPoint[nb.v0];
            const Point& b = mesh.meshPoint[nb.v1];
            const Point& c = mesh.meshPoint[nb.v2];
            if (inCircumcircle(a, b, c, p)) {
                q.push(k);
                visited[k] = true;
            }
        }
    }
    return cavity;
}

static std::vector<CavityEdge> cavityBoundary(const Mesh& mesh, const std::vector<int>& cavity) {
    std::vector<CavityEdge> cavityEdge;
    std::vector<bool> inCavity(mesh.meshTriangle.size(), false);

    for (const auto& c : cavity) {
        inCavity[c] = true;
    }
    
    for (const auto& c : cavity) {
        const Triangle& t = mesh.meshTriangle[c];

        if (t.n0 == -1 || !inCavity[t.n0]) {
            CavityEdge ce = {t.v1, t.v2, t.n0};
            cavityEdge.push_back(ce);
        }
        if (t.n1 == -1 || !inCavity[t.n1]) {
            CavityEdge ce = {t.v2, t.v0, t.n1};
            cavityEdge.push_back(ce);
        }
        if (t.n2 == -1 || !inCavity[t.n2]) {
            CavityEdge ce = {t.v0, t.v1, t.n2};
            cavityEdge.push_back(ce);
        }
    }

    return cavityEdge;
}

// Bowyer-Watson insertion using neighbours:
//   1. walk to the triangle containing p
//   2. BFS over neighbours to find the cavity (triangles whose circumcircle contains p)
//   3. replace the cavity by a fan of triangles {p, a, b}, one per outline edge,
//      reusing the cavity's slots so no other triangle index changes
// Requires valid neighbours on entry and keeps them valid. All triangles are CCW.
int insertPoint(Mesh& mesh, Point p) {
    int start = mesh.lastTriangle;
    if (start < 0 || start >= static_cast<int>(mesh.meshTriangle.size())) start = 0;
    int tIdx = walkToTriangle(mesh, start, p);
    if (tIdx == -1) {
        std::cerr << "insertPoint: walk did not find a triangle for (" << p.x << ", " << p.y << ")\n";
        std::abort();
    }
    int idx = static_cast<int>(mesh.addPoint(p));

    std::vector<int> cavity = findCavity(mesh, tIdx, p);
    std::vector<CavityEdge> cavEdge = cavityBoundary(mesh, cavity);

    // the new triangles: vertices and the neighbour across the outline edge
    std::vector<Triangle> newTri(cavEdge.size());
    for (std::size_t i = 0; i < newTri.size(); i++) {
        newTri[i].v0 = idx;
        newTri[i].v1 = cavEdge[i].a;
        newTri[i].v2 = cavEdge[i].b;
        newTri[i].n0 = cavEdge[i].outside;
    }

    // where each new triangle goes: reuse the cavity's slots, append the 2 extra ones
    std::vector<int> slot(newTri.size());
    for (std::size_t i = 0; i < newTri.size(); i++) {
        if (i < cavity.size()) {
            slot[i] = cavity[i];
        }
        else {
            mesh.meshTriangle.push_back(Triangle{});
            slot[i] = static_cast<int>(mesh.meshTriangle.size()) - 1;
        }
    }

    // link the fan: n1 = next triangle around p, n2 = previous one
    std::map<int,int> startAt, endAt;
    for (std::size_t i = 0; i < cavEdge.size(); i++) {
        startAt.insert({cavEdge[i].a, static_cast<int>(i)});
        endAt.insert({cavEdge[i].b, static_cast<int>(i)});
    }
    for (std::size_t i = 0; i < cavEdge.size(); i++) {
        int j = startAt.at(cavEdge[i].b);
        int k = endAt.at(cavEdge[i].a);
        newTri[i].n1 = slot[j];
        newTri[i].n2 = slot[k];
    }

    // write them into the mesh and point the outside triangles at them
    for (std::size_t i = 0; i < newTri.size(); i++) {
        mesh.meshTriangle[slot[i]] = newTri[i];
        if (cavEdge[i].outside == -1) {continue;}
        setNeighbour(mesh.meshTriangle[cavEdge[i].outside], cavEdge[i].a, cavEdge[i].b, slot[i]);
    }

    mesh.lastTriangle = slot[0];   // a triangle touching p: a good start for the next walk
    return idx;
}

void removeSuperTriangle(Mesh& mesh) { 
    std::vector<Triangle> result;

    for (const auto& tri : mesh.meshTriangle) {
        if (tri.v0 >= 3 && tri.v1 >= 3 && tri.v2 >= 3) {
            result.push_back(tri);
        }
    }

    for (auto& tri : result) {
        tri.v0 -= 3;
        tri.v1 -= 3;
        tri.v2 -= 3;
    }

    mesh.meshPoint.erase(mesh.meshPoint.begin(), mesh.meshPoint.begin() + 3);
    mesh.meshTriangle = result;
    mesh.lastTriangle = 0;   // the triangles were renumbered
    buildNeighbours(mesh);   // triangle indices changed, so the old neighbours are invalid
}

Mesh delaunay(const std::vector<Point>& points) {
    if (points.empty()) return Mesh{};
    Mesh mesh = makeSuperTriangle(points);
    for (const auto& point : points) {
        insertPoint(mesh, point);
    }
    removeSuperTriangle(mesh);
    return mesh;
}