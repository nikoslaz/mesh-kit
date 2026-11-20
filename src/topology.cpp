#include "topology.hpp"
#include "geometry.hpp"
#include "quality.hpp"
#include "delaunay.hpp"
#include <algorithm>
#include <array>
#include <set>
#include <utility>
#include <random>

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

std::map<std::array<int,3>, std::vector<int>> buildFaceMap(const TetMesh& mesh) {
    std::map<std::array<int,3>, std::vector<int>> faceMap;

    for (std::size_t i = 0; i < mesh.tets.size(); i++) {
        const Tet& tet = mesh.tets[i];
        std::array<int,3> key;
        // face 1: between v1 and v2 and v3
        key = { tet.v1, tet.v2, tet.v3 };
        std::sort(key.begin(), key.end());
        faceMap[key].push_back(static_cast<int>(i));

        // face 2: between v2 and v3 and v0
        key = { tet.v2, tet.v3, tet.v0 };
        std::sort(key.begin(), key.end());
        faceMap[key].push_back(static_cast<int>(i));

        // face 3: between v3 and v0 and v1
        key = { tet.v3, tet.v0, tet.v1 };
        std::sort(key.begin(), key.end());
        faceMap[key].push_back(static_cast<int>(i));

        // face 4: between v0 and v1 and v2
        key = { tet.v0, tet.v1, tet.v2 };
        std::sort(key.begin(), key.end());
        faceMap[key].push_back(static_cast<int>(i));
    }

    return faceMap;    
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //


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


std::vector<bool> insideTriangles(const Mesh& mesh, const std::vector<Segment>& segments, const std::vector<Point>& holes) {
    std::map<std::pair<int,int>, std::vector<int>> edgeMap = buildEdgeMap(mesh);

    std::set<std::pair<int,int>> walls;
    for (std::size_t ti = 0; ti < segments.size(); ti++) {
        const Segment& seg = segments[ti];
        std::pair<int,int> w;
        w = {std::min(seg.a,seg.b), std::max(seg.a,seg.b)};   
        walls.insert(w);
    }

    std::vector<bool> painted(mesh.meshTriangle.size(), false);
    std::vector<int> stack;

    for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
        const Triangle& t = mesh.meshTriangle[ti];
        bool touchesCorner = t.v0 < 3 || t.v1 < 3 || t.v2 < 3;
        if (touchesCorner) {
            painted[ti] = true;
            stack.push_back(static_cast<int>(ti));
        }
    }

    for (const auto& h : holes) {
        int t = findTriangleContaining(mesh, h);
        if (t == -1) {continue;}
        painted[t] = true;
        stack.push_back(t);
    }

    while (!stack.empty()) {
        int t = stack.back();
        stack.pop_back();
        Triangle tri = mesh.meshTriangle[t];
        std::pair<int,int> key;

        // edge 1: between v0 and v1
        key = { std::min(tri.v0, tri.v1), std::max(tri.v0, tri.v1) };   
        if (walls.count(key) == 0) {
            for (const auto& n : edgeMap[key]) {
                if (!painted[n]) {
                    painted[n] = true;
                    stack.push_back(n);
                }
            }
        }

        // edge 2: between v1 and v2
        key = { std::min(tri.v1, tri.v2), std::max(tri.v1, tri.v2) };   
        if (walls.count(key) == 0) {
            for (const auto& n : edgeMap[key]) {
                if (!painted[n]) {
                    painted[n] = true;
                    stack.push_back(n);
                }
            }
        }

        // edge 3: between v2 and v0
        key = { std::min(tri.v2, tri.v0), std::max(tri.v2, tri.v0) };
        if (walls.count(key) == 0) {
            for (const auto& n : edgeMap[key]) {
                if (!painted[n]) {
                    painted[n] = true;
                    stack.push_back(n);
                }
            }
        }
    }

    std::vector<bool> inside(mesh.meshTriangle.size());
    for (std::size_t t = 0; t < inside.size(); t++) inside[t] = !painted[t];
    return inside;
}

Mesh extractInside(const Mesh& mesh, const std::vector<bool>& inside) {
    std::vector<int> newIndex(mesh.meshPoint.size(), -1);
    Mesh result;

    for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
        if (!inside[ti]) {continue;}
        Triangle tri = mesh.meshTriangle[ti];
        int v[3] = {tri.v0, tri.v1, tri.v2};

        for (int k = 0; k < 3; k++) {
            if (newIndex[v[k]] == -1) {
                newIndex[v[k]] = (int) result.meshPoint.size();
                result.addPoint(mesh.meshPoint[v[k]]);
            }
            v[k] = newIndex[v[k]]; 
        }
        result.addTriangle({v[0], v[1], v[2]});
    }

    return result;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

void setNeighbour(Triangle& tri, int a, int b, int other) {
    if (tri.v0 != a && tri.v0 != b) {
        tri.n0 = other;
    }
    else if (tri.v1 != a && tri.v1 != b) {
        tri.n1 = other;
    }
    else if (tri.v2 != a && tri.v2 != b) {
        tri.n2 = other;
    }
}

void setNeighbour(Tet& tet, int a, int b, int c,int other) {
    if (tet.v0 != a && tet.v0 != b && tet.v0 != c) {
        tet.n0 = other;
    }
    else if (tet.v1 != a && tet.v1 != b && tet.v1 != c) {
        tet.n1 = other;
    }
    else if (tet.v2 != a && tet.v2 != b && tet.v2 != c) {
        tet.n2 = other;
    }
    else if (tet.v3 != a && tet.v3 != b && tet.v3 != c) {
        tet.n3 = other;
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //


void buildNeighbours(Mesh& mesh) {
    std::map<std::pair<int,int>, std::vector<int>> edgeMap = buildEdgeMap(mesh);

    for (auto& t : mesh.meshTriangle) {
        t.n0 = -1;
        t.n1 = -1;
        t.n2 = -1;
    }

    for (const auto& e : edgeMap) {
        if (e.second.size() != 2) continue;
        Triangle& t = mesh.meshTriangle[e.second[0]];
        Triangle& k = mesh.meshTriangle[e.second[1]];
        int a = e.first.first;
        int b = e.first.second;

        setNeighbour(t, a, b, e.second[1]);
        setNeighbour(k, a, b, e.second[0]);
    }
}

void buildNeighbours(TetMesh& mesh) {
    std::map<std::array<int,3>, std::vector<int>> faceMap = buildFaceMap(mesh);

    for (auto& t : mesh.tets) {
        t.n0 = -1;
        t.n1 = -1;
        t.n2 = -1;
        t.n3 = -1;
    }

    for (const auto& e : faceMap) {
        if (e.second.size() != 2) continue;
        Tet& t = mesh.tets[e.second[0]];
        Tet& k = mesh.tets[e.second[1]];
        int a = e.first[0];
        int b = e.first[1];
        int c = e.first[2];

        setNeighbour(t, a, b, c, e.second[1]);
        setNeighbour(k, a, b, c, e.second[0]);
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

int walkToTriangle(const Mesh& mesh, int start, const Point& p) {
    int t = start;
    int limit = static_cast<int>(mesh.meshTriangle.size()) + 10;
    for (int steps = 0; steps < limit; steps++) {
        const Triangle& tri = mesh.meshTriangle[t];

        Point a = mesh.meshPoint[tri.v0];
        Point b = mesh.meshPoint[tri.v1];
        Point c = mesh.meshPoint[tri.v2];

        int next;
        if (orientation(b, c, p) < 0) {
            next = tri.n0;
        }
        else if (orientation(c, a, p) < 0) {
            next = tri.n1;
        }
        else if (orientation(a, b, p) < 0) {
            next = tri.n2;
        }
        else {
            return t;
        }

        if (next == -1) return -1;
        t = next;
    }

    return -1;
}

int walkToTet(const TetMesh& mesh, int start, const Point3& p) {
    int t = start;
    int limit = static_cast<int>(mesh.tets.size()) + 10;
    for (int steps = 0; steps < limit; steps++) {
        const Tet& tet = mesh.tets[t];

        Point3 a = mesh.points[tet.v0];
        Point3 b = mesh.points[tet.v1];
        Point3 c = mesh.points[tet.v2];
        Point3 d = mesh.points[tet.v3];

        int next;
        if      (orientation(p, b, c, d) < 0) next = tet.n0;   // beyond the face opposite a
        else if (orientation(a, p, c, d) < 0) next = tet.n1;   // beyond the face opposite b
        else if (orientation(a, b, p, d) < 0) next = tet.n2;   // beyond the face opposite c
        else if (orientation(a, b, c, p) < 0) next = tet.n3;   // beyond the face opposite d
        else return t;                                         // inside (or on a face)

        if (next == -1) return -1;
        t = next;
    }

    return -1;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

std::vector<int> tetsAroundEdge(const TetMesh& mesh, const std::vector<std::vector<int>>& vt, int u, int v) {
    std::vector<int> result;

    for (const auto& i : vt[u]) {
        const Tet& t = mesh.tets[i];
        if (t.v0 == u || t.v1 == u || t.v2 == u || t.v3 == u) {
            if (t.v0 == v || t.v1 == v || t.v2 == v || t.v3 == v) {
                result.push_back(static_cast<int>(i));
            }
        }
    }

    return result;
}

bool tryFlipEdge(TetMesh& mesh, const std::vector<std::vector<int>>& vt, int u, int v) {
    std::vector<int> around = tetsAroundEdge(mesh, vt, u, v);

    if (around.size() == 3) {
        Tet old[3] = {mesh.tets[around[0]], mesh.tets[around[1]], mesh.tets[around[2]]};
        Tet out[2];
        bool result =  flip32Tets(mesh, around[0], around[1], around[2], out);
        if (result == false) { return false; }

        double oldQual = worstQuality(mesh, old, 3);
        double newQual = worstQuality(mesh, out, 2);
        if (newQual <= oldQual) { return false; }

        apply32(mesh, around[0], around[1], around[2]);

        return true;
    }

    return false;
}

// Tries the 2-3 flip on each face of tet t that has a neighbour (interior faces only, so the hull
// is not touched). Accepted only if the worst of the three new tets beats the worst of the two old.
bool tryFlipFace(TetMesh& mesh, int t) {
    const Tet t1 = mesh.tets[t];
    int ns[4] = {t1.n0, t1.n1, t1.n2, t1.n3};
    for (int i = 0; i < 4; i++) {
        int n = ns[i];
        if (n == -1) continue;

        Tet out[3];
        if (!flip23Tets(mesh, t, n, out)) continue;

        Tet old[2] = {mesh.tets[t], mesh.tets[n]};
        if (worstQuality(mesh, out, 3) <= worstQuality(mesh, old, 2)) continue;

        apply23(mesh, t, n);
        return true;
    }
    return false;
}
