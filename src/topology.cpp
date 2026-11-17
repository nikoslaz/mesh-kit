#include "topology.hpp"
#include "geometry.hpp"
#include <algorithm>
#include <set>
#include <utility>

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