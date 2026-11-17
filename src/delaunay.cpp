#include "delaunay.hpp"
#include <algorithm>
#include <map>
#include <utility>

Mesh delaunay(const std::vector<Point>& points) {
    Mesh mesh;
    if (points.empty()) return mesh;

    double maxX = points[0].x;
    double minX = points[0].x;
    double maxY = points[0].y;
    double minY = points[0].y;
    double width, height, m;
    double cx, cy;

    for (const auto& point: points) {
        mesh.addPoint(point);
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

    for (std::size_t i = 0; i < points.size(); i++) {
        std::vector<std::size_t> bad;

        for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
            const Triangle& tri = mesh.meshTriangle[ti];
            if (inCircumcircle(mesh.meshPoint[tri.v0],
               mesh.meshPoint[tri.v1],
               mesh.meshPoint[tri.v2],
               mesh.meshPoint[i]) ) {
                bad.push_back(ti);
            }
        }

        std::map<std::pair<int,int>, int> edgeCount;
        for (const auto& ti : bad) {
            const Triangle& tri = mesh.meshTriangle[ti];
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

        std::vector<Triangle> keep;
        for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
            if (std::find(bad.begin(), bad.end(), ti) == bad.end()) {
                keep.push_back(mesh.meshTriangle[ti]);
            }
        }

        for (const auto& entry : edgeCount) {
                if (entry.second == 1) {
                    Triangle newTriang = {static_cast<int>(i) ,entry.first.first, entry.first.second};
                    keep.push_back(newTriang);
                }

        }

        mesh.meshTriangle = keep;
    }

    // Cleanup, after ALL points are inserted: drop every triangle that uses a
    // super-triangle corner (those corners have indices >= n), then drop the
    // corners themselves (they are the last 3 points in the vector).
    const int n = static_cast<int>(points.size());
    std::vector<Triangle> result;
    for (const auto& tri : mesh.meshTriangle) {
        if (tri.v0 < n && tri.v1 < n && tri.v2 < n) {
            result.push_back(tri);
        }
    }
    mesh.meshTriangle = result;
    mesh.meshPoint.resize(n);

    return mesh;
}