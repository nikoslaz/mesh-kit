#include <iostream>
#include <set>

#include "delaunay.hpp"
#include "geometry.hpp"
#include "topology.hpp"

static int failures = 0;

static void checkBool(const char* name, bool got, bool expected) {
    bool ok = (got == expected);
    std::cout << (ok ? "PASS  " : "FAIL  ") << name << "\n";
    if (!ok) ++failures;
}

int main() {
    std::cout << "--- buildEdgeMap: unit square (2 triangles) ---\n";
    {
        Mesh mesh = delaunay({{0,0}, {1,0}, {1,1}, {0,1}});
        auto edges = buildEdgeMap(mesh);

        int border = 0, interior = 0, other = 0;
        for (const auto& e : edges) {
            if (e.second.size() == 1) ++border;
            else if (e.second.size() == 2) ++interior;
            else ++other;
        }
        checkBool("2 triangles", mesh.meshTriangle.size() == 2, true);
        checkBool("5 edges in total", edges.size() == 5, true);
        checkBool("4 border edges (1 triangle each)", border == 4, true);
        checkBool("1 interior edge (the diagonal, 2 triangles)", interior == 1, true);
        checkBool("no edge belongs to 3 or more triangles", other == 0, true);

        bool keysOrdered = true;
        for (const auto& e : edges) if (e.first.first >= e.first.second) keysOrdered = false;
        checkBool("keys are stored as (smaller, larger)", keysOrdered, true);
    }

    std::cout << "--- buildEdgeMap: general properties on a bigger mesh ---\n";
    {
        std::vector<Point> pts;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                pts.push_back({i + 0.01 * j * i, j + 0.013 * i});   // slightly jittered 5x5 grid
        Mesh mesh = delaunay(pts);
        auto edges = buildEdgeMap(mesh);

        std::size_t incidences = 0;
        bool valid = true;
        for (const auto& e : edges) {
            incidences += e.second.size();
            if (e.second.empty() || e.second.size() > 2) valid = false;
            for (int t : e.second)
                if (t < 0 || t >= static_cast<int>(mesh.meshTriangle.size())) valid = false;
        }
        checkBool("every edge belongs to 1 or 2 valid triangles", valid, true);
        checkBool("each triangle contributes exactly 3 edges",
                  incidences == 3 * mesh.meshTriangle.size(), true);

        // Euler's formula for a triangulated convex region: V - E + F = 1
        long V = static_cast<long>(mesh.meshPoint.size());
        long E = static_cast<long>(edges.size());
        long F = static_cast<long>(mesh.meshTriangle.size());
        checkBool("Euler: V - E + F = 1", V - E + F == 1, true);
    }

    std::cout << "--- findTriangleContaining ---\n";
    {
        // true if p lies in (or on the border of) triangle ti, checked independently
        auto contains = [](const Mesh& m, int ti, const Point& p) {
            const Triangle& t = m.meshTriangle[ti];
            double o1 = orientation(m.meshPoint[t.v0], m.meshPoint[t.v1], p);
            double o2 = orientation(m.meshPoint[t.v1], m.meshPoint[t.v2], p);
            double o3 = orientation(m.meshPoint[t.v2], m.meshPoint[t.v0], p);
            bool allNonNeg = o1 >= 0 && o2 >= 0 && o3 >= 0;
            bool allNonPos = o1 <= 0 && o2 <= 0 && o3 <= 0;
            return allNonNeg || allNonPos;
        };

        Mesh sq = delaunay({{0,0}, {1,0}, {1,1}, {0,1}});
        const Point inside[] = {{0.9, 0.2}, {0.2, 0.9}, {0.5, 0.3}, {0.3, 0.5}, {0.05, 0.02}};
        bool allFound = true;
        for (const auto& p : inside) {
            int ti = findTriangleContaining(sq, p);
            if (ti < 0 || !contains(sq, ti, p)) allFound = false;
        }
        checkBool("points inside the unit square are found in a triangle that contains them",
                  allFound, true);

        checkBool("point far outside (5,5) returns -1",   findTriangleContaining(sq, {5, 5})   == -1, true);
        checkBool("point just outside (1.1,0.5) returns -1", findTriangleContaining(sq, {1.1, 0.5}) == -1, true);
        checkBool("point just outside (-0.01,0.5) returns -1", findTriangleContaining(sq, {-0.01, 0.5}) == -1, true);
        checkBool("a corner (0,0) is found", findTriangleContaining(sq, {0, 0}) >= 0, true);
        checkBool("a point on the border (0.5,0) is found", findTriangleContaining(sq, {0.5, 0}) >= 0, true);

        // the centre lies on the diagonal, which is shared by two triangles: either one is fine
        int c = findTriangleContaining(sq, {0.5, 0.5});
        checkBool("centre on the diagonal is found in a triangle that contains it",
                  c >= 0 && contains(sq, c, {0.5, 0.5}), true);

        Mesh empty;
        checkBool("empty mesh returns -1", findTriangleContaining(empty, {0, 0}) == -1, true);

        // bigger mesh: every sample point inside the grid must be found
        std::vector<Point> pts;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) pts.push_back({double(i), double(j)});
        Mesh grid = delaunay(pts);
        bool gridOk = true;
        for (int k = 0; k < 100; ++k) {
            Point p{0.05 + 4.9 * ((k * 37) % 100) / 100.0, 0.05 + 4.9 * ((k * 91) % 100) / 100.0};
            int ti = findTriangleContaining(grid, p);
            if (ti < 0 || !contains(grid, ti, p)) gridOk = false;
        }
        checkBool("100 sample points in a 6x6 grid are all found", gridOk, true);
    }

    std::cout << "\n" << (failures == 0 ? "ALL TOPOLOGY TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
