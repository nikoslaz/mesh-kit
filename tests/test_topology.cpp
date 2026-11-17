#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

#include "delaunay.hpp"
#include "domain.hpp"
#include "geometry.hpp"
#include "quality.hpp"
#include "refine.hpp"
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

    std::cout << "--- insideTriangles: plate with a circular hole ---\n";
    {
        Domain d = makePlateWithHole(16);
        Mesh mesh = makeSuperTriangle(d.points);
        for (const auto& p : d.points) insertPoint(mesh, p);
        std::vector<Segment> segs;
        for (const auto& s : d.segments) segs.push_back({s.a + 3, s.b + 3});
        segs = splitEncroachedSegments(mesh, segs);

        std::vector<bool> inside = insideTriangles(mesh, segs, d.holes);
        checkBool("one answer per triangle", inside.size() == mesh.meshTriangle.size(), true);

        auto centroid = [&](const Triangle& t) {
            const Point& a = mesh.meshPoint[t.v0];
            const Point& b = mesh.meshPoint[t.v1];
            const Point& c = mesh.meshPoint[t.v2];
            return Point{(a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3};
        };

        std::size_t nIn = 0;
        bool noCorners = true, insideOk = true, outsideOk = true;
        double area = 0;
        for (std::size_t ti = 0; ti < inside.size(); ++ti) {
            const Triangle& t = mesh.meshTriangle[ti];
            bool touchesCorner = t.v0 < 3 || t.v1 < 3 || t.v2 < 3;
            Point c = centroid(t);
            double r = std::hypot(c.x - 5, c.y - 5);
            if (inside[ti]) {
                ++nIn;
                if (touchesCorner) noCorners = false;
                bool inPlate = c.x > 0 && c.x < 10 && c.y > 0 && c.y < 10;
                if (!inPlate || r < 1.9) insideOk = false;
                area += std::fabs(orientation(mesh.meshPoint[t.v0], mesh.meshPoint[t.v1],
                                              mesh.meshPoint[t.v2])) / 2;
            } else if (!touchesCorner) {
                // not inside and not around the plate: must be a triangle of the hole
                if (r > 2.0) outsideOk = false;
            }
        }
        checkBool("some triangles are inside", nIn > 0, true);
        checkBool("some triangles are not inside", nIn < inside.size(), true);
        checkBool("no inside triangle uses a super-triangle corner", noCorners, true);
        checkBool("inside triangles lie in the plate and outside the hole", insideOk, true);
        checkBool("every other non-inside triangle lies in the hole", outsideOk, true);

        int holeTri = findTriangleContaining(mesh, d.holes[0]);
        checkBool("triangle at the hole centre is not inside",
                  holeTri >= 0 && !inside[holeTri], true);

        // 10x10 plate minus the regular 16-gon of circumradius 2
        double expected = 100.0 - 0.5 * 16 * 2 * 2 * std::sin(2 * M_PI / 16);
        std::cout << "  (inside area " << area << ", expected " << expected << ")\n";
        checkBool("inside area = plate minus the hole polygon", std::fabs(area - expected) < 1e-9, true);
    }

    std::cout << "--- extractInside: small hand-made mesh ---\n";
    {
        // 6 points; triangle A uses 2,3,5 and B uses 3,5,4; triangle C (0,1,2) is not inside
        Mesh m;
        for (int i = 0; i < 6; ++i) m.addPoint({double(i), double(i * i)});
        m.addTriangle({2, 3, 5});
        m.addTriangle({3, 5, 4});
        m.addTriangle({0, 1, 2});

        Mesh r = extractInside(m, {true, true, false});
        checkBool("2 triangles kept", r.meshTriangle.size() == 2, true);
        checkBool("only the 4 used points are kept", r.meshPoint.size() == 4, true);
        // points are numbered in order of first use: 2->0, 3->1, 5->2, 4->3
        checkBool("triangle A becomes {0,1,2}",
                  r.meshTriangle[0].v0 == 0 && r.meshTriangle[0].v1 == 1 && r.meshTriangle[0].v2 == 2, true);
        checkBool("triangle B becomes {1,2,3}",
                  r.meshTriangle[1].v0 == 1 && r.meshTriangle[1].v1 == 2 && r.meshTriangle[1].v2 == 3, true);
        checkBool("new point 0 is old point 2",
                  r.meshPoint[0].x == 2 && r.meshPoint[0].y == 4, true);
        checkBool("new point 3 is old point 4",
                  r.meshPoint[3].x == 4 && r.meshPoint[3].y == 16, true);
        checkBool("the input mesh is not modified",
                  m.meshPoint.size() == 6 && m.meshTriangle.size() == 3, true);

        Mesh none = extractInside(m, {false, false, false});
        checkBool("nothing inside gives an empty mesh",
                  none.meshPoint.empty() && none.meshTriangle.empty(), true);

        Mesh all = extractInside(m, {true, true, true});
        checkBool("everything inside keeps 3 triangles and the 6 used points",
                  all.meshTriangle.size() == 3 && all.meshPoint.size() == 6, true);
    }

    std::cout << "--- extractInside: refined plate with a circular hole ---\n";
    {
        const double sqrt2 = std::sqrt(2.0);
        Domain d = makePlateWithHole(16);
        Mesh mesh = makeSuperTriangle(d.points);
        for (const auto& p : d.points) insertPoint(mesh, p);
        std::vector<Segment> segs;
        for (const auto& sg : d.segments) segs.push_back({sg.a + 3, sg.b + 3});
        refineMesh(mesh, segs, d.holes, sqrt2);
        std::vector<bool> inside = insideTriangles(mesh, segs, d.holes);

        std::size_t nInside = 0;
        for (bool b : inside) if (b) ++nInside;

        Mesh r = extractInside(mesh, inside);
        checkBool("one output triangle per inside triangle", r.meshTriangle.size() == nInside, true);
        checkBool("fewer points than before (corners and unused points are gone)",
                  r.meshPoint.size() < mesh.meshPoint.size(), true);

        std::vector<int> used(r.meshPoint.size(), 0);
        bool valid = true;
        for (const auto& t : r.meshTriangle)
            for (int v : {t.v0, t.v1, t.v2}) {
                if (v < 0 || v >= static_cast<int>(r.meshPoint.size())) valid = false;
                else ++used[v];
            }
        checkBool("all triangle indices are valid", valid, true);
        bool allUsed = true;
        for (int u : used) if (u == 0) allUsed = false;
        checkBool("every point is used by at least one triangle", allUsed, true);

        double area = 0, minAng = 180.0;
        for (const auto& t : r.meshTriangle) {
            const Point& a = r.meshPoint[t.v0];
            const Point& b = r.meshPoint[t.v1];
            const Point& c = r.meshPoint[t.v2];
            area += std::fabs(orientation(a, b, c)) / 2;
            minAng = std::min(minAng, minAngle(a, b, c) * 180.0 / M_PI);
        }
        double expected = 100.0 - 0.5 * 16 * 2 * 2 * std::sin(2 * M_PI / 16);
        std::cout << "  (final mesh: " << r.meshPoint.size() << " points, "
                  << r.meshTriangle.size() << " triangles, min angle " << minAng << " deg)\n";
        checkBool("area equals plate minus hole polygon", std::fabs(area - expected) < 1e-9, true);
        checkBool("minimum angle is at least 20 degrees", minAng >= 20.0, true);

        bool inPlate = true, outOfHole = true;
        for (const auto& p : r.meshPoint) {
            if (p.x < -1e-12 || p.x > 10 + 1e-12 || p.y < -1e-12 || p.y > 10 + 1e-12) inPlate = false;
            if (std::hypot(p.x - 5, p.y - 5) < 1.9) outOfHole = false;
        }
        checkBool("all points lie in the plate [0,10]x[0,10]", inPlate, true);
        checkBool("no point lies inside the hole", outOfHole, true);

        // topology of a plate with ONE hole: V - E + F = 0 (Euler characteristic of an annulus)
        auto edges = buildEdgeMap(r);
        bool edgesOk = true;
        std::size_t borderEdges = 0;
        for (const auto& e : edges) {
            if (e.second.size() > 2) edgesOk = false;
            if (e.second.size() == 1) ++borderEdges;
        }
        long V = static_cast<long>(r.meshPoint.size());
        long E = static_cast<long>(edges.size());
        long F = static_cast<long>(r.meshTriangle.size());
        checkBool("every edge belongs to 1 or 2 triangles", edgesOk, true);
        checkBool("Euler: V - E + F = 0 for a shape with one hole", V - E + F == 0, true);
        checkBool("border edges exist (outer wall and hole wall)", borderEdges > 0, true);
    }

    std::cout << "\n" << (failures == 0 ? "ALL TOPOLOGY TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
