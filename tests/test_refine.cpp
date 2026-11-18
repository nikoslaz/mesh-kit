#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
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

// true if some triangle of the mesh has both a and b as vertices
static bool isMeshEdge(const Mesh& mesh, int a, int b) {
    for (const auto& t : mesh.meshTriangle) {
        int hits = 0;
        for (int v : {t.v0, t.v1, t.v2}) if (v == a || v == b) ++hits;
        if (hits == 2) return true;
    }
    return false;
}

static bool anyEncroached(const Mesh& mesh, const std::vector<Segment>& segs) {
    for (const auto& s : segs) {
        for (std::size_t i = 3; i < mesh.meshPoint.size(); ++i) {
            int idx = static_cast<int>(i);
            if (idx == s.a || idx == s.b) continue;
            if (encroaches(mesh.meshPoint[s.a], mesh.meshPoint[s.b], mesh.meshPoint[i]))
                return true;
        }
    }
    return false;
}

// Build the mesh of a domain the way the caller must: super-triangle first,
// domain points inserted, segment indices shifted by 3.
static Mesh buildBase(const Domain& d, std::vector<Segment>& segs) {
    Mesh mesh = makeSuperTriangle(d.points);
    for (const auto& p : d.points) insertPoint(mesh, p);
    segs.clear();
    for (const auto& s : d.segments) segs.push_back({s.a + 3, s.b + 3});
    return mesh;
}

int main() {
    std::cout << "--- splitSegment ---\n";
    {
        Domain d;
        d.points = {{0, 0}, {4, 0}};
        d.segments = {{0, 1}};
        std::vector<Segment> segs;
        Mesh mesh = buildBase(d, segs);          // points: 0,1,2 = corners, 3 = (0,0), 4 = (4,0)
        std::size_t pointsBefore = mesh.meshPoint.size();

        int mi = splitSegment(mesh, segs, 0);

        checkBool("returns the index of the new point", mi == static_cast<int>(pointsBefore), true);
        checkBool("one point was added", mesh.meshPoint.size() == pointsBefore + 1, true);
        checkBool("new point is the midpoint (2,0)",
                  mesh.meshPoint[mi].x == 2.0 && mesh.meshPoint[mi].y == 0.0, true);
        checkBool("now 2 segments", segs.size() == 2, true);
        checkBool("first half is {3, mi}",  segs[0].a == 3  && segs[0].b == mi, true);
        checkBool("second half is {mi, 4}", segs[1].a == mi && segs[1].b == 4,  true);
        checkBool("both halves are mesh edges",
                  isMeshEdge(mesh, segs[0].a, segs[0].b) && isMeshEdge(mesh, segs[1].a, segs[1].b), true);
    }

    std::cout << "--- findEncroachedSegment ---\n";
    {
        // two segments: S0 = (0,0)-(4,0) and S1 = (10,0)-(14,0); mesh indices 3..6
        Domain d;
        d.points = {{0, 0}, {4, 0}, {10, 0}, {14, 0}};
        d.segments = {{0, 1}, {2, 3}};
        std::vector<Segment> segs;
        Mesh mesh = buildBase(d, segs);

        checkBool("point (2,1) is inside the circle of S0 -> 0",
                  findEncroachedSegment(mesh, segs, {2, 1}) == 0, true);
        checkBool("point (12,-1) is inside the circle of S1 -> 1",
                  findEncroachedSegment(mesh, segs, {12, -1}) == 1, true);
        checkBool("point (7,0) is between them, in neither circle -> -1",
                  findEncroachedSegment(mesh, segs, {7, 0}) == -1, true);
        checkBool("point (2,2) is exactly on the circle of S0 -> -1",
                  findEncroachedSegment(mesh, segs, {2, 2}) == -1, true);
        checkBool("point (2,5) far above -> -1",
                  findEncroachedSegment(mesh, segs, {2, 5}) == -1, true);
        checkBool("an endpoint of a segment does not encroach it -> not 0",
                  findEncroachedSegment(mesh, segs, {0, 0}) != 0, true);

        std::vector<Segment> none;
        checkBool("empty segment list -> -1",
                  findEncroachedSegment(mesh, none, {2, 1}) == -1, true);
    }

    std::cout << "--- findBadTriangle ---\n";
    {
        const double sqrt2 = std::sqrt(2.0);

        // hand-made mesh: triangle 0 is good (equilateral), triangle 1 is a thin sliver
        Mesh m;
        m.addPoint({0, 0});        // 0
        m.addPoint({2, 0});        // 1
        m.addPoint({1, std::sqrt(3.0)});   // 2
        m.addPoint({10, 0});       // 3
        m.addPoint({20, 0});       // 4
        m.addPoint({15, 0.1});     // 5
        m.addTriangle({0, 1, 2});
        m.addTriangle({3, 4, 5});

        checkBool("sliver is found when both triangles are inside",
                  findBadTriangle(m, {true, true}, sqrt2) == 1, true);
        checkBool("sliver is ignored when it is not inside",
                  findBadTriangle(m, {true, false}, sqrt2) == -1, true);
        checkBool("good triangle alone gives -1",
                  findBadTriangle(m, {true, true}, 100.0) == -1, true);
        checkBool("nothing inside gives -1",
                  findBadTriangle(m, {false, false}, sqrt2) == -1, true);
        checkBool("a stricter limit (0.5) flags the equilateral one first",
                  findBadTriangle(m, {true, true}, 0.5) == 0, true);
        checkBool("the limit is strict: ratio equal to the limit is not bad",
                  findBadTriangle(m, {true, false}, 1.0 / std::sqrt(3.0)) == -1, true);

        Mesh none;
        checkBool("empty mesh gives -1", findBadTriangle(none, {}, sqrt2) == -1, true);

        // on a real mesh: the refined segments of the plate still leave bad triangles
        Domain d = makePlateWithHole(16);
        std::vector<Segment> segs;
        Mesh plate = buildBase(d, segs);
        segs = splitEncroachedSegments(plate, segs);
        std::vector<bool> inside = insideTriangles(plate, segs, d.holes);
        int bad = findBadTriangle(plate, inside, sqrt2);
        checkBool("the plate has a bad inside triangle before refinement", bad >= 0, true);
        if (bad >= 0) {
            const Triangle& tr = plate.meshTriangle[bad];
            double ratio = radiusEdgeRatio(plate.meshPoint[tr.v0], plate.meshPoint[tr.v1],
                                           plate.meshPoint[tr.v2]);
            checkBool("the triangle it returns is inside", inside[bad], true);
            checkBool("the triangle it returns really has ratio > sqrt(2)", ratio > sqrt2, true);
            // and no earlier inside triangle is bad (it returns the FIRST one)
            bool firstOne = true;
            for (int k = 0; k < bad; ++k) {
                if (!inside[k]) continue;
                const Triangle& q = plate.meshTriangle[k];
                if (radiusEdgeRatio(plate.meshPoint[q.v0], plate.meshPoint[q.v1],
                                    plate.meshPoint[q.v2]) > sqrt2) firstOne = false;
            }
            checkBool("it is the first bad inside triangle", firstOne, true);
        }
    }

    std::cout << "--- splitEncroachedSegments: one segment, one close point ---\n";
    {
        // segment (0,0)-(4,0); the point (2,1) is inside its diametral circle (radius 2)
        Domain d;
        d.points = {{0, 0}, {4, 0}, {2, 1}};
        d.segments = {{0, 1}};
        std::vector<Segment> segs;
        Mesh mesh = buildBase(d, segs);
        checkBool("segment starts encroached", anyEncroached(mesh, segs), true);
        std::size_t pointsBefore = mesh.meshPoint.size();
        segs = splitEncroachedSegments(mesh, segs);
        checkBool("no segment is encroached afterwards", anyEncroached(mesh, segs), false);
        checkBool("at least one point was added", mesh.meshPoint.size() > pointsBefore, true);
        checkBool("segment list grew with the points",
                  segs.size() - 1 == mesh.meshPoint.size() - pointsBefore, true);
        bool allEdges = true;
        for (const auto& s : segs) if (!isMeshEdge(mesh, s.a, s.b)) allEdges = false;
        checkBool("every segment is an edge of the mesh", allEdges, true);
    }

    std::cout << "--- splitEncroachedSegments: plate with a circular hole ---\n";
    {
        Domain d = makePlateWithHole(16);
        std::vector<Segment> segs;
        Mesh mesh = buildBase(d, segs);
        std::size_t pointsBefore = mesh.meshPoint.size();
        std::size_t segsBefore = segs.size();
        checkBool("the plate starts with encroached segments", anyEncroached(mesh, segs), true);

        segs = splitEncroachedSegments(mesh, segs);

        checkBool("no segment is encroached afterwards", anyEncroached(mesh, segs), false);
        checkBool("segments were split", segs.size() > segsBefore, true);
        checkBool("each split added exactly one point and one segment",
                  segs.size() - segsBefore == mesh.meshPoint.size() - pointsBefore, true);

        bool allEdges = true;
        for (const auto& s : segs) if (!isMeshEdge(mesh, s.a, s.b)) allEdges = false;
        checkBool("every segment is an edge of the mesh", allEdges, true);

        bool valid = true;
        for (const auto& s : segs)
            if (s.a < 3 || s.b < 3 ||
                s.a >= static_cast<int>(mesh.meshPoint.size()) ||
                s.b >= static_cast<int>(mesh.meshPoint.size())) valid = false;
        checkBool("segment indices are valid", valid, true);

        std::cout << "  (segments: " << segsBefore << " -> " << segs.size() << ")\n";
    }

    std::cout << "--- refineMesh: plate with a circular hole ---\n";
    {
        const double sqrt2 = std::sqrt(2.0);
        Domain d = makePlateWithHole(16);
        std::vector<Segment> segs;
        Mesh mesh = buildBase(d, segs);

        std::vector<bool> before = insideTriangles(mesh, segs, d.holes);
        int badBefore = findBadTriangle(mesh, before, sqrt2);
        std::size_t pointsBefore = mesh.meshPoint.size();

        refineMesh(mesh, segs, d.holes, sqrt2);

        std::vector<bool> inside = insideTriangles(mesh, segs, d.holes);
        checkBool("there were bad triangles before", badBefore >= 0, true);
        checkBool("no bad triangle remains", findBadTriangle(mesh, inside, sqrt2) == -1, true);
        checkBool("no segment is encroached", anyEncroached(mesh, segs), false);

        bool allEdges = true;
        for (const auto& s : segs) if (!isMeshEdge(mesh, s.a, s.b)) allEdges = false;
        checkBool("every segment is still an edge of the mesh", allEdges, true);

        double minAng = 180.0, area = 0;
        std::size_t nIn = 0;
        for (std::size_t ti = 0; ti < inside.size(); ++ti) {
            if (!inside[ti]) continue;
            ++nIn;
            const Triangle& tr = mesh.meshTriangle[ti];
            const Point& a = mesh.meshPoint[tr.v0];
            const Point& b = mesh.meshPoint[tr.v1];
            const Point& c = mesh.meshPoint[tr.v2];
            minAng = std::min(minAng, minAngle(a, b, c) * 180.0 / M_PI);
            area += std::fabs(orientation(a, b, c)) / 2;
        }
        double expected = 100.0 - 0.5 * 16 * 2 * 2 * std::sin(2 * M_PI / 16);
        std::cout << "  (points " << pointsBefore << " -> " << mesh.meshPoint.size()
                  << ", inside triangles " << nIn << ", min angle " << minAng
                  << " deg, area " << area << ")\n";
        checkBool("minimum angle of the inside triangles is at least 20 degrees", minAng >= 20.0, true);
        checkBool("the area of the shape is unchanged", std::fabs(area - expected) < 1e-9, true);
        checkBool("the mesh stayed reasonably small (< 3000 points)", mesh.meshPoint.size() < 3000, true);

        // running it again must not change anything
        std::size_t pointsAfter = mesh.meshPoint.size();
        refineMesh(mesh, segs, d.holes, sqrt2);
        checkBool("refining an already refined mesh adds no points",
                  mesh.meshPoint.size() == pointsAfter, true);
    }

    std::cout << "--- meshDomain: the whole pipeline ---\n";
    {
        const double sqrt2 = std::sqrt(2.0);
        Domain d = makePlateWithHole(16);
        Domain copy = d;

        Mesh r = meshDomain(d, sqrt2);

        checkBool("the domain is not modified (points)",   d.points.size() == copy.points.size(), true);
        bool segsSame = d.segments.size() == copy.segments.size();
        for (std::size_t i = 0; segsSame && i < d.segments.size(); ++i)
            if (d.segments[i].a != copy.segments[i].a || d.segments[i].b != copy.segments[i].b) segsSame = false;
        checkBool("the domain is not modified (segments keep their indices)", segsSame, true);

        double area = 0, minAng = 180.0, maxRatio = 0;
        for (const auto& t : r.meshTriangle) {
            const Point& a = r.meshPoint[t.v0];
            const Point& b = r.meshPoint[t.v1];
            const Point& c = r.meshPoint[t.v2];
            area += std::fabs(orientation(a, b, c)) / 2;
            minAng = std::min(minAng, minAngle(a, b, c) * 180.0 / M_PI);
            maxRatio = std::max(maxRatio, radiusEdgeRatio(a, b, c));
        }
        double expected = 100.0 - 0.5 * 16 * 2 * 2 * std::sin(2 * M_PI / 16);
        std::cout << "  (" << r.meshPoint.size() << " points, " << r.meshTriangle.size()
                  << " triangles, min angle " << minAng << " deg)\n";
        // The exact point count depends on the order bad triangles are fixed in, so check the
        // topology instead: a plate with one hole has V - E + F = 0, and all 16 hole corners
        // and 4 plate corners must still be there.
        auto edges = buildEdgeMap(r);
        long V = static_cast<long>(r.meshPoint.size());
        long E = static_cast<long>(edges.size());
        long F = static_cast<long>(r.meshTriangle.size());
        checkBool("Euler: V - E + F = 0 (one hole)", V - E + F == 0, true);
        int corners = 0;
        for (const auto& q : d.points)
            for (const auto& p : r.meshPoint)
                if (p.x == q.x && p.y == q.y) { ++corners; break; }
        checkBool("all 20 input points are in the mesh", corners == 20, true);
        checkBool("area equals plate minus hole", std::fabs(area - expected) < 1e-9, true);
        checkBool("every triangle has ratio <= sqrt(2)", maxRatio <= sqrt2, true);
        checkBool("min angle >= 20 degrees", minAng >= 20.0, true);

        Mesh again = meshDomain(d, sqrt2);
        checkBool("calling it twice gives the same result",
                  again.meshPoint.size() == r.meshPoint.size() &&
                  again.meshTriangle.size() == r.meshTriangle.size(), true);

        Mesh finer = meshDomain(d, 1.0);
        checkBool("a stricter limit (1.0) gives a finer mesh",
                  finer.meshTriangle.size() > r.meshTriangle.size(), true);
        double fineMin = 180.0;
        for (const auto& t : finer.meshTriangle)
            fineMin = std::min(fineMin, minAngle(finer.meshPoint[t.v0], finer.meshPoint[t.v1],
                                                 finer.meshPoint[t.v2]) * 180.0 / M_PI);
        std::cout << "  (limit 1.0: " << finer.meshPoint.size() << " points, "
                  << finer.meshTriangle.size() << " triangles, min angle " << fineMin << " deg)\n";
        checkBool("a stricter limit gives a larger minimum angle", fineMin > minAng, true);

        Domain empty;
        Mesh e = meshDomain(empty, sqrt2);
        checkBool("an empty domain gives an empty mesh",
                  e.meshPoint.empty() && e.meshTriangle.empty(), true);
    }

    std::cout << "--- meshDomain: invalid input is rejected ---\n";
    {
        auto throwsInvalid = [](const Domain& d, double ratio) {
            try { meshDomain(d, ratio); }
            catch (const std::invalid_argument&) { return true; }
            catch (...) { return false; }
            return false;
        };
        Domain ok;
        ok.points = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
        ok.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

        checkBool("a valid domain does not throw", throwsInvalid(ok, 1.5), false);
        checkBool("ratio 0 is rejected",        throwsInvalid(ok, 0.0), true);
        checkBool("negative ratio is rejected", throwsInvalid(ok, -1.0), true);
        checkBool("NaN ratio is rejected",      throwsInvalid(ok, std::nan("")), true);
        checkBool("infinite ratio is rejected", throwsInvalid(ok, INFINITY), true);

        Domain tooBig = ok;   tooBig.segments.push_back({2, 4});
        checkBool("segment index past the last point is rejected", throwsInvalid(tooBig, 1.5), true);
        Domain negative = ok; negative.segments.push_back({-1, 2});
        checkBool("negative segment index is rejected", throwsInvalid(negative, 1.5), true);
        Domain loop = ok;     loop.segments.push_back({2, 2});
        checkBool("segment from a point to itself is rejected", throwsInvalid(loop, 1.5), true);
        Domain nanPoint = ok; nanPoint.points[1].x = std::nan("");
        checkBool("NaN coordinate is rejected", throwsInvalid(nanPoint, 1.5), true);
        Domain nanHole = ok;  nanHole.holes.push_back({std::nan(""), 0});
        checkBool("NaN hole marker is rejected", throwsInvalid(nanHole, 1.5), true);
        Domain empty;
        checkBool("an empty domain is valid", throwsInvalid(empty, 1.5), false);
    }

    std::cout << "\n" << (failures == 0 ? "ALL REFINE TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
