#include <iostream>
#include <vector>

#include "delaunay.hpp"
#include "domain.hpp"
#include "geometry.hpp"
#include "refine.hpp"

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

    std::cout << "\n" << (failures == 0 ? "ALL REFINE TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
