#include "refine.hpp"
#include "delaunay.hpp"
#include "geometry.hpp"

static bool isEncroached(const Mesh& mesh, const Segment& s) {
    for (std::size_t i = 3; i < mesh.meshPoint.size(); i++) {
        const int idx = static_cast<int>(i);
        if (idx == s.a || idx == s.b) continue;
        if (encroaches(mesh.meshPoint[s.a], mesh.meshPoint[s.b], mesh.meshPoint[i])) {
            return true;
        }
    }
    return false;
}

std::vector<Segment> splitEncroachedSegments(Mesh& mesh, std::vector<Segment> segments) {
    bool changed = true;

    while (changed) {
        changed = false;
        for (std::size_t i = 0; i < segments.size(); i++) {
            if (isEncroached(mesh, segments[i])) {
                Segment s = segments[i];
                const Point& A = mesh.meshPoint[s.a];
                const Point& B = mesh.meshPoint[s.b];
                Point m = { (A.x + B.x) / 2, (A.y + B.y) / 2 };
                int mi = insertPoint(mesh, m);
                segments[i] = {s.a, mi};
                segments.push_back({mi, s.b});
                changed = true;
                break;          
            }
        }
    }

    return segments;
}