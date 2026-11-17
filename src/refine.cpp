#include "refine.hpp"
#include <cmath>
#include <stdexcept>
#include <string>
#include "delaunay.hpp"
#include "geometry.hpp"
#include "quality.hpp"
#include "topology.hpp"

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

int splitSegment(Mesh& mesh, std::vector<Segment>& segments, std::size_t i) {
    Segment s = segments[i];
    const Point& A = mesh.meshPoint[s.a];
    const Point& B = mesh.meshPoint[s.b];
    Point m = { (A.x + B.x) / 2, (A.y + B.y) / 2 };
    int mi = insertPoint(mesh, m);
    segments[i] = {s.a, mi};
    segments.push_back({mi, s.b});
    return mi;
}

int findEncroachedSegment(const Mesh& mesh, const std::vector<Segment>& segments, const Point& p) {
    for (std::size_t i = 0; i < segments.size(); i++) {
        const Segment& s = segments[i];
        if (encroaches(mesh.meshPoint[s.a], mesh.meshPoint[s.b], p)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

int findBadTriangle(const Mesh& mesh, const std::vector<bool>& inside, double maxRatio) {
    for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
        if (!inside[ti]) continue;
        const Triangle& t = mesh.meshTriangle[ti];
        const Point& a = mesh.meshPoint[t.v0];
        const Point& b = mesh.meshPoint[t.v1];
        const Point& c = mesh.meshPoint[t.v2];
        if (radiusEdgeRatio(a, b, c) > maxRatio) {
            return static_cast<int>(ti);
        }
    }
    return -1;
}

std::vector<Segment> splitEncroachedSegments(Mesh& mesh, std::vector<Segment> segments) {
    bool changed = true;

    while (changed) {
        changed = false;
        for (std::size_t i = 0; i < segments.size(); i++) {
            if (isEncroached(mesh, segments[i])) {
                splitSegment(mesh, segments, i);
                changed = true;
                break;          
            }
        }
    }

    return segments;
}

void refineMesh(Mesh& mesh, std::vector<Segment>& segments, const std::vector<Point>& holes, double maxRatio) {
    segments = splitEncroachedSegments(mesh, segments);
    std::vector<bool> inside;

    for (int iter = 0; iter < 10000; iter++) {
        inside = insideTriangles(mesh, segments, holes);
        int bad = findBadTriangle(mesh, inside, maxRatio);
        if (bad == -1) break;      

        const Triangle& t = mesh.meshTriangle[bad];
        const Point& a = mesh.meshPoint[t.v0];
        const Point& b = mesh.meshPoint[t.v1];
        const Point& c = mesh.meshPoint[t.v2];
        Point m = circumcenter(a, b, c);
        
        int si = findEncroachedSegment(mesh, segments, m);
        if (si != -1) {
            splitSegment(mesh, segments, si);
            segments = splitEncroachedSegments(mesh, segments);
        }
        else {
            insertPoint(mesh, m);
        }
    }
}

Mesh meshDomain(const Domain& d, double maxRatio) {
    if (!(maxRatio > 0) || !std::isfinite(maxRatio))
        throw std::invalid_argument("meshDomain: maxRatio must be a positive number");
    for (std::size_t i = 0; i < d.points.size(); i++) {
        if (!std::isfinite(d.points[i].x) || !std::isfinite(d.points[i].y))
            throw std::invalid_argument("meshDomain: point " + std::to_string(i) + " has a non-finite coordinate");
    }
    const int n = static_cast<int>(d.points.size());
    for (std::size_t i = 0; i < d.segments.size(); i++) {
        const Segment& sg = d.segments[i];
        if (sg.a < 0 || sg.a >= n || sg.b < 0 || sg.b >= n)
            throw std::invalid_argument("meshDomain: segment " + std::to_string(i) + " refers to a point that does not exist");
        if (sg.a == sg.b)
            throw std::invalid_argument("meshDomain: segment " + std::to_string(i) + " has the same point at both ends");
    }
    for (std::size_t i = 0; i < d.holes.size(); i++) {
        if (!std::isfinite(d.holes[i].x) || !std::isfinite(d.holes[i].y))
            throw std::invalid_argument("meshDomain: hole " + std::to_string(i) + " has a non-finite coordinate");
    }

    Mesh mesh = makeSuperTriangle(d.points);

    for (const auto& p : d.points) {
        insertPoint(mesh, p);
    }

    std::vector<Segment> segs;
    for (const auto& s : d.segments) {
        segs.push_back({s.a + 3, s.b + 3});
    }

    refineMesh(mesh, segs, d.holes, maxRatio);
    std::vector<bool> inside = insideTriangles(mesh, segs, d.holes);

    return extractInside(mesh, inside);
}
