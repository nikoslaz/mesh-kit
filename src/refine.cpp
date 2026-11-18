#include "refine.hpp"
#include <cmath>
#include <algorithm>
#include <utility>
#include <map>
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

// ---------------------------------------------------------------------------------------
// Fast refinement. Instead of rescanning the whole mesh on every step, the refiner keeps:
//   - segOf:      which mesh edge is which segment
//   - inside:     inside/outside flag per triangle, updated from each insertion
//   - candidates: triangles that may be bad (only new triangles can become bad)
// and uses the report from insertPoint (InsertChange) to update them locally.
// ---------------------------------------------------------------------------------------

using EdgeKey = std::pair<int,int>;

static EdgeKey edgeKey(int a, int b) {
    return {std::min(a, b), std::max(a, b)};
}

// The triangle across edge (a, b) of tri, or -1.
static int neighbourAcross(const Triangle& tri, int a, int b) {
    if (tri.v0 != a && tri.v0 != b) return tri.n0;
    if (tri.v1 != a && tri.v1 != b) return tri.n1;
    return tri.n2;
}

static int oppositeVertex(const Triangle& tri, int a, int b) {
    if (tri.v0 != a && tri.v0 != b) return tri.v0;
    if (tri.v1 != a && tri.v1 != b) return tri.v1;
    return tri.v2;
}

static bool hasVertex(const Triangle& tri, int v) {
    return tri.v0 == v || tri.v1 == v || tri.v2 == v;
}

// Is segment s encroached? If s is an edge of the (Delaunay) mesh, a point encroaches it
// only if a vertex opposite s in one of its two triangles does, so two tests are enough.
// If s is missing from the mesh, fall back to checking every point.
static bool isEncroachedLocal(const Mesh& mesh, const Segment& s) {
    const Point& A = mesh.meshPoint[s.a];
    const Point& B = mesh.meshPoint[s.b];
    Point mid{(A.x + B.x) / 2, (A.y + B.y) / 2};
    int start = mesh.lastTriangle;
    if (start < 0 || start >= static_cast<int>(mesh.meshTriangle.size())) start = 0;
    int t = walkToTriangle(mesh, start, mid);
    if (t == -1 || !hasVertex(mesh.meshTriangle[t], s.a) || !hasVertex(mesh.meshTriangle[t], s.b)) {
        return isEncroached(mesh, s);
    }
    int sides[2] = {t, neighbourAcross(mesh.meshTriangle[t], s.a, s.b)};
    for (int side : sides) {
        if (side == -1) continue;
        int o = oppositeVertex(mesh.meshTriangle[side], s.a, s.b);
        if (o < 3) continue;   // super-triangle corners never count
        if (encroaches(A, B, mesh.meshPoint[o])) return true;
    }
    return false;
}

struct Refiner {
    Mesh& mesh;
    std::vector<Segment>& segments;
    const std::vector<Point>& holes;
    std::map<EdgeKey, int> segOf;
    std::vector<bool> inside;
    std::vector<int> candidates;

    Refiner(Mesh& m, std::vector<Segment>& s, const std::vector<Point>& h)
        : mesh(m), segments(s), holes(h) {
        for (std::size_t i = 0; i < segments.size(); i++)
            segOf[edgeKey(segments[i].a, segments[i].b)] = static_cast<int>(i);
    }

    // Segments with an edge in a removed or an added triangle: the only ones whose
    // encroachment can have changed.
    void touchedSegments(const InsertChange& ch, std::vector<int>& out) const {
        auto addEdges = [&](const Triangle& t) {
            const int v[3] = {t.v0, t.v1, t.v2};
            for (int k = 0; k < 3; k++) {
                auto it = segOf.find(edgeKey(v[k], v[(k + 1) % 3]));
                if (it != segOf.end()) out.push_back(it->second);
            }
        };
        for (const auto& t : ch.removed) addEdges(t);
        for (int t : ch.added) addEdges(mesh.meshTriangle[t]);
    }

    // Split segment i at its midpoint; returns what the insertion changed.
    InsertChange split(int i) {
        Segment s = segments[i];
        const Point& A = mesh.meshPoint[s.a];
        const Point& B = mesh.meshPoint[s.b];
        Point m = {(A.x + B.x) / 2, (A.y + B.y) / 2};
        InsertChange ch;
        int mi = insertPoint(mesh, m, &ch);
        segOf.erase(edgeKey(s.a, s.b));
        segments[i] = {s.a, mi};
        segments.push_back({mi, s.b});
        segOf[edgeKey(s.a, mi)] = i;
        segOf[edgeKey(mi, s.b)] = static_cast<int>(segments.size()) - 1;
        return ch;
    }

    // Split segments until none is encroached. `work` holds the segments to check.
    // Returns the triangles created along the way.
    std::vector<int> splitEncroached(std::vector<int> work) {
        std::vector<int> created;
        while (!work.empty()) {
            int i = work.back();
            work.pop_back();
            if (!isEncroachedLocal(mesh, segments[i])) continue;
            InsertChange ch = split(i);
            work.push_back(i);
            work.push_back(static_cast<int>(segments.size()) - 1);
            touchedSegments(ch, work);
            created.insert(created.end(), ch.added.begin(), ch.added.end());
        }
        return created;
    }

    // Inside/outside of every triangle from scratch: flood from the super-triangle and the
    // holes, stopping at segments. Only valid while every segment is an edge of the mesh.
    void floodInside() {
        const int T = static_cast<int>(mesh.meshTriangle.size());
        std::vector<bool> painted(T, false);
        std::vector<int> stack;
        for (int t = 0; t < T; t++) {
            const Triangle& tri = mesh.meshTriangle[t];
            if (tri.v0 < 3 || tri.v1 < 3 || tri.v2 < 3) {
                painted[t] = true;
                stack.push_back(t);
            }
        }
        for (const auto& h : holes) {
            int t = findTriangleContaining(mesh, h);
            if (t != -1 && !painted[t]) {
                painted[t] = true;
                stack.push_back(t);
            }
        }
        while (!stack.empty()) {
            int t = stack.back();
            stack.pop_back();
            const Triangle& tri = mesh.meshTriangle[t];
            const int v[3] = {tri.v0, tri.v1, tri.v2};
            const int n[3] = {tri.n2, tri.n0, tri.n1};   // neighbour across edge v[k] - v[k+1]
            for (int k = 0; k < 3; k++) {
                int nb = n[k];
                if (nb == -1 || painted[nb]) continue;
                if (segOf.count(edgeKey(v[k], v[(k + 1) % 3]))) continue;   // a wall
                painted[nb] = true;
                stack.push_back(nb);
            }
        }
        inside.assign(T, false);
        for (int t = 0; t < T; t++) inside[t] = !painted[t];
    }

    // New triangles take the inside flag of the cavity triangle they replaced. Correct when
    // no segment crosses the cavity, which holds for a circumcenter that encroaches nothing.
    void inheritInside(const InsertChange& ch) {
        std::vector<bool> from(ch.addedFrom.size());
        for (std::size_t i = 0; i < from.size(); i++) from[i] = inside[ch.addedFrom[i]];
        inside.resize(mesh.meshTriangle.size(), false);
        for (std::size_t i = 0; i < ch.added.size(); i++) inside[ch.added[i]] = from[i];
    }

    bool isBad(int t, double maxRatio) const {
        if (t < 0 || t >= static_cast<int>(mesh.meshTriangle.size()) || !inside[t]) return false;
        const Triangle& tri = mesh.meshTriangle[t];
        return radiusEdgeRatio(mesh.meshPoint[tri.v0], mesh.meshPoint[tri.v1],
                               mesh.meshPoint[tri.v2]) > maxRatio;
    }

    // A segment that point p (not yet inserted) would encroach, or -1. Only segments that are
    // edges of p's cavity can be encroached by p.
    int encroachedBy(const Point& p) const {
        std::vector<int> cavity = findCavityOf(mesh, p);
        if (cavity.empty()) return findEncroachedSegment(mesh, segments, p);
        int best = -1;
        for (int t : cavity) {
            const Triangle& tri = mesh.meshTriangle[t];
            const int v[3] = {tri.v0, tri.v1, tri.v2};
            for (int k = 0; k < 3; k++) {
                auto it = segOf.find(edgeKey(v[k], v[(k + 1) % 3]));
                if (it == segOf.end()) continue;
                const Segment& s = segments[it->second];
                if (encroaches(mesh.meshPoint[s.a], mesh.meshPoint[s.b], p) &&
                    (best == -1 || it->second < best))
                    best = it->second;
            }
        }
        return best;
    }

    void run(double maxRatio) {
        std::vector<int> all(segments.size());
        for (std::size_t i = 0; i < all.size(); i++) all[i] = static_cast<int>(i);
        splitEncroached(all);
        floodInside();
        for (int t = 0; t < static_cast<int>(mesh.meshTriangle.size()); t++) candidates.push_back(t);

        for (int iter = 0; iter < 10000 && !candidates.empty(); ) {
            int t = candidates.back();
            candidates.pop_back();
            if (!isBad(t, maxRatio)) continue;
            iter++;

            const Triangle& tri = mesh.meshTriangle[t];
            Point c = circumcenter(mesh.meshPoint[tri.v0], mesh.meshPoint[tri.v1], mesh.meshPoint[tri.v2]);
            int si = encroachedBy(c);
            if (si != -1) {
                // split the segment instead; t may still be bad afterwards, so keep it
                candidates.push_back(t);
                InsertChange ch = split(si);
                std::vector<int> work = {si, static_cast<int>(segments.size()) - 1};
                touchedSegments(ch, work);
                std::vector<int> created = splitEncroached(work);
                created.insert(created.end(), ch.added.begin(), ch.added.end());
                floodInside();
                candidates.insert(candidates.end(), created.begin(), created.end());
            }
            else {
                InsertChange ch;
                insertPoint(mesh, c, &ch);
                inheritInside(ch);
                candidates.insert(candidates.end(), ch.added.begin(), ch.added.end());
            }
        }
    }
};

std::vector<Segment> splitEncroachedSegments(Mesh& mesh, std::vector<Segment> segments) {
    static const std::vector<Point> noHoles;
    buildNeighbours(mesh);
    Refiner r(mesh, segments, noHoles);
    std::vector<int> all(segments.size());
    for (std::size_t i = 0; i < all.size(); i++) all[i] = static_cast<int>(i);
    r.splitEncroached(all);
    return segments;
}

void refineMesh(Mesh& mesh, std::vector<Segment>& segments, const std::vector<Point>& holes, double maxRatio) {
    buildNeighbours(mesh);
    Refiner r(mesh, segments, holes);
    r.run(maxRatio);
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
