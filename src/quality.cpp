#include "quality.hpp"
#include "geometry.hpp"
#include "topology.hpp"
#include "smoothing.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>

// Angle (in radians) at `corner` formed by the segments corner->p and corner->q.
double angleAt(const Point& corner, const Point& p, const Point& q) {
    double ux = p.x - corner.x;
    double uy = p.y - corner.y;
    double vx = q.x - corner.x;
    double vy = q.y - corner.y;

    double cross = ux*vy - uy*vx;
    double dot   = ux*vx + uy*vy;

    return std::atan2(std::fabs(cross), dot);
}

// Smallest interior angle (in radians) of triangle abc.
double minAngle(const Point& a, const Point& b, const Point& c) {
    return std::min({angleAt(a, b, c), angleAt(b, c, a), angleAt(c, a, b)});
}

std::vector<double> minAngles(const Mesh& mesh) {
    std::vector<double> result;
    for (std::size_t ti = 0; ti < mesh.meshTriangle.size(); ti++) {
        const Triangle& tri = mesh.meshTriangle[ti];
        double angle = minAngle(mesh.meshPoint[tri.v0], mesh.meshPoint[tri.v1], mesh.meshPoint[tri.v2]) * 180.0 / std::numbers::pi;
        result.push_back(angle);
    }
    return result;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

double radiusEdgeRatio(const Point& a, const Point& b, const Point& c) {
    Point p = circumcenter(a, b, c);
    double R = std::hypot(p.x - a.x, p.y - a.y);

    double ab = std::hypot(b.x - a.x, b.y - a.y);
    double bc = std::hypot(c.x - b.x, c.y - b.y);
    double ca = std::hypot(a.x - c.x, a.y - c.y);
    double shortest = std::min({ab, bc, ca});

    return R / shortest;
}

static Point3 sub(const Point3& p, const Point3& q) {
    return  Point3 {p.x - q.x, p.y - q.y, p.z - q.z};
}     

static double dot(const Point3& p, const Point3& q) {
    return (p.x*q.x + p.y*q.y + p.z*q.z);
}

static Point3 cross(const Point3& p, const Point3& q){
    return Point3{p.y*q.z - p.z*q.y, p.z*q.x - p.x*q.z, p.x*q.y - p.y*q.x};
}

static double norm(const Point3& p) {
    return std::sqrt(dot(p, p));
}

// Radius ratio 3r/R of a tetrahedron: 1 for a regular tet, 0 for a flat one.
double tetRadiusRatio(const Point3& a, const Point3& b, const Point3& c, const Point3& d) {
    Point3 u = sub(b, a);
    Point3 v = sub(c, a);
    Point3 w = sub(d, a);

    double vol6 = std::fabs(dot(u, cross(v, w)));   // 6 x volume
    if (vol6 == 0) return 0;                         // flat tet
    double V = vol6 / 6;

    double S = norm(cross(u, v)) / 2                  // face abc
             + norm(cross(u, w)) / 2                  // face abd
             + norm(cross(v, w)) / 2                  // face acd
             + norm(cross(sub(c, b), sub(d, b))) / 2; // face bcd
    double r = 3 * V / S;                             // inradius

    // circumradius: |circumcentre - a| = |top| / (2 * 6V)
    Point3 cu = cross(v, w), cv = cross(w, u), cw = cross(u, v);
    double uu = dot(u, u), vv = dot(v, v), ww = dot(w, w);
    Point3 top{uu * cu.x + vv * cv.x + ww * cw.x,
               uu * cu.y + vv * cv.y + ww * cw.y,
               uu * cu.z + vv * cv.z + ww * cw.z};
    double R = norm(top) / (2 * vol6);

    return 3 * r / R;
}

std::vector<double> tetQualities(const TetMesh& mesh) {
    std::vector<double> q;
    q.reserve(mesh.tets.size());
    for (const auto& t : mesh.tets)
        q.push_back(tetRadiusRatio(mesh.points[t.v0], mesh.points[t.v1], mesh.points[t.v2], mesh.points[t.v3]));
    return q;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

double worstQuality(const TetMesh& mesh, const Tet* tets, int count) {
    double minR = 1.0;

    for (int i = 0; i < count; i++) {
        const Tet& t = tets[i];
        double r = tetRadiusRatio(mesh.points[t.v0], mesh.points[t.v1], mesh.points[t.v2], mesh.points[t.v3]);
        if (r < minR) {
            minR = r;
        }
    }

    return minR;
}

std::vector<int> findSlivers(const TetMesh& mesh, double threshold) {
    std::vector<int> slivers;

    for (std::size_t i = 0; i < mesh.tets.size(); i++) {
        const Tet& t = mesh.tets[i];
        if (tetRadiusRatio(mesh.points[t.v0], mesh.points[t.v1], mesh.points[t.v2], mesh.points[t.v3]) < threshold ) {
            slivers.push_back(static_cast<int>(i));
        }
    }

    return slivers;
}

// Flips slivers away. Each pass tries every sliver once; a tet already fixed by an earlier flip in
// the same pass is skipped. Edges with a hull vertex are never flipped, so the boundary stays the same.
// Every accepted flip raises the worst quality of the tets it replaces. Returns the number of flips.
int removeSlivers(TetMesh& mesh, double threshold, int maxPasses) {
    int flips = 0;
    for (int pass = 0; pass < maxPasses; pass++) {
        std::vector<int> slivers = findSlivers(mesh, threshold);
        if (slivers.empty()) break;

        auto vt   = vertexTets(mesh);
        auto hull = hullVertices(mesh);   // flips never change the boundary, so this stays valid

        int flipsThisPass = 0;
        for (int s : slivers) {
            if (s >= static_cast<int>(mesh.tets.size())) continue;   // removed by an earlier 3-2 flip
            const Tet t = mesh.tets[s];
            if (worstQuality(mesh, &t, 1) >= threshold) continue;    // already replaced by a better tet

            bool flipped = false;
            int V[4] = {t.v0, t.v1, t.v2, t.v3};
            for (int a = 0; a < 4 && !flipped; a++) {
                for (int b = a + 1; b < 4 && !flipped; b++) {
                    if (hull[V[a]] || hull[V[b]]) continue;
                    if (tryFlipEdge(mesh, vt, V[a], V[b])) flipped = true;
                }
            }
            if (!flipped && tryFlipFace(mesh, s)) flipped = true;

            if (flipped) {
                flipsThisPass++;
                vt = vertexTets(mesh);   // tet indices changed
            }
        }

        if (flipsThisPass == 0) break;
        flips += flipsThisPass;
    }
    return flips;
}

// Alternates sliver removal with Laplacian smoothing (guarded, so no move makes a tet worse).
// Stops early when a round changes nothing. Returns the total number of flips.
int improveSlivers(TetMesh& mesh, double threshold, int rounds) {
    int flips = 0;
    for (int r = 0; r < rounds; r++) {
        int f = removeSlivers(mesh, threshold, 100000);
        int moves = laplacianSmooth3d(mesh, 3);
        flips += f;
        if (f == 0 && moves == 0) break;
    }
    return flips;
}
