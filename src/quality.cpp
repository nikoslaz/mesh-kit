#include "quality.hpp"
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
