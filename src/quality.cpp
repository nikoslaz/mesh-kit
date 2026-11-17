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

double radiusEdgeRatio(const Point& a, const Point& b, const Point& c) {
    Point p = circumcenter(a, b, c);
    double R = std::hypot(p.x - a.x, p.y - a.y);

    double ab = std::hypot(b.x - a.x, b.y - a.y);
    double bc = std::hypot(c.x - b.x, c.y - b.y);
    double ca = std::hypot(a.x - c.x, a.y - c.y);
    double shortest = std::min({ab, bc, ca});

    return R / shortest;
}
