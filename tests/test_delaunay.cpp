#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "delaunay.hpp"
#include "geometry.hpp"

static int failures = 0;

// ---------- small helpers: an independent convex hull to compare against ----------

static std::vector<Point> convexHull(std::vector<Point> p) {
    std::sort(p.begin(), p.end(), [](const Point& a, const Point& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    std::size_t n = p.size();
    if (n < 3) return p;
    std::vector<Point> h(2 * n);
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; i++) {                       // lower hull
        while (k >= 2 && orientation(h[k - 2], h[k - 1], p[i]) <= 0) k--;
        h[k++] = p[i];
    }
    for (std::size_t i = n - 1, t = k + 1; i > 0; i--) {        // upper hull
        while (k >= t && orientation(h[k - 2], h[k - 1], p[i - 1]) <= 0) k--;
        h[k++] = p[i - 1];
    }
    h.resize(k - 1);
    return h;                                                    // counter-clockwise
}

static double polygonArea(const std::vector<Point>& h) {
    double s = 0;
    for (std::size_t i = 0; i < h.size(); i++) {
        const Point& a = h[i];
        const Point& b = h[(i + 1) % h.size()];
        s += a.x * b.y - a.y * b.x;
    }
    return std::fabs(s) / 2;
}

static bool onSegment(const Point& a, const Point& b, const Point& p) {
    return orientation(a, b, p) == 0 &&
           p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) &&
           p.y >= std::min(a.y, b.y) && p.y <= std::max(a.y, b.y);
}

// number of input points lying on the boundary of the convex hull (corners + collinear ones)
static std::size_t boundaryCount(const std::vector<Point>& pts, const std::vector<Point>& hull) {
    std::size_t count = 0;
    for (const auto& p : pts) {
        for (std::size_t j = 0; j < hull.size(); j++) {
            if (onSegment(hull[j], hull[(j + 1) % hull.size()], p)) { count++; break; }
        }
    }
    return count;
}

// ---------- the checks ----------

static void report(const std::string& name, const std::vector<std::string>& problems,
                   const std::string& info) {
    if (problems.empty()) {
        std::cout << "PASS  " << name << "  (" << info << ")\n";
    } else {
        std::cout << "FAIL  " << name << "  (" << info << ")\n";
        for (const auto& p : problems) std::cout << "        - " << p << "\n";
        ++failures;
    }
}

// Full check of a triangulation of 3 or more points.
static void checkTriangulation(const std::string& name, const std::vector<Point>& pts) {
    Mesh m = delaunay(pts);
    const std::size_t n = pts.size();
    std::vector<std::string> problems;

    if (m.meshPoint.size() != n)
        problems.push_back("expected " + std::to_string(n) + " points, got " +
                           std::to_string(m.meshPoint.size()) + " (super-triangle points not removed?)");

    bool indicesOk = true;
    for (const auto& t : m.meshTriangle) {
        for (int v : {t.v0, t.v1, t.v2})
            if (v < 0 || static_cast<std::size_t>(v) >= m.meshPoint.size()) indicesOk = false;
    }
    if (!indicesOk) {
        problems.push_back("a triangle uses a point index that doesn't exist");
        report(name, problems, std::to_string(m.meshTriangle.size()) + " triangles");
        return;
    }

    double totalArea = 0;
    bool areasOk = true, delaunayOk = true;
    for (const auto& t : m.meshTriangle) {
        const Point& a = m.meshPoint[t.v0];
        const Point& b = m.meshPoint[t.v1];
        const Point& c = m.meshPoint[t.v2];
        double area = std::fabs(orientation(a, b, c)) / 2;
        if (area < 1e-12) areasOk = false;
        totalArea += area;
        for (std::size_t i = 0; i < n; i++)
            if (inCircumcircle(a, b, c, pts[i])) delaunayOk = false;
    }
    if (!areasOk)    problems.push_back("a triangle has zero area (degenerate)");
    if (!delaunayOk) problems.push_back("Delaunay property violated: a point is inside some circumcircle");

    std::vector<Point> hull = convexHull(pts);
    std::size_t expectedTris = 2 * n - 2 - boundaryCount(pts, hull);
    double expectedArea = polygonArea(hull);
    if (m.meshTriangle.size() != expectedTris)
        problems.push_back("expected " + std::to_string(expectedTris) + " triangles (2N-2-h), got " +
                           std::to_string(m.meshTriangle.size()));
    if (std::fabs(totalArea - expectedArea) > 1e-9 * std::max(1.0, expectedArea))
        problems.push_back("triangles cover area " + std::to_string(totalArea) +
                           " but the convex hull has area " + std::to_string(expectedArea));

    report(name, problems,
           std::to_string(n) + " points, " + std::to_string(m.meshTriangle.size()) + " triangles");
}

// Fewer than 3 points: no triangles, no crash.
static void checkTooFewPoints(const std::string& name, const std::vector<Point>& pts) {
    Mesh m = delaunay(pts);
    std::vector<std::string> problems;
    if (!m.meshTriangle.empty()) problems.push_back("expected 0 triangles");
    if (m.meshPoint.size() != pts.size()) problems.push_back("point count changed");
    report(name, problems, std::to_string(pts.size()) + " points, " +
                           std::to_string(m.meshTriangle.size()) + " triangles");
}

static std::vector<Point> randomPoints(unsigned seed, int count, double scale, double ox, double oy) {
    std::mt19937 rng(seed);
    std::vector<Point> pts;
    for (int i = 0; i < count; i++) {
        double x = rng() / 4294967296.0;   // [0,1), same on every platform
        double y = rng() / 4294967296.0;
        pts.push_back({ox + scale * x, oy + scale * y});
    }
    return pts;
}

static std::vector<Point> grid(int w, int h, double ox = 0, double oy = 0) {
    std::vector<Point> pts;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) pts.push_back({ox + i, oy + j});
    return pts;
}

int main() {
    std::cout << "--- tiny inputs ---\n";
    checkTooFewPoints("0 points", {});
    checkTooFewPoints("1 point",  {{3, 4}});
    checkTooFewPoints("2 points", {{0, 0}, {1, 1}});

    std::cout << "--- small hand-checkable inputs ---\n";
    checkTriangulation("single triangle",            {{0, 0}, {1, 0}, {0, 1}});
    checkTriangulation("unit square",                {{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    checkTriangulation("unit square, reversed order",{{0, 1}, {1, 1}, {1, 0}, {0, 0}});
    checkTriangulation("square + centre point",      {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}});
    checkTriangulation("centre point inserted first",{{0.5, 0.5}, {0, 0}, {1, 0}, {1, 1}, {0, 1}});
    checkTriangulation("thin triangle",              {{0, 0}, {10, 0}, {5, 0.5}});

    std::cout << "--- positions: negative / far from the origin ---\n";
    checkTriangulation("square at negative coords",  {{-6, -6}, {-5, -6}, {-5, -5}, {-6, -5}});
    checkTriangulation("square far from origin",     {{1000, 1000}, {1001, 1000}, {1001, 1001}, {1000, 1001}});

    std::cout << "--- grids (many cocircular and collinear points) ---\n";
    checkTriangulation("3x3 grid",  grid(3, 3));
    checkTriangulation("5x5 grid",  grid(5, 5));
    checkTriangulation("4x6 grid at (-2,-3)", grid(4, 6, -2, -3));

    std::cout << "--- random points (checked against an independent convex hull) ---\n";
    checkTriangulation("random 10 points",   randomPoints(1, 10, 1.0, 0, 0));
    checkTriangulation("random 50 points",   randomPoints(2, 50, 1.0, 0, 0));
    checkTriangulation("random 200 points",  randomPoints(3, 200, 100.0, -50, -50));
    checkTriangulation("random 500 points",  randomPoints(4, 500, 1000.0, 0, 0));
    checkTriangulation("random 100 points, far from origin", randomPoints(5, 100, 10.0, 5000, -5000));

    std::cout << "\n" << (failures == 0 ? "ALL DELAUNAY TESTS PASSED" : "SOME DELAUNAY TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
