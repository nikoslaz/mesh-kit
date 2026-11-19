#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <vector>

#include <fstream>
#include <string>

#include "delaunay3d.hpp"
#include "vtk_writer.hpp"
#include "geometry.hpp"
#include "quality.hpp"
#include "topology.hpp"

static int failures = 0;

static void checkBool(const char* name, bool got, bool expected) {
    bool ok = (got == expected);
    std::cout << (ok ? "PASS  " : "FAIL  ") << name << "\n";
    if (!ok) ++failures;
}

static void checkNear(const char* name, double got, double expected, double tol) {
    bool ok = std::fabs(got - expected) <= tol;
    std::cout << (ok ? "PASS  " : "FAIL  ") << name << "  (got " << got << ", expected " << expected << ")\n";
    if (!ok) ++failures;
}

static const Point3& P(const TetMesh& m, int i) { return m.points[i]; }

static int countNonPositive(const TetMesh& m) {
    int n = 0;
    for (const auto& t : m.tets)
        if (orientation(P(m, t.v0), P(m, t.v1), P(m, t.v2), P(m, t.v3)) <= 0) ++n;
    return n;
}

static double totalVolume(const TetMesh& m) {
    double v = 0;
    for (const auto& t : m.tets) {
        const Point3 &a = P(m, t.v0), &b = P(m, t.v1), &c = P(m, t.v2), &d = P(m, t.v3);
        double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
        double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
        double wx = d.x - a.x, wy = d.y - a.y, wz = d.z - a.z;
        v += std::fabs(ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)) / 6.0;
    }
    return v;
}

// every neighbour link must be symmetric and share exactly 3 vertices
static int countBadNeighbours(const TetMesh& m) {
    int bad = 0;
    for (int t = 0; t < (int)m.tets.size(); ++t) {
        const Tet& k = m.tets[t];
        int ns[4] = {k.n0, k.n1, k.n2, k.n3};
        int vs[4] = {k.v0, k.v1, k.v2, k.v3};
        for (int i = 0; i < 4; ++i) {
            if (ns[i] == -1) continue;
            const Tet& o = m.tets[ns[i]];
            if (o.n0 != t && o.n1 != t && o.n2 != t && o.n3 != t) ++bad;
            int ov[4] = {o.v0, o.v1, o.v2, o.v3};
            int shared = 0;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    if (vs[a] == ov[b]) ++shared;
            if (shared != 3) ++bad;
        }
    }
    return bad;
}

// brute force: no vertex strictly inside any tet's circumsphere
static int countNotDelaunay(const TetMesh& m) {
    int n = 0;
    for (const auto& t : m.tets) {
        for (int v = 0; v < (int)m.points.size(); ++v) {
            if (v == t.v0 || v == t.v1 || v == t.v2 || v == t.v3) continue;
            if (inCircumcircle(P(m, t.v0), P(m, t.v1), P(m, t.v2), P(m, t.v3), P(m, v))) { ++n; break; }
        }
    }
    return n;
}

static void checkMesh(const char* label, const TetMesh& m) {
    std::cout << "--- " << label << " (" << m.tets.size() << " tets) ---\n";
    checkBool("all tets positive", countNonPositive(m) == 0, true);
    checkBool("neighbours symmetric, 3 shared vertices", countBadNeighbours(m) == 0, true);
    checkBool("empty circumsphere (Delaunay)", countNotDelaunay(m) == 0, true);
}

int main() {
    std::vector<Point3> cube = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};

    {
        TetMesh m = delaunay3d(cube);
        checkMesh("cube corners", m);
        checkNear("volume = 1", totalVolume(m), 1.0, 1e-12);
    }
    {
        std::vector<Point3> pts = cube;
        pts.push_back({0.5, 0.5, 0.5});
        TetMesh m = delaunay3d(pts);
        checkMesh("cube corners + centre", m);
        checkNear("volume = 1", totalVolume(m), 1.0, 1e-12);
    }
    {
        std::mt19937 g(42);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        std::vector<Point3> pts = cube;
        for (int i = 0; i < 500; ++i) pts.push_back({u(g), u(g), u(g)});
        TetMesh m = delaunay3d(pts);
        checkMesh("cube corners + 500 random points", m);
        checkNear("volume = 1", totalVolume(m), 1.0, 1e-9);
        checkBool("all points used", m.points.size() == pts.size(), true);
    }

    {
        std::vector<Point3> pts = cube;
        pts.push_back({0.5, 0.5, 0.5});
        pts.push_back({0.5, 0.5, 0.5});   // exact duplicate of the centre
        pts.push_back({0, 0, 0});         // duplicate of a corner
        pts.push_back({0.25, 0.3, 0.7});
        TetMesh m = delaunay3d(pts);
        checkMesh("duplicate points", m);
        checkNear("volume = 1", totalVolume(m), 1.0, 1e-12);
        checkBool("input numbering kept (all points present)", m.points.size() == pts.size(), true);
        std::vector<bool> used(m.points.size(), false);
        for (const auto& t : m.tets) { used[t.v0] = used[t.v1] = used[t.v2] = used[t.v3] = true; }
        int unused = 0;
        for (bool u : used) if (!u) ++unused;
        checkBool("exactly the 2 duplicates are unused", unused == 2, true);
    }
    {
        std::cout << "--- VTK writer ---\n";
        std::vector<Point3> pts = cube;
        pts.push_back({0.5, 0.5, 0.5});
        TetMesh m = delaunay3d(pts);
        const std::string path = "test_delaunay3d_out.vtk";
        writeVTK(m, path);
        std::ifstream in(path);
        std::string line;
        std::size_t nPoints = 0, nCells = 0, nTypes = 0, total = 0;
        int cellLines = 0;
        while (std::getline(in, line)) {
            if (line.rfind("POINTS ", 0) == 0) nPoints = std::stoul(line.substr(7));
            else if (line.rfind("CELLS ", 0) == 0) {
                std::size_t sp = line.find(' ', 6);
                nCells = std::stoul(line.substr(6, sp - 6));
                total = std::stoul(line.substr(sp + 1));
            } else if (line.rfind("CELL_TYPES ", 0) == 0) nTypes = std::stoul(line.substr(11));
            else if (line.rfind("4 ", 0) == 0) ++cellLines;
        }
        checkBool("POINTS matches the vertex count", nPoints == m.points.size(), true);
        checkBool("CELLS count and size (5 numbers per tet)", nCells == m.tets.size() && total == 5 * nCells, true);
        checkBool("one '4 v0 v1 v2 v3' line per tet", cellLines == (int)m.tets.size(), true);
        checkBool("CELL_TYPES count matches", nTypes == m.tets.size(), true);
        std::remove(path.c_str());
    }
    {
        std::cout << "--- tetRadiusRatio ---\n";
        Point3 r0{1, 1, 1}, r1{1, -1, -1}, r2{-1, 1, -1}, r3{-1, -1, 1};   // regular tet
        checkNear("regular tet = 1", tetRadiusRatio(r0, r1, r2, r3), 1.0, 1e-12);
        Point3 o{0,0,0}, x{1,0,0}, y{0,1,0}, z{0,0,1};
        checkNear("corner tet = 0.732", tetRadiusRatio(o, x, y, z), 3.0 * (std::sqrt(3.0) - 1) / 3.0 * 1.0, 1e-12);
        checkNear("swapping two vertices gives the same value", tetRadiusRatio(o, y, x, z), tetRadiusRatio(o, x, y, z), 1e-15);
        checkNear("scaling does not change it", tetRadiusRatio({0,0,0},{5,0,0},{0,5,0},{0,0,5}), tetRadiusRatio(o, x, y, z), 1e-12);
        checkBool("nearly flat tet is close to 0", tetRadiusRatio(o, x, y, {0.3, 0.3, 1e-6}) < 1e-4, true);
        checkNear("flat tet = 0", tetRadiusRatio(o, x, y, {0.3, 0.3, 0}), 0.0, 0.0);
        TetMesh m = delaunay3d(cube);
        auto q = tetQualities(m);
        bool inRange = q.size() == m.tets.size();
        for (double v : q) if (!(v > 0 && v <= 1 + 1e-12)) inRange = false;
        checkBool("tetQualities: one value per tet, all in (0, 1]", inRange, true);
    }

    std::cout << (failures == 0 ? "\nAll tests passed.\n" : "\nSome tests FAILED.\n");
    return failures == 0 ? 0 : 1;
}
