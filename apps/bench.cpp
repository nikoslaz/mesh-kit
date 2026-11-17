// meshkit_bench: wall-clock timings of Delaunay triangulation and quality meshing.
//   meshkit_bench            default sizes
//   meshkit_bench <max_n>    Delaunay up to max_n random points (doubling from 500)
// Build with optimisation (cmake -DCMAKE_BUILD_TYPE=Release) for meaningful numbers.
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "delaunay.hpp"
#include "domain.hpp"
#include "refine.hpp"

using Clock = std::chrono::steady_clock;

static double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

int main(int argc, char** argv) {
    int maxN = argc > 1 ? std::atoi(argv[1]) : 4000;
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Delaunay triangulation of random points in the unit square\n";
    std::cout << "  points   seconds   triangles\n";
    for (int n = 500; n <= maxN; n *= 2) {
        std::mt19937 rng(1);
        std::vector<Point> pts;
        for (int i = 0; i < n; i++) pts.push_back({rng() / 4294967296.0, rng() / 4294967296.0});
        auto t = Clock::now();
        Mesh m = delaunay(pts);
        std::cout << std::setw(8) << n << std::setw(10) << secondsSince(t)
                  << std::setw(12) << m.meshTriangle.size() << "\n";
    }

    std::cout << "\nmeshDomain (max radius-edge ratio 1.2) on a 10x10 plate with a circular hole\n";
    std::cout << "  hole sides   seconds    points   triangles\n";
    for (int n : {16, 32, 64, 128, 256}) {
        Domain d = makePlateWithHole(n);
        auto t = Clock::now();
        Mesh m = meshDomain(d, 1.2);
        std::cout << std::setw(12) << n << std::setw(10) << secondsSince(t)
                  << std::setw(10) << m.meshPoint.size() << std::setw(12) << m.meshTriangle.size() << "\n";
    }
    return 0;
}
