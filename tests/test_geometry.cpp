#include <algorithm>
#include <cmath>
#include <iostream>

#include "geometry.hpp"
#include "quality.hpp"

static int failures = 0;

static void checkNear(const char* name, double got, double expected) {
    bool ok = std::fabs(got - expected) < 1e-12;
    std::cout << (ok ? "PASS  " : "FAIL  ") << name
              << "  (got " << got << ", expected " << expected << ")\n";
    if (!ok) ++failures;
}

static void checkBool(const char* name, bool got, bool expected) {
    bool ok = (got == expected);
    std::cout << (ok ? "PASS  " : "FAIL  ") << name
              << "  (got " << (got ? "true" : "false")
              << ", expected " << (expected ? "true" : "false") << ")\n";
    if (!ok) ++failures;
}

int main() {
    std::cout << "--- orientation ---\n";
    checkNear("ccw unit triangle",  orientation({0,0}, {1,0}, {0,1}),  1.0);
    checkNear("cw unit triangle",   orientation({0,0}, {0,1}, {1,0}), -1.0);
    checkNear("ccw bigger (2x2)",   orientation({0,0}, {2,0}, {0,2}),  4.0);
    checkNear("collinear points",   orientation({0,0}, {1,1}, {2,2}),  0.0);

    std::cout << "--- inCircumcircle, a,b,c counter-clockwise ---\n";
    // circle through (0,0),(1,0),(0,1): centre (0.5,0.5), radius sqrt(0.5)
    checkBool("ccw d=(0.5,0.5) centre",    inCircumcircle({0,0}, {1,0}, {0,1}, {0.5, 0.5}),  true);
    checkBool("ccw d=(0.25,0.25) inside",  inCircumcircle({0,0}, {1,0}, {0,1}, {0.25, 0.25}), true);
    checkBool("ccw d=(-0.1,-0.1) outside", inCircumcircle({0,0}, {1,0}, {0,1}, {-0.1, -0.1}), false);
    checkBool("ccw d=(2,2) outside",       inCircumcircle({0,0}, {1,0}, {0,1}, {2, 2}),      false);
    checkBool("ccw d=(1,1) on circle",     inCircumcircle({0,0}, {1,0}, {0,1}, {1, 1}),      false);
    // circle through (0,0),(4,0),(0,4): centre (2,2), radius sqrt(8)
    checkBool("ccw big d=(2,2) centre",    inCircumcircle({0,0}, {4,0}, {0,4}, {2, 2}),      true);
    checkBool("ccw big d=(4,4) on circle", inCircumcircle({0,0}, {4,0}, {0,4}, {4, 4}),      false);
    checkBool("ccw big d=(5,5) outside",   inCircumcircle({0,0}, {4,0}, {0,4}, {5, 5}),      false);

    std::cout << "--- inCircumcircle, a,b,c clockwise (needs the orientation check) ---\n";
    checkBool("cw d=(0.5,0.5) centre",     inCircumcircle({0,0}, {0,1}, {1,0}, {0.5, 0.5}),  true);
    checkBool("cw d=(0.25,0.25) inside",   inCircumcircle({0,0}, {0,1}, {1,0}, {0.25, 0.25}), true);
    checkBool("cw d=(-0.1,-0.1) outside",  inCircumcircle({0,0}, {0,1}, {1,0}, {-0.1, -0.1}), false);
    checkBool("cw d=(2,2) outside",        inCircumcircle({0,0}, {0,1}, {1,0}, {2, 2}),      false);
    checkBool("cw d=(1,1) on circle",      inCircumcircle({0,0}, {0,1}, {1,0}, {1, 1}),      false);

    std::cout << "--- circumcenter ---\n";
    {
        auto checkPoint = [](const char* name, Point got, Point expected) {
            bool ok = std::fabs(got.x - expected.x) < 1e-12 && std::fabs(got.y - expected.y) < 1e-12;
            std::cout << (ok ? "PASS  " : "FAIL  ") << name
                      << "  (got " << got.x << "," << got.y
                      << ", expected " << expected.x << "," << expected.y << ")\n";
            if (!ok) ++failures;
        };
        checkPoint("unit right triangle",        circumcenter({0,0}, {1,0}, {0,1}),  {0.5, 0.5});
        checkPoint("bigger right triangle",      circumcenter({0,0}, {4,0}, {0,4}),  {2, 2});
        checkPoint("vertex order b,c,a",         circumcenter({1,0}, {0,1}, {0,0}),  {0.5, 0.5});
        checkPoint("vertex order c,b,a (clockwise)", circumcenter({0,1}, {1,0}, {0,0}), {0.5, 0.5});
        checkPoint("centre at the origin",       circumcenter({-1,0}, {1,0}, {0,1}), {0, 0});
        // obtuse triangle: the circumcenter lies outside the triangle (below it)
        checkPoint("obtuse triangle, centre outside", circumcenter({0,0}, {4,0}, {2,0.5}), {2, -3.75});

        // property: the centre is equally far from all three vertices
        struct Tri { Point a, b, c; };
        const Tri tris[] = {
            {{0,0}, {1,0}, {0,1}},
            {{0,0}, {4,0}, {2,0.5}},
            {{-3,2}, {5,-1}, {1,7}},
            {{100,100}, {101,100}, {100.5,100.9}},
            {{0,0}, {10,0}, {5,1e-3}},
        };
        bool equidistant = true;
        for (const auto& t : tris) {
            Point o = circumcenter(t.a, t.b, t.c);
            double ra = std::hypot(o.x - t.a.x, o.y - t.a.y);
            double rb = std::hypot(o.x - t.b.x, o.y - t.b.y);
            double rc = std::hypot(o.x - t.c.x, o.y - t.c.y);
            double tol = 1e-9 * std::max({ra, rb, rc, 1.0});
            if (std::fabs(ra - rb) > tol || std::fabs(ra - rc) > tol) equidistant = false;
        }
        checkBool("centre is equidistant from a, b, c (5 triangles)", equidistant, true);
    }

    std::cout << "--- radiusEdgeRatio ---\n";
    {
        const double sqrt3 = std::sqrt(3.0);
        // equilateral triangle with side 2: R = 2/sqrt(3), shortest = 2  ->  1/sqrt(3)
        checkNear("equilateral = 1/sqrt(3) (the best possible)",
                  radiusEdgeRatio({0,0}, {2,0}, {1,sqrt3}), 1.0 / sqrt3);
        // right isosceles: R = sqrt(2)/2, shortest = 1
        checkNear("right isosceles = sqrt(2)/2",
                  radiusEdgeRatio({0,0}, {1,0}, {0,1}), std::sqrt(2.0) / 2.0);
        checkNear("same triangle, other vertex order",
                  radiusEdgeRatio({0,1}, {0,0}, {1,0}), std::sqrt(2.0) / 2.0);
        checkNear("scaling the triangle does not change the ratio",
                  radiusEdgeRatio({0,0}, {100,0}, {0,100}), std::sqrt(2.0) / 2.0);
        // 30-60-90 triangle (0,0),(2,0),(0,2*tan30)... use (0,0),(sqrt3,0),(0,1): R = 1, shortest = 1
        checkNear("30-60-90 triangle = 1",
                  radiusEdgeRatio({0,0}, {sqrt3,0}, {0,1}), 1.0);
        // thin sliver (0,0),(10,0),(5,0.1): R = 125.05, shortest = 5.001  ->  ratio ~ 25
        checkNear("thin sliver ratio ~ 25.005",
                  std::round(radiusEdgeRatio({0,0}, {10,0}, {5,0.1}) * 1000.0) / 1000.0, 25.005);
        checkBool("equilateral is below the Ruppert bound sqrt(2)",
                  radiusEdgeRatio({0,0}, {2,0}, {1,sqrt3}) < std::sqrt(2.0), true);
        checkBool("sliver is above the Ruppert bound sqrt(2)",
                  radiusEdgeRatio({0,0}, {10,0}, {5,0.1}) > std::sqrt(2.0), true);
    }

    std::cout << "--- encroaches (segment (0,0)-(2,0): circle centre (1,0), radius 1) ---\n";
    checkBool("inside, above the segment (1,0.5)",  encroaches({0,0}, {2,0}, {1, 0.5}),  true);
    checkBool("inside, below the segment (1,-0.5)", encroaches({0,0}, {2,0}, {1, -0.5}), true);
    checkBool("outside (1,1.5)",                    encroaches({0,0}, {2,0}, {1, 1.5}),  false);
    checkBool("exactly on the circle (1,1)",        encroaches({0,0}, {2,0}, {1, 1}),    false);
    checkBool("endpoint a is not encroaching",      encroaches({0,0}, {2,0}, {0, 0}),    false);
    checkBool("endpoint b is not encroaching",      encroaches({0,0}, {2,0}, {2, 0}),    false);
    checkBool("midpoint of the segment (1,0)",      encroaches({0,0}, {2,0}, {1, 0}),    true);
    checkBool("far away (10,10)",                   encroaches({0,0}, {2,0}, {10, 10}),  false);
    checkBool("order of a and b does not matter",   encroaches({2,0}, {0,0}, {1, 0.5}),  true);

    std::cout << "--- orientation, nearly collinear (robustness) ---\n";
    // p moves on a 256x256 grid of spacing 2^-53 around (0.5, 0.5); q and r lie on y = x.
    // The sign of orientation must not depend on which vertex the triangle starts from.
    {
        const Point q{12, 12};
        const Point r{24, 24};
        const double u = std::ldexp(1.0, -53);
        auto sgn = [](double v) { return (v > 0) - (v < 0); };
        int inconsistent = 0;
        for (int i = 0; i < 256; ++i) {
            for (int j = 0; j < 256; ++j) {
                const Point p{0.5 + i * u, 0.5 + j * u};
                int s1 = sgn(orientation(p, q, r));
                int s2 = sgn(orientation(q, r, p));
                int s3 = sgn(orientation(r, p, q));
                if (s1 != s2 || s2 != s3) ++inconsistent;
            }
        }
        std::cout << "inconsistent orientations: " << inconsistent << " of 65536\n";
        checkBool("orientation sign is consistent under rotation", inconsistent == 0, true);
    }

    std::cout << "\n" << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
