#include <cmath>
#include <iostream>

#include "geometry.hpp"

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
