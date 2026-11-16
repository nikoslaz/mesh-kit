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

    std::cout << "\n" << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
