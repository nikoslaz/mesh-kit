#include <cmath>
#include <iostream>

#include "domain.hpp"

static int failures = 0;

static void checkBool(const char* name, bool got, bool expected) {
    bool ok = (got == expected);
    std::cout << (ok ? "PASS  " : "FAIL  ") << name << "\n";
    if (!ok) ++failures;
}

int main() {
    std::cout << "--- makePlateWithHole(16) ---\n";
    Domain d = makePlateWithHole(16);

    checkBool("20 points",   d.points.size()   == 20, true);
    checkBool("20 segments", d.segments.size() == 20, true);
    checkBool("1 hole",      d.holes.size()    == 1,  true);

    bool onCircle = true;
    for (std::size_t i = 4; i < d.points.size(); ++i) {
        double r = std::hypot(d.points[i].x - 5.0, d.points[i].y - 5.0);
        if (std::fabs(r - 2.0) > 1e-12) onCircle = false;
    }
    checkBool("circle points are at distance 2 from (5,5)", onCircle, true);

    checkBool("last segment wraps around to {19, 4}",
              d.segments.back().a == 19 && d.segments.back().b == 4, true);
    checkBool("square closes: last square segment is {3, 0}",
              d.segments[3].a == 3 && d.segments[3].b == 0, true);

    bool indicesValid = true;
    for (const auto& s : d.segments) {
        if (s.a < 0 || s.b < 0 ||
            s.a >= static_cast<int>(d.points.size()) ||
            s.b >= static_cast<int>(d.points.size())) indicesValid = false;
    }
    checkBool("all segment indices point to existing points", indicesValid, true);

    std::cout << "\n" << (failures == 0 ? "ALL DOMAIN TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
