#include "domain.hpp"
#include <numbers>
#include <cmath>

Domain makePlateWithHole(int circleSegments) {
    Domain d;

    d.points.push_back({0,0});
    d.points.push_back({10,0});
    d.points.push_back({10,10});
    d.points.push_back({0,10});

    d.segments.push_back({0,1});
    d.segments.push_back({1,2});
    d.segments.push_back({2,3});
    d.segments.push_back({3,0});

    for (int k = 0; k < circleSegments; k++) {
        double angle = 2 * std::numbers::pi * k / circleSegments;
        d.points.push_back({5 + 2 * std::cos(angle), 5 + 2 * std::sin(angle)});
    }

    for (int k = 0; k < circleSegments; k++) {
        d.segments.push_back({4 + k, 4 + (k + 1) % circleSegments});
    }

    d.holes.push_back({5, 5});

    return d;
}