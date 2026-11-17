#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "mesh.hpp"

struct Segment { 
    int a; 
    int b; 
};  // indices into points

struct Domain {
    std::vector<Point> points;
    std::vector<Segment> segments;
    std::vector<Point> holes;       // one point inside each hole
};

Domain makePlateWithHole(int circleSegments);

#endif