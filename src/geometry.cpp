#include "geometry.hpp"
#include "predicate.hpp"

static void init() {
    static bool done = false;
    if (!done) { exactinit(); done = true; }
}

double orientation(const Point& a, const Point& b, const Point& c){ 
    init(); 

    double pa[2] = {a.x, a.y};
    double pb[2] = {b.x, b.y};
    double pc[2] = {c.x, c.y};

    return orient2d(pa, pb, pc);
}

bool inCircumcircle(const Point& a, const Point& b, const Point& c, const Point& d) {
    init(); 

    double pa[2] = {a.x, a.y};
    double pb[2] = {b.x, b.y};
    double pc[2] = {c.x, c.y};
    double pd[2] = {d.x, d.y};    

    double inc = incircle(pa, pb, pc, pd);
    double o = orientation(a, b, c);

    return (inc * o > 0);
}

bool encroaches(const Point& a, const Point& b, const Point& p) {
    return ((a.x - p.x)*(b.x - p.x) + (a.y - p.y)*(b.y - p.y) < 0);
}

Point circumcenter(const Point& a, const Point& b, const Point& c) { 
    double bx = b.x - a.x;
    double by = b.y - a.y;
    double cx = c.x - a.x;
    double cy = c.y - a.y;

    double d = 2 * (bx*cy - by*cx);

    double b2 = bx*bx + by*by;       // squared length of b
    double c2 = cx*cx + cy*cy;       // squared length of c

    double ux = (cy*b2 - by*c2) / d;
    double uy = (bx*c2 - cx*b2) / d;

    return {a.x + ux, a.y + uy};
}