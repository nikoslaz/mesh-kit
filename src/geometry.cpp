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
