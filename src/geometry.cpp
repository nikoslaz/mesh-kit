#include "geometry.hpp"

double orientation(const Point& a, const Point& b, const Point& c){ 
    double ux = (b.x - a.x);
    double uy = (b.y - a.y);
    double vx = (c.x - a.x);
    double vy = (c.y - a.y);
    return ux * vy - uy * vx;
}

bool inCircumcircle(const Point& a, const Point& b, const Point& c, const Point& d) {
    double p1 = a.x - d.x;
    double p2 = a.y - d.y;    
    double p3 = p1 * p1 + p2 * p2;
    double q1 = b.x - d.x;
    double q2 = b.y - d.y;    
    double q3 = q1 * q1 + q2 * q2;
    double r1 = c.x - d.x;
    double r2 = c.y - d.y;
    double r3 = r1 * r1 + r2 * r2;

    double out = p1*(q2*r3 - q3*r2)  -  p2*(q1*r3 - q3*r1)  +  p3*(q1*r2 - q2*r1);
    double o = orientation(a, b, c);

    if (out * o > 0) {
        return true;
    }
    else {
        return false;
    }
}
