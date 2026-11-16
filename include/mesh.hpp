#ifndef MESHSKIT
#define MESHSKIT

#include <vector>
#include <iostream>

struct Point {
    double x;
    double y;
};

struct Triangle {
    int v0;
    int v1;
    int v2;
};

class Mesh {
public:
    std::vector<Point> meshPoint;
    std::vector<Triangle> meshTriangle;

    std::size_t addPoint(Point point) {
        meshPoint.push_back(point);
        return meshPoint.size() - 1; // index
    }

    void addTriangle(Triangle triangle) {
        meshTriangle.push_back(triangle);
    }
};

#endif