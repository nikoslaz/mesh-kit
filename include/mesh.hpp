#ifndef MESHSKIT
#define MESHSKIT

#include <vector>
#include <iostream>

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

struct Point {
    double x;
    double y;
};

struct Triangle {
    int v0;
    int v1;
    int v2;
    int n0 = -1;
    int n1 = -1;
    int n2 = -1;
};

class Mesh {
public:
    std::vector<Point> meshPoint;
    std::vector<Triangle> meshTriangle;
    int lastTriangle = 0;

    std::size_t addPoint(Point point) {
        meshPoint.push_back(point);
        return meshPoint.size() - 1; // index
    }

    void addTriangle(Triangle triangle) {
        meshTriangle.push_back(triangle);
    }
};

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

struct Point3 {
    double x;
    double y;
    double z;
};

struct Tet {
    int v0;
    int v1;
    int v2;
    int v3;
    int n0 = -1;
    int n1 = -1;
    int n2 = -1;
    int n3 = -1;
};

class TetMesh {
public:
    std::vector<Point3> points;
    std::vector<Tet> tets;
    int lastTet = 0;

    std::size_t addPoint(Point3 point) {
        points.push_back(point);
        return points.size() - 1; // index
    }

    void addTet(Tet tet) {
        tets.push_back(tet);
    }
};

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif