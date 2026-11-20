#ifndef QUALITY
#define QUALITY

#include "mesh.hpp"
#include "geometry.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

double angleAt(const Point& corner, const Point& p, const Point& q);
double minAngle(const Point& a, const Point& b, const Point& c);
std::vector<double> minAngles(const Mesh& mesh);

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

double radiusEdgeRatio(const Point& a, const Point& b, const Point& c); // 2D
double tetRadiusRatio(const Point3& a, const Point3& b, const Point3& c, const Point3& d); // 3D, 3r/R in [0,1]
std::vector<double> tetQualities(const TetMesh& mesh);                                     // one radius ratio per tet

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

std::vector<int> findSlivers(const TetMesh& mesh, double threshold);
double worstQuality(const TetMesh& mesh, const Tet* tets, int count); // smallest radius ratio among count tets (1.0 if count is 0)
int removeSlivers(TetMesh& mesh, double threshold, int maxPasses);
int improveSlivers(TetMesh& mesh, double threshold, int rounds); // removeSlivers + smoothing, repeated

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

#endif