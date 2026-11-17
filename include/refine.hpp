#ifndef REFINE_HPP
#define REFINE_HPP

#include "domain.hpp"
#include "mesh.hpp"
#include <vector>

int findBadTriangle(const Mesh& mesh, const std::vector<bool>& inside, double maxRatio);
std::vector<Segment> splitEncroachedSegments(Mesh& mesh, std::vector<Segment> segments);
int splitSegment(Mesh& mesh, std::vector<Segment>& segments, std::size_t i);
int findEncroachedSegment(const Mesh& mesh, const std::vector<Segment>& segments, const Point& p);
void refineMesh(Mesh& mesh, std::vector<Segment>& segments,const std::vector<Point>& holes, double maxRatio);
Mesh meshDomain(const Domain& d, double maxRatio);

#endif