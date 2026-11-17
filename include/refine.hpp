#ifndef REFINE_HPP
#define REFINE_HPP

#include "domain.hpp"
#include "mesh.hpp"
#include <vector>

// Index of the first inside triangle whose radius-edge ratio exceeds maxRatio, or -1 if none.
int findBadTriangle(const Mesh& mesh, const std::vector<bool>& inside, double maxRatio);

std::vector<Segment> splitEncroachedSegments(Mesh& mesh, std::vector<Segment> segments);
int splitSegment(Mesh& mesh, std::vector<Segment>& segments, std::size_t i);
int findEncroachedSegment(const Mesh& mesh, const std::vector<Segment>& segments, const Point& p);

#endif