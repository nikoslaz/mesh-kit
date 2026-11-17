#ifndef REFINE_HPP
#define REFINE_HPP

#include "domain.hpp"
#include "mesh.hpp"
#include <vector>

std::vector<Segment> splitEncroachedSegments(Mesh& mesh, std::vector<Segment> segments);

#endif