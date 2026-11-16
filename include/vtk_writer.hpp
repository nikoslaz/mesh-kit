#ifndef VTK_WRITER
#define VTK_WRITER

#include "mesh.hpp"
#include "quality.hpp"
#include <string>

void writeVTK(const Mesh& mesh, const std::string& path);
void writeVTK(const Mesh& mesh, const std::string& path,
              const std::string& dataName, const std::vector<double>& cellValues);

#endif