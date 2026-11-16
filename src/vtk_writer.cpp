#include "vtk_writer.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>

void writeVTK(const Mesh& mesh, const std::string& path) {
    std::ofstream file;
    file.open(path);
    if (!file) {
        std::cout << "File didn't opened, exiting!" << std::endl;
        exit(1);
    }
    std::size_t sizePoints = mesh.meshPoint.size(); 
    std::size_t sizeTriangle = mesh.meshTriangle.size();

    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "MESH SKIT OUTPUT" << std::endl;
    file << "ASCII"<< std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << sizePoints << " double" << std::endl;

    for (const auto& point: mesh.meshPoint) {
        file << point.x << " " << point.y << " 0" << std::endl;
    }

    file << "CELLS " << sizeTriangle << " " << sizeTriangle * 4 << std::endl;

    for (const auto& triang: mesh.meshTriangle) {
        file << "3 " << triang.v0 << " " << triang.v1 << " " <<triang.v2 << " " << std::endl;
    }

    file << "CELL_TYPES " << sizeTriangle << std::endl;

    for (std::size_t i = 0; i < sizeTriangle; i++) {
        file << "5" << std::endl;
    }

    file.close();
}

void writeVTK(const Mesh& mesh, const std::string& path, const std::string& dataName, const std::vector<double>& cellValues) {
    std::ofstream file;
    file.open(path);
    if (!file) {
        std::cout << "File didn't opened, exiting!" << std::endl;
        exit(1);
    }
    std::size_t sizePoints = mesh.meshPoint.size(); 
    std::size_t sizeTriangle = mesh.meshTriangle.size();

    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "MESH SKIT OUTPUT" << std::endl;
    file << "ASCII"<< std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << sizePoints << " double" << std::endl;

    for (const auto& point: mesh.meshPoint) {
        file << point.x << " " << point.y << " 0" << std::endl;
    }

    file << "CELLS " << sizeTriangle << " " << sizeTriangle * 4 << std::endl;

    for (const auto& triang: mesh.meshTriangle) {
        file << "3 " << triang.v0 << " " << triang.v1 << " " <<triang.v2 << " " << std::endl;
    }

    file << "CELL_TYPES " << sizeTriangle << std::endl;

    for (std::size_t i = 0; i < sizeTriangle; i++) {
        file << "5" << std::endl;
    }

    file << "CELL_DATA " << mesh.meshTriangle.size() << std::endl;
    file << "SCALARS " << dataName << " double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;

    for (double v : cellValues) {
        file << v << '\n';
    }


    file.close();
}
