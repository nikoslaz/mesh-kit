#include "mesh_io.hpp"
#include <fstream>
#include <map>
#include <sstream>

void VtkWriter::write(const Mesh& mesh, const std::string& path) const {
    writeVTK(mesh, path);
}

std::unique_ptr<MeshWriter> makeWriter(const std::string& format) {
    if (format == "vtk") {
        return std::make_unique<VtkWriter>();
    }
    if (format == "nastran") {
        return std::make_unique<NastranWriter>();
    }
    throw std::invalid_argument("unknown mesh writer format: " + format);
}

Mesh GmshReader::read(const std::string& path) const {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open " + path);

    Mesh mesh;
    std::map<int,int> nodeIndex;         
    std::string line;

    while (std::getline(file, line)) {
        if (line == "$MeshFormat") {
            std::getline(file, line);
            if (line.rfind("2.2", 0) != 0) {
                throw std::runtime_error("Wrong versions, must be 2.2! " + path);
            }
        } else if (line == "$Nodes") {
            std::getline(file, line);
            int n = std::stoi(line);

            for (int k = 0; k < n; k++) {
                std::getline(file, line);
                std::istringstream ss(line);
                int id; double x, y, z;
                ss >> id >> x >> y >> z;

                std::size_t index = mesh.addPoint({x, y});
                nodeIndex[id] = static_cast<int>(index);
            }
        } else if (line == "$Elements") {
            std::getline(file, line);
            int m = std::stoi(line);

            for (int k = 0; k < m; k++) {
                std::getline(file, line);
                std::istringstream ss(line);
                int id, type, ntags;
                ss >> id >> type >> ntags;
                if (type != 2) continue;          
                int tag;
                for (int t = 0; t < ntags; t++) ss >> tag;   

                int a, b, c;
                ss >> a >> b >> c;                 
                mesh.addTriangle({nodeIndex.at(a), nodeIndex.at(b), nodeIndex.at(c)});
            }
        }
    }
    return mesh;
}

std::unique_ptr<MeshReader> makeReader(const std::string& format) {
    if (format == "gmsh") {
        return std::make_unique<GmshReader>();
    }
    throw std::invalid_argument("unknown mesh reader format: " + format);
}

void NastranWriter::write(const Mesh& mesh, const std::string& path) const {
    std::ofstream file(path);
    if (!file) throw std::runtime_error("cannot open " + path);

    file << "BEGIN BULK\n";

    for (std::size_t i = 0; i < mesh.meshPoint.size(); i++) {
        const Point& p = mesh.meshPoint[i];
        file << "GRID," << i + 1 << ",," << p.x << "," << p.y << ",0.0\n";
    }
    for (std::size_t i = 0; i < mesh.meshTriangle.size(); i++) {
        const Triangle& t = mesh.meshTriangle[i];
        file << "CTRIA3," << i + 1 << ",1," << t.v0 + 1 << "," << t.v1 + 1 << "," << t.v2 + 1 << "\n";
    }

    file << "ENDDATA\n";
}

static std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t\r");
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(" \t\r");
    return s.substr(a, b - a + 1);
}

std::map<std::string, std::string> readConfig(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open " + path);

    std::map<std::string, std::string> cfg;
    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;      

        std::size_t pos = line.find('=');
        if (pos == std::string::npos) continue;            

        std::string key   = trim(line.substr(0, pos));
        std::string value = trim(line.substr(pos + 1));
        cfg[key] = value;
    }
    return cfg;
}