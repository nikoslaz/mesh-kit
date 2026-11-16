// meshkit_cli
//   meshkit_cli <config.cfg>   read -> (smooth) -> write, as described by a key=value config file
//   meshkit_cli demo           random 60-point example, writes before/after VTK files to output/
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "delaunay.hpp"
#include "mesh_io.hpp"
#include "quality.hpp"
#include "smoothing.hpp"
#include "vtk_writer.hpp"

using Config = std::map<std::string, std::string>;

// Value of a required key, with a clear error instead of map::at's
static std::string need(const Config& cfg, const std::string& key) {
    auto it = cfg.find(key);
    if (it == cfg.end()) throw std::runtime_error("missing config key: " + key);
    return it->second;
}

// Prints "<label>_min X" and "<label>_mean Y" (easy for scripts to parse)
static void printStats(const std::string& label, const std::vector<double>& angles) {
    if (angles.empty()) return;
    double mn = *std::min_element(angles.begin(), angles.end());
    double mean = std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
    std::cout << label << "_min " << mn << "\n" << label << "_mean " << mean << "\n";
}

static void runConfig(const std::string& cfgPath) {
    Config cfg = readConfig(cfgPath);

    auto reader = makeReader(need(cfg, "reader"));
    auto writer = makeWriter(need(cfg, "writer"));
    Mesh mesh = reader->read(need(cfg, "input"));
    std::cout << "points " << mesh.meshPoint.size() << "\n"
              << "triangles " << mesh.meshTriangle.size() << "\n";

    std::vector<double> before = minAngles(mesh);

    int iterations = cfg.count("smooth") ? std::stoi(cfg.at("smooth")) : 0;
    if (iterations > 0) laplacianSmooth(mesh, iterations);
    std::vector<double> after = minAngles(mesh);

    printStats("min_angle_before", before);
    printStats("min_angle_after", after);

    std::filesystem::path out = need(cfg, "output");
    if (out.has_parent_path()) std::filesystem::create_directories(out.parent_path());
    writer->write(mesh, out.string());

    // Optional: per-triangle quality as CSV
    if (cfg.count("quality")) {
        std::ofstream q(cfg.at("quality"));
        if (!q) throw std::runtime_error("cannot open " + cfg.at("quality"));
        q << "triangle,min_angle_before,min_angle_after\n";
        for (std::size_t i = 0; i < before.size(); i++)
            q << i << "," << before[i] << "," << after[i] << "\n";
    }
}

static void demo() {
    std::filesystem::create_directories("output");
    std::mt19937 rng(42);
    std::vector<Point> pts;
    for (int i = 0; i < 60; i++)
        pts.push_back({rng() / 4294967296.0, rng() / 4294967296.0});

    Mesh mesh = delaunay(pts);
    auto before = minAngles(mesh);
    std::cout << "smallest angle before: " << *std::min_element(before.begin(), before.end()) << "\n";
    writeVTK(mesh, "output/random60_before.vtk", "min_angle", before);

    laplacianSmooth(mesh, 10);
    auto after = minAngles(mesh);
    std::cout << "smallest angle after:  " << *std::min_element(after.begin(), after.end()) << "\n";
    writeVTK(mesh, "output/random60_after.vtk", "min_angle", after);
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: meshkit_cli <config.cfg>\n"
                  << "       meshkit_cli demo\n";
        return 1;
    }
    try {
        std::string arg = argv[1];
        if (arg == "demo") demo();
        else runConfig(arg);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
