// Unit tests for mesh I/O: Gmsh reader, VTK and Nastran writers, factories, config parser.
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "mesh_io.hpp"

namespace fs = std::filesystem;

static int failures = 0;

static void check(const std::string& name, bool ok) {
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name << "\n";
    if (!ok) ++failures;
}

static fs::path tmpDir() {
    fs::path d = fs::temp_directory_path() / "meshkit_test_io";
    fs::create_directories(d);
    return d;
}

static std::string writeText(const std::string& name, const std::string& text) {
    fs::path p = tmpDir() / name;
    std::ofstream(p) << text;
    return p.string();
}

static std::string readText(const std::string& path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

static bool contains(const std::string& text, const std::string& part) {
    return text.find(part) != std::string::npos;
}

// Runs f and reports whether it threw an exception of type E
template <typename E, typename F>
static bool throws(F f) {
    try { f(); } catch (const E&) { return true; } catch (...) { return false; }
    return false;
}

static bool samePoint(const Point& p, double x, double y) { return p.x == x && p.y == y; }
static bool sameTri(const Triangle& t, int a, int b, int c) { return t.v0 == a && t.v1 == b && t.v2 == c; }

int main() {
    std::cout << "--- Gmsh reader ---\n";
    Mesh sq = makeReader("gmsh")->read("tests/data/square.msh");
    check("square.msh has 4 points", sq.meshPoint.size() == 4);
    check("square.msh has 2 triangles (points/lines ignored)", sq.meshTriangle.size() == 2);
    check("point coordinates",
          sq.meshPoint.size() == 4 &&
          samePoint(sq.meshPoint[0], 0, 0) && samePoint(sq.meshPoint[1], 1, 0) &&
          samePoint(sq.meshPoint[2], 1, 1) && samePoint(sq.meshPoint[3], 0, 1));
    check("triangle indices are 0-based",
          sq.meshTriangle.size() == 2 &&
          sameTri(sq.meshTriangle[0], 0, 1, 2) && sameTri(sq.meshTriangle[1], 0, 2, 3));

    // Node ids that are not 1..N, and an element with 3 tags instead of 2
    std::string odd = writeText("odd_ids.msh",
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$Nodes\n3\n10 0 0 0\n20 2 0 0\n30 0 2 0\n$EndNodes\n"
        "$Elements\n2\n1 1 2 0 1 10 20\n2 2 3 7 8 9 30 10 20\n$EndElements\n");
    Mesh m = makeReader("gmsh")->read(odd);
    check("non-sequential node ids are mapped", m.meshPoint.size() == 3 &&
          m.meshTriangle.size() == 1 && sameTri(m.meshTriangle[0], 2, 0, 1));

    std::string v4 = writeText("v4.msh", "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");
    check("MSH version 4 is rejected",
          throws<std::runtime_error>([&] { makeReader("gmsh")->read(v4); }));
    check("missing file throws",
          throws<std::runtime_error>([] { makeReader("gmsh")->read("does/not/exist.msh"); }));

    std::cout << "--- factories ---\n";
    check("makeReader(\"banana\") throws invalid_argument",
          throws<std::invalid_argument>([] { makeReader("banana"); }));
    check("makeWriter(\"banana\") throws invalid_argument",
          throws<std::invalid_argument>([] { makeWriter("banana"); }));
    check("makeWriter(\"vtk\") and (\"nastran\") succeed",
          makeWriter("vtk") != nullptr && makeWriter("nastran") != nullptr);

    std::cout << "--- Nastran writer ---\n";
    std::string bdf = (tmpDir() / "square.bdf").string();
    makeWriter("nastran")->write(sq, bdf);
    check("square.bdf content is exact", readText(bdf) ==
          "BEGIN BULK\n"
          "GRID,1,,0,0,0.0\nGRID,2,,1,0,0.0\nGRID,3,,1,1,0.0\nGRID,4,,0,1,0.0\n"
          "CTRIA3,1,1,1,2,3\nCTRIA3,2,1,1,3,4\n"
          "ENDDATA\n");

    std::cout << "--- VTK writer ---\n";
    std::string vtkPath = (tmpDir() / "square.vtk").string();
    makeWriter("vtk")->write(sq, vtkPath);
    std::string vtk = readText(vtkPath);
    check("VTK header",           contains(vtk, "DATASET UNSTRUCTURED_GRID"));
    check("VTK has 4 points",     contains(vtk, "POINTS 4 double"));
    check("VTK has 2 cells",      contains(vtk, "CELLS 2 8"));
    check("VTK triangle 0",       contains(vtk, "3 0 1 2"));
    check("VTK triangle 1",       contains(vtk, "3 0 2 3"));
    check("VTK cell types",       contains(vtk, "CELL_TYPES 2\n5\n5"));

    std::cout << "--- config parser ---\n";
    std::string cfgPath = writeText("test.cfg",
        "# a comment\n\n"
        "reader = gmsh\n"
        "   writer=nastran   \n"
        "this line has no equals sign\n"
        "  # indented comment\n"
        "output = a=b.bdf\n");
    auto cfg = readConfig(cfgPath);
    check("3 keys parsed", cfg.size() == 3);
    check("spaces trimmed", cfg.count("reader") && cfg.at("reader") == "gmsh");
    check("no spaces around '='", cfg.count("writer") && cfg.at("writer") == "nastran");
    check("value keeps later '='", cfg.count("output") && cfg.at("output") == "a=b.bdf");
    check("missing config file throws",
          throws<std::runtime_error>([] { readConfig("does/not/exist.cfg"); }));

    std::cout << "\n" << (failures == 0 ? "ALL IO TESTS PASSED" : "SOME IO TESTS FAILED")
              << " (" << failures << " failure(s))\n";
    return failures == 0 ? 0 : 1;
}
