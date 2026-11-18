cmake_minimum_required(VERSION 3.16)
project(meshkit LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Library: the actual mesh code
add_library(meshkit
    src/mesh.cpp
    src/vtk_writer.cpp
    src/geometry.cpp
    src/delaunay.cpp
    src/quality.cpp
    src/smoothing.cpp
    src/mesh_io.cpp
    src/predicate.c
    src/domain.cpp
    src/refine.cpp
    src/topology.cpp
    src/delaunay3d.cpp
)
# Shewchuk's predicates need every floating-point operation rounded separately: no fused
# multiply-add (-ffp-contract=off) and never -ffast-math. -w silences his old-style C.
set_source_files_properties(src/predicate.c PROPERTIES COMPILE_OPTIONS "-w;-ffp-contract=off")
target_include_directories(meshkit PUBLIC include)
target_compile_options(meshkit PRIVATE -Wall -Wextra -Wpedantic)

# Command-line program
add_executable(meshkit_cli apps/main.cpp)
target_link_libraries(meshkit_cli PRIVATE meshkit)
target_compile_options(meshkit_cli PRIVATE -Wall -Wextra -Wpedantic)

# Benchmark (use -DCMAKE_BUILD_TYPE=Release)
add_executable(meshkit_bench apps/bench.cpp)
target_link_libraries(meshkit_bench PRIVATE meshkit)

# Tests:  cmake --build build && ctest --test-dir build
enable_testing()
foreach(t test_geometry test_delaunay test_io test_domain test_refine test_topology)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE meshkit)
    add_test(NAME ${t} COMMAND ${t} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
