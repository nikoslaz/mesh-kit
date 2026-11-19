#include "smoothing.hpp"

// 3D guarded Laplacian smoothing (see smoothing.hpp).
int laplacianSmooth3d(TetMesh& mesh, int iterations) {
    // moving vertices does not change the connectivity, so these are built once
    auto vt   = vertexTets(mesh);
    auto nb   = vertexNeighbours(mesh, vt);
    auto hull = hullVertices(mesh);

    int count = 0;
    for (int iter = 0; iter < iterations; iter++) {
        for (std::size_t v = 0; v < vt.size(); v++) {
            if (hull[v])       continue;
            if (nb[v].empty()) continue;

            double sx = 0, sy = 0, sz = 0;
            for (int w : nb[v]) {
                sx += mesh.points[w].x;
                sy += mesh.points[w].y;
                sz += mesh.points[w].z;
            }
            const double n = static_cast<double>(nb[v].size());
            Point3 target{sx / n, sy / n, sz / n};

            if (tryMoveVertex(mesh, static_cast<int>(v), target, vt[v])) count++;
        }
    }
    return count;
}
