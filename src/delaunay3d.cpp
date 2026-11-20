#include "delaunay.hpp"
#include "topology.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <queue>
#include <utility>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <iostream>

struct CaviteFace {
    int a;
    int b;
    int c;
    int outside;
    int owner;     // the cavity tetra this face belongs to
    Tet ownerTet;  // copy of the owner (its slot may be reused before we build the new tets)
    int slot;      // the face is opposite ownerTet's vertex in this slot (0..3)
};

static bool inConflict(const TetMesh& mesh, int t, const Point3& p);   // defined with the ghost code below

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

static std::vector<int> findCavity(const TetMesh& mesh, int start, const Point3& p) {
    std::vector<int> cavity;
    std::vector<bool> visited(mesh.tets.size(), false);
    std::queue<int> q;

    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        int t = q.front();
        q.pop();
        cavity.push_back(t);

        const Tet& tet = mesh.tets[t];
        int nbrs[4] = {tet.n0, tet.n1, tet.n2, tet.n3};
        for (int k : nbrs) {
            if (k == -1 || visited[k]) continue;
            if (inConflict(mesh, k, p)) {
                q.push(k);
                visited[k] = true;
            }
        }
    }
    return cavity;
}

static std::vector<CaviteFace> cavityBoundary(const TetMesh& mesh, const std::vector<int>& cavity) {
    std::vector<CaviteFace> caviteFace;
    std::vector<bool> inCavity(mesh.tets.size(), false);

    for (const auto& c : cavity) {
        inCavity[c] = true;
    }
    
    for (const auto& c : cavity) {
        const Tet& t = mesh.tets[c];

        if (t.n0 == -1 || !inCavity[t.n0]) {
            CaviteFace ce = {t.v1, t.v2, t.v3, t.n0, c, t, 0};
            caviteFace.push_back(ce);
        }
        if (t.n1 == -1 || !inCavity[t.n1]) {
            CaviteFace ce = {t.v0, t.v2, t.v3, t.n1, c, t, 1};
            caviteFace.push_back(ce);
        }
        if (t.n2 == -1 || !inCavity[t.n2]) {
            CaviteFace ce = {t.v0, t.v1, t.v3, t.n2, c, t, 2};
            caviteFace.push_back(ce);
        }
        if (t.n3 == -1 || !inCavity[t.n3]) {
            CaviteFace ce = {t.v0, t.v1, t.v2, t.n3, c, t, 3};
            caviteFace.push_back(ce);
        }
    }

    return caviteFace;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

static std::vector<int> spatialOrder(const std::vector<Point3>& points) {
    if (points.empty()) return {};
    double minZ = points[0].z;
    double maxZ = points[0].z;
    double minY = points[0].y;
    double maxY = points[0].y;

    for (const auto& p : points) {
        minZ = std::min(minZ, p.z);
        maxZ = std::max(maxZ, p.z);
        minY = std::min(minY, p.y);
        maxY = std::max(maxY, p.y);
    }

    const int layers = std::max(1, static_cast<int>(std::cbrt(points.size() / 2.0)));
    const int rows   = layers;

    const double depth = maxZ - minZ;
    auto layerOf = [&](const Point3& p) {
        if (depth <= 0) return 0;
        return std::min(layers - 1, static_cast<int>((p.z - minZ) / depth * layers));
    };

    const double height = maxY - minY;
    auto rowOf = [&](const Point3& p) {
        if (height <= 0) return 0;
        return std::min(rows - 1, static_cast<int>((p.y - minY) / height * rows));
    };

    std::vector<int> order(points.size());
    for (std::size_t i = 0; i < order.size(); i++) order[i] = static_cast<int>(i);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
        int li = layerOf(points[i]), lj = layerOf(points[j]);
        if (li != lj) return li < lj;

        int ri = rowOf(points[i]), rj = rowOf(points[j]);
        if (ri != rj) return (li % 2 == 0) ? ri < rj : ri > rj;   // flip row direction on odd layers

        // same layer, same row: sort by x
        // the row's x-direction should flip each time you move to the next row.
        // "next row" parity depends on both layer and row: use (li + ri) % 2
        return ((li + ri) % 2 == 0) ? points[i].x < points[j].x : points[i].x > points[j].x;
    });
    
    return order;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

static bool collinear3(const Point3& a, const Point3& b, const Point3& c) {
    // project onto the three coordinate planes -> collinear in 3D <=> collinear in all three
    Point axy{a.x, a.y}, bxy{b.x, b.y}, cxy{c.x, c.y};
    Point ayz{a.y, a.z}, byz{b.y, b.z}, cyz{c.y, c.z};
    Point azx{a.z, a.x}, bzx{b.z, b.x}, czx{c.z, c.x};
    return orientation(axy, bxy, cxy) == 0 && orientation(ayz, byz, cyz) == 0 && orientation(azx, bzx, czx) == 0;
}

// positions in `order` of the first 4 non-coplanar points; all -1 if every point is coplanar
static std::array<int,4> findStartTet(const std::vector<Point3>& pts, const std::vector<int>& order) {
    const int n = static_cast<int>(order.size());
    std::array<int,4> none = {-1, -1, -1, -1};
    if (n < 4) return none;

    const Point3& a = pts[order[0]];

    // 1) first point different from a
    int i1 = -1;
    for (int k = 1; k < n; ++k) {
        const Point3& q = pts[order[k]];
        if (q.x != a.x || q.y != a.y || q.z != a.z) { 
            i1 = k; 
            break; 
        }
    }
    if (i1 == -1) return none;
    const Point3& b = pts[order[i1]];

    // 2) first point not collinear with a, b
    int i2 = -1;
    for (int k = i1 + 1; k < n; ++k) {
        if (!collinear3(a, b, pts[order[k]])) {
            i2 = k;
            break;
        }
    }
    if (i2 == -1) return none;
    const Point3& c = pts[order[i2]];

    // 3) first point not coplanar with a, b, c: use orientation(a, b, c, d) != 0
    int i3 = -1;
    for (int k = i2 + 1; k < n; ++k) {
        if (orientation(a, b, c, pts[order[k]]) != 0) {
            i3 = k;
            break;
        }
    }
    if (i3 == -1) return none;

    return {0, i1, i2, i3};
}

static constexpr int GHOST = 0;   // point index reserved for the vertex at infinity

// slot (0..3) of GHOST in t, or -1 if t is a real tet
static int ghostSlot(const Tet& t) {
    if (t.v0 == GHOST) return 0;
    if (t.v1 == GHOST) return 1;
    if (t.v2 == GHOST) return 2;
    if (t.v3 == GHOST) return 3;
    return -1;
}

// 1 real tet + 4 ghost tets, all neighbours set. points = [dummy, p0, p1, p2, p3]
static TetMesh makeGhostMesh(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3) {
    TetMesh mesh;

    mesh.addPoint({0, 0, 0});   // index 0 = GHOST
    mesh.addPoint(p0);          // indices 1..4
    mesh.addPoint(p1);
    mesh.addPoint(p2);
    mesh.addPoint(p3);

    if (orientation(p0, p1, p2, p3) >= 0) {
        Tet tet = {1, 2, 3, 4};
        mesh.addTet(tet);
    } 
    else {
        Tet tet = {1, 3, 2, 4};
        mesh.addTet(tet);
    }

    const Tet real = mesh.tets[0];      
    for (int i = 0; i < 4; i++) {
        Tet g = real;
        int* vs[4] = {&g.v0, &g.v1, &g.v2, &g.v3};
        *vs[i] = GHOST;   // the ghost sits outside the face opposite vertex i
        // GHOST is on the other side of that face than v_i was, so swap two of the other
        // vertices: then every tet is positive if GHOST is thought of as a point far outside
        if (i < 2) std::swap(g.v2, g.v3);
        else       std::swap(g.v0, g.v1);
        mesh.addTet(g);
    }

    buildNeighbours(mesh);
    mesh.lastTet = 0;

    return mesh;
}

// true if p lies inside the circumsphere of tet t (for a ghost tet: p "sees" its hull face)
static bool inConflict(const TetMesh& mesh, int t, const Point3& p) {
    const Tet& tet = mesh.tets[t];
    int vs[4] = {tet.v0, tet.v1, tet.v2, tet.v3};
    int ns[4] = {tet.n0, tet.n1, tet.n2, tet.n3};

    int g = ghostSlot(tet);   // -1 = not a ghost tet

    // real tet -> same vertex order as findCavity
    if (g == -1) {
        return inCircumcircle(mesh.points[vs[0]], mesh.points[vs[2]], mesh.points[vs[1]], mesh.points[vs[3]], p);
    }

    // ghost tet -> put p in the ghost's slot
    Point3 q[4];
    for (int i = 0; i < 4; i++){
        q[i] = (i == g) ? p : mesh.points[vs[i]];
    }
    double o = orientation(q[0], q[1], q[2], q[3]);

    if (o > 0) return true;
    if (o < 0) return false;
    return inConflict(mesh, ns[g], p);   
}

// Walk from `start` towards p. Returns a real tet containing p, or (p outside the hull) the ghost tet
// over the hull face that p sees. Either way the result is in conflict with p and seeds the cavity.
// Never calls orientation on a ghost tet: points[GHOST] is only a placeholder.
static int walkToConflict(const TetMesh& mesh, int start, const Point3& p) {
    int t = start;
    int g = ghostSlot(mesh.tets[t]);
    if (g != -1) {   // start from the real tet behind the hull face
        const Tet& s = mesh.tets[t];
        int ns[4] = {s.n0, s.n1, s.n2, s.n3};
        t = ns[g];
    }

    int limit = static_cast<int>(mesh.tets.size()) + 10;
    for (int steps = 0; steps < limit; steps++) {
        const Tet& tet = mesh.tets[t];
        const Point3& a = mesh.points[tet.v0];
        const Point3& b = mesh.points[tet.v1];
        const Point3& c = mesh.points[tet.v2];
        const Point3& d = mesh.points[tet.v3];

        int next;
        if      (orientation(p, b, c, d) < 0) next = tet.n0;   // beyond the face opposite a
        else if (orientation(a, p, c, d) < 0) next = tet.n1;   // beyond the face opposite b
        else if (orientation(a, b, p, d) < 0) next = tet.n2;   // beyond the face opposite c
        else if (orientation(a, b, c, p) < 0) next = tet.n3;   // beyond the face opposite d
        else return t;                                         // inside (or on a face)

        if (ghostSlot(mesh.tets[next]) != -1) return next;     // p is outside the hull, beyond this face
        t = next;
    }
    return -1;   // only if the neighbours are broken
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

int insertPoint3d(TetMesh& mesh, Point3 p) {
    int hint = (mesh.lastTet >= 0 && mesh.lastTet < static_cast<int>(mesh.tets.size())) ? mesh.lastTet : 0;
    int start = walkToConflict(mesh, hint, p);
    if (start == -1) return -1;
    int idx = mesh.addPoint(p);

    const Tet& s = mesh.tets[start];
    int vs[4] = {s.v0, s.v1, s.v2, s.v3};
    for (int v : vs) {
        if (v == GHOST) continue;   // placeholder point, not a real vertex
        const Point3& q = mesh.points[v];
        if (q.x == p.x && q.y == p.y && q.z == p.z) return idx;
    }

    std::vector<int> cavity = findCavity(mesh, start, p);
    std::vector<CaviteFace> caviteFace = cavityBoundary(mesh,cavity);
    std::vector<int> newTets;

    for (size_t i = 0; i < caviteFace.size(); i++) {
        // the owner tet with the vertex opposite this face replaced by p: p lies on the same side
        // of the face as that vertex did, so the new tet keeps the owner's (positive) orientation.
        // No orientation call, so this also works for faces that contain GHOST.
        Tet tet = caviteFace[i].ownerTet;
        int* tv[4] = {&tet.v0, &tet.v1, &tet.v2, &tet.v3};
        *tv[caviteFace[i].slot] = idx;
        tet.n0 = tet.n1 = tet.n2 = tet.n3 = -1;

        int slot;
        if (i < cavity.size()) { 
            slot = cavity[i]; 
            mesh.tets[slot] = tet; 
        }
        else { 
            mesh.addTet(tet); 
            slot = (int)mesh.tets.size() - 1; 
        }
        newTets.push_back(slot);
    }

    for (size_t i = 0; i < newTets.size(); i++) {
        if (caviteFace[i].outside != -1) {
            Tet& outsideTet = mesh.tets[caviteFace[i].outside];
            int a = caviteFace[i].a;
            int b = caviteFace[i].b;
            int c = caviteFace[i].c;
            setNeighbour(outsideTet, a, b, c, newTets[i]);
            setNeighbour(mesh.tets[newTets[i]], a, b, c, caviteFace[i].outside);
        }
    }

    std::map<std::pair<int,int>, std::vector<int>> edgeTets;
    for (size_t i = 0; i < newTets.size(); i++) {
        std::pair<int,int> key;
        int a = caviteFace[i].a;
        int b = caviteFace[i].b;
        int c = caviteFace[i].c;

        // edge 1: between a and b
        key = { std::min(a, b), std::max(a, b) };   
        edgeTets[key].push_back((newTets[i]));

        // edge 2: between b and c
        key = { std::min(b, c), std::max(b, c) };   
        edgeTets[key].push_back((newTets[i]));

        // edge 3: between c and a
        key = { std::min(c, a), std::max(c, a) };
        edgeTets[key].push_back((newTets[i]));
    }

    for (const auto& e : edgeTets) {
        assert(e.second.size() == 2);
        int x = e.first.first, y = e.first.second;
        int t1 = e.second[0], t2 = e.second[1];
        setNeighbour(mesh.tets[t1], x, y, idx, t2);
        setNeighbour(mesh.tets[t2], x, y, idx, t1);
    }

    // 3D can create fewer tets than the cavity held so we drop the unused cavity slots
    if (cavity.size() > newTets.size()) {
        std::vector<int> spare(cavity.begin() + newTets.size(), cavity.end());
        std::sort(spare.begin(), spare.end(), std::greater<int>());
        for (int s : spare) {
            int last = static_cast<int>(mesh.tets.size()) - 1;
            if (s != last) {
                mesh.tets[s] = mesh.tets[last];
                Tet& moved = mesh.tets[s];
                int* ns[4] = {&moved.n0, &moved.n1, &moved.n2, &moved.n3};
                for (int* n : ns) {
                    if (*n == -1) continue;
                    Tet& o = mesh.tets[*n];
                    int* os[4] = {&o.n0, &o.n1, &o.n2, &o.n3};
                    for (int* q : os) if (*q == last) *q = s;
                }
            }
            mesh.tets.pop_back();
        }
    }

    mesh.lastTet = static_cast<int>(mesh.tets.size()) - 1;   // a recent tet, valid after the compaction

    return idx;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

TetMesh delaunay3d(const std::vector<Point3>& points) {
    if (points.empty()) return TetMesh{};
    std::vector<int> order = spatialOrder(points);

    // the first 4 non-coplanar points (in insertion order) form the starting tet
    std::array<int,4> st = findStartTet(points, order);
    if (st[0] == -1) {   // every point is coplanar: no tets
        TetMesh flat;
        flat.points = points;
        return flat;
    }

    // insertion sequence: the 4 start points first, then the rest in spatial order
    std::vector<int> seq;
    seq.reserve(order.size());
    for (int k : st) seq.push_back(order[k]);
    for (int k = 0; k < static_cast<int>(order.size()); k++)
        if (k != st[0] && k != st[1] && k != st[2] && k != st[3]) seq.push_back(order[k]);

    // points[0] = GHOST, points[1..4] = start points, then one point per insertion
    TetMesh mesh = makeGhostMesh(points[seq[0]], points[seq[1]], points[seq[2]], points[seq[3]]);
    for (std::size_t k = 4; k < seq.size(); k++) {
        insertPoint3d(mesh, points[seq[k]]);
    }

    // drop the ghost tets and renumber: mesh vertex v (>= 1) is input point seq[v - 1]
    std::vector<Tet> real;
    for (const auto& t : mesh.tets) {
        if (ghostSlot(t) != -1) continue;
        Tet r;
        r.v0 = seq[t.v0 - 1];
        r.v1 = seq[t.v1 - 1];
        r.v2 = seq[t.v2 - 1];
        r.v3 = seq[t.v3 - 1];
        real.push_back(r);
    }

    TetMesh result;
    result.points = points;   // keep the input numbering
    result.tets = real;
    result.lastTet = 0;
    buildNeighbours(result);  // hull faces get -1 again
    return result;
}
