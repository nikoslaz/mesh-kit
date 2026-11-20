#include "geometry.hpp"
#include "predicate.hpp"
#include "topology.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <vector>

static void init() {
    static bool done = false;
    if (!done) { exactinit(); done = true; }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

double orientation(const Point& a, const Point& b, const Point& c){ 
    init(); 

    double pa[2] = {a.x, a.y};
    double pb[2] = {b.x, b.y};
    double pc[2] = {c.x, c.y};

    return orient2d(pa, pb, pc);
}

double orientation(const Point3& a, const Point3& b, const Point3& c, const Point3& d) {
    init(); 

    double pa[3] = {a.x, a.y, a.z};
    double pb[3] = {b.x, b.y, b.z};
    double pc[3] = {c.x, c.y, c.z};
    double pd[3] = {d.x, d.y, d.z};

    return orient3d(pa, pb, pc, pd);
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

bool encroaches(const Point& a, const Point& b, const Point& p) {
    return ((a.x - p.x)*(b.x - p.x) + (a.y - p.y)*(b.y - p.y) < 0);
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

bool inCircumcircle(const Point& a, const Point& b, const Point& c, const Point& d) {
    init(); 

    double pa[2] = {a.x, a.y};
    double pb[2] = {b.x, b.y};
    double pc[2] = {c.x, c.y};
    double pd[2] = {d.x, d.y};    

    double inc = incircle(pa, pb, pc, pd);
    double o = orientation(a, b, c);

    return (inc * o > 0);
}

bool inCircumcircle(const Point3& a, const Point3& b, const Point3& c, const Point3& d, const Point3& e) {
    init(); 

    double pa[3] = {a.x, a.y, a.z};
    double pb[3] = {b.x, b.y, b.z};
    double pc[3] = {c.x, c.y, c.z};
    double pd[3] = {d.x, d.y, d.z};
    double pe[3] = {e.x, e.y, e.z};      

    double inc = insphere(pa, pb, pc, pd, pe);
    double o = orientation(a, b, c, d);

    return (inc * o > 0);
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

Point circumcenter(const Point& a, const Point& b, const Point& c) { 
    double bx = b.x - a.x;
    double by = b.y - a.y;
    double cx = c.x - a.x;
    double cy = c.y - a.y;

    double d = 2 * (bx*cy - by*cx);

    double b2 = bx*bx + by*by;       // squared length of b
    double c2 = cx*cx + cy*cy;       // squared length of c

    double ux = (cy*b2 - by*c2) / d;
    double uy = (bx*c2 - cx*b2) / d;

    return {a.x + ux, a.y + uy};
}

// ------------------------------------------------------------------------------------------------------------------------------------------------ //

// t1 and t2 share a face. Fills out[3] with the three tets around edge (apex1, apex2).
// Returns false if any of them is not positive (edge de doesn't pass through face abc).
bool flip23Tets(const TetMesh& mesh, int t1, int t2, Tet out[3]) {
    Tet tet1 = mesh.tets[t1];
    Tet tet2 = mesh.tets[t2];

    int A[4]  = {tet1.v0, tet1.v1, tet1.v2, tet1.v3};
    int B[4]  = {tet2.v0, tet2.v1, tet2.v2, tet2.v3};
    int nA[4] = {tet1.n0, tet1.n1, tet1.n2, tet1.n3};

    bool found = false;
    int d = -1;
    int F[3];
    int n = 0;
    for (std::size_t i = 0; i < 4; i++) {
        if (nA[i] == t2) {
            found = true;
            d = A[i];
            for (std::size_t k = 0; k < 4; k++) {
                if (k != i) F[n++] = A[k];
            }
        }
    }
    if (found == false) return false;

    int e = -1;
    for (int i = 0; i < 4; i++) {
        if (B[i] != F[0] && B[i] != F[1] && B[i] != F[2]){
            e = B[i];
        } 
    }
    if (e == -1) return false;

    double o = orientation(mesh.points[F[0]], mesh.points[F[1]], mesh.points[F[2]], mesh.points[d]);
    if (o == 0) return false;
    if (o < 0)  std::swap(F[0], F[1]);

    int a = F[0];
    int b = F[1];
    int c = F[2];
    out[0] = Tet{a, b, e, d};
    out[1] = Tet{b, c, e, d};
    out[2] = Tet{c, a, e, d};

    for (std::size_t k = 0; k < 3; k++) {
        Point3 pa = mesh.points[out[k].v0];
        Point3 pb = mesh.points[out[k].v1];
        Point3 pc = mesh.points[out[k].v2];
        Point3 pd = mesh.points[out[k].v3];

        if (orientation(pa, pb, pc, pd) <= 0) {
            return false;
        }
    }

    return true;
}

// Neighbour slot i of t (the face opposite vertex i).
static int& nbr(Tet& t, int i) {
    switch (i) {
        case 0:  return t.n0;
        case 1:  return t.n1;
        case 2:  return t.n2;
        default: return t.n3;
    }
}

// The three vertices of t opposite slot i, sorted, so equal faces compare equal.
static std::array<int,3> faceKey(const Tet& t, int i) {
    int v[4] = {t.v0, t.v1, t.v2, t.v3};
    std::array<int,3> f{};
    int k = 0;
    for (int j = 0; j < 4; j++)
        if (j != i) f[k++] = v[j];
    std::sort(f.begin(), f.end());
    return f;
}

// The neighbours of the tets in `cluster` that are not in it (what a flip must re-attach to).
static std::vector<int> outerNeighbours(const TetMesh& mesh, const std::vector<int>& cluster) {
    std::vector<int> outer;
    for (int c : cluster) {
        const Tet& t = mesh.tets[c];
        for (int n : {t.n0, t.n1, t.n2, t.n3})
            if (n != -1 && std::find(cluster.begin(), cluster.end(), n) == cluster.end()) outer.push_back(n);
    }
    return outer;
}

// After a flip wrote new tets into `fresh`, links their faces to each other and to `outer`.
// Only these few tets are touched, instead of rebuilding the neighbours of the whole mesh.
static void relink(TetMesh& mesh, const std::vector<int>& fresh, const std::vector<int>& outer) {
    for (int f : fresh)
        for (int i = 0; i < 4; i++) nbr(mesh.tets[f], i) = -1;

    std::vector<int> all = fresh;
    all.insert(all.end(), outer.begin(), outer.end());

    for (int a : fresh) {
        for (int i = 0; i < 4; i++) {
            if (nbr(mesh.tets[a], i) != -1) continue;    // already linked from the other side
            const std::array<int,3> key = faceKey(mesh.tets[a], i);
            for (int b : all) {
                if (b == a) continue;
                for (int j = 0; j < 4; j++) {
                    if (faceKey(mesh.tets[b], j) == key) {
                        nbr(mesh.tets[a], i) = b;
                        nbr(mesh.tets[b], j) = a;
                    }
                }
            }
        }
    }
}

bool apply23(TetMesh& mesh, int t1, int t2) {
    Tet out[3];
    if(flip23Tets(mesh, t1, t2, out) == false) return false;
    std::vector<int> outer = outerNeighbours(mesh, {t1, t2});   // before the old tets are overwritten
    mesh.tets[t1] = out[0];
    mesh.tets[t2] = out[1];
    mesh.addTet(out[2]);
    relink(mesh, {t1, t2, static_cast<int>(mesh.tets.size()) - 1}, outer);
    return true;
}

bool flip32Tets(const TetMesh& mesh, int t1, int t2, int t3, Tet out[2]) {
    Tet tet1 = mesh.tets[t1];
    Tet tet2 = mesh.tets[t2];
    Tet tet3 = mesh.tets[t3];

    int A[4] = {tet1.v0, tet1.v1, tet1.v2, tet1.v3};
    int B[4] = {tet2.v0, tet2.v1, tet2.v2, tet2.v3};
    int C[4] = {tet3.v0, tet3.v1, tet3.v2, tet3.v3}; 

    std::map<int,int> count;
    for (std::size_t i = 0; i < 4; i++) {
        count[A[i]]++;
        count[B[i]]++;
        count[C[i]]++;
    }

    int nE = 0;
    int nF = 0;
    int E[4];
    int F[4];
    for (const auto& kv : count) {
        int vertex = kv.first;
        int times  = kv.second;
        
        if (times == 3)      { if (nE >= 2) return false; E[nE++] = vertex; }
        else if (times == 2) { if (nF >= 3) return false; F[nF++] = vertex; }
        else return false;
    }
    if (nE != 2 || nF != 3) return false;

    int d = E[0];
    int e = E[1];

    Point3 pa = mesh.points[F[0]];
    Point3 pb = mesh.points[F[1]];
    Point3 pc = mesh.points[F[2]];
    Point3 pd = mesh.points[d];
    
    if (orientation(pa, pb, pc, pd) == 0)     { return false; }
    else if (orientation(pa, pb, pc, pd) < 0) { std::swap(F[0], F[1]); }

    out[0] = Tet{F[0], F[1], F[2], d}; 
    out[1] = Tet{F[0], F[1], e, F[2]};

    for (std::size_t i = 0; i < 2; i ++) {
        pa = mesh.points[out[i].v0];
        pb = mesh.points[out[i].v1];
        pc = mesh.points[out[i].v2];
        pd = mesh.points[out[i].v3];
        if (orientation(pa, pb, pc, pd) <= 0) { return false; }
    }

    return true;
}

bool apply32(TetMesh& mesh, int t1, int t2, int t3) {
    Tet out[2];
    if(flip32Tets(mesh, t1, t2, t3, out) == false) return false;
    std::vector<int> outer = outerNeighbours(mesh, {t1, t2, t3});
    mesh.tets[t1] = out[0];
    mesh.tets[t2] = out[1];
    relink(mesh, {t1, t2}, outer);   // nothing points at t3 any more

    // remove t3 by moving the last tet into its slot (O(1), only the moved tet's index changes)
    const int last = static_cast<int>(mesh.tets.size()) - 1;
    if (t3 != last) {
        mesh.tets[t3] = mesh.tets[last];
        Tet& moved = mesh.tets[t3];
        for (int i = 0; i < 4; i++) {
            int n = nbr(moved, i);
            if (n == -1) continue;
            Tet& o = mesh.tets[n];
            for (int j = 0; j < 4; j++)
                if (nbr(o, j) == last) nbr(o, j) = t3;
        }
    }
    mesh.tets.pop_back();
    return true;
}
