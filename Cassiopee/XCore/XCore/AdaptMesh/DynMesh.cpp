#include "DynMesh.h"
#include "common/Karray.h"
#include "Array.h"
#include "Skin.h"
#include "Vec.h"
#include "Point.h"
#include "common/mem.h"
#include <stack>

DynMesh::DynMesh()
{}

DynMesh::DynMesh(Karray *karray)
{
    np = karray->npoints();
    nf = karray->nfaces();
    nc = karray->ncells();

    X.resize(np);
    Y.resize(np);
    Z.resize(np);
    memcpy(X.data(), karray->x, np * sizeof(E_Float));
    memcpy(Y.data(), karray->y, np * sizeof(E_Float));
    memcpy(Z.data(), karray->z, np * sizeof(E_Float));

    F.reserve(nf);

    for (E_Int i = 0; i < nf; i++) {
        E_Int np = -1;
        E_Int *pn = karray->get_face(i, np);
        std::vector<E_Int> points(np);
        for (E_Int j = 0; j < np; j++)
            points[j] = pn[j] - 1;
        F.push_back(points);
    }

    C.reserve(nc);
    for (E_Int i = 0; i < nc; i++) {
        E_Int nf = -1;
        E_Int *pf = karray->get_cell(i, nf);
        std::vector<E_Int> faces(nf);
        for (E_Int j = 0; j < nf; j++)
            faces[j] = pf[j] - 1;
        C.push_back(faces);
    }

    ftag.resize(nf, 0);

    owner.resize(nf, -1);
    neigh.resize(nf, -1);
    for (E_Int cid = 0; cid < nc; cid++) {
        const auto &pf = C[cid];
        for (E_Int fid : pf) {
            if (owner[fid] == -1) owner[fid] = cid;
            else neigh[fid] = cid;
        }
    }
}

void DynMesh::extract_points_from_ftag(ArrayI *pids)
{
    E_Int *ptag = (E_Int *)XMALLOC(np * sizeof(E_Int));
    memset(ptag, 0, np * sizeof(E_Int));
    pids->count = 0;

    for (E_Int fid = 0; fid < nf; fid++) {
        if (ftag[fid] != 1) continue;
        const auto &pn = F[fid];
        
        for (E_Int pid : pn) {
            pids->count += (ptag[pid] == 0);
            ptag[pid] = 1;
        }
    }

    pids->ptr = (E_Int *)XMALLOC(pids->count * sizeof(E_Int));
    E_Int *ptr = pids->ptr;

    for (E_Int pid = 0; pid < np; pid++) {
        if (ptag[pid] == 1)
            *ptr++ = pid;
    }

    XFREE(ptag);
}

void DynMesh::make_skin_connectivity(SkinGraph *skin_graph)
{
    // Count
    printf("skin count: %d\n", skin_graph->nf);
    skin_graph->xadj = (E_Int *)XMALLOC((skin_graph->nf+1) * sizeof(E_Int));
    E_Int *xadj = skin_graph->xadj;
    xadj[0] = 0;

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        xadj[i+1] = F[fid].size();
        xadj[i+1] += xadj[i];
    }

    skin_graph->fpts = (E_Int *)XMALLOC(xadj[skin_graph->nf] * sizeof(E_Int)); 

    // Populate
    E_Int *ptr = skin_graph->fpts;

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        const auto &pn = F[fid];
        for (E_Int p : pn) {
            *ptr++ = p;
        }
    }

    assert(ptr - skin_graph->fpts == xadj[skin_graph->nf]);
}

void DynMesh::make_face_centers(const E_Int NF, const E_Int *skin,
    Vec3f *fc)
{
    for (E_Int i = 0; i < NF; i++) {
        fc[i].x = fc[i].y = fc[i].z = 0.0;
        E_Int fid = skin[i];
        const auto &pn = F[fid];
        assert(ftag[skin[i]] == 1);
        for (E_Int p : pn) {
            fc[i].x += X[p];
            fc[i].y += Y[p];
            fc[i].z += Z[p];
        }
        fc[i].x /= pn.size(); fc[i].y /= pn.size(); fc[i].z /= pn.size();
    }
}

bool DynMesh::point_in_tri(const Point *p, E_Int tid) const
{
    const auto &pn = F[tid];
    E_Int A = pn[0], B = pn[1], C = pn[2];
    return Point_in_tri(p->x, p->y, p->z,
                        X[A], Y[A], Z[A],
                        X[B], Y[B], Z[B],
                        X[C], Y[C], Z[C]);
}

bool DynMesh::point_in_quad(const Point *p, E_Int qid) const
{
    // TODO(Imad): maybe compute face centers once in pre-pass
    // Star the quad into 4 triangles
    E_Float O[3] = {0.0, 0.0, 0.0};
    const auto &pn = F[qid];
    E_Int A = pn[0], B = pn[1], C = pn[2], D = pn[3];
    O[0] = (X[A] + X[B] + X[C] + X[D]) * 0.25;
    O[1] = (Y[A] + Y[B] + Y[C] + Y[D]) * 0.25;
    O[2] = (Z[A] + Z[B] + Z[C] + Z[D]) * 0.25;

    bool hit = false;

    // First triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[A], Y[A], Z[A],
                       X[B], Y[B], Z[B]);
    if (hit) return true;

    // Second triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[B], Y[B], Z[B],
                       X[C], Y[C], Z[C]);
    if (hit) return true;


    // Third triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[C], Y[C], Z[C],
                       X[D], Y[D], Z[D]);
    if (hit) return true;

    // Fourth triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[D], Y[D], Z[D],
                       X[A], Y[A], Z[A]);
    if (hit) return true;

    return false;
}

bool DynMesh::point_in_face(const Point *p, E_Int fid) const
{
    if (F[fid].size() == 4) return point_in_quad(p, fid);
    return point_in_tri(p, fid);
}

void DynMesh::prepare_for_refinement(ArrayI *ref_faces)
{
    ref_faces->count = 0;
    for (size_t i = 0; i < tri_graph.nf; i++) {
        ref_faces->count += (tri_graph.fdat[i] == 1);
    }
    ref_faces->ptr = (E_Int *)XMALLOC(ref_faces->count * sizeof(E_Int));
    E_Int *ptr = ref_faces->ptr;
    for (size_t i = 0; i < tri_graph.nf; i++) {
        if (tri_graph.fdat[i] == 1) {
            *ptr++ = i;
        }
    }

    // Resize data structures
    resize_point_data(ref_faces->count);
    resize_face_data(ref_faces->count);
}

void DynMesh::refine_faces(ArrayI *ref_faces)
{
    tri_graph.T.resize(tri_graph.nf + ref_faces->count*4, {-1, -1, -1});
    tri_graph.E.resize(tri_graph.nf + ref_faces->count*4, {-1, -1, -1});
    tri_graph.skin.resize(tri_graph.nf + ref_faces->count*4, -1);
    tri_graph.level.resize(tri_graph.nf + ref_faces->count*4, 0);

    printf("Faces before refinement: %d\n", nf);

    for (E_Int i = 0; i < ref_faces->count; i++) {
        E_Int fid = ref_faces->ptr[i];
        assert(tri_graph.level[fid] == 0);
        assert(face_is_tri(tri_graph.skin[fid]));
        refine_tri(fid);
    }
}

void DynMesh::resize_point_data(size_t nref_faces)
{
    size_t nnew_points = np + nref_faces * 3;
    X.resize(nnew_points);
    Y.resize(nnew_points);
    Z.resize(nnew_points);
}

void DynMesh::resize_face_data(size_t nref_faces)
{
    size_t nnew_faces = nf + nref_faces * 4; // 3 new triangles + 1 potential neighbor triangle
    F.resize(nnew_faces);
    owner.resize(nnew_faces);
    neigh.resize(nnew_faces);
    ftag.resize(nnew_faces, 0);
}

E_Int DynMesh::get_edge_index(E_Int nei, E_Int tri)
{
    for (E_Int i = 0; i < 3; i++) {
        if (tri_graph.E[nei][i] == tri)
            return i;
    }
    return -1;
}

void DynMesh::refine_tri(E_Int tri_idx)
{
    assert(tri_graph.level[tri_idx] == 0);
    puts("refining");
    fflush(stdout);
    E_Int tri = tri_graph.skin[tri_idx];

    // Refine the edges
    const auto &pn = F[tri];
    E_Int ec[3];

    assert(pn.size() == 3);

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        E_Int q = pn[(i+1)%pn.size()];
        uEdge edge(p, q);

        auto it = ecenter.find(edge);

        if (it == ecenter.end()) {
            X[np] = 0.5 * (X[p] + X[q]);
            Y[np] = 0.5 * (Y[p] + Y[q]);
            Z[np] = 0.5 * (Z[p] + Z[q]);
            ecenter[edge] = np;
            ec[i] = np;
            np++;
        } else {
            ec[i] = it->second;
        }
    }
    
    E_Int vn[3] = {pn[0], pn[1], pn[2]};

    // Create new triangles
    F[tri]  = {vn[0], ec[0], ec[2]};
    F[nf]   = {ec[0], vn[1], ec[1]};
    F[nf+1] = {ec[2], ec[1], vn[2]};
    F[nf+2] = {ec[0], ec[1], ec[2]};

    fchildren[tri] = {nf, nf+1, nf+2};
    ftag[nf] = ftag[nf+1] = ftag[nf+2] = ftag[tri];

    // Update owners and neighbours
    neigh[nf] = neigh[nf+1] = neigh[nf+2] = neigh[tri];
    owner[nf] = owner[nf+1] = owner[nf+2] = owner[tri];

    // Skin

    // Add the new triangles
    size_t NF = tri_graph.nf;
    
    auto &T = tri_graph.T;
    auto &E = tri_graph.E;
    
    T[tri_idx][0] = vn[0]; T[tri_idx][1] = ec[0]; T[tri_idx][2] = ec[2];
    assert(T[tri_idx][0] != T[tri_idx][1]);
    assert(T[tri_idx][0] != T[tri_idx][2]);
    assert(T[tri_idx][1] != T[tri_idx][2]);

    T[NF][0]   = ec[0]; T[NF][1]   = vn[1]; T[NF][2]   = ec[1];
    T[NF+1][0] = ec[2]; T[NF+1][1] = ec[1]; T[NF+1][2] = vn[2];
    T[NF+2][0] = ec[0]; T[NF+2][1] = ec[1]; T[NF+2][2] = ec[2];
    
    E_Int N0 = NF;
    E_Int N1 = NF+1;
    E_Int N2 = NF+2;

    E_Int A = tri_graph.E[tri_idx][0];
    E_Int B = tri_graph.E[tri_idx][1];
    E_Int C = tri_graph.E[tri_idx][2];
    
    E[tri_idx][0] = A; E[tri_idx][1] = N2; E[tri_idx][2] = C;
    E[N0][0] = A;  E[N0][1] = B;  E[N0][2] = N2;
    E[N1][0] = N2; E[N1][1] = B;  E[N1][2] = C;
    E[N2][0] = N0; E[N2][1] = N1; E[N2][2] = tri_idx;
 
    tri_graph.skin[NF] = nf;
    tri_graph.skin[NF+1] = nf+1;
    tri_graph.skin[NF+2] = nf+2;

    tri_graph.level[tri_idx]++;
    assert(tri_graph.level[tri_idx] == 1);
    tri_graph.level[NF] = tri_graph.level[NF+1] = tri_graph.level[NF+2] = 1;
    
    // Increment face count
    nf += 3; 
    tri_graph.nf += 3;

    NF = tri_graph.nf;

    bool cut_A = (A != -1) && (tri_graph.level[A] == 0 && tri_graph.fdat[A] == 0); 
    
    if (cut_A) {
        E_Int edge_idx = get_edge_index(A, tri_idx);
        assert(edge_idx != -1);
        E_Int P = T[A][(edge_idx+2)%3];

        E_Int D = E[A][(edge_idx+1)%3];
        E_Int e = E[A][(edge_idx+2)%3];

        // A points and neighbours
        T[A][0] = ec[0];   T[A][1] = vn[0]; T[A][2] = P;

        // NF points and neighbours
        T[NF][0] = ec[0]; T[NF][1] = P; T[NF][2] = vn[1];
        
        // e neighbours (if it exists)
        if (e != -1) {
            E_Int idx = get_edge_index(e, A);
            E[e][idx] = NF;
        }

        E[A][0] = tri_idx; E[A][1] = D; E[A][2] = NF;
        E[NF][0] = A;     E[NF][1] = e; E[NF][2] = N0;
        
        // N0 neighbours
        E[N0][0] = NF;// E[N1][0] = B; E[N1][1] = N2;
        
        tri_graph.skin[NF] = nf;

        assert(tri_graph.level[A] == 0);
        assert(tri_graph.level[NF] == 0);
        
        NF += 1;
        
        // Global stuff
        E_Int gA = tri_graph.skin[A];
        F[gA] = {ec[0], vn[0], P};
        fchildren[gA].push_back(nf);

        F[nf] = {ec[0], P, vn[1]};
        ftag[nf] = ftag[gA];
        owner[nf] = owner[gA];
        neigh[nf] = neigh[gA];

        nf += 1;
    }

    bool cut_B = (B != -1) && (tri_graph.level[B] == 0 && tri_graph.fdat[B] == 0); 
    
    if (cut_B) {
        E_Int edge_idx = get_edge_index(B, tri_idx);
        assert(edge_idx != -1);
        E_Int P = T[B][(edge_idx+2)%3];

        E_Int f = E[B][(edge_idx+1)%3];
        E_Int G = E[B][(edge_idx+2)%3];

        // B points and neighbours
        T[B][0] = ec[1]; T[B][1] = vn[1]; T[B][2] = P;

        // NF points and neighbours
        T[NF][0] = ec[1]; T[NF][1] = P; T[NF][2] = vn[2];
        
        // G neighbours (if it exists)
        if (G != -1) {
            E_Int idx = get_edge_index(G, B);
            E[G][idx] = NF;
        }

        E[B][0] = N0;    E[B][1] = f;     E[B][2] = NF;
        E[NF][0] = B;     E[NF][1] = G; E[NF][2] = N1;

        // N1 neighbours
        E[N1][1] = NF;
        
        assert(tri_graph.level[B] ==  0);
        assert(tri_graph.level[NF] == 0);
        
        tri_graph.skin[NF] = nf;
        NF += 1;
        
        // Global stuff
        E_Int gB = tri_graph.skin[B];
        F[gB] = {ec[1], vn[1], P};
        fchildren[gB].push_back(nf);

        F[nf] = {ec[1], P, vn[2]};
        ftag[nf] = ftag[gB];
        owner[nf] = owner[gB];
        neigh[nf] = neigh[gB];

        nf += 1;
    }

    bool cut_C = (C != -1) && (tri_graph.level[C] == 0 && tri_graph.fdat[C] == 0); 

    if (cut_C) {
        E_Int edge_idx = get_edge_index(C, tri_idx);
        assert(edge_idx != -1);
        E_Int P = T[C][(edge_idx+2)%3];

        E_Int I = E[C][(edge_idx+1)%3];
        E_Int H = E[C][(edge_idx+2)%3];
        
        // C points and neighbours
        T[C][0] = ec[2]; T[C][1] = P; T[C][2] = vn[0];

        // NF points and neighbours
        T[NF][0] = ec[2]; T[NF][1] = vn[2]; T[NF][2] = P;
        
        // I neighbours (if it exists)
        if (I != -1) {
            E_Int idx = get_edge_index(I, C);
            E[I][idx] = NF;
        }

        E[C][0] = NF;    E[C][1] = H; E[C][2] = tri_idx;
        E[NF][0] = N1;    E[NF][1] = I;     E[NF][2] = C;
        

        // N1 adjacency
        E[N1][2] = NF;
        
        assert(tri_graph.level[C] == 0);
        assert(tri_graph.level[NF] == 0);
        
        tri_graph.skin[NF] = nf;
        NF += 1;
        
        // Global stuff
        E_Int gC = tri_graph.skin[C];
        F[gC] = {ec[2], vn[2], P};
        fchildren[gC].push_back(nf);

        F[nf] = {ec[2], P, vn[0]};
        ftag[nf] = ftag[gC];
        owner[nf] = owner[gC];
        neigh[nf] = neigh[gC];

        nf += 1;
    }
}

DynMesh DynMesh::extract_conformized()
{
    // Keep all the points
    std::vector<E_Float> new_X(X), new_Y(Y), new_Z(Z);

    // Conformize the faces

    std::vector<std::vector<E_Int>> new_F(nf);

    for (E_Int fid = 0; fid < nf; fid++) {

        const auto &pn = F[fid];

        auto &new_face = new_F[fid];

        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];

            std::list<E_Int> epoints;

            extract_edge_points(p, q, epoints);

            epoints.pop_back();

            for (auto it = epoints.begin(); it != epoints.end(); it++)
                new_face.push_back(*it);
        }
    }

    // Update cell connectivity

    std::vector<std::vector<E_Int>> new_C(C.size());

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];

        auto &new_cell = new_C[i];

        for (E_Int fid : pf) {
            auto it = fchildren.find(fid);

            if (it == fchildren.end()) {
                new_cell.push_back(fid);
            } else {
                std::vector<E_Int> fleaves;
                get_fleaves(fid, fleaves);
                for (E_Int leaf : fleaves)
                    new_cell.push_back(leaf);
            }
        }
    }

    DynMesh new_M;
    new_M.np = np;
    new_M.X = X;
    new_M.Y = Y;
    new_M.Z = Z;
    new_M.nf = nf;
    new_M.F = new_F;
    new_M.nc = nc;
    new_M.C = new_C;
    new_M.ftag = ftag;


    return new_M;
}

void DynMesh::get_fleaves(E_Int face, std::vector<E_Int> &fleaves)
{
    fleaves.push_back(face);
    
    const auto it = fchildren.find(face);

    if (it == fchildren.end()) {
        return;
    }
    
    for (E_Int child : it->second)
        get_fleaves(child, fleaves);
}

PyObject *DynMesh::export_karray()
{
    E_Int sizeNGon = 0, sizeNFace = 0;

    for (const auto &pn : F) sizeNGon += (E_Int)pn.size();
    for (const auto &pf : C) sizeNFace += (E_Int)pf.size();

    const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

    PyObject *array = K_ARRAY::buildArray3(3, varString, np, nc, nf, "NGON",
        sizeNGon, sizeNFace, 3, false, 3);
    
    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(array, f, cn);

    E_Float *px = f->begin(1);
    for (E_Int i = 0; i < np; i++) px[i] = X[i];
    E_Float *py = f->begin(2);
    for (E_Int i = 0; i < np; i++) py[i] = Y[i];
    E_Float *pz = f->begin(3);
    for (E_Int i = 0; i < np; i++) pz[i] = Z[i];

    E_Int *indPG = cn->getIndPG();
    E_Int *ngon = cn->getNGon();
    E_Int *indPH = cn->getIndPH();
    E_Int *nface = cn->getNFace();

    indPG[0] = indPH[0] = 0;
    for (E_Int i = 0; i < nf; i++) indPG[i+1] = indPG[i] + (E_Int)F[i].size();
    for (E_Int i = 0; i < nc; i++) indPH[i+1] = indPH[i] + (E_Int)C[i].size();

    assert(indPG[nf] == sizeNGon);
    assert(indPH[nc] == sizeNFace);

    E_Int *ptr = ngon;

    for (E_Int i = 0; i < nf; i++) {
        const auto &pn = F[i];
        for (E_Int p : pn) *ptr++ = p+1;
    }

    ptr = nface;

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (E_Int f : pf) *ptr++ = f+1;
    }

    delete f;
    delete cn;

    return array;
}

void DynMesh::extract_edge_points(E_Int a, E_Int b, std::list<E_Int> &points)
{
    E_Int ref = 0;

    points.clear();
    points.push_back(a);
    points.push_back(b);

    do {
        ref = 0;

        assert(*std::prev(points.end()) == b);

        for (auto it = points.begin(); it != std::prev(points.end()); it++) {
            E_Int a = *it;
            E_Int b = *std::next(it);

            uEdge e(a, b);

            auto search = ecenter.find(e);

            if (search != ecenter.end()) {
                points.insert(std::next(it), search->second);
                ref = 1;
            }
        }
    } while (ref);
}

struct EdgeNode {
    E_Int p, q;
    E_Int fi, posi;
    mutable E_Int fj, posj;
    EdgeNode(E_Int p_, E_Int q_)
    {
        p = std::min(p_, q_);
        q = std::max(p_, q_);
        fi = posi = fj = posj = -1;
    }
    bool operator<(const EdgeNode &e) const
    {
        return (p < e.p) || (p == e.p && q < e.q);
    }
};

void DynMesh::make_tri_graph()
{
    tri_graph.nf = 0;
    
    // From tagged faces
    for (E_Int i = 0; i < nf; i++) {
        tri_graph.nf += (ftag[i] == 1);
    }
    tri_graph.skin.clear();
    tri_graph.skin.resize(tri_graph.nf);
    E_Int *skin = tri_graph.skin.data();
    E_Int *ptr = skin;
    for (E_Int i = 0; i < nf; i++) {
        if (ftag[i] == 1)
            *ptr++ = i;
    }

    tri_graph.T.clear();
    tri_graph.E.clear();
    tri_graph.level.clear();

    tri_graph.T.resize(tri_graph.nf, {-1, -1, -1});
    tri_graph.E.resize(tri_graph.nf, {-1, -1, -1});
    tri_graph.level.resize(tri_graph.nf, 0);

    std::set<EdgeNode> edges;

    for (size_t i = 0; i < tri_graph.nf; i++) {
        E_Int fid = tri_graph.skin[i];
        const auto &pn = F[fid];
        assert(pn.size() == 3);
        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];

            tri_graph.T[i][j] = p;

            E_Int q = pn[(j+1)%3];

            EdgeNode node(p, q);
            auto it = edges.find(node);
            if (it == edges.end()) {
                node.fi = i;
                node.posi = j;
                edges.insert(node);
            } else {
                assert(it->fi != -1);
                assert(it->posi != -1);
                assert(it->fj == -1);
                assert(it->posj == -1);
                it->fj = i;
                it->posj = j;
            }
        }
    }

    for (const auto &e : edges) {
        assert(e.fi != -1);
        tri_graph.E[e.fi][e.posi] = e.fj;
        if (e.fj != -1)
            tri_graph.E[e.fj][e.posj] = e.fi;
        else
            assert(e.posj == -1);
    }
}


void DynMesh::triangulate(const E_Int *faces, E_Int fcount)
{
    F.resize(nf + fcount);
    owner.resize(nf + fcount, -1);
    neigh.resize(nf + fcount, -1);
    ftag.resize(nf + fcount, 0);

    for (E_Int i = 0; i < fcount; i++) {
        E_Int fid = faces[i];

        auto &pn = F[fid];
        auto &tri = F[nf];
        
        tri = { pn[0], pn[2], pn[3] };
        pn = { pn[0], pn[1], pn[2] };

        auto &pown = C[owner[fid]];
        pown.push_back(nf);

        if (neigh[fid] != -1) {
            auto &pnei = C[neigh[fid]];
            pnei.push_back(nf);
        }

        owner[nf] = owner[fid];
        neigh[nf] = neigh[fid];
        ftag[nf] = ftag[fid];

        nf++;
    }
}
