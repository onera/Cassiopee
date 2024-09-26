#include "DynMesh.h"
#include "Karray.h"
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
    memcpy(X.data(), karray->X(), np * sizeof(E_Float));
    memcpy(Y.data(), karray->Y(), np * sizeof(E_Float));
    memcpy(Z.data(), karray->Z(), np * sizeof(E_Float));

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

    for (E_Int fid : factive) {
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

void DynMesh::extract_skin(E_Int *count, E_Int **skin)
{
    *count = 0;
    
    for (E_Int fid : factive) {
        *count += (neigh[fid] == -1);
    }

    *skin = (E_Int *)XMALLOC(*count * sizeof(E_Int));
    E_Int *ptr = *skin;

    for (E_Int fid : factive) {
        if (neigh[fid] == -1)
            *ptr++ = fid;
    }
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

void DynMesh::make_skin_graph(SkinGraph *skin_graph)
{
    extract_skin(&skin_graph->nf, &skin_graph->skin);
    make_skin_connectivity(skin_graph);
    SkinGraph_make_skin_neighbours(skin_graph);
}

void DynMesh::make_face_centers(const E_Int NF, const E_Int *skin,
    Vec3f *fc)
{
    for (E_Int i = 0; i < NF; i++) {
        fc[i].x = fc[i].y = fc[i].z = 0.0;
        E_Int fid = skin[i];
        assert(face_is_active(fid));
        const auto &pn = F[fid];
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

void DynMesh::init_adaptation_data(E_Int *tagged_faces, E_Int count)
{
    flevel.resize(nf, 0);

    for (E_Int i = 0; i < count; i++) factive.insert(tagged_faces[i]);
}

void DynMesh::init_adaptation_data()
{
    flevel.resize(nf, 0);

    for (E_Int i = 0; i < nf; i++) factive.insert(i);
}

void DynMesh::prepare_for_refinement(ArrayI *ref_faces)
{
    ref_faces->count = 0;
    for (E_Int fid : factive) {
        ref_faces->count += (fref[fid] == 1);
    }
    ref_faces->ptr = (E_Int *)XMALLOC(ref_faces->count * sizeof(E_Int));
    E_Int *ptr = ref_faces->ptr;
    for (E_Int fid = 0; fid < nf; fid++) {
        if (fref[fid] > 0) {
            assert(fref[fid] == 1);
            *ptr++ = fid;
        }
    }

    // Refine the lower-level faces first
    std::sort(ref_faces->ptr, ref_faces->ptr + ref_faces->count,
        [&] (E_Int i, E_Int j) { return flevel[i] < flevel[j]; });
    
    // Resize data structures
    resize_point_data(ref_faces->count);
    resize_face_data(ref_faces->count);
}

void DynMesh::refine_faces(ArrayI *ref_faces)
{
    for (E_Int i = 0; i < ref_faces->count; i++) {
        E_Int fid = ref_faces->ptr[i];
        if (face_is_tri(fid)) refine_tri(fid);
        else refine_quad(fid);
    }
}

void DynMesh::resize_point_data(size_t nref_faces)
{
    size_t nnew_points = np + nref_faces * 5;
    X.resize(nnew_points);
    Y.resize(nnew_points);
    Z.resize(nnew_points);
}

void DynMesh::resize_face_data(size_t nref_faces)
{
    size_t nnew_faces = nf + nref_faces * 4;
    F.resize(nnew_faces);
    flevel.resize(nnew_faces, -1);
    owner.resize(nnew_faces);
    neigh.resize(nnew_faces);
    ftag.resize(nnew_faces, 0);
}

void DynMesh::refine_tri(E_Int tri)
{
    // Refine the edges
    const auto &pn = F[tri];
    E_Int ec[3];

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

    // New face points
    E_Int nf0 = nf, nf1 = nf+1, nf2 = nf+2, nf3 = nf+3;

    F[nf0] = { pn[0], ec[0], ec[2] };
    F[nf1] = { ec[0], pn[1], ec[1] };
    F[nf2] = { ec[2], ec[1], pn[2] };
    F[nf3] = { ec[0], ec[1], ec[2] };

    // Disable quad and enable its children
    factive.erase(tri);
    factive.insert(nf0);
    factive.insert(nf1);
    factive.insert(nf2);
    factive.insert(nf3);

    // Set quad children pointers
    fchildren[tri] = { nf0, nf1, nf2, nf3 };
    flevel[nf0] = flevel[nf1] = flevel[nf2] = flevel[nf3] = flevel[tri] + 1;

    ftag[nf0] = ftag[nf1] = ftag[nf2] = ftag[nf3] = ftag[tri];

    neigh[nf0] = neigh[nf1] = neigh[nf2] = neigh[nf3] = neigh[tri];
    owner[nf0] = owner[nf1] = owner[nf2] = owner[nf3] = owner[tri];

    nf += 4;
}

void DynMesh::refine_quad(E_Int quad)
{
    // Refine the edges
    const auto &pn = F[quad];
    E_Int ec[4];

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

    // Face centroid

    X[np] = Y[np] = Z[np] = 0.0;
    for (E_Int i = 0; i < 4; i++) {
        E_Int p = pn[i];
        X[np] += X[p];
        Y[np] += Y[p];
        Z[np] += Z[p];
    }
    X[np] *= 0.25;
    Y[np] *= 0.25;
    Z[np] *= 0.25;

    // New face points
    E_Int nf0 = nf, nf1 = nf+1, nf2 = nf+2, nf3 = nf+3;

    F[nf0] = { pn[0], ec[0], np,    ec[3] };
    F[nf1] = { ec[0], pn[1], ec[1], np    };
    F[nf2] = { np   , ec[1], pn[2], ec[2] };
    F[nf3] = { ec[3], np   , ec[2], pn[3] };

    // Disable quad and enable its children
    factive.erase(quad);
    factive.insert(nf0);
    factive.insert(nf1);
    factive.insert(nf2);
    factive.insert(nf3);

    // Set quad children pointers
    fchildren[quad] = { nf0, nf1, nf2, nf3 };
    flevel[nf0] = flevel[nf1] = flevel[nf2] = flevel[nf3] = flevel[quad] + 1;

    ftag[nf0] = ftag[nf1] = ftag[nf2] = ftag[nf3] = ftag[quad];

    neigh[nf0] = neigh[nf1] = neigh[nf2] = neigh[nf3] = neigh[quad];
    owner[nf0] = owner[nf1] = owner[nf2] = owner[nf3] = owner[quad];

    np += 1;
    nf += 4;
}

DynMesh DynMesh::extract_conformized()
{
    // Keep all the points
    std::vector<E_Float> new_X(X), new_Y(Y), new_Z(Z);

    // Conformize the faces

    std::vector<std::vector<E_Int>> new_F(factive.size());

    E_Int new_nf = 0;
    
    std::map<E_Int, E_Int> new_fids;

    for (E_Int face : factive) {
        new_fids[face] = new_nf;

        const auto &pn = F[face];

        auto &new_face = new_F[new_nf];

        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];

            std::list<E_Int> epoints;

            extract_edge_points(p, q, epoints);

            epoints.pop_back();

            for (auto it = epoints.begin(); it != epoints.end(); it++)
                new_face.push_back(*it);
        }

        new_nf++;
    }

    // Update cell connectivity

    std::vector<std::vector<E_Int>> new_C(C.size());

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];

        auto &new_cell = new_C[i];

        for (E_Int face : pf) {

            if (face_is_active(face)) {
                new_cell.push_back(new_fids[face]);
            } else {
                std::vector<E_Int> fleaves;
                get_fleaves(face, fleaves);

                for (E_Int fleaf : fleaves)
                    new_cell.push_back(new_fids[fleaf]);
            }
        }
    }

    DynMesh new_M;
    new_M.np = np;
    new_M.X = X;
    new_M.Y = Y;
    new_M.Z = Z;
    new_M.nf = new_nf;
    new_M.F = new_F;
    new_M.nc = nc;
    new_M.C = new_C;

    for (E_Int face : factive) {
        new_M.factive.insert(new_fids[face]);
    }

    return new_M;
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

void DynMesh::get_fleaves(E_Int face, std::vector<E_Int> &fleaves)
{
    if (face_is_active(face)) {
        fleaves.push_back(face);
        return;
    }

    for (E_Int child : fchildren.at(face)) get_fleaves(child, fleaves);
}