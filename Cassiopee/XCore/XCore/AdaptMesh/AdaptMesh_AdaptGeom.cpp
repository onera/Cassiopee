#include "Mesh.h"
#include "common/mem.h"
#include "Array.h"
#include "Point.h"
#include "BVH.h"
#include "Point.h"
#include "Skin.h"
#include "Vec.h"
#include "DynMesh.h"
#include <stack>

// TODO(Imad): sorting routines

void locate_points_in_skin
(
    const Mesh *source,
    const ArrayI *points,
    const DynMesh *target,
    const E_Int *skin,
    const E_Int *indices,
    const BVH_node *bvh,
    PointFaces *ploc
)
{
    for (E_Int i = 0; i < points->count; i++) {
        E_Int pid = points->ptr[i];
        const Point p = {source->X[pid], source->Y[pid], source->Z[pid]};
        PointFaces *pfaces = &ploc[i];
        //assert(Point_in_Box3D(&p, &bvh->box));
        BVH_locate_point(bvh, target, skin, indices, &p, pfaces);
        /*
        if (pfaces->count == 0) {
            fprintf(stderr, "bvh_locate: failed at point index %d (%d)!\n",
                i, spid);
            point_write(&p);
            abort();
        }
        */
    }
}

void locate_points_in_skin
(
    const DynMesh *source,
    const ArrayI *points,
    const Mesh *target,
    const E_Int *skin,
    const E_Int *indices,
    const BVH_node *bvh,
    PointFaces *ploc
)
{
    for (E_Int i = 0; i < points->count; i++) {
        E_Int pid = points->ptr[i];
        const Point p = {source->X[pid], source->Y[pid], source->Z[pid]};
        PointFaces *pfaces = &ploc[i];
        //assert(Point_in_Box3D(&p, &bvh->box));
        BVH_locate_point(bvh, target, skin, indices, &p, pfaces);
        /*
        if (pfaces->count == 0) {
            fprintf(stderr, "bvh_locate: failed at point index %d (%d)!\n",
                i, spid);
            point_write(&p);
            abort();
        }
        */
    }
}

void init_skin_refinement_cells(const SkinGraph *skin_graph,
    const E_Int *fdat, Mesh *M, E_Int *rcount)
{
    E_Int nf = skin_graph->nf;
    E_Int count = 0;

    for (E_Int i = 0; i < nf; i++) {
        if (fdat[i] > 0) {
            E_Int fid = skin_graph->skin[i];
            E_Int own = M->owner[fid];
            // If cref[own] != 0, H18 is impossible, refinement has to be H27.
            assert(M->cref[own] == 0);
            M->cref[own] = 1;
            count++;
        }
    }

    *rcount = count;
}

void Mesh_set_face_as_cell_bottom(Mesh *M, E_Int fid, E_Int cid)
{
    // Make fid the bottom face
    E_Int *cell = Mesh_get_cell(M, cid);
    E_Int size = 4*M->cstride[cid];
    E_Int pos = Get_pos(fid, cell, size);
    assert(pos != -1);
    assert(pos % 4 == 0);
    Right_shift(cell, pos, size);
    assert(cell[0] == fid);
    E_Int *crange = Mesh_get_crange(M, cid);
    Right_shift(crange, pos/4, M->cstride[cid]);
}

void assign_face_refinement_data(Mesh *M)
{
    for (E_Int cid = 0; cid < M->nc; cid++) {
        if (M->cref[cid] == 0) continue;
        assert(M->cref[cid] == 1);

        E_Int *cell = Mesh_get_cell(M, cid);
        E_Int *crange = Mesh_get_crange(M, cid);
        E_Int cstride = M->cstride[cid];
        E_Int clvl = M->clevel[cid];

        for (E_Int i = 0; i < cstride; i++) {
            E_Int *pf = cell + 4*i;

            for (E_Int j = 0; j < crange[i]; j++) {
                E_Int face = pf[j];
                if (M->fref[face] == 1) continue;

                E_Int flvl = M->flevel[face];

                assert(flvl >= clvl);

                if (flvl == clvl) M->fref[face] = 1;
            }
        }
    }
}

// TODO(Imad): isotropic resizing for now

void Mesh_isolate_refinement_entities(Mesh *M, ArrayI *rcells, ArrayI *rfaces)
{
    rcells->count = 0;
    for (E_Int cid = 0; cid < M->nc; cid++) {
        rcells->count += (M->cref[cid] == 1);
    }
    rcells->ptr = (E_Int *)XMALLOC(rcells->count * sizeof(E_Int));
    E_Int *ptr = rcells->ptr;
    for (E_Int cid = 0; cid < M->nc; cid++) {
        if (M->cref[cid] > 0) {
            assert(M->cref[cid] == 1);
            *ptr++ = cid;
        }
    }

    rfaces->count = 0;
    for (E_Int fid = 0; fid < M->nf; fid++) {
        rfaces->count += (M->fref[fid] == 1);
    }
    rfaces->ptr = (E_Int *)XMALLOC(rfaces->count * sizeof(E_Int));
    ptr = rfaces->ptr;
    for (E_Int fid = 0; fid < M->nf; fid++) {
        if (M->fref[fid] > 0) {
            assert(M->fref[fid] == 1);
            *ptr++ = fid;
        }
    }
}

void Mesh_resize(Mesh *M, const ArrayI *rcells, const ArrayI *rfaces)
{
    // 7 new cells per refined cells
    E_Int cell_incr = rcells->count * 7;
    // 3 new faces per refined face + 12 new faces per refined cell
    E_Int face_incr = rfaces->count * 3 + rcells->count * 12;
    // Estimate point increase
    E_Int point_incr = 2 * face_incr;

    E_Int new_nc = M->nc + cell_incr;
    E_Int new_nf = M->nf + face_incr;
    E_Int new_np = M->np + point_incr;

    M->X = (E_Float *)XRESIZE(M->X, new_np * sizeof(E_Float));
    M->Y = (E_Float *)XRESIZE(M->Y, new_np * sizeof(E_Float));
    M->Z = (E_Float *)XRESIZE(M->Z, new_np * sizeof(E_Float));

    Mesh_resize_face_data(M, new_nf);
    Mesh_resize_cell_data(M, new_nc);
}

E_Int refine_quad_X(E_Int quad, Mesh *M)
{
    assert(M->fref[quad] == 1);
    E_Int NODES[8];

    E_Int *fpts = Mesh_get_face(M, quad);
    E_Int *frange = Mesh_get_frange(M, quad);

    // BOT
    NODES[0] = fpts[0];
    if (frange[0] == 2) {
        NODES[4] = fpts[1];
    } else {
        E_Int p = fpts[0];
        E_Int q = fpts[2];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[4]);
    }

    // RGT
    NODES[1] = fpts[2];
    NODES[5] = fpts[3];

    // TOP
    NODES[2] = fpts[4];
    if (frange[2] == 2) {
        NODES[6] = fpts[5];
    } else {
        E_Int p = fpts[4];
        E_Int q = fpts[6];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[6]);
    }

    // LFT
    NODES[3] = fpts[6];
    NODES[7] = fpts[7];


    // Second child
    fpts = Mesh_get_face(M, M->nf);
    fpts[0] = NODES[4];
    fpts[1] = -1;
    fpts[2] = NODES[1];
    fpts[3] = NODES[5];
    fpts[4] = NODES[2];
    fpts[5] = -1;
    fpts[6] = NODES[6];
    fpts[7] = -1;

    // First child replaces quad
    fpts = Mesh_get_face(M, quad);
    fpts[0] = NODES[0];
    fpts[1] = -1;
    fpts[2] = NODES[4];
    fpts[3] = -1;
    fpts[4] = NODES[6];
    fpts[5] = -1;
    fpts[6] = NODES[3];
    fpts[7] = NODES[7];

    // Update ranges and strides
    Mesh_update_face_range_and_stride(M, quad, M->nf, 1);

    // Conformize parent cells

    E_Int own = M->owner[quad];

    assert(M->clevel[own] == M->flevel[quad]);

    if (Mesh_conformize_cell_face(M, own, quad, M->nf, 2) != 0) return 1;

    E_Int nei = M->neigh[quad];

    if (nei != -1) {
        assert(M->clevel[nei] == M->flevel[quad]);

        if (Mesh_conformize_cell_face(M, nei, quad, M->nf, 2) != 0) return 1;
    }

    for (E_Int i = 0; i < 1; i++) {
        M->owner[M->nf+i] = own;
        M->neigh[M->nf+i] = nei;
    }

    // Update adaptation info
    M->flevel[quad]++;

    assert(M->fref[quad] == FACE_REFINED);

    for (E_Int i = 0; i < 1; i++) {
        E_Int fid = M->nf + i;

        M->flevel[fid] = M->flevel[quad];
        M->ftype[fid] = M->ftype[quad];

        M->fref[fid] = FACE_NEW;
    }

    M->fchildren[quad] = {quad, M->nf};

    M->fparent[M->nf] = quad;

    M->ftag[M->nf] = M->ftag[quad];

    M->fpattern[M->nf] = M->fpattern[quad];

    // Increment face/edge/point count
    M->nf += 1;

    return 0;
}

E_Int refine_quad_Y(E_Int quad, Mesh *M)
{
    assert(M->fref[quad] == 1);
    E_Int NODES[8];

    E_Int *fpts = Mesh_get_face(M, quad);
    E_Int *frange = Mesh_get_frange(M, quad);

    // BOT
    NODES[0] = fpts[0];
    NODES[4] = fpts[1];

    // RGT
    NODES[1] = fpts[2];
    if (frange[1] == 2) {
        NODES[5] = fpts[3];
    } else {
        E_Int p = fpts[2];
        E_Int q = fpts[4];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[5]);
    }

    // TOP
    NODES[2] = fpts[4];
    NODES[6] = fpts[5];
    
    // LFT
    NODES[3] = fpts[6];
    if (frange[3] == 2) {
        NODES[7] = fpts[7];
    } else {
        E_Int p = fpts[6];
        E_Int q = fpts[0];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[7]);
    }

    // Second child
    fpts = Mesh_get_face(M, M->nf);
    fpts[0] = NODES[7];
    fpts[1] = -1;
    fpts[2] = NODES[5];
    fpts[3] = -1;
    fpts[4] = NODES[2];
    fpts[5] = NODES[6];
    fpts[6] = NODES[3];
    fpts[7] = -1;

    // First child replaces quad
    fpts = Mesh_get_face(M, quad);
    fpts[0] = NODES[0];
    fpts[1] = NODES[4];
    fpts[2] = NODES[1];
    fpts[3] = -1;
    fpts[4] = NODES[5];
    fpts[5] = -1;
    fpts[6] = NODES[7];
    fpts[7] = -1;

    // Update ranges and strides
    Mesh_update_face_range_and_stride(M, quad, M->nf, 1);

    // Conformize parent cells

    E_Int own = M->owner[quad];

    assert(M->clevel[own] == M->flevel[quad]);

    if (Mesh_conformize_cell_face(M, own, quad, M->nf, 2) != 0) return 1;

    E_Int nei = M->neigh[quad];

    if (nei != -1) {
        assert(M->clevel[nei] == M->flevel[quad]);

        if (Mesh_conformize_cell_face(M, nei, quad, M->nf, 2) != 0) return 1;
    }

    for (E_Int i = 0; i < 1; i++) {
        M->owner[M->nf+i] = own;
        M->neigh[M->nf+i] = nei;
    }

    // Update adaptation info
    M->flevel[quad]++;

    assert(M->fref[quad] == FACE_REFINED);

    for (E_Int i = 0; i < 1; i++) {
        E_Int fid = M->nf + i;

        M->flevel[fid] = M->flevel[quad];
        M->ftype[fid] = M->ftype[quad];

        M->fref[fid] = FACE_NEW;
    }

    M->fchildren[quad] = {quad, M->nf};

    M->fparent[M->nf] = quad;

    M->ftag[M->nf] = M->ftag[quad];

    M->fpattern[M->nf] = M->fpattern[quad];

    // Increment face/edge/point count
    M->nf += 1;

    return 0;
}

E_Int refine_quad_dir(E_Int fid, Mesh *M)
{
    if (M->fpattern[fid] == DIR_X) return refine_quad_X(fid, M);
    if (M->fpattern[fid] == DIR_Y) return refine_quad_Y(fid, M);
    return Q9_refine(fid, M);
}

void Mesh_refine_dir(Mesh *M, ArrayI *ref_cells, ArrayI *ref_faces)
{
    std::set<E_Int> levelset;

    for (E_Int i = 0; i < ref_cells->count; i++) {
        levelset.insert(M->clevel[ref_cells->ptr[i]]);
    }

    for (E_Int i = 0; i < ref_faces->count; i++) {
        levelset.insert(M->flevel[ref_faces->ptr[i]]);
    }

    std::vector<E_Int> levels;
    for (E_Int level : levelset) levels.push_back(level);
    std::sort(levels.begin(), levels.end());

    //std::reverse(ref_cells->ptr, ref_cells->ptr + ref_cells->count);
    //std::reverse(ref_faces->ptr, ref_faces->ptr + ref_faces->count);

    E_Int cells_left = ref_cells->count-1;
    E_Int faces_left = ref_faces->count-1;

    for (E_Int level : levels) {
        while (faces_left >= 0 && M->flevel[ref_faces->ptr[faces_left]] == level) {
            E_Int face = ref_faces->ptr[faces_left];
            faces_left--;
            assert(M->fpattern[face] != -1);
            refine_quad_dir(face, M);
        }

        while (cells_left >= 0 && M->clevel[ref_cells->ptr[cells_left]] == level) {
            E_Int cell = ref_cells->ptr[cells_left];
            cells_left--;
            refine_cell_dir(cell, M);
        }
    }
}

#include "common/Karray.h"

PyObject *K_XCORE::AdaptMesh_AdaptGeom(PyObject *self, PyObject *args)
{
    PyObject *AMESH, *SLAVE, *TAGGED_FACES;

    if (!PYPARSETUPLE_(args, OOO_, &AMESH, &SLAVE, &TAGGED_FACES)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(AMESH, "AdaptMesh")) {
        RAISE("Bad first AdaptMesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(AMESH, "AdaptMesh");

    if (M->npc > 1) {
        RAISE("AdaptGeom is sequential.");
        return NULL;
    }

    Karray sarray;

    E_Int ret;

    ret = Karray_parse_ngon(SLAVE, sarray);

    if (ret != 0) {
        RAISE("Bad slave mesh.");
        return NULL;
    }

    puts("Preparing meshes for intersection...");

    // Check slave point tags
    E_Int *tagged_faces = NULL;
    E_Int tag_size = -1;
    ret = K_NUMPY::getFromNumpyArray(TAGGED_FACES, tagged_faces, tag_size, true);
    if (ret != 1) {
        Karray_free_ngon(sarray);
        RAISE("Bad slave points tag.");
        return NULL;
    }
 
    // Init slave mesh
    DynMesh S(&sarray);
    S.orient_skin(OUT);
    for (E_Int i = 0; i < tag_size; i++) {
        S.ftag[tagged_faces[i]] = 1;
    }

    S.triangulate(tagged_faces, tag_size);

    //S = S.extract_conformized();
    //return S.export_karray();

    // Refine M volumetric wrt to S tagged faces point cloud
    // Refine S surfacic wrt to M tagged faces point cloud

    E_Int iter = 0;
    E_Int max_iter = 10;
    E_Int ref_count_M, ref_count_S;

    do {
        iter++;
        
        /************************** M refinement **************************/
    
        printf("\niter: %d\n", iter);

        // Extract spoints from tagged sfaces
        ArrayI spoints;
        S.extract_points_from_ftag(&spoints);
        
        // We need the skin connectivity graph
        SkinGraph skin_graph = {0};
        Mesh_make_skin_graph(M, &skin_graph);
        printf("Skin: %d faces\n", skin_graph.nf);

        // BVH the skin
        Vec3f *skin_fc = (Vec3f *)XMALLOC(skin_graph.nf * sizeof(Vec3f));
        Mesh_make_face_centers(M, skin_graph.nf, skin_graph.skin, skin_fc);
        E_Int *indices = (E_Int *)XMALLOC(skin_graph.nf * sizeof(E_Int));
        for (E_Int i = 0; i < skin_graph.nf; i++) indices[i] = i;
        const Box3 huge = {-FLT_MAX, -FLT_MAX, -FLT_MAX,
                            FLT_MAX,  FLT_MAX,  FLT_MAX};
        BVH_node *bvh = BVH_make(M, skin_graph.skin, skin_fc, indices, 0,
            skin_graph.nf, &huge);
        puts("BVH constructed");

        
        // Locate spoints in skin
        PointFaces *sploc =
            (PointFaces *)XMALLOC(spoints.count * sizeof(PointFaces));
        memset(sploc, 0, spoints.count * sizeof(PointFaces));
        locate_points_in_skin
        (
            &S, &spoints,
            M, skin_graph.skin, indices,
            bvh,
            sploc
        );
        puts("Points located");

        if (iter == 1) {
            // We need the skin face normals
            Mesh_SkinGraph_compute_normals(M, &skin_graph);

            E_Int *xadj = skin_graph.xadj;
            E_Int *fnei = skin_graph.fnei;
            Vec3f *fnml = skin_graph.fnml;

            // Traverse the skin in breadth-first fashion
            bool *visited = (bool *)XMALLOC(skin_graph.nf * sizeof(bool));
            memset(visited, 0, skin_graph.nf * sizeof(bool));

            E_Int *fqueue = (E_Int *)XMALLOC(skin_graph.nf * sizeof(E_Int));
            memset(fqueue, -1, skin_graph.nf * sizeof(E_Int));

            // Start from the first mface
            E_Int fseed = indices[sploc[0].ptr[0]];

            E_Int front = 0, rear = 0;
            fqueue[rear++] = fseed;

            visited[fseed] = true;
            E_Int nvisited = 1;
            
            while (front != rear) {
                E_Int fid = fqueue[front++];

                E_Int start = xadj[fid];
                E_Int end = xadj[fid+1];
                E_Int stride = end - start;
                E_Int *pn = &fnei[start];

                Vec3f *fid_nml = &fnml[fid];

                for (E_Int i = 0; i < stride; i++) {
                    E_Int nei = pn[i];
                    if (visited[nei]) continue;

                    Vec3f *nei_nml = &fnml[nei];
                    E_Float dp = K_MATH::dot((E_Float *)fid_nml, (E_Float *)nei_nml, 3);
                    //assert(dp >= 0.0);
                    dp = std::max(dp, -1.0);
                    dp = std::min(dp, 1.0);
                    E_Float theta = acos(dp) * 180.0 / K_MATH::PI;
            
                    if (theta < 60.0 || theta > 120.0) {
                        visited[nei] = true;
                        fqueue[rear++] = nei;
                        nvisited++;
                    }
                }
            }

            // Set the base patch
            E_Int cvisit = 0;
            for (E_Int i = 0; i < skin_graph.nf; i++) {
                if (visited[i] == 0) continue;

                // Tag this face for later
                assert(M->ftag[skin_graph.skin[i]] == 0);
                M->ftag[skin_graph.skin[i]] = 1;

                E_Int fid = skin_graph.skin[i];
                E_Int cid = M->owner[fid];

                while (cid != -1) {
                    cvisit++;
                    Mesh_set_face_as_cell_bottom(M, fid, cid);
                    H18_reorder(cid, M);
                    assert(check_canon_hexa(cid, M) == 0);

                    E_Int *cell = Mesh_get_cell(M, cid);

                    E_Int map[4];
                    E_Int *face = Mesh_get_face(M, fid);
                    for (E_Int j = 0; j < 4; j++) map[j] = face[2*j];
                    E_Int reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[0]);
                    if (reorient) std::swap(map[1], map[3]); 

                    for (E_Int j = 0; j < 6; j++) {
                        fid = cell[4*j];
                        
                        if (M->fpattern[fid] != -1) continue;

                        if (j == 0 || j == 1) {
                            M->fpattern[fid] = DIR_ISO;
                        } else {
                            
                            E_Int i0;
                            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[j]);
                            face = Mesh_get_face(M, fid);

                            M->fpattern[fid] = DIR_X;
                                
                            if (j == 2) {
                                i0 = Get_pos(map[0], face, 8);                    
                            } else if (j == 3) {
                                i0 = Get_pos(map[1], face, 8);
                            } else if (j == 4) {
                                i0 = Get_pos(map[1], face, 8);
                            } else if (j == 5) {
                                i0 = Get_pos(map[2], face, 8);
                            }

                            assert(i0 != -1);
                            i0 /= 2;

                            if ((reorient == 0 && (i0 == 1 || i0 == 3)) ||
                                (reorient == 1 && (i0 == 0 || i0 == 2))) {
                                M->fpattern[fid] = DIR_Y;
                            }
                        }
                    }

                    // Step to top adjacent cell
                    fid = cell[4];
                    cid = Mesh_get_cnei(M, cid, fid);
                }
            }

            assert(cvisit == M->nc);

            XFREE(visited);
            XFREE(fqueue);
        }

        // Isolate faces that contain more than MAX_POINTS_PER_FACE spoints
        ArrayI rfaces;
        PointFaces_extract_by_threshold
        (
            sploc, spoints.count,
            skin_graph.skin, skin_graph.nf,
            1, // threshold
            &rfaces
        );
        puts("Refinement faces isolated");
        printf("Refinement faces: %d\n", rfaces.count);
        ref_count_M = rfaces.count;
        
        if (rfaces.count > 0) {
    
            // Smooth face refinement data
            E_Int *fdat = (E_Int *)XMALLOC(skin_graph.nf * sizeof(E_Int));
            memset(fdat, 0, skin_graph.nf * sizeof(E_Int));
            for (E_Int i = 0; i < rfaces.count; i++) {
                E_Int idx_in_skin = indices[rfaces.ptr[i]];
                fdat[idx_in_skin] = 1;
            }
            SkinGraph_smooth_ref_data(&skin_graph, fdat, M);
            puts("Face refinement smoothed out");
            
            E_Int smooth_nfref = 0;
            for (E_Int i = 0; i < skin_graph.nf; i++) {
                if (fdat[i] > 0) smooth_nfref++;
            }
            printf("Smooth refinement face count: %d\n", smooth_nfref);

            // Allocate
            M->cref = (E_Int *)XRESIZE(M->cref, M->nc * sizeof(E_Int));
            memset(M->cref, 0, M->nc * sizeof(E_Int));

            // Cells
            E_Int ref_cell_count;
            init_skin_refinement_cells(&skin_graph, fdat, M, &ref_cell_count);
            printf("Refinement cells: %d\n", ref_cell_count);
            if (ref_cell_count == 0) break;

            // Smooth cell refinement data
            Mesh_smooth_cell_refinement_data(M);
            puts("Cell refinement smoothed out");

            // Assign refinement data
            M->fref = (E_Int *)XRESIZE(M->fref, M->nf * sizeof(E_Int));
            memset(M->fref, 0, M->nf * sizeof(E_Int));
            assign_face_refinement_data(M);

            // Isolate cells/faces to be refined
            ArrayI ref_cells, ref_faces;
            Mesh_isolate_refinement_entities(M, &ref_cells, &ref_faces);
            printf("Refinement cells: %d\n", ref_cells.count);
            printf("Refinement faces: %d\n", ref_faces.count);
            
            // Resize for refinement
            Mesh_resize(M, &ref_cells, &ref_faces);
            puts("Mesh resized for refinement");

            // Sort entities by refinement level
            std::sort(ref_cells.ptr, ref_cells.ptr + ref_cells.count,
                [&] (E_Int i, E_Int j) { return M->clevel[i] > M->clevel[j]; });
            puts("Refinement cells sorted");
            std::sort(ref_faces.ptr, ref_faces.ptr + ref_faces.count,
                [&] (E_Int i, E_Int j) { return M->flevel[i] > M->flevel[j]; });
            puts("Refinement faces sorted");
            
            // Refine
            printf("Cells before refinement: %d\n", M->nc);
            printf("Faces before refinement: %d\n", M->nf);
            Mesh_refine_dir(M, &ref_cells, &ref_faces);

            printf("Cells after refinement: %d\n", M->nc);
            printf("Faces after refinement: %d\n", M->nf);
        
            Mesh_conformize_face_edge(M);

            ArrayI_free(&ref_cells);
            ArrayI_free(&ref_faces);
            XFREE(fdat);
        }

        BVH_free(bvh);
        ArrayI_free(&rfaces);
        ArrayI_free(&spoints);
        
        XFREE(sploc);
        XFREE(indices);
        XFREE(skin_fc);
        SkinGraph_free(&skin_graph);

        /************************** S refinement **************************/

        // Extract mpoints from tagged mfaces
        ArrayI mpoints;
        Mesh_extract_points_from_ftag(M, &mpoints);
        printf("mpoints: %d\n", mpoints.count);
        
        /*if (iter == 10)
        {
            std::vector<Point> points;
            for (E_Int i = 0; i < mpoints.count; i++) {
                E_Int p = mpoints.ptr[i];
                Point point = {M->X[p], M->Y[p], M->Z[p]};
                points.push_back(point);
            }
            points_write("mpoints", points);
        }
        */

        // We need the skin connectivity graph
        S.make_tri_graph();
        printf("Skin: %lu faces\n", S.tri_graph.nf);

        // BVH the skin
        skin_fc = (Vec3f *)XCALLOC(S.tri_graph.nf, sizeof(Vec3f));
        S.make_face_centers(S.tri_graph.nf, S.tri_graph.skin.data(), skin_fc);
        indices = (E_Int *)XMALLOC(S.tri_graph.nf * sizeof(E_Int));
        for (size_t i = 0; i < S.tri_graph.nf; i++) indices[i] = i;
        bvh = BVH_make(&S, S.tri_graph.skin.data(), skin_fc, indices, 0,
            S.tri_graph.nf, &huge);
        puts("BVH constructed");

        // Locate mpoints in skin
        PointFaces *mploc =
            (PointFaces *)XCALLOC(mpoints.count, sizeof(PointFaces));
        locate_points_in_skin
        (
            M, &mpoints,
            &S, S.tri_graph.skin.data(), indices,
            bvh,
            mploc
        );
        puts("Points located");

        // Isolate faces that contain more than MAX_POINTS_PER_FACE mpoints
        PointFaces_extract_by_threshold
        (
            mploc, mpoints.count,
            S.tri_graph.skin.data(), S.tri_graph.nf,
            1, // threshold
            &rfaces
        );
        puts("Refinement faces isolated");
        printf("Refinement faces: %d\n", rfaces.count);
        ref_count_S = rfaces.count;

        /*
        if (iter == 2)
        {
            npy_intp dims[2];

            // faces array
            dims[1] = 1;
            dims[0] = (npy_intp)rfaces.count;
        
            PyArrayObject *FACES = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
            E_Int *pf = (E_Int *)PyArray_DATA(FACES);
            E_Int *ptr = pf;

            for (E_Int j = 0; j < rfaces.count; j++) {
                *ptr++ = skin_graph.skin[indices[rfaces.ptr[j]]]+1;
            }

            return (PyObject *)FACES;
        }
        */

        if (rfaces.count > 0) {
    
            // Smooth face refinement data
            S.tri_graph.fdat.clear();
            S.tri_graph.fdat.resize(S.tri_graph.nf, 0);
            for (E_Int i = 0; i < rfaces.count; i++) {
                E_Int idx_in_skin = indices[rfaces.ptr[i]];
                S.tri_graph.fdat[idx_in_skin] = 1;
            }
            
            // Isolate refinement faces
            ArrayI ref_faces;
            S.prepare_for_refinement(&ref_faces);
            printf("Refinement faces: %d\n", ref_faces.count);

            // Refine
            S.refine_faces(&ref_faces);

            // Make sure to copy ftag data
            S = S.extract_conformized();

            return S.export_karray();
            
            ArrayI_free(&ref_faces);
        }

        // FREE
        S.tri_graph.clear();
        S.fchildren.clear();
        S.ecenter.clear();
        BVH_free(bvh);
        ArrayI_free(&rfaces);
        ArrayI_free(&mpoints);
        XFREE(mploc);
        XFREE(indices);
        XFREE(skin_fc);
        SkinGraph_free(&skin_graph);

    } while (iter < max_iter && (ref_count_M > 0 || ref_count_S > 0));


    PyObject *Sout = S.export_karray();

    Karray_free_ngon(sarray);

    return Sout;

    return Py_None;
}
