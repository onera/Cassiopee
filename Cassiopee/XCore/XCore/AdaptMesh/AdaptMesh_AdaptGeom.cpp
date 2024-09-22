#include "Mesh.h"
#include "common/mem.h"

void point_print(E_Float px, E_Float py, E_Float pz)
{
    FILE *fh = fopen("point", "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "1\n");
    fprintf(fh, "%f %f %f\n", px, py, pz);
    fclose(fh);
}

struct ArrayI {
    E_Int count;
    E_Int *ptr;
};

void Mesh_extract_faces_from_ftag(Mesh *M, ArrayI *fids)
{
    fids->count = 0;

    for (E_Int i = 0; i < M->nf; i++) {
        fids->count += (M->ftag[i] == 1);
    }

    fids->ptr = (E_Int *)XMALLOC(fids->count * sizeof(E_Int));
    E_Int *ptr = fids->ptr;

    for (E_Int i = 0; i < M->nf; i++) {
        if (M->ftag[i] == 1)
            *ptr++ = i;
    }
}

void Mesh_extract_points_from_ftag(Mesh *M, ArrayI *pids)
{
    // WARNING(Imad): ptag is reset here
    memset(M->ptag, 0, M->np * sizeof(E_Int));
    pids->count = 0;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        if (M->ftag[fid] != 1) continue;
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        for (E_Int i = 0; i < M->fstride[fid]; i++) {
            E_Int *pn = face + 2*i;
            for (E_Int j = 0; j < frange[i]; j++) {
                E_Int pid = pn[j];
                pids->count += (M->ptag[pid] == 0);
                M->ptag[pid] = 1;
            }
        }
    }

    pids->ptr = (E_Int *)XMALLOC(pids->count * sizeof(E_Int));
    E_Int *ptr = pids->ptr;

    for (E_Int pid = 0; pid < M->np; pid++) {
        if (M->ptag[pid] == 1)
            *ptr++ = pid;
    }
}

struct Box2D {
    E_Float xmin, ymin;
    E_Float xmax, ymax;
};

Box2D Mesh_box2D_faces(Mesh *M, ArrayI *fids)
{
    E_Float xmin, ymin, xmax, ymax;
    xmin = ymin = 0xFFFFFFFF;
    xmax = ymax = -xmin;

    for (E_Int i = 0; i < fids->count; i++) {
        E_Int fid = fids->ptr[i];
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            E_Int *pn = face + 2*j;
            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int pid = pn[k];
                if (M->X[pid] < xmin) xmin = M->X[pid];
                if (M->Y[pid] < ymin) ymin = M->Y[pid];
                if (M->X[pid] > xmax) xmax = M->X[pid];
                if (M->Y[pid] > ymax) ymax = M->Y[pid];
            }
        }
    }

    // Safety
    xmax = xmax + (xmax - xmin)*0.01;
    ymax = ymax + (ymax - ymin)*0.01;

    return { xmin, ymin, xmax, ymax };
}

#define MAX_FACES_PER_GRID 8

struct GridFaces {
    E_Int ptr[MAX_FACES_PER_GRID];
    E_Int count;
};

#define GRIDX 64
#define GRIDY 64

void Mesh_2D_grid_faces(Mesh *M, ArrayI *mfaces, Box2D box2D,
    GridFaces grid_faces[GRIDX][GRIDY])
{
    E_Float HX = (box2D.xmax - box2D.xmin) / GRIDX;
    E_Float HY = (box2D.ymax - box2D.ymin) / GRIDY;

    for (E_Int i = 0; i < mfaces->count; i++) {
        E_Int fid = mfaces->ptr[i];
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);

        E_Int minX = GRIDX, maxX = -1;
        E_Int minY = GRIDY, maxY = -1;
        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            E_Int *pn = face + 2*j;
            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int pid = pn[k];
                E_Int grid_x = (M->X[pid] - box2D.xmin) / HX;
                E_Int grid_y = (M->Y[pid] - box2D.ymin) / HY;
                assert(grid_x >= 0 && grid_x < GRIDX);
                assert(grid_y >= 0 && grid_y < GRIDY);
                if (grid_x < minX) minX = grid_x;
                if (grid_y < minY) minY = grid_y;
                if (grid_x > maxX) maxX = grid_x;
                if (grid_y > maxY) maxY = grid_y;
            }
        }

        for (E_Int grid_x = minX; grid_x <= maxX; grid_x++) {
            for (E_Int grid_y = minY; grid_y <= maxY; grid_y++) {
                GridFaces *grid = &grid_faces[grid_x][grid_y];
                if (grid->count >= MAX_FACES_PER_GRID) {
                    fprintf(stderr, "Grid is too coarse, refine it!\n");
                    abort();
                }
                grid->ptr[grid->count++] = i; // We need the index in mfaces, not the id
            }
        }
    }
}

struct TriBary {
    E_Float UX, UY;
    E_Float VX, VY;
    E_Float UU, VV, UV;
    E_Float inv_denom;
};

void Mesh_locate_spoints_in_mtris
(
    Mesh *S, ArrayI *spoints,
    Mesh *M, ArrayI *mtris,
    Box2D box2D, GridFaces grid_faces[GRIDX][GRIDY], E_Int *sploc)
{
    // Compute the mtris dot products once
    TriBary *tri_barys = (TriBary *)XMALLOC(mtris->count * sizeof(TriBary));
    for (E_Int i = 0; i < mtris->count; i++) {
        E_Int tid = mtris->ptr[i];
        assert(M->ftype[tid] == TRI);
        E_Int *tri = Mesh_get_face(M, tid);
        E_Int A = tri[0], B = tri[2], C = tri[4];

        E_Float ux = (M->X[B] - M->X[A]);
        E_Float uy = (M->Y[B] - M->Y[A]);
        E_Float vx = (M->X[C] - M->X[A]);
        E_Float vy = (M->Y[C] - M->Y[A]);

        tri_barys[i].UX = ux;
        tri_barys[i].UY = uy;

        tri_barys[i].VX = vx;
        tri_barys[i].VY = vy;

        tri_barys[i].UU = ux*ux + uy*uy;
        tri_barys[i].VV = vx*vx + vy*vy;
        tri_barys[i].UV = ux*vx + uy*vy;

        tri_barys[i].inv_denom = 1.0 / (tri_barys[i].UU * tri_barys[i].VV -
                                        tri_barys[i].UV * tri_barys[i].UV);
    }

    E_Float HX = (box2D.xmax - box2D.xmin) / GRIDX;
    E_Float HY = (box2D.ymax - box2D.ymin) / GRIDY;

    memset(sploc, -1, mtris->count * sizeof(E_Int));

    for (E_Int i = 0; i < spoints->count; i++) {
        E_Int spid = spoints->ptr[i];
        E_Float px = S->X[spid];
        E_Float py = S->Y[spid];
        E_Int grid_x = (px - box2D.xmin) / HX;
        E_Int grid_y = (py - box2D.ymin) / HY;
        assert(grid_x >= 0);
        assert(grid_x < GRIDX);
        assert(grid_y >= 0);
        assert(grid_y < GRIDY);
        GridFaces *grid = &grid_faces[grid_x][grid_y];

        if (grid->count == 0) {
            fprintf(stderr, "No faces found in point grid!\n");
            abort();
        };

        for (E_Int j = 0; j < grid->count; j++) {
            E_Int tidx = grid->ptr[j];
            E_Int tid = mtris->ptr[tidx];
            assert(M->ftype[tid] == TRI);

            E_Int *tri = Mesh_get_face(M, tid);
            E_Int A = tri[0];
            E_Float DX = px - M->X[A];
            E_Float DY = py - M->Y[A];

            TriBary *tb = &tri_barys[tidx];

            E_Float d20 = DX*tb->UX + DY*tb->UY;
            E_Float d21 = DX*tb->VX + DY*tb->VY;

            E_Float u = (tb->VV*d20 - tb->UV*d21) * tb->inv_denom;
            E_Float v = (tb->UU*d21 - tb->UV*d20) * tb->inv_denom;

            if (u >= 0.0 && v >= 0.0 && (u + v) <= 1.0) {
                sploc[i] = tid;
                break;
            }
        }

        if (sploc[i] == -1) {
            fprintf(stderr, "Failed to locate point %d!\n", spid);
            abort();
        }
    }



}

PyObject *K_XCORE::AdaptMesh_AdaptGeom(PyObject *self, PyObject *args)
{
    PyObject *AMESH, *SMESH;

    if (!PYPARSETUPLE_(args, OO_, &AMESH, &SMESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(AMESH, "AdaptMesh")) {
        RAISE("Bad first AdaptMesh hook.");
        return NULL;
    }

    if (!PyCapsule_IsValid(SMESH, "AdaptMesh")) {
        RAISE("Bad second AdaptMesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(AMESH, "AdaptMesh");
    Mesh *S = (Mesh *)PyCapsule_GetPointer(SMESH, "AdaptMesh");

    // Refine M volumetric wrt to S tagged faces point cloud
    // Refine S surfacic wrt to M tagged faces point cloud

    E_Int ref_M = 0, ref_S = 0;
    E_Int iter = 0;

    do {
        iter++;

        // Sfaces are tagged
        ArrayI spoints;
        Mesh_extract_points_from_ftag(S, &spoints);
        
        // Spoints must belong to one of the tagged Mfaces
        // All the tagged faces must be triangles
        ArrayI mtris;
        Mesh_extract_faces_from_ftag(M, &mtris);

        // Box the mfaces
        Box2D box2D = Mesh_box2D_faces(M, &mtris);
        
        // Grid the mfaces
        GridFaces grid_faces[GRIDX][GRIDY] = {0};
        Mesh_2D_grid_faces(M, &mtris, box2D, grid_faces);

        // Locate spoints in mfaces
        E_Int *sploc = (E_Int *)XMALLOC(spoints.count * sizeof(E_Int));
        Mesh_locate_spoints_in_mtris
        (
            S, &spoints,
            M, &mtris,
            box2D, grid_faces,
            sploc
        );





        XFREE(spoints.ptr);
        XFREE(mtris.ptr);
        XFREE(sploc);
    } while (ref_S);

    return Py_None;
}
