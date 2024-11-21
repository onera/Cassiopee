#include "Box.h"
#include "Mesh.h"
#include "FaceSort.h"
#include "DynMesh.h"

Box3 Box3_make
(
    const Mesh *M,
    const FaceSort *mfaces,
    E_Int start, E_Int end
)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = FLT_MAX;
    xmax = ymax = zmax = -FLT_MAX;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = mfaces[i].fid;
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        assert(M->fstride[fid] == 3);
        assert(M->ftype[fid] == TRI);
        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            E_Int *pn = face + 2*j;
            assert(frange[j] == 1);
            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int pid = pn[k];
                if (M->X[pid] < xmin) xmin = M->X[pid];
                if (M->Y[pid] < ymin) ymin = M->Y[pid];
                if (M->Z[pid] < zmin) zmin = M->Z[pid];

                if (M->X[pid] > xmax) xmax = M->X[pid];
                if (M->Y[pid] > ymax) ymax = M->Y[pid];
                if (M->Z[pid] > zmax) zmax = M->Z[pid];
            }
        }
    }

    // Safety
    E_Float dx = (xmax - xmin) * 0.01;
    E_Float dy = (ymax - ymin) * 0.01;
    E_Float dz = (zmax - zmin) * 0.01;
    xmin -= dx;
    ymin -= dy;
    zmin -= dz;
    xmax += dx;
    ymax += dy;
    zmax += dz;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}

void Box3_clamp(const Box3 *parent, Box3 *child)
{
    child->xmin = std::max(parent->xmin, child->xmin);
    child->ymin = std::max(parent->ymin, child->ymin);
    child->zmin = std::max(parent->zmin, child->zmin);
    child->xmax = std::min(parent->xmax, child->xmax);
    child->ymax = std::min(parent->ymax, child->ymax);
    child->zmax = std::min(parent->zmax, child->zmax);
}

Box3 Box3_make
(
    const Mesh *M,
    const E_Int *skin,
    const E_Int *indices,
    E_Int start, E_Int end
)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = FLT_MAX;
    xmax = ymax = zmax = -FLT_MAX;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = skin[indices[i]];
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);

        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            E_Int *pn = face + 2*j;
            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int pid = pn[k];
                if (M->X[pid] < xmin) xmin = M->X[pid];
                if (M->Y[pid] < ymin) ymin = M->Y[pid];
                if (M->Z[pid] < zmin) zmin = M->Z[pid];

                if (M->X[pid] > xmax) xmax = M->X[pid];
                if (M->Y[pid] > ymax) ymax = M->Y[pid];
                if (M->Z[pid] > zmax) zmax = M->Z[pid];
            }
        }
    }

    // Safety
    E_Float dx = (xmax - xmin) * 0.01;
    E_Float dy = (ymax - ymin) * 0.01;
    E_Float dz = (zmax - zmin) * 0.01;
    xmin -= dx;
    ymin -= dy;
    zmin -= dz;
    xmax += dx;
    ymax += dy;
    zmax += dz;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}

/* DynMesh */

Box3 Box3_make
(
    const DynMesh *M,
    const FaceSort *mfaces,
    E_Int start, E_Int end
)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = FLT_MAX;
    xmax = ymax = zmax = -FLT_MAX;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = mfaces[i].fid;
        const auto &pn = M->F[fid];
        
        for (E_Int pid : pn) {
            if (M->X[pid] < xmin) xmin = M->X[pid];
            if (M->Y[pid] < ymin) ymin = M->Y[pid];
            if (M->Z[pid] < zmin) zmin = M->Z[pid];
            if (M->X[pid] > xmax) xmax = M->X[pid];
            if (M->Y[pid] > ymax) ymax = M->Y[pid];
            if (M->Z[pid] > zmax) zmax = M->Z[pid];
        }
    }

    // Safety
    E_Float dx = (xmax - xmin) * 0.01;
    E_Float dy = (ymax - ymin) * 0.01;
    E_Float dz = (zmax - zmin) * 0.01;
    xmin -= dx;
    ymin -= dy;
    zmin -= dz;
    xmax += dx;
    ymax += dy;
    zmax += dz;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}

Box3 Box3_make
(
    const DynMesh *M,
    const E_Int *skin,
    const E_Int *indices,
    E_Int start, E_Int end
)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = FLT_MAX;
    xmax = ymax = zmax = -FLT_MAX;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = skin[indices[i]];
        const auto &pn = M->F[fid];

        for (E_Int pid : pn) {
            if (M->X[pid] < xmin) xmin = M->X[pid];
            if (M->Y[pid] < ymin) ymin = M->Y[pid];
            if (M->Z[pid] < zmin) zmin = M->Z[pid];
            if (M->X[pid] > xmax) xmax = M->X[pid];
            if (M->Y[pid] > ymax) ymax = M->Y[pid];
            if (M->Z[pid] > zmax) zmax = M->Z[pid];
        }
    }

    // Safety
    E_Float dx = (xmax - xmin) * 0.01;
    E_Float dy = (ymax - ymin) * 0.01;
    E_Float dz = (zmax - zmin) * 0.01;
    xmin -= dx;
    ymin -= dy;
    zmin -= dz;
    xmax += dx;
    ymax += dy;
    zmax += dz;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}