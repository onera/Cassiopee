#include "Point.h"
#include "constants.h"
#include "FaceSort.h"
#include "Box.h"
#include "Array.h"
#include "common/mem.h"

bool Point_in_tri(E_Float px, E_Float py, E_Float pz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz)
{
    // Normal vector to the plane
    E_Float Y[3] = {bx-ax, by-ay, bz-az};
    E_Float Z[3] = {cx-ax, cy-ay, cz-az};
    E_Float N[3];
    K_MATH::cross(Y, Z, N);

    E_Float X[3] = {px-ax, py-ay, pz-az};

    E_Float dp = K_MATH::dot(N, X, 3);
    
    // Is the point on the plane?
    if (dp < -TOL || dp > TOL) return 0;

    E_Float x1 = K_MATH::dot(X, Y, 3);
    E_Float y1 = K_MATH::dot(Y, Y, 3);
    E_Float z1 = K_MATH::dot(Z, Y, 3);
    E_Float x2 = K_MATH::dot(X, Z, 3);
    E_Float y2 = K_MATH::dot(Y, Z, 3);
    E_Float z2 = K_MATH::dot(Z, Z, 3);

    E_Float u = (x1*z2 - x2*z1) / (y1*z2 - y2*z1);
    if (u < -TOL || u > 1 + TOL) return false;

    E_Float v = (-x1*y2 + x2*y1) / (y1*z2 - y2*z1);
    if (v < -TOL || v > 1 + TOL) return false;

    E_Float w = 1 - u - v;
    if (w < -TOL || w > 1 + TOL) return false;

    return true;
}

bool Point_in_FaceSort(const Point *p, const FaceSort *f)
{
    E_Float DX = p->x - f->xa;
    E_Float DY = p->y - f->ya;
    E_Float DZ = p->z - f->za;
    
    E_Float d20 = DX*f->UX + DY*f->UY + DZ*f->UZ;
    E_Float d21 = DX*f->VX + DY*f->VY + DZ*f->VZ;

    E_Float u = (f->VV*d20 - f->UV*d21) * f->inv_denom;
    if (u < -TOL || u > 1.0 + TOL) return false;

    E_Float v = (f->UU*d21 - f->UV*d20) * f->inv_denom;
    if (v < -TOL || v > 1.0 + TOL) return false;

    E_Float w = 1.0 - (u + v);
    if (w < -TOL || w > 1.0 + TOL) return false;

    return true;
}

bool Point_in_Box3D(const Point *p, const Box3 *box)
{
    return (p->x >= box->xmin) && (p->x <= box->xmax) &&
           (p->y >= box->ymin) && (p->y <= box->ymax) &&
           (p->z >= box->zmin) && (p->z <= box->zmax);
}

void PointFaces_extract_by_threshold
(
    const PointFaces *sploc, E_Int spcount,
    const E_Int *skin, E_Int mcount,
    const E_Int threshold,
    ArrayI *faces
)
{
    E_Int *ftag = (E_Int *)XMALLOC(mcount * sizeof(E_Int));
    memset(ftag, 0, mcount * sizeof(E_Int));

    for (E_Int i = 0; i < spcount; i++) {
        const PointFaces *pfaces = &sploc[i];
        for (E_Int j = 0; j < pfaces->count; j++)
            ftag[pfaces->ptr[j]]++;
    }

    faces->count = 0;

    for (E_Int i = 0; i < mcount; i++) {
        if (ftag[i] > threshold)
            faces->count++;
    }

    faces->ptr = (E_Int *)XMALLOC(faces->count * sizeof(E_Int));
    E_Int *ptr = faces->ptr;

    for (E_Int i = 0; i < mcount; i++) {
        if (ftag[i] > threshold)
            *ptr++ = i;//skin[i];
    }

    XFREE(ftag);
}

void points_write(const char *fname, const std::vector<Point> &P)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", P.size());
    for (auto p : P) fprintf(fh, "%f %f %f\n", p.x, p.y, p.z);
    fclose(fh);
}

void point_write(E_Float px, E_Float py, E_Float pz)
{
    FILE *fh = fopen("point", "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "1\n");
    fprintf(fh, "%f %f %f\n", px, py, pz);
    fclose(fh);
}

void point_write(const Point *p)
{
    return point_write(p->x, p->y, p->z);
}