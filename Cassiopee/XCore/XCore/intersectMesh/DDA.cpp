#include "mesh.h"
#include "primitives.h"

AABB IMesh::AABB_face(const std::vector<E_Int> &pn) const
{
    AABB ret;
    for (E_Int p : pn) {
        if (X[p] > ret.xmax) ret.xmax = X[p];
        if (X[p] < ret.xmin) ret.xmin = X[p];
        if (Y[p] > ret.ymax) ret.ymax = Y[p];
        if (Y[p] < ret.ymin) ret.ymin = Y[p];
        if (Z[p] > ret.zmax) ret.zmax = Z[p];
        if (Z[p] < ret.zmin) ret.zmin = Z[p];
    }
    return ret;
}

void IMesh::hash_skin()
{
    NX = 100;
    NY = 100;
    NZ = 100;

    HX = (xmax - xmin) / NX;
    HY = (ymax - ymin) / NY;
    HZ = (zmax - zmin) / NZ;

    E_Int NXY = NX * NY;

    assert(bin_faces.empty());

    for (E_Int fid : skin) {
        const auto &pn = F[fid];

        AABB bbox = AABB_face(pn);

        E_Int imin = floor((bbox.xmin - xmin) / HX);
        E_Int imax = floor((bbox.xmax - xmin) / HX);
        E_Int jmin = floor((bbox.ymin - ymin) / HY);
        E_Int jmax = floor((bbox.ymax - ymin) / HY);
        E_Int kmin = floor((bbox.zmin - zmin) / HZ);
        E_Int kmax = floor((bbox.zmax - zmin) / HZ);

        for (E_Int k = kmin; k < kmax+1; k++) {
            for (E_Int j = jmin; j < jmax+1; j++) {
                for (E_Int i = imin; i < imax+1; i++) {
                    E_Int voxel = i + NX*j + NXY*k;
                    bin_faces[voxel].push_back(fid);
                }
            }
        }
    }
}