#include "mesh.h"
#include "primitives.h"

void IMesh::hash_skin()
{
    bin_faces.clear();

    bin_faces.resize(NXYZ);

    for (E_Int fid : skin) {
        const auto &pn = F[fid];

        E_Int Imin, Jmin, Kmin;
        E_Int Imax, Jmax, Kmax;

        Imin = Jmin = Kmin = NXYZ;
        Imax = Jmax = Kmax = -1;

        for (E_Int p : pn) {
            E_Float x = X[p];
            E_Float y = Y[p];
            E_Float z = Z[p];

            E_Int I = floor((x - xmin) / HX);
            E_Int J = floor((y - ymin) / HY);
            E_Int K = floor((z - zmin) / HZ);

            assert(I >= 0);
            assert(J >= 0);
            assert(K >= 0);

            assert(I < NX);
            assert(J < NY);
            assert(K < NZ);

            if (I < Imin) Imin = I;
            if (J < Jmin) Jmin = J;
            if (K < Kmin) Kmin = K;
            if (I > Imax) Imax = I;
            if (J > Jmax) Jmax = J;
            if (K > Kmax) Kmax = K;
        }

        for (E_Int I = Imin; I <= Imax; I++) {
            for (E_Int J = Jmin; J <= Jmax; J++) {
                for (E_Int K = Kmin; K <= Kmax; K++) {
                    E_Int voxel = get_voxel(I, J, K);
                    assert(voxel >= 0);
                    assert(voxel < NXYZ);
                    bin_faces[voxel].push_back(fid);
                }
            }
        }
    }
}

void IMesh::hash_patch()
{
    fmap.clear();

    for (E_Int fid : patch) {
        assert(face_is_active(fid));
    
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
                    fmap[voxel].push_back(fid);
                }
            }
        }
    }
}

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