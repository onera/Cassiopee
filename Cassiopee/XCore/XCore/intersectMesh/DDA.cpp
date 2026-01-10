/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
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
    bin_faces.clear();

    for (E_Int fid : patch) {
        assert(face_is_active(fid));
    
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
