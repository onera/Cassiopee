/*    
    Copyright 2013-2025 Onera.

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
#include "smesh.h"
#include "triangle.h"
#include "primitives.h"
#include "io.h"
#include "ray.h"

std::vector<PointLoc> Smesh::locate2(const Smesh &Sf) const
{
    std::vector<PointLoc> plocs(Sf.np);

    size_t on_vertex = 0, on_edge = 0;

    for (E_Int pid = 0; pid < Sf.np; pid++) {
        E_Float x = Sf.X[pid];
        E_Float y = Sf.Y[pid];
        E_Float z = Sf.Z[pid];

        E_Float dx = (E_Float) rand() / RAND_MAX;
        E_Float dy = (E_Float) rand() / RAND_MAX;
        E_Float dz = (E_Float) rand() / RAND_MAX; 
        Ray ray(x, y, z, dx, dy, dz, Ray::Policy::BOTH);

        HitData hit_data;
        intersect_ray(ray, 0, hit_data);

        /*
        if (pid == 326)
            intersect_ray(ray, 0, hit_data);
        else
            intersect_ray(ray, 0, hit_data);
        /*
        for (const auto &loc : hit_data.locs) {
            printf("%f ", loc.t);
        }
        puts("");

        if (pid == 326) {
            for (auto fid : hit_data.tested) {
                char name[64] = {0};
                sprintf(name, "fid%d.im", fid);
                write_face(name, fid);
            }
            double scale = 30;
            edge_write("ray.im", x, y, z, x+scale*dx, y+scale*dy, z+scale*dz);
            point_write("pid.im", x, y, z);
            int i = 0;
            for (const auto &loc : hit_data.locs) {
                printf("%f ", loc.t);
                char name[64] = {0};
                sprintf(name, "hit%d.im", i);
                point_write(name, loc.x, loc.y, loc.z);
                i++;
            }
            puts("");
            write_ngon("tested.im", hit_data.tested);
        }
        */

        if (hit_data.locs.empty()) {
            fprintf(stderr, "Could not locate point %d\n", pid);
            point_write("lost.im", x, y, z);
            //write_ngon("tested.im", hit_data.tested);
            assert(0);
        }

        auto &ploc = plocs[pid];
        E_Float t_abs_min = EFLOATMAX;
        for (const auto &loc : hit_data.locs) {
            if (fabs(loc.t) < t_abs_min) {
                t_abs_min = fabs(loc.t);
                ploc = loc;
            }
        }
    }

    /*
    std::vector<E_Int> fids;
    for (const auto &ploc : plocs) {
        E_Int fid = ploc.fid;
        fids.push_back(fid);
    }
    write_ngon("fids.im", fids);
    */

    return plocs;
}

void Smesh::make_bbox()
{
    NX = 100;
    NY = 100;
    NZ = 100;
    NXY = NX * NY;
    NXYZ = NXY * NZ;

    xmin = ymin = zmin = EFLOATMAX;
    xmax = ymax = zmax = EFLOATMIN;

    for (E_Int i = 0; i < np; i++) {
        if (X[i] < xmin) xmin = X[i];
        if (Y[i] < ymin) ymin = Y[i];
        if (Z[i] < zmin) zmin = Z[i];
        if (X[i] > xmax) xmax = X[i];
        if (Y[i] > ymax) ymax = Y[i];
        if (Z[i] > zmax) zmax = Z[i];
    }

    E_Float dx = xmax - xmin;
    E_Float dy = ymax - ymin;
    E_Float dz = zmax - zmin;

    //xmin = xmin - dx*0.01;
    //ymin = ymin - dy*0.01;
    //zmin = zmin - dz*0.01;
    //xmax = xmax + dx*0.01;
    //ymax = ymax + dy*0.01;
    //zmax = zmax + dz*0.01;

    //HX = (dx != 0) ? (xmax - xmin) / NX : 1;
    //HY = (dy != 0) ? (ymax - ymin) / NY : 1;
    //HZ = (dz != 0) ? (zmax - zmin) / NZ : 1;

    box.xmin = xmin, box.ymin = ymin, box.zmin = zmin;
    box.xmax = xmax, box.ymax = ymax, box.zmax = zmax;
}

inline
void Smesh::bin_face(E_Int fid)
{
    const auto &pn = Fc[fid];

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

void Smesh::hash_faces()
{
    bin_faces.clear();

    for (E_Int fid = 0; fid < nf; fid++) {
        bin_face(fid);
    }
}
