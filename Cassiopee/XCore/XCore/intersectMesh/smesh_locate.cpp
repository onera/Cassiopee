#include "smesh.h"
#include "triangle.h"
#include "primitives.h"
#include "io.h"

void Smesh::replace_by_projections(const std::vector<E_Int> &pids,
    const std::vector<PointLoc> &plocs)
{
    for (size_t i = 0; i < pids.size(); i++) {
        E_Int pid = pids[i];
        const auto &ploc = plocs[pid];
        X[pid] = ploc.x;
        Y[pid] = ploc.y;
        Z[pid] = ploc.z;
    }
}

void Smesh::correct_near_points_and_edges(Smesh &Sf,
    std::vector<PointLoc> &plocs)
{
    E_Int on_vertex = 0, on_edge = 0;
    for (size_t i = 0; i < plocs.size(); i++) {
        auto &ploc = plocs[i];

        E_Int fid = ploc.fid;
        assert(fid < nf);
        const auto &pn = Fc[fid];

        if (ploc.v_idx != -1) {
            on_vertex++;
            E_Int p = pn[ploc.v_idx];
            E_Float dx = X[p]-Sf.X[i];
            E_Float dy = Y[p]-Sf.Y[i];
            E_Float dz = Z[p]-Sf.Z[i];
            E_Float dist = dx*dx + dy*dy + dz*dz;
            if (dist >= Sf.min_pdist_squared) {
                fprintf(stderr, "Tight near-vertex situation!\n");
                point_write("mpoint", X[p], Y[p], Z[p]);
                point_write("spoint", Sf.X[i], Sf.Y[i], Sf.Z[i]);
                assert(0);
            } else {
                Sf.X[i] = X[p];
                Sf.Y[i] = Y[p];
                Sf.Z[i] = Z[p];
            }
        } else if (ploc.e_idx != -1) {
            on_edge++;
            E_Float u = ploc.bcrd[0];
            E_Float v = ploc.bcrd[1];
            E_Float w = ploc.bcrd[2];
            assert(Sign(w, NEAR_EDGE_TOL) == 0);
            u += w;
            assert(Sign(u+v-1) == 0);
            E_Int p = pn[ploc.e_idx];
            E_Int q = pn[(ploc.e_idx+1)%pn.size()];
            Sf.X[i] = u*X[p] + v*X[q];
            Sf.Y[i] = u*Y[p] + v*Y[q];
            Sf.Z[i] = u*Z[p] + v*Z[q];
        }
    }
    printf("on vertex: %d - on edge: %d\n", on_vertex, on_edge);
}

void Smesh::make_bbox()
{
    NX = 100;
    NY = 100;
    NZ = 100;
    NXY = NX * NY;
    NXYZ = NXY * NZ;

    xmin = ymin = zmin = std::numeric_limits<E_Float>::max();
    xmax = ymax = zmax = std::numeric_limits<E_Float>::min();

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

    xmin = xmin - dx*0.01;
    ymin = ymin - dy*0.01;
    zmin = zmin - dz*0.01;
    xmax = xmax + dx*0.01;
    ymax = ymax + dy*0.01;
    zmax = zmax + dz*0.01;

    HX = (xmax - xmin) / NX;
    HY = (ymax - ymin) / NY;
    HZ = (zmax - zmin) / NZ;
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

std::vector<PointLoc> Smesh::locate(const Smesh &Sf) const
{
    std::vector<PointLoc> ploc(Sf.np);

    for (E_Int pid = 0; pid < Sf.np; pid++) {

        E_Float x = Sf.X[pid];
        E_Float y = Sf.Y[pid];
        E_Float z = Sf.Z[pid];

        E_Int I = floor((x - xmin) / HX);
        E_Int J = floor((y - ymin) / HY);
        E_Int K = floor((z - zmin) / HZ);
        E_Int voxel = get_voxel(I, J, K);

        const auto &pf = bin_faces.at(voxel);

        bool found = false;

        auto &loc = ploc[pid];

        /*
        if (pid == 371)
            point_write("lost", x, y, z);
        */

        for (size_t i = 0; i < pf.size() && !found; i++) {
            E_Int fid = pf[i];
            const auto &pn = Fc[fid];
            const E_Float *fc = &fcenters[3*fid];

            /*
            if (pid == 371) {
                write_face("fid", fid);
                point_write("fc", fc[0], fc[1], fc[2]);
            }
            */

            for (size_t j = 0; j < pn.size(); j++) {
                E_Int p = pn[j];
                E_Int q = pn[(j+1)%pn.size()];

                /*
                if (pid == 371) {
                    point_write("p", X[p], Y[p], Z[p]);
                    point_write("q", X[q], Y[q], Z[q]);
                }
                */

                E_Float u, v, w;

                found = Triangle::is_point_inside(x, y, z,
                    X[p], Y[p], Z[p],
                    X[q], Y[q], Z[q],
                    fc[0], fc[1], fc[2],
                    u, v, w
                );

                if (found) {
                    loc.fid = fid;
                    loc.sub = j;
                    loc.bcrd[0] = u;
                    loc.bcrd[1] = v;
                    loc.bcrd[2] = w;

                    // on p
                    if      (Sign(1-u, NEAR_VERTEX_TOL) == 0)
                        loc.v_idx = j;
                    // on q
                    else if (Sign(1-v, NEAR_VERTEX_TOL) == 0)
                        loc.v_idx = (j+1)%pn.size();
                    // on edge {p, q}
                    else if (Sign(w, NEAR_EDGE_TOL) == 0)
                        loc.e_idx = j;

                    break;
                }
            }
        }

        if (!found) {
            point_write("lost", x, y, z);
            write_ngon("bin", pf);
        }

        assert(found);
    }

    return ploc;
}
