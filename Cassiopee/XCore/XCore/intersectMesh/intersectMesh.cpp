/*    
    Copyright 2013-2024 Onera.

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
#include "xcore.h"
#include "karray.h"
#include "common/common.h"
#include "mesh.h"
#include "smesh.h"
#include "dcel.h"
#include "vertex.h"
#include "face.h"
#include "hedge.h"
#include "io.h"
#include "cycle.h"
#include "triangle.h"
#include "primitives.h"

static
Int EdgeEdgeIntersect(Float ax, Float ay, Float az, Float bx, Float by,
    Float bz, Float px, Float py, Float pz, Float qx, Float qy, Float qz,
    Float &ix, Float &iy, Float &iz)
{
    // P1 = a
    // P2 = b
    // Q1 = p
    // Q2 = q


    Float d1[3] = {bx-ax, by-ay, bz-az};
    Float d2[3] = {qx-px, qy-py, qz-pz};
    Float r[3] = {px-ax, py-ay, pz-az};

    Float d1d2[3];
    K_MATH::cross(d1, d2, d1d2);
    Float denom = K_MATH::dot(d1d2, d1d2, 3);

    if (Sign(denom) == 0) {

        Float colli[3];
        K_MATH::cross(d1, r, colli);
        Float NORM = K_MATH::norm(colli, 3);
        if (Sign(NORM) == 0) {
            assert("collinear!" && 0);
        } else {
            return 0;
        }
    }

    Float tmp[3];
    K_MATH::cross(r, d2, tmp);
    Float t = K_MATH::dot(tmp, d1d2, 3);
    t /= denom;
    if (t < -TOL) return 0;

    K_MATH::cross(r, d1, tmp);
    Float u = K_MATH::dot(tmp, d1d2, 3);
    u /= denom;
    if (u < -TOL || u > 1 + TOL) return 0;

    ix = px + u*(qx - px);
    iy = py + u*(qy - py);
    iz = pz + u*(qz - pz);

    assert(Sign(u) != 0);
    assert(Sign(1-u) != 0);

    return 1;
}

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE, *MPATCH, *SPATCH;
  
    if (!PYPARSETUPLE_(args, OOOO_, &MASTER, &SLAVE, &MPATCH, &SPATCH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    Karray sarray;

    Int ret;

    ret = Karray_parse_ngon(MASTER, marray);

    if (ret != 0) return NULL;

    ret = Karray_parse_ngon(SLAVE, sarray);

    if (ret != 0) {
        Karray_free_ngon(marray);
        return NULL;
    }

    // Init and orient master/slave meshes
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);
    IMesh S(*sarray.cn, sarray.X, sarray.Y, sarray.Z, sarray.npts);

    // Check master intersection patch (zero-based)
    Int *mpatch = NULL;
    Int mpatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(MPATCH, mpatch, mpatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad master patch.");
        return NULL;
    }

    printf("Master patch: " SF_D_ " faces\n", mpatch_size);

    // Check slave intersection patch (zero-based)
    Int *spatch = NULL;
    Int spatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(SPATCH, spatch, spatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad slave patch.");
        return NULL;
    }

    printf("Slave patch: " SF_D_ " faces\n", spatch_size);

    for (Int i = 0; i < mpatch_size; i++) M.patch.insert(mpatch[i]);
    for (Int i = 0; i < spatch_size; i++) S.patch.insert(spatch[i]);

    M.orient_skin(OUT);
    S.orient_skin(IN);

    // Extract surface meshes
    Smesh Mf(M);
    Smesh Sf(S);

    Mf.write_ngon("Mf");
    Sf.write_ngon("Sf");

    Mf.make_point_edges();
    Sf.make_point_edges();

    Dcel D(Mf, Sf);

    D.locate_spoints(Mf, Sf);

    return Py_None;

    std::map<Int, std::vector<Int>> spoints_to_mfaces;

    for (Int sp = 0; sp < Sf.np; sp++) {
        
        for (Int mf = 0; mf < Mf.nf; mf++) {

            const auto &cn = Mf.F[mf];

            Float o[3] = {0, 0, 0};
            for (Int p : cn) {
                o[0] += Mf.X[p];
                o[1] += Mf.Y[p];
                o[2] += Mf.Z[p];
            }
            for (Int i = 0; i < 3; i++) o[i] /= cn.size();

            for (size_t i = 0; i < cn.size(); i++) {
                Int p = cn[i];
                Int q = cn[(i+1)%cn.size()];

                Float u, v, w;

                if (Triangle::is_point_inside(Sf.X[sp], Sf.Y[sp], Sf.Z[sp],
                    Mf.X[p], Mf.Y[p], Mf.Z[p],
                    Mf.X[q], Mf.Y[q], Mf.Z[q],
                    o[0], o[1], o[2],
                    u, v, w)) {
                    spoints_to_mfaces[sp].push_back(mf);
                    break;
                }
            }
        }

        assert(spoints_to_mfaces.find(sp) != spoints_to_mfaces.end());
    }

    std::map<Int, std::vector<Int>> mfaces_to_spoints;

    for (const auto &data : spoints_to_mfaces) {
        Int sp = data.first;
        const auto &mfaces = data.second;

        for (Int mf : mfaces) {
            mfaces_to_spoints[mf].push_back(sp);
        }
    }

    auto &mX = Mf.X;
    auto &mY = Mf.Y;
    auto &mZ = Mf.Z;
    auto &sX = Sf.X;
    auto &sY = Sf.Y;
    auto &sZ = Sf.Z;

    std::map<Int, std::vector<Int>> done_edges;

    std::vector<point> xpoints;

    for (const auto &data : mfaces_to_spoints) {
        Int mf = data.first;
        const auto &spts = data.second;

        // Compute face normal


        printf("doing mface %d\n", mf);
        fflush(stdout);

        for (Int sp : spts) {

            const auto &sedges = Sf.P2E[sp];

            for (Int se : sedges) {
                
                // Skip already visited edge
                if (done_edges.find(se) != done_edges.end()) continue;

                printf("    edge %d\n", se);



                // Edge starts from current point
                Int p = Sf.E[se].p;
                Int q = Sf.E[se].q;

                if (p != sp) std::swap(p, q);

                Float dir[3] = {sX[q]-sX[p], sY[q]-sY[p], sZ[q]-sZ[p]}; 
                
                Int current = mf;
                Int found = 0;
                Int max_walks = 10;
                Int walk = 0;

                done_edges[se].push_back(current);

                Float px = sX[p], py = sY[p], pz = sZ[p];

                Int last_hit_edge = -1;
                
                while (!found && walk < max_walks) {
                    printf("        current: %d\n", current);

                    // If q within current, stop
                    
                    auto it = mfaces_to_spoints.find(current);

                    if (it != mfaces_to_spoints.end()) {
                        const auto &pts = it->second;

                        for (Int pt : pts) {
                            if (pt == q) {
                                found = 1;
                                break;
                            }
                        }
                    }

                    if (found) break;

                    const auto &mpts = Mf.F[current];
                    Float o[3] = {0, 0, 0};
                    for (Int p : mpts) {
                        o[0] += Mf.X[p]; 
                        o[1] += Mf.Y[p]; 
                        o[2] += Mf.Z[p];
                    }
                    for (Int i = 0; i < 3; i++) o[i] /= mpts.size();
                    
                    Float N[3] = {};
                    Int a = mpts[0], b = mpts[1];
                    Float V[3] = {mX[a]-o[0], mY[a]-o[1], mZ[a]-o[2]};
                    Float U[3] = {mX[b]-o[0], mY[b]-o[1], mZ[b]-o[2]};
                    K_MATH::cross(V, U, N);
                    Float NORM = K_MATH::norm(N, 3);
                    for (Int i = 0; i < 3; i++) N[i] /= NORM;

                    // Project pq onto current face's plane

                    Float proj[3] = {};
                    Float dp = K_MATH::dot(dir, N, 3);
                    for (Int i = 0; i < 3; i++)
                        proj[i] = dir[i] - dp * N[i];
                    NORM = K_MATH::norm(proj, 3);
                    for (int i = 0; i < 3; i++) proj[i] /= NORM;

                    // Shoot a ray {p, proj}, get the intersected edge

                    Float qx = px + 2 * proj[0];
                    Float qy = py + 2 * proj[1];
                    Float qz = pz + 2 * proj[2];

                    Int hit_edge_idx = -1;
                    Float ix, iy, iz;

                    const auto &medges = Mf.F2E[current];
                    
                    for (size_t i = 0; i < medges.size(); i++) {
                        Int me = medges[i];
                        printf("            testing edge %d\n", me);

                        if (me == last_hit_edge) {
                            puts("          skipping last hit edge");
                            continue;
                        }

                        ix = 0, iy = 0, iz = 0;


                        Int a = Mf.E[me].p;
                        Int b = Mf.E[me].q;

                        Float ax = mX[a], ay = mY[a], az = mZ[a];
                        Float bx = mX[b], by = mY[b], bz = mZ[b];

                        Int hit = EdgeEdgeIntersect(
                            px, py, pz,
                            qx, qy, qz,
                            ax, ay, az,
                            bx, by, bz,
                            ix, iy, iz);
                        
                        printf("            %f %f %f\n", ix, iy, iz);

                        if (hit) {
                            hit_edge_idx = i;
                            last_hit_edge = me;
                            break;
                        }
                    }

                    assert(hit_edge_idx != -1); 

                    // Traverse edge
                    current = Mf.F2F[current][hit_edge_idx];
                    assert(current != -1);
                    done_edges[se].push_back(current);
                    walk++;
                    px = ix;
                    py = iy;
                    pz = iz;

                    xpoints.push_back(point(px, py, pz));
                }

                assert(found);
            }

        }
    }
    
    point_write("xpoints", xpoints);
   
    
    return Py_None;





    


































    return Py_None;
}
