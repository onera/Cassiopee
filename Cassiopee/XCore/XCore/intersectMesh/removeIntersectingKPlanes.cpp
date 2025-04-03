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
#include "xcore.h"
#include "common/common.h"
#include "common/Karray.h"
#include "mesh.h"
#include "ray.h"
#include "io.h"
#include "smesh.h"

static
void construct_faces(std::vector<std::vector<E_Int>> &faces,
    const std::vector<E_Int> &apoints, const std::vector<E_Int> &dpoints,
    E_Int &NF)
{
    const E_Int *min_ptr = apoints.data();
    const E_Int *max_ptr = dpoints.data();
    E_Int min_range = apoints.size();
    E_Int max_range = dpoints.size();
    
    if (apoints.size() > dpoints.size()) {
        std::swap(min_ptr, max_ptr);
        std::swap(min_range, max_range);
    }
    
    for (E_Int l = 0; l < min_range-1; l++) {
        faces.push_back({min_ptr[l],   max_ptr[l],
                         max_ptr[l+1], min_ptr[l+1]});
        NF++;
    }
    
    std::vector<E_Int> last_face;
    
    last_face.push_back(min_ptr[min_range-1]);
    
    for (E_Int l = min_range-1; l < max_range; l++)
        last_face.push_back(max_ptr[l]);
    
    faces.push_back(last_face);
    NF++;
}

PyObject *handle_slave(const IMesh *M, const Smesh &Mf, Karray& sarray)
{
    E_Int ni = sarray.ni;
    E_Int nj = sarray.nj;
    E_Int nk = sarray.nk;
    E_Int nij = ni * nj;

    E_Float *Xs = sarray.x;
    E_Float *Ys = sarray.y;
    E_Float *Zs = sarray.z;

    // Last k-plane outside of M
    std::vector<E_Int> kmax(nij, -1);

    // Indices of points to be projected
    std::vector<E_Int> proj_points;
    std::vector<E_Int> inside;

    for (E_Int j = 0; j < nj; j++) {
        for (E_Int i = 0; i < ni; i++) {

            E_Int base = i + ni*j;

            for (E_Int k = nk-1; k >= 0; k--) {

                E_Int idx = base + nij*k;

                E_Float px = Xs[idx];
                E_Float py = Ys[idx];
                E_Float pz = Zs[idx];

                if (!Mf.is_point_inside(px, py, pz)) {
                    assert(k >= 1);
                    kmax[base] = k-1;
                    //kmax[base] = k;

                    // Cache the point to be projected
                    E_Int proj_id = base + nij*kmax[base];
                    proj_points.push_back(proj_id);
                    inside.push_back(base + nij*k);

                    break;
                }
            }
        }
    }



    assert(proj_points.size() == (size_t)nij);

    //point_write("proj_points.im", Xs, Ys, Zs, proj_points);
    //point_write("inside.im", Xs, Ys, Zs, inside);
    //printf("points written!\n");
    //fflush(stdout);

    E_Int np = ni*nj*nk;

    // Project points onto marray surface
    std::unordered_map<E_Int, TriangleIntersection> point_hit_table;

    //std::vector<Point> projections;

    for (E_Int i = 0; i < nij; i++) {
        E_Int p = proj_points[i];
        E_Int q = p - nij;

        E_Float px = Xs[p];
        E_Float py = Ys[p];
        E_Float pz = Zs[p];

        E_Float dx = Xs[p] - Xs[q];
        E_Float dy = Ys[p] - Ys[q];
        E_Float dz = Zs[p] - Zs[q];

        E_Float NORM = sqrt(dx*dx + dy*dy + dz*dz);

        dx /= NORM, dy /= NORM, dz /= NORM;

        Ray ray(px, py, pz, dx, dy, dz, Ray::Policy::FORWARD);
        HitData hit_data;
        Mf.intersect_ray(ray, 0, hit_data);

        if (hit_data.locs.empty()) {
            point_write("lost.im", px, py, pz);
            fprintf(stderr, "WARNING: Failed to project point %d!\n", p);
            assert(0);
        }

        PointLoc ploc;
        E_Float min_abs_t = EFLOATMAX;
        for (const auto &loc : hit_data.locs) {
            if (fabs(loc.t) < min_abs_t) {
                min_abs_t = fabs(loc.t);
                ploc = loc;
            }
        }

        TriangleIntersection TI;
        TI.pid = p;
        TI.x = ploc.x, TI.y = ploc.y, TI.z = ploc.z;
        point_hit_table[p] = TI;

        //projections.push_back({ploc.x, ploc.y, ploc.z});
    }

    //point_write("projections.im", projections);

    // Construct the new faces and cells

    std::vector<std::vector<E_Int>> faces;
    std::vector<std::vector<E_Int>> cells;

    E_Int nf, nc;
    nf = nc = 0;

    for (E_Int j = 0; j < nj-1; j++) {
        for (E_Int i = 0; i < ni-1; i++) {

            // Ligne ij

            E_Int bases[4] = {
                i   + ni * j,
                i+1 + ni * j,
                i+1 + ni * (j+1),
                i   + ni * (j+1)
            };
            
            E_Int done = 0;
            
            for (E_Int k = 0; k <= kmax[bases[0]] && !done; k++) {

                E_Int a = bases[0] + k * nij;
                E_Int b = bases[1] + k * nij;
                E_Int c = bases[2] + k * nij;
                E_Int d = bases[3] + k * nij;

                E_Int A, B, C, D;
                A = B = C = D = -1;
                
                // Est-on sur un kmax ?

                if ((k == kmax[bases[0]]) ||
                    (k == kmax[bases[1]]) ||
                    (k == kmax[bases[2]]) ||
                    (k == kmax[bases[3]])) {
                    done = 1;
                }

                
                if (done) {

                    E_Int K;

                    // Accumulate the points
                    std::vector<E_Int> apoints;
                    assert(kmax[bases[0]] >= k);
                    K = k;
                    while (K <= kmax[bases[0]]) {
                        apoints.push_back(bases[0] + K*nij);
                        K++;
                    }
                    apoints.push_back(point_hit_table.at(apoints.back()).pid + np);

                    std::vector<E_Int> bpoints;
                    assert(kmax[bases[1]] >= k);
                    K = k;
                    while (K <= kmax[bases[1]]) {
                        bpoints.push_back(bases[1] + K*nij);
                        K++;
                    }
                    bpoints.push_back(point_hit_table.at(bpoints.back()).pid + np);

                    std::vector<E_Int> cpoints;
                    assert(kmax[bases[2]] >= k);
                    K = k;
                    while (K <= kmax[bases[2]]) {
                        cpoints.push_back(bases[2] + K*nij);
                        K++;
                    }
                    cpoints.push_back(point_hit_table.at(cpoints.back()).pid + np);

                    std::vector<E_Int> dpoints;
                    assert(kmax[bases[3]] >= k);
                    K = k;
                    while (K <= kmax[bases[3]]) {
                        dpoints.push_back(bases[3] + K*nij);
                        K++;
                    }
                    dpoints.push_back(point_hit_table.at(dpoints.back()).pid + np);

                    E_Int NF = 0;

                    // BOT
                    faces.push_back({a, b, c, d});
                    NF++;

                    // TOP
                    faces.push_back({apoints.back(),
                                     bpoints.back(),
                                     cpoints.back(),
                                     dpoints.back()});
                    NF++;

                    // LFT
                    construct_faces(faces, apoints, dpoints, NF);

                    // RGT
                    construct_faces(faces, bpoints, cpoints, NF);

                    // FRO
                    construct_faces(faces, bpoints, apoints, NF);

                    // BCK
                    construct_faces(faces, cpoints, dpoints, NF);


                    // Add the new cell and faces

                    std::vector<E_Int> cell(NF);
                    for (E_Int l = 0; l < NF; l++) cell[l] = nf+l;
                    cells.push_back(cell);

                    nf += NF;
                    nc += 1;
                }
                
                else {

                    A = a + nij;
                    B = b + nij;
                    C = c + nij;
                    D = d + nij;
                
                    faces.push_back({a, b, c, d});
                    faces.push_back({A, B, C, D});
                    faces.push_back({a, d, D, A});
                    faces.push_back({b, c, C, B});
                    faces.push_back({b, a, A, B});
                    faces.push_back({c, d, D, C});

                    cells.push_back({nf, nf+1, nf+2, nf+3, nf+4, nf+5});

                    nf += 6;
                    nc += 1;

                }


            }
        }
    }

    // Renumber points

    std::map<E_Int, E_Int> new_pids;
    std::map<E_Int, E_Int> old_pids;

    E_Int NP = 0;

    for (auto &pn : faces) {
        for (auto &p : pn) {
            auto it = new_pids.find(p);
            if (it  == new_pids.end()) {
                new_pids[p] = NP;
                old_pids[NP] = p;
                p = NP;
                NP++;
            } else {
                p = it->second; 
            }
        }
    }
    
    // Make coordinates

    IMesh new_M;

    new_M.F = faces;
    new_M.C = cells;

    auto &new_X = new_M.X;
    auto &new_Y = new_M.Y;
    auto &new_Z = new_M.Z;

    new_X.resize(NP);
    new_Y.resize(NP);
    new_Z.resize(NP);

    for (const auto &pdat : new_pids) {
        E_Int opid = pdat.first;
        E_Int npid = pdat.second;

        if (opid < np) {
            new_X[npid] = Xs[opid];
            new_Y[npid] = Ys[opid];
            new_Z[npid] = Zs[opid];
        } else {
            const auto &TI = point_hit_table.at(opid - np);
            new_X[npid] = TI.x;
            new_Y[npid] = TI.y;
            new_Z[npid] = TI.z;
        }
    }

    new_M.np = NP;
    new_M.nf = faces.size();
    new_M.nc = cells.size();

    PyObject *tpl = new_M.export_karray();

    // Tag the projected points
    npy_intp dims[2];
    dims[1] = 1;
    dims[0] = NP;
    PyArrayObject *tag = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    E_Float *ptag = (E_Float *)PyArray_DATA(tag);
    memset(ptag, 0, NP*sizeof(E_Float));

    for (const auto &pdat : new_pids) {
        E_Int opid = pdat.first;
        E_Int npid = pdat.second;

        if (opid >= np) {
            ptag[npid] = 1.0;
        }
    }

    PyObject *out = PyList_New(0);
    PyList_Append(out, tpl);
    PyList_Append(out, (PyObject *)tag);
    Py_DECREF(tpl);
    Py_DECREF(tag);
    
    return out;
}

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES;
  
    if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVES)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MASTER, "IntersectMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    IMesh *M = (IMesh *)PyCapsule_GetPointer(MASTER, "IntersectMesh");

    E_Int nslaves = PyList_Size(SLAVES);
    E_Int i, ret;

    std::vector<Karray> sarrays(nslaves);

    for (i = 0; i < nslaves; i++) {
        Karray sarray;
        PyObject *SLAVE = PyList_GetItem(SLAVES, i);
        ret = Karray_parse_structured(SLAVE, sarray);

        if (ret != 0) {
            break;
        }

        sarrays[i] = sarray;
    }

    if (ret != 0) {
        for (E_Int j = 0; j < i; j++)
            Karray_free_structured(sarrays[j]);
        return NULL;
    }

    PyObject *slaves_out = PyList_New(0);

    for (E_Int i = 0; i < nslaves; i++) {
        printf("Projecting %d / %d\n", i+1, nslaves);
        PyObject *st = handle_slave(M, M->Mf, sarrays[i]);
        PyList_Append(slaves_out, st);
        Py_DECREF(st);
        Karray_free_structured(sarrays[i]);
    }

    return slaves_out;
}
