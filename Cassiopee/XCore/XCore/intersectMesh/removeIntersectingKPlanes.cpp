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

#include <chrono>

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

    //using std::chrono::high_resolution_clock;
    //using std::chrono::duration_cast;
    //using std::chrono::milliseconds;
    //auto t1 = high_resolution_clock::now();

    for (E_Int j = 0; j < nj; j++) {
        for (E_Int i = 0; i < ni; i++) {

            E_Int base = i + ni*j;
            //bool was_inside = true;

            for (E_Int k = nk-1; k >= 0; k--) {

                E_Int idx = base + nij*k;

                E_Float px = Xs[idx];
                E_Float py = Ys[idx];
                E_Float pz = Zs[idx];

                if (!Mf.is_point_inside(px, py, pz)) {
                //if (M->is_point_inside(px, py, pz)) {
                    kmax[base] = k-1;

                    // Cache the point to be projected
                    E_Int proj_id = base + nij*kmax[base];
                    proj_points.push_back(proj_id);

                    //was_outside = false;
                    break;
                }
            }
            
            /*
            if (was_outside) {
                // i-j line completely outside of M
                // Projection points is the last point
                kmax[base] = nk-2;
                E_Int proj_id = base + nij*kmax[base];
                proj_points.push_back(proj_id);
            }
            */
        }
    }

    //auto t2 = high_resolution_clock::now();
    //auto ms_int = duration_cast<milliseconds>(t2-t1);
    //std::cout << ms_int.count() << "ms\n";

    assert(proj_points.size() == (size_t)nij);

    //point_write("proj_points.im", Xs, Ys, Zs, proj_points);
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

        std::vector<PointLoc> mlocs;
        Mf.ray_intersect_BVH(px, py, pz, dx, dy, dz, Mf.root_node_idx, mlocs);

        //if (mlocs.empty()) {
        //    point_write("lost.im", px, py, pz);
        //}

        assert(mlocs.size() > 0);
        PointLoc ploc;
        E_Float min_abs_t = EFLOATMAX;
        for (const auto &mloc : mlocs) {
            if (fabs(mloc.t) < min_abs_t) {
                min_abs_t = fabs(mloc.t);
                ploc = mloc;
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

E_Int get_kmax(IMesh *M, Karray& sarray)
{
    E_Int ni = sarray.ni;
    E_Int nj = sarray.nj;
    E_Int nk = sarray.nk;
    E_Int nij = ni * nj;

    E_Float *Xs = sarray.x;
    E_Float *Ys = sarray.y;
    E_Float *Zs = sarray.z;

    E_Int kmax = -1;

    for (E_Int k = 0; k < nk && kmax == -1; k++) {
        for (E_Int j = 0; j < nj && kmax == -1; j++) {
            for (E_Int i = 0; i < ni; i++) {
                E_Int idx = i + ni*j + nij*k;

                E_Float px = Xs[idx];
                E_Float py = Ys[idx];
                E_Float pz = Zs[idx];

                if (M->is_point_inside(px, py, pz)) {
                    kmax = k-1;
                    return kmax;
                }
            }
        }
    }

    assert(0);
    return kmax;
}

PyObject *handle_slave2(IMesh *M, Karray& sarray, E_Int kmax)
{
    E_Int ni = sarray.ni;
    E_Int nj = sarray.nj;
    E_Int nk = sarray.nk;
    E_Int nij = ni * nj;

    E_Float *Xs = sarray.x;
    E_Float *Ys = sarray.y;
    E_Float *Zs = sarray.z;

    // Indices of points to be projected
    std::vector<E_Int> proj_points;

    for (E_Int j = 0; j < nj; j++) {
        for (E_Int i = 0; i < ni; i++) {

            E_Int idx = i + ni*j + nij*kmax;

            proj_points.push_back(idx);
        }
    }

    //point_write("proj_points", Xs, Ys, Zs, proj_points);

    //E_Int np = ni*nj*nk;

    // Project points onto marray surface
    std::unordered_map<E_Int, TriangleIntersection> point_hit_table;

    for (E_Int i = 0; i < nij; i++) {
        E_Int p = proj_points[i];
        E_Int q = p + nij;

        E_Float px = Xs[p];
        E_Float py = Ys[p];
        E_Float pz = Zs[p];

        E_Float dx = Xs[q] - Xs[p];
        E_Float dy = Ys[q] - Ys[p];
        E_Float dz = Zs[q] - Zs[p];

        E_Float NORM = sqrt(dx*dx + dy*dy + dz*dz);

        dx /= NORM, dy /= NORM, dz /= NORM;

        TriangleIntersection TI;

        E_Int hit = M->project_point(px, py, pz, dx, dy, dz, TI, i);

        assert(hit);

        TI.pid = p;

        point_hit_table[p] = TI;
    }

    // Make out cartesian mesh
    PyObject *tpl;
    nk = kmax + 2;
    tpl = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, 3);

    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *c;
    char *varString, *eltType;
    K_ARRAY::getFromArray3(tpl, varString, f, ni, nj, nk, c, eltType);

    E_Float *xt = f->begin(1);
    E_Float *yt = f->begin(2);
    E_Float *zt = f->begin(3);

    // Copy all the points up to kmax
    for (E_Int k = 0; k < kmax+1; k++) {
        for (E_Int j = 0; j < nj; j++) {
            for (E_Int i = 0; i < ni; i++) {
                E_Int ind = i + j*ni + k*nij;
                xt[ind] = Xs[ind];
                yt[ind] = Ys[ind];
                zt[ind] = Zs[ind];
            }
        }
    }

    // Copy the projected points
    E_Int ind = nij*(kmax+1);
    for (E_Int i = 0; i < nij; i++) {
        E_Int p = proj_points[i];
        auto EH = point_hit_table[p];
        E_Float x = EH.x;
        E_Float y = EH.y;
        E_Float z = EH.z;
        xt[ind] = x;
        yt[ind] = y;
        zt[ind] = z;
        ind++;
    }

    // Tag the projected points
    npy_intp dims[2];
    dims[1] = 1;
    dims[0] = (npy_intp)ni*nj*nk;
    PyArrayObject *tag = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    E_Float *ptag = (E_Float *)PyArray_DATA(tag);
    for (E_Int i = 0; i < nij*(kmax+1); i++) ptag[i] = 0.0;
    for (E_Int i = nij*(kmax+1); i < nij*nk; i++) ptag[i] = 1.0;

    PyObject *out = PyList_New(0);
    PyList_Append(out, tpl);
    PyList_Append(out, (PyObject *)tag);
    RELEASESHAREDS(tpl, f);
    Py_DECREF(tpl);
    Py_DECREF(tag);

    return out;
}

#include "smesh.h"
//#include <random>
//#include "precise.h"

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

    /*
    E_Float a = 2.654351;
    E_Float b = 8.321359;
    E_Float c = 11.63157;
    E_Float d = 29.68484;
    auto E = difference_of_products(a, b, c, d);
    E.print();
    printf("%d\n", E.sign());
    exit(0);
    */

    IMesh *M = (IMesh *)PyCapsule_GetPointer(MASTER, "IntersectMesh");

    /*
    E_Float lo = -1.0;
    E_Float up = 1.0;
    std::uniform_real_distribution<E_Float> unif(lo, up);
    std::default_random_engine re;

    E_Float px = 0, py = 0, pz = 0;
    for (int i = 0; i < 1; i++) {
        //E_Float dx = unif(re);
        //E_Float dy = unif(re);
        //E_Float dz = unif(re);
        //E_Float n = sqrt(dx*dx + dy*dy + dz*dz);
        //dx /= n; dy /= n; dz /= n;
        E_Float dx = 1.0, dy = 0.0, dz = 0.0;
        Ray ray(px, py, pz, dx, dy, dz);
        bool inside = Mf.is_point_inside(ray);
        assert(inside);
    }

    return Py_None;
    */

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

    /*
    E_Int kmax = 10000000;

    for (E_Int i = 0; i < nslaves; i++) {
        E_Int k = get_kmax(M, sarrays[i]);
        printf("k: %d\n", k);
        if (k < kmax) kmax = k;
    }

    printf("kmax: %d\n", kmax);
    */

    PyObject *slaves_out = PyList_New(0);

    for (E_Int i = 0; i < nslaves; i++) {
        printf("Projecting %d / %d\n", i+1, nslaves);
        //PyObject *st = handle_slave2(M, sarrays[i], kmax);
        PyObject *st = handle_slave(M, M->Mf, sarrays[i]);
        PyList_Append(slaves_out, st);
        Py_DECREF(st);
        Karray_free_structured(sarrays[i]);
    }

    return slaves_out;
}
