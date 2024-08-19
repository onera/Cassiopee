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
#include <unordered_map>

#include "xcore.h"
#include "common/common.h"
#include "karray.h"
#include "mesh.h"
#include "ray.h"
#include "io.h"

static
PyObject *handle_slave(const IMesh &M, Karray& sarray, E_Int patch_size, E_Int *patch);

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES, *PATCH;
  
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVES, &PATCH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;

    E_Int ret;

    ret = Karray_parse_ngon(MASTER, marray);

    if (ret != 0) return NULL;

    E_Int nslaves = PyList_Size(SLAVES);
    E_Int i;

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
        Karray_free_ngon(marray);
        for (E_Int j = 0; j < i; j++)
            Karray_free_structured(sarrays[j]);
        return NULL;
    }
    
    // Check intersection patch
    E_Int *patch = NULL;
    E_Int patch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(PATCH, patch, patch_size, false);
    if (ret != 1) {
        Karray_free_ngon(marray);
        for (E_Int i = 0; i < nslaves; i++) Karray_free_structured(sarrays[i]);
        RAISE("Bad master patch.");
        return NULL;
    }

    // Zero-based
    for (E_Int i = 0; i < patch_size; i++) patch[i] -= 1;

    // Init and orient marray mesh
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);

    PyObject *out = PyList_New(0);

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *st = handle_slave(M, sarrays[i], patch_size, patch);
        PyList_Append(out, st);
        Py_DECREF(st);
        Karray_free_structured(sarrays[i]);
    }

    Karray_free_ngon(marray);

    delete [] patch;

    return out;
}

static
PyObject *handle_slave(const IMesh &M, Karray& sarray, E_Int patch_size, E_Int *patch)
{
    E_Float *Xs = sarray.X;
    E_Float *Ys = sarray.Y;
    E_Float *Zs = sarray.Z;

    E_Int ni = sarray.ni;
    E_Int nj = sarray.nj;
    E_Int nk = sarray.nk;
    E_Int nij = ni * nj;
    

    /**************************************************************************/

    // Detect at which k does the slave mesh intersect the marray mesh
    // Get projected points coordinates

    // Max plane index that doesn't intersection with marray bbox (zero-based)
    E_Int kmax = 0;

    std::vector<E_Int> inside_point;

    for (E_Int k = 0; k < nk; k++, kmax++) {
        E_Int inside = 0;

        E_Int p = nij*k;

        for (E_Int l = 0; l < nij; l++, p++) {

            if (M.is_point_inside(Xs[p], Ys[p], Zs[p])) {
                inside_point.push_back(p);
                inside = 1;
                break;
            }
            
        }

        if (inside) break;
    }

    //char fname2[128] = {};
    //sprintf(fname2, "inside%d", idx);
    //point_write(fname2, Xs, Ys, Zs, inside_point);

    printf("Intersection plane index k: %d\n", kmax);

    // points to be projected nij*(kmax-1) .. nij*kmax
    std::vector<E_Int> proj_points;
    E_Int ind = (kmax-1)*nij;

    for (E_Int l = 0; l < nij; l++) {
        proj_points.push_back(ind);
        ind++;
    }

    //char fname[128] = {};
    //sprintf(fname, "proj_points%d", idx);
    //point_write(fname, Xs, Ys, Zs, proj_points);

    /**************************************************************************/

    // Project points onto marray surface
    std::unordered_map<E_Int, TriangleIntersection> point_hit_table;

    for (E_Int i = 0; i < nij; i++) {
        E_Int p = proj_points[i];
        E_Int q = p + nij;

        E_Float px = Xs[p];
        E_Float py = Ys[p];
        E_Float pz = Zs[p];
        E_Float dx = Xs[q]-px;
        E_Float dy = Ys[q]-py;
        E_Float dz = Zs[q]-pz;

        const auto &X = M.X;
        const auto &Y = M.Y;
        const auto &Z = M.Z;

        TriangleIntersection TI;

        E_Int A, B, C;
        E_Int hit = -1;
        E_Float t_last = E_FLOAT_MAX;

        for (E_Int fid = 0; fid < patch_size; fid++) {
            E_Int face = patch[fid];
            const auto &pn = M.F[face];

            // TODO: handle any polygon
            assert(pn.size() == 4);

            A = pn[0]; B = pn[1]; C = pn[2];

            // First triangle ABC

            hit = MollerTrumbore(px, py, pz, dx, dy, dz,
                                 X[A], Y[A], Z[A],
                                 X[B], Y[B], Z[B],
                                 X[C], Y[C], Z[C],
                                 TI);

            if (hit && TI.t < t_last) {
                t_last = TI.t;
                TI.face = face;
                TI.tri = 0;
                point_hit_table[p] = TI;
                continue;
            }
            
            if (pn.size() == 3)
                continue;

            // Second triangle CDA
            A = pn[0]; B = pn[2]; C = pn[3];

            hit = MollerTrumbore(px, py, pz, dx, dy, dz,
                                 X[A], Y[A], Z[A],
                                 X[B], Y[B], Z[B],
                                 X[C], Y[C], Z[C],
                                 TI);

            if (hit && TI.t < t_last) {
                t_last = TI.t;
                TI.face = face;
                TI.tri = 1;
                point_hit_table[p] = TI;
            }
        }

        // point must hit!
        assert(point_hit_table.find(p) != point_hit_table.end());
        if (point_hit_table.find(p) == point_hit_table.end()) {
            puts("AIE!");
            abort();
        }
    }

    /*
    FILE *fh = fopen("hit_faces", "w");
    assert(fh);
    for (const auto hit : point_hit_table) {
        E_Int pt = hit.first;
        const auto &TI = hit.second;
        fprintf(fh, "%d %d\n", pt + nij, TI.face);
    }
    fclose(fh);
    */

    //char fname3[128] = {};
    //sprintf(fname3, "projection%d", idx);
    //edge_write(fname3, Xs, Ys, Zs, point_hit_table);

    /*************************************************************************/
    
    // planes: 0 .... kmax-1
    // plane kmax was removed
    // numbers of planes: kmax + 1

    // Make out cartesian mesh
    PyObject *tpl;
    nk = kmax + 1; 
    tpl = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, 3);

    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *c;
    char *varString, *eltType;
    K_ARRAY::getFromArray3(tpl, varString, f, ni, nj, nk, c, eltType);

    E_Float* xt = f->begin(1);
    E_Float* yt = f->begin(2);
    E_Float* zt = f->begin(3);

    // Copy all the points up to kmax
    for (E_Int k = 0; k < kmax; k++) {
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
    ind = nij*kmax;
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
    for (E_Int i = 0; i < nij*kmax; i++) ptag[i] = 0.0;
    for (E_Int i = nij*kmax; i < nij*nk; i++) ptag[i] = 1.0;

    PyObject *out = PyList_New(0);
    PyList_Append(out, tpl);
    PyList_Append(out, (PyObject *)tag);
    RELEASESHAREDS(tpl, f);
    Py_DECREF(tpl);
    Py_DECREF(tag);
    
    return out;
}
