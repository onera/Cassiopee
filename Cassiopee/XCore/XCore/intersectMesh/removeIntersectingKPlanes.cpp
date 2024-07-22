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

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE, *PATCH;
  
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVE, &PATCH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    Karray sarray;

    Int ret;

    ret = Karray_parse_ngon(MASTER, marray);

    if (ret != 0) return NULL;

    ret = Karray_parse_structured(SLAVE, sarray);

    if (ret != 0) {
        Karray_free_ngon(marray);
        return NULL;
    }

    // Check intersection patch
    Int *patch = NULL;
    Int patch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(PATCH, patch, patch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_structured(sarray);
        RAISE("Bad master patch.");
        return NULL;
    }

    // Zero-based
    for (Int i = 0; i < patch_size; i++) patch[i] -= 1;

    // Init and orient marray mesh
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);

    Float *Xs = sarray.X;
    Float *Ys = sarray.Y;
    Float *Zs = sarray.Z;

    Int ni = sarray.ni;
    Int nj = sarray.nj;
    Int nk = sarray.nk;
    Int nij = ni * nj;

    /**************************************************************************/

    // Detect at which k does the slave mesh intersect the marray mesh
    // Get projected points coordinates

    // Max plane index that doesn't intersection with marray bbox (zero-based)
    Int kmax = 0;

    for (Int k = 0; k < nk; k++) {
        Int inside = 0;

        for (Int j = 0; j < nj && !inside; j++) {
            for (Int i = 0; i < ni; i++) {
                Int p = i + ni*j + nij*k;

                if (M.is_point_inside(Xs[p], Ys[p], Zs[p])) {
                    inside = 1;
                    kmax = k;
                    break;
                }
            }
        }

        if (inside) break;
    }

    assert(kmax > 1);

    // points to be projected nij*(kmax-1) .. nij*kmax
    std::vector<Int> proj_points;
    Int kshift = (kmax-1)*nij;

    for (Int j = 0; j < nj; j++) {
        for (Int i = 0; i < ni; i++) {
            Int ind = i + j*ni + kshift;
            proj_points.push_back(ind);
        }
    }

    //point_write("proj_points", Xs, Ys, Zs, proj_points);

    /**************************************************************************/

    // Project points onto marray surface
    std::unordered_map<Int, TriangleIntersection> point_hit_table;

    for (Int i = 0; i < nij; i++) {
        Int p = proj_points[i];
        Int q = p + nij;

        Float px = Xs[p];
        Float py = Ys[p];
        Float pz = Zs[p];
        Float dx = Xs[q]-px;
        Float dy = Ys[q]-py;
        Float dz = Zs[q]-pz;

        const auto &X = M.X;
        const auto &Y = M.Y;
        const auto &Z = M.Z;

        TriangleIntersection TI;

        Int A, B, C;
        Int hit = -1;

        for (Int fid = 0; fid < patch_size; fid++) {
            Int face = patch[fid];
            const auto &pn = M.F[face];

            A = pn[0]; B = pn[1]; C = pn[2];

            // First triangle ABC

            hit = MollerTrumbore(px, py, pz, dx, dy, dz,
                                 X[A], Y[A], Z[A],
                                 X[B], Y[B], Z[B],
                                 X[C], Y[C], Z[C],
                                 TI);

            if (hit) {
                point_hit_table[p] = TI;
                break;
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

            if (hit) {
                point_hit_table[p] = TI;
                break;
            }
        }

        // point must hit!
        assert(hit == 1);
    }

    edge_write("projection", Xs, Ys, Zs, point_hit_table);

    /*************************************************************************/

    // Make out cartesian mesh
    PyObject *tpl;
    nk = kmax + 1; 
    tpl = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, 3);

    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *c;
    char *varString, *eltType;
    K_ARRAY::getFromArray3(tpl, varString, f, ni, nj, nk, c, eltType);

    Float* xt = f->begin(1);
    Float* yt = f->begin(2);
    Float* zt = f->begin(3);

    // Copy all the points up to kmax
    for (Int k = 0; k < kmax; k++) {
        for (Int j = 0; j < nj; j++) {
            for (Int i = 0; i < ni; i++) {
                Int ind = i + j*ni + k*nij;
                xt[ind] = Xs[ind];
                yt[ind] = Ys[ind];
                zt[ind] = Zs[ind];
            }
        }
    }

    // Copy the projected points
    for (Int i = 0; i < nij; i++) {
        Int p = proj_points[i];
        auto EH = point_hit_table[p];
        Float x = EH.x;
        Float y = EH.y;
        Float z = EH.z;
        xt[p+nij] = x;
        yt[p+nij] = y;
        zt[p+nij] = z;
    }

    // Tag the projected points
    npy_intp dims[2];
    dims[1] = 1;
    dims[0] = (npy_intp)ni*nj*nk;
    PyArrayObject *tag = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    Float *ptag = (Float *)PyArray_DATA(tag);
    for (Int i = 0; i < nij*kmax; i++) ptag[i] = 0.0;
    for (Int i = nij*kmax; i < nij*nk; i++) ptag[i] = 1.0;

    PyObject *out = PyList_New(0);
    PyList_Append(out, tpl);
    PyList_Append(out, (PyObject *)tag);
    RELEASESHAREDS(tpl, f);
    Py_DECREF(tpl);
    Py_DECREF(tag);
    Karray_free_ngon(marray);
    Karray_free_structured(sarray);

    return out;
}