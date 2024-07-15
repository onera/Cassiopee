#include <unordered_map>

#include "xcore.h"
#include "../common/common.h"
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

    E_Int ret;

    ret = Karray_parse_ngon(MASTER, marray);

    if (ret != 0) return NULL;

    ret = Karray_parse_structured(SLAVE, sarray);

    if (ret != 0) {
        Karray_free_ngon(marray);
        return NULL;
    }

    // Check intersection patch
    E_Int *patch = NULL;
    E_Int patch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(PATCH, patch, patch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_structured(sarray);
        RAISE("Bad master patch.");
        return NULL;
    }

    // Zero-based
    for (E_Int i = 0; i < patch_size; i++) patch[i] -= 1;

    // Init and orient marray mesh
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);

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

    for (E_Int k = 0; k < nk; k++) {
        E_Int inside = 0;

        for (E_Int j = 0; j < nj && !inside; j++) {
            for (E_Int i = 0; i < ni; i++) {
                E_Int p = i + ni*j + nij*k;

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
    std::vector<E_Int> proj_points;
    E_Int kshift = (kmax-1)*nij;

    for (E_Int j = 0; j < nj; j++) {
        for (E_Int i = 0; i < ni; i++) {
            E_Int ind = i + j*ni + kshift;
            proj_points.push_back(ind);
        }
    }

    //point_write("proj_points", Xs, Ys, Zs, proj_points);

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

        for (E_Int fid = 0; fid < patch_size; fid++) {
            E_Int face = patch[fid];
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
    for (E_Int i = 0; i < nij; i++) {
        E_Int p = proj_points[i];
        auto EH = point_hit_table[p];
        E_Float x = EH.x;
        E_Float y = EH.y;
        E_Float z = EH.z;
        xt[p+nij] = x;
        yt[p+nij] = y;
        zt[p+nij] = z;
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
    Karray_free_ngon(marray);
    Karray_free_structured(sarray);

    return out;
}