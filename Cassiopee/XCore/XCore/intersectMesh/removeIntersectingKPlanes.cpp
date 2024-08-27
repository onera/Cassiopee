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
PyObject *handle_slave(const IMesh &M, Karray& sarray, E_Int min_Kmax,
    std::set<E_Int> &faces_to_tri);

static
E_Int get_kmax(const IMesh &M, const Karray &sarray);

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES;
  
    if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVES)) {
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


    // Init and orient marray mesh
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);


    PyObject *slaves_out = PyList_New(0);

    std::set<E_Int> faces_to_tri;

    E_Int min_Kmax = std::numeric_limits<E_Int>::max();

    for (E_Int i = 0; i < nslaves; i++) {
        E_Int KMAX = get_kmax(M, sarrays[i]);
        if (KMAX < min_Kmax) min_Kmax = KMAX;
    }

    for (E_Int i = 0; i < nslaves; i++) {
        printf("Projecting %d / %d\n", i+1, nslaves);
        PyObject *st = handle_slave(M, sarrays[i], min_Kmax, faces_to_tri);
        PyList_Append(slaves_out, st);
        Py_DECREF(st);
        Karray_free_structured(sarrays[i]);
    }

    puts("Triangulating master projection faces...");

    E_Int nf = M.nf;

    E_Int face_incr = faces_to_tri.size();

    M.F.resize(nf + face_incr);

    std::vector<E_Int> owner(nf, -1), neigh(nf, -1);

    for (E_Int i = 0; i < M.nc; i++) {
        const auto &pf = M.C[i];
        for (E_Int fid : pf) {
            if (owner[fid] == -1) owner[fid] = i;
            else neigh[fid] = i;
        }
    }

    for (E_Int fid : faces_to_tri) {

        auto &pn = M.F[fid];

        assert(pn.size() == 4);

        std::vector<E_Int> tri0(3), tri1(3);

        tri0[0] = pn[0], tri0[1] = pn[1], tri0[2] = pn[2];

        tri1[0] = pn[0], tri1[1] = pn[2], tri1[2] = pn[3];

        pn = tri0;

        M.F[nf] = tri1;

        E_Int own = owner[fid];

        assert(own != -1);

        auto &pf = M.C[own];

        pf.push_back(nf);

        E_Int nei = neigh[fid];

        if (nei != -1) {
            auto &pf = M.C[nei];

            pf.push_back(nf);
        }

        nf++;
    }

    assert(nf == M.nf + face_incr);

    M.nf = nf;

    puts("Exporting projection meshes...");

    PyObject *master_out = M.export_karray();

    PyObject *out = PyList_New(0);
    PyList_Append(out, master_out);
    PyList_Append(out, slaves_out);

    Py_DECREF(master_out);
    Py_DECREF(slaves_out);

    Karray_free_ngon(marray);

    return out;
}

static
E_Int get_kmax(const IMesh &M, const Karray &sarray)
{
    const E_Float *Xs = sarray.X;
    const E_Float *Ys = sarray.Y;
    const E_Float *Zs = sarray.Z;

    const E_Int ni = sarray.ni;
    const E_Int nj = sarray.nj;
    const E_Int nk = sarray.nk;
    const E_Int nij = ni * nj;
    

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

    printf("Intersection plane index k: %d / %d\n", kmax, nk);

    return kmax;
}

static
PyObject *handle_slave(const IMesh &M, Karray& sarray, E_Int kmax,
    std::set<E_Int> &faces_to_tri)
{
    E_Float *Xs = sarray.X;
    E_Float *Ys = sarray.Y;
    E_Float *Zs = sarray.Z;

    E_Int ni = sarray.ni;
    E_Int nj = sarray.nj;
    E_Int nk = sarray.nk;
    E_Int nij = ni * nj;
    

    /**************************************************************************/

    // points to be projected nij*(kmax-1) .. nij*kmax
    std::vector<E_Int> proj_points;
    E_Int ind = (kmax-1)*nij;

    //std::vector<E_Int> upper_points;

    for (E_Int l = 0; l < nij; l++) {
        proj_points.push_back(ind);
        //upper_points.push_back(ind+nij);
        ind++;
    }

    //char fname[128] = {};
    //sprintf(fname, "proj_points%d", idx);
    //point_write("proj_points", Xs, Ys, Zs, proj_points);
    //point_write("upper_points", Xs, Ys, Zs, upper_points);

    /**************************************************************************/

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

        E_Int hit = M.project_point(px, py, pz, dx, dy, dz, TI);

        assert(hit);

        //printf("hit!\n");

        point_hit_table[p] = TI;
    }

    for (const auto &ploc : point_hit_table) {
        E_Int fid = ploc.second.face;
        const auto &pn = M.F[fid];
        if (pn.size() != 3) faces_to_tri.insert(fid);
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
