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
#include "common/common.h"
#include "karray.h"
#include "mesh.h"
#include "ray.h"
#include "io.h"

static
PyObject *handle_slave(IMesh *M, Karray& sarray, E_Int min_Kmax);

static
E_Int get_kmax(const IMesh *M, const Karray &sarray);

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES;
  
    if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVES)) {
        RAISE("Bad input.");
        return NULL;
    }

    puts("Parsing MASTER and SLAVE");

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

    E_Int min_Kmax = std::numeric_limits<E_Int>::max();

    puts("Getting min_Kmax...");

    for (E_Int i = 0; i < nslaves; i++) {
        E_Int KMAX = get_kmax(M, sarrays[i]);
        if (KMAX < min_Kmax) min_Kmax = KMAX;
    }

    //min_Kmax = 10;

    for (E_Int i = 0; i < nslaves; i++) {
        printf("Projecting %d / %d\n", i+1, nslaves);
        PyObject *st = handle_slave(M, sarrays[i], min_Kmax);
        PyList_Append(slaves_out, st);
        Py_DECREF(st);
        Karray_free_structured(sarrays[i]);
    }

    return slaves_out;
}

static
E_Int get_kmax(const IMesh *M, const Karray &sarray)
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

            if (M->is_point_inside(Xs[p], Ys[p], Zs[p])) {
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
PyObject *handle_slave(IMesh *M, Karray& sarray, E_Int kmax)
{
    printf("Projecting at kmax index: %d\n", kmax);
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

    std::vector<E_Int> upper_points;

    for (E_Int l = 0; l < nij; l++) {
        proj_points.push_back(ind);
        upper_points.push_back(ind+nij);
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

        //point_write("bad_point", px, py, pz);

        E_Float dx = Xs[q] - Xs[p];
        E_Float dy = Ys[q] - Ys[p];
        E_Float dz = Zs[q] - Zs[p];

        E_Float NORM = sqrt(dx*dx + dy*dy + dz*dz);

        dx /= NORM, dy /= NORM, dz /= NORM;

        TriangleIntersection TI;

        E_Int hit = M->project_point(px, py, pz, dx, dy, dz, TI, i);

        /*
        if (i == 8) {
            puts("WRITING BAD EDGE AND POINT");
            point_write("bad_point", px, py, pz);
            E_Float qx = px + 100 * dx;
            E_Float qy = py + 100 * dy;
            E_Float qz = pz + 100 * dz;
            edge_write("bad_edge", px, py, pz, qx, qy, qz);
        }
        */

        //printf("point %d: t = %f\n", i, TI.t);

        assert(hit);

        point_hit_table[p] = TI;
    }

    for (const auto &ploc : point_hit_table) {
        E_Int fid = ploc.second.face;
        const auto &pn = M->F[fid];
        assert(pn.size() == 4);
        M->faces_to_tri.insert(fid);
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
