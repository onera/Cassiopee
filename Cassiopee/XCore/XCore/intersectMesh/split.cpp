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
#include "common/Karray.h"
#include <stack>
#include <map>

PyObject *K_XCORE::split_connex(PyObject *self, PyObject *args)
{
    PyObject *MESH, *CTAG, *PTAG;
    if (!PYPARSETUPLE_(args, OOO_, &MESH, &CTAG, &PTAG)) 
    {
        RAISE("Bad input.");
        return NULL;
    }

    Karray array;
    if (Karray_parse_ngon(MESH, array) != 0) 
    {
        RAISE("Failed to parse NGon.");
        return NULL;
    }
    E_Int nc = array.ncells();
    E_Int np = array.npoints();

    E_Float *ctag = NULL;
    E_Int size = -1;
    E_Int ret = K_NUMPY::getFromNumpyArray(CTAG, ctag, size);
    if (ret != 1 || size != nc) {
        RAISE("Bad ctag array.");
        return NULL;
    }

    E_Float *ptag = NULL;
    ret = K_NUMPY::getFromNumpyArray(PTAG, ptag, size);
    if (ret != 1 || size != np) 
    {
        RAISE("Bad ptag array.");
        return NULL;
    }

    // Cell-cell connectivity
    E_Int nf = array.nfaces();
    std::vector<E_Int> owner(nf, -1), neigh(nf, -1);
    for (E_Int cid = 0; cid < nc; cid++) {
        E_Int stride = -1;
        E_Int *pf = array.get_cell(cid, stride);
        for (E_Int i = 0; i < stride; i++) {
            E_Int fid = pf[i]-1;
            if (owner[fid] == -1) owner[fid] = cid;
            else neigh[fid] = cid;
        }
    }

    std::vector<std::vector<E_Int>> clists;

    E_Int seed = 0;
    std::vector<E_Int> visited(nc, 0);

    while (1) {
        while (seed < nc && visited[seed]) seed++;
        if (seed == nc) break;
    
        std::stack<E_Int> stk;
        clists.push_back(std::vector<E_Int>());
        auto &cell_list = clists.back();
        assert(cell_list.empty());
        stk.push(seed);

        while (!stk.empty()) {
            E_Int cid = stk.top();
            stk.pop();

            if (visited[cid]) continue;

            visited[cid] = 1;

            cell_list.push_back(cid);

            E_Int stride = -1;
            E_Int *pf = array.get_cell(cid, stride);
            for (E_Int i = 0; i < stride; i++) {
                E_Int fid = pf[i]-1;
                E_Int nei = owner[fid] == cid ? neigh[fid] : owner[fid];
                if (nei == -1) continue;
                if (!visited[nei]) stk.push(nei);
            }
        }
    }

    assert(seed == nc);
    printf("Connex parts: %zu\n", clists.size());
    {
        E_Int size = 0;
        for (const auto &clist : clists) size += clist.size();
        assert(size == nc);
    }

    // Construct and export
    PyObject *out = PyList_New(0);

    E_Float *X = array.x;
    E_Float *Y = array.y;
    E_Float *Z = array.z;

    PyObject *alist = PyList_New(0);
    PyObject *ctlist = PyList_New(0);
    PyObject *ptlist = PyList_New(0);

    for (const auto &clist : clists) {
        // Build array from clist
        E_Int sizeNGon = 0, sizeNFace = 0;

        std::map<E_Int, E_Int> fmap;
        std::map<E_Int, E_Int> pmap;
        nf = 0;
        E_Int np = 0;
        nc = clist.size();

        for (E_Int cid : clist) {
            E_Int stride = -1;
            E_Int *pf = array.get_cell(cid, stride);
            sizeNFace += stride;

            for (E_Int i = 0; i < stride; i++) {
                E_Int fid = pf[i];

                if (fmap.find(fid) == fmap.end()) {
                    fmap[fid] = nf;
                    E_Int NP = -1;
                    E_Int *pn = array.get_face(fid-1, NP);
                    sizeNGon += NP;

                    for (E_Int j = 0; j < NP; j++) {
                        E_Int pid = pn[j];
                        if (pmap.find(pid) == pmap.end()) {
                            pmap[pid] = np;
                            np++;
                        }
                    }

                    nf++;
                }
            }
        }

        const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

        PyObject *arr = K_ARRAY::buildArray3(3, varString, np, nc, nf, "NGON",
            sizeNGon, sizeNFace, 3, false, 3);
        
        printf("np: %d\n", np);
        printf("nf: %d\n", nf);
        printf("nc: %d\n", nc);
    
        K_FLD::FldArrayF *f;
        K_FLD::FldArrayI *cn;
        K_ARRAY::getFromArray3(arr, f, cn);

        // Points
        E_Float *px = f->begin(1);
        E_Float *py = f->begin(2);
        E_Float *pz = f->begin(3);

        for (const auto &pids : pmap) {
            E_Int opid = pids.first-1;
            E_Int npid = pids.second;
            px[npid] = X[opid];
            py[npid] = Y[opid];
            pz[npid] = Z[opid];
        }

        // Faces
        E_Int *indPG = cn->getIndPG();
        E_Int *ngon = cn->getNGon();
        indPG[0] = 0;

        for (const auto &fids : fmap) 
        {
            E_Int ofid = fids.first-1;
            E_Int nfid = fids.second;
            assert(nfid < nf);
            E_Int NP = -1;
            array.get_face(ofid, NP);
            indPG[nfid+1] = NP;
        }
        for (E_Int fid = 0; fid < nf; fid++) indPG[fid+1] += indPG[fid];
        assert(indPG[nf] == sizeNGon);

        for (const auto &fids : fmap) 
        {
            E_Int ofid = fids.first-1;
            E_Int nfid = fids.second;
            E_Int NP = -1;
            E_Int* pn = array.get_face(ofid, NP);
            E_Int* ptr = &ngon[indPG[nfid]];
            for (E_Int i = 0; i < NP; i++) ptr[i] = pmap.at(pn[i])+1;
        }

        // Cells
        E_Int *indPH = cn->getIndPH();
        E_Int *nface = cn->getNFace();
        indPH[0] = 0;
        nc = 0;
        for (E_Int cid : clist) 
        {
            E_Int NF = -1;
            //E_Int *pf = 
            array.get_cell(cid, NF);
            indPH[nc+1] = NF;
            nc++;
        }
        for (E_Int cid = 0; cid < nc; cid++) indPH[cid+1] += indPH[cid];
        assert(indPH[nc] == sizeNFace);

        nc = 0;
        for (E_Int cid : clist) 
        {
            E_Int NF = -1;
            E_Int *pf = array.get_cell(cid, NF);
            E_Int *ptr = &nface[indPH[nc]];
            for (E_Int i = 0; i < NF; i++) ptr[i] = fmap.at(pf[i])+1 ;
            nc++;
        }

        npy_intp dims[2];
        dims[1] = 1;

        dims[0] = (npy_intp)nc;
        PyArrayObject *ct = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        E_Float *ptr = (E_Float *)PyArray_DATA(ct);
        for (E_Int i = 0; i < nc; i++) ptr[i] = ctag[clist[i]];

        dims[0] = (npy_intp)np;
        PyArrayObject *pt = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        ptr = (E_Float *)PyArray_DATA(pt);
        for (const auto &pids : pmap) {
            E_Int npid = pids.second;
            E_Int opid = pids.first-1;
            ptr[npid] = ptag[opid];
        }


        PyList_Append(alist, arr);
        Py_DECREF(arr);
        PyList_Append(ctlist, (PyObject *)ct);
        Py_DECREF(ct);
        PyList_Append(ptlist, (PyObject *)pt);
        Py_DECREF(pt);
    }

    PyList_Append(out, alist);
    PyList_Append(out, ctlist);
    PyList_Append(out, ptlist);
    Py_DECREF(alist);
    Py_DECREF(ctlist);
    Py_DECREF(ptlist);

    return out;
}