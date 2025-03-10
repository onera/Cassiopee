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
#include "mesh.h"
#include "common/Karray.h"

PyObject *K_XCORE::triangulate_skin(PyObject *self, PyObject *args)
{
    PyObject *MESH, *PTLISTS;
    if (!PYPARSETUPLE_(args, OO_, &MESH, &PTLISTS)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray array;
    if (Karray_parse_ngon(MESH, array) != 0) {
        RAISE("Mesh should be an NGon.");
        return NULL;
    }

    IMesh M(array);
    M.make_skin();

    std::vector<Py_BC> bcs_in;

    if (PTLISTS != Py_None) {
        E_Int nbcs = PyList_Size(PTLISTS);
        for (E_Int i = 0; i < nbcs; i++) {
            PyObject *PTLIST = PyList_GetItem(PTLISTS, i);
            Py_BC bc;
            E_Int ret = K_NUMPY::getFromNumpyArray(PTLIST, bc.ptr, bc.size, true);
            Py_DECREF(PTLIST);
            if (ret != 1) {
                RAISE("Couldn't extract pointlists.");
                Karray_free_ngon(array);
                return NULL;
            }
            bcs_in.push_back(bc);
        }

        std::unordered_map<E_Int, E_Int> fid_to_bc;

        for (size_t i = 0; i < bcs_in.size(); i++) {
            const auto &bc = bcs_in[i];
            for (E_Int j = 0; j < bc.size; j++) {
                E_Int fid = bc.ptr[j];
                fid_to_bc[fid] = i;
            }
        }

        M.triangulate_skin(fid_to_bc);

        std::vector<std::vector<E_Int>> new_bcs(nbcs);
        for (const auto &it : fid_to_bc) {
            new_bcs[it.second].push_back(it.first);
        }

        PyObject *out = PyList_New(0);
        PyList_Append(out, M.export_karray());
        PyObject *bcs_out = PyList_New(0);

        for (const auto &new_bc : new_bcs) {
            npy_intp dims[2];
            dims[0] = (npy_intp)new_bc.size();
            dims[1] = 1;
            PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
            E_Int *ptr = (E_Int *)PyArray_DATA(arr);
            memcpy(ptr, new_bc.data(), new_bc.size()*sizeof(E_Int));
            PyList_Append(bcs_out, (PyObject *)arr);
            Py_DECREF(arr);
        }

        PyList_Append(out, bcs_out);
        Py_DECREF(bcs_out);
        return out;
    } else {
        M.triangulate_skin();
        return M.export_karray();
    }
}

void IMesh::triangulate_skin(std::unordered_map<E_Int, E_Int> &fid_to_bc)
{
    E_Int to_tri = 0;

    for (E_Int fid : skin) {
        const auto &pn = F[fid];
        if (pn.size() == 3) continue;
        assert(pn.size() == 4);
        to_tri++;
    }

    F.resize(nf + 3*to_tri);
    X.resize(np + to_tri, 0);
    Y.resize(np + to_tri, 0);
    Z.resize(np + to_tri, 0);

    for (E_Int fid : skin) {
        const auto &pn = F[fid];
        if (pn.size() == 3) continue;
        for (E_Int pid : pn) {
            X[np] += X[pid];
            Y[np] += Y[pid];
            Z[np] += Z[pid];
        }
        X[np] *= 0.25;
        Y[np] *= 0.25;
        Z[np] *= 0.25;

        E_Int nodes[4] = {pn[0], pn[1], pn[2], pn[3]};

        F[fid] = {nodes[0], nodes[1], np};
        F[nf] = {nodes[1], nodes[2], np};
        F[nf+1] = {nodes[2], nodes[3], np};
        F[nf+2] = {nodes[3], nodes[0], np};

        E_Int own = owner[fid];
        auto &pf = C[own];
        size_t stride = pf.size();
        pf.resize(stride + 3);
        for (int i = 0; i < 3; i++) pf[stride+i] = nf+i;

        // One-based
        auto it = fid_to_bc.find(fid+1);
        if (it != fid_to_bc.end()) {
            for (int i = 0; i < 3; i++) fid_to_bc[nf+i+1] = it->second;
        }

        np += 1;
        nf += 3;
    }

    assert((size_t)nf == F.size());
    assert((size_t)np == X.size());
}

void IMesh::triangulate_skin()
{
    E_Int to_tri = 0;

    for (E_Int fid : skin) {
        const auto &pn = F[fid];
        if (pn.size() == 3) continue;
        assert(pn.size() == 4);
        to_tri++;
    }

    F.resize(nf + 3*to_tri);
    X.resize(np + to_tri, 0);
    Y.resize(np + to_tri, 0);
    Z.resize(np + to_tri, 0);

    for (E_Int fid : skin) {
        const auto &pn = F[fid];
        if (pn.size() == 3) continue;
        for (E_Int pid : pn) {
            X[np] += X[pid];
            Y[np] += Y[pid];
            Z[np] += Z[pid];
        }
        X[np] *= 0.25;
        Y[np] *= 0.25;
        Z[np] *= 0.25;

        E_Int nodes[4] = {pn[0], pn[1], pn[2], pn[3]};

        F[fid] = {nodes[0], nodes[1], np};
        F[nf] = {nodes[1], nodes[2], np};
        F[nf+1] = {nodes[2], nodes[3], np};
        F[nf+2] = {nodes[3], nodes[0], np};

        E_Int own = owner[fid];
        auto &pf = C[own];
        size_t stride = pf.size();
        pf.resize(stride + 3);
        for (int i = 0; i < 3; i++) pf[stride+i] = nf+i;

        np += 1;
        nf += 3;
    }

    assert((size_t)nf == F.size());
    assert((size_t)np == X.size());
}