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
    PyObject *MESH;
    if (!PYPARSETUPLE_(args, O_, &MESH)) {
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
    M.triangulate_skin();
    return M.export_karray();
}

/*
std::vector<PyArrayObject *> IMesh::triangulate_skin(const std::vector<Py_BC> &bcs_in,
    const std::unordered_map<E_Int, E_Int> &fid_to_bc)
{
    E_Int NF = nf;

    owner.resize(nf + skin.size(), -1);
    neigh.resize(nf + skin.size(), -1);

    std::vector<PyArrayObject *> new_bcs(bcs_in.size());
    std::vector<E_Int *> ptrs(bcs_in.size());
    std::vector<E_Int> counts(bcs_in.size(), 0);

    for (size_t i = 0; i < bcs_in.size(); i++) {
        E_Int new_size = bcs_in[i].size * 2;
        npy_intp dims[2];
        dims[0] = (npy_intp)new_size;
        dims[1] = 1;
        new_bcs[i] = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
        ptrs[i] = (E_Int *)PyArray_DATA(new_bcs[i]);
    }

    for (auto fid : skin) {
        assert(neigh[fid] == -1);

        const auto &pn = F[fid];
        if (pn.size() == 3) continue;

        assert(pn.size() == 4);

        E_Int nodes[4] = {pn[0], pn[1], pn[2], pn[3]};

        F.push_back({nodes[2], nodes[3], nodes[0]});
        F[fid] = {nodes[0], nodes[1], nodes[2]};
        assert(F[fid].size() == 3);

        E_Int own = owner[fid];
        auto &pf = C[own];
        pf.push_back(nf);

        owner[nf] = own;

        // Which bc (if any) does fid belong to?
        auto it = fid_to_bc.find(fid+1);
        if (it != fid_to_bc.end()) {
            E_Int bc_id = it->second;
            ptrs[bc_id][counts[bc_id]++] = fid+1;
            ptrs[bc_id][counts[bc_id]++] = nf+1;
        }

        nf++;
    }

    for (E_Int i = NF; i < nf; i++)
        skin.push_back(i);

    for (size_t i = 0; i < counts.size(); i++) {
        assert(counts[i] == 2*bcs_in[i].size);
    }

    return new_bcs;
}
*/

void IMesh::triangulate_skin()
{
    owner.clear();
    neigh.clear();
    owner.resize(nf, -1);
    neigh.resize(nf, -1);

    for (size_t cid = 0; cid < nc; cid++) {
        const auto &pf = C[cid];
        for (auto fid : pf) {
            if (owner[fid] == -1) owner[fid] = cid;
            else {
                assert(neigh[fid] == -1);
                neigh[fid] = cid;
            }
        }
    }

    skin.clear();
    for (size_t fid = 0; fid < nf; fid++) {
        if (neigh[fid] == -1) skin.push_back(fid);
    }

    // One new point per triangulated face
    E_Int NP = np + skin.size();

    // Three new faces per triangulated face
    E_Int NF = nf + 3*skin.size();

    for (auto fid : skin) {
        const auto &pn = F[fid];
        if (pn.size() == 3) continue;
        
        assert(pn.size() == 4);
        
        E_Int nodes[4] = {pn[0], pn[1], pn[2], pn[3]};

        E_Float cx, cy, cz;
        cx = cy = cz = 0;
        for (int i = 0; i < 4; i++) {
            cx += X[nodes[i]];
            cy += Y[nodes[i]];
            cz += Z[nodes[i]];
        }
        cx *= 0.25;
        cy *= 0.25;
        cz *= 0.25;

        F[fid] = {nodes[0], nodes[1], np};
        F.push_back({nodes[1], nodes[2], np});
        F.push_back({nodes[2], nodes[3], np});
        F.push_back({nodes[3], nodes[0], np});
        
        E_Int own = owner[fid];
        auto &pf = C[own];
        for (int i = 0; i < 3; i++) {
            pf.push_back(nf+i);
            owner.push_back(own);
            neigh.push_back(-1);
        }

        X.push_back(cx);
        Y.push_back(cy);
        Z.push_back(cz);

        nf += 3;
        np += 1;
    }

    assert(np == NP);
    assert(nf == NF);
}
