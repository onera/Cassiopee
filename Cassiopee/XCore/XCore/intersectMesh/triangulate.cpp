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

        std::vector<PyArrayObject *> new_bcs = M.triangulate_skin(bcs_in, fid_to_bc);

        PyObject *out = PyList_New(0);
        PyList_Append(out, M.export_karray());
        PyObject *bcs_out = PyList_New(0);
        for (const auto &new_bc : new_bcs) {
            PyList_Append(bcs_out, (PyObject *)new_bc);
            Py_DECREF(new_bc);
        }
        PyList_Append(out, bcs_out);
        Py_DECREF(bcs_out);
        return out;
    } else {
        M.triangulate_skin();
        return M.export_karray();
    }
}

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

void IMesh::triangulate_skin()
{
    E_Int NF = nf;

    owner.resize(nf + skin.size(), -1);
    neigh.resize(nf + skin.size(), -1);

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

        nf++;
    }

    for (E_Int i = NF; i < nf; i++)
        skin.push_back(i);
}