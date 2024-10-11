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

void IMesh::triangulate_skin()
{
    E_Int NF = nf;

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

        nf++;
    }

    for (E_Int i = NF; i < nf; i++)
        skin.push_back(i);
}
