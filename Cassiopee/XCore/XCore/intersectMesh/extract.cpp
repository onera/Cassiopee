#include "mesh.h"
#include "common/Karray.h"

PyObject *K_XCORE::extractCell(PyObject *self, PyObject *args)
{
    PyObject *MESH;
    E_Int cid;
  
    if (!PYPARSETUPLE_(args, O_ I_, &MESH, &cid)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;

    E_Int ret;

    ret = Karray_parse_ngon(MESH, marray);

    if (ret != 0) return NULL;

    IMesh M(marray);

    IMesh M_out;

    auto &C = M_out.C;
    auto &F = M_out.F;
    auto &X = M_out.X;
    auto &Y = M_out.Y;
    auto &Z = M_out.Z;

    E_Int NC = 1;

    E_Int NF = M.C[cid].size();
    
    std::vector<E_Int> cell(NF);

    for (E_Int i = 0; i < NF; i++) cell[i] = i;

    C.push_back(cell);

    for (E_Int fid : M.C[cid]) {
        F.push_back(M.F[fid]);
    }

    std::map<E_Int, E_Int> new_pids;

    E_Int NP = 0;

    for (auto &pn : F) {
        for (E_Int &p : pn) {
            auto it = new_pids.find(p);
            
            if (it == new_pids.end()) {
                new_pids[p] = NP;
                p = NP;
                NP++;
            } else {
                p = it->second;
            }
        }
    }

    X.resize(NP);
    Y.resize(NP);
    Z.resize(NP);

    for (const auto &pdat : new_pids) {
        E_Int opid = pdat.first;
        E_Int npid = pdat.second;
        X[npid] = M.X[opid];
        Y[npid] = M.Y[opid];
        Z[npid] = M.Z[opid];
    }

    M_out.nc = NC;
    M_out.nf = NF;
    M_out.np = NP;

    /*
    for (auto &cell : M_out.C) {
        for (auto f : cell) printf("%d ", f);
        puts("");
    }

    for (auto &face : M_out.F) {
        for (auto f : cell) printf("%d ", f);
        puts("");
    }
    */


    return M_out.export_karray();
}
