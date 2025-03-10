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
