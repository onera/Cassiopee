/*    
    Copyright 2013-2026 ONERA.

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
#include "smesh.h"

ierr Smesh_from_Karray(Smesh &Sf, const Karray &array)
{
    std::vector<E_Int> owner(array.nfaces(), -1);
    std::vector<E_Int> neigh(array.nfaces(), -1);
    E_Int nc = array.ncells();

    for (E_Int cid = 0; cid < nc; cid++) {
        E_Int nf = -1;
        E_Int *pf = array.get_cell(cid, nf);
        for (E_Int i = 0; i < nf; i++) {
            E_Int fid = pf[i]-1;
            if (owner[fid] == -1) {
                owner[fid] = cid;
            } else {
                if (neigh[fid] != -1) {
                    IERROR("Face %d is shared between multiple cells.\n", fid);
                    return IERR_NON_MANIFOLD_FACE;
                }
                neigh[fid] = cid;
            }
        }
    }

    std::vector<E_Int> skin;
    E_Int nf = array.nfaces();
    for (E_Int fid = 0; fid < nf; fid++) {
        if (neigh[fid] == -1)
            skin.push_back(fid);
    }

    Sf.np = 0;
    std::map<E_Int, E_Int> pmap;
    Sf.F.resize(skin.size());
    Sf.nf = skin.size();

    for (size_t i = 0; i < skin.size(); i++) {
        E_Int gfid = skin[i];
        E_Int np = -1;
        E_Int *pn = array.get_face(gfid, np);
        auto &F = Sf.F[i];
        F.resize(np, -1);
        for (E_Int j = 0; j < np; j++) {
            E_Int gpid = pn[j];
            auto it = pmap.find(gpid);
            if (it == pmap.end()) {
                pmap[gpid] = Sf.np;
                F[j] = Sf.np;
                Sf.np++;
            } else {
                F[j] = it->second;
            }
        }
    }

    Sf.P.resize(Sf.np);
    Sf.X.resize(Sf.np);
    Sf.Y.resize(Sf.np);
    Sf.Z.resize(Sf.np);
    for (const auto &pids : pmap) {
        E_Int gpid = pids.first - 1;
        E_Int lpid = pids.second;
        auto &P = Sf.P[lpid];
        P.x() = array.x[gpid];
        P.y() = array.y[gpid];
        P.z() = array.z[gpid];
        Sf.X[lpid] = array.x[gpid];
        Sf.Y[lpid] = array.y[gpid];
        Sf.Z[lpid] = array.z[gpid];
    }

    Sf.Fc = Sf.F;
    Sf.make_edge_data();

    return IERR_OK;
}

void Smesh::make_edge_data()
{
    F2E.resize(Fc.size());
    std::map<o_edge, E_Int, o_edge_cmp> edges;
    ne = 0;

    for (E_Int fid = 0; fid < nf; fid++) {
        auto &pn = Fc[fid];
        auto &pe = F2E[fid];
        pe.resize(pn.size());
        for (size_t i = 0; i < pn.size(); i++) {
            E_Int p = pn[i];
            E_Int q = pn[(i+1)%pn.size()];
            o_edge e(p, q);
            auto it = edges.find(e);
            if (it == edges.end()) {
                pe[i] = ne;
                edges[e] = ne;
                E.push_back(e);
                ne++;
            } else {
                pe[i] = it->second;
            }
        }
    }

    E2F.resize(ne);
    
    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pe = F2E[fid];
        
        for (E_Int eid : pe) {
            if (E2F[eid][0] == -1)
                E2F[eid][0] = fid;
            else
                E2F[eid][1] = fid;
        }
    }
}

ierr Smesh::check_manifold() const
{
    std::vector<int> count(ne, 0);
    for (const auto &pe : F2E) {
        for (auto eid : pe) {
            count[eid]++;
            if (count[eid] > 2) {
                IERROR("Edge %d is shared between more than two faces.\n", eid);
                return IERR_NON_MANIFOLD_EDGE;
            }
        }
    }
    return IERR_OK;
}

PyObject *K_XCORE::init_surf(PyObject *self, PyObject *args)
{
    PyObject *MESH;
    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray array;
    if (Karray_parse_ngon(MESH, array) != 0) {
        RAISE("Input mesh must be an NGon.");
        return NULL;
    }

    Smesh *Sf = new Smesh;
    if (Smesh_from_Karray(*Sf, array) != IERR_OK) {
        RAISE("Failed to create surface mesh from input NGon.");
        Karray_free_ngon(array);
        return NULL;
    }

    PyObject *hook = PyCapsule_New((void *)Sf, NULL, NULL);

    return hook;
}

PyObject *K_XCORE::check_surf_manifold(PyObject *self, PyObject *args)
{
    PyObject *SURF_HOOK;
    if (!PYPARSETUPLE_(args, O_, &SURF_HOOK)) {
        RAISE("Bad input.");
        return NULL;
    }

    Smesh *Sf = (Smesh *)PyCapsule_GetPointer(SURF_HOOK, NULL);

    bool is_manifold = Sf->check_manifold() != IERR_NON_MANIFOLD_EDGE;

    if (!is_manifold) {
        IERROR("Surface is non-manifold.\n");
    }

    return PyBool_FromLong(is_manifold);
}