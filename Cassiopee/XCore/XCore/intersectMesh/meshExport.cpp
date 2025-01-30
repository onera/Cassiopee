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

PyObject *IMesh::export_karray(E_Int remove_periodic) const
{
    if (remove_periodic) return export_karray_periodic();

    return export_karray_orig();
}

PyObject *IMesh::export_karray_periodic() const
{
    // Keep the cells whose tag is 1

    assert(ctag.size() == (size_t)nc);

    std::vector<E_Int> cells;

    for (E_Int cid = 0; cid < nc; cid++) {
        if (ctag[cid] == 1) cells.push_back(cid);
    }

    E_Int new_nc = (E_Int)cells.size();

    E_Int sizeNFace = 0;

    std::map<E_Int, E_Int> new_fids;
    E_Int new_nf = 0;

    for (E_Int i = 0; i < new_nc; i++) {
        E_Int cid = cells[i];
        const auto &pf = C[cid];

        sizeNFace += (E_Int)pf.size();

        for (const E_Int fid : pf) {
            if (new_fids.find(fid) == new_fids.end()) {
                new_fids[fid] = new_nf;
                new_nf++;
            }
        }
    }

    std::map<E_Int, E_Int> new_pids;
    E_Int new_np = 0;

    E_Int sizeNGon = 0;

    for (const auto &fdat : new_fids) {
        E_Int ofid = fdat.first;

        const auto &pn = F[ofid];

        sizeNGon += (E_Int)pn.size();

        for (const E_Int pid : pn) {
            if (new_pids.find(pid) == new_pids.end()) {
                new_pids[pid] = new_np;
                new_np++; 
            }
        }
    }

    const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

    PyObject *array = K_ARRAY::buildArray3(3, varString, new_np, new_nc,
        new_nf, "NGON", sizeNGon, sizeNFace, 3, false, 3);
    
    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(array, f, cn);

    // Make the coordinates

    E_Float *px = f->begin(1);
    E_Float *py = f->begin(2);
    E_Float *pz = f->begin(3);

    for (const auto &pdat : new_pids) {
        const E_Int opid = pdat.first;
        const E_Int npid = pdat.second;

        px[npid] = X[opid];
        py[npid] = Y[opid];
        pz[npid] = Z[opid];
    }

    // Make the element start offsets

    E_Int *indPG = cn->getIndPG();
    E_Int *ngon = cn->getNGon();
    E_Int *indPH = cn->getIndPH();
    E_Int *nface = cn->getNFace();

    indPG[0] = indPH[0] = 0;

    for (const auto &fdat : new_fids) {
        E_Int ofid = fdat.first;
        E_Int nfid = fdat.second;

        indPG[nfid+1] = (E_Int)F[ofid].size();
    }

    for (E_Int i = 0; i < new_nf; i++) indPG[i+1] += indPG[i];

    for (E_Int i = 0; i < new_nc; i++) {
        E_Int cid = cells[i];

        const auto &pf = C[cid];

        indPH[i+1] = (E_Int)pf.size();
    }

    for (E_Int i = 0; i < new_nc; i++) indPH[i+1] += indPH[i];

    // Make the element connectivities

    for (const auto &fdat : new_fids) {
        E_Int ofid = fdat.first;
        E_Int nfid = fdat.second;

        E_Int *ptr = &ngon[indPG[nfid]];

        const auto &pn = F[ofid];

        for (const E_Int opid : pn) *ptr++ = new_pids.at(opid) + 1;
    }

    for (E_Int i = 0; i < new_nc; i++) {
        E_Int cid = cells[i];

        E_Int *ptr = &nface[indPH[i]];
        
        const auto &pf = C[cid];

        for (const E_Int ofid : pf) *ptr++ = new_fids.at(ofid) + 1;
    }

    delete f;
    delete cn;

    return array;
}

PyObject *IMesh::export_karray_orig() const
{
    E_Int sizeNGon = 0, sizeNFace = 0;

    for (const auto &pn : F) sizeNGon += (E_Int)pn.size();
    for (const auto &pf : C) sizeNFace += (E_Int)pf.size();

    const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

    PyObject *array = K_ARRAY::buildArray3(3, varString, np, nc, nf, "NGON",
        sizeNGon, sizeNFace, 3, false, 3);
    
    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(array, f, cn);

    E_Float *px = f->begin(1);
    for (E_Int i = 0; i < np; i++) px[i] = X[i];
    E_Float *py = f->begin(2);
    for (E_Int i = 0; i < np; i++) py[i] = Y[i];
    E_Float *pz = f->begin(3);
    for (E_Int i = 0; i < np; i++) pz[i] = Z[i];

    E_Int *indPG = cn->getIndPG();
    E_Int *ngon = cn->getNGon();
    E_Int *indPH = cn->getIndPH();
    E_Int *nface = cn->getNFace();

    indPG[0] = indPH[0] = 0;
    for (E_Int i = 0; i < nf; i++) indPG[i+1] = indPG[i] + (E_Int)F[i].size();
    for (E_Int i = 0; i < nc; i++) indPH[i+1] = indPH[i] + (E_Int)C[i].size();

    assert(indPG[nf] == sizeNGon);
    assert(indPH[nc] == sizeNFace);

    E_Int *ptr = ngon;

    for (E_Int i = 0; i < nf; i++) {
        const auto &pn = F[i];
        for (E_Int p : pn) *ptr++ = p+1;
    }

    ptr = nface;

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (E_Int f : pf) *ptr++ = f+1;
    }

    delete f;
    delete cn;

    return array;
}

IMesh IMesh::extract_conformized()
{
    // Keep all the points
    std::vector<E_Float> new_X(X), new_Y(Y), new_Z(Z);

    // Conformize the faces

    std::vector<std::vector<E_Int>> new_F(factive.size());

    E_Int new_nf = 0;
    
    std::map<E_Int, E_Int> new_fids;

    for (E_Int face : factive) {
        new_fids[face] = new_nf;

        const auto &pn = F[face];

        auto &new_face = new_F[new_nf];

        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];

            std::list<E_Int> epoints;

            extract_edge_points(p, q, epoints);

            epoints.pop_back();

            for (auto it = epoints.begin(); it != epoints.end(); it++)
                new_face.push_back(*it);

            /*
            UEdge e(p, q);

            auto it = ecenter.find(e);

            if (it != ecenter.end()) {
                new_face.push_back(p);
                new_face.push_back(it->second);
            } else {
                new_face.push_back(p);
            }
            */
        }

        new_nf++;
    }

    // Update cell connectivity

    std::vector<std::vector<E_Int>> new_C(C.size());

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];

        auto &new_cell = new_C[i];

        for (E_Int face : pf) {

            if (face_is_active(face)) {
                new_cell.push_back(new_fids[face]);
            } else {
                std::vector<E_Int> fleaves;
                get_fleaves(face, fleaves);

                for (E_Int fleaf : fleaves)
                    new_cell.push_back(new_fids[fleaf]);
            }
        }
    }

    IMesh new_M;
    new_M.np = np;
    new_M.X = X;
    new_M.Y = Y;
    new_M.Z = Z;
    new_M.nf = new_nf;
    new_M.F = new_F;
    new_M.nc = nc;
    new_M.C = new_C;
    new_M.ctag = ctag;

    for (E_Int face : patch) {
        new_M.patch.insert(new_fids[face]);
        new_M.factive.insert(new_fids[face]);
    }

    return new_M;
}
