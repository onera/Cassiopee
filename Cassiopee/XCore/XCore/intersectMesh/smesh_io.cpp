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
#include "smesh.h"

void Smesh::write_edges(const char *fname, const std::set<E_Int> &eids) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", eids.size()*2);
    for (E_Int eid : eids) {
        const auto &e = E[eid];
        E_Int p = e.p;
        E_Int q = e.q;
        fprintf(fh, "%f %f %f\n", X[p], Y[p], Z[p]);
        fprintf(fh, "%f %f %f\n", X[q], Y[q], Z[q]);
    }
    fprintf(fh, "EDGES\n");
    fprintf(fh, "%zu\n", eids.size());
    for (size_t i = 0; i < 2*eids.size(); i++) {
        fprintf(fh, "%zu ", i);
    }
    fprintf(fh, "\n");

    fclose(fh);
}

void Smesh::write_points(const char *fname, const std::vector<E_Int> &pids) const
{
    FILE *fh= fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", pids.size());
    for (E_Int pid : pids){
        fprintf(fh, "%f %f %f\n", X[pid], Y[pid], Z[pid]);
    }

    fclose(fh);
}

void Smesh::write_points(const char *fname, const std::set<E_Int> &pids) const
{
    FILE *fh= fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", pids.size());
    for (E_Int pid : pids){
        fprintf(fh, "%f %f %f\n", X[pid], Y[pid], Z[pid]);
    }

    fclose(fh);
}

void Smesh::write_ngon(const char *fname, const std::set<E_Int> &fset) const
{
    std::vector<E_Int> flist;
    flist.reserve(fset.size());
    for (E_Int fid : fset) flist.push_back(fid);
    write_ngon(fname, flist);
}

void Smesh::write_ngon(const char *fname, const std::vector<E_Int> &faces) const
{
    std::vector<E_Int> INDPH(faces.size() + 1);
    INDPH[0] = 0;

    std::map<E_Int, E_Int> new_pids;
    std::map<E_Int, E_Int> new_eids;
    E_Int idx = 0;
    E_Int NP = 0;
    E_Int NE = 0;
    E_Int NF = (E_Int)faces.size();

    for (E_Int fid : faces) {
        const auto &pn = Fc[fid];
        const auto &pe = F2E[fid];
        
        INDPH[idx+1] = INDPH[idx] + (E_Int)pn.size();
        idx++;

        for (E_Int pid : pn) {
            if (new_pids.find(pid) == new_pids.end()) {
                new_pids[pid] = NP;
                NP++;
            }
        }

        for (E_Int eid : pe) {
            if (new_eids.find(eid) == new_eids.end()) {
                new_eids[eid] = NE;
                NE++;
            }
        }
    }

    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", new_pids.size());

    std::vector<E_Float> nX(NP), nY(NP), nZ(NP);
    for (const auto &pids : new_pids) {
        E_Int opid = pids.first;
        E_Int npid = pids.second;
        nX[npid] = X[opid];
        nY[npid] = Y[opid];
        nZ[npid] = Z[opid];
    }

    for (E_Int pid = 0; pid < NP; pid++) {
        fprintf(fh, "%f %f %f\n", nX[pid], nY[pid], nZ[pid]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, "%zu\n", new_eids.size()+1);
    E_Int sizeNGon = -2;
    for (size_t i = 0; i < new_eids.size() + 1; i++) {
        sizeNGon += 2;
        fprintf(fh, "%d ", sizeNGon);
    }
    fprintf(fh, "\n");
    assert(sizeNGon == 2*NE);

    std::vector<o_edge> nE(new_eids.size(), {-1, -1});
    for (const auto &eids : new_eids) {
        E_Int oeid = eids.first;
        E_Int neid = eids.second;
        nE[neid].p = new_pids[E[oeid].p];
        nE[neid].q = new_pids[E[oeid].q];
    }
    
    fprintf(fh, "NGON\n");
    fprintf(fh, "%d\n", 2*NE);
    for (const auto &e : nE) {
        fprintf(fh, "%d %d ", e.p, e.q);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, "%d\n", NF+1);
    for (E_Int i = 0; i < NF+1; i++)
        fprintf(fh, "%d ", INDPH[i]);
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, "%d\n", INDPH[NF]);
    for (size_t i = 0; i < faces.size(); i++) {
        const auto &pe = F2E[faces[i]];
        for (E_Int eid : pe) {
            fprintf(fh, "%d ", new_eids[eid]);
        }
    }
    fprintf(fh, "\n");

    fclose(fh);
}

void Smesh::write_face(const char *fname, E_Int fid) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    const auto &pn = Fc[fid];
    E_Int np = (size_t)pn.size();

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", np);
    for (E_Int p : pn) {
        fprintf(fh, "%f %f %f\n", X[p], Y[p], Z[p]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, "%d\n", np + 1);
    E_Int sizeNGon = -2;
    for (E_Int i = 0; i < np+1; i++) {
        sizeNGon += 2;
        fprintf(fh, "%d ", sizeNGon);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, "%d\n", 2*np);
    for (E_Int i = 0; i < np; i++) {
        fprintf(fh, "%d %d ", i, (i+1)%np);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, "%d\n", 2);
    fprintf(fh, "0 %d\n", np);

    fprintf(fh, "NFACE\n");
    fprintf(fh, "%d\n", np);
    for (E_Int i = 0; i < np; i++) fprintf(fh, "%d ", i);

    fclose(fh);
}

void Smesh::write_ngon(const char *fname) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (E_Int i = 0; i < np; i++) {
        fprintf(fh, "%f %f %f\n", X[i], Y[i], Z[i]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", ne+1);
    size_t sizeNGon = 0;
    fprintf(fh, "0 ");
    for (E_Int i = 0; i < ne; i++) {
        sizeNGon += 2;
        fprintf(fh, "%zu ", sizeNGon);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, "%zu\n", sizeNGon);
    for (E_Int i = 0; i < ne; i++) {
        fprintf(fh, SF_D_  " " SF_D_  " ", E[i].p, E[i].q);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    size_t sizeNFace = 0;
    fprintf(fh, "0 ");
    for (E_Int i = 0; i < nf; i++) {
        sizeNFace += F2E[i].size();
        fprintf(fh, "%zu ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, "%zu\n", sizeNFace);
    for (E_Int i = 0; i < nf; i++) {
        for (E_Int e : F2E[i]) fprintf(fh, SF_D_ " ", e);
    }
    fprintf(fh, "\n");

    fclose(fh);
}



