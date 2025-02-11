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
#include <stack>

#include "Mesh.h"
#include "common/mem.h"

void Mesh_get_cneis(Mesh *M, E_Int cid, E_Int &nn, E_Int neis[24])
{
    E_Int *cell = Mesh_get_cell(M, cid);
    E_Int *crange = Mesh_get_crange(M, cid);
    E_Int cstride = M->cstride[cid];

    nn = 0;

    for (E_Int i = 0; i < cstride; i++) {
        E_Int *pf = cell + 4*i;

        for (E_Int j = 0; j < crange[i]; j++) {
            E_Int face = pf[j];
            E_Int nei = Mesh_get_cnei(M, cid, face);
            if (nei != -1 && nei != neis[nn])
                neis[nn++] = nei;
        }
    }
}

static
E_Int Mesh_smooth_cref_local(Mesh *M)
{
    std::stack<E_Int> stk;
    for (E_Int i = 0; i < M->nc; i++) {
        if (M->cref[i] > 0) stk.push(i);
    }

    while (!stk.empty()) {
        E_Int cid = stk.top();
        stk.pop();

        E_Int nn, neis[24];
        Mesh_get_cneis(M, cid, nn, neis);

        //E_Int incr_cell = M->cref[cid] + M->clevel[cid];

        for (E_Int i = 0; i < nn; i++) {
            E_Int nei = neis[i];

            E_Int incr_cell = M->cref[cid] + M->clevel[cid];
            E_Int incr_nei = M->cref[nei] + M->clevel[nei];

            E_Int diff = abs(incr_nei - incr_cell);

            if (diff <= 1) continue;

            E_Int cell_to_mod = incr_cell > incr_nei ? nei : cid;

            M->cref[cell_to_mod] += diff-1;
            //M->cref[cell_to_mod] += 1;

            stk.push(cell_to_mod);
        }
    }
    return 0;
}

E_Int Mesh_smooth_cref(Mesh *M)
{
    E_Int exchange = 0, max_exchanges = 10;

    while (++exchange <= max_exchanges) {
        Mesh_smooth_cref_local(M);

        MPI_Barrier(MPI_COMM_WORLD);

        E_Int lstop = 0;

        for (E_Int i = 0; i < M->npp; i++) {
            PPatch *P = &M->pps[i];

            for (E_Int j = 0; j < P->nf; j++) {
                E_Int own = M->owner[P->pf[j]];
                assert(own >= 0 && own < M->nc);
                P->sbuf_i[j] = M->cref[own] + M->clevel[own];
            }

            MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid,
                MPI_COMM_WORLD, &M->reqs[M->nrq++]);
            MPI_Irecv(P->rbuf_i, P->nf, XMPI_INT, P->nei, P->nei,
                MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        }

        Mesh_comm_waitall(M);

        for (E_Int i = 0; i < M->npp; i++) {
            PPatch *P = &M->pps[i];

            for (E_Int j = 0; j < P->nf; j++) {
                E_Int own = M->owner[P->pf[j]];
                E_Int oval = M->cref[own] + M->clevel[own];
                E_Int nval = P->rbuf_i[j];

                if (nval > oval + 1) {
                    M->cref[own] += 1;
                    lstop = 1;
                }
            }
        }

        E_Int gstop = 0;
        MPI_Allreduce(&lstop, &gstop, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (gstop == 0) break;
    }


    if (exchange > max_exchanges) {
        merr("Smoothing exceeded max allowed iterations.");
        return 1;
    }

    return 0;
}

// TODO(Imad): exactly the same as Mesh_smooth_cref...
void Mesh_smooth_cell_refinement_data(Mesh *M)
{
    E_Int nc = M->nc;

    std::stack<E_Int> stk;

    for (E_Int cid = 0; cid < nc; cid++) {
        if (M->cref[cid] > 0)
            stk.push(cid);
    }

    while (!stk.empty()) {
        E_Int cid = stk.top();
        stk.pop();

        E_Int nn, neis[24];
        Mesh_get_cneis(M, cid, nn, neis);

        for (E_Int i = 0; i < nn; i++) {
            E_Int nei = neis[i];
            E_Int incr_nei = M->cref[nei] + M->clevel[nei];
            E_Int incr_cid = M->cref[cid] + M->clevel[cid];
            E_Int diff = abs(incr_nei - incr_cid);
            if (diff <= 1) continue;
            E_Int idx_to_modify = incr_cid > incr_nei ? nei : cid;
            M->cref[idx_to_modify] += diff-1;
            stk.push(idx_to_modify);
        }
    }
}
