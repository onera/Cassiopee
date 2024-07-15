#include <stack>

#include "Mesh.h"
#include "../common/mem.h"

inline
Int Mesh_get_cnei(Mesh *M, Int cid, Int fid)
{
    assert(cid == M->owner[fid] || cid == M->neigh[fid]);
    return (M->owner[fid] == cid) ? M->neigh[fid] : M->owner[fid];
}

void Mesh_get_cneis(Mesh *M, Int cid, Int &nn, Int neis[24])
{
    Int *cell = Mesh_get_cell(M, cid);
    Int *crange = Mesh_get_crange(M, cid);
    Int cstride = M->cstride[cid];

    for (Int i = 0; i < cstride; i++) {
        Int *pf = cell + 4*i;

        for (Int j = 0; j < crange[i]; j++) {
            Int face = pf[j];
            Int nei = Mesh_get_cnei(M, cid, face);
            if (nei != -1) neis[nn++] = nei;
        }
    }
}

static
Int Mesh_smooth_cref_local(Mesh *M)
{
    std::stack<Int> stk;
    for (Int i = 0; i < M->nc; i++) {
        if (M->cref[i] > 0) stk.push(i);
    }

    while (!stk.empty()) {
        Int cid = stk.top();
        stk.pop();

        Int nn = 0;
        Int neis[24];
        Mesh_get_cneis(M, cid, nn, neis);

        Int incr_cell = M->cref[cid] + M->clevel[cid];

        for (Int i = 0; i < nn; i++) {
            Int nei = neis[i];

            Int incr_nei = M->cref[nei] + M->clevel[nei];

            Int diff = abs(incr_nei - incr_cell);

            if (diff <= 1) continue;

            Int cell_to_mod = incr_cell > incr_nei ? nei : cid;

            M->cref[cell_to_mod] += 1;

            stk.push(cell_to_mod);
        }
    }
}

Int Mesh_smooth_cref(Mesh *M)
{
    Int exchange = 0, max_exchanges = 10;

    while (++exchange <= max_exchanges) {
        Mesh_smooth_cref_local(M);

        MPI_Barrier(MPI_COMM_WORLD);

        Int lstop = 0;

        for (Int i = 0; i < M->npp; i++) {
            PPatch *P = &M->pps[i];

            for (Int j = 0; j < P->nf; j++) {
                Int own = M->owner[P->pf[j]];
                assert(own >= 0 && own < M->nc);
                P->sbuf_i[j] = M->cref[own] + M->clevel[own];
            }

            MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid,
                MPI_COMM_WORLD, &M->reqs[M->nrq++]);
            MPI_Irecv(P->rbuf_i, P->nf, XMPI_INT, P->nei, P->nei,
                MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        }

        Mesh_comm_waitall(M);

        for (Int i = 0; i < M->npp; i++) {
            PPatch *P = &M->pps[i];

            for (Int j = 0; j < P->nf; j++) {
                Int own = M->owner[P->pf[j]];
                Int oval = M->cref[own] + M->clevel[own];
                Int nval = P->rbuf_i[j];

                if (nval > oval + 1) {
                    M->cref[own] += 1;
                    lstop = 1;
                }
            }
        }

        Int gstop;
        MPI_Allreduce(&lstop, &gstop, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (gstop == 0) break;
    }


    if (exchange > max_exchanges) {
        merr("Smoothing exceeded max allowed iterations.");
        return 1;
    }

    return 0;
}
