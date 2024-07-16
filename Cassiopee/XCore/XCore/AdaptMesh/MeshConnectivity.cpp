/*    
    Copyright 2013-2024 Onera.

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
#include "Mesh.h"
#include "common/mem.h"
#include "Edge.h"

static inline
Int register_edge(Mesh *M, Int p, Int q, std::map<UEdge, Int> &emap, Int *pe,
    Int pos)
{
    UEdge E(p, q);

    auto it = emap.find(E);

    Int ret;

    if (it == emap.end()) {
        M->edges[M->ne].p = E.p;
        M->edges[M->ne].q = E.q;
        pe[pos] = M->ne;
        ret = M->ne;
        emap[E] = M->ne++;
    } else {
        pe[pos] = it->second;
        ret = it->second;
    }

    return ret;
}

static inline
void count_edge(Mesh *M, Int p, Int q, std::map<UEdge, Int> &emap, Int &ne)
{
    UEdge E(p, q);

    auto it = emap.find(E);

    if (it == emap.end()) emap[E] = ne++;
}

void Mesh_make_edge_connectivity(Mesh *M)
{
    M->ne = 0;

    // Count edges

    std::map<UEdge, Int> emap;

    Int ne = 0;

    for (Int fid = 0; fid < M->nf; fid++) {
        Int *face = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);
        Int fstride = M->fstride[fid];
        Int size = 2*fstride; 

        for (Int i = 0; i < size; i += 2) {
            Int p = face[i];
            Int q = face[i+1];
            Int r = face[(i+2)%size];

            if (frange[i/2] == 1) {
                count_edge(M, p, r, emap, ne);
                
            } else {
                count_edge(M, p, q, emap, ne);
                count_edge(M, q, r, emap, ne);
            }
        }
    }

    M->ne = ne;

    M->edges = (UEdge *)XMALLOC(M->ne * sizeof(UEdge));

    assert(M->fedg == NULL);
    M->fedg = IntArray(8 * M->nf);
    memset(M->fedg, -1, 8 * M->nf * sizeof(Int));

    emap.clear();

    M->xedf = IntArray(M->ne+1);

    M->elevel = IntArray(M->ne);
    memset(M->elevel, -1, M->ne * sizeof(Int));
 

    M->ne = 0;

    for (Int fid = 0; fid < M->nf; fid++) {
        Int *face = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);
        Int fstride = M->fstride[fid];

        Int size = fstride*2;

        Int *pe = Mesh_get_fedges(M, fid);

        for (Int i = 0; i < size; i += 2) {
            Int p = face[i];
            Int q = face[i+1];
            Int r = face[(i+2)%size];
            Int eid;

            assert(p != -1);

            if (frange[i/2] == 1) {
                assert(r != -1);
                eid = register_edge(M, p, r, emap, pe, i);
                M->xedf[eid+1]++;

                M->elevel[eid] = M->flevel[fid];
            } else {
                assert(frange[i/2] == 2);
                assert(q != -1);
                eid = register_edge(M, p, q, emap, pe, i);
                M->xedf[eid+1]++;

                M->elevel[eid] = M->flevel[fid]+1;

                assert(r != -1);
                eid = register_edge(M, q, r, emap, pe, i+1);
                M->xedf[eid+1]++;


                M->elevel[eid] = M->flevel[fid]+1;
            }
        }
    }

    assert(ne == M->ne);

    for (Int i = 0; i < M->ne; i++) M->xedf[i+1] += M->xedf[i];

    assert(M->E2F == NULL);
    M->E2F = IntArray(M->xedf[M->ne]);
    memset(M->E2F, -1, M->xedf[M->ne]);
    std::vector<Int> count(M->ne, 0);

    for (Int fid = 0; fid < M->nf; fid++) {
        Int *fedg = Mesh_get_fedges(M, fid);
        Int *frange = Mesh_get_frange(M, fid);

        for (Int i = 0; i < M->fstride[fid]; i++) {
            Int *pe = fedg + 2*i;

            for (Int j = 0; j < frange[i]; j++) {
                Int eid = pe[j];
                M->E2F[M->xedf[eid] + count[eid]++] = fid;
            }
        }
    }
}


void Mesh_update_global_cell_ids(Mesh *M)
{
    Int new_cells = M->nc - M->nc_old;

    Int first_new_cell = 0;
    MPI_Scan(&new_cells, &first_new_cell, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    Int gnc_old = 0;
    MPI_Allreduce(&M->nc_old, &gnc_old, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    Int gnc = 0;
    MPI_Allreduce(&M->nc, &gnc, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (new_cells > 0) {
        M->l2gc = (Int *)XRESIZE(M->l2gc, M->nc * sizeof(Int));

        Int gstart = gnc_old + first_new_cell - new_cells;

        for (Int i = 0; i < new_cells; i++) {
            Int lcid = M->nc_old + i;
            Int gcid = gstart + i;
            assert(gcid < gnc);
            assert(M->g2lc.find(gcid) == M->g2lc.end());
            M->l2gc[lcid] = gcid;
            M->g2lc[gcid] = lcid; 
        }
    }

    //MPI_Barrier(MPI_COMM_WORLD);
}

void Mesh_update_ppatches(Mesh *M)
{
    M->face_to_ppatch.clear();

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        
        // Count

        Int new_nf = 0;

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];

            Int fstate = M->fref[face];

            assert(fstate != FACE_NEW);

            if (fstate == FACE_UNTOUCHED) new_nf += 1;
            else new_nf += M->fchildren.at(face).size();
        }

        // Allocate

        Int *pf = IntArray(new_nf);
        Int *pn = IntArray(new_nf);
        P->sbuf_i = (Int *)XRESIZE(P->sbuf_i, new_nf * sizeof(Int));

        // Fill

        Int *ptr = pf;

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];

            if (M->fref[face] == 0) {
                *ptr++ = face;
            } else {
                const auto &children = M->fchildren.at(face);
                for (Int child : children) *ptr++ = child;
            }
        }

        // Reassign

        assert(ptr - pf == new_nf);

        XFREE(P->pf);
        XFREE(P->pn);

        P->nf = new_nf;
        P->pf = pf;
        P->pn = pn;

        for (Int j = 0; j < P->nf; j++) M->face_to_ppatch[P->pf[j]] = i;
    }

    // Switch children order for the bigger proc for iso mode only
    if (!M->mode_2D) {
        for (Int i = 0; i < M->npp; i++) {
            PPatch *P = &M->pps[i];

            if (M->pid < P->nei) continue;

            Int j = 0;

            for (; j < P->nf; ) {
                Int face = P->pf[j++];

                Int fstate = M->fref[face];

                assert(fstate != FACE_NEW);

                if (fstate == FACE_UNTOUCHED) continue;

                assert(fstate == FACE_REFINED);

                std::swap(P->pf[j], P->pf[j+2]);
                j += 3;
            }

            assert(j == P->nf);
        }
    }

    // Exchange remote cell neighbours

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        Int *ptr = P->sbuf_i;

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];
            Int own = M->owner[face];
            Int gcid = M->l2gc[own];
            *ptr++ = gcid;
        }

        MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        MPI_Irecv(P->pn, P->nf, XMPI_INT, P->nei, P->nei,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
    }

    Mesh_comm_waitall(M);
}

void Mesh_update_bpatches(Mesh *M)
{
    M->face_to_bpatch.clear();

    for (Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];
        
        // Count

        Int new_nf = 0;

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];

            Int fstate = M->fref[face];

            assert(fstate != FACE_NEW);

            if (fstate == FACE_UNTOUCHED) new_nf += 1;
            else new_nf += M->fchildren.at(face).size();
        }

        // Allocate

        Int *pf = IntArray(new_nf);

        // Fill

        Int *ptr = pf;

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];

            if (M->fref[face] == FACE_UNTOUCHED) {
                *ptr++ = face;
            } else {
                assert(M->fref[face] == FACE_REFINED);
                const auto &children = M->fchildren.at(face);
                for (Int child : children) *ptr++ = child;
            }
        }

        // Reassign

        assert(ptr - pf == new_nf);

        XFREE(P->pf);

        P->nf = new_nf;
        P->pf = pf;

        for (Int j = 0; j < P->nf; j++) M->face_to_bpatch[P->pf[j]] = i;
    }
}

void Mesh_update_global_face_ids(Mesh *M)
{
    M->l2gf = (Int *)XRESIZE(M->l2gf, M->nf * sizeof(Int)); 
    for (Int i = M->nf_old; i < M->nf; i++) M->l2gf[i] = -1; 

    Int nfree = M->nf;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        nfree -= P->nf;
    }

    Int ndup = M->nf - nfree;

    Int gnfree;
    MPI_Allreduce(&nfree, &gnfree, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    Int gndup;
    MPI_Allreduce(&ndup, &gndup, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    gndup /= 2;

    Int gstart = M->gnf_old;

    // On veut commencer à partir du nombre de faces global effectif (sans doublons)
    // Nombre de faces total = Nombre de faces internal/boundary total + pfaces/2
    
    // Count the number of new faces that do not belong to a patch
    Int nibf = 0;
    for (Int i = M->nf_old; i < M->nf; i++) {
        assert(M->fref[i] == FACE_NEW);
        if (M->face_to_ppatch.find(i) == M->face_to_ppatch.end()) {
            nibf++;
        }
    }

    Int first_ibf;
    MPI_Scan(&nibf, &first_ibf, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    gstart = gstart + first_ibf - nibf;

    for (Int i = M->nf_old; i < M->nf; i++) {
        Int lfid = i;
        // Skip patch face
        if (M->face_to_ppatch.find(lfid) != M->face_to_ppatch.end()) continue;
        Int gfid = gstart++;
        assert(M->g2lf.find(gfid) == M->g2lf.end());
        assert(M->l2gf[lfid] == -1);
        M->l2gf[lfid] = gfid;
        M->g2lf[gfid] = lfid; 
    }

    // Toutes mes anciennes/nouvelles faces proc
    
    Int tnp_new = 0, tnp_old = 0;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int lfid = P->pf[j];
            Int fstate = M->fref[lfid];
            if (fstate == FACE_NEW) tnp_new++;
            else tnp_old++;
        }
    }

    // Count the pfaces that I control

    Int my_npf = 0;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid > P->nei) continue;

        for (Int j = 0; j < P->nf; j++) {
            Int lfid = P->pf[j];
            Int fstate = M->fref[lfid];
            if (fstate == FACE_NEW) my_npf++;
        }
    }

    Int Npf;
    MPI_Scan(&my_npf, &Npf, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    Int max_gstart;
    MPI_Allreduce(&gstart, &max_gstart, 1, XMPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    gstart = max_gstart + Npf - my_npf;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid > P->nei) continue;

        for (Int j = 0; j < P->nf; j++) {
            Int lfid = P->pf[j];
            Int fstate = M->fref[lfid];
            if (fstate == FACE_NEW) {
                Int gfid = gstart++;
                assert(M->g2lf.find(gfid) == M->g2lf.end());
                assert(M->l2gf[lfid] == -1);
                M->l2gf[lfid] = gfid;
                M->g2lf[gfid] = lfid;
            }
        }
    }

    // Allocate
    
    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        P->sbuf_i = (Int *)XRESIZE(P->sbuf_i, P->nf * sizeof(Int));
        P->rbuf_i = (Int *)XRESIZE(P->rbuf_i, P->nf * sizeof(Int));
    }

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid < P->nei) {

            for (Int j = 0; j < P->nf; j++) {
                P->sbuf_i[j] = M->l2gf[P->pf[j]];
            }

            MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid,
                MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        } else {
            MPI_Irecv(P->rbuf_i, P->nf, XMPI_INT, P->nei, P->nei,
                MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        }
    }

    Mesh_comm_waitall(M);


    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid < P->nei) continue;

        for (Int j = 0; j < P->nf; j++) {

            Int rfid = P->rbuf_i[j];

            Int lfid = P->pf[j];
            Int fstate = M->fref[lfid];

            if (fstate != FACE_NEW) {
                continue;
            }

            assert(M->g2lf.find(rfid) == M->g2lf.end());

            assert(M->l2gf[lfid] == -1);

            M->l2gf[lfid] = rfid;
            M->g2lf[rfid] = lfid;
        }
    }
}


void Mesh_make_cell_cells(Mesh *M)
{
    assert(M->owner);
    assert(M->neigh);

    assert(M->xneis == NULL);
    assert(M->cneis == NULL);

    M->xneis = IntArray(M->nc + 1);

    // Internal faces
    Int nif = 0;

    for (Int i = 0; i < M->nf; i++) {
        Int nei = M->neigh[i];
        if (nei == -1) continue;

        Int own = M->owner[i];
        M->xneis[own+1]++;
        M->xneis[nei+1]++;

        nif++;
    }

    // PPatch faces

    Int npf = 0;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        npf += P->nf;

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];
            assert(face >= 0 && face < M->nf);
            M->xneis[M->owner[face]+1]++;
        }
    }

    for (Int i = 0; i < M->nc; i++) M->xneis[i+1] += M->xneis[i];

    M->cneis = IntArray(2 * nif + npf);

    Int *count = IntArray(M->nc);

    // Local cell neighbours first

    for (Int i = 0; i < M->nf; i++) {
        Int nei = M->neigh[i];
        
        if (nei == -1) continue;

        Int own = M->owner[i];

        M->cneis[M->xneis[own] + count[own]++] = M->l2gc[nei];
        M->cneis[M->xneis[nei] + count[nei]++] = M->l2gc[own];
    }

    // Remote cell neighbours
    
    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int face =  P->pf[j];
            Int own = M->owner[face];
            M->cneis[M->xneis[own] + count[own]++] = P->pn[j];
        }
    }
    
    XFREE(count);
}
