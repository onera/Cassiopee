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
E_Int register_edge(Mesh *M, E_Int p, E_Int q, std::map<UEdge, E_Int> &emap, E_Int *pe,
    E_Int pos)
{
    UEdge E(p, q);

    auto it = emap.find(E);

    E_Int ret;

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
void count_edge(Mesh *M, E_Int p, E_Int q, std::map<UEdge, E_Int> &emap, E_Int &ne)
{
    UEdge E(p, q);

    auto it = emap.find(E);

    if (it == emap.end()) emap[E] = ne++;
}

void Mesh_make_edge_connectivity(Mesh *M)
{
    M->ne = 0;

    // Count edges

    std::map<UEdge, E_Int> emap;

    E_Int ne = 0;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        E_Int fstride = M->fstride[fid];
        E_Int size = 2*fstride; 

        for (E_Int i = 0; i < size; i += 2) {
            E_Int p = face[i];
            E_Int q = face[i+1];
            E_Int r = face[(i+2)%size];

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
    memset(M->fedg, -1, 8 * M->nf * sizeof(E_Int));

    emap.clear();

    M->xedf = IntArray(M->ne+1);

    M->elevel = IntArray(M->ne);
    memset(M->elevel, -1, M->ne * sizeof(E_Int));
 

    M->ne = 0;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        E_Int fstride = M->fstride[fid];

        E_Int size = fstride*2;

        E_Int *pe = Mesh_get_fedges(M, fid);

        for (E_Int i = 0; i < size; i += 2) {
            E_Int p = face[i];
            E_Int q = face[i+1];
            E_Int r = face[(i+2)%size];
            E_Int eid;

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

    for (E_Int i = 0; i < M->ne; i++) M->xedf[i+1] += M->xedf[i];

    assert(M->E2F == NULL);
    M->E2F = IntArray(M->xedf[M->ne]);
    memset(M->E2F, -1, M->xedf[M->ne]);
    std::vector<E_Int> count(M->ne, 0);

    for (E_Int fid = 0; fid < M->nf; fid++) {
        E_Int *fedg = Mesh_get_fedges(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);

        for (E_Int i = 0; i < M->fstride[fid]; i++) {
            E_Int *pe = fedg + 2*i;

            for (E_Int j = 0; j < frange[i]; j++) {
                E_Int eid = pe[j];
                M->E2F[M->xedf[eid] + count[eid]++] = fid;
            }
        }
    }
}

void Mesh_update_global_cell_ids(Mesh *M)
{
    if (M->npc == 1) return;

    E_Int new_cells = M->nc - M->nc_old;

    E_Int first_new_cell = new_cells;
    MPI_Scan(&new_cells, &first_new_cell, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    E_Int gnc_old = 0;
    MPI_Allreduce(&M->nc_old, &gnc_old, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    E_Int gnc = 0;
    MPI_Allreduce(&M->nc, &gnc, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (new_cells > 0) {
        M->l2gc = (E_Int *)XRESIZE(M->l2gc, M->nc * sizeof(E_Int));

        E_Int gstart = gnc_old + first_new_cell - new_cells;

        for (E_Int i = 0; i < new_cells; i++) {
            E_Int lcid = M->nc_old + i;
            E_Int gcid = gstart + i;
            assert(gcid < gnc);
            assert(M->g2lc.find(gcid) == M->g2lc.end());
            M->l2gc[lcid] = gcid;
            M->g2lc[gcid] = lcid; 
        }
    }
}

void Mesh_update_ppatches(Mesh *M)
{
    M->face_to_ppatch.clear();

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        
        // Count

        E_Int new_nf = 0;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face = P->pf[j];

            E_Int fstate = M->fref[face];

            assert(fstate != FACE_NEW);

            if (fstate == FACE_UNTOUCHED) new_nf += 1;
            else new_nf += M->fchildren.at(face).size();
        }

        // Allocate

        E_Int *pf = IntArray(new_nf);
        E_Int *pn = IntArray(new_nf);
        P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, new_nf * sizeof(E_Int));

        // Fill

        E_Int *ptr = pf;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face = P->pf[j];

            if (M->fref[face] == 0) {
                *ptr++ = face;
            } else {
                const auto &children = M->fchildren.at(face);
                for (E_Int child : children) *ptr++ = child;
            }
        }

        // Reassign

        assert(ptr - pf == new_nf);

        XFREE(P->pf);
        XFREE(P->pn);

        P->nf = new_nf;
        P->pf = pf;
        P->pn = pn;

        for (E_Int j = 0; j < P->nf; j++) M->face_to_ppatch[P->pf[j]] = i;
    }

    // Switch children order for the bigger proc for iso mode only
    if (!M->mode_2D) {
        for (E_Int i = 0; i < M->npp; i++) {
            PPatch *P = &M->pps[i];

            if (M->pid < P->nei) continue;

            E_Int j = 0;

            for (; j < P->nf; ) {
                E_Int face = P->pf[j++];

                E_Int fstate = M->fref[face];

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

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        E_Int *ptr = P->sbuf_i;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face = P->pf[j];
            E_Int own = M->owner[face];
            E_Int gcid = M->l2gc[own];
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

    for (E_Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];
        
        // Count

        E_Int new_nf = 0;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face = P->pf[j];

            E_Int fstate = M->fref[face];

            assert(fstate != FACE_NEW);

            if (fstate == FACE_UNTOUCHED) new_nf += 1;
            else new_nf += M->fchildren.at(face).size();
        }

        // Allocate

        E_Int *pf = IntArray(new_nf);

        // Fill

        E_Int *ptr = pf;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face = P->pf[j];

            if (M->fref[face] == FACE_UNTOUCHED) {
                *ptr++ = face;
            } else {
                assert(M->fref[face] == FACE_REFINED);
                const auto &children = M->fchildren.at(face);
                for (E_Int child : children) *ptr++ = child;
            }
        }

        // Reassign

        assert(ptr - pf == new_nf);

        XFREE(P->pf);

        P->nf = new_nf;
        P->pf = pf;

        for (E_Int j = 0; j < P->nf; j++) M->face_to_bpatch[P->pf[j]] = i;
    }
}

void Mesh_update_global_face_ids(Mesh *M)
{
    if (M->npc == 1) return;

    M->l2gf = (E_Int *)XRESIZE(M->l2gf, M->nf * sizeof(E_Int)); 
    for (E_Int i = M->nf_old; i < M->nf; i++) M->l2gf[i] = -1; 

    E_Int nfree = M->nf;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        nfree -= P->nf;
    }

    E_Int ndup = M->nf - nfree;

    E_Int gnfree = nfree;
    MPI_Allreduce(&nfree, &gnfree, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    E_Int gndup = ndup;
    MPI_Allreduce(&ndup, &gndup, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    gndup /= 2;

    E_Int gstart = M->gnf_old;

    // On veut commencer à partir du nombre de faces global effectif (sans doublons)
    // Nombre de faces total = Nombre de faces internal/boundary total + pfaces/2
    
    // Count the number of new faces that do not belong to a patch
    E_Int nibf = 0;
    for (E_Int i = M->nf_old; i < M->nf; i++) {
        assert(M->fref[i] == FACE_NEW);
        if (M->face_to_ppatch.find(i) == M->face_to_ppatch.end()) {
            nibf++;
        }
    }

    E_Int first_ibf;
    MPI_Scan(&nibf, &first_ibf, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    gstart = gstart + first_ibf - nibf;

    for (E_Int i = M->nf_old; i < M->nf; i++) {
        E_Int lfid = i;
        // Skip patch face
        if (M->face_to_ppatch.find(lfid) != M->face_to_ppatch.end()) continue;
        E_Int gfid = gstart++;
        assert(M->g2lf.find(gfid) == M->g2lf.end());
        assert(M->l2gf[lfid] == -1);
        M->l2gf[lfid] = gfid;
        M->g2lf[gfid] = lfid; 
    }

    // Toutes mes anciennes/nouvelles faces proc
    
    E_Int tnp_new = 0, tnp_old = 0;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            E_Int fstate = M->fref[lfid];
            if (fstate == FACE_NEW) tnp_new++;
            else tnp_old++;
        }
    }

    // Count the pfaces that I control

    E_Int my_npf = 0;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid > P->nei) continue;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            E_Int fstate = M->fref[lfid];
            if (fstate == FACE_NEW) my_npf++;
        }
    }

    E_Int Npf;
    MPI_Scan(&my_npf, &Npf, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    E_Int max_gstart;
    MPI_Allreduce(&gstart, &max_gstart, 1, XMPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    gstart = max_gstart + Npf - my_npf;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid > P->nei) continue;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            E_Int fstate = M->fref[lfid];
            if (fstate == FACE_NEW) {
                E_Int gfid = gstart++;
                assert(M->g2lf.find(gfid) == M->g2lf.end());
                assert(M->l2gf[lfid] == -1);
                M->l2gf[lfid] = gfid;
                M->g2lf[gfid] = lfid;
            }
        }
    }

    // Allocate
    
    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, P->nf * sizeof(E_Int));
        P->rbuf_i = (E_Int *)XRESIZE(P->rbuf_i, P->nf * sizeof(E_Int));
    }

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid < P->nei) {

            for (E_Int j = 0; j < P->nf; j++) {
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


    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        if (M->pid < P->nei) continue;

        for (E_Int j = 0; j < P->nf; j++) {

            E_Int rfid = P->rbuf_i[j];

            E_Int lfid = P->pf[j];
            E_Int fstate = M->fref[lfid];

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

    // E_Internal faces
    E_Int nif = 0;

    for (E_Int i = 0; i < M->nf; i++) {
        E_Int nei = M->neigh[i];
        if (nei == -1) continue;

        E_Int own = M->owner[i];
        M->xneis[own+1]++;
        M->xneis[nei+1]++;

        nif++;
    }

    // PPatch faces

    E_Int npf = 0;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        npf += P->nf;

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face = P->pf[j];
            assert(face >= 0 && face < M->nf);
            M->xneis[M->owner[face]+1]++;
        }
    }

    for (E_Int i = 0; i < M->nc; i++) M->xneis[i+1] += M->xneis[i];

    M->cneis = IntArray(2 * nif + npf);

    E_Int *count = IntArray(M->nc);

    // Local cell neighbours first

    for (E_Int i = 0; i < M->nf; i++) {
        E_Int nei = M->neigh[i];
        
        if (nei == -1) continue;

        E_Int own = M->owner[i];

        M->cneis[M->xneis[own] + count[own]++] = M->l2gc[nei];
        M->cneis[M->xneis[nei] + count[nei]++] = M->l2gc[own];
    }

    // Remote cell neighbours
    
    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int face =  P->pf[j];
            E_Int own = M->owner[face];
            M->cneis[M->xneis[own] + count[own]++] = P->pn[j];
        }
    }
    
    XFREE(count);
}

E_Int Mesh_get_global_face_count(Mesh *M)
{
    E_Int lcount = M->nf;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        if (M->pid < P->nei) continue;

        lcount -= P->nf;
    }

    E_Int gcount;
    MPI_Allreduce(&lcount, &gcount, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return gcount;
}