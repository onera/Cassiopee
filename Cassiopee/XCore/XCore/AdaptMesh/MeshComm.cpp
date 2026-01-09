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
#include <numeric>
#include <limits>
#include <unordered_map>

#include "Mesh.h"
#include "common/mem.h"
#include "scotch/ptscotch.h"

static
E_Int *compute_cell_weights(Mesh *M)
{
    E_Int *cwgts = IntArray(M->nc);

    E_Int gnc;
    MPI_Allreduce(&M->nc, &gnc, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    std::vector<E_Float> weights(M->nc, 0);

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int cval = M->cref[i];

        weights[i] = 1 << (3*cval);
    }

    E_Float min_weight = std::numeric_limits<E_Float>::max();

    for (E_Int i = 0; i < M->nc; i++) {
        if (min_weight > weights[i]) min_weight = weights[i];
    }

    E_Float gmin_weight;
    MPI_Allreduce(&min_weight, &gmin_weight, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    E_Float sum = 0.0;
    for (E_Int weight : weights) sum += weight;

    E_Float gsum;
    MPI_Allreduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    E_Float velotabsum = gsum / min_weight;

    E_Float upper_range = std::numeric_limits<SCOTCH_Num>::max() - 1;

    E_Float range_scale = 1;

    if (velotabsum > upper_range) {
        range_scale = 0.9 * upper_range / velotabsum;
        puts("Sum of weights overflows SCOTCH_Num");
    }

    for (E_Int i = 0; i < M->nc; i++) {
        cwgts[i] = ((weights[i]/min_weight - 1)*range_scale) + 1;
    }

    return cwgts;
}

static
E_Int *map_cell_graph(Mesh *M, E_Int *cwgts)
{
    assert(sizeof(SCOTCH_Num) <= sizeof(SCOTCH_Idx));
    assert(sizeof(SCOTCH_Idx) >= sizeof(void *));
    assert(sizeof(SCOTCH_Num) == sizeof(E_Int));
    assert(SCOTCH_numSizeof() == sizeof(SCOTCH_Num));

    int ret;

    SCOTCH_Dgraph graph;
    SCOTCH_dgraphInit(&graph, MPI_COMM_WORLD);

    if (M->pid == 0) puts("Building graph...");

    ret = SCOTCH_dgraphBuild(
        &graph,
        0,
        M->nc,
        M->nc,
        M->xneis,
        NULL,
        cwgts,
        M->l2gc,
        M->xneis[M->nc],
        M->xneis[M->nc],
        M->cneis,
        NULL,
        NULL
    );

    if (ret != 0) {
        merr("Failed to build Scotch graph.");
        return NULL;
    }

    if (M->pid == 0) puts("Checking Scotch graph...");
    
    ret = SCOTCH_dgraphCheck(&graph);

    if (ret != 0) {
        merr("Failed Scotch graph check.");
        return NULL;
    }

    if (M->pid == 0) puts("Building Scotch strategy...");

    SCOTCH_Strat strat;
    
    ret = SCOTCH_stratInit(&strat);

    if (ret != 0) {
        merr("Failed Scotch strat initialisation.");
        return NULL;
    }

    E_Int *cmap = IntArray(M->nc);

    if (M->pid == 0) puts("Partitioning graph...");

    ret = SCOTCH_dgraphPart(&graph, M->npc, &strat, cmap);
    if (ret != 0) {
        merr("Failed to map dual graph.");
        XFREE(cmap);
    }

    SCOTCH_dgraphExit(&graph);
    SCOTCH_stratExit(&strat);

    return cmap;
}

// TODO(Imad): point redistribution is not robust.

#define TOL 1e-10

int sign(E_Float x)
{
    if (x < -TOL) return -1;
    if (x > TOL) return 1;
    return 0;
}

struct xyz {
    E_Float x, y, z;

    bool operator<(const xyz &p) const {
        if (sign(x-p.x) < 0) return true;
        if (sign(x-p.x) == 0 && sign(y-p.y) < 0) return true;
        if (sign(x-p.x) == 0 && sign(y-p.y) == 0 && sign(z-p.z) < 0) return true;
        return false;
    }
};

static
E_Int Mesh_redistribute(Mesh *M, E_Int *cmap)
{
    /*
    if (M->pid == 0) puts("Verifying mesh integrity...");
    Mesh_set_orientation(M);
    if (M->mode_2D) {
        for (E_Int cid = 0; cid < M->nc; cid++) H18_reorder(cid, M);
    } else {
        for (E_Int cid = 0; cid < M->nc; cid++) H27_reorder(cid, M);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (M->pid == 0) {
        puts("Mesh integrity OK");
        fflush(stdout);
    }
    */

    // Allocate contiguous count arrays

    E_Int count_size = M->npc;

    int *counts = (int *)XCALLOC(19 * count_size, sizeof(int));

    int *cscount = counts;
    int *crcount = counts + 1  * count_size;
    int *scount  = counts + 2  * count_size;
    int *rcount  = counts + 3  * count_size;
    int *sfcount = counts + 4  * count_size;
    int *rfcount = counts + 5  * count_size;
    int *spcount = counts + 6  * count_size;
    int *rpcount = counts + 7  * count_size;
    int *sicount = counts + 8  * count_size;
    int *ricount = counts + 9  * count_size;
    int *sbcount = counts + 10 * count_size;
    int *rbcount = counts + 11 * count_size;
    int *idx     = counts + 12 * count_size;
    int *stagcount  = counts + 13 * count_size;
    int *rtagcount  = counts + 14 * count_size;
    int *snamecount = counts + 15 * count_size;
    int *rnamecount = counts + 16 * count_size;
    int *stypecount = counts + 17 * count_size;
    int *rtypecount = counts + 18 * count_size;

    // Allocate contiguous dists array

    E_Int dist_size = M->npc + 1;

    int *dists = (int *)XCALLOC(18 * dist_size, sizeof(int));

    int *csdist = dists;
    int *crdist = dists + 1  * dist_size;
    int *sdist  = dists + 2  * dist_size;
    int *rdist  = dists + 3  * dist_size;
    int *sfdist = dists + 4  * dist_size;
    int *rfdist = dists + 5  * dist_size;
    int *spdist = dists + 6  * dist_size;
    int *rpdist = dists + 7  * dist_size;
    int *sidist = dists + 8  * dist_size;
    int *sbdist = dists + 9  * dist_size;
    int *ridist = dists + 10 * dist_size;
    int *rbdist = dists + 11 * dist_size;
    int *stagdist  = dists + 12 * dist_size;
    int *rtagdist  = dists + 13 * dist_size;
    int *snamedist = dists + 14 * dist_size;
    int *rnamedist = dists + 15 * dist_size;
    int *stypedist = dists + 16 * dist_size;
    int *rtypedist = dists + 17 * dist_size;


    // Redistribute cells

    if (M->pid == 0) puts("Distributing cells...");

    for (E_Int i = 0; i < M->nc; i++) cscount[cmap[i]]++;

    MPI_Alltoall(cscount, 1, MPI_INT, crcount, 1, MPI_INT, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        csdist[i+1] = csdist[i] + cscount[i];
        crdist[i+1] = crdist[i] + crcount[i];
    }

    E_Int nc = crdist[M->npc];

    // Early return if bad cell map
    // TODO(Imad): could be made better...

    if (nc == 0) {
        merr("WARNING: Proc " SF_D_ " has 0 cells after mesh redistribution. "
             "This case is currently not handled.\n\n", M->pid);
    }

    E_Int min_nc = 0;
    MPI_Allreduce(&nc, &min_nc, 1, XMPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (min_nc == 0) {
        XFREE(counts);
        XFREE(dists);
        return 1;
    }

    // Send cell ids and strides

    E_Int *scids = IntArray(csdist[M->npc]);
    E_Int *rcids = IntArray(crdist[M->npc]);

    E_Int *scstride = IntArray(csdist[M->npc]);
    E_Int *rcstride = IntArray(crdist[M->npc]);
 
    for (E_Int i = 0; i < M->npc; i++) idx[i] = csdist[i];

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int dest = cmap[i];
        scids[idx[dest]] = M->l2gc[i];
        scstride[idx[dest]] = M->cstride[i];
        idx[dest]++;
    }

    MPI_Alltoallv(scids, cscount, csdist, XMPI_INT,
                  rcids, crcount, crdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(scstride, cscount, csdist, XMPI_INT,
                  rcstride, crcount, crdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    std::map<E_Int, E_Int> g2lc;

    for (E_Int i = 0; i < nc; i++) g2lc[rcids[i]] = i;

    // Send cell ranges

    for (E_Int i = 0; i < M->npc; i++) {
        scount[i] = 6 * cscount[i];
        rcount[i] = 6 * crcount[i];
        sdist[i+1] = sdist[i] + scount[i];
        rdist[i+1] = rdist[i] + rcount[i];
    }

    E_Int *scrange = IntArray(sdist[M->npc]);
    E_Int *rcrange = IntArray(rdist[M->npc]);

    for (E_Int i = 0; i < M->npc; i++) idx[i] = sdist[i];

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int where = cmap[i];
        E_Int *crange = Mesh_get_crange(M, i);
        for (E_Int j = 0; j < 6; j++) {
            scrange[idx[where]++] = crange[j];
        }
    }

    MPI_Alltoallv(scrange, scount, sdist, XMPI_INT,
                  rcrange, rcount, rdist, XMPI_INT,
                  MPI_COMM_WORLD);

    // Send cell faces

    E_Int gnf;
    MPI_Allreduce(&M->nf, &gnf, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        scount[i] = 24 * cscount[i];
        rcount[i] = 24 * crcount[i];
        sdist[i+1] = sdist[i] + scount[i];
        rdist[i+1] = rdist[i] + rcount[i];
    }

    E_Int *scells = IntArray(sdist[M->npc]);
    E_Int *rcells = IntArray(rdist[M->npc]);

    for (E_Int i = 0; i < M->npc; i++) idx[i] = sdist[i];

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int where = cmap[i];
        E_Int *cell = Mesh_get_cell(M, i);
        for (E_Int j = 0; j < 24; j++) {
            E_Int lfid = cell[j];

            scells[idx[where]] = (lfid == -1) ? -1 : M->l2gf[lfid];

            idx[where]++;
        }
    }

    MPI_Alltoallv(scells, scount, sdist, XMPI_INT,
                  rcells, rcount, rdist, XMPI_INT,
                  MPI_COMM_WORLD);

    // Faces

    if (M->pid == 0) puts("Distributing faces...");

    E_Int nf = 0;

    std::map<E_Int, E_Int> g2lf;
    g2lf[-1] = -1; // Trick to eliminate an if statement inside inner loop

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pc = &rcells[rdist[i]];

        for (E_Int j = 0; j < rcount[i]; j++) {
            E_Int gfid = pc[j];
            //if (gfid == -1) continue;

            if (g2lf.find(gfid) == g2lf.end()) {
                g2lf[gfid] = nf++;
                rfcount[i]++;
            }
        }
    }

    MPI_Alltoall(rfcount, 1, MPI_INT, sfcount, 1, MPI_INT, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        sfdist[i+1] = sfdist[i] + sfcount[i];
        rfdist[i+1] = rfdist[i] + rfcount[i];
    }

    E_Int *sfids = IntArray(sfdist[M->npc]);
    E_Int *rfids = IntArray(rfdist[M->npc]);

    nf = 0;
    g2lf.clear();
    g2lf[-1] = -1; // Same trick

    E_Int *ptr = rfids;

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pc = &rcells[rdist[i]];

        for (E_Int j = 0; j < rcount[i]; j++) {
            E_Int gfid = pc[j];
            //if (gfid == -1) continue;

            if (g2lf.find(gfid) == g2lf.end()) {
                g2lf[gfid] = nf++;
                *ptr++ = gfid;
            }
        }
    }

    assert(ptr - rfids == rfdist[M->npc]);

    MPI_Alltoallv(rfids, rfcount, rfdist, XMPI_INT,
                  sfids, sfcount, sfdist, XMPI_INT,
                  MPI_COMM_WORLD);

    // Send face strides and ranges 
    E_Int *sfstride = IntArray(sfdist[M->npc]);
    E_Int *rfstride = IntArray(rfdist[M->npc]);

    for (E_Int i = 0; i < M->npc; i++) {
        scount[i] = 4*sfcount[i];
        rcount[i] = 4*rfcount[i];
        sdist[i+1] = sdist[i] + scount[i];
        rdist[i+1] = rdist[i] + rcount[i];
    }

    E_Int *sfrange = IntArray(sdist[M->npc]);
    E_Int *rfrange = IntArray(rdist[M->npc]);

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &sfids[sfdist[i]];

        E_Int *pr = &sfrange[sdist[i]];
        E_Int *ps = &sfstride[sfdist[i]];

        for (E_Int j = 0; j < sfcount[i]; j++) {
            E_Int gfid = pf[j];
            E_Int lfid = M->g2lf.at(gfid);

            E_Int *frange = Mesh_get_frange(M, lfid);
            E_Int fstride = M->fstride[lfid];

            for (E_Int k = 0; k < 4; k++) *pr++ = frange[k];
            *ps++ = fstride;
        }
    }

    MPI_Alltoallv(sfstride, sfcount, sfdist, XMPI_INT,
                  rfstride, rfcount, rfdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(sfrange, scount, sdist, XMPI_INT,
                  rfrange, rcount, rdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    // Send faces and count points to be exchanged

    for (E_Int i = 0; i < M->npc; i++) {
        scount[i] = 8*sfcount[i];
        rcount[i] = 8*rfcount[i];
        sdist[i+1] = sdist[i] + scount[i];
        rdist[i+1] = rdist[i] + rcount[i];
    }

    E_Int *sfaces = IntArray(sdist[M->npc]);
    E_Int *rfaces = IntArray(rdist[M->npc]);

    std::map<E_Int, E_Int> pts;

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &sfids[sfdist[i]];

        E_Int *ptr = &sfaces[sdist[i]];

        // Send the points numbered from 0 to np-1
        // TODO(Imad): maybe cache this?
        pts.clear();
        pts[-1] = -1;
        E_Int count = 0;

        for (E_Int j = 0; j < sfcount[i]; j++) {
            E_Int gfid = pf[j];
            E_Int lfid = M->g2lf.at(gfid);

            E_Int *face = Mesh_get_face(M, lfid);

            for (E_Int k = 0; k < 8; k++) {
                E_Int lp = face[k];

                auto it = pts.find(lp);

                if (it == pts.end()) {
                    *ptr++ = count;
                    pts[lp] = count;
                    count++;
                } else {
                    *ptr++ = it->second;
                }
            }
        }

        spcount[i] = count;
    }
  
    MPI_Alltoallv(sfaces, scount, sdist, XMPI_INT,
                  rfaces, rcount, rdist, XMPI_INT,
                  MPI_COMM_WORLD);

    if (M->pid == 0) puts("Distributing points...");

    MPI_Alltoall(spcount, 1, MPI_INT, rpcount, 1, MPI_INT, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        spdist[i+1] = spdist[i] + spcount[i];
        rpdist[i+1] = rpdist[i] + rpcount[i];
    }

    E_Float *sx = FloatArray(spdist[M->npc]);
    E_Float *sy = FloatArray(spdist[M->npc]);
    E_Float *sz = FloatArray(spdist[M->npc]);
    E_Float *rx = FloatArray(rpdist[M->npc]);
    E_Float *ry = FloatArray(rpdist[M->npc]);
    E_Float *rz = FloatArray(rpdist[M->npc]);

    E_Float *px = sx;
    E_Float *py = sy;
    E_Float *pz = sz;

    for (E_Int i = 0; i < M->npc; i++) 
    {
        E_Int *pf = &sfids[sfdist[i]];

        //E_Int *ptr = &sfaces[sdist[i]];

        pts.clear();
        pts[-1] = -1;
        E_Int count = 0;

        for (E_Int j = 0; j < sfcount[i]; j++) 
        {
            E_Int gfid = pf[j];
            E_Int lfid = M->g2lf.at(gfid);

            E_Int *face = Mesh_get_face(M, lfid);

            for (E_Int k = 0; k < 8; k++) 
            {
                E_Int lp = face[k];

                auto it = pts.find(lp);

                if (it == pts.end()) 
                {
                    *px++ = M->X[lp];
                    *py++ = M->Y[lp];
                    *pz++ = M->Z[lp];
                    pts[lp] = count++;
                }
            }
        }

        assert(spcount[i] == count);
    }
    
    MPI_Alltoallv(sx, spcount, spdist, MPI_DOUBLE,
                  rx, rpcount, rpdist, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    MPI_Alltoallv(sy, spcount, spdist, MPI_DOUBLE,
                  ry, rpcount, rpdist, MPI_DOUBLE,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(sz, spcount, spdist, MPI_DOUBLE,
                  rz, rpcount, rpdist, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // Shift the point indices

    E_Int accum = 0;
    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pn = &rfaces[rdist[i]];

        for (E_Int j = 0; j < rcount[i]; j++) {
            if (pn[j] != -1) {
                pn[j] += accum;
            }
        }

        accum += rpcount[i];
    }

    if (M->pid == 0) puts("Eliminating duplicate points...");

    std::map<E_Int, E_Int> dups;
    E_Int np = 0;
    pts.clear();

    {
        std::map<xyz, E_Int> cloud;
        for (E_Int i = 0; i < rpdist[M->npc]; i++) {
            xyz crd = {rx[i], ry[i], rz[i]};
            auto it = cloud.find(crd);
            if (it == cloud.end()) {
                cloud[crd] = i;
                pts[i] = np;
                np++;
            } else {
                dups[i] = it->second;
            }
        }
    }

    if (M->pid == 0) puts("Renumbering points...");

    E_Float *X = FloatArray(np);
    E_Float *Y = FloatArray(np);
    E_Float *Z = FloatArray(np);

    pts[-1] = -1;
    dups[-1] = -1;

    for (E_Int i = 0; i < rdist[M->npc]; i++) {
        E_Int lp = rfaces[i];

        auto it = dups.find(lp);

        if (it == dups.end()) {

            E_Int new_pid = pts.at(lp);
            rfaces[i] = new_pid;
            X[new_pid] = rx[lp];
            Y[new_pid] = ry[lp];
            Z[new_pid] = rz[lp];
            
        } else {
            
            rfaces[i] = pts.at(it->second);
        }
    }

    // Init temp owner and neigh

    E_Int *owner = IntArray(nf);
    E_Int *neigh = IntArray(nf);
    for (E_Int i = 0; i < nf; i++) {
        owner[i] = -1;
        neigh[i] = -1;
    }

    for (E_Int i = 0; i < nc; i++) {
        E_Int *cell = &rcells[24*i];
        E_Int *crange = &rcrange[6*i];
        E_Int cstride = rcstride[i];

        for (E_Int j = 0; j < cstride; j++) {
            E_Int *pf = cell + 4*j;

            for (E_Int k = 0; k < crange[j]; k++) {
                E_Int gfid = pf[k];
                E_Int lfid = g2lf.at(gfid);
                assert(lfid >= -1 && lfid < nf);
                if (owner[lfid] == -1) owner[lfid] = i;
                else neigh[lfid] = i;
            }
        }
    }

    // Do comm patches

    if (M->pid == 0) puts("Creating comm patches...");

    memset(rcount, 0, M->npc * sizeof(int));

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &rfids[rfdist[i]];

        for (E_Int j = 0; j < rfcount[i]; j++) {
            E_Int gfid = pf[j];
            E_Int lfid = g2lf.at(gfid);

            if (neigh[lfid] != -1) continue;

            rcount[i]++;
        }
    }

    MPI_Alltoall(rcount, 1, MPI_INT, scount, 1, MPI_INT, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        rdist[i+1] = rdist[i] + rcount[i];
        sdist[i+1] = sdist[i] + scount[i];
    }

    // Request comm info

    E_Int *rbinfo = IntArray(rdist[M->npc]);
    E_Int *sbinfo = IntArray(sdist[M->npc]);

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &rfids[rfdist[i]];
        E_Int *pr = &rbinfo[rdist[i]];

        for (E_Int j = 0; j < rfcount[i]; j++) {
            E_Int gfid = pf[j];
            E_Int lfid = g2lf.at(gfid);

            if (neigh[lfid] != -1) continue;

            *pr++ = gfid;
        }
    }

    MPI_Alltoallv(rbinfo, rcount, rdist, XMPI_INT,
                  sbinfo, scount, sdist, XMPI_INT,
                  MPI_COMM_WORLD);


    memset(spcount, 0, M->npc * sizeof(int));
    memset(rpcount, 0, M->npc * sizeof(int));

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &sbinfo[sdist[i]];

        std::set<E_Int> bctag;

        for (E_Int j = 0; j < scount[i]; j++) {
            E_Int gfid = pf[j];
            assert(M->g2lf.find(gfid) != M->g2lf.end());

            E_Int lfid = M->g2lf[gfid];

            if (Mesh_face_is_iface(M, lfid)) {
                // Send cmap[lnei] + gfid + gnei
                sicount[i] += 3;
            }

            else if (Mesh_face_is_pface(M, lfid)) {
                E_Int patch_id = M->face_to_ppatch[lfid];
                E_Int proc = M->pps[patch_id].nei;
                assert(proc >= 0 && proc < M->npc);

                // We will be sending gfid and the requesting proc         
                rpcount[proc] += 2;
            }
            
            else {
                assert(Mesh_face_is_bface(M, lfid));
                
                // Send bc gid + gfid
                sbcount[i] += 2;

                // Send bcname to proc i
                bctag.insert(M->face_to_bpatch.at(lfid));
            }
        }

        for (E_Int tag : bctag) {
            stagcount[i] += 1;
            char *bcname = M->bps[tag].name;
            snamecount[i] += strlen(bcname) + 1;
            char *bctype = M->bps[tag].type;
            stypecount[i] += strlen(bctype) + 1;
        }
    }

    MPI_Alltoall(sicount, 1, MPI_INT, ricount, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(rpcount, 1, MPI_INT, spcount, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(sbcount, 1, MPI_INT, rbcount, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Alltoall(stagcount,  1, MPI_INT, rtagcount,  1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(snamecount, 1, MPI_INT, rnamecount, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(stypecount, 1, MPI_INT, rtypecount, 1, MPI_INT, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        sidist[i+1] = sidist[i] + sicount[i];
        spdist[i+1] = spdist[i] + spcount[i];
        sbdist[i+1] = sbdist[i] + sbcount[i];
        ridist[i+1] = ridist[i] + ricount[i];
        rpdist[i+1] = rpdist[i] + rpcount[i];
        rbdist[i+1] = rbdist[i] + rbcount[i];
        
        stagdist[i+1] = stagdist[i] + stagcount[i];
        rtagdist[i+1] = rtagdist[i] + rtagcount[i];
        
        snamedist[i+1] = snamedist[i] + snamecount[i];
        rnamedist[i+1] = rnamedist[i] + rnamecount[i];

        stypedist[i+1] = stypedist[i] + stypecount[i];
        rtypedist[i+1] = rtypedist[i] + rtypecount[i];
    }

    E_Int *sidata = IntArray(sidist[M->npc]);
    E_Int *spdata = IntArray(spdist[M->npc]);
    E_Int *sbdata = IntArray(sbdist[M->npc]);
    E_Int *ridata = IntArray(ridist[M->npc]);
    E_Int *rpdata = IntArray(rpdist[M->npc]);
    E_Int *rbdata = IntArray(rbdist[M->npc]);

    E_Int *stag = IntArray(stagdist[M->npc]);
    E_Int *rtag = IntArray(rtagdist[M->npc]);
    char *sname = (char *)XMALLOC(snamedist[M->npc] * sizeof(char));
    char *rname = (char *)XMALLOC(rnamedist[M->npc] * sizeof(char));
    char *stype = (char *)XMALLOC(stypedist[M->npc] * sizeof(char));
    char *rtype = (char *)XMALLOC(rtypedist[M->npc] * sizeof(char));

    for (E_Int i = 0; i < M->npc; i++) idx[i] = rpdist[i];

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &sbinfo[sdist[i]];
        
        E_Int *iptr = &sidata[sidist[i]];
        E_Int *bptr = &sbdata[sbdist[i]];

        std::set<E_Int> bctag;
        
        for (E_Int j = 0; j < scount[i]; j++) {
            E_Int gfid = pf[j];
            assert(M->g2lf.find(gfid) != M->g2lf.end());

            E_Int lfid = M->g2lf[gfid];

            if (Mesh_face_is_iface(M, lfid)) {
                E_Int own = M->owner[lfid];
                E_Int nei = M->neigh[lfid];

                E_Int send_id = cmap[own] == i ? nei : own;
                
                assert(cmap[send_id] != i);

                *iptr++ = cmap[send_id];
                *iptr++ = gfid;
                *iptr++ = M->l2gc[send_id];
            }
            
            else if (Mesh_face_is_pface(M, lfid)) {
                E_Int patch_id = M->face_to_ppatch[lfid];
                E_Int proc = M->pps[patch_id].nei;
                assert(proc >= 0 && proc < M->npc);    

                rpdata[idx[proc]++] = gfid;
                rpdata[idx[proc]++] = i;
            }
            
            else {
                E_Int lid = M->face_to_bpatch[lfid];
                E_Int gid = M->bps[lid].gid;
                
                *bptr++ = gid;
                *bptr++ = gfid;

                // Send bcname to proc i
                bctag.insert(M->face_to_bpatch.at(lfid));
            }
        }

        E_Int *tptr = &stag[stagdist[i]];
        char *nptr = &sname[snamedist[i]];
        char *Tptr = &stype[stypedist[i]];

        for (E_Int tag : bctag) {
            *tptr++ = M->bps[tag].gid;

            char *bcname = M->bps[tag].name;
            for (size_t j = 0; j < strlen(bcname); j++) *nptr++ = bcname[j];
            *nptr++ = '\0';

            char *bctype = M->bps[tag].type;
            for (size_t j = 0; j < strlen(bctype); j++) *Tptr++ = bctype[j];
            *Tptr++ = '\0';
        }
    }

    for (E_Int i = 0; i < M->npc; i++) assert(idx[i] == rpdist[i+1]);

    MPI_Alltoallv(sidata, sicount, sidist, XMPI_INT,
                  ridata, ricount, ridist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(rpdata, rpcount, rpdist, XMPI_INT,
                  spdata, spcount, spdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(sbdata, sbcount, sbdist, XMPI_INT,
                  rbdata, rbcount, rbdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(stag, stagcount, stagdist, XMPI_INT,
                  rtag, rtagcount, rtagdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(sname, snamecount, snamedist, MPI_CHAR,
                  rname, rnamecount, rnamedist, MPI_CHAR,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(stype, stypecount, stypedist, MPI_CHAR,
                  rtype, rtypecount, rtypedist, MPI_CHAR,
                  MPI_COMM_WORLD);

    // Parse boundary patches

    if (M->pid == 0) puts("Creating boundary patches...");

    E_Int nbp = 0;

    std::map<E_Int, E_Int> gbc_to_lbc;
    std::map<E_Int, E_Int> lbc_to_gbc;
    std::map<E_Int, E_Int> lbc_to_size;

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &rbdata[rbdist[i]];

        for (E_Int j = 0; j < rbcount[i]; ) {
            E_Int gid = pf[j++];

            // skip gfid
            j++;
            //E_Int gfid = pf[j++];
            //assert(g2lf.find(gfid) != g2lf.end());

            auto it = gbc_to_lbc.find(gid);

            if (it == gbc_to_lbc.end()) {
                gbc_to_lbc[gid] = nbp;
                lbc_to_gbc[nbp] = gid;
                lbc_to_size[nbp]++;
                nbp++;
            } else {
                lbc_to_size[it->second]++;
            }
        }
    }

    BPatch *bps = (BPatch *)XCALLOC(nbp, sizeof(BPatch));

    for (E_Int i = 0; i < nbp; i++) {
        BPatch *P = &bps[i];
        P->gid = lbc_to_gbc[i];
        P->nf = lbc_to_size[i];
        P->pf = IntArray(P->nf);
        P->nf = 0; // Reset for next loop.
        P->name = NULL;
        P->type = NULL;
    }

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *tptr = &rtag[rtagdist[i]];
        char *nptr = &rname[rnamedist[i]];
        char *Tptr = &rtype[rtypedist[i]];

        for (E_Int j = 0; j < rtagcount[i]; j++) {
            E_Int tag = tptr[j];

            // bcname
            char bcname[1024];
            char c;
            E_Int k = 0;
            while ((c = *nptr++)) {
                bcname[k++] = c;
            }
            bcname[k] = '\0';

            // bctype
            char bctype[1024];
            k = 0;
            while ((c = *Tptr++)) {
                bctype[k++] = c;
            }
            bctype[k] = '\0';

            E_Int ltag = gbc_to_lbc.at(tag);
            BPatch *P = &bps[ltag];
            if (P->name) continue;

            P->name = (char *)XMALLOC(strlen(bcname)+1);
            strcpy(P->name, bcname);
            assert(P->name[strlen(P->name)] == '\0');

            assert(P->type == NULL);

            P->type = (char *)XMALLOC(strlen(bctype)+1);
            strcpy(P->type, bctype);
            assert(P->type[strlen(P->type)] == '\0');
        }
    }

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &rbdata[rbdist[i]];

        for (E_Int j = 0; j < rbcount[i]; ) {
            E_Int gbc = pf[j++];
            E_Int gfid = pf[j++];

            E_Int lbc = gbc_to_lbc[gbc];

            BPatch *P = &bps[lbc];

            P->pf[P->nf] = g2lf[gfid];

            P->nf++;
        }
    }

    // Warning(Imad): reusing scount, rcount, sdist and rdist
    memset(scount, 0, M->npc * sizeof(int));

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &spdata[spdist[i]];

        for (E_Int j = 0; j < spcount[i]; ) {
            //E_Int gfid = pf[j++];
            j++; // skip gfid
            //if (M->g2lf.find(gfid) == M->g2lf.end()) abort();
            E_Int proc = pf[j++];

            // We will be sending cmap[gown], gfid and gown
            scount[proc] += 3;
        }
    }

    MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

    for (E_Int i = 0; i < M->npc; i++) {
        sdist[i+1] = sdist[i] + scount[i];
        rdist[i+1] = rdist[i] + rcount[i];
    }

    E_Int *spinfo = IntArray(sdist[M->npc]);
    E_Int *rpinfo = IntArray(rdist[M->npc]);

    for (E_Int i = 0; i < M->npc; i++) idx[i] = sdist[i];

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &spdata[spdist[i]];

        for (E_Int j = 0; j < spcount[i]; ) {
            E_Int gfid = pf[j++];
            E_Int proc = pf[j++];

            E_Int lfid = M->g2lf[gfid];
            assert(Mesh_face_is_pface(M, lfid));
            E_Int lown = M->owner[lfid];
            E_Int gown = M->l2gc[lown];

            spinfo[idx[proc]++] = cmap[lown];
            spinfo[idx[proc]++] = gfid;
            spinfo[idx[proc]++] = gown;
        }
    }

    for (E_Int i = 0; i < M->npc; i++) assert(idx[i] == sdist[i+1]);

    MPI_Alltoallv(spinfo, scount, sdist, XMPI_INT,
                  rpinfo, rcount, rdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    // Parse comm patches

    std::map<E_Int, E_Int> gproc_to_lproc;
    std::map<E_Int, E_Int> lproc_to_gproc;
    std::map<E_Int, E_Int> lproc_to_size;

    E_Int npp = 0;

    // Former internal faces contributions

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &ridata[ridist[i]];

        for (E_Int j = 0; j < ricount[i]; ) {
            E_Int proc = pf[j++];
            
            //E_Int gfid = pf[j++];
            j++; // skip gfid

            //E_Int gnei = pf[j++];
            j++; // skip gnei

            auto it = gproc_to_lproc.find(proc);

            if (it == gproc_to_lproc.end()) {
                gproc_to_lproc[proc] = npp;
                lproc_to_gproc[npp] = proc;
                lproc_to_size[npp]++;
                npp++;
            } else {
                lproc_to_size[it->second]++;
            }

            //assert(g2lf.find(gfid) != g2lf.end());

            //assert(g2lc.find(gnei) == g2lc.end());
        }
    }
    
    // Former ppatch faces contributions

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &rpinfo[rdist[i]];

        for (E_Int j = 0; j < rcount[i]; ) {
            E_Int proc = pf[j++];
            assert(proc != M->pid);

            auto it = gproc_to_lproc.find(proc);

            if (it == gproc_to_lproc.end()) {
                gproc_to_lproc[proc] = npp;
                lproc_to_gproc[npp] = proc;
                lproc_to_size[npp]++;
                npp++;
            } else {
                lproc_to_size[it->second]++;
            }

            // skip gfid;
            j++;

            // skip gnei;
            j++;

            //E_Int gfid = pf[j++];
            //assert(g2lf.find(gfid) != g2lf.end());

            //E_Int gnei = pf[j++];
            //assert(g2lc.find(gnei) == g2lc.end());
        }
    }

    PPatch *pps = (PPatch *)XCALLOC(npp, sizeof(PPatch));

    for (E_Int i = 0; i < npp; i++) {
        PPatch *P = &pps[i];
        P->nei = lproc_to_gproc[i];
        P->nf = lproc_to_size[i];
        P->pf = IntArray(P->nf);
        P->pn = IntArray(P->nf);
        P->nf = 0; // For the next loop.
    }

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &ridata[ridist[i]];

        for (E_Int j = 0; j < ricount[i]; ) {
            E_Int proc = pf[j++];
            E_Int gfid = pf[j++];
            E_Int gnei = pf[j++];

            PPatch *P = &pps[gproc_to_lproc[proc]];

            //P->pf[P->nf] = g2lf[gfid];
            P->pf[P->nf] = gfid;
            P->pn[P->nf] = gnei;

            P->nf++;
        }
    }

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &rpinfo[rdist[i]];

        for (E_Int j = 0; j < rcount[i]; ) {
            E_Int proc = pf[j++];
            E_Int gfid = pf[j++];
            E_Int gnei = pf[j++];

            PPatch *P = &pps[gproc_to_lproc[proc]];

            //P->pf[P->nf] = g2lf[gfid];
            P->pf[P->nf] = gfid;
            P->pn[P->nf] = gnei;

            P->nf++;
        }
    }

    // Sort pfaces/pneis
    
    for (E_Int i = 0; i < npp; i++) {
        PPatch *P = &pps[i];

        std::vector<E_Int> indices(P->nf);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
            [&] (E_Int i, E_Int j) { return P->pf[i] < P->pf[j]; });
        
        E_Int *pf = IntArray(P->nf);
        E_Int *pn = IntArray(P->nf);

        for (E_Int j = 0; j < P->nf; j++) {
            // Replace with local indices
            pf[j] = g2lf.at(P->pf[indices[j]]);
            pn[j] = P->pn[indices[j]];
        }

        XFREE(P->pf);
        XFREE(P->pn);

        P->pf = pf;
        P->pn = pn;
    }

    // Redistribute cref/clevel/ctype

    if (M->pid == 0) puts("Distributing cell adaptation data...");

    E_Int *scref = IntArray(csdist[M->npc]);
    E_Int *sclvl = IntArray(csdist[M->npc]);
    E_Int *sctyp = IntArray(csdist[M->npc]);
    
    E_Int *rcref = IntArray(nc);
    E_Int *rclvl = IntArray(nc);
    E_Int *rctyp = IntArray(nc);

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pc = &scids[csdist[i]];
        E_Int *pr = &scref[csdist[i]];
        E_Int *pl = &sclvl[csdist[i]];
        E_Int *pt = &sctyp[csdist[i]];
        
        for (E_Int j = 0; j < cscount[i]; j++) {
            E_Int gcid = pc[j];
            assert(M->g2lc.find(gcid) != M->g2lc.end());
            
            E_Int lcid = M->g2lc[gcid];
            assert(lcid >= 0 && lcid < M->nc);
            
            E_Int cval = M->cref[lcid];
            E_Int clvl = M->clevel[lcid];
            E_Int ctyp = M->ctype[lcid];
            
            *pr++ = cval;
            *pl++ = clvl;
            *pt++ = ctyp;
        }
    }

    MPI_Alltoallv(scref, cscount, csdist, XMPI_INT,
                  rcref, crcount, crdist, XMPI_INT,
                  MPI_COMM_WORLD);

    MPI_Alltoallv(sclvl, cscount, csdist, XMPI_INT,
                  rclvl, crcount, crdist, XMPI_INT,
                  MPI_COMM_WORLD);

    MPI_Alltoallv(sctyp, cscount, csdist, XMPI_INT,
                  rctyp, crcount, crdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    // Redistribute flevel/ftype/fref

    if (M->pid == 0) puts("Distributing face adaptation data...");
    
    E_Int *sflvl = IntArray(sfdist[M->npc]);
    E_Int *sftyp = IntArray(sfdist[M->npc]);
    E_Int *sfref = IntArray(sfdist[M->npc]);

    E_Int *rflvl = IntArray(nf);
    E_Int *rftyp = IntArray(nf);
    E_Int *rfref = IntArray(nf);

    for (E_Int i = 0; i < M->npc; i++) {
        E_Int *pf = &sfids[sfdist[i]];
        E_Int *pl = &sflvl[sfdist[i]];
        E_Int *pt = &sftyp[sfdist[i]];
        E_Int *pr = &sfref[sfdist[i]];
        
        for (E_Int j = 0; j < sfcount[i]; j++) {
            E_Int gfid = pf[j];
            E_Int lfid = M->g2lf.at(gfid);
            assert(lfid >= 0 && lfid < M->nf);
            
            E_Int flvl = M->flevel[lfid];
            E_Int ftyp = M->ftype[lfid];
            E_Int fref = M->fref[lfid];

            *pl++ = flvl;
            *pt++ = ftyp;
            *pr++ = fref;
        }
    }

    MPI_Alltoallv(sflvl, sfcount, sfdist, XMPI_INT,
                  rflvl, rfcount, rfdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(sftyp, sfcount, sfdist, XMPI_INT,
                  rftyp, rfcount, rfdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    MPI_Alltoallv(sfref, sfcount, sfdist, XMPI_INT,
                  rfref, rfcount, rfdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    // Make the new mesh

    if (M->pid == 0) puts("Creating the new mesh...");

    Mesh_reset_base_data(M);

    Mesh_reset_boundary_data(M);
    
    Mesh_reset_adaptation_data(M);
    
    Mesh_reset_comm_data(M);
    
    Mesh_reset_parallel_data(M);

    M->np = np;
    M->X = X;
    M->Y = Y;
    M->Z = Z;

    M->nf = nf;
    M->faces = rfaces;
    M->fstride = rfstride;
    M->frange = rfrange;

    M->nc = nc;
    M->cells = rcells;
    M->cstride = rcstride;
    M->crange = rcrange;

    M->nbp = nbp;
    M->bps = bps;

    assert(M->face_to_bpatch.empty());

    for (E_Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            assert(lfid >= 0 && lfid < M->nf);
            M->face_to_bpatch[lfid] = i;
        }
    }

    M->npp = npp;
    M->pps = pps;
    
    assert(M->face_to_ppatch.empty());

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            assert(lfid >= 0 && lfid < M->nf);
            M->face_to_ppatch[lfid] = i;
        }
    }

    M->g2lc = g2lc;
    M->g2lf = g2lf;

    M->l2gc = rcids;
    M->l2gf = rfids;

    // Set nface to local indices

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int *cell = Mesh_get_cell(M, i);
        E_Int *crange = Mesh_get_crange(M, i);
        E_Int cstride = M->cstride[i];

        for (E_Int j = 0; j < cstride; j++) {
            E_Int *pf = cell + 4*j;

            for (E_Int k = 0; k < crange[j]; k++) {
                pf[k] = M->g2lf[pf[k]];
                assert(pf[k] < M->nf);
            }
        }
    }

    // Adaptation

    M->cref = rcref;
    M->fref = rfref;

    M->clevel = rclvl;
    M->ctype = rctyp;

    M->flevel = rflvl;
    M->ftype = rftyp;

    /*
    if (M->pid == 0) puts("Verifying mesh integrity after balancing...");
    Mesh_set_orientation(M);
    if (M->mode_2D) {
        for (E_Int cid = 0; cid < M->nc; cid++) H18_reorder(cid, M);
    } else {
        for (E_Int cid = 0; cid < M->nc; cid++) H27_reorder(cid, M);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (M->pid == 0) puts("Mesh integrity OK.");
    */

    // Clean-up

    if (M->pid == 0) puts("Cleaning-up...");

    XFREE(stag);
    XFREE(rtag);
    XFREE(sname);
    XFREE(rname);
    XFREE(stype);
    XFREE(rtype);

    XFREE(rx);
    XFREE(ry);
    XFREE(rz);
    XFREE(sfref);
    XFREE(sftyp);
    XFREE(sflvl);
    XFREE(sctyp);
    XFREE(sclvl);
    XFREE(scref);
    
    XFREE(spinfo);
    XFREE(rpinfo);
    XFREE(rbdata);
    XFREE(rpdata);
    XFREE(ridata);
    XFREE(sbdata);
    XFREE(spdata);
    XFREE(sidata);
    XFREE(sbinfo);
    XFREE(rbinfo);
    XFREE(sx);
    XFREE(sy);
    XFREE(sz);
    
    XFREE(sfids);
    XFREE(sfaces);
    XFREE(sfstride);
    XFREE(sfrange);
    
    XFREE(scids);
    XFREE(scells);
    XFREE(scstride);
    XFREE(scrange);
    
    XFREE(owner);
    XFREE(neigh);

    XFREE(dists);
    XFREE(counts);

    if (M->pid == 0) puts("Done.");

    return 0;
}

E_Int Mesh_load_balance(Mesh *M)
{
    E_Int *cwgts = compute_cell_weights(M);

    E_Int *cmap = map_cell_graph(M, cwgts);

    XFREE(cwgts);

    if (cmap == NULL) return 1;

    E_Int ret = Mesh_redistribute(M, cmap);

    E_Int gret = 0;
    MPI_Allreduce(&ret, &gret, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    XFREE(cmap);

    E_Int gnf = 0;
    MPI_Allreduce(&M->nf, &gnf, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD); 

    return gret;
}
