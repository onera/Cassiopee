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

static
void PPatch_free(PPatch *pp)
{
    XFREE(pp->pf);
    XFREE(pp->pn);
    XFREE(pp->sbuf_i);
    XFREE(pp->rbuf_i);
}

static
void BPatch_free(BPatch *bp)
{
    XFREE(bp->pf);
    XFREE(bp->name);
    XFREE(bp->type);
}

void Mesh_reset_base_data(Mesh *M)
{
    M->np = 0;
    XFREE(M->X);
    XFREE(M->Y);
    XFREE(M->Z);

    M->ne = 0;
    XFREE(M->edges);
    XFREE(M->fedg);
    XFREE(M->xedf);
    XFREE(M->E2F);

    M->nf = 0;
    XFREE(M->faces);
    XFREE(M->frange);
    XFREE(M->fstride);

    M->nc = 0;
    XFREE(M->cells);
    XFREE(M->crange);
    XFREE(M->cstride);
    
    XFREE(M->owner);
    XFREE(M->neigh);
}

void Mesh_reset_boundary_data(Mesh *M)
{
    for (E_Int i = 0; i < M->nbp; i++) BPatch_free(&M->bps[i]);
    M->nbp = 0;
    XFREE(M->bps);

    M->face_to_bpatch.clear();
}

void Mesh_reset_adaptation_data(Mesh *M)
{
    M->ecenter.clear();

    XFREE(M->cref);
    XFREE(M->fref);
    XFREE(M->fpattern);

    XFREE(M->clevel);
    XFREE(M->flevel);
    XFREE(M->elevel);

    XFREE(M->ctype);
    XFREE(M->ftype);

    M->fchildren.clear();
    M->cchildren.clear();

    XFREE(M->fparent);
}

void Mesh_reset_comm_data(Mesh *M)
{
    for (E_Int i = 0; i < M->npp; i++) PPatch_free(&M->pps[i]);
    M->npp = 0;
    XFREE(M->pps);

    M->face_to_ppatch.clear();
}

void Mesh_reset_parallel_data(Mesh *M)
{
    XFREE(M->l2gc);
    XFREE(M->l2gf);

    M->g2lc.clear();
    M->g2lf.clear();

    XFREE(M->xneis);
    XFREE(M->cneis);
}

void Mesh_reset_tags(Mesh *M)
{
    XFREE(M->ctag);
    XFREE(M->ftag);
    XFREE(M->ptag);
}

void Mesh_free(Mesh *M)
{
    Mesh_reset_base_data(M);
    Mesh_reset_boundary_data(M);
    Mesh_reset_adaptation_data(M);
    Mesh_reset_comm_data(M);
    Mesh_reset_parallel_data(M);
    Mesh_reset_tags(M);
    delete M;
}
