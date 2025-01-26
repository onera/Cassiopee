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
#include "common/Karray.h"
#include "common/mem.h"

Mesh::Mesh()
{
    /* Base data */

    np = 0;
    X = Y = Z = NULL;

    ne = 0;
    edges = NULL;
    fedg = NULL;
    xedf = NULL;
    E2F = NULL;

    nf = 0;
    faces = NULL;
    fstride = NULL;
    frange = NULL;

    nc = 0;
    cells = NULL;
    cstride = NULL;
    crange = NULL;

    owner = neigh = NULL;

    nbp = 0;
    bps = NULL;

    /* Adaptation */

    mode_2D = NULL;

    cref = NULL;
    fref = NULL;
    fpattern = NULL;

    clevel = NULL;
    flevel = NULL;
    elevel = NULL;

    ctype = NULL;
    ftype = NULL;

    fparent = NULL;

    nc_old = 0;
    nf_old = 0;

    /* Parallel */
    pid = 0;
    npc = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &npc);
    
    nrq = 0;
    reqs = (MPI_Request *)XMALLOC(2 * npc * sizeof(MPI_Request));

    npp = 0;
    pps = NULL;

    l2gc = l2gf = NULL;

    xneis = cneis = NULL;
    
    ctag = NULL;
    ftag = NULL;
    ptag = NULL;
}

Mesh *Mesh_from_Karray(Karray *karray)
{    
    E_Int np = karray->npoints();
    E_Float *X = karray->x;
    E_Float *Y = karray->y;
    E_Float *Z = karray->z;

    E_Int nf = karray->nfaces();
    
    E_Int nc = karray->ncells();

    Mesh *M = new Mesh;

    M->np = np;
    M->X = FloatArray(np);
    M->Y = FloatArray(np);
    M->Z = FloatArray(np);

    memcpy(M->X, X, np * sizeof(E_Float));
    memcpy(M->Y, Y, np * sizeof(E_Float));
    memcpy(M->Z, Z, np * sizeof(E_Float));

    M->nf = nf;

    M->faces = IntArray(8 * M->nf);
    memset(M->faces, -1, 8 * M->nf * sizeof(E_Int));
    M->fstride = IntArray(M->nf);
    M->frange = IntArray(4 * M->nf);

    for (E_Int i = 0; i < M->nf; i++) {
        E_Int np = -1;
        E_Int *pn = karray->get_face(i, np);
        assert(np != -1);

        M->fstride[i] = np;
        
        E_Int *face = Mesh_get_face(M, i);
        E_Int *frange = Mesh_get_frange(M, i);

        for (E_Int j = 0; j < np; j++) {
            face[2*j] = pn[j] - 1;
            frange[j] = 1;
        }
    }

    M->nc = nc;

    M->cells = IntArray(24 * M->nc);
    memset(M->cells, -1, 24 * M->nc * sizeof(E_Int));
    M->cstride = IntArray(M->nc);
    M->crange = IntArray(6 * M->nc);

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int nf = -1;
        E_Int *pf = karray->get_cell(i, nf);
        assert(nf != -1);

        M->cstride[i] = nf;
        
        E_Int *cell = &M->cells[24*i];
        E_Int *crange = &M->crange[6*i];

        for (E_Int j = 0; j < nf; j++) {
            cell[4*j] = pf[j] - 1;
            crange[j] = 1;
        }
    }

    return M;
}
