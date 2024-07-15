#include "Mesh.h"
#include "Karray.h"
#include "../common/mem.h"

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

    clevel = NULL;
    flevel = NULL;
    elevel = NULL;

    ctype = NULL;
    ftype = NULL;

    fparent = NULL;

    nc_old = 0;
    nf_old = 0;

    /* Parallel */

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &npc);

    nrq = 0;
    reqs = (MPI_Request *)XMALLOC(2 * npc * sizeof(MPI_Request));

    npp = 0;
    pps = NULL;

    l2gc = l2gf = NULL;

    xneis = cneis = NULL;
}

Mesh *Mesh_from_Karray(Karray *karray)
{    
    Int np = karray->npoints();
    Float *X = karray->X();
    Float *Y = karray->Y();
    Float *Z = karray->Z();

    Int nf = karray->nfaces();
    //Int *indpg = karray->indpg();
    
    Int nc = karray->ncells();
    //Int *indph = karray->indph();

    Mesh *M = new Mesh;

    M->np = np;
    M->X = FloatArray(np);
    M->Y = FloatArray(np);
    M->Z = FloatArray(np);

    memcpy(M->X, X, np * sizeof(Float));
    memcpy(M->Y, Y, np * sizeof(Float));
    memcpy(M->Z, Z, np * sizeof(Float));

    M->nf = nf;

    M->faces = IntArray(8 * M->nf);
    memset(M->faces, -1, 8 * M->nf * sizeof(Int));
    M->fstride = IntArray(M->nf);
    M->frange = IntArray(4 * M->nf);

    for (Int i = 0; i < M->nf; i++) {
        Int np = -1;
        Int *pn = karray->get_face(i, np);
        assert(np != -1);

        M->fstride[i] = np;
        
        Int *face = Mesh_get_face(M, i);
        Int *frange = Mesh_get_frange(M, i);

        for (Int j = 0; j < np; j++) {
            face[2*j] = pn[j] - 1;
            frange[j] = 1;
        }
    }

    M->nc = nc;

    M->cells = IntArray(24 * M->nc);
    memset(M->cells, -1, 24 * M->nc * sizeof(Int));
    M->cstride = IntArray(M->nc);
    M->crange = IntArray(6 * M->nc);

    for (Int i = 0; i < M->nc; i++) {
        Int nf = -1;
        Int *pf = karray->get_cell(i, nf);
        assert(nf != -1);

        M->cstride[i] = nf;
        
        Int *cell = &M->cells[24*i];
        Int *crange = &M->crange[6*i];

        for (Int j = 0; j < nf; j++) {
            cell[4*j] = pf[j] - 1;
            crange[j] = 1;
        }
    }

    return M;
}
