#include "Proto.h"
#include "scotch/ptscotch.h"

void Comm_waitall(AMesh *M)
{
  MPI_Waitall(M->nrq, M->req, MPI_STATUSES_IGNORE);
  M->nrq = 0;
}

static
E_Int *make_dual_graph(AMesh *M, E_Int *XNEIS)
{
  // Internal faces
  E_Int nif = 0;
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;
    E_Int own = M->owner[i];
    XNEIS[own+1] += 1;
    XNEIS[nei+1] += 1;
    nif += 1;
  }
  assert(M->nif == nif);

  // How many proc faces do we have
  E_Int npf = 0;
  
  for (E_Int i = 0; i < M->npatches; i++) {
    npf += M->patches[i].nf;

    for (E_Int j = 0; j < M->patches[i].nf; j++) {
      E_Int face = M->patches[i].pf[j];
      XNEIS[M->owner[face]+1] += 1;
    }
  }

  for (E_Int i = 0; i < M->ncells; i++) XNEIS[i+1] += XNEIS[i];

  E_Int *CADJ = (E_Int *)XMALLOC((2*nif + npf) * sizeof(E_Int));

  // Second pass: fill global adjacency array
  E_Int *cneis = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  // Internal faces first
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;
    E_Int own = M->owner[i];

    CADJ[XNEIS[own] + cneis[own]++] = M->gcells[nei];
    CADJ[XNEIS[nei] + cneis[nei]++] = M->gcells[own];
  }

  // Proc faces
  for (E_Int i = 0; i < M->npatches; i++) {
    for (E_Int j = 0; j < M->patches[i].nf; j++) {
      E_Int face = M->patches[i].pf[j];
      E_Int own = M->owner[face];
      CADJ[XNEIS[own] + cneis[own]++] = M->patches[i].pn[j];
    }
  }

  XFREE(cneis);

  return CADJ;
}

static
E_Int *map_graph(AMesh *M, E_Int *XNEIS, E_Int *CADJ, E_Int *CWGT)
{
  SCOTCH_Dgraph graph;
  SCOTCH_dgraphInit(&graph, MPI_COMM_WORLD);

  E_Int ret;
  ret = SCOTCH_dgraphBuild(
    &graph,
    0,                // 0-based
    M->ncells,        // number of local cells
    M->ncells,        // compact graph
    XNEIS,            // local adjacency array
    NULL,             // compact graph
    CWGT,             // local cell weights
    M->gcells,        // local cell labels
    XNEIS[M->ncells], // local number of arcs (twice the number of edges)
    XNEIS[M->ncells], // compact graph
    CADJ,             // local adjacency array
    NULL,
    NULL
  );

  if (ret != 0) {
    fprintf(stderr, "AdaptMesh: Error in remapping dual graph.\n");
    assert(0);
    exit(1);
  }

  ret = SCOTCH_dgraphCheck(&graph);
  if (ret != 0) {
    fprintf(stderr, "AdaptMesh: Error in checking dual graph.\n");
    assert(0);
    exit(1);
  }

  SCOTCH_Strat strat;
  ret = SCOTCH_stratInit(&strat);
  if (ret != 0) {
    fprintf(stderr, "AdaptMesh: Failed to init strat\n");
    assert(0);
    exit(1);
  }

  E_Int *map = (E_Int *)XMALLOC(M->ncells * sizeof(E_Int));

  ret = SCOTCH_dgraphPart(&graph, M->npc, &strat, map);

  if (ret != 0) {
    fprintf(stderr, "AdaptMesh: Failed to map dual graph\n");
    assert(0);
    exit(1);
  }

  return map;
}

static
AMesh *reconstruct_mesh(AMesh *M, E_Int *cmap)
{
  AMesh *m = new AMesh;
  m->ecenter = new std::map<Edge, E_Int>;
  m->CT = new std::unordered_map<E_Int, E_Int>;
  m->FT = new std::unordered_map<E_Int, E_Int>;
  m->PT = new std::unordered_map<E_Int, E_Int>;

  E_Int npc = M->npc;

  auto &mCT = *(m->CT);
  auto &mFT = *(m->FT);
  auto &mPT = *(m->PT);

  auto &MCT = *(M->CT);
  auto &MFT = *(M->FT);
  auto &MPT = *(M->PT);

  /* CELLS */
  int *c_scount = (E_Int *)XCALLOC(npc, sizeof(int));
  int *c_rcount = (E_Int *)XMALLOC(npc * sizeof(int));

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = cmap[i];
    c_scount[where] += 1;
  }

  MPI_Alltoall(c_scount, 1, MPI_INT, c_rcount, 1, MPI_INT, MPI_COMM_WORLD);

  int *c_sdist = (E_Int *)XMALLOC((npc+1) * sizeof(int));
  int *c_rdist = (E_Int *)XMALLOC((npc+1) * sizeof(int));
  c_sdist[0] = c_rdist[0] = 0;
  for (E_Int i = 0; i < npc; i++) {
    c_sdist[i+1] = c_sdist[i] + c_scount[i];
    c_rdist[i+1] = c_rdist[i] + c_rcount[i];
  }

  m->ncells = c_rdist[npc];
  m->gcells = (E_Int *)XMALLOC(m->ncells * sizeof(E_Int));

  E_Int *scells = (E_Int *)XMALLOC(c_sdist[npc] * sizeof(E_Int));
  E_Int *c_stride = (E_Int *)XMALLOC(c_sdist[npc] * sizeof(E_Int));

  int *idx = (E_Int *)XMALLOC(npc * sizeof(int));
  for (E_Int i = 0; i < npc; i++) idx[i] = c_sdist[i];

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = cmap[i];
    scells[idx[where]] = M->gcells[i];
    c_stride[idx[where]] = M->indPH[i+1] - M->indPH[i];
    idx[where] += 1;
  }

  // Global cell indices
  MPI_Alltoallv(scells, c_scount, c_sdist, XMPI_INT,
                m->gcells, c_rcount, c_rdist, XMPI_INT,
                MPI_COMM_WORLD);
  
  for (E_Int i = 0; i < m->ncells; i++)
    mCT[m->gcells[i]] = i;
  
  m->indPH = (E_Int *)XMALLOC((m->ncells+1) * sizeof(E_Int));

  // indPH
  MPI_Alltoallv(c_stride, c_scount, c_sdist, XMPI_INT,
                m->indPH+1, c_rcount, c_rdist, XMPI_INT,
                MPI_COMM_WORLD);


  m->indPH[0] = 0;
  for (E_Int i = 0; i < m->ncells; i++) m->indPH[i+1] += m->indPH[i];

  // nface
  int *scount = (E_Int *)XCALLOC(npc, sizeof(int));
  int *rcount = (E_Int *)XMALLOC(npc * sizeof(int));

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = cmap[i];
    scount[where] += M->indPH[i+1] - M->indPH[i];
  }

  MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

  E_Int *sdist = (E_Int *)XMALLOC((npc+1) * sizeof(int));
  E_Int *rdist = (E_Int *)XMALLOC((npc+1) * sizeof(int));
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < npc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  E_Int *sdata = (E_Int *)XMALLOC(sdist[npc] * sizeof(E_Int));
  m->nface = (E_Int *)XCALLOC(rdist[npc], sizeof(E_Int));

  assert(rdist[npc] == m->indPH[m->ncells]);

  for (E_Int i = 0; i < npc; i++) idx[i] = sdist[i];

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = cmap[i];
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++)
      sdata[idx[where]++] = M->gfaces[M->nface[j]];
  }

  MPI_Alltoallv(sdata,    scount, sdist, XMPI_INT,
                m->nface, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  /* FACES */

  // Hash and request faces
  m->nfaces = 0;
  int *f_rcount = (int *)XCALLOC(npc, sizeof(int));
  int *f_scount = (int *)XMALLOC(npc * sizeof(int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pc = &m->gcells[c_rdist[i]];

    for (int j = 0; j < c_rcount[i]; j++) {
      E_Int cell = mCT[pc[j]];
      assert(cell >= 0 && cell < m->ncells);

      for (E_Int k = m->indPH[cell]; k < m->indPH[cell+1]; k++) {
        E_Int face = m->nface[k];

        if (mFT.find(face) == mFT.end()) {
          mFT[face] = m->nfaces++;
          f_rcount[i]++;
        }
      }
    }
  }

  MPI_Alltoall(f_rcount, 1, MPI_INT, f_scount, 1, MPI_INT, MPI_COMM_WORLD);

  int *f_rdist = (int *)XMALLOC((npc+1) * sizeof(int));
  int *f_sdist = (int *)XMALLOC((npc+1) * sizeof(int));
  f_rdist[0] = f_sdist[0] = 0;
  for (E_Int i = 0; i < npc; i++) {
    f_rdist[i+1] = f_rdist[i] + f_rcount[i];
    f_sdist[i+1] = f_sdist[i] + f_scount[i];
  }

  m->gfaces = (E_Int *)XMALLOC(m->nfaces * sizeof(E_Int));

  E_Int *sfaces = (E_Int *)XMALLOC(f_sdist[npc] * sizeof(E_Int));

  for (E_Int i = 0; i < npc; i++) idx[i] = f_rdist[i];

  int *vfaces = (E_Int *)XCALLOC(m->nfaces, sizeof(int));
  
  for (E_Int i = 0; i < npc; i++) {
    E_Int *pc = &m->gcells[c_rdist[i]];

    for (int j = 0; j < c_rcount[i]; j++) {
      E_Int cell = mCT[pc[j]];

      for (E_Int k = m->indPH[cell]; k < m->indPH[cell+1]; k++) {
        E_Int face = m->nface[k];

        if (!vfaces[mFT[face]]) {
          vfaces[mFT[face]]++;
          m->gfaces[idx[i]++] = face;
        }
      }
    }
  }

  // Global face indices       
  MPI_Alltoallv(m->gfaces, f_rcount, f_rdist, XMPI_INT,
                sfaces,    f_scount, f_sdist, XMPI_INT,
                MPI_COMM_WORLD);

  // Send face strides
  E_Int *f_stride = (E_Int *)XMALLOC(f_sdist[npc] * sizeof(E_Int));

  m->indPG = (E_Int *)XMALLOC((m->nfaces+1) * sizeof(E_Int));

  for (E_Int i = 0; i < npc; i++) scount[i] = 0;
  
  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *ps = &f_stride[f_sdist[i]];

    for (int j = 0; j < f_scount[i]; j++) {
      E_Int face = MFT[pf[j]];
      assert(face < M->nfaces);

      E_Int stride = M->indPG[face+1] - M->indPG[face];
      ps[j] = stride;
      scount[i] += stride;
    }
  }

  MPI_Alltoallv(f_stride,   f_scount, f_sdist, XMPI_INT,
                m->indPG+1, f_rcount, f_rdist, XMPI_INT,
                MPI_COMM_WORLD);

  m->indPG[0] = 0;
  for (E_Int i = 0; i < m->nfaces; i++) m->indPG[i+1] += m->indPG[i];

  MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

  // ngon
  for (E_Int i = 0; i < npc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }
  
  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));

  m->ngon = (E_Int *)XMALLOC(rdist[npc] * sizeof(E_Int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *pn = &sdata[sdist[i]];

    for (int j = 0; j < f_scount[i]; j++) {
      E_Int face = MFT[pf[j]];
      assert(face < M->nfaces);
      
      for (E_Int k = M->indPG[face]; k < M->indPG[face+1]; k++)
        *pn++ = M->gpoints[M->ngon[k]];
    }
  }

  MPI_Alltoallv(sdata  , scount, sdist, XMPI_INT,
                m->ngon, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  /* POINTS */

  // Request points
  m->npoints = 0;
  
  int *p_rcount = (int *)XCALLOC(npc, sizeof(int));
  int *p_scount = (int *)XMALLOC(npc * sizeof(int));
  
  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &m->gfaces[f_rdist[i]];

    for (int j = 0; j < f_rcount[i]; j++) {
      E_Int face = mFT[pf[j]];

      for (E_Int k = m->indPG[face]; k < m->indPG[face+1]; k++) {
        E_Int point = m->ngon[k];

        if (mPT.find(point) == mPT.end()) {
          mPT[point] = m->npoints++;
          p_rcount[i]++;
        }
      }
    }
  }

  MPI_Alltoall(p_rcount, 1, MPI_INT, p_scount, 1, MPI_INT, MPI_COMM_WORLD);

  int *p_sdist = (int *)XMALLOC((npc+1) * sizeof(int));
  int *p_rdist = (int *)XMALLOC((npc+1) * sizeof(int));
  p_sdist[0] = p_rdist[0] = 0;
  for (E_Int i = 0; i < npc; i++) {
    p_sdist[i+1] = p_sdist[i] + p_scount[i];
    p_rdist[i+1] = p_rdist[i] + p_rcount[i];
  }

  E_Int *spoints = (E_Int *)XMALLOC(p_sdist[npc] * sizeof(E_Int));

  m->gpoints = (E_Int *)XMALLOC(p_rdist[npc] * sizeof(E_Int));

  int *vpoints = (int *)XCALLOC(m->npoints, sizeof(int));

  for (E_Int i = 0; i < npc; i++) idx[i] = p_rdist[i];

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &m->gfaces[f_rdist[i]];

    for (int j = 0; j < f_rcount[i]; j++) {
      E_Int face = mFT[pf[j]];

      for (E_Int k = m->indPG[face]; k < m->indPG[face+1]; k++) {
        E_Int point = m->ngon[k];
        assert(mPT.find(point) != mPT.end());

        if (!vpoints[mPT[point]]) {
          vpoints[mPT[point]]++;
          m->gpoints[idx[i]++] = point;
        }
      }
    }
  }

  MPI_Alltoallv(m->gpoints, p_rcount, p_rdist, XMPI_INT,
                spoints,    p_scount, p_sdist, XMPI_INT,
                MPI_COMM_WORLD);

  // Send xyz
  m->x = (E_Float *)XMALLOC(m->npoints * sizeof(E_Float));
  m->y = (E_Float *)XMALLOC(m->npoints * sizeof(E_Float));
  m->z = (E_Float *)XMALLOC(m->npoints * sizeof(E_Float));

  E_Float *sx = (E_Float *)XMALLOC(p_sdist[npc] * sizeof(E_Float));
  E_Float *sy = (E_Float *)XMALLOC(p_sdist[npc] * sizeof(E_Float));
  E_Float *sz = (E_Float *)XMALLOC(p_sdist[npc] * sizeof(E_Float));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pp = &spoints[p_sdist[i]];
    E_Float *px = &sx[p_sdist[i]];
    E_Float *py = &sy[p_sdist[i]];
    E_Float *pz = &sz[p_sdist[i]];

    for (int j = 0; j < p_scount[i]; j++) {
      E_Int point = MPT[pp[j]];
      *px++ = M->x[point];
      *py++ = M->y[point];
      *pz++ = M->z[point];
    }
  }

  MPI_Alltoallv(sx,   p_scount, p_sdist, MPI_DOUBLE,
                m->x, p_rcount, p_rdist, MPI_DOUBLE,
                MPI_COMM_WORLD);
  
  MPI_Alltoallv(sy,   p_scount, p_sdist, MPI_DOUBLE,
                m->y, p_rcount, p_rdist, MPI_DOUBLE,
                MPI_COMM_WORLD);
  
  MPI_Alltoallv(sz,   p_scount, p_sdist, MPI_DOUBLE,
                m->z, p_rcount, p_rdist, MPI_DOUBLE,
                MPI_COMM_WORLD);

  
  for (E_Int i = 0; i < m->ncells; i++) {
    for (E_Int j = m->indPH[i]; j < m->indPH[i+1]; j++) {
      m->nface[j] = mFT[m->nface[j]];
    }
  }

  for (E_Int i = 0; i < m->nfaces; i++) {
    for (E_Int j = m->indPG[i]; j < m->indPG[i+1]; j++) {
      m->ngon[j] = mPT[m->ngon[j]];
    }
  }

  // Make parent elements
  m->owner = (E_Int *)XMALLOC(m->nfaces * sizeof(E_Int));
  m->neigh = (E_Int *)XMALLOC(m->nfaces * sizeof(E_Int));

  Orient_boundary(m);
  Build_own_nei(m);

  // Exchange ref_data
  m->ref_data = (int *)XMALLOC(m->ncells * sizeof(int));

  int *sref = (int *)XMALLOC(c_sdist[npc] * sizeof(int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pc = &scells[c_sdist[i]];
    E_Int *pr = &sref[c_sdist[i]];

    for (int j = 0; j < c_scount[i]; j++) {
      E_Int cell = MCT[pc[j]];
      *pr++ = M->ref_data[cell];
    }
  }
  
  MPI_Alltoallv(sref,        c_scount, c_sdist, MPI_INT,
                m->ref_data, c_rcount, c_rdist, MPI_INT,
                MPI_COMM_WORLD);

  /* BOUNDARY */

  // Note(Imad): We count all boundaries, much simpler to code!
  m->nbc = M->nbc;
  m->ptlists = (E_Int **)XMALLOC(m->nbc * sizeof(E_Int *));
  m->bcsizes = (E_Int *) XMALLOC(m->nbc * sizeof(E_Int));
  m->bcnames = (char **) XMALLOC(m->nbc * sizeof(char *));
  m->nbf = 0;

  for (E_Int i = 0; i < M->nbc; i++) {

    // Exchange ptlist by ptlist
    m->bcnames[i] = (char *)XMALLOC((strlen(M->bcnames[i])+1)*sizeof(char));
    strcpy(m->bcnames[i], M->bcnames[i]);

    E_Int *ptlist = M->ptlists[i];

    for (E_Int j = 0; j < npc; j++) scount[j] = 0;

    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      E_Int own = M->owner[face];
      E_Int where = cmap[own];
      scount[where] += 1;
    }

    MPI_Alltoall(scount, 1, MPI_INT,
                 rcount, 1, MPI_INT,
                 MPI_COMM_WORLD);

    for (E_Int j = 0; j < npc; j++) {
      sdist[j+1] = sdist[j] + scount[j];
      rdist[j+1] = rdist[j] + rcount[j];
    }

    m->bcsizes[i] = rdist[npc];
    m->nbf += m->bcsizes[i];

    m->ptlists[i] = (E_Int *)XMALLOC(m->bcsizes[i] * sizeof(E_Int));

    sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));

    for (E_Int j = 0; j < npc; j++) idx[j] = sdist[j];

    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      E_Int own = M->owner[face];
      E_Int where = cmap[own];
      sdata[idx[where]++] = M->gfaces[face];
    }

    MPI_Alltoallv(sdata,         scount, sdist, XMPI_INT,
                  m->ptlists[i], rcount, rdist, XMPI_INT,
                  MPI_COMM_WORLD);
    
    // Replace with local face indices
    ptlist = m->ptlists[i];

    for (E_Int j = 0; j < m->bcsizes[i]; j++) {
      assert(mFT.find(ptlist[j]) != mFT.end());
      ptlist[j] = mFT[ptlist[j]];
    }
  }

  /* COMM PATCHES */

  // Isolate pfaces
  memset(vfaces, 0, m->nfaces*sizeof(E_Int));

  int spcount = 0; // How many pfaces am i requesting

  // Eliminate boundary faces
  for (E_Int i = 0; i < m->nbc; i++) {
    E_Int *ptlist = m->ptlists[i];

    for (E_Int j = 0; j < m->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      assert(face >= 0 && face < m->nfaces);
      vfaces[face] = 1;
    }
  }

  // Eliminate internal faces
  for (E_Int i = 0; i < m->nfaces; i++) {
    // Skip boundary faces
    if (vfaces[i] == 1) continue;

    if (m->neigh[i] != -1) { // Skip internal faces
      vfaces[i] = 1;
    } else {
      spcount++;
    }
  }

  E_Int *pfaces = (E_Int *)XMALLOC(spcount * sizeof(E_Int));
  E_Int *f_it = pfaces;
  for (E_Int i = 0; i < m->nfaces; i++) {
    if (vfaces[i] == 0)
      *f_it++ = m->gfaces[i];
  }

  E_Int *RCOUNT = (E_Int *)XMALLOC(npc * sizeof(E_Int));
  MPI_Allgather(&spcount, 1, MPI_INT, RCOUNT, 1, MPI_INT, MPI_COMM_WORLD);

  E_Int *RDIST = (E_Int *)XMALLOC((npc+1) * sizeof(E_Int));
  RDIST[0] = 0;
  for (E_Int i = 0; i < npc; i++) {
    RDIST[i+1] = RDIST[i] + RCOUNT[i];
  }

  E_Int *RECV = (E_Int *)XMALLOC(RDIST[npc] * sizeof(E_Int));
  MPI_Allgatherv(pfaces, spcount, XMPI_INT,
                 RECV, RCOUNT, RDIST, XMPI_INT,
                 MPI_COMM_WORLD);
  
  memset(scount, 0, npc*sizeof(E_Int));
  
  // Look for requested faces in global face table
  for (E_Int i = 0; i < npc; i++) {
    E_Int *ptr = &RECV[RDIST[i]];

    if (i == m->pid) continue;
    
    for (E_Int j = 0; j < RCOUNT[i]; j++) {
      E_Int gface = ptr[j];
      auto search = mFT.find(gface);
      if (search != mFT.end()) {
        E_Int lface = search->second;
        assert(lface >= 0 && lface < m->nfaces);
        assert(m->neigh[lface] == -1);
        scount[i] += 2;
      }
    }
  }

  MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < npc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));
  E_Int *rdata = (E_Int *)XMALLOC(rdist[npc] * sizeof(E_Int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *ptr = &RECV[RDIST[i]];
    E_Int *ps = &sdata[sdist[i]];

    if (i == m->pid) continue;
    
    for (E_Int j = 0; j < RCOUNT[i]; j++) {
      E_Int gface = ptr[j];
      auto search = mFT.find(gface);
      if (search != mFT.end()) {
        E_Int lface = search->second;
        assert(lface >= 0 && lface < m->nfaces);
        *ps++ = gface;
        *ps++ = m->gcells[m->owner[lface]];
      }
    }
  }

  MPI_Alltoallv(sdata, scount, sdist, XMPI_INT,
                rdata, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  m->npf = rdist[npc]/2;
  m->npatches = 0;

  for (E_Int i = 0; i < npc; i++) {
    if (rcount[i] == 0) continue;
    m->npatches++;
  }

  m->patches = (Patch *)XCALLOC(m->npatches, sizeof(Patch));

  m->npatches = 0;

  for (E_Int i = 0; i < npc; i++) {
    if (rcount[i] == 0) continue;
    
    m->patches[m->npatches].nf = rcount[i]/2;
    m->patches[m->npatches].pf = (E_Int *)XMALLOC(rcount[i] * sizeof(E_Int));
    m->patches[m->npatches].pn = (E_Int *)XMALLOC(rcount[i] * sizeof(E_Int));
    m->patches[m->npatches].nei = i;

    E_Int *ptr = &rdata[rdist[i]];

    Patch *P = &m->patches[m->npatches];

    E_Int k = 0;

    for (E_Int j = 0; j < rcount[i];) {
      P->pf[k] = ptr[j++];
      P->pn[k] = ptr[j++];
      assert(mCT.find(P->pn[k]) == mCT.end());
      k++;
    }

    m->npatches++;
  }

  m->nif = m->nfaces - m->nbf - m->npf;

  for (E_Int i = 0; i < m->npatches; i++) {
    Patch *P = &m->patches[i];

    std::sort(P->pf, P->pf+P->nf);
    std::sort(P->pn, P->pn+P->nf, [&] (E_Int i, E_Int j) {
      return P->pf[i] < P->pf[j];
    });

    // Replace with local faces
    for (E_Int j = 0; j < P->nf; j++)
      P->pf[j] = mFT[P->pf[j]];
  }
  
  /* EDGE CENTER MAP */

  // Send edge+center is face was sent

  memset(scount, 0, npc * sizeof(int));

  Edge E;
  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];

    for (int j = 0; j < f_scount[i]; j++) {
      E_Int lface = MFT[pf[j]];

      E_Int np = -1;
      E_Int *pn = get_face(lface, np, M->ngon, M->indPG);

      for (E_Int k = 0; k < np; k++) {
        E_Int p = pn[k];
        E_Int q = pn[(k+1)%np];

        E.set(p, q);

        auto search = M->ecenter->find(E);

        if (search != M->ecenter->end())
          scount[i] += 3;
      }
    }
  }

  MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < npc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));
  rdata = (E_Int *)XRESIZE(rdata, rdist[npc] * sizeof(E_Int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *ptr = &sdata[sdist[i]];

    for (int j = 0; j < f_scount[i]; j++) {
      E_Int lface = MFT[pf[j]];

      E_Int np = -1;
      E_Int *pn = get_face(lface, np, M->ngon, M->indPG);

      for (E_Int k = 0; k < np; k++) {
        E_Int p = pn[k];
        E_Int q = pn[(k+1)%np];

        E.set(p, q);

        auto search = M->ecenter->find(E);

        if (search != M->ecenter->end()) {
          *ptr++ = M->gpoints[search->first.p0_];
          *ptr++ = M->gpoints[search->first.p1_];
          *ptr++ = M->gpoints[search->second];
        }
      }
    }
  }

  MPI_Alltoallv(sdata, scount, sdist, XMPI_INT,
                rdata, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  m->ecenter = new std::map<Edge, E_Int>;

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pp = &rdata[rdist[i]];

    int j = 0;
    for (; j < rcount[i];) {
      E_Int p = mPT[pp[j++]];
      E_Int q = mPT[pp[j++]];
      E_Int c = mPT[pp[j++]];
      
      E.set(p, q);
      m->ecenter->insert({E, c});
    }
  }

  /* CELL TREE */
  Tree *CT = M->cellTree;

  // First, send everything but children

  for (E_Int i = 0; i < npc; i++) {
    scount[i] = 4*c_scount[i];
    rcount[i] = 4*c_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));
  rdata = (E_Int *)XRESIZE(rdata, rdist[npc] * sizeof(E_Int));

  for (E_Int i = 0; i < npc; i++) idx[i] = sdist[i];

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pc = &scells[c_sdist[i]];
    E_Int *ptr = &sdata[sdist[i]];

    for (int j = 0; j < c_scount[i]; j++) {
      E_Int cell = MCT[pc[j]];
      
      *ptr++ = M->gcells[CT->parent_[cell]];
      *ptr++ = CT->level_[cell];
      *ptr++ = CT->type_[cell];
      *ptr++ = CT->state_[cell];
    }
  }

  MPI_Alltoallv(sdata, scount, sdist, XMPI_INT,
                rdata, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);
  
  m->cellTree = new Tree(m->ncells);
  Tree *ct = m->cellTree;

  // TODO(Imad): Minimize the ugly copying
  E_Int *ptr = rdata;
  for (E_Int i = 0; i < m->ncells; i++) {
    ct->parent_[i] = mCT[*ptr++];
    ct->level_[i] = *ptr++;
    ct->type_[i] = *ptr++;
    ct->state_[i] = *ptr++;
  }

  // Count children
  memset(scount, 0, npc*sizeof(int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pc = &scells[c_sdist[i]];

    for (int j = 0; j < c_scount[i]; j++) {
      E_Int cell = MCT[pc[j]];

      Children *cur = CT->children(cell);

      while (cur) {
        scount[i] += 2 + cur->n; // ngen + nchild per gen + children
        cur = cur->next;
      }
    }
  }

  MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < npc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));
  rdata = (E_Int *)XRESIZE(rdata, rdist[npc] * sizeof(E_Int));

  if (sdist[npc]) {
    for (E_Int i = 0; i < npc; i++) {
      E_Int *pc = &scells[c_sdist[i]];
      E_Int *ptr = &sdata[sdist[i]];

      for (int j = 0; j < c_scount[i]; j++) {
        E_Int cell = MCT[pc[j]];

        Children *cur = CT->children(cell);

        E_Int &ngen = *ptr++;
        E_Int &nchild = *ptr++;
        ngen = 0;
        nchild = 0;

        while (cur) {
          ngen += 1;
          nchild = cur->n;

          for (E_Int k = 0; k < cur->n; k++)
            *ptr++ = M->gcells[cur->pc[k]];

          cur = cur->next;
        }
      }
    }
  }

  MPI_Alltoallv(sdata, scount, sdist, XMPI_INT,
                rdata, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  if (rdist[npc]) {
    for (E_Int i = 0; i < npc; i++) {
      E_Int *pc = &m->gcells[c_rdist[i]];
      E_Int *ptr = &rdata[rdist[i]];

      for (E_Int j = 0; j < rcount[i];) {
        E_Int cell = MCT[pc[j]];

        E_Int ngen = *ptr++;
        E_Int nchild = *ptr++;

        if (ngen == 0) continue;

        Children *tmp = (Children *)XMALLOC(sizeof(Children));
        tmp->next = NULL;
        Children *cur = tmp->next;

        for (E_Int k = 0; k < ngen; j++) {
          Children *new_gen = (Children *)XMALLOC(sizeof(Children));
          new_gen->n = nchild;
          new_gen->next = NULL;

          for (E_Int l = 0; l < nchild; l++)
            new_gen->pc[l] = mCT[*ptr++];

          cur = new_gen;
          cur = cur->next;
        }

        ct->children_[cell] = tmp->next;
        XFREE(tmp);
      }
    }
  }

  /* FACE TREE */
  Tree *FT = M->faceTree;

  // First, send everything but children

  for (E_Int i = 0; i < npc; i++) {
    scount[i] = 4*f_scount[i];
    rcount[i] = 4*f_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));
  rdata = (E_Int *)XRESIZE(rdata, rdist[npc] * sizeof(E_Int));


  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *ptr = &sdata[sdist[i]];

    for (int j = 0; j < f_scount[i]; j++) {
      E_Int face = MFT[pf[j]];
      assert(face >= 0 && face < M->nfaces);
      *ptr++ = M->gfaces[FT->parent_[face]];
      *ptr++ = FT->level_[face];
      *ptr++ = FT->type_[face];
      *ptr++ = FT->state_[face];
    }
  }

  MPI_Alltoallv(sdata, scount, sdist, XMPI_INT,
                rdata, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  
  m->faceTree = new Tree(m->nfaces);
  Tree *ft = m->faceTree;

  // TODO(Imad): again with the ugly copying!
  ptr = rdata;
  for (E_Int i = 0; i < m->nfaces; i++) {
    ft->parent_[i] = mFT[*ptr++];
    ft->level_[i] = *ptr++;
    ft->type_[i] = *ptr++;
    ft->state_[i] = *ptr++;
  }

  // Count children
  memset(scount, 0, npc*sizeof(int));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];

    for (int j = 0; j < f_scount[i]; j++) {
      E_Int face = MFT[pf[j]];

      Children *cur = FT->children(face);

      while (cur) {
        scount[i] += 2 + cur->n; // ngen + nchild per gen + children
        cur = cur->next;
      }
    }
  }

  MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < npc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata = (E_Int *)XRESIZE(sdata, sdist[npc] * sizeof(E_Int));
  rdata = (E_Int *)XRESIZE(rdata, rdist[npc] * sizeof(E_Int));

  if (sdist[npc]) {
    for (E_Int i = 0; i < npc; i++) {
      E_Int *pf = &sfaces[f_sdist[i]];
      E_Int *ptr = &sdata[sdist[i]];

      for (int j = 0; j < f_scount[i]; j++) {
        E_Int face =  MFT[pf[j]];

        Children *cur = FT->children(face);

        E_Int &ngen = *ptr++;
        E_Int &nchild = *ptr++;
        ngen = 0;
        nchild = 0;

        while (cur) {
          ngen += 1;
          nchild = cur->n;

          for (E_Int k = 0; k < cur->n; k++)
            *ptr++ = M->gfaces[cur->pc[k]];

          cur = cur->next;
        }
      }
    }
  }

  MPI_Alltoallv(sdata, scount, sdist, XMPI_INT,
                rdata, rcount, rdist, XMPI_INT,
                MPI_COMM_WORLD);

  if (rdist[npc]) {
    for (E_Int i = 0; i < npc; i++) {
      E_Int *pf = &m->gfaces[f_rdist[i]];
      E_Int *ptr = &rdata[rdist[i]];

      for (E_Int j = 0; j < rcount[i];) {
        E_Int face = MFT[pf[j]];

        E_Int ngen = *ptr++;
        E_Int nchild = *ptr++;

        if (ngen == 0) continue;

        Children *tmp = (Children *)XMALLOC(sizeof(Children));
        tmp->next = NULL;
        Children *cur = tmp->next;

        for (E_Int k = 0; k < ngen; j++) {
          Children *new_gen = (Children *)XMALLOC(sizeof(Children));
          new_gen->n = nchild;
          new_gen->next = NULL;

          for (E_Int l = 0; l < nchild; l++)
            new_gen->pc[l] = mFT[*ptr++];

          cur = new_gen;
          cur = cur->next;
        }

        ft->children_[face] = tmp->next;
        XFREE(tmp);
      }
    }
  }

  //if (m->pid == 3)
  //  ft->print_children();

  /* FACE CENTERS */
  m->fc = (E_Float *)XMALLOC(3*m->nfaces * sizeof(E_Float));

  for (E_Int i = 0; i < npc; i++) {
    scount[i] = 3*f_scount[i];
    rcount[i] = 3*f_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  E_Float *sfc = (E_Float *)XMALLOC(sdist[npc] * sizeof(E_Float));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Float *ptr = &sfc[sdist[i]];

    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = MFT[pf[j]];
      E_Float *fc = &M->fc[3*face];

      for (E_Int k = 0; k < 3; k++)
        *ptr++ = *fc++;
    }
  }

  MPI_Alltoallv(sfc,   scount, sdist, MPI_DOUBLE,
                m->fc, rcount, rdist, MPI_DOUBLE,
                MPI_COMM_WORLD);
  
  /* CELL CENTERS */
  m->cx = (E_Float *)XMALLOC(m->ncells * sizeof(E_Float));
  m->cy = (E_Float *)XMALLOC(m->ncells * sizeof(E_Float));
  m->cz = (E_Float *)XMALLOC(m->ncells * sizeof(E_Float));

  E_Float *scx = (E_Float *)XMALLOC(c_sdist[npc] * sizeof(E_Float));
  E_Float *scy = (E_Float *)XMALLOC(c_sdist[npc] * sizeof(E_Float));
  E_Float *scz = (E_Float *)XMALLOC(c_sdist[npc] * sizeof(E_Float));

  for (E_Int i = 0; i < npc; i++) {
    E_Int *pc = &scells[c_sdist[i]];
    E_Float *px = &scx[c_sdist[i]];
    E_Float *py = &scy[c_sdist[i]];
    E_Float *pz = &scz[c_sdist[i]];
    
    for (E_Int j = 0; j < c_scount[i]; j++) {
      E_Int cell = MCT[pc[j]];

      *px++ = M->cx[cell];
      *py++ = M->cy[cell];
      *pz++ = M->cz[cell];
    }
  }

  MPI_Alltoallv(scx,  c_scount, c_sdist, MPI_DOUBLE,
               m->cx, c_rcount, c_rdist, MPI_DOUBLE,
               MPI_COMM_WORLD);
  
  MPI_Alltoallv(scy,  c_scount, c_sdist, MPI_DOUBLE,
               m->cy, c_rcount, c_rdist, MPI_DOUBLE,
               MPI_COMM_WORLD);
  
  MPI_Alltoallv(scz,  c_scount, c_sdist, MPI_DOUBLE,
               m->cz, c_rcount, c_rdist, MPI_DOUBLE,
               MPI_COMM_WORLD);

  // Copy parameters
  m->Tr = M->Tr;
  m->Tu = M->Tu;
  m->eps = M->eps;
  m->hmin = M->hmin;
  m->hmax = M->hmax;
  m->unrefine = M->unrefine;
  if (M->mode_2D) {
    m->mode_2D = (E_Float *)XMALLOC(3 * sizeof(E_Float));
    memcpy(m->mode_2D, M->mode_2D, 3*sizeof(E_Float));
  }

  // Free
  XFREE(c_scount);
  XFREE(c_rcount);
  XFREE(c_sdist);
  XFREE(c_rdist);
  XFREE(scells);
  XFREE(c_stride);
  XFREE(idx);
  XFREE(scount);
  XFREE(rcount);
  XFREE(sdist);
  XFREE(rdist);
  XFREE(sdata);
  XFREE(rdata);
  XFREE(f_rcount);
  XFREE(f_scount);
  XFREE(f_rdist);
  XFREE(f_sdist);
  XFREE(sfaces);
  XFREE(vfaces);
  XFREE(f_stride);
  XFREE(p_rcount);
  XFREE(p_scount);
  XFREE(p_rdist);
  XFREE(p_sdist);
  XFREE(spoints);
  XFREE(sx);
  XFREE(sy);
  XFREE(sz);
  XFREE(sref);
  XFREE(sfc);
  XFREE(scx);
  XFREE(scy);
  XFREE(scz);
  XFREE(RCOUNT);
  XFREE(RDIST);
  XFREE(RECV);
  XFREE(pfaces);
  mesh_drop(M);

  reorder_cells(m);

  // Check everything is in place
  assert(check_canon_cells(m));

  return m;
}

AMesh *Redistribute_mesh(AMesh *M)
{
  // Cell weights
  E_Int *CWGT = (E_Int *)XMALLOC(M->ncells * sizeof(E_Int));
  
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int wgt = std::max(M->ref_data[i], 0);
    CWGT[i] = 1 << wgt;
  }

  // Construct dual graph
  E_Int *XNEIS = (E_Int *)XCALLOC((M->ncells+1), sizeof(E_Int));

  E_Int *CADJ = make_dual_graph(M, XNEIS);

  // Get cell map
  E_Int *cmap = map_graph(M, XNEIS, CADJ, CWGT);

  // Produce a new mesh based on cell map
  AMesh *new_mesh = reconstruct_mesh(M, cmap);

  XFREE(CWGT);
  XFREE(XNEIS);
  XFREE(CADJ);
  XFREE(cmap);

  return new_mesh;
}

void Comm_interface_data_f(AMesh *M, E_Float *data, E_Int stride, E_Float **rbuf)
{
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    P->sbuf_f = (E_Float *)XRESIZE(P->sbuf_f, stride*P->nf*sizeof(E_Float));

    E_Int l = 0;
    for (E_Int j = 0; j < P->nf; j++) {
      E_Int own = M->owner[P->pf[j]];
      E_Float *ptr = &data[stride*own];
      for (E_Int k = 0; k < stride; k++)
        P->sbuf_f[l++] = ptr[k];
    }

    assert(l == stride*P->nf);

    MPI_Irecv(rbuf[i], stride*P->nf, MPI_DOUBLE, P->nei, 0, MPI_COMM_WORLD, &M->req[M->nrq++]);
    MPI_Isend(P->sbuf_f, stride*P->nf, MPI_DOUBLE, P->nei, 0, MPI_COMM_WORLD, &M->req[M->nrq++]);
    
    assert(M->nrq < 2*M->npc);
  }

  Comm_waitall(M);
}

void Comm_interface_data_i(AMesh *M, E_Int *data, E_Int stride, E_Int **rbuf)
{
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, stride*P->nf*sizeof(E_Int));

    E_Int l = 0;
    for (E_Int j = 0; j < P->nf; j++) {
      E_Int own = M->owner[P->pf[j]];
      E_Int *ptr = &data[stride*own];
      for (E_Int k = 0; k < stride; k++)
        P->sbuf_i[l++] = ptr[k];
    }

    assert(l == stride*P->nf);

    MPI_Irecv(rbuf[i], stride*P->nf, XMPI_INT, P->nei, 0, MPI_COMM_WORLD, &M->req[M->nrq++]);
    MPI_Isend(P->sbuf_i, stride*P->nf, XMPI_INT, P->nei, 0, MPI_COMM_WORLD, &M->req[M->nrq++]);
    
    assert(M->nrq < 2*M->npc);
  }

  Comm_waitall(M);
}
