#include "proto.h"
#include "scotch/ptscotch.h"

static
std::vector<E_Int> make_dual_graph(mesh *M, std::vector<E_Int> &xneis)
{
  // how many internal faces
  E_Int nif = 0;
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;
    E_Int own = M->owner[i];
    xneis[own+1]++;
    xneis[nei+1]++;
    nif++;
  }

  // how many proc faces
  E_Int npf = 0;
  for (E_Int i = 0; i < M->nppatches; i++) {
    npf += M->ppatches[i].nfaces;
    for (E_Int j = 0; j < M->ppatches[i].nfaces; j++) {
      E_Int face = M->ppatches[i].faces[j];
      E_Int own = M->owner[face];
      xneis[own+1]++;
    }
  }

  std::vector<E_Int> cadj(2*nif + npf);
  for (E_Int i = 0; i < M->ncells; i++)
    xneis[i+1] += xneis[i];

  assert(xneis[M->ncells] == (2*nif + npf));
 
  // second pass: fill
  std::vector<E_Int> cneis(M->ncells, 0);
  
  // internal faces
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;
    E_Int own = M->owner[i];

    E_Int start = xneis[own];
    E_Int *ps = &cadj[start];
    ps[cneis[own]++] = M->gcells[nei];

    start = xneis[nei];
    ps = &cadj[start];
    ps[cneis[nei]++] = M->gcells[own];

    assert(cneis[nei] <= 6);
    assert(cneis[own] <= 6);
  }

  // proc faces
  for (E_Int i = 0; i < M->nppatches; i++) {
    for (E_Int j = 0; j < M->ppatches[i].nfaces; j++) {
      E_Int face = M->ppatches[i].faces[j];
      E_Int own = M->owner[face];
      E_Int gnei = M->ppatches[i].gneis[j];
      E_Int start = xneis[own];
      E_Int *ps = &cadj[start];
      ps[cneis[own]++] = gnei;

      assert(cneis[own] <= 6);
    }
  }

  return cadj;
}

mesh *redistribute_mesh(mesh *M)
{
  // compute cell weights
  std::vector<E_Int> cwgt(M->ncells, 0);
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pr = &M->ref_data[3*i];
    E_Int pow = 0;
    for (E_Int j = 0; j < 3; j++)
      pow += pr[j];
    //cwgt[i] = 1 + (1<<pow);
    cwgt[i] = 1<<pow;
    //cwgt[i] = std::max(std::max(pr[0], pr[1]), pr[2]);
  }

  std::vector<E_Int> xadj(M->ncells+1, 0);

  std::vector<E_Int> cadj = make_dual_graph(M, xadj);

  SCOTCH_Dgraph graph;
  SCOTCH_dgraphInit(&graph, MPI_COMM_WORLD);

  E_Int ret;
  ret = SCOTCH_dgraphBuild(
    &graph,
    0,
    M->ncells,
    M->ncells,
    &xadj[0],
    NULL,
    &cwgt[0],
    &M->gcells[0],
    (E_Int)cadj.size(),    
    (E_Int)cadj.size(),
    &cadj[0],
    NULL,
    NULL);

  assert(ret == 0);

  ret = SCOTCH_dgraphCheck(&graph);

  assert(ret == 0);

  SCOTCH_Strat strat;
  ret = SCOTCH_stratInit(&strat);
  if (ret != 0) {
    fprintf(stderr, "SCOTCH_startInit(): Failed to init strat\n");
  }

  std::vector<E_Int> part(M->ncells, -1);
  ret = SCOTCH_dgraphPart(&graph, M->npc, &strat, &part[0]);

  if (ret != 0) {
    fprintf(stderr, "SCOTCH_dgraphPart(): Failed to map graph\n");
  }

  if (M->pid == 0)
    printf("Graph map OK\n");

  // exchange
  E_Int nproc = M->npc;
  mesh *m = new mesh;

  // send cell ids and strides
  std::vector<int> c_scount(nproc, 0);
  std::vector<int> c_rcount(nproc);
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = part[i];
    assert(where < m->npc && where >= 0);
    c_scount[where] += 1;
  }
  
  /*
  printf("%d -> ", m->pid);
  for (E_Int i = 0; i < nproc; i++)
    printf("%d ", c_scount[i]);
  puts("");
  */

  MPI_Alltoall(&c_scount[0], 1, MPI_INT, &c_rcount[0], 1, MPI_INT,
    MPI_COMM_WORLD);

  std::vector<int> c_sdist(nproc+1);
  std::vector<int> c_rdist(nproc+1);
  c_sdist[0] = c_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    c_sdist[i+1] = c_sdist[i] + c_scount[i];
    c_rdist[i+1] = c_rdist[i] + c_rcount[i];
  }

  E_Int nncells = c_rdist[nproc];
  m->ncells = nncells;

  m->gcells = (E_Int *)XCALLOC(nncells, sizeof(E_Int));
  std::vector<E_Int> scells(c_sdist[nproc]);
  
  std::vector<E_Int> c_stride(c_sdist[nproc]);
  std::vector<E_Int> idx(nproc);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = c_sdist[i];

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = part[i];
    scells[idx[where]] = M->gcells[i];
    c_stride[idx[where]] = M->xcells[i+1] - M->xcells[i];
    idx[where]++;
  }

  MPI_Alltoallv(&scells[0], &c_scount[0], &c_sdist[0], XMPI_INT,
                m->gcells, &c_rcount[0], &c_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  m->xcells = (E_Int *)XCALLOC((nncells+1), sizeof(E_Int));
  MPI_Alltoallv(&c_stride[0], &c_scount[0], &c_sdist[0], XMPI_INT,
                m->xcells+1, &c_rcount[0], &c_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  m->xcells[0] = 0;
  for (E_Int i = 0; i < nncells; i++)
    m->xcells[i+1] += m->xcells[i];

  assert(m->xcells[m->ncells] == 6*m->ncells);

  // send NFACE connectivity
  std::vector<int> scount(nproc, 0);
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = part[i];
    scount[where] += M->xcells[i+1] - M->xcells[i];
  }

  std::vector<int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> sdist(nproc+1);
  std::vector<int> rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  std::vector<E_Int> sdata(sdist[nproc]);
  m->NFACE = (E_Int *)XCALLOC(rdist[nproc], sizeof(E_Int));

  assert(rdist[nproc] == 6*m->ncells);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = part[i];
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++)
      sdata[idx[where]++] = M->gfaces[M->NFACE[j]];
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                m->NFACE, &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  // hash cells
  for (E_Int i = 0; i < m->ncells; i++)
    m->CT[m->gcells[i]] = i;

  // hash and request faces
  E_Int nnfaces = 0;
  std::vector<int> f_rcount(nproc, 0);
  std::vector<int> f_scount(nproc);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &m->gcells[c_rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int cell = m->CT[pc[j]];
      for (E_Int k = m->xcells[cell]; k < m->xcells[cell+1]; k++) {
        E_Int face = m->NFACE[k];
        if (m->FT.find(face) == m->FT.end()) {
          m->FT[face] = nnfaces++;
          f_rcount[i]++;
        }
      }
    }
  }

  MPI_Alltoall(&f_rcount[0], 1, MPI_INT, &f_scount[0], 1, MPI_INT,
    MPI_COMM_WORLD);

  std::vector<int> f_rdist(nproc+1);
  std::vector<int> f_sdist(nproc+1);
  f_rdist[0] = f_sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    f_rdist[i+1] = f_rdist[i] + f_rcount[i];
    f_sdist[i+1] = f_sdist[i] + f_scount[i];
  }

  assert(nnfaces == f_rdist[nproc]);
  m->nfaces = nnfaces;

  m->gfaces = (E_Int *)XCALLOC(m->nfaces, sizeof(E_Int));
  std::vector<E_Int> sfaces(f_sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = f_rdist[i];

  std::vector<E_Int> vfaces(nnfaces, 0);
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &m->gcells[c_rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int cell = m->CT[pc[j]];
      for (E_Int k = m->xcells[cell]; k < m->xcells[cell+1]; k++) {
        E_Int face = m->NFACE[k];
        assert(face > 0);
        if (!vfaces[m->FT[face]]) {
          vfaces[m->FT[face]]++;
          m->gfaces[idx[i]++] = face;
        }
      }
    }
  }
          
  MPI_Alltoallv(m->gfaces, &f_rcount[0], &f_rdist[0], XMPI_INT,
                &sfaces[0], &f_scount[0], &f_sdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  // send face strides
  std::vector<E_Int> f_stride(f_sdist[nproc]);
  m->xfaces = (E_Int *)XCALLOC((m->nfaces+1), sizeof(E_Int));
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *ps = &f_stride[f_sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = M->FT[pf[j]];
      assert(face < M->nfaces);
      E_Int stride = M->xfaces[face+1] - M->xfaces[face];
      ps[j] = stride;
      scount[i] += stride;
    }
  }

  MPI_Alltoallv(&f_stride[0], &f_scount[0], &f_sdist[0], XMPI_INT,
                m->xfaces+1, &f_rcount[0], &f_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  m->xfaces[0] = 0;
  for (E_Int i = 0; i < nnfaces; i++)
    m->xfaces[i+1] += m->xfaces[i];

  assert(m->xfaces[m->nfaces] == 4*m->nfaces);

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  // send NGON
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }
  sdata.resize(sdist[nproc]);
  assert(rdist[nproc] == 4*m->nfaces);
  m->NGON = (E_Int *)XCALLOC(rdist[nproc], sizeof(E_Int));

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *pn = &sdata[sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = M->FT[pf[j]];
      assert(face < M->nfaces);
      for (E_Int k = M->xfaces[face]; k < M->xfaces[face+1]; k++)
        *pn++ = M->gpoints[M->NGON[k]];
    }
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                m->NGON, &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  // request points
  E_Int nnpoints = 0;
  std::vector<int> p_rcount(nproc, 0);
  std::vector<int> p_scount(nproc);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &m->gfaces[f_rdist[i]];
    for (E_Int j = 0; j < f_rcount[i]; j++) {
      E_Int face = m->FT[pf[j]];
      for (E_Int k = m->xfaces[face]; k < m->xfaces[face+1]; k++) {
        E_Int point = m->NGON[k];
        if (m->PT.find(point) == m->PT.end()) {
          m->PT[point] = nnpoints++;
          p_rcount[i]++;
        }
      }
    }
  }

  MPI_Alltoall(&p_rcount[0], 1, MPI_INT, &p_scount[0], 1, MPI_INT,
    MPI_COMM_WORLD);

  std::vector<int> p_sdist(nproc+1);
  std::vector<int> p_rdist(nproc+1);
  p_sdist[0] = p_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    p_sdist[i+1] = p_sdist[i] + p_scount[i];
    p_rdist[i+1] = p_rdist[i] + p_rcount[i];
  }
  std::vector<E_Int> spoints(p_sdist[nproc]);
  m->gpoints = (E_Int *)XCALLOC(p_rdist[nproc], sizeof(E_Int));
  std::vector<E_Int> vpoints(nnpoints, 0);
  assert(nnpoints == p_rdist[nproc]);
  m->npoints = nnpoints;

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = p_rdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &m->gfaces[f_rdist[i]];
    for (E_Int j = 0; j < f_rcount[i]; j++) {
      E_Int face = m->FT[pf[j]];
      for (E_Int k = m->xfaces[face]; k < m->xfaces[face+1]; k++) {
        E_Int point = m->NGON[k];
        assert(m->PT.find(point) != m->PT.end());
        if (!vpoints[m->PT[point]]) {
          vpoints[m->PT[point]]++;
          m->gpoints[idx[i]++] = point;
        }
      }
    }
  }

  MPI_Alltoallv(m->gpoints, &p_rcount[0], &p_rdist[0], XMPI_INT,
                &spoints[0], &p_scount[0], &p_sdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  // send coordinates
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] = 3*p_scount[i];
    rcount[i] = 3*p_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }
  assert(sdist[nproc] == 3*p_sdist[nproc]);
  assert(rdist[nproc] == 3*p_rdist[nproc]);
  
  std::vector<E_Float> sxyz(sdist[nproc]);
  m->xyz = (E_Float *)XCALLOC(rdist[nproc], sizeof(E_Float));

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pp = &spoints[p_sdist[i]];
    E_Float *px = &sxyz[sdist[i]];
    for (E_Int j = 0; j < p_scount[i]; j++) {
      E_Int point = M->PT[pp[j]];
      E_Float *X = &M->xyz[3*point];
      for (E_Int k = 0; k < 3; k++)
        *px++ = X[k];
    }
  }

  MPI_Alltoallv(&sxyz[0], &scount[0], &sdist[0], MPI_DOUBLE,
                m->xyz, &rcount[0], &rdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

  // exchange ref_data
  m->ref_data = (E_Int *)XCALLOC(3*m->ncells, sizeof(E_Int));
  sdata.resize(3*c_sdist[nproc]);
  E_Int *p = &sdata[0];
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &scells[c_sdist[i]];
    for (E_Int j = 0; j < c_scount[i]; j++) {
      E_Int cell = M->CT[ptr[j]];
      E_Int *pr = &M->ref_data[3*cell];
      for (E_Int k = 0; k < 3; k++)
        *p++ = pr[k];
    }
  }

  for (E_Int i = 0; i < nproc; i++) {
    scount[i] = 3*c_scount[i];
    rcount[i] = 3*c_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  assert(rdist[nproc] == 3*m->ncells);

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                m->ref_data, &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  // init parent elements
  m->owner = (E_Int *)XCALLOC(m->nfaces, sizeof(E_Int));
  m->neigh = (E_Int *)XCALLOC(m->nfaces, sizeof(E_Int));
  for (E_Int i = 0; i < m->nfaces; i++) {
    m->owner[i] = -1;
    m->neigh[i] = -1;
  }
  
  for (E_Int i = 0; i < m->xcells[m->ncells]; i++)
    m->NFACE[i] = m->FT[m->NFACE[i]];

  for (E_Int i = 0; i < m->xfaces[m->nfaces]; i++)
    m->NGON[i] = m->PT[m->NGON[i]];
  
  
  for (E_Int i = 0; i < m->ncells; i++) {
    for (E_Int j = m->xcells[i]; j < m->xcells[i+1]; j++) {
      E_Int face = m->NFACE[j];
      if (m->owner[face] == -1) m->owner[face] = i;
      else m->neigh[face] = i;
    }
  }

  if (m->ncells != 0)
    topo_init_mesh(m);

  return m;
}
