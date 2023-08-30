#include "proto.h"
#include "scotch/ptscotch.h"

static
std::vector<E_Int> make_dual_graph(mesh *M, std::vector<E_Int> &xneis)
{
  for (auto& elem : xneis)
    elem = 0;

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
    }
  }

  return cadj;
}

mesh *redistribute_mesh(mesh *M)
{
  // compute cell weights
  // TODO(Imad): how to quantify cell weights?
  std::vector<E_Int> cwgt(M->ncells, 0);
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pr = &M->ref_data[3*i];
    for (E_Int j = 0; j < 3; j++)
      cwgt[i] += pr[j];
    //cwgt[i] = std::max(std::max(pr[0], pr[1]), pr[2]);
  }

  std::vector<E_Int> xadj(M->ncells+1);

  std::vector<E_Int> cadj = make_dual_graph(M, xadj);

  SCOTCH_Dgraph graph;
  SCOTCH_dgraphInit(&graph, MPI_COMM_WORLD);

  E_Int ret;
  ret = SCOTCH_dgraphBuild(
    &graph,
    0,
    M->ncells,
    M->ncells, // change ?
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

  std::vector<E_Int> part(M->ncells);
  ret = SCOTCH_dgraphPart(&graph, M->npc, &strat, &part[0]); 

  if (ret != 0) {
    fprintf(stderr, "SCOTCH_dgraphPart(): Failed to map graph\n");
  }

  if (M->pid == 0)
    printf("Graph map OK\n");

  // exchange
  E_Int nproc = M->npc;
  E_Int rank = M->pid;
  mesh *m = new mesh;

  // send cell ids and strides
  std::vector<E_Int> c_scount(nproc, 0);
  std::vector<E_Int> c_rcount(nproc);
  for (E_Int i = 0; i < M->ncells; i++)
    c_scount[part[i]]++;

  MPI_Alltoall(&c_scount[0], 1, MPI_INT, &c_rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> c_sdist(nproc+1);
  std::vector<E_Int> c_rdist(nproc+1);
  c_sdist[0] = c_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    c_sdist[i+1] = c_sdist[i] + c_scount[i];
    c_rdist[i+1] = c_rdist[i] + c_rcount[i];
  }

  E_Int nncells = c_rdist[nproc];

  m->gcells = (E_Int *)malloc(nncells * sizeof(E_Int));
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

  MPI_Alltoallv(&scells[0], &c_scount[0], &c_sdist[0], MPI_INT,
                m->gcells, &c_rcount[0], &c_rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  m->xcells = (E_Int *)malloc((nncells+1) * sizeof(E_Int));
  MPI_Alltoallv(&c_stride[0], &c_scount[0], &c_sdist[0], MPI_INT,
                m->xcells+1, &c_rcount[0], &c_rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  m->xcells[0] = 0;
  for (E_Int i = 0; i < nncells; i++)
    m->xcells[i+1] += m->xcells[i];

  // send NFACE connectivity
  std::vector<E_Int> scount(nproc, 0);
  for (E_Int i = 0; i < M->ncells; i++)
    scount[part[i]] += M->xcells[i+1] - M->xcells[i];

  std::vector<E_Int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sdist(nproc+1);
  std::vector<E_Int> rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  std::vector<E_Int> sdata(sdist[nproc]);
  m->NFACE = (E_Int *)malloc(rdist[nproc] * sizeof(E_Int));

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int where = part[i];
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++)
      sdata[idx[where]++] = M->gfaces[M->NFACE[j]];
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                m->NFACE, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  
  // hash cells
  for (E_Int i = 0; i < nncells; i++)
    m->CT[m->gcells[i]] = i;

  // hash and request faces
  E_Int nnfaces = 0;
  std::vector<E_Int> f_rcount(nproc, 0);
  std::vector<E_Int> f_scount(nproc);

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

  MPI_Alltoall(&f_rcount[0], 1, MPI_INT, &f_scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> f_rdist(nproc+1);
  std::vector<E_Int> f_sdist(nproc+1);
  f_rdist[0] = f_sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    f_rdist[i+1] = f_rdist[i] + f_rcount[i];
    f_sdist[i+1] = f_sdist[i] + f_scount[i];
  }

  assert(nnfaces == f_rdist[nproc]);

  m->gfaces = (E_Int *)malloc(nnfaces * sizeof(E_Int));
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
        if (!vfaces[m->FT[face]]) {
          vfaces[m->FT[face]]++;
          m->gfaces[idx[i]++] = face;
        }
      }
    }
  }
          
  MPI_Alltoallv(m->gfaces, &f_rcount[0], &f_rdist[0], MPI_INT,
                &sfaces[0], &f_scount[0], &f_sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // send face strides
  std::vector<E_Int> f_stride(f_sdist[nproc]);
  m->xfaces = (E_Int *)malloc((f_rdist[nproc]+1) * sizeof(E_Int));
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *ps = &f_stride[f_sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = M->FT[pf[j]];
      E_Int stride = M->xfaces[face+1] - M->xfaces[face];
      ps[j] = stride;
      scount[i] += stride;
    }
  }

  MPI_Alltoallv(&f_stride[0], &f_scount[0], &f_sdist[0], MPI_INT,
                m->xfaces+1, &f_rcount[0], &f_rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  m->xfaces[0] = 0;
  for (E_Int i = 0; i < nnfaces; i++)
    m->xfaces[i+1] += m->xfaces[i];

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  // send NGON
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }
  sdata.resize(sdist[nproc]);
  m->NGON = (E_Int *)malloc(rdist[nproc] * sizeof(E_Int));

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *pn = &sdata[sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = M->FT[pf[j]];
      for (E_Int k = M->xfaces[face]; k < M->xfaces[face+1]; k++)
        *pn++ = M->gpoints[M->NGON[k]];
    }
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                m->NGON, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);
  
  // request points
  E_Int nnpoints = 0;
  std::vector<E_Int> p_rcount(nproc, 0);
  std::vector<E_Int> p_scount(nproc);
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

  MPI_Alltoall(&p_rcount[0], 1, MPI_INT, &p_scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> p_sdist(nproc+1);
  std::vector<E_Int> p_rdist(nproc+1);
  p_sdist[0] = p_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    p_sdist[i+1] = p_sdist[i] + p_scount[i];
    p_rdist[i+1] = p_rdist[i] + p_rcount[i];
  }
  std::vector<E_Int> spoints(p_sdist[nproc]);
  m->gpoints = (E_Int *)malloc(p_rdist[nproc] * sizeof(E_Int));
  std::vector<E_Int> vpoints(nnpoints, 0);
  assert(nnpoints == p_rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = p_rdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &m->gfaces[f_rdist[i]];
    for (E_Int j = 0; j < f_rcount[i]; j++) {
      E_Int face = m->FT[pf[j]];
      for (E_Int k = m->xfaces[face]; k < m->xfaces[face+1]; k++) {
        E_Int point = m->NGON[k];
        if (!vpoints[m->PT[point]]) {
          vpoints[m->PT[point]]++;
          m->gpoints[idx[i]++] = point;
        }
      }
    }
  }

  MPI_Alltoallv(m->gpoints, &p_rcount[0], &p_rdist[0], MPI_INT,
                &spoints[0], &p_scount[0], &p_sdist[0], MPI_INT,
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
  m->xyz = (E_Float *)malloc(rdist[nproc] * sizeof(E_Float));

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

  M->ncells = nncells;
  M->nfaces = nnfaces;
  M->npoints = nnpoints;
  return m;
}
