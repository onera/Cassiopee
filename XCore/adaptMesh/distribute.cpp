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
      E_Int own = M->owner[j];
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

void redistribute_mesh(mesh *M, const std::vector<E_Int> &ref_data)
{
  // compute cell weights
  // TODO(Imad): for now, cell weight is max ref_data
  std::vector<E_Int> cwgt(M->ncells);
  for (E_Int i = 0; i < M->ncells; i++) {
    const E_Int *pr = &ref_data[3*i];
    for (E_Int j = 0; j < 3; j++)
      cwgt[i] = std::max(std::max(pr[0], pr[1]), pr[2]);
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
}
