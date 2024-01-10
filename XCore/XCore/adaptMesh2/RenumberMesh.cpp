#include "Proto.h"
#include <queue>
#include <numeric>
#include "Proto.h"

static
void renumber_facetree(Element **facetree, AMesh *M,
  const std::vector<E_Int> &new_faces)
{
  Element **new_facetree = (Element **)XMALLOC(M->nfaces * sizeof(Element *));

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int new_face = new_faces[i];
    new_facetree[new_face] = facetree[i];
    Element *elem = new_facetree[new_face];
    for (E_Int j = 0; j < elem->nchildren; j++)
      elem->children[j] = new_faces[elem->children[j]];
    elem->parent = new_faces[elem->parent];
  }

  XFREE(M->faceTree);
  M->faceTree = new_facetree;
}

static
void renumber_celltree(Element **celltree, AMesh *M,
  const std::vector<E_Int> &new_cells)
{
  Element **new_celltree = (Element **)XMALLOC(M->ncells * sizeof(Element *));

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int new_cell = new_cells[i];
    new_celltree[new_cell] = celltree[i];
    Element *elem = new_celltree[new_cell];
    for (E_Int j = 0; j < elem->nchildren; j++)
      elem->children[j] = new_cells[elem->children[j]];
    elem->parent = new_cells[elem->parent];
  }

  XFREE(M->cellTree);
  M->cellTree = new_celltree;
}

static
void renumber_own_nei(AMesh *M, const std::vector<E_Int> &new_faces)
{
  E_Int *new_owner = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  E_Int *new_neigh = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  //memset(new_owner, -1, M->nfaces*sizeof(E_Int));
  //memset(new_neigh, -1, M->nfaces*sizeof(E_Int));

  for (E_Int i = 0; i < M->nfaces; i++) {
    new_owner[new_faces[i]] = M->owner[i];
    new_neigh[new_faces[i]] = M->neigh[i];
  }

  XFREE(M->owner);
  XFREE(M->neigh);
  M->owner = new_owner;
  M->neigh = new_neigh;
}

// Shuffle face order: internal faces first
// Boundary faces are grouped together
void init_mesh_numbering(AMesh *M)
{
  M->nbf = 0;

  // Renumber boundary faces
  // Map[old bface] -> nif + i
  std::vector<E_Int> new_faces(M->nfaces, -1);

  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *bfaces = M->ptlists[i];
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int bface = bfaces[j];
      new_faces[bface] = M->nbf++;
    }
  }

  // Shift
  E_Int nif = M->nfaces - M->nbf;
  M->nif = 0;

  for (E_Int i = 0; i < M->nfaces; i++) {
    if (new_faces[i] == -1) new_faces[i] = M->nif++;
    else new_faces[i] += nif;
  }

  assert(M->nif + M->nbf == M->nfaces);

  // Renumber indPG
  E_Int *new_indPG = (E_Int *)XMALLOC((M->nfaces+1) * sizeof(E_Int));
  renumber_indPG(M, new_faces, new_indPG);

  // Renumber ngon
  E_Int *new_ngon = (E_Int *)XMALLOC(new_indPG[M->nfaces] * sizeof(E_Int));
  renumber_ngon(M, new_faces, new_indPG, new_ngon);

  XFREE(M->indPG);
  XFREE(M->ngon);
  M->indPG = new_indPG;
  M->ngon = new_ngon;

  // Renumber nface
  renumber_nface_shallow(M, new_faces);

  // Renumber owner and neigh
  renumber_own_nei(M, new_faces);

  // Renumber boundary pointlists
  renumber_boundary(M, new_faces);

  // Renumber comm patches
  renumber_comm_patches(M, new_faces);

  // Renumber adaptation trees
  renumber_facetree(M->faceTree, M, new_faces);
}

std::vector<E_Int> renumber_cells(AMesh *M)
{
  // We need cell-to-cell connecitivity
  std::vector<std::vector<E_Int>> cell_cells(M->ncells); 

  M->nif = 0;

  for (E_Int i = 0; i < M->nfaces; i++) {
    if (M->neigh[i] == -1) continue;
    cell_cells[M->owner[i]].push_back(M->neigh[i]);
    cell_cells[M->neigh[i]].push_back(M->owner[i]);
    M->nif++;
  }

  // Cuthill-McKee algorithm
  std::vector<E_Int> new_cells(M->ncells);
  std::queue<E_Int> next_cell;
  std::vector<E_Int> visited(M->ncells, 0);
  E_Int idx = 0; 
  std::vector<E_Int> neis;
  std::vector<E_Int> wgts;
  std::vector<E_Int> order;

  while (1) {
    // Look for least connected cell
    E_Int cur_cell = -1;
    E_Int min_wgt = E_IDX_NONE;

    for (E_Int i = 0; i < M->ncells; i++) {
      if (!visited[i] && cell_cells[i].size() < min_wgt) {
        min_wgt = cell_cells[i].size();
        cur_cell = i;
      }
    }

    if (cur_cell == -1) break;

    next_cell.push(cur_cell);

    while (!next_cell.empty()) {
      cur_cell = next_cell.front();
      next_cell.pop();

      if (visited[cur_cell]) continue;

      visited[cur_cell] = 1;

      new_cells[idx++] = cur_cell;

      const auto &neighbours = cell_cells[cur_cell];
      
      neis.clear();
      wgts.clear();

      for (size_t i = 0; i < neighbours.size(); i++) {
        E_Int nei = neighbours[i];
        if (!visited[nei]) {
          neis.push_back(nei);
          wgts.push_back(cell_cells[nei].size());
        }
      }

      order.clear();
      order.resize(wgts.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&] (E_Int i, E_Int j)
      {
        return wgts[i] < wgts[j];
      });

      for (size_t i = 0; i < order.size(); i++) {
        E_Int index = order[i];
        next_cell.push(neis[index]);
      }
    }
  }

  assert(idx == M->ncells);

  return new_cells;
}

static
std::vector<E_Int> invert(const std::vector<E_Int> &Map)
{
  std::vector<E_Int> inverse(Map.size());
  E_Int i = 0;
  for (const auto val : Map) {
    assert(val < Map.size());
    inverse[val] = i;
    i++;
  }
  return inverse;
}

// new_cells: cellOrder
// old_cells: reverseCellOrder
std::vector<E_Int> renumber_faces(AMesh *M,
  const std::vector<E_Int> &new_cells)
{
  std::vector<E_Int> new_faces(M->nfaces, -1);

  E_Int new_face = 0;

  std::vector<E_Int> neis, order;

  std::vector<E_Int> old_cells(invert(new_cells));

  for (E_Int new_cell = 0; new_cell < M->ncells; new_cell++) {
    E_Int old_cell = new_cells[new_cell];

    E_Int nf = -1;
    E_Int pf[24];

    get_full_cell(old_cell, M, nf, pf);

    neis.resize(nf);
    
    for (E_Int i = 0; i < nf; i++) {
      E_Int face = pf[i];

      E_Int old_nei = get_neighbour(old_cell, face, M);

      if (old_nei == -1) {
        assert(face >= M->nif);
        neis[i] = -1;
      } else {
        assert(face < M->nif);
        E_Int new_nei = old_cells[old_nei];
        assert(new_nei != new_cell);

        if (new_cell < new_nei) neis[i] = new_nei;
        else neis[i] = -1;
      }
    }

    order.clear();
    order.resize(neis.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&] (E_Int i, E_Int j) {
      return neis[i] < neis[j];
    });

    for (size_t i = 0; i < order.size(); i++) {
      E_Int idx = order[i];
      if (neis[idx] != -1)
        new_faces[pf[idx]] = new_face++;
    }
  }

  // Do not touch boundary faces
  for (E_Int i = new_face; i < M->nfaces; i++)
    new_faces[i] = i;

  for (E_Int i = 0; i < M->nfaces; i++)
    assert(new_faces[i] != -1);

  return new_faces;
}

void renumber_indPG(AMesh *M, const std::vector<E_Int> &new_faces,
  E_Int *new_indPG)
{
  new_indPG[0] = 0;
  for (E_Int i = 0; i < M->nfaces; i++)
    new_indPG[new_faces[i]+1] = get_stride(i, M->indPG);
  for (E_Int i = 0; i < M->nfaces; i++) new_indPG[i+1] += new_indPG[i];
  assert(new_indPG[M->nfaces] == M->indPG[M->nfaces]);
}

void renumber_ngon(AMesh *M, const std::vector<E_Int> &new_faces,
  E_Int *new_indPG, E_Int *new_ngon)
{
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int new_face = new_faces[i];
    E_Int *ptr = &new_ngon[new_indPG[new_face]];
    for (E_Int j = M->indPG[i]; j < M->indPG[i+1]; j++)
      *ptr++ = M->ngon[j];
  }
}

void renumber_indPH(AMesh *M, const std::vector<E_Int> &new_cells,
  E_Int *new_indPH)
{
  new_indPH[0] = 0;
  for (E_Int i = 0; i < M->ncells; i++)
    new_indPH[new_cells[i]+1] = get_stride(i, M->indPH);
  for (E_Int i = 0; i < M->ncells; i++) new_indPH[i+1] += new_indPH[i];
  assert(new_indPH[M->ncells] == M->indPH[M->ncells]);
}

void renumber_nface(AMesh *M, const std::vector<E_Int> &new_cells,
  E_Int *new_indPH, E_Int *new_nface)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int new_cell = new_cells[i];
    E_Int *ptr = &new_nface[new_indPH[new_cell]];
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++)
      *ptr++ = M->nface[j];
  }
}

void renumber_nface_shallow(AMesh *M, const std::vector<E_Int> &new_faces)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++)
      M->nface[j] = new_faces[M->nface[j]];
  }
}

void renumber_boundary(AMesh *M, const std::vector<E_Int> &new_faces)
{
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      ptlist[j] = new_faces[ptlist[j]];
    }
  }
}

void renumber_comm_patches(AMesh *M, const std::vector<E_Int> &new_faces)
{
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *p = &M->patches[i];
    for (E_Int j = 0; j < p->nfaces; j++)
      p->faces[j] = new_faces[p->faces[j]];
  }
}

void renumber_mesh(AMesh *M)
{
  std::vector<E_Int> new_cells = renumber_cells(M);

  std::vector<E_Int> new_faces = renumber_faces(M, new_cells);

  // Renumber indPG
  E_Int *new_indPG = (E_Int *)XMALLOC((M->nfaces+1) * sizeof(E_Int));
  renumber_indPG(M, new_faces, new_indPG);

  // Renumber ngon
  E_Int *new_ngon = (E_Int *)XMALLOC(new_indPG[M->nfaces] * sizeof(E_Int));
  renumber_ngon(M, new_faces, new_indPG, new_ngon);

  // Renumber indPH
  E_Int *new_indPH = (E_Int *)XMALLOC((M->ncells+1) * sizeof(E_Int));
  renumber_indPH(M, new_cells, new_indPH);

  // Renumber nface
  E_Int *new_nface = (E_Int *)XMALLOC(new_indPH[M->ncells] * sizeof(E_Int));
  renumber_nface(M, new_cells, new_indPH, new_nface);
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = new_indPH[i]; j < new_indPH[i+1]; j++)
      new_nface[j] = new_faces[new_nface[j]];
  }

  // Renumber owner and neigh
  E_Int *new_owner = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  E_Int *new_neigh = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));

  for (E_Int i = 0; i < M->nfaces; i++) {
    new_owner[new_faces[i]] = M->owner[i];
    new_neigh[new_faces[i]] = M->neigh[i];
  }

  for (E_Int i = 0; i < M->nfaces; i++) {
    new_owner[i] = new_cells[new_owner[i]];
    new_neigh[i] = new_cells[new_neigh[i]];
  }

  // Swap faces if necessary
  for (E_Int i = 0; i < M->nif; i++) {
    assert(new_neigh[i] != -1);
    E_Int own = new_owner[i];
    E_Int nei = new_neigh[i];

    if (nei < own) {
      E_Int np = -1;
      E_Int *pn = get_face(i, np, new_ngon, new_indPG);
      std::reverse(pn+1, pn+np);
      std::swap(new_owner[i], new_neigh[i]);
    }
  }
  XFREE(M->owner);
  XFREE(M->neigh);
  M->owner = new_owner;
  M->neigh = new_neigh;

  XFREE(M->indPG);
  XFREE(M->indPH);
  XFREE(M->ngon);
  XFREE(M->nface);
  M->indPG = new_indPG;
  M->indPH = new_indPH;
  M->ngon = new_ngon;
  M->nface = new_nface;

  // Renumber boundary pointlists
  renumber_boundary(M, new_faces);

  // Renumber comm patches
  renumber_comm_patches(M, new_faces);

  // Renumber cell/face tree
  renumber_facetree(M->faceTree, M, new_faces);
  renumber_celltree(M->cellTree, M, new_cells);
}
