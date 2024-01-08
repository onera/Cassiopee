#include "Proto.h"
#include <stack>
#include <numeric>

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
  std::stack<E_Int> next_cell;
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
      cur_cell = next_cell.top();
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

std::vector<E_Int> renumber_faces(AMesh *M,
  const std::vector<E_Int> &new_cells)
{
  std::vector<E_Int> new_faces(M->nfaces, -1);

  E_Int new_face = 0;

  std::vector<E_Int> neis, order;

  for (E_Int old_cell = 0; old_cell < M->ncells; old_cell++) {
    E_Int nf = -1;
    E_Int *pf = get_cell(old_cell, nf, M->nface, M->indPH);

    neis.resize(nf);
    
    for (E_Int i = 0; i < nf; i++) {
      E_Int face = pf[i];

      if (is_internal_face(face, M)) {
        E_Int old_nei = get_neighbour(old_cell, face, M);
        E_Int new_nei = new_cells[old_nei];
        E_Int new_cell = new_cells[old_cell];
        
        if (new_cell < new_nei) neis[i] = new_nei;
        else neis[i] = -1;
      } else {
        neis[i] = -1;
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

  return new_faces;
}

void renumber_mesh(AMesh *A, const std::vector<E_Int> &new_cells,
  const std::vector<E_Int> &new_faces)
{
  // Renumber 
}
