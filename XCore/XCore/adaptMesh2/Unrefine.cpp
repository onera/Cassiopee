#include "Proto.h"

void assign_pattern_to_adapt_faces(const std::vector<E_Int> &adapt_faces,
  std::vector<E_Int> &patterns, AMesh *M)
{
  if (!M->mode_2D) return;

  patterns.resize(adapt_faces.size());

  for (size_t i = 0; i < adapt_faces.size(); i++) {
    E_Int face = adapt_faces[i];
    E_Int own = M->owner[face];
    E_Int nf = -1;
    E_Int *pf = get_cell(own, nf, M->nface, M->indPH);

    E_Int fpos = -1;
    
    E_Int flvl = M->faceTree->level(face);
    E_Int olvl = M->cellTree->level(own);

    if (flvl == olvl)
      fpos = Get_pos(face, pf, nf);
    else
      fpos = Get_pos(master_face(face, M), pf, nf);
    
    assert(fpos != -1);

    // BOT and TOP: ISO
    // Other faces: DIR

    if (fpos == 0 || fpos == 1) patterns[i] = ISO;
    else patterns[i] = DIR;
  }
}

void update_external_pe_after_unref(E_Int cell, AMesh *M)
{
  Children *children = M->cellTree->children(cell);

  for (E_Int i = 0; i < children->n; i++) {
    E_Int child = children->pc[i];

    for (E_Int j = M->indPH[child]; j < M->indPH[child+1]; j++) {
      E_Int face = M->nface[j];

      if (M->owner[face] == child) {
        M->owner[face] = cell;
      } else if (M->neigh[face] == child) {
        M->neigh[face] = cell;
      } else {
        assert(M->owner[face] == -1);
        assert(M->neigh[face] == -1);
      }
    }
  }
}

void unrefine_cells(const std::vector<E_Int> &unref_cells, size_t start,
  size_t stop, AMesh *M)
{
  if (!M->mode_2D) {
    for (size_t i = start; i < stop; i++) {
      E_Int cell = unref_cells[i];
      E_Int type = M->cellTree->type(cell);
      switch (type) {
        case HEXA:
          H27_unrefine(cell, M);
          break;
        case TETRA:
          unrefine_tetra(cell, M);
          break;
        case PENTA:
          Pe18_unrefine(cell, M);
          break;
        case PYRA:
          unrefine_pyra(cell, M);
          break;
        default:
          assert(0);
          break;
      }
    }
  } else {
    for (size_t i = start; i < stop; i++) {
      E_Int cell = unref_cells[i];
      E_Int type = M->cellTree->type(cell);
      switch (type) {
        case HEXA:
          H18_unrefine(cell, M);
          break;
        case PENTA:
          Pe12_unrefine(cell, M);
          break;
        default:
          assert(0);
          break;
      }
    }
  }
}

void unrefine_faces(const std::vector<E_Int> &unref_faces,
  const std::vector<E_Int> &unref_patterns, size_t start,
  size_t stop, AMesh *M)
{
  if (!M->mode_2D) {
    for (size_t i = start; i < stop; i++) {
      E_Int face = unref_faces[i];
      E_Int type = M->faceTree->type(face);
      switch (type) {
        case    QUAD: Q9_unrefine(face, M); break;
        case    TRI : T6_unrefine(face, M); break;
        default     : assert(0); break;
      }
    }
  } else {
    for (size_t i = start; i < stop; i++) {
      E_Int face = unref_faces[i];
      E_Int type = M->faceTree->type(face);
      E_Int pattern = unref_patterns[i];

      switch (type) {
        case QUAD: {
          if      (pattern == ISO) Q9_unrefine(face, M);
          else if (pattern == DIR) Q6_unrefine(face, M);
          break;
        }
        case TRI:
          assert(pattern == ISO);
          T6_unrefine(face, M);
          break;
        default:
          assert(0);
          break;
      }
    }
  }
}

void unrefine_tetra(E_Int cell, AMesh *M)
{}
void unrefine_pyra(E_Int cell, AMesh *M)
{}

void get_unref_faces_and_cells(AMesh *M, std::vector<E_Int> &unref_faces,
  std::vector<E_Int> &unref_cells)
{
  unref_cells.clear();
  unref_faces.clear();

  std::set<E_Int> ucset;

  for (E_Int i = 0; i < M->ncells; i++) {
    if (M->ref_data[i] < 0) {
      ucset.insert(master_cell(i, M));
    }
  }

  for (auto ucell : ucset) unref_cells.push_back(ucell);

  // Analyse faces

  std::set<E_Int> ufset;
  std::set<E_Int> rfset;

  // Boundary
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];
    
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      
      E_Int own = M->owner[face];

      E_Int flvl = M->faceTree->level(face);
      E_Int olvl = M->cellTree->level(own);
      assert(olvl == flvl);
      E_Int oval = olvl + M->ref_data[own];

      if (oval < flvl) {
        ufset.insert(master_face(face, M));
      }
    }
  }

  for (E_Int i = 0; i < M->nfaces; i++) {
    if (M->neigh[i] == -1) continue;

    E_Int flvl = M->faceTree->level(i);

    E_Int own = M->owner[i];
    E_Int nei = M->neigh[i];

    E_Int lo = M->cellTree->level(own);
    E_Int ln = M->cellTree->level(nei);

    E_Int ro = M->ref_data[own];
    E_Int rn = M->ref_data[nei];

    E_Int oval = lo + ro;
    E_Int nval = ln + rn;

    // Refine if oval or nval is greater than flvl
    if ((oval == nval && ro == -1 && rn == 0) ||
        (oval == nval && rn == -1 && ro == 0) ||
        (oval < flvl && nval < flvl)) {
      ufset.insert(master_face(i, M));
    }
  }

  if (ucset.empty()) {
    assert(ufset.empty());
  }

  for (auto uface : ufset) unref_faces.push_back(uface);
}

void unrefine_mesh(AMesh *M, std::vector<E_Int> &unref_faces,
  std::vector<E_Int> &unref_cells)
{
  // Sort unref cells/faces by decreasing level
  std::sort(unref_faces.begin(), unref_faces.end(), [&] (E_Int i, E_Int j)
  {
    return M->faceTree->level(i) > M->faceTree->level(j);
  });

  std::sort(unref_cells.begin(), unref_cells.end(), [&] (E_Int i, E_Int j)
  {
    return M->cellTree->level(i) > M->cellTree->level(j);
  });

  std::vector<E_Int> face_distrib(1, 0), cell_distrib(1, 0);
  std::vector<E_Int> face_lvl, cell_lvl;

  E_Int max_lvl = M->faceTree->level(unref_faces[0]);

  face_lvl.push_back(max_lvl);

  for (size_t i = 1; i < unref_faces.size(); i++) {
    if (M->faceTree->level(unref_faces[i]) < max_lvl) {
      max_lvl = M->faceTree->level(unref_faces[i]);
      face_distrib.push_back(i);
      face_lvl.push_back(max_lvl);
    }
  }

  face_distrib.push_back(unref_faces.size());

  max_lvl = M->cellTree->level(unref_cells[0]);

  cell_lvl.push_back(max_lvl);

  for (size_t i = 1; i < unref_cells.size(); i++) {
    if (M->cellTree->level(unref_cells[i]) < max_lvl) {
      max_lvl = M->cellTree->level(unref_cells[i]);
      cell_distrib.push_back(i);
      cell_lvl.push_back(max_lvl);
    }
  }

  cell_distrib.push_back(unref_cells.size());

  assert(face_distrib.size() == cell_distrib.size());

  // Assign unref pattern to faces if 2D mode
  std::vector<E_Int> patterns;
  assign_pattern_to_adapt_faces(unref_faces, patterns, M);

  for (size_t i = 0; i < face_distrib.size()-1; i++) {
    unrefine_cells(unref_cells, cell_distrib[i], cell_distrib[i+1], M);
    
    unrefine_faces(unref_faces, patterns, face_distrib[i],
      face_distrib[i+1], M);
    
    for (E_Int j = cell_distrib[i]; j < cell_distrib[i+1]; j++) {
      E_Int ucell = unref_cells[j];
      Children *cur = M->cellTree->children(ucell);
      Children *prev = cur->next;
      M->cellTree->children_[ucell] = prev;

      // free
      XFREE(cur);
    }
    
    for (E_Int j = face_distrib[i]; j < face_distrib[i+1]; j++) {
      E_Int uface = unref_faces[j];
      Children *cur = M->faceTree->children(uface);
      if (cur == NULL) {
        assert(M->owner[uface] == -1);
        assert(M->neigh[uface] == -1);
        continue;
      }
      Children *prev = cur->next;
      M->faceTree->children_[uface] = prev;

      // free
      XFREE(cur);
    }
  }

  update_boundary_faces(M);

  // Eliminate GONE cells
  E_Int nc = 0;
  E_Int sizeNFace = 0;
  Tree *CT = M->cellTree;
  std::vector<E_Int> new_cells(M->ncells, -1);
  for (E_Int i = 0; i < M->ncells; i++) {
    if (CT->state(i) != GONE) {
      new_cells[i] = nc++;
      sizeNFace += get_stride(i, M->indPH);
    }
  }

  // Eliminate GONE faces
  E_Int nf = 0;
  E_Int sizeNGon = 0;
  Tree *FT = M->faceTree;
  std::vector<E_Int> new_faces(M->nfaces, -1);
  for (E_Int i = 0; i < M->nfaces; i++) {
    if (FT->state(i) != GONE) {
      new_faces[i] = nf++;
      sizeNGon += get_stride(i, M->indPG);
    }
  }
}