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
