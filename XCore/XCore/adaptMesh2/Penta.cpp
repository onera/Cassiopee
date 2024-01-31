#include "Proto.h"

void Pe18_refine(E_Int cell, AMesh *M)
{}
void Pe18_unrefine(E_Int cell, AMesh *M)
{}
void Pe12_refine(E_Int cell, AMesh *M)
{}
void Pe12_unrefine(E_Int cell, AMesh *M)
{}

void reorder_penta(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  // First tri is bottom
  E_Int bot = -1;
  for (E_Int i = 0; i < nf; i++) {
    if (get_stride(pf[i], M->indPG) == 3) {
      bot = pf[i];
      Right_shift(pf, i, 5);
      break;
    }
  }
  assert(bot != -1);

  E_Int common[3], map[3];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 3; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_Pe[0], M);
  if (reorient) std::swap(map[1], map[2]);

  E_Int top, lft, rgt, bck;
  top = lft = rgt = bck = -1;

  for (E_Int j = 1; j < 5; j++) {
    for (E_Int k = 0; k < 3; k++) common[k] = 0;

    E_Int face = pf[j];
    E_Int stride = get_stride(face, M->indPG);
    E_Int *pn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < stride; k++) {
      E_Int point = pn[k];

      for (E_Int l = 0; l < 3; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[2]) lft = face;
    else if (common[0] && common[1]) rgt = face;
    else if (common[1] && common[2]) bck = face;
    else                             top = face;
  }
  assert(top != -1);
  assert(lft != -1);
  assert(rgt != -1);
  assert(bck != -1);

  pf[1] = top;
  pf[2] = lft;
  pf[3] = rgt;
  pf[4] = bck;
}

E_Int check_canon_penta(E_Int cell, AMesh *M)
{
  E_Int NODES[6] = {-1, -1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_Pe[0], M);
  reorder_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT (in)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  assert(i0 != -1);
  reorient = get_reorient(face, cell, normalIn_Pe[2], M);
  reorder_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[2]);
  NODES[5] = local[2];
  NODES[3] = local[3];

  // RGT (out)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[3], M);
  reorder_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[3] == NODES[3]);
  NODES[4] = local[2];

  // BCK (in)
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[4], M);
  reorder_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);

  // TOP
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[3], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Pe[1], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[3]);
  assert(local[1] == NODES[4]);
  assert(local[2] == NODES[5]);

  return 1;
}