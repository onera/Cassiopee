#include "Proto.h"

void T6_refine(E_Int face, AMesh *)
{}

void T6_unrefine(E_Int face, AMesh *)
{}

void reorder_tri(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0)
{
  for (E_Int i = 0; i < 3; i++) local[i] = pn[i];
  Right_shift(local, i0, 3);
  if (reorient) std::swap(local[1], local[2]);
}