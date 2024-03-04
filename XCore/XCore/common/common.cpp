#include "common.h"

E_Int Get_pos(E_Int e, E_Int *pn, E_Int size)
{
  for (E_Int i = 0; i < size; i++) {
    if (pn[i] == e)
      return i;
  }
  return -1;
}

void Right_shift(E_Int *pn, E_Int pos, E_Int size)
{
  E_Int tmp[10];
  assert(size <= 10);
  for (E_Int i = 0; i < size; i++) tmp[i] = pn[i];
  for (E_Int i = 0; i < size; i++)
    pn[i] = tmp[(i+pos)%size];
}