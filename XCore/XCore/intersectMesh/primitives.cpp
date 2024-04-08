#include "proto.h"

point::point(E_Float X, E_Float Y, size_t ID, size_t HE)
: x(X), y(Y), id(ID), he(HE)
{}

segment::segment(point *P, point *Q, size_t ID, size_t HE)
{
  p = P;
  q = Q;
  id = ID;
  dx = q->x - p->x;
  dy = q->y - p->y;
  he = HE;
  assert(segment_is_lexico(this));
}
