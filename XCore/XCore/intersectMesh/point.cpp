#include "proto.h"

point::point()
{}

point::point(E_Float X, E_Float Y)
: x(X), y(Y), he(-1)
{}

int cmp_points(const point *p, const point *q)
{
  int s = geom_sign(p->x - q->x);
  return (s != 0) ? s : geom_sign(p->y - q->y);
}