#include "proto.h"

int get_sign(E_Float x)
{
  if (x > TOL) return 1;
  if (x < -TOL) return -1;
  return 0;
}

int f_eq(E_Float x, E_Float y)
{
    return get_sign(fabs(x-y)) == 0;
}

int get_orient(point *a, point *b, point *c)
{
    E_Float det = (b->x - a->x)*(c->y - a->y) - (b->y - a->y)*(c->x - a->x);
    return get_sign(det);
}

int segments_are_colli(segment *s1, segment *s2)
{
    return (get_orient(s1->p, s1->q, s2->p) == 0) &&
           (get_orient(s1->p, s1->q, s2->q) == 0);
}

int cmp_xyz(E_Float x1, E_Float y1, E_Float x2, E_Float y2)
{
    int sign_x = get_sign(x1 - x2);
    if (sign_x != 0) return sign_x;
    int sign_y = get_sign(y1 - y2);
    return sign_y;
}

int cmp_points(const point *p1, const point *p2)
{
    int sign = cmp_xyz(p1->x, p1->y, p2->x, p2->y);
    assert(sign == 0 || sign == -1 || sign == 1);
    if (sign == 0) assert(p1 == p2);
    return sign;
}

int cmp_segments_lexico(const segment *s1, const segment *s2)
{
  int cmp = cmp_points(s1->p, s2->p);
  if (cmp != 0) return cmp;
  cmp = s1->id - s2->id;
  if (cmp < 0) return -1;
  if (cmp > 0) return 1;
  assert(s1 == s2);
  return 0;
}

int segment_is_lexico(const segment *s)
{
  return cmp_points(s->p, s->q) <= 0;
}
