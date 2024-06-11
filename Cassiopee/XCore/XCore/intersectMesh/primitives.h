#pragma once

#include <vector>

#include "queue.h"
#include "status.h"

extern E_Float TOL;

void compute_intersection(Queue &Q, Snode *sit0, Snode *sit1,
    std::vector<Vertex *> &I);

E_Int compare(const Vertex &a, const Vertex &b);

E_Int compare(const Segment &s0, const Segment &s1, E_Float rx, E_Float ry);

E_Int cmp_mySeg(const Segment &s1, const Segment &s2);

E_Int cmp_points(E_Float x1, E_Float y1, E_Float x2, E_Float y2);

E_Float DifferenceOfProducts(E_Float a, E_Float b, E_Float c, E_Float d);

E_Float TwoDiff(E_Float a, E_Float b);

E_Int Sign(E_Float x);

E_Int orient3D(E_Float *A, E_Float *B, E_Float *C, E_Float *D);

E_Float dRand(E_Float dMin, E_Float dMax);
