#ifndef _INTERSECT_MESH_STRUCT_H
#define _INTERSECT_MESH_STRUCT_H

#include "../common/common.h"
#include "../common/mem.h"
#include "xcore.h"


#define INTERSECT_TOL 1e-6

struct Mesh {
  E_Int nc;
  E_Int nf;
  E_Int np;
  E_Int ne;

  E_Int *nface;
  E_Int *indPH;
  E_Int *ngon;
  E_Int *indPG;

  E_Float *x;
  E_Float *y;
  E_Float *z;

  E_Int *owner;
  E_Int *neigh;

  Mesh();
};

struct BBox {
  E_Float xmin;
  E_Float xmax;
  E_Float ymin;
  E_Float ymax;
  E_Float zmin;
  E_Float zmax;
  E_Float dmax;
};

struct Point {
  E_Float x;
  E_Float y;
  E_Float z;
};

struct Triangle {
  E_Int p0;
  E_Int p1;
  E_Int p2;
};

struct Edge_NO {
  E_Int p;
  E_Int q;

  Edge_NO(E_Int p, E_Int q);

  bool operator<(const Edge_NO &other) const;
};

struct Ray {
  E_Float p[3];
  E_Float q[3];
};

struct Edge_Hit {
  E_Float u;
  E_Float v;
  E_Float t;
  E_Int T; // Triangle id in master skin polyhedron

  E_Float x; // Coordinates of intersection points
  E_Float y;
  E_Float z;
};

#endif