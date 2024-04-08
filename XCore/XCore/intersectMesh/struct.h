#ifndef _INTERSECT_MESH_STRUCT_H
#define _INTERSECT_MESH_STRUCT_H

#include "../common/common.h"
#include "../common/mem.h"
#include "xcore.h"

#define INTERSECT_TOL 1e-6
#define TOL 1e-10
#define LEFT -1
#define RIGHT 1
#define ON 0

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

struct Triangle {
  E_Int A;
  E_Int B;
  E_Int C;
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
  E_Float w;

  E_Float t;

  E_Int F; // Quad intersected
  E_Int T; // First or second triangle of Q

  E_Float x; // Coordinates of intersection points
  E_Float y;
  E_Float z;

  E_Int is_vertex;
  E_Int on_edge;

  Edge_Hit() :
    u(-1), v(-1), w(-1), t(-1), F(-1), T(-1), x(0), y(0), z(0),
    is_vertex(-1), on_edge(-1)
  {}
};

/********************************************/

/* DCEL */

struct point {
  E_Float x;
  E_Float y;
  size_t id;
  size_t he;

  point(E_Float x, E_Float y, size_t id, size_t he);
};

struct half_edge {
  size_t orig;
  size_t twin;
  size_t next;
  size_t prev;
  size_t face;
};

struct face {
  size_t outer_component;
  std::vector<E_Int> inner_components;
};

struct dcel {
  std::vector<point *> points;
  std::vector<half_edge> half_edges;
  std::vector<face> faces;

  dcel(Mesh *);
  dcel(E_Int ni, E_Int nj, const std::vector<E_Float> &X,
    const std::vector<E_Float> &Y, const std::vector<E_Float> &Z);
  dcel(const dcel &, const dcel &);

  void resolve();
  int is_valid();
  int _check_duplicate_points();
};

struct segment {
  point *p;
  point *q;
  size_t id;
  E_Float dx;
  E_Float dy;
  size_t he; // corresponding half-edge

  segment(point *, point *, size_t, size_t);
};

struct event {
    point *p;
    segment *s;
    event *left;
    event *right;

    event(point *);
    void print();
};

struct status {
    segment *s;
    point *p;
    segment *c;
    status *left;
    status *right;

    status(segment *);
    void print();
};

#endif
