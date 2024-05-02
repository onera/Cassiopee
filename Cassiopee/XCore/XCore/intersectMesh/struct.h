#ifndef _INTERSECT_MESH_STRUCT_H
#define _INTERSECT_MESH_STRUCT_H

#include "../common/common.h"
#include "../common/mem.h"
#include "xcore.h"
#include <unordered_map>
#include <array>
#include <vector>

#define INTERSECT_TOL 1e-6
#define TOL 1e-10
#define RED 0
#define BLACK 1
#define INNER 0
#define OUTER 1
#define DEGEN -1
#define NO_IDEA 2

struct point {
  E_Float x, y, z;
};

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

struct vertex;
struct face;
struct cycle;

struct hedge {
  vertex *orig;
  hedge *twin;
  hedge *next;
  hedge *prev;
  face *left;
  int color;
  cycle *cycl;
  int s;

  hedge(vertex *);
};

struct face {
  hedge *rep;
  std::vector<hedge *> inner;
  int id;
  int ofp[2];

  face(); 
  void print_vertices();
  void print_parents();
};

struct vertex {
  E_Float x, y;
  E_Int nid;
  int oid[2]; // original point ids
  hedge *inc;
  hedge *left;

  vertex(E_Float, E_Float, int, bool=false);
  void print();
  void fprint(FILE *fh);
};

struct edge {
  vertex *p;
  vertex *q;
  int color;

  edge(vertex *P, vertex *Q);
};

struct segment {
  vertex *p;
  vertex *q;
  E_Int id;
  hedge *rep;

  segment(hedge *, E_Int);
  segment(edge *);
  segment(vertex *, vertex *, E_Int);
  inline int color() { return rep ? rep->color : NO_IDEA; }
};

/* SWEEP */
struct event {
  vertex *key;
  segment *inf;
  event *left;
  event *right;
  
  event(E_Float x, E_Float y, E_Int i);
  event(vertex *);
  void print();
  void inorder(std::vector<vertex *> &);
};

struct snode {
  segment *s;
  void *inf;
  snode *left;
  snode *right;

  snode(segment *);
  void print();
};

struct queue {
  event *root;
  int nelem;

  queue();

  void inorder(std::vector<vertex *> &);

  event *insert(vertex *);
  event *insert(E_Float, E_Float, E_Int);
  void erase(event *);
  void erase(vertex *);
  event *min();
  event *locate(E_Float, E_Float);
  vertex *locate_v(E_Float, E_Float);
  event *lookup(vertex *);
  int empty();

  event *_insert(event *&, vertex *);
  event *_insert(event *&, E_Float, E_Float, E_Int);
  event *_erase(event *, vertex *);
  event *_locate(event *, E_Float, E_Float);
  event *_lookup(event *, vertex *);

  void print();
};

struct status {
  snode *root;
  E_Float xs, ys;

  status();
  void print();

  snode *insert(segment *);
  void erase(segment *);
  snode *locate(segment *);
  snode *lookup(segment *);
  snode *pred(segment *);
  snode *pred(snode *);
  snode *succ(segment *);
  snode *succ(snode *);

  void update_sweep_position(E_Float, E_Float);

  snode *_insert(snode *&, segment *);
  snode *_erase(snode *, segment *);
  snode *_locate(snode *, segment *);
  snode *_lookup(snode *, segment *);
  void _pred(snode *, segment *, snode *&);
  void _succ(snode *, segment *, snode *&);

  int _segment_cmp(segment *, segment *);
  int _segment_cmp_sweep(segment *, segment *);
};

struct o_edge {
    E_Int p, q;

    o_edge(E_Int, E_Int);
};

struct smesh {
  std::vector<E_Float> X;
  std::vector<E_Float> Y;
  std::vector<std::vector<E_Int>> F;
  std::vector<o_edge> E;
  std::vector<std::vector<E_Int>> F2E;
  std::vector<std::array<E_Int, 2>> E2F;
  std::unordered_map<E_Int, E_Int> g2l_p; // mesh point id to local point id

  smesh(Mesh *, E_Int *, E_Int);
  void make_edges();
  void write_su2(const char *fname);
};

struct cycle {
  hedge *rep;
  int inout;
  vertex *left;
  cycle *prev;
  cycle *next;

  cycle(hedge *);
  void print();
  int get_size();
};

struct dcel {
  std::vector<vertex *> V;
  std::vector<hedge *> H;
  std::vector<face *> F;
  std::vector<cycle *> C;

  queue Q;

  dcel();
  dcel(const smesh &, const smesh &);
  void find_intersections();
  void set_color(int);
  void write_edges(const char *);
  void write_ngon(const char *);

  void _init_from_smesh(const smesh &, int color);
};

#endif
