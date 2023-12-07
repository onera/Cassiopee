#ifndef _STRUCT_H
#define _STRUCT_H

#include "xcore.h"
#include <map>
#include <mpi.h>

#define TETRA 0
#define PENTA 1
#define PYRA 2
#define HEXA 3

#define TRI 0
#define QUAD 1

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

extern const E_Int normalIn_T[4];
extern const E_Int normalIn_H[6];

struct pDirs {
  E_Float I[3];
  E_Float J[3];
  E_Float K[3];
  E_Float L[3];
};

struct Edge {
  E_Int p0_;
  E_Int p1_;

  Edge(E_Int p0, E_Int p1);

  void set(E_Int p0, E_Int p1);

  bool operator<(const Edge &e) const;
};

struct Patch {
  E_Int *faces;
  E_Int nfaces;
  E_Int nei_proc;
  E_Float *sbuf_d;
};

struct Element {
  E_Int *children;
  E_Int parent;
  uint8_t position;
  uint8_t type;
  uint8_t level;
};

struct AMesh {
  E_Int ncells;
  E_Int nfaces;
  E_Int npoints;

  E_Float *x;
  E_Float *y;
  E_Float *z;
  
  E_Int *owner;
  E_Int *neigh;
  
  E_Int *nface;
  E_Int *indPH;
  E_Int *ngon;
  E_Int *indPG;
  
  std::map<Edge, E_Int> ecenter;

  Patch *patches;
  E_Int npatches;

  int pid;
  int npc;
  int nreq;
  MPI_Request *req;

  Element **cellTree;
  Element **faceTree;

  E_Float *fc;
  E_Float *fa;
  E_Float *cx;
  E_Float *cy;
  E_Float *cz;
  
  E_Float **xn;
  E_Float **yn;
  E_Float **zn;
  E_Float **fn;

  E_Float *lsqG;
  E_Float *lsqGG;
  E_Float *lsqH;
  E_Float *lsqHH;

  short *ref_data;

  AMesh();
};


#endif
