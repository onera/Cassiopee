#ifndef _STRUCT_H
#define _STRUCT_H

#include "xcore.h"
#include <map>
#include <mpi.h>

#define TYPE_UNKNOWN -1
#define HEXA 0
#define TETRA 1
#define PENTA 2
#define PYRA 3

#define QUAD 0
#define TRI 1

#define MAXCHAR 1024

#define INTMAX E_IDX_NONE
#define INTMIN -(E_IDX_NONE-1)

#define RAISE(error) PyErr_SetString(PyExc_ValueError, (error))

#define GONE -2
#define UNREFINED -1
#define UNTOUCHED 0
#define REFINED 1

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

extern const E_Int normalIn_T[4];
extern const E_Int normalIn_H[6];
extern const E_Int normalIn_Pe[5];
extern const E_Int normalIn_Py[5];

struct pDirs {
  E_Float I[3];
  E_Float J[3];
  E_Float K[3];
  E_Float L[3];
};

struct Edge {
  E_Int p0_;
  E_Int p1_;

  Edge();

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

struct Children {
  E_Int n;
  E_Int pc[10];
  struct Children *next;
};

struct Tree {
  E_Int nelem_;
  E_Int *parent_;
  Children **children_;
  E_Int *level_;
  E_Int *type_;
  E_Int *state_;

  Tree(E_Int nelem);

  void compress(const std::vector<E_Int> &new_ids, E_Int new_size); 

  void renumber(const std::vector<E_Int> &new_ids);

  void resize(E_Int new_nelem);

  void print_children();

  void print_elem_histo(E_Int elem);

  void set_new_elem(E_Int elem, E_Int type, E_Int level);

  void set_parent_elem(E_Int elem, E_Int nchildren, E_Int pos);

  void set_child_elem(E_Int cpos, E_Int parent, E_Int type, E_Int level,
    E_Int pos);

  inline Children *children(E_Int elem)
  {
    return children_[elem];
  }

  inline E_Int level(E_Int elem)
  {
    return level_[elem];
  }

  inline E_Int type(E_Int elem)
  {
    return type_[elem];
  }

  inline E_Int parent(E_Int elem)
  {
    return parent_[elem];
  }

  inline E_Int state(E_Int elem)
  {
    return state_[elem];
  }
};


struct AMesh {
  E_Int ncells;
  E_Int nfaces;
  E_Int npoints;
  E_Int nif;
  E_Int nbf;

  E_Float *x;
  E_Float *y;
  E_Float *z;
   
  E_Int *nface;
  E_Int *indPH;
  E_Int *ngon;
  E_Int *indPG;

  E_Int *owner;
  E_Int *neigh;
  
  E_Int nbc;
  E_Int **ptlists;
  E_Int *bcsizes;
  char **bcnames;

  Patch *patches;
  E_Int npatches;

  int pid;
  int npc;
  int nreq;
  MPI_Request *req;

  std::map<Edge, E_Int> *ecenter;
  
  Tree *cellTree;
  Tree *faceTree;

  int *ref_data;
  E_Float Tr;
  E_Float Tu;

  E_Int *closed_indPG;
  E_Int *closed_ngon;

  E_Float *fc;
  E_Float *cx;
  E_Float *cy;
  E_Float *cz;

  AMesh();
};


#endif
