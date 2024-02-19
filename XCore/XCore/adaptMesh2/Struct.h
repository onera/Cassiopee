#ifndef _STRUCT_H
#define _STRUCT_H

#include "xcore.h"
#include <map>
#include <unordered_map>
#include <mpi.h>
#include "../common/common.h"

#define ISO 0
#define DIR 1

#define DIRX 0
#define DIRY 1

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
  E_Int nf;
  E_Int *pf;
  E_Int *pn;
  E_Int nei;
  E_Int *sbuf_i;
  E_Int *rbuf_i;
  E_Float *sbuf_f;
  E_Float *rbuf_f;
};

struct Children {
  E_Int n;
  E_Int pc[10];
  struct Children *next;
};

struct Tree {
  E_Int nelem_;
  E_Int *parent_;
  E_Int *level_;
  E_Int *type_;
  E_Int *state_;
  Children **children_;

  Tree(E_Int nelem);

  void drop();
  
  int check_numbering(E_Int N);

  void compress(const std::vector<E_Int> &new_ids, E_Int new_size); 

  void renumber(const std::vector<E_Int> &new_ids);

  void resize(E_Int new_nelem);

  void print_children();

  void print_face_types();

  void print_cell_types();

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

  /* Mesh */

  E_Int ncells;
  E_Int nfaces;
  E_Int npoints;
  E_Int nif;
  E_Int nbf;
  E_Int npf;

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

  /* Adaptation */
  
  E_Float Tr;
  E_Float Tu;
  E_Float eps;
  E_Float hmin;
  E_Float hmax;
  E_Int unrefine;
  E_Float *mode_2D;

  E_Int *ref_data;

  std::map<Edge, E_Int> *ecenter;

  Tree *cellTree;
  Tree *faceTree;

  E_Int prev_ncells;
  E_Int prev_nfaces;
  E_Int prev_npoints;

  E_Int onc;
  E_Int onf;
  E_Int onp;

  /* Parallel */

  int pid;
  int npc;
  int nrq;
  MPI_Request *req;

  E_Int *gcells;
  E_Int *gfaces;
  E_Int *gpoints;

  E_Int npatches;
  Patch *patches;

  std::unordered_map<E_Int, E_Int> *PT;
  std::unordered_map<E_Int, E_Int> *FT;
  std::unordered_map<E_Int, E_Int> *CT;

  E_Int *XNEIS;
  E_Int *CADJ;

  /* Export */

  E_Int *closed_indPG;
  E_Int *closed_ngon;
  E_Int *closed_indPH;
  E_Int *closed_nface;

  AMesh();
};


#endif
