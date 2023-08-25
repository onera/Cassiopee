#ifndef STRUCT_H
#define STRUCT_H

#include "intersector.h"
#include <mpi.h>
#include <unordered_map>
#include <map>

#define MAXNEI 6

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

extern const E_Int normalIn[6];
extern const double one_sixth;

struct proc_patch {
  E_Int *faces;
  E_Int nfaces;
  E_Int *gneis;
  E_Int nei_proc;
  E_Float *send_buf_d;
  E_Float *recv_buf_d;
  E_Int *send_buf_i;
  E_Int *recv_buf_i;
  E_Int *indices; // sorted in increasing gface value

  proc_patch()
  :
  faces(NULL), nfaces(-1), gneis(NULL), nei_proc(-1), send_buf_d(NULL), recv_buf_d(NULL),
  send_buf_i(NULL), recv_buf_i(NULL), indices(NULL)
  {}
};

struct edge {
  E_Int p0_;
  E_Int p1_;

  edge()
  {}

  edge(E_Int p0, E_Int p1)
  :
  p0_(std::min(p0, p1)),
  p1_(std::max(p0, p1))
  {}

  void set(E_Int p0, E_Int p1)
  {
    p0_ = std::min(p0, p1);
    p1_ = std::max(p0, p1);
  }

  bool operator<(const edge &a) const
  {
    return (p0_ < a.p0_) || (p0_ == a.p0_ && p1_ < a.p1_);
  }
};

struct mesh {
  E_Int ncells;
  E_Int nfaces;
  E_Int npoints;

  E_Float *xyz;

  E_Int *owner;
  E_Int *neigh;

  E_Int *NFACE;
  E_Int *xcells;
  E_Int *NGON;
  E_Int *xfaces;

  std::unordered_map<E_Int, E_Int> CT;
  std::unordered_map<E_Int, E_Int> FT;

  std::map<edge, E_Int> ET; // edge to center

  E_Float *cc;
  E_Float *fc;

  E_Float *lsqG; // matrix A for gradient computation
  E_Float *lsqGG; // matrix tA.A for gradient computation
  E_Float *lsqH; // matrix tA for hessian computation
  E_Float *lsqHH; // matrix tA.A for hessian computation

  proc_patch *ppatches;
  E_Int nppatches;

  E_Int pid;
  E_Int npc;
  E_Int nreq;
  MPI_Request *req;

  E_Float **pnei_coords;
  E_Float **pnei_flds;
  E_Float **pnei_grads;

  E_Int *gcells;
  E_Int *gfaces;

  mesh()
  :
  ncells(-1), nfaces(-1), npoints(-1), xyz(NULL), owner(NULL), neigh(NULL), NFACE(NULL),
  xcells(NULL), NGON(NULL), xfaces(NULL), CT(), FT(), ET(), cc(NULL), fc(NULL), lsqG(NULL),
  lsqGG(NULL), lsqH(NULL), lsqHH(NULL), ppatches(NULL), nppatches(-1), pid(-1), npc(-1), nreq(0),
  req(NULL), pnei_coords(NULL), pnei_flds(NULL), pnei_grads(NULL), gcells(NULL), gfaces(NULL)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &npc);
    req = (MPI_Request *)malloc(2*npc * sizeof(MPI_Request));
  }
};

struct tree {
  std::vector<E_Int> enabled;
  std::vector<E_Int> level;
  std::vector<E_Int> children;
  std::vector<E_Int> indir;
  std::vector<E_Int> parent;
  E_Int last; // fill next at children[last]
  E_Int size; // size of enabled, level, indir

  tree(E_Int, E_Int);
};

#endif
