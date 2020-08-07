#ifndef NUGA_VNGON_UNIT_H
#define NUGA_VNGON_UNIT_H

// for elements (polygons/polyhedra) storage
template <typename INT_type, int ALLOC>// 0 (for CPP), 1 (for C)
struct vngon_unit
{
  static constexpr int CALLOC = ALLOC;
  using INT_t = INT_type;

  INT_t*   elts;   // contiguous storage of "facets" (face or node) of elements (polyhedra or polygons)
  INT_t    nelts;  // nb of elts
  int      idx_start;
  INT_t*   range;  // facets for i-th elements are in range [ range[i]; range[i+1]-1 ]
  INT_t    nrange; // size of range array

  vngon_unit():elts(nullptr), nelts(0), idx_start(0), range(nullptr), nrange(0){}

  void release() { 
    NUGA::allocator<ALLOC>::deallocate(elts);  nelts = 0; 
    NUGA::allocator<ALLOC>::deallocate(range); nrange = 0;
  }

};

// interleaved 3D-coordinates structure (x0 y0 z0...xi yi zi...)
template <typename FLT_t, typename INT_t, int ALLOC>// 0 (for CPP), 1 (for C)
struct crd3D
{
  static constexpr int CALLOC = ALLOC;

  FLT_t * p;     // pointer to first x
  INT_t   n;     // nb of points

  inline FLT_t* get(INT_t i) { return (p + 3 * i); }
  inline INT_t  getSize()    { return n; }
  inline FLT_t* col(INT_t i) { return (p + 3 * i); }
  inline INT_t  cols() { return n; }
};

// polyhedral mesh
template <typename FLT_type, typename INT_type, int ALLOC>
struct phmesh
{
  using INT_t = INT_type;
  using FLT_t = FLT_type;
  static constexpr int CALLOC = ALLOC;

  using crd_t = crd3D<FLT_t, INT_t, ALLOC>;

  vngon_unit<INT_t, ALLOC>  pgs;
  vngon_unit<INT_t, ALLOC>  phs;
  crd_t                     crd;
};

using c_morse_t = vngon_unit<int, 1>;
using c_crd3D_t = crd3D<double, int, 1>;
using c_phmesh_t = phmesh<double, int, 1>;

#endif