/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_VNGON_UNIT_H
#define NUGA_VNGON_UNIT_H

#include "Nuga/include/allocator.hxx"

namespace NUGA
{

  // for elements (polygons/polyhedra) storage
  template <typename INT_type, int ALLOC>// 0 (for CPP), 1 (for C)
  struct vngon_unit
  {
    static constexpr int CALLOC = ALLOC;
    using INT_t = INT_type;

    INT_t*   elts;      // contiguous storage of "facets" (i.e. faces or nodes) of elements (polyhedra or polygons)
    int      idx_start; // 0-based or 1-based indices in elts ?
    INT_t*   range;     // facets for i-th elements are in range [ range[i]; range[i+1]-1 ]
    INT_t    nrange;    // size of range array (equals the nb of elts + 1)

    // Default constructor
    vngon_unit() :elts(nullptr), idx_start(0), range(nullptr), nrange(0) {}

    // Standard constructor
    vngon_unit(INT_t* es, int idstart, INT_t* rng, INT_t nrng)
      :elts(es), idx_start(idstart), range(rng), nrange(nrng) {}

    vngon_unit(const vngon_unit& rhs) :elts(rhs.elts), idx_start(rhs.idx_start), range(rhs.range), nrange(rhs.nrange) {}

    void release() {
      NUGA::allocator<ALLOC>::deallocate(elts);
      NUGA::allocator<ALLOC>::deallocate(range);
      nrange = 0;
    }

  };

  // interleaved 3D-coordinates structure (x0 y0 z0...xi yi zi...)
  template <typename FLT_t, typename INT_t, int ALLOC>// 0 (for CPP), 1 (for C)
  struct crd3D
  {
    static constexpr int CALLOC = ALLOC;

    FLT_t * p;     // pointer to first x
    INT_t   n;     // nb of points

    crd3D():p(nullptr),n(0){}
    crd3D(FLT_t * pts, INT_t npts) :p(pts), n(npts) {}

    inline FLT_t* get(INT_t i) { return (p + 3 * i); }
    inline INT_t  getSize() { return n; }
    inline FLT_t* col(INT_t i) { return (p + 3 * i); }
    inline INT_t  cols() { return n; }
    void release() { NUGA::allocator<ALLOC>::deallocate(p);  n = 0; }

  };

  // polyhedral mesh
  template <typename FLT_type, typename INT_type, int ALLOC>
  struct phmesh
  {
    using INT_t = INT_type;
    using FLT_t = FLT_type;
    using ngu_t = vngon_unit<INT_t, ALLOC>;
    static constexpr int CALLOC = ALLOC;

    using crd_t = crd3D<FLT_t, INT_t, ALLOC>;

    phmesh() {}
    phmesh(crd_t& coords, ngu_t& pg, ngu_t& ph):pgs(pg), phs(ph), crd(coords){}

    vngon_unit<INT_t, ALLOC>  pgs;
    vngon_unit<INT_t, ALLOC>  phs;
    crd_t                     crd;
  };

  using c_morse_t = vngon_unit<int, 1>;
  using c_crd3D_t = crd3D<E_Float, int, 1>;
  using c_phmesh_t = phmesh<E_Float, int, 1>;

}

#endif
