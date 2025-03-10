/*    
    Copyright 2013-2025 Onera.

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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef _GEOMALGO_H_
#define _GEOMALGO_H_

#include "Nuga/include/DynArray.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/macros.h"

namespace NUGA {

template <typename ElementType>
class GeomAlgo {
public:
  ///
  template <typename Coordinate_t, typename Connectivity_t, short DIM>
  static void neighboring (const K_FLD::ArrayAccessor<Coordinate_t>& coords,
                           const K_FLD::ArrayAccessor<Connectivity_t>& conn, Connectivity_t& neighbors);
  ///
  inline static void reversi_chimera_skin (const Vector_t<K_FLD::FloatArray*>& crds, const Vector_t<K_FLD::IntArray*>& cnts, bool outward=true);
  
  /// Based on quality
  inline static void get_swapE (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, 
                                const std::set<K_MESH::NO_Edge>& hE, E_Float tol, std::vector<std::pair<E_Int, E_Int> >& swapE);
  /// Based on quality (with input har edges)
  //inline static void get_swapE (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, 
   //                             E_Float tol, std::vector<std::pair<E_Int, E_Int> >& swapE);
  
  ///
  inline static E_Float angle_measure(const E_Float* ni, const E_Float* nj, const E_Float* E0, const E_Float* E1);

  ///
  inline static void min_quality(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, E_Float& minq, E_Int& imin);
  
private:

  GeomAlgo(void){}

  ~GeomAlgo(void){}
};

} // End namespace NUGA

#include "GeomAlgo.cxx"

#endif
