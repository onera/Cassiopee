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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_T6_HXX
#define NUGA_T6_HXX

#include "Nuga/include/Polygon.h"
#include "Nuga/include/subdivision.hxx"

namespace NUGA
{

  struct T6 /*: public K_MESH::Polygon*/ //risky as currently refining nodes are not interleaved
  {
    
    template <eSUBDIV_TYPE STYPE>
    static void split(ngon_unit& PGs, const E_Int* nodes, E_Int firstChild)
    {
      // set them in _ng.PGs
      E_Int* p = PGs.get_facets_ptr(firstChild);

      p[0] = nodes[0];  p[1] = nodes[3];  p[2] = nodes[5]; 
      p[4] = nodes[3];  p[5] = nodes[1];  p[6] = nodes[4]; 
      p[8] = nodes[5];  p[9] = nodes[4];  p[10] = nodes[2]; 
      p[12] = nodes[3]; p[13] = nodes[4]; p[14] = nodes[5]; 
    }

    ///
    static void retrieve_ordered_data
    (const ngon_unit& PGs, E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* T6nodes)
    {
      const E_Int* pN = PGs.get_facets_ptr(PGi);
      E_Int nb_edges = PGs.stride(PGi);

      for (int i = 0; i < nb_edges; i++)
        T6nodes[i] = pN[i];
 
      const E_Int* pNFils3 = PGs.get_facets_ptr(four_childrenPG[3]);

#ifdef DEBUG_HIERARCHICAL_MESH    
      const E_Int* pNFils1 = PGs.get_facets_ptr(four_childrenPG[1]);
      const E_Int* pNFils2 = PGs.get_facets_ptr(four_childrenPG[2]);
      const E_Int* pNFils0 = PGs.get_facets_ptr(four_childrenPG[0]);

      assert(pNFils0[2] == pNFils2[0]);
      assert(pNFils1[0] == pNFils0[1]);
      assert(pNFils1[2] == pNFils2[1]);
      assert(pNFils2[0] == pNFils3[2]);
      assert(pNFils0[2] == pNFils3[2]);
#endif

      T6nodes[3] = pNFils3[0];
      T6nodes[4] = pNFils3[1];
      T6nodes[5] = pNFils3[2];

      K_CONNECT::IdTool::right_shift<3>(&T6nodes[0], i0);
      K_CONNECT::IdTool::right_shift<3>(&T6nodes[3], i0);
      K_CONNECT::IdTool::right_shift<3>(&four_childrenPG[0], i0);

      if (reorient == true)
      {
        std::swap(T6nodes[1], T6nodes[2]);
        std::swap(T6nodes[3], T6nodes[5]);
        std::swap(four_childrenPG[1], four_childrenPG[2]);
      }
    }
  };

}

#endif
