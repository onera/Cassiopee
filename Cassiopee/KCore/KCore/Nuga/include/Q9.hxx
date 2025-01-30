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

#ifndef NUGA_Q9_HXX
#define NUGA_Q9_HXX

#include "Nuga/include/subdivision.hxx"

namespace NUGA
{
  struct Q9
  {
    template <eSUBDIV_TYPE STYPE>
    static void split(ngon_unit& PGs, const E_Int* nodes, E_Int firstChild)
    {
      // set them in PGs
      E_Int* q41 = PGs.get_facets_ptr(firstChild);
      E_Int* q42 = PGs.get_facets_ptr(firstChild+1);
      E_Int* q43 = PGs.get_facets_ptr(firstChild+2);
      E_Int* q44 = PGs.get_facets_ptr(firstChild+3);

      q41[0] = nodes[0]; q41[1] = nodes[4]; q41[2] = nodes[8]; q41[3] = nodes[7];
      q42[0] = nodes[4]; q42[1] = nodes[1]; q42[2] = nodes[5]; q42[3] = nodes[8];
      q43[0] = nodes[8]; q43[1] = nodes[5]; q43[2] = nodes[2]; q43[3] = nodes[6];
      q44[0] = nodes[7]; q44[1] = nodes[8]; q44[2] = nodes[6]; q44[3] = nodes[3];
    }

    ///
    static void retrieve_ordered_data
    (const ngon_unit& PGs, E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* Q9nodes)
    {
      const E_Int* nodes = PGs.get_facets_ptr(PGi);
      E_Int nb_nodes = PGs.stride(PGi);

      for (int i = 0; i < nb_nodes; i++) Q9nodes[i] = nodes[i];

      const E_Int* pNFils0 = PGs.get_facets_ptr(four_childrenPG[0]);
      const E_Int* pNFils2 = PGs.get_facets_ptr(four_childrenPG[2]);

#ifdef DEBUG_HIERARCHICAL_MESH    
      const E_Int* pNFils1 = PGs.get_facets_ptr(four_childrenPG[1]);
      const E_Int* pNFils3 = PGs.get_facets_ptr(four_childrenPG[3]);

      assert(pNFils0[2] == pNFils2[0]);
      assert(pNFils1[0] == pNFils0[1]);
      assert(pNFils1[2] == pNFils2[1]);
      assert(pNFils2[3] == pNFils3[2]);
      assert(pNFils0[3] == pNFils3[0]);
#endif

      // we got the 4 to 8 thanks to child 0 and child 2 (never swapped)
      Q9nodes[4] = pNFils0[1];
      Q9nodes[5] = pNFils2[1];
      Q9nodes[6] = pNFils2[3];
      Q9nodes[7] = pNFils0[3];
      Q9nodes[8] = pNFils0[2]; // centroid

      K_CONNECT::IdTool::right_shift<4>(&Q9nodes[0], i0);
      K_CONNECT::IdTool::right_shift<4>(&Q9nodes[4], i0);
      K_CONNECT::IdTool::right_shift<4>(&four_childrenPG[0], i0);

      if (reorient == true)
      {
        std::reverse(&Q9nodes[1], &Q9nodes[1] + 3);
        std::reverse(&Q9nodes[4], &Q9nodes[4] + 4);
        std::swap(four_childrenPG[1], four_childrenPG[3]);
      }

    }
  };
}

#endif
