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

#ifndef NUGA_PY13_HXX
#define NUGA_PY13_HXX

#include "splitting_t.hxx"
#include "Nuga/include/Pyramid.h"

namespace NUGA
{

  template <>
  class splitting_t<K_MESH::Pyramid, NUGA::XYZ, 1> : public splitting_base_t
  {
  public:
    E_Int FACES[20]; // BOT F1,F2,F3,F4
    E_Int nodes[14];

    template <typename arr_t>
    splitting_t(K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, E_Int centroidId,
                const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree)
    {
      // CONVENTION
      // BOT : INWARD
      // OTHERS : OUTWARD
      // N0 is first node of first face

      E_Int* BOT = FACES;
      E_Int* F1 = FACES + 4;
      E_Int* F2 = FACES + 8;
      E_Int* F3 = FACES + 12;
      E_Int* F4 = FACES + 16;

      // Loop throught the four faces and create the conventional TH10 based on tet's faces subdivision

      const E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
      E_Int t6[6];

      // 1. BOT
      E_Int PGi = pPGi[0] - 1;
      const E_Int* pN = ng.PGs.get_facets_ptr(PGi);

      nodes[0] = *pN;
      nodes[1] = *(pN + 1);
      nodes[2] = *(pN + 2);
      nodes[3] = *(pN + 3);

      // 1.1 convention : BOT is oriented inward
      bool reoriented(false);
      if (IS_OUTWARD(F2E, PGi, PHi))
      {
        std::swap(nodes[1], nodes[3]);
        reoriented = true;
      }

      // 1.2 retrieve BOT children according to the convention:
      // first node is the conventional first node
      const E_Int* pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)BOT[k] = *(pchild++);//init
      
      E_Int q9[9];
      Q9::retrieve_ordered_data(ng.PGs, PGi, 0/*i0*/, reoriented, BOT, q9);

      nodes[5] = q9[4];
      nodes[6] = q9[5];
      nodes[7] = q9[6];
      nodes[8] = q9[7];
      nodes[9] = q9[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(q9[0] == nodes[0]);
      assert(q9[1] == nodes[1]);
      assert(q9[2] == nodes[2]);
      assert(q9[3] == nodes[3]);
#endif

      // 2. F1
      PGi = pPGi[1] - 1;
      const E_Int* p = ng.PGs.get_facets_ptr(PGi);
      E_Int i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);

      // 1.1 convention : F1 is oriented outward
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      // 1.2 retrieve F1 children according to the convention:
      // first node is the conventional first node
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F1[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F1, t6);

      nodes[4] = t6[2];
      nodes[10] = t6[4];
      nodes[13] = t6[5];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[0]);
      assert(t6[1] == nodes[1]);
      assert(t6[3] == nodes[5]);
#endif

      // 3. F2
      PGi = pPGi[2] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[1]);

      // 1.1 convention : F1 is oriented outward
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      // 1.2 retrieve F1 children according to the convention:
      // first node is the conventional first node
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F2[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F2, t6);

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[1]);
      assert(t6[1] == nodes[2]);
      assert(t6[2] == nodes[4]);
      assert(t6[3] == nodes[6]);
      assert(t6[5] == nodes[10]);
#endif

      // 4. F3
      PGi = pPGi[3] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[2]);

      // 1.1 convention : F3 is oriented outward
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      // 1.2 retrieve F1 children according to the convention:
      // first node is the conventional first node
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F3[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F3, t6);

      nodes[11] = t6[5];
      nodes[12] = t6[4];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[2]);
      assert(t6[1] == nodes[3]);
      assert(t6[2] == nodes[4]);
      assert(t6[3] == nodes[7]);
#endif

      // 5. F4
      PGi = pPGi[4] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      i0 = K_CONNECT::IdTool::get_pos(p, 3, nodes[3]);

      // 1.1 convention : F1 is oriented outward
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      // 1.2 retrieve F1 children according to the convention:
      // first node is the conventional first node
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F4[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F4, t6);
    }

    ///
    template <typename pg_arr_t, typename ph_arr_t>
    void split(ngon_type& ng, E_Int PHi, tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree, K_FLD::IntArray& F2E,
               E_Int firstIntPG, E_Int firstPHChild)
    { 
      E_Int INT[13];
      // INTERNAL faces  
      for (size_t j = 0; j < 13; ++j)
        INT[j] = firstIntPG + j;

      E_Int* q41 = ng.PGs.get_facets_ptr(INT[0]);
      E_Int* q42 = ng.PGs.get_facets_ptr(INT[1]);
      E_Int* q43 = ng.PGs.get_facets_ptr(INT[2]);
      E_Int* q44 = ng.PGs.get_facets_ptr(INT[3]);
      E_Int* q45 = ng.PGs.get_facets_ptr(INT[4]);
      E_Int* q46 = ng.PGs.get_facets_ptr(INT[5]);
      E_Int* q47 = ng.PGs.get_facets_ptr(INT[6]);
      E_Int* q48 = ng.PGs.get_facets_ptr(INT[7]);
      E_Int* q49 = ng.PGs.get_facets_ptr(INT[8]);
      E_Int* q410 = ng.PGs.get_facets_ptr(INT[9]);
      E_Int* q411 = ng.PGs.get_facets_ptr(INT[10]);
      E_Int* q412 = ng.PGs.get_facets_ptr(INT[11]);
      E_Int* q413 = ng.PGs.get_facets_ptr(INT[12]);

      q41[0] = nodes[5]; q41[1] = nodes[9]; q41[2] = nodes[13];
      q42[0] = nodes[9]; q42[1] = nodes[5]; q42[2] = nodes[10];
      q43[0] = nodes[6]; q43[1] = nodes[9]; q43[2] = nodes[10];
      q44[0] = nodes[9]; q44[1] = nodes[6]; q44[2] = nodes[11];
      q45[0] = nodes[7]; q45[1] = nodes[9]; q45[2] = nodes[11];
      q46[0] = nodes[9]; q46[1] = nodes[7]; q46[2] = nodes[12];
      q47[0] = nodes[8]; q47[1] = nodes[9]; q47[2] = nodes[12];
      q48[0] = nodes[9]; q48[1] = nodes[8]; q48[2] = nodes[13];

      q49[0] = nodes[10]; q49[1] = nodes[13]; q49[2] = nodes[9];
      q410[0] = nodes[11]; q410[1] = nodes[10]; q410[2] = nodes[9];
      q411[0] = nodes[12]; q411[1] = nodes[11]; q411[2] = nodes[9];
      q412[0] = nodes[13]; q412[1] = nodes[12]; q412[2] = nodes[9];

      q413[0] = nodes[13]; q413[1] = nodes[10]; q413[2] = nodes[11]; q413[3] = nodes[12];

      // the 8 children of PH
      E_Int PHichildr[10];

      for (int j = 0; j < 10; ++j)
        PHichildr[j] = firstPHChild + j;

      E_Int* h271 = ng.PHs.get_facets_ptr(PHichildr[0]);
      E_Int* h272 = ng.PHs.get_facets_ptr(PHichildr[1]);
      E_Int* h273 = ng.PHs.get_facets_ptr(PHichildr[2]);
      E_Int* h274 = ng.PHs.get_facets_ptr(PHichildr[3]);
      E_Int* h275 = ng.PHs.get_facets_ptr(PHichildr[4]);
      E_Int* h276 = ng.PHs.get_facets_ptr(PHichildr[5]);
      E_Int* h277 = ng.PHs.get_facets_ptr(PHichildr[6]);
      E_Int* h278 = ng.PHs.get_facets_ptr(PHichildr[7]);
      E_Int* h279 = ng.PHs.get_facets_ptr(PHichildr[8]);
      E_Int* h2710 = ng.PHs.get_facets_ptr(PHichildr[9]);

      splitP13(INT, FACES, FACES + 4, FACES + 8, FACES + 12, FACES + 16, h271, h272, h273, h274, h275, h276, h277, h278, h279, h2710);

      // set them in the tree
      PHtree.set_children(PHi, PHichildr, 10);

      // update F2E
      update_F2E(ng, PHi, PHichildr, INT, PGtree, F2E);
    }

    ///
    static void splitP13
    (const E_Int* INT, const E_Int* BOT, const E_Int* F1, const E_Int* F2, const E_Int* F3, const E_Int* F4,
      E_Int* h271, E_Int* h272, E_Int* h273, E_Int* h274, E_Int* h275, E_Int* h276, E_Int* h277, E_Int* h278, E_Int* h279, E_Int* h2710)
    {
      h271[0] = BOT[0] + 1; h271[1] = F1[0] + 1; h271[2] = INT[0] + 1; h271[3] = INT[7] + 1; h271[4] = F4[1] + 1;
      h272[0] = BOT[1] + 1; h272[1] = F1[1] + 1; h272[2] = F2[0] + 1; h272[3] = INT[2] + 1; h272[4] = INT[1] + 1;
      h273[0] = BOT[2] + 1; h273[1] = INT[3] + 1; h273[2] = F2[1] + 1; h273[3] = F3[0] + 1; h273[4] = INT[4] + 1;
      h274[0] = BOT[3] + 1; h274[1] = INT[6] + 1; h274[2] = INT[5] + 1; h274[3] = F3[1] + 1; h274[4] = F4[0] + 1;

      h277[0] = INT[8] + 1; h277[1] = F1[3] + 1; h277[2] = INT[0] + 1; h277[3] = INT[1] + 1;
      h278[0] = INT[9] + 1; h278[1] = F2[3] + 1; h278[2] = INT[2] + 1; h278[3] = INT[3] + 1;
      h279[0] = INT[10] + 1; h279[1] = F3[3] + 1; h279[2] = INT[4] + 1; h279[3] = INT[5] + 1;
      h2710[0] = INT[11] + 1; h2710[1] = F4[3] + 1; h2710[2] = INT[6] + 1; h2710[3] = INT[7] + 1;

      h275[0] = INT[12] + 1; h275[1] = F1[2] + 1; h275[2] = F2[2] + 1; h275[3] = F3[2] + 1; h275[4] = F4[2] + 1;
      h276[0] = INT[12] + 1; h276[1] = INT[11] + 1; h276[2] = INT[10] + 1; h276[3] = INT[9] + 1; h276[4] = INT[8] + 1;
    }
    

    ///
    template<typename arr_t>
    void update_F2E
    (const ngon_type& ng, E_Int PHi, const E_Int *children, const E_Int* INT, const tree<arr_t>& PGtree, K_FLD::IntArray & F2E)
    {
      //
      const E_Int &nbc = subdiv_pol<K_MESH::Pyramid, ISO>::PHNBC;
      splitting_base_t::__update_outer_F2E(ng, PHi, children, nbc, PGtree, F2E);

      // INTERNAL faces

      F2E(0, INT[0]) = *(children);
      F2E(1, INT[0]) = *(children + 6);

      F2E(0, INT[1]) = *(children + 1);
      F2E(1, INT[1]) = *(children + 6);

      F2E(0, INT[2]) = *(children + 1);
      F2E(1, INT[2]) = *(children + 7);

      F2E(0, INT[3]) = *(children + 2);
      F2E(1, INT[3]) = *(children + 7);

      F2E(0, INT[4]) = *(children + 2);
      F2E(1, INT[4]) = *(children + 8);

      F2E(0, INT[5]) = *(children + 3);
      F2E(1, INT[5]) = *(children + 8);

      F2E(0, INT[6]) = *(children + 3);
      F2E(1, INT[6]) = *(children + 9);

      F2E(0, INT[7]) = *(children);
      F2E(1, INT[7]) = *(children + 9);

      F2E(0, INT[8]) = *(children + 5);
      F2E(1, INT[8]) = *(children + 6);

      F2E(0, INT[9]) = *(children + 5);
      F2E(1, INT[9]) = *(children + 7);

      F2E(0, INT[10]) = *(children + 5);
      F2E(1, INT[10]) = *(children + 8);

      F2E(0, INT[11]) = *(children + 5);
      F2E(1, INT[11]) = *(children + 9);

      F2E(0, INT[12]) = *(children + 5);
      F2E(1, INT[12]) = *(children + 4);
    }
        
   };

  using PY13 = splitting_t<K_MESH::Pyramid, NUGA::XYZ, 1>;
}

#endif
