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

#ifndef NUGA_PR18_HXX
#define NUGA_PR18_HXX

#include "splitting_t.hxx"

namespace NUGA
{

  template <>
  class splitting_t<K_MESH::Prism, NUGA::XYZ, 1> : public splitting_base_t
  {
  public:
    E_Int FACES[24]; // BOT00, BOTO1,...BOT03, TOP..., LEFT, ...RIGHT,...,FRONT...
    E_Int nodes[18];

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
      E_Int* TOP = FACES + 16;

      // Loop throught the four faces and create the conventional TH10 based on tet's faces subdivision

      const E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
      E_Int t6[6];

      // 1. BOT
      E_Int PGi = pPGi[0] - 1;
      const E_Int* pN = ng.PGs.get_facets_ptr(PGi);

      nodes[0] = *pN;
      nodes[1] = *(pN + 1);
      nodes[2] = *(pN + 2);

      // 1.1 convention : BOT is oriented inward
      bool reoriented(false);
      if (IS_OUTWARD(F2E, PGi, PHi))
      {
        std::swap(nodes[1], nodes[2]);
        reoriented = true;
      }

      // 1.2 retrieve BOT children according to the convention:
      // first node is the conventional first node
      const E_Int* pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)BOT[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, 0/*i0*/, reoriented, BOT, t6);

      nodes[6] = t6[3];
      nodes[7] = t6[4];
      nodes[8] = t6[5];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[0]);
      assert(t6[1] == nodes[1]);
      assert(t6[2] == nodes[2]);
#endif
      E_Int q9[9];

      // 2. F1
      PGi = pPGi[1] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);

      // 2.1 convention : F1 is oriented outward
      E_Int i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F1[k] = *(pchild++);//init

      Q9::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F1, q9);

      nodes[15] = q9[7];
      nodes[12] = q9[8];
      nodes[16] = q9[5];
      nodes[3] = q9[3];
      nodes[9] = q9[6];
      nodes[4] = q9[2];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(q9[0] == nodes[0]);
      assert(q9[4] == nodes[6]);
      assert(q9[1] == nodes[1]);
#endif

      // 3. F2
      PGi = pPGi[2] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);

      // 2.1 convention : F1 is oriented outward
      i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F2[k] = *(pchild++);//init

      Q9::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F2, q9);

      nodes[13] = q9[8];
      nodes[17] = q9[5];
      nodes[10] = q9[6];
      nodes[5] = q9[2];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(q9[0] == nodes[1]);
      assert(q9[4] == nodes[7]);
      assert(q9[1] == nodes[2]);
      assert(q9[7] == nodes[16]);
      assert(q9[3] == nodes[4]);
#endif

      // 4. F3
      PGi = pPGi[3] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);

      // 2.1 convention : F3 is oriented outward
      i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F3[k] = *(pchild++);//init

      Q9::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F3, q9);

      nodes[14] = q9[8];
      nodes[11] = q9[6];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(q9[0] == nodes[2]);
      assert(q9[4] == nodes[8]);
      assert(q9[1] == nodes[0]);
      assert(q9[7] == nodes[17]);
      assert(q9[5] == nodes[15]);
      assert(q9[3] == nodes[5]);
      assert(q9[2] == nodes[3]);
#endif

      // 2. TOP
      PGi = pPGi[4] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);

      // 2.1 convention : TOP is oriented outward
      i0 = K_CONNECT::IdTool::get_pos(pN, 3, nodes[3]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
        reoriented = true;

      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)TOP[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, TOP, t6);

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[3]);
      assert(t6[1] == nodes[4]);
      assert(t6[2] == nodes[5]);
      assert(t6[3] == nodes[9]);
      assert(t6[4] == nodes[10]);
      assert(t6[5] == nodes[11]);

      for (size_t i = 0; i < 18; ++i)
        assert(nodes[i] != IDX_NONE);
#endif
      
    }

    ///
    template <typename pg_arr_t, typename ph_arr_t>
    void split(ngon_type& ng, E_Int PHi, tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree, K_FLD::IntArray& F2E,
               E_Int firstIntPG, E_Int firstPHChild)
    { 
      const E_Int &nbc = subdiv_pol<K_MESH::Prism, ISO>::PHNBC;

      E_Int INT[10];
      for (size_t j = 0; j < 10; ++j)
        INT[j] = firstIntPG + j;

      ASSERT_IN_VECRANGE(ng.PGs, INT[0])
      E_Int* q41 = ng.PGs.get_facets_ptr(INT[0]);     // triangle

      ASSERT_IN_VECRANGE(ng.PGs, INT[1])
      E_Int* q42 = ng.PGs.get_facets_ptr(INT[1]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[2])
      E_Int* q43 = ng.PGs.get_facets_ptr(INT[2]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[3])
      E_Int* q44 = ng.PGs.get_facets_ptr(INT[3]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[4])
      E_Int* q45 = ng.PGs.get_facets_ptr(INT[4]);     // quad

      ASSERT_IN_VECRANGE(ng.PGs, INT[5])
      E_Int* q46 = ng.PGs.get_facets_ptr(INT[5]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[6])
      E_Int* q47 = ng.PGs.get_facets_ptr(INT[6]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[7])
      E_Int* q48 = ng.PGs.get_facets_ptr(INT[7]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[8])
      E_Int* q49 = ng.PGs.get_facets_ptr(INT[8]);

      ASSERT_IN_VECRANGE(ng.PGs, INT[9])
      E_Int* q410 = ng.PGs.get_facets_ptr(INT[9]);

      q41[0] = nodes[15]; q41[1] = nodes[12]; q41[2] = nodes[14];
      q42[0] = nodes[12]; q42[1] = nodes[16]; q42[2] = nodes[13];
      q43[0] = nodes[14]; q43[1] = nodes[13]; q43[2] = nodes[17];
      q44[0] = nodes[12]; q44[1] = nodes[13]; q44[2] = nodes[14];

      q45[0] = nodes[6]; q45[1] = nodes[8]; q45[2] = nodes[14]; q45[3] = nodes[12];  // �tage 1
      q46[0] = nodes[7]; q46[1] = nodes[6]; q46[2] = nodes[12]; q46[3] = nodes[13];
      q47[0] = nodes[8]; q47[1] = nodes[7]; q47[2] = nodes[13]; q47[3] = nodes[14];

      q48[0] = nodes[12]; q48[1] = nodes[14]; q48[2] = nodes[11]; q48[3] = nodes[9];  // �tage 2    
      q49[0] = nodes[13]; q49[1] = nodes[12]; q49[2] = nodes[9]; q49[3] = nodes[10];
      q410[0] = nodes[14]; q410[1] = nodes[13]; q410[2] = nodes[10]; q410[3] = nodes[11];

      // the 8 children of PH
      E_Int PHichildr[8];
      for (int j = 0; j < 8; ++j) {
        PHichildr[j] = firstPHChild + j;
      }

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[0])
      E_Int* h271 = ng.PHs.get_facets_ptr(PHichildr[0]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[1])
      E_Int* h272 = ng.PHs.get_facets_ptr(PHichildr[1]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[2])
      E_Int* h273 = ng.PHs.get_facets_ptr(PHichildr[2]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[3])
      E_Int* h274 = ng.PHs.get_facets_ptr(PHichildr[3]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[4])
      E_Int* h275 = ng.PHs.get_facets_ptr(PHichildr[4]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[5])
      E_Int* h276 = ng.PHs.get_facets_ptr(PHichildr[5]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[6])
      E_Int* h277 = ng.PHs.get_facets_ptr(PHichildr[6]);

      ASSERT_IN_VECRANGE(ng.PHs, PHichildr[7])
      E_Int* h278 = ng.PHs.get_facets_ptr(PHichildr[7]);

      splitPr18(INT, FACES, FACES + 4, FACES + 8, FACES + 12, FACES + 16, h271, h272, h273, h274, h275, h276, h277, h278);

      // set them in the tree
      PHtree.set_children(PHi, PHichildr, nbc);

      // update F2E
      update_F2E(ng, PHi, PHichildr, INT, PGtree, F2E);
    }
    ///
    static void splitPr18
    (const E_Int* INT, const E_Int* BOT, const E_Int* F1, const E_Int* F2, const E_Int* F3, const E_Int* TOP,
     E_Int* h271, E_Int* h272, E_Int* h273, E_Int* h274, E_Int* h275, E_Int* h276, E_Int* h277, E_Int* h278)
    {
      h271[0] = BOT[0]+1; h271[1] = F1[0]+1; h271[2] = INT[4]+1; h271[3] = F3[1]+1; h271[4] = INT[0]+1;
      h272[0] = BOT[1]+1; h272[1] = F1[1]+1; h272[2] = F2[0]+1; h272[3] = INT[5]+1; h272[4] = INT[1]+1;
      h273[0] = BOT[2]+1; h273[1] = INT[6]+1; h273[2] = F2[1]+1; h273[3] = F3[0]+1; h273[4] = INT[2]+1;
      h274[0] = INT[3]+1; h274[1] = INT[4]+1; h274[2] = INT[6]+1; h274[3] = INT[5]+1; h274[4] = BOT[3]+1;
 
      h275[0] = INT[0]+1; h275[1] = F1[3]+1; h275[2] = INT[7]+1; h275[3] = F3[2]+1; h275[4] = TOP[0]+1;
      h276[0] = INT[1]+1; h276[1] = F1[2]+1; h276[2] = F2[3]+1; h276[3] = INT[8]+1; h276[4] = TOP[1]+1;
      h277[0] = INT[2]+1; h277[1] = INT[9]+1; h277[2] = F2[2]+1; h277[3] = F3[3]+1; h277[4] = TOP[2]+1;
      h278[0] = INT[3]+1; h278[1] = INT[8]+1; h278[2] = INT[9]+1; h278[3] = INT[7]+1; h278[4] = TOP[3]+1; 
    } 

    ///
    template<typename arr_t>
    void update_F2E
    (const ngon_type& ng, E_Int PHi, const E_Int *children, const E_Int* INT, const tree<arr_t>& PGtree, K_FLD::IntArray & F2E)
    {
      //
      const E_Int &nbc = subdiv_pol<K_MESH::Prism, ISO>::PHNBC;
      splitting_base_t::__update_outer_F2E(ng, PHi, children, nbc, PGtree, F2E);

      // INTERNAL faces

      for (size_t k=0; k < 9; ++k)
      {
        assert(INT[k] > -1);
        assert(INT[k] < F2E.cols());
      }

      F2E(0, INT[0]) = *children;
      F2E(1, INT[0]) = *(children + 4);

      F2E(0, INT[1]) = *(children + 1);
      F2E(1, INT[1]) = *(children + 5);

      F2E(0, INT[2]) = *(children + 2);
      F2E(1, INT[2]) = *(children + 6);

      F2E(0, INT[3]) = *(children + 3);
      F2E(1, INT[3]) = *(children + 7);

      F2E(0, INT[4]) = *(children);
      F2E(1, INT[4]) = *(children + 3);

      F2E(0, INT[5]) = *(children + 1);
      F2E(1, INT[5]) = *(children + 3);

      F2E(0, INT[6]) = *(children + 2);
      F2E(1, INT[6]) = *(children + 3);

      F2E(0, INT[7]) = *(children + 4);
      F2E(1, INT[7]) = *(children + 7);

      F2E(0, INT[8]) = *(children + 5);
      F2E(1, INT[8]) = *(children + 7);

      F2E(0, INT[9]) = *(children + 6);
      F2E(1, INT[9]) = *(children + 7);
    }
        
   };

  using PR18 = splitting_t<K_MESH::Prism, NUGA::XYZ, 1>;
}

#endif
