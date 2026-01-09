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

#ifndef NUGA_TH10_HXX
#define NUGA_TH10_HXX

#include "splitting_t.hxx"
#include "Nuga/include/Tetrahedron.h"


namespace NUGA
{
  template <>
  class splitting_t<K_MESH::Tetrahedron, NUGA::XYZ, 1> : public splitting_base_t
  {
  public:
    E_Int FACES[16]; // BOT00, BOTO1,...BOT03, F1..., F2 ..., F3...
    E_Int nodes[10];
    E_Int ndiag;

    template <typename arr_t>
    splitting_t(K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, E_Int centroidId,
        const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree)
    {
      // CONVENTION
      // BOT : INWARD
      // F1, F2, F3 : OUTWARD
      // N0 is first node of first face

      E_Int* BOT = FACES;
      E_Int* F1 = FACES + 4;
      E_Int* F2 = FACES + 8;
      E_Int* F3 = FACES + 12;

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

      nodes[4] = t6[3];
      nodes[5] = t6[4];
      nodes[6] = t6[5];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[0]);
      assert(t6[1] == nodes[1]);
      assert(t6[2] == nodes[2]);
      //assert(t6[3] == nodes[3]);
#endif

      // 2. F1
      PGi = pPGi[1] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);
  
      // 2.1 convention : F1 is oriented outward
      E_Int i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
      {
        reoriented = true;
      }
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F1[k] = *(pchild++);//init
      
      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F1, t6);

      nodes[7] = t6[5];
      nodes[8] = t6[4];
      nodes[3] = t6[2];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[0]);
      assert(t6[1] == nodes[1]);
      assert(t6[3] == nodes[4]);
#endif

      // 2. F2
      PGi = pPGi[2] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);

      // 2.1 convention : F2 is oriented outward
      i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
      {
        reoriented = true;
      }
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F2[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F2, t6);

      nodes[9] = t6[4];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[1]);
      assert(t6[1] == nodes[2]);
      assert(t6[2] == nodes[3]);
      assert(t6[5] == nodes[8]);
      assert(t6[3] == nodes[5]);
#endif

      // 3. F3
      PGi = pPGi[3] - 1;
      pN = ng.PGs.get_facets_ptr(PGi);

      // 2.1 convention : F3 is oriented outward
      i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
      reoriented = false;
      if (IS_INWARD(F2E, PGi, PHi))
      {
        reoriented = true;
      }
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)F3[k] = *(pchild++);//init

      T6::retrieve_ordered_data(ng.PGs, PGi, i0, reoriented, F3, t6);

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(t6[0] == nodes[2]);
      assert(t6[1] == nodes[0]);
      assert(t6[3] == nodes[6]);
      assert(t6[2] == nodes[3]);
      assert(t6[4] == nodes[7]);
      assert(t6[5] == nodes[9]);
#endif

      // diagonal choice
      E_Float d1 = NUGA::sqrDistance(crd.col(nodes[8] - 1), crd.col(nodes[6] - 1), 3);
      E_Float d2 = NUGA::sqrDistance(crd.col(nodes[7] - 1), crd.col(nodes[5] - 1), 3);
      E_Float d3 = NUGA::sqrDistance(crd.col(nodes[9] - 1), crd.col(nodes[4] - 1), 3);

      ndiag = ((d1 <= d2) && (d1 <= d3)) ? 1 : ((d2 <= d1) && (d2 <= d3)) ? 2 : 3;

    }

    template <typename pg_arr_t, typename ph_arr_t>
    void split(ngon_type& ng, E_Int PHi, tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree, K_FLD::IntArray& F2E,
               E_Int firsIntPG, E_Int firstPHChild)
    {
      static constexpr E_Int nbc = subdiv_pol<K_MESH::Tetrahedron, ISO>::PHNBC;
      //static constexpr E_Int nbi = subdiv_pol<K_MESH::Tetrahedron, ISO>::NBI;

      E_Int* q41 = ng.PGs.get_facets_ptr(firsIntPG);
      assert(ng.PGs.stride(firsIntPG) == 3);
      E_Int* q42 = ng.PGs.get_facets_ptr(firsIntPG+1);
      assert(ng.PGs.stride(firsIntPG+1) == 3);
      E_Int* q43 = ng.PGs.get_facets_ptr(firsIntPG+2);
      assert(ng.PGs.stride(firsIntPG+2) == 3);
      E_Int* q44 = ng.PGs.get_facets_ptr(firsIntPG+3);
      assert(ng.PGs.stride(firsIntPG+3) == 3);
      E_Int* q45 = ng.PGs.get_facets_ptr(firsIntPG+4);
      assert(ng.PGs.stride(firsIntPG+4) == 3);
      E_Int* q46 = ng.PGs.get_facets_ptr(firsIntPG+5);
      assert(ng.PGs.stride(firsIntPG+5) == 3);
      E_Int* q47 = ng.PGs.get_facets_ptr(firsIntPG+6);
      assert(ng.PGs.stride(firsIntPG+6) == 3);
      E_Int* q48 = ng.PGs.get_facets_ptr(firsIntPG+7);
      assert(ng.PGs.stride(firsIntPG+7) == 3);

      q41[0] = nodes[4]; q41[1] = nodes[6]; q41[2] = nodes[7];
      q42[0] = nodes[5]; q42[1] = nodes[4]; q42[2] = nodes[8];
      q43[0] = nodes[6]; q43[1] = nodes[5]; q43[2] = nodes[9];
      q44[0] = nodes[8]; q44[1] = nodes[7]; q44[2] = nodes[9];

      if (ndiag == 1) {
        q45[0] = nodes[6]; q45[1] = nodes[8]; q45[2] = nodes[5];
        q46[0] = nodes[7]; q46[1] = nodes[8]; q46[2] = nodes[6];
        q47[0] = nodes[4]; q47[1] = nodes[6]; q47[2] = nodes[8];
        q48[0] = nodes[8]; q48[1] = nodes[6]; q48[2] = nodes[9];
      }
      else if (ndiag == 2) {
        q45[0] = nodes[5]; q45[1] = nodes[7]; q45[2] = nodes[8];
        q46[0] = nodes[6]; q46[1] = nodes[7]; q46[2] = nodes[5];
        q47[0] = nodes[4]; q47[1] = nodes[7]; q47[2] = nodes[5];
        q48[0] = nodes[5]; q48[1] = nodes[7]; q48[2] = nodes[9];
      }
      else {
        q45[0] = nodes[4]; q45[1] = nodes[9]; q45[2] = nodes[6];
        q46[0] = nodes[8]; q46[1] = nodes[9]; q46[2] = nodes[4];
        q47[0] = nodes[5]; q47[1] = nodes[9]; q47[2] = nodes[4];
        q48[0] = nodes[4]; q48[1] = nodes[9]; q48[2] = nodes[7];
      }


      // the 8 children of PH
      E_Int PHchildren[8];
      for (size_t k = 0; k < 8; ++k) PHchildren[k] = firstPHChild + k;

      E_Int* h271 = ng.PHs.get_facets_ptr(PHchildren[0]);
      assert(ng.PHs.stride(PHchildren[0]) == 4);
      E_Int* h272 = ng.PHs.get_facets_ptr(PHchildren[1]);
      assert(ng.PHs.stride(PHchildren[1]) == 4);
      E_Int* h273 = ng.PHs.get_facets_ptr(PHchildren[2]);
      assert(ng.PHs.stride(PHchildren[2]) == 4);
      E_Int* h274 = ng.PHs.get_facets_ptr(PHchildren[3]);
      assert(ng.PHs.stride(PHchildren[3]) == 4);
      E_Int* h275 = ng.PHs.get_facets_ptr(PHchildren[4]);
      assert(ng.PHs.stride(PHchildren[4]) == 4);
      E_Int* h276 = ng.PHs.get_facets_ptr(PHchildren[5]);
      assert(ng.PHs.stride(PHchildren[5]) == 4);
      E_Int* h277 = ng.PHs.get_facets_ptr(PHchildren[6]);
      assert(ng.PHs.stride(PHchildren[6]) == 4);
      E_Int* h278 = ng.PHs.get_facets_ptr(PHchildren[7]);
      assert(ng.PHs.stride(PHchildren[7]) == 4);

      E_Int INT[8];
      for (size_t k = 0; k < 8; ++k)
        INT[k] = firsIntPG + k;

      splitT10(INT, FACES, FACES + 4, FACES + 8, FACES + 12, h271, h272, h273, h274, h275, h276, h277, h278, ndiag);

      // set them in the tree
      PHtree.set_children(PHi, firstPHChild, nbc);

      update_F2E(ng, PHi, PHchildren, INT, PGtree, F2E, ndiag);
      
    }
   
   
    ///
    static void splitT10
    (const E_Int* INT, const E_Int* BOT, const E_Int* F1, const E_Int* F2, const E_Int* F3,
     E_Int* h271, E_Int* h272, E_Int* h273, E_Int* h274, E_Int* h275, E_Int* h276, E_Int* h277, E_Int* h278, E_Int ndiag)
    {
  
      if (ndiag==1)
      {
        h271[0] = BOT[0] + 1; h271[1] = F1[0] + 1; h271[2] = INT[0] + 1; h271[3] = F3[1] + 1;
        h272[0] = BOT[1] + 1; h272[1] = F1[1] + 1; h272[2] = F2[0] + 1; h272[3] = INT[1] + 1;
        h273[0] = BOT[2] + 1; h273[1] = INT[2] + 1; h273[2] = F2[1] + 1; h273[3] = F3[0] + 1;
        h274[0] = INT[3] + 1; h274[1] = F2[2] + 1; h274[2] = F3[2] + 1; h274[3] = F1[2] + 1;
          
        h275[0] = INT[0] + 1; h275[1] = INT[6] + 1; h275[2] = INT[5] + 1; h275[3] = F1[3] + 1;
        h276[0] = INT[1] + 1; h276[1] = BOT[3] + 1; h276[2] = INT[6] + 1; h276[3] = INT[4] + 1;;
        h277[0] = INT[2] + 1; h277[1] = INT[4] + 1; h277[2] = F2[3] + 1; h277[3] = INT[7] + 1;
        h278[0] = INT[3] + 1; h278[1] = INT[5] + 1; h278[2] = F3[3] + 1; h278[3] = INT[7] + 1;
      }
      else if (ndiag==2)
      {
        h271[0] = BOT[0] + 1; h271[1] = F1[0] + 1; h271[2] = INT[0] + 1; h271[3] = F3[1] + 1;
        h272[0] = BOT[1] + 1; h272[1] = F1[1] + 1; h272[2] = F2[0] + 1; h272[3] = INT[1] + 1;
        h273[0] = BOT[2] + 1; h273[1] = INT[2] + 1; h273[2] = F2[1] + 1; h273[3] = F3[0] + 1;
        h274[0] = INT[3] + 1; h274[1] = F2[2] + 1; h274[2] = F3[2] + 1; h274[3] = F1[2] + 1;;
          
        h275[0] = INT[0] + 1; h275[1] = BOT[3] + 1; h275[2] = INT[5] + 1; h275[3] = INT[6] + 1;
        h276[0] = INT[1] + 1; h276[1] = INT[6] + 1; h276[2] = F1[3] + 1; h276[3] = INT[4] + 1;
        h277[0] = INT[2] + 1; h277[1] = INT[5] + 1; h277[2] = INT[7] + 1; h277[3] = F3[3] + 1;
        h278[0] = INT[3] + 1; h278[1] = INT[4] + 1; h278[2] = INT[7] + 1; h278[3] = F2[3] + 1;
      }
      else
      {
        h271[0] = BOT[0] + 1; h271[1] = F1[0] + 1; h271[2] = INT[0] + 1; h271[3] = F3[1] + 1;
        h272[0] = BOT[1] + 1; h272[1] = F1[1] + 1; h272[2] = F2[0] + 1; h272[3] = INT[1] + 1;
        h273[0] = BOT[2] + 1; h273[1] = INT[2] + 1; h273[2] = F2[1] + 1; h273[3] = F3[0] + 1;
        h274[0] = INT[3] + 1; h274[1] = F2[2] + 1; h274[2] = F3[2] + 1; h274[3] = F1[2] + 1;
          
        h275[0] = INT[0] + 1; h275[1] = INT[4] + 1; h275[2] = F3[3] + 1; h275[3] = INT[7] + 1;
        h276[0] = INT[1] + 1; h276[1] = INT[6] + 1; h276[2] = INT[5] + 1; h276[3] = F2[3] + 1;
        h277[0] = INT[2] + 1; h277[1] = BOT[3] + 1; h277[2] = INT[6] + 1; h277[3] = INT[4] + 1;
        h278[0] = INT[3] + 1; h278[1] = F1[3] + 1; h278[2] = INT[7] + 1; h278[3] = INT[5] + 1;
      }
    }

    template <typename arr_t>
    void update_F2E
    (const ngon_type& ng, E_Int PHi, const E_Int *children, const E_Int* INT, const tree<arr_t>& PGtree, K_FLD::IntArray & F2E, E_Int ndiag)
    {
      //
      const E_Int &nbc = subdiv_pol<K_MESH::Tetrahedron, ISO>::PHNBC;
      __update_outer_F2E(ng, PHi, children, nbc, PGtree, F2E);

      // INTERNAL faces
      if (ndiag == 1)
      {
        F2E(0, INT[0]) = *(children);
        F2E(1, INT[0]) = *(children + 4);

        F2E(0, INT[1]) = *(children + 1);
        F2E(1, INT[1]) = *(children + 5);

        F2E(0, INT[2]) = *(children + 2);
        F2E(1, INT[2]) = *(children + 6);

        F2E(0, INT[3]) = *(children + 3);
        F2E(1, INT[3]) = *(children + 7);

        F2E(0, INT[4]) = *(children + 5);
        F2E(1, INT[4]) = *(children + 6);

        F2E(0, INT[5]) = *(children + 4);
        F2E(1, INT[5]) = *(children + 7);

        F2E(0, INT[6]) = *(children + 4);
        F2E(1, INT[6]) = *(children + 5);

        F2E(0, INT[7]) = *(children + 7);
        F2E(1, INT[7]) = *(children + 6);
      }
      else if (ndiag == 2)
      {
        F2E(0, INT[0]) = *(children);
        F2E(1, INT[0]) = *(children + 4);

        F2E(0, INT[1]) = *(children + 1);
        F2E(1, INT[1]) = *(children + 5);

        F2E(0, INT[2]) = *(children + 2);
        F2E(1, INT[2]) = *(children + 6);

        F2E(0, INT[3]) = *(children + 3);
        F2E(1, INT[3]) = *(children + 7);

        F2E(0, INT[4]) = *(children + 5);
        F2E(1, INT[4]) = *(children + 7);

        F2E(0, INT[5]) = *(children + 4);
        F2E(1, INT[5]) = *(children + 6);

        F2E(0, INT[6]) = *(children + 5);
        F2E(1, INT[6]) = *(children + 4);

        F2E(0, INT[7]) = *(children + 7);
        F2E(1, INT[7]) = *(children + 6);
      }
      else //diag == 3
      {
        F2E(0, INT[0]) = *(children);
        F2E(1, INT[0]) = *(children + 4);

        F2E(0, INT[1]) = *(children + 1);
        F2E(1, INT[1]) = *(children + 5);

        F2E(0, INT[2]) = *(children + 2);
        F2E(1, INT[2]) = *(children + 6);

        F2E(0, INT[3]) = *(children + 3);
        F2E(1, INT[3]) = *(children + 7);

        F2E(0, INT[4]) = *(children + 6);
        F2E(1, INT[4]) = *(children + 4);

        F2E(0, INT[5]) = *(children + 5);
        F2E(1, INT[5]) = *(children + 7);

        F2E(0, INT[6]) = *(children + 6);
        F2E(1, INT[6]) = *(children + 5);

        F2E(0, INT[7]) = *(children + 4);
        F2E(1, INT[7]) = *(children + 7);
      }
    }

        
   };

   using TH10 = splitting_t<K_MESH::Tetrahedron, NUGA::XYZ, 1>;
}

#endif
