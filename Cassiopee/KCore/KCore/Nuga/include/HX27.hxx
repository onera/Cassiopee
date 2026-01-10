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

#ifndef NUGA_H27_HXX
#define NUGA_H27_HXX

#include"splitting_t.hxx"
#include "Nuga/include/Hexahedron.h"

namespace NUGA
{
  template <>
  class splitting_t<K_MESH::Hexahedron, NUGA::XYZ, 1> : public splitting_base_t
  {
  public:
    E_Int FACES[24]; // BOT00, BOTO1,...BOT03, TOP..., LEFT, ...RIGHT,...,FRONT...
    E_Int nodes[27];

    ///
    template <typename arr_t>
    splitting_t(K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, E_Int centroidId, 
        const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree)
    {

      E_Int* BOT = FACES;
      E_Int* TOP = FACES + 4;
      E_Int* LEFT = FACES + 8;
      E_Int* RIGHT = FACES + 12;
      E_Int* FRONT = FACES + 16;
      E_Int* BACK = FACES + 20;

      const E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
      E_Int PGi = pPGi[0] - 1;
      const E_Int* pN = ng.PGs.get_facets_ptr(PGi);

      nodes[0] = *pN; // 0 -> PHi(0,0)

      if (F2E(1, PGi) == PHi) // for BOT, PH is the right element : well oriented
        for (int k = 1; k < 4; k++) nodes[k] = *(pN + k);

      else // otherwise : wrong orientation (swap 1 & 3)
      {
        nodes[1] = *(pN + 3);
        nodes[2] = *(pN + 2);
        nodes[3] = *(pN + 1);
      }
      // let's fill in nodes
      E_Int tmp[9];
      bool reorient;

      // BOT
      E_Int i0 = 0;
      reorient = need_a_reorient(pPGi[0] - 1, PHi, true, F2E);

      const E_Int* pchild = PGtree.children(pPGi[0] - 1);
      for (size_t k = 0; k < 4; ++k)BOT[k] = *(pchild++);//init

      /*if (reorient)*/ Q9::retrieve_ordered_data(ng.PGs, pPGi[0] - 1, i0, reorient, BOT, tmp);

      nodes[8] = tmp[4];
      nodes[9] = tmp[5];
      nodes[10] = tmp[6];
      nodes[11] = tmp[7];
      nodes[12] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmp[0] == nodes[0]);
      assert(tmp[1] == nodes[1]);
      assert(tmp[2] == nodes[2]);
      assert(tmp[3] == nodes[3]);
#endif

      // LEFT 
      const E_Int* p = ng.PGs.get_facets_ptr(pPGi[2] - 1);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
      pchild = PGtree.children(pPGi[2] - 1);
      for (size_t k = 0; k < 4; ++k)LEFT[k] = *(pchild++);//init
      reorient = need_a_reorient(pPGi[2] - 1, PHi, true, F2E);
      Q9::retrieve_ordered_data(ng.PGs, pPGi[2] - 1, i0, reorient, LEFT, tmp);

      nodes[21] = tmp[5];
      nodes[7] = tmp[2];
      nodes[16] = tmp[6];
      nodes[4] = tmp[3];
      nodes[18] = tmp[7];
      nodes[25] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmp[0] == nodes[0]);
      assert(tmp[1] == nodes[3]);
      assert(tmp[4] == nodes[11]);
#endif

      // RIGHT
      p = ng.PGs.get_facets_ptr(pPGi[3] - 1);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[1]);
      pchild = PGtree.children(pPGi[3] - 1);
      for (size_t k = 0; k < 4; ++k)RIGHT[k] = *(pchild++);//init

      reorient = need_a_reorient(pPGi[3] - 1, PHi, false, F2E);
      Q9::retrieve_ordered_data(ng.PGs, pPGi[3] - 1, i0, reorient, RIGHT, tmp);

      nodes[20] = tmp[5];
      nodes[6] = tmp[2];
      nodes[14] = tmp[6];
      nodes[5] = tmp[3];
      nodes[19] = tmp[7];
      nodes[23] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmp[0] == nodes[1]);
      assert(tmp[1] == nodes[2]);
      assert(tmp[4] == nodes[9]);
#endif

      // TOP
      p = ng.PGs.get_facets_ptr(pPGi[1] - 1);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[4]);
      pchild = PGtree.children(pPGi[1] - 1);
      for (size_t k = 0; k < 4; ++k)TOP[k] = *(pchild++);//init

      reorient = need_a_reorient(pPGi[1] - 1, PHi, false, F2E);
      Q9::retrieve_ordered_data(ng.PGs, pPGi[1] - 1, i0, reorient, TOP, tmp);

      nodes[13] = tmp[4];
      nodes[15] = tmp[6];
      nodes[17] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmp[0] == nodes[4]);
      assert(tmp[1] == nodes[5]);
      assert(tmp[2] == nodes[6]);
      assert(tmp[3] == nodes[7]);
      assert(tmp[5] == nodes[14]);
      assert(tmp[7] == nodes[16]);
#endif

      // FRONT
      p = ng.PGs.get_facets_ptr(pPGi[4] - 1);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[1]);
      pchild = PGtree.children(pPGi[4] - 1);
      for (size_t k = 0; k < 4; ++k)FRONT[k] = *(pchild++);//init

      reorient = need_a_reorient(pPGi[4] - 1, PHi, true, F2E);
      Q9::retrieve_ordered_data(ng.PGs, pPGi[4] - 1, i0, reorient, FRONT, tmp);

      nodes[22] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmp[0] == nodes[1]);
      assert(tmp[1] == nodes[0]);
      assert(tmp[2] == nodes[4]);
      assert(tmp[3] == nodes[5]);
      assert(tmp[4] == nodes[8]);
      assert(tmp[5] == nodes[18]);
      assert(tmp[6] == nodes[13]);
      assert(tmp[7] == nodes[19]);
#endif

      // BACK
      p = ng.PGs.get_facets_ptr(pPGi[5] - 1);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[2]);
      pchild = PGtree.children(pPGi[5] - 1);
      for (size_t k = 0; k < 4; ++k)BACK[k] = *(pchild++);//init

      reorient = need_a_reorient(pPGi[5] - 1, PHi, false, F2E);
      Q9::retrieve_ordered_data(ng.PGs, pPGi[5] - 1, i0, reorient, BACK, tmp);

      nodes[24] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmp[0] == nodes[2]);
      assert(tmp[1] == nodes[3]);
      assert(tmp[2] == nodes[7]);
      assert(tmp[3] == nodes[6]);
      assert(tmp[4] == nodes[10]);
      assert(tmp[5] == nodes[21]);
      assert(tmp[6] == nodes[15]);
      assert(tmp[7] == nodes[20]);
#endif

      K_MESH::Hexahedron::iso_barycenter(crd, nodes, 8, 1, crd.col(centroidId));
      //NUGA::refine_point_computer<K_MESH::Hexahedron>::compute_center(crd, nodes, 8, 1, crd.col(pos + i));

      nodes[26] = centroidId + 1;
    }

    ///
    template <typename pg_arr_t, typename ph_arr_t>
    void split(ngon_type& ng,E_Int PHi,  tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree, K_FLD::IntArray& F2E,
               E_Int firstIntPG, E_Int firstPHchild)
    {

      // INTERNAL faces creation
      E_Int p[27];
      p[0] = nodes[9];  p[1] = nodes[11];  p[2] = nodes[16];  p[3] = nodes[14];  p[4] = nodes[12];  p[5] = nodes[25];  p[6] = nodes[17];  p[7] = nodes[23];  p[8] = nodes[26];
      p[9] = nodes[18];  p[10] = nodes[19];  p[11] = nodes[20];  p[12] = nodes[21];  p[13] = nodes[22];  p[14] = nodes[23];  p[15] = nodes[24];  p[16] = nodes[25];  p[17] = nodes[26];
      p[18] = nodes[8];  p[19] = nodes[10];  p[20] = nodes[15];  p[21] = nodes[13];  p[22] = nodes[12];  p[23] = nodes[24];  p[24] = nodes[17];  p[25] = nodes[22];  p[26] = nodes[26];
      for (size_t j = 0; j < 3; ++j)
        NUGA::Q9::split<eSUBDIV_TYPE::ISO>(ng.PGs, &p[9 * j], firstIntPG + 4 * j);

      // HX8 creation
      E_Int *BOT(FACES), *TOP(FACES + 4), *LEFT(FACES + 8), *RIGHT(FACES + 12), *FRONT(FACES + 16), *BACK(FACES + 20);
      E_Int INT[12];
      for (size_t k = 0; k < 12; ++k)INT[k] = firstIntPG + k;

      E_Int* h271 = ng.PHs.get_facets_ptr(firstPHchild);
      E_Int* h272 = ng.PHs.get_facets_ptr(firstPHchild + 1);
      E_Int* h273 = ng.PHs.get_facets_ptr(firstPHchild + 2);
      E_Int* h274 = ng.PHs.get_facets_ptr(firstPHchild + 3);
      E_Int* h275 = ng.PHs.get_facets_ptr(firstPHchild + 4);
      E_Int* h276 = ng.PHs.get_facets_ptr(firstPHchild + 5);
      E_Int* h277 = ng.PHs.get_facets_ptr(firstPHchild + 6);
      E_Int* h278 = ng.PHs.get_facets_ptr(firstPHchild + 7);

      h271[0] = BOT[0] + 1; h271[1] = INT[4] + 1; h271[2] = LEFT[0] + 1; h271[3] = INT[8] + 1; h271[4] = FRONT[1] + 1; h271[5] = INT[1] + 1;
      h272[0] = BOT[1] + 1; h272[1] = INT[5] + 1; h272[2] = INT[8] + 1; h272[3] = RIGHT[0] + 1; h272[4] = FRONT[0] + 1; h272[5] = INT[0] + 1;
      h273[0] = BOT[2] + 1; h273[1] = INT[6] + 1; h273[2] = INT[9] + 1; h273[3] = RIGHT[1] + 1; h273[4] = INT[0] + 1; h273[5] = BACK[0] + 1;
      h274[0] = BOT[3] + 1; h274[1] = INT[7] + 1; h274[2] = LEFT[1] + 1; h274[3] = INT[9] + 1; h274[4] = INT[1] + 1; h274[5] = BACK[1] + 1;

      h275[0] = INT[4] + 1; h275[1] = TOP[0] + 1; h275[2] = LEFT[3] + 1; h275[3] = INT[11] + 1; h275[4] = FRONT[2] + 1; h275[5] = INT[2] + 1;
      h276[0] = INT[5] + 1; h276[1] = TOP[1] + 1; h276[2] = INT[11] + 1; h276[3] = RIGHT[3] + 1; h276[4] = FRONT[3] + 1; h276[5] = INT[3] + 1;
      h277[0] = INT[6] + 1; h277[1] = TOP[2] + 1; h277[2] = INT[10] + 1; h277[3] = RIGHT[2] + 1; h277[4] = INT[3] + 1; h277[5] = BACK[3] + 1;
      h278[0] = INT[7] + 1; h278[1] = TOP[3] + 1; h278[2] = LEFT[2] + 1; h278[3] = INT[10] + 1; h278[4] = INT[2] + 1; h278[5] = BACK[2] + 1;

      E_Int PHchildren[8];
      for (size_t k = 0; k < 8; ++k) PHchildren[k] = firstPHchild + k;

      // set them in the tree
      PHtree.set_children(PHi, PHchildren, 8);
      
      update_F2E(ng, PHi, PHchildren, INT, PGtree, F2E);
    }
   
    ///
    template <typename arr_t>
    void update_F2E
    (const ngon_type& ng, E_Int PHi, const E_Int *children, const E_Int* INT, const tree<arr_t>& PGtree, K_FLD::IntArray & F2E)
    {
      const E_Int &nbc = subdiv_pol<K_MESH::Hexahedron, ISO>::PHNBC;
      __update_outer_F2E(ng, PHi, children, nbc, PGtree, F2E);

      // INTERNAL faces
      F2E(0, INT[0]) = *(children + 1);
      F2E(1, INT[0]) = *(children + 2);

      F2E(0, INT[1]) = *(children);
      F2E(1, INT[1]) = *(children + 3);

      F2E(0, INT[2]) = *(children + 4);
      F2E(1, INT[2]) = *(children + 7);

      F2E(0, INT[3]) = *(children + 5);
      F2E(1, INT[3]) = *(children + 6);

      F2E(0, INT[4]) = *(children);
      F2E(1, INT[4]) = *(children + 4);

      F2E(0, INT[5]) = *(children + 1);
      F2E(1, INT[5]) = *(children + 5);

      F2E(0, INT[6]) = *(children + 2);
      F2E(1, INT[6]) = *(children + 6);

      F2E(0, INT[7]) = *(children + 3);
      F2E(1, INT[7]) = *(children + 7);

      F2E(0, INT[8]) = *(children);
      F2E(1, INT[8]) = *(children + 1);

      F2E(0, INT[9]) = *(children + 3);
      F2E(1, INT[9]) = *(children + 2);

      F2E(0, INT[10]) = *(children + 7);
      F2E(1, INT[10]) = *(children + 6);

      F2E(0, INT[11]) = *(children + 4);
      F2E(1, INT[11]) = *(children + 5);

    }
        
   };

   using HX27 = splitting_t<K_MESH::Hexahedron, NUGA::XYZ, 1>;
}

#endif
