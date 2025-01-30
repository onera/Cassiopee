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

#ifndef NUGA_H18_HXX
#define NUGA_H18_HXX

#include"splitting_t.hxx"
#include "Nuga/include/Hexahedron.h"

//#define DEBUG_HIERARCHICAL_MESH

namespace NUGA
{
  template <>
  class splitting_t<K_MESH::Hexahedron, NUGA::XY, 1> : public splitting_base_t
  {
  public:
    E_Int FACES[16]; // BOT00, BOTO1,...BOT03, TOP..., LEFT00, LEFT01, RIGHT00,...,FRONT...
    E_Int nodes[18];

    bool do_reorder;

    ///
    template <typename arr_t>
    splitting_t(K_FLD::FloatArray& crd, ngon_type& ng, E_Int PHi, E_Int centroidId,
      const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree, bool do_reord) : splitting_base_t(), do_reorder(do_reord)
    {

      E_Int* BOT = FACES;        //ISO
      E_Int* TOP = FACES + 4;    //ISO
      E_Int* LEFT = FACES + 8;   //DIR
      E_Int* RIGHT = FACES + 10; //DIR
      E_Int* FRONT = FACES + 12; //DIR
      E_Int* BACK = FACES + 14;  //DIR

#ifdef DEBUG_HIERARCHICAL_MESH
      if (PHi == 0)
      {
        const E_Int* faces = ng.PHs.get_facets_ptr(PHi);
        medith::write("BOT_0", crd, ng.PGs, faces[0]-1);
        medith::write("TOP_0", crd, ng.PGs, faces[1]-1);
        medith::write("LEFT_0", crd, ng.PGs, faces[2]-1);
        medith::write("RIGHT_0", crd, ng.PGs, faces[3]-1);
        medith::write("FRONT_0", crd, ng.PGs, faces[4]-1);
        medith::write("BACK_0", crd, ng.PGs, faces[5]-1);
      }
#endif

      if (do_reorder)
        reorder_as_XY(ng, PHi, PGtree, F2E); // swap faces to have bot/top ISO, remaining DIR (reorder here only for DIR_PROTO mode)
      
      const E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);

      if (PGtree.nb_children(pPGi[0] - 1) != 4 || PGtree.nb_children(pPGi[1] - 1) != 4 ||
          PGtree.nb_children(pPGi[2] - 1) != 2 || PGtree.nb_children(pPGi[3] - 1) != 2 ||
          PGtree.nb_children(pPGi[4] - 1) != 2 || PGtree.nb_children(pPGi[5] - 1) != 2) {
            ok_for_split = false;
            return;
          }

      // assert(PGtree.nb_children(pPGi[0] - 1) == 4);
      // assert(PGtree.nb_children(pPGi[1] - 1) == 4);
      // assert(PGtree.nb_children(pPGi[2] - 1) == 2);
      // assert(PGtree.nb_children(pPGi[3] - 1) == 2);
      // assert(PGtree.nb_children(pPGi[4] - 1) == 2);
      // assert(PGtree.nb_children(pPGi[5] - 1) == 2);

#ifdef DEBUG_HIERARCHICAL_MESH

      if (PHi == 127)
      {
        medith::write("elt", crd, ng, PHi);
        E_Int bot1 = PGtree.children(pPGi[0] - 1)[0];
        medith::write("BOT1", crd, ng.PGs, bot1);
        E_Int bot2 = PGtree.children(pPGi[0] - 1)[1];
        medith::write("BOT2", crd, ng.PGs, bot2);
        E_Int bot3 = PGtree.children(pPGi[0] - 1)[2];
        medith::write("BOT3", crd, ng.PGs, bot3);
        E_Int bot4 = PGtree.children(pPGi[0] - 1)[3];
        medith::write("BOT4", crd, ng.PGs, bot4);

        medith::write("TOP1", crd, ng.PGs, PGtree.children(pPGi[1] - 1)[0]);
        medith::write("TOP2", crd, ng.PGs, PGtree.children(pPGi[1] - 1)[1]);
        medith::write("TOP3", crd, ng.PGs, PGtree.children(pPGi[1] - 1)[2]);
        medith::write("TOP4", crd, ng.PGs, PGtree.children(pPGi[1] - 1)[3]);

        medith::write("LEFT1", crd, ng.PGs, PGtree.children(pPGi[2] - 1)[0]);
        medith::write("LEFT2", crd, ng.PGs, PGtree.children(pPGi[2] - 1)[1]);
        medith::write("RIGHT1", crd, ng.PGs, PGtree.children(pPGi[3] - 1)[0]);
        medith::write("RIGHT2", crd, ng.PGs, PGtree.children(pPGi[3] - 1)[1]);

        medith::write("FRONT1", crd, ng.PGs, PGtree.children(pPGi[4] - 1)[0]);
        medith::write("FRONT2", crd, ng.PGs, PGtree.children(pPGi[4] - 1)[1]);
        medith::write("BACK1", crd, ng.PGs, PGtree.children(pPGi[5] - 1)[0]);
        medith::write("BACK2", crd, ng.PGs, PGtree.children(pPGi[5] - 1)[1]);
    }
#endif

      E_Int PGi = pPGi[0]-1;
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
      E_Int tmpQ9[9];
      bool reorient;

      // BOT
      E_Int i0 = 0;
      reorient = need_a_reorient(PGi, PHi, true, F2E);

      const E_Int* pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)BOT[k] = *(pchild++);//init

      /*if (reorient)*/ Q9::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, BOT, tmpQ9);

      nodes[8] = tmpQ9[4];
      nodes[9] = tmpQ9[5];
      nodes[10] = tmpQ9[6];
      nodes[11] = tmpQ9[7];
      nodes[12] = tmpQ9[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmpQ9[0] == nodes[0]);
      assert(tmpQ9[1] == nodes[1]);
      assert(tmpQ9[2] == nodes[2]);
      assert(tmpQ9[3] == nodes[3]);
#endif

      // LEFT 
      E_Int tmpQ6[6];
      PGi = pPGi[2] - 1;
      const E_Int* p = ng.PGs.get_facets_ptr(PGi);
      //std::cout << p[0] << "/" << p[1] << "/" << p[2] << "/" << p[3] << std::endl;
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
      assert(i0 != -1);
      pchild = PGtree.children(PGi);
      assert(PGtree.nb_children(PGi) == 2);
      for (size_t k = 0; k < 2; ++k)LEFT[k] = *(pchild++);//init

      reorient = need_a_reorient(PGi, PHi, true, F2E);
      Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, LEFT, tmpQ6);

      nodes[4] = tmpQ6[3];
      nodes[7] = tmpQ6[2];
      nodes[16] = tmpQ6[5];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmpQ6[0] == nodes[0]);
      assert(tmpQ6[1] == nodes[3]);
      assert(tmpQ6[4] == nodes[11]);
#endif

      // RIGHT
      PGi = pPGi[3] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      //std::cout << p[0] << "/" << p[1] << "/" << p[2] << "/" << p[3] << std::endl;
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[1]);
      assert(i0 != -1);
      pchild = PGtree.children(PGi);
      assert(PGtree.nb_children(PGi) == 2);
      for (size_t k = 0; k < 2; ++k)RIGHT[k] = *(pchild++);//init

      reorient = need_a_reorient(PGi, PHi, false, F2E);
      Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, RIGHT, tmpQ6);

      nodes[5] = tmpQ6[3];
      nodes[6] = tmpQ6[2];
      nodes[14] = tmpQ6[5];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmpQ6[0] == nodes[1]);
      assert(tmpQ6[1] == nodes[2]);
      assert(tmpQ6[4] == nodes[9]);
#endif 

      // TOP
      PGi = pPGi[1] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[4]);
      assert(i0 != -1);
      pchild = PGtree.children(PGi);
      for (size_t k = 0; k < 4; ++k)TOP[k] = *(pchild++);//init

      reorient = need_a_reorient(PGi, PHi, false, F2E);
      Q9::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, TOP, tmpQ9);

      nodes[13] = tmpQ9[4];
      nodes[15] = tmpQ9[6];
      nodes[17] = tmpQ9[8];

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(nodes[14] == tmpQ9[5]);
      assert(nodes[16] == tmpQ9[7]);
      assert(nodes[5] == tmpQ9[1]);
      assert(nodes[6] == tmpQ9[2]);
      assert(nodes[4] == tmpQ9[0]);
      assert(nodes[7] == tmpQ9[3]);
#endif 

      // FRONT
      PGi = pPGi[4] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[1]);
      assert(i0 != -1);
      pchild = PGtree.children(PGi);
      assert(PGtree.nb_children(PGi) == 2);
      for (size_t k = 0; k < 2; ++k)FRONT[k] = *(pchild++);//init

      reorient = need_a_reorient(PGi, PHi, true, F2E);
      Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, FRONT, tmpQ6);

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmpQ6[0] == nodes[1]);
      assert(tmpQ6[1] == nodes[0]);
      assert(tmpQ6[2] == nodes[4]);
      assert(tmpQ6[3] == nodes[5]);
      assert(tmpQ6[4] == nodes[8]);
      assert(tmpQ6[5] == nodes[13]);
#endif

      // BACK
      PGi = pPGi[5] - 1;
      p = ng.PGs.get_facets_ptr(PGi);
      i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[2]);
      assert(i0 != -1);
      pchild = PGtree.children(PGi);
      assert(PGtree.nb_children(PGi) == 2);
      for (size_t k = 0; k < 2; ++k)BACK[k] = *(pchild++);//init

      reorient = need_a_reorient(PGi, PHi, false, F2E);
      Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, BACK, tmpQ6);

#ifdef DEBUG_HIERARCHICAL_MESH
      assert(tmpQ6[0] == nodes[2]);
      assert(tmpQ6[1] == nodes[3]);
      assert(tmpQ6[2] == nodes[7]);
      assert(tmpQ6[3] == nodes[6]);
      assert(tmpQ6[4] == nodes[10]);
      assert(tmpQ6[5] == nodes[15]);
#endif

    }

    template <typename arr_t>
    void reorder_as_XY(ngon_type& ng, E_Int PHi, const tree<arr_t> & PGtree, const K_FLD::IntArray& F2E)
    {
      E_Int* faces = ng.PHs.get_facets_ptr(PHi);
      E_Int nfaces = ng.PHs.stride(PHi);

      std::vector<E_Int> reord;
      reord.reserve(6);
      reord.insert(reord.end(), faces, faces + nfaces);

      E_Int nbc = PGtree.nb_children(faces[0] - 1);
      if (nbc == 4)
      {
        //swap enventually to be relevant with view regardin PH(0,0)
        /*const E_Int * nodes = ng.PGs.get_facets_ptr(reord[0] - 1);
        const E_Int* p = ng.PGs.get_facets_ptr(reord[2] - 1);
        E_Int i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
        if (i0 == -1)
        std::swap(reord[2], reord[3]);//left/right
        p = ng.PGs.get_facets_ptr(reord[4] - 1);
        i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
        if (i0 == -1)
        std::swap(reord[4], reord[5]);//front/back*/
      }

      nbc = PGtree.nb_children(faces[2] - 1);
      if (nbc == 4)
      {
        // left  -> bot
        reord[0] = faces[2];
        // right -> top
        reord[1] = faces[3];
        // top   -> left
        reord[2] = faces[1];
        // bot   -> right
        reord[3] = faces[0];

        //front/back unchanged
        reord[4] = faces[4];
        reord[5] = faces[5];

        for (size_t k = 0; k < 6; ++k)faces[k] = reord[k];
        K_MESH::Hexahedron::reorder_pgs(ng, F2E, PHi);

        //swap enventually to be relevant with view regardin PH(0,0)
        /*const E_Int * nodes = ng.PGs.get_facets_ptr(reord[0] - 1);
        const E_Int* p = ng.PGs.get_facets_ptr(reord[2] - 1);
        E_Int i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
        if (i0 == -1)
        std::swap(reord[2], reord[3]);//left/right
        p = ng.PGs.get_facets_ptr(reord[4] - 1);
        i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
        if (i0 == -1)
        std::swap(reord[4], reord[5]);//front/back*/

        return;
      }

      nbc = PGtree.nb_children(faces[4] - 1);
      if (nbc == 4)
      {
        // front -> bot
        reord[0] = faces[4];
        // back  -> top
        reord[1] = faces[5];
        // top   -> front
        reord[4] = faces[1];
        // bot   -> back
        reord[5] = faces[0];

        //left/right unchanged
        reord[2] = faces[2];
        reord[3] = faces[3];

        for (size_t k = 0; k < 6; ++k)faces[k] = reord[k];
        K_MESH::Hexahedron::reorder_pgs(ng, F2E, PHi);

        //swap enventually to be relevant with view regardin PH(0,0)
        /*const E_Int * nodes = ng.PGs.get_facets_ptr(reord[0] - 1);
        const E_Int* p = ng.PGs.get_facets_ptr(reord[2] - 1);
        E_Int i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
        if (i0 == -1)
        std::swap(reord[2], reord[3]);//left/right
        p = ng.PGs.get_facets_ptr(reord[4] - 1);
        i0 = K_CONNECT::IdTool::get_pos(p, 4, nodes[0]);
        if (i0 == -1)
        std::swap(reord[4], reord[5]);//front/back*/

        return;
      }
    }

    ///
    template <typename pg_arr_t, typename ph_arr_t>
    void split(ngon_type& ng,E_Int PHi,  tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree, K_FLD::IntArray& F2E,
               E_Int firstIntPG, E_Int firstPHchild)
    {

      // INTERNAL faces creation
      E_Int p[12];
      p[0] = nodes[9];  p[1] = nodes[11];  p[2] = nodes[16];  p[3] = nodes[14];  p[4] = nodes[12];  p[5] = nodes[17];
      p[6] = nodes[8];  p[7] = nodes[10];  p[8] = nodes[15];  p[9] = nodes[13];  p[10] = nodes[12];  p[11] = nodes[17];
      
      for (size_t j = 0; j < 2; ++j)
        NUGA::Q6::split(ng.PGs, &p[6 * j], Xd, firstIntPG + 2 * j);

      // HX8 creation
      E_Int *BOT(FACES), *TOP(FACES + 4), *LEFT(FACES + 8), *RIGHT(FACES + 10), *FRONT(FACES + 12), *BACK(FACES + 14);
      E_Int INT[4];
      for (size_t k = 0; k < 4; ++k)INT[k] = firstIntPG + k;

      E_Int* h181 = ng.PHs.get_facets_ptr(firstPHchild);
      E_Int* h182 = ng.PHs.get_facets_ptr(firstPHchild + 1);
      E_Int* h183 = ng.PHs.get_facets_ptr(firstPHchild + 2);
      E_Int* h184 = ng.PHs.get_facets_ptr(firstPHchild + 3);
     

      h181[0] = BOT[0] + 1; h181[1] = TOP[0] + 1; h181[2] = LEFT[0] + 1; h181[3] = INT[2] + 1  ; h181[4] = FRONT[1] + 1; h181[5] = INT[1] + 1;
      h182[0] = BOT[1] + 1; h182[1] = TOP[1] + 1; h182[2] = INT[2] + 1 ; h182[3] = RIGHT[0] + 1; h182[4] = FRONT[0] + 1; h182[5] = INT[0] + 1;
      h183[0] = BOT[2] + 1; h183[1] = TOP[2] + 1; h183[2] = INT[3] + 1 ; h183[3] = RIGHT[1] + 1; h183[4] = INT[0] + 1  ; h183[5] = BACK[0] + 1;
      h184[0] = BOT[3] + 1; h184[1] = TOP[3] + 1; h184[2] = LEFT[1] + 1; h184[3] = INT[3] + 1  ; h184[4] = INT[1] + 1  ; h184[5] = BACK[1] + 1;

      E_Int PHchildren[4];
      for (size_t k = 0; k < 4; ++k) PHchildren[k] = firstPHchild + k;

      // set them in the tree
      PHtree.set_children(PHi, PHchildren, 4);
      
      update_F2E(ng, PHi, PHchildren, INT, PGtree, F2E);
    }
   
    ///
    template <typename arr_t>
    void update_F2E
    (const ngon_type& ng, E_Int PHi, const E_Int *children, const E_Int* INT, const tree<arr_t>& PGtree, K_FLD::IntArray & F2E)
    {
      //const E_Int &nbc = subdiv_pol<K_MESH::Hexahedron, DIR>::PHNBC;
      __update_outer_F2E(ng, PHi, children, 4/*nbc*/, PGtree, F2E);

      // INTERNAL faces
      F2E(0, INT[0]) = *(children + 1);
      F2E(1, INT[0]) = *(children + 2);

      F2E(0, INT[1]) = *(children);
      F2E(1, INT[1]) = *(children + 3);

      F2E(0, INT[2]) = *(children);
      F2E(1, INT[2]) = *(children + 1);

      F2E(0, INT[3]) = *(children + 3);
      F2E(1, INT[3]) = *(children + 2);
    }
        
   };

   using HX18 = splitting_t<K_MESH::Hexahedron, NUGA::XY, 1>;
}

#endif
