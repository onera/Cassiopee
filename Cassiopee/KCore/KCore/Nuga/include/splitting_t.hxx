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

#ifndef NUGA_SPLITTING_T_HXX
#define NUGA_SPLITTING_T_HXX

#include "subdivision.hxx"

static bool need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R, const K_FLD::IntArray & F2E)
{
  if (F2E(1, PGi) == PHi && oriented_if_R == true) return false;
  else if (F2E(0, PGi) == PHi && oriented_if_R == false) return false;
  else return true;
}

#define IS_OUTWARD(F2E, PGi, PHi) (F2E(0,PGi)==PHi)
#define IS_INWARD(F2E, PGi, PHi) (F2E(1,PGi)==PHi)
#define IS_DIFF_FROM_2(x,a1,a2) ((x!=a1) && (x!=a2))
#define IS_DIFF_FROM_3(x, a1, a2, a3) ((x!=a1) && (x!=a2) && (x!=a3))

///////////////////////////////////////////////////////////////////////////////////////

namespace NUGA
{
  struct splitting_base_t
  {
    bool ok_for_split;

    splitting_base_t() : ok_for_split(true) {}

    template <typename arr_t>
    static void __update_outer_F2E(const ngon_type& ng, E_Int parentPHi, const E_Int* childrenPHi, E_Int nchildren, const tree<arr_t>& PGtree, K_FLD::IntArray& F2E);

    template <typename arr_t>
    static void __propagate_PH_neighbour_in_PG_descendance(E_Int PGi, E_Int side, E_Int PH, const tree<arr_t>& PGtree, K_FLD::IntArray& F2E);
  };

  template <typename ELT_t, NUGA::eDIR dir, short ORDER>
  class splitting_t;

  ///
  template <typename arr_t>
  void splitting_base_t::__update_outer_F2E
  (const ngon_type& ng, E_Int parentPHi, const E_Int* childrenPHi, E_Int nchildren, const tree<arr_t>& PGtree, K_FLD::IntArray& F2E)
  {
    //E_Int count(0);

    for (E_Int c = 0; c < nchildren; ++c)
    {
      E_Int PHi = childrenPHi[c];
      const E_Int* faces = ng.PHs.get_facets_ptr(PHi);
      E_Int nfaces = ng.PHs.stride(PHi);

      for (E_Int i = 0; i < nfaces; ++i)
      {
        E_Int PGi = faces[i] - 1;
        assert(PGi < F2E.cols());

        //E_Int parentPGi = PGtree.parent(PGi);

        if (F2E(0, PGi) == IDX_NONE && F2E(1, PGi) == IDX_NONE) // INNER FACE
          continue;

        E_Int side = IS_OUTWARD(F2E, PGi, parentPHi) ? 0 : 1; //child PH is on same side of child PG as parent PH regarding parent PG 
        /*if (side == 1)
        {
          if (!IS_INWARD(F2E, parentPGi, parentPHi)) {
            std::cout << "PHi/parentPHi : " << PHi << "/" << parentPHi << std::endl;
            std::cout << "PGi/parentPGi : " << PGi << "/" << parentPGi << std::endl;
            std::cout << "F2E(0, PGi) : " << F2E(0, PGi) << std::endl;
            std::cout << "F2E(1, PGi) : " << F2E(1, PGi) << std::endl;
            std::cout << "F2E(0, parentPGi) : " << F2E(0, parentPGi) << std::endl;
            std::cout << "F2E(1, parentPGi) : " << F2E(1, parentPGi) << std::endl;
          }
        }*/

        // replace parentPHi and propagate this PHi to PGi entire descendance
        __propagate_PH_neighbour_in_PG_descendance(PGi, side, PHi, PGtree, F2E);

      }
    }
  }

  ///
  template <typename arr_t>
  void splitting_base_t::__propagate_PH_neighbour_in_PG_descendance(E_Int PGi, E_Int side, E_Int PH, const tree<arr_t>& PGtree, K_FLD::IntArray& F2E)
  {
    F2E(side, PGi) = PH;

    E_Int nbc = PGtree.nb_children(PGi);
    if (nbc == 0) return;

    const E_Int* children = PGtree.children(PGi);
    for (E_Int n = 0; n < nbc; ++n)
      __propagate_PH_neighbour_in_PG_descendance(children[n], side, PH, PGtree, F2E);
  }

}

#endif
