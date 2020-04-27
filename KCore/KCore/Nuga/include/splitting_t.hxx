/*
 
 
 
              NUGA 
 
 
 
 */
#ifndef NUGA_SPLITTING_T_HXX
#define NUGA_SPLITTING_T_HXX

#include "subdivision.hxx"
 // CONVENTION UTILS //////////////////////////////////////////////////////////////////
inline bool need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R, const K_FLD::IntArray & F2E)
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
    template <typename arr_t>
    static void __update_outer_F2E(const ngon_type& ng, E_Int parentPHi, const E_Int* childrenPHi, E_Int nchildren, const tree<arr_t>& PGtree, K_FLD::IntArray& F2E);

  };

  template <typename ELT_t, eSUBDIV_TYPE STYPE, short ORDER>
  struct splitting_t;

  ///
  template <typename arr_t>
  void splitting_base_t::__update_outer_F2E
  (const ngon_type& ng, E_Int parentPHi, const E_Int* childrenPHi, E_Int nchildren, const tree<arr_t>& PGtree, K_FLD::IntArray& F2E)
  {
    E_Int count(0);

    for (E_Int c = 0; c < nchildren; ++c)
    {
      E_Int PHi = childrenPHi[c];
      const E_Int* faces = ng.PHs.get_facets_ptr(PHi);
      E_Int nfaces = ng.PHs.stride(PHi);

      for (E_Int i = 0; i < nfaces; ++i)
      {
        E_Int PGi = faces[i] - 1;
        E_Int parentPGi = PGtree.parent(PGi);

        if (parentPGi == E_IDX_NONE) continue; // inner face

        assert(PGi < F2E.cols());

        if (F2E(0, PGi) == E_IDX_NONE && F2E(1, PGi) == E_IDX_NONE) // init with parent data for outer faces
        {
          std::cout << "toto" << std::endl;
          F2E(0, PGi) = F2E(0, parentPGi);
          F2E(1, PGi) = F2E(1, parentPGi);
        }

        E_Int side = IS_OUTWARD(F2E, parentPGi, parentPHi) ? 0 : 1; //child PH is on same side of child PG as parent PH regarding parent PG 
        if (side == 1)
        {
          if (!IS_INWARD(F2E, parentPGi, parentPHi)) {
            std::cout << "PHi/parentPHi : " << PHi << "/" << parentPHi << std::endl;
            std::cout << "PGi/parentPGi : " << PGi << "/" << parentPGi << std::endl;
            std::cout << "F2E(0, PGi) : " << F2E(0, PGi) << std::endl;
            std::cout << "F2E(1, PGi) : " << F2E(1, PGi) << std::endl;
            std::cout << "F2E(0, parentPGi) : " << F2E(0, parentPGi) << std::endl;
            std::cout << "F2E(1, parentPGi) : " << F2E(1, parentPGi) << std::endl;
          }
        }

        F2E(side, PGi) = PHi;

        //fixme : to be really genric, the following must be recursive for the entire descendance of PGi
        E_Int nb_child_PGi = PGtree.nb_children(PGi);
        if (nb_child_PGi == 0) continue;

        const E_Int* childPGi = PGtree.children(PGi);
        for (E_Int n = 0; n < nb_child_PGi; ++n)
          F2E(side, childPGi[n]) = PHi;
      }

    }
  }

}

#endif
