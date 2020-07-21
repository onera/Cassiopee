/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_SELECTOR_HXX
#define NUGA_SELECTOR_HXX

#include "Nuga/include/BbTree.h"

namespace NUGA
{

  namespace selector
  {
    ///
    template<typename mesh_t>
    static void reduce_to_box(mesh_t & m, K_SEARCH::BBox3D const & ibox, bool brute_force)
    {
      int nbcells = m.ncells();
      std::vector<bool> keep(nbcells, false);

      if (brute_force)
      {
        for (int i = 0; i< nbcells; ++i)
        {
          K_SEARCH::BBox3D bxi;
          m.bbox(i, bxi);

          if (K_SEARCH::BbTree<3>::boxesAreOverlapping(&ibox, &bxi, E_EPSILON))
            keep[i] = true;
        }
      }
      else
      {
         auto loc = m.get_localizer(); // build a bbtree
         std::vector<E_Int> cands;
         loc->get_tree()->getOverlappingBoxes(ibox.minB, ibox.maxB, cands);
         for (size_t i = 0; i < cands.size(); ++i) keep[cands[i]] = true;
      }

      m.compress(keep);
    }

  };

}

#endif
