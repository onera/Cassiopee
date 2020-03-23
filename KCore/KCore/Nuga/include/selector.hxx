/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_SELECTOR_HXX
#define NUGA_SELECTOR_HXX

#include "Search/BbTree.h"

namespace NUGA
{

  namespace selector
  {
    ///
    template<typename mesh_t>
    static void reduce_to_box(mesh_t & m, K_SEARCH::BBox3D const & ibox)
    {
      int nbcells = m.ncells();
      std::vector<bool> keep(nbcells, false);

      for (int i = 0; i< nbcells; ++i)
      {
        K_SEARCH::BBox3D bxi;
        m.bbox(i, bxi);

        if (K_SEARCH::BbTree<3>::boxesAreOverlapping(&ibox, &bxi, E_EPSILON))
          keep[i] = true;
      }

      m.compress(keep);
    }

  };

}

#endif
