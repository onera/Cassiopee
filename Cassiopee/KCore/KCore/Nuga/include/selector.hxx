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

          if (K_SEARCH::BbTree<3>::boxesAreOverlapping(&ibox, &bxi, EPSILON))
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

  }

}

#endif
