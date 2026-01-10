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

#ifndef NUGA_Q6_HXX
#define NUGA_Q6_HXX

#include "Nuga/include/subdivision.hxx"

namespace NUGA
{
  struct Q6
  {
    static void split(ngon_unit& PGs, const E_Int* nodes, NUGA::eDIR d, E_Int firstChild)
    {
      // set them in PGs
      E_Int* q41 = PGs.get_facets_ptr(firstChild);
      E_Int* q42 = PGs.get_facets_ptr(firstChild+1);

      if (d == Xd)
      {
        q41[0] = nodes[0]; q41[1] = nodes[4]; q41[2] = nodes[5]; q41[3] = nodes[3];
        q42[0] = nodes[4]; q42[1] = nodes[1]; q42[2] = nodes[2]; q42[3] = nodes[5];
      }
      else // d == Y
      {
        q41[0] = nodes[0]; q41[1] = nodes[1]; q41[2] = nodes[4]; q41[3] = nodes[5];
        q42[0] = nodes[5]; q42[1] = nodes[4]; q42[2] = nodes[2]; q42[3] = nodes[3];
      }
    }

    ///
    static void retrieve_ordered_data
    (const ngon_unit& PGs, E_Int PGi, E_Int i0, bool reorient, E_Int* two_childrenPG, E_Int* Q6nodes)
    {
      const E_Int* nodes = PGs.get_facets_ptr(PGi);
      E_Int nb_nodes = PGs.stride(PGi);

      for (int i = 0; i < nb_nodes; i++) Q6nodes[i] = nodes[i];

      const E_Int* pNFils0 = PGs.get_facets_ptr(two_childrenPG[0]);
      //const E_Int* pNFils1 = PGs.get_facets_ptr(two_childrenPG[1]);

      NUGA::eDIR dir = (pNFils0[1] == nodes[1]) ? Y : Xd;

      K_CONNECT::IdTool::right_shift<4>(Q6nodes, i0);

      if (reorient) std::swap(Q6nodes[1], Q6nodes[3]);

      //NUGA::eDIR dout;

      // 1. setting the last 2 points, the refining point fetched from the first PG child.
      // 2. eventually swapping the children to fit to the convention
      // 3. output the dir based on 1. the local dir (the stored one), 2. the orientation and 3. the first node id i0
      if (dir == Xd)
      {
        if (!reorient)
        {
          if (i0 == 0)
          {
            Q6nodes[4] = pNFils0[1];
            Q6nodes[5] = pNFils0[2];
            //dout = Xd;
          }
          else if (i0 == 1)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[1];
            //dout = Y;
            std::swap(two_childrenPG[0], two_childrenPG[1]);

          }
          else if (i0 == 2)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[1];
            //dout = Xd;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
          else if (i0 == 3)
          {
            Q6nodes[4] = pNFils0[1];
            Q6nodes[5] = pNFils0[2];
            //dout = Y;
          }
        }
        else // opposite orientation stored
        {
          if (i0 == 0)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[1];
            //dout = Y;
          }
          else if (i0 == 1)
          {
            Q6nodes[4] = pNFils0[1];
            Q6nodes[5] = pNFils0[2];
            //dout = Xd;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
          else if (i0 == 2)
          {
            Q6nodes[4] = pNFils0[1];
            Q6nodes[5] = pNFils0[2];
            //dout = Y;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
          else if (i0 == 3)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[1];
            //dout = Xd;
          }
        }
      }
      else // dir == Y
      {
        if (!reorient)
        {
          if (i0 == 0)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[3];
            //dout = Y;  
          }
          else if (i0 == 1)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[3];
            //dout = Xd;
          }
          else if (i0 == 2)
          {
            Q6nodes[4] = pNFils0[3];
            Q6nodes[5] = pNFils0[2];
            //dout = Y;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
          else if (i0 == 3)
          {
            Q6nodes[4] = pNFils0[3];
            Q6nodes[5] = pNFils0[2];
            //dout = Xd;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
        }
        else // opposite orientation stored
        {
          if (i0 == 0)
          {
            Q6nodes[4] = pNFils0[3];
            Q6nodes[5] = pNFils0[2];
            //dout = Xd;
          }
          else if (i0 == 1)
          {
            Q6nodes[4] = pNFils0[3];
            Q6nodes[5] = pNFils0[2];
            //dout = Y;
          }
          else if (i0 == 2)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[3];
            //dout = Xd;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
          else if (i0 == 3)
          {
            Q6nodes[4] = pNFils0[2];
            Q6nodes[5] = pNFils0[3];
            //dout = Y;
            std::swap(two_childrenPG[0], two_childrenPG[1]);
          }
        }
      }
    }
  };
}

#endif
