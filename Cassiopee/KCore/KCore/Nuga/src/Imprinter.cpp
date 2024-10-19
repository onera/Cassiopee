/*    
    Copyright 2013-2024 Onera.

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

#include "Nuga/include/Imprinter.h"
#include "Nuga/include/macros.h"
#include <algorithm>

DELAUNAY::Imprinter::Imprinter
(const K_FLD::FloatArray& posS, const K_FLD::IntArray& connectT3):
_posS(posS), _connectT3(connectT3)
{
}

///
void
DELAUNAY::Imprinter::run(const K_FLD::FloatArray& posB0, const K_FLD::IntArray& connectB0,
                         K_FLD::FloatArray& posUV, int_vector_type& cell_indices, std::vector<E_Int>* hard_nodes)
{
  // Get the contout node indices;
  std::vector<E_Int>              B0_nodes;
  connectB0.uniqueVals(B0_nodes);

  if (hard_nodes)
    B0_nodes.insert(B0_nodes.end(), hard_nodes->begin(), hard_nodes->end());

  int_vector_type::const_iterator max = std::max_element(ALL(B0_nodes));

  cell_indices.clear();
  cell_indices.resize(*max + 1, IDX_NONE);

  posUV.clear();
  posUV.resize(2, *max + 1);

  E_Int nb_tris = _connectT3.cols();
  E_Int Ni;
  
  K_FLD::FloatArray::const_iterator pNi;
  
  E_Float UV[2];
  E_Float u,v, dmin, di;
  bool inside;

  for (size_t i = 0; i < B0_nodes.size(); ++i)
  {
    Ni  = B0_nodes[i];
    pNi = posB0.col(Ni);
    dmin = NUGA::FLOAT_MAX;

    for (size_type Si = 0; (Si < nb_tris); ++Si)
    {
      di = element_type::minDistanceToPoint(_posS, _connectT3.col(Si), pNi, UV, inside);
      u = UV[0];
      v = UV[1];

      if ((u == NUGA::FLOAT_MAX) || (v == NUGA::FLOAT_MAX))
        //return; // Numerical error : fixme
        continue;

      if (di < dmin)
      {
        dmin = di;
        posUV(0, Ni) = u;
        posUV(1, Ni) = v;
        cell_indices[Ni] = Si;
      }
    }
  }
}
