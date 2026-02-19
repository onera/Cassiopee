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

#ifndef _CONTOURSPLITTER_CXX_
#define _CONTOURSPLITTER_CXX_

#include "Nuga/include/ContourSplitter.h"

///
template <typename ElementType, typename BoundaryType>
void
ContourSplitter<ElementType, BoundaryType>::splitConnectivity
(const K_FLD::IntArray& connectBin, const std::set<BoundaryType>& cuttingEntities,
 std::vector<K_FLD::IntArray> & connectBout)
{
  typename NUGA::EltAlgo<ElementType>::NeighbourType    neighbors;
  std::vector<E_Int>                                              colors;
  E_Int                                                           maxcol, NELT(connectBin.cols());
  K_FLD::IntArray::const_iterator                                 pS;

  connectBout.clear();

  if (connectBin.cols() == 0)
    return;
  
  // Get the neighbouring of each element.
  NUGA::EltAlgo<ElementType>::getNeighbours(connectBin, neighbors);

  if (neighbors.empty())
    return;

  // Color the elements.
  NUGA::EltAlgo<ElementType>::coloring(connectBin, neighbors, cuttingEntities, colors);

  // Reserve memory to avoid reallocation.
  maxcol = *std::max_element(colors.begin(), colors.end());
  connectBout.resize(maxcol+1);
  for (size_t k = 0; k < connectBout.size(); ++k)
    connectBout[k].reserve(ElementType::NB_NODES, NELT);

  // Split by color.
  for (E_Int c = 0; c < NELT; ++c)
  {
    pS = connectBin.col(c);
    connectBout[colors[c]].pushBack(pS, pS + ElementType::NB_NODES);
  }
}

///
template <typename ElementType, typename BoundaryType>
void
ContourSplitter<ElementType, BoundaryType>::splitConnectivity
(const K_FLD::IntArray& connectBin, const std::set<BoundaryType>& cuttingEntities,
 std::vector< std::vector<E_Int> > & scolors)
{
  typename NUGA::EltAlgo<ElementType>::NeighbourType  neighbors;
  std::vector<E_Int>                                                     colors;
  E_Int                                                                  maxcol, NELT(connectBin.cols());
  //K_FLD::IntArray::const_iterator                                        pS;

  scolors.clear();

  if (connectBin.cols() == 0)
    return;
  
  // Get the neighbouring of each element.
  NUGA::EltAlgo<ElementType>::getNeighbours(connectBin, neighbors);

  if (neighbors.empty())
    return;

  // Color the elements.
  NUGA::EltAlgo<ElementType>::coloring(connectBin, neighbors, cuttingEntities, colors);

  // Reserve memory to avoid reallocation.
  maxcol = *std::max_element(colors.begin(), colors.end());
  scolors.resize(maxcol+1);
  
  // Split by color.
  for (E_Int c = 0; c < NELT; ++c)
  {
    //pS = connectBin.col(c);
    scolors[colors[c]].push_back(c);
  }
}

///
template <typename ElementType, typename BoundaryType>
void
ContourSplitter<ElementType, BoundaryType>::splitCoordinates
(std::vector<K_FLD::IntArray>& connects, const K_FLD::FloatArray& pos, std::vector<K_FLD::FloatArray>& poss)
{
  std::vector<bool>      flag0(pos.cols(), false), flag;
  std::vector<E_Int>     indices, newIds;
  K_FLD::FloatArray      Posi;
  size_t                 nbConnects(connects.size());                  

  // Fast return
  if (nbConnects < 1)    return; // Error


  poss.reserve(nbConnects);
 
  for (size_t i = 0; i < nbConnects; ++i)
  {
    flag = flag0;
    Posi = pos;

    if (nbConnects > 1)
    {
      connects[i].uniqueVals(indices);

      for (size_t j = 0; j < indices.size(); ++j)           // Flag the mesh nodes.
        flag[indices[j]] = true;

      K_FLD::FloatArray::compact (Posi, flag, newIds);     // Remove any node not in the mesh.
      K_FLD::IntArray::changeIndices(connects[i], newIds); // Change mesh indices accordingly.
    }

    poss.push_back(Posi);
    //connects[i].resize(connects[i].rows(), connects[i].cols());
  }
}

#endif
