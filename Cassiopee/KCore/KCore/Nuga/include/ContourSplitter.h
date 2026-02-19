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

#ifndef _CONTOURSPLITTER_H_
#define _CONTOURSPLITTER_H_

#include "Nuga/include/EltAlgo.h"


template <typename ElementType, typename BoundaryType>
class ContourSplitter
{
public:
  /// Splits the connectivity (only) with the edge N0N1 and stores the bits into separate containers.
  /** WARNING : it is the responsibility of the caller to destroy connectBout.*/
  static void splitConnectivity(const K_FLD::IntArray& connectBin, const std::set<BoundaryType>& cuttingEntities,
                     std::vector<K_FLD::IntArray> & connectBout);

  ///
  static void splitConnectivity(const K_FLD::IntArray& connectBin, const std::set<BoundaryType>& cuttingEntities,
                                std::vector< std::vector<E_Int> > & scolors);

  /// Splits and compacts the coordinates for each connectivity bits.
  static void splitCoordinates(std::vector<K_FLD::IntArray>& connects, const K_FLD::FloatArray& pos,
                               std::vector<K_FLD::FloatArray>& poss);

private:

  ContourSplitter(void){}

  ~ContourSplitter(void){}
};

#include "ContourSplitter.cxx"

#endif
