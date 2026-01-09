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

#include "Nuga/include/TSSplitter.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/ContourSplitter.h"

// Splits the connectivity (only) with the edge N0N1 and stores the bits 
// into separate containers.
void TSSplitter::split
(const K_FLD::IntArray& BAR, const K_FLD::IntArray& polyLine, std::vector<K_FLD::IntArray> & connectBout)
{
  std::set<K_MESH::NO_Edge> cuttingEdges;
  for (E_Int Si = 0; Si < polyLine.cols(); ++Si)
    cuttingEdges.insert(polyLine.col(Si));

  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(BAR, cuttingEdges, connectBout);
}

// Splits the connectivity with the edge N0N1 and stores the connectivity 
// bits and corresponding coordinates into separate containers.
void TSSplitter::split
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectBin, const K_FLD::IntArray& polyLine,
 std::vector<K_FLD::FloatArray> & posOut, std::vector<K_FLD::IntArray> & connectBout)
{
  split(connectBin, polyLine, connectBout);

  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitCoordinates(connectBout, pos, posOut);
}
