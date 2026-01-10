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

#pragma once
#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/DefContainers.h"
#include "Nuga/include/KdTree.h"
#include <vector>

class PatchMaker
{
public:
  ///
  static void run (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectS, 
                   const NUGA::int_vector_type & pairs, E_Float angle_tolerance,
                   std::vector<K_FLD::IntArray> & connectBout);

private:

  ///
  static E_Int __update_normals(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
                                const std::vector<E_Int> &pairs, K_FLD::FloatArray& normals);

  ///
  static void
    __flag_critical_nodes(const K_FLD::FloatArray& normals, E_Float critical_angle,
                          const std::vector<E_Int> & pairs, const std::vector<E_Int> & colors,
                          std::vector< std::vector<E_Int> > & sorted_nodes,
                          std::vector<E_Int> & critical_nodes);

  ///
  static E_Int
    __flag_critical_nodes_on_contour(const std::vector<E_Int>& nodes, E_Int N0, const K_FLD::FloatArray& normals,
                                     E_Float critical_angle, const std::vector<E_Int> & pairs,
                                     const std::vector<E_Int>& colors, std::vector<E_Int> & critical_nodes);

  static E_Float __getAngle(const E_Float* n1, const E_Float* n2);

  static void __build_cutting_edges(const std::vector<E_Int>& pairs,
                                    const std::vector<E_Int>& colors,
                                    const std::vector<E_Int>& critical_nodes,
                                    K_FLD::IntArray& cuttingEdges,
                                    K_FLD::IntArray& periodicEdges);

  PatchMaker(void){}
  ~PatchMaker(void){}
};
