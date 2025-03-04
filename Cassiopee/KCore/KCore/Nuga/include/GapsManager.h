/*    
    Copyright 2013-2025 Onera.

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
#include "Nuga/include/Triangle.h"
#include <vector>

class NodeAssociator;

class GapsManager
{
public:
  enum eCase {POST_NODAL = 0, POST_CENTER = 1, UNDEFINED = 2};
  enum eMode {SURFACIC = 0, PLANAR = 1};
public:
  ///
  static void run (const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components,
                   std::vector<K_FLD::FloatArray>& posFs, std::vector<K_FLD::IntArray>& connectFs,
                   eCase pcase, E_Bool refine = false, eMode mode = SURFACIC);
private:
  typedef K_MESH::Triangle      element_type;
private:
  static NodeAssociator* __buildAssociator(eCase pcase);

  static bool __isFreeContour(const K_FLD::IntArray& connectB, const std::vector<E_Int> & nmates);
  
  static void __get_external_contour(const K_FLD::FloatArray& coord, const std::vector<K_FLD::IntArray>& connectBs, E_Int& i0);

  ///
  GapsManager(void){}
  ~GapsManager(void){}
};
