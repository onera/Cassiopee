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

#ifndef NUGA_ADAPT_CELLS_H
#define NUGA_ADAPT_CELLS_H

#include "Nuga/include/cdatastruct.hxx"
//#include <vector>

namespace NUGA
{
  // adapts m with respect to input source points
  // Returns a polyhedral conformal mesh for polygons and polyhedra.
  // Just the leaves : all intermediary levels are removed
  int adapt_cells(c_phmesh_t& m,             // i/o mesh
                  const c_crd3D_t& src_pts);  // source points

  //// adapts m with respect to input source points
  //// Returns a polyhedral conformal mesh with ascendent history for polygons and polyhedra.
  //// Just the leaves : all intermediary levels are removed
  //int adapt_cells(c_phmesh_t& m,             // i/o mesh
  //  const c_crd3D_t& src_pts,                // source points
  //  std::vector<int>& pgoids,                // polygon original id in m before adaptation
  //  std::vector<int>& phoids);               // polyhedron original id in m before adaptation

  // adapts m with respect to input source points
  // PERSISTANT : returns the entire genealogy
  //int adapt_cells(const phmesh_t& m, const crd3D_t& src_pts, hmesh& hm);
}

#endif