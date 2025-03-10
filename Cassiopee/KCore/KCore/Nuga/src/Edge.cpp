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

#include "Nuga/include/Edge.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/maths.hxx"
#include <math.h>

//=============================================================================
void
K_MESH::Edge::getBoundary(const Edge& E1, const Edge& E2, boundary_type& b)
{
  b = IDX_NONE;
  const E_Int& n10 = E1._nodes[0];
  const E_Int& n11 = E1._nodes[1];
  const E_Int& n20 = E2._nodes[0];
  const E_Int& n21 = E2._nodes[1];

  if (n11 == n20) b = n11;
  else if (n10 == n21) b = n10;
  else if (n10 == n20) b = n10;
  else if (n11 == n21) b = n11;
}
