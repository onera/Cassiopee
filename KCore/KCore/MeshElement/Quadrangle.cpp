/*    
    Copyright 2013-2018 Onera.

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

#include "Quadrangle.h"
#include "Triangle.h"

const E_Int K_MESH::Quadrangle::NB_NODES = 4;

namespace K_MESH
{
  void Quadrangle::triangulate(const E_Int* nodes, E_Int* target)
  {    
    E_Int t[3];
    t[0] = nodes[0]; t[1] = nodes[1]; t[2] = nodes[3];
    std::copy(&t[0], &t[0]+3, target);
    std::copy(&nodes[0]+1, &nodes[0]+4, target+3);
  }
}
