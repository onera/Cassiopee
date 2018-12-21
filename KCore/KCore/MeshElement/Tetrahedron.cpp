/*    
    Copyright 2013-2019 Onera.

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

#include "Tetrahedron.h"

namespace K_MESH
{
  
const E_Int K_MESH::Tetrahedron::NB_NODES=4;
const E_Int K_MESH::Tetrahedron::NB_TRIS=4;
const E_Int K_MESH::Tetrahedron::NB_BOUNDS=4;

bool Tetrahedron::is_inside
(const E_Float* Ni, const E_Float* Nj, const E_Float* Nk, const E_Float* Nl, const E_Float* pt)
{
  if (K_FUNC::zzdet4(Ni, Nj, Nl, pt) > E_EPSILON)
    return false;
  if (K_FUNC::zzdet4(Nj, Nk, Nl, pt) > E_EPSILON)
    return false;
  if (K_FUNC::zzdet4(Nk, Ni, Nl, pt) > E_EPSILON)
    return false;
  if (K_FUNC::zzdet4(Ni, Nk, Nj, pt) > E_EPSILON)
    return false;
  return true;
}

void Tetrahedron::triangulate(E_Int* target)
{  
  target[0]=_nodes[0]; target[1]=_nodes[1]; target[2]=_nodes[3];
  target[3]=_nodes[0]; target[4]=_nodes[3]; target[5]=_nodes[2];
  target[6]=_nodes[0]; target[7]=_nodes[2]; target[8]=_nodes[1];
  target[9]=_nodes[1]; target[10]=_nodes[2]; target[11]=_nodes[3];
}

}


