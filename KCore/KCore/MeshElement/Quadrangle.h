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

#ifndef __K_MESH_QUADRANGLE_H__
#define __K_MESH_QUADRANGLE_H__

#include "Def/DefTypes.h"
#include "Fld/DynArray.h"

namespace K_MESH
{

  class Quadrangle
  {
  public:
    static const E_Int NB_NODES;
  public:
    void setNodes(E_Int n0, E_Int n1, E_Int n2, E_Int n3){_nodes[0]=n0;_nodes[1]=n1;_nodes[2]=n2;_nodes[3]=n3;}
    //void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
    inline E_Int nnodes() { return 4;}
    
    ///
    inline const E_Int* nodes(){return _nodes;}
    
    /// Gets the i-th node.
    inline const E_Int& node(const E_Int& i) const {return _nodes[i];}
    
    static void triangulate(const E_Int* nodes, E_Int* target);//WARNING : triangulation is put at target address contiguously
    
  private:
    E_Int _nodes[4];
  };
}

#endif
