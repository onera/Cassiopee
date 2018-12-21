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
#ifndef __K_MESH_PYRAMID_H__
#define __K_MESH_PYRAMID_H__

#include "Def/DefTypes.h"
#include "Fld/DynArray.h"


namespace K_MESH
{

class Pyramid {
  
  public:
    static const E_Int NB_NODES = 5;
    static const E_Int NB_TRIS = 6;
    static const E_Int NB_BOUNDS = 5;
      
  public:
    Pyramid():_shift(0){}
    ~Pyramid(){}
    
    Pyramid(const E_Int* nodes, E_Int shift=0):_shift(shift){ for (size_t i = 0; i< NB_NODES; ++i)_nodes[i]=*(nodes++) + shift;}
    
    inline E_Int node(E_Int i){return _nodes[i]+_shift;}
    
    E_Int* nodes() { return _nodes;}
    // const E_Int* nodes() const { return _nodes;}
    
    // E_Int nb_nodes() const {return NB_NODES;}
    E_Int nb_tris() const {return NB_TRIS;}
    
    //void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
    
    
    
    ///
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& acrd) {} // dummy since it is a basic element
        
    inline void triangle(E_Int i, E_Int* target)
    {
      assert (i >= 0 && i < NB_TRIS);
      
      switch (i)
      {
        case 0 : target[0] = _nodes[0]; target[1] = _nodes[1]; target[2] = _nodes[4]; break;  //014
        case 1 : target[0] = _nodes[1]; target[1] = _nodes[2]; target[2] = _nodes[4]; break;  //124
        case 2 : target[0] = _nodes[2]; target[1] = _nodes[3]; target[2] = _nodes[4]; break;  //234
        case 3 : target[0] = _nodes[3]; target[1] = _nodes[0]; target[2] = _nodes[4]; break;  //304
        
        case 4 : target[0] = _nodes[2]; target[1] = _nodes[1]; target[2] = _nodes[0]; break;  //210
        case 5 : target[0] = _nodes[2]; target[1] = _nodes[0]; target[2] = _nodes[3]; break;  //203
        
        default:break;
      }
    }
      
      
   ///  
   template <typename CoordAcc> inline
    void iso_barycenter(const CoordAcc& coord, E_Float* G)
    { 
      //
      for (size_t d=0; d < 3; ++d) G[d]=0.;
    
      for (E_Int i=0; i < NB_NODES; ++i)
      {
        for (size_t d=0; d < 3; ++d)
        {
          //std::cout << "v : " << coord.getVal(node(i), d) << std::endl;
          G[d] += coord.getVal(node(i), d);
        }
      }
    
      E_Float k = 1./(E_Float)NB_NODES;
    
      for (size_t i = 0; i < 3; ++i) G[i] *= k;
      //std::cout << "G : " << G[0] << "/" << G[1] << "/" << G[2] << std::endl;
    
    }
    
    ///
    template<typename box_t, typename CoordAcc>
    void bbox(const CoordAcc& acrd, box_t&bb) const
    {
      for (E_Int i = 0; i < 3; ++i)
        {bb.minB[i] = K_CONST::E_MAX_FLOAT; bb.maxB[i] = -K_CONST::E_MAX_FLOAT;}

      bb.compute(acrd, _nodes, NB_NODES, _shift/*idx start*/);
    }
    
private:
  
    Pyramid(const Pyramid& orig);
  
private:
    E_Int _shift;
    E_Int _nodes[5];
};

}
#endif	/* __K_MESH_TETRAHEDRON_H__ */

