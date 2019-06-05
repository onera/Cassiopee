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
    static const E_Int NB_TRIS;
    
    Quadrangle(){};
    
    template <typename NodeIterator>
    explicit Quadrangle(const NodeIterator pN)
    {_nodes[0] = *pN; _nodes[1] = *(pN+1); _nodes[2] = *(pN+2); _nodes[3] = *(pN+3);}

    void setNodes(E_Int n0, E_Int n1, E_Int n2, E_Int n3){_nodes[0]=n0;_nodes[1]=n1;_nodes[2]=n2;_nodes[3]=n3;}
    //void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
       
    E_Int nb_nodes() const {return NB_NODES;}
    E_Int nb_tris() const { return NB_TRIS;}
    E_Int nbounds() const { return 4;}
    
    ///
    inline E_Int* nodes(){return _nodes;} 
    inline const E_Int* nodes() const {return _nodes;}
    
    /// Gets the i-th node.
    inline const E_Int& node(const E_Int& i) const {return _nodes[i];}
    
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& crd){}
    
    inline void triangle(E_Int i, E_Int* target)
    {
      assert (i >= 0 && i < 2);
      
      if (i==0)
      {
        target[0]=_nodes[0]; target[1]=_nodes[1]; target[2]=_nodes[2];
      }
      else
      {
        target[0]=_nodes[0]; target[1]=_nodes[2]; target[2]=_nodes[3];
      }
    }
    
    static inline void triangulate(const E_Int* nodes, E_Int* target)//WARNING : triangulation is put at target address contiguously
    {
      E_Int t[3];
      t[0] = nodes[0]; t[1] = nodes[1]; t[2] = nodes[3];
      std::copy(&t[0], &t[0]+3, target);
      std::copy(&nodes[0]+1, &nodes[0]+4, target+3);
    }
    
    template<typename box_t, typename CoordAcc>
    void bbox(const CoordAcc& acrd, box_t&bb) const
    {
      for (E_Int i = 0; i < 3; ++i)
      {bb.minB[i] = K_CONST::E_MAX_FLOAT; bb.maxB[i] = -K_CONST::E_MAX_FLOAT;}

      box_t b;
      b.compute(acrd, _nodes, NB_NODES, 0/*idx start*/);
      for (E_Int i = 0; i < NB_NODES; ++i)
      {
        bb.minB[i] = std::min(bb.minB[i], b.minB[i]);
        bb.maxB[i] = std::max(bb.maxB[i], b.maxB[i]);
      }
    }
    
    template <typename CoordAcc>
    void iso_barycenter(const CoordAcc& coord, E_Float* G)
    { 
      //
      for (size_t d=0; d < 3; ++d) G[d]=0.;

      for (E_Int i=0; i < NB_NODES; ++i)
      {
        for (size_t d=0; d < 3; ++d)
        {
          //std::cout << "v : " << coord.getVal(node(i), d) << std::endl;
          G[d] += coord.getVal(_nodes[i], d);
        }
      }
    }
    
  private:
    E_Int _nodes[4];
  };
}

#endif
