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
#ifndef __K_MESH_HEXAHEDRON_H__
#define __K_MESH_HEXAHEDRON_H__

#include "MeshElement/Quadrangle.h"
#include "Fld/DynArray.h"
#include "Fld/ArrayAccessor.h"
#include "Fld/ngon_t.hxx"
#include "Def/DefTypes.h"

namespace K_MESH
{

class Hexahedron {
  
public:
  static const E_Int NB_NODES;
  static const E_Int NB_TRIS;
  
  typedef K_MESH::Quadrangle boundary_type;
    
public:
  Hexahedron(){}
  Hexahedron(E_Int* nodes, E_Int shift=0):_shift(shift){ for (size_t i = 0; i< 8; ++i)_nodes[i]=*(nodes++);}
 ~Hexahedron(){}
  
  void setNodes(E_Int* nodes){for (size_t i = 0; i< 8; ++i)_nodes[i]=*(nodes++);}
  
  inline E_Int node(E_Int i){return _nodes[i]+_shift;}
  
  inline E_Int nbounds() { return 6;}
  
  template <typename Connectivity_t>
  inline void set(const K_FLD::ArrayAccessor<Connectivity_t>& connect, E_Int K)
  {connect.getEntry(K, _nodes);}
  
  void triangulate(E_Int* target);//WARNING : connectT3 is Apended (not cleared upon entry)
  
  static void reorder_pgs(ngon_type& ng, const K_FLD::IntArray& F2E, E_Int i);
  
  static void get_internal(E_Int* nodes, E_Int* p);
  
  static void get_orient(const ngon_type& ng, const K_FLD::IntArray& F2E, E_Int PHi, E_Int* PHi_orient);
  
  static bool pt_is_inside(const ngon_type& ng, const K_FLD::FloatArray& crd, E_Int PHi, const E_Int* PHi_orient, const E_Float* pt);
  
  inline void getBoundary(E_Int n, boundary_type& b) const {
    
    switch (n)
    {
      case 0: b.setNodes(_nodes[0], _nodes[1], _nodes[2], _nodes[3]);break;
      case 1: b.setNodes(_nodes[4], _nodes[5], _nodes[6], _nodes[7]);break;
      case 2: b.setNodes(_nodes[0], _nodes[4], _nodes[7], _nodes[3]);break;
      case 3: b.setNodes(_nodes[1], _nodes[2], _nodes[6], _nodes[5]);break;
      case 4: b.setNodes(_nodes[0], _nodes[1], _nodes[5], _nodes[4]);break;
      case 5: b.setNodes(_nodes[2], _nodes[3], _nodes[7], _nodes[6]);break;
      default : break;
    }
  }
  
  inline void getBoundary(E_Int n, E_Int* nodes) const {
    
    switch (n)
    {
      case 0: nodes[0]=_nodes[0]; nodes[1]=_nodes[1]; nodes[2]=_nodes[2]; nodes[3]=_nodes[3];break;
      case 1: nodes[0]=_nodes[4]; nodes[1]=_nodes[5]; nodes[2]=_nodes[6]; nodes[3]=_nodes[7];break;
      case 2: nodes[0]=_nodes[0]; nodes[1]=_nodes[4]; nodes[2]=_nodes[7]; nodes[3]=_nodes[3];break;
      case 3: nodes[0]=_nodes[1]; nodes[1]=_nodes[2]; nodes[2]=_nodes[6]; nodes[3]=_nodes[5];break;
      case 4: nodes[0]=_nodes[0]; nodes[1]=_nodes[1]; nodes[2]=_nodes[5]; nodes[3]=_nodes[4];break;
      case 5: nodes[0]=_nodes[2]; nodes[1]=_nodes[3]; nodes[2]=_nodes[7]; nodes[3]=_nodes[6];break;
      default : break;
    }
  }
  
  ///
  template <typename CoordAcc>
  inline void iso_barycenter(const CoordAcc& coord, E_Float* G);
  
private:
  
  Hexahedron(const Hexahedron& orig);
 
private:
    E_Int _shift;
    E_Int _nodes[8];

};

template <typename CoordAcc> inline
void Hexahedron::iso_barycenter(const CoordAcc& coord, E_Float* G)
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

}
#endif	/* __K_MESH_HEXAHEDRON_H__ */
