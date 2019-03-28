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
  static const E_Int NB_BOUNDS;
  static const E_Int NB_EDGES;
  
  typedef K_MESH::Quadrangle boundary_type;
    
public:
  Hexahedron(){}
  Hexahedron(const E_Int* nodes, E_Int shift=0):_shift(shift){ for (size_t i = 0; i< 8; ++i)_nodes[i]=*(nodes++);}
 ~Hexahedron(){}
 
  E_Int* nodes() { return _nodes;}
  const E_Int* nodes() const { return _nodes;}
  
  E_Int nb_nodes() const {return NB_NODES;}
  E_Int nb_tris() const {return NB_TRIS;}
  
  void setNodes(E_Int* nodes){for (size_t i = 0; i< 8; ++i)_nodes[i]=*(nodes++);}
  
  inline E_Int node(E_Int i){return _nodes[i]+_shift;}
  
  inline E_Int nbounds() { return 6;}
  
  template <typename Connectivity_t>
  inline void set(const K_FLD::ArrayAccessor<Connectivity_t>& connect, E_Int K)
  {connect.getEntry(K, _nodes);}
  
  void triangulate(E_Int* target);//WARNING : connectT3 is Apended (not cleared upon entry)
  
  inline void triangle(E_Int i, E_Int* target)
  {
    assert (i >= 0 && i < NB_TRIS);
    
    switch (i)
    {
      case 0 : target[0] = _nodes[0]; target[1] = _nodes[3]; target[2] = _nodes[1]; break;  //031 BOTTOM
      case 1 : target[0] = _nodes[1]; target[1] = _nodes[3]; target[2] = _nodes[2]; break;  //132
      case 2 : target[0] = _nodes[4]; target[1] = _nodes[5]; target[2] = _nodes[7]; break;  //457 TOP
      case 3 : target[0] = _nodes[7]; target[1] = _nodes[5]; target[2] = _nodes[6]; break;  //756
      
      case 4 : target[0] = _nodes[3]; target[1] = _nodes[0]; target[2] = _nodes[7]; break;  //307 LEFT
      case 5 : target[0] = _nodes[7]; target[1] = _nodes[0]; target[2] = _nodes[4]; break;  //704
      case 6 : target[0] = _nodes[1]; target[1] = _nodes[2]; target[2] = _nodes[6]; break;  //126 RIGHT
      case 7 : target[0] = _nodes[1]; target[1] = _nodes[6]; target[2] = _nodes[5]; break;  //165
      case 8 : target[0] = _nodes[0]; target[1] = _nodes[1]; target[2] = _nodes[5]; break;  //015 FONT
      case 9 : target[0] = _nodes[0]; target[1] = _nodes[5]; target[2] = _nodes[4]; break;  //054
      case 10 : target[0] = _nodes[2]; target[1] = _nodes[3]; target[2] = _nodes[6]; break; //236 BACK
      case 11 : target[0] = _nodes[3]; target[1] = _nodes[7]; target[2] = _nodes[6]; break; //376
      default:break;
    }
  }
  
  ///
  template <typename TriangulatorType, typename acrd_t>
  void triangulate (const TriangulatorType& dt, const acrd_t& acrd) {} //dummy : for genericity
    
  static void get_internal(E_Int* nodes, E_Int* p);
  
  static void get_edges(const E_Int* nodes, Vector_t<K_MESH::NO_Edge>& edges);
  
  static bool cross(const ngon_type& ng, K_FLD::FloatArray& crd, E_Int* face, E_Int nb_faces, K_FLD::FloatArray& data, E_Float* P0, E_Float* P1, E_Float& lambda0, E_Float& lambda1, E_Float tolerance);
  
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
  
  template<typename box_t, typename CoordAcc>
  void bbox(const CoordAcc& acrd, box_t&bb) const
  {
    for (E_Int i = 0; i < 3; ++i)
      {bb.minB[i] = K_CONST::E_MAX_FLOAT; bb.maxB[i] = -K_CONST::E_MAX_FLOAT;}

    bb.compute(acrd, _nodes, NB_NODES, 0/*idx start*/);
  }
  
  template< typename ngo_t>
  static void reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i);
  
  
  ///
  template <typename CoordAcc>
  inline void iso_barycenter(const CoordAcc& coord, E_Float* G);
  
  template <typename ngunit_t>
  static inline void iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G);
  
private:
  
  Hexahedron(const Hexahedron& orig);
 
private:
    E_Int _shift;
    E_Int _nodes[8];

};

template< typename ngo_t>
void Hexahedron::reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i) // bot, top, left, right, front, back
{
  std::map<E_Int,E_Int> glmap; // crd1 to 0-26 indexes
  E_Int nb_faces = ng.PHs.stride(i); 
  E_Int* faces = ng.PHs.get_facets_ptr(i);
  E_Int PGi = faces[0] - 1;
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);
  
  // but convention, first face is bottom, first node is 0 in local numbering (0 to 26)

  glmap[*pN] = 0; // PHi(0,0) -> 0  
  glmap[*(pN+1)] = 1;
  glmap[*(pN+2)] = 2;
  glmap[*(pN+3)] = 3;

  if (F2E(1,PGi) != i) // for BOT, PH is the right element. if not, wrong orientation => swap of 1 and 3
  { 
    glmap[*(pN+3)] = 1;
    glmap[*(pN+1)] = 3;
  }
  E_Int TopId(E_IDX_NONE),LeftId(E_IDX_NONE),RightId(E_IDX_NONE),FrontId(E_IDX_NONE),BackId(E_IDX_NONE);

  for (int k = 1; k < 6; ++k)
  {
    int count = 0;
    Vector_t<bool> commonNodes(4,false);
    E_Int testedPG = faces[k]-1;
    E_Int* pNode = ng.PGs.get_facets_ptr(testedPG);

    for (int j = 0; j < 4; ++j)
    {
      auto it = glmap.find(pNode[j]);
      if (it != glmap.end())
      {
        // found
        count++;
        commonNodes[it->second] = true;
      }
    }
    if (count == 0) // no common point, the ith PG is the TOP
      TopId = k;
    else if (commonNodes[0] && commonNodes[1])
      FrontId = k;
    else if (commonNodes[1] && commonNodes[2])
      RightId = k;
    else if (commonNodes[2] && commonNodes[3])
      BackId = k;
    else if (commonNodes[0] && commonNodes[3])
      LeftId = k;
  }
  
  E_Int mol[6];

  mol[0] = faces[0];
  mol[1] = faces[TopId];
  mol[2] = faces[LeftId];
  mol[3] = faces[RightId];
  mol[4] = faces[FrontId];
  mol[5] = faces[BackId];

  for (int i = 0; i < nb_faces; ++i)
    faces[i] = mol[i];
}

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

  template <typename ngunit_t>
  inline void Hexahedron::iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G)
  {
    //WARNING : assuming reodrederd pgs : first is bottom, second is top
    
    E_Int new_bary[8];

    for (int i = 0; i < 2; ++i) // 8 points : bot and top nodes
    {
      const E_Int* nodes = PGs.get_facets_ptr(first_pg[i]-index_start);
      E_Int nb_nodes = PGs.stride(first_pg[i]-index_start);
      
      for (int k = 0; k  < nb_nodes; ++k)
        new_bary[nb_nodes*i+k] = nodes[k];   
    }
    
    K_MESH::Polyhedron<STAR_SHAPED>::iso_barycenter(crd, new_bary, 8, 1, G);
  }

}
#endif	/* __K_MESH_HEXAHEDRON_H__ */
