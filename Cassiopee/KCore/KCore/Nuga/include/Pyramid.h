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

#ifndef __K_MESH_PYRAMID_H__
#define __K_MESH_PYRAMID_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/ArrayAccessor.h"


namespace K_MESH
{

class Pyramid {

  public:
    static constexpr E_Int NB_NODES = 5;
    static constexpr E_Int NB_TRIS = 6;
    static constexpr E_Int NB_BOUNDS = 5;
      
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
    
    template <typename CoordAcc> inline
    void iso_barycenter(const CoordAcc& coord, E_Float* G);
    template <typename ngunit_t> inline
    static void iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G);
    
    template< typename ngo_t>
    static void reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i);

    template <typename ngunit_t>
    static bool is_of_type(const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs)
    {
      if (nb_pgs != 5) return false;
      E_Int s1(0), s2(0);

      for (int i = 0; i<5; i++)
      {
        if (PGs.stride(*(first_pg + i) - 1) == 3) ++s1;
        else if (PGs.stride(*(first_pg + i) - 1) == 4) ++s2;
        else return false;
      }

      return ((s1 == 4) && (s2 == 1));
    }
    
    ///
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& acrd) {} // dummy since it is a basic element
    ///
    template <typename acrd_t>
    E_Int cvx_triangulate (const acrd_t& acrd) {return 0;}
        
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
    template<typename box_t, typename CoordAcc>
    void bbox(const CoordAcc& acrd, box_t&bb) const
    {
      for (E_Int i = 0; i < 3; ++i)
        {bb.minB[i] = NUGA::FLOAT_MAX; bb.maxB[i] = -NUGA::FLOAT_MAX;}

      bb.compute(acrd, _nodes, NB_NODES, _shift/*idx start*/);
    }

private:
    E_Int _shift;
    E_Int _nodes[5];
};

///  
template <typename CoordAcc> inline
void Pyramid::iso_barycenter(const CoordAcc& coord, E_Float* G)
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
template <typename ngunit_t>
inline void Pyramid::iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G)
{
  //WARNING : assuming reodrederd pgs : first is bottom, second is top

  E_Int new_bary[5];

  // face BOT
  const E_Int* nodes = PGs.get_facets_ptr(first_pg[0]-index_start);
  E_Int nb_nodes = PGs.stride(first_pg[0]-index_start);

  for (int k = 0; k  < nb_nodes; ++k)
    new_bary[k] = nodes[k];   

  // Summit
  const E_Int* F1 = PGs.get_facets_ptr(first_pg[1]-index_start); // face F1
  E_Int nb_F1 = PGs.stride(first_pg[1]-index_start);
  
  for (int i=0; i< nb_F1; i++){
      E_Int n(0);
      for (int j=0; j < nb_nodes; j++){
          if ( *(F1+i) != *(nodes+j) ){
              n++;
          }
      }
      if ( n == (nb_nodes-1))
        new_bary[4] = *(F1+i);
  }
  
  K_MESH::Polyhedron<STAR_SHAPED>::iso_barycenter(crd, new_bary, 5, 1, G);
}

///
template< typename ngo_t>
void Pyramid::reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i)
{
  std::map<E_Int,E_Int> glmap; 
  
  assert (ng.PHs.stride(i) == 5);
  E_Int* faces = ng.PHs.get_facets_ptr(i);
  E_Int BOT = IDX_NONE;

  E_Int ibot(0);
  for (int i=0; i< 5; i++)
  {
    if (ng.PGs.stride(faces[i] - 1)!=4) // BOTTOM is the QUAD
      continue;
    
    assert (BOT == IDX_NONE); // one QUAD only
    
    BOT = faces[i] - 1;
    ibot = i;
  }
  
  assert (BOT != IDX_NONE);
  
  E_Int* pN = ng.PGs.get_facets_ptr(BOT);

  
  glmap[*pN] = 0;
  glmap[*(pN+1)] = 1;
  glmap[*(pN+2)] = 2;
  glmap[*(pN+3)] = 3;


  if (F2E(1,BOT) != i) // for BOT, PH is the right element. if not, wrong orientation => swap of 1 and 3
  { 
    glmap[*(pN+1)] = 3;
    glmap[*(pN+3)] = 1;
  }

  E_Int F1Id(IDX_NONE), F2Id(IDX_NONE), F3Id(IDX_NONE), F4Id(IDX_NONE);

  bool commonNodes[4];
  for (int k = 1; k < 5; ++k)
  {
    commonNodes[0] = commonNodes[1] = commonNodes[2] = commonNodes[3] = false;
    int ki = (ibot + k) % 5;
    
    E_Int testedPG = faces[ki]-1;
    
    E_Int* pNode = ng.PGs.get_facets_ptr(testedPG);
    assert(ng.PGs.stride(testedPG) == 3);

    for (int j = 0; j < 3; ++j)
    {
      auto it = glmap.find(pNode[j]);
      if (it != glmap.end())
        commonNodes[it->second] = true;
    }
    if (commonNodes[0] && commonNodes[1])
      F1Id = ki;
    else if (commonNodes[1] && commonNodes[2])
      F2Id = ki;
    else if (commonNodes[2] && commonNodes[3])
      F3Id = ki;
    else if (commonNodes[3] && commonNodes[0])
      F4Id = ki;
  }
  
  assert (F1Id != IDX_NONE);
  assert (F2Id != IDX_NONE);
  assert (F3Id != IDX_NONE);
  assert (F4Id != IDX_NONE);

  assert (F1Id != F2Id && F1Id != F3Id && F1Id != F4Id);
  assert (F2Id != F1Id && F2Id != F3Id && F2Id != F4Id);
  assert (F3Id != F1Id && F3Id != F2Id && F3Id != F4Id);
  assert (F4Id != F1Id && F4Id != F2Id && F4Id != F3Id);

  E_Int mol[] = {faces[ibot], faces[F1Id], faces[F2Id], faces[F3Id], faces[F4Id]};

  for (int i = 0; i < 5; ++i)
    faces[i] = mol[i];
}

}
#endif	/* __K_MESH_PYRAMID_H__ */

