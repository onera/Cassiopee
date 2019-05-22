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
#ifndef __K_MESH_TETRAHEDRON_H__
#define __K_MESH_TETRAHEDRON_H__

#include "Def/DefTypes.h"
#include "Fld/DynArray.h"
#include "MeshElement/Triangle.h"
#define Vector_t std::vector


static const E_Float aa = 0.25*(1. - 1./::sqrt(5.));
static const E_Float bb = 1. -3.*aa;//(5. + 3.*::sqrt(5.)) / 20.;

namespace K_MESH
{

class Tetrahedron {
  
  public:
    static const E_Int NB_NODES;
    static const E_Int NB_TRIS;
    static const E_Int NB_BOUNDS;
    static const E_Int NB_EDGES;

    typedef K_MESH::Triangle       boundary_type;
  
  public:
    Tetrahedron():_shift(0){}
    ~Tetrahedron(){}
    
    Tetrahedron(const E_Int* nodes, E_Int shift=0):_shift(shift){ for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++) + shift;}
    
    // Constructor from a NGON. non-oriented version
    template <typename ngunit_t>
    Tetrahedron(const ngunit_t & PGs, const E_Int* first_pg):_shift(0)
    {
      // WARNING : non-oriented _nodes upon exit
      // WARNING : assume a index start at 1 on node upon entry
      
      // work on the first 2 triangle to find out the nodes
      
      const E_Int* nodes = PGs.get_facets_ptr(first_pg[0]-1);
       
      
#ifdef DEBUG_MESH_ELEMENT
      E_Int nb_nodes = PGs.stride(first_pg[0]-1);
      assert (nb_nodes == 3);
#endif
      _nodes[0] = nodes[0]-1;
      _nodes[1] = nodes[1]-1;
      _nodes[2] = nodes[2]-1;
          
      nodes = PGs.get_facets_ptr(first_pg[1]-1);
      

#ifdef DEBUG_MESH_ELEMENT
      nb_nodes = PGs.stride(first_pg[1]-1); 
      assert (nb_nodes == 3);
#endif

      if ((nodes[0]-1 != _nodes[0]) && (nodes[0]-1 != _nodes[1]) && (nodes[0]-1 != _nodes[2])) _nodes[3] = nodes[0]-1;
      else if ((nodes[1]-1 != _nodes[0]) && (nodes[1]-1 != _nodes[1]) && (nodes[1]-1 != _nodes[2])) _nodes[3] = nodes[1]-1;
      else if ((nodes[2]-1 != _nodes[0]) && (nodes[2]-1 != _nodes[1]) && (nodes[2]-1 != _nodes[2])) //_nodes[3] = nodes[1];
        _nodes[3] = nodes[2]-1;
    }
    
    inline E_Int node(E_Int i){return _nodes[i]+_shift;}
    
    E_Int* nodes() { return _nodes;}
    const E_Int* nodes() const { return _nodes;}
    
    E_Int nb_nodes() const {return NB_NODES;}
    E_Int nb_tris() const {return NB_TRIS;}
    
    void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
    
    static bool is_inside(const E_Float* Ni, const E_Float* Nj, const E_Float* Nk, const E_Float* Nl, const E_Float* pt);
    
    void triangulate(E_Int* target);
    
    E_Float quality(const K_FLD::FloatArray& crd, E_Float* Vol);
    
    void edge_length_extrema(const K_FLD::FloatArray& crd, E_Float& Lmin2, E_Float& Lmax2);
    
    ///
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& acrd) {} //dummy : for genericity
        
    static void get_edges(const E_Int* nodes, Vector_t<K_MESH::NO_Edge>& edges);

    
    inline void triangle(E_Int i, E_Int* target)
    {
      assert (i >= 0 && i < 4);
      
      switch (i)
      {
        case 0 : target[0]=_nodes[0]; target[1]=_nodes[1]; target[2]=_nodes[3]; break;
        case 1 : target[0]=_nodes[0]; target[1]=_nodes[3]; target[2]=_nodes[2]; break;
        case 2 : target[0]=_nodes[0]; target[1]=_nodes[2]; target[2]=_nodes[1]; break;
        case 3 : target[0]=_nodes[1]; target[1]=_nodes[2]; target[2]=_nodes[3]; break;
        default:break;
      }
    }
    
    template< typename ngo_t>
    static void reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i);
    
    
    
    ///
    static inline E_Float volume(const E_Float* p1, const E_Float* p2, const E_Float* p3, const E_Float* p4)
    {return (K_FUNC::zzdet4(p1, p2, p3, p4)/6.);}
    
    inline void getBoundary(E_Int n, boundary_type& b) const {
    
    switch (n)
    {
      case 0: b.setNodes(_nodes[1], _nodes[2], _nodes[3]);break;
      case 1: b.setNodes(_nodes[0], _nodes[3], _nodes[2]);break;
      case 2: b.setNodes(_nodes[0], _nodes[1], _nodes[3]);break;
      case 3: b.setNodes(_nodes[0], _nodes[2], _nodes[1]);break;
      default : break;
    }
  }
    
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
  
  template <typename ngunit_t>
  static inline void iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G)
  {    
    Tetrahedron t(PGs, first_pg);
    using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
    acrd_t acrd(crd);
    t.iso_barycenter<acrd_t>(acrd, G);
  }
  
  template<typename box_t, typename CoordAcc>
  void bbox(const CoordAcc& acrd, box_t&bb) const
  {
    for (E_Int i = 0; i < 3; ++i)
      {bb.minB[i] = K_CONST::E_MAX_FLOAT; bb.maxB[i] = -K_CONST::E_MAX_FLOAT;}

    bb.compute(acrd, _nodes, NB_NODES, _shift/*idx start*/);
  }
    
    /*void split8(K_FLD::FloatArray& crd, K_FLD::IntArray& cntTH4)
    {
      
      E_Int& v0 = _nodes[0];
      E_Int& v1 = _nodes[1];
      E_Int& v2 = _nodes[3];
      E_Int& v3 = _nodes[4];
      
      E_Int v4= crd.cols();
      E_Int v5=v4+1;
      E_Int v6=v4+2;
      E_Int v7=v4+3;
      E_Int v8=v4+4;
      E_Int v9=v4+5;
      
      crd.resize(3, v4+6, 0.);
      
      K_FUNC::sum<3>(0.5, crd.col(v0), 0.5, crd.col(v2), crd.col(v4));
      K_FUNC::sum<3>(0.5, crd.col(v1), 0.5, crd.col(v2), crd.col(v5));
      K_FUNC::sum<3>(0.5, crd.col(v0), 0.5, crd.col(v1), crd.col(v6));
      
      K_FUNC::sum<3>(0.5, crd.col(v2), 0.5, crd.col(v3), crd.col(v7));
      K_FUNC::sum<3>(0.5, crd.col(v1), 0.5, crd.col(v3), crd.col(v8));
      K_FUNC::sum<3>(0.5, crd.col(v0), 0.5, crd.col(v3), crd.col(v9));
      
      E_Int TH4[4];

      TH4[0]=v1;
      TH4[1]=v5;
      TH4[2]=v6;
      TH4[3]=v8;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v5;
      TH4[1]=v2;
      TH4[2]=v4;
      TH4[3]=v7;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v6;
      TH4[1]=v4;
      TH4[2]=v0;
      TH4[3]=v9;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v8;
      TH4[1]=v7;
      TH4[2]=v9;
      TH4[3]=v3;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v6;
      TH4[1]=v5;
      TH4[2]=v4;
      TH4[3]=v8;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v9;
      TH4[1]=v6;
      TH4[2]=v4;
      TH4[3]=v8;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v8;
      TH4[1]=v7;
      TH4[2]=v5;
      TH4[3]=v4;
      cntTH4.pushBack(TH4, TH4+4);
      
      TH4[0]=v7;
      TH4[1]=v8;
      TH4[2]=v9;
      TH4[3]=v4;
      cntTH4.pushBack(TH4, TH4+4);
      
    }*/
        
    static void eval(const E_Float* p0, const E_Float* p1, const E_Float* p2, const E_Float* p3,
                     E_Float u, E_Float v, E_Float w, E_Float t, E_Float* X)
    {            
      X[0] = u*p0[0] + v*p1[0] + w*p2[0] + t*p3[0];
      X[1] = u*p0[1] + v*p1[1] + w*p2[1] + t*p3[1];
      X[2] = u*p0[2] + v*p1[2] + w*p2[2] + t*p3[2];
    }
    

    
    static void integ_O2(const E_Float* p0, const E_Float* p1, const E_Float* p2, const E_Float* p3,
                         const E_Float* X0, E_Float f0, const E_Float* grad0, E_Float& val)
    {
      val = 0.;
      
      E_Float X[3], X0X[3], gradf, f;;
      
      eval(p0,p1,p2,p3, aa, aa, aa, bb, X);
      K_FUNC::diff<3>(X, X0, X0X);
      gradf = K_FUNC::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      
      eval(p0,p1,p2,p3, bb, aa, aa, aa, X);
      K_FUNC::diff<3>(X, X0, X0X);
      gradf = K_FUNC::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      
      eval(p0,p1,p2,p3, aa, bb, aa, aa, X);
      K_FUNC::diff<3>(X, X0, X0X);
      gradf = K_FUNC::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      
      eval(p0,p1,p2,p3, aa, aa, bb, aa, X);
      K_FUNC::diff<3>(X, X0, X0X);
      gradf = K_FUNC::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      val /= 4.;
      
    }
    
    typedef E_Float (*pFunc)(E_Float x, E_Float y, E_Float z);
    
    static void integ_O2(const E_Float* p0, const E_Float* p1, const E_Float* p2, const E_Float* p3,
                               pFunc func, E_Float& val)
    {
      val = 0.;
      
      E_Float X[3];
      
      eval(p0,p1,p2,p3, aa, aa, aa, bb, X); 
      val +=func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, bb, aa, aa, aa, X);
      val +=func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, aa, bb, aa, aa, X);
      val +=func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, aa, aa, bb, aa, X);
      val +=func(X[0], X[1], X[2]);

      val /= 4.;
      
    }
    
    static void integ_O3(const E_Float* p0, const E_Float* p1, const E_Float* p2, const E_Float* p3,
                               pFunc func, E_Float& val)
    {
      val = 0.;
      
      E_Float X[3];
      eval(p0,p1,p2,p3, 0.25, 0.25, 0.25, 0.25, X); 
      val +=(-4./5.)*func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, 1./6., 1./6., 1./6., 0.5, X); 
      val +=(9./20.)*func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, 1./6., 1./6., 0.5, 1./6., X); 
      val +=(9./20.)*func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, 1./6., 0.5, 1./6., 1./6., X); 
      val +=(9./20.)*func(X[0], X[1], X[2]);
      
      eval(p0,p1,p2,p3, 0.5, 1./6., 1./6., 1./6., X); 
      val +=(9./20.)*func(X[0], X[1], X[2]);
      
    }
  
private:
  
  Tetrahedron(const Tetrahedron& orig);
  
private:
    E_Int _shift;
    E_Int _nodes[4];
};

template< typename ngo_t>
void Tetrahedron::reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i)
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

  if (F2E(1,PGi) != i) // for BOT, PH is the right element. if not, wrong orientation => swap of 1 and 3
  { 
    glmap[*(pN+1)] = 2;
    glmap[*(pN+2)] = 1;
  }
  E_Int F1Id(E_IDX_NONE), F2Id(E_IDX_NONE), F3Id(E_IDX_NONE);

  for (int k = 1; k < 4; ++k)
  {
    int count = 0;
    std::vector<bool> commonNodes(3,false);
    E_Int testedPG = faces[k]-1;
    E_Int* pNode = ng.PGs.get_facets_ptr(testedPG);

    for (int j = 0; j < 3; ++j)
    {
      auto it = glmap.find(pNode[j]);
      if (it != glmap.end())
      {
        // found
        count++;
        commonNodes[it->second] = true;
      }
    }
    if (commonNodes[0] && commonNodes[1])
      F1Id = k;
    else if (commonNodes[1] && commonNodes[2])
      F2Id = k;
    else if (commonNodes[2] && commonNodes[0])
      F3Id = k;
    }

  E_Int mol[4];

  mol[0] = faces[0];
  mol[1] = faces[F1Id];
  mol[2] = faces[F2Id];
  mol[3] = faces[F3Id];

  for (int i = 0; i < nb_faces; ++i)
    faces[i] = mol[i];
}

}
#endif	/* __K_MESH_TETRAHEDRON_H__ */

