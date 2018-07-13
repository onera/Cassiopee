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
#ifndef __K_MESH_TETRAHEDRON_H__
#define __K_MESH_TETRAHEDRON_H__

#include "Def/DefTypes.h"
#include "Fld/DynArray.h"
#include "MeshElement/Triangle.h"

static const E_Float aa = 0.25*(1. - 1./::sqrt(5.));
static const E_Float bb = 1. -3.*aa;//(5. + 3.*::sqrt(5.)) / 20.;

namespace K_MESH
{

class Tetrahedron {
  
  public:
    static const E_Int NB_NODES;
    static const E_Int NB_TRIS;
    typedef K_MESH::Triangle       boundary_type;
  
  public:
    Tetrahedron(){}
    ~Tetrahedron(){}
    
    void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
    
    static bool is_inside(const E_Float* Ni, const E_Float* Nj, const E_Float* Nk, const E_Float* Nl, const E_Float* pt);
    
    void triangulate(E_Int* target);
    
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
    E_Int _nodes[4];
};

}
#endif	/* __K_MESH_TETRAHEDRON_H__ */

