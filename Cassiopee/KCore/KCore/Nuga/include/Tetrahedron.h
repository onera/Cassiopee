/*    
    Copyright 2013-2026 ONERA.

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

#ifndef __K_MESH_TETRAHEDRON_H__
#define __K_MESH_TETRAHEDRON_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/IdTool.h"

#define Vector_t std::vector


static const E_Float aa = 0.25*(1. - 1./sqrt(5.));
static const E_Float bb = 1. -3.*aa;//(5. + 3.*sqrt(5.)) / 20.;

namespace K_MESH
{

class Tetrahedron 
{  
  public:
    static constexpr E_Int NB_NODES = 4;
    static constexpr E_Int NB_TRIS = 4;
    static constexpr E_Int NB_BOUNDS = 4;
    static constexpr E_Int NB_EDGES = 6;

    typedef K_MESH::Triangle       boundary_type;

    enum eShapeType { REGULAR, SPIKE, SLICE1, SLICE2, KNIFE1, KNIFE2, DELTA, UMBRELLA};
  
  public:
    Tetrahedron():_shift(0){}
    ~Tetrahedron(){}
    
    Tetrahedron(const E_Int* nodes, E_Int shift=0):_shift(shift){ for (size_t i = 0; i< 4; ++i)_nodes[i]=nodes[i] + shift;}
    
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
      _nodes[3] = -1;
          
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
    
    template <typename ngunit_t>
    static bool is_of_type(const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs) {
      if (nb_pgs != 4) return false;

      for (int i = 0; i<4; i++)
        if (PGs.stride(*(first_pg + i) - 1) != 3) return false;

      return true;
    }

    eShapeType shape_type(const K_FLD::FloatArray& crd)
    {
      E_Int ns[4], nods[4][4];
      K_MESH::Triangle::eDegenType ftype[4];
      E_Float FACTOR = 3;
      
      nods[0][0] = _nodes[0];
      nods[0][1] = _nodes[1];
      nods[0][2] = _nodes[2];
      nods[0][3] = _nodes[3];
      ftype[0] = K_MESH::Triangle::degen_type_angular(crd, nods[0][0], nods[0][1], nods[0][2], FACTOR, ns[0]);

      nods[1][0] = _nodes[0];
      nods[1][1] = _nodes[1];
      nods[1][2] = _nodes[3];
      nods[1][3] = _nodes[2];
      ftype[1] = K_MESH::Triangle::degen_type_angular(crd, nods[1][0], nods[1][1], nods[1][2], FACTOR, ns[1]);

      nods[2][0] = _nodes[1];
      nods[2][1] = _nodes[2];
      nods[2][2] = _nodes[3];
      nods[2][3] = _nodes[0];
      ftype[2] = K_MESH::Triangle::degen_type_angular(crd, nods[2][0], nods[2][1], nods[2][2], FACTOR, ns[2]);

      nods[3][0] = _nodes[2];
      nods[3][1] = _nodes[0];
      nods[3][2] = _nodes[3];
      nods[3][3] = _nodes[1];
      ftype[3] = K_MESH::Triangle::degen_type_angular(crd, nods[3][0], nods[3][1], nods[3][2], FACTOR, ns[3]);

      E_Int NB_OK{ 0 }, NB_SPIKES{ 0 }, NB_HATS{ 0 };

      for (size_t k = 0; k < 4; ++k)
      {
        if (ftype[k] == K_MESH::Triangle::eDegenType::OK)
          ++NB_OK;
        else if (ftype[k] == K_MESH::Triangle::eDegenType::SPIKE)
          ++NB_SPIKES;
        else if (ftype[k] == K_MESH::Triangle::eDegenType::HAT)
          ++NB_HATS;
      }

      if (NB_OK == 4)                           return REGULAR;
      if (NB_HATS == 2 && NB_SPIKES == 2)
      {
        // reorder _nodes to have the first two with the mel => WARNING : UNORIENTED
        E_Int n[4], count(0);
        if (ftype[0] == K_MESH::Triangle::eDegenType::HAT)
        {
          n[count++] = nods[0][ns[0]];
        }
        if (ftype[1] == K_MESH::Triangle::eDegenType::HAT)
        {
          n[count++] = nods[1][ns[1]];
        }
        if (ftype[2] == K_MESH::Triangle::eDegenType::HAT)
        {
          n[count++] = nods[2][ns[2]]; 
        }
        if (ftype[3] == K_MESH::Triangle::eDegenType::HAT)
        {
          n[count++] = nods[3][ns[3]];
        }
        assert(count == 2);

        // fills remaining
        for (size_t k = 0; k < 4; ++k)
          if (_nodes[k] != n[0] && _nodes[k] != n[1])
            n[count++] = _nodes[k];

        assert(count == 4);

        for (size_t k = 0; k < 4; ++k)
          _nodes[k] = n[k];

        return KNIFE2;
      }
      if (NB_HATS == 1 /*&& NB_SPIKES == 0*/)
      {
        if (ftype[0] == K_MESH::Triangle::eDegenType::HAT)
        {
          _nodes[0] = nods[0][ns[0]];           // reflex node
          _nodes[1] = nods[0][(ns[0] + 1) % 3]; // first node of opposite edge to split
          _nodes[2] = nods[0][(ns[0] + 2) % 3]; // second node of opposite edge to split
          _nodes[3] = nods[0][3];               // last tet node
        }
        else if (ftype[1] == K_MESH::Triangle::eDegenType::HAT)
        {
          _nodes[0] = nods[1][ns[1]];           // reflex node
          _nodes[1] = nods[1][(ns[1] + 1) % 3]; // first node of opposite edge to split
          _nodes[2] = nods[1][(ns[1] + 2) % 3]; // second node of opposite edge to split
          _nodes[3] = nods[1][3];               // last tet node
        }
        else if (ftype[2] == K_MESH::Triangle::eDegenType::HAT)
        {
          _nodes[0] = nods[2][ns[2]];           // reflex node
          _nodes[1] = nods[2][(ns[2] + 1) % 3]; // first node of opposite edge to split
          _nodes[2] = nods[2][(ns[2] + 2) % 3]; // second node of opposite edge to split
          _nodes[3] = nods[2][3];               // last tet node
        }
        else if (ftype[3] == K_MESH::Triangle::eDegenType::HAT)
        {
          _nodes[0] = nods[3][ns[3]];           // reflex node
          _nodes[1] = nods[3][(ns[3] + 1) % 3]; // first node of opposite edge to split
          _nodes[2] = nods[3][(ns[3] + 2) % 3]; // second node of opposite edge to split
          _nodes[3] = nods[3][3];               // last tet node
        }

        return DELTA;
      }
      if (NB_HATS >= 3)
      {
        //todo

        return KNIFE1;
      }

      if (NB_SPIKES >= 2) // SLICE1, SLICE2 or SPIKE ?
      {
        const E_Float* P0 = crd.col(_nodes[0]);
        const E_Float* P1 = crd.col(_nodes[1]);
        const E_Float* P2 = crd.col(_nodes[2]);
        const E_Float* P3 = crd.col(_nodes[3]);

        std::pair<E_Float, int> palma[4];

        palma[0] = std::make_pair(K_MESH::Triangle::surface<3>(P0, P1, P2), 1);
        palma[1] = std::make_pair(K_MESH::Triangle::surface<3>(P0, P1, P3), 2);
        palma[2] = std::make_pair(K_MESH::Triangle::surface<3>(P2, P0, P3), 4);
        palma[3] = std::make_pair(K_MESH::Triangle::surface<3>(P1, P2, P3), 8);

        std::sort(palma, palma + 4);

        if (FACTOR * palma[0].first < palma[1].first) // smin << 3 others
        {
          // reorder _nodes to have the first three to collapse => WARNING : UNORIENTED
          if (palma[0].second == 2) std::swap(_nodes[2], _nodes[3]);
          else if (palma[0].second == 4) std::swap(_nodes[1], _nodes[3]);
          else if (palma[0].second == 8) std::swap(_nodes[0], _nodes[3]);

          return SPIKE;
        }

        bool is_slice1 = (NB_OK == 2 && NB_SPIKES == 2) || (FACTOR * palma[1].first < palma[2].first); // smin1 & smin2 << 2 others

        if (is_slice1) 
        {
          // reorder _nodes to have the first two to collapse => WARNING : UNORIENTED

          if (palma[0].second + palma[1].second == 5) // N0 & N2 to collapse
          {
            std::swap(_nodes[1], _nodes[2]);
          }
          else if (palma[0].second + palma[1].second == 9) // N1 & N2 to collapse
          {
            std::swap(_nodes[0], _nodes[2]);
          }
          else if (palma[0].second + palma[1].second == 6) // N0 & N3 to collapse
          {
            std::swap(_nodes[1], _nodes[3]);
          }
          else if (palma[0].second + palma[1].second == 10) // N1 & N3 to collapse
          {
            std::swap(_nodes[0], _nodes[3]);
          }
          else if (palma[0].second + palma[1].second == 12) // N2 & N3 to collapse
          {
            K_CONNECT::IdTool::right_shift<4>(_nodes, 2);
          }

          return SLICE1;
        }

        // SLICE 2
        // reorder _nodes to have the first two to collapse , the third and fourth to collapse too => WARNING : UNORIENTED
        if (NB_SPIKES != 4)
          std::cout << "LOGIC ERROR" << std::endl;

        int count = 0; //fist 2 spike should be enough
        int pairs[2][2];
        pairs[0][0] = -1;
        pairs[0][1] = -1;
        pairs[1][0] = -1;
        pairs[1][1] = -1;
        if (ftype[0] == K_MESH::Triangle::eDegenType::SPIKE)
        {
          pairs[count][0] = nods[0][(ns[0] + 1) % 3]; // the two nodes but ns[0]
          pairs[count++][1] = nods[0][(ns[0] + 2) % 3];
        }
        if (ftype[1] == K_MESH::Triangle::eDegenType::SPIKE)
        {
          pairs[count][0] = nods[1][(ns[1] + 1) % 3]; // the two nodes but ns[0]
          pairs[count++][1] = nods[1][(ns[1] + 2) % 3];
        }
        if (ftype[2] == K_MESH::Triangle::eDegenType::SPIKE && count < 2)
        {
          pairs[count][0] = nods[2][(ns[2] + 1) % 3]; // the two nodes but ns[0]
          pairs[count++][1] = nods[2][(ns[2] + 2) % 3];
        }
        if (ftype[3] == K_MESH::Triangle::eDegenType::SPIKE && count < 2)
        {
          pairs[count][0] = nods[3][(ns[3] + 1) % 3]; // the two nodes but ns[0]
          pairs[count++][1] = nods[3][(ns[3] + 2) % 3];
        }

        assert(count == 2);
        assert(pairs[0][0] != pairs[1][0]);
        assert(pairs[0][0] != pairs[1][1]);
        assert(pairs[0][1] != pairs[1][0]);
        assert(pairs[0][1] != pairs[1][1]);

        _nodes[0] = pairs[0][0];
        _nodes[1] = pairs[0][1];
        _nodes[2] = pairs[1][0];
        _nodes[3] = pairs[1][1];

        return SLICE2;
      }

      //if (NB_SPIKES == 4) return SLICE2;
      //if (NB_SPIKES >= 3) return SPIKE;
      //if (NB_SPIKES == 2) return SLICE1;
      
      return REGULAR;
    }

    template <typename ngunit_t>
    inline static E_Int get_opposite_face_to_node(const ngunit_t & PGs, const E_Int* first_pg, E_Int N)
    {
      // Fopp is the face not having N

      for (size_t f = 0; f < 4; ++f)
      {
        E_Int Fi = first_pg[f] - 1;
        const E_Int* nodes = PGs.get_facets_ptr(Fi);
        //E_Int nnodes = PGs.stride(Fi);
        bool hasN = false;
        for (size_t n = 0; n < 3; ++n)
        {
          if (N == (nodes[n] - 1))
          {
            hasN = true;
            break;
          }
        }

        if (!hasN) return Fi;
      }

      return IDX_NONE;
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
    
    //dummies for generecity
    ///
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& acrd) {} //dummy : for genericity
    template <typename acrd_t>
    E_Int cvx_triangulate (const acrd_t& acrd) {return 0;}
        
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
    {return (NUGA::zzdet4(p1, p2, p3, p4)/6.);}
    
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
      {bb.minB[i] = NUGA::FLOAT_MAX; bb.maxB[i] = -NUGA::FLOAT_MAX;}

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
      
      NUGA::sum<3>(0.5, crd.col(v0), 0.5, crd.col(v2), crd.col(v4));
      NUGA::sum<3>(0.5, crd.col(v1), 0.5, crd.col(v2), crd.col(v5));
      NUGA::sum<3>(0.5, crd.col(v0), 0.5, crd.col(v1), crd.col(v6));
      
      NUGA::sum<3>(0.5, crd.col(v2), 0.5, crd.col(v3), crd.col(v7));
      NUGA::sum<3>(0.5, crd.col(v1), 0.5, crd.col(v3), crd.col(v8));
      NUGA::sum<3>(0.5, crd.col(v0), 0.5, crd.col(v3), crd.col(v9));
      
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
      NUGA::diff<3>(X, X0, X0X);
      gradf = NUGA::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      
      eval(p0,p1,p2,p3, bb, aa, aa, aa, X);
      NUGA::diff<3>(X, X0, X0X);
      gradf = NUGA::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      
      eval(p0,p1,p2,p3, aa, bb, aa, aa, X);
      NUGA::diff<3>(X, X0, X0X);
      gradf = NUGA::dot<3>(X0X, grad0);
      f = f0 + gradf;
      
      val +=f;
      
      eval(p0,p1,p2,p3, aa, aa, bb, aa, X);
      NUGA::diff<3>(X, X0, X0X);
      gradf = NUGA::dot<3>(X0X, grad0);
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
  E_Int F1Id(IDX_NONE), F2Id(IDX_NONE), F3Id(IDX_NONE);

  bool commonNodes[3];

  for (int k = 1; k < 4; ++k)
  {
    int count = 0;
    commonNodes[0] = commonNodes[1] = commonNodes[2] = false;

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

  assert (F1Id != IDX_NONE && F1Id != 0 && F1Id != F2Id && F1Id != F3Id);
  assert (F2Id != IDX_NONE && F2Id != 0 && F2Id != F1Id && F2Id != F3Id);
  assert (F3Id != IDX_NONE && F3Id != 0 && F3Id != F1Id && F3Id != F2Id);

  for (int i = 0; i < nb_faces; ++i)
    faces[i] = mol[i];
}

}
#endif	/* __K_MESH_TETRAHEDRON_H__ */

