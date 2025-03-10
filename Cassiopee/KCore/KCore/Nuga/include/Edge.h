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

#ifndef __K_MESH_EDGE_H__
#define __K_MESH_EDGE_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/maths.hxx"
#include "Nuga/include/DynArray.h"
#include <vector>
#include <algorithm>

namespace K_MESH
{

class Edge
{
public:
  typedef E_Int         size_type;
  typedef E_Int         node_type;
  typedef E_Int         boundary_type; // An edge is bounded by 2 nodes.
  typedef Edge          self_type;

public:
  static constexpr E_Int NB_NODES = 2;

public: /* Constructors, Destructor and operators*/

  ///
  explicit Edge(node_type n0, node_type n1){_nodes[0] = n0; _nodes[1] = n1;}
  ///
  explicit Edge(const E_Int* pN){_nodes[0] = *pN; _nodes[1] = *(pN+1);}
  ///
  explicit Edge(const K_FLD::IntArray& cnt, E_Int i){_nodes[0] = cnt(0,i); _nodes[1] = cnt(1,i);}
  ///
  explicit Edge(void){_nodes[0] = _nodes[1] = IDX_NONE;}
  ///
  virtual ~Edge(void){}

  /// Assignment operator.
  //Edge& operator=(const self_type& E){_nodes[0] = E._nodes[0]; _nodes[1] = E._nodes[1]; return *this;}

  /// Equality operator.
  inline bool operator==(const self_type& E)
  {return ((_nodes[0] == E._nodes[0]) && (_nodes[1] == E._nodes[1]));}

  /// Inequality operator.
  inline bool operator <(const self_type & E) const
  {return (_nodes[0] < E._nodes[0]) || (!(E._nodes[0] < _nodes[0]) && (_nodes[1] < E._nodes[1]));}

public: /* Set and Get methods */
  
  /// Sets the edge nodes.
  template <typename NodeIterator>
  inline void setNodes(const NodeIterator pN){_nodes[0] = *pN; _nodes[1] = *(pN+1);}
  
  /// Set the edge nodes.
  virtual inline void setNodes(node_type N0, node_type N1){_nodes[0] = N0; _nodes[1] = N1;}

  /// Gets the i-th node.
  inline const node_type& node(const size_type& i) const {return _nodes[i];}
  
  ///
  inline E_Int* begin() { return &_nodes[0];}
  inline E_Int* end() { return &_nodes[0]+NB_NODES;}
  inline const E_Int* begin() const { return &_nodes[0];}
  inline const E_Int* end() const { return &_nodes[0]+NB_NODES;}

  template <typename box_t, typename crd_t>
  void bbox (const crd_t&crd, box_t& box) const
  {
    box = box_t(crd, _nodes, NB_NODES);
  }

  double Lref2(const K_FLD::FloatArray& crd) const
  {
    return NUGA::sqrDistance(crd.col(_nodes[0]), crd.col(_nodes[1]), 3);
  }
  double Lref2(const std::vector<E_Float>& nodal_tol2) const
  {
    return std::min(nodal_tol2[_nodes[0]], nodal_tol2[_nodes[1]]);
  }

  /// Gets the i-th node (required to abstract some algorithms to be valid both in 2D and 3D.
  inline void getBoundary(size_type n, boundary_type& b){b = _nodes[n];}
  static inline void getBoundary(const E_Int* pN, size_type n, boundary_type& b){b = *(pN+n);}

public : /* Static functions related to edges */

  ///
  static void getBoundary(const self_type& E1, const self_type& E2, boundary_type& b);

  ///Computes the minimum distance between the lines supported by P0P1 and Q0Q1.
  /** u0 and u1 are the parameter along each line to the points closest to each other.
      in case of intersection both parameters allow to compute the intersection point and min_distance = 0.*/
  template <short DIM>
  static void lineLineMinDistance
              (const E_Float* P0, const E_Float* P1, const E_Float* Q0, const E_Float* Q1,
               E_Float& u0, E_Float& u1,
               E_Float tolerance, E_Bool& parallel, E_Bool& coincident, E_Float& min_distance);
                               
  ///
  template <short DIM>
  static bool intersect (const E_Float* P0, const E_Float* P1, const E_Float* Q0, const E_Float* Q1,
                         E_Float tol, E_Bool tol_is_absolute,
                         E_Float& u00, E_Float& u01, E_Float& u10, E_Float& u11, E_Bool& overlap);

  ///
  template <short DIM>
  static bool intersect 
    (const K_FLD::FloatArray& pos, E_Int N0, E_Int N1, E_Int M0, E_Int M1, E_Float tol, E_Bool tol_is_absolute,
    E_Float& u00, E_Float& u01, E_Float& u10, E_Float& u11, E_Bool& overlap, bool enforceCoplanarity=false);
  
  ///
  template <short DIM>
  static E_Float
  linePointMinDistance(const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda) ;
  
  ///
  template <short DIM>
  static E_Float
  linePointMinDistance2(const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda) ;

  ///
  template <short DIM>
  static E_Float
  edgePointMinDistance (const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda);
  
  ///
  template <short DIM>
  inline static E_Float
  edgePointMinDistance2 (const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda);

public : 
  static bool _enforce_coplanarity;
protected:
  E_Int _nodes[2];
};

////////////////

class NO_Edge : public Edge // Non oriented edge.
{
public:
  NO_Edge(E_Int n0, E_Int n1): Edge(n0 < n1 ? n0 : n1, n0 < n1 ? n1 : n0){}
  NO_Edge(const E_Int* pN): Edge(*pN < *(pN+1) ? *pN : *(pN+1), *pN < *(pN+1) ? *(pN+1) : *pN){}
  NO_Edge(void): Edge(){}
  inline void setNodes(E_Int N0, E_Int N1){_nodes[0] = (N0 < N1) ? N0 : N1; _nodes[1] = (N0 < N1) ? N1 : N0;}
  template <typename NodeIterator>
  inline void setNodes(const NodeIterator pN){_nodes[0] = (*pN < *(pN+1)) ? *pN : *(pN+1); _nodes[1] = (*pN < *(pN+1)) ? *(pN+1) : *pN;}
  ~NO_Edge(){}
};

////////////////
struct aEdge : public Edge
{
  E_Float m_Lref2;
  
  aEdge(const Edge& e, const K_FLD::FloatArray& crd):m_Lref2(-1.)
  {
    v1[0] = crd(0, e.node(0));
    v1[1] = crd(1, e.node(0));
    v1[2] = crd(2, e.node(0));
    
    v2[0] = crd(0, e.node(1));
    v2[1] = crd(1, e.node(1));
    v2[2] = crd(2, e.node(1));
  } // from "mesh" to autonmous
  
  aEdge(const Edge& e, const K_FLD::FloatArray& crd, E_Float L2r):aEdge(e, crd){m_Lref2 = L2r;}

  aEdge(const K_FLD::IntArray& cnt, E_Int i, const K_FLD::FloatArray& crd)
  {
    const E_Float* pt1 = crd.col(cnt(0, i));
    const E_Float* pt2 = crd.col(cnt(1, i));

    v1[0] = pt1[0];
    v1[1] = pt1[1];
    v1[2] = pt1[2];

    v2[0] = pt2[0];
    v2[1] = pt2[1];
    v2[2] = pt2[2];

    m_Lref2 = Lref2();
  }

  double Lref2() const { return (m_Lref2 > 0.) ? m_Lref2 : NUGA::sqrDistance(v1, v2, 3);}
 
  E_Float v1[3], v2[3];
};

} // End namespace K_MESH

//// Template implementation //////////
template <short DIM>
bool
K_MESH::Edge::intersect 
(const E_Float* P0, const E_Float* P1, const E_Float* Q0, const E_Float* Q1,
 E_Float tol, E_Bool tol_is_absolute,
 E_Float& u00, E_Float& u01, E_Float& u10, E_Float& u11, E_Bool& overlap) 
{
  E_Float u0, u1, min_d, tol0(tol), tol1(tol);
  E_Bool  parallel, coincident;
  
  overlap = false;
  u00 = u01 = u10 = u11 = NUGA::FLOAT_MAX;

  lineLineMinDistance<DIM>(P0, P1, Q0, Q1, u0, u1, tol, parallel, coincident, min_d);

  E_Float E0[DIM], E1[DIM];
  NUGA::diff<DIM>(P1, P0, E0);
  E_Float L0 = NUGA::sqrNorm<DIM>(E0);
  L0 = 1. / L0;
  NUGA::diff<DIM>(Q1, Q0, E1);
  E_Float L1 = NUGA::sqrNorm<DIM>(E1);
  L1 = 1. / L1;
  E_Float l0 = ::sqrt(L0);
  E_Float l1 = ::sqrt(L1);

  if (tol_is_absolute)
  {
    tol0 /= l0;
    tol1 /= l1;
  }

  if (min_d > std::max(tol0 * l0, tol1 * l1))      // Agonic lines
    return false;
  
  if (!coincident) // The lines have exactly one point of intersection.
  {
    u00 = u0;
    u10 = u1;

    return ((u0 >= -tol0) && (u0 <= (1. + tol0)) && (u1 >= -tol1) && (u1 <= (1. + tol1)));
  }
  else                  // The support line is the same for the 2 edges : up to 2 point of intersection.
  {
   E_Float V00[DIM], V01[DIM];   
   
   NUGA::diff<DIM>(Q0, P0, V00);
   NUGA::diff<DIM>(Q1, P0, V01);

   E_Float s00 = NUGA::dot<DIM>(E0, V00) * L0; // fixme FIXME FIXME a revoir
   E_Float s01 = NUGA::dot<DIM>(E0, V01) * L0;
   if (s01 < s00)
     std::swap(s00, s01);

   if ( (s01 < -tol0) || (s00 > (1. + tol0)) ) //  x----------x    x------x    
     return false;
   
   NUGA::diff<DIM>(P1, Q0, V01);

   E_Float s10 = - NUGA::dot<DIM>(E1, V00) * L1;
   E_Float s11 = NUGA::dot<DIM>(E1, V01) * L1;
   if (s11 < s10)
     std::swap(s10, s11);

   u00 = (s00 < (1. - tol0)) ? s00 : 1.; // return in range [0.;1.]
   u00 = (u00 > tol0) ? u00 : 0.;        //
   u01 = (s01 < (1. - tol0)) ? s01 : 1.; //
   u01 = (u01 > tol0) ? u01 : 0.;        //
   u10 = (s10 < (1. - tol1)) ? s10 : 1.; //
   u10 = (u10 > tol1) ? u10 : 0.;        //
   u11 = (s11 < (1. - tol1)) ? s11 : 1.; //
   u11 = (u11 > tol1) ? u11 : 0.;        //

   if (u00 == u01)                      //  x----------x------x
   {
     //u00 = (::fabs(s00) < ::fabs(s01)) ? s00 : s01;
     u01 = NUGA::FLOAT_MAX;
     //u10 = (::fabs(s10) < ::fabs(s11)) ? s10 : s11;
     u11 = NUGA::FLOAT_MAX;
   }
   else                                 //  x-------|--x      or     x-----x     or  x-------x
     overlap = true;                    //          x--|---x      x--|-----|--x      x-------x

   return true;
  }
}

//// Template implementation //////////
template <short DIM>
bool
K_MESH::Edge::intersect 
(const K_FLD::FloatArray& pos, E_Int N0, E_Int N1, E_Int M0, E_Int M1, E_Float tol, E_Bool tol_is_absolute,
 E_Float& u00, E_Float& u01, E_Float& u10, E_Float& u11, E_Bool& overlap, bool enforceCoplanarity) 
{
  E_Float u0, u1, min_d, tol0(tol), tol1(tol);
  E_Bool  parallel, coincident, sharing_a_node;

  enforceCoplanarity |= (DIM == 2); // 2D => true

  const E_Float* P0 = pos.col(N0);
  const E_Float* P1 = pos.col(N1);
  const E_Float* Q0 = pos.col(M0);
  const E_Float* Q1 = pos.col(M1);
  
  overlap = false;
  u00 = u01 = u10 = u11 = NUGA::FLOAT_MAX;

  lineLineMinDistance<DIM>(P0, P1, Q0, Q1, u0, u1, tol, parallel, coincident, min_d);

  sharing_a_node = ((N0 == M0) || (N0 == M1) || (N1 == M0) || (N1 == M1));

  //special case : sharing a node
  {
    if (sharing_a_node)
    {
      min_d = 0.;
      coincident = parallel;
    }
  }
  
  if (enforceCoplanarity && !parallel)
    min_d=0.;

  E_Float E0[DIM], E1[DIM];
  NUGA::diff<DIM>(P1, P0, E0);
  E_Float L0 = NUGA::sqrNorm<DIM>(E0);
  L0 = 1. / L0;
  NUGA::diff<DIM>(Q1, Q0, E1);
  E_Float L1 = NUGA::sqrNorm<DIM>(E1);
  L1 = 1. / L1;
  E_Float l0 = ::sqrt(L0);
  E_Float l1 = ::sqrt(L1);

  if (tol_is_absolute)
  {
    tol0 *= l0; // * because l0 is an inverse distance
    tol1 *= l1;
  }
  
  // tol0 & tol1 are now relative

  if (min_d > std::max(tol0 * l0, tol1 * l1))      // Agonic lines
    return false;

  if (!coincident) // The lines have exactly one point of intersection.
  {
    if (sharing_a_node)
    {
      if (N0 == M0)
      {
        u0=u1=0.;
      }
      else if (N0 == M1)
      {
        u0=0.; u1=1.;
      }
      else if (N1 == M0)
      {
        u0=1.; u1=0.;
      }
      else if (N1 == M1)
      {
        u0=1.; u1=1.;
      }
    }

    u00 = u0;
    u10 = u1;

    return ((u0 >= -tol0) && (u0 <= (1. + tol0)) && (u1 >= -tol1) && (u1 <= (1. + tol1)));
  }
  else                  // The support line is the same for the 2 edges : up to 2 point of intersection.
  {
    E_Float V00[DIM], V01[DIM];   
   
    NUGA::diff<DIM>(Q0, P0, V00);
    NUGA::diff<DIM>(Q1, P0, V01);

    E_Float s00 = NUGA::dot<DIM>(E0, V00) * L0; //pos de Q0 sur E0
    E_Float s01 = NUGA::dot<DIM>(E0, V01) * L0; //pos de Q1 sur E0
    if (s01 < s00)
      std::swap(s00, s01);

    if ( (s01 < -tol0) || (s00 > (1. + tol0)) ) //  x----------x    x------x    
      return false;

    NUGA::diff<DIM>(P1, Q0, V01);

    E_Float s10 = -NUGA::dot<DIM>(E1, V00) * L1; //pos de P0 sur E1
    E_Float s11 = NUGA::dot<DIM>(E1, V01) * L1;  //pos de P1 sur E1
    if (s11 < s10)
      std::swap(s10, s11);

    u00 = (s00 < (1. - tol0)) ? s00 : 1.; // return in range [0.;1.]
    u00 = (u00 > tol0) ? u00 : 0.;        //
    u01 = (s01 < (1. - tol0)) ? s01 : 1.; //
    u01 = (u01 > tol0) ? u01 : 0.;        //
    u10 = (s10 < (1. - tol1)) ? s10 : 1.; //
    u10 = (u10 > tol1) ? u10 : 0.;        //
    u11 = (s11 < (1. - tol1)) ? s11 : 1.; //
    u11 = (u11 > tol1) ? u11 : 0.;        //

    if (u00 == u01)                      //  x----------x------x
    {
      //u00 = (::fabs(s00) < ::fabs(s01)) ? s00 : s01;
      u01 = NUGA::FLOAT_MAX;
      //u10 = (::fabs(s10) < ::fabs(s11)) ? s10 : s11;
      u11 = NUGA::FLOAT_MAX;
    }
    else                                 //  x-------|--x      or     x-----x     or  x-------x
      overlap = true;                    //          x--|---x      x--|-----|--x      x-------x

    return true;
  }
}

#define MAX(a,b) ((a<b)?b:a)
#define MIN(a,b) ((a<b)?a:b)
///
template <short DIM>
void
K_MESH::Edge::lineLineMinDistance
(const E_Float* P0, const E_Float* P1, const E_Float* Q0, const E_Float* Q1,
 E_Float& u0, E_Float& u1,
 E_Float abstol, E_Bool& parallel, E_Bool& coincident, E_Float& min_distance)
{
  E_Float E0[DIM], E1[DIM], V01[DIM];
  E_Float L02, L12, tol2(abstol*abstol);
  
  u0 = u1 = min_distance = NUGA::FLOAT_MAX;
  coincident = parallel = false;

  NUGA::diff<DIM>(P1, P0, E0);
  NUGA::diff<DIM>(Q1, Q0, E1);
  NUGA::diff<DIM>(P0, Q0, V01);
  L02 = NUGA::sqrNorm<DIM>(E0);
  L12 = NUGA::sqrNorm<DIM>(E1);
  
  if ((L02 < tol2) || (L12 < tol2)) // skip degen
    return;

  E_Float l1, l2;
  E_Float d1 = linePointMinDistance2<DIM>(P0,P1, Q0, l1);
  E_Float d2 = linePointMinDistance2<DIM>(P0,P1, Q1, l2);
  E_Float dmax12 = MAX(d1, d2);
  //E_Float dmin12 = MIN(d1, d2);
    
  if (dmax12 < tol2)
  {
    coincident=parallel=true;
    min_distance=0.;
    return;
  }
    
  d1 = linePointMinDistance2<DIM>(Q0,Q1, P0, l1);
  d2 = linePointMinDistance2<DIM>(Q0,Q1, P1, l2);
  E_Float dmax22 = MAX(d1, d2);
  //E_Float dmin22 = MIN(d1, d2);
  
  if (dmax22 < tol2)
  {
    coincident=parallel=true;
    min_distance=0.;
    return;
  }
      
  //At this stage, computing the determinant should be safe, ie. no parallelism
  E_Float a, b, c, det2, u,v=0; //u,v are internal temporaries
  
  det2 = NUGA::sqrCross<DIM>(E0, E1);
  c = NUGA::dot<DIM>(E1, V01); 
  parallel = (::fabs(det2) < tol2 * L02 * L12);
  
  if (DIM == 2) // more testing for robsutness
  {
  	E_Float normal_Q0Q1[] = { -E1[1], E1[0] };
    double pow2D_P0_Q0Q1 = NUGA::dot<2>(V01, normal_Q0Q1);
  
    E_Float V02[DIM];
    NUGA::diff<DIM>(P1, Q0, V02);
    E_Float pow2D_P1_Q0Q1 = NUGA::dot<2>(V02, normal_Q0Q1);

    parallel &= ((pow2D_P0_Q0Q1 * pow2D_P1_Q0Q1) > 0.);
  }
  

  if (!parallel)
  {
    if (DIM == 3)
    {
      a = NUGA::dot<DIM>(E0, E1);
      b = NUGA::dot<DIM>(E0, V01);
      det2 = 1. / det2;
      u = det2 * (a*c - L12 * b);
      v = det2 * (L02*c - a * b);
      E_Float IJ[DIM];
      for (E_Int i = 0; i < DIM; ++i)
        IJ[i] = V01[i] + u * E0[i] - v * E1[i];
      min_distance = NUGA::sqrNorm<DIM>(IJ);

      if (min_distance > MIN(dmax12, dmax22)) // numerical error catch (nearly //). fixme : SHOULD BE dmin12 and dmin22. cf. XXX
      {
        parallel = true;
        u = 0.;
        v = c / L12;
      }
      else // params are good
      {
        u0 = u; u1 = v;
      }
    }
    else // 2D
    {
      E_Float normal_Q0Q1[] = { -E1[1], E1[0] };
      NUGA::normalize<2>(normal_Q0Q1);
      double a = -NUGA::dot<2>(E0, normal_Q0Q1);
      double pow2D_P0_Q0Q1 = NUGA::dot<2>(V01, normal_Q0Q1);
      double K = pow2D_P0_Q0Q1 / a;

      double X[2];
      for (E_Int i = 0; i < 2; ++i)
        X[i] = P0[i] + K * E0[i];

      d1 = linePointMinDistance2<DIM>(P0, P1, X, u0);
      d2 = linePointMinDistance2<DIM>(Q0, Q1, X, u1);
    }
  }
  else //parallel
  {
    u = 0.;
    v = c / L12;
  }

  if (parallel || (DIM != 2))
  {
    E_Float IJ[DIM];
    for (E_Int i = 0; i < DIM; ++i)
      IJ[i] = V01[i] + u * E0[i] - v * E1[i];
    min_distance = ::sqrt(NUGA::sqrNorm<DIM>(IJ));
  }
  else
    min_distance = 0.;

  coincident = (parallel && (min_distance < abstol));
}

//=============================================================================
template <short DIM>
E_Float
K_MESH::Edge::linePointMinDistance2
(const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda)
{
  E_Float       V0[DIM], V1[DIM], L;
 
  NUGA::diff<DIM> (P1, P0, V0);
  NUGA::diff<DIM> (P, P0, V1);
  lambda = 0.;

  L = NUGA::normalize<DIM>(V0);
  lambda = NUGA::dot<DIM> (V0, V1) / L;
  return NUGA::sqrCross<DIM>(V0, V1);
}

//=============================================================================
template <short DIM>
E_Float
K_MESH::Edge::linePointMinDistance
(const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda) 
{
  return ::sqrt(linePointMinDistance2<DIM>(P0,P1,P,lambda));
}

//=============================================================================
template <short DIM>
E_Float
K_MESH::Edge::edgePointMinDistance
(const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda) 
{
  return ::sqrt(edgePointMinDistance2<DIM>(P0,P1,P,lambda));
}

//=============================================================================
template <short DIM>
E_Float
K_MESH::Edge::edgePointMinDistance2
(const E_Float* P0, const E_Float* P1, const E_Float* P, E_Float& lambda) 
{
  E_Float       V0[DIM], V1[DIM], L, d2;
 
  NUGA::diff<DIM> (P1, P0, V0);
  NUGA::diff<DIM> (P, P0, V1);
  lambda = 0.;

  L = NUGA::normalize<DIM>(V0);
  lambda = NUGA::dot<DIM> (V0, V1) / L;

  if (lambda < 0.)      // Closest point is P0
    d2 = NUGA::sqrDistance(P0, P, DIM);
  else if (lambda > 1.) // Closest point is P1
    d2 = NUGA::sqrDistance(P1, P, DIM);
  else                  // The projection of P is inside the edge so d is the distance point-edge.
    d2 = NUGA::sqrCross<DIM>(V0, V1);

  return d2;
}

#endif
