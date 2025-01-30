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

#include "Nuga/include/Triangle.h"
#include <limits>

const E_Int K_MESH::Triangle::NB_NODES = 3;
const E_Int K_MESH::Triangle::NB_TRIS = 1;
//=============================================================================
E_Float K_MESH::Triangle::surface(const E_Float* p1, const E_Float* p2, 
                                  const E_Float* p3, size_type dim)
{
  if (dim == 2)
    return surface<2>(p1, p2, p3);
  if (dim == 3)
    return surface<3>(p1, p2, p3);

  return NUGA::FLOAT_MAX;
}

//=============================================================================
NUGA::size_type
K_MESH::Triangle::getOppLocalNodeId
(NUGA::size_type K, NUGA::size_type n,
 const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors)
{
  assert (n < NB_NODES);
  assert (K < neighbors.cols());

  size_type Kadj = neighbors(n, K);

  if (Kadj == IDX_NONE) return IDX_NONE;

  K_FLD::IntArray::const_iterator pK = connect.col(K);
  size_type N = *(pK+(n+1)%NB_NODES);
  K_FLD::IntArray::const_iterator pKadj = connect.col(Kadj);
  size_type a = getLocalNodeId(pKadj, N);

  return (a+1)%NB_NODES;
}

//=============================================================================
void K_MESH::Triangle::getBoundary
(const Triangle&  T1, const Triangle&  T2, K_MESH::NO_Edge& b)
{
  K_MESH::NO_Edge Ei, Ej;
  for (E_Int i = 0; i < NB_NODES; ++i)
  {
    T1.getBoundary(i, Ei);
    for (E_Int j = 0; j < NB_NODES; ++j)
    {
      T2.getBoundary(j, Ej);
      if (Ei == Ej)
      {
        b = Ei; return;
      }
    }
  }
}

//=============================================================================
void K_MESH::Triangle::getBoundary
(const Triangle&  T1, const Triangle&  T2, E_Int& i1, E_Int& i2)
{
  K_MESH::NO_Edge Ei, Ej;
  for (i1=0; i1 < NB_NODES; ++i1)
  {
    T1.getBoundary(i1, Ei);
    for (i2=0; i2 < NB_NODES; ++i2)
    {
      T2.getBoundary(i2, Ej);
      if (Ei == Ej)
        return;
    }
  }
}

//=============================================================================
E_Int K_MESH::Triangle::getOrientation(const Triangle&  T, 
                                       const E_Int& Ni, const E_Int& Nj, 
                                       E_Bool& same_orient)
{
  same_orient = false;
 
  for (E_Int n = 0; n < NB_NODES; ++n)
  {
    if ((T._nodes[n] == Ni) && (T._nodes[(n+1)%NB_NODES] == Nj))
    {
      same_orient = true;
      return 0;
    }
    if ((T._nodes[n] == Nj) && (T._nodes[(n+1)%NB_NODES] == Ni))
      return 0;
  }

  return -1;
}

K_MESH::Triangle::eDegenType K_MESH::Triangle::degen_type(const K_FLD::FloatArray& crd, E_Int N0, E_Int N1, E_Int N2, E_Float tol2, E_Float lambdac, E_Int& ns)
{
  ns = IDX_NONE;// ns for "special node" : the pick for a spike, the hat node for a hat
  //
  E_Float normal[3];
  K_MESH::Triangle::normal(crd.col(N0), crd.col(N1), crd.col(N2), normal);
  E_Float l2 = ::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  //
  if (::fabs(l2 - 1.) < EPSILON)
    return OK;

  // we have a degen !

  std::pair<E_Float, E_Int> palma[3];
  E_Int N[] = { N0, N1, N2 };

  palma[0] = std::make_pair(NUGA::sqrDistance(crd.col(N0), crd.col(N1), 3), 2);
  palma[1] = std::make_pair(NUGA::sqrDistance(crd.col(N0), crd.col(N2), 3), 1);
  palma[2] = std::make_pair(NUGA::sqrDistance(crd.col(N1), crd.col(N2), 3), 0);
  
  std::sort(&palma[0], &palma[0] + 3);

  // the last is the biggest edge and therefore is the base of the triangle
  // Ntop is the node to project on the base
  E_Int n = palma[2].second;
  E_Int& Ntop = N[n];
  E_Float& Lbase = palma[2].first;

  //if the biggest is too small, all of them are and this element need to be collapsed
  // the factor 4 is for implicit small : if the projected point cut in 2 pieces smaller than tol <=> Lbase < 2* tol
  if (Lbase*Lbase < 4.*tol2) return SMALL; 

  E_Float lambda;
  //E_Float d = 
  K_MESH::Edge::edgePointMinDistance<3>(crd.col(N[(n + 1) % 3]), crd.col(N[(n + 2) % 3]), crd.col(Ntop), lambda);

  //assert(lambda >= -1.e-15 && lambda <= 1. + 1.e-15);

  if ( (lambda < lambdac) || (lambda*lambda*Lbase*Lbase < tol2) )
  {
    ns = (n + 2) % 3;
    return SPIKE;
  }

  if ((lambda > 1. - lambdac) || ((1. - lambda) * (1. - lambda) * Lbase * Lbase < tol2))
  {
    ns = (n + 1) % 3;
    return SPIKE;
  }
    
  ns = n;
  return HAT;
}

K_MESH::Triangle::eDegenType K_MESH::Triangle::degen_type_angular(const K_FLD::FloatArray& crd, E_Int N0, E_Int N1, E_Int N2, E_Float FACTOR, E_Int& ns)
{
  eDegenType type = OK;
  ns = IDX_NONE;

  const double* P0 = crd.col(N0);
  const double* P1 = crd.col(N1);
  const double* P2 = crd.col(N2);
  
  double P0P1[3], P1P2[3], P2P0[3];
  NUGA::diff<3>(P1, P0, P0P1);
  NUGA::diff<3>(P2, P1, P1P2);
  NUGA::diff<3>(P0, P2, P2P0);

  double a0 = NUGA::PI - NUGA::normals_angle(P2P0, P0P1);
  double a1 = NUGA::PI - NUGA::normals_angle(P0P1, P1P2);
  double a2 = NUGA::PI - NUGA::normals_angle(P1P2, P2P0);

  std::pair<double, int> palma[3];
  
  palma[0] = std::make_pair(a0, 0);
  palma[1] = std::make_pair(a1, 1);
  palma[2] = std::make_pair(a2, 2);

  std::sort(&palma[0], &palma[0] + 3);

  double ratio1 = palma[1].first / palma[0].first;
  double ratio2 = palma[2].first / palma[1].first;

  if (FACTOR * ratio1 < ratio2) // ratio1 << ratio2 : the greatest angle is far from the 2 others
  {
    ns = palma[2].second;
    return HAT;
  }
  else if (FACTOR * ratio2 < ratio1) //the smallest angle is far from the 2 others
  {
    ns = palma[0].second;
    return SPIKE;
  }
  return type;
}
