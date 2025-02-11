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

//#define DEBUG_SWAPPER

#include "Nuga/include/SwapperT3.h"
#include <map>
#include <set>
#include "Nuga/include/Edge.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/IdTool.h"
#ifdef DEBUG_SWAPPER
#include "Nuga/include/medit.hxx"
std::string wdir = "";
#endif

//
SwapperT3::SwapperT3(){}
SwapperT3::~SwapperT3(){}

//
SwapperT3::eDegenType SwapperT3::degen_type2(const K_FLD::FloatArray& crd, E_Int N0, E_Int N1, E_Int N2, E_Float tol2, E_Float lambdac, E_Int& ns)
{
  ns = IDX_NONE;// ns for "special node" : the pick for a spike, the hat node for a hat
  //
  E_Float normal[3], MINQUAL{ ZERO_M }, BADQUAL_RTOL{ 0.25 }; //bad qual => we need to make it vanish => has to fall into spike/hat/small so big tol..
  K_MESH::Triangle::normal(crd.col(N0), crd.col(N1), crd.col(N2), normal);
  E_Float l2 = ::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  //
  bool normal_failure = !(::fabs(l2 - 1.) < EPSILON);
  
  double q = K_MESH::Triangle::qualityG<3>(crd.col(N0), crd.col(N1), crd.col(N2));
  bool good_qual = (q > MINQUAL);

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
  E_Float& Lbase2 = palma[2].first;
  //E_Float& Lmin2 = palma[0].first;

  /*std::cout << "Lmin : " << ::sqrt(palma[0].first) << std::endl;
  std::cout << "Lmax : " << ::sqrt(palma[2].first) << std::endl;
  std::cout << " ratio : " << ::sqrt(palma[0].first / palma[2].first) << std::endl;*/

  E_Float lambda;
  E_Float d = K_MESH::Edge::edgePointMinDistance<3>(crd.col(N[(n + 1) % 3]), crd.col(N[(n + 2) % 3]), crd.col(Ntop), lambda);
  if (!normal_failure && (d*d > tol2) && good_qual) return OK;

  if (!good_qual) // increase the tolerance
    tol2 = BADQUAL_RTOL*Lbase2;

  //if the biggest is too small, all of them are and this element need to be collapsed
  // the factor 4 is for implicit small : if the projected point cut in 2 pieces smaller than tol <=> Lbase < 2* tol
  if (Lbase2 < 4.*tol2)
    return SMALL;

  //assert(lambda >= -1.e-15 && lambda <= 1. + 1.e-15);
  //std::cout << "lambda : " << lambda << std::endl;

  if ( (lambda < lambdac) || (lambda*lambda*Lbase2 < tol2) )
  {
    ns = (n + 2) % 3;
    return SPIKE;
  }

  if ((lambda > 1. - lambdac) || ((1. - lambda) * (1. - lambda) * Lbase2 < tol2))
  {
    ns = (n + 1) % 3;
    return SPIKE;
  }
    
  ns = n;
  return HAT;
}

//
E_Int SwapperT3::clean(const K_FLD::FloatArray& coord, E_Float tol, K_FLD::IntArray& connect, std::vector<E_Int> & coids, std::vector<E_Int>& cnids)
{

  E_Float tol2(tol*tol);
  E_Int nb_tris0;

  K_CONNECT::IdTool::init_inc(cnids, connect.cols());

  do
  {
    K_FLD::IntArray::const_iterator   pS;
    E_Int                             Si;
    std::vector<E_Int> oidc;

    nb_tris0 = connect.cols();

    std::vector<E_Int> nodnids;
    K_CONNECT::IdTool::init_inc(nodnids, coord.cols());

    for (Si = 0; Si < nb_tris0; ++Si)
    {
      //std::cout << Si << std::endl;

      pS = connect.col(Si);

      E_Int N0 = *pS;
      E_Int N1 = *(pS + 1);
      E_Int N2 = *(pS + 2);

      E_Int ni;
      eDegenType type = degen_type2(coord, N0, N1, N2, tol2, 0., ni);

      if (type == SMALL)
      {
        E_Int Ni = std::min(N0, std::min(N1, N2));
        nodnids[N0] = nodnids[N1] = nodnids[N2] = Ni;
      }
      else if (type == SPIKE)
      {
        E_Int N1 = *(pS + (ni + 1) % 3);
        E_Int N2 = *(pS + (ni + 2) % 3);
        E_Int Ni = std::min(N1, N2);
        nodnids[N1] = nodnids[N2] = Ni;
      }
    }

    K_FLD::IntArray::changeIndices(connect, nodnids); //node ids

    std::vector<E_Int> newIDs;
    if (remove_degen(connect, newIDs))// 
    {
      K_CONNECT::valid pred(newIDs);
      K_CONNECT::IdTool::compress(coids, pred);
      K_CONNECT::IdTool::propagate(newIDs, cnids);
    }
  }
  while (connect.cols() < nb_tris0);

  return 0;
}

bool SwapperT3::has_degen(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, E_Float tol2)
{
  for (E_Int i = 0; i < connect.cols(); ++i)
  {
    K_FLD::IntArray::const_iterator pS = connect.col(i);

    E_Int N0 = *pS;
    E_Int N1 = *(pS + 1);
    E_Int N2 = *(pS + 2);

    E_Int ni;
    eDegenType type = degen_type2(coord, N0, N1, N2, tol2, 0., ni);
    if (type != OK)
      return true;
  }
  return false;

}

bool SwapperT3::has_hat(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, E_Float tol2)
{
  for (E_Int i = 0; i < connect.cols(); ++i)
  {
    K_FLD::IntArray::const_iterator pS = connect.col(i);

    E_Int N0 = *pS;
    E_Int N1 = *(pS + 1);
    E_Int N2 = *(pS + 2);

    E_Int ni;
    eDegenType type = degen_type2(coord, N0, N1, N2, tol2, 0., ni);
    if (type == HAT)
      return true;
  }
  return false;

}

///
E_Int SwapperT3::remove_degen
(K_FLD::IntArray& connect, NUGA::int_vector_type& newIDs)
{
  E_Int                       Si, COLS(connect.cols()), ROWS(connect.rows());
  K_FLD::IntArray::iterator   pS;
  K_FLD::IntArray             connectOut;
  
  connectOut.reserve(ROWS, COLS);

  newIDs.clear();
  newIDs.resize(COLS, IDX_NONE);

  for (Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);

    if ((*pS != *(pS+1)) && (*pS != *(pS+2)) && (*(pS+1) != *(pS+2)))
    {
      connectOut.pushBack(pS, pS + ROWS);
      newIDs[Si] = connectOut.cols() - 1;
    }
  }

  connect = connectOut;

  return (COLS - connect.cols());
}

//
E_Int SwapperT3::run (const K_FLD::FloatArray& coord, E_Float tol, K_FLD::IntArray& connect, std::vector<E_Int> & oids)
{
  K_FLD::IntArray::const_iterator   pS;
  E_Int                             Si, nb_tris0(connect.cols()), n, N0, N1, N2;
  //Ni, NB_NODES(3);
  //E_Float                           d2, l0, d0, lambda, Pk[3], dMax;
  E_Float tol2(tol*tol), lambdac=0.0;
  K_MESH::NO_Edge                   Ei;

  typedef std::map<K_MESH::NO_Edge, std::vector<E_Int> > e_to_pts_t;
  e_to_pts_t E_to_Pts;
  e_to_pts_t::iterator itE;
  //E_Float tolx = -1.;
  std::set<K_MESH::NO_Edge> frozenE;

#ifdef DEBUG_SWAPPER
  K_FLD::IntArray connectSpike, connectHat, connectSmall, connectUncomp;
#endif

  // 1st pass : detect bad triangles and store splitting points per edges :
  for (Si = 0; Si < nb_tris0; ++Si)
  {
    pS = connect.col(Si);

    N0 = *pS;
    N1 = *(pS + 1);
    N2 = *(pS + 2);

    E_Int ni;
    eDegenType type = degen_type2(coord, N0, N1, N2, tol2, lambdac, ni);

#ifdef DEBUG_SWAPPER
    if (type == OK) continue;
    else if (type == HAT)
      connectHat.pushBack(pS, pS + 3);
    else if (type == SPIKE)
      connectSpike.pushBack(pS, pS + 3);
    else if (type == SMALL)
      connectSmall.pushBack(pS, pS + 3);
#endif

    //assert(type == OK || type == HAT);

    if (type != HAT) continue;

    // reassign nodes to have N0 as hat
    N0 = *(pS + ni);
    N1 = *(pS + (ni + 1) % 3);
    N2 = *(pS + (ni + 2) % 3);

    Ei.setNodes(N1, N2);
    itE = E_to_Pts.find(Ei);

    if (itE != E_to_Pts.end()) continue; // current rule one point at a time

    //add N0 to NiNJ
    E_to_Pts[Ei].push_back(N0);
  }

#ifdef DEBUG_SWAPPER
  if (connectSpike.cols()) medith::write("spike_triangles.mesh", coord, connectSpike, "TRI");
  if (connectHat.cols()) medith::write("hat_triangles.mesh", coord, connectHat, "TRI");
  if (connectSmall.cols()) medith::write("small_triangles.mesh", coord, connectSmall, "TRI");
  if (connectUncomp.cols()) medith::write("uncomp_triangles.mesh", coord, connectUncomp, "TRI");
#endif

  if (E_to_Pts.empty())
    return 0;

  // sort the points on edges
  //for (itE = E_to_Pts.begin(); itE != E_to_Pts.end(); ++itE)
  //__tidy_edge(coord, itE->second);

  // 2nd pass : do the split
  E_Int nb_split, nb_tris;
  //E_Int n1;
  //size_t sz;
  K_FLD::IntArray new_connect;
  //E_Int itermax = 10, iter = 1;
  do
  {
    frozenE.clear();
    nb_split = 0;
    nb_tris = connect.cols();
    new_connect = connect;
    std::vector<bool> keep(connect.cols(), true);

#ifdef DEBUG_SWAPPER
    K_FLD::IntArray connectSplit, connectInit;
    K_FLD::IntArray::iterator pSi;
#endif

    for (Si = 0; Si < nb_tris; ++Si)
    {
      pS = connect.col(Si);

      //
      for (n = 0; n < 3; ++n)
      {
        E_Int Ni = *(pS + n);
        E_Int n1 = (n + 1) % 3;
        E_Int Nip1 = *(pS + n1);
        E_Int Nip2 = *(pS + (n + 2) % 3);

        Ei.setNodes(Ni, Nip1);

        itE = E_to_Pts.find(Ei);
        if (itE == E_to_Pts.end())
          continue;

        if (frozenE.find(Ei) != frozenE.end()) continue;

        frozenE.insert(K_MESH::NO_Edge(Nip1, Nip2));
        frozenE.insert(K_MESH::NO_Edge(Nip2, Ni));

        E_Int N = itE->second[0];

        if (Nip2 == N) {
          keep[Si] = false;
          continue; //discard the guilty
        }

        //split into 2 triangles
        new_connect(n1, Si) = N;
        new_connect.pushBack(pS, pS + 3);
        E_Int last = new_connect.cols() - 1;
        new_connect(n, last) = N;

#ifdef DEBUG_SWAPPER
        connectSplit.pushBack(new_connect.col(Si), new_connect.col(Si) + 3);
        connectSplit.pushBack(new_connect.col(last), new_connect.col(last) + 3);
        connectInit.pushBack(connect.col(Si), connect.col(Si) + 3);
#endif
        oids.push_back(oids[Si]);
      }
    }

#ifdef DEBUG_SWAPPER
    medith::write("split_triangles.mesh", coord, connectSplit, "TRI");
    medith::write("init_triangles.mesh", coord, connectInit, "TRI");
#endif

    keep.resize(new_connect.cols(), true);

    // now sync for degeneracies removed in between
    K_CONNECT::keep<bool> pred(keep);
    K_CONNECT::IdTool::compress(new_connect, pred);
    K_CONNECT::IdTool::compress(oids, pred);

    nb_split = new_connect.cols() - connect.cols();
    connect = new_connect;

  } while (nb_split);

  return (connect.cols() - nb_tris0);
}

void SwapperT3::edit_T3_caracs(const K_FLD::FloatArray& crd, E_Int* pN)
{
  E_Float hmin, Lmin, lambda_min;
  E_Int himin=-1, limin=-1;
  Lmin = hmin = lambda_min = NUGA::FLOAT_MAX;
  E_Float lambda, h;
  
  for (E_Int n=0; n < 3; ++n)
  {
    h = K_MESH::Edge::edgePointMinDistance<3>(crd.col(pN[n]), crd.col(pN[(n + 1) % 3]), crd.col(pN[(n + 2) % 3]), lambda);
    if (h < hmin)
    {
      hmin = h;
      lambda_min = lambda;
      himin = n;
    }
    
    E_Float d = ::sqrt(NUGA::sqrDistance(crd.col(pN[n]), crd.col(pN[(n + 1) % 3]), 3));
    
    if (d < Lmin)
    {
      limin = n;
      Lmin = d;
    }
    
  }
  
  std::cout << "worst h/lambda is : " << hmin << "/" << lambda_min << " reached at i :" << himin << std::endl;
  std::cout << "worst L is : " << Lmin << " reached at i :" << limin << std::endl;
  
}
