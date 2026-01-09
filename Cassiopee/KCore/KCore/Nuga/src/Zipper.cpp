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

#include "Nuga/include/Zipper.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/BARSplitter.h"
#include "Nuga/include/merge.h"
#ifdef WIN32
#ifdef E_TIME
#include "Nuga/include/chrono.h"
#endif
#endif
#ifdef WIN32
#ifdef E_DEBUG
#include "meshIO/meshIO.h"
#endif
#endif

///
Zipper::Zipper(bool reorient):_reorient(reorient)
{
}

///
Zipper::~Zipper(void)
{
}

///
void
Zipper::setMates
(const K_FLD::FloatArray& posT3, const K_FLD::IntArray& connectT3,
 const K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates)
{
#ifdef E_TIME
  NUGA::chrono c;
  c.start();
#endif
  nmates.resize(posT3.cols(), UNKNO);

  //std::cout << "zipping : get boundary" << std::endl;

  E_Int nbE2 = connectE2.cols();
  if (nbE2 == 0)
    return;

#ifdef E_TIME
  std::cout << c.elapsed() << std::endl;
  c.start();
  //std::cout << "zipping : set bottoms" << std::endl;
#endif

  // For each boundary node, get its bottom mates.
  K_FLD::IntArray bottomE2;
  __setBottoms(posT3, connectT3, connectE2, bottomE2);
  
#ifdef E_TIME
  std::cout << c.elapsed() << std::endl;
  c.start();
  //std::cout << "zipping : virtual edges" << std::endl;
#endif

  // For each edge compute its virtual edge.
  K_FLD::FloatArray vpos;
  K_FLD::IntArray vE2;
  std::vector<E_Int> realID;
  __computeVirtualEdges(posT3, connectE2, bottomE2, vpos, vE2, realID);

#ifdef E_TIME
  std::cout << c.elapsed() << std::endl;
  c.start();
  //std::cout << "zipping : nodes mates" << std::endl;
#endif

  // Set the node pairs (zippable) and free nodes (external contours).
  __setNodeMates(vpos, vE2, realID, nmates);
  

#ifdef WIN32
#ifdef E_DEBUG
  /*
  K_FLD::IntArray cc;
  for (size_t i = 0; i < nmates.size(); ++i)
  {
    if ((nmates[i] < 0) || (nmates[i] == IDX_NONE))
      continue;
    E_Int E[] = {i, nmates[i]};
    cc.pushBack(E, E+2);
  }
  meshIO::write("links.mesh", posT3, cc);
  */

  K_FLD::IntArray ccc;
  for (E_Int c = 0; c < connectE2.cols(); ++c)
  {
    if ((nmates[connectE2(0,c)] == FREE)||(nmates[connectE2(1,c)] == FREE))
      ccc.pushBack(connectE2.col(c), connectE2.col(c)+2);
  }
  meshIO::write("free.mesh", posT3, ccc);
  return;
#endif
#endif

}

void 
Zipper::__setBottoms
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectT3,
 const K_FLD::IntArray& connectE2, K_FLD::IntArray& bottomE2)
{
  K_FLD::IntArray connectAT3, neighbors;

  bottomE2.clear();
  E_Int none = IDX_NONE;
  bottomE2.resize(connectE2.rows(), connectE2.cols(), &none);


  NUGA::MeshTool::getAttachedElements(connectE2, connectT3, connectAT3);
  NUGA::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectAT3, neighbors);

  std::map<K_MESH::Edge, E_Int> edge_to_col;
  std::map<K_MESH::Edge, E_Int>::const_iterator itE;
  E_Int nbE2 = connectE2.cols();
  for (E_Int i = 0; i < nbE2; ++i)
    edge_to_col.insert(std::make_pair(K_MESH::Edge(connectE2(0, i), connectE2(1, i)), i));

  E_Int nbelts = connectAT3.cols(), S, Sn, bn;
  K_FLD::IntArray::const_iterator pS, pSn;
  K_MESH::Triangle t;
  K_MESH::Edge b;
  E_Int Ni, Nj, Nk;
  for (S = 0; S < nbelts; ++S)
  {
    pS = connectAT3.col(S);

    for (E_Int n = 0; n < 3; ++n)
    {
      Ni = *(pS+(n+1)%3);
      Nj = *(pS+(n+2)%3);

      itE = edge_to_col.find(K_MESH::Edge(Ni, Nj));

      if (itE == edge_to_col.end())
        continue;

      Nk = *(pS+n);

      if (NUGA::sqrDistance(pos.col(Ni), pos.col(Nk), 3) < NUGA::sqrDistance(pos.col(Nj), pos.col(Nk), 3))
      {
        bottomE2(0, itE->second) = Nk;
        Sn = neighbors((n+1)%3, S);
        if (Sn == IDX_NONE)
          continue;
        pSn = connectAT3.col(Sn);
        bn = K_MESH::Triangle::getOppLocalNodeId(S, (n+1)%3, connectAT3, neighbors);
        Nk = *(pSn + bn);
        bottomE2(1, itE->second) = Nk;
      }
      else
      {
        bottomE2(1, itE->second) = Nk;
        E_Int x = K_MESH::Triangle::getLocalNodeId(pS, Nj);
        Sn = neighbors(x, S);
        if (Sn == IDX_NONE)
          continue;
        pSn = connectAT3.col(Sn);
        bn = K_MESH::Triangle::getOppLocalNodeId(S, x, connectAT3, neighbors);
        Nk = *(pSn + bn);
        bottomE2(0, itE->second) = Nk;
      }
    }
  }
}

///
void
Zipper::__computeVirtualEdges
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2,
 const K_FLD::IntArray& bottomE2,
 K_FLD::FloatArray& vpos, K_FLD::IntArray& vE2,
 std::vector<E_Int>& realID)
{
  vpos.clear();
  vE2.clear();
  vE2.reserve(2, connectE2.cols());
  std::vector<E_Int> nodesE2;
  connectE2.uniqueVals(nodesE2);
  vpos.reserve(3, nodesE2.size());

  K_FLD::IntArray::const_iterator pE2, pBot;
  E_Int Bi, Bj, ve[2];
  E_Float Vi[3], Vj[3], Pi[3], Pj[3];

  E_Int Ni, Nj, nbE2(connectE2.cols());
  for (E_Int i = 0; i < nbE2; ++i)
  {
    pE2 = connectE2.col(i);
    pBot = bottomE2.col(i);
    Ni = *pE2;
    Nj = *(pE2+1);
    Bi = *pBot;
    if (Bi == IDX_NONE)
      continue;
    Bj = *(pBot+1);
    if (Bj == IDX_NONE)
      continue;

    //virtual points.
    NUGA::diff<3>(pos.col(Ni), pos.col(Bi), Vi);
    NUGA::diff<3>(pos.col(Nj), pos.col(Bj), Vj);

    for (E_Int k = 0; k < 3; ++k)
    {
      Vi[k] *= 0.5;//fixme adjust
      Vj[k] *= 0.5;

      Pi[k] = pos(k, Ni) + Vi[k];
      Pj[k] = pos(k, Nj) + Vj[k];
    }

    vpos.pushBack(Pi, Pi+3);
    vpos.pushBack(Pj, Pj+3);

    realID.push_back(Ni);
    realID.push_back(Nj);

    ve[0] = vpos.cols() -2;
    ve[1] = vpos.cols() -1;
    vE2.pushBack(ve, ve+2);
  }

  std::vector<E_Int> newIds;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> posAcc(vpos);
  ::merge(posAcc, EPSILON, newIds);
  K_FLD::IntArray::changeIndices(vE2, newIds);

#ifdef WIN32
#ifdef E_DEBUG
  //meshIO::write("virtual.mesh", vpos, vE2);
#endif
#endif
}


///
void
Zipper::__setNodeMates
(const K_FLD::FloatArray& vpos, const K_FLD::IntArray& vE2, const std::vector<E_Int>& realID,
 std::vector<E_Int>& nmates)
{
  if (vE2.cols() == 0)
    return;
  // Compute the radius for each node to detect mates.
  //std::vector<E_Float> R2;
  //__computeRadius(vpos, vE2, posT3, realID, R2);

  // Get the nodes neighbouring.
  NUGA::int_pair_vector_type node_to_nodes;
  if (BARSplitter::getNodesNeighBouring(vE2, node_to_nodes))
    return;// Error

  // Put the contour in a hset.
  NUGA::non_oriented_edge_set_type hE2;
  NUGA::non_oriented_edge_set_type::const_iterator itE;
  E_Int nbE2 = vE2.cols();
  for (E_Int i = 0; i < nbE2; ++i)
    hE2.insert(K_MESH::NO_Edge(vE2(0, i), vE2(1, i)));

  // Set the virtual edge mates (with 2 nodes attached to 2 node of another virtual edge)
  std::vector<E_Int> nodesvE2;;
  vE2.uniqueVals(nodesvE2);

  K_FLD::ArrayAccessor<K_FLD::FloatArray> pA(vpos);
  K_SEARCH::KdTree<> tree(pA, nodesvE2);

  E_Int ni0, nj0, Ni0, Nj0, ni1, nj1, Ni1, Nj1;
  K_MESH::NO_Edge E;
  K_FLD::IntArray::const_iterator pB;
  bool unknown = false, free;
  for (E_Int i = 0; i < vE2.cols(); ++i)
  {
    pB = vE2.col(i);
    ni0 = *pB;
    nj0 = *(pB+1);
    Ni0 = realID[ni0];
    Nj0 = realID[nj0];
    unknown = false;

    ni1 = tree.getClosest(ni0);
    assert(ni1 != IDX_NONE);
    Ni1 = realID[ni1];
    free = (ni1 == node_to_nodes[ni0].first);
    free |= (ni1 == node_to_nodes[ni0].second);
    if (free)
    {
      nmates[Ni0] = FREE;
      unknown = true;
    }
    else if (nmates[Ni1] == FREE) // Reset a "non true" free node
      nmates[Ni1] = UNKNO;

    nj1 = tree.getClosest(nj0);
    assert(nj1 != IDX_NONE);
    Nj1 = realID[nj1];
    free = (nj1 == node_to_nodes[nj0].first);
    free |= (nj1 == node_to_nodes[nj0].second);
    if (free)
    {
      nmates[Nj0] = FREE;
      unknown = true;
    }
    else if (nmates[Nj1] == FREE) // Reset a "non true" free node
      nmates[Nj1] = UNKNO;

    if (unknown)
      continue;

    if ((ni0 != tree.getClosest(ni1)) || (nj0 != tree.getClosest(nj1)))
      continue;

    E.setNodes(ni1, nj1);
    if (hE2.find(E) == hE2.end())
      continue;

    nmates[Ni0] = Ni1;
    nmates[Ni1] = Ni0;
    nmates[Nj0] = Nj1;
    nmates[Nj1] = Nj0;
  }
}

///
void
Zipper::zip
(const K_FLD::IntArray& connectE2,
 const std::vector<E_Int>& nmates, K_FLD::IntArray& zipped)
{
  zipped.clear();
  
  E_Int nbE2 = connectE2.cols();
  if (nbE2 == 0)
    return;

  if (nmates.empty())
    return;

  NUGA::oriented_edge_set_type hE2;
  for (E_Int i = 0; i < nbE2; ++i)
    hE2.insert(K_MESH::Edge(connectE2(0, i), connectE2(1, i)));

  //
  E_Int Ni, Nj, Nio, Njo, T1[3], T2[3];
  for (E_Int i = 0; i < nbE2; ++i)
  {
     Ni = connectE2(0,i);
     Nj = connectE2(1,i);
     Nio = nmates[Ni];
     Njo = nmates[Nj];

     if ((Nio < 0)|| (Nio == IDX_NONE))
       continue;
     if ((Njo < 0)|| (Njo == IDX_NONE))
       continue;

     if ((Nio == Nj) || (Njo == Ni))
       continue;

     if (hE2.find(K_MESH::Edge(Njo, Nio)) == hE2.end())
       continue;

     //assert(nmates[Nio] == Ni);
     //assert(nmates[Njo] == Nj);

     // zip
     T1[0] = Nj;
     T1[1] = Ni;
     T1[2] = Nio;
     T2[0] = Nio;
     T2[1] = Njo;
     T2[2] = Nj;

     zipped.pushBack(T1, T1+3);
     zipped.pushBack(T2, T2+3);

     hE2.erase(K_MESH::Edge(Ni, Nj));
     hE2.erase(K_MESH::Edge(Njo, Nio));
  }

/*
///
void
Zipper::__computeRadius
(const K_FLD::FloatArray& vpos, const K_FLD::IntArray& vE2, const K_FLD::FloatArray& pos,
 const std::vector<E_Int>& realID, std::vector<E_Float>& R2)
{
  R2.clear();
  R2.resize(vpos.cols(), NUGA::FLOAT_MAX);

  K_FLD::FloatArray L2;
  NUGA::MeshTool::computeEdgesSqrLengths<3>(vpos, vE2, L2);

  E_Float d2;
  E_Int   n, Ni, Nj, NBPOINTS(vpos.cols()), NBEDGE(vE2.cols());
  K_FLD::IntArray::const_iterator pB;

  for (n = 0; n < NBPOINTS; ++n)
  {
    d2 = NUGA::sqrDistance(vpos.col(n), pos.col(realID[n]), 3);
    R2[n] = d2;
  }

  for (n = 0; n < NBEDGE; ++n)
  {
    d2 = L2(0,n);
    pB = vE2.col(n);
    Ni = *pB;
    Nj = *(pB+1);
    R2[Ni] = std::min(R2[Ni], d2);
    R2[Nj] = std::min(R2[Nj], d2);
  }
}
*/

}

