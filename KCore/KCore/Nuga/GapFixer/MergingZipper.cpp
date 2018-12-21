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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "MergingZipper.h"
#include "Connect/MeshTool.h"
#include "Connect/merge.h"
#ifdef WIN32
#ifdef E_TIME
#include "../Delaunay/chrono.h"
#endif
#endif
#ifdef WIN32
#ifdef E_DEBUG
#include "meshIO/meshIO.h"
#endif
#endif

///
MergingZipper::MergingZipper(bool reorient):Zipper(reorient)
{
}

///
MergingZipper::~MergingZipper(void)
{
}

///
void
MergingZipper::setMates
(const K_FLD::FloatArray& posT3, const K_FLD::IntArray& connectT3,
 const K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates)
{
#ifdef E_TIME
  DELAUNAY::chrono c;
  c.start();
#endif

  if (!_reorient) // Merging mates are only useful if we need to reorient.
    return;

  nmates.resize(posT3.cols(), UNKNO);

  E_Int nbE2 = connectE2.cols();
  if (nbE2 == 0)
    return;

  std::vector<E_Int> realID(nmates.size(), E_IDX_NONE);
  for (size_t i = 0; i < realID.size(); ++i)
    realID[i] = i;

#ifdef E_TIME
  std::cout << c.elapsed() << std::endl;
  c.start();
  std::cout << "zipping : nodes mates" << std::endl;
#endif

  // Set the node pairs (mergeable) and free nodes (external contours).
  __setNodeMates(posT3, connectE2, realID, nmates);
}


void
MergingZipper::merge
(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
 K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates)
{
  std::vector<E_Int>  nodes, new_IDs;
  //E_Int nb_bound;
  E_Int Mi;

  connectE2.uniqueVals(nodes); // Component's boundary nodes.

  //nb_bound = nodes.size();

  // Merge the remaining boundary nodes based on the min boundary edge length
  // (e.g. several patches coincident at the same node).

  E_Float tol = __computeTolerance(pos, connectE2, nodes);

  K_FLD::ArrayAccessor<K_FLD::FloatArray> posAcc(pos);
  std::sort(nodes.begin(), nodes.end());

  E_Int nb_merges = ::merge(posAcc, tol, nodes, nodes, new_IDs); // Do the merge.

  if (nb_merges)
    K_FLD::IntArray::changeIndices(connect, new_IDs);

  //E_Int duplis = 
  __removeDuplicates(connect); //fixme : necessary ?

  if (_reorient)
  {
    //unmate
    for (size_t i = 0; i < new_IDs.size(); ++i)
    {
      if ((new_IDs[i] != (E_Int)i) && (new_IDs[i] != E_IDX_NONE)) // i.e a merge occured
      {
        Mi = nmates[i];
        nmates[i] = E_IDX_NONE; // unmate
        if ((Mi != E_IDX_NONE) && (Mi > -1))
          nmates[Mi] = E_IDX_NONE; // unmate
        if (new_IDs[i] > -1)
        {
          Mi = nmates[new_IDs[i]];
          nmates[new_IDs[i]] = E_IDX_NONE; // unmate
          if ((Mi != E_IDX_NONE) && (Mi > -1))
            nmates[Mi] = E_IDX_NONE; // unmate
        }
      }
    }
  }
}

///
E_Float
MergingZipper::__computeTolerance
(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
 const std::vector<E_Int>& nodes)
{
  E_Float           tol = K_CONST::E_MAX_FLOAT, max_d;
  K_FLD::FloatArray L2;
  K_FLD::IntArray   connectA;

  K_CONNECT::MeshTool::computeMinMaxEdgeSqrLength<3>(pos, connect, tol, max_d);
  return std::max((::sqrt(tol) - E_EPSILON), E_EPSILON);
}

E_Int
MergingZipper::__removeDuplicates(K_FLD::IntArray& connect)
{
  std::vector<K_CONT_DEF::non_oriented_int_pair_set_type> node_to_edges;
  std::vector<E_Int>              nodes;
  K_FLD::IntArray                 clean_connect;
  K_FLD::IntArray::const_iterator pS;
  E_Int                           Si, n, NB_COLS(connect.cols()), n0/*min id*/, N0, n1, n2, ret;

  n0 = E_IDX_NONE;
  connect.uniqueVals(nodes);
  node_to_edges.resize(*std::max_element(nodes.begin(), nodes.end())+1);

  for (Si = 0; Si < NB_COLS; ++Si)
  {
    N0 = E_IDX_NONE;
    pS = connect.col(Si);

    for (n = 0; (n < 3); ++n)
    {
      if (*(pS+n) < N0)
      {
        N0 = *(pS+n);
        n0 = n;
      }
    }
    n1 = (n0+1)%3;
    n2 = (n1+1)%3;

    if (node_to_edges[N0].insert(K_MESH::NO_Edge(*(pS+n1), *(pS+n2))).second)
      clean_connect.pushBack(pS, pS+3);
  }

  ret = connect.cols() - clean_connect.cols();
  if (ret) connect = clean_connect;
  return ret;
}
