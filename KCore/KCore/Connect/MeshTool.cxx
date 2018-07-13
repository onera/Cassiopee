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

#ifndef __GENERATOR_K_CONNECT_MESHTOOL_CXX__
#define __GENERATOR_K_CONNECT_MESHTOOL_CXX__

#include "MeshTool.h"

namespace K_CONNECT
{
  
template <typename OutputIterator>
E_Int
MeshTool::getContainingElements
(const E_Float* point, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
 const int_vector_type& ancestors, const bool_vector_type& mask, const tree_type& tree, OutputIterator out, size_type & N0) const
{
  size_type         K(E_IDX_NONE), b, count(0), Ni, Nj, Kneigh(E_IDX_NONE), A0;
  int_vector_type   anc;
  E_Float           q;

  // Fast search
  N0 = tree.getClose(point);
  
  //assert (N0 != E_IDX_NONE);
  if (N0 == E_IDX_NONE)
    return -1;

  A0 = ancestors[N0];
  //assert(mask[A0] == true);
  if (!mask[A0])
    return -1;
  if (A0 != E_IDX_NONE)
     K = __getContainingElement(point, pos, connect, neighbors, mask, A0);

  // Deeper search : loop over all the closest node ancestors.
  if (K == E_IDX_NONE)
  {
    N0 = tree.getClosest(point);
    getAncestors(N0, connect, ancestors, neighbors, std::back_inserter(anc));
    size_type sz = (size_type)anc.size();
    for (size_type i = 0; (i < sz) && (K == E_IDX_NONE); ++i)
    {
      if (anc[i] != E_IDX_NONE)
        K = __getContainingElement(point, pos, connect, neighbors, mask, anc[i]);
    }
  }

  // Brute force search : loop over all elements.
  if (K == E_IDX_NONE)
    K = __getContainingElement(point, pos, connect, neighbors, mask);

  if (K == E_IDX_NONE) // Error.
    return -1;

  // Check if point is matching any K node or is falling on an edge.
  K_FLD::IntArray::const_iterator pK = connect.col(K);
  for (b = 0; (b < K_MESH::Triangle::NB_NODES) && (count < 2); ++b)
  {
    Ni = *(pK + (b+1) % K_MESH::Triangle::NB_NODES);
    Nj = *(pK + (b+2) % K_MESH::Triangle::NB_NODES);
    assert(Ni != Nj);
    q = K_MESH::Triangle::qualityG<2>(point, pos.col(Ni), pos.col(Nj));
    if (q</*_tolerance*/E_EPSILON)
    {
      ++count;
      Kneigh = neighbors (b, K);
    }
  }
 
  if (count == 2) // Matching a node.
  {
    // Find out what node it is.
    N0 = (K_FUNC::sqrDistance(point, pos.col(*(pK+1)), 2) < K_FUNC::sqrDistance(point, pos.col(*(pK+2)), 2)) ? *(pK+1) : *(pK+2);
    N0 = (K_FUNC::sqrDistance(point, pos.col(N0), 2) < K_FUNC::sqrDistance(point, pos.col(*pK), 2)) ? N0 : *pK;
    
    getAncestors(N0, connect, ancestors, neighbors, out);
    return 2;
  }

  *(out++) = K;

  if ((count == 1) && (Kneigh != E_IDX_NONE)) // On an inner edge.
     *(out++) = Kneigh;

  return count;
}

template <typename OutputIterator>
void
MeshTool::getAncestors
(size_type N, const K_FLD::IntArray& connectM, const int_vector_type& ancestors, const K_FLD::IntArray& neighbors, OutputIterator out) const
{
  // Fast returns.
  if (N == E_IDX_NONE)              return;
  if (ancestors[N] == E_IDX_NONE)   return;

  size_type Kseed(ancestors[N]), K(Kseed), n;

  do
  {
    *(out++) = K;
    n = K_MESH::Triangle::getLocalNodeId(connectM.col(K), N);
    K = neighbors((n+1)%element_type::NB_NODES, K);
  }
  while ((K != Kseed) && (K != E_IDX_NONE));
  //fixme : algo pas complet : voir si necessite de recuperer tous les elements. autre choix, empecher de chosir un point de la box...
}

/*
template <typename KPredicate, typename BPredicate>
void
MeshTool::getConnexSet
(size_type K, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
 int_set_type& set, int_pair_set_type& bset,
 const KPredicate& kpred, const BPredicate& bpred) const
{
  if (!set.insert(K).second) // already in
    return;

  K_FLD::IntArray::const_iterator pK = connect.col(K);
  size_type Ni, Nj, K1;
  std::pair<size_type, size_type> Ei, rEi;

  for (size_type b = 0; b < K_MESH::Triangle::NB_NODES; ++b)
  {
    rEi.second = Ei.first = Ni = *(pK + (b+1) % K_MESH::Triangle::NB_NODES);
    rEi.first = Ei.second = Nj = *(pK + (b+2) % K_MESH::Triangle::NB_NODES);

    K1 = neighbors(b, K);

    if ((K1 == E_IDX_NONE) || bpred(Ei) || bpred(rEi))
    {
      bset.insert(int_pair_type(K,b));
      continue;
    }

    size_type j = K_MESH::Triangle::getOppLocalNodeId(K, b, connect, neighbors);
    
    if (bset.find(int_pair_type(K1,j)) != bset.end())
      continue;
    
    if (kpred(K1))
      getConnexSet(K1, connect, neighbors, set, bset, kpred, bpred);
    else
      bset.insert(int_pair_type(K,b));
  }
}
*/

template <typename KPredicate, typename BPredicate>
void
MeshTool::getConnexSet1
(size_type S, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
 int_set_type& set, int_pair_set_type& bset,
 const KPredicate& kpred, const BPredicate& bpred) const
{
  if (!set.insert(S).second) //already in
    return;

  size_type Sn, Ni(E_IDX_NONE), Nj(E_IDX_NONE);
  K_MESH::NO_Edge Ei;

  for (size_type b = 0; b < K_MESH::Triangle::NB_NODES; ++b)
  {
    Sn = neighbors(b, S);

    if ((Sn == E_IDX_NONE) || (_inval.find(Sn) != _inval.end()))
      continue;

    Ni = connect((b+1) % K_MESH::Triangle::NB_NODES, S);
    Nj = connect((b+2) % K_MESH::Triangle::NB_NODES, S);
    Ei.setNodes(Ni,Nj);

    if (bpred(Ei))
      _inval.insert(Sn);
    else if (kpred(Sn))
      getConnexSet1(Sn, connect, neighbors, set, bset, kpred, bpred);
  }
}
/*
///
template <typename KPredicate, typename BPredicate>
void
MeshTool::getConnexSet2
(size_type S, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
 int_set_type& set, int_pair_set_type& bset,
 const KPredicate& kpred, const BPredicate& bpred) const
{
 _pool.clear();
 _pool.push_back(S);
 _inval.clear();

 size_type Sn, Ni, Nj;
 K_MESH::NO_Edge E;

  while (!_pool.empty())
  {
    S = _pool[0];

    for (size_type b = 0; b < K_MESH::Triangle::NB_NODES; ++b)
    {
        Sn = neighbors(b, S);
        
        if ((Sn == E_IDX_NONE) || (_inval.find(Sn) != _inval.end()))
          continue;

        Ni = connect((b+1) % K_MESH::Triangle::NB_NODES, S);
        Nj = connect((b+2) % K_MESH::Triangle::NB_NODES, S);
        E.setNodes(Ni, Nj);
        
        if (bpred(E))
          _inval.insert(Sn);
        else if ((set.find(Sn) == set.end()) && kpred(Sn))
          _pool.push_back(Sn);
    }

    set.insert(S);
    _pool.pop_front();
  }
}
*/

///
template <E_Int DIM>
void
MeshTool::computeMinMaxEdgeSqrLength
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Float& min_d, E_Float& max_d)
{
  E_Float d;
  E_Int                           COLS(connect.cols()), NB_NODES(connect.rows()), seg, c, n;//, Ni, Nj;
  K_FLD::IntArray::const_iterator pS;

  min_d = K_CONST::E_MAX_FLOAT;
  max_d = - K_CONST::E_MAX_FLOAT;

  // Fast returns
  if (NB_NODES < 2 || COLS == 0)
    return;

  seg = (NB_NODES > 2) ? NB_NODES : 1; 

  for (c = 0; c < COLS; ++c)
  {
    pS = connect.col(c);

    for (n = 0; n < seg; ++n)
    {
      d = K_FUNC::sqrDistance(pos.col(*(pS+n)), pos.col(*(pS+(n+1)%NB_NODES)), DIM);
      /*if (d < min_d)
      {
        Ni = *(pS+n);
        Nj = *(pS+(n+1)%NB_NODES);
      }*/
      min_d = (d < min_d) ? d : min_d;
      max_d = (d > max_d) ? d : max_d;
    }
  }
}

template <E_Int DIM>
void
MeshTool::computeEdgesSqrLengths
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& L)
{
  E_Int COLS(connect.cols()), NB_NODES(connect.rows()), seg, c, n;
  K_FLD::IntArray::const_iterator pS;
  E_Float minL, maxL, l;

  seg = (NB_NODES > 2) ? NB_NODES : 1; 

  L.clear();
  L.resize(2, connect.cols());

  for (c = 0; c < COLS; ++c)
  {
    minL  = K_CONST::E_MAX_FLOAT;
    maxL  = -K_CONST::E_MAX_FLOAT;
    
    pS = connect.col(c);
    for (n = 0; n < seg; ++n)
    {
      l = K_FUNC::sqrDistance(pos.col(*(pS+n)), pos.col(*(pS+(n+1)%NB_NODES)), DIM);
      minL = (l < minL) ? l : minL;
      maxL = (l > maxL) ? l : maxL;
    }
    
    L(0, c) = minL;
    L(1, c) = maxL;
  }
}

///
template <E_Int DIM>
void
MeshTool::computeIncidentEdgesSqrLengths
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& L)
{
  E_Int                           COLS(connect.cols()), NB_NODES(connect.rows());
  E_Int                           Ni, Nj, seg, c, n;
  K_FLD::IntArray::const_iterator pS;
  E_Float                         l, min(-K_CONST::E_MAX_FLOAT), max(K_CONST::E_MAX_FLOAT);

  seg = (NB_NODES > 2) ? NB_NODES : 1; 

  L.clear();
  L.resize(1, pos.cols(), &max);  // init mins to DBL_MAX.
  L.resize(2, pos.cols(), &min);  // init maxs to -DBL_MAX.

  for (c = 0; c < COLS; ++c)
  {
    pS = connect.col(c);
    for (n = 0; n < seg; ++n)
    {
      Ni = *(pS+n);
      Nj = *(pS+(n+1)%NB_NODES);
      l = K_FUNC::sqrDistance(pos.col(Ni), pos.col(Nj), DIM);
      // update Ni.
      min = L(0, Ni);
      max = L(1, Ni);
      L(0, Ni) = (l < min) ? l : min;
      L(1, Ni) = (l > max) ? l : max;
      // update Nj.
      min = L(0, Nj);
      max = L(1, Nj);
      L(0, Nj) = (l < min) ? l : min;
      L(1, Nj) = (l > max) ? l : max;
    }
  }
}


///
template <typename Element_t>
bool
MeshTool::detectDuplicated
(const K_FLD::IntArray& connect, std::vector<E_Int>& dupIds)
{
  dupIds.clear();

  typedef std::map<Element_t, E_Int> ELT_to_Col_map;
  ELT_to_Col_map elt_to_col;
  typename ELT_to_Col_map::iterator it;
  E_Int Si, COLS(connect.cols());
  Element_t Ei;
  K_FLD::IntArray::const_iterator pS;
  bool has_duplicates=false;

  dupIds.resize(COLS);

  for (Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);
    Ei.setNodes(pS);
    it = elt_to_col.find(Ei);
    if (it == elt_to_col.end())
    {
      dupIds[Si] = Si;
      elt_to_col[Ei] = Si;
    }
    else
    {
      dupIds[Si] = it->second;
      has_duplicates=true;
    }
  }
  return has_duplicates;
}

///
template <typename Element_t>
E_Int
MeshTool::removeDuplicated
(K_FLD::IntArray& connect, std::vector<E_Int>& dupIds)
{
  K_FLD::IntArray                 cOut;
  E_Int                           Si, COLS(connect.cols());
  K_FLD::IntArray::const_iterator pS;
  
  //
  bool has_duplis = detectDuplicated<Element_t>(connect, dupIds);

  if (!has_duplis)
    return 0;
  //
  for (Si = 0; Si < COLS; ++Si)
  {
    if (dupIds[Si] == Si)
    {
      pS = connect.col(Si);
      cOut.pushBack(pS, pS+Element_t::NB_NODES);
    }
  }
  E_Int nb_removed = connect.cols() - cOut.cols();
  connect = cOut;
  return nb_removed;  
}

///
template <typename EdgeType>
void
MeshTool::extractEdges
(K_FLD::IntArray& connect, std::set<EdgeType>& edges)
{
  //
  size_type                       i, k, COLS(connect.cols()), ROWS(connect.rows());
  K_FLD::IntArray::const_iterator pS(connect.begin());

  edges.clear();

  //
  for (i = 0; i < COLS; ++i, pS = pS+ROWS)
  {
    for (k = 0; k < ROWS; ++k)
      edges.insert( EdgeType(*(pS+k), *(pS+(k+1)%ROWS)) );
  }
}

}

#endif
