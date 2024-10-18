/*    
    Copyright 2013-2024 Onera.

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

#ifndef _NUGA_MESHTOOL_CXX__
#define _NUGA_MESHTOOL_CXX__

#include "Nuga/include/MeshTool.h"
#include "Nuga/include/DynArray.h"

namespace NUGA
{
  
template <typename OutputIterator>
E_Int
MeshTool::getContainingElements
(const E_Float* point, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
 const int_vector_type& ancestors, const bool_vector_type& mask, const tree_type& tree, OutputIterator out, size_type & N0) const
{
  size_type         K(IDX_NONE), b, count(0), Ni, Nj, Kneigh(IDX_NONE), A0;
  int_vector_type   anc;
  E_Float           q;

  // Fast search
  N0 = tree.getClose(point);
  
  //assert (N0 != IDX_NONE);
  if (N0 == IDX_NONE)
    return -1;

  A0 = ancestors[N0];
  //assert(mask[A0] == true);
  if (!mask[A0])
    return -1;
  if (A0 != IDX_NONE)
     K = __getContainingElement(point, pos, connect, neighbors, mask, A0);

  // Deeper search : loop over all the closest node ancestors.
  if (K == IDX_NONE)
  {
    N0 = tree.getClosest(point);
    getAncestors(N0, connect, ancestors, neighbors, std::back_inserter(anc));
    size_type sz = (size_type)anc.size();
    for (size_type i = 0; (i < sz) && (K == IDX_NONE); ++i)
    {
      if (anc[i] != IDX_NONE)
        K = __getContainingElement(point, pos, connect, neighbors, mask, anc[i]);
    }
  }

  // Brute force search : loop over all elements.
  if (K == IDX_NONE)
    K = __getContainingElement(point, pos, connect, neighbors, mask);

  if (K == IDX_NONE) // Error.
    return -1;

  // Check if point is matching any K node or is falling on an edge.
  K_FLD::IntArray::const_iterator pK = connect.col(K);
  for (b = 0; (b < K_MESH::Triangle::NB_NODES) && (count < 2); ++b)
  {
    Ni = *(pK + (b+1) % K_MESH::Triangle::NB_NODES);
    Nj = *(pK + (b+2) % K_MESH::Triangle::NB_NODES);
    assert(Ni != Nj);
    q = K_MESH::Triangle::qualityG<2>(point, pos.col(Ni), pos.col(Nj));
    if (q</*_tolerance*/EPSILON)
    {
      ++count;
      Kneigh = neighbors (b, K);
    }
  }
 
  if (count == 2) // Matching a node.
  {
    // Find out what node it is.
    N0 = (NUGA::sqrDistance(point, pos.col(*(pK+1)), 2) < NUGA::sqrDistance(point, pos.col(*(pK+2)), 2)) ? *(pK+1) : *(pK+2);
    N0 = (NUGA::sqrDistance(point, pos.col(N0), 2) < NUGA::sqrDistance(point, pos.col(*pK), 2)) ? N0 : *pK;
    
    getAncestors(N0, connect, ancestors, neighbors, out);
    return 2;
  }

  *(out++) = K;

  if ((count == 1) && (Kneigh != IDX_NONE)) // On an inner edge.
     *(out++) = Kneigh;

  return count;
}

template <typename OutputIterator>
void
MeshTool::getAncestors
(size_type N, const K_FLD::IntArray& connectM, const int_vector_type& ancestors, const K_FLD::IntArray& neighbors, OutputIterator out) const
{
  // Fast returns.
  if (N == IDX_NONE)              return;
  if (ancestors[N] == IDX_NONE)   return;

  size_type Kseed(ancestors[N]), K(Kseed), n;

  do
  {
    *(out++) = K;
    n = K_MESH::Triangle::getLocalNodeId(connectM.col(K), N);
    K = neighbors((n+1)%element_type::NB_NODES, K);
  }
  while ((K != Kseed) && (K != IDX_NONE));
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

    if ((K1 == IDX_NONE) || bpred(Ei) || bpred(rEi))
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

  size_type Sn, Ni(IDX_NONE), Nj(IDX_NONE);
  K_MESH::NO_Edge Ei;

  for (size_type b = 0; b < K_MESH::Triangle::NB_NODES; ++b)
  {
    Sn = neighbors(b, S);

    if ((Sn == IDX_NONE) || (_inval.find(Sn) != _inval.end()))
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
        
        if ((Sn == IDX_NONE) || (_inval.find(Sn) != _inval.end()))
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

  min_d = NUGA::FLOAT_MAX;
  max_d = - NUGA::FLOAT_MAX;

  // Fast returns
  if (NB_NODES < 2 || COLS == 0)
    return;

  seg = (NB_NODES > 2) ? NB_NODES : 1; 

  for (c = 0; c < COLS; ++c)
  {
    pS = connect.col(c);

    for (n = 0; n < seg; ++n)
    {
      d = NUGA::sqrDistance(pos.col(*(pS+n)), pos.col(*(pS+(n+1)%NB_NODES)), DIM);
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
    minL  = NUGA::FLOAT_MAX;
    maxL  = -NUGA::FLOAT_MAX;
    
    pS = connect.col(c);
    for (n = 0; n < seg; ++n)
    {
      l = NUGA::sqrDistance(pos.col(*(pS+n)), pos.col(*(pS+(n+1)%NB_NODES)), DIM);
      minL = (l < minL) ? l : minL;
      maxL = (l > maxL) ? l : maxL;
    }
    
    L(0, c) = minL;
    L(1, c) = maxL;
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

template <typename Triangulator>
void MeshTool::refine_T3s(const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, 
                          std::map<K_MESH::NO_Edge, Vector_t<E_Int>>& edge_to_refined_edge, std::vector<E_Int>& oids)
{
  oids.clear();
  K_FLD::IntArray new_cnt;
  
  //first sort the points
  std::vector<std::pair<E_Float, E_Int> > sorted_nodes;
  std::vector<E_Int> usorted_nodes;
  for (auto& i : edge_to_refined_edge)
  {
    E_Int N1 = i.first.node(0);
    E_Int N2 = i.first.node(1);
    
    sorted_nodes.clear();
    sorted_nodes.push_back(std::make_pair(-NUGA::FLOAT_MAX,N1));
    sorted_nodes.push_back(std::make_pair(NUGA::FLOAT_MAX, N2));

    E_Float dN1N22 = NUGA::sqrDistance(coord.col(N1), coord.col(N2), 3);

    for (auto& N : i.second)
    {
      if (N == N1 || N == N2) continue;
      
      E_Float dN1n2 = NUGA::sqrDistance(coord.col(N1), coord.col(N), 3);
      sorted_nodes.push_back(std::make_pair(dN1n2/dN1N22, N));
    }

    std::sort(sorted_nodes.begin(), sorted_nodes.end());
    usorted_nodes.clear();
    
    E_Float prev_d = 0.;
    for (size_t n = 0; n < sorted_nodes.size(); ++n)
    {
      E_Float d = sorted_nodes[n].first;
      if (d != prev_d)
      {
        usorted_nodes.push_back(sorted_nodes[n].second);
        prev_d = d;
      }
    }
    
    edge_to_refined_edge[i.first] = usorted_nodes; // now sorted and unic vals
  }
  
  K_MESH::NO_Edge E;
  std::vector<E_Int> pg_molec;
  Triangulator dt;

  for (E_Int i = 0; i < connectT3.cols(); ++i)
  {
    pg_molec.clear();

    for (E_Int n = 0; n < 3; ++n)
    {
      E_Int Nn = connectT3(n, i);
      E_Int Np1 = connectT3((n+1)%3, i);

      E.setNodes(Nn, Np1);

      auto it = edge_to_refined_edge.find(E);
      if (it != edge_to_refined_edge.end())
      {
        auto nodes = it->second;
        if (Nn == E.node(0)) //same orientation
          for (size_t j = 0; j< nodes.size() - 1; ++j)
            pg_molec.push_back(nodes[j]);
        else
          for (E_Int j = nodes.size() - 1; j>0; --j)
            pg_molec.push_back(nodes[j]);
      }
      else
        pg_molec.push_back(Nn);
    }

    K_MESH::Polygon::triangulate(dt, coord, &pg_molec[0], pg_molec.size(), 0, new_cnt, false, true);
    oids.resize(new_cnt.cols(), i);
  }
  
  connectT3 = new_cnt;
}

void MeshTool::append_all_topo_paths
(const K_FLD::FloatArray& crd, E_Int Nstart, E_Int Nend, const id_to_ids_t& nodes_graph, idlists_t& paths)
{
  
  paths.clear();
  
  std::deque<E_Int> path;
  
  id_to_ids_t::const_iterator itb = nodes_graph.find(Nstart);
  std::cout << "before assert..." << std::endl;
  assert(itb != nodes_graph.end());

  std::cout << "ONE LINK MISSING..." << std::endl;
  // ONE LINK MISSING
  id_to_ids_t::const_iterator ite = nodes_graph.find(Nend);
  std::cout << "before assert..." << std::endl;
  assert(ite != nodes_graph.end());
  {
    std::cout << "before set_intersection..." << std::endl;
    std::vector<E_Int> commonN;
    std::set_intersection(itb->second.begin(), itb->second.end(), ite->second.begin(), ite->second.end(), std::back_inserter(commonN));
    
#ifdef DEBUG_SPLITTER
    for (E_Int i=0; i < commonN.size(); ++i) std::cout << commonN[i] << " ";
     std::cout << std::endl;
#endif
     
    for (auto&n : commonN)
    {
      if (n == Nstart || n == Nend) continue;
      
      path.clear();
      
      path.push_back(Nstart);
      path.push_back(n);
      path.push_back(Nend);
      
      paths.push_back(path);
    }
  }
  std::cout << "TWO LINK MISSING..." << std::endl;
  // TWO LINKS MISSING
  for (auto nb = itb->second.begin(); nb != itb->second.end(); ++nb)
  {
    E_Int Ni = *nb;
    
    if (Ni == Nend)
      continue;
    
    for (auto ne = ite->second.begin(); ne != ite->second.end(); ++ne)
    {
      E_Int Nj = *ne;
      
      if (Nj == Nstart)
        continue;

      auto itNj = nodes_graph.find(Nj);
      assert(itNj != nodes_graph.end());

      for (auto ne2 = itNj->second.begin(); ne2 != itNj->second.end(); ++ne2)
      {
        E_Int Nk = *ne2;
        if (Nk != Ni)
          continue;
        
        path.clear();
      
        path.push_back(Nstart);
        path.push_back(Ni);
        path.push_back(Nj);
        path.push_back(Nend);
      
        paths.push_back(path);
      }
    }
  }
  
  //if (!paths.empty()) return;
  std::cout << "THREE LINK MISSING..." << std::endl;
  // 3 LINKS MISSING
  for (auto nb = itb->second.begin(); nb != itb->second.end(); ++nb)
  {
    E_Int Ni = *nb;
    if (Ni == Nend)
      continue;
    
    auto itb1 = nodes_graph.find(Ni);
    
    for (auto ne = ite->second.begin(); ne != ite->second.end(); ++ne)
    {
      E_Int Nj = *ne;
     
      if (Nj == Nstart)
        continue;
      
      auto ite1 = nodes_graph.find(Nj);

      std::vector<E_Int> commonN;
      std::set_intersection(itb1->second.begin(), itb1->second.end(), ite1->second.begin(), ite1->second.end(), std::back_inserter(commonN));

      for (auto& n :commonN)
      {
        if (n == Nstart || n == Nend) continue;
        
        path.clear();
      
        path.push_back(Nstart);
        path.push_back(Ni);
        path.push_back(n);
        path.push_back(Nj);
        path.push_back(Nend);
      
        paths.push_back(path);
      }
    }
  }
  
}

template <typename IntCONT>
inline void MeshTool::get_farthest_point_to_edge 
(const K_FLD::FloatArray& crd, E_Int Ni, E_Int Nj, const IntCONT& list, E_Int& Nf, E_Float& d2)
{
  d2 = 0.;
  Nf = IDX_NONE;
  
  //std::cout << "Edge : " << Ni << "/" << Nj << std::endl;
  //std::cout << "dim : " << crd.rows() << std::endl;
  
  E_Float d22 = 0., lambda;
  for (size_t k = 0; k < list.size(); ++k)
  {
    //std::cout << "checked point : " << list[k] << std::endl;
    d22 = K_MESH::Edge::edgePointMinDistance2<2>(crd.col(Ni), crd.col(Nj), crd.col(list[k]), lambda);
    
    //std::cout << "distance found : " << d22 << std::endl;
    //std::cout << "lambda ? : " << lambda << std::endl;
    
    if (lambda < -EPSILON || lambda > 1 + EPSILON)
    {// do not consider this path
      d2 = NUGA::FLOAT_MAX;
      Nf = IDX_NONE;
      break;
    }
    
    if (d22 > d2)
    {
      Nf = list[k];
      d2 = d22;
    }
  }
  
  //std::cout << "farthest point :" << Nf << std::endl;
  //std::cout << "farthest dist :" << d2 << std::endl;
}

template <typename IntCont, short DIM >
inline void
NUGA::MeshTool::reorder_nodes_on_edge
(const K_FLD::FloatArray& pos, IntCont& nodes, int idx_start, std::vector<std::pair<E_Float, E_Int> >& sorter)
{
  E_Float L;
  E_Int   E[2];

  size_t nb_nodes = nodes.size();

  const E_Float *P0 = pos.col(nodes[0] - idx_start);
  L = NUGA::sqrDistance(pos.col(nodes[1] - idx_start), P0, DIM);

  sorter.clear();
  sorter.push_back(std::make_pair(0., nodes[0]));

  for (size_t k = 1; k < nb_nodes; ++k)
  {
    L = NUGA::sqrDistance(pos.col(nodes[k] - idx_start), P0, DIM);
    sorter.push_back(std::make_pair(L, nodes[k]));
  }

  std::sort(sorter.begin(), sorter.end());

  nodes.clear();
  E[0] = sorter[0].second;
  nodes.push_back(E[0]);

  for (size_t s = 0; s < sorter.size() - 1; ++s)
  {
    E[1] = sorter[s + 1].second;

    if (E[1] != E[0])
    {
      nodes.push_back(E[1]);
      E[0] = E[1];
    }
  }
}
  
}

#endif
