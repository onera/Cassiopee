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

#include "Nuga/include/BARSplitter.h"
#include "Nuga/include/ContourSplitter.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/KdTree.h"
#ifdef DEBUG_BARSPLITTER
#include <iostream>
#endif
//=============================================================================
// Splits the connectivity (only) with the edge N0N1 and stores the bits into 
// separate containers.
//=============================================================================
void BARSplitter::split(const K_FLD::FloatArray& pos, E_Int dim, 
                        const K_FLD::IntArray& BAR, 
                        E_Int N0, E_Int N1,
                        std::vector<K_FLD::IntArray>& connectBout, 
                        E_Bool enclose, E_Float tolerance)
{
  connectBout.clear();

  // Fast returns
  // N0 or N1 are not in the mesh
  // or N0N1 is an existing edge
  if (__fastReturn(BAR, N0, N1)) return;

  __split(pos, dim, BAR, N0, N1, connectBout, enclose, tolerance);
}

//=============================================================================
/* Split the connectivity with the edge N0N1 and stores the connectivity bits 
   and corresponding coordinates into separate containers. */
//=============================================================================
void BARSplitter::split(
  const K_FLD::FloatArray& pos, E_Int dim, 
  const K_FLD::IntArray& connectBin, E_Int N0, E_Int N1,
  std::vector<K_FLD::FloatArray> & posOut, 
  std::vector<K_FLD::IntArray> & connectBout,
  E_Bool enclose, E_Float tolerance)
{
  split(pos, dim, connectBin, N0, N1, connectBout, enclose, tolerance);

  ContourSplitter<K_MESH::Edge, E_Int>::splitCoordinates(connectBout, pos, posOut);
}

//=============================================================================
void 
BARSplitter::split
(const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& connectBin, const K_FLD::IntArray& cuttingEdges,
 std::vector<K_FLD::IntArray> & connectBout, E_Float tolerance)
{
  E_Int nb_edges = cuttingEdges.cols(), N0, N1;
  std::vector<K_FLD::IntArray> cOut;

  // Clear output.
  connectBout.clear();
  connectBout.push_back(connectBin);

  for (E_Int e = 0; e < nb_edges; ++e)
  {
    N0 = cuttingEdges(0, e);
    N1 = cuttingEdges(1, e);

    for (size_t c = 0; c < connectBout.size(); ++c)
    {
      split(pos, dim, connectBout[c], N0, N1, cOut, 1/*enclose*/, tolerance);
      if (cOut.size() == 2)
      {
        connectBout[c] = cOut[0];
        connectBout.push_back(cOut[1]);
        break;
      }
    }
  }
}

//=============================================================================
void 
BARSplitter::split
(const K_FLD::FloatArray& pos, E_Int dim, const std::vector<K_FLD::IntArray>& connectIns, const K_FLD::IntArray& cuttingEdges,
 std::vector<K_FLD::IntArray> & connectBout, E_Float tolerance)
{
  std::vector<K_FLD::IntArray> cOut;
  for (size_t i = 0; i < connectIns.size(); ++i)
  {
    split(pos, dim, connectIns[i], cuttingEdges, cOut, tolerance);
    connectBout.insert(connectBout.end(), cOut.begin(), cOut.end());
  }
}

//=============================================================================
void BARSplitter::__split
(const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& BAR, E_Int N0, E_Int N1,
 std::vector<K_FLD::IntArray> & connectBout, E_Bool enclose, E_Float tolerance)
{
  connectBout.clear();

  int_set_type cuttingNodes;
  cuttingNodes.insert(N0);
  cuttingNodes.insert(N1);

  ContourSplitter<K_MESH::Edge, E_Int>::splitConnectivity(BAR, cuttingNodes, connectBout);

  // Close the contour if required and possible (Ensure that N0N1 doesn't intersect any other edge).
  enclose = enclose && __canBeClosed(pos, dim, BAR, N0, N1, tolerance);
  if (enclose)// Get orientation and close.
  {
    E_Int N0N1[] = {N0, N1};
    E_Int N1N0[] = {N1, N0};
    E_Int orient;
    for (size_t i = 0; i < connectBout.size(); ++i)
    {
      orient = __getOrientation(connectBout[i], N0);
      if (orient == 1)
        connectBout[i].pushBack(N0N1, N0N1+2);
      else if (orient == 0)
        connectBout[i].pushBack(N1N0, N1N0+2);
    }
  }
}

//=============================================================================
E_Bool
BARSplitter::__fastReturn(const K_FLD::IntArray& connectB, E_Int N0, E_Int N1)
{
  if (!__isBARNode(connectB, N0)) return true;
  if (!__isBARNode(connectB, N1)) return true;
  if (__isBAREdge(connectB, N0, N1)) return true;
  return false;
}

//=============================================================================
E_Bool
BARSplitter::__isBARNode(const K_FLD::IntArray& connectB, E_Int N)
{
  K_FLD::IntArray::const_iterator pS;
  E_Int NB_ELTS = connectB.cols();
  E_Bool is_in = false;

  for (E_Int Si = 0; (Si < NB_ELTS) && (is_in == false); ++Si)
  {
    pS = connectB.col(Si);
    is_in = ((*pS == N) || (*(pS+1) == N));
  }

  return is_in;
}

//=============================================================================
E_Bool
BARSplitter::__isBAREdge(const K_FLD::IntArray& connectB, E_Int N0, E_Int N1)
{
  K_FLD::IntArray::const_iterator pS;
  E_Bool is_in = false;

  for (E_Int Si = 0; (Si < connectB.cols()) && (is_in == false); ++Si)
  {
    pS = connectB.col(Si);
    is_in = ((*pS == N0) && (*(pS+1) == N1)) || ((*pS == N1) && (*(pS+1) == N0));
  }

  return is_in;
}

//=============================================================================
E_Bool
BARSplitter::__canBeClosed
(const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& connectBin, E_Int N0, E_Int N1, E_Float tolerance)
{
  bool                              doable(true);
  E_Bool                            overlap;
  K_FLD::FloatArray::const_iterator P0 = pos.col(N0);
  K_FLD::FloatArray::const_iterator P1 = pos.col(N1);
  E_Int                             Ni, Nj, NBOUND(connectBin.cols());
  E_Float                           u00, u01, u10, u11;

  typedef bool (*pInterFunc)(const E_Float* P0, const E_Float* P1, const E_Float* Q0, const E_Float* Q1,
                             E_Float tolerance, E_Bool tol_is_abolute, 
                             E_Float& u00, E_Float& u01, E_Float& u10,
                             E_Float& u11, E_Bool& overlap);

  pInterFunc intersect;

  if (dim == 2)
    intersect = &K_MESH::Edge::intersect<2>;
  else
    intersect = &K_MESH::Edge::intersect<3>;

  for (E_Int Si = 0; (Si < NBOUND) && doable; ++Si)
  {
    Ni = connectBin(0,Si);
    Nj = connectBin(1,Si);

    if ((Ni == N0) || (Ni == N1) || (Nj == N0) || (Nj == N1))
      continue;

    doable = !intersect(P0, P1, pos.col(Ni), pos.col(Nj), tolerance, true/*absolute*/, u00, u01, u10, u11, overlap);
  }

  return doable;
}

//=============================================================================
E_Int
BARSplitter::__getOrientation(const K_FLD::IntArray& connectB, E_Int N0)
{
  K_FLD::IntArray::const_iterator pS;
  E_Int orient = -1;

  for (E_Int Si = 0; (Si < connectB.cols()) && (orient == -1); ++Si)
  {
    pS = connectB.col(Si);
    for (E_Int n = 0; (n < 2) && (orient == -1); ++n)
      orient = (*(pS+n) == N0) ? n : -1;
  }

  return orient;
}

//=============================================================================
E_Int
BARSplitter::get_node_to_nodes(const K_FLD::IntArray& connectE2, std::map< E_Int, std::vector<E_Int> >& node_to_nodes)
{
  K_FLD::IntArray::const_iterator pS;
  E_Int                           NBELTS(connectE2.cols()), e;
  
  node_to_nodes.clear();

  for (e = 0; e < NBELTS; ++e)
  {
    pS = connectE2.col(e);
    
    node_to_nodes[*pS].push_back(*(pS+1));
    node_to_nodes[*(pS+1)].push_back(*pS);
  }
  return 0;
}

//=============================================================================
// Manifold version of the above function : assume manifoldness to give a right result
E_Int
BARSplitter::get_node_to_nodes
(const K_FLD::IntArray& connectE2, std::map< E_Int, std::pair<E_Int, E_Int> >& node_to_nodes)
{
  // WARNING : do not detect non-manifoldness. We get 2 connections (max) per node.

  K_FLD::IntArray::const_iterator       pS;
  E_Int                                 NBELTS(connectE2.cols()), e;
  
  std::map< E_Int, std::pair<E_Int, E_Int> >::iterator it;
  
  node_to_nodes.clear();

  for (e = 0; e < NBELTS; ++e)
  {
    pS = connectE2.col(e);
    
    it = node_to_nodes.find(*pS);
    if (it == node_to_nodes.end())
    {
      node_to_nodes[*pS]; 
      it = node_to_nodes.find(*pS);
      it->second.first=it->second.second=IDX_NONE;
    }

    if (it->second.second == IDX_NONE)
    {it->second.second=*(pS+1);}
    else
      it->second.first=*(pS+1);
    
    it = node_to_nodes.find(*(pS+1));
    if (it == node_to_nodes.end())
    {
      node_to_nodes[*(pS+1)]; 
      it = node_to_nodes.find(*(pS+1));
      it->second.first=it->second.second=IDX_NONE;
    }

    if (it->second.first == IDX_NONE)
      node_to_nodes[*(pS+1)].first=*pS;
    else
      node_to_nodes[*(pS+1)].second=*pS;
  }
  return 0;
}

//=============================================================================
// Manifold version of the above function : fails if non-manifold contour is provided upon entry.
E_Int
BARSplitter::getNodesNeighBouring
(const K_FLD::IntArray& connectE2, NUGA::int_pair_vector_type& node_to_nodes)
{
  NUGA::int_vector_type           nodes, count;
  K_FLD::IntArray::const_iterator       pS;
  E_Int                                 NBELTS(connectE2.cols()), e, n;
  std::pair<E_Int, E_Int>               pnone(IDX_NONE, IDX_NONE);
  
  node_to_nodes.clear();

  connectE2.uniqueVals(nodes);
  node_to_nodes.resize(*std::max_element(nodes.begin(), nodes.end())+1, pnone);
  count.resize(node_to_nodes.size(), 0);

  for (e = 0; e < NBELTS; ++e)
  {
    pS = connectE2.col(e);
    for (n = 0; n < 2; ++n)
    {
      if (++count[*(pS+n)] > 2)// Deal only with manifold nodes.
        return 1;
    }
    node_to_nodes[*pS].second = *(pS+1);
    node_to_nodes[*(pS+1)].first = *pS;
  }
  return 0;
}

//=============================================================================
void
BARSplitter::split_periodic
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& c1, const K_FLD::IntArray& c2,
 E_Int N00, E_Int N01, E_Int N10, E_Int N11, std::vector<K_FLD::IntArray>& cs)
{
  E_Int n00(IDX_NONE), n01(IDX_NONE), n10(IDX_NONE), n11(IDX_NONE);

  assert (pos.rows() >= 3); // fixme : should be dim instead of rows.

  // Clear output.
  cs.clear();

  // Find out which nodes belong to which contour.

  if (__isBARNode(c1, N00))
  {
    n00 = N00;
    n01 = N01;
  }
  if ((n00 == IDX_NONE) && __isBARNode(c1, N01))
  {
    n00 = N01;
    n01 = N00;
  }
  if (__isBARNode(c1, N10))
  {
    n10 = N10;
    n11 = N11;
  }
  if ((n10 == IDX_NONE) && __isBARNode(c1, N11))
  {
    n10 = N11;
    n11 = N10;
  }

  bool done = false;

  done =  ( (n00 == IDX_NONE) || (n01 == IDX_NONE) || (n10 == IDX_NONE) || (n11 == IDX_NONE)) 
    || !__isBARNode(c2, n01) || !__isBARNode(c2, n11);

  if (done)
    return;

  // Do the split.

  std::vector<K_FLD::IntArray> cOut1, cOut2;

  __split(pos, 3, c1, n00, n10, cOut1, false/*enclose*/);
  if (cOut1.size() != 2)
    return;// failed to split into 2 bits
  
  assert (cOut1.size() == 2);
  __split(pos, 3, c2, n01, n11, cOut2, false/*enclose*/); 
  if (cOut2.size() != 2)
    return;// failed to split into 2 bits

  // Join the bits.

  E_Int orient00 = __getOrientation(cOut1[0], n00);
  E_Int orient01 = __getOrientation(cOut2[0], n01);

  assert ((orient00 != -1) && (orient01 != -1));

  E_Int E1[]  = {n00, n01};
  E_Int rE1[] = {n01, n00};
  E_Int E2[]  = {n10, n11};
  E_Int rE2[] = {n11, n10};

  cs.resize(2);
  cs[0] = cOut1[0];
  cs[1] = cOut1[1];

  if (orient00 == 0) // cOut1[0] starts from n00 and ends to n10. Opposite situation for cOut1[1].
  {
    cs[0].pushBack(rE1, rE1+2);
    cs[0].pushBack(E2, E2+2);
    cs[1].pushBack(E1, E1+2);
    cs[1].pushBack(rE2, rE2+2);
    if (orient01 == 0) // cOut2[0] starts from n01 (and ends to n11). Opposite situation for cOut2[1].
    {
      cs[0].pushBack(cOut2[1]);
      cs[1].pushBack(cOut2[0]);
    }
    else
    {
      cs[0].pushBack(cOut2[0]);
      cs[1].pushBack(cOut2[1]);
    }
  }
  else 
  {
    cs[0].pushBack(E1, E1+2);
    cs[0].pushBack(rE2, rE2+2);
    cs[1].pushBack(rE1, rE1+2);
    cs[1].pushBack(E2, E2+2);
    if (orient01 == 1) //
    {
      cs[0].pushBack(cOut2[1]);
      cs[1].pushBack(cOut2[0]);
    }
    else
    {
      cs[0].pushBack(cOut2[0]);
      cs[1].pushBack(cOut2[1]);
    }
  }
}


//=============================================================================
void
BARSplitter::split_periodic
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& ics,
 const K_FLD::IntArray& cuttingEdges, const std::vector<E_Int>& colors, 
 std::vector<K_FLD::IntArray>& ocs, K_FLD::IntArray& unusedEdges)
{
  bool carry_on = true, restart = false;
  
  std::vector<K_FLD::IntArray> cOut;
  NUGA::non_oriented_edge_set_type edges;
  for (E_Int i = 0; i < cuttingEdges.cols(); ++i)
    edges.insert(K_MESH::NO_Edge(cuttingEdges(0,i), cuttingEdges(1,i)));

  std::vector<E_Int> current_color = colors;

  ocs = ics;

  E_Int Ns[4], count = 0;

  while (carry_on)
  {
    restart = false;
    E_Int ocsSize = ocs.size();
    for (E_Int i = 0; i < ocsSize; ++i)
    {
      K_FLD::IntArray& c1 = ocs[i];

      for (E_Int j = i+1; j < ocsSize; ++j)
      {
        K_FLD::IntArray& c2 = ocs[j];

        count = 0;
        for (NUGA::non_oriented_edge_set_type::const_iterator it = edges.begin(); (it != edges.end()) && (count < 4); ++it)
        {
          if ((current_color[it->node(0)] == i) && (current_color[it->node(1)] == j))
          {
            Ns[count++] = it->node(0);
            Ns[count++] = it->node(1);
          }
          else if ((current_color[it->node(1)] == i) && (current_color[it->node(0)] == j))
          {
            Ns[count++] = it->node(1);
            Ns[count++] = it->node(0);
          }
        }

        if (count != 4) continue;

        //std::cout << " found" << std::endl;

        split_periodic(pos, c1, c2, Ns[0], Ns[1], Ns[2], Ns[3], cOut);

        restart = !cOut.empty();
        if (restart)
        {
          ocs[i] = cOut[0];
          ocs[j] = cOut[1];
          //update colors
          for (E_Int k = 0; k < ocs[i].cols(); ++k)
            current_color[ocs[i](0,k)] = current_color[ocs[i](1,k)] = i;
          for (E_Int k = 0; k < ocs[j].cols(); ++k)
            current_color[ocs[j](0,k)] = current_color[ocs[j](1,k)] = j;

          edges.erase(K_MESH::NO_Edge(Ns[0], Ns[1]));
          edges.erase(K_MESH::NO_Edge(Ns[2], Ns[3]));

          break;
        }
      }
      if (restart)
        break;
    }
    carry_on = restart && !edges.empty();
  }
  if (!edges.empty())
  {
    E_Int E[2];
    for (NUGA::non_oriented_edge_set_type::iterator it = edges.begin(); it != edges.end(); ++it)
    {
      E[0] = it->node(0);
      E[1] = it->node(1);
      unusedEdges.pushBack(E, E+2);
    }
  }
}

//=============================================================================
E_Float BARSplitter::__getAngle(const E_Float* n1, const E_Float* n2)
{
  //WARNING : n1 and n2 are asummed to be unitary
  E_Float n[3], c, s;
  NUGA::crossProduct<3>(n1, n2, n);
  s = NUGA::normalize<3>(n);
  c = NUGA::dot<3>(n1, n2);
  return ::atan2(s, c);
}

//=============================================================================
void BARSplitter::split_loops
(const K_FLD::IntArray& connectB, const K_FLD::FloatArray& coord, std::vector<K_FLD::IntArray> & loops)
{
  //WARNING : assume that loops have sharped angle at concentration node (to be able to associate them).
  //WARNING : designed for ngon_t::close_phs use or similar config.
  
  loops.clear();
  
  std::set<K_MESH::NO_Edge> edges;
  for (E_Int i = 0; i < connectB.cols(); ++i)
  {
    const E_Int& N0=connectB(0,i);
    const E_Int& N1=connectB(1,i);
    edges.insert(K_MESH::NO_Edge(N0,N1));
  }
  
  typedef std::map<E_Int, std::vector<E_Int> > ntn_t;
  ntn_t node_to_nodes;
  BARSplitter::get_node_to_nodes(connectB, node_to_nodes);
  
#ifdef DEBUG_BARSPLITTER
  for (ntn_t::const_iterator it=node_to_nodes.begin(); it != node_to_nodes.end(); ++it)
  {
    assert(it->second.size()!=0 && it->second.size()%2==0);
    std::cout << it->first << ": ";
    for (size_t i=0;i<it->second.size();++i) std::cout << it->second[i] << " ";
    std::cout << std::endl;
  }
#endif
  
  E_Int c=-1, N, Nnext=0, Nprev;
  std::set<K_MESH::NO_Edge>::iterator itE;
  K_MESH::Edge E;
  K_MESH::NO_Edge noE;
  E_Float n1[3], n2[3], mina, alpha;
  
  while (!edges.empty())
  {
    ++c;
    loops.resize(c+1);
    mina = NUGA::FLOAT_MAX;
    
    itE = edges.begin();
    
    E_Int N0=itE->node(0);
    E_Int N1=itE->node(1);
    
    loops[c].pushBack(itE->begin(), itE->end());
    edges.erase(*itE);
    
    N = N1;
    Nprev = N0;
    
    while (N != N0)
    {
      E_Int sz = node_to_nodes[N].size();
      assert (sz !=0 && sz%2==0);
      
      if (sz == 2)
        Nnext = (node_to_nodes[N][0] != Nprev) ? node_to_nodes[N][0] : node_to_nodes[N][1];
      else //GEOM distinction
      {
        NUGA::diff<3>(coord.col(Nprev), coord.col(N), n1);
        NUGA::normalize<3>(n1);
        for (E_Int i = 0; i < sz; ++i)
        {
          if (node_to_nodes[N][i] == Nprev) continue;
          
          NUGA::diff<3>(coord.col(node_to_nodes[N][i]), coord.col(N), n2);
          NUGA::normalize<3>(n2);
          
          alpha = ::fabs(__getAngle(n1, n2));
          if (alpha < mina)
          {
            mina = alpha;
            Nnext = node_to_nodes[N][i];
          }
        }
      }
      
      E.setNodes(N, Nnext);
      loops[c].pushBack(E.begin(), E.end());
      noE.setNodes(N, Nnext);
      edges.erase(noE);
      Nprev = N;
      N = Nnext;
      
    }
  }
}

E_Int BARSplitter::min_angle_node(const std::vector<E_Int>& sorted_nodes, const K_FLD::FloatArray& coord)
{
  E_Int Nmin(IDX_NONE);
  E_Int sz(sorted_nodes.size());
  E_Int iprev=sz-1;
  E_Float n1[3], n2[3], mina(NUGA::FLOAT_MAX), alpha;
  
  for (E_Int i = 0; i < sz; ++i)
  {
    const E_Int& Nprev = sorted_nodes[iprev];
    const E_Int& N = sorted_nodes[i];
    const E_Int& Nnext = sorted_nodes[(i+1)%sz];
    
    NUGA::diff<3>(coord.col(Nprev), coord.col(N), n1);
    NUGA::diff<3>(coord.col(Nnext), coord.col(N), n2);
    
    alpha = ::fabs(__getAngle(n1, n2));
    if (alpha < mina)
    {
      mina = alpha;
      Nmin = N;
    }
    
    iprev=i;
  }
  
  return Nmin;
  
}
