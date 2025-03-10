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

#ifndef _BARSPLITTER_H_
#define _BARSPLITTER_H_

#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/EltAlgo.h"
#include <deque>

class BARSplitter
{
public:
  typedef NUGA::int_set_type   int_set_type;
public:
  // Splits the connectivity (only) with the edge N0N1 and stores the bits into separate containers.
  /* WARNING: it is the responsibility of the caller to destroy connectBout.*/
  static void split (const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& connectBin, E_Int N0, E_Int N1,
                     std::vector<K_FLD::IntArray> & connectBout, E_Bool enclose = 1, E_Float tolerance = EPSILON);

  /* Splits the connectivity with the edge N0N1 and stores the connectivity bits and corresponding coordinates
      into separate containers.*/
   /* WARNING: it is the responsibility of the caller to destroy posOut and connectBout.*/
  static void split (const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& connectBin, E_Int N0, E_Int N1,
                     std::vector<K_FLD::FloatArray> & posOut, std::vector<K_FLD::IntArray> & connectBout,
                     E_Bool enclose = 1, E_Float tolerance = EPSILON);

  // Splits the connectivity with the input edges list and stores the bits into separate containers.
  /* WARNING: it is the responsibility of the caller to destroy connectBout. */
  static void split (const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& connectBin, const K_FLD::IntArray& cuttingEdges,
                     std::vector<K_FLD::IntArray> & connectBout, E_Float tolerance = EPSILON);

  // Splits the connectivity with the input edges list and stores the bits into separate containers.
  /* WARNING: it is the responsibility of the caller to destroy connectBout. */
  static void split (const K_FLD::FloatArray& pos, E_Int dim, const std::vector<K_FLD::IntArray>& connectBin,
                     const K_FLD::IntArray& cuttingEdges, std::vector<K_FLD::IntArray> & connectBout,
                     E_Float tolerance = EPSILON);
  // Separate single loops (TOPO algo)
  static void split_loops (const K_FLD::IntArray& connectBin, const K_FLD::FloatArray& coord, std::vector<K_FLD::IntArray> & connectBout);


  static void split_periodic (const K_FLD::FloatArray& pos, const K_FLD::IntArray& ic1, const K_FLD::IntArray& ic2,
                              E_Int N00, E_Int N01, E_Int N10, E_Int N11,
                              std::vector<K_FLD::IntArray>& cs);

  static void split_periodic (const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& ics,
                              const K_FLD::IntArray& cuttingEdges, const std::vector<E_Int>& colors,
                              std::vector<K_FLD::IntArray>& cs,
                              K_FLD::IntArray& unusedEdges);

  static E_Int getNodesNeighBouring(const K_FLD::IntArray& connectE2, NUGA::int_pair_vector_type& node_to_nodes);
  
  static E_Int get_node_to_nodes(const K_FLD::IntArray& connectE2, std::map< E_Int, std::vector<E_Int> >& node_to_nodes);

  // for closed contours
  template <typename Vect> 
  static E_Int getSortedNodes(const K_FLD::IntArray& connectB, Vect& nodes);
  
  // for both closed/open lines
  template <typename Vect>
  static E_Int getSortedChainNodes(const K_FLD::IntArray& connectB, Vect& nodes);
  
  ///
  static E_Int get_node_to_nodes(const K_FLD::IntArray& connectE2, std::map< E_Int, std::pair<E_Int, E_Int> >& node_to_nodes);

  ///
  template<typename EdgeCONT> static void split_eset_into_manifold_chains(const EdgeCONT& edges, std::vector<std::deque<E_Int>>& chains, std::vector<E_Int>& edgecol);
  
  ///
  template<typename EdgeCONT, typename IdCONT>
  static void compute_refinement(const EdgeCONT& free_edges, const IdCONT &chain,
                          std::map < K_MESH::NO_Edge, IdCONT>& edge_to_refine_nodes);
   
  // Returs the node having the smallest angle (absolute value).
  static E_Int min_angle_node(const std::vector<E_Int>& sorted_nodes, const K_FLD::FloatArray& coord);

private:
  BARSplitter(void){}
  ~BARSplitter(void){}

  ///
public://fixme!!!
  static void __split (const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& BAR, E_Int N0, E_Int N1,
                       std::vector<K_FLD::IntArray> & connectBout, E_Bool enclose, E_Float tolerance = EPSILON);
private:
  static E_Bool __fastReturn(const K_FLD::IntArray& connectBin, E_Int N0, E_Int N1);

  static E_Bool __isBARNode(const K_FLD::IntArray& connectB, E_Int N);

  static E_Bool __isBAREdge(const K_FLD::IntArray& connectB, E_Int N0, E_Int N1);

  static E_Bool __canBeClosed(const K_FLD::FloatArray& pos, E_Int dim, const K_FLD::IntArray& connectB,
                              E_Int N0, E_Int N1, E_Float tolerance);

  static E_Int __getOrientation(const K_FLD::IntArray& connectB, E_Int N0);
  
  static E_Float __getAngle(const E_Float* n1, const E_Float* n2);
  
  
  
};

//=============================================================================
template <typename Vect>
E_Int
BARSplitter::getSortedNodes // assume a closed polyline : DO NOT ADD TWICE A NODE TO CLOSE
(const K_FLD::IntArray& connectB, Vect& nodes)
{
  std::map< E_Int, std::pair<E_Int, E_Int> >    node_to_nodes;
  E_Int  Ni, Nnext(IDX_NONE);
  E_Int  Nprev, err(0);

  nodes.clear();

  err = get_node_to_nodes(connectB, node_to_nodes);
  if (err) return err;
  
  // starting as connectB orientation : if connectB is consistently oriented, it will be the same orientation upon exit
  nodes.push_back(connectB(0, 0));
  nodes.push_back(connectB(1, 0));
  Nprev = nodes[0];
  //
  for (E_Int i = 1; (i < connectB.cols()-1); ++i)
  {
    Ni = nodes[i];
    Nnext = node_to_nodes[Ni].second;
    if (Nnext==Nprev)
      Nnext = node_to_nodes[Ni].first;
    nodes.push_back(Nnext);
    Nprev=Ni;
  }
  return 0;
}

//=============================================================================
template <typename Vect>
E_Int
BARSplitter::getSortedChainNodes // open polyline : ADD TWICE THE NODE TO DISTINGUISH CLOSED/OPEN
(const K_FLD::IntArray& connectB, Vect& nodes)
{
  std::map< E_Int, std::pair<E_Int, E_Int> >    node_to_nodes;
  E_Int  Ni, Nstart(IDX_NONE), Nnext(IDX_NONE);
  E_Int  Nprev(IDX_NONE), err(0);

  nodes.clear();

  err = get_node_to_nodes(connectB, node_to_nodes);
  if (err)
    return err;

  Nstart = connectB(0,0);
  
  //nodes.reserve(connectB.cols());
  nodes.push_back(Nstart);
  
  E_Int nb_nodes = node_to_nodes.size();
  //
  for (E_Int i = 0; (i < nb_nodes); ++i)
  {
    Ni = nodes[i];
    Nnext = node_to_nodes[Ni].second;
    if (Nnext==Nprev)
      Nnext = node_to_nodes[Ni].first;
    if (Nnext == IDX_NONE)
      break;
    nodes.push_back(Nnext);
    Nprev=Ni;
  }
  
  if (Nnext != IDX_NONE) // therefore close
  {
    assert (Nnext == nodes[0]);
    return 0;
  }
  
  Nprev=nodes[1];
  for (E_Int i = 0; (i < nb_nodes); ++i)
  {
    Ni = nodes[0];
    Nnext = node_to_nodes[Ni].second;
    if (Nnext==Nprev)
      Nnext = node_to_nodes[Ni].first;
    if (Nnext == IDX_NONE)
      break;
    nodes.push_front(Nnext);
    Nprev=Ni;
  }
  
  //nodes.insert(nodes.end(), tmp.begin(), tmp.end());
  
  return 0;
}

template<typename EdgeCONT>
void BARSplitter::split_eset_into_manifold_chains(const EdgeCONT& edges, std::vector<std::deque<E_Int>>& chains, std::vector<E_Int>& edgecol)
{
  K_FLD::IntArray cntE, neighborsE;
  std::vector<E_Int> colorsE;
  for (auto it = edges.begin(); it != edges.end(); ++it)
    cntE.pushBack(it->begin(), it->end());

  NUGA::template EltAlgo<K_MESH::Edge>::getManifoldNeighbours(cntE, neighborsE, false);
  NUGA::template EltAlgo<K_MESH::Edge>::coloring(neighborsE, colorsE);

  // 2.b. classify chains : type 0 (closed), type 1 (open polyline)
  E_Int nb_chains = *std::max_element(colorsE.begin(), colorsE.end()) + 1;
  chains.resize(nb_chains);
  std::vector<K_FLD::IntArray> cnt_chains(nb_chains);

  edgecol.clear();
  edgecol.resize(edges.size());

  //2.c. dispacth by color and create the poly lines
  for (E_Int i = 0; i < cntE.cols(); ++i)
  {
    cnt_chains[colorsE[i]].pushBack(cntE.col(i), cntE.col(i) + 2);
    edgecol[i] = colorsE[i];
  }
  for (E_Int c = 0; c < nb_chains; ++c)
    BARSplitter::getSortedChainNodes(cnt_chains[c], chains[c]);
}

///
template<typename EdgeCONT, typename IdCONT>
void BARSplitter::compute_refinement
(const EdgeCONT& edges, const IdCONT &chain,
 std::map < K_MESH::NO_Edge, IdCONT>& edge_to_refine_nodes)
{
  edge_to_refine_nodes.clear();

  for (auto& e : edges)
  {
    E_Int I[] = {0, 0};
    E_Int count(-1);

    for (size_t i = 0; i < chain.size(); ++i)
    {
      if (chain[i] != e.node(0) && chain[i] != e.node(1)) continue;
      I[++count] = i;
    }

    assert(count == 1);

    for (E_Int j = I[0] + 1; j < I[1]; ++j)
      edge_to_refine_nodes[e].push_back(chain[j]);
  }
}

#endif
