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

#include "Nuga/include/PatchMaker.h"
#include "Nuga/include/Zipper.h"
#include "Nuga/include/ContourSplitter.h"
#include "Nuga/include/BARSplitter.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/MeshTool.h"

#ifdef E_TIME
#include "Nuga/include/chrono.h"
#endif
#ifdef WIN32
#ifdef E_DEBUG
#include "meshIO/meshIO.h"
#include "meshIO/medit.h"
#include <sstream>
#endif
#endif

#ifdef DEBUG_PATCHMAKER
#include "IO/DynArrayIO.h"
#endif

using namespace NUGA;

void
PatchMaker::run
(const K_FLD::FloatArray& posS, const K_FLD::IntArray& connectS, const int_vector_type & pairs,
 E_Float angle_tolerance, std::vector<K_FLD::IntArray> & connectBout)
{
#ifdef E_TIME
  NUGA::chrono c;
  c.start();
#endif
  K_FLD::IntArray connectBin;
  /* Split the connectivity into connex bits. */
  // Set a color for each surface.

  std::vector<K_FLD::IntArray> cOutS, cOut0, temp;
  int_vector_type surface_colors;
  non_oriented_edge_set_type dummyS;
  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(connectS, dummyS, cOutS);
  
#ifdef DEBUG_PATCHMAKER
  {
    std::vector<E_Int> colz;
    K_FLD::IntArray conn;
    for (size_t i=0; i< cOutS.size(); ++i)
    {
      conn.pushBack(cOutS[i]);
      colz.resize(conn.cols(), i);
    }
    K_CONVERTER::DynArrayIO::write("colored_surfaces.mesh", posS, conn, "TRI", 0, &colz);
  }
#endif

#ifdef E_TIME
  std::cout << "split surface : " << c.elapsed() << std::endl;
  std::cout << "split surface : nb surfaces " << cOutS.size() << std::endl;
  c.start();
#endif

  int_set_type dummy;
  K_FLD::IntArray connectb;
  for (size_t scol = 0; scol < cOutS.size(); ++scol)
  {
    NUGA::MeshTool::getBoundary(cOutS[scol], connectb);
    connectBin.pushBack(connectb);
    
#ifdef DEBUG_PATCHMAKER
    //K_CONVERTER::DynArrayIO::write("connectB.mesh", posS, connectb, "BAR");
#endif
    
    ContourSplitter<K_MESH::Edge, E_Int>::splitConnectivity(connectb, dummy, temp);
    cOut0.insert(cOut0.end(), temp.begin(), temp.end());
    // Set the surface color.
    surface_colors.resize(surface_colors.size()+temp.size(), scol);
  }
#ifdef E_TIME
  std::cout << "split contours : " << c.elapsed() << std::endl;
  c.start();
#endif

  /* Sort the nodes for each contours and set the contours color.*/
  E_Int nb_contours = cOut0.size();
  // Resize the colors vector.
  int_vector_type colors, nodes;
  connectBin.uniqueVals(nodes);
  if (!nodes.empty())
    colors.resize(*std::max_element(nodes.begin(), nodes.end())+1, IDX_NONE);
  // Do the sorting and coloring.
  std::vector<int_vector_type> sorted_nodes(nb_contours); 

  for (E_Int c = 0; c < nb_contours; ++c)
  {
    
    BARSplitter::getSortedNodes(cOut0[c], sorted_nodes[c]);
   
    for (size_t n = 0; n < sorted_nodes[c].size(); ++n) // assigns colors.
      colors[sorted_nodes[c][n]] = c;
  }
#ifdef E_TIME
  std::cout << "get sorted nodes : " << c.elapsed() << std::endl;
  c.start();
#endif

  /* Compute normals. */
  K_FLD::FloatArray normals;
  NUGA::MeshTool::computeNodeNormals(posS, connectS, normals);
  
  // Update the normals with the pairs info.
  //E_Int err = __update_normals(posS, connectBin, pairs, normals);
  //if (err)
    //return;
  
#ifdef E_TIME
  std::cout << "normals : " << c.elapsed() << std::endl;
  c.start();
#endif


  /* Flag the nodes where the cuts needs to be done.*/
  int_vector_type critical_nodes;
  __flag_critical_nodes(normals, angle_tolerance, pairs, colors, sorted_nodes, critical_nodes);

#ifdef E_TIME
  std::cout << "flag critical nodes : " << c.elapsed() << std::endl;
  c.start();
#endif


  /* Build the cuttingEdges.*/
  K_FLD::IntArray cuttingEdges, periodicEdges;
  __build_cutting_edges(pairs, colors, critical_nodes, cuttingEdges, periodicEdges);
  
#ifdef E_TIME
  std::cout << "build cutting edges : " << c.elapsed() << std::endl;
#endif


#ifdef WIN32
#ifdef E_DEBUG
  meshIO::write("periodicEdges.mesh", posS, periodicEdges);
  meshIO::write("cuttingedges.mesh", posS, cuttingEdges);
#endif
#endif

#ifdef E_TIME
  c.start();
#endif

  /* Build the patches. */
  std::vector<K_FLD::IntArray> cOut1;
  K_FLD::IntArray unusedEdges;
  // Join the pairs of periodic contour.
  
  BARSplitter::split_periodic(posS, cOut0, periodicEdges, colors, cOut1, unusedEdges);
  cuttingEdges.pushBack(unusedEdges); // Append the unused nedges to the cutting edges.

#ifdef E_TIME
  std::cout << "split periodic : " << c.elapsed() << std::endl;
  c.start();
#endif

  // Split the simple contours.
  BARSplitter::split(posS, 3, cOut1, cuttingEdges, connectBout);

#ifdef E_TIME
  std::cout << "split simple : " << c.elapsed() << std::endl;
  c.start();
#endif

#ifdef WIN32
#ifdef E_DEBUG

  K_FLD::IntArray cc;
  std::vector<E_Int> colorz;
  E_Int col = 0;
  for (size_t i = 0; i < connectBout.size(); ++i)
  {
    cc.pushBack(connectBout[i]);
    std::ostringstream o;
    ++col;
    if (col == meshIO::medit::RED)
      ++col;
    colorz.resize(cc.cols(), col);
    //o << "patch_" << i << ".mesh";
    //meshIO::write(o.str().c_str(), posS, connectBout[i]);
  }

  // std::cout << "number of patches " << connectBout.size() << std::endl;

  // Add the remaining mates. (in red)
  E_Int E[2], mate, cols(cc.cols());
  K_FLD::IntArray cc1 = cc;
  E_Int psize = E_Int(pairs.size());
  for (E_Int i = 0; i < cols; ++i)
  {
    for (E_Int k = 0; k < 2; ++k)
    {
      if ((cc(k,i) < 0) || (cc(k,i) >= psize))
        continue;
      mate = pairs[cc(k,i)];
      if ((mate < 0) || (mate == IDX_NONE))
        continue;

      E[0] = cc(k,i);
      E[1] = mate;
      cc1.pushBack(E, E+2);
    }
  }
  
  colorz.resize(cc1.cols(), meshIO::medit::RED);
   
  meshIO::write("patches.mesh", posS, cc1, colorz);
  
#endif
#endif
}

///
E_Int
PatchMaker::__update_normals
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
 const std::vector<E_Int> &pairs, K_FLD::FloatArray& normals)
{
  E_Int                           Ni, Nj, B0, B1, err;
  size_t                          i, nb_nodes, nb_bound(connectB.cols());
  K_FLD::IntArray::const_iterator pS;
  E_Float                         B[3], E[3], W[3], L, *p;
  NUGA::int_vector_type     nodesS;
  NUGA::int_set_type        hN;
  K_FLD::FloatArray               new_normals(normals);

  // Fast returns
  if (pos.cols() == 0)      return 1;
  if (connectB.cols() == 0) return 2;
  if (pairs.empty())        return 3;
  if (normals.cols() == 0)  return 4;

  // Get node neighbors
  NUGA::int_pair_vector_type     node_to_nodes;
  err = BARSplitter::getNodesNeighBouring(connectB, node_to_nodes);
  if (err)
    return err;
  connectB.uniqueVals(nodesS);

  // Get the nodes for which the normal can be recomputed.
  // i.e nodes with a mate.
  nb_nodes = nodesS.size();
  for (i = 0; i < nb_nodes; ++i)
  {
    Ni = nodesS[i];
    Nj = pairs[Ni];
    if ((Nj < 0) || (Nj == IDX_NONE))
      continue;
    if (Nj == node_to_nodes[Ni].first)
      continue;
    if (Nj == node_to_nodes[Ni].second)
      continue;
    hN.insert(Ni);
    for (E_Int k = 0; k < 3; ++k) // reset the normal.
      new_normals(k, Ni) = 0.;
  }

  // Add the new contributions
  for (size_t i = 0; i < nb_bound; ++i)
  {
    pS = connectB.col(i);
    B0 = *pS;
    B1 = *(pS+1);
    NUGA::diff<3>(pos.col(B1), pos.col(B0), B);

    for (E_Int n = 0; n < 2; ++n)
    {
      Ni = *(pS+n);
      if (hN.find(Ni) == hN.end())
        continue;
      Nj = pairs[Ni];
      NUGA::diff<3>(pos.col(Nj), pos.col(Ni), E);

      NUGA::crossProduct<3>(E, B, W);

      L = NUGA::normalize<3>(W);
      if (L < EPSILON)
      {
        for (E_Int k = 0; k < 3; ++k)
          new_normals(k, Ni) = 0.;
        hN.erase(Ni);
      }
      else
      {
        for (E_Int k = 0; k < 3; ++k)
          new_normals(k, Ni) += W[k];
      }
    }
  }

  // Normalize and assign
  for (std::set<E_Int>::const_iterator i = hN.begin(); i != hN.end(); ++i)
  {
    Ni = *i;
    p = new_normals.col(Ni);
    L  = NUGA::normalize<3>(p); 
    if (L != 0.)
     for (E_Int k = 0; k < 3; ++k)
       normals(k, Ni) = new_normals(k, Ni);
  }

  return 0;
}

///
E_Float PatchMaker::__getAngle(const E_Float* n1, const E_Float* n2)
{
  E_Float n[3], c, s;
  NUGA::crossProduct<3>(n1, n2, n);
  s = NUGA::normalize<3>(n);
  c = NUGA::dot<3>(n1, n2);

  return atan2(s, c);
}

///
void
PatchMaker::__flag_critical_nodes
(const K_FLD::FloatArray& normals,
 E_Float critical_angle,
 const std::vector<E_Int> & pairs,
 const std::vector<E_Int> & colors,
 std::vector< std::vector<E_Int> > & sorted_nodes,
 std::vector<E_Int> & critical_nodes)
{
  // Stamp the critical nodes on each contour and synchronize on others.
  E_Int nb_contours = sorted_nodes.size();
  NUGA::bool_vector_type::iterator itC;
  E_Int N0, prev, next, Ni, Np, Nn, c;
  size_t nb_nodes, n;

  critical_nodes.resize(pairs.size(), 0);

  // Stamp border nodes.
  for (c = 0; c < nb_contours; ++c)
  {
    nb_nodes = sorted_nodes[c].size();
    if (nb_nodes == 0)
      continue; // fixme : error ?

    for (n = 0; (n < nb_nodes); ++n)
    {
      Ni = sorted_nodes[c][n];
      
      if ((pairs[Ni] < 0) || (pairs[Ni] == IDX_NONE)) // border nodes have a mate.
        continue;
      prev = (n != 0) ? n-1 : nb_nodes-1;
      Np = sorted_nodes[c][prev];
      next = (n+1)%nb_nodes;
      Nn = sorted_nodes[c][next];
      //if ((pairs[Np] < 0) || (pairs[Np] == IDX_NONE))
      if (pairs[Np] == Zipper::FREE)
        critical_nodes[Ni] = 1;
      //else if ((pairs[Nn] < 0) || (pairs[Nn] == IDX_NONE))
      else if (pairs[Nn] == Zipper::FREE)
        critical_nodes[Ni] = 1;
    }
  }

  for (c = 0; c < nb_contours; ++c)
  {
    N0 = IDX_NONE;
    nb_nodes = sorted_nodes[c].size();
    if (nb_nodes == 0)
      continue; // fixme : error ?
    
    for (n = 0; (n < nb_nodes) && (N0 == IDX_NONE); ++n)
    {
      if (critical_nodes[sorted_nodes[c][n]] == 1)
        N0 = sorted_nodes[c][n];
    }

    for (n = 0; (n < nb_nodes) && (N0 == IDX_NONE); ++n)
    {
      if ((pairs[sorted_nodes[c][n]] != IDX_NONE) && (pairs[sorted_nodes[c][n]] >= 0))
        N0 = sorted_nodes[c][n];
    }
   
    if (N0 == IDX_NONE)
      N0 = sorted_nodes[c][0];

    critical_nodes[N0] = 1;

    __flag_critical_nodes_on_contour(sorted_nodes[c], N0, normals, critical_angle, pairs, colors, critical_nodes);

    // Propagate the critical nodes.
    nb_nodes = critical_nodes.size();
    for (n = 0; n < nb_nodes; ++n)
    {
      if ((critical_nodes[n] == 1) && (pairs[n] != IDX_NONE) && (pairs[n] >= 0))
        critical_nodes[pairs[n]] = 1;
    }
  }
}

///
E_Int
PatchMaker::__flag_critical_nodes_on_contour
(const std::vector<E_Int>& nodes, E_Int N0, const K_FLD::FloatArray& normals,
 E_Float critical_angle, const std::vector<E_Int>& pairs, const std::vector<E_Int>& colors,
 std::vector<E_Int> & flag)
{
  E_Int istart = 0, Ni, Nj, Nstart = N0;
  E_Float alpha;

  std::map<E_Int, E_Int> colors_to_join;
  std::map<E_Int, E_Int>::iterator it;

  while (nodes[istart] != N0){++istart;}

  size_t k, sz = nodes.size();
  for (k = 1; k < sz; ++k)
  {
    const E_Float* D = normals.col(Nstart);

    Ni = nodes[(istart+k)%sz];

    if (Ni == N0) // Done
      break;

    Nj = pairs[Ni];

    if ((Nj == IDX_NONE) || (Nj < 0))
      continue;

    if (flag[Ni] == false)
    {
      alpha = __getAngle(D, normals.col(Ni));
      alpha = (alpha < 0.) ? -alpha : alpha;
      if ( alpha > critical_angle)
        flag[Ni] = true;
    }

    if (flag[Ni] == true)
    {
      Nstart = Ni;
      // Ensure to have 2 attach for each pair of peridodic boundary.
      it = colors_to_join.find(colors[Nj]);
      if (it == colors_to_join.end()) 
        colors_to_join[colors[Nj]] = 1;
      else
        ++colors_to_join[colors[Nj]];
    }
  }

  // Ensure to have 2 attach for each pair of peridodic boundary.
  for (k = 0; k < sz; ++k)
  {
    Ni = nodes[k];

    if (flag[Ni] == true)
      continue;

    Nj = pairs[Ni];

    if ((Nj == IDX_NONE) || (Nj < 0))
      continue;

    it = colors_to_join.find(colors[Nj]);
    if ((it != colors_to_join.end()) && (it->second > 1))
      continue;

    flag[Ni] = true;
    if (it == colors_to_join.end()) 
      colors_to_join[colors[Nj]] = 1;
    else
      ++colors_to_join[colors[Nj]];
  }

  return 0;
}

///
void
PatchMaker::__build_cutting_edges
(const std::vector<E_Int>& pairs, const std::vector<E_Int>& colors,
 const std::vector<E_Int>& critical_nodes,
 K_FLD::IntArray& cuttingEdges, K_FLD::IntArray& periodicEdges)
{
  std::map<K_MESH::NO_Edge, E_Int> counter;
  std::map<K_MESH::NO_Edge, E_Int>::iterator it;
  NUGA::non_oriented_edge_set_type edges;

  K_MESH::NO_Edge E;
  E_Int Ei[2];
  E_Int Ni,Nj;
  for (size_t i = 0; i < critical_nodes.size(); ++i)
  {
    if (critical_nodes[i] == 0)
      continue;
    Ni = i;
    Nj = pairs[Ni];
    if ((Nj == IDX_NONE) || (Nj < 0))//fixme
      continue;
    //assert (Nj != IDX_NONE);
    E.setNodes(Ni, Nj);
    Ei[0] = E.node(0);
    Ei[1] = E.node(1);

    if (!edges.insert(E).second) // already taken into account.
      continue;

    E.setNodes(colors[Ei[0]], colors[Ei[1]]);
    it = counter.find(E);
    if (it == counter.end())
    {
      counter[E] = 1;
      periodicEdges.pushBack(Ei, Ei+2);
    }
    else if (it->second < 2)
    {
      ++it->second;
      periodicEdges.pushBack(Ei, Ei+2);
    }
    else
      cuttingEdges.pushBack(Ei, Ei+2);
  }
}

