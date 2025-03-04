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

#include "Nuga/include/NodeAssociator.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/ContourSplitter.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/Zipper.h"
#include "Nuga/include/BARSplitter.h"

#ifdef E_TIME
#include "Nuga/include/chrono.h"
#endif

#ifdef DEBUG_NODEASSOCIATOR
#include "IO/DynArrayIO.h"
#include <sstream>
#endif

bool NodeAssociator::reorient = false;


///
NodeAssociator::NodeAssociator(void)
{
}

///
NodeAssociator::~NodeAssociator(void)
{
}

///
void
NodeAssociator::make_pairs
(const K_FLD::FloatArray& pos, const std::vector< K_FLD::IntArray* >& components,
 std::vector<E_Int> &nmates, K_FLD::IntArray& OneSurface)
{
  nmates.clear();
  nmates.resize(pos.cols(), IDX_NONE);

  std::vector<E_Int> ncolors(pos.cols(), IDX_NONE);

  OneSurface.clear();

#ifdef E_TIME
  DELAUNAY::chrono c, glob;
  c.start();
  glob.start();
#endif

   // Intersections
  std::vector<K_FLD::IntArray> clean_comps;
  __removeOverlaps(pos, components, nmates, clean_comps);

  for (size_t comp = 0; comp < clean_comps.size(); ++comp)
    OneSurface.pushBack(clean_comps[comp]);

#ifdef WIN32
#ifdef E_DEBUG  
  meshIO::write("cleaned.mesh", pos, OneSurface);
#endif
#endif

  std::vector< std::vector<E_Int> > cnodes;
  std::vector< E_Int> scolors;
  std::vector<K_FLD::IntArray> surfaces;
  NUGA::non_oriented_edge_set_type dummyS;

  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(OneSurface, dummyS, surfaces);

  K_FLD::IntArray connectb;
  std::vector<K_FLD::IntArray> connectbs;
  NUGA::int_set_type dummy;
  std::vector<E_Int> nodes;
  E_Int color = 0;
  for (size_t s = 0; s < surfaces.size(); ++s)
  {
    NUGA::MeshTool::getBoundary(surfaces[s], connectb);
    connectbs.clear();
    ContourSplitter<K_MESH::Edge, E_Int>::splitConnectivity(connectb, dummy, connectbs);

    for (size_t c = 0; c < connectbs.size(); ++c)
    {
      connectbs[c].uniqueVals(nodes);
      cnodes.push_back(nodes);
      // Set the node color (contour color).
      for (size_t i = 0; i < nodes.size(); ++i)
        ncolors[nodes[i]] = color;
      ++color;
    }

    scolors.resize(scolors.size()+connectbs.size(), s);
  }

  __make_pairs(pos, cnodes, scolors, nmates);
}

///
void
NodeAssociator::__make_pairs
(const K_FLD::FloatArray& posB, const std::vector< std::vector<E_Int> >& contour_nodes,
 const std::vector<E_Int>& surface_colors, std::vector<E_Int> &pairs)
{
  // Fast returns.
  if (posB.cols() == 0) return;

  std::vector<std::pair<E_Float, std::pair<E_Int, E_Int> > > palmares;
  std::vector<E_Int> closest_color;
  E_Int nb_contours = contour_nodes.size();
  const E_Float* Pi;
  E_Int Ni, Nj;
  E_Float d, dmin;
  //pairs.clear();
  std::vector<E_Int> nodes;
  std::vector<K_SEARCH::KdTree<>* > trees;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> posAcc(posB);

#ifdef E_TIME
  DELAUNAY::chrono cc;
  cc.start();
  std::cout << "building trees" << std::endl;
#endif

  E_Int maxID = -1, id;
  for (E_Int i = 0; i < nb_contours; ++i)
  {
    //std::cout << i  << " " << sorted_nodes[i].size() << std::endl;
    trees.push_back(new K_SEARCH::KdTree<>(posAcc, contour_nodes[i]));
    id = *std::max_element(contour_nodes[i].begin(), contour_nodes[i].end());
    maxID = std::max(maxID, id);
  }

#ifdef E_TIME
  cc.start();
#endif

  //pairs.resize(maxID+1, IDX_NONE);
  closest_color.resize(maxID+1, IDX_NONE);
  NUGA::non_oriented_int_pair_set_type hset;
  K_MESH::NO_Edge Ei;

  for (E_Int i = 0; i < nb_contours; ++i)
  {
    const std::vector<E_Int>& LISTi = contour_nodes[i];

    for (size_t k = 0; k < LISTi.size(); ++k)
    {
      Ni = LISTi[k];
      Pi = posB.col(Ni);
      dmin = NUGA::FLOAT_MAX;

      for (E_Int j = 0; j < nb_contours; ++j)
      {
        if (i == j)
          continue;
        if (surface_colors[j] == surface_colors[i])
          continue;

        Nj = trees[j]->getClosest(Pi);  
        d = NUGA::sqrDistance(Pi, posB.col(Nj), 3);
        if (d < dmin)
        {
          dmin = d;
          closest_color[Ni] = j;
        }

        Ei.setNodes(Ni, Nj);

        if (!hset.insert(Ei).second)
        {
          if (closest_color[Ni] != j)
            continue;
          if (closest_color[Nj] != i)
            continue; 

          palmares.push_back(std::make_pair(d, std::make_pair(Ni, Nj)));
        }
      }
    }
  }

#ifdef E_TIME
  std::cout << "pairs : getting the pairs " << palmares.size() << " in "<< cc.elapsed() << std::endl;
  cc.start();
#endif

  std::sort(palmares.begin(), palmares.end());

  E_Int count = 0;
  for (size_t i = 0; i < palmares.size(); ++i)
  {
    const std::pair<E_Int, E_Int>& p = palmares[i].second;
    Ni = p.first;
    Nj = p.second;

    if ((pairs[Ni] != IDX_NONE) && (pairs[Ni] >= 0))
      continue;
    if ((pairs[Nj] != IDX_NONE) && (pairs[Nj] >= 0))
      continue;

    pairs[Nj] = Ni;
    pairs[Ni] = Nj;
    count++;
  }

  // Final cleaning
  for (E_Int i = 0; i < nb_contours; ++i)
    delete trees[i];

#ifdef WIN32
#ifdef E_DEBUG
  K_FLD::IntArray ccc;
  for (size_t i = 0; i < pairs.size(); ++i)
  {
    E_Int E[2];
    if (pairs[i] != IDX_NONE)
    {
      E[0] = i;
      E[1] = pairs[i];
      ccc.pushBack(E, E+2);
    }
  }
  meshIO::write("pairs.mesh", posB, ccc);
#endif
#endif
}

///
void
NodeAssociator::__removeOverlaps
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, std::vector<E_Int>& nmates,
 std::vector<K_FLD::IntArray>& componentsOut)
{
  // Get the overlapping pairs
  
#ifdef DEBUG_NODEASSOCIATOR
  {
    K_FLD::IntArray conB;
    for (size_t c = 0; c < components.size(); ++c)
    {
      NUGA::MeshTool::getBoundary(*components[c], conB);
    
      std::ostringstream o;
      o << "compB_" << c << ".mesh";
      K_CONVERTER::DynArrayIO::write(o.str().c_str(), pos, conB, "BAR");
    }
  }
#endif

  //std::vector<std::vector<K_FLD::IntArray> > parts(components.size());
  //NUGA::non_oriented_edge_set_type dum;
  //for (E_Int i = 0; i < components.size(); ++i)
    //ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(pos, *components[i], dum, parts[i]);

  std::vector<XPair> pairs;
  Intersector::getXPairs(pos, components, pairs);
  
  // Make the priorization
  //fixme
  std::vector<std::vector<E_Int> > componentsKeep, componentsRemove, componentsUnchanged;
  __priorize_X_zones(pos, components, pairs, componentsKeep, componentsRemove, componentsUnchanged);

  // Update mates for overlapping parts (-1)
  K_FLD::IntArray connectB0, connectBk, connectBr, connectBu, Sk, Sr, Su;
  std::vector<E_Int> nodes0, nodesk, nodesr, nodes, nodesu;
  for (size_t c = 0; c < components.size(); ++c)
  {
    nodes.clear();
    nodes0.clear();
    nodesk.clear();
    nodesr.clear();
    nodesu.clear();
    Sk.clear();
    Sr.clear();
    Su.clear();
    
    for (size_t i = 0; i < componentsKeep[c].size(); ++i)
      Sk.pushBack(components[c]->col(componentsKeep[c][i]), components[c]->col(componentsKeep[c][i])+3);
    for (size_t i = 0; i < componentsRemove[c].size(); ++i)
      Sr.pushBack(components[c]->col(componentsRemove[c][i]), components[c]->col(componentsRemove[c][i])+3);
    for (size_t i = 0; i < componentsUnchanged[c].size(); ++i)
      Su.pushBack(components[c]->col(componentsUnchanged[c][i]), components[c]->col(componentsUnchanged[c][i])+3);

    NUGA::MeshTool::getBoundary(*components[c], connectB0);
    connectB0.uniqueVals(nodes0);
    
    NUGA::MeshTool::getBoundary(Sk, connectBk);
    connectBk.uniqueVals(nodesk);
    NUGA::MeshTool::getBoundary(Sr, connectBr);
    connectBr.uniqueVals(nodesr);
    NUGA::MeshTool::getBoundary(Su, connectBu);
    connectBu.uniqueVals(nodesu);
    
#ifdef DEBUG_NODEASSOCIATOR
    std::ostringstream o;
    o << "Sk_" << c << ".mesh";
    K_CONVERTER::DynArrayIO::write(o.str().c_str(), pos, Sk, "BAR");
    o.str("");
    o << "Sr_" << c << ".mesh";
    K_CONVERTER::DynArrayIO::write(o.str().c_str(), pos, Sr, "BAR");
    o.str("");
    o << "Su_" << c << ".mesh";
    K_CONVERTER::DynArrayIO::write(o.str().c_str(), pos, Su, "BAR");
#endif
    
    std::sort(nodes0.begin(), nodes0.end());
    std::sort(nodesk.begin(), nodesk.end());
    std::sort(nodesr.begin(), nodesr.end());
    nodes = nodesr;
    std::set_intersection(nodes0.begin(), nodes0.end(), nodesk.begin(), nodesk.end(), std::back_inserter(nodes));
    //std::set_difference(nodesr.begin(), nodesr.end(), nodes0.begin(), nodes0.end(), std::back_inserter(nodes)); 

    // update overlapping nodes.
    for (size_t n = 0; n < nodes.size(); ++n)
      nmates[nodes[n]] = Zipper::OVERLAP;
  }

  // Update the components for exit.
  componentsOut.resize(components.size());
  std::set<E_Int> rem;
  std::vector<E_Int> newId;
  
  for (size_t c = 0; c < components.size(); ++c)
  {
    rem.clear();
    newId.clear();
    rem.insert(componentsRemove[c].begin(), componentsRemove[c].end());
    newId.resize(components[c]->cols());
    for (size_t k = 0; k < newId.size(); ++k)
      newId[k] = k;
    componentsOut[c].clear();
    for (E_Int i = 0; i < components[c]->cols(); ++i)
    {
      if (rem.find(i) == rem.end())
      {
        componentsOut[c].pushBack(components[c]->col(i), components[c]->col(i)+3);
        newId[i] = componentsOut[c].cols()-1;
      }
    }
    
    // update pairs ids accordingly.
    E_Int cc = (E_Int)c;
    for (size_t j = 0; j < pairs.size(); ++j)
    {
      if (pairs[j].S0 == cc)
        pairs[j].C0 = newId[pairs[j].C0];
      if (pairs[j].S1 == cc)
        pairs[j].C1 = newId[pairs[j].C1];
    }
  }

  // Reorient
  if (reorient)
    __reorientComponents(pos, componentsOut, pairs);
}

void
NodeAssociator::__priorize_X_zones
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, const std::vector<XPair>& pairs, 
 std::vector<std::vector<E_Int> >& componentsKeep, std::vector<std::vector<E_Int> >& componentsRemove,
 std::vector<std::vector<E_Int> >& componentsUnchanged)
{
  std::map<K_MESH::NO_Edge, std::pair<K_FLD::IntArray, K_FLD::IntArray> > overlaps;
  std::map<K_MESH::NO_Edge, std::pair<K_FLD::IntArray, K_FLD::IntArray> >::const_iterator it;
  std::map<K_MESH::NO_Edge, K_MESH::Edge> prior;
  std::map<K_MESH::NO_Edge, K_MESH::Edge>::const_iterator it1;
  K_FLD::IntArray::const_iterator pS1, pS2;

  K_FLD::IntArray iCs1,oCs1;
  componentsKeep.resize(components.size());
  componentsRemove.resize(components.size());
  componentsUnchanged.resize(components.size());
    
  for (size_t p = 0; p < pairs.size(); ++p)
  {
    E_Int c1 = pairs[p].C0;
    E_Int s1 = pairs[p].S0;
    E_Int c2 = pairs[p].C1;
    E_Int s2 = pairs[p].S1;

    pS1 = components[s1]->col(c1);
    pS2 = components[s2]->col(c2);
    
    if (s1 < s2)
    {
      K_MESH::NO_Edge e(s1, s2);
      overlaps[e].first.pushBack(pS1, pS1+3);
      overlaps[e].second.pushBack(pS2, pS2+3);
    }
    else
    {
      K_MESH::NO_Edge e(s2, s1);
      overlaps[e].first.pushBack(pS2, pS2+3);
      overlaps[e].second.pushBack(pS1, pS1+3);
    }
  }

  std::vector<E_Int> nodes1, nodes2;
  for (it = overlaps.begin(); it != overlaps.end(); ++it)
  {
    const K_MESH::NO_Edge& e = it->first;
    
    it->second.first.uniqueVals(nodes1);
    it->second.second.uniqueVals(nodes2);
    prior[e] = K_MESH::Edge(nodes1.size(), nodes2.size());
  }


  std::vector <std::vector<E_Int> > remains(components.size());
  std::vector <std::vector<bool> > xing(components.size());

  for (size_t s = 0; s < components.size(); ++s)
  {
    remains[s].resize(components[s]->cols(), -1);
    xing[s].resize(components[s]->cols(), false);
  }

  for (size_t p = 0; p < pairs.size(); ++p)
  {
    E_Int c1 = pairs[p].C0;
    E_Int s1 = pairs[p].S0;
    E_Int c2 = pairs[p].C1;
    E_Int s2 = pairs[p].S1;

    pS1 = components[s1]->col(c1);
    pS2 = components[s2]->col(c2);

    xing[s1][c1] = true;
    xing[s2][c2] = true;

    K_MESH::NO_Edge e(s1, s2);
    it1 = prior.find(e);
    E_Int ns1 = it1->second.node(0);
    E_Int ns2 = it1->second.node(1);
    E_Int s;
    if (ns1 > ns2)
      s = it1->first.node(1);
    else
      s = it1->first.node(0);

    if (s == s1)
    {
      remains[s1][c1] = 0;
      if (remains[s2][c2] != 0)
        remains[s2][c2] = 1;
    }
    else
    {
      remains[s2][c2] = 0;
      if (remains[s1][c1] != 0)
        remains[s1][c1] = 1;
    }
  }

  std::set<E_Int> hN, hNKeep;
  std::vector<E_Int> colors;
  for (size_t s = 0; s < components.size(); ++s)
  {
    for (E_Int c = 0; c < components[s]->cols(); ++c)
    {
      if (remains[s][c] == 1)
      {
        componentsKeep[s].push_back(c);
        hNKeep.insert((*components[s])(0,c));
        hNKeep.insert((*components[s])(1,c));
        hNKeep.insert((*components[s])(2,c));
      }
      else if (remains[s][c] == 0)
      {
        componentsRemove[s].push_back(c);
        hN.insert((*components[s])(0,c));
        hN.insert((*components[s])(1,c));
        hN.insert((*components[s])(2,c));
      }
      else if (remains[s][c] == -1)
      {
        componentsUnchanged[s].push_back(c);
        hN.insert((*components[s])(0,c));
        hN.insert((*components[s])(1,c));
        hN.insert((*components[s])(2,c));
      }

      if (xing[s][c])
      {
        iCs1.pushBack(components[s]->col(c), components[s]->col(c)+3);
        colors.push_back(s);
      }
      else
        oCs1.pushBack(components[s]->col(c), components[s]->col(c)+3);
    }
  }

  //_hard_nodes.clear();
  //std::set_difference(hN.begin(), hN.end(), hNKeep.begin(), hNKeep.end(), std::back_inserter(_hard_nodes));
#ifdef WIN32
#ifdef E_DEBUG
  meshIO::write("intersecting.mesh", pos, iCs1, colors);
  meshIO::write("out_of_x.mesh", pos, oCs1);
#endif
#endif
}

void
NodeAssociator::__priorize_X_zones2
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, const std::vector<XPair>& pairs, 
 std::vector<std::vector<E_Int> >& componentsKeep, std::vector<std::vector<E_Int> >& componentsRemove)
{
  std::map<K_MESH::NO_Edge, std::pair<K_FLD::IntArray, K_FLD::IntArray> > overlaps;
  std::map<K_MESH::NO_Edge, std::pair<K_FLD::IntArray, K_FLD::IntArray> >::const_iterator it;
  
  K_FLD::IntArray::const_iterator pS1, pS2;

  K_FLD::IntArray iCs1,oCs1;
  componentsKeep.resize(components.size());
  componentsRemove.resize(components.size());
    
  std::vector <std::vector<E_Int> > remains(components.size());
  std::vector <std::vector<bool> > xing(components.size());

  for (size_t s = 0; s < components.size(); ++s)
  {
    remains[s].resize(components[s]->cols(), -1);
    xing[s].resize(components[s]->cols(), false);
  }

  for (size_t p = 0; p < pairs.size(); ++p)
  {
    E_Int c1 = pairs[p].C0;
    E_Int s1 = pairs[p].S0;
    E_Int c2 = pairs[p].C1;
    E_Int s2 = pairs[p].S1;

    pS1 = components[s1]->col(c1);
    pS2 = components[s2]->col(c2);

    xing[s1][c1] = true;
    xing[s2][c2] = true;

    E_Float surf1 = K_MESH::Triangle::surface(pos.col(*pS1), pos.col(*(pS1+1)), pos.col(*(pS1+2)), 3);
    E_Float surf2 = K_MESH::Triangle::surface(pos.col(*pS2), pos.col(*(pS2+1)), pos.col(*(pS2+2)), 3);

    remains[s1][c1] = (surf1 < surf2);
    remains[s2][c2] = !remains[s1][c1];
  }

  std::set<E_Int> hN, hNKeep;
  std::vector<E_Int> colors;
  for (size_t s = 0; s < components.size(); ++s)
  {
    for (E_Int c = 0; c < components[s]->cols(); ++c)
    {
      if (remains[s][c] == 1)
      {
        componentsKeep[s].push_back(c);
        hNKeep.insert((*components[s])(0,c));
        hNKeep.insert((*components[s])(1,c));
        hNKeep.insert((*components[s])(2,c));
      }
      else if (remains[s][c] == 0)
      {
        componentsRemove[s].push_back(c);
        hN.insert((*components[s])(0,c));
        hN.insert((*components[s])(1,c));
        hN.insert((*components[s])(2,c));
      }

      if (xing[s][c])
      {
        iCs1.pushBack(components[s]->col(c), components[s]->col(c)+3);
        colors.push_back(s);
      }
      else
        oCs1.pushBack(components[s]->col(c), components[s]->col(c)+3);
    }
  }

  //_hard_nodes.clear();
  //std::set_difference(hN.begin(), hN.end(), hNKeep.begin(), hNKeep.end(), std::back_inserter(_hard_nodes));
#ifdef WIN32
#ifdef E_DEBUG
  meshIO::write("intersecting.mesh", pos, iCs1, colors);
  meshIO::write("out_of_x.mesh", pos, oCs1);
#endif
#endif
}

void gather_group(E_Int g, const std::map<E_Int, std::set<E_Int> >& groups, std::vector<bool>& inspected, std::set<E_Int>& gp)
{
  std::map<E_Int, std::set<E_Int> >::const_iterator itG = groups.find(g);
  if (itG == groups.end())
    return; // Error.
  std::set<E_Int>::const_iterator it = itG->second.begin();

  for (; it != itG->second.end(); ++it)
  {
    if (inspected[*it])
      continue;
    gp.insert(*it);
    inspected[*it] = true;
    gather_group(*it, groups, inspected, gp);
  }
}

///
void
NodeAssociator::__getRelativeOrient
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& cS, const std::vector<E_Int>& nmates,
 std::map<K_MESH::NO_Edge, E_Int>& surfpair_to_orient)
{
  std::vector<K_FLD::IntArray> cB;
  E_Int nb_surf;
  std::vector<E_Int> ncolor(pos.cols(), IDX_NONE);
  std::vector<std::vector<E_Int> > sorted_nodes;
  std::vector<NUGA::oriented_edge_set_type> hE2s;
  std::map<K_MESH::NO_Edge, E_Int>::const_iterator it;

  surfpair_to_orient.clear();
  
  nb_surf = cS.size();

  cB.resize(nb_surf);
  sorted_nodes.resize(nb_surf);
  hE2s.resize(nb_surf);
  
  // Get the associated boundaries
  std::vector<E_Int> nodes;
  for (E_Int s = 0; s < nb_surf; ++s)
  {
    NUGA::MeshTool::getBoundary(cS[s], cB[s]);
    BARSplitter::getSortedNodes(cB[s], sorted_nodes[s]);
    for (E_Int i = 0; i < cB[s].cols(); ++i)
      hE2s[s].insert(K_MESH::Edge(cB[s](0,i), cB[s](1,i)));
    for (size_t i = 0; i < sorted_nodes[s].size(); ++i)
      ncolor[sorted_nodes[s][i]] = s;
  }

  // Set the relative orientations
  E_Int Ni, Nprev, Nnext, next, prev, Nj, Nio, Njo, sio, sjo;
  K_MESH::Edge E;
  K_MESH::NO_Edge cont_pair;
  K_FLD::IntArray odd;
  for (E_Int s = 0; s < nb_surf; ++s)
  {
    E_Int nb_nodes = sorted_nodes[s].size();
    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      Ni = sorted_nodes[s][n];
      Nj = IDX_NONE;
      Nio = nmates[Ni];
      if ((Nio == IDX_NONE) || (Nio <0))
        continue;
      next = (n+1)%nb_nodes;
      prev = (n != 0) ? n-1 : nb_nodes-1;
      Nnext = sorted_nodes[s][next];
      Nprev = sorted_nodes[s][prev];
      if ((nmates[Nnext] >= 0)&&(nmates[Nnext] != IDX_NONE))
        Nj = Nnext;
      else if ((nmates[Nprev] >= 0)&&(nmates[Nprev] != IDX_NONE))
        Nj = Nprev;
      if (Nj == IDX_NONE)
        continue;
      Njo = nmates[Nj];
      if ((Njo == IDX_NONE) || (Njo <0))
        continue;

      sio = ncolor[Nio]; 
      sjo = ncolor[Njo]; 

      if (sio != sjo)
        continue;

      if (sio == s)
        continue;

      cont_pair.setNodes(s, sio);
      it = surfpair_to_orient.find(cont_pair);
      if (it != surfpair_to_orient.end())
        continue;

      if (hE2s[s].find(K_MESH::Edge(Ni, Nj)) == hE2s[s].end())
        E.setNodes(Nio, Njo);
      else
        E.setNodes(Njo, Nio);

      surfpair_to_orient[cont_pair] = hE2s[sio].find(E) != hE2s[sio].end();
      /*
      if (surfpair_to_orient[cont_pair] == 0)
      {
        E_Int Ed[] = {Ni, Nj};
        if (E.node(0) == Nio)
          std::swap(Ed[0], Ed[1]);
        
        odd.pushBack(Ed, Ed+2);
        Ed[0] = E.node(0);
        Ed[1] = E.node(1);
        odd.pushBack(Ed, Ed+2);
      }
      */
    }
  }
}

///
void
NodeAssociator::__reorientComponent
(const K_FLD::FloatArray& pos, K_FLD::IntArray& component, const std::vector<E_Int>& nmates)
{

  // Split the surfaces.
  NUGA::non_oriented_edge_set_type dummy;
  std::vector<K_FLD::IntArray> cS;
  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(component, dummy, cS);
  E_Int nb_surf = cS.size();
  
  // Get the relative orientations.
  std::map<K_MESH::NO_Edge, E_Int> surfpair_to_orient;
  __getRelativeOrient(pos, cS, nmates, surfpair_to_orient);

  std::map<K_MESH::NO_Edge, E_Int>::iterator it;
  
  std::vector<bool> flipGrp(nb_surf, false);

  std::vector<std::set<E_Int> > attached(nb_surf);
  for (std::map<K_MESH::NO_Edge, E_Int>::const_iterator it = surfpair_to_orient.begin(); it != surfpair_to_orient.end(); ++it)
  {
    attached[it->first.node(0)].insert(it->first.node(1));
    attached[it->first.node(1)].insert(it->first.node(0));
  }

  bool need_to_reorient;
  do 
  {
    need_to_reorient = false;
    E_Int surf = -1;
    for ( it = surfpair_to_orient.begin(); (it != surfpair_to_orient.end());++it)
    {
      if (it->second == 0)
      {
        need_to_reorient = true;
        if (attached[it->first.node(0)].size() < attached[it->first.node(1)].size())
          surf = it->first.node(0);
        else
          surf = it->first.node(1);
        break;
      }
    }

    if (!need_to_reorient) // Orientation is consistent.
      break;

    if (surf == -1)
      break; // error fixme

    flipGrp[surf] = !flipGrp[surf];
    
    // propagate the change
    for ( it = surfpair_to_orient.begin(); (it != surfpair_to_orient.end());++it)
    {
      if ((it->first.node(0) == surf) || (it->first.node(1) == surf))
        it->second = !it->second;
      //std::cout << "rel orient : " << it->first.node(0) << " " << it->first.node(1) << " : " << it->second << std::endl;
    }
    //std::cout << "###############################################" << std::endl << std::endl;
  }
  while (need_to_reorient);

  E_Int count = 0;
  for (E_Int s = 0; s < nb_surf; ++s)
  {
    //std::cout << " need a flip ? " << s << " : " << flipGrp[s] << std::endl;
    if (flipGrp[s])
      ++count;
    else
      --count;
  }
  if (::abs(count) == nb_surf)
    return;

  bool to_reverse = (count < 0) ? true : false;
  for (E_Int s = 0; s < nb_surf; ++s)
  {
    if (flipGrp[s] == to_reverse)
      NUGA::MeshTool::flipT3(cS[s]);
  }

  component.clear();
  for (E_Int s = 0; s < nb_surf; ++s)
    component.pushBack(cS[s]);

}

void
NodeAssociator::__getRelativeOrient
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& comps, const std::vector<XPair>& pairs,
 std::map<K_MESH::NO_Edge, E_Int>& comppair_to_orient)
{
  K_MESH::NO_Edge cont_pair;
  K_FLD::IntArray::const_iterator pS1, pS2;
  std::map<K_MESH::NO_Edge, E_Int>::const_iterator it;
  E_Float Norm1[3], Norm2[3];
  E_Int c0, c1, s0, s1;
  for (size_t xp = 0; xp <  pairs.size(); ++xp)
  {
    cont_pair.setNodes(pairs[xp].S0, pairs[xp].S1);
    it = comppair_to_orient.find(cont_pair);
    if (it != comppair_to_orient.end())
      continue;

    s0 = pairs[xp].S0;
    c0 = pairs[xp].C0;

    if (c0 >= comps[s0].cols()) //obsolete pair
      continue;

    s1 = pairs[xp].S1;
    c1 = pairs[xp].C1;

    if (c1 >= comps[s1].cols()) //obsolete pair
      continue;

    pS1 = comps[pairs[xp].S0].col(pairs[xp].C0);
    pS2 = comps[pairs[xp].S1].col(pairs[xp].C1);
    K_MESH::Triangle::normal(pos.col(*pS1), pos.col(*(pS1+1)), pos.col(*(pS1+2)), Norm1);
    K_MESH::Triangle::normal(pos.col(*pS2), pos.col(*(pS2+1)), pos.col(*(pS2+2)), Norm2);

    comppair_to_orient[cont_pair] = (NUGA::dot<3>(Norm1, Norm2) > 0.);
  }
}

void
NodeAssociator::__reorientComponents
(const K_FLD::FloatArray& pos, std::vector<K_FLD::IntArray>& components, const std::vector<XPair>& pairs)
{
  
  E_Int nb_comps = components.size();

  //std::cout << " in reorient comps" << std::endl;
  
  // Get the relative orientations.
  std::map<K_MESH::NO_Edge, E_Int> comppair_to_orient;
  __getRelativeOrient(pos, components, pairs, comppair_to_orient);

  //std::cout << " relative orient done" << std::endl;

  std::map<K_MESH::NO_Edge, E_Int>::iterator it;
  
  std::vector<bool> flipGrp(nb_comps, false);

  std::vector<std::set<E_Int> > attached(components.size());
  for (std::map<K_MESH::NO_Edge, E_Int>::const_iterator it = comppair_to_orient.begin(); it != comppair_to_orient.end(); ++it)
  {
    attached[it->first.node(0)].insert(it->first.node(1));
    attached[it->first.node(1)].insert(it->first.node(0));
  }

  std::set<E_Int> tried;

  bool need_to_reorient, found1, found2;
  do 
  {
    need_to_reorient = false;
    E_Int comp(-1), comp1, comp2;
    for ( it = comppair_to_orient.begin(); (it != comppair_to_orient.end());++it)
    {
      if (it->second == 0)
      {
        need_to_reorient = true;
        comp1 = it->first.node(0);
        comp2 = it->first.node(1);
        found1 = tried.find(comp1) != tried.end();
        found2 = tried.find(comp2) != tried.end();

        if (attached[comp1].size() < attached[comp2].size())
        {
          if (!found1 || found2)
            comp = comp1;
          else
            comp = comp2;
        }
        else
        {
          if (!found2 || found1)
            comp = comp2;
          else
            comp = comp1;
        }
        break;
      }
    }

    if (!need_to_reorient) // Orientation is consistent.
      break;
    if (comp == -1)
      break;// error fixme

    flipGrp[comp] = !flipGrp[comp];
    tried.insert(comp);
    
    // propagate the change
    for ( it = comppair_to_orient.begin(); (it != comppair_to_orient.end());++it)
    {
      if ((it->first.node(0) == comp) || (it->first.node(1) == comp))
        it->second = !it->second;
      //std::cout << "rel orient : " << it->first.node(0) << " " << it->first.node(1) << " : " << it->second << std::endl;
    }
    //std::cout << "###############################################" << std::endl << std::endl;
    //return;
  }
  while (need_to_reorient);

  //std::cout << " reorient done" << std::endl;

  E_Int count = 0;
  for (E_Int s = 0; s < nb_comps; ++s)
  {
    //std::cout << " need a flip ? " << s << " : " << flipGrp[s] << std::endl;
    if (flipGrp[s])
      ++count;
    else
      --count;
  }
  if (::abs(count) == nb_comps)
    return;

  bool to_reverse = (count < 0) ? true : false;
  for (E_Int s = 0; s < nb_comps; ++s)
  {
    if (flipGrp[s] == to_reverse)
      NUGA::MeshTool::flipT3(components[s]);
  }

  //std::cout << " flip done" << std::endl;
}

