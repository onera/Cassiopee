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

#ifndef __TRI_OVERLAPFIX_CXX__
#define __TRI_OVERLAPFIX_CXX__

#include "Nuga/include/TRI_Conformizer.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/BbTree.h"
#include "Nuga/include/BARSplitter.h"

using namespace K_SEARCH;

#ifdef FLAG_STEP
#include "Nuga/include/chrono.h"
#endif

namespace NUGA
{

///////////////////////////////////////////////////////////////////////////////
class Contour
{
public:
  Contour(const std::vector<E_Int>& cnodes){_nodes.insert(cnodes.begin(), cnodes.end());}
  Contour(const Contour& c){_nodes = c._nodes;}

  bool operator<(const Contour& c) const
  {
    if (_nodes.size() != c._nodes.size())
      return (_nodes.size() < c._nodes.size());
    std::set<E_Int>::const_iterator it(_nodes.begin()), itc(c._nodes.begin()), itEnd(_nodes.end());
    while ((it != itEnd) && (*it == *itc)){++it;++itc;}
    return ((it != itEnd) && (*it < *itc));
  }

private:
  std::set<E_Int> _nodes;
};
///////////////////////////////////////////////////////////////////////////////

///
template <short DIM>
void
TRI_Conformizer<DIM>::__run_correction_beta
(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
 NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc,
 E_Float tolerance)
{
  K_FLD::FloatArray                   normals;
  K_FLD::IntArray::const_iterator     pS;
  std::vector<E_Int>                  anodes, inodes, dupIds;
  std::vector<std::vector<E_Int> >    zones, contour_to_zones;
  algo_type::BoundToEltType           E_to_T;
  std::set<K_MESH::NO_Edge>           common_edges;
  std::vector<bool>                   patho_elts, free_elts;

  // Get the 'edge-to-triangles' map.
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  NUGA::EltAlgo<K_MESH::Triangle>::getBoundToElements(actv, E_to_T);

  //
  NUGA::MeshTool::detectDuplicated(connect, dupIds, false /*strict orient*/);

  // Common edges.
  __get_common_edges(pos, connect, dupIds, E_to_T, common_edges, patho_elts, free_elts);

  // Zones (delimited by common edges)
  // warning : zones must be up to date ans consistent with connect.
  __get_zones(pos, connect, E_to_T, common_edges, zones);  

  E_Int                               nb_zones(zones.size()), sz;
  std::vector<BbTree3D*>              treeZs(nb_zones);
  std::vector<std::vector<BBox3D*> >  boxeZs(nb_zones);
  std::vector<K_FLD::IntArray>        connectZs(nb_zones);
  std::vector<std::vector<E_Int> >    xcTemp(nb_zones), ancTemp(nb_zones);

  // Compute node normals.
  NUGA::MeshTool::computeNodeNormals(pos, connect, normals);

  // Create zone connectivities and associated trees.  
  for (E_Int z = 0; z < nb_zones; ++z)
  {
    std::vector<E_Int>& zone = zones[z];
    sz = zone.size();
    boxeZs[z].resize(sz);
    xcTemp[z].resize(sz);
    ancTemp[z].resize(sz);

    for (E_Int i = 0; i < sz; ++i)
    {
      E_Int c = zone[i];
      pS = connect.col(c);
      connectZs[z].pushBack(pS, pS+3);
      boxeZs[z][i] = new BBox3D(pos, pS, 3);
      xcTemp[z][i] = xc[c];
      ancTemp[z][i] = ancestors[c];
    }

    treeZs[z] = new BbTree3D(boxeZs[z]);
  }

  // Sort surfaces supported by the same contour.
  __buildSameSupportSurfacePairs(pos, connectZs, contour_to_zones);

  //
  std::set<E_Int> overlapz;
  E_Int Zi, Zj, nb_ti, nb_tj;;
  //
  for (size_t g = 0; g < contour_to_zones.size(); ++g)
  {
    overlapz.clear();
    std::vector<E_Int>& group = contour_to_zones[g];

    if (group.size() < 2)
      continue;

    //
    for (size_t i = 0; i < group.size(); ++i)
    {
      Zi = group[i];
      nb_ti = connectZs[Zi].cols();
      assert (nb_ti != 0);

      //DynArrayIO::write("zi.mesh", pos, connectZs[Zi]);

      for (size_t j = i+1; j < group.size(); ++j)
      {
        Zj = group[j];
        nb_tj = connectZs[Zj].cols();
        assert (nb_tj != 0);

        //DynArrayIO::write("zj.mesh", pos, connectZs[Zj]);

        if ((nb_tj == 1) && (nb_ti == 1))
          continue;

        if (__areOverlapping(pos, tolerance, normals, connectZs[Zi], *treeZs[Zi], connectZs[Zj], *treeZs[Zj]))
        {
          overlapz.insert(Zi);
          overlapz.insert(Zj);
        }
      }
    }

    if (overlapz.size() < 2)
      continue;

    // Now find the one with the greatest number of nodes and duplicate it.
    E_Int Zt(IDX_NONE);
    E_Float nt[DIM], ni[DIM];
    std::vector<std::vector<E_Int> > nodesZs(nb_zones);
    bool same_orient;
    size_t nb_elts(0);

    for (std::set<E_Int>::const_iterator it = overlapz.begin(); it != overlapz.end(); ++it)
    {
      Zi = *it;
      connectZs[Zi].uniqueVals(nodesZs[Zi]);
      if (nodesZs[Zi].size() > nb_elts)
      {
        Zt = Zi;
        nb_elts = nodesZs[Zi].size();
      }
    }

    overlapz.erase(Zt);
    __computeAveragedNormal(normals, nodesZs[Zt], nt);

    for (std::set<E_Int>::const_iterator it = overlapz.begin(); it != overlapz.end(); ++it)
    {
      Zi = *it;
      __computeAveragedNormal(normals, nodesZs[Zi], ni);
      same_orient = (NUGA::dot<3>(nt, ni) > EPSILON);

      connectZs[Zi] = connectZs[Zt];
      xcTemp[Zi] = xcTemp[Zt];
      ancTemp[Zi].resize(ancTemp[Zt].size(), ancTemp[Zi][0]);

      if (!same_orient)
        NUGA::MeshTool::flipT3(connectZs[Zi]);
    }
  }

  // update connect and xc before exit.
  connect.clear();
  xc.clear();
  ancestors.clear();

  for (E_Int z = 0; z < nb_zones; ++z)
  {
    connect.pushBack(connectZs[z]);
    xc.insert(xc.end(), xcTemp[z].begin(), xcTemp[z].end());
    ancestors.insert(ancestors.end(), ancTemp[z].begin(), ancTemp[z].end());
  }

  // destroy objects.
  for (E_Int z = 0; z < nb_zones; ++z)
  {
    for (size_t i = 0; i < boxeZs[z].size(); ++i)
      delete boxeZs[z][i];
    delete treeZs[z];
  }
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__run_correction_gamma
(const std::set<K_MESH::NO_Edge>&xpairs, NUGA::int_vector_type&colors,
 NUGA::int_vector_type&xr, K_FLD::IntArray& connect,
 NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc, const K_FLD::FloatArray& pos, NUGA::int_vector_type* priority)
{

  if (xpairs.empty())
  {
    if (priority) priority->clear();
    return;
  }
  
  if (priority == 0) return;

  std::vector<bool> keep(connect.cols(), true);
  std::set<E_Int> oids;
  std::vector<E_Int> priorityOut;

  for (std::set<K_MESH::NO_Edge>::const_iterator it = xpairs.begin(); it != xpairs.end(); ++it)
  {
    oids.insert(it->node(0));
    oids.insert(it->node(1));
  }

  std::map< E_Int, std::map<E_Int, std::vector<E_Int> > > T_to_col_to_sortedC;
  std::map<E_Int, K_FLD::IntArray > col_to_connectC;
  std::map<E_Int, K_FLD::IntArray >::const_iterator it1;
  std::map<E_Int, std::vector<E_Int> >::iterator it21, it22;

  E_Int nb_t3s;
  K_FLD::IntArray connectB;
  //
  for (std::set<E_Int>::const_iterator it = oids.begin(); it != oids.end(); ++it)
  {
    const E_Int& oTi = *it;
    nb_t3s = xr[oTi + 1] - xr[oTi];

    col_to_connectC.clear();

    // separate each sub element per color
    for (E_Int i = 0; i < nb_t3s; ++i)
    {
      const E_Int& Ki = xr[oTi] + i;
      col_to_connectC[colors[Ki]].pushBack(connect.col(Ki), connect.col(Ki) + 3);
    }

    // get boundary of each color
    for (it1 = col_to_connectC.begin(); it1 != col_to_connectC.end(); ++it1)
    {
      connectB.clear();
      //const K_FLD::IntArray& conecC = it1->second;
      const E_Int& color = it1->first;
      NUGA::MeshTool::getBoundary(it1->second, connectB);
      BARSplitter::getSortedNodes(connectB, T_to_col_to_sortedC[oTi][color]);
    }
  }

//#define MIN(a,b) ((a<b)? a: b)
//#define MAX(a,b) ((a>b)? a: b)
  
#ifdef DEBUG_TRI_CONFORMIZER
  std::set<E_Int> good_cols;
#endif
  
  // invalidate colors
  std::set<E_Int> bad_cols;
  E_Int c1, c2, *p1, *p2, sz1, sz2;
  for (std::set<K_MESH::NO_Edge>::const_iterator it = xpairs.begin(); it != xpairs.end(); ++it)
  {
    const E_Int& T1 = it->node(0);
    const E_Int& T2 = it->node(1);

    std::map<E_Int, std::vector<E_Int> >& c_to_sC1 = T_to_col_to_sortedC[T1];
    std::map<E_Int, std::vector<E_Int> >& c_to_sC2 = T_to_col_to_sortedC[T2];

    for (it21 = c_to_sC1.begin(); it21 != c_to_sC1.end(); ++it21)
    {
      c1 = it21->first;
      p1 = &(it21->second[0]);
      sz1 = it21->second.size();

      for (it22 = c_to_sC2.begin(); it22 != c_to_sC2.end(); ++it22)
      {
        c2 = it22->first;
        p2 = &(it22->second[0]);
        sz2 = it22->second.size();

        if (sz2 != sz1) continue;

        if (K_CONNECT::IdTool::equal(p1, p2, sz1, true/*permut accepted*/, false/*strict orient*/))
        {
          // invalidate one subdom
          
          bad_cols.insert( ((*priority)[T2] < (*priority)[T1])? c2 : c1);
          const E_Int& kpT = ((*priority)[T2] < (*priority)[T1])? T1 : T2;
          const E_Int& kpc = ((*priority)[T2] < (*priority)[T1])? c1 : c2;
          const E_Int& removed_priority = ((*priority)[T2] < (*priority)[T1]) ?(*priority)[T2] : (*priority)[T1];
          
          nb_t3s = xr[kpT + 1] - xr[kpT];
          for (E_Int i = 0; i < nb_t3s; ++i)
          {
            const E_Int& Ki = xr[kpT] + i;
            //if (colors[Ki]==kpc)ancestors[Ki]=-ancestors[Ki];  // mark it as having a duplicate
            //if (colors[Ki]==kpc)(*priority)[ancestors[kpT]]=-(*priority)[ancestors[kpT]];
            if (colors[Ki] == kpc && removed_priority > 0)priorityOut.push_back(Ki);
          }
#ifdef DEBUG_TRI_CONFORMIZER
          good_cols.insert(((*priority)[T2] < (*priority)[T1])? c1 : c2);
#endif
        }
      }
    }
  }

  for (E_Int i = 0; i < connect.cols(); ++i)
  {
    keep[i] = (bad_cols.find(colors[i]) == bad_cols.end());
  }


#ifdef DEBUG_TRI_CONFORMIZER
  std::vector<bool> ukeep(keep);
  K_CONNECT::IdTool::negative(ukeep);
  K_FLD::IntArray tmp(connect);
  K_CONNECT::keep<bool> upred(ukeep);
  K_CONNECT::IdTool::compress(tmp, upred);
  medith::write("Oremoved.mesh", pos, tmp, "TRI");
  TRI_debug::write_wired("OremovedN.mesh", pos, tmp, true);
  std::vector<bool> keep2(connect.cols(), false);
  for (size_t i = 0; i < connect.cols(); ++i)
  {
    keep2[i] = (good_cols.find(colors[i]) != good_cols.end());
  }
  K_FLD::IntArray tmp2(connect);
  K_CONNECT::keep<bool> gpred(keep2);
  K_CONNECT::IdTool::compress(tmp2, gpred);
  medith::write("Okept.mesh", pos, tmp2, "TRI");
  TRI_debug::write_wired("OkeptN.mesh", pos, tmp2, true);
#endif
  
  K_CONNECT::keep<bool> pred(keep);
  K_CONNECT::IdTool::compress(connect, pred);
  K_CONNECT::IdTool::compress(ancestors, pred);
  std::vector<E_Int> nids;
  K_CONNECT::IdTool::compress(xc, pred, nids);
 
  *priority = priorityOut;
  
  for (size_t i = 0; i < priority->size(); ++i)
    (*priority)[i] = nids[(*priority)[i]];
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__buildSameSupportSurfacePairs
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray >& connectZs,
 std::vector<std::vector<E_Int> >& grouped_zones)
{
  K_FLD::IntArray connectb, tmp;
  std::vector<E_Int> nodes, dids;
  std::map<Contour, std::vector<E_Int> > contour_to_zones;

  grouped_zones.clear();

  for (size_t z = 0; z < connectZs.size(); ++z)
  {
    tmp = connectZs[z];
 
    NUGA::MeshTool::removeDuplicated(tmp, dids, false);

    NUGA::MeshTool::getBoundary(tmp, connectb);
    connectb.uniqueVals(nodes);
    Contour Ci(nodes);
    contour_to_zones[Ci].push_back(z);
  }

  for (std::map<Contour, std::vector<E_Int> >::iterator it = contour_to_zones.begin(); it != contour_to_zones.end(); ++it)
    grouped_zones.push_back(it->second);
}

///
template <short DIM>
E_Bool
TRI_Conformizer<DIM>::__areOverlapping
(const K_FLD::FloatArray& pos, E_Float tolerance, const K_FLD::FloatArray& normals,
 const K_FLD::IntArray& connectZ1, const K_SEARCH::BbTree3D& treeZ1,
 const K_FLD::IntArray& connectZ2, const K_SEARCH::BbTree3D& treeZ2)
{

  // warning : we assume here that the 2 input surfaces have the same support.

  std::vector<E_Int>  anodes1, bnodes, inodes1;
  std::vector<E_Int>  anodes2, inodes2;
  K_FLD::IntArray     connectb;
  E_Int               Ni;
  E_Bool              isFar;

  connectZ1.uniqueVals(anodes1);
  connectZ2.uniqueVals(anodes2);

  NUGA::MeshTool::getBoundary(connectZ1, connectb);

  connectb.uniqueVals(bnodes);
  std::sort(bnodes.begin(), bnodes.end());

  std::sort(anodes1.begin(), anodes1.end());
  std::set_difference(anodes1.begin(), anodes1.end(), bnodes.begin(), bnodes.end(), std::back_inserter(inodes1));
  std::sort(anodes2.begin(), anodes2.end());
  std::set_difference(anodes2.begin(), anodes2.end(), bnodes.begin(), bnodes.end(), std::back_inserter(inodes2));

  if (inodes1.empty() && inodes2.empty())
    return true;

  for (size_t i = 0; i < inodes1.size(); ++i)
  {
    Ni = inodes1[i];
    isFar = __IsNodeFarFromSurface(pos, normals, Ni, connectZ2, treeZ2, tolerance);
    if (isFar == true)
      return false;
  }
  
  for (size_t i = 0; i < inodes2.size(); ++i)
  {
    Ni = inodes2[i];
    isFar = __IsNodeFarFromSurface(pos, normals, Ni, connectZ1, treeZ1, tolerance);
    if (isFar == true)
      return false;
  }

  return true;
}

///
template <short DIM>
E_Bool
TRI_Conformizer<DIM>::__IsNodeFarFromSurface
(const K_FLD::FloatArray& pos, const K_FLD::FloatArray& normals,
 E_Int N, const K_FLD::IntArray& connectZ, const K_SEARCH::BbTree3D& treeZ, E_Float tolerance)
{
  E_Float                           P[3], dMax(-1.), d, UV[2];
  std::vector<E_Int>                boxes;
  K_FLD::IntArray::const_iterator   pS;
  bool                              inside;

  for (E_Int k = 0; k < DIM; ++k)
    P[k] = pos(k, N) + normals(k, N);

  boxes.clear();
  treeZ.getIntersectingBoxes(pos.col(N), P, boxes, tolerance);

  for (size_t i=0; (i < boxes.size()) && (dMax < tolerance); ++i)
  {
    pS = connectZ.col(boxes[i]);
    d = K_MESH::Triangle::minDistanceToPoint(pos, pS, pos.col(N), UV, inside);
    if (inside)
      dMax = std::max(dMax, d);
  }

  if (dMax < 0.)
    return true;

  return (dMax >= tolerance);
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__computeAveragedNormal
(const K_FLD::FloatArray& normals, const std::vector<E_Int>& nodes, E_Float* ni)
{
  ni[0] = ni[1] = ni[2] = 0.;

  E_Int sz = nodes.size();
  
  if (sz == 0)
    return;

  for (size_t i = 0; i < nodes.size(); ++i)
  {
    for (E_Int k = 0; k < DIM; ++k)
      ni[k] += normals(k, nodes[i]);
  }

  NUGA::normalize<DIM>(ni);
}

///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__get_common_edges
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& dupIds,
 const NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType& bound_to_elts,
 std::set<K_MESH::NO_Edge>& common_edges,
 std::vector<bool>& patho_elts, std::vector<bool>& free_elts)
{
  // Common edges are intersecting edges AND boundary of overlapping zones edges.
  // We assume here that the input triangles are all manifold : i.e : any edge is shared by 2 triangles (can be an open surface so can be 1).
  // After intersection, intersectiong edge are those which are shared by more than 2 triangles.

  NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType::const_iterator it;
  E_Int                           COLS(connect.cols());
  size_t                          Si, nb_elts;
  std::set<E_Int>                 unique_elts;

  common_edges.clear();
  patho_elts.clear();
  free_elts.clear();
  patho_elts.resize(COLS, false);
  free_elts.resize(COLS, false);

  for (it = bound_to_elts.begin(); it != bound_to_elts.end(); ++it)
  {
    const std::vector<E_Int>& elts = it->second; 

    unique_elts.clear();
    nb_elts = elts.size();

    for (Si = 0; Si < nb_elts; ++Si)
      unique_elts.insert(dupIds[elts[Si]]);

    nb_elts = unique_elts.size();

    if (nb_elts == 2)
      continue;
    else if (nb_elts == 1)  // Free.
    {
      nb_elts = elts.size();
      for (Si = 0; Si < nb_elts; ++Si)
        free_elts[elts[Si]] = true;
    }
    else                    // 3 or more : non manifold : considered as common edge.             
    {
      common_edges.insert(it->first);

      if ((nb_elts > 4))    // Pathology : 5 or more.
      {
        nb_elts = elts.size();
        for (Si = 0; Si < nb_elts; ++Si)
          patho_elts[elts[Si]] = true;
      }
    }
  }

#ifdef DEBUG_TRI_CONFORMIZER
  K_FLD::IntArray com;
  for (std::set<K_MESH::NO_Edge>::iterator it = common_edges.begin(); it != common_edges.end(); ++it)
  {
    E_Int E[]={it->node(0), it->node(1)};
    com.pushBack(E, E+2);
  }
  std::ostringstream o;
  o << "common_" << parent_type::_iter << ".mesh";
  medith::write(o.str().c_str(), pos, com, "BAR");
#endif

  return 0;
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__get_zones
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 const algo_type::BoundToEltType& E_to_T,
 const std::set<K_MESH::NO_Edge>& common_edges, std::vector< std::vector<E_Int> >& zones)
{
  //K_FLD::IntArray::const_iterator           pS;
  E_Int                                     c;
  K_MESH::NO_Edge                           Ei;
  K_MESH::Triangle                          Ti;
  algo_type::BoundToEltType::const_iterator itE;
  size_t                                    sz;
  algo_type::NeighbourType                  neighbors;

  zones.clear();

#ifdef FLAG_STEP
  NUGA::chrono tt;
  tt.start();
#endif

  // Get the neighbouring of each element.   
  neighbors.resize(connect.cols());
  E_Int Si, Sj;
  //
  for (itE = E_to_T.begin(); itE != E_to_T.end(); ++itE)
  {
    const std::vector<E_Int>& elts = itE->second;
    sz = elts.size();
    for (size_t i = 0; i < sz; ++i)
    {
      Si = elts[i];     
      for (size_t j = i+1; j < sz; ++j)
      {
        Sj = elts[j];
        neighbors[Si].push_back(Sj);
        neighbors[Sj].push_back(Si);
      }
    }
  }

  if (neighbors.empty())
    return;

  // Color the elements.
  std::vector<E_Int> colors;
  NUGA::EltAlgo<K_MESH::Triangle>::coloring(connect, neighbors, common_edges, colors);

  // Color by zone.
  E_Int maxcol = *std::max_element(colors.begin(), colors.end());
  zones.resize(maxcol+1);

  for (c = 0; c < connect.cols(); ++c)
  {
    //pS = connect.col(c);
    zones[colors[c]].push_back(c);
  }

#ifdef FLAG_STEP
  std::cout << "inter : build zones " << tt.elapsed() << std::endl;
  tt.start();
#endif
}

/*
///
template <short DIM>
void
TRI_Conformizer<DIM>::__run_correction_alpha1
(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
 NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc)
{
  std::set<K_MESH::NO_Edge>::const_iterator pEi;
  K_FLD::IntArray::const_iterator           pS;
  E_Int                                     c;
  K_MESH::NO_Edge                           Ei;
  K_FLD::IntArray                           cOut, sconnect, connectb, odd, uconnect;
  NUGA::bool_vector_type              xcOut, pathological, free;
  NUGA::int_vector_type               ancOut, odd_colors, ids, dupIds;
  bool                                      has_free_elt, has_patho_elt, keep;
  std::vector<std::vector<E_Int> >          zonesOut;
  algo_type::BoundToEltType                 E_to_T;
  std::vector<std::vector<E_Int> >          zones;
  std::set<K_MESH::NO_Edge>                 common_edges;
  
#ifdef FLAG_STEP
  NUGA::chrono tt;
  tt.start();
#endif

  // Get the 'edge-to-triangles' map.
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  NUGA::EltAlgo<K_MESH::Triangle>::getBoundToElements(actv, E_to_T);
  // Detect the duplicated elements.
  DELAUNAY::MeshTool::detectDuplicated(connect, dupIds, false);//strict orient

#ifdef FLAG_STEP
  std::cout << "inter : bound to elts " << tt.elapsed() << std::endl;
  tt.start();
#endif

  // Common edges.
  __get_common_edges(pos, connect, dupIds, E_to_T, common_edges, pathological, free);

#ifdef FLAG_STEP
  std::cout << "inter : common edges " << tt.elapsed() << std::endl;
  tt.start();
#endif

  // Zones (delimited by common edges).
  __get_zones(pos, connect, E_to_T, common_edges, zones);

#ifdef FLAG_STEP
  std::cout << "inter : zones " << tt.elapsed() << std::endl;
  tt.start();
#endif

  std::set<E_Int> tmp;

  // Remove any pathological bit of mesh.
  for (size_t z = 0; z < zones.size(); ++z)
  {
    std::vector<E_Int>& zone = zones[z];
    
    sconnect.clear();
    sconnect.reserve(3, zone.size());
    uconnect.clear();

    tmp.clear();

    has_free_elt = has_patho_elt = false;
    keep = true;

    for (size_t i = 0; (i < zone.size()) && keep; ++i)
    {
      c = zone[i];
      pS = connect.col(c);
      sconnect.pushBack(pS, pS+3);

      tmp.insert(dupIds[c]);

      has_free_elt |= free[c];
      has_patho_elt |= pathological[c];
      keep = !(has_free_elt && has_patho_elt);
    }

    if (keep) // another criterium : parasite closed surface attached to the intersection line.
    {
      if (tmp.size() < sconnect.cols()) // getBoundary works properly on a mesh without duplicated elements.
      {
        for (std::set<E_Int>::const_iterator jj = tmp.begin(); jj != tmp.end(); ++jj)
          uconnect.pushBack(connect.col(*jj), connect.col(*jj)+3);
        DELAUNAY::MeshTool::getBoundary(uconnect, connectb);
      }
      else
        DELAUNAY::MeshTool::getBoundary(sconnect, connectb);

      keep = (connectb.cols() != 0); // not a closed zone.
    }

#ifdef DEBUG_TRI_CONFORMIZER
    //std::ostringstream oo;
    //oo << "bit_" << _iter << "_" << z << ".mesh";
    //DynArrayIO::write(oo.str().c_str(), pos, sconnect);
#endif

    if (keep)
    {
      E_Int N0 = cOut.cols();
      cOut.pushBack(sconnect);
      E_Int N1 = cOut.cols();

      std::vector<E_Int> elts;
      for (E_Int i = N0; i < N1; ++i)
        elts.push_back(i);

      zonesOut.push_back(elts);

      for (size_t i = 0; i < zone.size(); ++i)
      {
        c = zone[i];
        pS = connect.col(c);
        xcOut.push_back(xc[c]);
        ancOut.push_back(ancestors[c]);
      }
    }
#ifdef DEBUG_TRI_CONFORMIZER
    else
    {
      odd.pushBack(sconnect);

      if (has_free_elt && (connectb.cols() != 0))
        odd_colors.resize(odd_colors.size() + sconnect.cols(), 0);
      else if (!has_free_elt && (connectb.cols() == 0))
        odd_colors.resize(odd_colors.size() + sconnect.cols(), 1);
      else
        odd_colors.resize(odd_colors.size() + sconnect.cols(), 2);
    }
#endif
  }

#ifdef FLAG_STEP
  //std::cout << "inter : correction 1 : main loop " << tt.elapsed() << std::endl;
  //tt.start();
#endif

#ifdef DEBUG_TRI_CONFORMIZER
  if (odd.cols() > 0)
  {
    std::ostringstream o;
    o << "odd.mesh";
    medith::write(o.str().c_str(), pos, odd, "TRI", 0, &odd_colors);
  }
#endif
  
  connect = cOut;
  xc = xcOut;
  ancestors = ancOut;
  zones = zonesOut;
}
*/
/*
///
template <short DIM>
void
TRI_Conformizer<DIM>::__compute_splitting_points
(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 const NUGA::bool_vector_type& xc,
 const std::set<K_MESH::NO_Edge>& commone_no_edges,
 E_Float tol, std::map<K_MESH::NO_Edge, std::vector<E_Int> >& edge_to_point)
{
  K_FLD::IntArray::const_iterator   pS;
  E_Int                             Si, nb_tris(connect.cols()), n, N0, Ni, Nj, NB_NODES(3);
  E_Float                           d2, l0, d0, lambda, eps(EPSILON), *Pi, *Pj, Pk[DIM], dMax, tol2(tol*tol);
  K_MESH::NO_Edge                   Ei;
  
  for (Si = 0; Si < nb_tris; ++Si)
  {
    pS = connect.col(Si);
    
    if (!xc[Si]) // Skip if the triangle is not coming from an intersection.
      continue;

    for (n = 0; n < NB_NODES; ++n)
    {
      N0 = *(pS + n);
      Ni = *(pS + (n+1)%NB_NODES);
      Nj = *(pS + (n+2)%NB_NODES);

      Ei.setNodes(Ni, Nj);
   
      if (commone_no_edges.find(Ei) == commone_no_edges.end())
        continue;

      Pi = pos.col(Ni);
      Pj = pos.col(Nj);

      d2 = NUGA::sqrDistance(pos.col(N0), Pi, 3);
      if (d2 < tol2)
        continue;
      d2 = NUGA::sqrDistance(pos.col(N0), Pj, 3);
      if (d2 < tol2)
        continue;
      d2 = NUGA::sqrDistance(Pi, Pj, 3);
      if (d2 < tol2)
        continue;

      d0 = K_MESH::Edge::edgePointMinDistance<DIM>(Pi, Pj, pos.col(N0), lambda);

      if ((lambda < eps) || (lambda > 1. - eps) || (d0 > tol))
        continue;

      l0 = (lambda < 0.5) ? lambda : 1. - lambda;
      dMax = (l0*l0*d2 > tol2) ? tol : EPSILON;

      if (d0 < dMax)
      {
        for (E_Int i = 0; i < DIM; ++i) // Create the new point.
          Pk[i] = (1. - lambda) * Pi[i] + lambda * Pj[i];

        pos.pushBack(Pk, Pk+3);
        Ni = pos.cols() - 1;
        edge_to_point[Ei].push_back(Ni); 

        break;
      }
    }
  }
}
*/
///
/*template <short DIM>
E_Int
TRI_Conformizer<DIM>::__fix_bad_triangles
(K_FLD::FloatArray& coord, E_Float tol, K_FLD::IntArray& connect,
 NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc)
{
  K_FLD::IntArray::const_iterator   pS;
  E_Int                             Si, nb_tris0(connect.cols()), n, N0, N1, N2, Ni, NB_NODES(3);
  E_Float                           d2, l0, d0, lambda, eps(EPSILON), *Pi, *Pj, Pk[DIM], dMax, tol2(tol*tol);
  K_MESH::NO_Edge                   Ei;
  
  typedef std::map<K_MESH::NO_Edge, std::vector<E_Int> > e_to_pts_t;
  e_to_pts_t E_to_Pts;
  e_to_pts_t::iterator itE;
  E_Float tolx = -1.;
  
#ifdef DEBUG_TRI_CONFORMIZER
  K_FLD::IntArray connectBad;
#endif
  
  // 1st pass : detect bad triangles and store splitting points per edges :
  for (Si = 0; Si < nb_tris0; ++Si)
  {
    pS = connect.col(Si);
    
    if (!xc[Si]) // Skip if the triangle is not coming from an intersection.
      continue;
    
    for (n = 0; n < NB_NODES; ++n)
    {
      N0 = *(pS + n);
      N1 = *(pS + (n+1)%NB_NODES);
      N2 = *(pS + (n+2)%NB_NODES);

      Ei.setNodes(N1, N2);

      Pi = coord.col(N1);
      Pj = coord.col(N2);

      d2 = NUGA::sqrDistance(coord.col(N0), Pi, 3);
      if (d2 < tol2)
        continue;
      d2 = NUGA::sqrDistance(coord.col(N0), Pj, 3);
      if (d2 < tol2)
        continue;
      d2 = NUGA::sqrDistance(Pi, Pj, 3);
      if (d2 < tol2)
        continue;

      d0 = K_MESH::Edge::edgePointMinDistance<DIM>(Pi, Pj, coord.col(N0), lambda);

      if ((lambda < eps) || (lambda > 1. - eps) || (d0 > tol-eps))
        continue;
      
      l0 = (lambda < 0.5) ? lambda : 1. - lambda;
      dMax = (l0*l0*d2 > tol2) ? tol : EPSILON;

      //if (d0 < dMax)
      {
        for (E_Int i = 0; i < DIM; ++i) // Create the new point.
          Pk[i] = (1. - lambda) * Pi[i] + lambda * Pj[i];

        coord.pushBack(Pk, Pk+3);
        Ni = coord.cols() - 1;
        
        tolx = std::max(tolx, d0);
        
        itE = E_to_Pts.find(Ei);
        if (itE == E_to_Pts.end())
        {
          E_to_Pts[Ei].push_back(N1);
          E_to_Pts[Ei].push_back(N2);
        }
        E_to_Pts[Ei].push_back(Ni);
        
#ifdef DEBUG_TRI_CONFORMIZER
        connectBad.pushBack(pS, pS+3);
#endif
        //if (N0==10816 || N1 == 10816 || N2 == 10816)
          //std::cout << "caught " << Si << std::endl;
        
        break;
      }
    }
  }
  
  if (E_to_Pts.empty())
    return 0;
  
#ifdef DEBUG_TRI_CONFORMIZER
  medith::write("bad_triangles.mesh", coord, connectBad, "TRI");
  connectBad.clear();
#endif
  
  // sort the points on edges
  for (itE = E_to_Pts.begin(); itE != E_to_Pts.end(); ++itE)
    __tidy_edge(coord, itE->second);
  
  // 2nd pass : do the split
  K_FLD::IntArray new_connect(connect);
  
  E_Int nb_split, nb_tris, nb_new_tris, n1;
  size_t sz;
  
  do
  {
  nb_split = 0;
  nb_new_tris=nb_tris = connect.cols();
  
#ifdef DEBUG_TRI_CONFORMIZER
  //K_FLD::IntArray connectSplit;
  //K_FLD::IntArray::iterator pSi;
#endif
  
  for (Si = 0; Si < nb_tris; ++Si)
  {
    pS = connect.col(Si);
    
    if (!xc[Si]) // Skip if the triangle is not coming from an intersection.
      continue;
    
    for (n = 0; n < NB_NODES; ++n)
    {
      N0 = *(pS + n);
      n1 = (n+1)%NB_NODES;
      N1 = *(pS + n1);
      N2 = *(pS + (n+2)%NB_NODES);
      
      Ei.setNodes(N0, N1);
      
      itE = E_to_Pts.find(Ei);
      if (itE == E_to_Pts.end())
        continue;
      
      nb_split++;
      
#ifdef DEBUG_TRI_CONFORMIZER
      //connectBad.pushBack(pS, pS+3);
      //E_Int nbtris0=nb_new_tris;
#endif
      
      const std::vector<E_Int>& nodes = itE->second;
      sz = nodes.size();
        
      if (nodes[0] == N0) // same orientation
      {
        new_connect(n1,Si)=nodes[1];
        for (size_t i=1; i < sz-1; ++i)
        {
          new_connect.pushBack(pS, pS+3);
          new_connect(n,nb_new_tris) = nodes[i];
          new_connect(n1,nb_new_tris) = nodes[i+1];
          nb_new_tris++;
        }
      }
      else //oposite orientation
      {
        new_connect(n1,Si)=nodes[sz-1];
        for (size_t i=sz-1; i > 0; --i)
        {
          new_connect.pushBack(pS, pS+3);
          new_connect(n,nb_new_tris) = nodes[i];
          new_connect(n1,nb_new_tris) = nodes[i-1];
          nb_new_tris++;
        }
      }
      
#ifdef DEBUG_TRI_CONFORMIZER
      //pSi = new_connect.col(Si);
      //connectSplit.pushBack(pSi, pSi+3);
#endif
      
      xc.resize(nb_new_tris, true);
      ancestors.resize(nb_new_tris, ancestors[Si]);
      
      break;
    }
  }
  
#ifdef DEBUG_TRI_CONFORMIZER
  //medith::write("triangles_to_split.mesh", coord, connectBad, "TRI");
  //for (size_t i=nb_tris; i < new_connect.cols(); ++i)
  //{
  //  pSi = new_connect.col(i);
  //  connectSplit.pushBack(pSi, pSi+3);
  //}
  //medith::write("triangles_splitted.mesh", coord, connectSplit, "TRI");
#endif
  
  connect = new_connect;
  }
  while (nb_split);
  
  //last clean before exit
  std::vector<E_Int> nids;
  if (parent_type::__removeDegenerated(connect, nids))// returns the new ids of connect.
  {
    IdTool::compress_vector(ancestors, nids);
    IdTool::compress_vector(xc, nids);
  }
  
  return (connect.cols() - nb_tris0);
}*/

// UNUSED FUNCIONS  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
#define GOOD_MERGE(E0,E1) ((nids[E0] == nids[E1]) || (commonEdges.find(K_MESH::NO_Edge(nids[E0], nids[E1])) != commonEdges.end()))
///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__swap_edges
(std::vector<std::pair<E_Int, E_Int> >& swapE, K_FLD::IntArray& connect, std::vector<E_Int> ancestors, K_FLD::IntArray& neighbors,
 const std::vector<E_Int>& nids, std::set<K_MESH::NO_Edge>& commonEdges)
{
  E_Int S, b, Sn, bn, Ni, Nj, Nk, Nl, S1, S2, b1, b2, count(1), Sremain;
  K_FLD::IntArray::iterator pS, pSn;
  bool do_swap;
  
  while (count && !swapE.empty())
  {
    count = 0;
    
    E_Int Xnb = swapE.size();
         
   for (size_t i=0; i < Xnb; ++i)
   {
      std::pair<E_Int, E_Int>& E = swapE[i];

      S = E.first;
      b = E.second;
      Sremain = IDX_NONE;

      pS = connect.col(S);
      Ni = *(pS + b);
      Nj = *(pS + (b+1) % 3);
      Nl = *(pS + (b+2) % 3);

      Sn = neighbors(b, S);
      if (Sn == IDX_NONE)
        continue;

      pSn = connect.col(Sn);
      bn = element_type::getOppLocalNodeId(S, b, connect, neighbors);

      Nk = *(pSn + bn);
      
     
      do_swap=!GOOD_MERGE(Nj, Nl); // if the edge won't merge properly
      
      /// DO THE SWAP //////////////////////////////////////////////////////////
      if (do_swap)
      {
      // Neighbors to modify
      S1 = neighbors((b+1) % 3, S);
      S2 = neighbors((bn+1) % 3, Sn);
      b1 = element_type::getOppLocalNodeId(S, (b+1) % 3, connect, neighbors);
      b2 = element_type::getOppLocalNodeId(Sn, (bn+1) % 3, connect, neighbors);

      // Update elements S and Sn (connect)
      *(pS  + (b+2)  % 3) = Nk;
      *(pSn + (bn+2) % 3) = Ni;
    
      // Update the neighboring (neighbors)
      neighbors((b+1)  % 3, S)  = Sn;
      neighbors((bn+1) % 3, Sn) = S;
      if ((S1 != IDX_NONE) && (b1 != IDX_NONE))
        neighbors(b1, S1)              = Sn;
      neighbors(bn, Sn)                = S1;
      if ((S2 != IDX_NONE) && (b2 != IDX_NONE))
        neighbors(b2, S2)              = S;
      neighbors(b, S)                  = S2;
      
      // update the ancestor (keep surviving color after merge
      if (GOOD_MERGE(Nj, Ni) || GOOD_MERGE(Ni, Nl)) // S will be crunched so Sn will survive
        Sremain = Sn;
      else if (GOOD_MERGE(Nj, Nk) || GOOD_MERGE(Nk, Nl)) //versa
        Sremain = S;
      assert(Sremain != IDX_NONE);
      ancestors[S]=ancestors[Sn]=Sremain;
      //////////////////////////////////////////////////////////////////////////
      
      // Update the remaining swap edges
      // (Sn, (bn+1)) -> (S,b)
      // (S, (b+1))   -> (Sn, bn)
      for (E_Int j = 0; j < Xnb; ++j)
      {
        //if (i==j)
        //  continue;
        
        std::pair<E_Int, E_Int>& E = swapE[j];
        if ((E.first == Sn) && (E.second == (bn+1) % 3))
        {
          E.first = S;
          E.second = b;
        }
        else if ((E.first == S) && (E.second == (b+1) % 3))
        {
          E.first = Sn;
          E.second = bn;
        }
      }
      }
      
            
      if (!do_swap || (nids[Nk]==Nk) && (nids[Ni] == Ni)) // ! do_swap means that it will merge properly
      {
        ++count;
        //Remove the edge.
         E = swapE.back(); // compacting : put the last in the current for the next pass
         swapE.pop_back(); // and remove the last.
      }
    }
  }

  return 0;
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__clean_topology
(K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, NUGA::bool_vector_type& xc, 
 const std::vector<E_Int>& source_nodes_vec, std::vector<E_Int>& common_nodes_vec,
 const K_FLD::IntArray& connectBNM, const std::vector<E_Int>& nids
#ifdef DEBUG_TRI_CONFORMIZER
        ,const K_FLD::FloatArray& coord
#endif
)
{
  // Get the working connect : any element of connect attached to a source node
  K_FLD::IntArray connectW;
  std::vector<E_Int> oids; // to apply the changes to connect
  std::vector<E_Int> ancW;
  
  std::set<E_Int> source_nodes_set(source_nodes_vec.begin(), source_nodes_vec.end());
  std::set<E_Int> common_nodes_set(common_nodes_vec.begin(), common_nodes_vec.end());
  
  std::set<K_MESH::NO_Edge> commonEdges;
  for (size_t i=0; i < connectBNM.cols(); ++i) commonEdges.insert(K_MESH::NO_Edge(connectBNM(0,i), connectBNM(1,i)));
  
  for (size_t i=0; i < connect.cols(); ++i)
  {
    for (size_t j=0; j < 3; ++j)
    {
      if (source_nodes_set.find(connect(j,i)) != source_nodes_set.end())
      {
        oids.push_back(i);
        connectW.pushBack(connect.col(i), connect.col(i)+3);
        ancW.push_back(ancestors[i]);
        break;
      }
    }
  }
    
#ifdef DEBUG_TRI_CONFORMIZER
  medith::write("toswap0.mesh", coord, connectW, "TRI");
#endif
    
 #ifdef DEBUG_TRI_CONFORMIZER   
  K_FLD::IntArray ee;
  //std::cout << 709 << "->" << nids[709] << std::endl;
  //std::cout << 710 << "->" << nids[710] << std::endl;
  std::vector<E_Int> colors;
#endif
  
  // Get the edges to swap (Kb format)
  std::vector<std::pair<E_Int, E_Int> > swapKb;
  std::set<K_MESH::NO_Edge> for_unicity;
  K_MESH::NO_Edge E; 
  for (size_t i=0; i < connectW.cols(); ++i)
  {
#ifdef DEBUG_TRI_CONFORMIZER
    E_Int sszz = ee.cols();
#endif
    
    for (size_t j=0; j < 3; ++j)
    {
      const E_Int& Ni = connectW(j,i);
      const E_Int& Nj = connectW((j+1)%3,i);
      const E_Int& nNi = nids[Ni];
      const E_Int& nNj = nids[Nj];
      
      if (Ni == nNi && Nj == nNj) // not a potential faulty edge
        continue;
      
      if (GOOD_MERGE(Ni,Nj)) //will merge properly
        continue;
      
      // potentital faulty have both nodes merging on the contour.
      if (common_nodes_set.find(nNi) == common_nodes_set.end())
        continue;
      if (common_nodes_set.find(nNj) == common_nodes_set.end())
        continue;
      
      //faulty edge : need a swap
      if (for_unicity.insert(K_MESH::NO_Edge(Ni, Nj)).second)//not already in
      {
        swapKb.push_back(std::make_pair(i, (j+2)%3));
#ifdef DEBUG_TRI_CONFORMIZER
        ee.pushBack(E.begin(), E.end());
#endif
      }
    }
#ifdef DEBUG_TRI_CONFORMIZER
    colors.push_back(sszz < ee.cols() ? 1 : 0);
#endif
  }
#ifdef DEBUG_TRI_CONFORMIZER
  medith::write("swapE.mesh", coord, ee, "BAR");
  medith::write("toswapM0.mesh", coord, connectW, "TRI", 0, &colors);
#endif
  
  if (!swapKb.empty())
  {
    K_FLD::IntArray neighbors;
    NUGA::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours (connectW, neighbors, false);//false : accept non manifold
    
    //the ones to target for a swap : connecW nodes - any going to be merged (remove also contou nodes to reduce the set, but they are par of the swapping
    std::vector<E_Int> swapnodes, result; 
    connectW.uniqueVals(swapnodes);
    std::sort(swapnodes.begin(), swapnodes.end());
    std::set_difference(swapnodes.begin(), swapnodes.end(), common_nodes_vec.begin(), common_nodes_vec.end(), std::back_inserter(result));
    
    E_Int err = __swap_edges(swapKb, connectW, ancW, neighbors, nids, commonEdges);
    
#ifdef DEBUG_TRI_CONFORMIZER
    std::vector<E_Int> colz(connectW.cols(), 0);
    for (size_t i=0; i < swapKb.size(); ++i)
    {
      E_Int K = swapKb[i].first;
      E_Int b = swapKb[i].second;
      
      colz[K]=1;
    }   
    medith::write("toswap1.mesh", coord, connectW, "TRI", 0, &colz);
#endif
  
    //Apply swapping to connect : now we should have a coorect topology to prevent faulty elements after merging
    for (size_t i=0; i < connectW.cols(); ++i)
      for (size_t j=0; j < 3; ++j)
        connect(j, oids[i]) = connectW(j,i);
    // update ancesors
    for (size_t i=0; i < ancW.size(); ++i)
      ancestors[oids[i]]=ancW[i];
  }

  //Apply (supposedly safely) the merge
  parent_type::__clean(nids, connect, ancestors, &xc);
  
#ifdef DEBUG_TRI_CONFORMIZER
  {
    std::vector<E_Int> colz(connect.cols(), 0);
    for (size_t i=0; i < swapKb.size(); ++i)
    {
      E_Int K = swapKb[i].first;
      E_Int b = swapKb[i].second;
      
      colz[oids[K]]=1;
    }   
    medith::write("toswap1m.mesh", coord, connect, "TRI", 0, &colz);
  }
#endif
} 

///
template <short DIM>
void
TRI_Conformizer<DIM>::__multi_swap
(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, NUGA::bool_vector_type& xc, std::vector<E_Int>& ancestors)
{
  // swap all the triangles connected to any flat triangle.
    
#ifdef FLAG_STEP
    
    NUGA::chrono tt;
    tt.start();
#endif
    
    std::vector<std::pair<E_Int, E_Int> > swapE;
    __get_swapE(pos, connect, swapE);

#ifdef FLAG_STEP
  if (NUGA::chrono::verbose > 1) std::cout << "Conformizer : correction 1 : finding swapE " << tt.elapsed() << std::endl;
  tt.start();
#endif
  
  if (swapE.empty())
    return;
  
  K_FLD::ArrayAccessor<K_FLD::IntArray> acT3(connect);
  typedef NUGA::EltAlgo<K_MESH::Triangle> algoT3;
  algoT3::BoundToEltType noE_to_oTs; // non oriented edge to oriented triangles (1 to n).
  
  algoT3::getBoundToElements(acT3, noE_to_oTs);
  
#ifdef FLAG_STEP
  if (NUGA::chrono::verbose > 1) std::cout << "Conformizer : correction 1 : getBoundToElements " << tt.elapsed() << std::endl;
  tt.start();
#endif
  
  std::vector<std::vector<E_Int> > T_to_Swap;
  __get_triangles_to_swap(pos, connect, parent_type::_connect0, ancestors, swapE, T_to_Swap);
  
#ifdef FLAG_STEP
  if (NUGA::chrono::verbose > 1) std::cout << "Conformizer : correction 1 : __get_triangles_to_swap " << tt.elapsed() << std::endl;
  tt.start();
#endif
  
  algoT3::BoundToEltType::iterator it;
  K_MESH::NO_Edge noE;
  Vector_t<bool> keep(connect.cols(), true);
  K_FLD::IntArray splitT3s;
  E_Int *pS, *pSk, T[3], j;
  Vector_t<E_Int> new_xc;
#ifdef DEBUG_CONFORMIZER
  Vector_t<E_Int> tmp;
#endif

  E_Int swsz=swapE.size();
  for (size_t i=0; i < swsz; ++i)
  {
    pS = connect.col(swapE[i].first);
    j = swapE[i].second;
    noE.setNodes(*(pS+(j+1)%3), *(pS+(j+2)%3));
    it =  noE_to_oTs.find(noE);
    assert (it != noE_to_oTs.end());
    
    Vector_t<E_Int> & T3s = it->second;
    
#ifdef DEBUG_CONFORMIZER
    tmp.insert(tmp.end(), T3s.begin(), T3s.end());
#endif
    
    size_t sz = it->second.size();
                                    //    x T[2]        x
                                    // ------    =>   / | \
                                    //  \  /          \ | /
    for (size_t k=0; k <sz; ++k)    //   \/            \|/
    {                               //   T[0]
      keep[T3s[k]] = false;
      
      if (swapE[i].first == T3s[k])
        continue;
      
      pSk = connect.col(T3s[k]);
      E_Int n = IDX_NONE;
      
      //get the edge rank
      for (size_t u=0; u<3; ++u)
      {
        if ( (*(pSk+u) == noE.node(0) && *(pSk+(u+1)%3) == noE.node(1)) || (*(pSk+u) == noE.node(1) && *(pSk+(u+1)%3) == noE.node(0)) )
        {n=u;break;}  
      }
      assert (n != IDX_NONE);
      
      //std::cout << "inital T and n : " << *(pSk) << "," << *(pSk+1) << "," << *(pSk+2) << " for n : " << n << std::endl;
	  
      T[0] = *(pSk+(n+2)%3);
      T[1]=*(pSk+n);
	  T[2] = *(pS + j);
	  splitT3s.pushBack(T, T+3);
      
      //std::cout << "T1 : " << T[0] << "," << T[1] << "," << T[2] << std::endl;
      
      T[1]=T[2];
      T[2]=*(pSk+(n+1)%3);
      splitT3s.pushBack(T, T+3);
      
      //std::cout << "T2 : " << T[0] << "," << T[1] << "," << T[2] << std::endl;
      
      new_xc.push_back(xc[T3s[k]]);
      new_xc.push_back(xc[T3s[k]]);
    }
  }
  
#ifdef FLAG_STEP
  if (NUGA::chrono::verbose > 1) std::cout << "Conformizer : correction 1 : splitting " << tt.elapsed() << std::endl;
  tt.start();
#endif
  
#ifdef DEBUG_CONFORMIZER
//  K_FLD::IntArray ctmp;
//  for (size_t i=0; i< tmp.size(); ++i) ctmp.pushBack(connect.col(tmp[i]), connect.col(tmp[i])+3);
//  medith::write("beforesplit.mesh", pos, ctmp, "TRI");
//  medith::write("aftersplit.mesh", pos, splitT3s, "TRI");
#endif
  
  {
    Vector_t<E_Int> nids;
    K_FLD::IntArray::compact(connect, keep, nids);
    K_CONNECT::IdTool::compact(xc, nids);
    connect.pushBack(splitT3s);
    xc.insert(xc.end(), new_xc.begin(), new_xc.end());
  }
  
  ///////////////////////////////////////////////////////////////////////
}
*/


/*
///
template <short DIM>
void
TRI_Conformizer<DIM>::__get_triangles_to_swap
 (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& connect0,
  const std::vector<E_Int>& ancestors, const std::vector<std::pair<E_Int, E_Int> >& swapE,
  std::vector<std::vector<E_Int> >& T_to_Swap)
{
  E_Int Si,j;
  size_t sz = swapE.size(), szt;
  std::vector<E_Int> oT3s;
  K_MESH::NO_Edge noE, noEt;
  const E_Int* pS;
    
  T_to_Swap.clear();
  T_to_Swap.resize(sz);
    
  //Filtrate first the potential triangles : work with ancestors to use the kdtree
  for (size_t e=0; e < sz; ++e)
  {
    Si=swapE[e].first;
    j=swapE[e].second;
    
    pS = connect.col(Si);
    noE.setNodes(*(pS+(j+1)%3), *(pS+(j+2)%3));
    
    BBox3D b(pos, connect0.col(ancestors[Si]), 3);
    
    parent_type::_tree->getOverlappingBoxes(&b.minB[0], &b.maxB[0], oT3s);
  }
  
  std::set<E_Int> uancs;
  uancs.insert(oT3s.begin(), oT3s.end());
  
  for (size_t e=0; e < sz; ++e)
  {
      //to do
      
  }
    
//    szt=T3s.size();
//    for (size_t t=0; t<szt; ++t)
//    {
//      if (T3s[t]==Si)
//        continue;
//      pS=connect.col(T3s[t]);
//      
//      for (size_t n=0; n<3; ++n)
//      {
//        noEt.setNodes(*(pS+n), *(pS+(n+1)%3));
//        if (noE == noEt)
//          T_to_Swap[e].push_back(T3s[t]);
//      }
//    }
//  }
}*/
}

#endif
