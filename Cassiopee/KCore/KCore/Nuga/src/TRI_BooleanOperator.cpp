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

#include "Nuga/include/TRI_BooleanOperator.h"
#include "Nuga/include/TRI_Conformizer.h"

#include "Nuga/include/defs.h"
#include "Nuga/include/TSSplitter.h"
#include "Nuga/include/ContourSplitter.h"
#include "Nuga/include/BbTree.h"

#ifdef DEBUG_TRI_BOOLEAN
#include "IO/io.h"
#include <sstream>
#include <iostream>
#endif

namespace NUGA
{

int TRI_BooleanOperator::case_nb = 0;

#define C_IN 1
#define C_OUT 0
#define C_OVERLAP 2
#define C_NONE -1

#ifdef DEBUG_CONFORMIZER
E_Int NUGA::ConformizerRoot::xtest_counter = 0;
E_Int NUGA::ConformizerRoot::fastdiscard_counter = 0;
E_Int NUGA::ConformizerRoot::split_counter = 0;
E_Int NUGA::ConformizerRoot::split_fastdiscard_counter = 0;
E_Int NUGA::ConformizerRoot::degen_counter = 0;
#endif

//=============================================================================
TRI_BooleanOperator::TRI_BooleanOperator
(const K_FLD::FloatArray& coord1, const K_FLD::IntArray& connect1,
 const K_FLD::FloatArray& coord2, const K_FLD::IntArray& connect2,
 E_Float tolerance, int itermax):parent_type(coord1, connect1, coord2, connect2, tolerance, new TRI_Conformizer<3>(), itermax)
{
  ++case_nb;
  //std::cout << "CASE " << case_nb << std::endl;
}

///
TRI_BooleanOperator::~TRI_BooleanOperator(void)
{
}

///
E_Int TRI_BooleanOperator::get_1_minus_2
(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  colors.clear();

  if (!initialized()) return 1;

  // = _connect_1_out_2 + swapped (_connect_2_in_1)

  coord = _coord;
  connect = _connect_2_in_1;
  colors.resize(connect.cols(), 1);

  for (E_Int i = 0; i < connect.cols(); ++i)
    std::swap(connect(1,i), connect(2,i));

  connect.pushBack(_connect_1_out_2);
  colors.resize(connect.cols(), 0);

  connect.pushBack(_connect_inter_opcoordite_1);
  colors.resize(connect.cols(), 0);

  return 0;
}

///
E_Int TRI_BooleanOperator::get_2_minus_1(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  colors.clear();

  if (!initialized()) return 1;

  // = _connect_2_out_1 + swapped (_connect_1_in_2)

  coord = _coord;
  connect = _connect_1_in_2;
  colors.resize(connect.cols(), 0);

  for (E_Int i = 0; i < connect.cols(); ++i)
    std::swap(connect(1,i), connect(2,i));

  connect.pushBack(_connect_2_out_1);
  colors.resize(connect.cols(), 1);

  K_FLD::IntArray temp = _connect_inter_opcoordite_1;
  NUGA::MeshTool::flipT3(temp);
  
  connect.pushBack(temp);
  colors.resize(connect.cols(), 1);

  return 0;
}

///
E_Int TRI_BooleanOperator::getUnion
(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;

  coord = _coord;
  connect = _connect_1_out_2;
  connect.pushBack(_connect_2_out_1);

  colors.clear();
  colors.resize(_connect_1_out_2.cols(), 0);
  colors.resize(colors.size() + _connect_2_out_1.cols(), 1);

  connect.pushBack(_connect_inter_same);
  colors.resize(connect.cols(), 0);

  return 0;
}

///
E_Int TRI_BooleanOperator::getIntersection
(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;

  coord = _coord;
  connect = _connect_1_in_2;
  connect.pushBack(_connect_2_in_1);

  colors.clear();
  colors.resize(_connect_1_in_2.cols(), 0);
  colors.resize(colors.size() + _connect_2_in_1.cols(), 1);

  connect.pushBack(_connect_inter_same);
  colors.resize(connect.cols(), 1);

  return 0;
}

///
E_Int TRI_BooleanOperator::getIntersectionBorder
(K_FLD::FloatArray& coord, K_FLD::IntArray& connectB)
{
  connectB.clear();
  
  if (!initialized()) return 1;

  coord = _coord;
  
  E_Int Ei[2];
  for (std::set<K_MESH::Edge>::const_iterator i = _hXCO.begin(); i != _hXCO.end(); ++i)
  {
    Ei[0] = i->node(0);
    Ei[1] = i->node(1);
    connectB.pushBack(Ei, Ei+2);
  }
  return 0;
}

///
E_Int TRI_BooleanOperator::__compute_x_contour(std::set<K_MESH::NO_Edge>& hXC)
{
  std::vector<K_MESH::NO_Edge> hE1, hE2, hXE, hXE1, hOv, hC;

  hXC.clear();

  __get_mesh_edges(_connects[0], hE1);
  __get_mesh_edges(_connects[1], hE2);

  // hEX contains all the intersecting edges.
  std::sort(hE1.begin(), hE1.end());
  std::sort(hE2.begin(), hE2.end());
  std::set_intersection(hE1.begin(), hE1.end(), hE2.begin(), hE2.end(), std::back_inserter(hXE));

  // Now get the overlap triangles.
  K_FLD::IntArray connectC, connectCb;
  __getCommonTriangles(connectC);
  __get_mesh_edges(connectC, hOv);

  //
  std::sort(hXE.begin(), hXE.end());
  std::sort(hOv.begin(), hOv.end());
  std::set_difference(hXE.begin(), hXE.end(), hOv.begin(), hOv.end(), std::back_inserter(hXE1));

  hXE = hXE1;

  NUGA::MeshTool::getBoundary(connectC, connectCb);
  //DynArrayIO::write("common.mesh", _coord, connectC);
  for (E_Int c = 0; c < connectCb.cols(); ++c)
    hC.push_back(K_MESH::NO_Edge(connectCb(0,c), connectCb(1,c)));

  // In fine, hXC = hXE + hC
  hXC.insert(hXE.begin(), hXE.end());
  hXC.insert(hC.begin(), hC.end()); // total : X + overlap
  
#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
  /*{
    E_Int Ei[2];
    K_FLD::IntArray connectXB;
    for (std::set<K_MESH::NO_Edge>::const_iterator i = hXC.begin(); i != hXC.end(); ++i)
    {
      Ei[0] = i->node(0);
      Ei[1] = i->node(1);
      connectXB.pushBack(Ei, Ei+2);
    }
    medith::write("xcontour.mesh", _coord, connectXB, "BAR");
  }*/
#endif

  return 0;
}

// Check if surfaces are closed
E_Int TRI_BooleanOperator::check_sanity()
{
  bool is_closed = isClosed(_coord, _connects[0]);
  if (!is_closed)
  {
    printf("Warning: Boolean: surface 0 is not closed.\n");
    return 1;
  } 
  is_closed = isClosed(_coord, _connects[1]);
  if (!is_closed) 
  {
    printf("Warning: Boolean: surface 1 is not closed.\n");
    return 1;
  }
  return 0;
}

///
E_Int TRI_BooleanOperator::compute_zones()
{
  E_Int ret(0), S0, i=0;
  std::set<K_MESH::NO_Edge>  hXC;
  
  _hXCO.clear();

  // Computes the common boundary (intersection edges and overlap zones boundaries).
  ret = __compute_x_contour(hXC);
  if (ret) 
  {
    printf("Warning: Boolean: fail to compute intersection boundary.\n");
    return ret;
  }

  // Get an external element to initialize colors (to 0).
  ret = __find_external_element(_coord, _connects, S0, i);
  if (ret) 
  {
    printf("Warning: Boolean: fail to find initial element for coloring.\n");
    return ret;
  }
  //
  ret = __alternate_coloring(_connects[i], S0, hXC, _colors[i]);
  if (ret) 
  {
    printf("Warning: Boolean: fail to find alternate element for coloring.\n");
    return ret;
  }

  // Fill the relevant containers.
  K_FLD::IntArray *pCout, *pCin, *pCinter;
  K_FLD::IntArray connect_Inter_1, connect_Inter_2;

  if (i == 0)
  {
    pCin  = &_connect_1_in_2;
    pCout = &_connect_1_out_2;
    pCinter = &connect_Inter_1;
  }
  else
  {
    pCin  = &_connect_2_in_1;
    pCout = &_connect_2_out_1;
    pCinter = &connect_Inter_1;
  }

  __fillZoneContainers(_colors[i], _connects[i], pCin, pCout, pCinter);
  
  ret = __reorient_x_contour(*pCout, hXC, _hXCO);
  if (ret) 
  {
    printf("Warning: Boolean: fail to reorient intersection boundary.\n");
    return ret;
  }

  ret = __color_other_part(_connects[(i+1)%2], hXC, _hXCO, _colors[(i+1)%2]);
  if (ret) 
  {
    printf("Warning: Boolean: fail to color zone.\n");
    return ret;
  }

  // Fill the relevant containers.
  if (i == 0)
  {
    pCin = &_connect_2_in_1;
    pCout = &_connect_2_out_1;
    pCinter = &connect_Inter_2;
  }
  else
  {
    pCin = &_connect_1_in_2;
    pCout = &_connect_1_out_2;
    pCinter = &connect_Inter_1;
  }

  i = (i+1)%2;
  __fillZoneContainers(_colors[i], _connects[i], pCin, pCout, pCinter);

  connect_Inter_1.pushBack(connect_Inter_2);
  __splitOverlapByOrientation(connect_Inter_1, _connect_inter_same, _connect_inter_opcoordite_1);

  return ret;
}

///
void
TRI_BooleanOperator::__fillZoneContainers
(const std::vector<E_Int> & colors, const K_FLD::IntArray& connect, 
 K_FLD::IntArray* pCin, K_FLD::IntArray* pCout, K_FLD::IntArray* pCinter)
{
  K_FLD::IntArray::const_iterator pS;
  size_t                          sz(colors.size());

  for (size_t c = 0; c < sz; ++c)
  {
    pS = connect.col(c);
    if (colors[c] == C_OUT) // i.e OUT or OVERLAP.
      pCout->pushBack(pS, pS+3);
    else if (colors[c] == C_IN)
      pCin->pushBack(pS, pS+3);
    else if (colors[c] == C_OVERLAP)
      pCinter->pushBack(pS, pS+3);
  }
}

///
E_Int
TRI_BooleanOperator::__get_mesh_edges
(const K_FLD::IntArray& connect, std::vector<K_MESH::NO_Edge>& no_edges)
{
  NUGA::non_oriented_edge_set_type  hE;
  K_FLD::IntArray::const_iterator         pS;
  E_Int                                   c, n, cols(connect.cols());
  K_MESH::NO_Edge                         Ei;

  no_edges.clear();

  for (c = 0; c < cols; ++c)
  {
    pS = connect.col(c);

    for (n = 0; n < K_MESH::Triangle::NB_NODES; ++n)
    {
      const E_Int& Ni = *(pS+n);
      const E_Int& Nj = *(pS + (n+1) % K_MESH::Triangle::NB_NODES);

      Ei.setNodes(Ni, Nj);

      if (hE.insert(Ei).second)
        no_edges.push_back(Ei);
    }
  }
  return 0;
}

///
E_Int
TRI_BooleanOperator::__find_external_element
(const K_FLD::FloatArray& coord, K_FLD::IntArray* connects, E_Int& S0, E_Int& j)
{
  E_Int                           ret(0), N0, N10, NB_POINTS(coord.cols());
  K_FLD::IntArray::const_iterator pS;
  E_Int                           c, n, nb_tris;

  S0 = IDX_NONE;

  // minN(maxN) will contains the 3 node indices contributing to the min(max) corner.
  // e.g. minN[1] is the index of the node having the lowest y-coordinate.
  {
    E_Int minN[3], maxN[3], minN1[3], maxN1[3];
    NUGA::MeshTool::boundingBox(coord, connects[0], minN, maxN);
    NUGA::MeshTool::boundingBox(coord, connects[1], minN1, maxN1);
 /*  
    for (E_Int k = 0; k < 3; ++k)
    {
      minN[k] = (minN1[k] < minN[k]) ? minN1[k] : minN[k];
      maxN[k] = (maxN[k] < maxN1[k]) ? maxN1[k] : maxN[k];
    }
*/
    N0  = minN[0];  // could be any of the 6 indices of minN and maxN.
    N10  = minN1[0];
    if ((N0 >= NB_POINTS) || (N10 >= NB_POINTS))
      return 1;
    if (coord(0, N0) > coord(0, N10))
      N0 = N10; 
  }

  for (E_Int s = 0; (s < 2) ; ++s)
  {
    K_FLD::IntArray& connect = connects[s];
    nb_tris = connect.cols();

    for (c = 0; (c < nb_tris) ; ++c)
    {
      pS = connect.col(c);
      for (n = 0; n < 3; ++n)
      {
        if (*(pS+n) == N0)
        {
          S0 = c;
          j = s;
          return 0;
        }
      }
    }
  }
  
  return ret;
}

///
E_Int
TRI_BooleanOperator::__alternate_coloring
(const K_FLD::IntArray& connect, E_Int S0, const std::set<K_MESH::NO_Edge>& hXC, std::vector<E_Int>& colors)
{

  // Initialize colors by zone (keeping overlap color).
  std::vector<E_Int>      zones;
  std::vector<E_Int>       zcolors;
  std::set<E_Int>         pool;
  E_Int                   c0, z0;

  // Get the neighbouring of each element.
  NUGA::EltAlgo<K_MESH::Triangle>::NeighbourType    neighbors;
  NUGA::EltAlgo<K_MESH::Triangle>::getNeighbours(connect, neighbors);
  if (neighbors.empty())
    return 1;

  // Get the zones.
  NUGA::EltAlgo<K_MESH::Triangle>::coloring(connect, neighbors, hXC, zones);
  assert(colors.size() == zones.size());

  // Get the zone neighbouring.
  std::map<E_Int, std::set<E_Int> > zone_to_zones;
  __getZoneNeighbouring(connect, zones, zone_to_zones);
  
  // Initialize zone colors and the pool.
  zcolors.resize(*std::max_element(zones.begin(), zones.end())+1, IDX_NONE);
  for (size_t Si = 0; Si < colors.size(); ++Si) // Set overlap color.
  {
    if (colors[Si] == C_OVERLAP)
      zcolors[zones[Si]] = C_OVERLAP;
  }
  z0 = zones[S0];
  zcolors[z0] = C_OUT;
  pool.insert(z0);

  // Coloring loop.
  std::map<E_Int, std::set<E_Int> >::iterator it;
  std::set<E_Int>::iterator itZ;

  while (!pool.empty())
  {
    z0 = *pool.begin();
    c0 = (E_Int)zcolors[z0];
    pool.erase(z0);
    it = zone_to_zones.find(z0);

    if (it == zone_to_zones.end())
      continue;

    for (itZ = it->second.begin(); itZ != it->second.end(); ++itZ)
    {
      z0 = *itZ;
      if (zcolors[z0] == IDX_NONE)
      {
        zcolors[z0] = (E_Int)((c0+1)%2);
        pool.insert(z0);
      }
    }
  }

  // Set colors by element.
  size_t sz = colors.size();
  for (size_t Si = 0; Si < sz; ++Si)
    colors[Si] = zcolors[zones[Si]];

#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
  medith::write("alternate_col.mesh", _coord, connect, "TRI", 0, &colors);
#endif


  return 0;
}

///
void
TRI_BooleanOperator::__getZoneNeighbouring
(const K_FLD::IntArray& connect, 
 std::vector<E_Int>& zones, std::map<E_Int, std::set<E_Int> >& zone_to_zones)
{
  // Set the edge_to_triangles structure.
  NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType bound_to_elts;
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  NUGA::EltAlgo<K_MESH::Triangle>::getBoundToElements(actv, bound_to_elts);
  NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType::iterator it, itEnd(bound_to_elts.end());

  std::set<E_Int> cols;
  std::set<E_Int>::const_iterator itZ, itZ2;
  E_Int i, sz;

  for (it = bound_to_elts.begin(); it != itEnd; ++it)
  {
    cols.clear();
    sz = it->second.size();
    for (i = 0; i < sz; ++i)
    {cols.insert(zones[it->second[i]]);}
    if (cols.size() > 1)
    {
      for (itZ = cols.begin(); itZ != cols.end(); ++itZ)
      {
        itZ2 = itZ;
        ++itZ2;
        for (; itZ2 != cols.end(); ++itZ2)
        {
          zone_to_zones[*itZ].insert(*itZ2);
          zone_to_zones[*itZ2].insert(*itZ);
        }
      }
    }
  }
}

///
E_Int
TRI_BooleanOperator::__reorient_x_contour
(const K_FLD::IntArray& connect, const std::set<K_MESH::NO_Edge>& hXC, std::set<K_MESH::Edge>& hXCO)
{
  E_Int                               ret(0);
  K_FLD::IntArray                     boundary;
  
  K_MESH::Edge                        E;
  std::set<K_MESH::Edge>              temp;

  hXCO.clear();

  if (hXC.empty())
    return 0;

  hXCO.insert(hXC.begin(), hXC.end());

  // Overlap contours included.
  K_FLD::IntArray connectU(connect), connectC;
  __getCommonTriangles(connectC);
  connectU.pushBack(connectC);
  NUGA::MeshTool::getBoundary(connectU, boundary);
    
  E_Int cols = boundary.cols();
  for (E_Int c = 0; c < cols; ++c)
  {
    E.setNodes(boundary.col(c));
    if (hXCO.erase(E))
      temp.insert(K_MESH::Edge(E.node(1), E.node(0)));
  }

  hXCO.insert(temp.begin(), temp.end());

#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
  {
    E_Int Ei[2];
    K_FLD::IntArray connectXB;
    std::set<K_MESH::NO_Edge> tmp;
    for (std::set<K_MESH::Edge>::const_iterator i = hXCO.begin(); i != hXCO.end(); ++i)
    {
      Ei[0] = i->node(0);
      Ei[1] = i->node(1);
      
      if (Ei[0] == Ei[1])
        std::cout << "degenrated" << std::endl;
      
      tmp.insert(K_MESH::NO_Edge(i->node(0), i->node(1)));
      
      connectXB.pushBack(Ei, Ei+2);
    }

    std::cout << "number of edges : " << hXCO.size() << std::endl;
    std::cout << "number of oriented edges : " << tmp.size() << std::endl;
    
    medith::write("reoriented_xcontour.mesh", _coord, connectXB, "BAR");
  }
#endif

  return ret;
}

///
E_Int
TRI_BooleanOperator::__color_other_part
(K_FLD::IntArray& connect, const NUGA::non_oriented_edge_set_type& hXC,
 const NUGA::oriented_edge_set_type& hXCO, std::vector<E_Int>& colors)
{
  E_Int                         ret(0), c, col;
  K_FLD::IntArray               connectbo;
  bool                          same_orient, opp_orient;
  K_MESH::Edge                  E, rE;
  K_FLD::IntArray               connectOut, Ci;
  std::vector<E_Int>            colorsOut;

  // Fast return.
  if (hXC.empty() || hXCO.empty())
  {
    colors.resize(connect.cols(), 0);
    return 0;
  }
  
#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
  {
    K_FLD::IntArray tmp;
    E_Int Ei[2];
    for (NUGA::non_oriented_edge_set_type::const_iterator i = hXC.begin(); i != hXC.end(); ++i)
    {
      Ei[0] = i->node(0);
      Ei[1] = i->node(1);
      tmp.pushBack(Ei, Ei+2);
    }
    medith::write("hXC.mesh", _coord, tmp);
  }
#endif

  // Detect the different zones.
  std::vector< std::vector<E_Int> > scolors;
  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(connect, hXC, scolors);
  
  // Loop through the zone and compare the boundary orientation with connectXB's.
  E_Int nb_zones = scolors.size();
  for (E_Int z = 0; z < nb_zones; ++z)
  {
    
    Ci.clear();
    for (size_t i = 0; i < scolors[z].size(); ++i)
      Ci.pushBack(connect.col(scolors[z][i]), connect.col(scolors[z][i])+3);
    
#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
    std::ostringstream o;
    o << "m_" << z << ".mesh";
    medith::write(o.str().c_str(), _coord, Ci);
#endif

    connectOut.pushBack(Ci);

    col = C_OVERLAP;
    if (colors[scolors[z][0]] != C_OVERLAP) // Overlap zone (already colored).
    {
      //
      NUGA::MeshTool::getBoundary(Ci, connectbo);
      same_orient = opp_orient = false;
      col = C_NONE;
      
#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
      {
        std::ostringstream o;
        o << "zb_" << z << ".mesh";
        medith::write(o.str().c_str(), _coord, connectbo);
      }
#endif

      for (c = 0; (c<connectbo.cols()) && (col == C_NONE); ++c)
      {
        E.setNodes(connectbo.col(c));
        same_orient = (hXCO.find(E) != hXCO.end());
        if (!same_orient)
        {
          rE.setNodes(E.node(1), E.node(0));
          opp_orient = (hXCO.find(rE) != hXCO.end());
        }

        if (same_orient || opp_orient) // found in any orientation.
          col = same_orient ? C_OUT : C_IN;// if same orientation -> out
      }
      assert (col != C_NONE);
    }

    colorsOut.resize(connectOut.cols(), col);
  }

  connect = connectOut;
  colors = colorsOut;

#ifdef E_DEBUG_TRI_BOOLEANOPERATOR
  medith::write("other_part_col.mesh", _coord, connect, "TRI", 0, &colors);
#endif

  return ret;
}

///
E_Int
TRI_BooleanOperator::__getCommonTriangles(K_FLD::IntArray& connectC)
{
  K_FLD::IntArray       connectG(_connects[0]);
  std::vector<E_Int>    dupIds;
  E_Int                 did, p1, p2, Si, Sj, c0(_connects[0].cols());

  connectC.clear();
  connectG.pushBack(_connects[1]);

  _colors[0].resize(_connects[0].cols(), C_NONE);
  _colors[1].resize(_connects[1].cols(), C_NONE);

  NUGA::MeshTool::detectDuplicated(connectG, dupIds, false /*strict orient*/);

  E_Int sz = (E_Int)dupIds.size();
  for (E_Int i = 0; i < sz; ++i)
  {
    did = dupIds[i];
    if (did != i)
    {
      connectC.pushBack(connectG.col(did), connectG.col(did)+3);

      p1 = (i < c0) ? 0 : 1;
      p2 = (did < c0) ? 0 : 1; 

      Si = (i < c0) ? i : i - c0;
      Sj = (did < c0) ? did : did - c0;

      _colors[p1][Si] = _colors[p2][Sj] = C_OVERLAP;
    }
  }
  
  return connectC.cols();
}

///
void
TRI_BooleanOperator:: __splitOverlapByOrientation
(const K_FLD::IntArray& connectInter, K_FLD::IntArray& connectIS, K_FLD::IntArray& connectIO)
{
  std::vector<E_Int> nids;

  NUGA::MeshTool::detectDuplicated(connectInter, nids, true/*strict orient*/);
  E_Int end = (E_Int)nids.size(), start = end / 2;
  for (E_Int i = start; i < end; ++i)
  {
    if (nids[i] != i)
      connectIS.pushBack(connectInter.col(i), connectInter.col(i)+3);
    else
      connectIO.pushBack(connectInter.col(i), connectInter.col(i)+3);
  }
  NUGA::MeshTool::removeDuplicated(connectIO, nids, false/*strict orient*/);

  assert (2 *(connectIO.cols() + connectIS.cols()) == connectInter.cols());

}

bool
TRI_BooleanOperator::isClosed(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect)
{
  NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType bound_to_elts;
  NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType::const_iterator it;
  
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  //maxNeigh = 
  NUGA::EltAlgo<K_MESH::Triangle>::getBoundToElements(actv, bound_to_elts);
  //if (maxNeigh != 2) return false;
  //minNeigh = maxNeigh;

  bool is_not_close = false;

  std::vector<E_Int> color_pathos(connect.cols(), 0);
  std::vector<bool> mask(connect.cols(), false);

  for (it = bound_to_elts.begin(); it != bound_to_elts.end(); ++it)
  {
    //minNeigh = std::min(E_Int(it->second.size()), minNeigh);
    if (it->second.size() == 1)
    {
      color_pathos[*it->second.begin()] = 1;
      mask[*it->second.begin()] = true;
      is_not_close = true;
    }
  }
#ifdef DEBUG_TRI_BOOLEAN
  if (is_not_close)
  {
    static int s_file_count = 0;
    std::ostringstream o;
    o << "unclosed_" << s_file_count << ".mesh";
    medith::write(o.str().c_str(), coord, connect, "TRI", 0, &color_pathos);
    
    o.str("");
    o << "only_pathos_" << s_file_count++ << ".mesh";
    medith::write(o.str().c_str(), coord, connect, "TRI", &mask);
    std::cout << "UNCLOSED !!!! " << std::endl << std::endl;
    
    // attached to pathos
    K_FLD::IntArray::const_iterator pS;
    K_MESH::NO_Edge E;
    std::vector<bool> new_mask(mask.size(), false);
    for (size_t i=0; i < connect.cols(); ++i)
    {
      if (!mask[i])
        continue;
      
      pS = connect.col(i);
      for (size_t n=0; n < 3; ++n)
      {
        E.setNodes(*(pS+n), *(pS+(n+1)%3));
        it = bound_to_elts.find(E);
        if (it != bound_to_elts.end())
          for (size_t j=0; j < it->second.size(); ++j)
            new_mask[it->second[j]]=true;
      }
    }
    o.str("");
    o << "pathos_and_Co_" << s_file_count++ << ".mesh";
    medith::write(o.str().c_str(), coord, connect, "TRI", &new_mask);
  }
#endif
  return !is_not_close;
}


}

