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

#if defined(DEBUG_CONFORMIZER) || defined(DEBUG_BOOLEAN)

#include "Nuga/include/TRI_debug.h"
#include "Nuga/include/medit.hxx"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/IdTool.h"
#include "Nuga/include/EltAlgo.h"
#include <iostream>
///

std::string medith::wdir = "";


void TRI_debug::draw_connected_to_node_T3s
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni)
{
  K_FLD::IntArray tmp;
  for (size_t i = 0; i < connectT3.cols(); ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      if (connectT3(j,i) == Ni)
      {
        tmp.pushBack(connectT3.col(i), connectT3.col(i)+3);
        break;
      }
    }
  }
  std::ostringstream o;
  o << "T3_connected_to_" << Ni << "_toT3.mesh";
  medith::write(o.str().c_str(), coord, tmp);
}

///
void TRI_debug::get_connected_to_T3
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ti, std::vector<E_Int> & cIds)
{
  if (Ti >= connectT3.cols())
    return;
  
  cIds.clear();
  
  E_Int Ni = connectT3(0,Ti);
  E_Int Nj = connectT3(1,Ti);
  E_Int Nk = connectT3(2,Ti);
  
  std::set<K_MESH::NO_Edge> edges;
  edges.insert(K_MESH::NO_Edge(Ni, Nj));
  edges.insert(K_MESH::NO_Edge(Nj, Nk));
  edges.insert(K_MESH::NO_Edge(Ni, Nk));
  
  K_MESH::NO_Edge E;
  for (size_t i = 0; i < connectT3.cols(); ++i)
  {
    if (i == Ti)
      continue;
    for (size_t j = 0; j < 3; ++j)
    {
      E.setNodes(connectT3(j, i), connectT3((j+1)%3, i));
      if (edges.find(E) != edges.end())
      {
        cIds.push_back(i);
        break;
      }
    }
  }
}
///
void TRI_debug::draw_connected_to_T3
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ti)
{
  if (Ti >= connectT3.cols())
    return;
  
  K_FLD::IntArray tmp;
  tmp.pushBack(connectT3.col(Ti), connectT3.col(Ti)+3);
  std::vector<E_Int> colors(1, COL_RED), cIds;
  
  TRI_debug::get_connected_to_T3(coord, connectT3, Ti, cIds);
  
  K_MESH::NO_Edge E;
  for (size_t i = 0; i < cIds.size(); ++i)
  {
    tmp.pushBack(connectT3.col(cIds[i]), connectT3.col(cIds[i])+3);
    colors.resize(colors.size()+1, COL_DEFAULT);
  }
  
  //
  std::ostringstream o;
  o << "connected_to_T3_" << Ti << ".mesh";
  medith::write(o.str().c_str(), coord, tmp, "TRI", 0, &colors);
}

///

void TRI_debug::draw_connected_to_T3
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk/*, bool shrink*/)
{
  std::vector<E_Int> colors;
  K_FLD::IntArray connOut;
  K_FLD::FloatArray coordOut;
  connected_to_T3(coord, connectT3, Ni, Nj, Nk, /*shrink,*/ coordOut, connOut, colors);
  //
  std::ostringstream o;
  o << "connected_to_T3.mesh" ;
  medith::write(o.str().c_str(), coordOut, connOut, "TRI", 0, &colors);
}

///

void TRI_debug::connected_to_T3
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk/*, bool shrink*/,
       K_FLD::FloatArray& coordOut, K_FLD::IntArray& connOut, std::vector<E_Int>& colors)
{
  std::set<K_MESH::NO_Edge> edges;
  edges.insert(K_MESH::NO_Edge(Ni, Nj));
  edges.insert(K_MESH::NO_Edge(Nj, Nk));
  edges.insert(K_MESH::NO_Edge(Ni, Nk));
  
  connOut.clear();
  colors.clear();
  E_Int Ti = IDX_NONE;
  
  //const K_FLD::FloatArray *pC = &coord;
  coordOut=coord;
    
  K_MESH::NO_Edge E, E2, E3;
  for (size_t i = 0; i < connectT3.cols(); ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      E.setNodes(connectT3(j, i), connectT3((j+1)%3, i));
      if (edges.find(E) != edges.end())
      {
        connOut.pushBack(connectT3.col(i), connectT3.col(i)+3);
        E2.setNodes(connectT3((j+1)%3, i), connectT3((j+2)%3, i));
        E3.setNodes(connectT3((j+2)%3, i), connectT3(j%3, i));
        
        if ((edges.find(E2) != edges.end()) && (edges.find(E3) != edges.end()))
        {
          colors.resize(colors.size()+1, COL_RED);
          Ti = connOut.cols()-1;
        }
        else
          colors.resize(colors.size()+1, COL_DEFAULT);
        break;
      }
    }
  }
  
  /*if (shrink)
  {
    coordCpy = coord;
    pC=&coordCpy;
    
    E_Float L = 0., Lmax=0.;
    K_MESH::NO_Edge Emax;
    for (std::set<K_MESH::NO_Edge>::const_iterator it = edges.begin(); it != edges.end(); ++it)
    {
      L = std::max(L, NUGA::sqrDistance(coord.col(it->node(0)), coord.col(it->node(1)), 3));
      if (Lmax < L)
      {
        Lmax = L;
        Emax = *it;
      }
    }
    L = sqrt(L);
    
    for (size_t i = 0; i < tmp.cols(); ++i)
    {
      if (Ti == i)
        continue;
      size_t j=0;
      for (; j<3; ++j)
      {
        E.setNodes(tmp(j, i), tmp((j+1)%3, i));
        if (Ti != i)
        {
          if (edges.find(E) != edges.end())
            break;
        }
        else
        {
          if (E == Emax)
          {
            break;
          }
        }
      }
      
      const E_Int& N0 = E.node(0);
      const E_Int& N1 = E.node(1);
      const E_Int& N2 = tmp((j+2)%3, i);
      
      E_Float N0N1[3];
      NUGA::diff<3>(coordCpy.col(N1), coordCpy.col(N0), N0N1);
      //E_Float l=sqrt(NUGA::sqrNorm<3>(N0N1));
      E_Float Ph[3];
      NUGA::sum<3>(0.5, N0N1, coordCpy.col(N0), Ph);
      
      E_Float PhN2[3];
      NUGA::diff<3>(coordCpy.col(N2), Ph, PhN2);
      NUGA::normalize<3>(PhN2);
      
      NUGA::sum<3>(L, PhN2, Ph, coordCpy.col(N2));
    }
  }*/
}

void TRI_debug::draw_connected_to_E2
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj)
{
  K_FLD::IntArray tmp;
  for (size_t i=0; i < connectT3.cols(); ++i)
  {
    if ( (connectT3(0,i) == Ni || connectT3(1,i) == Ni || connectT3(2,i) == Ni) &&
         (connectT3(0,i) == Nj || connectT3(1,i) == Nj || connectT3(2,i) == Nj) )
      tmp.pushBack(connectT3.col(i), connectT3.col(i)+3); 
  }
  medith::write("attached_toE.mesh", coord, tmp);
}

///

E_Int TRI_debug::get_T3_index
(const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk)
{
  std::set<E_Int> nodes;
  nodes.insert(Ni); nodes.insert(Nj); nodes.insert(Nk);
  E_Int id = IDX_NONE;
  for (size_t i=0; i < connectT3.cols(); ++i)
  {
    if ((nodes.find(connectT3(0, i)) != nodes.end()) && (nodes.find(connectT3(1, i)) != nodes.end()) && (nodes.find(connectT3(2, i)) != nodes.end()))
      id = i;
  }
  
  return id;
  
}

///

void TRI_debug::remove_T3
(const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, E_Int Ti)
{
  K_FLD::IntArray tmp;
  for (size_t i=0; i < connectT3.cols(); ++i)
  {
    if (Ti != i)
      tmp.pushBack(connectT3.col(i), connectT3.col(i)+3);
  }
  
  connectT3=tmp;
}

///

void TRI_debug::remove_T3
(const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk)
{
  E_Int Ti = get_T3_index(connectT3, Ni, Nj, Nk);
  remove_T3(coord, connectT3, Ti);
}

///
void TRI_debug::get_T3_neighbors
(const char* fname, E_Int Ti, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, bool both_orient)
{
  K_FLD::IntArray connect;
  std::vector<E_Int> colors;
  
  connect.pushBack(connectT3.col(Ti), connectT3.col(Ti)+3);
  colors.resize(1,COL_RED);
  
  for (size_t i=0; i < 3; ++i)
  {
    E_Int Tngh =  neighbors(i, Ti);
    connect.pushBack(connectT3.col(Tngh), connectT3.col(Tngh)+3);
    colors.resize(colors.size()+1,COL_GREEN);
  }
  
  if (both_orient)
  {
    E_Int stride = connectT3.cols()/2;
    E_Int Topp = (Ti <stride) ? Ti + stride : Ti -stride;
    connect.pushBack(connectT3.col(Topp), connectT3.col(Topp)+3);
    colors.resize(colors.size()+1,COL_NEXT);
    
    for (size_t i=0; i < 3; ++i)
    {
      E_Int Tngh =  neighbors(i, Topp);
      connect.pushBack(connectT3.col(Tngh), connectT3.col(Tngh)+3);
      colors.resize(colors.size()+1,COL_NEXT2);
    }
  }
  medith::write(fname, coord, connect, "TRI", 0, &colors);
}

void TRI_debug::coloring_frames
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, E_Int Kseed)
{
  size_t cols(neighbors.size()), K, Kn, color(1), sz;
  std::vector<E_Int> cpool;
  std::vector<E_Int>::const_iterator  itC;
  K_FLD::IntArray neighs;
  E_Int count(0);
  
  E_Int NONE_COL = 0;
  
  std::vector<E_Int> colors(connectT3.cols(), NONE_COL);
  std::vector<bool> keep(connectT3.cols(), false);
  
  std::ostringstream o;

 
  cpool.push_back(Kseed);

  while (!cpool.empty())
  {
    K = cpool.back();
    cpool.pop_back();
    
    if (colors[K] != NONE_COL)
      continue;
    
    
    colors[K] = color;
    keep[K]=true;
    
    o.str("");
    o << "frame_";
    if (count < 10) o << "0";
    o << count++ << ".mesh";
    medith::write(o.str().c_str(), coord, connectT3, "TRI", &keep/*, &colors*/);
    
    for (size_t i = 0; i < 3; ++i)
    {
      Kn = neighbors(i, K);

      if ((Kn != IDX_NONE) && (colors[Kn] == NONE_COL)) // Not colored.
        cpool.push_back(Kn);
    }
  } 
}

void TRI_debug::get_same_ancestor_T3s
(E_Int Ti, const std::vector<E_Int>& ancestors, std::vector<E_Int>& oids)
{  
  oids.clear();
  
  for (size_t i = 0; i < ancestors.size(); ++i)
  {
    if (ancestors[i] == Ti)
      oids.push_back(i);
  }
}

void TRI_debug::draw_same_ancestor_T3s
(E_Int Ti, E_Int iter, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::vector<E_Int>& ancestors)
{
  std::ostringstream o;
  o << "PG_" << Ti << "_it" << iter << ".mesh";
  
  std::vector<E_Int> oids;
  TRI_debug::get_same_ancestor_T3s(Ti, ancestors, oids);
  
  K_FLD::IntArray tmp;
  tmp.append_selection(connectT3, oids);
  
  medith::write(o.str().c_str(), coord, tmp, "TRI");
  
}

///
void TRI_debug::write_wired(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, bool normal, const std::vector<E_Int>* T3colors, const std::vector<bool>* keep)
{
  K_FLD::FloatArray crd(coord);
  K_FLD::IntArray connectE;
  K_MESH::Edge E;
  E_Int n0, n1;
  E_Float Ni[3], P0[3], P1[3], FACTOR, Lmin, L2;
  std::vector<E_Int> E2colors;
  for (size_t i = 0; i < connectT3.cols(); ++i)
  {
    if (keep && (*keep)[i] == false)
      continue;
    
    Lmin = NUGA::FLOAT_MAX;
    for (size_t j=0; j<3;++j)
    {
      E.setNodes(connectT3(j,i), connectT3((j+1)%3,i));
      connectE.pushBack(E.begin(), E.begin()+2);
      L2 = NUGA::sqrDistance(coord.col(E.node(0)), coord.col(E.node(1)), 3);
      Lmin = (L2 < Lmin) ? L2 : Lmin;
    }
    
    Lmin = 0.5*sqrt(Lmin);
    
    if (T3colors)
      E2colors.resize(E2colors.size()+3, (*T3colors)[i]);
      
    if (normal)
    {    
      K_MESH::Triangle::isoG(coord, connectT3.col(i), P0); // P0 is the center of Ti
      K_MESH::Triangle::normal(coord, connectT3.col(i), Ni);
      NUGA::normalize<3>(Ni);
      //NUGA::sum<3>(P0, Ni, P1);
      NUGA::sum<3>(1., P0, Lmin, Ni, P1);
      crd.pushBack(P0, P0+3);
      n0=crd.cols()-1;
      crd.pushBack(P1, P1+3);
      n1=crd.cols()-1;
      E.setNodes(n0,n1);
      connectE.pushBack(E.begin(), E.begin()+2);
      
      if (T3colors)
        E2colors.resize(E2colors.size()+1, 0); 
    }
  }
    
  std::vector<E_Int>* colors = 0;
  if (T3colors)
    colors=&E2colors;
  
  medith::write(fname, crd, connectE, "BAR", 0, colors);
}

///
void TRI_debug::write_wired(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const K_FLD::FloatArray& normals, const std::vector<E_Int>* T3colors, const std::vector<bool>* keep, bool localid)
{
  K_FLD::FloatArray crd(coord);
  K_FLD::IntArray connectE;
  K_MESH::Edge E;
  E_Int n0, n1;
  E_Float P0[3], P1[3];
  const E_Float* Ni;
  std::vector<E_Int> E2colors;
  for (size_t i = 0; i < connectT3.cols(); ++i)
  {
    if (keep && (*keep)[i] == false)
      continue;
    for (size_t j=0; j<3;++j)
    {
      E.setNodes(connectT3(j,i), connectT3((j+1)%3,i));
      connectE.pushBack(E.begin(), E.begin()+2);
    }
    
    if (T3colors)
      E2colors.resize(E2colors.size()+3, (*T3colors)[i]);
      
    K_MESH::Triangle::isoG(coord, connectT3.col(i), P0); // P0 is the center of Ti
    Ni = normals.col(i);
    NUGA::sum<3>(1., P0, 0.005, Ni, P1);
    crd.pushBack(P0, P0+3);
    n0=crd.cols()-1;
    crd.pushBack(P1, P1+3);
    n1=crd.cols()-1;
    E.setNodes(n0,n1);
    connectE.pushBack(E.begin(), E.begin()+2);
      
    if (T3colors)
      E2colors.resize(E2colors.size()+1, 0); 
  }
    
  std::vector<E_Int>* colors = 0;
  if (T3colors)
    colors=&E2colors;
  
  std::vector<E_Int> nids;
  K_FLD::FloatArray *pCrd(&crd), tmpCrd;
  if (localid){
    tmpCrd=crd;
    NUGA::MeshTool::compact_to_mesh(tmpCrd, connectE, nids);
    pCrd=&tmpCrd;
    if (colors) K_CONNECT::IdTool::compact(*colors, nids);
  }
  
  medith::write(fname, *pCrd, connectE, "BAR", 0, colors);
}

///
bool TRI_debug::analyze_T3_set(E_Int setid, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const std::vector<E_Int>& T3s, bool& is_manifold, bool& is_closed, std::vector<E_Int>* colors)
{
  bool healthy = true;
  is_manifold = is_closed = true;

  std::map<E_Int, K_FLD::IntArray> col_to_cnt;
  std::map<E_Int, K_FLD::IntArray>::const_iterator it;

  std::ostringstream o;

  // write out this set
  {
    K_FLD::IntArray tmp;
    tmp.append_selection(cnt, T3s);
    if (tmp.cols())
    {
      o << "faulty_parasite_" << setid << ".mesh";
      medith::write(o.str().c_str(), crd, tmp, "TRI", 0, colors);
    }
    else
      return healthy;

    if (!colors)
      col_to_cnt[0] = tmp;
  }

  if (colors)
  {
    for (size_t i = 0; i < T3s.size(); ++i)
    {
      const E_Int& Ti = T3s[i];
      E_Int Ci = (i < colors->size()) ? (*colors)[i] : 0;
      Ci = (Ci == IDX_NONE || Ci < 0) ? 0 : Ci;
      col_to_cnt[Ci].pushBack(cnt.col(Ti), cnt.col(Ti) + 3);
    }
  }

  typedef  NUGA::EltAlgo<K_MESH::Triangle> algo_t;
  typedef K_FLD::ArrayAccessor<K_FLD::IntArray> acnt_t;

  K_FLD::IntArray cB_Free, cB_NMnfld;
  algo_t::BoundToEltType b_to_e;
  algo_t::BoundToEltType::const_iterator itE;

  for (it = col_to_cnt.begin(); it != col_to_cnt.end(); ++it)
  {
    const E_Int& color = it->first;
    const K_FLD::IntArray& cn = it->second;
    
    acnt_t acnt(cn);
    b_to_e.clear();
    
    algo_t::getBoundToElements(acnt, b_to_e);

    cB_Free.clear();
    cB_NMnfld.clear();

    // get non manifold edges and freee edges
    for (itE = b_to_e.begin(); itE != b_to_e.end(); ++itE)
    {
      const K_MESH::NO_Edge& edge = itE->first;
      const std::vector<E_Int> & t3set = itE->second;

      if (t3set.size() == 1)
        cB_Free.pushBack(edge.begin(), edge.end());
      else if (t3set.size() > 2)
        cB_NMnfld.pushBack(edge.begin(), edge.end());
    }
    
    std::vector<E_Int> nids;

    //output free edges
    o.str("");
    o << "frees_for_col_" << color << ".mesh";
    K_FLD::FloatArray tmpCrd=crd;
    NUGA::MeshTool::compact_to_mesh(tmpCrd, cB_Free, nids);
    medith::write(o.str().c_str(), tmpCrd, cB_Free, "BAR");
    
    is_closed = (cB_Free.cols()==0);
 
    //output non manifold edges
    o.str("");
    o << "non_mnfld_for_col_" << color << ".mesh";
    tmpCrd=crd;
    NUGA::MeshTool::compact_to_mesh(tmpCrd, cB_NMnfld, nids);
    medith::write(o.str().c_str(), tmpCrd, cB_NMnfld, "BAR");
    
    is_manifold = (cB_NMnfld.cols()==0);

  }
  
  healthy = is_manifold && is_closed;

  return healthy;
}

#endif
