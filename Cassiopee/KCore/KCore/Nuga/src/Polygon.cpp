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

#include "Nuga/include/Polygon.h"
#include "Nuga/include/IdTool.h"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/GeomAlgo.h"
#include <vector>
#define Vector_t std::vector

#ifdef DEBUG_POLYGON
#include "IO/DynArrayIO.h"
#include <iostream>
#endif

namespace K_MESH
{ 

//=============================================================================
/*void Polygon::getBoundary
(const Polygon&  T1, const Polygon&  T2, K_MESH::NO_Edge& b)
{
  K_MESH::NO_Edge Ei, Ej;
  for (E_Int i = 0; i < NB_NODES; ++i)
  {
    T1.getBoundary(i, Ei);
    for (E_Int j = 0; j < NB_NODES; ++j)
    {
      T2.getBoundary(j, Ej);
      if (Ei == Ej)
      {
        b = Ei; return;
      }
    }
  }
}
*/
//=============================================================================
void Polygon::getBoundary
(const Polygon&  T1, const Polygon&  T2, E_Int& i1, E_Int& i2)
{
  K_MESH::NO_Edge Ei, Ej;
  for (i1=0; i1 < T1._nb_nodes; ++i1)
  {
    T1.getBoundary(i1, Ei);
    for (i2=0; i2 < T2._nb_nodes; ++i2)
    {
      T2.getBoundary(i2, Ej);
      if (Ei == Ej)
        return;
    }
  }
}

//=============================================================================
E_Int Polygon::getOrientation
(const Polygon&  PG, const E_Int& Ni, const E_Int& Nj, E_Bool& same_orient)
{
  same_orient = false;
 
  for (E_Int n = 0; n < PG._nb_nodes; ++n)
  {
    if ((PG._nodes[n] == Ni) && (PG._nodes[(n+1)%PG._nb_nodes] == Nj))
    {
      same_orient = true;
      return 0;
    }
    if ((PG._nodes[n] == Nj) && (PG._nodes[(n+1)%PG._nb_nodes] == Ni))
      return 0;
  }

  return -1;
}

// ///
// bool Polygon:: is_convex
// (const K_FLD::FloatArray& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start,
// const E_Float* normal, E_Float PI_ratio, E_Int& iworst)
// {
//   //
//   E_Float Ei[3], Ej[3];
//   bool convex = true;
//   iworst = IDX_NONE;
  
//   E_Float angle_max = NUGA::PI* PI_ratio; // a fraction betwen 0 and Pi
//   E_Float cos_min = ::cos(angle_max); // cos is decreasing on [0; Pi]
//   E_Float Z[3];
//   for (E_Int i = 1; i < nb_nodes + 1; ++i)
//   {
//     E_Int ei = nodes[i%nb_nodes] - index_start;
//     E_Int eim1 = nodes[i - 1] - index_start;
//     E_Int eip1 = nodes[(i + 1) % nb_nodes] - index_start;

//     NUGA::diff<3>(coord.col(ei), coord.col(eim1), &Ei[0]); //fixme : no normalization ??
//     NUGA::diff<3>(coord.col(eip1), coord.col(ei), &Ej[0]);
//     NUGA::sum<3>(normal, coord.col(ei), Z);

//     E_Float det = NUGA::zzdet4(coord.col(eim1), coord.col(ei), coord.col(eip1), Z);
    
//     if (det >= 0.) continue; // convex
    
//     NUGA::normalize<3>(Ei);
//     NUGA::normalize<3>(Ej);

//     E_Float c = NUGA::dot<3>(Ei, Ej);
    
//     if (c < cos_min) // angle > anle max
//     {
//       convex = false;
//       iworst = i%nb_nodes;
//       cos_min = c;
//     }
//   }

//   return convex;
// }

///
bool Polygon:: is_convex
(const K_FLD::FloatArray& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start,
const E_Float* normal, E_Float convexity_tol, E_Int& iworst, E_Int& ibest)
{
  //
  E_Float Ei[3], Ej[3];
  bool convex = true;
  ibest = iworst = IDX_NONE;
  
  E_Float Z[3], det_min(-convexity_tol), det_max(0.);
  for (E_Int i = 1; i < nb_nodes + 1; ++i)
  {
    E_Int ei = nodes[i%nb_nodes] - index_start;
    E_Int eim1 = nodes[i - 1] - index_start;
    E_Int eip1 = nodes[(i + 1) % nb_nodes] - index_start;

    NUGA::diff<3>(coord.col(ei), coord.col(eim1), &Ei[0]); //fixme : no normalization ??
    NUGA::diff<3>(coord.col(eip1), coord.col(ei), &Ej[0]);
    NUGA::sum<3>(normal, coord.col(ei), Z);

    E_Float det = NUGA::zzdet4(coord.col(eim1), coord.col(ei), coord.col(eip1), Z);

    det /= (NUGA::normalize<3>(Ei)*NUGA::normalize<3>(Ej)); //normalization to really have a angular-based test.
    //std::cout << det << std::endl;
    if (det < det_min)
    {
      convex = false;
      iworst = i%nb_nodes;
      det_min = det;
    }
    
    if (det > det_max)
    {
      ibest = i%nb_nodes;
      det_max = det;
    }
  }

  return convex;
}

/// Predicate
bool Polygon::is_convex
(const K_FLD::FloatArray& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, const E_Float *normal, E_Float convexity_tol)
{
  E_Float Z[3], Ei[3], Ej[3];
  for (E_Int i = 1; i < nb_nodes + 1; ++i)
  {
    E_Int ei = nodes[i%nb_nodes] - index_start;
    E_Int eim1 = nodes[i - 1] - index_start;
    E_Int eip1 = nodes[(i + 1) % nb_nodes] - index_start;

    NUGA::diff<3>(coord.col(ei), coord.col(eim1), &Ei[0]);
    NUGA::diff<3>(coord.col(eip1), coord.col(ei), &Ej[0]);

    /*NUGA::normalize<3>(Ei);
    NUGA::normalize<3>(Ej);
    NUGA::crossProduct<3>(Ei, Ej, Ek);

    ps = NUGA::dot<3>(Ek, normals.col(globalId[K]));
    */
    NUGA::sum<3>(normal, coord.col(ei), Z);

    E_Float uu = NUGA::zzdet4(coord.col(eim1), coord.col(ei), coord.col(eip1), Z);
    /*
    if ((ps*uu < 0.) && (::fabs(ps) > EPSILON) && (::fabs(uu) > EPSILON))
    {
    assert(false);
    }*/

    if (uu < -convexity_tol) //concavity : fixme need also a better criterion base on angles
      return false;
  }

  return true;
}

///
E_Int Polygon::build_pg_neighborhood
(const ngon_unit& PGS, ngon_unit& neighbor,
const E_Int* first_pg, E_Int nb_pgs,
const std::set<K_MESH::NO_Edge>* wall/*extra specified walls*/)
{
  neighbor.clear();
  PGS.updateFacets();

  if (nb_pgs == 0 || first_pg == 0) // consider all pgs
    neighbor = PGS; //same molecules as the connectivity
  else
  {
    neighbor.clear();
    neighbor._NGON.resize(2, 0);
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      E_Int nb_nodes = PGS.stride(PGi);
      const E_Int* begin = PGS.get_facets_ptr(PGi) - 1;
      const E_Int* end = begin + nb_nodes + 1;
      neighbor._NGON.insert(neighbor._NGON.end(), begin, end);
    }
    neighbor._NGON[0] = nb_pgs;
    neighbor._NGON[1] = neighbor._NGON.size() - 2;
    assert(neighbor.size() == nb_pgs);
  }

  neighbor.reset_facets();
  

  std::map<K_MESH::NO_Edge, std::pair<K_MESH::Edge, K_MESH::Edge> > noe_to_bound;//oriented edge to its position in the PG
  K_MESH::NO_Edge E;

  std::map<K_MESH::NO_Edge, std::pair<K_MESH::Edge, K_MESH::Edge> >::iterator it;

  if (nb_pgs && first_pg != 0)
  {
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      E_Int nb_nodes = PGS.stride(PGi);
      const E_Int* nodes = PGS.get_facets_ptr(PGi);

      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        E.setNodes(*(nodes + j), *(nodes + (j + 1) % nb_nodes));
        it = noe_to_bound.find(E);
        if (it == noe_to_bound.end())
        {
          //it = noe_to_bound[E];
          noe_to_bound[E].first = K_MESH::Edge(i, j);
          noe_to_bound[E].second = K_MESH::Edge(-1, -1);
        }
        else
        {
          if (it->second.second.node(0) == -1) //first time
          {
            it->second.second = K_MESH::Edge(i, j);
          }
          else //non-manifold so treat as a cut
            it->second.second = K_MESH::Edge(IDX_NONE, IDX_NONE);
        }
      }
    }
  }
  else
  {
    nb_pgs = PGS.size();
    for (E_Int PGi = 0; PGi < nb_pgs; ++PGi)
    {
      E_Int nb_nodes = PGS.stride(PGi);
      const E_Int* nodes = PGS.get_facets_ptr(PGi);

      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        E.setNodes(*(nodes + j), *(nodes + (j + 1) % nb_nodes));
        it = noe_to_bound.find(E);
        if (it == noe_to_bound.end())
        {
          //it = noe_to_bound[E];
          noe_to_bound[E].first = K_MESH::Edge(PGi, j);
          noe_to_bound[E].second = K_MESH::Edge(-1, -1);
        }
        else
        {
          if (it->second.second.node(0) == -1) //first time
          {
            it->second.second = K_MESH::Edge(PGi, j);
          }
          else //non-manifold so treat as a cut
            it->second.second = K_MESH::Edge(IDX_NONE, IDX_NONE);
        }
      }
    }
  }

  if (wall)
    for (std::set<K_MESH::NO_Edge>::const_iterator i = wall->begin(); i != wall->end(); ++i) noe_to_bound.erase(*i);

  std::map<K_MESH::NO_Edge, std::pair<K_MESH::Edge, K_MESH::Edge> >::iterator itE = noe_to_bound.end();
  for (it = noe_to_bound.begin(); it != itE; ++it)
  {
    if (it->second.second.node(0) == IDX_NONE) // non-manifold
      continue;
    if (it->second.second.node(0) == -1)         // free-edge
      continue;

    const E_Int& pg1 = it->second.first.node(0);
    const E_Int& n1 = it->second.first.node(1);
    const E_Int& pg2 = it->second.second.node(0);
    const E_Int& n2 = it->second.second.node(1);

    neighbor.get_facet(pg1, n1) = pg2;
    neighbor.get_facet(pg2, n2) = pg1;
  }

  return 0;
}

///
E_Int Polygon::get_sharp_edges
(const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient, const ngon_unit& lneighbors,
 std::set<K_MESH::NO_Edge>& sharp_edges, E_Float angular_threshold, const E_Float** normals)
{
  // WARNING : node are 1-based uupon exit (in reflex_edge and convex_edges))

  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);

  E_Float ni[3], nj[3];

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int PGi = *(first_pg + i) - 1;
    const E_Int* pNi = PGS.get_facets_ptr(PGi);
    E_Int  nb_nodes = PGS.stride(PGi);
    const E_Int* pKn = lneighbors.get_facets_ptr(i);

    bool reversed = (orient[i] == -1);
    E_Int er = get_oriented_normal(crd, PGS, PGi, reversed, ni, normals);
    if (er) continue; // degen element

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int e0 = *(pNi + n);
      E_Int e1 = *(pNi + (n + 1) % nb_nodes);

      if (reversed)
        std::swap(e0, e1);

      const E_Float* E0 = crd.col(e0 - 1);
      const E_Float* E1 = crd.col(e1 - 1);

      E_Int j = *(pKn + n);
      if (j == IDX_NONE)
        continue;
      E_Int PGj = *(first_pg + j) - 1;
      

      //const E_Int* pNj = PGS.get_facets_ptr(PGj);
      //E_Int nb_nodsj = PGS.stride(PGj);

      bool rev = (orient[j] == -1);
      E_Int er = get_oriented_normal(crd, PGS, PGj, rev, nj, normals);
      if (er) continue; // degen element so we consider its edges as permeable

      // Concave or not ?
      E_Float alpha = NUGA::angle_measure(ni, nj, E0, E1);
      alpha = ::fabs(NUGA::PI - alpha);

      if (alpha >= angular_threshold)
        sharp_edges.insert(K_MESH::NO_Edge(e0, e1));
    }
  }

  return 0;
}

///
E_Int Polygon::update_neighbor_with_sharp_edges
(const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient, ngon_unit& lneighbors,
 E_Float angular_threshold, const E_Float** normals)
{
  // WARNING : node are 1-based uupon exit (in reflex_edge and convex_edges))

  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);

  E_Float ni[3], nj[3];

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int PGi = *(first_pg + i) - 1;
    const E_Int* pNi = PGS.get_facets_ptr(PGi);
    E_Int  nb_nodes = PGS.stride(PGi);
    E_Int* pKn = lneighbors.get_facets_ptr(i);

    bool reversed = (orient[i] == -1);
    E_Int er = get_oriented_normal(crd, PGS, PGi, reversed, ni, normals);
    if (er) continue; // degen element

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int e0 = *(pNi + n);
      E_Int e1 = *(pNi + (n + 1) % nb_nodes);

      if (reversed)
        std::swap(e0, e1);

      const E_Float* E0 = crd.col(e0 - 1);
      const E_Float* E1 = crd.col(e1 - 1);

      E_Int j = *(pKn + n);
      if (j == IDX_NONE)
        continue;
      E_Int PGj = *(first_pg + j) - 1;

      //const E_Int* pNj = PGS.get_facets_ptr(PGj);
      //E_Int nb_nodsj = PGS.stride(PGj);

      bool rev = (orient[j] == -1);
      E_Int er = get_oriented_normal(crd, PGS, PGj, rev, nj, normals);
      if (er) continue; // degen element so we consider its edges as permeable

      // Concave or not ?
      E_Float alpha = NUGA::angle_measure(ni, nj, E0, E1);
      alpha = ::fabs(NUGA::PI - alpha);

      if (alpha < angular_threshold)
        continue;

      *(pKn + n) = IDX_NONE; //set a cut in lneighbors
    }
  }

  return 0;
}


///
E_Int Polygon::get_oriented_normal(const K_FLD::FloatArray& crd, const ngon_unit& pgs, E_Int PGi/*glob id : sync with normals*/, bool reverse, E_Float* Normi, const E_Float** normals)
{
  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);

  const E_Int* pNi = pgs.get_facets_ptr(PGi);
  E_Int nb_nodes = pgs.stride(PGi);

  if (normals == 0)
    K_MESH::Polygon::normal<acrd_t, 3>(acrd, pNi, nb_nodes, 1, Normi);
  else
  {
    Normi[0] = normals[0][PGi];//fixme
    Normi[1] = normals[1][PGi];
    Normi[2] = normals[2][PGi];
  }

  if (reverse)
  {
    Normi[0] = -Normi[0];
    Normi[1] = -Normi[1];
    Normi[2] = -Normi[2];
  }

  E_Float l2 = sqrt(Normi[0] * Normi[0] + Normi[1] * Normi[1] + Normi[2] * Normi[2]);
  if (::fabs(l2 - 1.) < EPSILON) // NOT DEGEN
    return 0;
  return 1;
}

E_Int Polygon::full_agglomerate
(const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Float angular_tol, const E_Int* orient, ngon_unit& agg_pgs, std::vector<E_Int>& nids, const E_Float** normals)
{
  
  // Some face will be agglomerated, some not so a PG have 3 states : unchanged/deleted/created
  // so nids is sized a initial nb_pgs. But any id > nb_pgs refer to an agglomerate in agg_pgs.
  // nids[i] = i : UNCHANGED
  // nids[i] = IDX_NONE DELETED
  // nids[i] = i + nb_pgs : i-th PG in agg_pgs
  
//  typedef std::set<K_MESH::NO_Edge> eset_t;

  agg_pgs.clear();
  K_CONNECT::IdTool::init_inc(nids, nb_pgs);

  ngon_unit lneighbors;
  build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);

  // 
  update_neighbor_with_sharp_edges(crd, PGS, first_pg, nb_pgs, orient, lneighbors, angular_tol, normals);

  std::vector<E_Int> colors;
  NUGA::EltAlgo<K_MESH::Polygon>::coloring(lneighbors, colors);

  E_Int nb_connex = *std::max_element(colors.begin(), colors.end()) + 1;
  if (nb_connex == nb_pgs) // they are all unconnected polygons so nothing todo. The caller must check nids/agglo_pgs emptyness
    return 0;

  std::vector<E_Int> ids, ori;
  std::set<K_MESH::Edge> w_oe_set;//fixme hpc
  std::vector<std::deque<E_Int>> PG;
  std::map<E_Int, E_Int> w_n_map;//fixme hpc
  for (E_Int c = 0; c < nb_connex; ++c)
  {
    ids.clear();
    ori.clear();
    E_Int id0 = IDX_NONE;
    
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      if (colors[i] != c) continue;
        
      if (ids.empty()) id0 = i;// we store the first id for each color to hold the aggregate, other agglomerated are ignored since kept to IDX_NONE
        
      ids.push_back(*(first_pg+i)-1);
      ori.push_back(orient[i]);
    }
    
#ifdef DEBUG_AGGLOMERATOR
    assert (id0 > -1 && id0 < nids.size());
#endif
    
    if (ids.size() == 1) // i.e. unchanged
      continue;

    E_Int err = get_boundary(crd, PGS, ids, PG, ori, w_oe_set, w_n_map);
    if (err) continue;

    size_t sz = PG[0].size();
    PG[0].push_front(sz);
    agg_pgs.add(PG[0]);
    E_Int shft = agg_pgs.size()-1;

    // we set the first to hold the aggregate, other agglomerated are discarded by setting them to NONE
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      if (colors[i] != c)
        continue;
      if (i == id0) nids[i] = shft + nb_pgs;
      else nids[i] = IDX_NONE;
    }
  }
    
  return 0;
}

///
E_Int Polygon::shuffle_triangulation()
{
  
  if (_triangles == nullptr) return 0;
  
  E_Int ntris = nb_tris();
  if (ntris == 1) return 0;
  
  // Get inner edges
  std::vector<std::pair<E_Int, E_Int> > wpair_set;
  {
    std::set<K_MESH::NO_Edge> tmp;
    K_MESH::NO_Edge E;
    for (E_Int i=0; i < ntris; ++i)
    {
      for (size_t n=0; n < 3; ++n)
      {
        E.setNodes(_triangles[3*i + n], _triangles[3*i + (n+1)%3]);
        if (!tmp.insert(E).second)//already in => inner edge
          wpair_set.push_back(std::make_pair(i, (n+2)%3));
      }
    }
  }
  
  if (wpair_set.empty()) return 0; //Triangle
  
  // _triangle => IntArray (for fast_swap_edges)
  K_FLD::IntArray cT3(3, ntris);
  for (E_Int i=0; i < ntris; ++i)
    for (size_t n=0; n < 3; ++n) cT3(n,i) = _triangles[3*i + n];
  
  K_FLD::IntArray neighbors;
  NUGA::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(cT3, neighbors);
  NUGA::EltAlgo<K_MESH::Triangle>::fast_swap_edges(wpair_set, cT3, neighbors);
  
  // IntArray => _triangles
  for (E_Int i=0; i < ntris; ++i)
    for (size_t n=0; n < 3; ++n) _triangles[3*i + n] = cT3(n,i);
  
  return 0;
}

} //namespace
