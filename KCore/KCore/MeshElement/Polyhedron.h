/*    
    Copyright 2013-2016 Onera.

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
#ifndef __K_MESH_POLYHEDRON_H__
#define __K_MESH_POLYHEDRON_H__

#include "Fld/DynArray.h"
#include "Fld/ArrayAccessor.h"
#include "MeshElement/Polygon.h"
#include "MeshElement/Tetrahedron.h"
#include "Connect/GeomAlgo.h"
#include "Connect/BARSplitter.h"

#define dPATHO_PH_NONE 0
#define dCENTROID_NOT_STAR 2
#define dISO_BARY_NOT_STAR 3
#define dOPEN_PHS 4
#define dCONCAVITY_TO_SPLIT 5
#define dDELAUNAY_FAILURE 7

#ifdef DEBUG1_POLYHEDRON
#include "Nuga/Boolean/NGON_debug.h"
#endif

namespace K_MESH
{
//static int level=0;

//enum TopoShape { UNKNOWN, STAR_SHAPED, CONVEX, CONCAVE};
#define UNKNOWN 0
#define STAR_SHAPED 7

#define ZETOL 2.e-2

template <int TopoShape>
class Polyhedron
{   
public:
  const ngon_unit* _pgs;
  const E_Int* _faces;
  E_Int _nb_faces;
  E_Int* _triangles;
    
public:
  
  Polyhedron():_pgs(nullptr), _faces(nullptr), _triangles(nullptr){}
  
  Polyhedron(const ngon_unit* pgs, const E_Int* faces):_pgs(pgs), _faces(faces), _triangles(nullptr){}
  
  ~Polyhedron(){ if (_triangles != nullptr) delete [] _triangles;}
  
  ///
  template <typename CoordAcc>
  static inline void iso_barycenter(const CoordAcc& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, E_Float* G)
  {
    //
    for (size_t d = 0; d < 3; ++d) G[d] = 0.;

    for (E_Int i = 0; i < nb_nodes; ++i)
    {
      for (size_t d = 0; d < 3; ++d)
      {
        //std::cout << "v : " << coord.getVal(node(i), d) << std::endl;
        G[d] += coord.getVal(nodes[i]-index_start, d);
      }
    }

    E_Float k = 1. / (E_Float)nb_nodes;

    for (size_t i = 0; i < 3; ++i) G[i] *= k;
    //std::cout << "G : " << G[0] << "/" << G[1] << "/" << G[2] << std::endl;
  }
  
  static inline void iso_barycenter(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, E_Float* G)
  {
    //
    for (size_t d = 0; d < 3; ++d) G[d] = 0.;

    for (E_Int i = 0; i < nb_nodes; ++i)
    {
      const E_Float* p = crd.col(nodes[i]-index_start);
      for (size_t d = 0; d < 3; ++d) G[d] += p[d];
    }

    E_Float k = 1. / (E_Float)nb_nodes;

    for (size_t i = 0; i < 3; ++i) G[i] *= k;
    //std::cout << "G : " << G[0] << "/" << G[1] << "/" << G[2] << std::endl;
  }
  
  template<typename box_t, typename CoordAcc>
  void iso_barycenter(const CoordAcc& acrd, E_Float* G, E_Int index_start) const
  {
    for (size_t d = 0; d < 3; ++d) G[d] = 0.;
    E_Int counter = 0;
    
    for (E_Int i=0; i < _nb_faces; ++i)
    {
      E_Int PGi = *(_faces + i) - 1;
      const E_Int* nodes = _pgs->get_facets_ptr(PGi);
      E_Int nb_nodes =  _pgs->stride(PGi);
      
      for (E_Int j=0; j < nb_nodes; ++j)
      {
        ++counter;
        for (size_t d = 0; d < 3; ++d)
          G[d] += acrd.getVal(nodes[j]-index_start, d);
      }
    }
    
    E_Float k = 1. / (E_Float)counter;

    for (size_t i = 0; i < 3; ++i) G[i] *= k;
  }
  
  template<typename box_t, typename CoordAcc>
  void bbox(const CoordAcc& acrd, box_t&bb) const
  {
    for (E_Int i = 0; i < 3; ++i)
    {bb.minB[i] = K_CONST::E_MAX_FLOAT; bb.maxB[i] = -K_CONST::E_MAX_FLOAT;}
    
    box_t b;
    for (E_Int i=0; i < _nb_faces; ++i)
    {
      E_Int PGi = *(_faces + i) - 1;
      const E_Int* nodes = _pgs->get_facets_ptr(PGi);
      E_Int nb_nodes =  _pgs->stride(PGi);
      
      b.compute(acrd, nodes, nb_nodes, 1/*idx start*/);
      
      for (E_Int i = 0; i < 3; ++i)
      {
        bb.minB[i] = std::min(bb.minB[i], b.minB[i]);
        bb.maxB[i] = std::max(bb.maxB[i], b.maxB[i]);
      }
    }
  }
  
  static inline void reorient(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, ngon_unit& opgs)
  {
    // build PG neighborhood for that given PH
    opgs.clear();
    
    ngon_unit neighbors;
    std::vector<E_Int> oids, orient(nb_pgs, 1);
    PGS.extract(first_pg, nb_pgs, opgs, oids);
    K_MESH::Polygon::build_pg_neighborhood(opgs, neighbors); 
    K_CONNECT::EltAlgo<K_MESH::Polygon>::reversi_connex(opgs, neighbors, 0/*ref*/, orient);
    
    //Apply new orientation
    for (E_Int i = 0; i < opgs.size(); ++i)
    {
      if (orient[i] == -1)
      {
        E_Int s = opgs.stride(i);
        E_Int* p = opgs.get_facets_ptr(i);
        std::reverse(p, p + s);
      }
    }
  }
  
  
  static inline void reorient_PHT3(K_FLD::IntArray& connectT3)
  {
    //Dual graph
    K_FLD::IntArray neighbors;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectT3, neighbors);
  
    // Make facets orientation consistent
    std::vector<E_Int> orient;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::reversi_connex(connectT3, neighbors, 0/*Kseed*/, orient);
    //permut the third and second nodes
    for (size_t k = 0; k < orient.size(); ++k)
      if (orient[k] == -1)
        std::swap(connectT3(1,k), connectT3(2,k));
    return;
  }
  
  static inline void reorient_PHT3(K_FLD::IntArray& connectT3,  K_FLD::IntArray& neighbors)
  {
    //Dual graph
    neighbors.clear();
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectT3, neighbors);
  
    // Make facets orientation consistent
    std::vector<E_Int> orient;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::reversi_connex(connectT3, neighbors, 0/*Kseed*/, orient);
    //permut the third and second nodes
    for (size_t k = 0; k < orient.size(); ++k)
      if (orient[k] == -1)
      {
        std::swap(connectT3(1,k), connectT3(2,k));
        std::swap(neighbors(1,k), neighbors(2,k));
      }
    return;
  }

  /// Predicate on given element in a mesh
  static E_Int is_concave
    (const K_FLD::FloatArray& crd, const ngon_unit& PGS, 
     const E_Int* first_pg, E_Int nb_pgs, bool open, const E_Int* orient, 
     bool & concave, E_Float threshold = E_EPSILON, const E_Float** normals = 0)
  {
    concave = false;

    if (nb_pgs == 4) // TH4
      return 0;
    
    // angular criterion
    E_Float angle_threshold = K_CONST::E_PI*(1. - threshold);
    angle_threshold = std::min(K_CONST::E_PI, angle_threshold);
    angle_threshold = std::max(angle_threshold, E_EPSILON);
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
    
    // build PG neighborhood for that given PH
    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);
    
    E_Float ni[3], nj[3];

    for (E_Int i=0; i <nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      //std::cout << PGi << std::endl;
      const E_Int* pNi = PGS.get_facets_ptr(PGi);
      E_Int  nb_nodes = PGS.stride(PGi);
      E_Int* pKn = lneighbors.get_facets_ptr(i);//i because lneighbor is local : sized as nb_pgs

      bool reversed = (orient[i] == -1);
      E_Int er = Polygon::get_oriented_normal(crd, PGS, PGi, reversed, ni, normals);
      if (er) continue; // degen element

      for (E_Int n = 0; n < nb_nodes; ++n)
      {
        E_Int e0 = *(pNi + n);
        E_Int e1 = *(pNi + (n + 1) % nb_nodes);

        if (reversed)
          std::swap(e0, e1);

        const E_Float* E0 = crd.col(e0 - 1);
        const E_Float* E1 = crd.col(e1 - 1);

        E_Int j = *(pKn + n);//0-based and is between 0 and nb_pgs (referring to lpgs : but the same position in ng)

        //assert (PGj >= 0 && PGj < nb_pgs);
        if (j == E_IDX_NONE)
        {
          if (open) continue;
          return 1; // inconsistency : the PH is supposed closed
        }

        E_Int PGj = *(first_pg + j) - 1;

        //const E_Int* pNj = PGS.get_facets_ptr(PGj);
        //E_Int nb_nodsj = PGS.stride(PGj);

        bool rev = (orient[j] == -1);
        E_Int er = Polygon::get_oriented_normal(crd, PGS, PGj, rev, nj, normals);
        if (er) continue; // degen element

        // Concave or not ?
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);

        if (alpha < angle_threshold)
        {
          concave=true;
#ifdef DEBUG1_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGj);
          
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG1o.mesh", crd, PGS, PGi, ni);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG2o.mesh", crd, PGS, PGj, nj);
#endif
          return 0;
        }
      }
    }
  
    return 0;
  }
  
  /// 
  static E_Int min_max_angles
    (const K_FLD::FloatArray& crd, const ngon_unit& PGS, 
     const E_Int* first_pg, E_Int nb_pgs, bool open, const E_Int* orient, E_Float &minA, E_Float& maxA, const E_Float** normals = 0)
  {
        
    minA = 7.;// just a value more than 2Pi
    maxA = -1.;//
        
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
    
    // build PG neighborhood for that given PH
    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);
    
    E_Float ni[3], nj[3];

    for (E_Int i=0; i <nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      //std::cout << PGi << std::endl;
      const E_Int* pNi = PGS.get_facets_ptr(PGi);
      E_Int  nb_nodes = PGS.stride(PGi);
      E_Int* pKn = lneighbors.get_facets_ptr(i);//i because lneighbor is local : sized as nb_pgs

      bool reversed = (orient[i] == -1);
      E_Int er = Polygon::get_oriented_normal(crd, PGS, PGi, reversed, ni, normals);
      if (er) continue; // degen element

      for (E_Int n = 0; n < nb_nodes; ++n)
      {
        E_Int e0 = *(pNi + n);
        E_Int e1 = *(pNi + (n + 1) % nb_nodes);

        if (reversed)
          std::swap(e0, e1);

        const E_Float* E0 = crd.col(e0 - 1);
        const E_Float* E1 = crd.col(e1 - 1);

        E_Int j = *(pKn + n);//0-based and is between 0 and nb_pgs (referring to lpgs : but the same position in ng)

        //assert (PGj >= 0 && PGj < nb_pgs);
        if (j == E_IDX_NONE)
        {
          if (open) continue;
          return 1; // inconsistency : the PH is supposed closed
        }

        E_Int PGj = *(first_pg + j) - 1;

        //const E_Int* pNj = PGS.get_facets_ptr(PGj);
        //E_Int nb_nodsj = PGS.stride(PGj);

        bool rev = (orient[j] == -1);
        E_Int er = Polygon::get_oriented_normal(crd, PGS, PGj, rev, nj, normals);
        if (er) continue; // degen element

        // Concave or not ?
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);
        
#ifdef DEBUG1_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGj);
          
          std::ostringstream o;
          o << "PG_" << PGi << ".mesh";
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG(o.str().c_str(), crd, PGS, PGi, ni);
          o.str("");
          o << "PG_" << PGj << ".mesh";
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG(o.str().c_str(), crd, PGS, PGj, nj);
#endif
        
        if (alpha < minA)
        {
          minA = alpha;
#ifdef DEBUG1_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGj);
          
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG1o.mesh", crd, PGS, PGi, ni);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG2o.mesh", crd, PGS, PGj, nj);
#endif
        }
        if (alpha > maxA)
        {
          maxA = alpha;
        }
      }
    }
  
    return 0;
  }

  /// Predicate on given pair of neighbor element in a mesh
  static E_Int is_concave_pair
    (const K_FLD::FloatArray& crd, const ngon_unit& PGS, 
     const E_Int* pgsi, E_Int nb_pgsi, const E_Int* pgsj, E_Int nb_pgsj,
     const E_Int* orienti, const E_Int* orientj, 
     bool & concave, E_Float & commonSurf, E_Float concave_threshold = E_EPSILON, const E_Float** normals = 0)
  {
    concave = false;
    commonSurf = 0.;

    if (nb_pgsi + nb_pgsj == 8) // TH4-TH4
      return 0;
   
    // angular criterion
    E_Float angle_threshold = K_CONST::E_PI*(1. - concave_threshold);
    angle_threshold = std::min(K_CONST::E_PI, angle_threshold);
    angle_threshold = std::max(angle_threshold, E_EPSILON);

    // aggregate
    std::vector<E_Int> maski(nb_pgsi, 1), maskj(nb_pgsj, 1);
    {
      std::set<E_Int> tmp, commonPG;
      tmp.insert(pgsi, pgsi + nb_pgsi);
      for (E_Int i = 0; i < nb_pgsj; ++i)
      {
        if (tmp.find(*(pgsj + i)) != tmp.end())
        {
          maskj[i] = 0;
          commonPG.insert(*(pgsj + i));
        }
      }

      if (commonPG.empty()) // not adjacent !
        return 1;

      for (E_Int i = 0; i < nb_pgsi; ++i)
      {
        E_Int PGi = *(pgsi + i);
        if (commonPG.find(PGi) != commonPG.end())
        {
          E_Int nb_nodes = PGS.stride(PGi-1);
          const E_Int* nodes = PGS.get_facets_ptr(PGi-1);
          commonSurf += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, nodes, nb_nodes, 1);
          maski[i] = 0;
        }
      }
    }

    // build the aggregate and its neighboring 
    Vector_t<E_Int> lorient, oids;
    ngon_unit lpgs;

    for (E_Int i = 0; i < nb_pgsi; ++i)
    {
      if (maski[i] == 0)
        continue;
      E_Int PGi = *(pgsi + i) - 1;
      const E_Int* nodes = PGS.get_facets_ptr(PGi);
      E_Int nb_nodes = PGS.stride(PGi);
      lpgs.add(nb_nodes, nodes);
      lorient.push_back(orienti[i]);
      oids.push_back(PGi);
    }
    for (E_Int i = 0; i < nb_pgsj; ++i)
    {
      if (maskj[i] == 0)
        continue;
      E_Int PGi = *(pgsj + i) - 1;
      const E_Int* nodes = PGS.get_facets_ptr(PGi);
      E_Int nb_nodes = PGS.stride(PGi);
      lpgs.add(nb_nodes, nodes);
      lorient.push_back(orientj[i]);
      oids.push_back(PGi);
    }

    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(lpgs, lneighbors);

    E_Float ni[3], nj[3];
    E_Int nb_tot_pgs = lpgs.size();
    for (E_Int i = 0; i <nb_tot_pgs; ++i)
    {
      //std::cout << PGi << std::endl;
      E_Int* nodes = lpgs.get_facets_ptr(i);
      E_Int  nb_nodes = lpgs.stride(i);
      E_Int* pKn = lneighbors.get_facets_ptr(i);

      bool reversed = (lorient[i] == -1);
      E_Int er = K_MESH::Polygon::get_oriented_normal(crd, lpgs, i, reversed, ni, normals);
      if (er) continue; // degen element

      for (E_Int n = 0; n < nb_nodes; ++n)
      {
        E_Int e0 = *(nodes + n);
        E_Int e1 = *(nodes + (n + 1) % nb_nodes);

        if (reversed)
          std::swap(e0, e1);

        const E_Float* E0 = crd.col(e0 - 1);
        const E_Float* E1 = crd.col(e1 - 1);

        E_Int j = *(pKn + n);//0-based
        if (j == E_IDX_NONE)
          return 1;

        const E_Int* nodesj = lpgs.get_facets_ptr(j);
        E_Int nb_nodsj = lpgs.stride(j);

        bool rev = (lorient[j] == -1);
        E_Int er = K_MESH::Polygon::get_oriented_normal(crd, lpgs, j, rev, nj, normals);
        if (er) continue; // degen element

        // Concave or not ?
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);

        if (alpha < angle_threshold)
        {
          concave = true;
#ifdef DEBUG1_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, lpgs, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, lpgs, PGj);

          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG1o.mesh", crd, lpgs, PGi, ni);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG2o.mesh", crd, lpgs, PGj, nj);
#endif
          return 0;
        }
      }
    }

    return 0;
  }

  ///
  static E_Int is_concave
  (const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, bool open, const E_Int* orient, bool & concave, 
   std::map<K_MESH::NO_Edge, E_Float>& relex_edges, std::set<K_MESH::NO_Edge>& convex_edges,
   E_Float concave_threshold = E_EPSILON, E_Float convex_threshold = E_EPSILON,
   const E_Float** normals=0)
  {
    // WARNING : node are 1-based uupon exit (in reflex_edge and convex_edges))
    
    concave = false;

    if (nb_pgs == 4) // TH4
      return 0;
    
    // angular criterions
    E_Float angle_concave = K_CONST::E_PI*(1. - concave_threshold);
    angle_concave = std::min(K_CONST::E_PI, angle_concave);
    angle_concave = std::max(angle_concave, E_EPSILON);
    E_Float angle_convex = K_CONST::E_PI*(1. + convex_threshold);
    angle_convex = std::max(K_CONST::E_PI, angle_convex);
    angle_convex = std::min(angle_convex, 2.*K_CONST::E_PI - E_EPSILON);
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
  
    // build PG neighborhood for that given PH
    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);
    
    E_Float ni[3], nj[3];

    for (E_Int i=0; i <nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      const E_Int* pNi      = PGS.get_facets_ptr(PGi);
      E_Int  nb_nodes = PGS.stride(PGi);
      E_Int* pKn      = lneighbors.get_facets_ptr(i);
      
      bool reversed=(orient[i] == -1);
      E_Int er = Polygon::get_oriented_normal(crd, PGS, PGi, reversed, ni, normals);
      if (er) continue; // degen element
      
      for (E_Int n=0; n < nb_nodes; ++n)
      {
        E_Int e0 = *(pNi+n);
        E_Int e1 = *(pNi+(n+1)%nb_nodes);
        
        if (reversed)
          std::swap(e0,e1);
        
        const E_Float* E0 = crd.col(e0-1);
        const E_Float* E1 = crd.col(e1-1);

        E_Int j = *(pKn + n);
        if (j == E_IDX_NONE)
        {
          if (open) continue;
          return -1; // inconsitency : the PH is supposed closed
        }

        E_Int PGj = *(first_pg + j) - 1;
        
        //const E_Int* pNj = PGS.get_facets_ptr(PGj);
        //E_Int nb_nodsj = PGS.stride(PGj);

        bool rev=(orient[j] == -1);
        E_Int er = Polygon::get_oriented_normal(crd, PGS, PGj, rev, nj, normals);
        //if (er) continue; // degen element

        // Concave or not ?
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);
        
        // adding er test for concave but not convex is a hack : 
        // doesn't hurt to add an extra convex but it does hurt to add a concave edge as we might end up with a falsely-non-manifold discarded chain
        if (!er && alpha < angle_concave)
        {
          concave=true;
          relex_edges.insert(std::make_pair(K_MESH::NO_Edge(e0,e1), alpha));

#ifdef DEBUG1_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGj);
          
          std::ostringstream o;
          o << "PG_" << PGi << ".mesh";
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG(o.str().c_str(), crd, PGS, PGi, ni);
          o.str("");
          o << "PG_" << PGj << ".mesh";
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG(o.str().c_str(), crd, PGS, PGj, nj);
#endif
        }
        else if (alpha > angle_convex || er)
          convex_edges.insert(K_MESH::NO_Edge(e0, e1));
      }
    }
  
    return 0;
  }
  
  ///
  static E_Int is_concave
  (const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, bool open, const E_Int* orient, bool & concave, 
   E_Int& nb_reflex_edges,
   E_Float concave_threshold = E_EPSILON, E_Float convex_threshold = E_EPSILON,
   const E_Float** normals=0)
  {
    // WARNING : node are 1-based uupon exit (in reflex_edge and convex_edges))
    
    concave = false;
    nb_reflex_edges = 0;
    std::set<K_MESH::NO_Edge> reflex_edges;
   
    if (nb_pgs == 4) // TH4
      return 0;
    
    // angular criterions
    E_Float angle_concave = K_CONST::E_PI*(1. - concave_threshold);
    angle_concave = std::min(K_CONST::E_PI, angle_concave);
    angle_concave = std::max(angle_concave, E_EPSILON);
    E_Float angle_convex = K_CONST::E_PI*(1. + convex_threshold);
    angle_convex = std::max(K_CONST::E_PI, angle_convex);
    angle_convex = std::min(angle_convex, 2.*K_CONST::E_PI - E_EPSILON);
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
  
    // build PG neighborhood for that given PH
    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);
    
    E_Float ni[3], nj[3];

    for (E_Int i=0; i <nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      const E_Int* pNi      = PGS.get_facets_ptr(PGi);
      E_Int  nb_nodes = PGS.stride(PGi);
      E_Int* pKn      = lneighbors.get_facets_ptr(i);
      
      bool reversed=(orient[i] == -1);
      E_Int er = Polygon::get_oriented_normal(crd, PGS, PGi, reversed, ni, normals);
      if (er) continue; // degen element
      
      for (E_Int n=0; n < nb_nodes; ++n)
      {
        E_Int e0 = *(pNi+n);
        E_Int e1 = *(pNi+(n+1)%nb_nodes);
        
        if (reversed)
          std::swap(e0,e1);
        
        const E_Float* E0 = crd.col(e0-1);
        const E_Float* E1 = crd.col(e1-1);

        E_Int j = *(pKn + n);
        if (j == E_IDX_NONE)
        {
          if (open) continue;
          return -1; // inconsitency : the PH is supposed closed
        }

        E_Int PGj = *(first_pg + j) - 1;
        
        //const E_Int* pNj = PGS.get_facets_ptr(PGj);
        //E_Int nb_nodsj = PGS.stride(PGj);

        bool rev=(orient[j] == -1);
        E_Int er = Polygon::get_oriented_normal(crd, PGS, PGj, rev, nj, normals);
        //if (er) continue; // degen element

        // Concave or not ?
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);
        
        // adding er test for concave but not convex is a hack : 
        // doesn't hurt to add an extra convex but it does hurt to add a concave edge as we might end up with a falsely-non-manifold discarded chain
        if (!er && alpha < angle_concave)
        {
          concave=true;
          reflex_edges.insert(K_MESH::NO_Edge(e0,e1));

#ifdef DEBUG1_POLYHEDRON
//          K_FLD::FloatArray crdt(crd);
//          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGi);
//          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, PGS, PGj);
//          
//          std::ostringstream o;
//          o << "PG_" << PGi << ".mesh";
//          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG(o.str().c_str(), crd, PGS, PGi, ni);
//          o.str("");
//          o << "PG_" << PGj << ".mesh";
//          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG(o.str().c_str(), crd, PGS, PGj, nj);
#endif
        }
      }
    }
    
    nb_reflex_edges = reflex_edges.size();
    return 0;
  }
  
  #define SIGN(a) ((a < -E_EPSILON) ? -1 : ((a > E_EPSILON) ? 1 : 0))  
  
  /// is the polyhedron star_shaped regarding the input point ?
  static bool is_star_shaping
  (const E_Float* point, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient, E_Int& faultyFacet, const E_Float** normals = 0)
  {
    // Idea : the point must be bellow all the PGs pseudo-plans => approximate star-shaping test
    // WARNING : faultyFacet is a loc id
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
    
    E_Float Normi[3], P0Pt[3];
    E_Int sgn(0);
        
    //
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      const E_Int* pNi = PGS.get_facets_ptr(PGi);
      //E_Int nb_nodes = PGS.stride(PGi);
      
      bool reverse = orient[i] == -1;
      E_Int er = Polygon::get_oriented_normal(crd, PGS, PGi, reverse, Normi, normals);
      if (er) continue; // degen element

      E_Int p0 = *(pNi) - 1; // first PG point
      K_FUNC::diff<3>(point, crd.col(p0), P0Pt);
      
      E_Float d = K_FUNC::dot<3>(P0Pt, Normi);
      
      if (i==0)
        sgn = SIGN(d);
      
      E_Float s = SIGN(d);
      
      if (s*sgn != 1) //not all giving a consistent answer
      {
        faultyFacet=i; //loc id
        //std::cout << i << std::endl;
        return false;
      }
    }
    
    return true;
  }
  
  /// is the polyhedron star_shaped regarding the input point ?
  template <typename TriangulatorType>
  static bool is_star_shaping_PTH3
  (const E_Float* point, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Int& faultyFacet)
  {
    // Idea : the point must be bellow all the T3-plans => exact bu not on the real shape
    // WARNING : faultyFacet is a loc id
    
    ngon_unit opgs;
    reorient(PGS, first_pg, nb_pgs, opgs);
    
    TriangulatorType t;
    K_FLD::IntArray connectT3;
    E_Int err = triangulate(t, opgs, crd, connectT3, true, false); // PH -> PHT3
    if (connectT3.cols() == 0 || err)
    {
      std::cout << "could not triangulate properly" << std::endl;
      return true; // silent meshing failures to avoid to treat it as a concavity to split
    }

    //MIO::write("PHT3.mesh", crd, connectT3, "TRI");
    //
        
    E_Int nb_t3s = connectT3.cols();
    K_FLD::IntArray::const_iterator pK;
    E_Float P0Pt[3], Normi[3];
    E_Int sgn(0);
    //
    for (E_Int i=0; i < nb_t3s; ++i)
    {
      pK = connectT3.col(i); 
      
      K_MESH::Triangle::normal(crd, pK, Normi);
      E_Float l2 = ::sqrt(Normi[0]*Normi[0]+Normi[1]*Normi[1]+Normi[2]*Normi[2]);
    
      if (::fabs(l2 - 1.) >= E_EPSILON) continue;  // DEGEN : not a good quality triangulation

      E_Int p0 = *pK; // first T3 point
      K_FUNC::diff<3>(point, crd.col(p0), P0Pt);
      
      E_Float d = K_FUNC::dot<3>(P0Pt, Normi);
      
      if (i==0)
        sgn = SIGN(d);
      
      E_Float s = SIGN(d);
      
      if (s*sgn != 1) //not all giving a consistent answer
      {
        faultyFacet=i; //loc id
        //std::cout << i << std::endl;
        return false;
      }
    }
    
    return true;
  }
  
  //
  static E_Int get_star_shaping_pt
  (E_Float* point, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient, const E_Float** normals = 0)
  {
    // Idea : first guess : isobary, then move the point to try to make it staring.
    E_Float Normi[3], P0Pt[3];
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
    
    std::vector<E_Int> phnodes;
    unique_nodes(PGS, first_pg, nb_pgs, phnodes);
    iso_barycenter(acrd, &phnodes[0], phnodes.size(), 1, point);
    
    while (nb_pgs--) // number of attempts
    {
      E_Int i = E_IDX_NONE;//faulty face
      if (is_star_shaping(point, crd, PGS, first_pg, nb_pgs, orient, i, normals))
        return 0;
      
      E_Int PGi = *(first_pg+i)-1;
      const E_Int* pNi = PGS.get_facets_ptr(PGi);
      
      // move the point
      E_Int p0 = *(pNi) - 1; // first PG point
      K_FUNC::diff<3>(point, crd.col(p0), P0Pt);
      
      bool reverse = orient[i] == -1;
      E_Int er = Polygon::get_oriented_normal(crd, PGS, PGi, reverse, Normi, normals);
      if (er) continue; // degen element
      
      E_Float d = K_FUNC::dot<3>(P0Pt, Normi);
      assert(d >= -E_EPSILON); //is indeed faulty
      
      K_FUNC::sum<3>(1., point, -d-10*E_EPSILON/*-d*(1.01)*/, Normi, point); //fixme !!
      
    };
    return 1;
  }

// is the polyhedron star_shaped regarding its centroid ?
  template<typename TriangulatorType>
  static E_Int is_centroid_star_shaped
  (const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient, const E_Float** normals = 0)
  {
    E_Float centroid[3], v;
    E_Int err = metrics2<TriangulatorType>(t, crd, PGS, first_pg, nb_pgs, v, centroid, false/*not all cvx*/);
    if (err) return dDELAUNAY_FAILURE; // cannot tell as the triangulation failed
    
    E_Int faultyFacet;
    if (is_star_shaping(centroid, crd, PGS, first_pg, nb_pgs, orient, faultyFacet, normals))
      return 0; // OK : PATHO_PH_NONE

    return dCENTROID_NOT_STAR;  
  }
  
  // is the polyhedron star_shaped regarding its centroid ?
  template<typename TriangulatorType>
  static E_Int is_centroid_star_shaped_PTH3
  (const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Float** normals = 0)
  {
    //centroid is tested against PHT3 to be consistent as it is computed on the PHT3 rep
    
    E_Float centroid[3], V;
    E_Int err = metrics2(t, crd, PGS, first_pg, nb_pgs, V, centroid, false/*not all cvx*/);
    if (err) return dDELAUNAY_FAILURE; // cannot tell as the triangulation failed

    E_Int faultyFacet;
    if (is_star_shaping_PTH3<TriangulatorType>(centroid, crd, PGS, first_pg, nb_pgs, faultyFacet))
      return 0; // OK : PATHO_PH_NONE

    return dCENTROID_NOT_STAR;  
  }
  
  // is the polyhedron star_shaped regarding its centroid ?
  template<typename TriangulatorType>
  static E_Int is_iso_bary_star_shaped
  (const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient, const E_Float** normals = 0)
  {
    // isoG tar-shaping is teste against real PH (instaed of PHT3)
    
    E_Float isoG[3];
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
    acrd_t acrd(crd);
  
    std::vector<E_Int> phnodes;
    unique_nodes(PGS, first_pg, nb_pgs, phnodes);
    iso_barycenter(acrd, &phnodes[0], phnodes.size(), 1, isoG);

    E_Int faultyFacet;
    if (is_star_shaping(isoG, crd, PGS, first_pg, nb_pgs, orient, faultyFacet, normals))
      return 0; // OK : PATHO_PH_NONE

    return dISO_BARY_NOT_STAR; 
  }

  // is the polyhedron closed ?
  static bool is_closed(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
  {
    // build PG neighborhood for that given PH
    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);

    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      const E_Int* p = lneighbors.get_facets_ptr(i);
      E_Int stride = lneighbors.stride(i);

      for (E_Int j = 0; j < stride; ++j)
        if (*(p + j) == E_IDX_NONE) 
        {
#ifdef DEBUG1_POLYHEDRON
          const E_Int* nodes = PGS.get_facets_ptr(*(first_pg + i) - 1);
          std::cout << j+1 << "-th edge of PG (0-based) " << *(first_pg+i)-1 << " is free !" << std::endl;

          std::cout << "edge nodes (0-based) are : " << *(nodes + j)-1 << " and " << *(nodes + (j+1)%stride)-1 << std::endl;
#endif
          return false;
        }
    }

    return true;
  }
  
  template <typename Triangulator_t>
  static E_Int is_pathological(const Triangulator_t& dt,  const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient)
  {
    //CURRENT RULE :
      // If OPEN => BUG !
      // If convex : OK
      // else      :              
      //     if (PHT3-)star-shaped for at least one of 2 points (centroid, iso-bary) ==> OK  
      //     else : need a split
      
    
    bool concave;
    E_Int err = K_MESH::Polyhedron<UNKNOWN>::is_concave(crd, PGS, first_pg, nb_pgs, false /*i.e. logic for supposed close PHs*/, orient, concave, 0.1);
  
    if (err)
      return dOPEN_PHS;
    
    if (!concave) return dPATHO_PH_NONE;
    
    E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_iso_bary_star_shaped(dt, crd, PGS, first_pg, nb_pgs, orient);
    if (res == 0) return dPATHO_PH_NONE;
    assert (res == dISO_BARY_NOT_STAR);
    
    // DO WE ALSO ADD THE PHT3 VERSION FOR ISO BARY ?
      
    res = K_MESH::Polyhedron<UNKNOWN>::is_centroid_star_shaped(dt, crd, PGS, first_pg, nb_pgs, orient);
    if (res == 0) return dPATHO_PH_NONE;
    if (res == dDELAUNAY_FAILURE) return dDELAUNAY_FAILURE;
    assert (res == dCENTROID_NOT_STAR);
        
    res = K_MESH::Polyhedron<UNKNOWN>::is_centroid_star_shaped_PTH3(dt, crd, PGS, first_pg, nb_pgs);
    if (res == 0) return dPATHO_PH_NONE;
    if (res == dDELAUNAY_FAILURE) return dDELAUNAY_FAILURE;
    assert (res == dCENTROID_NOT_STAR);
        
    // here : concave and computable but all star test failed => need a split
    return dCONCAVITY_TO_SPLIT;
  }
  
  // more output : reflex edges
  template <typename Triangulator_t>
  static E_Int is_pathological(const Triangulator_t& dt,  const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient,
                               std::map<K_MESH::NO_Edge, E_Float>& reflex_edges, std::set<K_MESH::NO_Edge>& convex_edges,
                               E_Float concave_threshold = E_EPSILON, E_Float convex_threshold = E_EPSILON)
  {
    //CURRENT RULE :
      // If OPEN => BUG !
      // If convex : OK
      // else      :              
      //     if PHT3-star-shaped for at least one of 2 points (centroid, iso-bary) ==> OK  
      //     else : need a split
      
    
    bool concave;
    E_Int err = K_MESH::Polyhedron<UNKNOWN>::is_concave(crd,PGS, first_pg, nb_pgs, false, orient, concave, 
                                                        reflex_edges, convex_edges, concave_threshold, convex_threshold);
  
    if (err)
      return dOPEN_PHS;
    
    if (!concave) return dPATHO_PH_NONE;
    
    E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_iso_bary_star_shaped(dt, crd, PGS, first_pg, nb_pgs, orient);
    if (res == 0) return dPATHO_PH_NONE;
    assert (res == dISO_BARY_NOT_STAR);
    
    // DO WE ALSO ADD THE PHT3 VERSION FOR ISO BARY ?
      
    res = K_MESH::Polyhedron<UNKNOWN>::is_centroid_star_shaped(dt, crd, PGS, first_pg, nb_pgs, orient);
    if (res == 0) return dPATHO_PH_NONE;
    if (res == dDELAUNAY_FAILURE) return dDELAUNAY_FAILURE;
    assert (res == dCENTROID_NOT_STAR);
        
    res = K_MESH::Polyhedron<UNKNOWN>::is_centroid_star_shaped_PTH3(dt, crd, PGS, first_pg, nb_pgs);
    if (res == 0) return dPATHO_PH_NONE;
    if (res == dDELAUNAY_FAILURE) return dDELAUNAY_FAILURE;
    assert (res == dCENTROID_NOT_STAR);
        
    // here : concave and computable but all star test failed => need a split
    return dCONCAVITY_TO_SPLIT;
  }
  
  ///
  template<typename TriangulatorType>
  static E_Int metrics
  (const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Float&V, E_Float* G, E_Float* pStar = 0)
  {
    G[0]=G[1]=G[2]=0.;

    V=K_CONST::E_MAX_FLOAT;
    E_Float p4[3], v;
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);

    std::vector<E_Int> phnodes;
    bool is_th4=false;

    if (nb_pgs == 4) // TH4
    {
      unique_nodes(PGS, first_pg, nb_pgs, phnodes);
      const E_Int& nb_nodes = phnodes.size();
      if (nb_nodes == 4)
        is_th4=true;
    }
  
    if (is_th4)
    {
      const E_Float* p1 = crd.col(phnodes[0]-1);
      const E_Float* p2 = crd.col(phnodes[1]-1);
      const E_Float* p3 = crd.col(phnodes[2]-1);
      const E_Float* p4 = crd.col(phnodes[3]-1);

      V=::fabs(::K_MESH::Tetrahedron::volume(p1, p2, p3, p4));
  
      G[0] = 0.25 * (p1[0]+p2[0]+p3[0]+p4[0]);
      G[1] = 0.25 * (p1[1]+p2[1]+p3[1]+p4[1]);
      G[2] = 0.25 * (p1[2]+p2[2]+p3[2]+p4[2]);

      return 0;
    }
    else if (phnodes.empty())
      unique_nodes(PGS, first_pg, nb_pgs, phnodes);

    if (!pStar) //assume the iso bary is..bof bof
      K_MESH::Polyhedron<TopoShape>::iso_barycenter(acrd, &phnodes[0], phnodes.size(), 1, p4);
    else
    {
      p4[0]=pStar[0]; p4[1]=pStar[1]; p4[2]=pStar[2];
    }
    
    ngon_unit opgs;
    if (TopoShape != STAR_SHAPED)
      reorient(PGS, first_pg, nb_pgs, opgs);
    else
    {
      std::vector<E_Int> oids;
      PGS.extract(first_pg, nb_pgs, opgs, oids);
    }

    K_FLD::IntArray connectT3;
    E_Int err = triangulate(t, opgs, crd, connectT3, true, false); // PH -> PHT3
    if (connectT3.cols() == 0 || err)
    {
      std::cout << "could not triangulate properly" << std::endl;
      return 1;
    }

    //MIO::write("PHT3.mesh", crd, connectT3, "TRI");


    V=0.;
    for (E_Int i=0; i < connectT3.cols(); ++i)
    {
      const K_FLD::IntArray::const_iterator pS = connectT3.col(i);
      const E_Float* p1 = crd.col(*pS);
      const E_Float* p2 = crd.col(*(pS+1));
      const E_Float* p3 = crd.col(*(pS+2));

      v=::fabs(Tetrahedron::volume(p1,p2,p3,p4));//must be absolute value when not reorienting as for the star-shaped case

      V += v;
  
      G[0] += 0.25 * v * (p1[0]+p2[0]+p3[0]+p4[0]);
      G[1] += 0.25 * v * (p1[1]+p2[1]+p3[1]+p4[1]);
      G[2] += 0.25 * v * (p1[2]+p2[2]+p3[2]+p4[2]);
    }
  
    G[0] /= V;
    G[1] /= V;
    G[2] /= V;
    
    return 0;
  }

  /// is the polyhedron a prism ?
  static bool is_prism(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Int*& generators, E_Int*& HX6opposites)
  {
    // returns true iff it'a prism (including hexa)
    // in case of hexa, fill in HX6opposites which gives the opposite sides : i is opposed to HX6opposites[i] and vice and versa
    // in case of prism, fill in the 2 generators

    std::vector<E_Int> nb_nods(nb_pgs);
    //
    generators[0] = generators[1] = E_IDX_NONE;

    E_Int nb_quads = 0;
    E_Int nb_gens = 0;
    E_Int nb_gen_nodes = 0;
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      nb_nods[i] = PGS.stride(PGi);

      if (nb_nods[i] == 4)
        ++nb_quads;
      else
        ++nb_gens;

      if (nb_gens > 2) return false;
      if (nb_gen_nodes != 0 && nb_gen_nodes != nb_nods[i]) return false;

      nb_gen_nodes = nb_nods[i];
    }
    bool is_hexa = (nb_quads == nb_pgs) && (nb_pgs == 6);
    if ((nb_quads != nb_pgs - 2) && !is_hexa) return false;

    // at this point : either hexa, classical prism or unknown (e.g an hexa with a quad split into 2 T3)

    // build PG neighborhood for that given PH
    ngon_unit lneighbors;
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs);

    if (is_hexa) // compute each pair of opposite ids
    {
      for (E_Int i = 0; i < nb_pgs; ++i) HX6opposites[i] = E_IDX_NONE;
      //
      const E_Int* neigh = lneighbors.get_facets_ptr(0);

      for (E_Int n = 0; n < 4; ++n)
      {
        E_Int j = *(neigh + n);
        if (j == E_IDX_NONE) return false; //open cells !
        E_Int jopp = *(neigh + (n + 2) % 4);
        if (jopp == E_IDX_NONE) return false; //open cells !

        HX6opposites[j] = jopp;
        HX6opposites[jopp] = j;
      }
      //O is opposite to ?
      for (E_Int i = 1; i < nb_pgs; ++i)
        if (HX6opposites[i] == E_IDX_NONE){
          HX6opposites[0] = HX6opposites[i];
          HX6opposites[i] = HX6opposites[0];
        break;
      }

      return true;
    }
    else // CLASSIC PRISM ?
    {
      //find generators
      for (E_Int p = 0; p < nb_pgs; ++p)
      {
        if (nb_nods[p] != 4)
        {
          E_Int k = (generators[0] == E_IDX_NONE) ? 0 : 1;
          generators[k] = p;
        }
      }
        
      // check that all the generator neighbors are quads (check only for one of them)
      const E_Int& first_gen = generators[0];
      const E_Int* neigh = lneighbors.get_facets_ptr(first_gen);
      for (E_Int n = 0; n < nb_nods[first_gen]; ++n)
      {
        E_Int j = *(neigh + n);
        if (j == E_IDX_NONE) return false; //open cells !
        if (nb_nods[j] != 4) return false;
      }

      return true;
    }
  }

  ///
  static bool has_baffles
  (const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, std::map<K_MESH::NO_Edge, E_Int>& w_map, E_Int* buffer_baffle_ids)
  {
    w_map.clear();
    
    std::map<K_MESH::NO_Edge, E_Int>::iterator it;
    K_MESH::NO_Edge E;
 
    // Count for each eadge how many times it is referenced
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      const E_Int* nodes = PGS.get_facets_ptr(PGi);
      E_Int nb_nodes = PGS.stride(PGi);
      
      for (E_Int n=0; n < nb_nodes; ++n)
      {
        E.setNodes(*(nodes + n), *(nodes + (n + 1) % nb_nodes));
        it = w_map.find(E);
        if (it == w_map.end()) w_map[E] = 1;
        else ++it->second;
      }
    }
    
    //
    bool has_bfl=false;
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      const E_Int* nodes = PGS.get_facets_ptr(PGi);
      E_Int nb_nodes = PGS.stride(PGi);
      
      buffer_baffle_ids[i]=i;
      
      for (E_Int n=0; n < nb_nodes; ++n)
      {
        E.setNodes(*(nodes + n), *(nodes + (n + 1) % nb_nodes));
        it = w_map.find(E);
        
        if (it->second > 1) continue; //manifold or not
        
        //baffle !
        buffer_baffle_ids[i]=E_IDX_NONE; //set to NONE (for removing) any baffle
        has_bfl=true;
        break;
      }
    }
    
    return has_bfl;
  }

  ///
  template<typename TriangulatorType>
  static E_Int metrics2
  (const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Float&V, E_Float* G, bool all_pgs_are_cvx)
  {
    G[0]=G[1]=G[2]=0.;

    V=K_CONST::E_MAX_FLOAT;
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);

    bool is_th4=false;

    if (nb_pgs == 4) // TH4
    {
      std::vector<E_Int> phnodes;
      unique_nodes(PGS, first_pg, nb_pgs, phnodes);
      const E_Int& nb_nodes = phnodes.size();
      if (nb_nodes == 4)
        is_th4=true;
    
      if (is_th4)
      {
        const E_Float* p1 = crd.col(phnodes[0]-1);
        const E_Float* p2 = crd.col(phnodes[1]-1);
        const E_Float* p3 = crd.col(phnodes[2]-1);
        const E_Float* p4 = crd.col(phnodes[3]-1);

        V=::fabs(::K_MESH::Tetrahedron::volume(p1, p2, p3, p4));
  
        G[0] = 0.25 * (p1[0]+p2[0]+p3[0]+p4[0]);
        G[1] = 0.25 * (p1[1]+p2[1]+p3[1]+p4[1]);
        G[2] = 0.25 * (p1[2]+p2[2]+p3[2]+p4[2]);

        return 0;
      }
    }
    
    ngon_unit opgs;
    reorient(PGS, first_pg, nb_pgs, opgs);
    
    K_FLD::IntArray connectT3;

    E_Int err(0);

    if (all_pgs_are_cvx)
    {
      typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
      acrd_t acrd(crd);
              
      for (E_Int i = 0; (i<nb_pgs); ++i)
      {
        const E_Int* nodes = opgs.get_facets_ptr(i);
        E_Int nb_nodes = opgs.stride(i);
        E_Float normal[3];
        K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, normal);
        E_Int iworst, ibest;
        K_MESH::Polygon:: is_convex(crd, nodes, nb_nodes, 1, normal, 1.e-8/*convexity_tol*/, iworst, ibest);
        K_MESH::Polygon::cvx_triangulate(crd, nodes, nb_nodes, ibest, 1, connectT3);
      }
      
      //MIO::write("PHT3.mesh", crd, connectT3, "TRI");    
    }
    else
    {
      err = triangulate(t, opgs, crd, connectT3, false, false); // PH -> PHT3
      //MIO::write("PHT30.mesh", crd, connectT3, "TRI");
    }
    
    err &= (connectT3.cols() != 0);

    if (!err) metrics(crd, connectT3, V/*acc*/, G/*Cacc*/);

    return err;
  }

  ///
  template <typename TriangulatorType>
  static E_Int triangulate
    (const TriangulatorType& dt, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const K_FLD::FloatArray& crd, K_FLD::IntArray& connectT3, std::vector<E_Int>& colors, bool do_not_shuffle, bool improve_quality)
  {
    connectT3.clear();
    colors.clear();
    PGS.updateFacets();

    E_Int err(0);
    for (E_Int i = 0; (i<nb_pgs) && !err; ++i)
    {
      E_Int PGi = *(first_pg + i) - 1;
      err = K_MESH::Polygon::triangulate(dt, crd, PGS.get_facets_ptr(PGi), PGS.stride(PGi), 1/*index start*/, connectT3, do_not_shuffle, improve_quality);
      colors.resize(connectT3.cols(), i); // assign color i to the appendaned T3s
    }
    return err;
  }
  
  
   ///
  template <typename TriangulatorType>
  static E_Int triangulate
    (const TriangulatorType& dt, const ngon_unit& lpgs, const K_FLD::FloatArray& crd, K_FLD::IntArray& connectT3, bool do_not_shuffle, bool improve_quality)
  {
    connectT3.clear();
    lpgs.updateFacets();
    
  
    E_Int nb_pgs = lpgs.size();

    E_Int err(0);
    for (E_Int i = 0; (i<nb_pgs) && !err; ++i)
    {
      err = K_MESH::Polygon::triangulate(dt, crd, lpgs.get_facets_ptr(i), lpgs.stride(i), 1/*index start*/, connectT3, do_not_shuffle, improve_quality);
    }
    return err;
  }
  
  E_Int nb_tris()
  {
    E_Int ntris=0;
    for (E_Int i=0; i < _nb_faces; ++i)
    {
      E_Int PGi = *(_faces + i) - 1;
      ntris += (_pgs->stride(PGi)-2);
    }
    return ntris;
  }
  
   ///
  template <typename TriangulatorType, typename acrd_t>
  E_Int triangulate
    (const TriangulatorType& dt, const acrd_t& acrd)
  {
    _pgs->updateFacets();

    E_Int ntris = nb_tris();
    _triangles = new E_Int[ntris*3];
    
    //to get local
    K_FLD::FloatArray crd;
    std::vector<E_Int> oids, lpgs;

    E_Int err(0);
    E_Float Pt[3];
    E_Int * pstart = _triangles;
    for (E_Int i = 0; (i<_nb_faces) && !err; ++i)
    {
      E_Int PGi = *(_faces + i) - 1;
      E_Int stride = _pgs->stride(PGi);
      const E_Int* nodes = _pgs->get_facets_ptr(PGi);
      
      crd.clear();
      lpgs.clear();
      lpgs.resize(stride);
      oids.resize(stride);
      for (E_Int n=0; n < stride; ++n){
        lpgs[n]=n; oids[n]=nodes[n]-1;
        acrd.getEntry(nodes[n]-1, Pt);
        crd.pushBack(Pt, Pt+3);
      }
      
      //std::cout << crd  << std::endl;
      //MIO::write("crd.mesh", crd, K_FLD::IntArray(), "TRI");
      
      err = K_MESH::Polygon::triangulate_inplace(dt, crd, &lpgs[0], stride, 0/*index start*/, pstart, true/*do_not_shuffle*/, false/*improve_quality*/);
      E_Int pgntris = stride - 2;
      for (E_Int k=0; k < pgntris * 3; ++k, ++pstart) *pstart = oids[*pstart]; //get back to global but starting at 0
    }
    
    return err;
  }
  
  inline void triangle(E_Int i, E_Int* target)
  {
    assert (_triangles != nullptr);
    const E_Int* p = &_triangles[i*3];
    target[0] = *(p++);
    target[1] = *(p++);
    target[2] = *p;
    
  }
  
///
  static E_Int merge_two_phs
    (const K_FLD::FloatArray& crd, const ngon_unit& PGS, 
    const E_Int* pgs1, E_Int nb_pgs1, const E_Int* pgs2, E_Int nb_pgs2, 
    const E_Int* orient1, const E_Int* orient2, ngon_unit& mPH, ngon_unit& oorient)
{
  
  std::set<E_Int> tmp, commonPG;
  tmp.insert(pgs1, pgs1+nb_pgs1);
  
  for (E_Int i=0; i < nb_pgs2; ++i)
    if (tmp.find(*(pgs2+i)) != tmp.end())
      commonPG.insert(*(pgs2+i));
  
  std::vector<E_Int> molecPH, molecO;
  E_Int sz = nb_pgs1 + nb_pgs2 - 2*commonPG.size();
  molecPH.reserve(sz);
  molecO.reserve(sz);
  
  for (E_Int i=0; i < nb_pgs1; ++i)
    if (commonPG.find(*(pgs1+i)) == commonPG.end())
    {
      molecPH.push_back(*(pgs1+i));
      molecO.push_back(*(orient1+i));
    }
  for (E_Int i=0; i < nb_pgs2; ++i)
    if (commonPG.find(*(pgs2+i)) == commonPG.end())
    {
      molecPH.push_back(*(pgs2+i));
      molecO.push_back(*(orient2+i));
    }
  
  mPH.add(molecPH.size(), &molecPH[0]);
  mPH._type.push_back(INNER);//fixme
  oorient.add(molecO.size(), &molecO[0]);
  
  return 0;
}

 ///
  static E_Int merge_two_phs
  (const K_FLD::FloatArray& crd, const ngon_unit& PGS,
   const E_Int* pgs1, E_Int nb_pgs1, const E_Int* pgs2, E_Int nb_pgs2,
   ngon_unit& mPH)
  {
    mPH.clear();
    if (nb_pgs2*nb_pgs1 == 0) return 1;

    std::set<E_Int> tmp, commonPG;
    tmp.insert(pgs1, pgs1+nb_pgs1);
  
    for (E_Int i=0; i < nb_pgs2; ++i)
      if (tmp.find(*(pgs2+i)) != tmp.end())
        commonPG.insert(*(pgs2+i));
  
    std::vector<E_Int> molecPH;
    
    if (commonPG.size() == nb_pgs1 && nb_pgs1 == nb_pgs2) // PH1 is the same as PH2 => by convention the merge is one of them
    { 
      molecPH.reserve(nb_pgs1);
      for (E_Int i=0; i < nb_pgs1; ++i)
        molecPH.push_back(*(pgs1+i));
    }
    else
    {
      E_Int sz = nb_pgs1 + nb_pgs2 - 2*commonPG.size();
      molecPH.reserve(sz);
  
      for (E_Int i=0; i < nb_pgs1; ++i)
        if (commonPG.find(*(pgs1+i)) == commonPG.end())
          molecPH.push_back(*(pgs1+i));
      for (E_Int i=0; i < nb_pgs2; ++i)
        if (commonPG.find(*(pgs2+i)) == commonPG.end())
          molecPH.push_back(*(pgs2+i));
    }

    mPH.add(molecPH.size(), &molecPH[0]);
    mPH._type.push_back(INNER);//fixme

  return 0;
}

static void unique_nodes(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, std::vector<E_Int>& unodes)
{
  unodes.clear();
  std::set<E_Int> tmp;

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int pgi = *(first_pg+i) - 1;
    E_Int nb_nodes = PGS.stride(pgi);
    const E_Int* nodes = PGS.get_facets_ptr(pgi);

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int N = *(nodes + n);
      if (tmp.insert(N).second) //not already in
        unodes.push_back(N);
    }
  }
}

static E_Int unique_nodes(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
{
  std::set<E_Int> tmp;
  E_Int nb_unodes(0);

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int pgi = *(first_pg+i) - 1;
    E_Int nb_nodes = PGS.stride(pgi);
    const E_Int* nodes = PGS.get_facets_ptr(pgi);

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int N = *(nodes + n);
      if (tmp.insert(N).second) //not already in
        ++nb_unodes;
    }
  }
  return nb_unodes;
}

inline static void expressions (E_Float w0 , E_Float w1 , E_Float w2 , E_Float &f1 , E_Float& f2 , E_Float& f3 , E_Float& g0 , E_Float& g1 , E_Float& g2 )
{
  E_Float temp0 = w0 + w1 ;   f1 = temp0 + w2 ;
  E_Float temp1 = w0 * w0 ;   
  E_Float temp2 = temp1 + w1 * temp0 ;
 
  f2 = temp2 + w2 * f1 ;   
  f3 = w0 * temp1 + w1 * temp2 + w2 * f2 ;
  g0 = f2 + w0 * ( f1 + w0 ) ;
  g1 = f2 + w1 * ( f1 + w1 ) ;
  g2 = f2 + w2 * ( f1 + w2 ) ;
}
  
inline static void metrics(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cT3, E_Float & volume, E_Float* centroid)
{
  const E_Float  mult[10] = {1./6. ,1./24. ,1./24. ,1./24. ,1./60. ,1./60. ,1./60. ,1./120. ,1./120. ,1./120.};
  E_Float  intg[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};//  order :  1 ,  x ,  y ,  z ,  x ^2 ,  y ^2 ,  z ^2 ,  xy ,  yz ,  zx
  
  K_FLD::IntArray::const_iterator pK;
  E_Int i0, i1, i2;
  E_Float f1x, f2x, f3x, g0x, g1x, g2x, f1y, f2y, f3y, g0y, g1y, g2y, f1z, f2z, f3z, g0z, g1z, g2z;
  
  for (E_Int i=0; i < cT3.cols(); ++i)
  {
    pK = cT3.col(i);
    i0=*(pK); i1 = *(pK+1); i2 = *(pK+2);
    
    const E_Float *p0 = crd.col(i0);
    const E_Float *p1 = crd.col(i1);
    const E_Float *p2 = crd.col(i2);
    
    E_Float a1 = p1[0] - p0[0];
    E_Float b1 = p1[1] - p0[1];
    E_Float c1 = p1[2] - p0[2];
    
    E_Float a2 = p2[0] - p0[0];
    E_Float b2 = p2[1] - p0[1];
    E_Float c2 = p2[2] - p0[2];
    
    E_Float d0 = b1 * c2 - b2 * c1;
    E_Float d1 = a2 * c1 - a1 * c2;
    E_Float d2 = a1 * b2 - a2 * b1;
    
    expressions(p0[0],p1[0],p2[0], f1x, f2x, f3x, g0x, g1x, g2x);
    expressions(p0[1],p1[1],p2[1], f1y, f2y, f3y, g0y, g1y, g2y);
    expressions(p0[2],p1[2],p2[2], f1z, f2z, f3z, g0z, g1z, g2z);
    
    intg[0] += d0*f1x;
    intg[1] += d0*f2x;
    intg[2] += d1*f2y;
    intg[3] += d2*f2z;
    intg[4] += d0*f3x;
    intg[5] += d1*f3y;
    intg[6] += d2*f3z;
    intg[7] += d0 * (p0[1] * g0x + p1[1] * g1x + p2[1] * g2x);
    intg[8] += d1 * (p0[2] * g0y + p1[2] * g1y + p2[2] * g2y);
    intg[9] += d2 * (p0[0] * g0z + p1[0] * g1z + p2[0] * g2z);;
  }
  
  for (E_Int i=0; i < 10; ++i) intg[i] *= mult[i];
  
  volume = intg[0];
  
  E_Float k=1./volume;
  
  centroid[0] = k * intg[1];
  centroid[1] = k * intg[2];
  centroid[2] = k * intg[3];
  
}

//typedef E_Float (*pFunc)(E_Float x, E_Float y, E_Float z);
//
//template<typename TriangulatorType>
//static void integ_O3(const K_FLD::FloatArray& crd, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, pFunc func, E_Float& val)
//{
//  val = 0.;
//  
//  E_Float FG[3], CG[3], vol;
//  TriangulatorType dt;
//  
//
//  K_MESH::Polyhedron<STAR_SHAPED>::metrics2<TriangulatorType>(dt, crd, PGS, first_pg, nb_pgs, vol, CG);
//  
//  const E_Float* p3 = CG;
//
//  for (size_t f = 0; f < nb_pgs; ++f)
//  {
//    E_Int fi = *(first_pg+f) - 1;
//    const E_Int* nodes = PGS.get_facets_ptr(fi);
//    E_Int nb_nodes = PGS.stride(fi);
//    
//    K_MESH::Polygon::centroid<3>(crd, nodes, nb_nodes, 1, FG);
//    
//     const E_Float* p0 = FG;
//    
//    for (E_Int n=0; n < nb_nodes; ++n)
//    {
//      E_Float v;
//      const E_Float* p1 = crd.col(*(nodes+n) - 1);
//      const E_Float* p2 = crd.col(*(nodes+(n+1)%nb_nodes) - 1);
//    
//      K_MESH::Tetrahedron::integ_O3(p0, p1, p2, p3, func, v);
//     
//      val +=v;
//    }
//  } 
//}

static void get_orient(const ngon_unit& PHs, E_Int PHi, const K_FLD::IntArray& F2E, E_Int* PHi_orient)
{    
  E_Int nb_pgs = PHs.stride(PHi);
  const E_Int* p = PHs.get_facets_ptr(PHi);
  
  for (int i = 0; i < nb_pgs; ++i)
    PHi_orient[i] = (F2E(1,p[i]-1) == PHi) ? -1 : 1;
}

static bool pt_is_inside(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const K_FLD::FloatArray& crd, const E_Int* PHi_orient, const E_Float* pt, E_Float tolerance)
{
  
  for (int i = 0; i < nb_pgs; i++)
  {
    const E_Int* pN = PGS.get_facets_ptr(first_pg[i]-1);
    E_Float det = K_FUNC::zzdet4(crd.col(pN[0]-1), crd.col(pN[1]-1), crd.col(pN[2]-1), pt); // approx : consider the first 3 nodes of the PG for defining the plane
        
    E_Int s = zSIGN(det,ZERO_M);
        
    if ( s == 0 ) // if det is in the random area, calculate the distance h
    {
      E_Float u[3];
      E_Float v[3];
      for (int j = 0; j < 3; ++j)
      {
        u[j] = crd.col(pN[1]-1)[j] - crd.col(pN[0]-1)[j];
        v[j] = crd.col(pN[2]-1)[j] - crd.col(pN[0]-1)[j];
      }
        
      E_Float norm = ::sqrt(K_FUNC::sqrCross<3>(u,v));
      E_Float h = det / norm;
        
      s = zSIGN(h,tolerance);
    }
    // position will be known comparing the sign to the orientation of the PG
    if ((PHi_orient[i] == 1) && (s == 1) ) return false;
    else if ((PHi_orient[i] == -1) && (s == - 1) ) return false;      

  }
  return true;
}

static bool is_HX8(const ngon_unit& PGs, const E_Int* firstPG, E_Int nb_pgs)
{    
  if (nb_pgs != 6) return false;
  
  for (int i=0; i<6; i++)
    if (PGs.stride(*(firstPG+i)-1) != 4)
      return false;

  return true;
}

static bool is_TH4(const ngon_unit& PGs, const E_Int* firstPG, E_Int nb_pgs)
{    
  if (nb_pgs != 4) return false;
  
  for (int i=0; i<4; i++)
    if (PGs.stride(*(firstPG+i)-1) != 3) return false; 

  return true;
}

static bool is_PY5(const ngon_unit& PGs, const E_Int* firstPG, E_Int nb_pgs)
{    
  E_Int s1(0), s2(0); 
    
  if (nb_pgs != 5) return false;
  
  for (int i=0; i<5; i++)
  {
    if (PGs.stride(*(firstPG+i)-1) == 3) ++s1;
    else if (PGs.stride(*(firstPG+i)-1) == 4) ++s2;
    else return false;
  }

  return ((s1 == 4) && (s2 == 1));
}

static bool is_PR6(const ngon_unit& PGs, const E_Int* firstPG, E_Int nb_pgs)
{    
  E_Int s1(0), s2(0);

  if (nb_pgs != 6) return false;
    
  for (int i=0; i<6; i++)
  {
    if (PGs.stride(*(firstPG+i)-1) == 3) ++s1;
    else if (PGs.stride(*(firstPG+i)-1) == 4) ++s2;
    else return false;
  }
  
  return ((s1==2) && (s2==3));
 
}

};


}
#endif	/* __K_MESH_HEXAHEDRON_H__ */
