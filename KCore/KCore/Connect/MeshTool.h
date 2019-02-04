/*    
    Copyright 2013-2019 Onera.

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

#ifndef K_CONNECT_MESHTOOL_H_
#define K_CONNECT_MESHTOOL_H_

#include "Search/KdTree.h"
#include "Def/DefContainers.h"
#include "MeshElement/Triangle.h"
#include "MeshElement/Polygon.h"
#include "Fld/ngon_unit.h"

#include <deque>

#define CONVEXITY_TOL 1.e-8 //temporoary to be consistent with BooleanNgon

namespace K_CONNECT{

class MeshTool
{
  public:

    typedef K_CONT_DEF::size_type         size_type;
    typedef K_CONT_DEF::int_vector_type   int_vector_type;
    typedef K_CONT_DEF::int_set_type      int_set_type;
    typedef K_CONT_DEF::int_pair_type     int_pair_type;
    typedef K_CONT_DEF::int_pair_set_type int_pair_set_type;
    typedef K_CONT_DEF::int_pair_set_type boundary_set_type;
    typedef K_CONT_DEF::bool_vector_type bool_vector_type;
    typedef K_MESH::Triangle              element_type;
    typedef K_SEARCH::KdTree<>            tree_type;
    
    typedef std::map<E_Int, int_set_type >          id_to_ids_t;
    typedef std::vector<std::deque<E_Int> >         idlists_t;
 
  public:

    ///
    MeshTool(const tree_type& tree, E_Float tolerance = E_EPSILON);
    ///
    ~MeshTool(void);

    /// Returns the elements containing the input point.
    template <typename OutputIterator>
    E_Int getContainingElements
      (const E_Float* point, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
       const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, const bool_vector_type& mask, 
       const tree_type& tree, OutputIterator out, size_type & N0) const;

    template <typename MetricType>
    void getIntersectionNodes(size_type Ni, size_type Nj, int_vector_type& Qis, std::vector<MetricType>& Mis){/*fixme*/Qis.push_back(Ni);Qis.push_back(Nj);}

    template <typename OutputIterator>
    void getAncestors(size_type N, const K_FLD::IntArray& connectM, const int_vector_type& ancestors, const K_FLD::IntArray& neighbors, OutputIterator out) const;
/*
    template <typename KPredicate, typename BPredicate>
    void getConnexSet(size_type Kseed, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
                                 int_set_type& set, boundary_set_type& bset,
                                 const KPredicate& kpred, const BPredicate& bpred) const;
*/
    template <typename KPredicate,typename BPredicate>
    void getConnexSet1(size_type Kseed, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
                                 int_set_type& set, boundary_set_type& bset,
                                 const KPredicate& kpred, const BPredicate& bpred) const;
/*
    template <typename KPredicate,typename BPredicate>
    void getConnexSet2(size_type Kseed, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors,
                                 int_set_type& set, boundary_set_type& bset,
                                 const KPredicate& kpred, const BPredicate& bpred) const;
*/
    static void getBoundary(const K_FLD::IntArray& connect, K_FLD::IntArray& connectB);
    static void getBoundaryT3Mesh(const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, K_CONT_DEF::int_pair_vector_type& boundaries);
    static void getBoundaryT3Mesh(const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, K_FLD::IntArray& cB);
    static void getNonManifoldEdges(const K_FLD::IntArray& connect, K_FLD::IntArray& connectB);

    static void getAttachedElements(const K_FLD::IntArray& connectB, const K_FLD::IntArray& connectS, 
                                    K_FLD::IntArray& connectElts);
    
    static void build_node_arity(const std::set<K_MESH::NO_Edge>& edges, std::map<E_Int, E_Int>& node_to_count);

    static void computeNodeNormals(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& normals);
    static E_Int computeNodeNormals(const K_FLD::FloatArray& crd, const ngon_unit& pgs, K_FLD::FloatArray& normals, E_Int smooth_iters = 0);
    
    static E_Int computeNodeNormalsFromPGNormals (const K_FLD::FloatArray& PG_normals, const ngon_unit& pgs, K_FLD::FloatArray& node_normals);
    
    static E_Int smoothNodeNormals(const ngon_unit& pgs, K_FLD::FloatArray& normals, E_Int smooth_iters=1);

    static void compact_to_mesh(K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& new_IDs, E_Int* N0 = 0);
    
    static void compact_to_mesh(const K_FLD::FloatArray& pos0, const K_FLD::IntArray& connect0,
                                K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
                                std::vector<E_Int>& oids);

    /// L stores the min, and the max length for each element (column-wise). The first row is for mins.
    template <E_Int DIM>
    static void computeEdgesSqrLengths(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& L);

    ///
    template <E_Int DIM>
    static void computeMinMaxEdgeSqrLength(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Float& min_d, E_Float& max_d);

    /// L stores the min, and the max incident edge's length for each node (column-wise). The first row is for mins.
    template <E_Int DIM>
    static void computeIncidentEdgesSqrLengths(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& L);
    
    static void computeIncidentEdgesSqrLengths(const K_FLD::FloatArray& crd, const ngon_unit& pgs, K_FLD::FloatArray& L);
    
    static void computeMinMaxIndices (const K_FLD::IntArray& connect, E_Int& min_i, E_Int& max_i);

    ///
    static void boundingBox (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2, 
                             E_Float* minB, E_Float* maxB);
    ///
    static void boundingBox (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, 
                             E_Int* minB, E_Int* maxB);

    static void flipT3(K_FLD::IntArray& connectT3);

    static void metricAnisoE22D(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2, double kn, K_FLD::FloatArray& metric);
    ///
    static bool detectDuplicated(const K_FLD::IntArray& connect, std::vector<E_Int>& dupIds, bool strict_orient = true);
    ///
    static E_Int removeDuplicated(K_FLD::IntArray& connect, std::vector<E_Int>& dupIds, bool strict_orient = true);

    ///
    template <typename EdgeType>
    static void extractEdges(K_FLD::IntArray& connect, std::set<EdgeType>& edges);

    /** accessors*/
    inline const tree_type& getTree() const {return _tree;} 

    inline E_Float getTolerance() const { return _tolerance;}

    static E_Int aggregate_convex(const K_FLD::FloatArray&crd, const K_FLD::IntArray& connectT3, const E_Float* normal, ngon_unit& agg_pgs, E_Float convexity_tol = CONVEXITY_TOL);

    static E_Int get_polygonal_boundary
      (const K_FLD::IntArray& cT3, const std::vector<E_Int>& ids, std::deque<E_Int>& PGi,
      std::set<K_MESH::Edge>& w_oe_set, std::map<E_Int, E_Int>& w_n_map);
    static E_Int starify_from_node(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, E_Int istar, K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, const E_Float* normal, E_Float convexity_tol);
    
    static E_Int cut_mesh_convex_line
      (const K_FLD::FloatArray& crd, const K_FLD::IntArray& connectT3, const E_Float* normal, const K_FLD::IntArray& connectB,
      E_Int K0, E_Int n0, K_FLD::IntArray& neighbors,
      K_FLD::IntArray& connectT31, K_FLD::IntArray& connectT32);

    static E_Int get_edges_lying_on_plane(const K_FLD::FloatArray& crd, E_Int index_start, const std::set<K_MESH::NO_Edge>& edges, E_Int Np, const E_Float* normal, E_Float tol_rel, std::set<K_MESH::NO_Edge>& lyingEs);
    static void burn_free_branches(std::set<K_MESH::NO_Edge>& edges, std::map<E_Int, E_Int>& node_to_count);
    static E_Float get_max_deviation (const K_FLD::FloatArray& crd, const K_FLD::IntArray& cT3, const K_FLD::IntArray& neighT3);
    
    template <typename Triangulator>
  static void refine_T3s(const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, 
                         std::map<K_MESH::NO_Edge, Vector_t<E_Int> >& edge_to_refined_edge, std::vector<E_Int>& oids);
    
    static inline void append_all_topo_paths (const K_FLD::FloatArray& crd, E_Int Nstart, E_Int Nend, const id_to_ids_t& nodes_graph, idlists_t& paths);
    
    template <typename IntCONT>
    static inline void get_farthest_point_to_edge (const K_FLD::FloatArray& crd, E_Int Ni, E_Int Nj, const IntCONT& list, E_Int& Nf, E_Float& d2);
    
private:

  E_Int __getContainingElement(const E_Float* point, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                               const K_FLD::IntArray& neighbors, const bool_vector_type& mask, size_type Kseed = E_IDX_NONE) const;

  E_Int __searchDirection(size_type K, const E_Float* point, const K_FLD::FloatArray& pos,
                          const K_FLD::IntArray& connect, E_Bool random) const;

  ///
  template <typename StoringTriangleType>
  static bool detectDuplicated(const K_FLD::IntArray& connect, std::vector<E_Int>& dupIds);

  ///
  template <typename StoringTriangleType>
  static E_Int removeDuplicated(K_FLD::IntArray& connect, std::vector<E_Int>& dupIds);



private:

  const tree_type&                  _tree;
  E_Float                           _tolerance;
  mutable std::deque<size_type>     _pool;
public://fixme
  mutable int_set_type              _inval;
};

} // End namespace K_CONNECT

#include "MeshTool.cxx"

#endif /* K_CONNECT_MESHTOOL_H_ */
