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

#ifndef __TRI_CONFORMIZER_H__
#define __TRI_CONFORMIZER_H__

#include "Conformizer.h"

#include <vector>
#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/T3Mesher.h"

namespace NUGA
{

template<short DIM>
class TRI_Conformizer : public Conformizer<DIM, K_MESH::Triangle>
{

public :
  typedef Conformizer<DIM, K_MESH::Triangle>   parent_type;
  typedef K_MESH::Triangle                     element_type;
  typedef std::vector< std::vector<E_Int> >    edge_container_type;
  typedef NUGA::EltAlgo<K_MESH::Triangle> algo_type;

public:
  /// Constructor.
  TRI_Conformizer(bool wnh = false/*with node history*/);
  /// Destructor.
  virtual ~TRI_Conformizer(void){}
  
  // Overridden Methods ///////////////////////////////////////////////////////

#ifndef DEBUG_TRI_CONFORMIZER
protected:
#else
public:
#endif
  ///
  void __set_tolerances(E_Float Lmin, E_Float Lmax, E_Float  user_tolerance);
  ///
  void __prepare_data(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect);
  ///
  E_Int __intersect(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                     T3& t1, T3& t2, E_Float tol);
  ///
  void __update_data(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& newIDs);
  
  /// Splits the triangles by triangulation.
  E_Int __split_Elements(const K_FLD::FloatArray& pos, K_FLD::IntArray & connect,
                        NUGA::bool_vector_type& xc,
                        NUGA::int_vector_type& ancestors);
  
  ///Hook to do a merge toward the intersection line (TRI)
  E_Int __simplify_and_clean(const K_FLD::FloatArray& pos, E_Float tolerance, K_FLD::IntArray& connect,
                             NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc);
  ///Hook to do a merge toward the intersection line (TRI)
  E_Int __simplify_and_clean2(const K_FLD::FloatArray& pos, E_Float tolerance, K_FLD::IntArray& connect,
    NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc);
  /// when duplicates are baffles
  void __process_duplicates(K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, NUGA::bool_vector_type& xc);
  
  ///Hook to manage overlapping zones
  void __run_correction_beta(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                            NUGA::int_vector_type& ancestors,
                            NUGA::bool_vector_type& xc,
                            E_Float tolerance);

  void __run_correction_gamma(const std::set<K_MESH::NO_Edge>&xpairs, NUGA::int_vector_type&colors,
                              NUGA::int_vector_type&xr, K_FLD::IntArray& connect,
                              NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc
                              , const K_FLD::FloatArray& pos, Vector_t<E_Int>* priority);
  
///
/*void __clean_topology(K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, NUGA::bool_vector_type& xc, 
                       const std::vector<E_Int>& source_nodes_vec, std::vector<E_Int>& common_nodes_vec,
                       const K_FLD::IntArray& connectBNM, const std::vector<E_Int>& nids
#ifdef DEBUG_TRI_CONFORMIZER
        ,const K_FLD::FloatArray& coord
#endif
);
///
E_Int __swap_edges(std::vector<std::pair<E_Int, E_Int> >& sapwE, K_FLD::IntArray& connectM, 
                   std::vector<E_Int> ancestors, K_FLD::IntArray& neighbors,
                   const std::vector<E_Int>& nids, std::set<K_MESH::NO_Edge>& commonEdges);
*/

#ifdef DEBUG_TRI_CONFORMIZER
  ///
  void drawElements(const char* fname, const char* filefmt, const K_FLD::FloatArray& coord,
                    const K_FLD::IntArray& connect, const std::vector<T3> & elts, bool localid = false, std::vector<E_Int>* colors=0);
  ///
  void draw(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Float tol, E_Int what, K_FLD::IntArray& out);

#endif
  
  /////////////////////////////////////////////////////////////////////////////

#ifndef DEBUG_TRI_CONFORMIZER
private:
#else
public:
#endif
  
  ///
  E_Int __iterative_run(DELAUNAY::T3Mesher<E_Float>& mesher, K_FLD::FloatArray& crd, K_FLD::IntArray& cB, 
                        NUGA::int_vector_type& hnodes, DELAUNAY::MeshData& data, std::vector<E_Int>& nids,
                        bool ignore_unforceable_edge, bool silent_last_iter);
  ///
  E_Int __get_mesh_data(const K_FLD::FloatArray & pos, const K_FLD::IntArray & connect, const T3& t, edge_container_type& Edges,
                        K_FLD::FloatArray& p, K_FLD::IntArray& c, std::vector<E_Int>& oids);
  ///
  void __improve_triangulation_quality(const T3& tri, const std::vector<E_Int>& revIDs,
                                       std::vector<std::pair<E_Int, E_Int> >& wpair_set,
                                       std::set<K_MESH::NO_Edge>& wedge_set,
                                       DELAUNAY::MeshData& data);
  ///
  E_Bool __intersect(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, T3& t,
                     edge_container_type& edges, E_Int idxE, E_Float tol, std::vector<E_Int>& nodes, E_Int* pidx, E_Bool& coplanar);
  
  ///
  bool __fast_discard(const K_FLD::FloatArray& pos, const E_Int* T1, const E_Int* T2, E_Float tol);
  
  ///
  void __tidy_edges(const K_FLD::FloatArray& pos, edge_container_type& edges);
  ///
  void __tidy_edge (const K_FLD::FloatArray& coord, std::vector<E_Int>& nodes);
  ///
  void __reorder_nodes_on_edge(const K_FLD::FloatArray& pos, std::vector<E_Int>& nodes);
  
  ///
  E_Int __get_connectB2(const K_FLD::FloatArray & pos, const K_FLD::IntArray & connect,
                       const T3& t, edge_container_type& edges, std::set<K_MESH::Edge>& hBO,
                       std::vector<E_Int>& hNodes);
  ///
  inline bool __get_B0_edges(const K_FLD::IntArray & connect,
                      const T3& t, edge_container_type& Edges,
                      std::set<K_MESH::Edge>& hBO, std::set<E_Int>& Nodes0);
  ///
  inline void __get_Inner_edges(const K_FLD::FloatArray & pos, const K_FLD::IntArray & connect,
                         const T3& t, edge_container_type& Edges,
                         std::set<K_MESH::Edge>& hBO, std::set<E_Int>& Nodes0);
  ///
  inline void __get_Imprint_edges(const T3& t, std::set<K_MESH::Edge>& hBO, std::set<E_Int>& Nodes0);
  ///
  E_Bool is_inside(const K_FLD::IntArray& connect, const T3& t, const K_FLD::FloatArray& pos, E_Int Ni, E_Float tol);
  ///
  void __compact_to_mesh(const K_FLD::FloatArray& posIn, const K_FLD::IntArray& connectIn, 
                         K_FLD::FloatArray& posOut, K_FLD::IntArray& connectOut, 
                         std::vector<E_Int>& revIDs, std::vector<E_Int> * hard_nodes = 0);
  ///
  void __transform(const E_Float* P0, const E_Float* P1, const E_Float* P2, K_FLD::FloatArray& pos);
  ///
  void __transform(K_FLD::FloatArray& pos, const K_FLD::FloatArray& t);


  private: // correction specific
    ///
    void __get_zones(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                     const algo_type::BoundToEltType& E_to_T,
                     const std::set<K_MESH::NO_Edge>& common_edges, std::vector< std::vector<E_Int> >& zones);
    ///
    E_Int __get_common_edges(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& dupIds,
                             const NUGA::EltAlgo<K_MESH::Triangle>::BoundToEltType& bound_to_elts,
                             std::set<K_MESH::NO_Edge>& common_edges,
                             std::vector<bool>& patho_elts, std::vector<bool>& free_elts);
    ///
    void __compute_splitting_points(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                                    const NUGA::bool_vector_type& xc,
                                    const std::set<K_MESH::NO_Edge>& commone_no_edges,
                                    E_Float tol, std::map<K_MESH::NO_Edge, std::vector<E_Int> >& edge_to_point);
    ///
    void __buildSameSupportSurfacePairs(const K_FLD::FloatArray& pos, 
                                        const std::vector<K_FLD::IntArray >& connectZs, std::vector<std::vector<E_Int> >& grouped_zones);
    ///
    E_Bool __areOverlapping(const K_FLD::FloatArray& pos, E_Float tolerance, const K_FLD::FloatArray& normals,
                            const K_FLD::IntArray& connectZ1, const K_SEARCH::BbTree3D& treeZ1,
                            const K_FLD::IntArray& connectZ2, const K_SEARCH::BbTree3D& treeZ2);
    ///
    E_Bool __IsNodeFarFromSurface(const K_FLD::FloatArray& pos, const K_FLD::FloatArray& normals,
                                  E_Int N, const K_FLD::IntArray& connectZ, const K_SEARCH::BbTree3D& treeZ, E_Float tolerance);
    ///
    void __computeAveragedNormal(const K_FLD::FloatArray& normals, const std::vector<E_Int>& nodes, E_Float* ni);

#ifdef DEBUG_TRI_CONFORMIZER
    ///
    void drawTandE(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS, E_Int Ni, E_Int Nj);
    ///
    void T3nodes (const K_FLD::IntArray& connect, const T3& t, std::set<E_Int>& nodes);
    ///
    void drawT3(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Int ith_elts, bool compact=false);
    ///
    bool detect_duplis_and_baffles(const K_FLD::IntArray& connect);
#endif

private:

  std::vector<std::vector<E_Int> > _edges;
  std::vector<E_Int> _wnodes_vec, _wnodes_vecclean;
  std::set<E_Int>    _wnodes_set;

  K_FLD::FloatArray _P,_iP;
  E_Float _U1[3], _U2[3], _U3[3];

#ifdef FLAG_STEP
  E_Float tcomp, ttra, trun, tcon, tcon0, tcon1, tcon2, tnot;
#endif

};
}

#include "TRI_Conformizer.cxx"
#include "TRI_correction.cxx"

#endif
