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
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef __GENERATOR_DELAUNAY_KERNEL_H__
#define __GENERATOR_DELAUNAY_KERNEL_H__

#include "Predicates.h"
#include "MeshData.h"
#include "Connect/MeshTool.h"
#include "macros.h"

#include <map>

namespace DELAUNAY
{

typedef   float                                    CONSTRAINED ;
typedef   int                                      UNCONSTRAINED ;

template <typename T>
class Kernel
{

 public: /** Typedefs */

    typedef   K_CONT_DEF::size_type                    size_type;
    typedef   K_CONT_DEF::int_vector_type              int_vector_type;
    typedef   K_CONT_DEF::int_set_type                 int_set_type;
    typedef   K_CONT_DEF::int_pair_type                int_pair_type;
    typedef   K_CONT_DEF::int_pair_vector_type         int_pair_vector_type;
    typedef   K_CONT_DEF::int_pair_set_type            int_pair_set_type;
    typedef   K_CONT_DEF::bool_vector_type             bool_vector_type;
    typedef   K_CONT_DEF::non_oriented_edge_set_type   non_oriented_edge_set_type;

    typedef   Kernel                                   self_type;
    typedef   K_MESH::Triangle                         element_type;

    typedef   no_predicate<K_MESH::NO_Edge>            unconstrained_predicate;
    typedef   HardEdgeCriterion                        constrained_predicate;
        
public:

  /// Constructor.
  explicit Kernel(MeshData& data, const K_CONNECT::MeshTool& tool);
 
         
  ///
  ~Kernel(void);

  ///
  template <typename ConstraintType>
  E_Int insertNode(size_type N, const T& m, const ConstraintType& dummy);

  /// sets the hard edges
  void setConstraint(const K_CONT_DEF::non_oriented_edge_set_type& hard_edges);

private:

  /// Gets the cavity associated to the input node N regarding the metric value m.
  /** Here m is evaluated at N but a more sophisticated approach could be to take a mean of the triangle summits.*/
  template <typename isConstrained>
  E_Int __getCavity(size_type N, const T& m, const K_FLD::FloatArray& pos,
                    const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, int_set_type& cavity,
                    int_pair_vector_type& cboundary);

  /// Remeshes the cavity
  E_Int __remeshCavity(size_type N, K_FLD::IntArray & connect, K_FLD::IntArray& neighbors,
                      int_vector_type& ancestors, const int_set_type& cavity,
                      const int_pair_vector_type& cboundary);

  /// Invalidates old elements of the remeshed cavity.
  void __invalidCavityElements(const int_set_type& cavity, const K_FLD::IntArray& connect,
                               bool_vector_type& mask);

  ///
  template <typename isConstrained>
  E_Int __getInitialCavity(size_type N, 
                           const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, 
                           const K_FLD::IntArray& neighbors,
                           const int_vector_type& ancestors, int_set_type& base, int_set_type& cavity,
                           int_pair_set_type& cboundary);

  ///
  template <typename isConstrained>
  inline
  void __appendCavity(const K_FLD::IntArray& neighbors, const int_set_type& base,
                      int_set_type& cavity, int_pair_set_type& cboundary);

  ///
  E_Int __fixCavity(size_type N, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors,
                    const int_vector_type& ancestors, const int_set_type& base, int_set_type& cavity,
                    int_pair_set_type& cboundary);

  ///
  E_Int __ensureEmptyCavity(const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, 
                            const int_set_type& base, int_set_type& cavity, int_pair_set_type& cboundary);

  ///
  E_Int __ensureStarShape(size_type N, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectM, 
                          const K_FLD::IntArray& neighbors, const int_set_type& base,
                          int_set_type& cavity, int_pair_set_type& cboundary);

  ///
  void __getSortedBoundary(const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors, size_type Ki, size_type b0,
                           const int_set_type& cavity, int_pair_set_type& cboundary,
                           int_set_type& visitedK, int_pair_vector_type& sorted_boundary); 

  E_Int __getSortedBoundary2(const K_FLD::IntArray& connect, int_pair_set_type& in, int_pair_vector_type& out);

  ///
  void __getBoundary(int_set_type& cavity, const K_FLD::IntArray& neighbors, int_pair_set_type& cboundary);

private:

  MeshData& _data;
  const K_CONNECT::MeshTool& _tool;
  DelaunayCriterion<T>              _Ball_pred;
  unconstrained_predicate           _unconstrained_pred;
  constrained_predicate             *_constrained_pred;

  // Temporary containers (put here for speed).
  int_set_type          _cavity;
  int_pair_vector_type  _cboundary;
  int_set_type           _base;
  int_pair_set_type     _sbound, _real_cboundary;
  int_set_type          _visited;
  int_set_type          inodes;
  int_vector_type       Ancs;
  int_vector_type       elements;
  std::map<size_type, int_pair_set_type::iterator> _node_to_rightS;
  
public:
 size_type _Nmatch;
#ifdef E_TIME
public:
  double inval_time;
  double remesh_time;
  double cavity_time;
  double init_cavity_time;
  double sorting_bound_time;
  double fix_cavity_time;
  double _base_time;
  double _append_time;
#endif

};

}

#include "Kernel.cxx"

#endif

