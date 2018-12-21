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

#ifndef __CONFORMIZER_H__
#define	__CONFORMIZER_H__

#include <vector>
#include "Def/DefTypes.h"
#include "Fld/DynArray.h"
#include "Def/DefContainers.h"
#include "Search/BbTree.h"

#ifdef DEBUG_CONFORMIZER
#define TRUPLE 533111,533401,533110
static E_Int zTi = 436;

#endif

namespace NUGA
{
  
template<typename Element_t>
struct isDegenerated : public std::unary_function <const K_FLD::IntArray::iterator&, bool>
{
  explicit isDegenerated(){}

  inline bool operator() (const  K_FLD::IntArray::iterator& pS) const
  {return ((*pS == *(pS+1)) || (*pS == *(pS+2)) || (*(pS+1) == *(pS+2)));}
};

template<>
struct isDegenerated<K_MESH::Edge> : public std::unary_function <const K_FLD::IntArray::iterator&, bool>
{
  explicit isDegenerated(){}

  inline bool operator() (const  K_FLD::IntArray::iterator& pS) const
  {return (*pS == *(pS+1));}
};

template <typename Element_t>
struct Struc
{
typedef void Type;
};

struct T3
{
  T3(E_Int si):id(si){}
  E_Int id; // Triangle index.
  std::vector<E_Int> edges;
  K_FLD::IntArray    Xedges;
};

template <>
struct Struc<K_MESH::Triangle>
{
  typedef T3 Type;
};

struct E2
{
  E2(E_Int e):id(e){}
  E_Int id; // Edge index.
  std::vector<E_Int> nodes;
};

template <>
struct Struc<K_MESH::Edge>
{
  typedef E2 Type;
};

class ConformizerRoot{
public:
  virtual E_Int run(K_FLD::FloatArray& coord, K_FLD::IntArray& cBAR, std::vector<E_Int>& ancestor, std::vector<E_Int>* ovlp_priority=0, 
            E_Float tolerance = -1., E_Int X0=0, E_Int itermax=10) = 0;
  virtual ~ConformizerRoot(){}
protected:
  ConformizerRoot(){}
  
#ifdef DEBUG_CONFORMIZER
  static E_Int fastdiscard_counter, xtest_counter;
  static E_Int split_fastdiscard_counter, split_counter, degen_counter;
#endif
  
};


template <E_Int DIM, typename Element_t>
class Conformizer : public ConformizerRoot {

public:
  /// Main method.
  E_Int run(K_FLD::FloatArray& coord, K_FLD::IntArray& cBAR, std::vector<E_Int>& ancestor, std::vector<E_Int>* ovlp_priority=0,
            E_Float tolerance = -1., E_Int X0 = 0, E_Int itermax=10);
  ///
  const std::vector<E_Int>& get_node_history(){return _node_history;}
  ///
  virtual ~Conformizer(){}
  
public://fixme
  bool _split_swap_afterwards; // TRI specific

protected:
  ///
  Conformizer(bool wnh = false): _split_swap_afterwards(false), _absolute_tol(true)/*fixme : no choice yet*/, _with_node_history(wnh), _X0(0) {}
  
  // Methods to override : interface to implement
protected:
  ///
  virtual void __set_tolerances(E_Float Lmin, E_Float Lmax, E_Float  user_tolerance) = 0;
  ///
  virtual void __prepare_data(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect);

#ifndef DEBUG_CONFORMIZER
private:
#else
public:
#endif
  ///
  virtual E_Int __intersect(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                             typename Struc<Element_t>::Type& t1, typename Struc<Element_t>::Type& t2, E_Float tol) = 0;
  ///
  virtual void __update_data(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& newIDs) = 0;
  ///
  virtual E_Int __split_Elements(const K_FLD::FloatArray& pos, K_FLD::IntArray & connect,
                                K_CONT_DEF::bool_vector_type& xc,
                                K_CONT_DEF::int_vector_type& ancestors) = 0;
  
  /// Hook inside __run to merge toward the intersection line (TRI)
  virtual E_Int __simplify_and_clean(const K_FLD::FloatArray& pos, E_Float tolerance, K_FLD::IntArray& connect,
                                     K_CONT_DEF::int_vector_type& ancestors, K_CONT_DEF::bool_vector_type& xc){return 0;}
  /// Hook inside __run to merge toward the intersection line (TRI)
  virtual E_Int __simplify_and_clean2(const K_FLD::FloatArray& pos, E_Float tolerance, K_FLD::IntArray& connect,
    K_CONT_DEF::int_vector_type& ancestors, K_CONT_DEF::bool_vector_type& xc){return 0;}
  /// Hook inside run after the loop to manage overlapping zones
  virtual void __run_correction_beta(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                                     K_CONT_DEF::int_vector_type& ancestors, 
                                     K_CONT_DEF::bool_vector_type& xc,
                                     E_Float tolerance){}
  /// Hook inside run after the loop to manage overlapping zones
  virtual void __run_correction_gamma(const std::set<K_MESH::NO_Edge>&xpairs, K_CONT_DEF::int_vector_type&colors,
                                      K_CONT_DEF::int_vector_type&xr, K_FLD::IntArray& connect,
                                      K_CONT_DEF::int_vector_type& ancestors, K_CONT_DEF::bool_vector_type& xc
                                      , const K_FLD::FloatArray& pos, K_CONT_DEF::int_vector_type* priority){}
   
  /////////////////////////////////////////////////////////////////////////////
  
  // Common Factorized Methods
   
#ifndef DEBUG_CONFORMIZER
protected:
#else
public:
#endif
  
 ///
 E_Int __merge_clean (E_Float tol, const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, 
               std::vector<E_Int>& newIDs, K_CONT_DEF::bool_vector_type* xc, E_Int fromIdx);
 ///
 E_Int __merge_clean
 (E_Float tol, const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, 
 std::vector<E_Int>& new_IDs, K_CONT_DEF::bool_vector_type* xc, 
 const K_CONT_DEF::int_vector_type& source, const K_CONT_DEF::int_vector_type& target);
 
 ///
 void __clean
 (const std::vector<E_Int>& new_IDs, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, K_CONT_DEF::bool_vector_type* xc);

  ///
  E_Int __removeDegenerated(K_FLD::IntArray& connect, K_CONT_DEF::int_vector_type& newIDs);
  
#ifdef DEBUG_CONFORMIZER
  virtual void drawElements(const char* fname, const char* filefmt, const K_FLD::FloatArray& coord,
                            const K_FLD::IntArray& connect, 
                            const std::vector<typename Struc<Element_t>::Type> & elts, bool localid = false, std::vector<E_Int>* colors=0) = 0;
  
  virtual void drawT3(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Int ith_elts, bool compact=false){}
  E_Int _nbX;
  virtual bool detect_duplis_and_baffles(const K_FLD::IntArray& connect){return false;}
#endif
   
#ifndef DEBUG_CONFORMIZER
private:
#else
public:
#endif
  ///
  E_Int __merge_duplicated_nodes(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                                        E_Float tolerance, E_Int fromIdx = 0);
  ///
  void __initialize(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, E_Float tolerance,
                    std::vector<E_Int>& ancestors, K_CONT_DEF::bool_vector_type& xc, E_Int X0);
  ///
  void __finalize();
  ///
  E_Int __run(K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
              std::vector<E_Int>& ancestors, K_CONT_DEF::bool_vector_type& xc, E_Float tolerance = -1.); 
  ///
  E_Int __compute_intersections(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                                        K_CONT_DEF::bool_vector_type& xc, E_Float tolerance);
  ///
  void __compute_min_edge_length (const K_FLD::FloatArray& pos,
                                  const K_FLD::IntArray& connect, const K_CONT_DEF::bool_vector_type& xc, E_Float &Lmin, E_Float& Lmax);
  
  ///
  void __update_tolerance(const K_FLD::FloatArray& coord, K_CONT_DEF::bool_vector_type& xc, const std::vector<E_Int>& ancestors);
  

#ifndef DEBUG_CONFORMIZER
protected:
#else
public:
#endif
  
  std::vector<K_SEARCH::BoundingBox<DIM>*> _boxes, bxtmp;
  K_SEARCH::BoundingBox<DIM>* _pool;
  K_SEARCH::BbTree<DIM> *_tree;
  std::vector<typename Struc<Element_t>::Type> _elements;
  E_Float _tolerance, _tol_x, _tol_clean;
  E_Int _iter;
  bool _absolute_tol;
  std::vector<std::pair<E_Float, E_Int> > _sorterFI;
  E_Int _initialrows;
  bool _with_node_history;
  std::vector<E_Int> _node_history;
  E_Int _N0;
  K_FLD::IntArray _connect0;
  E_Int _X0; // start testing from that id (when 2 valid input mesh, avoid self-X tests)
  E_Int _itermax;
  bool _needs_another_iter;
  bool _one_pass_mode;
  
public:
  std::set<K_MESH::NO_Edge> _xpairs;
  std::vector<E_Int> _colors, _xr, *_priority;


};

}
#include "Conformizer.cxx"

#endif	/* CONFORMIZER_H */

