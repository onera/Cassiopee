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

//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __CONFORMIZER_CXX__
#define __CONFORMIZER_CXX__

#include "Nuga/include/DynArray.h"


#include "Conformizer.h"
#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/IdTool.h"

#ifdef DEBUG_CONFORMIZER
#include <sstream>
#include "IO/io.h"
#include <iostream>
static int xtesttotal=0;
static int xtestreal=0;
#endif

#ifdef DEBUG_TRI_CONFORMIZER
#include "TRI_debug.h"
#endif

#ifdef FLAG_STEP
#include "Nuga/include/chrono.h"
#include <iostream>
#endif

#define MAX_PERCENT 0.05
#define ZERO_MACHINE EPSILON

namespace NUGA
{
  
//
template <short DIM, typename Element_t>
void
Conformizer<DIM, Element_t>::__prepare_data
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect)
{
  E_Int nb_elts = connect.cols();
  
  _elements.clear();
  _elements.reserve(nb_elts);
  K_CONNECT::IdTool::init_inc(_elements, nb_elts);
  
  if (_brute_force) return;

  _boxes.clear();
  _boxes.resize(nb_elts);

  _pool = new K_SEARCH::BoundingBox<DIM>[nb_elts];

  for (E_Int i = 0; i < nb_elts; ++i)
  {
    _boxes[i] = &_pool[i];
    _boxes[i]->compute(pos, connect.col(i), Element_t::NB_NODES);
  }
  if (_iter == 1 && _X0 > 0)
  {
    bxtmp = _boxes;
    bxtmp.resize(_X0);
    _tree = new K_SEARCH::BbTree<DIM>(bxtmp);
  }
  else
    _tree = new K_SEARCH::BbTree<DIM>(_boxes);
}
  
///
template <short DIM, typename Element_t>
E_Int Conformizer<DIM, Element_t>::__removeDegenerated
(K_FLD::IntArray& connect, NUGA::int_vector_type& newIDs)
{
  E_Int                       Si, COLS(connect.cols()), ROWS(connect.rows());
  K_FLD::IntArray::iterator   pS;
  K_FLD::IntArray             connectOut;
  isDegenerated<Element_t>    isDegen;
  
  connectOut.reserve(ROWS, COLS);

  newIDs.clear();
  newIDs.resize(COLS, IDX_NONE);

  for (Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);
    
    //if ((*pS != *(pS+1)) && (*pS != *(pS+2)) && (*(pS+1) != *(pS+2)))
    if (!isDegen(pS))
    {
      connectOut.pushBack(pS, pS+ ROWS);
      newIDs[Si] = connectOut.cols() - 1;
    }
  }

  connect = connectOut;
 
  return (COLS - connect.cols());
}
  
///
template <short DIM, typename Element_t>
E_Int Conformizer<DIM, Element_t>::run
(K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, Vector_t<E_Int>* priority, E_Float tolerance, E_Int X0, E_Int itermax)
{
  E_Int                         iter(1), xnb;
  bool                          carry_on(false);
  NUGA::bool_vector_type  xc;

  if (connect.cols() == 0) return 0;
  
  _itermax = itermax;
  _one_pass_mode = (_itermax == 1);
  _priority = priority;

  //assert(pos.rows() == DIM);
  _initialrows = -1;
  if (DIM < pos.rows()) //fixme : right place to make the problem 2D ?
  {
    _initialrows = pos.rows();
    pos.resize(DIM, pos.cols());
  }
  
#ifdef FLAG_STEP
     std::cout << "Conformizer::__initialize: X0 "        << X0 << std::endl;
     std::cout << "Conformizer::__initialize: tolerance " << tolerance << std::endl;
     std::cout << "Conformizer::__initialize: nb elts "   << connect.cols() << std::endl;
     std::cout << "Conformizer::__initialize: xc size "   << xc.size() << std::endl;
#endif

  this->__initialize(pos, connect, tolerance, ancestors, xc, X0);

  do
  {
#ifdef FLAG_STEP
    std::cout << "Conformizer::run: ITER " << iter << std::endl;
    std::cout << "Conformizer::run: Initial Nb Points : " << pos.cols() << std::endl;
#endif

    ++_iter;
    _needs_another_iter=false;

#ifdef FLAG_STEP
    chrono chr,c;
    chr.start();c.start();
#endif
        
    // Initialize data structures and create the BST.
    this->__prepare_data(pos, connect);
    
#ifdef FLAG_STEP
    std::cout << "Conformizer::run : __prepare_data : " << c.elapsed() << std::endl;
    c.start();
#endif

    //
    xnb = this->__run(pos, connect, ancestors, xc, _tolerance);
    
    if (_needs_another_iter)
    { ++_itermax; _itermax = std::min((E_Int)10, _itermax);} // max thershold set to 10.
    
#ifdef FLAG_STEP
    std::cout << "Conformizer::run : __run : " << c.elapsed() << std::endl;
    c.start();
#endif

    // Delete the BST.
    this->__finalize();
    
    if (xnb < 0) break;
    
#ifdef FLAG_STEP
    std::cout << "Conformizer::run : __finalize : " << c.elapsed() << std::endl;
#endif
    
    // If not specified by the user, update (if bigger) the tolerance based only on the intersecting ORIGINAL triangles set.
   if ( (_iter == 1) && (tolerance == 0.) && !_one_pass_mode ) // First pass.
   {
     //E_Float oldtol=_tolerance;
     __update_tolerance(pos, xc, ancestors);
     //if (_tolerance > oldtol) fixme : doesn't work. But something might be done here to clear small elements...
       //this->__simplify_and_clean(pos, _tolerance, connect, ancestors, xc);
   }

#ifdef FLAG_STEP
    std::cout << "Conformizer::run : Final Nb Points : " << pos.cols() << std::endl;
#endif

#ifdef FLAG_STEP
     std::cout << "Conformizer::run : iter time : " << chr.elapsed() << std::endl;
#endif
    
    _connect0.release(); //no required anymore

    carry_on = (xnb > 0) && (iter++ < _itermax);
  }
  while (carry_on);
  
  if ((xnb && !_one_pass_mode) || (_one_pass_mode && xnb<0))  // could not converge properly or the T3Mesher failed for an element (-1)
    return 1;
  
  // for TRI only
  if (!_one_pass_mode) //do no use when attempt a one pass procees (e.g. for the boolean NGON)
  {
#ifdef FLAG_STEP
    std::cout << "__run_correction_beta " << pos.cols() << std::endl;
#endif
    this->__run_correction_beta(pos, connect, ancestors, xc, _tolerance);
  }
  else // process wrong overlap splits
  {
    this->__run_correction_gamma(_xpairs, _colors, _xr, connect, ancestors, xc, pos, priority);
  }

  if (_initialrows > 0) //fixme : can we exit with a real 2D coordinate matrix ?
  {
    E_Float zero = 0.;
    pos.resize(_initialrows, pos.cols(), &zero);
  }
  
//#ifdef E_DEBUG
  if (_with_history)
  {
    //check histo validity with distance
  }
//#endif

  if (connect.cols() == 0)
    return 1;
  return 0;
}

///
template <short DIM, typename Element_t>
void
Conformizer<DIM, Element_t>::__initialize
(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, E_Float tolerance,
 std::vector<E_Int>& ancestors, NUGA::bool_vector_type& xc, E_Int X0)
{
  _X0 = X0;
  _X1 = 0;
  if (X0 < 0)
  {
    _X0 = -X0;
    _X1 = pos.cols(); 
  }
  
  // initialize T3s ancestors
  ancestors.resize(connect.cols(), 0);
  for (size_t i = 0; i < ancestors.size(); ++i)ancestors[i] = i;
  
  // initialize node history
  if (_with_history)
  {
    _node_history.resize(pos.cols());
    for (size_t i = 0; i < _node_history.size(); ++i)_node_history[i] = i;
  }
  
  _connect0=connect;

  // Initial cleaning
  std::vector<E_Int> nids;
  //E_Int nb_merges = 
  __merge_clean(EPSILON, pos, connect, ancestors, nids, 0, 0);
  
  xc.clear();
  xc.resize(connect.cols(), true); // done to ensure that tolerance is computed taking into account everything

  // Initialize tolerances
  E_Float Lmin, Lmax;
  __compute_min_edge_length(pos, connect, xc, Lmin, Lmax);
    
  this->__set_tolerances(Lmin, Lmax, tolerance);
  
  // now set xc properly taking into account _X0
  if (_X0 != 0) // avoid self-X tests
  {
    xc.clear();
    xc.resize(_X0, false);
    xc.resize(connect.cols(), true);
  }

//#ifdef DEBUG_CONFORMIZER
#ifdef FLAG_STEP
  {
    std::cout << "input tol       : "  << tolerance << std::endl;
    std::cout << "min edge length : "  << Lmin << std::endl;
    std::cout << "_tolerance      : " << _tolerance << std::endl;
  }
#endif
  //medith::write("cleaned.tp", pos, connect, "BAR");
//#endif

}

///
template <short DIM, typename Element_t>
E_Int
Conformizer<DIM, Element_t>::__merge_clean
(E_Float tol, const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, 
 std::vector<E_Int>& new_IDs, NUGA::bool_vector_type* xc, E_Int fromIdx)
{
  K_FLD::ArrayAccessor<K_FLD::FloatArray> posAcc(pos);
  NUGA::int_vector_type nodes;
     
  // Do The merge.
  if (fromIdx==0)
  {
    nodes.resize(posAcc.size());
    for (size_t i = 0; i < nodes.size(); ++i)nodes[i]=i;
    
    return __merge_clean(tol, pos, connect, ancestors, new_IDs, xc, nodes, nodes);
  }
  else if (fromIdx > 0)
  {
    std::vector<E_Int> source, target;
    connect.uniqueVals(target);
    for (E_Int i = fromIdx; i < pos.cols(); ++i) source.push_back(i);
    target.insert(target.end(), source.begin(), source.end());
    
    return __merge_clean(tol, pos, connect, ancestors, new_IDs, xc, source, target);
  }
  else // (fromIdx < 0) => _X1 is set : rule for target : [0, ... X0, X1+1, ...pos.cols(), X0+1, ... X1]. The less prior are those in the tail of input pos 
  {
    std::vector<E_Int> source, target;
    for (E_Int i = 0; i < _X0; ++i) target.push_back(i);
    for (E_Int i = _X1; i < pos.cols(); ++i) target.push_back(i);
    for (E_Int i = _X0; i < _X1; ++i) target.push_back(i);
    
    source = target;
    std::reverse(ALL(source));

    return __merge_clean(tol, pos, connect, ancestors, new_IDs, xc, source, target);
  }
}

///
template <short DIM, typename Element_t>
E_Int
Conformizer<DIM, Element_t>::__merge_clean
(E_Float tol, const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, 
 std::vector<E_Int>& new_IDs, NUGA::bool_vector_type* xc, const NUGA::int_vector_type& source, const NUGA::int_vector_type& target)
{
    K_FLD::ArrayAccessor<K_FLD::FloatArray> posAcc(pos);   
    E_Int nb_merges(0);
    
    new_IDs.clear();
    
    // Do The merge.
    nb_merges = merge(posAcc, tol, source, target, new_IDs);

    if (nb_merges > 0)
      __clean(new_IDs, connect, ancestors, xc);
    else
      new_IDs.clear();
      
    return nb_merges;
}

///
template <short DIM, typename Element_t>
void
Conformizer<DIM, Element_t>::__clean
(const std::vector<E_Int>& new_IDs, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, NUGA::bool_vector_type* xc)
{  
  K_FLD::IntArray::changeIndices(connect, new_IDs);
  K_FLD::IntArray::changeIndices(_connect0, new_IDs);
      
  if (_with_history)
    K_CONNECT::IdTool::propagate(new_IDs, _node_history);
      
  std::vector<E_Int> nids;
  if (__removeDegenerated(connect, nids))// returns the new ids of connect.
  {
    K_CONNECT::valid pred(nids);
    K_CONNECT::IdTool::compress(ancestors, pred);
    if (xc)
      K_CONNECT::IdTool::compress(*xc, pred);
    for (auto& e : _elements)
      e.id = nids[e.id];
    K_CONNECT::IdTool::compress(_elements, pred);
  }
}

///
template<short DIM, typename Element_t>
E_Int Conformizer<DIM, Element_t>::__run
(K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& ancestors,
 NUGA::bool_vector_type& xc, E_Float tolerance)
{

#ifdef FLAG_STEP
  chrono tt;
  tt.start();
#endif
  
#ifdef DEBUG_CONFORMIZER
  if (detect_duplis_and_baffles(connect))
    std::cout << "oddities : after __prepare_data : FOUND !!!!" << std::endl;
#endif

#ifdef FLAG_STEP
   std::cout << "Conformizer : init : " << tt.elapsed() << std::endl;
  tt.start();
#endif

#ifdef DEBUG_TRI_CONFORMIZER
  {
  std::ostringstream o;
  o << "beforeX_" << _iter << ".mesh";
  medith::write(o.str().c_str(), pos, connect);

  E_Int Ti = TRI_debug::get_T3_index(connect, TRUPLE);
  E_Int a = (Ti != IDX_NONE) ? ancestors[Ti] : IDX_NONE;
  }
#endif
  
#ifdef DEBUG_CONFORMIZER
  NUGA::ConformizerRoot::fastdiscard_counter = NUGA::ConformizerRoot::xtest_counter = 0;
  NUGA::ConformizerRoot::split_counter = NUGA::ConformizerRoot::split_fastdiscard_counter = NUGA::ConformizerRoot::degen_counter = 0;
#endif

#ifdef FLAG_STEP
   std::cout << "Conformizer : get x pairs : " << tt.elapsed() << std::endl;
   std::cout << "Conformizer : nb pts : " << pos.cols() << std::endl;
   std::cout << "Conformizer : nb_elts : " << connect.cols() << std::endl;
   std::cout << "Conformizer : xc size : " << xc.size() << std::endl;
   std::cout << "Conformizer : tol_x :  " << _tol_x << std::endl;

  tt.start();
#endif
  
  // Computes the intersection nodes and edges.
  _N0 = pos.cols(); // First index of new points.
  E_Int nbX = 0;
  if (!_brute_force)
    nbX=this->__compute_intersections_w_localizer(pos, connect, xc, _tol_x/*EPSILON for now*/);
  else
    nbX = this->__compute_intersections_brute(pos, connect, xc, _tol_x/*EPSILON for now*/);

#ifdef FLAG_STEP
   std::cout << "Conformizer : compute intersections : " << tt.elapsed() << std::endl;
  tt.start();
#endif
  
#ifdef DEBUG_CONFORMIZER
  if (detect_duplis_and_baffles(connect))
    std::cout << "oddities : after __compute_intersections : FOUND !!!!" << std::endl;
#endif
  
#ifdef DEBUG_CONFORMIZER
  _nbX = nbX;
  std::ostringstream o;
  o << "allX_" << _iter << ".mesh";
  //medith::write(o.str().c_str(), pos, connect, 0/*elt type*/, &xc);
  medith::write(o.str().c_str(), pos, connect, "TRI");
#endif
  

  
  if (nbX == 0) // Done.
    return 0;
  
#ifdef FLAG_STEP
  tt.start();
#endif

  // Merges the new nodes.
  std::vector<E_Int> nids;
  // fixme : tol here is set to 1.-9 : small to ensure convergence of the process (compute X / merge) but not too small to prevent meshing (splitElements) to fail...

#ifdef FLAG_STEP
  E_Int nb_merges =
#endif
    int n0 = (_X1 != 0) ? -_N0 : _N0; 
    __merge_clean(_tol_clean/*1.e-9 for now*/, pos, connect, ancestors, nids, &xc, n0);

#ifdef DEBUG_CONFORMIZER
  {
    std::ostringstream o;
    o << "after_merge_" << _iter << ".mesh";
    medith::write(o.str().c_str(), pos, connect);
}
#endif

#ifdef FLAG_STEP
    std::cout << "inter : __merge_clean (clean after X - nb merged : " << nb_merges << " ) :" << tt.elapsed() << std::endl;
  tt.start();
#endif
  
#ifdef DEBUG_CONFORMIZER
  if (detect_duplis_and_baffles(connect))
    std::cout << "oddities : after __merge_clean : FOUND !!!!" << std::endl;
#endif
  
  // Propagate the new indices in the data AND TIDY EDGES.
  {
    this->__update_data(pos, connect, nids);

#ifdef FLAG_STEP
  std::cout << "inter: __update_data after cleaning:" << tt.elapsed() << std::endl;
  tt.start();
#endif
  }
  
#ifdef DEBUG_CONFORMIZER
  if (detect_duplis_and_baffles(connect))
    std::cout << "oddities: after __update_data: FOUND !!!!" << std::endl;
#endif

#ifdef E_DEBUG
  {
    std::ostringstream o;
    o << "BefSplit_" << _iter << ".mesh";
    K_FLD::IntArray xing;
    for (size_t i = 0; i < xc.size(); ++i)
      //if (xc[i])
        xing.pushBack(connect.col(i), connect.col(i)+3);
    medith::write(o.str().c_str(), pos, xing);
  }
#endif
  
#ifdef FLAG_STEP
  {
    tt.start();
     std::cout << "Conformizer : splitting " << connect.cols() << " elements ... " << std::endl;
    /*std::cout << " pos sz : " << pos.cols() << std::endl;
    std::cout << " connect sz : " << connect.cols() << std::endl;
    std::cout << " xc sz : " << xc.size() << std::endl;
    std::cout << " ancestors sz : " << ancestors.size() << std::endl;*/
#ifdef DEBUG_CONFORMIZER
    E_Int Ti = TRI_debug::get_T3_index(connect, TRUPLE);
    E_Int a = (Ti != IDX_NONE) ? ancestors[Ti] : IDX_NONE;
#endif
  }
#endif
#ifdef DEBUG_CONFORMIZER
  {
    std::ostringstream o;
    o << "before_splitX_" << _iter << ".mesh";
    medith::write(o.str().c_str(), pos, connect);
  }
#endif

  // Split the triangles.  
  E_Int err = __split_Elements (pos, connect, xc, ancestors); 
  if (err)
  {
#ifdef DEBUG_CONFORMIZER
    std::cout << "failed on split" << std::endl;
#endif
    return -1;
  }
  
#ifdef FLAG_STEP
  std::cout << "Conformizer : splitting elements : " << tt.elapsed() << std::endl;
  std::cout << "Conformizer : new connect size : " << connect.cols() << std::endl;
#endif
    
  
#ifdef DEBUG_CONFORMIZER
  if (detect_duplis_and_baffles(connect))
    std::cout << "oddities : after __split_Elements : FOUND !!!!" << std::endl;
#endif

#ifdef DEBUG_CONFORMIZER
  if (nbX)
  {
    std::ostringstream o;
    o << "splitX_" << _iter << ".mesh";
    medith::write(o.str().c_str(), pos, connect);
  }
#endif

#ifdef FLAG_STEP
  tt.start();
#endif
  
#ifdef DEBUG_TRI_CONFORMIZER
  E_Int Ti = TRI_debug::get_T3_index(connect, TRUPLE);
  E_Int a = (Ti != IDX_NONE) ? ancestors[Ti] : IDX_NONE;
  //E_Int b = ancestors[5416];
  //std::cout << "before : " << 5416 << "->" << ancestors[5416] << std::endl;
  //std::cout << "before : " << 5595 << "->" << ancestors[5595] << std::endl;
#endif
  
  //fixme : improvement ?
  // might be a better place to update the tolerance but require a finer tolerance
  //based on hard node concept and local metric preservation to avoid triangles to get too big when some of their nodes are merged on the contour
  //if (_iter == 1 && tolerance == 0.) // First pass.
  //   __update_tolerance(pos, xc, ancestors);

#ifdef DEBUG_CONFORMIZER
  E_Int nb_merges2 = 0;
#endif
  if (_one_pass_mode && _split_swap_afterwards)
#ifdef DEBUG_CONFORMIZER
    nb_merges2 = 
#endif
    this->__simplify_and_clean2(pos, EPSILON, connect, ancestors, xc);
  else if (!_one_pass_mode)
#ifdef DEBUG_CONFORMIZER
    nb_merges2 = 
#endif
    this->__simplify_and_clean(pos, _tolerance, connect, ancestors, xc);

#ifdef DEBUG_CONFORMIZER
  if (detect_duplis_and_baffles(connect))
    std::cout << "oddities : after __simplify_and_clean : FOUND !!!!" << std::endl;
#endif

  if (connect.cols() == 0) // Error.
    nbX = 0;
  
#ifdef DEBUG_CONFORMIZER
  if (nb_merges2)
  {
    std::ostringstream o;
    o << "spe_mergeX_" << _iter << ".mesh";
    medith::write(o.str().c_str(), pos, connect, 0/*elt type*/, &xc);
  }
#endif

  // Clean unused nodes. fixme : see if it can be only once at the end.
  if (nbX)
  {
    std::vector<E_Int> nids;
    NUGA::MeshTool::compact_to_mesh(pos, connect, nids);

    if (_with_history)
    {
      K_CONNECT::IdTool::propagate(nids, _node_history);
      //update also parent_t::_elements to get x history (from which edge this x point commes from ?)
      this->__update_data(pos, connect, nids);
    }
    K_FLD::IntArray::changeIndices(_connect0, nids);
  }

#ifdef FLAG_STEP
   std::cout << "Conformizer : compact : " << tt.elapsed() << std::endl;
  tt.start();
#endif

#ifdef DEBUG_CONFORMIZER
  {
    std::vector<E_Int> colors;
    for (E_Int c = 0; c < connect.cols(); ++c)
      colors.push_back(ancestors[c]);
    
    std::ostringstream o;
    o << "connecto_" << _iter << ".mesh";
    medith::write(o.str().c_str(), pos, connect);
  }
#endif
  
#ifdef DEBUG_TRI_CONFORMIZER
  if (nbX)
  {
    std::ostringstream o;
    o << "solvedX_" << _iter << ".mesh";
    medith::write(o.str().c_str(), pos, connect, 0/*elt type*/, &xc);
    E_Int Ti = TRI_debug::get_T3_index(connect,TRUPLE);
    E_Int a = (Ti != IDX_NONE) ? ancestors[Ti] : IDX_NONE;
  }
#endif

  return nbX; //should be 0 if OK
}

///
template <short DIM, typename Element_t>
void
Conformizer<DIM, Element_t>::__finalize()
{
  if (_pool != nullptr) 
  {
    delete [] _pool;
    _pool = nullptr;
  }

  _boxes.clear();

  if (_tree != nullptr)
  {
    delete _tree;
    _tree = nullptr;
  } 
}

///
template <short DIM, typename Element_t>
void
Conformizer<DIM, Element_t>::__update_tolerance(const K_FLD::FloatArray& coord, NUGA::bool_vector_type& xc, const std::vector<E_Int>& ancestors)
{
  NUGA::bool_vector_type xc0(_connect0.cols()+1, false);
  for (size_t i=0; i < xc.size(); ++i)
  {
    xc0[ancestors[i]] = xc[i];
  }
  //medith::write("X.mesh", pos, _connect0, 0/*elt type*/, &xc0);
  E_Float Lmin, Lmax;
  __compute_min_edge_length(coord, _connect0, xc0, Lmin, Lmax);
  E_Float new_tol = MAX_PERCENT*Lmin;
  _tolerance = std::max(new_tol, _tolerance);
  
#ifdef DEBUG_CONFORMIZER
  std::cout << "new tol: "  << _tolerance << std::endl;
#endif
  
}

///
template <short DIM, typename Element_t>
void
Conformizer<DIM, Element_t>::__compute_min_edge_length
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const NUGA::bool_vector_type& xc, E_Float &Lmin, E_Float& Lmax)
{
  E_Float  min_d(EPSILON), max_d;
  K_FLD::IntArray relevant_elts;
  E_Int ROWS(connect.rows());

  for (E_Int c = 0; c < connect.cols(); ++c)
  {
    if (xc[c])
      relevant_elts.pushBack(connect.col(c), connect.col(c)+ROWS);
  }
  
  NUGA::MeshTool::computeMinMaxEdgeSqrLength<DIM>(pos, relevant_elts, min_d, max_d);

  Lmin= sqrt(min_d);
  Lmax = sqrt(max_d);
}

///
template <short DIM, typename Element_t>
E_Int
Conformizer<DIM, Element_t>::__compute_intersections_w_localizer
(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, NUGA::bool_vector_type& xc, E_Float tolerance)
{
  E_Int i, j, nb_boxes, nb_elts(xc.size()), nbX(0)/*, v1(-1), v2(-1), nbO(0)*/;
  std::vector<E_Int> bs;
  E_Int x;
  std::set<K_MESH::NO_Edge> processed;
  NUGA::bool_vector_type new_xc(xc.size(), false);
  
  typedef typename Struc<Element_t>::Type DS_Type;
  
#ifdef DEBUG_TRI_CONFORMIZER
  /*E_Int counter(0);
  E_Int Ti = TRI_debug::get_T3_index(connect, 17911, 17912, 18156);
  std::cout << "Ti : " << Ti << std::endl;
  E_Int Ti2 = TRI_debug::get_T3_index(connect,17882, 17911, 17849);
  std::cout << "Ti2 : " << Ti2 << std::endl;*/
#endif
  
#ifdef DEBUG_CONFORMIZER
  std::vector<DS_Type> xelts, aelts;
  E_Int iTi=0;
  std::ostringstream o;
  E_Int count=0;
#endif

  //
  for (i = 0; i < nb_elts; ++i)
  {
    if (!xc[i])
      continue;

    DS_Type& e1 = _elements[i];

    // Get surface2's boxes that are overlapping e1.
    bs.clear();
    _tree->getOverlappingBoxes(_boxes[i]->minB, _boxes[i]->maxB, bs);
    nb_boxes = bs.size();
    // Loop through the caught triangles and do the intersection test.
    for (j = 0; j < nb_boxes; ++j)
    {
      if (bs[j] == i)
        continue;
      
      if (!processed.insert(K_MESH::NO_Edge(i, bs[j])).second)
        continue;
      
      DS_Type& e2 = _elements[bs[j]];
      
#ifdef DEBUG_CONFORMIZER
      if (e1.id == zTi)iTi=i;
      else if (e2.id == zTi)iTi=bs[j];
      E_Int jj = bs[j];
      
      if (e1.id == zTi)
        aelts.push_back(e2);
      else if (e2.id == zTi)
        aelts.push_back(e1);
      
      //if ((i == Ti && j == 11) || (i == 82624 && j == 25))
      if (e1.id == zTi || e2.id == zTi)
      {
        std::vector<DS_Type> selec;
        selec.push_back(e1);
        selec.push_back(e2);
        o.str("");
        o << "selec_" << count++ << ".mesh";
        this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, selec, false);
      }
      
#endif

#ifdef DEBUG_CONFORMIZER
    NUGA::ConformizerRoot::xtest_counter++;
#endif

    x = this->__intersect(pos, connect, e1, e2, tolerance); // intersection test and trace storage

    if (x)
    {
      ++nbX;
      new_xc[i] = new_xc[bs[j]] = true;
      //      v1 = e1.id;
      //      v2 = e2.id;

#ifdef DEBUG_CONFORMIZER
      if (e1.id == zTi || e2.id == zTi)
      {
        if (e1.id != zTi) xelts.push_back(e1);
        if (e2.id != zTi) xelts.push_back(e2);
      }
#endif

      if (x == 2)
      {
        if (!_one_pass_mode)_needs_another_iter = true; // overlap might require a second pass as each triangle split is done separtely so triangulation can be different.
        if (_one_pass_mode && _iter == 1)_xpairs.insert(K_MESH::NO_Edge(i, bs[j]));
      }
    }
    }
  }

  xc = new_xc;
  
#ifdef FLAG_STEP
   std::cout << "cloud size (nb of T3s)   : " << connect.cols() << std::endl;
#ifdef DEBUG_CONFORMIZER
   std::cout << "total nb of X tests done : " << NUGA::ConformizerRoot::xtest_counter << std::endl;
   std::cout << "nb of X quicly discarded : " << NUGA::ConformizerRoot::fastdiscard_counter << std::endl;
   std::cout << " REDUCING NB OF TESTS ??? : " << xtestreal << " instead of " << xtesttotal << std::endl;
#endif
   std::cout << "nb of X detected         : " << nbX << std::endl;
#endif

#ifdef DEBUG_CONFORMIZER
  
  o.str("");
  o << "Xmolecule_" << zTi << ".mesh";
  xelts.push_back(iTi);
  xelts.push_back(zTi);
  this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, xelts);
  drawT3(pos, connect, zTi, true);
  
  o.str("");
  o << "molecule_" << zTi << ".mesh";
  std::vector<E_Int> colors(aelts.size(), 0);
  colors.push_back(1);
  aelts.push_back(zTi);
  if (aelts.size() > 1)
    this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, aelts, 0, &colors);
  
#endif
  
  return nbX;
}

///
template <short DIM, typename Element_t>
E_Int
Conformizer<DIM, Element_t>::__compute_intersections_brute
(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, NUGA::bool_vector_type& xc, E_Float tolerance)
{
  NUGA::bool_vector_type new_xc(xc.size(), false);

  typedef typename Struc<Element_t>::Type DS_Type;

#ifdef DEBUG_TRI_CONFORMIZER
  /*E_Int counter(0);
  E_Int Ti = TRI_debug::get_T3_index(connect, 17911, 17912, 18156);
  std::cout << "Ti : " << Ti << std::endl;
  E_Int Ti2 = TRI_debug::get_T3_index(connect,17882, 17911, 17849);
  std::cout << "Ti2 : " << Ti2 << std::endl;*/
#endif

#ifdef DEBUG_CONFORMIZER
  std::vector<DS_Type> xelts, aelts;
  E_Int iTi = 0;
  std::ostringstream o;
  E_Int count = 0;
#endif

  E_Int nbelts = _elements.size();
  E_Int nbX(0);
  //
  for (int i = 0; i < _X0; ++i)
  {
    DS_Type& e1 = _elements[i];

    // Loop through the second part and do the intersection test.
    for (int j = _X0; j < nbelts; ++j)
    {
      DS_Type& e2 = _elements[j];

#ifdef DEBUG_CONFORMIZER
      /*if (e1.id == zTi)iTi = i;
      else if (e2.id == zTi)iTi = j;

      if (e1.id == zTi)
        aelts.push_back(e2);
      else if (e2.id == zTi)
        aelts.push_back(e1);*/

      //if ((i == Ti && j == 11) || (i == 82624 && j == 25))
      //if (e1.id == zTi || e2.id == zTi)
      {
        std::vector<DS_Type> selec;
        selec.push_back(e1);
        selec.push_back(e2);
        
        o.str("");
        o << "selec_" << count++ << ".mesh";
        this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, selec, false);
      }

#endif

#ifdef DEBUG_CONFORMIZER
      NUGA::ConformizerRoot::xtest_counter++;
#endif

      E_Int x = this->__intersect(pos, connect, e1, e2, tolerance); // intersection test and trace storage

      if (x)
      {
        ++nbX;
        new_xc[i] = new_xc[j] = true;
        //      v1 = e1.id;
        //      v2 = e2.id;

#ifdef DEBUG_CONFORMIZER
        if (e1.id == zTi || e2.id == zTi)
        {
          if (e1.id != zTi) xelts.push_back(e1);
          if (e2.id != zTi) xelts.push_back(e2);
        }
#endif

        if (x == 2)
        {
          if (!_one_pass_mode)_needs_another_iter = true; // overlap might require a second pass as each triangle split is done separtely so triangulation can be different.
          if (_one_pass_mode && _iter == 1)_xpairs.insert(K_MESH::NO_Edge(i, j));
        }
      }
    }
  }

  xc = new_xc;

#ifdef FLAG_STEP
   std::cout << "cloud size (nb of T3s)   : " << connect.cols() << std::endl;
#ifdef DEBUG_CONFORMIZER
   std::cout << "total nb of X tests done : " << NUGA::ConformizerRoot::xtest_counter << std::endl;
   std::cout << "nb of X quicly discarded : " << NUGA::ConformizerRoot::fastdiscard_counter << std::endl;
   std::cout << " REDUCING NB OF TESTS ??? : " << xtestreal << " instead of " << xtesttotal << std::endl;
#endif
   std::cout << "nb of X detected         : " << nbX << std::endl;
#endif

#ifdef DEBUG_CONFORMIZER

  /*o.str("");
  o << "Xmolecule_" << zTi << ".mesh";
  xelts.push_back(iTi);
  xelts.push_back(zTi);
  this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, xelts);
  drawT3(pos, connect, zTi, true);

  o.str("");
  o << "molecule_" << zTi << ".mesh";
  std::vector<E_Int> colors(aelts.size(), 0);
  colors.push_back(1);
  aelts.push_back(zTi);
  if (aelts.size() > 1)
    this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, aelts, 0, &colors);*/

#endif

  return nbX;
}

}

#endif
