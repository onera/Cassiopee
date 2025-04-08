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

#ifndef __KCORE_SEARCH_BBTREE_CXX__
#define __KCORE_SEARCH_BBTREE_CXX__

#include "Nuga/include/BbTree.h"
#include <algorithm>
#include "Nuga/include/ngon_t.hxx"

using ngon_type = ngon_t<K_FLD::IntArray>;

/// Builds a tree and inserts the boxes (externally built) from begin to end.
template <short DIM, typename BBoxType>
K_SEARCH::BbTree<DIM, BBoxType>::BbTree
(const std::vector<BBoxType*> & boxes, E_Float tolerance)
:_boxes(boxes), _tolerance(tolerance), _owes_boxes(false), _pool(NULL)
{

  size_type none = IDX_NONE, size = boxes.size();

  _root_id = (size > 1) ? size : 0;
  
  //
  _tree.resize(BBTREE_ROWS, size, &none); //fixme 
  _tree.reserve(BBTREE_ROWS, 2*size); //fixme
  
  std::vector<size_type> indices(size);
  for (size_type i = 0; i < size; ++i) indices[i] = i;

  __insert(indices.begin(), indices.end());

#ifdef E_TIME1
//_append_tree = _get_box_boxes = _get_longest_axis = _build_vectors = _fill_vectors = 0.;
#endif
}

namespace K_SEARCH
{
/// Builds a tree and inserts the boxes (externally built) from begin to end.
template <>
template <typename array_t>
BbTree<3>::BbTree
(const K_FLD::FloatArray& crd, const array_t& ng, E_Float tolerance)
: _tolerance(tolerance), _owes_boxes(true), _pool(NULL)
{
    
  E_Int nb_elts = ng.PHs.size();
  _pool = new K_SEARCH::BBox3D[nb_elts];

  _boxes.reserve(nb_elts);
  
  Vector_t<E_Int> nodes;

  for (int i = 0; i < nb_elts; i++){
    ng.nodes_ph(i, nodes, true);
    K_SEARCH::BBox3D* box = &_pool[i];
    box->compute(crd, nodes);
    _boxes.push_back(box);
  }
  
  size_type none = IDX_NONE, size = _boxes.size();

  _root_id = (size > 1) ? size : 0;
  
  //
  _tree.resize(BBTREE_ROWS, size, &none); //fixme 
  _tree.reserve(BBTREE_ROWS, 2*size); //fixme
  
  std::vector<size_type> indices(size);
  for (size_type i = 0; i < size; ++i) indices[i] = i;

  __insert(indices.begin(), indices.end());
}

template <> inline
BbTree<3>::BbTree
(const K_FLD::FloatArray& crd, const ngon_unit& PGs, E_Float tolerance)
: _tolerance(tolerance), _owes_boxes(true), _pool(NULL)
{
    
  E_Int nb_elts = PGs.size();
  _pool = new K_SEARCH::BBox3D[nb_elts];

  _boxes.reserve(nb_elts);

  for (int i = 0; i < nb_elts; i++){
    
    const E_Int* nodes = PGs.get_facets_ptr(i);
    E_Int stride = PGs.stride(i);
    K_SEARCH::BBox3D* box = &_pool[i];
    box->compute(crd, nodes, stride, 1);
    _boxes.push_back(box);
  }
  
  size_type none = IDX_NONE, size = _boxes.size();

  _root_id = (size > 1) ? size : 0;
  
  //
  _tree.resize(BBTREE_ROWS, size, &none); //fixme 
  _tree.reserve(BBTREE_ROWS, 2*size); //fixme
  
  std::vector<size_type> indices(size);
  for (size_type i = 0; i < size; ++i) indices[i] = i;

  __insert(indices.begin(), indices.end());
}

template <> inline
BbTree<3>::BbTree
(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt)
: _tolerance(EPSILON), _owes_boxes(true), _pool(NULL)
{
  E_Int nb_elts = cnt.cols();
  _pool = new K_SEARCH::BBox3D[nb_elts];

  _boxes.reserve(nb_elts);

  E_Int stride = cnt.rows();

  for (int i = 0; i < nb_elts; i++){
    
    const E_Int* nodes = cnt.col(i);
    
    K_SEARCH::BBox3D* box = &_pool[i];
    box->compute(crd, nodes, stride, 0);
    _boxes.push_back(box);
  }
  
  size_type none = IDX_NONE, size = _boxes.size();

  _root_id = (size > 1) ? size : 0;
  
  //
  _tree.resize(BBTREE_ROWS, size, &none); //fixme 
  _tree.reserve(BBTREE_ROWS, 2*size); //fixme
  
  std::vector<size_type> indices(size);
  for (size_type i = 0; i < size; ++i) indices[i] = i;

  __insert(indices.begin(), indices.end());
}

}

///
template <short DIM, typename BBoxType>
void
K_SEARCH::BbTree<DIM, BBoxType>::getOverlappingBoxes
(const E_Float* minB, const E_Float* maxB, std::vector<size_type>& out) const
{
  BBoxType box(minB, maxB);
  __getOverlappingBoxes(&box, _root_id, out);
}

///
template <short DIM, typename BBoxType>
void
K_SEARCH::BbTree<DIM, BBoxType>::__getOverlappingBoxes
(BBoxType* box, E_Int node, std::vector<size_type>& out) const
{
  E_Int leftc(_tree(0, node)), rightc(_tree(1, node));
  
  // Leaf node
  if (leftc == IDX_NONE)
  {
    if (boxesAreOverlapping(_boxes[node], box, _tolerance))
      out.push_back(node);
    return;
  }

  // Left child
  if (boxesAreOverlapping(_boxes[leftc], box, _tolerance))
    __getOverlappingBoxes(box, leftc, out);

  // Right child
  if (boxesAreOverlapping(_boxes[rightc], box, _tolerance))
    __getOverlappingBoxes(box, rightc, out);
}

///============================================================================
template <short DIM, typename BBoxType>
bool
K_SEARCH::BbTree<DIM, BBoxType>::hasAnOverlappingBox
(const E_Float* minB, const E_Float* maxB) const
{
  BBoxType box(minB, maxB);
  return __hasAnOverlappingBox(&box, _root_id);
}

///
template <short DIM, typename BBoxType>
bool
K_SEARCH::BbTree<DIM, BBoxType>::__hasAnOverlappingBox
(BBoxType* box, E_Int node) const
{
  E_Int leftc(_tree(0, node)), rightc(_tree(1, node));
  
  // Leaf node
  if (leftc == IDX_NONE)
    return boxesAreOverlapping(_boxes[node], box, _tolerance);

  // Left child
  bool found = false;
  if (boxesAreOverlapping(_boxes[leftc], box, _tolerance))
    found = __hasAnOverlappingBox(box, leftc);

  // Right child
  if (!found && boxesAreOverlapping(_boxes[rightc], box, _tolerance))
    found = __hasAnOverlappingBox(box, rightc);
  
  return found;
}

//=============================================================================
template <short DIM, typename BBoxType>
void
K_SEARCH::BbTree<DIM, BBoxType>::getIntersectingBoxes
(const E_Float* P0, const E_Float* P1, std::vector<size_type>& out, E_Float tolerance, bool strict) const
{
  if (!strict)
    __getIntersectingBoxes(P0, P1, _root_id, tolerance, out);
  else
  {
    BBox3D boxSeg;
    for (size_t i=0; i < DIM; ++i)
    {
      boxSeg.minB[i]=std::min(P0[i], P1[i]);
      boxSeg.maxB[i]=std::max(P0[i], P1[i]);
    }
    __getIntersectingBoxesStrict(P0, P1, _root_id, tolerance, out, boxSeg);
  }
}

//=============================================================================
template <short DIM, typename BBoxType>
void
K_SEARCH::BbTree<DIM, BBoxType>::__getIntersectingBoxes
(const E_Float* P0, const E_Float* P1, E_Int node, const E_Float& abstol, std::vector<size_type>& out) const
{
  E_Int leftc(_tree(0, node)), rightc(_tree(1, node));

  // Leaf node
  if (leftc == IDX_NONE)
  {
    if (__boxIntersectRay(*_boxes[node], P0, P1, abstol))
      out.push_back(node);
    return;
  }
  
  // Left child
  if (__boxIntersectRay(*_boxes[leftc], P0, P1, abstol))
    __getIntersectingBoxes(P0, P1, leftc, abstol, out);

  // Right child
  if (__boxIntersectRay(*_boxes[rightc], P0, P1, abstol))
    __getIntersectingBoxes(P0, P1, rightc, abstol, out);
}

//=============================================================================
template <short DIM, typename BBoxType>
void
K_SEARCH::BbTree<DIM, BBoxType>::__getIntersectingBoxesStrict
(const E_Float* P0, const E_Float* P1, E_Int node, const E_Float& abstol, std::vector<size_type>& out, const BBoxType & boxSeg) const
{
  E_Int leftc(_tree(0, node)), rightc(_tree(1, node));

  // Leaf node
  if (leftc == IDX_NONE)
  {
    if (__boxIntersectSeg(*_boxes[node], P0, P1, abstol, boxSeg))
      out.push_back(node);
    return;
  }
  
  // Left child
  if (__boxIntersectSeg(*_boxes[leftc], P0, P1, abstol, boxSeg))
    __getIntersectingBoxesStrict(P0, P1, leftc, abstol, out, boxSeg);

  // Right child
  if (__boxIntersectSeg(*_boxes[rightc], P0, P1, abstol, boxSeg))
    __getIntersectingBoxesStrict(P0, P1, rightc, abstol, out, boxSeg);
}


///
template <short DIM, typename BBoxType>
E_Int
K_SEARCH::BbTree<DIM, BBoxType>::__insert
(const std::vector<size_type>::const_iterator& begin,
 const std::vector<size_type>::const_iterator& end)
{

#ifdef E_TIME1
  NUGA::chrono c1;
#endif

  size_type size = end - begin;
  if (size == 1)
    return *begin;

#ifdef E_TIME1
  c1.start();
#endif

  size_type col = _tree.cols();
  _tree.pushBack(_newL, _newL + BBTREE_ROWS);

#ifdef E_TIME1
  _append_tree += c1.elapsed();
  c1.start();
#endif

  BBoxType* box = new BBoxType;
  _boxes.push_back(box);
  __getBoxesBox(_boxes, begin, end, *box);

#ifdef E_TIME1
  _get_box_boxes += c1.elapsed();
#endif

  if (size == 2)
  {
    _tree(0, col) = *begin;
    _tree(1, col) = *(begin+1);
    return col;
  }


#ifdef E_TIME1
  c1.start();
#endif

  size_type axis = __getLongestSideAxis(*box);

#ifdef E_TIME1
  _get_longest_axis += c1.elapsed();
#endif

  E_Float      x = (box->getMinB(axis) + box->getMaxB(axis));

#ifdef E_TIME1
  c1.start();
#endif

  std::vector<size_type> left, right;
  left.reserve(size);
  right.reserve(size);
  std::vector<size_type>::const_iterator it;

#ifdef E_TIME1
  _build_vectors += c1.elapsed();
  c1.start();
#endif

  for (it = begin; it != end; ++it)
  {
    const E_Int& id = *it;
    const BBoxType* box = _boxes[id];

    if ((box->getMinB(axis) + box->getMaxB(axis)) < x)
      left.push_back(id);
    else
      right.push_back(id);
  }

  std::vector<size_type>::const_iterator r_end(right.end()), l_end(left.end()), r_beg(right.begin()), l_beg(left.begin());

  if (r_end == r_beg)
  {
    size_type median = ((l_end - l_beg)/2);
    //it = l_beg + median;
    l_end = r_beg = l_beg + median;
    //l_beg = left.begin();
    r_end = left.end();
  }
  else if (l_end == l_beg) 
  {
    size_type median = (r_end - r_beg) / 2;
    //it = r_beg + median;
    l_end = r_beg = r_beg + median;
    l_beg = right.begin();
    //r_end = right.end();
  }

#ifdef E_TIME1
  _fill_vectors += c1.elapsed();
#endif

  _tree(0, col) = __insert(l_beg, l_end);
  _tree(1, col) = __insert(r_beg, r_end);

  return col;
}

///
template <short DIM, typename BBoxType>
void
K_SEARCH::BbTree<DIM, BBoxType>::__getBoxesBox
(const std::vector<BBoxType*>& boxes, const std::vector<size_type>::const_iterator& begin,
 const std::vector<size_type>::const_iterator& end, BBoxType& box) const
{
  for (std::vector<size_type>::const_iterator it = begin; it != end; ++it)
  {
    const BBoxType* Bi = boxes[*it];

    for (E_Int i = 0; i < DIM; ++i)
    {
      box.minB[i] = (box.minB[i] > Bi->getMinB(i)) ? Bi->getMinB(i) : box.minB[i];
      box.maxB[i] = (box.maxB[i] < Bi->getMaxB(i)) ? Bi->getMaxB(i) : box.maxB[i];
    }
  }
}

///
template <short DIM, typename BBoxType>
E_Int
K_SEARCH::BbTree<DIM, BBoxType>::__getLongestSideAxis
(const BBoxType& box) const
{
  E_Int axis(0);
  E_Float L = -NUGA::FLOAT_MAX, Li;
  for (E_Int i = 0; i < DIM; ++i)
  {
    Li = box.maxB[i] - box.minB[i];
    if (L < Li)
    {
      axis = i;
      L = Li;
    }
  }
  return axis;
}

template <short DIM, typename BBoxType>
E_Bool
K_SEARCH::BbTree<DIM, BBoxType>::boxesAreOverlapping (const BBoxType * bb1, const BBoxType * bb2, const E_Float& tol) 
{
  if (bb1 == bb2) return true;

  for (E_Int i = 0; i < DIM; ++i)
  {
    if (bb1->minB[i] > bb2->maxB[i] + tol) return false;
    if (bb2->minB[i] > bb1->maxB[i] + tol) return false;
  }  
  return true;
}

template <short DIM, typename BBoxType>
E_Bool
K_SEARCH::BbTree<DIM, BBoxType>::box1IsIncludedinbox2 (const BBoxType * bb1, const BBoxType * bb2, const E_Float& tol) 
{
  if (bb1 == bb2) return false;

  for (E_Int i = 0; i < DIM; ++i)
  {
    if (bb1->minB[i] <= bb2->minB[i] - tol) return false;
    if (bb2->maxB[i] <= bb1->maxB[i] - tol) return false;
  }  
  return true;
}

namespace K_SEARCH
{
///
template <> inline
E_Bool
BbTree<2>::__boxIntersectRay
(const BBox2D & box, const E_Float* P0, const E_Float* P1, const E_Float& abstol) const
{
  return false;//fixme
}

///
inline void __compute_x0_x1(const E_Float* P0, const E_Float* P0P1, const E_Float* mB, E_Int ind1, E_Int ind2, E_Float lambda, E_Float& x0, E_Float& x1)
{
  x0 = (P0[ind1] + lambda * P0P1[ind1]) - mB[ind1];
  x1 = (P0[ind2] + lambda * P0P1[ind2]) - mB[ind2];
}

///
template <> inline
E_Bool
BbTree<3>::__boxIntersectRay
(const BBox3D & box, const E_Float* P0, const E_Float* P1, const E_Float& abstol) const
{
  E_Float P0P1[3];
  NUGA::diff<3>(P1, P0, P0P1);
  
  E_Float denom, lambda, x0, x1;
  E_Float dx(box.maxB[0] - box.minB[0]), dy(box.maxB[1] - box.minB[1]), dz(box.maxB[2] - box.minB[2]);

  // Front and back faces
  denom = P0P1[1];
  if ( (EPSILON < denom) || (denom < - EPSILON) ) 
  {
    denom = 1./denom;
    
    lambda = (box.minB[1] - P0[1]) * denom;
    __compute_x0_x1(P0, P0P1, box.minB, 0, 2, lambda, x0, x1);
    
    if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
      return true;

    
    lambda = (box.maxB[1] - P0[1]) * denom;
    __compute_x0_x1(P0, P0P1, box.minB, 0, 2, lambda, x0, x1);
    
    if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
      return true;
  }

  // Right and left faces
  denom = P0P1[0];
  if ( (EPSILON < denom) || (denom < - EPSILON) ) 
  {
    denom = 1./denom;
    
    lambda = (box.minB[0] - P0[0]) * denom;
    __compute_x0_x1(P0, P0P1, box.minB, 1, 2, lambda, x0, x1);
    
    if ((x0 > -abstol) && (x0 < dy + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
      return true;

    lambda = (box.maxB[0] - P0[0]) * denom;
    __compute_x0_x1(P0, P0P1, box.minB, 1, 2, lambda, x0, x1);
    
    if ((x0 > -abstol) && (x0 < dy + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
      return true;
  }

  // Top and bottom faces
  denom = P0P1[2];
  if ( (EPSILON < denom) || (denom < - EPSILON) ) 
  {
    denom = 1./denom;

    lambda = (box.minB[2] - P0[2]) * denom;
    __compute_x0_x1(P0, P0P1, box.minB, 0, 1, lambda, x0, x1);

    if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dy + abstol))
      return true;

    lambda = (box.maxB[2] - P0[2]) * denom;
    __compute_x0_x1(P0, P0P1, box.minB, 0, 1, lambda, x0, x1);

    if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dy + abstol))
      return true;
  }

  return false;
}

///
template <> inline
E_Bool/*fixme*/
BbTree<3>::__boxIntersectSeg
(const BBox3D & box, const E_Float* P0, const E_Float* P1, const E_Float& abstol, const BBox3D & boxSeg) const
{
  // Fast return : box are not overlapping
  if (!boxesAreOverlapping(&box, &boxSeg, abstol))
    return false;
  
  // Fast return : seg fully inside box
  bool ins=true;
  for (E_Int i = 0; i < 3; ++i)
  {
    if (box.minB[i] > boxSeg.minB[i] + abstol) {ins=false;break;}
    if (boxSeg.maxB[i] > box.maxB[i] + abstol) {ins=false;break;}
  }  
  if (ins)
    return 1;
  
  E_Float P0P1[3];
  NUGA::diff<3>(P1, P0, P0P1);
  
  E_Float denom, lambda, x0, x1;
  E_Float dx(box.maxB[0] - box.minB[0]), dy(box.maxB[1] - box.minB[1]), dz(box.maxB[2] - box.minB[2]);

  // Front and back faces
  denom = P0P1[1];
  if ((EPSILON < denom) || (denom < -EPSILON)) // fixme : iconsistency with abstol (several similar test below)
  {
    denom = 1./ denom;
    
    lambda = (box.minB[1] - P0[1]) * denom;
    if ( (-EPSILON < lambda) && (lambda < 1. + EPSILON) )
    {
      __compute_x0_x1(P0, P0P1, box.minB, 0, 2, lambda, x0, x1);
      
      if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
        return 1;
    }

    lambda = (box.maxB[1] - P0[1]) * denom;
    if ( (-EPSILON < lambda) && (lambda < 1. + EPSILON) )
    {
      __compute_x0_x1(P0, P0P1, box.minB, 0, 2, lambda, x0, x1);

      if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
        return 1;
    }
  }

  // Right and left faces
  denom = P0P1[0];
  if ((EPSILON < denom) || (denom < -EPSILON))
  {
    denom = 1./ denom;
    
    lambda = (box.minB[0] - P0[0]) * denom;
    if ( (-EPSILON < lambda) && (lambda < 1. + EPSILON) )
    {
      __compute_x0_x1(P0, P0P1, box.minB, 1, 2, lambda, x0, x1);
      
      if ((x0 > -abstol) && (x0 < dy + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
        return 1;
    }

    lambda = (box.maxB[0] - P0[0])*denom;
    if ( (-EPSILON < lambda) && (lambda < 1. + EPSILON) )
    {
      __compute_x0_x1(P0, P0P1, box.minB, 1, 2, lambda, x0, x1);
      
      if ((x0 > -abstol) && (x0 < dy + abstol) && (x1 > -abstol) && (x1 < dz + abstol))
        return 1;
    }
  }

  // Top and bottom faces
  denom = P0P1[2];
  if ((EPSILON < denom) || (denom < -EPSILON))
  {
    denom = 1./ denom;
    lambda = (box.minB[2] - P0[2]) * denom;
    if ( (-EPSILON < lambda) && (lambda < 1. + EPSILON) )
    {
      __compute_x0_x1(P0, P0P1, box.minB, 0, 1, lambda, x0, x1);
      
      if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dy + abstol))
        return 1;
    }

    lambda = (box.maxB[2] - P0[2]) * denom;
    if ( (-EPSILON < lambda) && (lambda < 1. + EPSILON) )
    {
      __compute_x0_x1(P0, P0P1, box.minB, 0, 1, lambda, x0, x1);
      
      if ((x0 > -abstol) && (x0 < dx + abstol) && (x1 > -abstol) && (x1 < dy + abstol))
        return 1;
    }
  }

  return 0;
}
}


#endif
