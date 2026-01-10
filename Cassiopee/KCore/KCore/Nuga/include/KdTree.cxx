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

#ifndef __KCORE_SEARCH_KDTREE_CXX__
#define __KCORE_SEARCH_KDTREE_CXX__

#include "Nuga/include/KdTree.h"

// ============================================================================
/// Builds a tree and inserts all the nodes of the input coordinate array pos.
// ============================================================================
template <typename CoordArrayType>
K_SEARCH::KdTree<CoordArrayType>::KdTree(const coord_access_type& posAcc,
                                         E_Float tolerance, bool do_omp)
:_posAcc(posAcc), _tree_sz(0), _dim(posAcc.stride()), _tolerance(tolerance*tolerance)
{
  size_type none = IDX_NONE;
  _tree.resize(3, posAcc.size(), &none);

  typedef std::vector<size_type> Vector;

  size_type val = 0;
  Vector indices(posAcc.size());
  Vector::iterator it(indices.begin()), itEnd(indices.end());

  while (it != itEnd) *(it++) = val++;

  if (!do_omp){
    //std::cout << "MERGE SEQ 1" << std::endl;
    __insert(indices.begin(), itEnd, 0/*depth*/);
  }
  else
  {
    //std::cout << "MERGE OMP1" << std::endl;
    //omp_set_nested(1);
#pragma omp parallel
    {
      #pragma omp single
      {
        __insert_omp(indices.begin(), itEnd, 0/*depth*/);
        __set_tree(indices.begin(), indices.end());
      }
    } 
  }
}

// ============================================================================
/// Builds a tree and inserts all the points of the input coordinates array.
// ============================================================================
template <typename CoordArrayType>
K_SEARCH::KdTree<CoordArrayType>::KdTree(const coord_access_type& posAcc, 
                                         std::vector<size_type> indices/*passed by value*/,
                                         E_Float tolerance, bool do_omp)
 :_posAcc(posAcc), _tree_sz(0), _dim(posAcc.stride()), _tolerance(tolerance*tolerance)
{
  size_type none = IDX_NONE;
  _tree.resize(3, _tree_sz + indices.size(), &none);

  if (!do_omp){
    //std::cout << "MERGE SEQ 2" << std::endl;
    __insert(indices.begin(), indices.end(), 0/*depth*/);
  }
  else
  {
    //std::cout << "MERGE OMP 2" << std::endl;
#pragma omp parallel
    {

      #pragma omp single
      {
        __insert_omp(indices.begin(), indices.end(), 0/*depth*/);
        __set_tree(indices.begin(), indices.end());
      }
    }
  }
}

// ============================================================================
/// Reset.
// ============================================================================
template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::clear()
{
  _tree_sz   = 0;
  _dim       = _posAcc.stride();

  //size_type none = IDX_NONE;
  _tree.clear();
}

// ============================================================================
/// Refill an existing tree
// ============================================================================
template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::build(std::vector<size_type>* indices, E_Float tolerance)
{
  _tolerance = tolerance*tolerance;

  clear();

  size_type none = IDX_NONE;

  if (indices != nullptr)
  {
    _tree.resize(3, _tree_sz + indices->size(), &none);
    __insert(indices->begin(), indices->end(), 0/*depth*/);
  }
  else
  {
    _tree.resize(3, _posAcc.size(), &none);
    
    std::vector<size_type> indices(_posAcc.size());
    std::vector<size_type>::iterator it(indices.begin()), itEnd(indices.end());
    size_type val = 0;
    while (it != itEnd) *(it++) = val++;

    __insert(indices.begin(), itEnd, 0/*depth*/);
  }
}

//
template <typename CoordArrayType>
template <typename InputIterator>
E_Int
K_SEARCH::KdTree<CoordArrayType>::__insert
(InputIterator begin, InputIterator end, size_type depth)
{
  if (begin == end) return IDX_NONE;

  InputIterator                    it(begin + (end - begin)/2);
  size_type                        cols(_tree_sz);
  
  //if (_posAcc.isOutOfRange(*it)) return IDX_NONE;

  if (it != begin)
  {
    size_type d = depth % _dim;
    std::nth_element(begin, it, end, [this,d](E_Int Ni, E_Int Nj) {return  this->_posAcc.getVal(Ni, d) <  this->_posAcc.getVal(Nj, d);});
  }

  _tree(0, _tree_sz++) = *it; // Set N as the parent, increase tree size by one.

   if (begin != it)
     _tree(1,cols) = __insert (begin, it, depth+1);

   if (++it != end)
     _tree(2,cols) = __insert (it, end, depth+1);

   return cols;
}

//
template <typename CoordArrayType>
template <typename InputIterator>
void
K_SEARCH::KdTree<CoordArrayType>::__insert_omp
(InputIterator begin, InputIterator end, size_type depth)
{
  if (begin == end) return;

  InputIterator it(begin + (end - begin)/2);
  auto itn=it; ++itn;
  size_type dn = depth+1;
    
  if (it != begin)
  {
    size_type d = depth % _dim;
    std::nth_element(begin, it, end, [this,d](E_Int Ni, E_Int Nj) {return  this->_posAcc.getVal(Ni, d) <  this->_posAcc.getVal(Nj, d); });
  }

  #pragma omp task
  {
    __insert_omp (begin, it, dn);
  }

  //#pragma omp task : save a thread by using the current thread
  {
    __insert_omp (itn, end, dn);
  }

  #pragma omp taskwait
}

//
template <typename CoordArrayType>
template <typename InputIterator>
E_Int
K_SEARCH::KdTree<CoordArrayType>::__set_tree
(InputIterator begin, InputIterator end)
{
  if (begin == end) return IDX_NONE;

  InputIterator                    it(begin + (end - begin)/2);
  size_type                        cols(_tree_sz);

  _tree(0, _tree_sz++) = *it; // Set N as the parent, increase tree size by one.

   if (begin != it)
     _tree(1,cols) = __set_tree (begin, it);

   if (++it != end)
     _tree(2,cols) = __set_tree (it, end);

   return cols;
}

// ============================================================================
/// Inserts a node in the tree.
// ============================================================================
template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::insert (size_type N)
{
  size_type none = IDX_NONE;
  _tree.resize(3, _tree_sz+1, &none);
  __insert(N);
}

// ============================================================================
/// Returns a close node. It is not necessarily the closest one.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClose(const E_Float* point) const
{
  size_type m(IDX_NONE);
  E_Float dist2;
  __getClosest_through_path(point, m, dist2);
  return m;
}

// ============================================================================
/// Returns a close node and the square of the distance to it. It is not necessarily the closest one. 
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClose(const E_Float* point, E_Float& dist2) const
{
  size_type m(IDX_NONE);
  __getClosest_through_path(point, m, dist2);
  return m;
}

// ============================================================================
// Returns a close node to node of index N. It is not necessarily the 
//closest one.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClose (E_Int n) const
{
  size_type m(IDX_NONE);
  E_Float dist2;
  __getClosest_through_path(n, m, dist2);
  return m;
}

// ============================================================================
// Returns a close node to node N (it is not necessarily the closest one) and the distance between them.
// ============================================================================
template <typename CoordArrayType>
E_Int K_SEARCH::KdTree<CoordArrayType>::getClose (E_Int n, E_Float& dist2) const
{
  size_type m(IDX_NONE);
  E_Float Xn[3];
  _posAcc.getEntry(n, Xn);
  __getClosest_through_path(n, Xn, m, dist2);
  return m;
}

// ============================================================================
/// Returns the closest node.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClosest(const E_Float* point) const
{
  E_Float dist2;
  size_type m = getClose(point, dist2);
  if (_tolerance < dist2)
    __seek_closest(point, 0/*root col*/, 0/*axis*/, dist2, m);
  
  return m;
}

// ============================================================================
/// Returns the closest node and the square distance to it.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClosest(const E_Float* point, E_Float& dist2) const
{
  size_type m = getClose(point, dist2);
  if (_tolerance < dist2)
    __seek_closest(point, 0/*root col*/, 0/*axis*/, dist2, m);
  
  return m;
}

// ============================================================================
/// Returns the closest node.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClosest(E_Int n) const 
{
  E_Float dist2;
  size_type m = getClose(n, dist2);
  
  if (_tolerance < dist2)
  {
    E_Float Xn[3];
    _posAcc.getEntry(n, Xn);
    __seek_closest(n, Xn, 0/*root col*/, 0/*axis*/, dist2, m);
  }

  return m;
}

// ============================================================================
/// Returns the closest node and the square distance to it.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClosest(E_Int n, E_Float& dist2) const 
{
  size_type m = getClose(n, dist2);
  
  if (_tolerance < dist2)
  {
    E_Float Xn[3];
    _posAcc.getEntry(n, Xn);
    __seek_closest(n, Xn, 0/*root col*/, 0/*axis*/, dist2, m);
  }

  return m;
}

// ============================================================================
/// Returns the closest node and the square distance to it using an input guessed distance.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClosest(E_Int n, const E_Float& guessed_d2, E_Float& dist2) const 
{
  size_type m = IDX_NONE;  
  dist2 = guessed_d2;

  if (_tolerance < dist2){
    E_Float Xn[3];
    _posAcc.getEntry(n, Xn);
    __seek_closest(n, Xn, 0/*root col*/, 0/*axis*/, dist2, m);
  }

  return m;
}

// ============================================================================
/// Returns the closest node and the square distance to it using an input guessed distance.
// ============================================================================
template <typename CoordArrayType>
E_Int
K_SEARCH::KdTree<CoordArrayType>::getClosest(const E_Float* point, const E_Float& guessed_d2, E_Float& dist2) const 
{
  size_type m = IDX_NONE;  
  dist2 = guessed_d2;

  if (_tolerance*_tolerance < dist2)
    __seek_closest(point, 0/*root col*/, 0/*axis*/, dist2, m);

  return m;
}

// ============================================================================
/// Returns all the nodes in the input box by appending the vector 'out'.
// ============================================================================
template <typename CoordArrayType>
void
K_SEARCH::KdTree<CoordArrayType>::getInBox
(const E_Float* minB, const E_Float* maxB, std::vector<size_type>& out) const
{
  // Check that the box is not ill posed
  bool ok = true;
  for (size_type i = 0; i < _dim; ++i)
    ok &= (minB[i] < maxB[i]);
  if (!ok) return;

  //
  __getInBox(0/*root*/, 0/*axis*/, minB, maxB, out);
}

// ============================================================================
/// Returns all the nodes in the input sphere centered on C by appending the vector 'out'.
// ============================================================================
template <typename CoordArrayType>
void
K_SEARCH::KdTree<CoordArrayType>::getInSphere
(const E_Float* C, E_Float radius, std::vector<size_type>& out) const
{
  E_Float mB[3], MB[3];
  for (E_Int i = 0; i < _dim; ++i)
  {
    MB[i] = C[i]+radius;
    mB[i] = C[i]-radius;
  }

  //
  __getInBox(0/*root*/, 0/*axis*/, mB, MB, out);

  //
  E_Float d2;
  size_t sz = out.size();
  for (size_t i = 0; i < sz;)
  {
    size_type & Ni = out[i];      
    d2 = _posAcc.dist2(C, Ni);
    if (d2 < radius*radius)
    {
      ++i;
    }
    else Ni=out[--sz];
  }
  out.resize(sz);
}

// ============================================================================
/// Returns all the nodes in the input sphere centered on C (and their distance to it)
/// by appending the vector 'out'.
// ============================================================================
template <typename CoordArrayType>
void
K_SEARCH::KdTree<CoordArrayType>::getInSphere
(const E_Float* C, E_Float radius, std::vector<size_type>& out, std::vector<E_Float>& dist2) const
{
  E_Float mB[3], MB[3];
  for (size_t i = 0; i < _dim; ++i)
  {
    MB[i] = C[i]+radius;
    mB[i] = C[i]-radius;
  }

  //
  __getInBox(0/*root*/, 0/*axis*/, mB, MB, out);

  // Reduce caught nodes to the open sphere (C, radius)
  E_Float d2, r2(radius*radius);
  size_t sz = out.size();
  dist2.reserve(sz);
  for (size_t i = 0; i < sz;)
  {
    size_type & Ni = out[i];      
    d2 = _posAcc.dist2(C, Ni);
    if (d2 < r2)
    {
      dist2.push_back(d2);
      ++i;
    }
    else Ni = out[--sz];
  }
  out.resize(sz);
}

// ============================================================================
/// Returns all the nodes in the input sphere centered on node N by appending the vector 'out'.
// ============================================================================
template <typename CoordArrayType>
void
K_SEARCH::KdTree<CoordArrayType>::getInSphere
(E_Int N, E_Float radius, std::vector<size_type>& out) const
{

  E_Float Xn[3], mB[3], MB[3];
  _posAcc.getEntry(N, Xn);

  for (size_t i = 0; i < _dim; ++i)
  {
    MB[i] = Xn[i]+radius;
    mB[i] = Xn[i]-radius;
  }

  //
  __getInBox(0/*root*/, 0/*axis*/, mB, MB, out);

  // Reduce caught nodes to the open sphere (C, radius)
  E_Float d2, r2(radius*radius);
  size_t sz = out.size();
  for (size_t i = 0; i < sz;)
  {
    size_type & Ni = out[i];      
    d2 = _posAcc.dist2(N, Ni);
    if ((N != Ni) && (d2 < r2)) ++i;
    else Ni = out[--sz];
  }
  out.resize(sz);
}

// ============================================================================
/// Returns all the nodes in the input sphere centered on node N (and their distance to it)
/// by appending the vector 'out'.
// ============================================================================
template <typename CoordArrayType>
void
K_SEARCH::KdTree<CoordArrayType>::getInSphere
(E_Int N, E_Float radius, std::vector<size_type>& out, std::vector<E_Float>& dist2) const
{
  E_Float Xn[3], mB[3], MB[3];
  _posAcc.getEntry(N, Xn);

  for (E_Int i = 0; i < _dim; ++i)
  {
    MB[i] = Xn[i]+radius;
    mB[i] = Xn[i]-radius;
  }

  //
  __getInBox(0/*root*/, 0/*axis*/, mB, MB, out);

  // Reduce caught nodes to the open sphere (C, radius)
  E_Float d2;
  size_t sz = out.size();
  dist2.reserve(sz);
  for (size_t i = 0; i < sz;)
  {
    size_type & Ni = out[i];      
    d2 = _posAcc.dist2(Xn, Ni);
    if ((Ni != N) && (d2 < radius*radius))
    {
      dist2.push_back(d2);
      ++i;
    }
    else
      Ni=out[--sz];
  }
  out.resize(sz);
}

/*****************  private ***********************************/

template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::__insert(size_type n)
{
  size_type i(0), axis(0), parent(0), child(1), *tbegin(_tree.begin());
  const size_type *pi;

  _posAcc.getEntry(n,_Xn);

  while (i < _tree_sz)
  {
    pi = tbegin + 3*i; // pointer to the ith _tree's node.
    parent = i;
    child = (_Xn[axis] < _posAcc.getVal(*pi, axis)) ? 1: 2;
    axis = (axis+1) % _dim;
    i = *(pi + child); //_tree(child, i)
  }

  *(tbegin + (3 * parent) + child) = _tree_sz; //_tree(child, parent) = _tree_sz
  *(tbegin + (3 * _tree_sz)) = n; //_tree(0, _tree_sz) = n;
  ++_tree_sz;
}

template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::__getClosest_through_path(const E_Float* pt, size_type& m, E_Float& d2) const
{
  size_type i(0), axis(0);
  const size_type * Ni, *tbegin(_tree.begin());
  E_Float d2tmp, D;

  m = IDX_NONE;
  d2 = NUGA::FLOAT_MAX;

  while ((i < _tree_sz) && (_tolerance < d2))
  {
    Ni = tbegin + 3*i; // pointer to the ith _tree's node.
    D = (pt[axis] - _posAcc.getVal(*Ni, axis));

    if (D*D < d2)
    {
      d2tmp = _posAcc.dist2(pt, *Ni);

      if (d2tmp < d2)
      {
        d2 = d2tmp;
        m = *Ni;
      }
    }
    i = (D < 0) ? *(Ni+1) : *(Ni+2);
    axis = (axis+1) % _dim;
  }
}

template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::__getClosest_through_path(size_type n, const E_Float *Xn, size_type& m, E_Float& d2) const
{
  size_type i(0), axis(0);
  const size_type *Ni, *tbegin(_tree.begin());

  E_Float d2tmp;

  m = IDX_NONE;
  d2 = NUGA::FLOAT_MAX;
  
  while ((i < _tree_sz) && (_tolerance < d2))
  {
    Ni = tbegin + 3*i; // pointer to the ith _tree's node.
    d2tmp = _posAcc.dist2(n, *Ni);

    if ((d2tmp < d2) && (n != *Ni))
    {
      d2 = d2tmp;
      m = *Ni;
    }
    i= (Xn[axis] < _posAcc.getVal(*Ni, axis)) ? *(Ni+1) : *(Ni+2);
    axis = (axis+1) % _dim;
  }
}

template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::__seek_closest
(size_type n, const E_Float *Xn, size_type ci, size_type axis, E_Float& d2, size_type& m) const
{
  if (ci == IDX_NONE) return;
  const size_type* Ni = _tree.begin() + 3 * ci;// pointer to the ith _tree's node.
  
  E_Float D = (Xn[axis] - _posAcc.getVal(*Ni, axis));

  if ((D*D) >= d2)
  {
    if (D < 0)
      __seek_closest(n, Xn, *(Ni+1), (axis+1)%_dim, d2, m); //left traversal
    else
      __seek_closest(n, Xn, *(Ni+2), (axis+1)%_dim, d2, m); //right traversal
  }
  else
  {
    E_Float di2 = _posAcc.dist2(n, *Ni);
    if ((di2 < d2) && ( n != *Ni))
    {
		  d2 = di2;
		  m = *Ni;
	  }
    __seek_closest(n, Xn, *(Ni+1), (axis+1)%_dim, d2, m); //left traversal
    __seek_closest(n, Xn, *(Ni+2), (axis+1)%_dim, d2, m); //right traversal
  }
}

template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::__seek_closest
(const E_Float *pt, size_type ci, size_type axis, E_Float& d2, size_type& m) const
{
	if (ci == IDX_NONE) return;

	const size_type* Ni = _tree.begin() + 3 * ci;// pointer to the ith _tree's node.
  E_Float D = (pt[axis] - _posAcc.getVal(*Ni, axis));
 
  if ((D*D) >= d2)
  {
    if (D < 0)
      __seek_closest(pt, *(Ni+1), (axis+1)%_dim, d2, m); //left traversal
    else
      __seek_closest(pt, *(Ni+2), (axis+1)%_dim, d2, m); //right traversal
  }
  else
  {
    E_Float di2 = _posAcc.dist2(pt, *Ni);
    if (di2 < d2)
    {
		  d2 = di2;
		  m = *Ni;
	  }
    __seek_closest(pt, *(Ni+1), (axis+1)%_dim, d2, m); //left traversal
    __seek_closest(pt, *(Ni+2), (axis+1)%_dim, d2, m); //right traversal
  }
}

template <typename CoordArrayType>
void K_SEARCH::KdTree<CoordArrayType>::__getInBox
(size_type ci, size_type axis, const E_Float* mBox, const E_Float* MBox, std::vector<size_type>& out) const
{
  if (ci == IDX_NONE) return;

	const size_type* Ni = _tree.begin() + 3 * ci;// pointer to the ith _tree's node.

  E_Float xi = _posAcc.getVal(*Ni, axis);

  bool do_left = (mBox[axis] <= xi);
  bool do_right = (xi <= MBox[axis]);

  if (do_left && do_right) // might be in so do the complete checking
  {
    bool is_in=true;
    size_type ax = (axis+1)%_dim;
    for (E_Int j = 0; j < _dim-1; ++j)
    {
      is_in &= ((mBox[ax] <= _posAcc.getVal(*Ni, ax)) && (_posAcc.getVal(*Ni, ax) <= MBox[ax]));
      ax = (ax+1)%_dim;
    }
    if (is_in)
      out.push_back(*Ni);
  }

  if (do_left) // left search
    __getInBox(*(Ni+1), (axis+1)%_dim, mBox, MBox, out);
  if (do_right) //right search
    __getInBox(*(Ni+2), (axis+1)%_dim, mBox, MBox, out);
}
#endif
