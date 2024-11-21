/*    
    Copyright 2013-2024 Onera.

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

#ifndef __GENERATOR_DELAUNAY_KERNEL_CXX__
#define __GENERATOR_DELAUNAY_KERNEL_CXX__

#include "Kernel.h"
#ifdef E_DEBUG2
#include "IO/io.h"
#endif

namespace DELAUNAY
{

/// Constructor
template <typename T>
Kernel<T>::Kernel(MeshData& data, const NUGA::MeshTool& tool)
    : _data(&data), _tool(&tool), _Ball_pred(*data.pos, _data->connectM), _constrained_pred(nullptr), _Nmatch(IDX_NONE)
{
  data.mask.resize(_data->connectM.cols(), true);

#ifdef E_TIME
  inval_time = remesh_time = cavity_time = 0.;
  sorting_bound_time = fix_cavity_time = init_cavity_time = 0.;
  _append_time = _base_time = 0.;
#endif
}

template <typename T>
void Kernel<T>::set(MeshData& data, const NUGA::MeshTool& tool)
{
  clear();

  _data = &data;
  _tool = &tool;
  if (_constrained_pred) delete _constrained_pred;
  _constrained_pred = nullptr;
  _Nmatch = IDX_NONE;
  // , _Ball_pred(*data.pos, _data->connectM), _constrained_pred(nullptr)

  _data->mask.resize(_data->connectM.cols(), true);

}

template <typename T>
void Kernel<T>::clear()
{
  _cavity.clear();
  _cboundary.clear();
  _base.clear();
  _sbound.clear();
  _real_cboundary.clear();
  _visited.clear();
  inodes.clear();
  Ancs.clear();
  elements.clear();
  _node_to_rightS.clear();
}

/// Destructor
template <typename T>
Kernel<T>::~Kernel(void)
{
  if (_constrained_pred) delete _constrained_pred;
}

///
template <typename T>
template <typename ConstraintType>
E_Int
Kernel<T>::insertNode(size_type N, const T& m, const ConstraintType& dummy)
{
#ifdef E_TIME
  chrono c;
  c.start();
#endif

  //
  E_Int ret = __getCavity<ConstraintType> (N, m, *_data->pos, _data->neighbors, _data->ancestors, _cavity, _cboundary);

  if (ret == -1 || ret == 2)
    return ret;

#ifdef E_TIME
  cavity_time += c.elapsed();
  c.start();
#endif

  //
  ret = __remeshCavity (N, _data->connectM, _data->neighbors, _data->ancestors, _cavity, _cboundary);
  if (ret)
    return ret;

#ifdef E_TIME
  remesh_time += c.elapsed();
  c.start();
#endif

  //
  __invalidCavityElements (_cavity, _data->connectM, _data->mask);

#ifdef E_TIME
  inval_time += c.elapsed();
#endif

  return ret;
}

template <typename T>
void
Kernel<T>::setConstraint(const NUGA::non_oriented_edge_set_type& hard_edges)
{
  if (_constrained_pred) delete _constrained_pred;
  _constrained_pred = new constrained_predicate(hard_edges);
}

///
template <typename T>
template <typename ConstraintType>
E_Int
Kernel<T>::__getCavity
(size_type N, const T& m,
 const K_FLD::FloatArray& pos, const K_FLD::IntArray& neighbors,
 const int_vector_type& ancestors, int_set_type& cavity,
 int_pair_vector_type& cboundary)
{

#ifdef E_TIME
  chrono c;
  c.start();
#endif

  _Ball_pred.setPoint(pos.col(N), m); // Set the cavity "center" and associated metric.

  E_Int ret = __getInitialCavity<ConstraintType>(N, pos, _data->connectM, neighbors,
                                 ancestors, _base, cavity, _sbound);

#ifdef E_TIME
  init_cavity_time += c.elapsed();
  c.start();
#endif

#ifdef E_DEBUG2
  if (ret == 2)
  {
    K_FLD::IntArray cc;
    for (int_set_type::const_iterator it = cavity.begin(); it != cavity.end(); ++it)
    {
      K_FLD::IntArray::const_iterator pK = _data->connectM.col(*it);
      cc.pushBack(pK, pK+3);
    }
    medith::write ("init_cav.mesh", *_data->pos, cc, "TRI");
  }
#endif
  
  if (ret == 2)
    return ret;

  if (ret == -1) // Error
    return ret;;

  if ((_base.size() > 2) || (cavity.size() != _base.size()))
    ret = __fixCavity(N, pos, _data->connectM, neighbors, ancestors, _base, cavity, _sbound);


  if (ret == -1) // Error
    return ret;

#ifdef E_DEBUG2 
  if (ret == 2)
  {
  K_FLD::IntArray cc;
  for (int_set_type::const_iterator it = cavity.begin(); it != cavity.end(); ++it)
  {
  K_FLD::IntArray::const_iterator pK = _data->connectM.col(*it);
  cc.pushBack(pK, pK+3);
  }

  medith::write ("fixed_cav.mesh", *_data->pos, cc, "TRI");  
  }
#endif

#ifdef E_TIME
  fix_cavity_time += c.elapsed();
  c.start();
#endif

  cboundary.clear();
  //_visited.clear();
  //__getSortedBoundary(_data->connectM, neighbors, *_base.begin(), 0, cavity, _sbound, _visited, cboundary);
  ret = __getSortedBoundary2(_data->connectM, _sbound, cboundary);
  _sbound.clear();

#ifdef E_TIME
  sorting_bound_time += c.elapsed();
#endif

  return ret;
}

///
template <typename T>
template <typename ConstraintType>
E_Int
Kernel<T>::__getInitialCavity
(size_type N, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  int_vector_type ba;
  size_type Nmatch;

  cavity.clear();
  cboundary.clear();
  base.clear();

#ifdef E_TIME
  chrono c;
  c.start();
#endif

  E_Int ret = _tool->getContainingElements (pos.col(N), pos, connect, neighbors, ancestors, _data->mask,
                                           _tool->getTree(), std::back_inserter(ba), Nmatch);

  size_type sz = (size_type)ba.size();
  for (size_type i = 0; i < sz; ++i)//fixme
    base.insert(ba[i]);

#ifdef E_TIME
  _base_time += c.elapsed();
  c.start();
#endif

  if (ret == 1) // On an edge.
  {
    //fixme
  }
  else if (ret == 2) // On an existing node.
  {
    this->_Nmatch = Nmatch;
    return ret;
  }
  else if (ret == -1) // Error
    return ret;

#ifdef E_DEBUG2 
  
  if (ret == 2)
  {  
  K_FLD::IntArray cc;
  for (int_set_type::iterator i = base.begin(); i != base.end(); ++i)
  {
  const E_Int* p = connect.col(*i);
  cc.pushBack(p,p+3);
  }
  medith::write("base.mesh", pos, cc);
  }
#endif

  __appendCavity<ConstraintType> (neighbors, base, cavity, cboundary);

#ifdef E_TIME
  _append_time += c.elapsed();
#endif

#ifdef E_DEBUG2
  if (ret == 2)
  {
  K_FLD::IntArray cc;
  for (int_set_type::iterator i = cavity.begin(); i != cavity.end(); ++i)
  {
  const E_Int* p = connect.col(*i);
  cc.pushBack(p,p+3);
  }
  medith::write("cavity.mesh", pos, cc);
  }
#endif

  return ret;

}

/// Default implementation : constrained iso/aniso
template <typename T>
template <typename ConstraintType>
void
Kernel<T>::__appendCavity
(const K_FLD::IntArray& neighbors, const int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  assert(_constrained_pred);
  for (int_set_type::const_iterator it = base.begin(); it != base.end(); ++it)
  {
    _tool->_inval.clear();
    _tool->getConnexSet1(*it, _data->connectM/*fixme*/, neighbors, cavity, 
                        cboundary, _Ball_pred, *_constrained_pred);
  }

  __getBoundary(cavity, neighbors, cboundary);
}

// Specialization : unconstrained iso.
template <> 
template <> inline
void
Kernel<E_Float>::__appendCavity<UNCONSTRAINED>
(const K_FLD::IntArray& neighbors, const int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  for (int_set_type::const_iterator it = base.begin(); it != base.end(); ++it)
  {
    _tool->_inval.clear();
    _tool->getConnexSet1(*it, _data->connectM/*fixme*/, neighbors, cavity, 
                        cboundary, _Ball_pred, _unconstrained_pred);
  }

  __getBoundary(cavity, neighbors, cboundary);
}

// Specialization unconstrained aniso.
template <> 
template <> inline
void
Kernel<Aniso2D>::__appendCavity<UNCONSTRAINED>
(const K_FLD::IntArray& neighbors, const int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  for (int_set_type::const_iterator it = base.begin(); it != base.end(); ++it)
  {
    _tool->_inval.clear();
    _tool->getConnexSet1(*it, _data->connectM/*fixme*/, neighbors, cavity, 
                        cboundary, _Ball_pred, _unconstrained_pred);
  }

  __getBoundary(cavity, neighbors, cboundary);
}

///
template <typename T>
E_Int
  Kernel<T>::__remeshCavity
(size_type N, K_FLD::IntArray & connect, K_FLD::IntArray& neighbors,
 int_vector_type& ancestors, const int_set_type& cavity,
 const int_pair_vector_type& cboundary)
{
  // Fast return.
  if (cavity.empty())       return -1; // Error
  if (cboundary.empty())    return -1; // Error

  size_type K0, Kprev(IDX_NONE), K1(connect.cols()), Kstart(connect.cols()), Kadj,j,jadj;
  K_FLD::IntArray::const_iterator pK0;
  size_type                   newN[] = {IDX_NONE, IDX_NONE, IDX_NONE};

  int_pair_type Bi;
  size_type triangle[3], N0, N1;

  triangle[2] = N;

  bool set_color = !_data->colors.empty();

  size_type color = -1;
  if (set_color)
    color = _data->colors[*cavity.begin()];

  size_type bsize = (size_type)cboundary.size();
  for (size_type i = 0; i < bsize; ++i)
  {
    Bi = cboundary[i];
    K0 = Bi.first;
    j = Bi.second;
    pK0 = connect.col(K0);

    N0 = triangle[0] = *(pK0 + (j+1) % element_type::NB_NODES);
    N1 = triangle[1] = *(pK0 + (j+2) % element_type::NB_NODES);

    connect.pushBack (triangle, triangle + 3);
    if (color != -1)
      _data->colors.push_back(color);

    Kadj = neighbors(j, K0);
    neighbors.pushBack (newN, newN + 3);
    neighbors(2, K1) = Kadj; // External adj

    if (Kadj != IDX_NONE)
    {
      jadj = element_type::getOppLocalNodeId(K0, j, connect, neighbors);
      assert (neighbors (jadj, Kadj) == K0);
      neighbors(jadj, Kadj) = K1; // Reciprocal
    }

    ancestors[N0] = ancestors[N1] = K1; // ancestors update

    if (Kprev != IDX_NONE) // Internal adj
    {
      neighbors(0, Kprev) = K1;
      neighbors(1, K1) = Kprev;
    }
    Kprev = K1++;
  }

  neighbors(0, Kprev) = Kstart;
  neighbors(1, Kstart) = Kprev;
  ancestors[N] = Kstart;

  return 0;
}

///
template <typename T>
void
  Kernel<T>::__invalidCavityElements
(const int_set_type& cavity, const K_FLD::IntArray& connect,
 bool_vector_type& mask)
{
  mask.resize(connect.cols(), true);
  for (int_set_type::const_iterator i = cavity.begin(); i != cavity.end(); ++i)
    mask[*i] = false;
}


///
template <typename T>
E_Int
  Kernel<T>::__fixCavity
(size_type N, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors,
 const int_vector_type& ancestors, const int_set_type& base, int_set_type& cavity,
 int_pair_set_type& cboundary)
{
  E_Int err = __ensureEmptyCavity(connectM, neighbors, ancestors, base, cavity, cboundary);
  if (err != -1)
    err = __ensureStarShape(N, pos, connectM, neighbors, base, cavity, cboundary);
  return err;
}

///
template <typename T>
void
  Kernel<T>::__getSortedBoundary
(const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors, size_type Ki, size_type b0,
 const int_set_type& cavity, int_pair_set_type& cboundary, int_set_type& visitedK, int_pair_vector_type& sorted_boundary)
{
  size_type                       Nj, Kj, b, b1;
  int_pair_type                   Bi;
  K_FLD::IntArray::const_iterator  pKi, pKj;

  if (cboundary.empty())
    return;

  visitedK.insert(Ki);

  for (size_type i = 0; i < element_type::NB_NODES; ++i)
  {
    b         = (b0+i) % element_type::NB_NODES;
    Bi.first  = Ki;
    Bi.second = b;

    Kj  = neighbors (b, Ki);

#ifdef E_DEBUG
    if (Ki != IDX_NONE)
      pKi = connectM.col(Ki);
    if (Kj != IDX_NONE)
      pKj = connectM.col(Kj);
#endif

    if (IS_IN(cboundary, Bi) || (Kj == IDX_NONE))//fixme : est-ce que la seconde condition est necessaire ?
    {
      sorted_boundary.push_back(Bi);
      cboundary.erase(Bi);
      continue;
    }

    if (IS_IN(visitedK,Kj)) // already visited;
      continue;

    pKi = connectM.col(Ki);
    pKj = connectM.col(Kj);
    Nj  = *(pKi + (b+2)% element_type::NB_NODES);
    b1  = element_type::getLocalNodeId (pKj, Nj);

    __getSortedBoundary(connectM, neighbors, Kj, b1, cavity, cboundary, visitedK, sorted_boundary);
  }
}

///
template <typename T>
E_Int
  Kernel<T>::__getSortedBoundary2
(const K_FLD::IntArray& connect, int_pair_set_type& in, int_pair_vector_type& out)
{
  size_type sz(in.size());
  if (sz == 0)
    return 1;
  
  out.resize(sz);
  _node_to_rightS.clear();
  size_type S,b,Ni;
  for (int_pair_set_type::iterator it = in.begin(); it != in.end(); ++it)
  {
    S = it->first;
    b = it->second;
    Ni = connect((b+1)%3, S);
    _node_to_rightS[Ni] = it;
  }

  out[0] = *in.begin();

  for (size_type i = 0; i < sz-1; ++i)
  {
    S = out[i].first;
    b = out[i].second;
    Ni = connect((b+2)%3, S);
    out[i+1] = *_node_to_rightS[Ni];
  }
  return 0;
}

/* vtune reports the uncommented is better but not very obvious looking at CPU...
///
template <typename T>
E_Int
  Kernel<T>::__ensureEmptyCavity
(const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, const int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  K_FLD::IntArray::const_iterator pK;
  int_pair_set_type::const_iterator it(cboundary.begin()), itEnd(cboundary.end());
  int_set_type::const_iterator itC(cavity.begin()), itCend(cavity.end());

  inodes.clear();

  // Store the cavity nodes.
  for (; itC != itCend; ++itC)
  {
    pK = connectM.col(*itC);
    inodes.insert (*pK);
    inodes.insert (*(pK + 1));
    inodes.insert (*(pK + 2));
  }

  // Store the boundary nodes.
  for (; it != itEnd; ++it)
  {
    const int_pair_type & Bi = *it;
    const size_type& K = Bi.first;
    const size_type& b = Bi.second;
    pK = connectM.col(K);

    inodes.erase (*(pK + (b+1) % element_type::NB_NODES));
    //inodes.erase (*(pK + (b+2) % element_type::NB_NODES));
  }

  if (inodes.empty())
    return 0;
  
  int_set_type removable = cavity;
  for (int_set_type::const_iterator itB = base.begin(); itB != base.end(); ++itB)
    removable.erase(*itB);
  
  while (!inodes.empty())
  {
    Ancs.clear();
    _tool->getAncestors (*inodes.begin(), ancestors, neighbors, std::back_inserter(Ancs));

    elements.clear();
    std::sort(ALL(Ancs));
    std::set_intersection (ALL(Ancs), ALL(removable), std::back_inserter(elements));

    if (elements.empty())
      return -1;
    
    const size_type& K = *elements.begin();
    pK = connectM.col(K);

    cavity.erase (K);
    removable.erase(K);

    for (size_type i = 0; i < element_type::NB_NODES; ++i)
    {
      inodes.erase(*(pK+i));
      int_pair_type B(K,i);
      if (!cboundary.insert(B).second) // if already in
        cboundary.erase(B);
    }
  }
  return 0;
}*/

///
template <typename T>
E_Int
  Kernel<T>::__ensureEmptyCavity
(const K_FLD::IntArray& connectM, const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, const int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  K_FLD::IntArray::const_iterator pK;
  int_pair_set_type::const_iterator it(cboundary.begin()), itEnd(cboundary.end());
  int_set_type::const_iterator itC(cavity.begin()), itCend(cavity.end());

  inodes.clear();
  
  // Store the inner nodes
  for (; itC != itCend; ++itC)
  {
    pK = connectM.col(*itC);
    inodes.insert(*pK);
    inodes.insert(*(pK + 1));
    inodes.insert(*(pK + 2));
  }
  
  // Store the boundary nodes.
  for (; it != itEnd; ++it)
  {
    const int_pair_type & Bi = *it;
    const size_type& K = Bi.first;
    const size_type& b = Bi.second;
    pK = connectM.col(K);

    inodes.erase (*(pK + (b+1) % element_type::NB_NODES));
    //inodes.insert (*(pK + (b+2) % element_type::NB_NODES));
  }
        
  if (inodes.empty())
    return 0;
  
  int_set_type removable = cavity;
  for (int_set_type::const_iterator itB = base.begin(); itB != base.end(); ++itB)
    removable.erase(*itB);
  
  while (!inodes.empty())
  {
    Ancs.clear();
    _tool->getAncestors (*inodes.begin(), connectM, ancestors, neighbors, std::back_inserter(Ancs));

    elements.clear();
    std::sort(ALL(Ancs));
    std::set_intersection (ALL(Ancs), ALL(removable), std::back_inserter(elements));

    if (elements.empty())
      return -1;
    
    const size_type& K = *elements.begin();
    pK = connectM.col(K);

    cavity.erase (K);
    removable.erase(K);

    for (size_type i = 0; i < element_type::NB_NODES; ++i)
    {
      inodes.erase(*(pK+i));
      int_pair_type B(K,i);
      if (!cboundary.insert(B).second) // if already in
        cboundary.erase(B);
    }
  }
  return 0;
}

///
template <typename T>
E_Int
  Kernel<T>::__ensureStarShape
(size_type N, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectM, 
 const K_FLD::IntArray& neighbors, const int_set_type& base,
 int_set_type& cavity, int_pair_set_type& cboundary)
{
  bool carry_on;
  size_type S,b, Ni, Nj;
  K_FLD::IntArray::const_iterator pS;
  E_Float tolerance(_tool->getTolerance()), q;

  do 
  {
    carry_on = false;
    for (int_pair_set_type::iterator it = cboundary.begin(); it != cboundary.end(); ++it)
    {
      const int_pair_type& B = *it;
      S = B.first;
      b = B.second;

      if (base.find(S) != base.end())
        continue;

      pS = connectM.col(S);

      Ni = *(pS + (b+1)% element_type::NB_NODES);
      Nj = *(pS + (b+2)% element_type::NB_NODES);

      q = K_MESH::Triangle::qualityG<2>(pos.col(N), pos.col(Ni), pos.col(Nj));

      if (q > tolerance)
        continue;

      cavity.erase(S);
      carry_on = true;
    }

    if (carry_on)
      __getBoundary(cavity, neighbors, cboundary);
  }
  while(carry_on);

  return 0;
}
  
///
template <typename T>
void
  Kernel<T>::__getBoundary
(int_set_type& cavity, const K_FLD::IntArray& neighbors, int_pair_set_type& cboundary)
{
  cboundary.clear();
  for (int_set_type::const_iterator it = cavity.begin(); it != cavity.end(); ++it)
  {
    size_type S = *it;
    for (size_type i = 0; i < 3; ++i)
    {
      size_type Sn = neighbors(i, S);
      if (cavity.find(Sn) == cavity.end())
        cboundary.insert(int_pair_type(S,i));
    }
  }
}



}

#endif
