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

#ifndef __K_CONNECT_IDTOOL_H__
#define __K_CONNECT_IDTOOL_H__

#include "Nuga/include/defs.h"
#include <vector>
#include<map>
#define Vector_t std::vector

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_unit.h"


namespace K_CONNECT
{

class IdTool {
  
public:
  /// bijective
  static void reverse_indirection
  (const Vector_t<E_Int> & assoc, Vector_t<E_Int>& reverse_assoc);
  /// non-bijective : convert a n-to-one vector (tipically an oids) to a map
  static void reverse_indirection(const E_Int*oids, E_Int sz, std::map<E_Int, std::vector<E_Int>>& split_map);
  /// non-bijective : convert a n-to-one vector (tipically an oids) to a ngon_unit
  static void reverse_indirection(E_Int nb_pgs, const E_Int*oids, E_Int sz, ngon_unit& split_graph);
  /// 
  template < short LIST_sZ >
  static void right_shift(E_Int* list, int shift);
  // variable stride version
  static void right_shift(E_Int* list, E_Int sz, int shift);
  ///
  static void reverse_sorting(Vector_t<E_Int> & vec);
  ///
  static void negative (Vector_t<bool>& flag);
  static void negative(Vector_t<E_Int>& flag);
  ///
  static void propagate(const Vector_t<E_Int>& nids, Vector_t<E_Int>& oids);
  /// Convert a flag array to a corresponding compacted indirection old_to_new
  static void build_indir(const std::vector<bool>& keep, std::vector<E_Int> & nids);
  /// Convert a flag array to both corresponding compacted indirections old_to_new and new_to_old
  static void build_indir(const std::vector<bool>& keep, std::vector<E_Int> & nids, std::vector<E_Int> & oids);
  ///
  template < typename T>
  static E_Int compress_unic(std::vector<T>& vec);
  ///
  template < typename T, typename Predicate_t>
  static E_Int compress(std::vector<T>& vec, const Predicate_t& P);
  ///
  template < typename T, typename Predicate_t>
  static E_Int compress(std::vector<T>& vec, const Predicate_t& P, std::vector<E_Int>& nids);
  ///
  template < typename T, typename Predicate_t>
  static E_Int compress(K_FLD::DynArray<T>& arr, const Predicate_t& P);
  ///
  template < typename T, typename Predicate_t>
  static E_Int compress(K_FLD::DynArray<T>& arr, const Predicate_t& P, std::vector<E_Int>& nids);
  ///
  template < typename T>
  static E_Int compress(K_FLD::DynArray<T>& arr, const std::vector<E_Int>& keepids, E_Int idx_start);
  ///
  template < typename T>
  static K_FLD::DynArray<T> compress_(K_FLD::DynArray<T> const & arr, const std::vector<E_Int>& keepids, E_Int idx_start);
 
  ///
  template < typename T, typename Predicate_t>
  static E_Int compress(K_FLD::FldArray<T>& arr, const Predicate_t& P);
  ///
  static E_Int max(const K_FLD::IntArray& connect);
  static E_Int max(const Vector_t<E_Int>& vec);
  ///
  static E_Int min(const K_FLD::IntArray& connect);
  static E_Int min(const Vector_t<E_Int>& vec);
  ///
  static E_Int sum(const Vector_t<E_Int>& vec){E_Int s = 0; for (size_t i = 0; i < vec.size(); ++i) { s += vec[i]; } return s;}
  ///
  static void compact(std::vector<E_Int>& vec, const std::vector<bool> & flag);
  ///
  template <typename T>
  static void compact(std::vector<T>& vec, const std::vector<E_Int> & nids);
  ///
  template <typename T>
  static std::vector<T> compact_(std::vector<T> const & vec, std::vector<E_Int> const & nids);
  ///
  template <typename T>
  static void compact(K_FLD::DynArray<T>& arr, const std::vector<E_Int> & nids);
  ///
  template <typename T>
  static E_Int compact (E_Int* target, E_Int& szio, const T* keep/*, Vector2& new_Ids*/);
  
  ///
  template <typename T>
  static bool equal_vec(const Vector_t<T>& a, const Vector_t<T>& b);
  
  ///
  static bool equal(const E_Int* p, const E_Int* q, E_Int n, bool permut_accepted, bool strict_orient);

  static bool check_in_range(const E_Int*p, E_Int sz, E_Int minval, E_Int maxval)
  {
    for (E_Int k = 0; k < sz; ++k)
    {
      if (p[k] < minval) return false;
      if (p[k] > maxval) return false;
    }
    return true;
  }
  
  ///
  template<typename VecDec>
  static void shift (VecDec & vec, E_Int shift){for (size_t i=0;i<vec.size(); ++i)if (vec[i] != IDX_NONE) vec[i] += shift;}
  template<typename VecDec>
  static void shift (VecDec & vec, E_Int from, E_Int shift){for (size_t i=from;i<vec.size(); ++i) if (vec[i] != IDX_NONE) vec[i]+=shift;}

  template<typename VecDec>
  static void init_inc(VecDec & vec, E_Int sz){ vec.clear(); vec.resize(sz, IDX_NONE); for (E_Int i = 0; i < sz; ++i) vec[i] = i; }
  
  template<typename VecDec>
  static void init_inc(VecDec & vec, E_Int sz, E_Int shift){ vec.clear(); vec.resize(sz, IDX_NONE); for (E_Int i = 0; i < sz; ++i) vec[i] = i+shift; }

  template <typename T>
  static void init_inc(K_FLD::DynArray<T>& arr, E_Int row, E_Int sz){ 
    arr.resize(arr.rows(), sz, IDX_NONE);  
    for (E_Int i = 0; i < sz; ++i) arr(row, i) = i; }

  template <typename IntCont, typename Edge_t>
  static void loc2glob(std::map<Edge_t, IntCont>& m, E_Int idx_start/*in m*/, const std::vector<E_Int>& oids);


 static inline E_Int get_pos(const E_Int* vec, E_Int n, E_Int val)
  {
    for (E_Int i = 0; i < n; i++)
      if (vec[i] == val) return i;
    assert(false);
    return -1;
  }

};

/// 
template <short LIST_SZ>
void IdTool::right_shift(E_Int* list, int shift)
{
  if (shift == 0) return;

  E_Int tmp[LIST_SZ];
  for (int i =0; i < LIST_SZ; ++i)
    tmp[i] = list[(i+shift)% LIST_SZ];

  for (int i =0; i < LIST_SZ; ++i)
    list[i] = tmp[i];    
}

///
template<>
inline void IdTool::shift<K_FLD::IntArray> (K_FLD::IntArray & arr, E_Int from, E_Int shift)
{
  for (E_Int i=from;i<arr.cols(); ++i) 
  {
    for (E_Int j = 0; j < arr.rows(); ++j)
      if (arr(j,i) != IDX_NONE) arr(j,i) +=shift;
  }
}

///
template <typename T>
bool IdTool::equal_vec(const Vector_t<T>& a, const Vector_t<T>& b)
{
  if (a.size() != b.size())
    return false;
  
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[i])
      return false;
  
  return true;
}

///
template < typename T>
E_Int IdTool::compress_unic(std::vector<T>& vec)
{
  std::set<T> uset;
  std::vector<T> new_vec;

  new_vec.reserve(vec.size());

  for (size_t i=0; i < vec.size(); ++i)
    if (uset.insert(vec[i]).second) // not already in
      new_vec.push_back(vec[i]);

  vec = new_vec;
  return 0;
}

///
struct unchanged // : public std::unary_function <E_Int, bool>
{
  unchanged(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i]==i);
  }

  const std::vector<E_Int>& _indir;
};

///
struct valid // : public std::unary_function <E_Int, bool>
{
  valid(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i] != IDX_NONE);
  }

  const std::vector<E_Int>& _indir;
};

///
struct invalid // : public std::unary_function <E_Int, bool>
{
  invalid(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i] == IDX_NONE);
  }

  const std::vector<E_Int>& _indir;
};

///
template< typename T = bool>
struct keep // : public std::unary_function <E_Int, bool>
{
  keep(const std::vector<T>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i]);
  }

  const std::vector<T>& _indir;
};
template< typename T = bool>
struct keep2 // : public std::unary_function <E_Int, bool>
{
  keep2(const T* indir, E_Int n):_indir(indir), _n(n){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i]);
  }

  const T* _indir;
  E_Int _n;
};

struct strictly_positive // : public std::unary_function <E_Int, bool>
{
  strictly_positive(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i] > 0);
  }

  const std::vector<E_Int>& _indir;
};

///
template < typename T, typename Predicate_t>
E_Int IdTool::compress(std::vector<T>& vec, const Predicate_t& P)
{
  //assert (vec.size() == P._indir.size());
  size_t         i, cols(vec.size());
  std::vector<T> new_vec;
  for (i = 0; i < cols; ++i)
  {
    if (P(i)) new_vec.push_back(vec[i]);
  }

  E_Int ret = vec.size() - new_vec.size();
  vec = new_vec;
  return ret;
}

///
template < typename T, typename Predicate_t>
E_Int IdTool::compress(std::vector<T>& vec, const Predicate_t& P, std::vector<E_Int>& nids)
{
  assert(vec.size() == P._indir.size());
  size_t         i, cols(vec.size()), count(0);
  std::vector<T> new_vec;

  nids.clear();
  nids.resize(vec.size(), IDX_NONE);

  for (i = 0; i < cols; ++i)
  {
    if (P(i))
    {
      new_vec.push_back(vec[i]);
      nids[i] = count++;
    }
  }

  E_Int ret = vec.size() - new_vec.size();
  vec = new_vec;
  return ret;
}

///
template < typename T, typename Predicate_t>
E_Int IdTool::compress(K_FLD::DynArray<T>& arr, const Predicate_t& P)
{
  if (arr.cols() == 0)
    return 0;
  //assert (arr.cols() == P._indir.size());
  size_t         i, cols(arr.cols()), stride(arr.rows());
  K_FLD::DynArray<T> new_arr;
  for (i = 0; i < cols; ++i)
  {
    if (P(i)) new_arr.pushBack(arr.col(i), arr.col(i)+stride);
  }

  E_Int ret = arr.cols() - new_arr.cols();
  arr = new_arr;
  return ret;
}

///
template < typename T, typename Predicate_t>
E_Int IdTool::compress(K_FLD::DynArray<T>& arr, const Predicate_t& P, std::vector<E_Int>& nids)
{
  if (arr.cols() == 0)
    return 0;
  //assert (arr.cols() == P._indir.size());
  size_t         i, cols(arr.cols()), stride(arr.rows()), count(0);
  
  nids.clear();
  nids.resize(cols, IDX_NONE);
  
  K_FLD::DynArray<T> new_arr;
  for (i = 0; i < cols; ++i)
  {
    if (P(i))
    {
      new_arr.pushBack(arr.col(i), arr.col(i)+stride);
      nids[i] = count++;
    }
  }

  E_Int ret = arr.cols() - new_arr.cols();
  arr = new_arr;
  return ret;
}

template < typename T>
E_Int IdTool::compress(K_FLD::DynArray<T>& arr, const std::vector<E_Int>& keepids, E_Int idx_start)
{
  E_Int ret = arr.cols();
  arr = std::move(compress_(arr, keepids, idx_start));
  ret -= arr.cols();
 
  return ret;
}

template < typename T>
K_FLD::DynArray<T> IdTool::compress_(K_FLD::DynArray<T> const & arr, const std::vector<E_Int>& keepids, E_Int idx_start)
{
  K_FLD::DynArray<T> new_arr;
  E_Int stride(arr.rows());
  
  new_arr.reserve(stride, keepids.size());
  for (size_t i = 0; i < keepids.size(); ++i)
  {
    E_Int id = keepids[i] - idx_start;
    new_arr.pushBack(arr.col(id), arr.col(id) + stride);
  }
  return new_arr;
}

/// compact (with reordering) of a vector using an indirection
template <typename T>
void IdTool::compact(std::vector<T>& vec, const std::vector<E_Int> & nids)
{
  if (vec.empty() || nids.empty()) return;
  vec = compact_(vec, nids);
}

/// compact (with reordering) of a vector using an indirection
template <typename T>
std::vector<T> IdTool::compact_(std::vector<T> const & vec, std::vector<E_Int> const & nids)
{
  std::vector<T> tmp;

  E_Int s = max(nids) + 1;
  tmp.resize(s);

  for (size_t i = 0; i < vec.size(); ++i)
  {
    const E_Int& ni = nids[i];
    if (ni != IDX_NONE)
      tmp[ni] = vec[i];
  }
  return tmp;
}

/// compact (with reordering) of a vector using a indirection from compacting an IntArray
template <typename T>
void IdTool::compact(K_FLD::DynArray<T>& arr, const std::vector<E_Int> & nids)
{
  size_t COLS(arr.cols());
  if (COLS==0) return;
  
  K_FLD::DynArray<T> tmp;
    
  size_t s = max(nids)+1, ROWS(arr.rows());
  tmp.resize(ROWS, s);
    
  for (size_t i=0; i< COLS; ++i )
  {
    const E_Int& ni=nids[i];
    if (ni != IDX_NONE)
    {
      for (size_t k=0; k < ROWS; ++k)
        tmp(k,ni)=arr(k,i);
    }
  }
  arr=tmp;
}

///
template <typename T>
E_Int IdTool::compact (E_Int* target, E_Int& szio, const T* keep/*, Vector2& new_Ids*/)
{

  bool  carry_on(false);
  E_Int i1(0), i2(szio-1);

  //new_Ids.clear();
  //new_Ids.resize(cols, IDX_NONE);

  do{

    while ((i1 <= i2) && keep[i1]){/*new_Ids[i1] = i1;*/ ++i1;}  // Get the first empty column.
    while ((i1 <= i2) && !keep[i2]){--i2;} // Get the first column to move from the tail.

    carry_on = (i1 < i2);

    if (carry_on)
    { // Copy column i2 in column i1.
      //new_Ids[i2] = i1;
      E_Int& v2 = target[i2--];
      E_Int& v1 = target[i1++];
      v1 = v2;
    }
  }
  while (carry_on);

  E_Int szi = szio;
  szio = i1;

  return (szi - szio);
}

///
template <typename IntCont, typename Edge_t>
void IdTool::loc2glob(std::map<Edge_t, IntCont>& m, E_Int idx_start /* in m*/, const std::vector<E_Int>& oids)
{
  std::map<Edge_t, IntCont> gm;
  Edge_t e;

  for (auto& it : m)
  {
    e.setNodes(oids[it.first.node(0) - idx_start], oids[it.first.node(1) - idx_start]);
    for (auto & n : it.second) n = oids[n - idx_start];
    gm[e] = it.second;
  }

  m = gm;
}

}

#endif
