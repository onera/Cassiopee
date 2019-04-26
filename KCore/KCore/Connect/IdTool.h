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

#ifndef __K_CONNECT_IDTOOL_H__
#define __K_CONNECT_IDTOOL_H__

#include "Def/DefTypes.h"
#include <vector>
#define Vector_t std::vector

#include "Fld/DynArray.h"
#include "Fld/ngon_unit.h"


namespace K_CONNECT
{

class IdTool {
  
public:
  /// bijective
  static void reverse_indirection
  (const Vector_t<E_Int> & assoc, Vector_t<E_Int>& reverse_assoc);
  /// non-bijective : convert a n-to-one vector (tipically an oids) to a ngon_unit
  static void reverse_indirection(E_Int nb_pgs, const E_Int*oids, E_Int sz, ngon_unit& split_graph);
  /// 
  template < E_Int S >
  static void right_shift(E_Int* list, E_Int sz);
  ///
  static void reverse_sorting(Vector_t<E_Int> & vec);
  ///
  static void negative (Vector_t<bool>& flag);
  ///
  static void propagate(const Vector_t<E_Int>& nids, Vector_t<E_Int>& oids);
  /// Convert a flag array to a corresponding compacted indirection old_to_new
  static void build_indir(const std::vector<bool>& keep, std::vector<E_Int> & nids);
  /// Convert a flag array to both corresponding compacted indirections old_to_new and new_to_old
  static void build_indir(const std::vector<bool>& keep, std::vector<E_Int> & nids, std::vector<E_Int> & oids);
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
  static E_Int compress(K_FLD::FldArray<T>& arr, const Predicate_t& P);
  ///
  static E_Int max(const K_FLD::IntArray& connect);
  static E_Int max(const Vector_t<E_Int>& vec);
  ///
  static E_Int min(const K_FLD::IntArray& connect);
  static E_Int min(const Vector_t<E_Int>& vec);
  ///
  static void compact(std::vector<E_Int>& vec, const std::vector<bool> & flag);
  template <typename T>
  static void compact(std::vector<T>& vec, const std::vector<E_Int> & nids);
  template <typename T>
  static void compact(K_FLD::DynArray<T>& arr, const std::vector<E_Int> & nids);
  
  ///
  template <typename T>
  static bool equal_vec(Vector_t<T>& a, Vector_t<T>& b);
  
  ///
  static bool equal(const E_Int* p, const E_Int* q, E_Int n, bool permut_accepted, bool strict_orient);
  
  ///
  template<typename VecDec>
  static void shift (VecDec & vec, E_Int shift){for (size_t i=0;i<vec.size(); ++i)if (vec[i] != E_IDX_NONE) vec[i] += shift;}
  template<typename VecDec>
  static void shift (VecDec & vec, E_Int from, E_Int shift){for (size_t i=from;i<vec.size(); ++i) if (vec[i] != E_IDX_NONE) vec[i]+=shift;}

  template<typename VecDec>
  static void init_inc(VecDec & vec, E_Int sz){ vec.clear(); vec.resize(sz, E_IDX_NONE); for (E_Int i = 0; i < sz; ++i) vec[i] = i; }

  template <typename T>
  static void init_inc(K_FLD::DynArray<T>& arr, E_Int row, E_Int sz){ 
    arr.resize(arr.rows(), sz, E_IDX_NONE);  
    for (E_Int i = 0; i < sz; ++i) arr(row, i) = i; }

};

/// 
template <E_Int S>
void IdTool::right_shift(E_Int* list, E_Int sz)
{
    if (sz == 0) return; 

      E_Int tmp[S];

      for (int i =0; i < S; ++i){
          tmp[i] = list[(i+sz)%S];
      }

      for (int i =0; i < S; ++i){
          list[i] = tmp[i];
      }    
}

///
template <typename T>
bool IdTool::equal_vec(Vector_t<T>& a, Vector_t<T>& b)
{
  if (a.size() != b.size())
    return false;
  
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[i])
      return false;
  
  return true;
}

///
struct unchanged : public std::unary_function <E_Int, bool>
{
  unchanged(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i]==i);
  }

  const std::vector<E_Int>& _indir;
};

///
struct valid : public std::unary_function <E_Int, bool>
{
  valid(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i] != E_IDX_NONE);
  }

  const std::vector<E_Int>& _indir;
};

///
struct invalid : public std::unary_function <E_Int, bool>
{
  invalid(const std::vector<E_Int>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i] == E_IDX_NONE);
  }

  const std::vector<E_Int>& _indir;
};

///
template< typename T = bool>
struct keep : public std::unary_function <E_Int, bool>
{
  keep(const std::vector<T>& indir):_indir(indir){}
  inline bool operator() (E_Int i ) const
  {
    return (_indir[i]);
  }

  const std::vector<T>& _indir;
};

struct strictly_positive : public std::unary_function <E_Int, bool>
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
  assert (vec.size() == P._indir.size());
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
  nids.resize(vec.size(), E_IDX_NONE);

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
  assert (arr.cols() == P._indir.size());
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

///fixme : mem corruption
/*template < typename T, typename Predicate_t>
E_Int IdTool::compress(K_FLD::FldArray<T>& arr, const Predicate_t& P)
{
  //assert (arr.getSize() == P._indir.size());
  size_t         i, cols(arr.getSize());
  K_FLD::FldArray<T> new_arr(1, arr.getNfld());
  new_arr.resize(0,new_arr.getNfld());
  for (i = 0; i < cols; ++i)
  {
    if (P(i)) new_arr.pushBack(arr, i);
  }
  E_Int ret = arr.getSize() - new_arr.getSize();
  arr = new_arr;
  return ret;
}*/

/// compact (with reordering) of a vector using a indirection from compacting an IntArray
template <typename T>
void IdTool::compact(std::vector<T>& vec, const std::vector<E_Int> & nids)
{
  
  E_Int sz = vec.size();
  if (sz == 0) return;
  
  std::vector<T> tmp;
    
  E_Int s = max(nids)+1;
  tmp.resize(s);
    
  for (E_Int i=0; i<sz; ++i )
  {
    const E_Int& ni=nids[i];
    if (ni != E_IDX_NONE)
      tmp[ni]=vec[i];
  }
  vec=tmp;
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
    if (ni != E_IDX_NONE)
    {
      for (size_t k=0; k < ROWS; ++k)
        tmp(k,ni)=arr(k,i);
    }
  }
  arr=tmp;
}

}

#endif
