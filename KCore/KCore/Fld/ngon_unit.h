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

#ifndef __NGON_UNIT_H__
#define	__NGON_UNIT_H__

#include <vector>
#include <set>
#include "Def/DefTypes.h"
#include <assert.h>

#include "Fld/DynArray.h"

#define Vector_t std::vector

template <typename T> struct ngon_t;

#define INNER 0
#define INITIAL_SKIN 1
#define CONNEXION_SKIN 3
#define UNUSED -1 // for PH in _ng1 or _ng2 not involved

#define IS_EXTERNAL(a) (a==INITIAL_SKIN || a==CONNEXION_SKIN)

class ngon_unit
{
  public:
    
    friend struct ngon_t<K_FLD::IntArray>;
    
    ///
    ngon_unit():_dirty(true){};
    ///
    explicit ngon_unit (const E_Int* begin);//Cassiopee format
    ///
    ngon_unit(const ngon_unit& ngin);
    /// 
    ngon_unit(const E_Int* begin, E_Int sz, E_Int nbe);
        
    /// ngon_unit& operator=(const Vector_t<E_Int>& vNGON);
    ngon_unit& operator=(const ngon_unit& ngin);
    ///
    E_Float operator[](const E_Int i){return _NGON[i];}
    /// Refresh the facets starts.
    void updateFacets() const;
    
    /// memory
    void reserve (E_Int sz){_NGON.reserve(sz); _facet.reserve(sz); _type.reserve(sz); _ancEs.reserve(2, sz);} //fixme : wrong for _NGON
    
    void expand_n_fixed_stride (E_Int n, E_Int strd){
      E_Int sz0 = _facet.size();
      _NGON.resize(_NGON.size() + n*(strd+1), strd) ; _facet.resize(sz0 + n) ; 
      if (!_type.empty()) _type.resize(_type.size() + n) ;
      if (_ancEs.cols()!=0) _ancEs.resize(2, _ancEs.cols() + n);
//      E_Int pos = _facet[sz0-1] + strd+1;
//      // set the facets pos
//      for (E_Int i=0; i< n; ++i, pos += strd+1)
//        _facet[sz0 + i] = pos;

      _NGON[0] += n;
      _NGON[1] = _NGON.size() - 2;
      // the following update has been added since the commented block above doesn't work in Basic element mode.
      _dirty=true;
      updateFacets(); 
    }
    E_Int capacity(){return _NGON.capacity();}
    // Interrogations
    /// Returns the number of facets for the i-th PH(PG)  (ArrayAccessor interface)
    inline E_Int stride(E_Int i) const {return _NGON[_facet[i]];}
    /// Returns the j-th element.  (ArrayAccessor interface)
    inline void getEntry(const E_Int& j, E_Int* pE) const
    {for (E_Int k=0;k<_NGON[_facet[j]]/*stride(j)*/;++k) pE[k]=_NGON[_facet[j]+1+k];/*get_facet(j,k)*/}
    /// Returns the j-th facet of the i-th element
    inline E_Int get_facet(E_Int i, E_Int j) const {return _NGON[_facet[i]+1+j];}
    inline E_Int& get_facet(E_Int i, E_Int j) {return _NGON[_facet[i]+1+j];}
    /// Returns the j-th facet of the i-th element
    inline const E_Int* get_facets_ptr(E_Int i) const {return &_NGON[_facet[i]+1];}
    inline E_Int* get_facets_ptr(E_Int i) {return &_NGON[_facet[i]+1];}
    /// Returns the number of entities
    inline E_Int size() const {return !_NGON.empty() ? _NGON[0] : 0;}
    ///
    inline const E_Int* begin() const {return &_NGON[0];}
    ///
    E_Int get_facets_max_id() const ;
    ///
    void get_indices_of_type (E_Int FLAG, Vector_t<E_Int>& indices) const;

    void get_stride_extrema(E_Int& mins, E_Int& maxs) const;
    ///
    void flag_indices_of_type (E_Int FLAG, Vector_t<bool>& flag) const;
    ///
    void find_elts_with_facets(const std::set<E_Int>& facets, std::vector<E_Int> & elts) const;
    ///
    void unique_indices(std::vector<E_Int>& indices) const;
    ///
    void extract (const Vector_t<E_Int>& indices, ngon_unit& ng_out, Vector_t<E_Int>& oldIds) const;
    void extract (const E_Int* ptr, E_Int n, ngon_unit& ng_out, Vector_t<E_Int>& oids) const;
    ///
    void extract_of_type (E_Int FLAG, ngon_unit& ng_out, Vector_t<E_Int>& oldIds) const;
    /// 
    template <typename Predicate_t>
    void extract_by_predicate (const Predicate_t& P, ngon_unit& ngo, Vector_t<E_Int>& oids, Vector_t<E_Int>& nids) const;
    
    // Transformations
    ///
    void append(const ngon_unit& ng);
    ///
    void append(const Vector_t<E_Int>& NGON);
    ///
    void append(const ngon_unit& cngon_in, const E_Int* first, E_Int nb_elts);
    ///
    void clear() { _NGON.clear(); _facet.clear(); _type.clear(); _ancEs.clear(); _dirty = true; }
    ///
    void change_indices (const Vector_t<E_Int>& nIds, bool zerobased);
    ///
    void get_degenerated(Vector_t<E_Int>& indices);
    ///
    void get_duplicated(Vector_t<E_Int>& indices);
    ///
    E_Int remove_duplicated();
    ///
    E_Int remove_consecutive_duplicated();
    ///
    E_Int remove_entities (const Vector_t<E_Int>& to_remove, Vector_t<E_Int>& nids);
    ///
    E_Int remove_facets(const Vector_t<E_Int>& facet_ids, Vector_t<E_Int>& nids, E_Int min_facets_nb=0); //0-based nids
    
    ///
    void shift(E_Int shift);
    
    /// 
    void sort_by_type(Vector_t<E_Int>& nids, Vector_t<E_Int>& oids) const;
    
    ///
    template <typename T>
    static void shift(Vector_t<T>& vec, const T& val){for (size_t i = 0; i < vec.size(); ++i)vec[i]+=val;}
    ///
    void reset_facets();//set all facet values to E_IDX_NONE (useful for building a neighbor table)
    
    /// warning : need a call to updateFacets afterwards
    template <template<typename, typename> class Container_t, typename Element,typename Allocator_t>
    void add (const Container_t<Element, Allocator_t>& molecule);

    /// warning : need a call to updateFacets afterwards
    void add(E_Int n, const E_Int* facet_ptr);

    //check
    bool is_fixed_stride(E_Int& stride) const;
    bool attributes_are_consistent() const ;
    
    /// Conversions

    /// convert a fixed-stride storage to a ngon_unit storage : work only for polygons and tetrahedra
    static void convert_fixed_stride_to_ngon_unit(const K_FLD::IntArray&cnt, E_Int shift, ngon_unit& nguo);
    /// reciprocal 
    static void convert_ngon_unit_to_fixed_stride(const ngon_unit& ngui, E_Int shift, K_FLD::IntArray& cnto);

    // tranfer attributes for a reduced set (tipically an agglomeration)
    void compact_attributes(const ngon_unit& ngi, const Vector_t<E_Int>& nids);
    void compact_attributes(const ngon_unit& ngi, const Vector_t<bool>& keep);
    // tranfer attributes for a larger set (tipically a split)
    void spread_attributes(const ngon_unit& ngi, const Vector_t<E_Int>& oids);
    
private:
    /// warning : need a call to updateFacets afterwards
    void __add (const ngon_unit& ng, E_Int ith);
    ///
    
     
public:
    Vector_t<E_Int> _NGON;
    mutable Vector_t<E_Int> _facet;// _facet[i] is the index in _data of the sub entity of i-th element.
    Vector_t<E_Int>  _type;
    K_FLD::IntArray  _ancEs;
    mutable bool _dirty; 
};


///
template <template<typename, typename> class Container_t, typename Element, typename Allocator_t>
void ngon_unit::add (const Container_t<Element, Allocator_t>& molecule)
{
  // molecule here is one PH or PG
  if (_NGON.empty())
  {
    _NGON.resize(2,0);
    _dirty=false;
  }
  
  _NGON[0]++;
  _NGON.insert(_NGON.end(), molecule.begin(), molecule.end());
  _NGON[1]=_NGON.size()-2;
  
  if (!_dirty )
    _facet.push_back(_NGON.size()-molecule.size());
}

///
template <typename Predicate_t>
void ngon_unit::extract_by_predicate (const Predicate_t& P, ngon_unit& ngo, Vector_t<E_Int>& oids, Vector_t<E_Int>& nids) const
{
  // 0-based indices
  
  size_t sz = size();
  assert (_type.size() == sz);
  
  oids.clear();
  
  updateFacets();
  
  nids.clear();
  nids.resize(sz, E_IDX_NONE);
  
  ngo.clear();
  
  // resize
  E_Int nb_ents = 0;
  E_Int cumul_stride = 0;
  for (size_t i = 0; i < sz; ++i){
    if (P(i)){
      ++nb_ents;
      cumul_stride += stride(i);
    }
  }
  
  if (nb_ents == 0) return;
  
  oids.resize(nb_ents);
  
  ngo._NGON.resize(2 + cumul_stride + nb_ents, 0);
  ngo._facet.resize(nb_ents);
  
  bool update_types = !_type.empty();
  if (update_types) ngo._type.resize(nb_ents);
  bool update_ancs = (_ancEs.cols() != 0);
  if (update_ancs) ngo._ancEs.resize(2, nb_ents);
  ngo._dirty = false;
  
  // fill it
  E_Int pos(2);
  nb_ents = 0;
  for (size_t i = 0; i < _type.size(); ++i)
  {
    if (!P(i)) continue;
    
    oids[nb_ents] = i;
    nids[i] = nb_ents;
    
    ngo._facet[nb_ents] = pos;

    E_Int s = stride(i);
    ngo._NGON[pos++] = s;

    const E_Int* ptr = get_facets_ptr(i);
    std::copy(ptr, ptr + s, &ngo._NGON[pos]);

    pos += s;

    if (update_types) ngo._type[nb_ents]=_type[i];
    if (update_ancs) ngo._ancEs.pushBack(_ancEs.col(i), _ancEs.col(i) + 2);

    ngo._NGON[0] = ++nb_ents; 
  }
  
  ngo._NGON[1] = ngo._NGON.size() - 2;
  
}

#ifdef DEBUG_BOOLEAN
#include <iostream>
inline std::ostream &operator<<(std::ostream& out, const ngon_unit& ng)
{
  out << "####################################" << std::endl;
/*
  // Print out the matrix.
  size_t nb_ents = ng._NGON[0];
  size_t sz = ng._NGON[1];
  for (E_Int i = 2; i < sz; ++i)
  {
    std::cout << ng._NGON[i] << ": ";
    
    for (size_t j=i+1; j < i+1+ng._NGON[i]; ++j)
      std::cout << ng._NGON[j] << " ";
    
    i += ng._NGON[i];
    
    std::cout << std::endl;
  }
*/
  
  out << "############ NGON ####################" << std::endl;
  for (size_t i = 0; i < ng._NGON.size(); ++i)
    out << ng._NGON[i] << " " << std::endl;
  out << std::endl;
  out << "####################################" << std::endl;
  
  out << "############ FACETS ###################" << std::endl;
  for (size_t i = 0; i < ng._facet.size(); ++i)
    out << ng._facet[i] << " "<< std::endl;
  out << std::endl;
  out << "####################################" << std::endl;
  
  return out;
}
#endif

#endif	/* __NGON_UNIT_H__ */

