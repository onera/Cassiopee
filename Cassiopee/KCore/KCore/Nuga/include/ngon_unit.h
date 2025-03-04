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

#ifndef __NGON_UNIT_H__
#define	__NGON_UNIT_H__

#include <vector>
#include <set>
#include <assert.h>

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/mem_chunk.hxx"
#include "Nuga/include/macros.h"

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
    explicit ngon_unit (const E_Int* begin);//begin is stored in Cassiopee format
    explicit ngon_unit(const K_FLD::IntArray& cass_arr);
    ///
    ngon_unit(const ngon_unit& ngin);
    ngon_unit(ngon_unit&& ngin);
    /// 
    ngon_unit(const E_Int* begin, E_Int sz, E_Int nbe);//begin is stored in Cassiopee format

    /// morse to ngon_unit
    ngon_unit(const E_Int* pngon, const E_Int* prange, E_Int rangesz);
        
    /// ngon_unit& operator=(const Vector_t<E_Int>& vNGON);
    ngon_unit& operator=(const ngon_unit& ngin);

    /// export functions
    // into morse
    template <typename vngu>
    void relay_mem(vngu& ngo)
    {
      ngo.release(); // delete current mem

      using INT_t = typename vngu::INT_t;
      INT_t nbe = size();
      INT_t sz = _NGON.size() - 2 - nbe; //removing cartouche and facets_nb entries

      // allocate memory for the morse data structure
      using memchunk_t = NUGA::mem_chunk<INT_t, NUGA::allocator<ngo.CALLOC>>;
      memchunk_t ngon(sz);
      memchunk_t range(nbe + 1); //one-passed-the-end

      // plug
      ngo.elts = ngon.data;
      //ngo.nelts = sz;
      ngo.idx_start = 1;
      ngo.range = range.data;
      ngo.nrange = range.alloc_sz;

      //_NGON/_facets => morse
      INT_t k{ 0 };
      for (INT_t i = 0; i < nbe; ++i)
      {
        const INT_t* facets = this->get_facets_ptr(i);
        const INT_t nbf = this->stride(i);

        ngo.range[i] = k;

        for (INT_t j = 0; j < nbf; ++j)
          ngo.elts[k++] = facets[j];
      }
      ngo.range[nbe] = sz;
    }

    ///
    E_Float operator[](const E_Int i){return _NGON[i];}
    /// Refresh the facets starts.
    void updateFacets() const;
    
    /// memory
    void reserve (E_Int sz){_NGON.reserve(sz); _facet.reserve(sz); _type.reserve(sz); _ancEs.reserve(2, sz);} //fixme : wrong for _NGON
    
    void expand_n_fixed_stride (E_Int n, E_Int strd){

      if (_NGON.empty())
      {
        _NGON.resize(2, 0);
        _dirty = true;
      }

      E_Int sz0 = (E_Int)_facet.size();
      _NGON.resize(_NGON.size() + n*(strd+1), strd) ; _facet.resize(sz0 + n) ; 
      if (!_type.empty()) _type.resize(_type.size() + n) ;
      if (_ancEs.cols()!=0) _ancEs.resize(2, _ancEs.cols() + n);
//      E_Int pos = _facet[sz0-1] + strd+1;
//      // set the facets pos
//      for (E_Int i=0; i< n; ++i, pos += strd+1)
//        _facet[sz0 + i] = pos;

      _NGON[0] += n;
      _NGON[1] = (E_Int)_NGON.size() - 2;
      // the following update has been added since the commented block above doesn't work in Basic element mode.
      _dirty=true;
      updateFacets();
    }

    void expand_n_ab_fixed_stride(E_Int nn, E_Int na, E_Int astrd, E_Int nb, E_Int bstrd) {

      /// compound  na*astrd, nb*bstrd...nn times

      if (_NGON.empty())
      {
        _NGON.resize(2, 0);
        _dirty = true;
      }

      E_Int sz_tot = (na * (astrd+1) + nb * (bstrd+1)) * nn;
      E_Int sz0    = (E_Int)_facet.size();
      E_Int ng_sz0 = (E_Int)_NGON.size();

      _NGON.resize(ng_sz0 + sz_tot); _facet.resize(sz0 + nn * (na + nb));
      
      // set the coorrect strides in NGON
      for (E_Int i = 0; i < nn; ++i)
      {
        E_Int pos = ng_sz0 + i * (na * (astrd + 1) + nb * (bstrd + 1));
        for (E_Int j = 0; j < na*(astrd + 1); ++j)
        {
          assert ((pos + j) > -1 && (pos + j) < E_Int(_NGON.size()));
          _NGON[pos + j] = astrd;
        }
        for (E_Int j = 0; j < nb*(bstrd + 1); ++j)
        {
          int v = pos + j + na * (astrd + 1);
          assert (v > -1 && v < E_Int(_NGON.size()));
          _NGON[v] = bstrd;
        }
      }
      
      if (!_type.empty()) _type.resize(_type.size() + nn * (na + nb));
      if (_ancEs.cols() != 0) _ancEs.resize(2, _ancEs.cols() + nn * (na + nb));
      
      _NGON[0] += nn * (na + nb);
      _NGON[1] = _NGON.size() - 2;
      
      _dirty = true;
      updateFacets();
    }

    // warning : pregnant has non-zero vals
    void expand_variable_stride(E_Int n, const E_Int* pregnant)
    {
      if (_NGON.empty())
      {
        _NGON.resize(2, 0);
        _dirty = true;
      }

      E_Int sz0 = _facet.size();

      E_Int nb_tot_child = 0;
      for (E_Int i = 0; i < n; ++i) nb_tot_child += pregnant[i];
      
      E_Int pos = _NGON.size();
      _NGON.resize(pos + nb_tot_child + n, IDX_NONE); 
      _facet.resize(sz0 + n);

      //put the nb_facets
      for (E_Int i = 0; i < n; ++i)
      {
        ASSERT_IN_VECRANGE(_NGON, pos)
        _NGON[pos] = pregnant[i];
        pos += pregnant[i] + 1;
      }

      if (!_type.empty()) _type.resize(_type.size() + n);
      if (_ancEs.cols() != 0) _ancEs.resize(2, _ancEs.cols() + n);

      _NGON[0] += n;
      _NGON[1] = _NGON.size() - 2;

      _dirty = true;
      updateFacets();

    }

    E_Int capacity(){return _NGON.capacity();}
    // Interrogations
    /// Returns the number of facets for the i-th PH(PG)  (ArrayAccessor interface)
    inline E_Int stride(E_Int i) const
    {
      ASSERT_IN_VECRANGE(_facet, i)
      ASSERT_IN_VECRANGE(_NGON, _facet[i])
      return _NGON[_facet[i]];
    }
    /// Returns the j-th element.  (ArrayAccessor interface)
    inline void getEntry(const E_Int& j, E_Int* pE) const
    {
      ASSERT_IN_VECRANGE(_facet, j)
      ASSERT_IN_VECRANGE(_NGON, _facet[j])

      for (E_Int k=0;k<_NGON[_facet[j]]/*stride(j)*/;++k)
      {
        ASSERT_IN_VECRANGE(_NGON, _facet[j]+1+k)
        pE[k]=_NGON[_facet[j]+1+k];/*get_facet(j,k)*/
      }
    }
    /// Returns the j-th facet of the i-th element
    inline E_Int get_facet(E_Int i, E_Int j) const 
    {
      ASSERT_IN_VECRANGE(_facet, i)
      ASSERT_IN_VECRANGE(_NGON, _facet[i]+1+j)
      return _NGON[_facet[i]+1+j];
    }
    
    inline E_Int& get_facet(E_Int i, E_Int j)
    {
      ASSERT_IN_VECRANGE(_facet, i)
      ASSERT_IN_VECRANGE(_NGON, _facet[i]+1+j)
      return _NGON[_facet[i]+1+j];
    }
    
    /// Returns the j-th facet of the i-th element
    inline const E_Int* get_facets_ptr(E_Int i) const
    {
      ASSERT_IN_VECRANGE(_facet, i)
      ASSERT_IN_VECRANGE(_NGON, _facet[i]+1)
      return &_NGON[_facet[i]+1];
    }

    inline E_Int* get_facets_ptr(E_Int i)
    {
      ASSERT_IN_VECRANGE(_facet, i)
      ASSERT_IN_VECRANGE(_NGON, _facet[i]+1)
      return &_NGON[_facet[i]+1];
    }
    
    inline const E_Int* begin(E_Int i) const
    {
      ASSERT_IN_VECRANGE(_facet, i)
      ASSERT_IN_VECRANGE(_NGON, _facet[i]+1)
      return &_NGON[_facet[i]+1];
    }//SYNOMYM
    
    /// Returns the number of entities
    inline E_Int size() const {return !_NGON.empty() ? _NGON[0] : 0;}
    inline E_Int getSize() { return size(); } // to match DynArray interface
    ///
    inline const E_Int* begin() const {return &_NGON[0];}
    ///
    E_Int get_facets_max_id() const ;
    ///
    E_Int get_facets_min_id() const;
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
    void extract (const Vector_t<E_Int>& indices, ngon_unit& ng_out, Vector_t<E_Int>& oldIds, E_Int idx_start=0) const;
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
    /// WARNING : nids is 0-based, "this" is kept 1-based
    void change_indices (const Vector_t<E_Int>& nIds, E_Int idx_start=1);
    ///
    void get_degenerated(E_Int min_nb_facets, Vector_t<E_Int>& indices);
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
    void shift(E_Int shift, E_Int from=0);
    
    /// 
    void sort_by_type(Vector_t<E_Int>& nids, Vector_t<E_Int>& oids) const;
    
    ///
    template <typename T>
    static void shift(Vector_t<T>& vec, const T& val){for (size_t i = 0; i < vec.size(); ++i)vec[i]+=val;}
    ///
    void reset_facets();//set all facet values to IDX_NONE (useful for building a neighbor table)
    
    /// warning : need a call to updateFacets afterwards
    template <template<typename, typename> class Container_t, typename Element,typename Allocator_t>
    void add (const Container_t<Element, Allocator_t>& molecule);

    /// warning : need a call to updateFacets afterwards
    void add(int n, const E_Int* facet_ptr, int shift = 0);

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
    mutable Vector_t<E_Int>  _type;
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
  assert (_type.empty() || (_type.size() == sz));
  
  oids.clear();
  
  updateFacets();
  
  nids.clear();
  nids.resize(sz, IDX_NONE);
  
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
  for (size_t i = 0; i < sz; ++i)
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


#endif	/* __NGON_UNIT_H__ */

