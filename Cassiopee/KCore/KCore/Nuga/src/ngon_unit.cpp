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

#include "Nuga/include/ngon_unit.h"
#include "Nuga/include/IdTool.h"
#include "Nuga/include/Edge.h"
#include <assert.h>
#include <algorithm>
#include <map>

#ifdef DEBUG_NGON_UNIT
#include <iostream>
#endif

/// Methods Definitions
ngon_unit::ngon_unit(const E_Int* begin):_dirty(true)
{
  if (begin == nullptr) return;
  
  E_Int l = *(begin+1)+2;
  const E_Int* end = begin+l;
  _NGON.clear();
  _NGON.insert(_NGON.end(), begin, end);
  updateFacets();
}

///
ngon_unit::ngon_unit(const K_FLD::IntArray& cass_arr):ngon_unit(cass_arr.begin()){}

///
ngon_unit::ngon_unit(const E_Int* begin, E_Int sz, E_Int nbe):_dirty(true)
{
  const E_Int* end = begin+sz;
  _NGON.clear();
  _NGON.resize(2, 0);
  _NGON.insert(_NGON.end(), begin, end);
  _NGON[0]=nbe;
  _NGON[1]=sz;  
  updateFacets();
}

/// morse to ngon_unit
ngon_unit::ngon_unit(const E_Int* pngon, const E_Int* prange, E_Int rangesz)
{
  E_Int nbe = rangesz - 1;

  E_Int nbftot{ 0 };
  for (E_Int i = 0; i < nbe; ++i)
    nbftot += prange[i + 1] - prange[i];

  _NGON.resize(2 + nbftot + nbe, 0); // Like Cassiopee storage

  E_Int k{ 2 };
  for (E_Int i = 0; i < nbe; ++i)
  {
    E_Int nbf = prange[i + 1] - prange[i];

    _NGON[k++] = nbf;
    for (E_Int u = 0; u < nbf; ++u)
      _NGON[k++] = pngon[prange[i] + u];
  }

  _NGON[0] = nbe;
  _NGON[1] = _NGON.size() - 2;
  _dirty = true;
  updateFacets();
}

template <typename Container>
E_Int getPosFacets(const E_Int* data, E_Int facets_start, E_Int nb_facets/*nb of or PG or nodes*/, Container& posFacets)
{
  posFacets.clear();
  posFacets.resize(nb_facets, 1);
  for (E_Int i = 0; i < nb_facets; i++)
  {
    posFacets[i] = facets_start;
    facets_start += data[facets_start] + 1;
  }
  return 1;
}

///
void ngon_unit::updateFacets() const
{
  if (_dirty && !_NGON.empty())
  {
    getPosFacets(&_NGON[0], 2, _NGON[0], _facet);
    _dirty=false;
  }
}

///
ngon_unit& ngon_unit::operator=(const ngon_unit& ng)
{
  _NGON.clear();
  _facet.clear();
  _NGON = ng._NGON;
  _type = ng._type;
  _ancEs = ng._ancEs;
  _dirty=ng._dirty;
  if (!_dirty)
    _facet = ng._facet;
  return *this;
}

///
ngon_unit::ngon_unit(const ngon_unit& ngin)
{
  *this = ngin;
}

//move version
ngon_unit::ngon_unit(ngon_unit&& ngin):
  _NGON(std::move(ngin._NGON)), _facet(std::move(ngin._facet)),
  _type(std::move(ngin._type)), _ancEs(std::move(ngin._ancEs)), _dirty(ngin._dirty)
{
  ngin._NGON.clear();
  ngin._facet.clear();
  ngin._type.clear();
  ngin._ancEs.clear();
}

///
bool ngon_unit::attributes_are_consistent() const
{
  if (!_type.empty() && _type.size() != _facet.size())
    return false;
  if (_ancEs.cols() != 0 && _ancEs.cols() != (E_Int)_facet.size())
    return false;
  return true;
}

void ngon_unit::get_stride_extrema(E_Int& mins, E_Int& maxs) const
{
  updateFacets();

  E_Int nb_pg(size());
  maxs = -1;
  mins = -1;
  // Loop through the elements and increment face_count
  for (E_Int i = 0; i < nb_pg; ++i)
  {
    E_Int s = stride(i);
    mins = (mins == -1 || s < mins) ? s : mins;
    maxs = (maxs < s) ? s : maxs;
  }
}


///
void ngon_unit::add(int n, const E_Int* facet_ptr, int shift)
{
  // molecule here is one PH or PG
  if (_NGON.empty())
  {
    _NGON.resize(2, 0);
    _dirty = false;
  }

  _NGON[0]+=1;
  _NGON.push_back(n);
  _NGON.insert(_NGON.end(), facet_ptr, facet_ptr + n);
  _NGON[1] = _NGON.size() - 2;

  if (shift != 0)
  {
    E_Int from = _NGON.size() - n;
    K_CONNECT::IdTool::shift(_NGON, from, shift);
  }

  if (!_dirty)
    _facet.push_back(_NGON.size() - n - 1);
}

///
E_Int ngon_unit::get_facets_max_id() const
{
  updateFacets();
  
  E_Int nb_nodes, nb_pg(size()), id;  
  E_Int maxId=-1;
  // Loop through the elements and increment face_count
  for (E_Int i = 0; i < nb_pg; ++i)
  {
    nb_nodes = stride(i);
    for (E_Int j = 0; j < nb_nodes; ++j)
    {
      id = get_facet(i,j);
      maxId = (maxId < id) ? id : maxId;
    }
  }
  return maxId;
}

E_Int ngon_unit::get_facets_min_id() const
{
  updateFacets();

  E_Int nb_nodes, nb_pg(size()), id;
  E_Int minId = IDX_NONE;
  // Loop through the elements and increment face_count
  for (E_Int i = 0; i < nb_pg; ++i)
  {
    nb_nodes = stride(i);
    for (E_Int j = 0; j < nb_nodes; ++j)
    {
      id = get_facet(i, j);
      minId = (id < minId) ? id : minId;
    }
  }
  return minId;
}
///
void ngon_unit::append(const ngon_unit& cngon_in)
{
  if (cngon_in.size() == 0)
    return;

  if (_NGON.empty())
    *this=cngon_in;
  else
  {
    E_Int sz1 = _NGON.size();
    _NGON[0] += cngon_in._NGON[0];
    _NGON.insert(_NGON.end(), cngon_in._NGON.begin()+2, cngon_in._NGON.end());
    _NGON[1]=_NGON.size()-2;
    if (!_type.empty()) //transfer the info if it is already computed for the left hand side
      _type.insert(_type.end(), cngon_in._type.begin(), cngon_in._type.end());
    if (_ancEs.cols()) //transfer the info if it is already computed for the left hand side
      _ancEs.pushBack(cngon_in._ancEs);
      
    if (!_dirty && !cngon_in._dirty)
    {
      E_Int sz2 = _facet.size();
      _facet.insert(_facet.end(), cngon_in._facet.begin(), cngon_in._facet.end());
      K_CONNECT::IdTool::shift(_facet, sz2, sz1-2); //tail shift
    }
    else
      _dirty=true;
  }
 
#ifdef DEBUG_NGON_UNIT
  assert(this->attributes_are_consistent());
#endif
}

void ngon_unit::append(const ngon_unit& cngon_in, const E_Int* first, E_Int nb_elts)
{
  if (!first)
    return;
  if (cngon_in.size() == 0)
    return;
  
  for (E_Int i=0; i < nb_elts; ++i)
    __add(cngon_in, *(first+i)-1);
  
#ifdef DEBUG_NGON_UNIT
  assert(this->attributes_are_consistent());
#endif
  
}

///
void ngon_unit::append(const Vector_t<E_Int>& NGON)
{
  if (_NGON.empty())
    _NGON=NGON;
  else
  {
    _NGON[0] += NGON[0];
    _NGON.insert(_NGON.end(), NGON.begin()+2, NGON.end());
    _NGON[1]=_NGON.size()-2;
  }
  _dirty=true;
}

///
void ngon_unit::__add (const ngon_unit& ng, E_Int ith)
{
  if (_NGON.empty())
  {
    _NGON.resize(2,0);
    _dirty=false;
  }
  
  _NGON[0]++;
  E_Int sz = ng.stride(ith);
  const E_Int* ptr = ng.get_facets_ptr(ith);
  _NGON.push_back(sz);
  _NGON.insert(_NGON.end(), ptr, ptr+sz);
  _NGON[1]=_NGON.size()-2;
  if (!_dirty )
    _facet.push_back(_NGON.size()-1-sz);
  if (!ng._type.empty()) //transfer externality if the info is already computed for the left hand side. fixme hpc
    _type.push_back(ng._type[ith]);
  if (ng._ancEs.cols()) //transfer ancestors if the info is already computed for the left hand side. fixme hpc
    _ancEs.pushBack(ng._ancEs.col(ith), ng._ancEs.col(ith)+2);
}

///
E_Int ngon_unit::remove_entities (const Vector_t<E_Int>& to_remove, Vector_t<E_Int>& nids) //0-based nids
{
  nids.clear();

  if (to_remove.empty())
    return 0;

  E_Int sz(size()), count(0);
  std::set<E_Int> remset(to_remove.begin(), to_remove.end());//fixme
  std::set<E_Int>::const_iterator remend(remset.end());

  nids.resize(sz, IDX_NONE);

  updateFacets();

  ngon_unit ngu;

  for (E_Int i = 0; (i < sz); ++i)
  {
    if (remset.find(i) != remend)
      continue;

    ngu.__add(*this, i);
    nids[i]=count++;
  }
  
  //transfer externality
  ngu.compact_attributes(*this, nids);

  *this=ngu;
  updateFacets();
  
#ifdef DEBUG_NGON_UNIT
  assert(this->attributes_are_consistent());
#endif
  
  return sz-size();
}

///
E_Int ngon_unit::remove_facets(const Vector_t<E_Int>& nfacids, Vector_t<E_Int>& nids, E_Int min_facet_nb) //0-based nids
{

  ngon_unit ngu;
  Vector_t<E_Int> molecule;
  
  updateFacets();
  
  E_Int nb_elts = size();
  
  nids.clear(); //new element ids
  nids.resize(nb_elts, IDX_NONE);

  E_Int count(-1), degen(0);
  //std::cout << "nb elts : " << PHs.size() << std::endl;
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    E_Int nb_facets = this->stride(i);
    if (nb_facets == 0) //fixme : required ?
    {
#ifdef DEBUG_NGON_UNIT
      std::cout << i << "-th entity is stored with 0 facet !" << std::endl;
#endif
      continue;
    }

    const E_Int* facets = this->get_facets_ptr(i);
    
    molecule.clear();
    molecule.push_back(nb_facets);
    
    for (E_Int f = 0; f < nb_facets; ++f)
    {
      E_Int Fi = *(facets+f) - 1;

#ifdef DEBUG_NGON_UNIT
      if (Fi < 0 || Fi >= nfacids.size())
      {
        std::cout << "elt : " <<  f << std::endl;
        std::cout << "elt stride : " << nb_facets << std::endl;
        std::cout << "BOUM : " << Fi << std::endl;
        return 0;
      }
#endif
      if (nfacids[Fi] == IDX_NONE)
        continue;
      molecule.push_back(nfacids[Fi]+1);
    }
    molecule[0] = molecule.size() - 1;
    if (molecule[0]>min_facet_nb)//i.e has at least one valid facet (a PH has at least one PG, a PG has at least 3 nodes).
    {
      ngu.add(molecule);
      nids[i]=++count;
    }
    else if (min_facet_nb > 0)
    {
      ++degen;

#ifdef DEBUG_NGON_UNIT
      std::cout << "degen ! " <<  i << std::endl;
#endif
    }
#ifdef DEBUG_NGON_UNIT
    else
    {
      std::cout << "warning ! " <<  i << std::endl;
    }
#endif
  }
  
  //transfer externality
  ngu.compact_attributes(*this, nids);

  *this = ngu;
  
#ifdef DEBUG_NGON_UNIT
  assert(this->attributes_are_consistent());
#endif

  return (nb_elts - ngu.size());// the diff tells if some there wer some degen
}
  
///
void ngon_unit::shift(E_Int val, E_Int from)
{
  if (val == 0) return;

  updateFacets();

  size_t sz(_facet.size()), nb_nodes;  
  for (size_t i = from; i < sz; ++i)
  {
    nb_nodes = stride(i);
    E_Int* ptr = get_facets_ptr(i);
    for (size_t j = 0; j < nb_nodes; ++j)
    {
      E_Int& fj = *(ptr+j);
      if (fj != IDX_NONE) fj += val;
    }
  }
}

///set all facet values to IDX_NONE (useful for building a neighbor table)
void
ngon_unit::reset_facets()
{
  updateFacets();
  
  size_t sz(_facet.size()), nb_nodes;  
  for (size_t i = 0; i < sz; ++i)
  {
    nb_nodes = stride(i);
    E_Int* ptr = get_facets_ptr(i);
    for (size_t j = 0; j < nb_nodes; ++j)
      *(ptr+j) = IDX_NONE;
  }
}

///
void ngon_unit::extract
(const Vector_t<E_Int>& indices, ngon_unit& ng_out, Vector_t<E_Int>& oldIds, E_Int idx_start) const
{
  // 0-based indices
  ng_out.clear();
  oldIds.clear();
  oldIds.reserve(indices.size());
  updateFacets();
  
  for (size_t i = 0; i < indices.size(); ++i)
  {
    ng_out.__add(*this, indices[i] - idx_start);
    oldIds.push_back(indices[i]);
  }

  ng_out.updateFacets();
}

///
void ngon_unit::extract
(const E_Int* ptr, E_Int n, ngon_unit& ng_out, Vector_t<E_Int>& oids) const
{
  // 0-based indices
  ng_out.clear();
  oids.clear();
  oids.reserve(n);
  updateFacets();
  
  for (E_Int i = 0; i < n; ++i)
  {
    ng_out.__add(*this, *(ptr+i)-1);
    oids.push_back(*(ptr+i)-1);
  }
  ng_out.updateFacets();
}

///
void ngon_unit::extract_of_type (E_Int FLAG, ngon_unit& ng_out, Vector_t<E_Int>& oldIds) const
{
#ifdef DEBUG_NGON_UNIT
  assert ( attributes_are_consistent());
#endif
  
  Vector_t<E_Int> indices;
  get_indices_of_type(FLAG, indices);
  ng_out.clear();
  extract(indices, ng_out, oldIds);//get E32 skin faces
}

///
void ngon_unit::get_indices_of_type(E_Int FLAG, Vector_t<E_Int>& indices) const
{
  // 0-based indices
  for (size_t i = 0; i < _type.size(); ++i)
    if (_type[i] == FLAG)indices.push_back(i);
}

///
void ngon_unit::flag_indices_of_type (E_Int FLAG, Vector_t<bool>& flag) const
{

  flag.clear();
  flag.resize(_type.size(), false);
  for (size_t i = 0; i < _type.size(); ++i) flag[i] = (_type[i] == FLAG);
}

///
void ngon_unit::find_elts_with_facets(const std::set<E_Int>& facets, std::vector<E_Int> & elts) const
{
  elts.clear();
  
  updateFacets();
  
  E_Int nb_facets, nb_elts(size()), id;
  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      if (facets.find(id) != facets.end())
      {
        elts.push_back(i);
        break;
      }
    }
  }  
}

///
void ngon_unit::unique_indices(std::vector<E_Int>& indices) const
{
  std::set<E_Int> tmp;
  E_Int nb_facets, nb_elts(size()), id;
  
  indices.clear();
  
  updateFacets();
  
  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      tmp.insert(id);
    }
  }
  
  indices.insert(indices.end(), tmp.begin(), tmp.end());
}

/// NOT CONSECUTIVE : (N0, N1, N1, N2, N3, N3, N3, N4, N1, N7...) => (N0, N1, N2, N3, N4, N7...)
//any doublon occuring is removed => not ok when facets are node (PG destructive), ok when facet is PG
E_Int ngon_unit::remove_duplicated()
{
  std::set<E_Int> stmp;
  E_Int nb_facets, nbf, nb_elts(size()), id, nbe(0);
  ngon_unit ngtmp;
  
  updateFacets();
  
  ngtmp.clear(); 
  ngtmp._NGON.resize(2, 0);

  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    
    if (nb_facets == 1)
    {
      ngtmp._NGON.push_back(1);
      ngtmp._NGON.push_back(get_facet(i,0));
      ++nbe;
      continue;
    }
    
    stmp.clear();
    nbf=0;
    ngtmp._NGON.push_back(0);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      //std::cout << "i / id " << i << " : " << id << std::endl;
      if (stmp.insert(id).second)
      {
        ngtmp._NGON.push_back(id);
        ++nbf;
      } 
    }
    
    if (nbf)
    {
      ++nbe;
      ngtmp._NGON[ngtmp._NGON.size()-1-nbf]=nbf;//set the nb of facets
    }
    else
      ngtmp._NGON.pop_back();
  }
  
  ngtmp._NGON[0]=nbe;
  ngtmp._NGON[1]=ngtmp._NGON.size()-2;
  
  E_Int prev_sz = this->_NGON[1];
  
  // hack to preserve info (type and history) : must work because this function does not change the number of facets
  if (!_type.empty())
    ngtmp._type=this->_type;
  if (_ancEs.cols())
    ngtmp._ancEs=this->_ancEs;
  
  *this=ngtmp;
  
#ifdef DEBUG_NGON_UNIT
  assert(this->attributes_are_consistent());
#endif
  
  updateFacets();
  
  return this->_NGON[1] - prev_sz;
}

/// CONSECUTIVE (N0, N1, N1, N2, N3, N3, N3, N4, N1, N7...) => (N0, N1, N2, N3, N4, N1, N7...) :
//=> ok when facets are node (PG safe), not ok when facet is PG because some duplicates remain
E_Int ngon_unit::remove_consecutive_duplicated()
{
  E_Int nb_facets, nbf, nb_elts(size()), id, prev_id, nbe(0);
  ngon_unit ngtmp;
  
  updateFacets();
  
  ngtmp.clear(); 
  ngtmp._NGON.resize(2, 0);
  
  std::vector<E_Int> nids(nb_elts, IDX_NONE);
  
  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    
    if (nb_facets == 1)
    {
      ngtmp._NGON.push_back(1);
      ngtmp._NGON.push_back(get_facet(i,0));
      nids[i] = nbe++;
      continue;
    }
    
    nbf=0;
    ngtmp._NGON.push_back(0);
    prev_id=get_facet(i,nb_facets-1);//initialized with the last one to see if not equal to the first one
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      //std::cout << "i / id " << i << " : " << id << std::endl;
      if (id != prev_id)
      {
        ngtmp._NGON.push_back(id);
        prev_id=id;
        ++nbf;
      } 
    }
    if (nbf)
    {
      nids[i] = nbe++;
      ngtmp._NGON[ngtmp._NGON.size()-1-nbf]=nbf;//set the nb of facets
    }
    else
      ngtmp._NGON.pop_back();
  }
  
  ngtmp._NGON[0]=nbe;
  ngtmp._NGON[1]=ngtmp._NGON.size()-2;
  
  E_Int prev_sz = (this->_NGON[1]);
  
  //transfer attributes
  ngtmp.compact_attributes(*this, nids);
  
  *this=ngtmp;
  
#ifdef DEBUG_NGON_UNIT
  assert(this->attributes_are_consistent());
#endif
  
  return ngtmp._NGON[1] - prev_sz;
}

///
/*
if PGS : min_nb_facets = 1 for lineic, 3 otherwise
if PHs : min_nb_facets = 2 for lineic, 4 otherwise
*/
//=============================================================================
void ngon_unit::get_degenerated(E_Int min_nb_facets, Vector_t<E_Int>& indices) 
{
  updateFacets();
  
  E_Int ngon_dim = min_nb_facets - 1;
  E_Int s, nb_elts(size());
  std::set<E_Int> unic;
  indices.clear();
  
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    s = stride(i);
    
    if (s < min_nb_facets)
    {
      indices.push_back(i);
      continue;
    }

    E_Int* facets = get_facets_ptr(i);
    
    // remove duplicates
    unic.clear();
    unic.insert(get_facets_ptr(i), get_facets_ptr(i)+s);
    E_Int news=unic.size();

    if (news < min_nb_facets)
      indices.push_back(i);
    else if ( (ngon_dim == 2) && (news < s) ) //has some duplicates : just some pinches or has some duplicated egdes ?
    {
      std::set<K_MESH::NO_Edge> sedges;
      
      for (E_Int n = 0; n < s; ++n)
      {
        int Fi     = facets[n];
        int Fip1   = facets[(n + 1) % s];

        sedges.insert(K_MESH::NO_Edge(Fi, Fip1));
      }

      if (sedges.size() < 3)
        indices.push_back(i);
    }
  }
}

///
void ngon_unit::change_indices (const Vector_t<E_Int>& nIds, E_Int idx_start)
{
  if (nIds.empty()) return;

  updateFacets();
  
  E_Int nb_facets, nb_elts(size());
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      E_Int& id = get_facet(i,j);
      if (id == IDX_NONE) continue;
      E_Int nid = nIds[id-idx_start];
      id = (nid != IDX_NONE) ? nid + idx_start : IDX_NONE;
    }
  }
}

bool ngon_unit::is_fixed_stride(E_Int& strd) const
{
  updateFacets();
  
  strd=stride(0);
  
  E_Int nb_elts(size());
  for (E_Int i = 1; i < nb_elts; ++i)
  {
    if (strd != stride(i))
      return false;
  }
  return true;
}

/// convert a fixed-stride storage to a ngon_unit storage
void ngon_unit::convert_fixed_stride_to_ngon_unit(const K_FLD::IntArray&cnt, E_Int shift, ngon_unit& nguo)
{
  //
  E_Int nb_elts = cnt.cols();
  E_Int stride = cnt.rows();

  E_Int sz = (stride + 1)*nb_elts + 2;
  nguo._NGON.resize(sz);

  nguo._NGON[0] = nb_elts;
  nguo._NGON[1] = sz - 2;

  E_Int pos = 2;
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    const E_Int* p = cnt.col(i);
    nguo._NGON[pos++] = stride;
    for (E_Int j = 0; j < stride; ++j)
      nguo._NGON[pos++] = *(p + j) + shift;
  }

  nguo.updateFacets();
}

void ngon_unit::convert_ngon_unit_to_fixed_stride(const ngon_unit& ngui, E_Int shift, K_FLD::IntArray& cnto)
{
  if (ngui.size() == 0) return;
  //
  E_Int stride;
  if (!ngui.is_fixed_stride(stride))
    return;
  
  E_Int nb_elts = ngui.size();
  
  cnto.clear();
  cnto.resize(stride, nb_elts);

  for (E_Int i = 0; i < nb_elts; ++i)
  {
    const E_Int* p = ngui.get_facets_ptr(i);
    
    for (E_Int j=0; j < stride; ++j)
      cnto(j,i)=*(p+j)-shift;
  }
}

/// transfer attribute (reduction)
void ngon_unit::compact_attributes(const ngon_unit& ngi, const Vector_t<E_Int>& nids)
{
  if (ngi._type.empty() && (ngi._ancEs.cols() == 0)) return;
  
  //transfer externality
  if (!ngi._type.empty())
  {
    _type = ngi._type;
    K_CONNECT::IdTool::compact(_type, nids);
  }
  if (ngi._ancEs.cols() != 0)
  {
    _ancEs = ngi._ancEs;
    K_CONNECT::IdTool::compact(_ancEs, nids);
  }
  
  this->updateFacets();

#ifdef DEBUG_NGON_UNIT
    assert(attributes_are_consistent());
#endif
}

/// transfer attribute (reduction)
void ngon_unit::compact_attributes(const ngon_unit& ngi, const Vector_t<bool>& keep)
{
  if (ngi._type.empty() && (ngi._ancEs.cols() == 0)) return;
  
  //transfer externality
  K_CONNECT::keep<bool> pred(keep);
  
  if (!ngi._type.empty())
  {
    _type = ngi._type;
    K_CONNECT::IdTool::compress(_type, pred);
  }
  if (ngi._ancEs.cols() != 0)
  {
    _ancEs = ngi._ancEs;
    K_CONNECT::IdTool::compress(_ancEs, pred);
  }
  
  this->updateFacets();
  
#ifdef DEBUG_NGON_UNIT
    assert(attributes_are_consistent());
#endif
}

/// transfer attribute (spread)
void ngon_unit::spread_attributes(const ngon_unit& ngi, const Vector_t<E_Int>& oids)
{
  if (ngi._type.empty() && (ngi._ancEs.cols() == 0)) return;

  size_t sz = oids.size();

  if (sz == 0) return;

#ifdef DEBUG_NGON_UNIT
  assert(ngi.attributes_are_consistent());
  assert (attributes_are_consistent());
  assert(_facet.size() == sz);
#endif

  //transfer externality
  if (!ngi._type.empty())
  {
    _type.resize(oids.size(), IDX_NONE);
    for (size_t i = 0; i < sz; ++i)
      _type[i] = ngi._type[oids[i]];
  }
  if (ngi._ancEs.cols() != 0)
  {
    _ancEs.resize(2, sz, IDX_NONE);
    for (size_t i = 0; i < sz; ++i)
    {
      _ancEs(0, i) = ngi._ancEs(0, oids[i]);
      _ancEs(1, i) = ngi._ancEs(1, oids[i]);
    }
  }

  this->updateFacets();

#ifdef DEBUG_NGON_UNIT
  assert(attributes_are_consistent());
#endif
}

void ngon_unit::sort_by_type(Vector_t<E_Int>& nids, Vector_t<E_Int>& oids) const
{
  using type_to_oid_t = std::pair<E_Int, E_Int>;
  
  size_t sz = _type.size();
  
  Vector_t<type_to_oid_t> type_to_oid(sz);
  
  for (size_t i=0; i < sz; ++i)
  {
    type_to_oid[i].first = _type[i];
    type_to_oid[i].second = i;
  }
  
  std::sort(type_to_oid.begin(), type_to_oid.end());
  
  nids.clear();
  nids.resize(sz, IDX_NONE);
  oids.clear();
  oids.resize(sz, IDX_NONE);
  
  for (size_t i=0; i < sz; ++i){
    nids[type_to_oid[i].second] = i;
    oids[i] = type_to_oid[i].second;
  }
}
