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

#ifndef ZONE_T_HXX
#define ZONE_T_HXX

#include <map>
#include <vector>

#ifdef DEBUG_ZONE_T
#include "Nuga/include/medit.hxx"
#endif


#define Vector_t std::vector

#define PH_INNER_COL -1
#define UNSET_COL     0
#define PH_GHOST      99

#define PG_COL_RANGE   10000 // <=> sum(n_bcs(zi)) over i zones must be <= 1000

#define PG_INNER_COL   0  * PG_COL_RANGE
#define PG_JOIN_COL    1  * PG_COL_RANGE
#define PG_LAY1_IN_COL 2  * PG_COL_RANGE
#define PG_LAY1_BC_COL 3  * PG_COL_RANGE
#define PG_BC          4  * PG_COL_RANGE
#define PG_LAY2_IN_COL 5  * PG_COL_RANGE
#define PG_LAY2_BC_COL 6  * PG_COL_RANGE
#define OTHER_LAYS_COL 7  * PG_COL_RANGE  //above 2nd lauyer , anything is unsorted for now
#define PG_GHOST       99 * PG_COL_RANGE

#define NEIGHBOR1(PH0, F2E, shift, Fi) ((F2E[Fi] == PH0 ) ? F2E[Fi+shift] : F2E[Fi])

inline bool is_a_bc(E_Int type, E_Int nlayers)
{
  if ((PG_LAY1_BC_COL <= type && type < PG_LAY1_BC_COL + PG_COL_RANGE) && (nlayers >= 2))
    return true;
  if ((PG_LAY2_BC_COL <= type && type < PG_LAY2_BC_COL + PG_COL_RANGE) && (nlayers >=3))
    return true;
  if (PG_BC <= type && type < PG_BC + PG_COL_RANGE)
    return true;
  return false;
}

namespace NUGA
{
  

///
template <typename crd_t, typename ngo_t>
class zone_t
{
  public:
  E_Int _id;
  
  public:
    
    enum eEntType { INTERNAL = -1, JOIN = 0, BC = 1};
    using join_t = std::map<E_Int, std::pair<zone_t*, Vector_t<E_Int> > >;
    using bc_t    = std::map<E_Int, Vector_t<E_Int> >; 
    
    zone_t():_id(IDX_NONE),_F2E_NONE(0){}
    zone_t(E_Int id, crd_t& crd, ngo_t& ng, E_Int* f2e, E_Int F2E_NONE):_id(id), _crd(crd), _ng(ng), _F2E_NONE(F2E_NONE){
      
      E_Int nb_pgs = ng.PGs.size();
      _F2Es.insert(_F2Es.begin(), f2e, f2e + 2*nb_pgs);
    }
    
    zone_t(const zone_t& z):_id(z._id), _crd(z._crd), _ng(z._ng), _F2Es(z._F2Es), _F2E_NONE(z._F2E_NONE)
    {   
      copy_boundaries(z); 
    }
    
    zone_t& operator=(const zone_t& z)
    {
      _id = z._id;
      _crd = z._crd;
      _ng = z._ng;
      _F2Es = z._F2Es;
      _F2E_NONE = z._F2E_NONE;
            
      copy_boundaries(z);

      return *this;
      
    }
    
    static void join_zones(zone_t& z1, E_Int* r12, zone_t& z2, E_Int* r21, E_Int nbj, E_Int ojid)
    {
      // find available id : the lowest minus 1 avail in each zone
      E_Int jid1 = (z1._joins.empty()) ? 0 : z1._joins.begin()->first;
      E_Int jid2 = (z2._joins.empty()) ? 0 : z2._joins.begin()->first;
      E_Int jid = - (std::max(fabs(1.*jid1), fabs(1.*jid2)) + 1) ; // negative and 1-based to not overwrite any join when changing colors in init_pgs_color
      
      Vector_t<E_Int> tmp;
      tmp.insert(tmp.begin(), r12, r12 + nbj);
      
      z1._joins[jid] = std::make_pair(&z2, tmp);

      tmp.clear();
      tmp.insert(tmp.begin(), r21, r21 + nbj);
      
      z2._joins[jid] = std::make_pair(&z1, tmp);

      z1._rac_inId2outId[jid] = z2._rac_inId2outId[jid] = ojid;
    }
    
    //E_Int nb_joins() { return _zone_to_join.size();}
    
    bool empty() const {return (_ng.PHs.size() == 0);}
    
    bool is_disconnected_from(const zone_t* z, Vector_t<E_Int>& racids){ 
      racids.clear();
      for (const auto& itrz : _joins)
        if (itrz.second.first == z) { racids.push_back(itrz.first); }
      return racids.empty();
    }
    
    
    void append(zone_t& z, bool keep_joins);
    
    bool set_ph_type(E_Int PHi, E_Int val) { 
      if (_ng.PHs._type[PHi] == UNSET_COL){
        _ng.PHs._type[PHi] = val;
        return true;
      }
      return false;
    }
        
    const Vector_t<E_Int>& get_join_pg_list(E_Int jid) const { 
           
      auto it = _joins.find(jid);
      assert (it != _joins.end());
     
      return it->second.second;
    }
    
    Vector_t<E_Int>& get_f2e() { return _F2Es;}
    
    ngon_type& get_ng() { return _ng;}
    
    void change_joins_zone(zone_t* from, zone_t* to) 
    {
      for (auto& r2z : _joins)
      {
        zone_t* zj = r2z.second.first;
        if (zj == from) //that_z and not ztmp because join pointers are not update with operator=
          r2z.second.first = to;
      }
    }
    
    void compress_join(E_Int jid, const Vector_t<E_Int>& nids) 
    {

      auto j = _joins.find(jid);
      assert (j != _joins.end());

      K_CONNECT::valid pred(nids);
      K_CONNECT::IdTool::compress(j->second.second, pred);
      if (j->second.second.empty())
        _joins.erase(jid);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    //void extract_positive_types(zone_t<crd_t, ngo_t, true>& z);
    E_Int reduce_to_positive_types();
    
    void sort_by_type() ;
    void insert_ghosts_on_bcs(E_Int type, E_Int nb_layers);
    
    void init_pg_types(E_Int& color);
    void set_join_types(E_Int& color);
    void set_BC_types(E_Int& color);
    void reset_ph_types(E_Int t);
        
    void reorient_f2e_boundaries();
    
    void set_f2e_across_joins();
    
    void flag_connected_to_type(E_Int val, E_Int connec_val, Vector_t<E_Int> & zones4nextpass);
    
    void copy_boundaries(const zone_t& z);
    
    void update_boundaries(const Vector_t<E_Int>& pgnids, E_Int idx_start);
    
    void add_boundary(E_Int obcid, const E_Int* ids, E_Int n);

    void sort_pointLists();
    
    void set_pg_colors();

    void __distinguish_bcs();
    
    void color_ranges(std::vector<E_Int>& PGcolors, std::vector<E_Int>& PHcolors);
    
#ifdef DEBUG_ZONE_T
void draw_boundaries()
{
  Vector_t<E_Int> ids;
  E_Int nb_pgs = _ng.PGs.size();
  for (size_t i=0; i < nb_pgs; ++i)
  {
    if (_F2Es[i+nb_pgs] == _F2E_NONE)ids.push_back(i);
  }
  
  //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs("bound", _crd, _ng.PGs, ids, false);
}
#endif
    
//private:
    // Mesh
    K_FLD::FloatArray _crd;
    ngon_type         _ng;
    
    // F2E
    Vector_t<E_Int>                      _F2Es;        // separated-fields array
    E_Int _F2E_NONE;

    // BCs and Joins
    join_t                            _joins;
    bc_t _bcs;      // Boundaries
    
    std::map <E_Int, E_Int>           _rac_inId2outId; // Internal id to Original id 
    std::map <E_Int, E_Int>           _bc_inId2outId;  // Internal id to Original id 
    
    
    
};

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::reset_ph_types(E_Int t)
{
  // init all PH to 0 for now
  E_Int nb_phs = _ng.PHs.size();
  _ng.PHs._type.clear();
  _ng.PHs._type.resize(nb_phs, t);
}

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::set_join_types(E_Int& color)
{
  join_t new_joins = _joins;

  Vector_t<E_Int> to_remove;
  for (auto itbi : _joins)
  {
    E_Int ojid = itbi.first;
        
    if (ojid >= 0) {//need only to set the color to the PGs (already partly processed by its join mate)
    
      E_Int nbj = itbi.second.second.size();
      
      for (E_Int i = 0; i < nbj; ++i) 
      {
        E_Int id = itbi.second.second[i] - 1;
        assert (id > -1 && id < _ng.PGs._type.size());
        _ng.PGs._type[id] = ojid;
      }
      continue;
    }
    
    auto itrz = _joins.find(ojid);
    assert (itrz != _joins.end());
    zone_t& opp_z = *itrz->second.first;

    to_remove.push_back(ojid);
    
    // change the key in the map in both sides
    const Vector_t<E_Int>& ids = itbi.second.second;
    
    E_Int nbj = ids.size();
    
    for (E_Int i = 0; i < nbj; ++i)
    {
      E_Int id = ids[i] - 1;
      assert (id > -1 && id < _ng.PGs._type.size());
      _ng.PGs._type[id] = color;
    }
    
    new_joins[color] = std::make_pair(&opp_z,ids);

    _rac_inId2outId[color] = _rac_inId2outId[ojid];
    _rac_inId2outId.erase(ojid);
        
    auto oitrz = opp_z._joins.find(ojid);
    assert (oitrz != _joins.end());
    
    const Vector_t<E_Int>& opp_ids = oitrz->second.second;
  
    opp_z._joins[color] = std::make_pair(this, opp_ids);
    opp_z._joins.erase(ojid);
    
    ++color;
  }
  
  if (!to_remove.empty())
  {
    _joins = new_joins;
    for (size_t i = 0; i < to_remove.size(); ++i)
      _joins.erase(to_remove[i]);
  }
}

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::set_BC_types(E_Int& color)
{
  bc_t new_bc_to_ids;
  std::map<E_Int, E_Int> new_bc_inId2outId;
  //std::cout << "set_BC_types" << std::endl;
  for (auto& itbi : _bcs)
  {
    const Vector_t<E_Int>& bids = itbi.second;
    for (size_t i = 0; i < bids.size(); ++i) _ng.PGs._type[bids[i] - 1] = color;
    new_bc_to_ids[color] = bids;
    new_bc_inId2outId[color] = _bc_inId2outId[itbi.first];
    ++color;
  }
  
  _bcs = new_bc_to_ids;
  _bc_inId2outId = new_bc_inId2outId;
  
}

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::init_pg_types(E_Int& color)
{    
  //std::cout << "init_pg_types 1 : color : " << color << std::endl;
  // init PGs type based on joins and BCs : INTERNAL < JOIN <  IN_LAYER1 < BC_LAYER1 < BC < IN_LAYER2 < BC_LAYER2 < ...
  size_t nb_pgs = _ng.PGs.size();
  _ng.PGs._type.clear();
  _ng.PGs._type.resize(nb_pgs, color++);
  
  set_join_types(color); 
  set_BC_types(color);
}

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::reorient_f2e_boundaries()
{
  E_Int sz = _F2Es.size()/2;
  for (E_Int i = 0; i < sz; ++i)
    if (_F2Es[i] == _F2E_NONE)
      std::swap(_F2Es[i], _F2Es[i+sz]);
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::set_f2e_across_joins()
{
  
  for (auto itri : _joins)
  {
    E_Int id = itri.first;
    zone_t* that_z = itri.second.first;       // joined zone
    auto & this_ids = itri.second.second;      // list of join PG on "this" side
    
    if (this_ids.size() == 0) continue;
    
    auto & this_f2e = this->_F2Es;    // F2E on "this" side
    
    const Vector_t<E_Int>& that_ids = that_z->get_join_pg_list(id);
    Vector_t<E_Int>& that_f2e = that_z->get_f2e();
    E_Int f2esz = that_f2e.size() / 2;
    
    for (size_t i = 0; i < this_ids.size(); ++i)
    {
      E_Int thisPGi = this_ids[i] - 1;
      
      E_Int this_leftPH = this_f2e[thisPGi]; // index start must be 1
      E_Int lid = this_leftPH;
      
      E_Int thatPGi = that_ids[i] - 1;
      E_Int& that_rightPH = that_f2e[thatPGi + f2esz]; // index start must be 1
      assert (that_rightPH == this->_F2E_NONE); // mandatory previous call to reorient_f2e_boundaries
      
      that_rightPH = - lid;
    }
  }
}

///
template <typename crd_t, typename ngo_t>
E_Int zone_t<crd_t, ngo_t>::reduce_to_positive_types()
{
  /// 1. Reduce PHs
  Vector_t<E_Int> phnids;
  {
    K_CONNECT::strictly_positive pred(_ng.PHs._type);
    Vector_t<E_Int> phoids;
    ngon_unit phs;
    _ng.PHs.extract_by_predicate(pred, phs, phoids, phnids);
    _ng.PHs = phs;
  }
  
  if (_ng.PHs.size() == 0) return 0;
  
  /// 2. Reduce PGs => pgnids
  
  // 2.1 flag right PGs 
  Vector_t<E_Int> keepPG;
  _ng.flag_facets_of_elts(1, keepPG);

  // 2.3 compress PGs with pgnids
  ngon_unit pgs;
  Vector_t<E_Int> pgnids, pgoids;
  K_CONNECT::keep<E_Int> pred(keepPG);

  _ng.PGs.extract_by_predicate(pred, pgs, pgoids, pgnids);
  _ng.PGs = pgs;
  _ng.PHs.change_indices(pgnids);
  
      
  // 2.4 compress and update accordingly F2E
  {
    Vector_t<E_Int> tmp(keepPG);
    keepPG.insert(keepPG.end(), tmp.begin(), tmp.end());//double it for F2E
    K_CONNECT::keep<E_Int> pred1(keepPG); 
    K_CONNECT::IdTool::compress(_F2Es, pred1);
    
    //E_Int nb_pgs = _F2Es.size() / 2;
    
    for (size_t i=0; i < _F2Es.size(); ++i) 
    {
      E_Int PHi = _F2Es[i];
      if (PHi == _F2E_NONE || PHi < 0) continue; // BC or JOIN
      _F2Es[i] = (phnids[PHi - 1] ==  IDX_NONE) ? _F2E_NONE :  phnids[PHi - 1] + 1;
      if (_F2Es[i] == IDX_NONE) _F2Es[i] = _F2E_NONE;
    }
  }
    
  /// 3. joins and BCs with pgnids
  update_boundaries(pgnids, 1);
  
  /// 4. compact nodes
  ngon_type::compact_to_used_nodes(_ng.PGs, _crd);
  
  return _ng.PHs.size();
  
}

///
/// \param val
/// \param connec_val
/// \param traversed_zones
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::flag_connected_to_type
(E_Int val, E_Int connec_val, Vector_t<E_Int> & zones4nextpass)
{
  size_t nb_phs = _ng.PHs.size();
  
  size_t nb_pgs = _F2Es.size() / 2;
    
  for (size_t i=0; i < nb_phs; ++i)
  {
    if (_ng.PHs._type[i] != connec_val ) continue; // consider only elements attached to those having type "connec_val"

    E_Int PHi = i+1;

    const E_Int * pPGi = _ng.PHs.get_facets_ptr(i);
    E_Int nb_faces = _ng.PHs.stride(i);
    //std::cout << "nb faces : " << nb_faces << std::endl;
    for (E_Int n=0; n < nb_faces;++n)
    {
      E_Int PGi = *(pPGi + n) - 1;
      E_Int neighPH = NEIGHBOR1(PHi, _F2Es, nb_pgs, PGi);
      
      if (neighPH == this->_F2E_NONE) continue; // BC
      
      else if (neighPH > 0)                // INNER
      {
        E_Int& typ = _ng.PHs._type[neighPH - 1];
        
        if (typ != UNSET_COL) continue;     // already colored
        typ = val;                  // set it !
        zones4nextpass[_id] = true; // flag this zone for next layer
      }
      else /*if (is_join_PG(PGi))*/         // JOIN
      {
        E_Int racid = _ng.PGs._type[PGi];
                
        auto itrz = _joins.find(racid);
        assert (itrz != _joins.end());
        
        zone_t* rz = itrz->second.first;
        
        if (rz->set_ph_type(-neighPH-1, val))
          zones4nextpass[rz->_id] = true;   
      }
    }
  } 
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::set_pg_colors()
{
  
#ifdef DEBUG_ZONE_T
  {
    ngon_unit PHex;
    Vector_t<E_Int> oIds;
    _ng.PHs.extract_of_type(1, PHex, oIds);
    
    ngon_type ngt(_ng.PGs, PHex);
    
    K_FLD::IntArray cnto;
    ngt.export_to_array(cnto);
    medith::write("layer1.plt", _crd, cnto, "NGON");
  }
  {
    ngon_unit PHex;
    Vector_t<E_Int> oIds;
    _ng.PHs.extract_of_type(2, PHex, oIds);
    
    ngon_type ngt(_ng.PGs, PHex);
    
    K_FLD::IntArray cnto;
    ngt.export_to_array(cnto);
    medith::write("layer2.plt", _crd, cnto, "NGON");
  }
#endif
  E_Int nb_pgs = _ng.PGs.size();
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int eL    = _F2Es[i];
    E_Int eR    = _F2Es[i+nb_pgs];
#ifdef DEBUG_ZONE_T
    assert(eL > 0);
#endif
    E_Int typeL = _ng.PHs._type[eL-1];
    E_Int typeR = (eR != _F2E_NONE) ? _ng.PHs._type[eR-1] : IDX_NONE;
    
    if (typeL == typeR) // PG_INNER_COL, PG_LAY1_IN_COL, PG_LAY2_IN_COL
    {
      if (typeL == PH_INNER_COL)
        _ng.PGs._type[i] = PG_INNER_COL;
      else if (typeL == 1)
        _ng.PGs._type[i] = PG_LAY1_IN_COL;
      else if (typeL == 2)
        _ng.PGs._type[i] = PG_LAY2_IN_COL;
      else 
        _ng.PGs._type[i] = OTHER_LAYS_COL;
    }
    else if (typeR == IDX_NONE) // BC
    {
      if (typeL == PH_INNER_COL)
        _ng.PGs._type[i] = PG_BC;
      else if (typeL == 1)
        _ng.PGs._type[i] = PG_LAY1_BC_COL;
      else if (typeL == 2)
        _ng.PGs._type[i] = PG_LAY2_BC_COL;
      else 
        _ng.PGs._type[i] = OTHER_LAYS_COL;
    }
    else // JOIN or INNER in layers
    {
      E_Int m = K_FUNC::E_min(typeL, typeR);
      if (m == PH_INNER_COL) // JOIN (form inner to 1st layer)
        _ng.PGs._type[i] = PG_JOIN_COL;
      else if (m == 1) // 1st layer
        _ng.PGs._type[i] = PG_LAY1_IN_COL;
      else if (m == 2)
        _ng.PGs._type[i] = PG_LAY2_IN_COL;
      else
        _ng.PGs._type[i] = OTHER_LAYS_COL;
    }
  }

#ifdef DEBUG_ZONE_T
  E_Int pg_inner(0),pg_join(0), lay1in(0), lay1bc(0), bc(0), lay2in(0), lay2bc(0), other(0);
          
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    if ( _ng.PGs._type[i] == PG_INNER_COL)++pg_inner;
    else if ( _ng.PGs._type[i] == PG_JOIN_COL)++pg_join;
    else if ( _ng.PGs._type[i] == PG_BC)++bc;
    else if ( _ng.PGs._type[i] == PG_LAY1_IN_COL)++lay1in;
    else if ( _ng.PGs._type[i] == PG_LAY1_BC_COL)++lay1bc;
    else if ( _ng.PGs._type[i] == PG_LAY2_IN_COL)++lay2in;
    else if ( _ng.PGs._type[i] == PG_LAY2_BC_COL)++lay2bc;
    else if ( _ng.PGs._type[i] == OTHER_LAYS_COL)++other;
  }
  
  std::cout << " PG_INNER_COL : " << pg_inner << std::endl;
  std::cout << " PG_JOIN_COL : " << pg_join << std::endl;
  std::cout << " PG_BC : " << bc << std::endl;
  std::cout << " PG_LAY1_IN_COL : " << lay1in << std::endl;
  std::cout << " PG_LAY1_BC_COL : " << lay1bc << std::endl;
  std::cout << " PG_LAY2_IN_COL : " << lay2in << std::endl;
  std::cout << " PG_LAY2_BC_COL : " << lay2bc << std::endl;
  std::cout << " OTHER_LAYS_COL : " << other << std::endl;
#endif

  __distinguish_bcs();
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::__distinguish_bcs()
{
  for (auto b=_bcs.begin(); b != _bcs.end(); ++b)
  {
    E_Int bcid = b->first;
    Vector_t<E_Int>& ids = b->second;
    size_t nbc = ids.size();
    for (size_t i=0; i < nbc; ++i)
    {
      E_Int PGi = ids[i] - 1;
      if (_ng.PGs._type[PGi] != PG_BC && _ng.PGs._type[PGi] != PG_LAY1_BC_COL && _ng.PGs._type[PGi] != PG_LAY2_BC_COL) continue;
      _ng.PGs._type[PGi] += bcid; // first bc will have 1001, second 1002...
    }
  }
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::add_boundary(E_Int obcid, const E_Int* ids, E_Int n)
{
  E_Int newId = (_bcs.empty() ? 0 : _bcs.rbegin()->first) + 1;
  Vector_t<E_Int>& vids = _bcs[newId];
  vids.insert(vids.begin(), ids, ids+n);
  //std::cout << "adding in map : " << newId <<"/"<<obcid << std::endl;
  _bc_inId2outId[newId] = obcid;
}


/// 
/// \param z
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::copy_boundaries(const zone_t& z)
{  
  for (auto b : z._bcs)
    _bcs[b.first] = b.second;
  
  for (const auto& i : z._joins)
    _joins[i.first] = i.second;
  
  _bc_inId2outId = z._bc_inId2outId;
  _rac_inId2outId = z._rac_inId2outId;
}

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::update_boundaries(const Vector_t<E_Int>& pgnids, E_Int idx_start)
{
  Vector_t<E_Int> vtmp, to_remove, nids;

  if (!_bcs.empty())
  {
    for (auto b=_bcs.begin(); b != _bcs.end(); ++b){
      
      E_Int bcid = b->first;
      size_t nbc = b->second.size();
      
      if (nbc == 0) { 
        to_remove.push_back(bcid);
        continue;
      }
      
      vtmp.clear();
      nids.clear();
      
      vtmp.reserve(nbc);
      nids.resize(nbc, IDX_NONE);
      E_Int count(0);
      for (size_t i = 0; i < nbc; ++i)
      {
        E_Int PGi = b->second[i] - idx_start;
        if (pgnids[PGi] == IDX_NONE) continue;
        vtmp.push_back(pgnids[PGi] + idx_start);
        nids[i] = count++;
      }

      if (!vtmp.empty())
        _bcs[bcid] = vtmp;
      else
        to_remove.push_back(bcid);
    }

    for (size_t i=0; i < to_remove.size(); ++i)
    {
      _bcs.erase(to_remove[i]);
      _bc_inId2outId.erase(to_remove[i]);
    }
  }

  to_remove.clear();
  
  for (auto j : _joins){
      
    E_Int jid = j.first;
    size_t nbj = j.second.second.size();
    
    if (j.second.second.empty()) { 
      to_remove.push_back(jid);
      continue;
    }
    
    vtmp.clear();
    nids.clear();

    vtmp.reserve(nbj);
    nids.resize(nbj, IDX_NONE);
    E_Int count(0);
    for (size_t i = 0; i < nbj; ++i)
    {
      E_Int PGi = j.second.second[i] - idx_start;
      if (pgnids[PGi] == IDX_NONE) continue;
      vtmp.push_back(pgnids[PGi] + idx_start);
      nids[i] = count++;
    }
    
    zone_t* opp_z = _joins[jid].first;
    
    if (!vtmp.empty())
    {
      _joins[jid].second = vtmp;
      // update opposite accordingly
      if (vtmp.size() < nbj) opp_z->compress_join(jid, nids);
    }
    else
    {
      to_remove.push_back(jid);
      opp_z->_joins.erase(jid);
    }
  }
  
  for (size_t i=0; i < to_remove.size(); ++i)
    _joins.erase(to_remove[i]);
}

///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::sort_pointLists()
{
  for (auto& b : _bcs)
  {
     std::vector<E_Int>& ids = b.second;
     std::sort(ids.begin(), ids.end());

  }
  for (auto j : _joins)
  {
    std::vector<E_Int>& ids = j.second.second;
    std::sort(ids.begin(), ids.end());
  }
}


///
template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::append( zone_t& that_z, bool keep_joins)
{
  if (that_z.empty()) return;
  
  // 1. If disconnected (i.e. no join) => just concatenate
  Vector_t<E_Int> jids ;
  
  bool discnt = is_disconnected_from(&that_z, jids);
  
  if (discnt)
  {
    zone_t ztmp(that_z); // fixme : its a clone, the join are still pointing to that_z
    E_Int shft = _crd.cols();
    ngon_type ngtmp;
    ngtmp.PGs = ztmp._ng.PGs;
    ngtmp.PHs = that_z._ng.PHs;
    ngtmp.PGs.shift(shft);
    _crd.pushBack(ztmp._crd); // concatenate points
    
    _ng.append(ngtmp);
    
    _bcs.insert(that_z._bcs.begin(), that_z._bcs.end());
    _joins.insert(that_z._joins.begin(), that_z._joins.end());
    
    size_t sz0 = _F2Es.size() /2;
    Vector_t<E_Int> tmp;
    tmp.insert(tmp.end(), &_F2Es[sz0], &_F2Es[sz0] + sz0);
    _F2Es.resize(sz0);
    
    size_t sz1 = ztmp._F2Es.size() /2;
    Vector_t<E_Int> tmp1;
    tmp1.insert(tmp1.end(), &ztmp._F2Es[sz1], &ztmp._F2Es[sz1] + sz1);
    ztmp._F2Es.resize(sz1);

    _F2Es.insert(_F2Es.end(), ztmp._F2Es.begin(), ztmp._F2Es.end());
    _F2Es.insert(_F2Es.end(), tmp.begin(), tmp.end());
    _F2Es.insert(_F2Es.end(), tmp1.begin(), tmp1.end());
  }
  else // has joins
  {
    assert (!jids.empty());
   
    E_Int shftPG0 = _ng.PGs.size();
    E_Int shftPH0 = _ng.PHs.size();
    
    // start the appending
        
    E_Int shft = _crd.cols();
    ngon_type ngtmp;
    ngtmp.PGs = that_z._ng.PGs;
    ngtmp.PHs = that_z._ng.PHs;
    ngtmp.PGs.shift(shft);
    
    _crd.pushBack(that_z._crd); // concatenate points
    
    // reindex and append pgs not in any join, update this F2E
    
    E_Int that_nb_pgs = ngtmp.PGs.size();
    Vector_t<E_Int> pgnids(that_nb_pgs, IDX_NONE);
    
    for (size_t jj = 0; jj < jids.size(); ++jj)
    {
      auto itri = _joins.find(jids[jj]);
      assert (itri != _joins.end());
      auto & this_ids = itri->second.second;
    
      const Vector_t<E_Int>&  that_ids = that_z.get_join_pg_list(jids[jj]);
      
      E_Int nb_pgsj = this_ids.size();
      assert (nb_pgsj == that_ids.size());

      // replace join ids
      for (E_Int i = 0; i < nb_pgsj; ++i)
      {
        E_Int that_id = that_ids[i] - 1;
        assert (that_id < that_nb_pgs && that_id >=0);
        
        E_Int this_id = this_ids[i] - 1;
        
        pgnids[that_id] = this_id;
        
        E_Int this_right = this_id+shftPG0;
        if (_F2Es[this_right] != _F2E_NONE && _F2Es[this_right] >= 0) this_right = this_id; // take left instead
        assert (_F2Es[this_right] == _F2E_NONE || _F2Es[this_right] < 0);
        
        E_Int that_left = that_id;
        E_Int that_right = that_left+that_nb_pgs;
        if (that_z._F2Es[that_right] != _F2E_NONE && that_z._F2Es[that_right] >= 0) std::swap(that_left, that_right);
        assert (that_z._F2Es[that_right] == _F2E_NONE || that_z._F2Es[that_right] < 0);
        
        _F2Es[this_right] = that_z._F2Es[that_left] + shftPH0; // assign that left to be this right
      }
    }

    // Now any NONE in pgnids is referring to an inner PG : those to be extracted
    {
      ngon_unit pgs;
      K_CONNECT::invalid pred(pgnids); // i.e inner PGs
      Vector_t<E_Int> poids, pnids;
      ngtmp.PGs.extract_by_predicate(pred, pgs, poids, pnids);
      
      _ng.PGs.append(pgs);

      E_Int nb_pgs1 = _ng.PGs.size();
      
      // F2E    
      //reshape accordingly F2E     
//      _F2Es.resize(2*szn, -1);
//      for (E_Int i = 0; i < sz0; ++i)
//        _F2Es[i+szn] = _F2Es[i+sz0];
//      for (E_Int i=sz0; i < szn; ++i)_F2Es[i]=-1;
      {
        Vector_t<E_Int> newF2E(2*nb_pgs1, _F2E_NONE);
        for (E_Int i = 0; i < shftPG0; ++i)
        {
          newF2E[i] = _F2Es[i];
          newF2E[i+nb_pgs1] = _F2Es[i+shftPG0];
        }
        
        _F2Es = newF2E;  
      }

      for (E_Int i = 0; i < that_nb_pgs; ++i)
      {
        if (pgnids[i] != IDX_NONE) continue;
        pgnids[i] = pnids[i] + shftPG0;
      }

#ifdef DEBUG_ZONE_T
      this->draw_boundaries();
#endif
      
      size_t sz1 = that_z._F2Es.size() / 2;
      
      // shift the value considering the on going appending
      for (size_t i=0; i < that_z._F2Es.size(); ++i)
        that_z._F2Es[i] = (that_z._F2Es[i] != that_z._F2E_NONE && that_z._F2Es[i] > 0 ) ? that_z._F2Es[i] + shftPH0 : _F2E_NONE;
      
      // append the value
      
      for (E_Int i = 0; i < pgs.size(); ++i)
      {
        E_Int left = that_z._F2Es[poids[i]];
        E_Int right = that_z._F2Es[poids[i] + sz1];
        if (left == _F2E_NONE)std::swap(left,right);
        _F2Es[i + shftPG0]           = left;
        _F2Es[i + shftPG0 + nb_pgs1] = right;
      }
    }
    
#ifdef DEBUG_ZONE_T
      this->draw_boundaries();
#endif
    
    ngtmp.PHs.change_indices(pgnids); 

    _ng.PHs.append(ngtmp.PHs);

    that_z.update_boundaries(pgnids, 1);
    
    // all zone joined to that_z must point now on "this" now
    for (auto& r2z : that_z._joins)
    {
      zone_t* zi = r2z.second.first;
      zi->change_joins_zone(&that_z, this);
    }
    
    _bcs.insert(that_z._bcs.begin(), that_z._bcs.end());
    _bc_inId2outId.insert(that_z._bc_inId2outId.begin(), that_z._bc_inId2outId.end());
    
    if (keep_joins) // preserve common joins and keep initial ids (from this_z point of view (PointList) 
    {
      for (const auto& it_that : that_z._joins)
      {
        E_Int jid = it_that.first;
        auto it_this = this->_joins.find(jid);
        
        if (it_this == this->_joins.end()) //only append if new join (i.e join with a third zone)
          this->_joins[jid] = it_that.second;
      }
    }
    else // append the new one, delete the join between these 2 zones
    {
      _joins.insert(that_z._joins.begin(), that_z._joins.end());
      _rac_inId2outId.insert(that_z._rac_inId2outId.begin(), that_z._rac_inId2outId.end());
      
      for (size_t jj = 0; jj < jids.size(); ++jj)
      {
        _joins.erase(jids[jj]);
        _rac_inId2outId.erase(jids[jj]);
      }
    }
  }
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::sort_by_type()
{
  Vector_t<E_Int> nids, oids;
  _ng.PGs.sort_by_type(nids, oids);

  // do the change
  {
    size_t nb_pgs = _ng.PGs.size();
    ngon_unit pgs;
    pgs.reserve(nb_pgs);
    pgs._NGON.reserve(_ng.PGs._NGON.size());
    
    Vector_t<E_Int> typ(nb_pgs);

    for (size_t i=0; i < nb_pgs; ++i){
      pgs.add(_ng.PGs.stride(oids[i]), _ng.PGs.get_facets_ptr(oids[i]));
      typ[i] = _ng.PGs._type[oids[i]];
    }

    _ng.PGs = pgs;
    _ng.PGs._type = typ;
  }
  
  update_boundaries(nids, 1/*index_start*/); //will only refresh joins and bc ids (as nids is a permutation, no reduction
  
#ifdef DEBUG_ZONE_T
  {
    std::set<E_Int> colors( _ng.PGs._type.begin(),  _ng.PGs._type.end());
    E_Int count(0);
    for (auto c = colors.begin(); c != colors.end(); ++c)
    {
      std::stringstream o;
      o << "col_" << count++ << ".plt";
      ngon_unit pgs;
      Vector_t<E_Int> oids;
      _ng.PGs.extract_of_type(*c, pgs, oids);
      
      ngon_type ng(pgs, true);
      K_FLD::IntArray cnto;
      ng.export_to_array(cnto);
      medith::write(o.str().c_str(), _crd, cnto, "NGON");
      
    }
  }
#endif

  _ng.PHs.change_indices(nids);
  
  _ng.PHs.sort_by_type(nids, oids);
  
  // do the change
  {
    size_t nb_phs = _ng.PHs.size();
    ngon_unit phs;
    phs.reserve(nb_phs);
    phs._NGON.reserve(_ng.PHs._NGON.size());

    Vector_t<E_Int> typ(nb_phs);
    
    for (size_t i=0; i < nb_phs; ++i){
      phs.add(_ng.PHs.stride(oids[i]), _ng.PHs.get_facets_ptr(oids[i]));
      typ[i] = _ng.PHs._type[oids[i]];
    }

    _ng.PHs = phs;
    _ng.PHs._type = typ;
  }
  
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::insert_ghosts_on_bcs(E_Int type, E_Int nb_layers) //type 1 : single-PG host, 2: degen ghost, 3 : degen ghost with top creation (rather than bot duplication)
{

  // first shift types to make room for BC ghost (just before 2nd layer)
  size_t nb_phs = _ng.PHs.size();
  for (size_t i=0; i < nb_phs; ++i)
  {
    if (_ng.PHs._type[i] > 1 && _ng.PHs._type[i] != IDX_NONE) ++_ng.PHs._type[i];
  }
  // Now add one cell per BC pg
  size_t nb_pgs = _ng.PGs.size();
  if (type == 1)
  {
    // Now add one cell per BC pg
    E_Int pg;

    // RULE
    // for any nb of layers, add PG_BC
    // for nlayers >=2 , add also PG_LAY1_BC_COL
    // for nlayer n>=3 , add also PG_LAY2_BC_COL
    // above not handled

    for (size_t i=0; i < nb_pgs; ++i)
    {
      bool ad = is_a_bc(_ng.PGs._type[i], nb_layers);

      if (ad)
      {
        pg = i+1;
        _ng.PHs.add(1, &pg);
      }
    }
    _ng.PHs.updateFacets();
    nb_phs = _ng.PHs.size();
    _ng.PHs._type.resize(nb_phs, 2); // set GHOST cells color
  }
  else if (type == 2 || type == 3)
  {
    std::vector<E_Int> pglist;
    for (size_t i=0; i < nb_pgs; ++i)
    {
      bool ad = is_a_bc(_ng.PGs._type[i], nb_layers);

      if (ad)
        pglist.push_back(i);
    }
    ngon_type::add_flat_ghosts(_ng, pglist, (type == 3));
    //replace temporarily to put ghost at the right position
    nb_phs = _ng.PHs.size();
    for (size_t i = 0; i < nb_phs; ++i)
      if (_ng.PHs._type[i] == PH_GHOST)_ng.PHs._type[i] = 2; 
  }
  
  
  
  // put the PHs of type 2 at the right place
  Vector_t<E_Int> nids, oids;
  _ng.PHs.sort_by_type(nids, oids);
  
  // do the change
  {
    size_t nb_phs = _ng.PHs.size();
    ngon_unit phs;
    phs.reserve(nb_phs);
    phs._NGON.reserve(_ng.PHs._NGON.size());

    Vector_t<E_Int> typ(nb_phs);
    
    for (size_t i=0; i < nb_phs; ++i){
      phs.add(_ng.PHs.stride(oids[i]), _ng.PHs.get_facets_ptr(oids[i]));
      typ[i] = _ng.PHs._type[oids[i]];
    }

    _ng.PHs = phs;
    _ng.PHs._type = typ;
  }
  
  // Now get back to previous types
  nb_phs = _ng.PHs.size();
  for (size_t i=0; i < nb_phs; ++i)
  {
    if (_ng.PHs._type[i] == 2) _ng.PHs._type[i] = PH_GHOST;
    else if (_ng.PHs._type[i] > 1 && _ng.PHs._type[i] != IDX_NONE) --_ng.PHs._type[i];
  }
}

template <typename crd_t, typename ngo_t>
void zone_t<crd_t, ngo_t>::color_ranges
(std::vector<E_Int>& PGcolors, std::vector<E_Int>& PHcolors)
{
  PGcolors.clear(); PHcolors.clear();
  E_Int pgColCur=_ng.PGs._type[0];
  
  PGcolors.push_back(0);
  E_Int i=0;

  while (i<_ng.PGs.size())
  {
    while (i<_ng.PGs.size() && _ng.PGs._type[i] ==  pgColCur)++i;
    pgColCur = _ng.PGs._type[i];
    PGcolors.push_back(i);
    ++i;
  }
  
  E_Int phColCur=_ng.PHs._type[0];
  PHcolors.push_back(0);
  i=0;

  while (i<_ng.PHs.size())
  {
    while (i<_ng.PHs.size() && _ng.PHs._type[i] ==  phColCur)++i;
    phColCur = _ng.PHs._type[i];
    PHcolors.push_back(i);
    ++i;
  }
}

} // NUGA

#endif /* ZONE_T_HXX */

