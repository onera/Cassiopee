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
//Authors : Sam Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)
 
#ifndef NUGA_TREE_HXX
#define NUGA_TREE_HXX

#include "Nuga/include/Hexahedron.h"
#include "Nuga/include/Polyhedron.h"

using ngon_type = ngon_t<K_FLD::IntArray>;

namespace NUGA
{

//
template <typename array>
class array_trait;

// children_array : ngon_unit for PGs/PHs/HybridBasic, IntArray for MonoBasic (T3, Q4, HX6...) 
template <typename children_array>
class tree
{  
  private:
    ngon_unit *       _entities; // Facets(PG) or Cells (PHs)
    children_array    _children; // 
    Vector_t<E_Int>   _parent; //sized as entities
    Vector_t<E_Int>   _indir; //sized as entities
    Vector_t<E_Int>   _level; //sized as entities
    Vector_t<bool>    _enabled; //sized as entities
    
  public:
    explicit tree(ngon_unit & entities):_entities(&entities){ resize_hierarchy(entities.size());}
    
    void set_entities(ngon_unit& ngu) { _entities = &ngu ;}//relocate for hook

    const Vector_t<E_Int>& level() const {return _level;}
        
    // to make sizes consistent : need to be called when refining the mesh
    void resize_hierarchy(size_t nb_ent)
    {
      _parent.resize(nb_ent, IDX_NONE);
      _indir.resize(nb_ent, IDX_NONE);
      _level.resize(nb_ent, 0);
      _enabled.resize(nb_ent, true);
    }

    void reset_enabled(E_Int sz, bool val) // for PG enabling in hmesh
    {
      _enabled.clear();
      _enabled.resize(sz, val);
    }

    const Vector_t<bool>& get_enabled() { return _enabled;}
    
    // expand a tree with fixed-stride-children
    void resize(const Vector_t<E_Int>& ids, E_Int stride)
    {
      // first available local id : one passed-the-end before appending : important to get it before resizing _children
      E_Int locid = array_trait<children_array>::size(_children); // nb_child
      
      // get the total nb of children to add
      E_Int nb_new_children = stride * ids.size();
      
      // expand the children array
      array_trait<children_array>::resize_children(_children, stride, nb_new_children);

      // expand remaining attributes
      E_Int current_sz = _parent.size();
      resize_hierarchy(current_sz + nb_new_children);

      // set the local id of each entity promoted for refinement.
      E_Int n = ids.size();
      for (E_Int i=0; i < n; ++i)
        _indir[ids[i]] = locid++;
    }

    // expand a tree with variable-stride-children
    // pregnant[i] is the nb of children for ids[i]
    void resize(const Vector_t<E_Int>& ids, const Vector_t<E_Int>& pregnant)
    {
      if (pregnant.empty()) return;

      assert(ids.size() == pregnant.size());

      //first available local id : one passed-the-end before appending : important to get it before resizing _children
      E_Int locid = array_trait<children_array>::size(_children); // nb_child
      
      // get the total nb of new entities after refining
      E_Int nb_new_children = K_CONNECT::IdTool::sum(pregnant);

      // expand the children array
      _children.expand_variable_stride(pregnant.size(), &pregnant[0]);

      // expand remaining attributes
      E_Int current_sz = _parent.size();
      resize_hierarchy(current_sz + nb_new_children); // _espace parent 

      // set the local id of each entity promoted for refinement.
      E_Int n = ids.size();
      for (E_Int i=0; i < n; ++i)
      {
        ASSERT_IN_VECRANGE(_indir, ids[i])
        _indir[ids[i]] = locid++;
      }
    }

    inline E_Int size() const
    {
      size_t sz = _indir.size();

      if (sz != _parent.size()) return -1;
      if (sz != _level.size())  return -1;
      if (sz != _indir.size()) return -1;
      if (sz != _enabled.size()) return -1;

      return E_Int(sz);
    }
    
    inline const E_Int& get_level(E_Int i /*zero based*/) const { ASSERT_IN_VECRANGE(_level, i);  return _level[i];}
    
    inline void set_level(E_Int i /*zero based*/, E_Int level) { ASSERT_IN_VECRANGE(_level, i);  _level[i] = level;}
    
    inline E_Int parent(E_Int i /*zero based*/) const { ASSERT_IN_VECRANGE(_parent, i);  return _parent[i];}
    
    void get_oids(std::vector<E_Int>& oids) const ;

    E_Int get_root_parent(E_Int i) const;

    //WRONG (enablin logic) and unused
    // void add_children(E_Int i/*zero based*/, const E_Int* children, E_Int n){
     
    //   assert(i < _entities->size());
    //   _indir[i] = array_trait<children_array>::size(_children);// size of _children

    //   array_trait<children_array>::add_children(_children, children, n);
      
    //   for (size_t c=0; c<n; ++c) _parent[children[c]] = i;
      
    //   _level.resize(_level.size()+n, _level[i]+1);
    //   // enable the children, disable himself
    //   _enabled.resize(_level.size()+n, true); // WRONG : use method enable
    //   _enabled[i] = false;                    // WRONG
    // }
    
    // set n random children ids
    void set_children(E_Int i/*zero based*/, const E_Int* childr, E_Int n){
     
      E_Int* there = children(i);
      assert(there != NULL);
      
      std::copy(childr, childr+n, there);
      
      for (E_Int c=0; c<n; ++c) _parent[childr[c]] = i;
      
    }

    // set n contiguous children id
    void set_children(E_Int i/*zero based*/, E_Int firstChild, E_Int n) {

      E_Int* there = children(i);
      assert(there != NULL);

      for (E_Int c = 0; c < n; ++c)
      {
        *(there++) = firstChild + c;
        _parent[firstChild + c] = i;
      }
    }

    //
    E_Int nb_children(E_Int i /*zero based*/) const {
      ASSERT_IN_VECRANGE(_indir, i)
      if (_indir[i] == IDX_NONE) return 0;
      return array_trait<children_array>::nb_children(_children, _indir[i]);}
    
    //
    const E_Int* children(E_Int i /*zero based*/) const {
      ASSERT_IN_VECRANGE(_indir, i)
      if (_indir[i] == IDX_NONE) return nullptr;
      return array_trait<children_array>::children(_children, _indir[i]);
    }
    
    E_Int* children(E_Int i /*zero based*/) {
      ASSERT_IN_VECRANGE(_indir, i)
      if (_indir[i] == IDX_NONE) return nullptr;
      return array_trait<children_array>::children(_children, _indir[i]);
    }
    
    void enable(E_Int i /*zero based*/, bool act_on_genealogy = true)
    {
       ASSERT_IN_VECRANGE(_enabled, i)

       _enabled[i] = true;

       if (!act_on_genealogy) return;
       
       // disable its parent
       E_Int pid = parent(i);
       if ( pid !=  IDX_NONE )_enabled[pid] = false;

       // disable its children
       E_Int nbc = nb_children(i);
       const E_Int* childr = children(i);
       for (E_Int n = 0; n < nbc; ++n)
         _enabled[*(childr + n)] = false;
    }
      
    inline bool is_enabled(E_Int i /*zero based*/) const { ASSERT_IN_VECRANGE(_enabled, i);  return _enabled[i]; }

    E_Int get_enabled_parent(E_Int i, E_Int& parent) const;

    E_Int get_enabled_descendants(E_Int i, std::vector<E_Int>& ids) const;

    void get_enabled_relatives(E_Int i, std::vector<E_Int>& ids) const;
    
    void get_enabling_histo
    (
      const std::vector<E_Int> & h_ids,
      std::map<E_Int, std::vector<E_Int>>& hph_oids,
      std::map<E_Int, std::vector<E_Int>>& hph_nids
    );

    void __get_enabled_descendants(E_Int i, std::vector<E_Int>& ids) const;

};

/// 
template <typename children_array>
void tree<children_array>::get_oids(std::vector<E_Int>& oids) const
{
  size_t nb_ents(_parent.size());
  
  oids.clear();
  K_CONNECT::IdTool::init_inc(oids, nb_ents);
  
  for (size_t i=0; i < nb_ents; ++i)
  {
    E_Int pid = _parent[i];
    while (pid != IDX_NONE) //get back to root _parent
    {
      oids[i] = pid;
      ASSERT_IN_VECRANGE(_parent, pid);
      pid = _parent[pid];
    };
  }
}

/// 
template <typename children_array>
E_Int tree<children_array>::get_root_parent(E_Int i) const
{
  E_Int pid = _parent[i];
  while (pid != IDX_NONE) //get back to root _parent
  {
    ASSERT_IN_VECRANGE(_parent, pid);

    E_Int npid = _parent[pid];
    if (npid == IDX_NONE) break;
    pid = npid;
  };
  return pid;
}

///
template <typename children_array>
E_Int tree<children_array>::get_enabled_parent(E_Int i, E_Int& pid) const
{
  pid = _parent[i];

  while (pid != IDX_NONE && !_enabled[pid])
  {
    ASSERT_IN_VECRANGE(_parent, pid);
    pid = _parent[pid];
  }

  return 0;
}

///
template <typename children_array>
E_Int tree<children_array>::get_enabled_descendants(E_Int i, std::vector<E_Int>& ids) const
{
  ids.clear();
  if (_enabled[i]) return 0;

  E_Int nbc = nb_children(i);
  if (nbc == 0) return 0;
  const E_Int* childr = children(i);

  for (E_Int c=0; c <nbc; ++c)
    __get_enabled_descendants(childr[c], ids);

  return 0;
}

/// For a given node in thhe tree, retrieves enabled ids (either itself, a parent or all enabled descendants)
template <typename children_array>
void tree<children_array>::get_enabled_relatives(E_Int i, std::vector<E_Int>& ids) const
{
  // rules : itself enabled, just return (leaving ids empty)
  //       : returns enabled parent in ids. WARNNG : assume hmesh cannot have single child descendant
  //       
  ASSERT_IN_VECRANGE(_enabled, i);
  ids.clear();
  if (_enabled[i]) return;  // still enabled

  get_enabled_descendants(i, ids);
  if (!ids.empty()) return; // descendants are enabled

  E_Int pid{ IDX_NONE };
  get_enabled_parent(i, pid);
  assert(pid != IDX_NONE);

  ids.push_back(pid);
}

///
template <typename children_array>
void tree<children_array>::get_enabling_histo
(
  const std::vector<E_Int> & h_ids,
  std::map<E_Int, std::vector<E_Int>>& hph_oids,
  std::map<E_Int, std::vector<E_Int>>& hph_nids
 )
{
  //
  E_Int nb_phs = h_ids.size();

  hph_oids.clear();
  hph_nids.clear();
  
  Vector_t<E_Int> ids;
  //
  for (E_Int i = 0; i < nb_phs; ++i) // previous enabled hmesh ids
  {
    E_Int hphid = h_ids[i];

    get_enabled_relatives(hphid, ids);

    if (ids.empty()) // hphid is still enabled
    {
      hph_nids[hphid].push_back(hphid);
      hph_oids[hphid].push_back(hphid);
    }
    else if (ids.size() == 1) // parent => agglom : store as neg val to distinguish from intact case
    {
      hph_oids[ids[0]].push_back(hphid); // will contain all the bits that have been agglom into ids[0]
      hph_nids[hphid].push_back(-(ids[0] + 1)); // parent origin store as neg
    }
    else
    {
      for (size_t c = 0; c < ids.size(); ++c) // if one entry => parent id, otherwise descendants ids
      {
        hph_oids[ids[c]].push_back(-(hphid + 1)); // parent origin store as neg
        hph_nids[hphid].push_back(ids[c]);
      }
    }
  }
}

template <typename children_array>
void tree<children_array>::__get_enabled_descendants(E_Int i, std::vector<E_Int>& ids) const
{
  if (_enabled[i])
  {
    ids.push_back(i);
    return;
  }

  E_Int nbc = nb_children(i);
  const E_Int* childr = children(i);

  for (E_Int c=0; c <nbc; ++c)
    __get_enabled_descendants(childr[c], ids);
}

//////////////  ARRAY TRAITS : how to get size, expand an array whether it is fixed stride (IntArray) or not (ngon_unit) 

///
template<>
class array_trait<ngon_unit>
{
  public:
    
  static E_Int size(const ngon_unit& arr) { return arr.size();}
  
  static void add_children(ngon_unit& arr, const E_Int* children, E_Int n){
      arr.add(n, children);//alexis : set _type for children
      arr.updateFacets();
  }

  static E_Int nb_children(const ngon_unit& arr, E_Int loci){  
      return arr.stride(loci);
  }
    
  static const E_Int* children(const ngon_unit& arr, E_Int loci) {
      return arr.get_facets_ptr(loci);
  }
  
  static E_Int* children(ngon_unit& arr, E_Int loci) {
      return arr.get_facets_ptr(loci);
  }
 
  //
  static void resize_children(ngon_unit& arr, E_Int stride, E_Int nb_new_children) //fixed stride
  {
    arr.expand_n_fixed_stride(nb_new_children, stride);
  }

};

// FIXED STRIDE CHILDREN ARRAY : BASIC ELEMENTS (2D and 3D)
template<>
class array_trait<K_FLD::IntArray>
{
  public:
    
  static E_Int size(const K_FLD::IntArray& arr) { return arr.cols();}

  static void add_children(K_FLD::IntArray& arr, const E_Int* children, E_Int n){
      arr.pushBack(children, children+n);
  }
    
  static E_Int nb_children(const K_FLD::IntArray& arr, E_Int loci){  
      return arr.rows();
  }
    
  static const E_Int* children(const K_FLD::IntArray& arr, E_Int loci) {
      return arr.col(loci);
  }
  
  static E_Int* children(K_FLD::IntArray& arr, E_Int loci) {
      return arr.col(loci);
  }
  
  static void resize_children(K_FLD::IntArray& arr, E_Int stride, E_Int nb_children){
    arr.resize(stride, arr.cols() + nb_children/stride, IDX_NONE);
  }

};


}


#endif
