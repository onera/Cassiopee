/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_TREE_HXX
#define NUGA_TREE_HXX

#include "MeshElement/Hexahedron.h"
#include "MeshElement/Polyhedron.h"

#define Vector_t std::vector
using ngon_type = ngon_t<K_FLD::IntArray>;

namespace NUGA
{
  
// 
template <typename ELT>
class tree;
//
template <typename array>
class children_trait;

// children_array : ngon_unit for PGs/PHs/HybridBasic, IntArray for MonoBasic (T3, Q4, HX6...) 
template <typename children_array>
class ent_tree
{
  template <typename ELT> friend struct tree;
  
  private:
    ngon_unit &       _entities;
    children_array    _children;
    E_Int             _nb_children;
    Vector_t<E_Int>   _parent; //sized as entities
    Vector_t<E_Int>   _indir; //sized as entities
    Vector_t<E_Int>   _level; //sized as entities
    Vector_t<bool>    _enabled; //sized as entities
    
  public:
    explicit ent_tree(ngon_unit & entities, E_Int nbc):_entities(entities), _nb_children(nbc){ resize_hierarchy(entities.size());}
    
    const Vector_t<E_Int>& level() const {return _level;}
        
    // to make sizes consistent : need to be called when refining the mesh
    void resize_hierarchy(size_t nb_ent)
    {
      _parent.resize(nb_ent, E_IDX_NONE);
      _indir.resize(nb_ent, E_IDX_NONE);
      _level.resize(nb_ent, 0);
      _enabled.resize(nb_ent, true);
    }
    
    inline const E_Int& get_level(E_Int i /*zero based*/) const {return _level[i];}
    
    inline void set_level(E_Int i /*zero based*/, E_Int level) {_level[i] = level;}
    
    inline E_Int get_parent_size() {return (E_Int)_parent.size();}
    
    inline E_Int parent(E_Int i /*zero based*/){ return _parent[i];}
    
    //
    void add_children(E_Int i/*zero based*/, const E_Int* children, E_Int n){
     
      assert(i < _entities.size());
      _indir[i] = children_trait<children_array>::size(_children);// size of _children

      children_trait<children_array>::add_children(_children, children, n);
      
      for (size_t c=0; c<n; ++c) _parent[children[c]] = i;
      
      _level.resize(_level.size()+n, _level[i]+1);
      // enable the children, disable himself
      _enabled.resize(_level.size()+n, true);
      _enabled[i] = false;
    }
    
    void set_children(E_Int i/*zero based*/, const E_Int* childr, E_Int n){
     
      E_Int* there = children(i);
      assert(there != NULL);
      
      std::copy(childr, childr+n, there);
      
      for (size_t c=0; c<n; ++c) _parent[childr[c]] = i;
      
    }

    //
    E_Int nb_children(E_Int i /*zero based*/){
      if (_indir[i] == E_IDX_NONE) return 0;
      return children_trait<children_array>::nb_children(_children, _indir[i]);}
    
    //
    const E_Int* children(E_Int i /*zero based*/) const {
      if (_indir[i] == E_IDX_NONE) return nullptr;
      return children_trait<children_array>::children(_children, _indir[i]);
    }
    
    E_Int* children(E_Int i /*zero based*/) {
      if (_indir[i] == E_IDX_NONE) return nullptr;
      return children_trait<children_array>::children(_children, _indir[i]);
    }
    
    void enable(E_Int i /*zero based*/)
    {
       agglomerate(i);
       
       // disable its parent
       _enabled[parent(i)] = false;
    }
    
    inline void agglomerate(E_Int i /*zero based*/)
    {
       _enabled[i] = true;
       
       // disable its children
       E_Int nbc = nb_children(i);
       const E_Int* childr = children(i);
       for (size_t n = 0; n < nbc; ++n) 
         _enabled[*(childr+n)] = false;
    }
    
    inline bool is_enabled(E_Int i /*zero based*/){ return _enabled[i];}
    
    void disable_one_elt(E_Int i /*zero based*/)
    {
        _enabled[i] = false;
    }
    
    void enable_one_elt(E_Int i /*zero based*/)
    {
        _enabled[i] = true;
    }
};



template<>
struct tree<K_MESH::Hexahedron>
{
  ent_tree<K_FLD::IntArray> _PGtree, _PHtree;
  
  public:
    tree(ngon_type& ng): _PGtree(ent_tree<K_FLD::IntArray>(ng.PGs, 4)), _PHtree(ent_tree<K_FLD::IntArray>(ng.PHs, 8)){}
    
    void resize(ent_tree<K_FLD::IntArray>& tree, const Vector_t<E_Int>& PGref)
    {
      E_Int n = PGref.size();

      // expand children
      E_Int NB_CHILDREN = tree._nb_children;
      E_Int children_sz = tree._children.cols();

      tree._children.resize(NB_CHILDREN, children_sz + n, E_IDX_NONE);

      // expand hierarchy data
      E_Int current_sz = tree._parent.size();
      tree.resize_hierarchy(current_sz + n*NB_CHILDREN);

      // set the starting-child-index for each input PG to refine
      E_Int pos = children_sz;
      for (E_Int i=0; i < n; ++i)
        tree._indir[PGref[i]] = pos++;
    }

};

template<>
struct tree<K_MESH::Polyhedron<UNKNOWN> >
{
    ent_tree<ngon_unit> _PGtree, _PHtree;

  public:
    tree(ngon_type& ng): _PGtree(ent_tree<ngon_unit>(ng.PGs, -1)), _PHtree(ent_tree<ngon_unit>(ng.PHs, -1)){}
  
  
  
};



template<>
class children_trait<ngon_unit>
{
  public:
    
  static E_Int size(const ngon_unit& arr) { return arr.size();}
  
  static void add_children(ngon_unit& arr, const E_Int* children, E_Int n){
      arr.add(n, children);
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
};

template<>
class children_trait<K_FLD::IntArray>
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
};


}


#endif
