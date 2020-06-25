/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)
 
#ifndef NUGA_TREE_HXX
#define NUGA_TREE_HXX

#include "MeshElement/Hexahedron.h"
#include "MeshElement/Polyhedron.h"

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
    explicit tree(ngon_unit & entities, E_Int nbc/*dummy*/):_entities(&entities){ resize_hierarchy(entities.size());}
    
    void set_entities(ngon_unit& ngu) { _entities = &ngu ;}//relocate for hook

    const Vector_t<E_Int>& level() const {return _level;}
        
    // to make sizes consistent : need to be called when refining the mesh
    void resize_hierarchy(size_t nb_ent)
    {
      _parent.resize(nb_ent, E_IDX_NONE);
      _indir.resize(nb_ent, E_IDX_NONE);
      _level.resize(nb_ent, 0);
      _enabled.resize(nb_ent, true);
    }
    
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
        _indir[ids[i]] = locid++;
    }

    inline E_Int size() const
    {
      E_Int sz = _indir.size();

      if (sz != _parent.size()) return -1;
      if (sz != _level.size())  return -1;
      if (sz != _indir.size()) return -1;
      if (sz != _enabled.size()) return -1;

      return sz;
    }
    
    inline const E_Int& get_level(E_Int i /*zero based*/) const {return _level[i];}
    
    inline void set_level(E_Int i /*zero based*/, E_Int level) {_level[i] = level;}
    
    inline E_Int parent(E_Int i /*zero based*/) const { return _parent[i];}
    
    void get_oids(std::vector<E_Int>& oids) const ; //WARNING : NOT VALID AFTER CONFOMIZE
    
    
    //
    void add_children(E_Int i/*zero based*/, const E_Int* children, E_Int n){
     
      assert(i < _entities->size());
      _indir[i] = array_trait<children_array>::size(_children);// size of _children

      array_trait<children_array>::add_children(_children, children, n);
      
      for (size_t c=0; c<n; ++c) _parent[children[c]] = i;
      
      _level.resize(_level.size()+n, _level[i]+1);
      // enable the children, disable himself
      _enabled.resize(_level.size()+n, true);
      _enabled[i] = false;
    }
    
    // set n random children ids
    void set_children(E_Int i/*zero based*/, const E_Int* childr, E_Int n){
     
      E_Int* there = children(i);
      assert(there != NULL);
      
      std::copy(childr, childr+n, there);
      
      for (size_t c=0; c<n; ++c) _parent[childr[c]] = i;
      
    }

    // set n contiguous children id
    void set_children(E_Int i/*zero based*/, E_Int firstChild, E_Int n) {

      E_Int* there = children(i);
      assert(there != NULL);

      for (size_t c = 0; c < n; ++c)
      {
        *(there++) = firstChild + c;
        _parent[firstChild + c] = i;
      }
    }

    //
    E_Int nb_children(E_Int i /*zero based*/) const {
      if (_indir[i] == E_IDX_NONE) return 0;
      return array_trait<children_array>::nb_children(_children, _indir[i]);}
    
    //
    const E_Int* children(E_Int i /*zero based*/) const {
      if (_indir[i] == E_IDX_NONE) return nullptr;
      return array_trait<children_array>::children(_children, _indir[i]);
    }
    
    E_Int* children(E_Int i /*zero based*/) {
      if (_indir[i] == E_IDX_NONE) return nullptr;
      return array_trait<children_array>::children(_children, _indir[i]);
    }
    
    void enable(E_Int i /*zero based*/)
    {
       _enabled[i] = true;
       
       // disable its parent
       if ( parent(i) !=  E_IDX_NONE )_enabled[parent(i)] = false;

       // disable its children
       E_Int nbc = nb_children(i);
       const E_Int* childr = children(i);
       for (E_Int n = 0; n < nbc; ++n)
         _enabled[*(childr + n)] = false;
    }
      
    inline bool is_enabled(E_Int i /*zero based*/) const { return _enabled[i];}
};

/// WARNING : true while the tree is alive (NOT AFTER CONFORMIZE))
template <typename children_array>
void tree<children_array>::get_oids(std::vector<E_Int>& oids) const
{
  E_Int nb_ents(_parent.size());
  
  oids.clear();
  K_CONNECT::IdTool::init_inc(oids, nb_ents);
  
  for (size_t i=0; i < nb_ents; ++i)
  {
    E_Int pid = _parent[i];
    while (pid != E_IDX_NONE) //get back to root _parent
    {
      oids[i] = pid;
      pid = _parent[pid];
    };
  }
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
    arr.resize(stride, arr.cols() + nb_children/stride, E_IDX_NONE);
  }
  
};


}


#endif
