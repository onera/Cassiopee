/*
 
 
 
              NUGA 
 
 
 
 */
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef LOCALIZER_HXX
#define LOCALIZER_HXX

namespace NUGA
{

template <typename Tree_t, typename acrd_t, typename acnt_t>
class localizer
{
  public:
    
    using box_t = typename Tree_t::box_type;
  
  public:
  
  localizer(const Tree_t& tree, E_Float tolerance):_tree(&tree), _owner(false), _tolerance(tolerance){}
  
  localizer(const acrd_t& coords, const acnt_t& connect, E_Float tolerance):_owner(true), _tolerance(tolerance){__create_tree(coords, connect);}
  
  ~localizer(){ 
    if (_owner) __destroy_tree();
  }
  
  template<typename ELT>
  void get_candidates(const ELT& e, const acrd_t crde, std::vector<E_Int>& candidates) const ;
  
  private:
    void __create_tree(const acrd_t& crd, const acnt_t& cnt);
    void __destroy_tree();
  
  private:
    const Tree_t* _tree;
    std::vector<box_t> _boxes;
    bool _owner;
    E_Float _tolerance;
};

///
template <typename Tree_t, typename acrd_t, typename acnt_t>
template<typename ELT>
void localizer<Tree_t, acrd_t, acnt_t>::get_candidates(const ELT& e, const acrd_t crde, std::vector<E_Int>& candidates) const { 
  
  box_t bb;
  e.bbox(crde, bb);
  
  candidates.clear();
  _tree->getOverlappingBoxes(bb.minB, bb.maxB, candidates);
}



///
template <typename Tree_t, typename acrd_t, typename acnt_t>
void localizer<Tree_t, acrd_t, acnt_t>::__create_tree(const acrd_t& coords, const acnt_t& connect)
{

#ifdef DEBUG_COLLIDER
  chrono c;
  c.start();
#endif
  
  K_FLD::IntArray e;
  E_Int s, sz(connect.size());
 
#ifdef DEBUG_COLLIDER 
  std::cout << "create " << sz << " boxes.." <<  std::endl;
#endif

  for (E_Int i = 0; i < sz; ++i)
  {
    s = connect.stride(i);
    e.reserve(1, s);
    connect.getEntry(i, e.begin());
    
    box_t* bb = new box_t(coords, e.begin(), s);
    _boxes.push_back(bb);
  }
  
#ifdef DEBUG_COLLIDER
  std::cout << "boxes done : " << c.elapsed() << std::endl;
  
  std::cout << "create tree.." <<  std::endl;
  c.start();
#endif
  
  _tree = new Tree_t(_boxes, _tolerance);
  
#ifdef DEBUG_COLLIDER
  std::cout << "tree done : " << c.elapsed() << std::endl;
#endif
}

///
template <typename Tree_t, typename acrd_t, typename acnt_t>
void localizer<Tree_t, acrd_t, acnt_t>::__destroy_tree()
{
  _boxes.clear();
  delete _tree;
  _tree=nullptr;
}

}

#endif /* LOCALIZER_HXX */

