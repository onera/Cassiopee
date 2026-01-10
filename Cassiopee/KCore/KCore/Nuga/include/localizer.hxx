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

#ifndef LOCALIZER_HXX
#define LOCALIZER_HXX

#ifdef DEBUG_COLLIDER
#include "Nuga/include/chrono.h"
#endif

namespace NUGA
{

template <typename Tree_t>
class localizer
{
  public:
    
    using box_t = typename Tree_t::box_type;
  
  public:

  localizer() = delete;
  localizer(const localizer&) = delete; // to avoid to copy a loc owned by mesh_t that will destroy it
  
  localizer(Tree_t& tree, E_Float tolerance):_tree(&tree), _owner(false), _tolerance(tolerance){_tree->_tolerance = _tolerance;}

  // tree owner version
  localizer(Tree_t* tree, E_Float tolerance):_tree(tree), _owner(true), _tolerance(tolerance){_tree->_tolerance = _tolerance;}
  
  template <typename acrd_t, typename acnt_t>
  localizer(const acrd_t& coords, const acnt_t& connect, E_Float tolerance):_owner(true), _tolerance(tolerance){__create_tree(coords, connect);}
  
  ~localizer() { if (_owner) __destroy_tree(); }
  
  template<typename ELT, typename crd_t>
  void get_candidates(const ELT& e, const crd_t& crde, std::vector<E_Int>& candidates, int idx_start = 0, E_Float RTOL = 0.) const ;

  void get_candidates(const box_t& bb, std::vector<E_Int>& candidates, int idx_start = 0, E_Float RTOL = 0.) const ;
  
  const Tree_t* get_tree() const { return _tree;}
  
  private:
    template <typename acrd_t, typename acnt_t>
    void __create_tree(const acrd_t& crd, const acnt_t& cnt);

    void __destroy_tree();
  
  private:
    Tree_t* _tree;
    std::vector<box_t> _boxes;
    bool _owner;
    E_Float _tolerance;
};

///
template <typename Tree_t>
template<typename ELT, typename crd_t>
void localizer<Tree_t>::get_candidates(const ELT& e, const crd_t& crde, std::vector<E_Int>& candidates, int idx_start, E_Float RTOL) const {
  
  box_t bb;
  e.bbox(crde, bb);

  if (RTOL > 0.) bb.enlarge(RTOL);

  bb.minB[0] -= ZERO_M;
  bb.minB[1] -= ZERO_M;
  bb.minB[2] -= ZERO_M;

  bb.maxB[0] += ZERO_M;
  bb.maxB[1] += ZERO_M;
  bb.maxB[2] += ZERO_M;
  
  get_candidates(bb, candidates, idx_start);
}

///
template <typename Tree_t>
void localizer<Tree_t>::get_candidates(const box_t& bb, std::vector<E_Int>& candidates, int idx_start, E_Float RTOL) const {

  candidates.clear();
  _tree->getOverlappingBoxes(bb.minB, bb.maxB, candidates);

  if (idx_start != 0)
    K_CONNECT::IdTool::shift(candidates, idx_start);
}



///
template <typename Tree_t>
template <typename acrd_t, typename acnt_t>
void localizer<Tree_t>::__create_tree(const acrd_t& coords, const acnt_t& connect)
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
template <typename Tree_t>
void localizer<Tree_t>::__destroy_tree()
{
  _boxes.clear();
  delete _tree;
  _tree=nullptr;
}

}

#endif /* LOCALIZER_HXX */

