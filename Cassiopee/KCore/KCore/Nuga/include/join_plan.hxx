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

#ifndef NUGA_JOIN_PLAN_HXX
#define NUGA_JOIN_PLAN_HXX

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_unit.h"
#include "Nuga/include/IdTool.h"
#include "Nuga/include/tree.hxx"

#define NO_CHILDREN -1
#define NO_GRAND_CHILDREN -2

namespace NUGA
{

//
template <typename arr_t>
class join_plan
{
public:

  inline static E_Int extract_compact_enabled_tree(const tree<arr_t>& PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, K_FLD::IntArray& plan);
  
  inline static void extract_sub_plans(const tree<arr_t>& PGtree, E_Int PGi, const K_FLD::IntArray& plan, std::map<E_Int, K_FLD::IntArray>& odata);

  inline static bool one_child_requires(const K_FLD::IntArray& plan);

private:
  ///
  inline static void __extract_subplan(E_Int i, const  K_FLD::IntArray& plan, K_FLD::IntArray& sub_plan);
  inline static E_Int __extract_compact_enabled_tree(const tree<arr_t> & PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, K_FLD::IntArray& plan);
};

///
template <typename arr_t> inline
E_Int join_plan<arr_t>::__extract_compact_enabled_tree
(const tree <arr_t> & PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, K_FLD::IntArray& plan)
{
  //std::cout << "__extract_compact_tree : begin" << std::endl;

  if (PGtree.is_enabled(PGi)) return NO_CHILDREN; //stop at enabled entities

  E_Int nbc = PGtree.nb_children(PGi);
  //std::cout << "nbc : " << nbc << std::endl;

  const E_Int* pchild = PGtree.children(PGi);
  //std::cout << "pchild : " << pchild << std::endl;

  STACK_ARRAY(E_Int, nbc, children);
  for (E_Int i = 0; i < nbc; ++i) children[i] = pchild[i];

  //std::cout << "before F : children.get/reverse/i0" << children.get() << "/" << reverse << "/" << i0 << std::endl;
  // to put in receiver ref frame
  F(children.get(), nbc, reverse, i0);
  //std::cout << "after F" << std::endl;

  //to compact : dont create empty columns for all-leaves children
  E_Int nb_tot_enabled_grand_children{ 0 };
  for (E_Int i = 0; i < nbc; ++i)
  {
    // since PGi is disabled, if one of it children is disabled, it means that it has grand children
    
    //std::cout << "children[i] : " << children[i] << std::endl;
    nb_tot_enabled_grand_children += E_Int(!PGtree.is_enabled(children[i]));
  }

  E_Int cols = plan.cols();

  if (nb_tot_enabled_grand_children == 0)
  {
    plan.pushBack(IDX_NONE, 4); // for a directional face, the last two elements in the column will be empty
    for (E_Int uu = 0; uu < nbc; ++uu) plan(uu, cols) = NO_CHILDREN; //set for real children
    
    return  cols;
  }
 
  //std::cout << "nb_tot_enabled_grand_children : " << nb_tot_enabled_grand_children << std::endl;

  plan.pushBack(IDX_NONE, 4); // for a directional face, the last two elements in the column will be empty
  //std::cout << plan << std::endl;

  for (E_Int i = 0; i < nbc; ++i)
  {
    E_Int val = __extract_compact_enabled_tree(PGtree, children[i], F, reverse, i0, plan);
    //std::cout << "plan sz / cols /val : " << plan.cols() << "/" << cols<< "/" << val << std::endl;
    plan(i, cols) = val;
  }

  //std::cout << "__extract_compact_tree : end" << std::endl;

  return cols;
}

///
template <typename arr_t> inline
void join_plan<arr_t>::__extract_subplan
(E_Int i, const K_FLD::IntArray& plan, K_FLD::IntArray& sub_plan)
{
  //
  sub_plan.clear();

  if (i == NO_CHILDREN) return;

  if (i == NO_GRAND_CHILDREN)
  {
    sub_plan.resize(plan.rows(), 1, NO_CHILDREN);
    return;
  }

  if (i == IDX_NONE) return; //less childen than 4 plan rows

  std::vector<E_Int> nodepool;
  nodepool.push_back(i);

  std::vector<bool> keep(plan.cols(), false);
  keep[i] = true;

  while (!nodepool.empty())
  {
    E_Int c = nodepool.back();
    nodepool.pop_back();

    for (E_Int j = 0; j < plan.rows(); ++j)
    {
      E_Int v = plan(j, c);
      if (v <= 0 || v == IDX_NONE) continue;
      nodepool.push_back(v);
      keep[v] = true;
    }
  }
  //std::cout << plan << std::endl;
  sub_plan = plan;
  std::vector<E_Int> nids;
  K_CONNECT::keep<bool> pred_keep(keep);
  K_CONNECT::IdTool::compress(sub_plan, pred_keep, nids);
  //std::cout << sub_plan << std::endl;
  // now update "pointers" (column references)
  for (E_Int ii = 0; ii < sub_plan.cols(); ++ii)
  {
    for (E_Int j = 0; j < sub_plan.rows(); ++j)
    {
      E_Int& p = sub_plan(j, ii);
      if (p <= 0 || p == IDX_NONE) continue;
      assert(p > 1); // first column is gone and is the only one to point on 2nd so must not encounter p==1
      p = nids[p];
      assert(p > 0); // must at least point on the second column
    }
  }
  //std::cout << sub_plan << std::endl;
}

///
template <typename arr_t>
E_Int join_plan<arr_t>::extract_compact_enabled_tree(const tree<arr_t>& PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, K_FLD::IntArray& plan)
{
  plan.clear();
  E_Int ret = __extract_compact_enabled_tree(PGtree, PGi, F, reverse, i0, plan);
  //std::cout << plan << std::endl;
  return ret;
}

///
template <typename arr_t> inline
void join_plan<arr_t>::extract_sub_plans(const tree<arr_t>& PGtree, E_Int PGi, const K_FLD::IntArray& plan, std::map<E_Int, K_FLD::IntArray>& odata)
{
  if (plan.cols() == 0) // no plan
    return;
  if (plan.rows() == 1) // no grand children
    return;

  const E_Int* children = PGtree.children(PGi);
  
#ifdef DEBU_JOIN_PLAN
  E_Int nbchildren = PGtree.nb_children(PGi);
  assert(nbchildren != 0);    //must have been subdivided or had children already
  //assert(nbchildren == 4);    // WARNING : ISO ASSUMPTION  : NBC = 4
#endif
  
  // at least one column (sized as NBC)
 
#ifdef DEBU_JOIN_PLAN
  E_Int nnodes = plan.cols();
  assert(nnodes >= 1);
  //assert(plan.rows() == nbchildren);
#endif
  
  K_FLD::IntArray sub_plan;

  for (E_Int n = 0; n <plan.rows(); ++n)
  {
    __extract_subplan(plan(n, 0), plan, sub_plan);
    if (sub_plan.cols() != 0) odata[children[n]] = sub_plan;
  }
}

///
template<typename arr_t> inline
bool join_plan<arr_t>::one_child_requires(const K_FLD::IntArray& plan)
{
  if (plan.cols() == 0) return false; // nothing planned

  bool require{ false };

  for (E_Int k = 0; (k < plan.rows()) && !require; ++k)
  {
    const E_Int& plan_for_childk = plan(k, 0);
    //require if grand children are planned : at least one of PGi's children is planned for having children
    require |= (plan_for_childk != NO_CHILDREN) && (plan_for_childk != IDX_NONE);
  }

  return require;
}

}
#endif
