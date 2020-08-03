/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_JOIN_PLAN_HXX
#define NUGA_JOIN_PLAN_HXX

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_unit.h"
#include "Nuga/include/IdTool.h"

#define NO_CHILDREN -1
#define NO_GRAND_CHILDREN -2

namespace NUGA
{

//
template <typename arr_t>
class join_plan
{
public:

  static E_Int extract_compact_tree(const tree<arr_t>& PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, arr_t& plan);
  
  static void extract_sub_plans(const tree<arr_t>& PGtree, E_Int PGi, const arr_t& plan, std::map<E_Int, arr_t>& odata);

  static bool one_child_requires(const arr_t& plan);

private:
  ///
  static void __extract_subplan(E_Int i, const arr_t& plan, arr_t& sub_plan);
  static E_Int __extract_compact_tree(const tree<arr_t> & PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, arr_t& plan);  
};

///
template <>
E_Int join_plan<K_FLD::IntArray>::__extract_compact_tree
(const tree < K_FLD::IntArray > & PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, K_FLD::IntArray& plan)
{
  //std::cout << "__extract_compact_tree : begin" << std::endl;
  E_Int nbc = PGtree.nb_children(PGi);
  if (nbc == 0) return NO_CHILDREN;

  //std::cout << "nbc : " << nbc << std::endl;

  const E_Int* pchild = PGtree.children(PGi);

  //std::cout << "pchild : " << pchild << std::endl;

  STACK_ARRAY(E_Int, nbc, children);
  for (E_Int i = 0; i < nbc; ++i) children[i] = pchild[i];

  //std::cout << "before F : children.get/reverse/i0" << children.get() << "/" << reverse << "/" << i0 << std::endl;
  // to put in receiver ref frame
  F(children.get(), reverse, i0);
  //std::cout << "after F" << std::endl;

  //to compact : dont create empty columns for all-leaves children
  E_Int nb_tot_grand_children{ 0 };
  for (E_Int i = 0; i < nbc; ++i)
  {
    //std::cout << "children[i] : " << children[i] << std::endl;
    nb_tot_grand_children += PGtree.nb_children(children[i]);
  }

  E_Int cols = plan.cols();

  if (nb_tot_grand_children == 0)
  {
    plan.pushBack(NO_CHILDREN, nbc);
    return  cols;
  }
 
  //std::cout << "nb_tot_grand_children : " << nb_tot_grand_children << std::endl;

  plan.pushBack(E_IDX_NONE, nbc);
  //std::cout << plan << std::endl;

  for (E_Int i = 0; i < nbc; ++i)
  {
    E_Int val = __extract_compact_tree(PGtree, children[i], F, reverse, i0, plan);
    //std::cout << "plan sz / cols /val : " << plan.cols() << "/" << cols<< "/" << val << std::endl;
    plan(i, cols) = val;
  }

  //std::cout << "__extract_compact_tree : end" << std::endl;

  return cols;
}

///
template <>
void join_plan<K_FLD::IntArray>::__extract_subplan
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
      if (v <= 0) continue;
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
  for (E_Int i = 0; i < sub_plan.cols(); ++i)
  {
    for (E_Int j = 0; j < sub_plan.rows(); ++j)
    {
      E_Int& p = sub_plan(j, i);
      if (p <= 0) continue;
      assert(p > 1); // first column is gone and is the only one to point on 2nd so must not encounter p==1
      p = nids[p];
      assert(p > 0); // must at least point on the second column
    }
  }
  //std::cout << sub_plan << std::endl;
}

///
template <>
E_Int join_plan<K_FLD::IntArray>::extract_compact_tree(const tree<K_FLD::IntArray>& PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, K_FLD::IntArray& plan)
{
  plan.clear();
  E_Int ret = __extract_compact_tree(PGtree, PGi, F, reverse, i0, plan);
  //std::cout << plan << std::endl;
  return ret;
}

///
template <>
void join_plan<K_FLD::IntArray>::extract_sub_plans(const tree<K_FLD::IntArray>& PGtree, E_Int PGi, const K_FLD::IntArray& plan, std::map<E_Int, K_FLD::IntArray>& odata)
{
  if (plan.cols() == 0) // no plan
    return;
  if (plan.rows() == 1) // no grand children
    return;

  const E_Int* children = PGtree.children(PGi);
  
#ifdef DEBU_JOIN_PLAN
  E_Int nbchildren =
#endif
    PGtree.nb_children(PGi);

#ifdef DEBU_JOIN_PLAN
  assert(nbchildren != 0);    //must have been subdivided or had children already
  assert(nbchildren == 4);    // WARNING : ISO ASSUMPTION  : NBC = 4
#endif
  // at least one column (sized as NBC)
 
#ifdef DEBU_JOIN_PLAN
  E_Int nnodes = plan.cols();
  assert(nnodes >= 1);
  assert(plan.rows() == nbchildren);
#endif
  
  K_FLD::IntArray sub_plan;

  for (E_Int n = 0; n <plan.rows(); ++n)
  {
    __extract_subplan(plan(n, 0), plan, sub_plan);
    if (sub_plan.cols() != 0) odata[children[n]] = sub_plan;
  }
}

///
template<>
bool join_plan<K_FLD::IntArray>::one_child_requires(const K_FLD::IntArray& plan)
{
  if (plan.cols() == 0) return false; // nothing planned

  bool require{ false };

  for (E_Int k = 0; (k < plan.rows()) && !require; ++k)
  {
    const E_Int& plan_for_childk = plan(k, 0);
    //require if grand children are planned : at least one of PGi's children is planned for having children
    require |= (plan_for_childk != NO_CHILDREN);
  }

  return require;
}

//////////////// ngon_unit impl //////////////////////////////
///
template <>
E_Int join_plan<ngon_unit>::extract_compact_tree(const tree<ngon_unit>& PGtree, E_Int PGi, reordering_func F, bool reverse, E_Int i0, ngon_unit& plan)
{
  //todo
  return 1;
}

///
template <>
void join_plan<ngon_unit>::extract_sub_plans(const tree<ngon_unit>& PGtree, E_Int PGi, const ngon_unit& plan, std::map<E_Int, ngon_unit>& odata)
{
  //todo
  return;
}

///
template <>
bool join_plan<ngon_unit>::one_child_requires(const ngon_unit& plan)
{
  //todo
  return false;
}

}
#endif
