/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __BAR_CONFORMIZER_CXX__
#define __BAR_CONFORMIZER_CXX__

#include "BAR_Conformizer.h"

namespace NUGA
{

///  
template <E_Int DIM>
std::vector<int> BAR_Conformizer<DIM>::get_x_history()
{
  std::vector<int> xhis;
  size_t nb_edges = parent_type::_elements.size();
  int maxid = -1;
  for (size_t i = 0; i < nb_edges; ++i)
  {
    E2& e = parent_type::_elements[i];
    size_t sz = e.nodes.size();
    if (sz == 0)
      continue;
    maxid = std::max(maxid, *std::max_element(ALL(e.nodes)));
  }
  
  xhis.resize(maxid+1, IDX_NONE);

  for (size_t i = 0; i < nb_edges; ++i)
  {
    E2& e = parent_type::_elements[i];
    size_t sz = e.nodes.size();
    if (sz == 0)
      continue;
    for (size_t j = 0; j < sz; ++j)
    {
      int xid = e.nodes[j];
      if (xhis[xid] == IDX_NONE) xhis[xid] = e.id; // keep first (prior to lower ids)
    }
  }

  return xhis;
}

///
template <E_Int DIM>
void
BAR_Conformizer<DIM>::__set_tolerances(E_Float Lmin, E_Float Lmax, E_Float  user_tolerance)
{  
    
  // min edge length (MEL) : Lmin
  
  // The current accpeted range for tolerance is between 1e-10% and 1.% of MEL
  E_Float tol_min = EPSILON*Lmin;
  E_Float tol_max = 0.01*Lmin;
  
  if (user_tolerance < 0.)
    user_tolerance = -user_tolerance*Lmin;
  
  if (user_tolerance > 0.)
  {
    parent_type::_tolerance = std::min(user_tolerance, tol_max);
    parent_type::_tolerance = std::max(user_tolerance, tol_min);
  }
  else
    parent_type::_tolerance = tol_max;
    
  // In case of bar it is supposedly easier to make the tolerances consistent so hopefully the following will work :
  parent_type::_tol_x=parent_type::_tol_clean=parent_type::_tolerance;
}
  
///
template <E_Int DIM>
E_Int
BAR_Conformizer<DIM>::__intersect
(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E2& e1, E2& e2, E_Float tolerance)
{
  
  K_FLD::IntArray::const_iterator pS1 = connect.col(e1.id);
  K_FLD::IntArray::const_iterator pS2 = connect.col(e2.id);
  
  E_Bool overlap;
  E_Float u[4];
  E_Int N[]={*pS1, *(pS1+1), *pS2, *(pS2+1)};
  E2* EE[] = {&e1,&e2};
  
  bool intersect = K_MESH::Edge::intersect<DIM>(pos, N[0],N[1], N[2], N[3], 
                               tolerance, parent_type::_absolute_tol, u[0]/*u00*/, u[1]/*u01*/, u[2]/*u10*/, u[3]/*u11*/, overlap); 
  if (!intersect)
    return false;
  
  E_Float /*Edge[DIM][2], */IP[DIM], Edg1[DIM], Edg2[DIM];
  NUGA::diff<DIM>(pos.col(*(pS1+1)), pos.col(*pS1), /*Edge[0]*/Edg1);
  NUGA::diff<DIM>(pos.col(*(pS2+1)), pos.col(*pS2), /*Edge[1]*/Edg2);
  
  E_Float tol[] ={tolerance, tolerance};
  if (parent_type::_absolute_tol)
  {
    tol[0] /= ::sqrt(NUGA::sqrNorm<DIM>(/*Edge[0]*/Edg1));
    tol[1] /= ::sqrt(NUGA::sqrNorm<DIM>(/*Edge[1]*/Edg2));
  }
  
  E_Int n1, Ni, ret = false;
  for (E_Int n = 0; n < 4; ++n)
  {
    n1 = (n < 2) ? 0 : 1;
    
    if ((u[n] < tol[n1]) || (u[n] > (1. - tol[n1])))
      continue;
    
    ret = true;
    
    if (n1 == 0)
    {
    
    for (E_Int k = 0; k < DIM; ++k)
      IP[k] = pos(k, N[n]) + u[n]*Edg1[k];
    }
    else
    {
      for (E_Int k = 0; k < DIM; ++k)
      IP[k] = pos(k, N[n]) + u[n]*Edg2[k];
    }
    
    pos.pushBack(IP, IP+DIM);//addPoint(IP, pos);
    Ni = pos.cols()-1;
        
    EE[n1]->nodes.push_back(Ni);
  }
  
  return ret;
}

///
template <E_Int DIM>
void
BAR_Conformizer<DIM>::__update_data
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& newIDs)
{
  if (newIDs.empty()) return;

  std::set<E_Int> pool;
  E_Int nb_edges = parent_type::_elements.size();
  for (E_Int i = 0; i < nb_edges; ++i)
  {
    E2& e = parent_type::_elements[i];
    pool.clear();

    size_t sz = e.nodes.size();
    if (sz == 0)
      continue;
    for (size_t j = 0; j < sz; ++j)
      pool.insert(newIDs[e.nodes[j]]);
    e.nodes.clear();
    e.nodes.insert(e.nodes.begin(), pool.begin(), pool.end());
    const E_Float *P0 = pos.col(connect(0,e.id));
    __reorder_nodes_on_edge(pos, e.nodes, P0);
  }
}

template <E_Int DIM>
void
BAR_Conformizer<DIM>::__reorder_nodes_on_edge
(const K_FLD::FloatArray& pos, std::vector<E_Int>& nodes, const E_Float *P0)
{
  E_Float L;
  E_Int   E[2];

  size_t nb_nodes = nodes.size();

  parent_type::_sorterFI.clear();
  
  for (size_t k = 0; k < nb_nodes; ++k)
  {
    L = NUGA::sqrDistance(pos.col(nodes[k]), P0, DIM);
    parent_type::_sorterFI.push_back(std::make_pair(L, nodes[k]));
  }

  std::sort(parent_type::_sorterFI.begin(), parent_type::_sorterFI.end());

  nodes.clear();
  E[0] = parent_type::_sorterFI[0].second;
  nodes.push_back(E[0]);

  for (size_t s = 0; s < parent_type::_sorterFI.size()-1; ++s)
  {
    E[1] = parent_type::_sorterFI[s+1].second;

    if (E[1] != E[0])
    {
      nodes.push_back(E[1]);
      E[0] = E[1];
    }
  }
}

///  
template <E_Int DIM>
E_Int
BAR_Conformizer<DIM>::__split_Elements
(const K_FLD::FloatArray& pos, K_FLD::IntArray & connect,
 NUGA::bool_vector_type& xc,
 NUGA::int_vector_type& ancestors)
{
  K_FLD::IntArray::const_iterator pS;
  K_FLD::IntArray connectOut;
  NUGA::int_vector_type ancOut;
  NUGA::bool_vector_type xcOut;
  E_Int E[2];
  
#define ADD_ONE_EDGE(pS, anci, xci) \
{connectOut.pushBack(pS, pS+2);ancOut.push_back(anci);xcOut.push_back(xci);}
  
  E_Int nb_edges = parent_type::_elements.size();
  for (E_Int i = 0; i < nb_edges; ++i)
  {
    E2& e = parent_type::_elements[i];
    pS = connect.col(e.id);
    
    if (e.nodes.empty())//nothing to do
      ADD_ONE_EDGE(pS, ancestors[i], false)
    else
    {
      E[0]=*pS;
      E[1]=e.nodes[0];
      ADD_ONE_EDGE(E, ancestors[i], true)
      
      size_t sz = e.nodes.size();
      for (size_t j = 0; j < (sz-1); ++j)
      {
        E[0]=e.nodes[j];
        E[1]=e.nodes[j+1];
        ADD_ONE_EDGE(E, ancestors[i], true)
      }
      
      E[0]=e.nodes[sz-1];
      E[1]=*(pS+1);
      ADD_ONE_EDGE(E, ancestors[i], true)      
    } 
  }
  connect = connectOut;
  xc = xcOut;
  ancestors = ancOut;
  
  return 0;
}

}

#endif
