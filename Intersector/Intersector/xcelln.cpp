
/*    
    Copyright 2013-2020 Onera.

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

#ifndef VISUAL
#include "intersector.h"
#endif

#include <map>
#include "Nuga/include/macros.h"
#include "Fld/ngon_t.hxx"
#include "Nuga/include/mesh_t.hxx"
//#include "Nuga/include/classifyer.hxx"
#include "Nuga/include/masker.hxx"
#include "Nuga/include/xcelln.hxx"
// #include <vector>
// #include "Fld/DynArray.h"
// 
// #include "MeshElement/Hexahedron.h"
// #include "MeshElement/Polyhedron.h"
// #include "Search/BbTree.h"
// #include "Nuga/include/localizer.hxx"
// #include "Nuga/Delaunay/Triangulator.h"
// #include "Nuga/include/ph_clipper.hxx"
// #include "Nuga/include/polyhedron.hxx"

// //#define FLAG_STEP
// //#define OUTPUT_XCELLN

// #ifdef FLAG_STEP
// #include "chrono.h"
// int chrono::verbose = 2;
// #endif

// #ifdef NETBEANSZ
// //#define DEBUG_XCELLN
// #endif

// #ifdef OUTPUT_XCELLN
// #include <sstream>
// #include "Nuga/include/medit.hxx"
// #endif

// #ifdef DEBUG_XCELLN
// #include "Nuga/include/medit.hxx"
// #include "Nuga/Boolean/NGON_debug.h"
// using NGDBG  = NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>;
// #endif


 using IntVec = std::vector<E_Int>;
 using prior_t = std::map< E_Int, IntVec>;

 using ii_pair_t = std::pair<E_Int, E_Int>;
 using no_ii_pair_t = K_MESH::NO_Edge;

struct aIsLessThanb : std::binary_function<int, int, bool>
{
  
  aIsLessThanb(const prior_t p):_p(p){}
  
  bool operator()(int a, int b) const { 
    
    auto it = _p.find(a);
    if (it == _p.end()) return false;
    
    for (size_t i=0; i < it->second.size(); ++i)
      if (it->second[i] == b) return true;
    
    return false;
  }
  
  prior_t _p;
};
   
inline void comp_priorities(const std::vector<ii_pair_t> & priority,  prior_t & decrease_prior_per_comp, IntVec& rank_wnps)
{
  // WARNING DECREASING PRIORITY UPON EXIT
  
  // 0. get the max comp id
  E_Int max_comp_id = -1;
  for (auto it = priority.begin(); it != priority.end(); ++it)
  {
    max_comp_id = std::max(max_comp_id, it->first);
    max_comp_id = std::max(max_comp_id, it->second);
  }
  
  E_Int nb_comps = max_comp_id + 1;

  decrease_prior_per_comp.clear();
  rank_wnps.resize(nb_comps, 0);

  // // 1. CLEANING INPUT (discard redundancy & bilateral)
  std::set<no_ii_pair_t> upairs;
  for (auto it = priority.begin(); it != priority.end(); ++it) // make pairs unique and unilateral
    upairs.insert(no_ii_pair_t(it->first, it->second));
  // // in case of bilateral input keep the first
  std::vector<ii_pair_t> clean_priority;
  for (auto it = priority.begin(); it != priority.end(); ++it)
  {
    no_ii_pair_t p(it->first, it->second);
    if (upairs.find(p) == upairs.end()) continue; // p is redundant here and first occ has been remove in upairs so discard p
    upairs.erase(p);
    clean_priority.push_back(*it);
  }
  
  // // 2. 
  // // store the priors per comp
  for (size_t i = 0; i < clean_priority.size(); ++i)
  {
    E_Int Lcompid = clean_priority[i].first;
    E_Int Rcompid = clean_priority[i].second;

    //std::cout << " pair : " << Lcompid << "/" << Rcompid << std::endl;
    
    decrease_prior_per_comp[Lcompid].push_back(Rcompid);
    ++rank_wnps[Lcompid];
    
    decrease_prior_per_comp[Rcompid].push_back(Lcompid); //opposite pair : for considering WNPs
  }

  // IntVec all_comps;
  // for (auto it = decrease_prior_per_comp.begin(); it != decrease_prior_per_comp.end(); ++it)
  //   all_comps.push_back(it->first);

  aIsLessThanb pred(decrease_prior_per_comp);
  rank_wnps.resize(1+decrease_prior_per_comp.rbegin()->first, 0); //sized as the max comp num

  E_Int i=0;
  for (auto it = decrease_prior_per_comp.begin(); it != decrease_prior_per_comp.end(); ++it, ++i)
  {
    E_Int compid = it->first;
    std::sort(ALL(it->second), pred);
    std::reverse(ALL(it->second)); // WARNING DECREASING PRIORITY DONE HERE
    //std::cout << "WNP rank for comp " << it->first << " is " << rank_wnps[compid] << std::endl;
  }
}

using zmesh_t = NUGA::pg_smesh_t;
using bmesh_t = NUGA::edge_mesh_t;
//using zmesh_t = NUGA::ph_mesh_t;
//using bmesh_t = NUGA::pg_smesh_t;

void compute_zone(K_FLD::FloatArray& z_crd, K_FLD::IntArray& z_cnt,
                  const IntVec& z_priorities, E_Int rank_wnp,
                  const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
                  std::vector< std::vector<E_Int>> &mask_wall_ids, 
                  std::vector<E_Float>& z_xcelln, bool binary_mode, E_Float col_X, E_Float RTOL)
{
  zmesh_t  z_mesh(z_crd, z_cnt);  // polygonal surface mesh
  std::vector<bmesh_t*> mask_meshes;
  
  if (binary_mode)
  {
    NUGA::masker<zmesh_t, bmesh_t> classs(RTOL, col_X);
    classs.prepare(z_mesh, mask_crds, mask_cnts, mask_wall_ids, z_priorities, rank_wnp, mask_meshes);
    classs.compute(z_mesh, mask_meshes, z_xcelln);
  }
  else
  {
    NUGA::xcellnv<zmesh_t, bmesh_t> classs(RTOL);
    classs.prepare(z_mesh, mask_crds, mask_cnts, mask_wall_ids, z_priorities, rank_wnp, mask_meshes);
    classs.compute(z_mesh, mask_meshes, z_xcelln);
  }
  
  
#ifdef DEBUG_XCELLN
  std::cout << "TERMINATE" << std::endl;
#endif

  for (size_t u=0; u < mask_meshes.size(); ++u)
    delete mask_meshes[u];
}

void MOVLP_XcellN(const std::vector<K_FLD::FloatArray*> &crds, const std::vector<K_FLD::IntArray*>& cnts,
                  const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority, 
                  const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
                  std::vector< std::vector<E_Int>> &mask_wall_ids, 
                  std::vector< std::vector<E_Float>>& xcelln, bool binary_mode, E_Float col_X, E_Float RTOL)
{
  //priority map
  prior_t sorted_comps_per_comp;
  IntVec rank_wnps;
  comp_priorities(priority, sorted_comps_per_comp, rank_wnps);

  E_Int nb_zones = crds.size();

  xcelln.resize(nb_zones);
  for (E_Int u=0; u < nb_zones; ++u)
    if(cnts[u] != nullptr) xcelln[u].resize((*cnts[u])[0], 1.);

#ifndef NETBEANSZ
//#pragma omp parallel for
#endif
  for (E_Int z=0; z < nb_zones; ++z)
  { 
    //std::cout << "processing zone : " << z << " from comp " << comp_id[z] << std::endl;

    auto it = sorted_comps_per_comp.find(comp_id[z]);
    if (it == sorted_comps_per_comp.end()) continue;

    E_Int z_rank_wnp = rank_wnps[it->first];
    IntVec& z_priorities = it->second;

    //if (z != 26) continue;
    //std::cout << "calling compute_zone for zone : " << z << " over " << nb_zones << std::endl;
    
    compute_zone(*crds[z], *cnts[z], z_priorities, z_rank_wnp, mask_crds, mask_cnts, mask_wall_ids, xcelln[z], binary_mode, col_X, RTOL);

    //fixme : hack for fully inside non prior
    // inferior but uncolored => assume it means fully in so IN
    bool full_out = (*std::min_element(ALL(xcelln[z])) == 1.);
    bool is_inferior = (rank_wnps[comp_id[z]] == z_priorities.size() && !z_priorities.empty());
    
    if (full_out && is_inferior) 
    {
      //std::cout << "full OUT rank : " << rank_wnps[comp_id[z]] << std::endl;
      E_Int sz = xcelln[z].size();
      xcelln[z].clear();
      xcelln[z].resize(sz, NUGA::IN);
    }
  }
}


#if !  defined(NETBEANS) && ! defined(VISUAL)
PyObject* K_INTERSECTOR::XcellNSurf(PyObject* self, PyObject* args)
{
  PyObject *zones, *base_num, *masks, *priorities, *wall_ids;
  E_Int binary_mode(1);
  E_Float col_X(0.5), RTOL(0.05);

  if (!PyArg_ParseTuple(args, "OOOOOldd", &zones, &base_num, &masks, &wall_ids, &priorities, &binary_mode, &col_X, &RTOL)) return NULL;

  E_Int nb_zones = PyList_Size(zones);
  E_Int nb_basenum = PyList_Size(base_num);
  E_Int nb_masks = PyList_Size(masks);
  E_Int nb_priority_pairs = PyList_Size(priorities);
  //E_Int nb_wall_sets = PyList_Size(wall_ids);

  //std::cout << "zones/masks/wallids/prior_pairs/base_num/binary/colX : " << nb_zones << "/" << nb_masks << "/" << nb_wall_sets << "/" << nb_priority_pairs << "/" << nb_basenum << "/" << binary_mode << "/" << col_X << std::endl;

  if (nb_zones != nb_basenum)
  {
    std::cout << "nb zones vs nb_basenum : " << nb_zones << "/" << nb_basenum << std::endl;
    PyErr_SetString(PyExc_ValueError,
       "XcellN: must have as many base ids as zones.");
      return NULL;
  }
  
  std::vector<K_FLD::FloatArray*> crds(nb_zones, nullptr), mask_crds(nb_masks, nullptr);
  std::vector<K_FLD::IntArray*>   cnts(nb_zones, nullptr), mask_cnts(nb_masks, nullptr);
  std::vector<std::vector<E_Int>>                          mask_wall_ids(nb_masks);

  std::vector<std::pair<E_Int, E_Int>> priority;
  std::vector<E_Int> comp_id;
  std::vector< std::vector<E_Float>> xcelln(nb_zones);

  E_Int err(0);
 
  // get the zones
  char* z_varString, *z_eltType;
  for (E_Int i=0; (i < nb_zones) && !err; ++i)
  {
    //std::cout << "getting zone in list : " << i << std::endl;
    PyObject* py_zone = PyList_GetItem(zones, i);

    err = check_is_NGON(py_zone, crds[i], cnts[i], z_varString, z_eltType);
  }

  // get the masks
  char* msk_varString, *msk_eltType;
  for (E_Int i=0; (i < nb_masks) && !err; ++i)
  {
    PyObject* py_mask = PyList_GetItem(masks, i);

    err = check_is_BAR(py_mask, mask_crds[i], mask_cnts[i], msk_varString, msk_eltType);
    //std::cout << "mask sizes : " << mask_crds[i]->cols() << " points" << std::endl;
    //std::cout << "mask sizes : " << mask_cnts[i]->cols() << " cells" << std::endl;
  }

  // get the wall ids
  for (E_Int i=0; (i < nb_masks) && !err; ++i)
  {
    PyObject* py_wall_ids = PyList_GetItem(wall_ids, i);
    if (py_wall_ids != Py_None)
    {
      E_Int nfld, sz, *data;
      /*E_Int ok =  */K_NUMPY::getFromNumpyArray(py_wall_ids, data, sz, nfld, 1/*shared*/, 0 /*inverse*/);

      mask_wall_ids[i].insert(mask_wall_ids[i].end(), data, data+sz);

      //std::cout << "wall retrieved : ok ? " << ok <<"/" << nfld << "/" << sz << std::endl;
    }
  }

  if (err)
  {
    for (E_Int i=0; i < nb_zones; ++i)
    {
      delete crds[i];
      delete cnts[i];
    }
    for (E_Int i=0; i < nb_masks; ++i)
    {
      delete mask_crds[i];
      delete mask_cnts[i];
    }

    return nullptr;
  }

  // get the priority pairs
  for (E_Int i=0; i < nb_priority_pairs; ++i)
  {
    //std::cout << "getting priority pairs in list : " << i << std::endl;
    PyObject* py_priority = PyList_GetItem(priorities, i);
    if (!PyTuple_Check(py_priority)) continue;

    PyObject* s = PyTuple_GetItem(py_priority, 0);
    E_Int p0 = PyLong_AsLong(s);
    s = PyTuple_GetItem(py_priority, 1);
    E_Int p1 = PyLong_AsLong(s);

    priority.push_back(std::make_pair(p0,p1));
  }

  // get the base_num (zone 's component id')
  for (E_Int i=0; i < nb_zones; ++i)
  {
    //std::cout << "getting priority pairs in list : " << i << std::endl;
    PyObject* py_comp_id = PyList_GetItem(base_num, i);
    E_Int c_id = PyLong_AsLong(py_comp_id);
    comp_id.push_back(c_id);
    //std::cout << "zone/comp : " << i << "/" << c_id << std::endl;
  }

  //std::cout << "calling MOVLP_XcellN..." << std::endl;

  ///
  MOVLP_XcellN (crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xcelln, binary_mode, col_X, RTOL);

  //std::cout << "ng of xcelln upon exit : " << xcelln.size() << std::endl;

  PyObject *l(PyList_New(0)), *tpl;

  for (E_Int i = 0; i < nb_zones; ++i)
  {
    E_Int sz = xcelln[i].size();
    //std::cout << "size of xcelln for zone " << i << ": " << sz << std::endl;
    K_FLD::FloatArray xcellno(1,sz);
    for (E_Int j = 0; j < sz; ++j)xcellno(0,j)=xcelln[i][j];
    tpl = K_ARRAY::buildArray(xcellno, "xcelln", *cnts[i], -1, z_eltType, false);
    //tpl = K_NUMPY::buildNumpyArray(&xcelln[i][0], xcelln[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  for (E_Int i=0; i < nb_zones; ++i)
  {
    delete crds[i];
    delete cnts[i];
  }
  for (E_Int i=0; i < nb_masks; ++i)
  {
    delete mask_crds[i];
    delete mask_cnts[i];
  }

  return l;

}
#endif

