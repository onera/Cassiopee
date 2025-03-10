
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

//#define FLAG_STEP
//#define DEBUG_XCELLN
//#define CLASSIFYER_DBG

#ifndef VISUAL
#include "intersector.h"
#endif

#include <map>
#include "Nuga/include/macros.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/masker.hxx"
#include "Nuga/include/xcelln.hxx"

#if defined(DEBUG_XCELLN) || defined(CLASSIFYER_DBG)
#include "Nuga/include/medit.hxx"
// #include "Nuga/include/NGON_debug.h"
// using NGDBG  = NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>;
std::string medith::wdir = "./";
#endif

using namespace NUGA;

#if ! defined(NETBEANS) && ! defined(VISUAL)

/// Default impl : for both output_type 0 (mask) and 1 (xcelln) for SURAFCE and VOLUME
template<typename classifyer_t>
void
pyMOVLP_XcellN
(const std::vector<K_FLD::FloatArray> &crds, const std::vector<K_FLD::IntArray>& cnts,
  const std::vector< std::vector<E_Int>> &zone_wall_ids,
  const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
  const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
  std::vector< std::vector<E_Int>> &mask_wall_ids,
  E_Float RTOL, const char* varString, char* eltType, PyObject* l)
{
  using outdata_t = typename classifyer_t::outdata_t;
  std::vector<outdata_t> xcelln;
  ///
  NUGA::MOVLP_xcelln_zones<classifyer_t>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xcelln, RTOL);

  //std::cout << "ng of xcelln upon exit : " << xcelln.size() << std::endl;
  E_Int nb_zones = crds.size();
  for (E_Int i = 0; i < nb_zones; ++i)
  {
    E_Int sz = xcelln[i].size();
    //std::cout << "size of xcelln for zone " << i << ": " << sz << std::endl;
    K_FLD::FloatArray xcellno(1, sz);
    for (E_Int j = 0; j < sz; ++j)xcellno(0, j) = xcelln[i][j];
    PyObject* tpl = K_ARRAY::buildArray(xcellno, varString, cnts[i], -1, eltType, false);
    //tpl = K_NUMPY::buildNumpyArray(&xcelln[i][0], xcelln[i].size(), 1);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
}

/// polyclip SURFACE impl (output type 2  : xcellno)
template<>
void
pyMOVLP_XcellN<NUGA::xcellno<pg_smesh_t, edge_mesh_t>>
(const std::vector<K_FLD::FloatArray> &crds, const std::vector<K_FLD::IntArray>& cnts,
  const std::vector< std::vector<E_Int>> &zone_wall_ids,
  const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
  const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
  std::vector< std::vector<E_Int>> &mask_wall_ids,
  E_Float RTOL, const char* varString, char* eltType, PyObject* l)
{
  using classifyer_t = NUGA::xcellno<pg_smesh_t, edge_mesh_t>;
  std::vector<typename classifyer_t::outdata_t> xmesh;
  ///
  NUGA::MOVLP_xcelln_zones<classifyer_t>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xmesh, RTOL);

  E_Int nb_zones = crds.size();
  for (E_Int i = 0; i < nb_zones; ++i)
  {
    //std::cout << "ng of cells upon exit : " << xmesh[i].mesh.ncells() << std::endl;
    //std::cout << "pushing out the mesh" << std::endl;
    // 1. ngon_unit => CASSIOPEE surface NGON
    K_FLD::IntArray cnto;
    ngon_type::export_surfacic_FN_to_EFN(xmesh[i].mesh.cnt, cnto);

    //std::cout << "buildarray mesh" << std::endl;
    PyObject *tpl = K_ARRAY::buildArray(xmesh[i].mesh.crd, varString, cnto, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // 2. history
    //std::cout << "history : " << xmesh[i].mesh.flag.size()  << std::endl;
    std::vector<E_Int> oids(xmesh[i].mesh.ncells(), E_IDX_NONE);
    for (size_t j = 0; j < oids.size(); ++j) oids[j] = xmesh[i].mesh.flag[j];

    //E_Int minf = *std::min_element(ALL(xmesh[i].mesh.flag));
    //E_Int maxf = *std::max_element(ALL(xmesh[i].mesh.flag));
    
    //std::cout << "buildArray histo : " << minf << "/" << maxf << std::endl;
    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
}

/// polyclip VOLUME impl (output type 2  : xcellno)
template<>
void
pyMOVLP_XcellN<NUGA::xcellno<ph_mesh_t, pg_smesh_t>>
(const std::vector<K_FLD::FloatArray> &crds, const std::vector<K_FLD::IntArray>& cnts,
  const std::vector< std::vector<E_Int>> &zone_wall_ids,
  const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
  const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
  std::vector< std::vector<E_Int>> &mask_wall_ids,
  E_Float RTOL, const char* varString, char* eltType, PyObject* l)
{
  using classifyer_t = NUGA::xcellno<ph_mesh_t, pg_smesh_t>;
  std::vector<typename classifyer_t::outdata_t> xmesh;
  ///
  NUGA::MOVLP_xcelln_zones<classifyer_t>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xmesh, RTOL);

  E_Int nb_zones = crds.size();
  for (E_Int i = 0; i < nb_zones; ++i)
  {
    //std::cout << "ng of cells upon exit : " << xmesh[i].mesh.ncells() << std::endl;
    //std::cout << "pushing out the mesh" << std::endl;
    // 1. ngon_unit => CASSIOPEE surface NGON
    K_FLD::IntArray cnto;
    xmesh[i].mesh.cnt.export_to_array(cnto);

    //std::cout << "buildarray mesh" << std::endl;
    PyObject *tpl = K_ARRAY::buildArray(xmesh[i].mesh.crd, varString, cnto, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // 2. history
    //std::cout << "history : " << xmesh[i].mesh.flag.size()  << std::endl;
    std::vector<E_Int> oids(xmesh[i].mesh.ncells(), E_IDX_NONE);
    for (size_t j = 0; j < oids.size(); ++j) oids[j] = xmesh[i].mesh.flag[j];

    //E_Int minf = *std::min_element(ALL(xmesh[i].mesh.flag));
    //E_Int maxf = *std::max_element(ALL(xmesh[i].mesh.flag));
    
    //std::cout << "buildArray histo : " << minf << "/" << maxf << std::endl;
    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
}

PyObject* K_INTERSECTOR::XcellN(PyObject* self, PyObject* args)
{
  PyObject *zones, *zwall_ids, *base_num, *masks, *priorities, *wall_ids;
  E_Int output_type(1);
  E_Float RTOL(0.05);

  if (!PYPARSETUPLE_(args, OOOO_ OO_ I_ R_, &zones, &zwall_ids, &base_num, &masks, 
    &wall_ids, &priorities, &output_type, &RTOL)) return NULL;

  E_Int nb_zones = PyList_Size(zones);
  //E_Int nb_wall_zones = PyList_Size(zwall_ids);
  E_Int nb_basenum = PyList_Size(base_num);
  E_Int nb_masks = PyList_Size(masks);
  E_Int nb_priority_pairs = PyList_Size(priorities);
  //E_Int nb_wall_sets = PyList_Size(wall_ids);

  //std::cout << "nb_zones/nb_wall_zones/nb_basenum/nb_masks/nb_priority_pairs/output_type : " << nb_zones << "/" << nb_wall_zones << "/" << nb_basenum << "/" << nb_masks << "/" << nb_priority_pairs << "/" << output_type << std::endl;
  
  if (nb_zones != nb_basenum)
  {
    std::cout << "xCellN: Info: nb zones / nb_basenum: " << nb_zones << "/" << nb_basenum << std::endl;
    PyErr_SetString(PyExc_ValueError,
       "XcellN: must have as many base ids as zones.");
      return NULL;
  }
  
  std::vector<K_FLD::FloatArray> crds(nb_zones), mask_crds(nb_masks);
  std::vector<K_FLD::IntArray>   cnts(nb_zones), mask_cnts(nb_masks);
  std::vector<std::vector<E_Int>>                mask_wall_ids(nb_masks), zone_wall_ids(nb_zones);

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

    err = getFromNGON(py_zone, crds[i], false, cnts[i], z_varString, z_eltType);
  }
  if (err)
  {
    PyErr_SetString(PyExc_TypeError,
                    "xCellN: at least one input NGON is invalid.");
    return NULL;
  }

  // get the masks (BAR in 2D, nuga NGON in 3D)
  char* msk_varString, *msk_eltType;
  std::vector<std::string> types;
  std::string s1("BAR");
  std::string s2("NGON");
  types.push_back(s1);
  types.push_back(s2);

  for (E_Int i=0; (i < nb_masks) && !err; ++i)
  {
    PyObject* py_mask = PyList_GetItem(masks, i);
    if (py_mask == Py_None) continue;

    err = get_of_type(types, py_mask, mask_crds[i], false, mask_cnts[i], msk_varString, msk_eltType);
    // std::cout << "mask sizes : " << mask_crds[i].cols() << " points" << std::endl;
    // if (strcmp("NGON", msk_eltType) == 0)
    // {
    //   ngon_type ng(mask_cnts[i]);
    //   std::cout << "mask sizes : "<< ng.PGs.size() << " cells" << std::endl;
    // }
    // else
    //   std::cout << "mask sizes : "<< mask_cnts[i].cols() << " cells" << std::endl;
  }
  if (err)
  {
    PyErr_SetString(PyExc_TypeError,
                    "xCellN: one mask is invalid.");
    return NULL;
  }

  bool DIM3 = (strcmp("NGON", msk_eltType) == 0);
  //std::cout << "DIM3 ? " << DIM3 << std::endl;

  // get the zone wall ids
  for (E_Int i=0; i < nb_zones; ++i)
  {
    PyObject* py_zwall_ids = PyList_GetItem(zwall_ids, i);
    //std::cout << "py_zwall_ids" << py_zwall_ids << std::endl;
    if (py_zwall_ids != Py_None)
    {
      E_Int nfld, sz, *data;
      E_Int ok =  K_NUMPY::getFromNumpyArray(py_zwall_ids, data, sz, nfld, 1/*shared*/);
      if (ok == 0) continue;

      zone_wall_ids[i].insert(zone_wall_ids[i].end(), data, data+sz);

      //std::cout << "zone wall retrieved : ok ? " << ok <<"/" << nfld << "/" << sz << std::endl;
    }
  }

  // get the mask wall ids
  assert(wall_ids != Py_None);
  for (E_Int i=0; i < nb_masks; ++i)
  {
    PyObject* py_wall_ids = PyList_GetItem(wall_ids, i);
    if (py_wall_ids != Py_None)
    {
      E_Int nfld, sz, *data;

#ifdef DEBUG_XCELLN
      E_Int ok =  
#endif
      K_NUMPY::getFromNumpyArray(py_wall_ids, data, sz, nfld, 1/*shared*/);

      mask_wall_ids[i].insert(mask_wall_ids[i].end(), data, data+sz);

      //std::cout << "wall retrieved : ok ? " << ok <<"/" << nfld << "/" << sz << std::endl;
    }
  }

  if (err)
  {
    PyErr_SetString(PyExc_TypeError,
                    "xCellN: one error occured in input.");
    return NULL;
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

  PyObject *l(PyList_New(0));

  if (!DIM3)
  {
    if (output_type == 0)
      pyMOVLP_XcellN<NUGA::masker<pg_smesh_t, edge_mesh_t>>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, "ratio", z_eltType, l);
    else if (output_type == 1)
      pyMOVLP_XcellN<NUGA::xcellnv<pg_smesh_t, edge_mesh_t>>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, "ratio", z_eltType, l);
    else if (output_type == 2)    
      pyMOVLP_XcellN<NUGA::xcellno<pg_smesh_t, edge_mesh_t>>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, z_varString, z_eltType, l);
  }
  else
  {
    if (output_type == 0)
      pyMOVLP_XcellN<NUGA::masker<ph_mesh_t, pg_smesh_t>>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, "ratio", z_eltType, l);
    else if (output_type == 1)
      pyMOVLP_XcellN<NUGA::xcellnv<ph_mesh_t, pg_smesh_t>>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, "ratio", z_eltType, l);
    else if (output_type == 2)    
      pyMOVLP_XcellN<NUGA::xcellno<ph_mesh_t, pg_smesh_t>>(crds, cnts, zone_wall_ids, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, z_varString, z_eltType, l);
  }

  return l;

}
#endif

