
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
#include "Nuga/include/masker.hxx"
#include "Nuga/include/xcelln.hxx"

// //#define FLAG_STEP
// #define DEBUG_XCELLN

#ifdef DEBUG_XCELLN
#include "Nuga/include/medit.hxx"
// #include "Nuga/Boolean/NGON_debug.h"
// using NGDBG  = NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>;
std::string medith::wdir = "";
#endif

#if !  defined(NETBEANS) && ! defined(VISUAL)

using zmesh_t = NUGA::pg_smesh_t;
using bmesh_t = NUGA::edge_mesh_t;

/// Default impl : for both output_type 0 (mask) and 1 (xcelln)
template<typename classifyer_t>
void
pyMOVLP_XcellNSurf
(const std::vector<K_FLD::FloatArray*> &crds, const std::vector<K_FLD::IntArray*>& cnts,
  const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
  const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
  std::vector< std::vector<E_Int>> &mask_wall_ids,
  E_Float RTOL, char* varString, char* eltType, PyObject* l)
{
  using outdata_t = typename classifyer_t::outdata_t;
  std::vector<outdata_t> xcelln;
  ///
  NUGA::MOVLP_xcelln_zones<classifyer_t>(crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xcelln, RTOL);

  //std::cout << "ng of xcelln upon exit : " << xcelln.size() << std::endl;
  E_Int nb_zones = crds.size();
  for (E_Int i = 0; i < nb_zones; ++i)
  {
    E_Int sz = xcelln[i].size();
    //std::cout << "size of xcelln for zone " << i << ": " << sz << std::endl;
    K_FLD::FloatArray xcellno(1, sz);
    for (E_Int j = 0; j < sz; ++j)xcellno(0, j) = xcelln[i][j];
    PyObject* tpl = K_ARRAY::buildArray(xcellno, varString, *cnts[i], -1, eltType, false);
    //tpl = K_NUMPY::buildNumpyArray(&xcelln[i][0], xcelln[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
}

/// polyclip impl (output type 2  : xcellno
template<>
void
pyMOVLP_XcellNSurf<NUGA::xcellno<zmesh_t, bmesh_t>>
(const std::vector<K_FLD::FloatArray*> &crds, const std::vector<K_FLD::IntArray*>& cnts,
  const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
  const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
  std::vector< std::vector<E_Int>> &mask_wall_ids,
  E_Float RTOL, char* varString, char* eltType, PyObject* l)
{
  using classifyer_t = NUGA::xcellno<zmesh_t, bmesh_t>;
  std::vector<typename classifyer_t::outdata_t> xmesh;
  ///
  NUGA::MOVLP_xcelln_zones<classifyer_t>(crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xmesh, RTOL);

  //std::cout << "ng of xcelln upon exit : " << xcelln.size() << std::endl;
  E_Int nb_zones = crds.size();
  for (E_Int i = 0; i < nb_zones; ++i)
  {
    // pushing out the mesh
    // 1. ngon_unit => CASSIOPEE surface NGON
    K_FLD::IntArray cnto;
    ngon_type::export_surfacic_FN_to_EFN(xmesh[i].mesh.cnt, cnto);

    // 2. history
    std::vector<E_Int> oids(xmesh[i].mesh.cnt.size(), E_IDX_NONE);
    for (E_Int j = 0; j < oids.size(); ++j) oids[j] = xmesh[i].mesh.cnt._ancEs(0, j);
    
    PyObject *tpl = K_ARRAY::buildArray(xmesh[i].mesh.crd, varString, cnto, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
}

PyObject* K_INTERSECTOR::XcellNSurf(PyObject* self, PyObject* args)
{
  PyObject *zones, *base_num, *masks, *priorities, *wall_ids;
  E_Int output_type(1);
  E_Float RTOL(0.05);

  if (!PyArg_ParseTuple(args, "OOOOOld", &zones, &base_num, &masks, &wall_ids, &priorities, &output_type, &RTOL)) return NULL;

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

  PyObject *l(PyList_New(0));

  if (output_type == 0)
    pyMOVLP_XcellNSurf<NUGA::masker<zmesh_t, bmesh_t>>(crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, "xcelln", z_eltType, l);
  else if (output_type == 1)
    pyMOVLP_XcellNSurf<NUGA::xcellnv<zmesh_t, bmesh_t>>(crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, "xcelln", z_eltType, l);
  else if (output_type == 2)
    pyMOVLP_XcellNSurf<NUGA::xcellno<zmesh_t, bmesh_t>>(crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, RTOL, z_varString, z_eltType, l);

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

