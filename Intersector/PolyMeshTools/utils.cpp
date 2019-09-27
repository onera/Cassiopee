/*    
    Copyright 2013-2019 Onera.

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

# include "intersector.h"
# include <string>
# include <sstream> 

//#define FLAG_STEP
#ifdef FLAG_STEP
#include "chrono.h"
E_Int chrono::verbose=1;
#endif

# include "Fld/ngon_t.hxx"
# include "Nuga/Delaunay/Triangulator.h"
#include "Nuga/include/localizer.hxx"
#include "Nuga/include/collider.hxx"

//#include <iostream>
#include <memory>

using namespace std;
using namespace K_FLD;


//=============================================================================
/* Creates 4 zones : 1) uncomputable polygons 2) uncomputable polyhedra 
   3) uncomputable polyhedra & neighbors 4) complementary of 3) */
//=============================================================================
PyObject* K_INTERSECTOR::extractUncomputables(PyObject* self, PyObject* args)
{
  E_Int neigh_level(1);
  PyObject *arr;

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &neigh_level)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), uphs/*uncomputable phs*/, uphs_wv1/*uncomputable phs with neighbors*/, remaining;
  ngon_unit upgs; //uncomputable polygons

  err = ngon_type::extract_uncomputables<DELAUNAY::Triangulator>(crd, ngi, neigh_level, upgs, uphs, uphs_wv1, remaining);
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "extract_uncomputables failed.");
    delete f; delete cn;
    return NULL;
  } 
  
  PyObject *l(PyList_New(0)), *tpl;

  if (upgs.size() == 0)
  {
    std::cout << "OK : there are no uncomputable polygons." << std::endl;
    PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnt, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else
  {
    // zone 1 : uncomputable pgs
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(upgs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      ngon_type ngo(upgs, false/*one ph per pg*/);
      ngo.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    // zone 2 : uncomputable phs
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(uphs.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      uphs.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    // zone 3 : uncomputable phs and first neighborhood
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(uphs_wv1.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      uphs_wv1.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    // zone 4 : complementary of above selection
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(remaining.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      remaining.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }
  
  
  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX*/
//=============================================================================
PyObject* K_INTERSECTOR::extractPathologicalCells(PyObject* self, PyObject* args)
{
  E_Int neigh_level(2);
  PyObject *arr;

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &neigh_level)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), neigh_phs/*patho neighbors*/, remaining_phs;
  std::vector<ngon_type> phsv;

  //std::cout << "neigh level : " <<  neigh_level << std::endl;
 
  err = ngon_type::extract_pathological_PHs<DELAUNAY::Triangulator>(crd, ngi, neigh_level, 
                                                                    phsv, neigh_phs, remaining_phs);
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "extract_pathological_phs failed.");
    delete f; delete cn;
    return NULL;
  } 
  
  PyObject *l(PyList_New(0)), *tpl;

  if (phsv.empty())
  {
    std::cout << "OK : all the cells are star-shaped regarding there centroids." << std::endl;
    tpl = K_ARRAY::buildArray(crd, varString, cnt, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else
  {
    E_Int nphs = phsv[0].PHs.size();
    //if (nphs)
    // zone 1 : 
    {
      if (nphs) std::cout << "there are " << nphs << " open cells (bug somewhere)." << std::endl;

      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(phsv[0].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[0].export_to_array(cnto);
      //patho_name = "open_cells";
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    nphs = phsv[1].PHs.size();
    //if (nphs)
    // zone 2 : 
    {
      if (nphs) std::cout << "there are " << nphs << " cells with degenerated polygons (showstopper?)." << std::endl;

      K_FLD::FloatArray crdtmp(crd);

      ngon_type::compact_to_used_nodes(phsv[1].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[1].export_to_array(cnto);

      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

    nphs = phsv[2].PHs.size();
    //if (nphs)
    // zone 3 : 
    {
      if (nphs) std::cout << "there are " << phsv[2].PHs.size() << " cells with some  delaunay-failure polygons." << std::endl;
      K_FLD::FloatArray crdtmp(crd);

      ngon_type::compact_to_used_nodes(phsv[2].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[2].export_to_array(cnto);

      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    nphs = phsv[3].PHs.size();
    if (nphs) std::cout << "there are " << phsv[3].PHs.size() << " non-centroid-star-shaped cells that can be split." << std::endl;
    // zone 4 : 
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(phsv[3].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[3].export_to_array(cnto);

      //std::cout << "pg : pg " << phsv[3].PGs.size() << " versus " << cnto[0] << std::endl;

      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

    // zone 5 : neighbors
    // {
    //   K_FLD::FloatArray crdtmp(crd);
    //   ngon_type::compact_to_used_nodes(neigh_phs.PGs, crdtmp); //reduce points
    //   //export to numpy
    //   K_FLD::IntArray cnto;
    //   neigh_phs.export_to_array(cnto);
    //   tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
    //   PyList_Append(l, tpl);
    //   Py_DECREF(tpl);
    // }

    // // zone 6 : complementary of above selections
    // if (remaining_phs.PHs.size())
    // {
    //   K_FLD::FloatArray crdtmp(crd);
    //   ngon_type::compact_to_used_nodes(remaining_phs.PGs, crdtmp); //reduce points
    //   //export to numpy
    //   K_FLD::IntArray cnto;
    //   remaining_phs.export_to_array(cnto);
    //   tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
    //   PyList_Append(l, tpl);
    //   Py_DECREF(tpl);
    // }
  }
  
  
  delete f; delete cn;
  return l;
}


//=============================================================================
/* Creates 2 zones : 1) outerlayer with firt neighborhoo 2) complementary */
//=============================================================================
PyObject* K_INTERSECTOR::extractOuterLayers(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int N(1), discard_external(0);

  if (!PYPARSETUPLEI(args, "Oll", "Oii", &arr, &N, &discard_external)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), outer, remaining;
 
  err = ngon_type::extract_n_outer_layers(crd, ngi, N, outer, remaining, discard_external);
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "extract_outer_layers failed.");
    delete f; delete cn;
    return NULL;
  } 
  
  PyObject *l(PyList_New(0)), *tpl;

  {
    //std::cout << "there are " << outer.PHs.size() << " outer cells detected." << std::endl;
    // zone 1 : outer
    if (outer.PHs.size())
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(outer.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      outer.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    // zone 2 : complementary 
    if (remaining.PHs.size())
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(remaining.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      remaining.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }
  
  
  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractNthCell(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &nth)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (nth >= ngi.PHs.size())
  {
    std::cout << "ERROR " << nth << " : is out of range" << std::endl;
    return NULL;
  }
 
  
  ngon_unit ph;
  ph.add(ngi.PHs.stride(nth), ngi.PHs.get_facets_ptr(nth));

  PyObject *l(PyList_New(0)), *tpl;
  
  // Extract the cell
  {
    ngon_type one_ph(ngi.PGs, ph);
    std::vector<E_Int> pgnids, phnids;
    one_ph.remove_unreferenced_pgs(pgnids, phnids);
    K_FLD::FloatArray crdtmp(crd);//we need crd for the neighbor zone
    ngon_type::compact_to_used_nodes(one_ph.PGs, crdtmp);

    K_FLD::IntArray cnto;
    one_ph.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  // Exract also its neighbors
  {
    std::vector<bool> keepPG(ngi.PGs.size(), false), keepPH;
    E_Int nb_pgs = ngi.PHs.stride(nth);
    E_Int* pgs = ngi.PHs.get_facets_ptr(nth);
    for (E_Int p=0; p < nb_pgs; ++p) keepPG[*(pgs+p)-1] = true;
    ngi.flag_PHs_having_PGs(keepPG, keepPH);

    keepPH[nth] = false;
    for (size_t i=0; i < keepPH.size(); ++i)
      if (keepPH[i]) std::cout << "neighbor : " << i << std::endl;

    ngon_type ngo;
    std::vector<E_Int> pgnids, phnids;
    ngon_type::select_phs(ngi, keepPH, pgnids, ngo);

    ngo.remove_unreferenced_pgs(pgnids, phnids);
    ngon_type::compact_to_used_nodes(ngo.PGs, crd);

    for (E_Int i=0; i < ngo.PHs.size(); ++i)
    {
      ngon_type one_ph;
      one_ph.addPH(ngo, i, true);

      K_FLD::IntArray cnto;
      one_ph.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }

  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::removeNthCell(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PyArg_ParseTuple(args, "Ol", &arr, &nth)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (nth >= ngi.PHs.size())
  {
    std::cout << "ERROR " << nth << " : is out of range" << std::endl;
    return NULL;
  }
  
  ngon_unit phs;
  for (E_Int i = 0; i < ngi.PHs.size(); ++i)
  {
    if (i == nth) continue;
    phs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
  }
  
  ngon_type ng(ngi.PGs, phs);

  std::vector<E_Int> pgnids, phnids;
  ng.remove_unreferenced_pgs(pgnids, phnids);
  ngon_type::compact_to_used_nodes(ng.PGs, crd);
  
  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractNthFace(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &nth)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (nth >= ngi.PGs.size())
  {
    std::cout << "ERROR " << nth << " : is out of range" << std::endl;
    return NULL;
  }

  PyObject *l(PyList_New(0)), *tpl;

  // Extract the face
  {
    ngon_unit pg;
    pg.add(ngi.PGs.stride(nth), ngi.PGs.get_facets_ptr(nth));

    ngon_type one_ph(pg, true);

    K_FLD::FloatArray crdtmp(crd); //we need crd for the parent elements zone
    ngon_type::compact_to_used_nodes(one_ph.PGs, crdtmp);

    K_FLD::IntArray cnto;
    one_ph.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  // Extract the Elements sharing this face
  {
    std::vector<bool> keepPG(ngi.PGs.size(), false);
    keepPG[nth]=true;
    std::vector<bool> keepPH;
    ngi.flag_PHs_having_PGs(keepPG, keepPH);

    for (size_t i=0; i < keepPH.size(); ++i)
      if (keepPH[i]) std::cout << "sharing element : " << i << std::endl;

    ngon_type ngo;
    std::vector<E_Int> pgnids, phnids;
    ngon_type::select_phs(ngi, keepPH, pgnids, ngo);

    ngo.remove_unreferenced_pgs(pgnids, phnids);
    ngon_type::compact_to_used_nodes(ngo.PGs, crd);

    for (E_Int i=0; i < ngo.PHs.size(); ++i)
    {
      ngon_type one_ph;
      one_ph.addPH(ngo, i, true);

      K_FLD::IntArray cnto;
      one_ph.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }

  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::checkCellsClosure(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  //K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  err = ngon_type::check_phs_closure(ngi);

  delete f; delete cn;

#ifdef E_DOUBLEINT
    return Py_BuildValue("l", long(err));
#else
    return Py_BuildValue("i", err);
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::removeBaffles(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  E_Int nb_baffles = ngi.remove_baffles();
  std::cout << "number of baffles removed : " << nb_baffles << std::endl;

  K_FLD::IntArray cnto;
  ngi.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::checkForDegenCells(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  //K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  E_Int degen = 0, imin(0), imax(0), mins(4), maxs(0);
  for (E_Int i=0; i< ngi.PHs.size(); ++i)
  {
    E_Int s = ngi.PHs.stride(i);
    if (s < mins)
    {
      mins=s;
      imin=i;
    }
    if (maxs < s)
    {
      maxs=s;
      imax=i;
    }
    if (s < 4) ++degen;
  }

  if (degen == 0) std::cout << "OK : There are no cell with less than 4 faces." << std::endl;
  else
    {
      std::cout << "ERROR : There are " << degen << " cells with less than 4 faces." << std::endl;
      std::cout << "the min of " << mins << " is reached first at : " << imin << std::endl;
      std::cout << "the max of " << maxs << " is reached first at : " << imax << std::endl;
    }

  delete f; delete cn;

#ifdef E_DOUBLEINT
    return Py_BuildValue("l", long(err));
#else
    return Py_BuildValue("i", err);
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::edgeLengthExtrema(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);
  E_Float Lmin,Lmax;
  E_Int imin, imax;

  ngon_type::edge_length_extrema(ngi.PGs, crd, Lmin, imin, Lmax, imax);

  std::cout << "Minimum Edge Length : " << Lmin << " reached at PG : " << imin << std::endl;
  std::cout << "Maximum Edge Length : " << Lmax << " reached at PG : " << imax << std::endl;

  delete f; delete cn;

#ifdef E_DOUBLEREAL
    return Py_BuildValue("d", double(Lmin));
#else
    return Py_BuildValue("f", float(Lmin));
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::computeAspectRatio(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Float vmin(0.);

  if (!PYPARSETUPLEF(args, "Od", "Of", &arr, &vmin)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  ngon_unit neighborsi;
  ngi.build_ph_neighborhood(neighborsi);

  Vector_t<E_Float> aspect_ratio;
  ngon_type::stats_bad_volumes<DELAUNAY::Triangulator>(crd, ngi, neighborsi, vmin, aspect_ratio);
  
  size_t sz = aspect_ratio.size();
  FloatArray ar(1, sz);
  for (size_t i = 0; i < sz; ++i) ar[i] = aspect_ratio[i];

  PyObject* tpl = K_ARRAY::buildArray(ar, "aspect_ratio", *cn, -1, "NGON", true);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeUserDefinedBC(PyObject* self, PyObject* args)
{
  PyObject *arr, *pgs;
  E_Float height(0.25);
  E_Int   strategy(0);   // 0 : CST_ABS , 1 : CST_REL_MEAN, 2 : CST_REL_MIN, 3 : VAR_REL_MEAN, 4 : VAR_REL_MIN
  E_Int create_ghost(true);
 
  if (!PYPARSETUPLE(args, "OOdll", "OOdii", "OOfll", "OOfii", &arr, &pgs, &height, &strategy, &create_ghost)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  std::vector<E_Int> PGlist;
  // Passing the specified wall pgs to the boolean to ignore cells that fall inside bodies
  {

    FldArrayI* inds=NULL;
    E_Int res=0;
    if (pgs != Py_None)
      res = K_NUMPY::getFromNumpyArray(pgs, inds, true);

    std::unique_ptr<FldArrayI> pL(inds); // to avoid to call explicit delete at several places in the code.
  
    //std::cout << "result for NUMPY is : " << res << std::endl;
    if ((res == 1) && (inds != NULL)  && (inds->getSize() != 0))
    {
      size_t nb_special_pgs = (size_t)inds->getSize();
      //E_Int minid(INT_MAX), maxid(-1);
      PGlist.resize(nb_special_pgs);
      for (size_t i = 0; i < nb_special_pgs; ++i) 
      {
        PGlist[i]=(*inds)[i]-1;
        //std::cout << pgsList[i] << std::endl;
        //minid = std::min(minid, pgsList[i]);
        //maxid = std::max(maxid, pgsList[i]);
      }

      //std::cout << "min/max : " << minid << "/" << maxid << std::endl;
    }
  }

  ngi.flag_externals(INITIAL_SKIN);
  bool has_been_reversed;
  DELAUNAY::Triangulator dt;
  err = ngon_type::reorient_skins(dt, crd, ngi, has_been_reversed);
  //std::cout << "reversed ? " << has_been_reversed << std::endl;
  if (!err)
  {
    ngon_type::eExtrudeStrategy strat = (ngon_type::eExtrudeStrategy)strategy;
    err = ngon_type::extrude_faces(crd, ngi, PGlist, height, bool(create_ghost), strat);
    //std::cout << "extrude_faces status : " << err << std::endl;
  }

  PyObject* tpl = NULL;

  if (!err)
  {
    K_FLD::IntArray cnto;
    ngi.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }
  
  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::reorientExternalFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngio(cnt), outer, remaining;

  ngio.flag_externals(INITIAL_SKIN);

  bool has_been_reversed;
  DELAUNAY::Triangulator t;
  err = ngon_type::reorient_skins(t, crd, ngio, has_been_reversed);

  // if (has_been_reversed)
  //   std::cout << "reorientExternalFaces : external faces has been reoriented" << std::endl;
  // else
  //   std::cout << "reorientExternalFaces : external faces orientation is correct" << std::endl;

  PyObject* tpl = NULL;

  if (!err)
  {
    K_FLD::IntArray cnto;
    ngio.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }
  
  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* reorient specified polygons. */
//=============================================================================
PyObject* K_INTERSECTOR::reorientSpecifiedFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;
  E_Int dir(1); //1 : outward -1 : inward

  if (!PYPARSETUPLEI(args, "OOl", "OOi", &arr, &py_pgs, &dir)) return NULL;

  if (dir != -1 && dir != 1) dir = 1;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  std::cout << "before numpy" << std::endl;

  E_Int res=0;
  E_Int* pgsList=NULL;
  E_Int size, nfld;
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, true/*shared*/, 0);

  std::cout << "after numpy" << std::endl;

  std::cout << res << std::endl;

  if (res != 1) return NULL;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngio(cnt);

  std::cout << "after ngio construct" << std::endl;

  std::vector<E_Int> plist, oids;
  plist.insert(plist.end(), pgsList, pgsList+size);

  std::cout << "after insert" << std::endl;
  //K_CONNECT::IdTool::shift(plist, -1);

  std::cout << "min pg specified : " << *std::min_element(pgsList, pgsList+size) << std::endl;
  std::cout << "max pg specified : " << *std::max_element(pgsList, pgsList+size) << std::endl;

  ngon_unit pgs;
  ngio.PGs.extract(plist, pgs, oids);

  std::vector<E_Int> orient;
  ngon_type::reorient_connex_PGs(pgs, (dir==-1), orient);

  // replace reverted polygons
  E_Int count(0);
  E_Int nb_pgs = pgs.size();
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    if (orient[i] == 1) continue;

    ++count;
    
    E_Int PGi = oids[i];

    E_Int s = ngio.PGs.stride(PGi);
    E_Int* p = ngio.PGs.get_facets_ptr(PGi);
    std::reverse(p, p + s);
  }

  std::cout << "nb of reoriented : "  << count  << " over " << size << " in pglist"<< std::endl;
    
  K_FLD::IntArray cnto;
  ngio.export_to_array(cnto);
  
  // pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);   
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::diffMesh(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;

  if (!PyArg_ParseTuple(args, "OO", &arr1, &arr2)) return NULL;

  K_FLD::FloatArray *f1(0), *f2(0);
  K_FLD::IntArray *cn1(0), *cn2(0);
  char *varString1, *varString2, *eltType1, *eltType2;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f1, cn1, varString1, eltType1);
  if (err) return NULL;
  // Check array # 2
  err = check_is_NGON(arr2, f2, cn2, varString2, eltType2);
  if (err) return NULL;

  std::unique_ptr<K_FLD::FloatArray> pf1(f1), pf2(f2);   //for memory cleaning
  std::unique_ptr<K_FLD::IntArray> pcn1(cn1), pcn2(cn2); //for memory cleaning

  K_FLD::FloatArray & crd = *f1;
  K_FLD::IntArray & cnt = *cn1;
  K_FLD::FloatArray & crd2 = *f2;
  K_FLD::IntArray & cnt2 = *cn2;

  // std::cout << "crd1 : " << crd.cols() << "/" << crd.rows() << std::endl;
  // std::cout << "cnt1 : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  // std::cout << "crd2 : " << crd2.cols() << "/" << crd2.rows() << std::endl;
  // std::cout << "cnt2 : " << cnt2.cols() << "/" << cnt2.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ng(cnt), ng1(cnt), ng2(cnt2);

  size_t nb_cells1 = (size_t)ng.PHs.size();
  size_t nb_cells2 = (size_t)ng2.PHs.size();

  // concatenate the meshes
  E_Int shft = crd.cols();
  ng2.PGs.shift(shft);
  ng.append(ng2);
  crd.pushBack(crd2);

  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);

  // detect identical cells
  std::vector<E_Int> nids; // cell ids in the concatenated ng
  E_Int nb_match = ng.detect_phs_with_same_centroid<acrd_t> (acrd, nids);

  //std::cout << "detect_phs_with_same_centroid : " << nb_match << std::endl;
  //
  if (nb_match == 0)
    std::cout << "the meshes are totally unmatching" << std::endl;

  std::vector<bool> keep1(nb_cells1, true);
  std::vector<bool> keep2(nb_cells2, true);
  //std::vector<bool> *pKa, *pKb;

  //
  size_t sz = nids.size();
  for (size_t i=0; i < sz; ++i)
  {
    size_t nid = nids[i];
    if (nid == i) continue;

    // 2 cells are matching
    if (i < nb_cells1)
      keep1[i] = false;
    else
      keep2[i-nb_cells1]=false;

    if (nid < nb_cells1)
      keep1[nid] = false;
    else
      keep2[nid-nb_cells1]=false;
  }

  PyObject *l(PyList_New(0)), *tpl;

  {
    ngon_type ng1o;
    std::vector<E_Int> pgnids, phnids;
    ng1.select_phs(ng1, keep1, phnids, ng1o);
    ng1o.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdtmp(crd);
    ngon_type::compact_to_used_nodes(ng1o.PGs, crdtmp); //reduce points
    //export to numpy
    K_FLD::IntArray cnto;
    ng1o.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString1, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  {
    ngon_type ng2o;
    std::vector<E_Int> pgnids, phnids;
    ng2.select_phs(ng2, keep2, phnids, ng2o);
    ng2o.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdtmp(crd);
    ngon_type::compact_to_used_nodes(ng2o.PGs, crdtmp); //reduce points
    //export to numpy
    K_FLD::IntArray cnto;
    ng2o.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString1, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::statsUncomputableFaces(PyObject* self, PyObject* args)
{

  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  std::vector<ngon_type::ePathoPG> flags;
  ngon_type::detect_uncomputable_pgs<DELAUNAY::Triangulator>(crd, ngi.PGs, flags);
  
  delete f; delete cn;

#ifdef E_DOUBLEINT
  return Py_BuildValue("l", long(0));
#else
  return Py_BuildValue("i", 0);
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::statsSize(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int comp_metrics(1);

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &comp_metrics)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  // Span of the mesh
  // Create the box
  K_SEARCH::BBox3D box;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd(crd);
  box.compute(acrd);
  // Box center and deltas
  E_Float dMax=0.;
  for (int i=0; i < 3; ++i)
   dMax = std::max(dMax, box.maxB[i] - box.minB[i]);
  //std::cout << "the span is : " << dMax << std::endl;

  E_Float smin(0.), smax(0.);
  E_Float vmin, vmax(-1.);

  if (comp_metrics == 1)
  {
    //
    E_Int imin, imax;
	  
	  ngon_type::surface_extrema(ngi.PGs, crd, smin, imin, smax, imax);
	  //std::cout << "the " << imin << "-th face has the smallest surface : " << smin << std::endl;
	  //std::cout << "the " << imax << "-th face has the biggest surface : " << smax << std::endl;
	  //
	  
	  ngon_type::volume_extrema<DELAUNAY::Triangulator>(ngi, crd, vmin, imin, vmax, imax);
	  //std::cout << "the " << imin << "-th cells has the smallest volume : " << vmin << std::endl;
    // if not a single cell
	  //if (imax != E_IDX_NONE) std::cout << "the " << imax << "-th cells has the biggest volume : " << vmax << std::endl;
  }

  if (ngi.PGs.size() == 1) smax = -1.;
  if (ngi.PHs.size() == 1) vmax = -1.;

  PyObject *l(PyList_New(0)), *tpl;
  
  delete f; delete cn;

#ifdef E_DOUBLEINT
  tpl = Py_BuildValue("d", long(dMax));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", dMax);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(smin));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", smin);
  PyList_Append(l, tpl);
#endif
#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(smax));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", smax);
  PyList_Append(l, tpl);
#endif
#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(vmin));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", vmin);
  PyList_Append(l, tpl);
#endif
#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(vmax));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", vmax);
  PyList_Append(l, tpl);
#endif

  return l;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::convert2Polyhedron(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);
  K_FLD::IntArray cnto;
  ngi.export_surfacic_view(cnto);
  ngon_unit pgs(cnto.begin());
  ngon_type ngo(pgs, true);

  cnto.clear();
  ngo.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::oneZonePerCell(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  PyObject *l(PyList_New(0)), *tpl;

  E_Int nb_phs = ngi.PHs.size();
  std::cout << "nb phs : " << nb_phs << std::endl;
  for (E_Int i=0; i < nb_phs; ++i)
  {
    ngon_type ngo;
    ngo.addPH(ngi, i, true);

    K_FLD::FloatArray crdtmp(crd);
    ngon_type::compact_to_used_nodes(ngo.PGs, crdtmp); //reduce points
    //export to numpy
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::closeOctalCells(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  ngon_type::close_phs(ngi, crd);

  K_FLD::IntArray cnto;
  ngi.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Converts a surfacic NGON from Cassiopee format to nuga format*/
//=============================================================================
PyObject* K_INTERSECTOR::convertNGON2DToNGON3D(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), outer, remaining;

  ngi.export_surfacic_view(cnt);

  ngon_unit pgs(&cnt[0]);

  PyObject* tpl = NULL;

  if (cnt.cols() != 0)
  {
  	ngon_type ngo(pgs, true);
  	K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }

  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* remove any cell contributing to a non-manifold boundary */
//=============================================================================
PyObject* K_INTERSECTOR::removeNonManifoldExternalCells(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  //std::cout << "nb initial phs : " << ngi.PHs.size() << std::endl;

  ngi.flag_external_pgs(INITIAL_SKIN);
  Vector_t<E_Int> oids;
  ngon_unit pg_ext;
  ngi.PGs.extract_of_type (INITIAL_SKIN, pg_ext, oids);

  //std::cout << "nb of external pgs : " << pg_ext.size() << std::endl;
    
  std::set<E_Int> pgids;
  ngon_type::get_pgs_with_non_manifold_edges(pg_ext, pgids);

  if (pgids.empty()) std::cout << "removeNonManifoldExternalCells : surface is clean" << std::endl;

  //for (auto i=pgids.begin(); i != pgids.end(); ++i) std::cout << "loc non manif id : " << *i << std::endl;

  std::set<E_Int> npgids;
  for (auto &i : pgids) npgids.insert(oids[i-1]);

  //for (auto i=npgids.begin(); i != npgids.end(); ++i) std::cout << "glob non manif id : " << *i << std::endl;
    
  std::set<E_Int> PHlist;
  ngi.get_PHs_having_PGs(npgids, 0, PHlist);

  //for (auto i=PHlist.begin(); i != PHlist.end(); ++i) std::cout << "PH id : " << *i << std::endl;

  ngi.remove_phs(PHlist);

  //std::cout << "nb final phs : " << ngi.PHs.size() << std::endl;

  K_FLD::IntArray cnto;
  ngi.export_to_array(cnto);

  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Computes centroids*/
//=============================================================================
PyObject* K_INTERSECTOR::centroids(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  K_FLD::FloatArray centroids;
  ngon_type::centroids<DELAUNAY::Triangulator>(ngi, crd, centroids);

  K_FLD::IntArray cnto;

  PyObject* tpl = K_ARRAY::buildArray(centroids, varString, cnto, 0, "NODE", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* retrieves any polygon that are overlapping */
//=============================================================================
PyObject* K_INTERSECTOR::getOverlappingFaces(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  E_Float RTOL(0.1), PS_MIN(0.95);
  E_Float d2[3], *pdir2(nullptr);

  if (!PYPARSETUPLEF(args, "OOdd(ddd)", "OOff(fff)", &arr1, &arr2, &RTOL, &PS_MIN, &d2[0], &d2[1], &d2[2])) return NULL;

  if (d2[0] != 0. || d2[1] != 0. || d2[2] != 0.)
    pdir2 = d2;

  K_FLD::FloatArray *f1(0), *f2(0);
  K_FLD::IntArray *cn1(0), *cn2(0);
  char *varString1, *varString2, *eltType1, *eltType2;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f1, cn1, varString1, eltType1);
  if (err) return NULL;

  // Check array # 2
  err = check_is_NGON(arr2, f2, cn2, varString2, eltType2);
  if (err) return NULL;

  std::unique_ptr<K_FLD::FloatArray> pf1(f1), pf2(f2);   //for memory cleaning
  std::unique_ptr<K_FLD::IntArray> pcn1(cn1), pcn2(cn2); //for memory cleaning

  K_FLD::FloatArray & crd1 = *f1;
  K_FLD::IntArray & cnt1 = *cn1;
  K_FLD::FloatArray & crd2 = *f2;
  K_FLD::IntArray & cnt2 = *cn2;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ng1(cnt1), ng2(cnt2);

  std::vector<E_Int> isx1, isx2;
    
  using tree_t = K_SEARCH::BbTree3D;
  using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
  using acnt_t = K_FLD::ArrayAccessor<K_FLD::IntArray>;
  
  tree_t tree(crd2, ng2.PGs);

  // compute a tolerance for the Tree
  E_Float tol,Lmax;
  E_Int imin, imax;
  ngon_type::edge_length_extrema(ng1.PGs, crd1, tol, imin, Lmax, imax);

  tol *= RTOL;
    
  NUGA::localizer<tree_t, acrd_t, acnt_t> localiz(tree, tol);
    
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
  NUGA::COLLIDE::compute_overlap<K_MESH::Polygon, K_MESH::Polygon>(crd1, ng1.PGs, crd2, ng2.PGs, localiz, PS_MIN, isx1, isx2, RTOL, true/*shuffle triangulation*/, pdir2);

#ifdef FLAG_STEP
  std::cout << "v0 : " << c.elapsed() << std::endl;
  std::cout << "nb x : " << std::count(isx1.begin(), isx1.end(), 1) << std::endl;
#endif

  PyObject *l(PyList_New(0)), *tpl;

  std::vector<E_Int> pgids1, pgids2;

  for (size_t i=0; i < isx1.size(); ++i)
    if (isx1[i]) pgids1.push_back(i);
  for (size_t i=0; i < isx2.size(); ++i)
    if (isx2[i]) pgids2.push_back(i);
    
  tpl = K_NUMPY::buildNumpyArray(&pgids1[0], pgids1.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  tpl = K_NUMPY::buildNumpyArray(&pgids2[0], pgids2.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  return l;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::merge(PyObject* self, PyObject* args)
{

  PyObject *arr1, *arr2;
  E_Float tolerance(1.e-15);
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &arr1, &arr2, &tolerance)) return NULL;

  char *varString, *eltType;
  E_Int ni, nj, nk;
  K_FLD::FloatArray *f1(0), *f2(0);
  K_FLD::IntArray *cn1(0);
  
  E_Int res = K_ARRAY::getFromArray(arr1, varString, f1, ni, nj, nk, cn1, eltType);

  // if (strcmp(eltType, "NODE") != 0)
  // {
  //   PyErr_SetString(PyExc_TypeError, "input error : invalid array, must be a NODE array.");
  //   delete f1; delete cn1;
  //   return nullptr;
  // }

  res = K_ARRAY::getFromArray(arr2, varString, f2, ni, nj, nk, cn1, eltType);

  // if (strcmp(eltType, "NODE") != 0)
  // {
  //   PyErr_SetString(PyExc_TypeError, "input error : invalid array, must be a NODE array.");
  //   delete f1, delete f2; delete cn1;
  //   return nullptr;
  // }
 
  K_FLD::FloatArray crd = *f1;
  crd.pushBack(*f2);

  K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(crd);
  Vector_t<E_Int> nids;
  E_Int nb_merges = ::merge(ca, tolerance, nids);

  E_Int sz = f2->cols();
  Vector_t<E_Int> nids_for_2(sz);
  for (E_Int i=0; i < sz; ++i)
  {
    nids_for_2[i] = nids[i+sz];
  }

  PyObject*tpl = K_NUMPY::buildNumpyArray(&nids_for_2[0], sz, 1, 0);
  delete f1, delete f2; delete cn1;

  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::oneph(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  char *varString, *eltType;
  E_Int ni, nj, nk;
  K_FLD::FloatArray *f(0);
  K_FLD::IntArray *cn(0);
  
  E_Int res = K_ARRAY::getFromArray(arr, varString, f, ni, nj, nk, cn, eltType);
  if ( (strcmp(eltType, "TRI") != 0) && (strcmp(eltType, "QUAD") != 0) )
  {
    PyErr_SetString(PyExc_TypeError, "input error : invalid array, must be a TRI or QUAD array.");
    delete f; delete cn;
    return nullptr;
  }

  //K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  K_FLD::FloatArray& crd = *f;

  ngon_unit pgs;
  ngon_unit::convert_fixed_stride_to_ngon_unit(cnt, 1, pgs);
 
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ng(pgs, 1);
  //ngon_type::clean_connectivity(ng, crd);

  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=======================  Intersector/PolyMeshTools/utils.cpp ====================
