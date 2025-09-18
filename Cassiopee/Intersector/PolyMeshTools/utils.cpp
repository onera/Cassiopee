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

# include "intersector.h"
# include <string>
# include <sstream> 

//#define FLAG_STEP
#ifdef FLAG_STEP
#include "chrono.h"
E_Int chrono::verbose=1;
#endif

//#include "Nuga/include/medit.hxx"
//std::string medith::wdir = "./";

# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"
#include "Nuga/include/localizer.hxx"
#include "Nuga/include/collider.hxx"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/displacement.hxx"
#include "Nuga/include/estimator.hxx"
#include "Nuga/include/Prism.h"
#include "Nuga/include/Hexahedron.h"
//#include <iostream>
#include <memory>
#include "dico_to_stl.h"
#ifdef _MPI
// close_cells use mpi calls
#include "Nuga/include/close_cells.hxx"
#endif

using namespace std;
using namespace K_FLD;
using namespace NUGA;

//=============================================================================
/* Creates 4 zones : 1) uncomputable polygons 2) uncomputable polyhedra 
   3) uncomputable polyhedra & neighbors 4) complementary of 3) */
//=============================================================================
PyObject* K_INTERSECTOR::extractUncomputables(PyObject* self, PyObject* args)
{
  E_Int neigh_level(1);
  PyObject *arr;

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &neigh_level)) return NULL;

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

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &neigh_level)) return NULL;

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

  if (!PYPARSETUPLE_(args, O_ II_, &arr, &N, &discard_external)) return NULL;

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

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &nth)) return NULL;

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
  // Exract also its shell
  {
    std::vector<bool> wprocessed/*externalized to not reallocate it each time*/;
    std::vector<E_Int> shellPHs, boundPGs;
    std::vector<bool> keepPH;
    keepPH.resize(ngi.PHs.size(), false);

    ngon_unit neighborsi;
    ngi.build_ph_neighborhood(neighborsi);

    ngon_type::ph_shell(ngi, nth, neighborsi, shellPHs, boundPGs, wprocessed);


    // extract shell
    
    for (size_t u = 0; u < shellPHs.size(); ++u)keepPH[shellPHs[u]] = true;

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
/* extractBiggestCell */
//=============================================================================
PyObject* K_INTERSECTOR::extractBiggestCell(PyObject* self, PyObject* args)
{

  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd(crd);

  E_Int nth=-1;
  E_Float dm(0.);
  for (E_Int i=0; i < ngi.PHs.size(); ++i)
  {
    K_SEARCH::BBox3D box;
    K_MESH::Polyhedron<0> PH(&ngi.PGs, ngi.PHs.get_facets_ptr(i), ngi.PHs.stride(i));
    PH.bbox(acrd, box);

    E_Float dx0 = box.maxB[0]-box.minB[0];
    E_Float dx1 = box.maxB[1]-box.minB[1];
    E_Float dx2 = box.maxB[2]-box.minB[2];

    if ( (dm < dx0) || (dm < dx1) || (dm < dx2))
    {
      dm = std::max(dx0, std::max(dx1,dx2));
      nth = i;
    }

  }
 
  std::cout << "biggest cell is : " << nth << std::endl;
  
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

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &nth)) return NULL;

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
PyObject* K_INTERSECTOR::removeNthFace(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &nth)) return NULL;

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
  
  ngon_unit pgs;
  std::vector<E_Int> nids(ngi.PGs.size(), IDX_NONE);
  E_Int count=0;
  for (E_Int i = 0; i < ngi.PGs.size(); ++i)
  {
    if (i == nth) continue;
    pgs.add(ngi.PGs.stride(i), ngi.PGs.get_facets_ptr(i));
    nids[i]=count++;
  }

  ngon_type ng;
  ng.PGs = pgs;
  ng.PHs = ngi.PHs;
  std::vector<E_Int> nfids;

  ng.PHs.remove_facets(nids, nfids);
  ng.updateFacets();

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
PyObject* K_INTERSECTOR::detectIdenticalCells(PyObject* self, PyObject* args)
{

  //std::cout << "detectIdenticalCells : begin" << std::endl;

  PyObject *arr;
  E_Int clean(0);
  E_Float tol(1.e-15);

  if (!PYPARSETUPLE_(args, O_ R_ I_, &arr, &tol, &clean)) return NULL;

  //std::cout << "detectIdenticalCells : after parse" << std::endl;

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

  Vector_t<E_Int> nids;
  ngi.detect_phs_with_same_centroid (crd, nids);
  bool found=false;

  E_Int nb_phs = ngi.PHs.size();
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (nids[i] != i)
    {
      std::cout << "detectIdenticalCells : " << i << " is identical to " << nids[i] << std::endl;
      found=true;
    }
  }

  if (!found)
    std::cout << "detectIdenticalCells : OK. No duplicates found." << std::endl;

  if (!clean || !found)
  {
    delete f; delete cn;
    return arr;
  }

  //std::cout << "detectIdenticalCells : clean" << std::endl;

  ngon_unit phs;
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (nids[i] != i) continue;
    phs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
  }

  //std::cout << "detectIdenticalCells : create output" << std::endl;
  
  ngon_type ng(ngi.PGs, phs);

  std::vector<E_Int> pgnids, phnids;
  ng.remove_unreferenced_pgs(pgnids, phnids);
  ngon_type::compact_to_used_nodes(ng.PGs, crd);

  //std::cout << "detectIdenticalCells : build array" << std::endl;
  
  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);

  //std::cout << "detectIdenticalCells : end" << std::endl;
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::detectOverConnectedFaces(PyObject* self, PyObject* args)
{

  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  std::vector<E_Int> pgs_occurrences;
  ngi.get_facets_connexion_nb(pgs_occurrences);

  bool error = false;
  for (size_t i=0; i < pgs_occurrences.size(); ++i)
  {
    if (pgs_occurrences[i] > 2)
    {
      std::cout << "multi PG : " << i << std::endl;
      error = true;
    }
  }

  if (!error)
    std::cout << "OK : no multiple PGs" << std::endl;
  
  K_FLD::IntArray cnto;
  ngi.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::collapseSmallEdges(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float edge_ratio(-1.), Lmax(-1.);

  if (!PYPARSETUPLE_(args, O_ RR_, &arr, &edge_ratio, &Lmax)) return NULL;

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
  ngon_unit pgso;

  std::vector<E_Int> nids;
  bool carry_on{false};

  //
  do
  {
    //double npgs0 = ngi.PGs.size();
    //std::cout << "collapse iter : " << iter++ << std::endl;

    //std::cout << "collapse_micro_edge" << std::endl;
    ngi.collapse_micro_edge(crd, edge_ratio, Lmax, nids);

    //validate/invalidate moves by flux
    ngon_unit neighborsi;
    ngi.build_ph_neighborhood(neighborsi);
    std::vector<E_Int> PHlist;
    K_CONNECT::IdTool::init_inc(PHlist, ngi.PHs.size());
    //std::cout << "validate_moves_by_fluxes" << std::endl;
    E_Int nb_valid_moves = ngon_type::validate_moves_by_fluxes<DELAUNAY::Triangulator>(nids, crd, ngi, neighborsi, PHlist);
    //std::cout << "nb valid moves : " << nb_valid_moves << std::endl;
    carry_on=(nb_valid_moves > 0);
    //std::cout << "clean_connectivity" << std::endl;
    ngi.PGs.change_indices(nids);
    ngon_type::clean_connectivity(ngi, crd, 3/*ngon_dim*/, 0./*tolerance*/, true/*remove_dup_phs*/, false/*do_omp*/);
    //std::cout << "nb pgs : " << ngi.PGs.size() << std::endl;
    //std::cout << "nb pts : " << crd.cols() << std::endl;
    //std::cout << "nb phs : " << ngi.PHs.size() << std::endl;
    
    // //if (carry_on)
    // {
      
    //   carry_on = (ngi.PGs.size() < npgs0);
    //   if (carry_on)npgs0 = ngi.PGs.size();
    // }
    //std::cout << "carry on ? : " << carry_on << std::endl;

  } while (carry_on);  

  //std::cout << "collapseMicroRegions : build array" << std::endl;
  ngon_type ngo = ngi;//(pgso, true);
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);

  //std::cout << "collapseMicroRegions : end" << std::endl;
  
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

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &nth)) return NULL;

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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

PyObject* K_INTERSECTOR::checkCellsFlux(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;
  if (!PYPARSETUPLE_(args, OO_, &arr, &PE)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  K_NUMPY::getFromNumpyArray(PE, cFE);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElment ne correpsond pas au nb de pgs" << std::endl;
    RELEASESHAREDN(PE, cFE);
    delete f; delete cn;
    return nullptr;
  }

  std::vector<E_Int> orient;
  E_Int imax=-1, uimax=-1, imax2=-1;
  E_Float fluxmax = -1., ufluxmax=-1., fluxmax2=-1.;
  for (E_Int i=0; i < ngi.PHs.size(); ++i)
  {
    orient.clear();

    const E_Int* pF = ngi.PHs.get_facets_ptr(i);
    E_Int nbf = ngi.PHs.stride(i);
    orient.resize(nbf, 1);

    for (E_Int f = 0; f < nbf; ++f)
    {
      E_Int PGi = *(pF+f) - 1;
      //std::cout << "PGi bef wwong :" << PGi << std::endl;
      if ((*cFE)(PGi, 1) != i+1) orient[f] = -1;
      assert (((*cFE)(PGi, 1) == i+1) || ((*cFE)(PGi, 2) == i+1) );
    }

    //std::cout << "computing flux for PH : " << i << std::endl;
    K_MESH::Polyhedron<0> PH(ngi, i);
    E_Float flxVec[3];
    PH.flux(crd, &orient[0], flxVec);

    E_Float flux = ::sqrt(K_FUNC::sqrNorm<3>(flxVec));
    E_Float s = PH.surface(crd);

    if (flux > ufluxmax)
    {
      uimax = i;
      ufluxmax = flux;
    }

    flux /= s; // normalizing

    if (flux > fluxmax)
    {
      imax = i;
      fluxmax = flux;
    }

    double Lref = ::sqrt(PH.Lref2(crd, NUGA::ISO_MEAN));

    flux /= Lref; // normalizing #2

    if (flux > fluxmax2)
    {
      imax2 = i;
      fluxmax2 = flux;
    }


  }

  std::cout << "normalized max flux is : " << fluxmax << " reached at cell : " << imax << std::endl;
  std::cout << "unormalized max flux is : " << ufluxmax << " reached at cell : " << uimax << std::endl;
  std::cout << "normalized (vol) max flux is : " << fluxmax2 << " reached at cell : " << imax2 << std::endl;

  delete f; delete cn;

  PyObject *l(PyList_New(0)), *tpl;

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(imax));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", imax);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(fluxmax));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(fluxmax);
  PyList_Append(l, tpl);
#endif
  
  RELEASESHAREDN(PE, cFE);
  return l;
}

PyObject* K_INTERSECTOR::checkAngularExtrema(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;
  if (!PYPARSETUPLE_(args, OO_, &arr, &PE)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  K_NUMPY::getFromNumpyArray(PE, cFE);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElment ne correpsond pas au nb de pgs" << std::endl;
    delete f; delete cn;
    return nullptr;
  }

  std::vector<E_Int> orient;
  E_Int imax{-1}, imin{-1};
  E_Int PG1max = -1;
  E_Int PG2max = -1;
  E_Float minA{7.}, maxA{-1.};
  
  //
  for (E_Int i=0; i < ngi.PHs.size(); ++i)
  {
    orient.clear();

    const E_Int* pF = ngi.PHs.get_facets_ptr(i);
    E_Int nbf = ngi.PHs.stride(i);
    orient.resize(nbf, 1);

    for (E_Int f = 0; f < nbf; ++f)
    {
      E_Int PGi = *(pF+f) - 1;
      //std::cout << "PGi bef wwong :" << PGi << std::endl;
      if ((*cFE)(PGi, 1) != i+1) orient[f] = -1;
      assert (((*cFE)(PGi, 1) == i+1) || ((*cFE)(PGi, 2) == i+1) );
    }

    E_Float mA, MA;
    E_Int maxAPG1 = -1, maxAPG2 = -1;
    err = K_MESH::Polyhedron<UNKNOWN>::min_max_angles(crd, ngi.PGs, pF, nbf, false, &orient[0], mA, MA, maxAPG1, maxAPG2);

    if (MA > maxA)
    {
      imax = i;
      maxA = MA;
      PG1max = maxAPG1;
      PG2max = maxAPG2;
    }
    if (mA < minA)
    {
      imin = i;
      minA = mA;
    }
  }

  std::cout << "minimal dihedral angle: " << minA << " reached at cell: " << imin << std::endl;
  std::cout << "maximal dihedral angle: " << maxA << " reached at cell: " << imax << std::endl;
  std::cout << "maximal dihedral angle for PGs: " << PG1max << " and " << PG2max << std::endl;

  delete f; delete cn;

  PyObject *l(PyList_New(0)), *tpl;

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(imax));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", imax);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(maxA));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(maxA);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(imin));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", imin);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(minA));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(minA);
  PyList_Append(l, tpl);
#endif
  
  return l;
}

int comp_vol(const K_FLD::FloatArray& crd, const ngon_type& ngi, const FldArrayI* cFE, std::vector<E_Int>& orient, E_Int i, DELAUNAY::Triangulator & dt, double &v)
{
  //std::cout << "PH : " << i << std::endl;
  orient.clear();

  const E_Int* pF = ngi.PHs.get_facets_ptr(i);
  E_Int nbf = ngi.PHs.stride(i);
  orient.resize(nbf, 1);

  for (E_Int j = 0; j < nbf; ++j)
  {
    E_Int PGi = *(pF+j) - 1;
    //std::cout << "PGi bef wwong :" << PGi << std::endl;
    if ((*cFE)(PGi, 1) != i+1) orient[j] = -1;
    //assert (((*cFE)(PGi, 1) == i+1) || ((*cFE)(PGi, 2) == i+1) );
  }

  //std::cout << "computing flux for PH : " << i << std::endl;
  K_MESH::Polyhedron<0> PH(ngi, i);
  
  E_Int err = PH.volume<DELAUNAY::Triangulator>(crd, &orient[0], v, dt);

  return err;
}

PyObject* K_INTERSECTOR::checkCellsVolume(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;
  if (!PYPARSETUPLE_(args, OO_, &arr, &PE)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  K_NUMPY::getFromNumpyArray(PE, cFE);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElment ne correpsond pas au nb de pgs" << std::endl;
    delete f; delete cn;
    return nullptr;
  }

  E_Int imin=-1;
  E_Float vmin = NUGA::FLOAT_MAX;
  E_Int imax=-1;
  E_Float vmax = 0;

  E_Int nb_max_threads = __NUMTHREADS__;
  //std::cout << "nb threads max : " << nb_max_threads << std::endl;
    
  std::vector<E_Int> im(nb_max_threads, IDX_NONE);
  std::vector<E_Float> vm(nb_max_threads, NUGA::FLOAT_MAX);
  std::vector<std::vector<E_Int>> orient(nb_max_threads);

  std::vector<E_Int> iM(nb_max_threads, IDX_NONE);
  std::vector<E_Float> vM(nb_max_threads, 0.);

  E_Int i, id{0};
  DELAUNAY::Triangulator dt;

#pragma omp parallel shared(vm, im, vM, iM, ngi, crd, cFE, orient) private (i, id, dt) default(none)
{
  id = __CURRENT_THREAD__;
  //std::cout << "before loop thread : " << id  << std::endl;
#pragma omp for //schedule(dynamic)
  for (i=0; i < ngi.PHs.size(); ++i)
  {
    double v;
    E_Int err = comp_vol(crd, ngi, cFE, orient[id], i, dt, v);
    if (!err && v < vm[id]) // min for current thread
    {
      im[id] = i;
      vm[id] = v;
    }
    if (!err && v > vM[id]) // min for current thread
    {
      iM[id] = i;
      vM[id] = v;
    }
  }
}

  for (E_Int i=0; i < nb_max_threads; ++i)
  {
    if (vm[i] < vmin)
    {
      imin = im[i];
      vmin = vm[i];
    }
    if (vM[i] > vmax)
    {
      imax = iM[i];
      vmax = vM[i];
    }
  }

  std::cout << "min vol is : " << vmin << " reached at cell : " << imin << std::endl;
  std::cout << "max vol is : " << vmax << " reached at cell : " << imax << std::endl;

  delete f; delete cn;

  PyObject *l(PyList_New(0)), *tpl;

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(imin));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", imin);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(vmin));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(vmin);
  PyList_Append(l, tpl);
#endif
  
  return l;
}

PyObject* K_INTERSECTOR::checkCellsVolumeAndGrowthRatio(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;
  if (!PYPARSETUPLE_(args, OO_, &arr, &PE)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  K_NUMPY::getFromNumpyArray(PE, cFE);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElement ne correpsond pas au nb de pgs" << std::endl;
    delete f; delete cn;
    return nullptr;
  }

  E_Int nphs = ngi.PHs.size();

  std::vector<double> vols(nphs, NUGA::FLOAT_MAX);

  E_Int nb_max_threads = __NUMTHREADS__;
  //std::cout << "nb threads max : " << nb_max_threads << std::endl;
    
  std::vector<std::vector<E_Int>> orient(nb_max_threads);

  E_Int i, id{0};
  DELAUNAY::Triangulator dt;

#pragma omp parallel shared(ngi, crd, cFE, orient, vols) private (i, id, dt) default(none)
  {
    id = __CURRENT_THREAD__;
    //std::cout << "before loop thread : " << id  << std::endl;
#pragma omp for //schedule(dynamic)
    for (i=0; i < ngi.PHs.size(); ++i)
    {
      double v;
      E_Int err = comp_vol(crd, ngi, cFE, orient[id], i, dt, v);
      if (!err)
        vols[i] = v;
    }
  }

  //
  ngon_unit neighborsi;
  ngi.build_ph_neighborhood(neighborsi);

  Vector_t<E_Float> growth_ratio;
  ngon_type::stats_bad_volumes<DELAUNAY::Triangulator>(crd, ngi, neighborsi, vols, -1., growth_ratio);

  E_Int ivolmin = -1,igrmin = -1;
  E_Float volmin{NUGA::FLOAT_MAX}, grmin{NUGA::FLOAT_MAX};

  for (size_t i=0; i < vols.size(); ++i)
  {
    if (vols[i] < volmin)
    {
      volmin = vols[i];
      ivolmin = i;
    }
  }

  for (size_t i=0; i < growth_ratio.size(); ++i)
  {
    if (growth_ratio[i] < grmin)
    {
      grmin = growth_ratio[i];
      igrmin = i;
    }
  }
 
  delete f; delete cn;
 
  PyObject *l(PyList_New(0)), *tpl;

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(ivolmin));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", ivolmin);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(volmin));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(volmin);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(igrmin));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", igrmin);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(grmin));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(grmin);
  PyList_Append(l, tpl);
#endif
  
  return l;
}

PyObject* K_INTERSECTOR::extractBadVolCells(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;
  double aratio{0.125}, vmin{0.};
  int nneighs{1};

  if (!PYPARSETUPLE_(args, OO_ RR_ I_, &arr, &PE, &aratio, &vmin, &nneighs)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  K_NUMPY::getFromNumpyArray(PE, cFE);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElement ne correpsond pas au nb de pgs" << std::endl;
    delete f; delete cn;
    return nullptr;
  }

  E_Int nphs = ngi.PHs.size();

  std::vector<E_Int> orient;
  std::vector<double> vols(nphs, NUGA::FLOAT_MAX);
  
  //compute volumes using input orientation 
  for (E_Int i=0; i < ngi.PHs.size(); ++i)
  {
    orient.clear();

    const E_Int* pF = ngi.PHs.get_facets_ptr(i);
    E_Int nbf = ngi.PHs.stride(i);
    orient.resize(nbf, 1);

    for (E_Int f = 0; f < nbf; ++f)
    {
      E_Int PGi = *(pF+f) - 1;
      //std::cout << "PGi bef wwong :" << PGi << std::endl;
      if ((*cFE)(PGi, 1) != i+1) orient[f] = -1;
      assert (((*cFE)(PGi, 1) == i+1) || ((*cFE)(PGi, 2) == i+1) );
    }

    //std::cout << "computing flux for PH : " << i << std::endl;
    K_MESH::Polyhedron<0> PH(ngi, i);
    double v;
    DELAUNAY::Triangulator dt;
    E_Int err = PH.volume<DELAUNAY::Triangulator>(crd, &orient[0], v, dt);

    if (!err)
      vols[i] = v;
    else
    {
      //std::cout << "error to triangulate cell " << i << "at face : " << err-1 << std::endl;
      //medith::write("badcell", crd, ngi, i);
      //medith::write("faultyPG", crd, ngi.PGs.get_facets_ptr(err-1), ngi.PGs.stride(err-1), 1);
    }
  }

  ngon_unit neighborsi;
  ngi.build_ph_neighborhood(neighborsi);

  Vector_t<E_Float> growth_ratio;
  ngon_type::stats_bad_volumes<DELAUNAY::Triangulator>(crd, ngi, neighborsi, vols, -1., growth_ratio);

  std::vector<bool> keep(nphs, false);
  E_Int badcount=0;
  for (size_t i=0; i < (size_t)nphs; ++i)
  {
    if ( (growth_ratio[i] < aratio) || (vols[i] < vmin) ) {
      ++badcount;
      keep[i]=true;
    }
  }

  //std::cout << "nb of bad cells found : " << badcount << " (over " << nphs << ")" << std::endl;

  // extend with second neighborhood and separate from non-involved polyhedra
  for (E_Int j=0; j< nneighs; ++j)
    ngon_type::flag_neighbors(ngi, keep);
    
  ngon_type ngo;  
  {
    Vector_t<E_Int> ngpids;
    ngi.select_phs(ngi, keep, ngpids, ngo);
  }

  ngon_type::compact_to_used_nodes(ngo.PGs, crd); //reduce points

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

PyObject* K_INTERSECTOR::extractOverConnectedCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  int nneighs{1};

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &nneighs)) return NULL;

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

  std::vector<E_Int> nfconnects;
  ngi.get_facets_connexion_nb(nfconnects);

  std::vector<bool> keep(ngi.PHs.size(), false);

  for (size_t k = 0; k < (size_t)ngi.PHs.size(); ++k)
  {
    size_t nf = ngi.PHs.stride(k);
    const E_Int* faces = ngi.PHs.get_facets_ptr(k);
    size_t nb_overcon = 0;
    for (size_t f = 0; f < nf; ++f)
    {
      E_Int PGi = faces[f] - 1;
      if (nfconnects[PGi] > 2) ++nb_overcon;
    }
    if (nb_overcon == nf) // => fully over connected
    {
      keep[k]= true;
    }
  }

  // extend with second neighborhood and separate from non-involved polyhedra
  for (E_Int j=0; j< nneighs; ++j)
    ngon_type::flag_neighbors(ngi, keep);
    
  ngon_type ngo;  
  {
    Vector_t<E_Int> ngpids;
    ngi.select_phs(ngi, keep, ngpids, ngo);
  }

  ngon_type::compact_to_used_nodes(ngo.PGs, crd); //reduce points

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

///
PyObject* K_INTERSECTOR::volume(PyObject* self, PyObject* args)
{
  PyObject *arr, *axcelln;
  if (!PYPARSETUPLE_(args, OO_, &arr, &axcelln)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (xcelln)
  bool use_xcelln = false;

  K_FLD::FloatArray* xcelln(nullptr);
  K_FLD::IntArray *cn1(0);
  char *varString1;
  E_Int ni, nj, nk;
  E_Int res = 0;
  if (axcelln != Py_None) res = K_ARRAY::getFromArray(axcelln, varString1, xcelln, ni, nj, nk, cn1, eltType);
  if (res == 1) use_xcelln = true;

  //std::cout << py_xcelln << std::endl;

  // E_Int res = 0;
  // if (py_xcelln != Py_None)
  // {
  //   std::cout << "get numpy " << std::endl;

  //   res = K_NUMPY::getFromNumpyArray(py_xcelln, xcelln);
  //   std::cout << xcelln << std::endl;
  //   if (res == 1) use_xcelln = true;
  // }

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  // std::cout << "use_xcelln ? " << use_xcelln << std::endl;
  // std::cout << "xcelln ? " << xcelln << std::endl;
  // std::cout << "res : " << res << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (use_xcelln && (xcelln != nullptr) && (ngi.PHs.size() != xcelln->getSize()))
  {
    std::cout << "le champ xcelln ne correpsond pas au nb de polyedres => pas pris en compte" << std::endl;
    std::cout << "nb phs : " << ngi.PHs.size() << std::endl;
    std::cout << "taille xcelln : " << xcelln->getSize() << std::endl;
    use_xcelln = false;
  }

  std::vector<E_Float> vols;
  ngon_type::volumes<DELAUNAY::Triangulator>(crd, ngi, vols, false/*not all cvx*/, true/*new algo*/);

  //std::cout << "use_xcelln ?" << use_xcelln << std::endl;
  E_Float V = 0.;
  if (use_xcelln)
  {
    for (size_t i = 0; i < vols.size(); ++i)
      V += vols[i] * (*xcelln)[i];
  }
  else
    for (size_t i = 0; i < vols.size(); ++i)
      V += vols[i];

  delete f; delete cn;

#ifdef E_DOUBLEREAL
  return Py_BuildValue("d", V);
#else
    return Py_BuildValue("f", float(V));
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::removeBaffles(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  E_Int nb_cells_w_baffles = ngi.remove_baffles();
  if (nb_cells_w_baffles) std::cout << "number of cells with baffles found : " << nb_cells_w_baffles << std::endl;

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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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
PyObject* K_INTERSECTOR::checkForBigCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int N{8};

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &N)) return NULL;

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
  E_Int maxf{0};
  E_Int idm{IDX_NONE};
  E_Int count{0};
  std::vector<bool> keep(ngi.PHs.size(), false);
  for (E_Int i=0; i< ngi.PHs.size(); ++i)
  {
    E_Int s = ngi.PHs.stride(i);
    if (s > N)
    {
      ++count;
      keep[i]=true;
    }
    if (s > maxf)
    {
      maxf = s;
      idm = i;
    }
    maxf = std::max(maxf, s);
    
  }

  PyObject* tpl{nullptr};

  if (idm != IDX_NONE) std::cout << idm << " is the biggest with " << maxf << " faces" << std::endl;

  if (count > 0)
  {
    std::cout << count << " cells over the specified number of faces have been found." << std::endl;
    
    std::vector<E_Int> pgnids, phnids;
    ngon_type ngo;
    ngi.select_phs(ngi, keep, phnids, ngo);
    ngo.remove_unreferenced_pgs(pgnids, phnids);
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }
  else
    std::cout << "No cells over the specified number of faces have been found." << std::endl;

  delete f; delete cn;
  return tpl;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::edgeLengthExtrema(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  //std::cout << "Minimum Edge Length : " << Lmin << " reached at PG : " << imin << std::endl;
  //std::cout << "Maximum Edge Length : " << Lmax << " reached at PG : " << imax << std::endl;

  delete f; delete cn;

#ifdef E_DOUBLEREAL
    return Py_BuildValue("d", double(Lmin));
#else
    return Py_BuildValue("f", float(Lmin));
#endif
}

PyObject* K_INTERSECTOR::edgeLengthMax(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  //std::cout << "Minimum Edge Length : " << Lmin << " reached at PG : " << imin << std::endl;
  //std::cout << "Maximum Edge Length : " << Lmax << " reached at PG : " << imax << std::endl;

  delete f; delete cn;

#ifdef E_DOUBLEREAL
    return Py_BuildValue("d", double(Lmax));
#else
    return Py_BuildValue("f", float(Lmax));
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::computeGrowthRatio(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Float vmin(0.);

  if (!PYPARSETUPLE_(args, O_ R_, &arr, &vmin)) return NULL;

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

  std::vector<E_Float> vols;
  ngon_type::volumes<DELAUNAY::Triangulator>(crd, ngi, vols, false/*not all cvx*/, false/* ! new algo*/);


  Vector_t<E_Float> growth_ratio;
  ngon_type::stats_bad_volumes<DELAUNAY::Triangulator>(crd, ngi, neighborsi, vols, vmin, growth_ratio);
  
  size_t sz = growth_ratio.size();
  FloatArray ar(1, sz);
  for (size_t i = 0; i < sz; ++i) ar[i] = growth_ratio[i];

  PyObject* tpl = K_ARRAY::buildArray(ar, "growth_ratio", *cn, -1, "NGON", true);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeBC(PyObject* self, PyObject* args)
{
  PyObject *arr, *pgs;
  E_Float height(0.25);
  E_Int   strategy(0);   // 0 : CST_ABS , 1 : CST_REL_MEAN, 2 : CST_REL_MIN, 3 : VAR_REL_MEAN, 4 : VAR_REL_MIN
  E_Int create_ghost(true);
 
  if (!PYPARSETUPLE_(args, OO_ R_ II_, &arr, &pgs, &height, &strategy, &create_ghost)) return NULL;

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
      res = K_NUMPY::getFromNumpyArray(pgs, inds);

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
/* XXX : encore un beau commentaire explicite de Sam!! */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeSurf(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float layer_height(0.25);
  E_Int   strategy(0), nlayers(1);   // 0 : CST_ABS , 1 : CST_REL_MEAN, 2 : CST_REL_MIN, 3 : VAR_REL_MEAN, 4 : VAR_REL_MIN
 
  if (!PYPARSETUPLE_(args, O_ R_ II_, &arr, &layer_height, &nlayers, &strategy)) return NULL;

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
  ngon_type ngio(cnt);

  ngio.PHs.clear();

  std::vector<E_Int> pglist;
  K_CONNECT::IdTool::init_inc(pglist, ngio.PGs.size());

  ngon_type::eExtrudeStrategy strat = (ngon_type::eExtrudeStrategy)strategy;

  std::vector<E_Int> tops;
  int smooth_iters = 0;

  for (int l = 0; l < nlayers; ++l)
  {
    tops.clear();
    ngon_type::extrude_faces(crd, ngio, pglist, layer_height, true, strat, smooth_iters, &tops);
    pglist.clear();
    pglist.insert(pglist.end(), tops.begin(), tops.end());
  }

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
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeRevolSurf(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float pt[3], dir[3];
  E_Int   nlayers(1);   // 0 : CST_ABS , 1 : CST_REL_MEAN, 2 : CST_REL_MIN, 3 : VAR_REL_MEAN, 4 : VAR_REL_MIN

  //std::cout << "extrudeRevolSurf : 1" << std::endl;
 
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ I_, &arr, &pt[0], &pt[1], &pt[2], &dir[0], &dir[1], &dir[2], &nlayers)) return NULL;

  //std::cout << "extrudeRevolSurf : 2" << std::endl;

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
  ngon_type ngio(cnt);

  ngio.PHs.clear();

  std::vector<E_Int> pglist;
  K_CONNECT::IdTool::init_inc(pglist, ngio.PGs.size());

  std::vector<E_Int> tops;

  E_Float angle = K_FUNC::normalize<3>(dir);

  //std::cout << "extrudeRevolSurf : 3 : angle : " << angle << std::endl;

  for (int l = 0; l < nlayers; ++l)
  {
    //std::cout << "extrudeRevolSurf : layer : " << l << std::endl;
    tops.clear();
    ngon_type::extrude_revol_faces(crd, ngio, pglist, dir, pt, angle, &tops);

    pglist.clear();
    pglist.insert(pglist.end(), tops.begin(), tops.end());
  }

  //std::cout << "extrudeRevolSurf : 4" << std::endl;

  PyObject* tpl = NULL;

  if (!err)
  {
    K_FLD::IntArray cnto;
    ngio.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }

  //std::cout << "extrudeRevolSurf : 5" << std::endl;
  
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

  if (!PYPARSETUPLE_(args, OO_ I_, &arr, &py_pgs, &dir)) return NULL;

  if (dir != -1 && dir != 1) dir = 1;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "before numpy" << std::endl;

  E_Int res=0;
  E_Int* pgsList=NULL;
  E_Int size, nfld;
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld);

  //std::cout << "after numpy" << std::endl;

  std::cout << res << std::endl;

  if (res != 1) return NULL;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngio(cnt);

  //std::cout << "after ngio construct" << std::endl;

  std::vector<E_Int> plist, oids;
  plist.insert(plist.end(), pgsList, pgsList+size);

  //std::cout << "after insert" << std::endl;
  //K_CONNECT::IdTool::shift(plist, -1);

  //std::cout << "min pg specified : " << *std::min_element(pgsList, pgsList+size) << std::endl;
  //std::cout << "max pg specified : " << *std::max_element(pgsList, pgsList+size) << std::endl;

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

  //std::cout << "nb of reoriented : "  << count  << " over " << size << " in pglist"<< std::endl;
    
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
PyObject* K_INTERSECTOR::reorient(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int dir(1); //1 : outward -1 : inward

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &dir)) return NULL;

  if (dir != -1 && dir != 1) dir = 1;

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
  
  ngon_type ngio(cnt);

  ngio.flag_externals(INITIAL_SKIN);
  bool has_been_reversed;
  DELAUNAY::Triangulator t;
  err = ngon_type::reorient_skins(t, crd, ngio, has_been_reversed);

  if (dir == -1)
  {
    Vector_t<E_Int> oids;
    ngon_unit pg_ext;
    ngio.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);

    E_Int nb_pgs = pg_ext.size();
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      E_Int PGi = oids[i];

      E_Int s = ngio.PGs.stride(PGi);
      E_Int* p = ngio.PGs.get_facets_ptr(PGi);
      std::reverse(p, p + s);
    }
  }
    
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
PyObject* K_INTERSECTOR::externalFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_skipids;
  E_Int geo_dim(-1);

  if (!PYPARSETUPLE_(args, OO_ I_, &arr, &py_skipids, &geo_dim)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return nullptr;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type::eGEODIM geodim = (geo_dim == -1) ? ngon_type::get_ngon_geodim(cnt) : (ngon_type::eGEODIM)geo_dim;

  //std::cout << "GEO dIM ? " << geodim << std::endl;

  if (geodim == ngon_type::eGEODIM::ERROR)
  {
    std::cout << "externalFaces : Input Error : mesh is corrupted." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::MIXED)
  {
    std::cout << "externalFaces : Input Error : mesh mixed elt types (lineic and/or surfacic and /or volumic." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::LINEIC)
  {
    std::cout << "externalFaces : Unsupported : lineic NGON are not handled." << std::endl;
    return nullptr;
  }

  // so either SURFACIC, SURFACIC_CASSIOPEE or VOLUMIC

  if (geodim == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnt1;
    ng.export_surfacic_view(cnt1);
    //std::cout << "exported" << std::endl;
    geodim = ngon_type::eGEODIM::SURFACIC;
    cnt=cnt1;
  }

  PyObject *l(PyList_New(0));
  
  // SURFACIC OR VOLUMIC ?

  if (geodim == ngon_type::eGEODIM::SURFACIC)
  {
    //std::cout << "mesh object..." << std::endl;
    NUGA::pg_smesh_t mesh(crd, cnt);
    //std::cout << "getting boundary..." << std::endl;
    NUGA::edge_mesh_t ef;
    mesh.get_boundary(ef);
    // pushing out the BAR
    //std::cout << "pushing out" << std::endl;
    PyObject *tpl = K_ARRAY::buildArray(ef.crd, varString, ef.cnt, -1, "BAR", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    PyList_Append(l, Py_None); //no edge history
  }
  else if (geodim == ngon_type::eGEODIM::VOLUMIC)
  {
    NUGA::ph_mesh_t mesh(crd, cnt);
    NUGA::pg_smesh_t ef;
    std::vector<E_Int> oids, ancestors;
    mesh.get_boundary(ef, oids, ancestors);

    //std::cout << "ef.ncells vs oids sz : " << ef.ncells() << " / " << oids.size() << std::endl;

    // get PG ids to skip 
    std::vector<bool> keep;
    std::set<E_Int> discard_ids;

    //std::cout << "py_skipids : " << py_skipids << std::endl;
    if (py_skipids != Py_None)
    {
      keep.resize(ef.ncells(), true);
      E_Int *ptL, size, nfld;
      K_NUMPY::getFromNumpyArray(py_skipids, ptL, size, nfld);

      discard_ids.insert(ptL, ptL+size);

      for (size_t k=0; k < (size_t)ef.ncells(); ++k)
        if (discard_ids.find(oids[k]) != discard_ids.end())
          keep[k] = false;

      // remove joins if any specified
      ef.compress(keep);
      K_CONNECT::IdTool::compact(oids, keep);
    }

    // pushing out the NGON SURFACE
    K_FLD::IntArray cnto;
    ngon_type ng(ef.cnt, true);
    ng.export_to_array(cnto);
    PyObject *tpl = K_ARRAY::buildArray(ef.crd, varString, cnto, -1, "NGON", false);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out history  
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  
  delete f; delete cn;
  return l;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::diffMesh(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  if (!PYPARSETUPLE_(args, OO_, &arr1, &arr2)) return NULL;

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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  if (!PYPARSETUPLE_(args, O_ I_, &arr, &comp_metrics)) return NULL;

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
	  std::cout << "the " << imin << "-th cells has the smallest volume : " << vmin << std::endl;
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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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
PyObject* K_INTERSECTOR::oneZonePerFace(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  E_Int nb_pgs = ngi.PGs.size();
  std::cout << "nb pgs : " << nb_pgs << std::endl;
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    ngon_unit ngu;
    ngu.add(ngi.PGs.stride(i), ngi.PGs.get_facets_ptr(i));
    
    ngon_type ngo(ngu, true);

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
PyObject* K_INTERSECTOR::immerseNodes(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  double TOL{0.};

  if (!PYPARSETUPLE_(args, OO_ R_, &arr1, &arr2, &TOL)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  K_FLD::FloatArray* sf(0);
  K_FLD::IntArray* scn(0);
  char* svarString, *seltType;
  // Check array # 2
  err = check_is_NGON(arr2, sf, scn, svarString, seltType);
  if (err) return NULL;

  K_FLD::FloatArray & surf_crd = *sf;
  K_FLD::IntArray & surf_cnt = *scn;
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  //////// IMMERSING NODES

  using zmesh_t = NUGA::ph_mesh_t;
  using bmesh_t = NUGA::pg_smesh_t;

  // mesh on which some nodes will move
  zmesh_t m1(crd, cnt);
  // surface to test
  bmesh_t s2(surf_crd, surf_cnt);

  // reorient outward
  //std::cout << " S2 ORIENT ? " << s2.oriented << std::endl;
  if (s2.oriented == 0)
  {
    ngon_type ng2(s2.cnt, true/*one ph for all*/);
    bool has_been_reversed;
    DELAUNAY::Triangulator dt;
    ngon_type::reorient_skins(dt, s2.crd, ng2, has_been_reversed);
    s2.cnt = ng2.PGs;
    s2.oriented = 1;
  }
  else if (s2.oriented == -1)
    s2.reverse_orient();

  auto opp_dirs = NUGA::immerse_nodes(m1, s2, TOL);

  // // grabb involved points
  // auto crd = m1.crd;
  // std::vector<bool> keep(crd.cols(), false);
  // for (size_t i = 0; i < opp_dirs.size(); ++i)
  //   keep[opp_dirs[i].flag] = true;

  // K_CONNECT::keep<bool> pred(keep);
  // K_CONNECT::IdTool::compress(crd, pred);

  ////////////////////////
  
  K_FLD::IntArray cnto;
  m1.cnt.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(m1.crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::closeCells(PyObject* self, PyObject* args)
{
  PyObject *py_arrs(nullptr);
  PyObject *py_zids(nullptr), *py_zid_to_rid_to_list(nullptr);
  PyObject *py_rid_to_zones(nullptr);
  if (!PYPARSETUPLE_(args, OOOO_, &py_arrs, &py_zids, &py_zid_to_rid_to_list, &py_rid_to_zones)) return NULL;

  // 1. GET MESHES 
  E_Int nb_meshes{0};
  if (PyList_Check(py_arrs))
    nb_meshes = PyList_Size(py_arrs);

  //std::cout << "nb_meshes : " << nb_meshes << std::endl;

  if (nb_meshes == 0) return nullptr;

  std::vector<K_FLD::FloatArray*> f(nb_meshes);
  std::vector<K_FLD::IntArray*> cn(nb_meshes);
  char* varString, *eltType;

  std::vector<NUGA::ph_mesh_t> ph_meshes(nb_meshes);
  std::vector<NUGA::ph_mesh_t*> ptr_ph_meshes(nb_meshes);

  using ngon_type = ngon_t<K_FLD::IntArray>;
  
  for (E_Int m = 0; m < nb_meshes; ++m)
  {
    PyObject* arr = PyList_GetItem(py_arrs, m);
    // Check array
    E_Int err = check_is_NGON(arr, f[m], cn[m], varString, eltType);
    if (err) return nullptr;
    
    // conversion to the generic mesh interface
    ph_meshes[m].crd = *(f[m]);
    ph_meshes[m].cnt = std::move(ngon_type(*(cn[m])));

    ptr_ph_meshes[m] = &ph_meshes[m];
  }

  // 2. GET ZIDS 
  std::vector<int> zids(nb_meshes);
  assert (nb_meshes == PyList_Size(py_zids));
  for (int m = 0; m < nb_meshes; ++m)
  {
    PyObject* pyz = PyList_GetItem(py_zids, m);
    E_Int zid = PyInt_AsLong(pyz);
    zids[m] = zid;
  }

  // 3. GET  GET RID TO LIST MAP  MAP 
  std::map<int, std::map<int, std::vector<E_Int>>> zid_to_rid_to_list;
  convert_dico_to_map___int_int_vecint(py_zid_to_rid_to_list, zid_to_rid_to_list);
  //assert (zid_to_rid_to_list.size() == nb_meshes);

  // 4. GET RID_TO_ZONES MAP 
  //todo VD : py_rid_to_zones => rid_to_zones
  std::map<int, std::pair<int,int>> rid_to_zones;
  convert_dico_to_map__int_pairint(py_rid_to_zones, rid_to_zones);

  // 5. CLOSE
#ifdef _MPI
  // CB:  close cells is removed from compilation without mpi
  using para_algo_t = NUGA::omp_algo<NUGA::ph_mesh_t, E_Float>; // SEQ or multi-SEQ
  using closecell_t = NUGA::close_cells< para_algo_t, NUGA::ph_mesh_t>;

  closecell_t cc;
  cc.run(ptr_ph_meshes, zids, zid_to_rid_to_list, rid_to_zones);
#endif

  // pushing out the result : the set of closed meshes
  PyObject *l(PyList_New(0));
  for (size_t i=0; i < (size_t)nb_meshes; ++i)
  {
    K_FLD::IntArray cnto;
    ptr_ph_meshes[i]->cnt.export_to_array(cnto);

    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(ptr_ph_meshes[i]->crd, varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  for (int m = 0; m < nb_meshes; ++m)
  {
    delete f[m]; delete cn[m];
  }

  return l;
}

//=============================================================================
/* Converts a surfacic NGON from Cassiopee format to nuga format*/
//=============================================================================
PyObject* K_INTERSECTOR::convertNGON2DToNGON3D(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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
/* Converts a surfacic Basic to NGON nuga format*/
//=============================================================================
PyObject* K_INTERSECTOR::convertBasic2NGONFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_BASICF(arr, f, cn, varString, eltType);
  if (err)
  {
    std::cout << "convertBasic2NGONFaces : ERROR : " << err << std::endl;
    return nullptr;
  }

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  ngon_unit pgs;
  ngon_unit::convert_fixed_stride_to_ngon_unit(cnt, 1, pgs);

  PyObject* tpl = NULL;

  if (pgs.size() != 0)
  {
    ngon_type wNG(pgs, true);
    K_FLD::IntArray cnto;
    wNG.export_to_array(cnto);
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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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
/* Computes centroids */
//=============================================================================
PyObject* K_INTERSECTOR::centroids(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type::eGEODIM geodim = ngon_type::get_ngon_geodim(cnt);

  //std::cout << "GEO dIM ? " << geodim << std::endl;

  if (geodim == ngon_type::eGEODIM::ERROR)
  {
    std::cout << "centroids : Input Error : mesh is corrupted." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::MIXED)
  {
    std::cout << "centroids : Input Error : mesh mixed elt types (lineic and/or surfacic and /or volumic." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::LINEIC)
  {
    std::cout << "centroids : Unsupported : lineic NGON are not handled." << std::endl;
    return nullptr;
  }

  // so either SURFACIC, SURFACIC_CASSIOPEE or VOLUMIC

  if (geodim == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnt1;
    ng.export_surfacic_view(cnt1);
    //std::cout << "exported" << std::endl;
    geodim = ngon_type::eGEODIM::SURFACIC;
    cnt=cnt1;
  }

  ngon_type ngi(cnt);

  K_FLD::FloatArray cents;

  if (geodim == ngon_type::eGEODIM::SURFACIC)
    ngon_type::centroids(ngi.PGs, crd, cents);
  else // (geodim == ngon_type::eGEODIM::VOLUMIC)
    ngon_type::centroids<DELAUNAY::Triangulator>(ngi, crd, cents);

  K_FLD::IntArray cnto;

  PyObject* tpl = K_ARRAY::buildArray(cents, varString, cnto, 0, "NODE", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Computes volumes */
//=============================================================================
PyObject* K_INTERSECTOR::volumes(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int algo(0);
  E_Int all_pgs_cvx(0);

  if (!PYPARSETUPLE_(args, O_ II_, &arr, &algo, &all_pgs_cvx)) return NULL;

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

  std::vector<E_Float> vols;
  ngon_type::volumes<DELAUNAY::Triangulator>(crd, ngi, vols, (all_pgs_cvx == 1), (algo == 1));

  size_t sz = vols.size();
  FloatArray ar(1, sz);
  for (size_t i = 0; i < sz; ++i) 
    ar[i] = vols[i];

  PyObject* tpl = K_ARRAY::buildArray(ar, "volumes", *cn, -1, "NGON", true);

  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* retrieves any polygon that are overlapping */
//=============================================================================
PyObject* K_INTERSECTOR::getOverlappingFaces(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  E_Float RTOL(0.1), AMAX(0.1);
  E_Float d2[3], *pdir2(nullptr);

  if (!PYPARSETUPLE_(args, OO_ RR_ TRRR_, &arr1, &arr2, &RTOL, &AMAX, &d2[0], &d2[1], &d2[2])) return NULL;

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
    
  tree_t tree(crd2, ng2.PGs);

  // compute a tolerance for the Tree
  E_Float tol,Lmax;
  E_Int imin, imax;
  ngon_type::edge_length_extrema(ng1.PGs, crd1, tol, imin, Lmax, imax);

  tol *= RTOL;
    
  NUGA::localizer<tree_t> localiz(tree, tol);
    
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif

  AMAX = std::max(0., AMAX);
  AMAX = std::min(K_CONST::E_PI, AMAX);
  double PS_MIN = ::cos(AMAX);

  NUGA::COLLIDE::compute_overlap<K_MESH::Polygon, K_MESH::Polygon>(crd1, ng1.PGs, crd2, ng2.PGs, localiz, isx1, isx2, RTOL, PS_MIN, true/*shuffle triangulation*/, pdir2);

#ifdef FLAG_STEP
  std::cout << "v0 : " << c.elapsed() << std::endl;
  std::cout << "nb x : " << std::count(isx1.begin(), isx1.end(), 1) << std::endl;
#endif

  PyObject *l(PyList_New(0)), *tpl;

  std::vector<E_Int> pgids1, pgids2;

  for (size_t i=0; i < isx1.size(); ++i)
    if (isx1[i] != IDX_NONE) pgids1.push_back(i);
  for (size_t i=0; i < isx2.size(); ++i)
    if (isx2[i] != IDX_NONE) pgids2.push_back(i);
    
  tpl = K_NUMPY::buildNumpyArray(&pgids1[0], pgids1.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  tpl = K_NUMPY::buildNumpyArray(&pgids2[0], pgids2.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  return l;

}

//=============================================================================

//=============================================================================
PyObject* K_INTERSECTOR::getCollidingTopFaces(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  E_Float RTOL(1.e-12);

  if (!PYPARSETUPLE_(args, OO_ I_, &arr1, &arr2, &RTOL)) return NULL;

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

  ngon_type::eGEODIM dim1 = ngon_type::get_ngon_geodim(cnt1);
  ngon_type::eGEODIM dim2 = ngon_type::get_ngon_geodim(cnt2);
 
  if (dim1 != ngon_type::eGEODIM::VOLUMIC)
  {
    std::cout << "getCollidingTopFaces : INPUT ERROR : input mesh colidee mesh must be volumic" << std::endl;
    return nullptr;
  }
  
  if (dim2 == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt2);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnttmp;
    ng.export_surfacic_view(cnttmp);
    //std::cout << "exported" << std::endl;
    dim2 = ngon_type::eGEODIM::SURFACIC;
    cnt2=cnttmp;
  }

  if (dim2 != ngon_type::eGEODIM::SURFACIC)
  {
    std::cout << "getCollidingTopFaces : INPUT ERROR : input mesh colider mesh must be surfacic" << std::endl;
    return nullptr;
  }

  ph_mesh_t  m1(crd1, cnt1);
  pg_smesh_t m2(crd2, cnt2);
  
  std::vector<bool> isx1, isx2;
  NUGA::COLLIDE::compute(m1,m2, isx1, isx2, RTOL);

  auto neighbors = m1.get_neighbors();

  std::vector<E_Int> freezePG(m1.cnt.PGs.size(), false);
  for (size_t i = 0; i < isx1.size(); ++i)
  {
    if (!isx1[i]) continue;

    bool is_hexa = K_MESH::Hexahedron::is_of_type(m1.cnt.PGs, m1.cnt.PHs.get_facets_ptr(i), m1.cnt.PHs.stride(i));
    if (!is_hexa) continue;

    size_t nneighs = neighbors->stride(i);

    const E_Int* pneighs = neighbors->get_facets_ptr(i);

    const E_Int* pfaces = m1.cnt.PHs.get_facets_ptr(i);

    //freeze top of boundary elements
    for (size_t k = 0; k< nneighs; ++k)
    {
      if (pneighs[k] != IDX_NONE) continue;

      int top = K_MESH::Hexahedron::get_opposite(m1.cnt.PGs, pfaces, k);
      assert (top != IDX_NONE);

      freezePG[pfaces[top]-1]=true;
      //std::cout << "freezing " << pfaces[top]-1 << std::endl;
    }    
  }

  std::set<E_Int> pgids;
  for (size_t i = 0; i < isx1.size(); ++i)
  {
    if (!isx1[i]) continue;

    bool is_prism{false}, is_hexa{false};
    is_hexa = K_MESH::Hexahedron::is_of_type(m1.cnt.PGs, m1.cnt.PHs.get_facets_ptr(i), m1.cnt.PHs.stride(i));
    if (!is_hexa)
      is_prism = K_MESH::Prism::is_of_type(m1.cnt.PGs, m1.cnt.PHs.get_facets_ptr(i), m1.cnt.PHs.stride(i));

    if (!is_prism && !is_hexa) continue;

    size_t nneighs = neighbors->stride(i);

    const E_Int* pneighs = neighbors->get_facets_ptr(i);

    const E_Int* pfaces = m1.cnt.PHs.get_facets_ptr(i);

    E_Int top=IDX_NONE, topOK=IDX_NONE;

    // seek for colliging regions boundary
    for (size_t k = 0; k< nneighs; ++k)
    {
      if (pneighs[k] == IDX_NONE) continue;  // boundary
      if (pneighs[k] != IDX_NONE && isx1[pneighs[k]]) continue; // inner face

      E_Int PGi = pfaces[k]-1;

      if (freezePG[PGi]) continue;

      if (is_prism && m1.cnt.PGs.stride(PGi) != 3) // bot must be triangle
        continue;

      //
      if (is_hexa)
        top = K_MESH::Hexahedron::get_opposite(m1.cnt.PGs, pfaces, k);
      else // prism
        top = K_MESH::Prism::get_opposite(m1.cnt.PGs, pfaces, k);

      if (top == IDX_NONE) continue;
      E_Int PGtop = pfaces[top]-1;
      if (freezePG[PGtop]) continue;

      if (pneighs[top] == IDX_NONE) continue;
      if (!isx1[pneighs[top]]) continue;

      if (topOK != IDX_NONE)
      {
        topOK = IDX_NONE;
        break;
      }

      //bot = k;
      topOK=top;
    }

    if (topOK == IDX_NONE) continue;

    pgids.insert(pfaces[topOK]-1);
  }

  //std::cout << "nb of tops : " << pgids.size() << std::endl;
    
  PyObject* tpl = nullptr;
  if (!pgids.empty())
  {
    std::vector<E_Int> tmp(ALL(pgids));
    tpl = K_NUMPY::buildNumpyArray(&tmp[0], tmp.size(), 1, 0);
  }
  return tpl;

}

//=============================================================================
/* retrieves any cells that are colliding */
//=============================================================================
PyObject* K_INTERSECTOR::getCollidingCells(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  E_Float RTOL(1.e-12);
  E_Int only_externals{0};

  if (!PYPARSETUPLE_(args, OO_ R_ I_, &arr1, &arr2, &RTOL, &only_externals)) return NULL;

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

  ngon_type::eGEODIM dim1 = ngon_type::get_ngon_geodim(cnt1);
  ngon_type::eGEODIM dim2 = ngon_type::get_ngon_geodim(cnt2);
 
  if (dim1 == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt1);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnttmp;
    ng.export_surfacic_view(cnttmp);
    //std::cout << "exported" << std::endl;
    dim1 = ngon_type::eGEODIM::SURFACIC;
    cnt1=cnttmp;
  }

  std::vector<bool> isx1, isx2;
  
  if (dim2 == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt2);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnttmp;
    ng.export_surfacic_view(cnttmp);
    //std::cout << "exported" << std::endl;
    dim2 = ngon_type::eGEODIM::SURFACIC;
    cnt2=cnttmp;
  }

  if (dim1 == ngon_type::eGEODIM::SURFACIC && dim2 == ngon_type::eGEODIM::LINEIC) // S vs BAR
  {
    pg_smesh_t m1(crd1, cnt1);
    edge_mesh_t m2(crd2, cnt2);
    NUGA::COLLIDE::compute(m1,m2, isx1, isx2, RTOL); 
  }
  /*else if (dim1 == ngon_type::eGEODIM::SURFACIC && dim2 == ngon_type::eGEODIM::SURFACIC) // S vs S
  {
    pg_smesh_t m1(crd1, cnt1);
    pg_smesh_t m2(crd2, cnt2);
    NUGA::COLLIDE::compute(m1,m2, isx1, isx2, RTOL); 
  }*/
  else if (dim1 == ngon_type::eGEODIM::VOLUMIC && dim2 == ngon_type::eGEODIM::SURFACIC) // V vs S
  {
    ph_mesh_t m1(crd1, cnt1);
    pg_smesh_t m2(crd2, cnt2);
    NUGA::COLLIDE::compute(m1,m2, isx1, isx2, RTOL); 
  }
  else if (dim1 == ngon_type::eGEODIM::VOLUMIC && dim2 == ngon_type::eGEODIM::VOLUMIC) // V vs V (take only PGs into account => pg_smesh_t)
  {
    ph_mesh_t m1(crd1, cnt1);
    pg_smesh_t m2(crd2, cnt2);
    NUGA::COLLIDE::compute(m1,m2, isx1, isx2, RTOL);

    if (only_externals)
    {
      m1.cnt.flag_externals(1);
      ngon_type::discard_holes_by_box(m1.crd, m1.cnt);

      for (size_t i=0; i < isx1.size(); ++i)
      {
       if (!isx1[i]) continue;
       if (m1.cnt.PHs._type[i] != 1) isx1[i] = 0;
      }
    }
  }
  else
  {
    std::cout << "getCollidingCells : INPUT ERROR : input mesh dimensiosn combination is not handled yet." << std::endl;
    return nullptr;
  }

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
/* Computes cell sensor data from the metric in a mesh */
//=============================================================================
PyObject* K_INTERSECTOR::estimateAdapReq(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  E_Float RTOL(1.e-12);
  E_Int minv{0}, maxv{5}, metric_policy{2}; // 0:MIN ; 1:MEAN ; 2:MAX

  if (!PYPARSETUPLE_(args, OO_ I_ R_ II_, &arr1, &arr2, &metric_policy, &RTOL, &minv, &maxv)) return NULL;

  K_FLD::FloatArray *f1(0), *f2(0);
  K_FLD::IntArray *cn1(0), *cn2(0);
  char *varString1, *varString2, *eltType1, *eltType2;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f1, cn1, varString1, eltType1);
  if (err) return NULL;

  // Check array # 2
  err = check_is_NGON(arr2, f2, cn2, varString2, eltType2);
  if (err) return NULL;

  std::cout << "metric_policy : " << metric_policy << std::endl;
  NUGA::eMetricType M = (NUGA::eMetricType)metric_policy;
 //std::cout << "M : " << M << std::endl;

  std::unique_ptr<K_FLD::FloatArray> pf1(f1), pf2(f2);   //for memory cleaning
  std::unique_ptr<K_FLD::IntArray> pcn1(cn1), pcn2(cn2); //for memory cleaning

  K_FLD::FloatArray & crd1 = *f1;
  K_FLD::IntArray & cnt1 = *cn1;
  K_FLD::FloatArray & crd2 = *f2;
  K_FLD::IntArray & cnt2 = *cn2;

  ngon_type::eGEODIM dim1 = ngon_type::get_ngon_geodim(cnt1);
  ngon_type::eGEODIM dim2 = ngon_type::get_ngon_geodim(cnt2);
 
  if (dim1 == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt1);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnttmp;
    ng.export_surfacic_view(cnttmp);
    //std::cout << "exported" << std::endl;
    dim1 = ngon_type::eGEODIM::SURFACIC;
    cnt1=cnttmp;
  }

  std::vector<E_Int> data;
  
  if (dim2 == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt2);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnttmp;
    ng.export_surfacic_view(cnttmp);
    //std::cout << "exported" << std::endl;
    dim2 = ngon_type::eGEODIM::SURFACIC;
    cnt2=cnttmp;
  }

  if (dim1 == ngon_type::eGEODIM::SURFACIC && dim2 == ngon_type::eGEODIM::LINEIC) // S vs BAR
  {
    pg_smesh_t m1(crd1, cnt1);
    edge_mesh_t m2(crd2, cnt2);
    NUGA::estimate_adap_req(m1, m2, M, RTOL, data, minv, maxv);
  }
  else if (dim1 == ngon_type::eGEODIM::SURFACIC && dim2 == ngon_type::eGEODIM::SURFACIC) // S vs S
  {
    pg_smesh_t m1(crd1, cnt1);
    edge_mesh_t m2(crd2, cnt2);
    NUGA::estimate_adap_req(m1, m2, M, RTOL, data, minv, maxv);
  }
  else if (dim1 == ngon_type::eGEODIM::VOLUMIC && dim2 == ngon_type::eGEODIM::SURFACIC) // V vs S
  {
    ph_mesh_t m1(crd1, cnt1);
    pg_smesh_t m2(crd2, cnt2);
    NUGA::estimate_adap_req(m1, m2, M, RTOL, data, minv, maxv);
  }
  else if (dim1 == ngon_type::eGEODIM::VOLUMIC && dim2 == ngon_type::eGEODIM::VOLUMIC) // V vs V (take only PGs into account => pg_smesh_t)
  {
    ph_mesh_t m1(crd1, cnt1);
    pg_smesh_t m2(crd2, cnt2);
    NUGA::estimate_adap_req(m1, m2, M, RTOL, data, minv, maxv);
  }
  else
  {
    std::cout << "estimateAdapReq : INPUT ERROR : input mesh dimensiosn combination is not handled yet." << std::endl;
    return nullptr;
  }

    
  PyObject* tpl = K_NUMPY::buildNumpyArray(&data[0], data.size(), 1, 0);
  return tpl;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getNthNeighborhood(PyObject* self, PyObject* args)
{
  PyObject *arr1, *py_ids;
  E_Int N{0};

  if (!PYPARSETUPLE_(args, O_ I_ O_, &arr1, &N, &py_ids)) return NULL;

  K_FLD::FloatArray *f1(0);
  K_FLD::IntArray *cn1(0);
  char *varString1, *eltType1;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f1, cn1, varString1, eltType1);
  if (err) return NULL;

  // Check numpy array
  FldArrayI* fids;
  K_NUMPY::getFromNumpyArray(py_ids, fids);

  std::vector<E_Int> ids(fids->getSize());
  for (size_t i=0; i < ids.size();++i)
    ids[i] = (*fids)[i];
  delete fids;

  std::unique_ptr<K_FLD::FloatArray> pf1(f1);   //for memory cleaning
  std::unique_ptr<K_FLD::IntArray> pcn1(cn1); //for memory cleaning

  K_FLD::FloatArray & crd1 = *f1;
  K_FLD::IntArray & cnt1 = *cn1;
  

  ngon_type::eGEODIM dim1 = ngon_type::get_ngon_geodim(cnt1);
  
  if (dim1 == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt1);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnttmp;
    ng.export_surfacic_view(cnttmp);
    //std::cout << "exported" << std::endl;
    dim1 = ngon_type::eGEODIM::SURFACIC;
    cnt1=cnttmp;
  }

  std::vector<E_Int> neighs;

  /*if (dim1 == ngon_type::eGEODIM::SURFACIC )
  {
    pg_smesh_t m1(crd1, cnt1);
    m1.get_nth_neighborhood(N, ids, neighs);
  }
  else*/ if (dim1 == ngon_type::eGEODIM::VOLUMIC)
  {
    ph_mesh_t m1(crd1, cnt1);
    m1.get_nth_neighborhood(N, ids, neighs);
  }
  else
  {
    std::cout << "getNthNeighborhood : INPUT ERROR : input mesh dimension combination is not handled yet." << std::endl;
    return nullptr;
  }

  PyObject *tpl = K_NUMPY::buildNumpyArray(&neighs[0], neighs.size(), 1, 0);

  return tpl;

}

//=============================================================================
/* retrieves any polygon that are connecting 2 aniso HEXA */
//=============================================================================
PyObject* K_INTERSECTOR::getAnisoInnerFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float aniso_ratio(0.05);
  
  if (!PYPARSETUPLE_(args, O_ R_, &arr, &aniso_ratio)) return NULL;

  K_FLD::FloatArray *f(0);
  K_FLD::IntArray *cn(0);
  char *varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  std::unique_ptr<K_FLD::FloatArray> pf(f);   //for memory cleaning
  std::unique_ptr<K_FLD::IntArray> pcn(cn); //for memory cleaning

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ng(cnt);

  E_Int nb_pgs = ng.PGs.size();
  std::vector<E_Int> flag(nb_pgs, 0);
  E_Int nb_phs = ng.PHs.size();
  E_Int bot, top;
  for (E_Int i=0; i < nb_phs; ++i)
  {
    const E_Int* faces = ng.PHs.get_facets_ptr(i);
    if (K_MESH::Polyhedron<0>::is_aniso_HX8(crd, ng.PGs, faces, ng.PHs.stride(i), aniso_ratio, bot, top))
    {
      ++flag[*(faces+bot)-1];
      ++flag[*(faces+top)-1];
    }
  }
  std::vector<E_Int> pgids;
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    if (flag[i] == 2) //inner
      pgids.push_back(i);
  }

  std::cout << "getAnisoInnerFaces : " << pgids.size() << std::endl;

  PyObject* tpl = K_NUMPY::buildNumpyArray(&pgids[0], pgids.size(), 1, 0);

  return tpl;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::merge(PyObject* self, PyObject* args)
{

  PyObject *arr1, *arr2;
  E_Float tolerance(1.e-15);
  if (!PYPARSETUPLE_(args, OO_ R_, &arr1, &arr2, &tolerance)) return NULL;

  char *varString, *eltType;
  E_Int ni, nj, nk;
  K_FLD::FloatArray *f1(0), *f2(0);
  K_FLD::IntArray *cn1(0);
  
  /*E_Int res = */K_ARRAY::getFromArray(arr1, varString, f1, ni, nj, nk, cn1, eltType);

  // if (strcmp(eltType, "NODE") != 0)
  // {
  //   PyErr_SetString(PyExc_TypeError, "input error : invalid array, must be a NODE array.");
  //   delete f1; delete cn1;
  //   return nullptr;
  // }

  /*res = */K_ARRAY::getFromArray(arr2, varString, f2, ni, nj, nk, cn1, eltType);

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
  /*E_Int nb_merges = */::merge(ca, tolerance, nids, true /*do_omp*/);

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
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

  char *varString, *eltType;
  E_Int ni, nj, nk;
  K_FLD::FloatArray *f(0);
  K_FLD::IntArray *cn(0);
  
  /*E_Int res = */K_ARRAY::getFromArray(arr, varString, f, ni, nj, nk, cn, eltType);
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

  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::concatenate(PyObject* self, PyObject* args)
{

  E_Float tol(1.e-15);
  PyObject* arrs;

  if (!PYPARSETUPLE_(args, O_ R_, &arrs, &tol))
  {
    PyErr_SetString(PyExc_TypeError, "concatenate : wrong args");
    return NULL;
  }

  E_Int nb_zones = PyList_Size(arrs);
  

  std::vector<K_FLD::FloatArray*> crds(nb_zones, nullptr);
  std::vector<K_FLD::IntArray*>   cnts(nb_zones, nullptr);
  char* varString, *eltType;

  ngon_type::eGEODIM geodim = ngon_type::eGEODIM::UNSET;
  E_Int err = 0;

  // get the zones
  for (E_Int i=0; (i < nb_zones) && !err; ++i)
  {
    //std::cout << "getting zone in list : " << i << std::endl;
    PyObject* py_zone = PyList_GetItem(arrs, i);
    
    err = check_is_NGON(py_zone, crds[i], cnts[i], varString, eltType);
    if (err) break;

    ngon_type::eGEODIM gd = ngon_type::get_ngon_geodim(*(cnts[i]));
    //std::cout << "gd : " << gd << std::endl;

    if (geodim == ngon_type::eGEODIM::UNSET)
      geodim = gd;

    if (geodim != gd)
    {
      //
      PyErr_SetString(PyExc_TypeError, "concatenate : mixed type inputs not handled");
      err = 1; break;
    }
    else if (geodim == ngon_type::eGEODIM::ERROR)
    {
      //
      std::ostringstream o;
      o << "concatenate : input << " << i << " is not handled";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1; break;
    }
    // std::cout << "zone sizes : " << crds[i]->cols() << " points" << std::endl;
    // std::cout << "zone sizes : " << cnts[i]->cols() << " cells" << std::endl;
  }

  if (err)
  {
    for (E_Int i=0; i < nb_zones; ++i)
    {
      delete crds[i];
      delete cnts[i];
    }
    return NULL;
  }

  // join and close


  // Creation z_pgnids et z_phnids par zone
  // z_pgnids[i][k] = j <=> la keme face de la zone i a pour indice j dans ng (concatenation) 
  std::vector<std::vector<E_Int>> z_pgnids(nb_zones) ;
  std::vector<std::vector<E_Int>> z_phnids(nb_zones) ;
  
  ngon_type ng;
  K_FLD::FloatArray crd; 
  std::vector<E_Float> Lmin2, *pLmin2(nullptr);
  bool tol_is_relative = (tol < 0.);
  
  for (size_t i=0; i < cnts.size(); ++i)
  {
    //std::cout << "appending" << std::endl;
    ngon_type ngt(*cnts[i]);
    ngt.PGs.shift(crd.cols());
    
    K_CONNECT::IdTool::init_inc(z_pgnids[i], ngt.PGs.size()); 
    K_CONNECT::IdTool::init_inc(z_phnids[i], ngt.PHs.size());

    // Shift pour indir. dans tab. globaux  
    for ( size_t j=0; j<(size_t)ngt.PGs.size(); ++j){ z_pgnids[i][j] +=  ng.PGs.size(); } 
    for ( size_t j=0; j<(size_t)ngt.PHs.size(); ++j){ z_phnids[i][j] +=  ng.PHs.size(); }

    ng.append(ngt);
    crd.pushBack(*crds[i]);

    if (tol_is_relative)
    {
      // compute nodal metric for this zone
      NUGA::ph_mesh_t mesh(*crds[i], *cnts[i]);
      auto lm2 = mesh.get_nodal_metric2(); // importtant to be ISO_MIN to detect closest nodes in warpes meshes

      Lmin2.insert(Lmin2.end(), ALL(lm2));
    }
  }

  if (tol_is_relative) pLmin2 = &Lmin2;

  // Creation z_pgnids et z_phnids global
  std::vector<E_Int> glo_pgnids;
  std::vector<E_Int> glo_phnids;

  // std::cout << "before clean : nb_phs/phs/crd : " << ng.PHs.size() << "/" << ng.PGs.size() << "/" << crd.cols() << std::endl;

  ngon_type::clean_connectivity(ng, crd, -1/*ngon_dim*/, tol/*tolerance*/, true/*remove_dup_phs*/, true/*do_omp*/, &glo_pgnids, &glo_phnids, pLmin2 );

  //Propagation des nouveaux ids dans z_pgnids/z_phnids
  for (size_t i=0; i < (size_t)nb_zones; ++i)
  {
    // propagate pgnids
    for ( size_t j=0; j< z_pgnids[i].size(); ++j) 
    {
      E_Int iids     = z_pgnids[i][j];
      if ( glo_pgnids[iids] == E_IDX_NONE ){ glo_pgnids[iids] = -1 ;} //convert IDX_NONE value into -1
      z_pgnids[i][j] = glo_pgnids[iids];
    }
    // propagate phnids
    for ( size_t j=0; j< z_phnids[i].size(); ++j) 
    {
      E_Int iids     = z_phnids[i][j];
      if ( glo_phnids[iids] == E_IDX_NONE ) { glo_phnids[iids] = -1 ;} //convert IDX_NONE value into -1
      z_phnids[i][j] = glo_phnids[iids];
    }
  }

  if (geodim == ngon_type::eGEODIM::SURFACIC)
  //if (geodim == SURFACIC) // NUGA SURF
    ng = ngon_type(ng.PGs, true);

  //std::cout << "after clean : nb_phs/phs/crd : " << ng.PHs.size() << "/" << ng.PGs.size() << "/" << crd.cols() << std::endl;
 

  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  
  PyObject* l   = PyList_New(0);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  PyList_Append(l, tpl);
 
  //Sortie indirection pgnids pour chaque zone
  for (size_t i=0; i < (size_t)nb_zones; ++i)
  {
    tpl = K_NUMPY::buildNumpyArray(&z_pgnids[i][0], z_pgnids[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  //Sortie indirection phnids pour chaque zone
  for (size_t i=0; i < (size_t)nb_zones; ++i)
  {
    tpl = K_NUMPY::buildNumpyArray(&z_phnids[i][0], z_phnids[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }  

  
  for (E_Int i=0; i < nb_zones; ++i)
  {
    delete crds[i];
    delete cnts[i];
  }
  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::drawOrientation(PyObject* self, PyObject* args)
{

  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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

  K_FLD::IntArray cntE;
  K_FLD::FloatArray crdE;
  for (E_Int i=0; i < ngi.PGs.size(); ++i)
  {
    
    const E_Int* nodes = ngi.PGs.get_facets_ptr(i);
    E_Int nb_nodes = ngi.PGs.stride(i);

    E_Float c[3], n[3], top[3];
    K_MESH::Polygon::centroid<3>(crd, nodes, nb_nodes, 1, c);
    K_MESH::Polygon::normal<K_FLD::FloatArray, 3>(crd, nodes, nb_nodes, 1, n);

    K_MESH::Polygon pg(nodes, nb_nodes, -1);
    double Lref = ::sqrt(pg.Lref2(crd));

    K_FUNC::sum<3>(1., c, Lref, n, top);

    E_Int id = crdE.cols();
    crdE.pushBack(c, c+3);
    crdE.pushBack(top, top+3);

    E_Int e[] = {id, id+1};
    cntE.pushBack(e, e+2);

  }
  
  PyObject* tpl = K_ARRAY::buildArray(crdE, varString, cntE, -1, "BAR", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getFaceIdsWithCentroids(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  if (!PYPARSETUPLE_(args, OO_, &arr1, &arr2)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f, cn, varString, eltType);
  if (err) return NULL;

  E_Int ni, nj, nk;
  K_FLD::FloatArray* f2(nullptr);
  K_FLD::IntArray* cn2(nullptr);
  char* varString2, *eltType2;
  K_ARRAY::getFromArray(arr2, varString2, f2, ni, nj, nk, cn2, eltType2);

  K_FLD::FloatArray& crd = *f;
  K_FLD::IntArray& cnt = *cn;

  K_FLD::FloatArray& crd2 = *f2;

  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type::eGEODIM geodim = ngon_type::get_ngon_geodim(cnt);

  //std::cout << "GEO dIM ? " << geodim << std::endl;

  if (geodim == ngon_type::eGEODIM::ERROR)
  {
    std::cout << "externalFaces : Input Error : mesh is corrupted." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::MIXED)
  {
    std::cout << "externalFaces : Input Error : mesh mixed elt types (lineic and/or surfacic and /or volumic." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::LINEIC)
  {
    std::cout << "externalFaces : Unsupported : lineic NGON are not handled." << std::endl;
    return nullptr;
  }

  // so either SURFACIC, SURFACIC_CASSIOPEE or VOLUMIC

  if (geodim == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnt1;
    ng.export_surfacic_view(cnt1);
    //std::cout << "exported" << std::endl;
    geodim = ngon_type::eGEODIM::SURFACIC;
    cnt=cnt1;
  }

  PyObject *tpl = nullptr;
  
  // SURFACIC OR VOLUMIC ?

  ngon_type ngi(cnt);
  ngon_unit& PGS = ngi.PGs;

  K_FLD::FloatArray cents;
  ngon_type::centroids(PGS, crd, cents);

  E_Int nb_pts2 = crd2.cols();
  //std::cout << "nb pts : " << nb_pts2 << std::endl;

  using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
  acrd_t acrd(cents);
  K_SEARCH::KdTree<> tree(acrd);

  std::vector<E_Int> ids;

  for (E_Int i=0; i < nb_pts2; ++i)
  {
    E_Int N = tree.getClosest(crd2.col(i));
    if (N != IDX_NONE)
      ids.push_back(N+1);
  }

  tpl = K_NUMPY::buildNumpyArray(&ids[0], ids.size(), 1, 0);

  delete f; delete cn;
  delete f2; delete cn2;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getFaceIdsCollidingVertex(PyObject* self, PyObject* args)
{
  PyObject *arr1;
  E_Float vtx[3];
  if (!PYPARSETUPLE_(args, O_ TRRR_, &arr1, &vtx[0], &vtx[1], &vtx[2])) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f, cn, varString, eltType);
  if (err) return NULL;


  K_FLD::FloatArray& crd = *f;
  K_FLD::IntArray& cnt = *cn;

  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type::eGEODIM geodim = ngon_type::get_ngon_geodim(cnt);

  //std::cout << "GEO dIM ? " << geodim << std::endl;

  if (geodim == ngon_type::eGEODIM::ERROR)
  {
    std::cout << "externalFaces : Input Error : mesh is corrupted." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::MIXED)
  {
    std::cout << "externalFaces : Input Error : mesh mixed elt types (lineic and/or surfacic and /or volumic." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::LINEIC)
  {
    std::cout << "externalFaces : Unsupported : lineic NGON are not handled." << std::endl;
    return nullptr;
  }

  // so either SURFACIC, SURFACIC_CASSIOPEE or VOLUMIC

  if (geodim == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnt1;
    ng.export_surfacic_view(cnt1);
    //std::cout << "exported" << std::endl;
    geodim = ngon_type::eGEODIM::SURFACIC;
    cnt=cnt1;
  }

  PyObject *tpl = nullptr;
  
  // SURFACIC OR VOLUMIC ?

  ngon_type ngi(cnt);

  //std::cout << "nb pts : " << nb_pts2 << std::endl;

  std::vector<E_Int> cands;
  std::set<E_Int> ids;

  double TOL = 1.e-6;

  if (geodim == ngon_type::eGEODIM::SURFACIC)
  {
    //std::cout << "mesh object..." << std::endl;
    NUGA::pg_smesh_t mesh(crd, cnt);
    //std::cout << "getting boundary..." << std::endl;
    auto loc = mesh.get_localizer();

    E_Float minB[] = { vtx[0] - TOL, vtx[1] - TOL, vtx[2] - TOL};
    E_Float maxB[] = { vtx[0] + TOL, vtx[1] + TOL, vtx[2] + TOL};
      
    K_SEARCH::BBox3D box(minB, maxB);

    cands.clear();
    loc->get_candidates(box, cands, 0);
    ids.insert(ALL(cands));
  }
  else if (geodim == ngon_type::eGEODIM::VOLUMIC)
  {
    //std::cout << "mesh object..." << std::endl;
    ngon_type ng(cnt);
    NUGA::pg_smesh_t mesh(crd, cnt);
    //std::cout << "getting boundary..." << std::endl;
    auto loc = mesh.get_localizer();

    E_Float minB[] = { vtx[0] - TOL, vtx[1] - TOL, vtx[2] - TOL};
    E_Float maxB[] = { vtx[0] + TOL, vtx[1] + TOL, vtx[2] + TOL};

    std::cout << "min : " << minB[0] << "/" << minB[1] << "/" << minB[2] << std::endl;
    std::cout << "max : " << maxB[0] << "/" << maxB[1] << "/" << maxB[2] << std::endl;
      
    K_SEARCH::BBox3D box(minB, maxB);

    cands.clear();
    loc->get_candidates(box, cands, 0);
    ids.insert(ALL(cands));
    
  }

  cands.clear();
  cands.insert(cands.end(), ALL(ids));

  tpl = K_NUMPY::buildNumpyArray(&cands[0], cands.size(), 1, 0);


  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getCells(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  E_Int is_face_id{1};

  if (!PYPARSETUPLE_(args, OO_ I_, &arr1, &arr2, &is_face_id)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f, cn, varString, eltType);
  if (err) return NULL;

  E_Int* ids{nullptr};
  E_Int size, nfld;
  if (arr2 != Py_None)
      K_NUMPY::getFromNumpyArray(arr2, ids, size, nfld);

  if (ids == nullptr) return NULL;

  K_FLD::FloatArray& crd = *f;
  K_FLD::IntArray& cnt = *cn;

  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type ng(cnt), ngo;

  std::vector<bool> keep;

  if (is_face_id)
  {
    std::set<E_Int> sids(ids, ids+size);    
    std::vector<E_Int> elts;
    ng.PHs.find_elts_with_facets(sids, elts);

    keep.resize(ng.PHs.size(), false);
    for (size_t i=0; i < elts.size(); ++i) keep[elts[i]]=true;
  }
  else // cell id
  {
    keep.resize(ng.PHs.size(), false);
    for (size_t i=0; i < (size_t)size; ++i) keep[ids[i]]=true;
  }

  std::vector<E_Int> oids;
  ngo.PGs = ng.PGs;
  for (size_t i=0; i < keep.size(); ++i)
  {
    if (keep[i] == false) continue;
    ngo.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
    oids.push_back(i+1);
  }

  ngo.PHs.updateFacets();
  std::vector<E_Int> pgnids, phnids;
  ngo.remove_unreferenced_pgs(pgnids, phnids);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  PyObject *l(PyList_New(0));
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);

  PyList_Append(l, tpl);
  Py_DECREF(tpl);


  //std::cout << "oids sz : ?" << oids.size() << std::endl;
  //std::cout << "ng PHs sz : ?" << ngo.PHs.size() << std::endl;
  
  FloatArray foids(1, oids.size());
  for (size_t i = 0; i < oids.size(); ++i) foids[i] = oids[i];

  //tpl = K_ARRAY::buildArray(foids, "oid", *cn, -1, "NGON*", true);
  tpl = K_ARRAY::buildArray(foids, "oid", cnto, -1, "NGON", false);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  delete f; delete cn;
  return l;

}

PyObject* K_INTERSECTOR::getFaces(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;
  if (!PYPARSETUPLE_(args, OO_, &arr1, &arr2)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f, cn, varString, eltType);
  if (err) return NULL;

  E_Int* pgids=NULL;
  E_Int size, nfld;
  if (arr2 != Py_None)
      K_NUMPY::getFromNumpyArray(arr2, pgids, size, nfld);

  K_FLD::FloatArray& crd = *f;
  K_FLD::IntArray& cnt = *cn;

  ngon_type ng(cnt);

  ngon_unit nguo;
  std::vector<E_Int> tmp(size);
  for (E_Int i=0; i < size; ++i)tmp[i] = pgids[i];

  std::vector< E_Int> oids;
  ng.PGs.extract(tmp, nguo, oids);

  ngon_type ngo(nguo, true);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);

  delete f; delete cn;
  return tpl;

}

//=======================  Intersector/PolyMeshTools/utils.cpp ====================
