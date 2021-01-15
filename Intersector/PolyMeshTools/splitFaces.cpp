/*    
    Copyright 2013-2021 Onera.

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


# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"
# include "Nuga/include/Splitter.h"
# include <memory>
//#include <iostream>

using namespace std;
using namespace NUGA;


// update the point list according to a split expressed by oids
PyObject* K_INTERSECTOR::updatePointLists(PyObject* self, PyObject* args)
{
  PyObject *py_oids, *py_ptLists;

  if (!PyArg_ParseTuple(args, "OO", &py_oids, &py_ptLists)) return NULL;

  E_Int nb_bcs = PyList_Size(py_ptLists);

  E_Int sz, r;
  E_Int* oids;
  E_Int res = K_NUMPY::getFromNumpyArray(py_oids, oids, sz, r, true/*shared*/);
  if (res != 1) return NULL;
  
  // WARNING : oids has no E_IDX_NONE (created entities are self referring) and is 0-based 
  E_Int nb_pgs = *std::max_element(oids, oids+sz) + 1;

  ngon_unit split_graph;
  K_CONNECT::IdTool::reverse_indirection(nb_pgs,  oids, sz, split_graph);

  PyObject *l(PyList_New(0)), *tpl;
  std::vector<E_Int> new_ptl;

  for (E_Int i=0; i < nb_bcs; ++i)
  {
    new_ptl.clear();
    PyObject* o = PyList_GetItem(py_ptLists, i);

    E_Int ptl_sz;
    FldArrayI out;
    K_ARRAY::getFromList(o, out);
    ptl_sz = out.getSize();

    for (E_Int j=0; j < ptl_sz; ++j)
    {
      E_Int oid = out(j,1)-1;
            
      E_Int nb_bits = split_graph.stride(oid);
      const E_Int* pbits = split_graph.get_facets_ptr(oid);

      if (nb_bits == 1 && pbits[0] == E_IDX_NONE)  // untouched
        new_ptl.push_back(oid+1);
      else
        for (E_Int u=0; u<nb_bits; ++u )
          new_ptl.push_back(pbits[u]+1);
    }

    tpl = K_NUMPY::buildNumpyArray(&new_ptl[0], new_ptl.size(), 1, 0);
    
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

   return l;
  
}


//=============================================================================
/* Triangulates exterior faces (any Polygon). */
//=============================================================================
PyObject* K_INTERSECTOR::triangulateExteriorFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int int_or_ext(2); //0 : internals only, 1: external only, 2: both
  E_Int improve_qual(0);

  if (!PYPARSETUPLEI(args, "Oll", "Oii", &arr, &int_or_ext, &improve_qual)) return NULL;

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
  
  ngon_type ngi(cnt), ngo;
  NUGA::transfo_t qual_param;
  qual_param.improve_qual=improve_qual ? true : false;

  Splitter::triangulate_external_pgs<DELAUNAY::Triangulator>(crd, ngi, int_or_ext, qual_param, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Triangulates BC. */
//=============================================================================
PyObject* K_INTERSECTOR::triangulateSpecifiedFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;
  E_Int /*int_or_ext(2), */improve_qual(1); //0 : internals only, 1: external only, 2: both

  if (!PYPARSETUPLEI(args, "OOl", "OOi", &arr, &py_pgs, &improve_qual)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  E_Int res=0;
  E_Int* pgsList=NULL;
  E_Int size, nfld;
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, true/*shared*/);

  if (res != 1) return NULL;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngi(cnt), ngo;

  // enable history
  ngi.PGs._ancEs.resize((E_Int)2, (E_Int)ngi.PGs.size(), (E_Int)0);
  for (E_Int i=0; i < ngi.PGs.size(); ++i) ngi.PGs._ancEs(0,i)=i;

  NUGA::transfo_t qual_param;
  qual_param.improve_qual=improve_qual ? true : false;

  Splitter::triangulate_specified_pgs<DELAUNAY::Triangulator>(crd, ngi, pgsList, size, qual_param, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  PyObject *l(PyList_New(0)), *tpl;
  
  // pushing out the mesh
  tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out history  
  E_Int sz = ngo.PGs.size();

  std::vector<E_Int> oids(sz);
  for (E_Int i = 0; i < sz; ++i) oids[i] = ngo.PGs._ancEs(0,i);

  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);
  
  
  delete f; delete cn;
  return l;
}

PyObject* K_INTERSECTOR::triangulateNFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;
  E_Int improve_qual(1), min_nvertices(4);

  if (!PYPARSETUPLEI(args, "OllO", "OiiO", &arr, &improve_qual, &min_nvertices, &py_pgs)) return NULL;

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
  
  E_Int* pgsList=NULL;
  E_Int size, nfld, res(1);
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, true/*shared*/);

  if (res != 1) return NULL;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngi(cnt), ngo;

  // enable history
  E_Int nb_pgs=(E_Int)ngi.PGs.size();
  ngi.PGs._ancEs.resize((E_Int)2, nb_pgs, (E_Int)0);
  for (E_Int i=0; i < nb_pgs; ++i) ngi.PGs._ancEs(0,i)=i;

  // discard pgs belonging to given BC
  std::vector<bool> addit(nb_pgs, true);
  if (pgsList != NULL)
    for (E_Int k = 0; k < size; ++k)
    {
     E_Int pg=pgsList[k];
      addit[pg] = false;
    }

  std::vector<E_Int> pgs;
  for (E_Int i=0; i < ngi.PGs.size(); ++i) if (ngi.PGs.stride(i) >= min_nvertices && addit[i]) pgs.push_back(i);

  NUGA::transfo_t qual_param;
  qual_param.improve_qual=improve_qual ? true : false;

  //std::cout << "nb of pgs to proceed over total : " << pgs.size() << "/" << ngi.PGs.size() << std::endl;

  Splitter::triangulate_specified_pgs<DELAUNAY::Triangulator>(crd, ngi, &pgs[0], pgs.size(), qual_param, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  PyObject *l(PyList_New(0)), *tpl;

  //std::cout << "nb pgs before/after : " << ngi.PGs.size() << "/" << ngo.PGs.size() << std::endl;
  
  // pushing out the mesh
  tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out history  
  E_Int sz = ngo.PGs.size();

  std::vector<E_Int> oids(sz);
  for (E_Int i = 0; i < sz; ++i) oids[i] = ngo.PGs._ancEs(0,i);

  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);
  
  
  delete f; delete cn;
  return l;
}

//=============================================================================
/* Convexify any concave polygon found in the input polyhedral mesh. */
//=============================================================================
PyObject* K_INTERSECTOR::convexifyFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float convexity_tol(1.e-8);

  if (!PYPARSETUPLEF(args, "Od", "Of", &arr, &convexity_tol)) return NULL;

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
  
  ngon_type ngi(cnt), ngo;
  Splitter::split_pgs(Splitter::convexify_pgs<DELAUNAY::Triangulator>, crd, ngi, convexity_tol, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  
  delete f; delete cn;
  return tpl;
}

//============================================================================================================
/* Split (convexify, starify) some targeted polygons on targeted cells 
 * (typically bad polyhedra -concaves, non-centroid-star-shaped-)
 * to prepare the split of those bad cells.*/
//============================================================================================================
PyObject* K_INTERSECTOR::prepareCellsSplit(PyObject* self, PyObject* args)
{
  PyObject *arr(0);
  E_Float PH_conc_threshold(1./3.);
  E_Float PH_cvx_threshold(0.05);
  E_Float PG_cvx_threshold(1.e-8);
  E_Int PH_set(0); // 0 for concave cells or 1 for non-centroid-star_shaped cells
  E_Int split_policy (0); // 0 : convexify concave pgs on PH set. 1 : starify concave pgs on PH set. 2 : starify any pgs at concave-chains ends.

  if (!PYPARSETUPLE(args, 
                    "Ollddd", "Oiiddd", "Ollfff", "Oiifff",
                    &arr, &PH_set, &split_policy, &PH_conc_threshold, &PH_cvx_threshold, &PG_cvx_threshold)) return NULL;

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
  
  ngon_type ngi(cnt), ngo;
  
  transfo_t::ePHset phset = transfo_t::ePHset(PH_set);
  transfo_t::eSplitPolicy policy = transfo_t::eSplitPolicy(split_policy);
  
  err = Splitter::prepareCellsSplit<DELAUNAY::Triangulator>(crd, ngi, phset, policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold, ngo);

  PyObject* tpl = NULL;
  if (err)
  {
    if (err == 1)
      PyErr_SetString(PyExc_TypeError, "prepareCellsSplit : could not orient properly.");
    else if (err == 2)
      PyErr_SetString(PyExc_TypeError, "prepareCellsSplit : could not dectect correctly pathologies.");
    else if (err == 3)
      PyErr_SetString(PyExc_TypeError, "prepareCellsSplit : could not split polygons correctly.");
  }
  else
  {
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
  
    tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  }

  delete f; delete cn;
  return tpl;
}

//=======================  Intersector/PolyMeshTools/split_faces.cpp ====================
