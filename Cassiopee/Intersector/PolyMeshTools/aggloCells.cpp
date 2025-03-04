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

# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"
# include "Nuga/include/Agglomerator.h"


using namespace std;
using namespace NUGA;


//=============================================================================
/* Agglomerate cells with a too high aspect ratio */
//=============================================================================
PyObject* K_INTERSECTOR::agglomerateSmallCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float vmin(0.), vratio(0.01), angle_threshold{1.e-12};
  E_Int method(0);

  if (!PYPARSETUPLE_(args, O_ RRR_ I_,
                     &arr, &vmin, &vratio, &angle_threshold, &method)) return NULL;

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

  E_Int nb_aggs(0);
  NUGA::Agglomerator::agglomerate_small_phs<DELAUNAY::Triangulator>(crd, ngi, vmin, vratio, ngo, nb_aggs, angle_threshold, method);

  PyObject *l(PyList_New(0)), *tpl;

  {
    // zone 1 : mesh
    {
      K_FLD::IntArray cnto;
      ngo.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, "NGON", false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }

  delete f; delete cn;
  return l;
}

//=============================================================================
/* Eradicate by shell-agglomeration cells with a too high aspect ratio */
//=============================================================================
PyObject* K_INTERSECTOR::shellAgglomerateSmallCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float vmin(0.), vratio(1000.);

  if (!PYPARSETUPLE_(args, O_ RR_, &arr, &vmin, &vratio)) return NULL;

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

  E_Int nb_aggs(0);
  NUGA::Agglomerator::shell_agglomerate_small_phs<DELAUNAY::Triangulator>(crd, ngi, vmin, vratio, ngo, nb_aggs);

  PyObject *l(PyList_New(0)), *tpl;

  {
    // zone 1 : mesh
    {
      K_FLD::IntArray cnto;
      ngo.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, "NGON", false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

  }

  delete f; delete cn;
  return l;
}

//=============================================================================
/* Agglomerate cells with a too high aspect ratio */
//=============================================================================
PyObject* K_INTERSECTOR::agglomerateNonStarCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int debug=0;
  double angle_threshold{1.e-12};

  if (!PYPARSETUPLE_(args, O_ R_, &arr, &angle_threshold)) return NULL;

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

  E_Int nb_aggs(0);
  NUGA::Agglomerator::agglomerate_non_star_phs<DELAUNAY::Triangulator>(crd, ngi, ngo, nb_aggs, angle_threshold);

  PyObject *l(PyList_New(0)), *tpl;

  {
    // zone 1 : mesh
    {
      K_FLD::IntArray cnto;
      ngo.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, "NGON", false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

    if (debug)
    {
      for (E_Int i=0; i < nb_aggs; ++i)
      {
        ngon_unit ph;
        ph.add(ngo.PHs.stride(i), ngo.PHs.get_facets_ptr(i));

        ngon_type one_ph(ngo.PGs, ph);
        std::vector<E_Int> pgnids, phnids;
        one_ph.remove_unreferenced_pgs(pgnids, phnids);
  
        K_FLD::FloatArray cr(crd);
        ngon_type::compact_to_used_nodes(one_ph.PGs, cr);

        K_FLD::IntArray cnto;
        one_ph.export_to_array(cnto);

        tpl = K_ARRAY::buildArray(cr, varString, cnto, -1, "NGON", false);
        PyList_Append(l, tpl);
        Py_DECREF(tpl);
      }
    }
  }

  delete f; delete cn;
  return l;
}

//=============================================================================
/* Agglomerate cells where polygons are specified */
//=============================================================================
PyObject* K_INTERSECTOR::agglomerateCellsWithSpecifiedFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;

  if (!PyArg_ParseTuple(args, "OO", &arr, &py_pgs)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  E_Int* pgsList=NULL;
  E_Int size=0, nfld=1;
  if (py_pgs != Py_None)
    K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, 1);

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngio(cnt);

  std::vector<E_Int> pgnids;
  NUGA::Agglomerator::agglomerate_phs_having_pgs(crd, ngio, pgsList, size, pgnids);

  K_FLD::IntArray cnto;
  ngio.export_to_array(cnto);

  PyObject *l(PyList_New(0)), *tpl;

  tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, "NGON", false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  //convert IDX_NONE value into -1
  for (size_t k=0; k < pgnids.size(); ++k)
    if (pgnids[k] == IDX_NONE)
      pgnids[k]=-1;

  tpl = K_NUMPY::buildNumpyArray(&pgnids[0], pgnids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);
  
  delete f; delete cn;
  return l;

}
//=======================  Intersector/PolyMeshTools/aggloFaces.cpp ====================
