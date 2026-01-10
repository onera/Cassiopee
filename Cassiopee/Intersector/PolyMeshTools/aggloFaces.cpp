/*    
    Copyright 2013-2026 ONERA.

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

//#include <iostream>

using namespace std;
using namespace NUGA;

//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::simplifyCells(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_skipids;
  E_Float angular_threshold(0.);
  E_Int treat_externals(1);

  if (!PYPARSETUPLE_(args, O_ I_ R_ O_, 
                    &arr, &treat_externals, &angular_threshold, &py_skipids)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // get PG ids to skip
  std::vector<E_Int>* skipPGids = nullptr;
  if (py_skipids != Py_None)
  {
    skipPGids = new std::vector<E_Int>;
    E_Int *ptL, size, nfld;
    K_NUMPY::getFromNumpyArray(py_skipids, ptL, size, nfld);
    skipPGids->insert(skipPGids->end(), ptL, ptL + size);
  }
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), ngo;
  
  ngon_unit orienti, oriento, neighborsi, neighborso;
  ngi.build_ph_neighborhood(neighborsi);
  ngi.build_orientation_ngu<DELAUNAY::Triangulator>(crd, ngi, orienti);
  
  NUGA::Agglomerator::simplify_phs (crd, ngi, orienti, neighborsi, angular_threshold, bool(treat_externals), ngo, oriento, neighborso, nullptr/*PHlist*/, skipPGids);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  if (skipPGids != nullptr) delete skipPGids;
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Agglomerate superfuous nodes (overdefined polygons) */
//=============================================================================
PyObject* K_INTERSECTOR::simplifyFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) return nullptr;

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
  
  ngon_type::simplify_pgs(ngi, crd);

  K_FLD::IntArray cnto;
  ngi.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::simplifySurf(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float angular_threshold(0.);

  if (!PYPARSETUPLE_(args, O_ R_, &arr, &angular_threshold)) return NULL;

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
  std::vector<E_Int> nids, orient(ngi.PGs.size(), 1);
  ngon_unit agglo_pgs;

  K_MESH::Polygon::full_agglomerate(crd, ngi.PGs, ngi.PGs.get_facets_ptr(0), ngi.PGs.size(), angular_threshold, &orient[0], agglo_pgs, nids);
  
  ngon_type ngo(agglo_pgs, true);
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;

  delete f; delete cn;
  return tpl;
}

PyObject* K_INTERSECTOR::collapseUncomputableFaces(PyObject* self, PyObject* args)
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
  
  std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngio(cnt);

  NUGA::Agglomerator::collapse_uncomputable_pgs<DELAUNAY::Triangulator>(crd, ngio);

  K_FLD::IntArray cnto;
  ngio.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);

  delete f; delete cn;
  return tpl;
}

PyObject* K_INTERSECTOR::collapseSmallCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float vmin{0.}, grmin{0.};

  if (!PYPARSETUPLE_(args, O_ RR_, &arr, &vmin, &grmin)) return NULL;

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

  //double vmin = 3.e-15;//itermax 3
  //double vmin  = 1.e-16; //itermax 4
  //double vratio = -1;
  NUGA::Agglomerator::collapse_small_tetras2<DELAUNAY::Triangulator>(crd, ngio, vmin, grmin);

  K_FLD::IntArray cnto;
  ngio.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);

  delete f; delete cn;
  return tpl;
}


//=======================  Intersector/PolyMeshTools/aggloFaces.cpp ====================
