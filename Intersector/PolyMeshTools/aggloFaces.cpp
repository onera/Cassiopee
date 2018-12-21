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


# include <string>
# include <sstream> 
# include "intersector.h"
# include "Fld/ngon_t.hxx"
# include "Nuga/Delaunay/Triangulator.h"
# include "Nuga/Boolean/Agglomerator.h"

//#include <iostream>

using namespace std;
using namespace NUGA;


//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::simplifyCells(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float angular_threshold(0.);
  E_Int treat_externals(1);

  if (!PYPARSETUPLE(args, "Old", "Oid", "Olf", "Oif", 
                    &arr, &treat_externals, &angular_threshold)) return NULL;

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
  
  ngon_unit orienti, oriento, neighborsi, neighborso;
  ngi.build_ph_neighborhood(neighborsi);
  ngi.build_orientation_ngu<DELAUNAY::Triangulator>(crd, ngi, orienti);
  
  NUGA::Agglomerator::simplify_phs (crd, ngi, orienti, neighborsi, angular_threshold, bool(treat_externals), ngo, oriento, neighborso);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  
  delete f; delete cn;
  return tpl;
}

PyObject* K_INTERSECTOR::collapseUncomputableFaces(PyObject* self, PyObject* args)
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


//=======================  Intersector/PolyMeshTools/aggloFaces.cpp ====================
