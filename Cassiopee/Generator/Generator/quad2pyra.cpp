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

# include "generator.h"
# include <string>
# include <sstream> 

#include "Nuga/include/Polygon.h" 

//#include <iostream>
#include <memory>


E_Int check_is_of_type(const char* type, PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  E_Int ni, nj, nk;
  
  E_Int res = K_ARRAY::getFromArray(arr, varString, f1, ni, nj, nk,
                                    cn1, eltType);
     
  bool err = (res !=2);
  err |= (strcmp(eltType, type) != 0);
  if (err)
  {
    //std::cout << "input error : err => " << err << std::endl;
    //std::cout << "input error : eltType => " << eltType << std::endl;
    std::ostringstream o;
    o << "input error : invalid array, must be a unstructured " << type << " array.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
    delete f1; delete cn1;
    return 1;
  }

  // Check coordinates.
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    PyErr_SetString(PyExc_TypeError, "input error : can't find coordinates in array.");//fixme  conformUnstr
    delete f1; delete cn1;
    return 1;
  }
  
  return 0;
}

void quad_to_pyra(K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Float q)
{
  E_Int nb_q4 = cnt.cols();
  E_Int nb_pts = crd.cols();
  
  cnt.resize(5, cnt.cols());
  crd.resize(3, nb_pts + nb_q4);
  
  E_Float G[3]/*iso_bar*/, n[3]/*normal*/;
  
#pragma omp parallel for private(G, n)
  for (E_Int i=0; i< nb_q4; ++i)
  {
    // approx centroid : iso bary
    K_MESH::Polygon::iso_barycenter<K_FLD::FloatArray, 3>(crd, cnt.col(i), 4, 0, G);
    // normal
    K_MESH::Polygon::ndS<K_FLD::FloatArray, 3>(crd, cnt.col(i), 4, 0, n);
    K_FUNC::normalize<3>(n);
    
    E_Float d = ::sqrt(K_FUNC::sqrDistance(crd.col(cnt(0,i)), G, 3)); // half diagonal of the quad
    E_Float h = q * d;
    
    K_FUNC::sum<3>(1., G, h, n, crd.col(nb_pts + i));
    cnt(4,i) = nb_pts + i; 
  }
  
}

//=============================================================================
/* Creates a pyramid for each input quad */
//=============================================================================
PyObject* K_GENERATOR::quad2Pyra(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float hratio(0.5); 

  if (!PYPARSETUPLE_(args, O_ R_, &arr, &hratio)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_of_type("QUAD", arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  quad_to_pyra(crd, cnt, hratio);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnt, -1, "PYRA", false);

  delete f; delete cn;
  return tpl;
}

//=======================  Generator/quad2pyra.cpp ====================
