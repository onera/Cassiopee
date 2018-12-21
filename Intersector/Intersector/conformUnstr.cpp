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

// Making a Triangles soup conformal

# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/Boolean/TRI_Conformizer.h"
# include "Nuga/Boolean/BAR_Conformizer.h"
//#include <iostream>

using namespace std;
using namespace NUGA;

E_Int check_args(PyObject* z1arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  E_Int ni, nj, nk;
  
  E_Int res = K_ARRAY::getFromArray(z1arr, varString, f1, ni, nj, nk,
                                    cn1, eltType);
     
  bool err = (res !=2);
  err |= (strcmp(eltType, "TRI") != 0) && (strcmp(eltType, "BAR") != 0);
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "conformUnstr : invalid array, must be a unstructured TRI or BAR array.");
    delete f1; delete cn1;
    return 1;
  }

  // Check coordinates.
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    PyErr_SetString(PyExc_TypeError, "conformUnstr : can't find coordinates in array.");
    delete f1; delete cn1;
    return 1;
  }
  
  return 0;
}

//=============================================================================
/* Intersection de 2 surfaces fermees. */
//=============================================================================
PyObject* K_INTERSECTOR::conformUnstr(PyObject* self, PyObject* args)
{
  PyObject *z1arr, *z2arr;
  E_Float tolerance=0.;
  E_Int left_or_right_or_both=0, itermax=10;
  if (!PYPARSETUPLE(args, 
                    "OOdll", "OOdii", 
                    "OOfll", "OOfii", 
                    &z1arr, &z2arr, &tolerance, &left_or_right_or_both, 
                    &itermax)) { return NULL; }

  K_FLD::FloatArray* f1(0), *f2(0);
  K_FLD::IntArray* cn1(0), *cn2(0);
  char* varString, *eltType1, *eltType2;
  // Check array # 1
  E_Int err = check_args(z1arr, f1, cn1, varString, eltType1);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f1;
  K_FLD::IntArray & cnt = *cn1;
  E_Int min_z2_id = cnt.cols();
  E_Int X0 = 0;
    
  if (z2arr != Py_None)
  {
    // Check array # 2
    err = check_args(z2arr, f2, cn2, varString, eltType2);
    if (err)
      return NULL;
    
    if (strcmp(eltType1, eltType2) != 0) // do not handle mixed type yet*
    {
      PyErr_SetString(PyExc_TypeError, "conformUnstr : invalid arrays, both must be of same type, unstructured TRI or BAR arrays.");
      return NULL;
    }
    
    cn2->shift(crd.cols());
    crd.pushBack(*f2);
    cnt.pushBack(*cn2);
    
    X0 = min_z2_id;
  }
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  //~ std::cout << "colors : " << colors.size() << std::endl;
  //~ std::cout << "tolerance : " << tolerance << std::endl;
  //~ std::cout << "X0 : " << X0 << std::endl;
  
  vector<E_Int> colors;
  
  if (strcmp(eltType1, "TRI") == 0)
  {
    TRI_Conformizer<3> conformizer;
    err = conformizer.run(crd, cnt, colors, 0, tolerance, X0, itermax);
    if (err) PyErr_SetString(PyExc_TypeError, "conformUnstr : conformizer failed.");
  }
  else // BAR
  {
    BAR_Conformizer<3> conformizer;
    err = conformizer.run(crd, cnt, colors, 0, tolerance, X0, itermax);
    if (err) PyErr_SetString(PyExc_TypeError, "conformUnstr : conformizer failed.");
  }
  
  PyObject* tpl = NULL;
  if (!err)
  {
    vector<E_Int> nids;
    K_FLD::IntArray* cntOut = 0;
    E_Int rows = cnt.rows();
    
    if (left_or_right_or_both == 2) // both
    {
      cntOut = &cnt;
    }
    else if (left_or_right_or_both == 1) // right
    {
      cntOut = new K_FLD::IntArray;
      //
      for (E_Int i=0; i < cnt.cols(); ++i)
      {
        if (colors[i] >= min_z2_id)
          cntOut->pushBack(cnt.col(i), cnt.col(i)+rows);
      }
    }
    else //left
    {
      //
      cntOut = new K_FLD::IntArray;
      for (E_Int i=0; i < cnt.cols(); ++i)
      {
        if (colors[i] < min_z2_id)
          cntOut->pushBack(cnt.col(i), cnt.col(i)+rows);
      }
    }
    
    //
    K_CONNECT::MeshTool::compact_to_mesh(crd, *cntOut, nids);
    //std::cout << "nb elts : " << cntOut->cols() << " and type (nb rows) : " << cntOut->rows() << std::endl; 
    tpl = K_ARRAY::buildArray(crd, varString, *cntOut, -1, eltType1, false);
    
    if (left_or_right_or_both != 2) delete cntOut;
  }
  
  delete f1; delete f2; delete cn1; delete cn2;
  return tpl;
}



//=======================  Generator/conforrmTri.cpp ====================
