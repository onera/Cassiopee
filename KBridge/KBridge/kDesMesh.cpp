/*    
    Copyright 2013-2018 Onera.

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

// conversion to Cassiopee DesMesh

#include "kbridge.h"

# define E_DOUBLEREAL
//# define E_DOUBLEINT 
# define _DEF_USE_ISO_
# define _DEF_USE_USING_
# define __STD std
# define _ISO_LIST_
# define _ISO_VECTOR_
# define _E_USE_STANDARD_IOSTREAM_
# define _DEF_USE_ISO_

# include <vector>
#include <list>
#include <map>	
using namespace std;

# include "Descp/Base/DesApi.h"

//=============================================================================
// Convert an array to DesMesh
//=============================================================================
PyObject* K_KBRIDGE::array2KDesMesh(PyObject* self, PyObject* args)
{
  char* swigName;
  PyObject* array;
	
  if (!PyArg_ParseTuple(args, "Os", &array, &swigName)) return NULL;
	
  DesMesh* desMesh;
  BlkMesh* mesh;
  E_Int ret = DesApi::getInstance()->hackMesh(swigName, desMesh, mesh);
  if (ret != 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "array2KDesMesh: second arg is not a mesh.");
    return NULL;
  }
  if (mesh == NULL)
  {
    desMesh->submit();
    mesh = desMesh->getMesh();
  }
	
  // Check array
  E_Int im, jm, km;
  K_FLD::FldArrayF* f;
  K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
	
  res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, 
                              cn, eltType); 
	
  if (res != 1)
  {
    if (res == 2)  {delete f; delete cn;}
    PyErr_SetString(PyExc_TypeError,
                    "array2KDesMesh: mesh must be initialised from a structured array.");
    return NULL;
  }
	
  // Check: coordinates present ?
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
	 
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "array2KDesMesh: coordinates not found in array.");
    return NULL;
  }
  posx++; posy++; posz++;
	
  // Check: nfld = n varString ?
  E_Int nfld = f->getNfld();
  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int sizevars = vars.size(); 
  if (sizevars != nfld)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "array2KDesMesh: numpy array and varString are not coherent.");
    return NULL;
  }
	 
  // Build the node array
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Int npts =  im*jm*km;
  FldNodeF coord(npts, 3);
  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);
  for (E_Int i = 0; i < npts; i++)
  {
    coordx[i] = xt[i];
    coordy[i] = yt[i];
    coordz[i] = zt[i];
  }
	 
  // Set everything
  mesh->setNewCoord(coord, im, jm, km);
  mesh->findIfDirectOrNot();
  desMesh->setI(KEY_IM, im);
  desMesh->setI(KEY_JM, jm);
  desMesh->setI(KEY_KM, km);
	   
  // Free created data
  delete f;
  for (E_Int i = 0; i < sizevars; i++) delete [] vars[i];
	 
  Py_INCREF(Py_None);
  return Py_None;
}
