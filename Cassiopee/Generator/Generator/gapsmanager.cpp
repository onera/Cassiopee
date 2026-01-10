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
# include "Nuga/include/GapsManager.h"
using namespace std;

//=============================================================================
/* Maillage d'un trou composé de plusieurs contours à connecter. */
//=============================================================================
PyObject* K_GENERATOR::gapsmanager(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int probtype=0 /* for post, nodal mesh */;
  E_Int refine(0), coplanar(0);
  if (!PYPARSETUPLE_(args, O_ III_,
                    &arrays, &probtype, &refine, &coplanar)) return NULL;
  
  // Check every arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, "gapsmanager: arrays argument must be a list.");
    return NULL;
  }
  
  E_Int n = PyList_Size(arrays);
  PyObject* tpl;
  std::vector<K_FLD::IntArray*> components(n);
  K_FLD::FloatArray *pPos;
  E_Int res;
  char* eltType; char* varString;
  E_Int ni, nj, nk;

  // Retrieve the meshes and gather all in a unique coordinates array.
  K_FLD::FloatArray pos;
  for (E_Int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(arrays, i);
    res = K_ARRAY::getFromArray(tpl, varString, pPos, ni, nj, nk, components[i], eltType);
    
    // Coordinates
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    
    if (res != 2)
    {
      printf("Warning: gapsmanager: zone " SF_D_ " must be unstructured. Skipped.\n", i);
    }
    else 
    {
      components[i]->shift(pos.cols()); // shift indices.
      
      E_Int fields[]={posx, posy, posz};
      pos.pushBack(*pPos, fields, 3);
      
      delete pPos;
    }
  }

  // Run the gap manager.
  std::vector<K_FLD::FloatArray> posFs;
  std::vector<K_FLD::IntArray> connectFs;
  GapsManager::run(pos, components, posFs, connectFs, GapsManager::eCase(probtype), refine, GapsManager::eMode(coplanar));

  // Cleaning.
  for (E_Int i = 0; i < n; i++) delete components[i];

  // Formation des array de sortie
  PyObject* l = PyList_New(0);
  
  for (size_t i = 0; i < connectFs.size(); ++i)
  {
    tpl = K_ARRAY::buildArray(posFs[i], varString, connectFs[i], -1, eltType);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}
