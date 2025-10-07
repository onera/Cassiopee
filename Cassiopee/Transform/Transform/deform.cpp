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

# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Deform mesh by moving surface of a given vector. */
// ============================================================================
PyObject* K_TRANSFORM::deform(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vectorNames;
  if (!PYPARSETUPLE_(args, OO_, &array, &vectorNames))
      return NULL;

  if (PyList_Check(vectorNames) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "deform: vector names must be a list.");
    return NULL;
  }
  E_Int nvectors = PyList_Size(vectorNames);
  if (nvectors != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "deform: vector names must be a list of 3 strings.");
    return NULL;
  }

  // Check array
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(array, varString1, f1,
                                      im1, jm1, km1, cn1, eltType1);

  // Vecteur et array valides (structure ou non structure) ?
  if (res1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "deform: 1st argument is invalid.");
    return NULL;
  }
  vector<char*> vars;
  for (E_Int nov = 0; nov < 3; nov++)
  {
    char* vects;
    PyObject* tpl0 = PyList_GetItem(vectorNames, nov);
    if (PyString_Check(tpl0))
    {
      vects = PyString_AsString(tpl0);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0))
    {
      vects = (char*)PyUnicode_AsUTF8(tpl0);
    }
#endif
    else
    {
      RELEASESHAREDB(res1, array,f1,cn1);
      PyErr_SetString(PyExc_TypeError,
                      "deform: vector fields must be strings.");
      return NULL;
    }
    vars.push_back(vects);
  }
  // Presence des coordonnees ?
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    PyErr_SetString(PyExc_ValueError,
                    "deform: coordinates not found in 1st argument.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int posdx = K_ARRAY::isNamePresent(vars[0], varString1);
  E_Int posdy = K_ARRAY::isNamePresent(vars[1], varString1);
  E_Int posdz = K_ARRAY::isNamePresent(vars[2], varString1);
  if ( posdx == -1 || posdy == -1 || posdz == -1)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    PyErr_SetString(PyExc_ValueError,
                    "deform: vector not found in zone.");
    return NULL;
  }
  posdx++; posdy++; posdz++;

  // Construit l'array resultat et l'initialise par copie
  E_Int api = f1->getApi();
  E_Int npts = f1->getSize();
  PyObject* tpl;
  FldArrayF* newF;

  if (res1 == 1) // structured
  {
    tpl = K_ARRAY::buildArray3(*f1, varString1, im1, jm1, km1, api);
  }
  else // unstructured
  {
    tpl = K_ARRAY::buildArray3(*f1, varString1, *cn1, eltType1, api);
  }
  K_ARRAY::getFromArray3(tpl, newF);

  // Deformation
  E_Float* newFx = newF->begin(posx);
  E_Float* newFy = newF->begin(posy);
  E_Float* newFz = newF->begin(posz);
  E_Float* f1x = f1->begin(posx); E_Float* f2x = f1->begin(posdx);
  E_Float* f1y = f1->begin(posy); E_Float* f2y = f1->begin(posdy);
  E_Float* f1z = f1->begin(posz); E_Float* f2z = f1->begin(posdz);

  #pragma omp parallel for
  for (E_Int i = 0; i < npts; i++)
  {
    newFx[i] = f1x[i] + f2x[i];
    newFy[i] = f1y[i] + f2y[i];
    newFz[i] = f1z[i] + f2z[i];
  }

  RELEASESHAREDS(tpl, newF);
  RELEASESHAREDB(res1, array, f1, cn1);
  return tpl;
}
