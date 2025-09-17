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
# include <string.h>
# include "post.h"

//=============================================================================
/* Compute the divergence of a set of vector fields given at mesh nodes.
   The divergence is given on cell centers. */
//=============================================================================
PyObject* K_POST::computeDiv(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  if (!PYPARSETUPLE_(args, OO_, &array, &vars0)) return NULL;

  // extract variables constituting components of the vector whose div is calculated
  std::vector<char*> vars;
  if (PyList_Check(vars0) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: a list of 3 variables for div computation must be defined.");
    return NULL;
  }
  Py_ssize_t nvars = PyList_Size(vars0);
  if (nvars != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: 3 variables must be defined to extract the div.");
    return NULL;
  }
  for (Py_ssize_t i = 0; i < nvars; i++)
  {
    PyObject* tpl0 = PyList_GetItem(vars0, i);
    if (PyString_Check(tpl0))
    {
      char* str = PyString_AsString(tpl0);
      vars.push_back(str);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0)) 
    {
      char* str = (char*)PyUnicode_AsUTF8(tpl0);
      vars.push_back(str);
    }
#endif
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeDiv: varname must be a string.");
      return NULL;
    }
  }

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: invalid array.");
    return NULL;
  }
  E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
  E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
  E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
  if (posu == -1 || posv == -1 || posw == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: at least one variable was not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posu++; posv++; posw++;

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld < 6)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  if (strncmp(vars[0], vars[1], strlen(vars[0])-1) != 0 ||
      strncmp(vars[1], vars[2], strlen(vars[1])-1) != 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: invalid names for vector component fields.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  char* sv0 = vars[0]; char* sv1 = vars[1]; char* sv2 = vars[2];
  char s0 = sv0[strlen(sv0)-1];
  char s1 = sv1[strlen(sv1)-1];
  char s2 = sv2[strlen(sv2)-1];
  if (s0 != 'X' || s1 != 'Y' || s2 != 'Z')
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: error with the order of given scalar fields.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  char* varStringOut = new char [strlen(vars[0])+3]; // +3 for "div", -1 for trailing 'X'
  strcpy(varStringOut, "div");
  char* pt = varStringOut;
  char* v = vars[0];
  pt += 3;
  for (size_t j = 0; j < strlen(v)-1; j++)
  {
    *pt = v[j]; pt++;
  }
  *pt = '\0';

  PyObject* tpl = NULL;
  FldArrayF* f2;
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1, ni-1);
    E_Int nj1 = K_FUNC::E_max(1, nj-1);
    E_Int nk1 = K_FUNC::E_max(1, nk-1);
    tpl = K_ARRAY::buildArray3(1, varStringOut, ni1, nj1, nk1);
    K_ARRAY::getFromArray3(tpl, f2);
    E_Int ierr = computeDivStruct(
      ni, nj, nk,
      f->begin(posx), f->begin(posy), f->begin(posz),
      f->begin(posu), f->begin(posv), f->begin(posw),
      f2->begin()
    );
    if (ierr == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeDiv: Not valid for 1D.");
      delete [] varStringOut;
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
  }
  else if (res == 2)
  {
    E_Int npts = f->getSize();
    E_Int api = f->getApi();
    E_Bool center = true;
    E_Bool copyConnect = true;

    if (strcmp(eltType, "NGON") == 0)
    {
      tpl = K_ARRAY::buildArray3(
        1, varStringOut, npts, *cn, "NGON",
        center, api, copyConnect
      );
      K_ARRAY::getFromArray3(tpl, f2);
      E_Int ierr = computeDivNGon(
        f->begin(posx), f->begin(posy), f->begin(posz),
        f->begin(posu), f->begin(posv), f->begin(posw), *cn,
        f2->begin()
      );
      if (ierr == 1)
      {
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv: divergence can only be computed for 3D NGONs.");
        delete [] varStringOut;
        RELEASESHAREDS(tpl, f2);
        RELEASESHAREDB(res, array, f, cn); return NULL;
      }
    }
    else  // ME
    {
      tpl = K_ARRAY::buildArray3(
        1, varStringOut, npts, *cn, eltType,
        center, api, copyConnect
      );
      K_ARRAY::getFromArray3(tpl, f2);
      computeDivUnstruct(
        *cn, eltType,
        f->begin(posx), f->begin(posy), f->begin(posz),
        f->begin(posu), f->begin(posv), f->begin(posw),
        f2->begin()
      );
    }
  }
  delete [] varStringOut;
  RELEASESHAREDS(tpl, f2);
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
