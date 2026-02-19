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
# include <string.h>
# include "post.h"

//=============================================================================
/* Calcul du rotationnel d'un champ defini par un vecteur (u,v,w) en noeuds. 
   Le rotationnel est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeCurl(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  if (!PYPARSETUPLE_(args, OO_, &array, &vars0)) return NULL;
  
  //extraction des variables constituant le vecteur dont le rot est calcule
  std::vector<char*> vars;
  if (PyList_Check(vars0) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: a list of 3 variables for curl computation must be defined.");
    return NULL; 
  }
  if (PyList_Size(vars0) != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: 3 variables must be defined to extract the curl.");
    return NULL;
  }
  for (Py_ssize_t i = 0; i < PyList_Size(vars0); i++)
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
                      "computeCurl: varname must be a string.");
      return NULL;
    }
  }
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
  E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
  E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
  if (posu == -1 || posv == -1 || posw == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: one variable was not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posu++; posv++; posw++;
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  posx++; posy++; posz++;
  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();
  if (nfld < 6)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  // calcul du rotationnel 
  PyObject* tpl = NULL;
  FldArrayF* f2;
  char* varStringOut = new char[15];
  strcpy(varStringOut, "rotx,roty,rotz"); 
  
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1, ni-1);
    E_Int nj1 = K_FUNC::E_max(1, nj-1);
    E_Int nk1 = K_FUNC::E_max(1, nk-1);
    tpl = K_ARRAY::buildArray3(3, varStringOut, ni1, nj1, nk1, api);
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* rotx = f2->begin(1);
    E_Float* roty = f2->begin(2);
    E_Float* rotz = f2->begin(3);

    computeCurlStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz),
                      f->begin(posu), f->begin(posv), f->begin(posw),
                      rotx, roty, rotz);
  }
  else if (res == 2)  // non structure
  {
    E_Int npts = f->getSize();
    E_Int api = f->getApi();
    E_Bool center = true;
    E_Bool copyConnect = true;

    if (strcmp(eltType, "NGON") == 0)
    {
      tpl = K_ARRAY::buildArray3(
        3, varStringOut, npts, *cn, "NGON",
        center, api, copyConnect);
      K_ARRAY::getFromArray3(tpl, f2);
      E_Float* rotx = f2->begin(1);
      E_Float* roty = f2->begin(2);
      E_Float* rotz = f2->begin(3);

      E_Int ierr = computeCurlNGon(
        f->begin(posx), f->begin(posy), f->begin(posz),
        f->begin(posu), f->begin(posv), f->begin(posw), *cn, 
        rotx, roty, rotz);    

      if (ierr == 1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "computeCurl: curl can only be computed for 3D NGONs.");
        delete [] varStringOut;
        RELEASESHAREDS(tpl, f2);
        RELEASESHAREDB(res, array, f, cn); return NULL;         
      }      
    }
    else  // ME
    {
      tpl = K_ARRAY::buildArray3(
        3, varStringOut, npts, *cn, eltType,
        center, api, copyConnect);
      K_ARRAY::getFromArray3(tpl, f2);
      E_Float* rotx = f2->begin(1);
      E_Float* roty = f2->begin(2);
      E_Float* rotz = f2->begin(3);
      
      computeCurlUnstruct(
        *cn, eltType,
        f->begin(posx), f->begin(posy), f->begin(posz),
        f->begin(posu), f->begin(posv), f->begin(posw),
        rotx, roty, rotz
      );              
    }   
  }
  delete [] varStringOut;
  RELEASESHAREDS(tpl, f2);
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
