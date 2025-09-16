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
/* Calcul du gradient d un ensemble de champs definis en noeuds. 
   Le gradient est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeGrad(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* varname;
  if (!PYPARSETUPLE_(args, OO_, &array, &varname)) return NULL;
  
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, 
                                     eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: invalid array.");
    return NULL;
  }
  // check varname
  char* var = NULL;
  if (PyString_Check(varname)) var = PyString_AsString(varname);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(varname)) var = (char*)PyUnicode_AsUTF8(varname);
#endif
  else
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: varname must be a string.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  } 
  E_Int posv = K_ARRAY::isNamePresent(var, varString);
  if (posv == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: variable not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posv++;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld == 3)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  PyObject* tpl = NULL;
  FldArrayF* f2;
  char* varStringOut;
  computeGradVarsString(var, varStringOut);

  if (res == 1) 
  {
    E_Int ni1 = K_FUNC::E_max(1, ni-1);
    E_Int nj1 = K_FUNC::E_max(1, nj-1);
    E_Int nk1 = K_FUNC::E_max(1, nk-1);
    tpl = K_ARRAY::buildArray3(3, varStringOut, ni1, nj1, nk1);
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* gradx = f2->begin(1);
    E_Float* grady = f2->begin(2);
    E_Float* gradz = f2->begin(3);

    computeGradStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      f->begin(posv),
                      gradx, grady, gradz);
  }
  else if (res == 2) 
  {
    E_Int npts = f->getSize();
    E_Int api = f->getApi();
    E_Bool center = true;
    E_Bool copyConnect = true;

    if (strcmp(eltType, "NGON") == 0)
    {
      // Build unstructured NGON array from existing connectivity & empty fields
      tpl = K_ARRAY::buildArray3(
        3, varStringOut, npts, *cn, "NGON",
        center, api, copyConnect
      );
      K_ARRAY::getFromArray3(tpl, f2);
      E_Float* gradx = f2->begin(1);
      E_Float* grady = f2->begin(2);
      E_Float* gradz = f2->begin(3);
      
      E_Int ierr = computeGradNGon(
        f->begin(posx), f->begin(posy), f->begin(posz), 
        f->begin(posv), *cn, 
        gradx, grady, gradz);

      if (ierr == 1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "computeGrad: gradient can only be computed for 3D NGONs.");
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
      E_Float* gradx = f2->begin(1);
      E_Float* grady = f2->begin(2);
      E_Float* gradz = f2->begin(3);
      computeGradUnstruct(
        f->begin(posx), f->begin(posy), f->begin(posz), *cn, eltType,
        f->begin(posv), gradx, grady, gradz
      );
    }
  }
  delete [] varStringOut;
  RELEASESHAREDS(tpl, f2);
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}

//=============================================================================
/* A partir de la chaine de variables initiale: (x,y,z,var1,var2,...)
   Cree la chaine (gradxvar1,gradyvar1,gradzvar1, gradxvar2, ....) 
   Cette routine alloue varStringOut */
//=============================================================================
void K_POST::computeGradVarsString(char* varString, char*& varStringOut)
{
  std::vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int c = -1;
  E_Int varsSize = vars.size();
  E_Int sizeVarStringOut = 0;
  for (E_Int v = 0; v < varsSize; v++)
  {
    E_Int vsize = strlen(vars[v]);
    sizeVarStringOut += vsize + 6;
  }
  varStringOut = new char [3*sizeVarStringOut];

  for (E_Int v = 0; v < varsSize; v++)
  {
    char*& var0 = vars[v];
    if (strcmp(var0, "x") != 0 && 
        strcmp(var0, "y") != 0 && 
        strcmp(var0, "z") != 0)
    {
      if (c == -1)
      {
        strcpy(varStringOut, "gradx");
        c = 1;
      }
      else strcat(varStringOut, ",gradx");
      strcat(varStringOut, var0);
      strcat(varStringOut, ",grady");
      strcat(varStringOut, var0);
      strcat(varStringOut, ",gradz");
      strcat(varStringOut, var0);
    }
  } 
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];
}
