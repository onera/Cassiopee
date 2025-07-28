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

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Calcul de la norme du gradient d'un ensemble de champs definis en noeuds. 
   La norme du gradient est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeNormGrad(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* varname;
  if (!PyArg_ParseTuple(args, "OO", &array, &varname)) return NULL;
  
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: invalid array.");
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
                    "computeNormGrad: varname must be a string.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  } 
  E_Int posv = K_ARRAY::isNamePresent(var, varString);
  if (posv == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: variable not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posv++;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: coordinates not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld == 3)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: no field to compute.");
    RELEASESHAREDB(res,array, f, cn); return NULL;
  }
  PyObject* tpl;
  E_Int sizeVarStringOut=strlen(var)+5;//grad var
  char* varStringOut = new char[sizeVarStringOut];
  strcpy(varStringOut, "grad"); strcat(varStringOut,var);
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1,ni-1);
    E_Int nj1 = K_FUNC::E_max(1,nj-1);
    E_Int nk1 = K_FUNC::E_max(1,nk-1);    
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray(1, varStringOut, ni1, nj1, nk1);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    E_Float* fyp = new E_Float[ncells]; E_Float* fzp = new E_Float[ncells];
    computeGradStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      f->begin(posv),
                      fnp, fyp, fzp);
    // stockage de la norme
#pragma omp parallel for shared(fnp, ncells) if (ncells > 50)
    for (E_Int i = 0; i < ncells; i++) fnp[i] = sqrt(fnp[i]*fnp[i]+fyp[i]*fyp[i]+fzp[i]*fzp[i]);

    delete [] fyp; delete [] fzp;
  }
  else // non structure 
  {
    if (strcmp(eltType, "BAR")!= 0 &&
        strcmp(eltType, "TRI")!= 0 &&
        strcmp(eltType, "QUAD")!= 0 &&
        strcmp(eltType, "TETRA")!= 0 &&
        strcmp(eltType, "HEXA")!= 0 &&
        strcmp(eltType, "PENTA")!= 0 )        
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeNormGrad: not a valid element type.");
      RELEASESHAREDU(array,f, cn); return NULL;
    }
    E_Int npts = f->getSize();
    tpl = K_ARRAY::buildArray(1, varStringOut, npts, cn->getSize(), -1, eltType, true, cn->getSize()*cn->getNfld());
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    E_Int nelts = cn->getSize();
    FldArrayF fp(nelts, 1, fnp, true);
    E_Float* fyp = new E_Float[nelts]; E_Float* fzp = new E_Float[nelts];

    computeGradNS(eltType, npts, *cn, 
                  f->begin(posx), f->begin(posy), f->begin(posz), 
                  f->begin(posv), fp.begin(1), fyp, fzp);    

    // stockage de la norme
#pragma omp parallel for shared(fnp, nelts) if (nelts > 50)
    for (E_Int i = 0; i < nelts; i++) fnp[i] = sqrt(fnp[i]*fnp[i]+fyp[i]*fyp[i]+fzp[i]*fzp[i]);

    delete [] fyp; delete [] fzp;
  }
  
  RELEASESHAREDB(res,array,f,cn);
  delete [] varStringOut;
  return tpl;
}
