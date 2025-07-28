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
/* Calcul du rotationnel d'un champ defini par un vecteur (u,v,w) en noeuds. 
   La norme du rotationnel est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeNormCurl(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  if (!PyArg_ParseTuple(args, "OO", &array, &vars0)) return NULL;
  
  //extraction des variables constituant le vecteur dont le rot est calcule
  vector<char*> vars;
  if (PyList_Check(vars0) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeNormCurl: a list of 3 variables for curl computation must be defined.");
    return NULL; 
  }
  if (PyList_Size(vars0) != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeNormCurl: 3 variables must be defined to extract the curl.");
    return NULL;
  }
  for (int i = 0; i < PyList_Size(vars0); i++)
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
                      "computeNormCurl: varname must be a string.");
      return NULL;
    }
  }
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormCurl: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
  E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
  E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
  if (posu == -1 || posv == -1 || posw == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeNormCurl: one variable was not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posu++; posv++; posw++;

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormCurl: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld < 6)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormCurl: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  PyObject* tpl;
  char* varStringOut = new char[4];
  strcpy(varStringOut,"rot"); 
  
  // calcul du rotationnel 
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1,ni-1);
    E_Int nj1 = K_FUNC::E_max(1,nj-1);
    E_Int nk1 = K_FUNC::E_max(1,nk-1);
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray(1, varStringOut, ni1, nj1, nk1);   
    E_Float* rotnp = K_ARRAY::getFieldPtr(tpl);
    E_Float* rotyp = new E_Float[ncells]; E_Float* rotzp = new E_Float[ncells];
    computeCurlStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz),
                      f->begin(posu), f->begin(posv), f->begin(posw),
                      rotnp, rotyp, rotzp);
             
    // stockage de la norme
#pragma omp parallel for shared(rotnp, ncells) if (ncells > 50)
    for (E_Int i = 0; i < ncells; i++) rotnp[i] = sqrt(rotnp[i]*rotnp[i]+rotyp[i]*rotyp[i]+rotzp[i]*rotzp[i]);
    delete [] rotyp; delete [] rotzp;
  }
  else // non structure 
  {
    if (strcmp(eltType, "TRI")   != 0 &&
        strcmp(eltType, "QUAD")  != 0 &&
        strcmp(eltType, "TETRA") != 0 &&
        strcmp( eltType, "HEXA") != 0 &&
        strcmp(eltType, "PENTA") != 0 ) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeNormCurl: not a valid element type.");
      RELEASESHAREDU(array,f, cn); return NULL;
    }
    
    E_Int npts = f->getSize();
    tpl = K_ARRAY::buildArray(1, varStringOut, npts, cn->getSize(), -1, eltType, true, cn->getSize()*cn->getNfld());
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    E_Float* rotnp = K_ARRAY::getFieldPtr(tpl);
    E_Int nelts = cn->getSize();
    E_Float* rotyp = new E_Float[nelts]; E_Float* rotzp = new E_Float[nelts];
  
    // calcul du rotationnel aux centres des elements
    computeCurlNS(eltType, npts, *cn, 
                  f->begin(posx), f->begin(posy), f->begin(posz),
                  f->begin(posu), f->begin(posv), f->begin(posw),
                  rotnp, rotyp, rotzp);
                            
    // stockage de la norme
#pragma omp parallel for shared(rotnp, nelts) if (nelts > 50)
    for (E_Int i = 0; i < nelts; i++) rotnp[i] = sqrt(rotnp[i]*rotnp[i]+rotyp[i]*rotyp[i]+rotzp[i]*rotzp[i]);
    delete [] rotyp; delete [] rotzp;
  } 
  delete [] varStringOut;
  RELEASESHAREDB(res,array,f,cn);
  return tpl;
}
