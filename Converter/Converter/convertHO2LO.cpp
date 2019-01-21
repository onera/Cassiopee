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
#include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert HO mesh to LO mesh */
// ============================================================================
PyObject* K_CONVERTER::convertHO2LO(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int mode;
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &array, &mode)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray2(array, varString, 
                               f, ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2)
  {
     PyErr_SetString(PyExc_TypeError, 
                     "convertHO2LO: array is invalid.");
     return NULL;
  }
  if (res == 1)
  {   
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertHO2LO: array must be unstructured.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "convertHO2LO: array must not be NGON.");
    return NULL;
  }
  
  // Caracteristiques de l'array input
  E_Int nelts = cn->getSize();
  E_Int nfld = f->getNfld();
  E_Int nvertex = f->getSize();
  E_Int api = f->getApi();

  PyObject* o = NULL;
  if (mode == 0) // Coarse
  {
    // Caracteristique de l'array de sortie
    char outEltType[128]; E_Int d;
    K_ARRAY::eltString2TypeId(eltType, outEltType, d, d, d);

    // directement buildArray2
    o = K_ARRAY::buildArray2(nfld, varString, nvertex, nelts, -1, outEltType, false, 0, 0, 0, api);

    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);

    // Ne pas utiliser (*fo) = (*f); peut reallouer 
    fo->copy(*f, 1, nfld); // copie les champs est ok

    K_CONNECT::connectHO2LO(eltType, *cn, *co, 0);
    RELEASESHAREDU(o, fo, co);    
  }
  else if (mode == 1) // fine
  {
    if (K_STRING::cmp(eltType, 5, "BAR_3") == 0)
    {
      E_Int neltsF = nelts*2;
      o = K_ARRAY::buildArray2(nfld, varString, nvertex, neltsF, -1, "BAR", false, 0, 0, 0, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray2(o, fo, co);
      fo->copy(*f, 1, nfld);       
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(2*i,1) = (*cn)(i,1);
        (*co)(2*i,2) = (*cn)(i,3);
        (*co)(2*i+1,1) = (*cn)(i,3);
        (*co)(2*i+1,2) = (*cn)(i,2);
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 5, "TRI_6") == 0)
    {
      E_Int neltsF = nelts*4;
      o = K_ARRAY::buildArray2(nfld, varString, nvertex, neltsF, -1, "TRI", false, 0, 0, 0, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray2(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(4*i,1) = (*cn)(i,1);
        (*co)(4*i,2) = (*cn)(i,4);
        (*co)(4*i,3) = (*cn)(i,6);

        (*co)(4*i+1,1) = (*cn)(i,4);
        (*co)(4*i+1,2) = (*cn)(i,2);
        (*co)(4*i+1,3) = (*cn)(i,5);
        
        (*co)(4*i+2,1) = (*cn)(i,6);
        (*co)(4*i+2,2) = (*cn)(i,4);
        (*co)(4*i+2,3) = (*cn)(i,5);
        
        (*co)(4*i+3,1) = (*cn)(i,3);
        (*co)(4*i+3,2) = (*cn)(i,6);
        (*co)(4*i+3,3) = (*cn)(i,5);
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 6, "QUAD_8") == 0)
    {
      E_Int neltsF = nelts*6;
      o = K_ARRAY::buildArray2(nfld, varString, nvertex, neltsF, -1, "TRI", false, 0, 0, 0, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray2(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(6*i,1) = (*cn)(i,1);
        (*co)(6*i,2) = (*cn)(i,5);
        (*co)(6*i,3) = (*cn)(i,8);

        (*co)(6*i+1,1) = (*cn)(i,2);
        (*co)(6*i+1,2) = (*cn)(i,6);
        (*co)(6*i+1,3) = (*cn)(i,5);
        
        (*co)(6*i+2,1) = (*cn)(i,3);
        (*co)(6*i+2,2) = (*cn)(i,7);
        (*co)(6*i+2,3) = (*cn)(i,6);
        
        (*co)(6*i+3,1) = (*cn)(i,4);
        (*co)(6*i+3,2) = (*cn)(i,8);
        (*co)(6*i+3,3) = (*cn)(i,7);

        (*co)(6*i+4,1) = (*cn)(i,7);
        (*co)(6*i+4,2) = (*cn)(i,5);
        (*co)(6*i+4,3) = (*cn)(i,6);

        (*co)(6*i+5,1) = (*cn)(i,7);
        (*co)(6*i+5,2) = (*cn)(i,8);
        (*co)(6*i+5,3) = (*cn)(i,5);
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 6, "QUAD_9") == 0)
    {
      E_Int neltsF = nelts*8;
      o = K_ARRAY::buildArray2(nfld, varString, nvertex, neltsF, -1, "TRI", false, 0, 0, 0, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray2(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(8*i,1) = (*cn)(i,1);
        (*co)(8*i,2) = (*cn)(i,5);
        (*co)(8*i,3) = (*cn)(i,8);

        (*co)(8*i+1,1) = (*cn)(i,8);
        (*co)(8*i+1,2) = (*cn)(i,5);
        (*co)(8*i+1,3) = (*cn)(i,9);
        
        (*co)(8*i+2,1) = (*cn)(i,9);
        (*co)(8*i+2,2) = (*cn)(i,5);
        (*co)(8*i+2,3) = (*cn)(i,2);
        
        (*co)(8*i+3,1) = (*cn)(i,9);
        (*co)(8*i+3,2) = (*cn)(i,2);
        (*co)(8*i+3,3) = (*cn)(i,6);

        (*co)(8*i+4,1) = (*cn)(i,7);
        (*co)(8*i+4,2) = (*cn)(i,9);
        (*co)(8*i+4,3) = (*cn)(i,6);

        (*co)(8*i+5,1) = (*cn)(i,7);
        (*co)(8*i+5,2) = (*cn)(i,6);
        (*co)(8*i+5,3) = (*cn)(i,3);

        (*co)(8*i+6,1) = (*cn)(i,4);
        (*co)(8*i+6,2) = (*cn)(i,8);
        (*co)(8*i+6,3) = (*cn)(i,9);

        (*co)(8*i+7,1) = (*cn)(i,4);
        (*co)(8*i+7,2) = (*cn)(i,9);
        (*co)(8*i+7,3) = (*cn)(i,7);
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 8, "TETRA_10") == 0)
    {
      E_Int N = 8;
      E_Int neltsF = nelts*N;
      o = K_ARRAY::buildArray2(nfld, varString, nvertex, neltsF, -1, "TETRA", false, 0, 0, 0, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray2(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(N*i,1) = (*cn)(i,1);
        (*co)(N*i,2) = (*cn)(i,5);
        (*co)(N*i,3) = (*cn)(i,7);
        (*co)(N*i,4) = (*cn)(i,8);

        (*co)(N*i+1,1) = (*cn)(i,2);
        (*co)(N*i+1,2) = (*cn)(i,6);
        (*co)(N*i+1,3) = (*cn)(i,9);
        (*co)(N*i+1,4) = (*cn)(i,5);
        
        (*co)(N*i+2,1) = (*cn)(i,5);
        (*co)(N*i+2,2) = (*cn)(i,6);
        (*co)(N*i+2,3) = (*cn)(i,7);
        (*co)(N*i+2,4) = (*cn)(i,9);
        
        (*co)(N*i+3,1) = (*cn)(i,6);
        (*co)(N*i+3,2) = (*cn)(i,3);
        (*co)(N*i+3,3) = (*cn)(i,10);
        (*co)(N*i+3,4) = (*cn)(i,7);
        
        (*co)(N*i+4,1) = (*cn)(i,10);
        (*co)(N*i+4,2) = (*cn)(i,9);
        (*co)(N*i+4,3) = (*cn)(i,8);
        (*co)(N*i+4,4) = (*cn)(i,7);
        
        (*co)(N*i+5,1) = (*cn)(i,5);
        (*co)(N*i+5,2) = (*cn)(i,8);
        (*co)(N*i+5,3) = (*cn)(i,7);
        (*co)(N*i+5,4) = (*cn)(i,9);

        (*co)(N*i+6,1) = (*cn)(i,6);
        (*co)(N*i+6,2) = (*cn)(i,10);
        (*co)(N*i+6,3) = (*cn)(i,9);
        (*co)(N*i+6,4) = (*cn)(i,7);
      
        (*co)(N*i+7,1) = (*cn)(i,4);
        (*co)(N*i+7,2) = (*cn)(i,8);
        (*co)(N*i+7,3) = (*cn)(i,9);
        (*co)(N*i+7,4) = (*cn)(i,10);
      }
      RELEASESHAREDU(o, fo, co);   
    }
  }
  RELEASESHAREDB(res, array, f, cn);
  return o;
}
