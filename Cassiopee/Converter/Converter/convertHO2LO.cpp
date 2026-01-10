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
  if (!PYPARSETUPLE_(args, O_ I_, &array, &mode)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, 
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

    // directement buildArray3
    o = K_ARRAY::buildArray3(nfld, varString, nvertex, nelts, outEltType, false, api);

    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray3(o, fo, co);

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
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "BAR", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
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
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "TRI", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
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
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "TRI", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
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
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "TRI", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
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
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "TETRA", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
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
    else if (K_STRING::cmp(eltType, 7, "HEXA_20") == 0)
    {
      E_Int N = 22;
      E_Int neltsF = nelts*N;
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "TETRA", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(N*i,1) = (*cn)(i,1);
        (*co)(N*i,2) = (*cn)(i,9);
        (*co)(N*i,3) = (*cn)(i,12);
        (*co)(N*i,4) = (*cn)(i,13);

        (*co)(N*i+1,1) = (*cn)(i,2);
        (*co)(N*i+1,2) = (*cn)(i,10);
        (*co)(N*i+1,3) = (*cn)(i,9);
        (*co)(N*i+1,4) = (*cn)(i,14);
 
        (*co)(N*i+2,1) = (*cn)(i,3);
        (*co)(N*i+2,2) = (*cn)(i,11);
        (*co)(N*i+2,3) = (*cn)(i,10);
        (*co)(N*i+2,4) = (*cn)(i,15);

        (*co)(N*i+3,1) = (*cn)(i,11);
        (*co)(N*i+3,2) = (*cn)(i,4);
        (*co)(N*i+3,3) = (*cn)(i,12);
        (*co)(N*i+3,4) = (*cn)(i,16);

        (*co)(N*i+4,1) = (*cn)(i,15);
        (*co)(N*i+4,2) = (*cn)(i,16);
        (*co)(N*i+4,3) = (*cn)(i,11);
        (*co)(N*i+4,4) = (*cn)(i,12);

        (*co)(N*i+5,1) = (*cn)(i,10);
        (*co)(N*i+5,2) = (*cn)(i,11);
        (*co)(N*i+5,3) = (*cn)(i,12);
        (*co)(N*i+5,4) = (*cn)(i,15);
        
        (*co)(N*i+6,1) = (*cn)(i,9);
        (*co)(N*i+6,2) = (*cn)(i,14);
        (*co)(N*i+6,3) = (*cn)(i,13);
        (*co)(N*i+6,4) = (*cn)(i,12);
      
        (*co)(N*i+7,1) = (*cn)(i,9);
        (*co)(N*i+7,2) = (*cn)(i,10);
        (*co)(N*i+7,3) = (*cn)(i,14);
        (*co)(N*i+7,4) = (*cn)(i,12);
       
        (*co)(N*i+8,1) = (*cn)(i,15);
        (*co)(N*i+8,2) = (*cn)(i,12);
        (*co)(N*i+8,3) = (*cn)(i,16);
        (*co)(N*i+8,4) = (*cn)(i,20);
        
        (*co)(N*i+9,1) = (*cn)(i,10);
        (*co)(N*i+9,2) = (*cn)(i,12);
        (*co)(N*i+9,3) = (*cn)(i,15);
        (*co)(N*i+9,4) = (*cn)(i,20);
        
        
        (*co)(N*i+10,1) = (*cn)(i,10);
        (*co)(N*i+10,2) = (*cn)(i,18);
        (*co)(N*i+10,3) = (*cn)(i,20);
        (*co)(N*i+10,4) = (*cn)(i,15);
        
        
        (*co)(N*i+11,1) = (*cn)(i,16);
        (*co)(N*i+11,2) = (*cn)(i,8);
        (*co)(N*i+11,3) = (*cn)(i,19);
        (*co)(N*i+11,4) = (*cn)(i,20);
        
        (*co)(N*i+12,1) = (*cn)(i,7);
        (*co)(N*i+12,2) = (*cn)(i,19);
        (*co)(N*i+12,3) = (*cn)(i,15);
        (*co)(N*i+12,4) = (*cn)(i,18);
        
        (*co)(N*i+13,1) = (*cn)(i,15);
        (*co)(N*i+13,2) = (*cn)(i,19);
        (*co)(N*i+13,3) = (*cn)(i,16);
        (*co)(N*i+13,4) = (*cn)(i,20);
        
        (*co)(N*i+14,1) = (*cn)(i,18);
        (*co)(N*i+14,2) = (*cn)(i,19);
        (*co)(N*i+14,3) = (*cn)(i,20);
        (*co)(N*i+14,4) = (*cn)(i,15);
        
        (*co)(N*i+15,1) = (*cn)(i,5);
        (*co)(N*i+15,2) = (*cn)(i,17);
        (*co)(N*i+15,3) = (*cn)(i,20);
        (*co)(N*i+15,4) = (*cn)(i,13);
        
        (*co)(N*i+16,1) = (*cn)(i,6);
        (*co)(N*i+16,2) = (*cn)(i,18);
        (*co)(N*i+16,3) = (*cn)(i,17);
        (*co)(N*i+16,4) = (*cn)(i,14);
        
        
        (*co)(N*i+17,1) = (*cn)(i,13);
        (*co)(N*i+17,2) = (*cn)(i,17);
        (*co)(N*i+17,3) = (*cn)(i,14);
        (*co)(N*i+17,4) = (*cn)(i,12);
        
        (*co)(N*i+18,1) = (*cn)(i,20);
        (*co)(N*i+18,2) = (*cn)(i,18);
        (*co)(N*i+18,3) = (*cn)(i,17);
        (*co)(N*i+18,4) = (*cn)(i,10);
        
        (*co)(N*i+19,1) = (*cn)(i,14);
        (*co)(N*i+19,2) = (*cn)(i,18);
        (*co)(N*i+19,3) = (*cn)(i,10);
        (*co)(N*i+19,4) = (*cn)(i,17);
        
        (*co)(N*i+20,1) = (*cn)(i,12);
        (*co)(N*i+20,2) = (*cn)(i,20);
        (*co)(N*i+20,3) = (*cn)(i,13);
        (*co)(N*i+20,4) = (*cn)(i,17);
        
        (*co)(N*i+21,1) = (*cn)(i,10);
        (*co)(N*i+21,2) = (*cn)(i,20);
        (*co)(N*i+21,3) = (*cn)(i,12);
        (*co)(N*i+21,4) = (*cn)(i,17); 
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 7, "HEXA_27") == 0)
    {
      E_Int N = 8;
      E_Int neltsF = nelts*N;
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "HEXA", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(N*i,1) = (*cn)(i,4);
        (*co)(N*i,2) = (*cn)(i,12);
        (*co)(N*i,3) = (*cn)(i,21);
        (*co)(N*i,4) = (*cn)(i,11);
        (*co)(N*i,5) = (*cn)(i,16);
        (*co)(N*i,6) = (*cn)(i,25);
        (*co)(N*i,7) = (*cn)(i,27);
        (*co)(N*i,8) = (*cn)(i,24);

        (*co)(N*i+1,1) = (*cn)(i,12);
        (*co)(N*i+1,2) = (*cn)(i,1);
        (*co)(N*i+1,3) = (*cn)(i,9);
        (*co)(N*i+1,4) = (*cn)(i,21);
        (*co)(N*i+1,5) = (*cn)(i,25);
        (*co)(N*i+1,6) = (*cn)(i,13);
        (*co)(N*i+1,7) = (*cn)(i,22);
        (*co)(N*i+1,8) = (*cn)(i,27);
 
      
        (*co)(N*i+2,1) = (*cn)(i,11);
        (*co)(N*i+2,2) = (*cn)(i,21);
        (*co)(N*i+2,3) = (*cn)(i,10);
        (*co)(N*i+2,4) = (*cn)(i,3);
        (*co)(N*i+2,5) = (*cn)(i,24);
        (*co)(N*i+2,6) = (*cn)(i,27);
        (*co)(N*i+2,7) = (*cn)(i,23);
        (*co)(N*i+2,8) = (*cn)(i,15);
        
        (*co)(N*i+3,1) = (*cn)(i,21);
        (*co)(N*i+3,2) = (*cn)(i,9);
        (*co)(N*i+3,3) = (*cn)(i,2);
        (*co)(N*i+3,4) = (*cn)(i,10);
        (*co)(N*i+3,5) = (*cn)(i,27);
        (*co)(N*i+3,6) = (*cn)(i,22);
        (*co)(N*i+3,7) = (*cn)(i,14);
        (*co)(N*i+3,8) = (*cn)(i,23);
        
        (*co)(N*i+4,1) = (*cn)(i,16);
        (*co)(N*i+4,2) = (*cn)(i,25);
        (*co)(N*i+4,3) = (*cn)(i,27);
        (*co)(N*i+4,4) = (*cn)(i,24);
        (*co)(N*i+4,5) = (*cn)(i,8);
        (*co)(N*i+4,6) = (*cn)(i,20);
        (*co)(N*i+4,7) = (*cn)(i,26);
        (*co)(N*i+4,8) = (*cn)(i,19);

        (*co)(N*i+5,1) = (*cn)(i,25);
        (*co)(N*i+5,2) = (*cn)(i,13);
        (*co)(N*i+5,3) = (*cn)(i,22);
        (*co)(N*i+5,4) = (*cn)(i,27);
        (*co)(N*i+5,5) = (*cn)(i,20);
        (*co)(N*i+5,6) = (*cn)(i,5);
        (*co)(N*i+5,7) = (*cn)(i,17);
        (*co)(N*i+5,8) = (*cn)(i,26);
        
        (*co)(N*i+6,1) = (*cn)(i,24);
        (*co)(N*i+6,2) = (*cn)(i,27);
        (*co)(N*i+6,3) = (*cn)(i,23);
        (*co)(N*i+6,4) = (*cn)(i,15);
        (*co)(N*i+6,5) = (*cn)(i,19);
        (*co)(N*i+6,6) = (*cn)(i,26);
        (*co)(N*i+6,7) = (*cn)(i,18);
        (*co)(N*i+6,8) = (*cn)(i,7);
      
        (*co)(N*i+7,1) = (*cn)(i,27);
        (*co)(N*i+7,2) = (*cn)(i,22);
        (*co)(N*i+7,3) = (*cn)(i,14);
        (*co)(N*i+7,4) = (*cn)(i,23);
        (*co)(N*i+7,5) = (*cn)(i,26);
        (*co)(N*i+7,6) = (*cn)(i,17);
        (*co)(N*i+7,7) = (*cn)(i,6);
        (*co)(N*i+7,8) = (*cn)(i,18);
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 8, "PENTA_18") == 0)
    {
      E_Int N = 8;
      E_Int neltsF = nelts*N;
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "PENTA", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(N*i,1) = (*cn)(i,1);
        (*co)(N*i,2) = (*cn)(i,7);
        (*co)(N*i,3) = (*cn)(i,9);
        (*co)(N*i,4) = (*cn)(i,10);
        (*co)(N*i,5) = (*cn)(i,16);
        (*co)(N*i,6) = (*cn)(i,18);
  
        (*co)(N*i+1,1) = (*cn)(i,7);
        (*co)(N*i+1,2) = (*cn)(i,2);
        (*co)(N*i+1,3) = (*cn)(i,8);
        (*co)(N*i+1,4) = (*cn)(i,16);
        (*co)(N*i+1,5) = (*cn)(i,11);
        (*co)(N*i+1,6) = (*cn)(i,17);
        
      
        (*co)(N*i+2,1) = (*cn)(i,7);
        (*co)(N*i+2,2) = (*cn)(i,8);
        (*co)(N*i+2,3) = (*cn)(i,9);
        (*co)(N*i+2,4) = (*cn)(i,16);
        (*co)(N*i+2,5) = (*cn)(i,17);
        (*co)(N*i+2,6) = (*cn)(i,18);
        
        (*co)(N*i+3,1) = (*cn)(i,9);
        (*co)(N*i+3,2) = (*cn)(i,8);
        (*co)(N*i+3,3) = (*cn)(i,3);
        (*co)(N*i+3,4) = (*cn)(i,18);
        (*co)(N*i+3,5) = (*cn)(i,17);
        (*co)(N*i+3,6) = (*cn)(i,12);
        
        (*co)(N*i+4,1) = (*cn)(i,10);
        (*co)(N*i+4,2) = (*cn)(i,16);
        (*co)(N*i+4,3) = (*cn)(i,18);
        (*co)(N*i+4,4) = (*cn)(i,4);
        (*co)(N*i+4,5) = (*cn)(i,13);
        (*co)(N*i+4,6) = (*cn)(i,15);
        
        (*co)(N*i+5,1) = (*cn)(i,16);
        (*co)(N*i+5,2) = (*cn)(i,11);
        (*co)(N*i+5,3) = (*cn)(i,17);
        (*co)(N*i+5,4) = (*cn)(i,13);
        (*co)(N*i+5,5) = (*cn)(i,5);
        (*co)(N*i+5,6) = (*cn)(i,14);
        
        (*co)(N*i+6,1) = (*cn)(i,16);
        (*co)(N*i+6,2) = (*cn)(i,17);
        (*co)(N*i+6,3) = (*cn)(i,18);
        (*co)(N*i+6,4) = (*cn)(i,13);
        (*co)(N*i+6,5) = (*cn)(i,14);
        (*co)(N*i+6,6) = (*cn)(i,15);
        
        (*co)(N*i+7,1) = (*cn)(i,18);
        (*co)(N*i+7,2) = (*cn)(i,17);
        (*co)(N*i+7,3) = (*cn)(i,12);
        (*co)(N*i+7,4) = (*cn)(i,15);
        (*co)(N*i+7,5) = (*cn)(i,14);
        (*co)(N*i+7,6) = (*cn)(i,6);      
      }
      RELEASESHAREDU(o, fo, co);   
    }
    else if (K_STRING::cmp(eltType, 7, "PYRA_14") == 0)
    {
      E_Int N = 16;
      E_Int neltsF = nelts*N;
      o = K_ARRAY::buildArray3(nfld, varString, nvertex, neltsF, "TETRA", false, api);
      FldArrayF* fo; FldArrayI* co;
      K_ARRAY::getFromArray3(o, fo, co);
      fo->copy(*f, 1, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        (*co)(N*i,1) = (*cn)(i,1);
        (*co)(N*i,2) = (*cn)(i,6);
        (*co)(N*i,3) = (*cn)(i,9);
        (*co)(N*i,4) = (*cn)(i,10);
        
        (*co)(N*i+1,1) = (*cn)(i,6);
        (*co)(N*i+1,2) = (*cn)(i,2);
        (*co)(N*i+1,3) = (*cn)(i,7);
        (*co)(N*i+1,4) = (*cn)(i,11);
      
        (*co)(N*i+2,1) = (*cn)(i,7);
        (*co)(N*i+2,2) = (*cn)(i,3);
        (*co)(N*i+2,3) = (*cn)(i,8);
        (*co)(N*i+2,4) = (*cn)(i,12);
        
        (*co)(N*i+3,1) = (*cn)(i,8);
        (*co)(N*i+3,2) = (*cn)(i,4);
        (*co)(N*i+3,3) = (*cn)(i,9);
        (*co)(N*i+3,4) = (*cn)(i,13);
        

        (*co)(N*i+4,1) = (*cn)(i,6);
        (*co)(N*i+4,2) = (*cn)(i,14);
        (*co)(N*i+4,3) = (*cn)(i,9);
        (*co)(N*i+4,4) = (*cn)(i,10);
        
        (*co)(N*i+5,1) = (*cn)(i,6);
        (*co)(N*i+5,2) = (*cn)(i,7);
        (*co)(N*i+5,3) = (*cn)(i,14);
        (*co)(N*i+5,4) = (*cn)(i,11);
        
        (*co)(N*i+6,1) = (*cn)(i,7);
        (*co)(N*i+6,2) = (*cn)(i,8);
        (*co)(N*i+6,3) = (*cn)(i,14);
        (*co)(N*i+6,4) = (*cn)(i,12);
        

        (*co)(N*i+7,1) = (*cn)(i,8);
        (*co)(N*i+7,2) = (*cn)(i,9);
        (*co)(N*i+7,3) = (*cn)(i,14);
        (*co)(N*i+7,4) = (*cn)(i,13);

        (*co)(N*i+8,1) = (*cn)(i,6);
        (*co)(N*i+8,2) = (*cn)(i,14);
        (*co)(N*i+8,3) = (*cn)(i,10);
        (*co)(N*i+8,4) = (*cn)(i,11);

        (*co)(N*i+9,1) = (*cn)(i,7);
        (*co)(N*i+9,2) = (*cn)(i,14);
        (*co)(N*i+9,3) = (*cn)(i,11);
        (*co)(N*i+9,4) = (*cn)(i,12);

        (*co)(N*i+10,1) = (*cn)(i,8);
        (*co)(N*i+10,2) = (*cn)(i,14);
        (*co)(N*i+10,3) = (*cn)(i,12);
        (*co)(N*i+10,4) = (*cn)(i,13);

        (*co)(N*i+11,1) = (*cn)(i,9);
        (*co)(N*i+11,2) = (*cn)(i,14);
        (*co)(N*i+11,3) = (*cn)(i,10);
        (*co)(N*i+11,4) = (*cn)(i,13);

        
        (*co)(N*i+12,1) = (*cn)(i,14);
        (*co)(N*i+12,2) = (*cn)(i,10);
        (*co)(N*i+12,3) = (*cn)(i,11);
        (*co)(N*i+12,4) = (*cn)(i,12);

        (*co)(N*i+13,1) = (*cn)(i,14);
        (*co)(N*i+13,2) = (*cn)(i,10);
        (*co)(N*i+13,3) = (*cn)(i,12);
        (*co)(N*i+13,4) = (*cn)(i,13);

        (*co)(N*i+14,1) = (*cn)(i,10);
        (*co)(N*i+14,2) = (*cn)(i,11);
        (*co)(N*i+14,3) = (*cn)(i,12);
        (*co)(N*i+14,4) = (*cn)(i,5);

        (*co)(N*i+15,1) = (*cn)(i,10);
        (*co)(N*i+15,2) = (*cn)(i,12);
        (*co)(N*i+15,3) = (*cn)(i,13);
        (*co)(N*i+15,4) = (*cn)(i,5);

      
      }
      RELEASESHAREDU(o, fo, co);   
    }

  }

  RELEASESHAREDB(res, array, f, cn);
  return o;
}
