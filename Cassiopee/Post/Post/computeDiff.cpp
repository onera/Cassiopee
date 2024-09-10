/*    
    Copyright 2013-2024 Onera.

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
# include "post.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Compute diff*/
//=============================================================================
PyObject* K_POST::computeDiff(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* varname;
  if (!PyArg_ParseTuple(args, "OO", &array, &varname)) return NULL;
   
  // Check array
  char* varString; char* eltType; 
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, 
                                    eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeDiff: invalid array.");
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
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "computeDiff: varname must be a string.");
    return NULL;
  }
   
  E_Int posv = K_ARRAY::isNamePresent(var, varString);
  if (posv == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "computeDiff: variable not found in array.");
    return NULL;
  }
  posv++;
  // check if cellN is in 
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString); posc++;

  E_Float* cellN = NULL;
  if (posc > 0) cellN = f->begin(posc);

  /*-----------------------------*/
  E_Float* field = f->begin(posv);
  E_Int dim = 3;

  E_Int ind, indm1, indp1;
  E_Float vali1, vali2, valj1, valj2, valk1, valk2;
  E_Int npts = f->getSize();
  E_Float find;
  E_Int ninj = ni*nj;

  if (res == 1) 
  {
    FldArrayF* diff = new FldArrayF(npts,1); diff->setAllValuesAtNull();
    E_Float* difft = diff->begin();
    if (ni < 3 || nj < 3 || nk < 3) dim = 2;
    if (posc < 1)
    {  
      switch (dim)
      {
        case 2:
          if (nk < 3) 
          {
            for (E_Int j = 0; j < nj; j++)
              for (E_Int i = 0; i < ni; i++)
              {
                ind = i + j * ni;
                find = field[ind];
                // direction i
                indm1 = ind-1; indp1 = ind+1;
                if (i == 0) indm1 = ind; 
                else if (i == ni-1) indp1 = ind;
                vali1 = K_FUNC::E_abs(find-field[indm1]);
                vali2 = K_FUNC::E_abs(find-field[indp1]);
                //direction j
                indm1 = ind-ni; indp1 = ind+ni;
                if (j == 0) indm1 = ind; 
                else if (j == nj-1) indp1 = ind;
                valj1 = K_FUNC::E_abs(find-field[indm1]);
                valj2 = K_FUNC::E_abs(find-field[indp1]);
                //diff
                difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                if (nk == 2) difft[ind+ninj] = difft[ind];
              }
          }
          else if (nj < 3) 
          {
            for (E_Int k = 0; k < nk; k++)
              for (E_Int i = 0; i < ni; i++)
              {
                ind = i + k*ni;
                find = field[ind];
                // direction i
                indm1 = ind-1; indp1 = ind+1;
                if (i == 0) indm1 = ind; 
                else if (i == ni-1) indp1 = ind;
                vali1 = K_FUNC::E_abs(find-field[indm1]);
                vali2 = K_FUNC::E_abs(find-field[indp1]);
                //direction k
                indm1 = ind-ni; indp1 = ind+ni;
                if (k == 0) indm1 = ind; 
                else if (k == nk-1) indp1 = ind;
                valj1 = K_FUNC::E_abs(find-field[indm1]);
                valj2 = K_FUNC::E_abs(find-field[indp1]);
                //diff
                difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                if (nj == 2) difft[ind+ni] = difft[ind];
              }
          }
          else if (ni < 3) 
          {
            for (E_Int k = 0; k < nk; k++)
              for (E_Int j = 0; j < nj; j++)
              {
                ind = j + k*ni;
                find = field[ind];
                // direction j
                indm1 = ind-1; indp1 = ind+1;
                if (j == 0) indm1 = ind; 
                else if (j == nj-1) indp1 = ind;
                vali1 = K_FUNC::E_abs(find-field[indm1]);
                vali2 = K_FUNC::E_abs(find-field[indp1]);
                //direction k
                indm1 = ind-ni; indp1 = ind+ni;
                if (k == 0) indm1 = ind; 
                else if (k == nk-1) indp1 = ind;
                valj1 = K_FUNC::E_abs(find-field[indm1]);
                valj2 = K_FUNC::E_abs(find-field[indp1]);
                //diff
                difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                if (ni == 2) difft[ind+1] = difft[ind];
              }
          }
          break;
        case 3:
          for (E_Int k = 0; k < nk; k++)
            for (E_Int j = 0; j < nj; j++)
              for (E_Int i = 0; i < ni; i++)            
              {
                ind = i + j*ni + k*ninj;
                find = field[ind];
                // direction i
                indm1 = ind-1; indp1 = ind+1;
                if (i == 0) indm1 = ind; 
                else if (i == ni-1) indp1 = ind;
                vali1 = K_FUNC::E_abs(find-field[indm1]);
                vali2 = K_FUNC::E_abs(find-field[indp1]);
                //direction j 
                indm1 = ind-ni; indp1 = ind+ni;
                if (j == 0) indm1 = ind; 
                else if (j == nj-1) indp1 = ind;
                valj1 = K_FUNC::E_abs(find-field[indm1]);
                valj2 = K_FUNC::E_abs(find-field[indp1]);
                //direction k
                indm1 = ind-ninj; indp1 = ind+ninj;
                if (k == 0) indm1 = ind; 
                else if (k == nk-1) indp1 = ind;
                valk1 = K_FUNC::E_abs(find-field[indm1]);
                valk2 = K_FUNC::E_abs(find-field[indp1]);
                difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                difft[ind] = K_FUNC::E_max(difft[ind], valk1, valk2);
              }
          break;
      }
    } // pas de cellN
    else // cellN
    {
      switch (dim)
      {
        case 2:
          if ( nk < 3 ) 
          {
            for (E_Int j = 0; j < nj; j++)
              for (E_Int i = 0; i < ni; i++)
              {
                ind = i + j * ni;
                find = field[ind];
                if (cellN[ind] == 1.) 
                {
                  // direction i
                  indm1 = ind-1; indp1 = ind+1;
                  if (i == 0) indm1 = ind; 
                  else if (i == ni-1) indp1 = ind;
                  vali1 = K_FUNC::E_abs(find-field[indm1]);
                  vali2 = K_FUNC::E_abs(find-field[indp1]);
                  //direction j
                  indm1 = ind-ni; indp1 = ind+ni;
                  if (j == 0) indm1 = ind; 
                  else if (j == nj-1) indp1 = ind;
                  valj1 = K_FUNC::E_abs(find-field[indm1]);
                  valj2 = K_FUNC::E_abs(find-field[indp1]);
                  //diff
                  difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                }
                if (nk == 2) difft[ind+ninj] = difft[ind];
              }
          }
          else if (nj < 3) 
          {
            for (E_Int k = 0; k < nk; k++)
              for (E_Int i = 0; i < ni; i++)
              {
                ind = i + k*ni;
                find = field[ind];
                if (cellN[ind] == 1.) 
                {
                  // direction i
                  indm1 = ind-1; indp1 = ind+1;
                  if (i == 0) indm1 = ind; 
                  else if (i == ni-1) indp1 = ind;
                  vali1 = K_FUNC::E_abs(find-field[indm1]);
                  vali2 = K_FUNC::E_abs(find-field[indp1]);
                  //direction k
                  indm1 = ind-ni; indp1 = ind+ni;
                  if (k == 0) indm1 = ind; 
                  else if (k == nk-1) indp1 = ind;
                  valj1 = K_FUNC::E_abs(find-field[indm1]);
                  valj2 = K_FUNC::E_abs(find-field[indp1]);
                  //diff
                  difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                }
                if (nj == 2) difft[ind+ni] = difft[ind];
              }
          }
          else if (ni < 3)
          {
            for (E_Int k = 0; k < nk; k++)
              for (E_Int j = 0; j < nj; j++)
              {
                ind = j + k*ni;
                find = field[ind];
                if (cellN[ind] == 1.) 
                {
                  // direction j
                  indm1 = ind-1; indp1 = ind+1;
                  if (j == 0) indm1 = ind; 
                  else if (j == nj-1) indp1 = ind;
                  vali1 = K_FUNC::E_abs(find-field[indm1]);
                  vali2 = K_FUNC::E_abs(find-field[indp1]);
                  //direction k
                  indm1 = ind-ni; indp1 = ind+ni;
                  if (k == 0) indm1 = ind; 
                  else if (k == nk-1) indp1 = ind;
                  valj1 = K_FUNC::E_abs(find-field[indm1]);
                  valj2 = K_FUNC::E_abs(find-field[indp1]);
                  //diff
                  difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                }
                if (ni == 2) difft[ind+1] = difft[ind];
              }
          }
          break;
        case 3:          
          for (E_Int k = 0; k < nk; k++)
            for (E_Int j = 0; j < nj; j++)
              for (E_Int i = 0; i < ni; i++)            
              {
                ind = i + j*ni + k*ninj;
                find = field[ind];
                if (cellN[ind] == 1.) 
                {
                  // direction i
                  indm1 = ind-1; indp1 = ind+1;
                  if (i == 0) indm1 = ind; 
                  else if (i == ni-1) indp1 = ind;
                  vali1 = K_FUNC::E_abs(find-field[indm1]);
                  vali2 = K_FUNC::E_abs(find-field[indp1]);
                  //direction j 
                  indm1 = ind-ni; indp1 = ind+ni;
                  if (j == 0) indm1 = ind; 
                  else if (j == nj-1) indp1 = ind;
                  valj1 = K_FUNC::E_abs(find-field[indm1]);
                  valj2 = K_FUNC::E_abs(find-field[indp1]);
                  //direction k
                  indm1 = ind-ninj; indp1 = ind+ninj;
                  if (k == 0) indm1 = ind; 
                  else if (k == nk-1) indp1 = ind;
                  valk1 = K_FUNC::E_abs(find-field[indm1]);
                  valk2 = K_FUNC::E_abs(find-field[indp1]);
                  difft[ind] = K_FUNC::E_max(vali1, vali2, valj1, valj2);
                  difft[ind] = K_FUNC::E_max(difft[ind], valk1, valk2);
                }
              }
          break;
      }      
    }
    RELEASESHAREDB(res, array, f, cn);
    PyObject* tpl = K_ARRAY::buildArray(*diff, var, ni, nj, nk);
    delete diff;
    return tpl;
  } // fin structure
  else // non structure
  {
    E_Int ok = 0; 
    FldArrayF* diff;
    // centres ou noeuds ?
    if (strcmp(eltType,"BAR") == 0 || strcmp(eltType,"TRI") == 0 || 
        strcmp(eltType,"QUAD") == 0 || strcmp(eltType,"HEXA") == 0 ||
        strcmp(eltType,"TETRA") == 0 || strcmp(eltType,"PENTA") == 0) 
    {
      diff = new FldArrayF(npts,1); diff->setAllValuesAtNull();
      if (posc > 0)
        computeDiffForVertexWithCellN(npts, *cn, field, cellN, diff->begin(1));
      else 
        computeDiffForVertex(npts, *cn, field, diff->begin(1));
      ok = 1;
    }    
    if (ok == 0) 
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError, 
                      "computeDiff: unknown element type for unstructured array.");
      return NULL;
    }
    FldArrayI* cn2 = new FldArrayI(cn->getSize(), cn->getNfld());
    *cn2 = *cn; //copy

    RELEASESHAREDU(array, f, cn);
    PyObject* tpl = K_ARRAY::buildArray(*diff, var, *cn2, -1, eltType);
    delete diff; delete cn2;
    return tpl;
  }
}
//=============================================================================
/* Calcul pour les elements - loc = vertex */
//=============================================================================
void K_POST::computeDiffForVertex(E_Int npts, FldArrayI& cn, 
                                  E_Float* field, E_Float* difft)
{
  vector< vector<E_Int> > cVN(npts);
  K_CONNECT::connectEV2VNbrs(cn, cVN);                          
  E_Float valmax; E_Float diff;
  for (E_Int ind = 0; ind < npts; ind++)
  {
    E_Float val = field[ind]; valmax = 0.;
    vector<E_Int>& voisins = cVN[ind];//vertices voisins
    E_Int nvoisins = voisins.size();
    for (E_Int nov = 0; nov < nvoisins; nov++)
    {
      E_Int v = voisins[nov]-1;
      diff = K_FUNC::E_abs(field[v]-val);
      valmax = K_FUNC::E_max(valmax, diff);
    }
    difft[ind] = valmax;
  }
}
//=============================================================================
/* Meme chose mais avec prise en compte du cellN */
//=============================================================================
void K_POST::computeDiffForVertexWithCellN(
  E_Int npts, FldArrayI& cn, E_Float* field, 
  E_Float* cellN, E_Float* difft)
{
  vector< vector<E_Int> > cVN(npts);
  K_CONNECT::connectEV2VNbrs(cn, cVN);                          
  E_Float valmax; E_Float diff;
  for (E_Int ind = 0; ind < npts; ind++)
  {
    if (cellN[ind] == 1.) 
    {
      E_Float val = field[ind]; valmax = 0.;
      vector<E_Int>& voisins = cVN[ind];//vertices voisins
      E_Int nvoisins = voisins.size();
      for (E_Int nov = 0; nov < nvoisins; nov++)
      {
        E_Int v = voisins[nov]-1;
        diff = K_FUNC::E_abs(field[v]-val);
        valmax = K_FUNC::E_max(valmax, diff);
      }
      difft[ind] = valmax;
    }
  }
}
