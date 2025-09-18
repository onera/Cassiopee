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
#include "converter.h"
//# include <vector>
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Convertit un tableau defini en noeuds en un tableau defini en centres */
// ============================================================================
PyObject* K_CONVERTER::center2Node2(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int type;
  PyObject* BCFields; // optional indR+fields on BCs
  if (!PYPARSETUPLE_(args, O_ I_ O_, &array, &type, &BCFields)) return NULL;

  // Check BCFields
  if (BCFields != Py_None)
  {
    if (PyList_Check(BCFields) == false) 
    {
      PyErr_SetString(PyExc_TypeError, 
                       "center2Node: BCFields must be a list of indR+fields.");
      return NULL;
    }
  }

  // printf("Converter: center2Node: Api 3\n");

  E_Int ni, nj, nk;
  FldArrayI* c; FldArrayF* FCenter;
  char* eltType; char* varString;
  E_Int res; 
  res = K_ARRAY::getFromArray3(array, varString, FCenter, 
                               ni, nj, nk, c, eltType);
  if (res != 1 && res != 2) return NULL;

  /* Essaie de trouver la variables cellN. Les traitements sont un peu
     differents pour cette variable. */
  E_Int cellN = -1;
  E_Int ret = K_ARRAY::isCellNatureField2Present(varString);
  if (ret != -1) cellN = ret+1;

  // Retourne le mode de sortie du champ cellnaturefield (1 (0,1), 2 (0,1,2) 
  // ou 3 (0, 1, -interpolationblock))
  E_Int mod = 0; E_Float nature;
  if (cellN != -1)
  {
    mod = 1;
    E_Float* cellNat = FCenter->begin(cellN);
    for (E_Int ind = 0; ind < FCenter->getSize(); ind++) 
    {
      nature = cellNat[ind];
      if (K_FUNC::fEqualZero(nature-2.) == true) { mod = 2; break; }
      else if (nature < -0.2) { mod = 3; break; }
    }
  }

  E_Int api = FCenter->getApi();

  if (res == 1)
  {
    E_Int nin, njn, nkn;
    E_Int nfld = FCenter->getNfld();
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    
    if (ni == 1)
    {
      if ((nj != 1)&&(nk != 1))
      { nin = 1; njn = nj+1; nkn = nk+1; }
      else if (nj == 1)
      { nin = 1; njn = 1; nkn = nk+1; }
      else //if (nk == 1)
      { nin = 1; njn = 1; nkn = nk+1; }
    }
    else if (nj == 1)
    {
      if ((ni != 1)&&(nk != 1))
      { nin = ni+1; njn = 1; nkn = nk+1; }
      else if (ni == 1)
      { nin = 1; njn = 1; nkn = nk+1; }
      else // if (nk == 1)
      { nin = ni+1; njn = 1; nkn = 1; }
    }
    else if (nk == 1)
    {
      if ((ni != 1)&&(nj != 1))
      { nin = ni+1; njn = nj+1; nkn = 1; }
      else if (ni == 1)
      { nin = 1; njn = nj+1; nkn = 1; }
      else //if (nj == 1)
      { nin = ni+1; njn = 1; nkn = 1; }
    }
    else
    { nin = ni+1; njn = nj+1; nkn = nk+1; }
    PyObject*tpl = K_ARRAY::buildArray3(nfld, varString, nin, njn, nkn);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF FNode(nin*njn*nkn, nfld, fnp, true);
    ret = K_LOC::center2nodeStruct2(*FCenter, ni, nj, nk, cellN, mod, 
                                   posx, posy, posz, FNode, nin, njn, nkn,
                                   type);
    // Boundary corrections
    if (BCFields != Py_None)
    {
      //PyObject* indR = PyList_GetItem(BCFields, 0);
      //PyObject* fields = PyList_GetItem(BCFields, 1);
      //E_Int res = K_ARRAY::getFromArray3(fields, varString, FCenter, 
      //                                   ni, nj, nk, c, eltType);
      //center2NodeStructBorder(FNode, nin, njn, nkn);
      //RELEASESHAREDB(res, fields, );
    }

    RELEASESHAREDS(array, FCenter);
    if (ret == 0) return NULL;
    return tpl;
  }
  else if (res == 2)
  {
    E_Int nb = 0; E_Int nfld = FCenter->getNfld();
    if (// K_STRING::cmp(eltType, "NGON") == 0 || 
        K_STRING::cmp(eltType, "NGON*") == 0)
    {
      E_Int nelts = FCenter->getSize();
      vector< vector<E_Int> > cEV(nelts);
      K_CONNECT::connectNG2EV(*c, cEV);

      /* Calcul de Nb en eliminant les vertex non references
      E_Int nptsmax = 0;
      for (E_Int et = 0; et < nelts; et++) nptsmax += cEV[et].size();
      FldArrayI count(nptsmax); count.setAllValuesAtNull();
      E_Int* countp = count.begin();
      for (E_Int et = 0; et < nelts; et++)
      {
        vector<E_Int>& vertices = cEV[et]; E_Int nvert = vertices.size();   
        for (E_Int nov = 0; nov < nvert; nov++)
        {
          E_Int indv = vertices[nov]-1;         
          if (countp[indv] == 0) {countp[indv]++; nb++;}
        }
      }
      count.malloc(0);
      */

      /* Calcul de Nb en prenant le numero max du vertex dans les faces */
      E_Int nb = 0;
      for (E_Int et = 0; et < nelts; et++)
      {
        vector<E_Int>& vertices = cEV[et]; E_Int nvert = vertices.size();   
        for (E_Int nov = 0; nov < nvert; nov++)
        {
          E_Int indv = vertices[nov];
          nb = K_FUNC::E_max(nb, indv);
        }
      }

      FldArrayF* FNode = new FldArrayF(nb, nfld);
      ret = K_LOC::center2nodeNGon2(*FCenter, *c, cEV, *FNode, cellN, mod, type);
    
      for (E_Int et = 0; et < nelts; et++) cEV[et].clear();
      cEV.clear();
      PyObject* tpl = K_ARRAY::buildArray3(*FNode, varString, *c, eltType, api);
      delete FNode;
      RELEASESHAREDU(array, FCenter, c);
      if (ret == 0) return NULL;
      return tpl;
    }
    else // elements basiques 
    {    
      FldArrayF* FNode;
      if (// K_STRING::cmp(eltType, "NODE") == 0 || 
          K_STRING::cmp(eltType, "NODE*") == 0)
      {
        nb = FCenter->getSize();
        FNode = new FldArrayF(nb, nfld);
        E_Float* fn = FNode->begin();
        E_Float* fc = FCenter->begin();
        E_Int nfld = FCenter->getNfld();
        for (E_Int i = 0; i < nb*nfld; i++) fn[i] = fc[i]; 
        ret = 1;
      }
      else if (K_STRING::cmp(eltType, "BAR*") == 0   || 
               K_STRING::cmp(eltType, "TRI*") == 0   ||
               K_STRING::cmp(eltType, "QUAD*") == 0  ||
               K_STRING::cmp(eltType, "TETRA*") == 0 ||
               K_STRING::cmp(eltType, "HEXA*") == 0  ||
               K_STRING::cmp(eltType, "PENTA*") == 0 ||
               K_STRING::cmp(eltType, "PYRA*") == 0 ) 
      {
        E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
        E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
        E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
        E_Int ne = c->getSize(); // nombre de centres = nombre d'elements
        E_Int nt = c->getNfld();
        for (E_Int n = 1; n <= nt; n++)
        {
          E_Int* cn = c->begin(n);
          for (E_Int e = 0; e < ne; e++) nb = K_FUNC::E_max(nb, cn[e]);
        }
        FNode = new FldArrayF(nb, nfld);      
        ret = K_LOC::center2nodeUnstruct2(*FCenter, *c, cellN, mod, posx, posy, posz, *FNode, type);
      }
      else // autres
      {
        PyErr_SetString(PyExc_TypeError, 
                        "center2Node: unstructured array must be eltType*.");
        RELEASESHAREDU(array, FCenter, c); return NULL;
      }
      PyObject* tpl = K_ARRAY::buildArray3(*FNode, varString, *c, eltType, api);
      delete FNode;
      RELEASESHAREDU(array, FCenter, c);
      if (ret == 0) return NULL;
      return tpl;
    }
  }
  else return NULL;
}

// center2Node for border
