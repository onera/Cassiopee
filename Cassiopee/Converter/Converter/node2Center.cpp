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

// conversion noeuds/centres

#include <stdio.h>
#include <string.h>
#include "converter.h"
#include "kcore.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Convertit un tableau defini en noeuds en un tableau defini en centres */
// ============================================================================
PyObject* K_CONVERTER::node2Center(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int sorted;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &sorted)) return NULL;

  PyObject* tpl;
  E_Int ni, nj, nk;
  char* varString; char* eltType;
  E_Int res; 
  FldArrayF* FNode; FldArrayI* c;
  res = K_ARRAY::getFromArray3(array, varString, FNode, ni, nj, nk,
                               c, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "node2Center: array is invalid.");
    return NULL;
  }

  /* Essaie de trouver la variable cellN. Les traitements sont un peu
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
    E_Float* cellNat = FNode->begin(cellN);
    for (E_Int ind = 0; ind < FNode->getSize(); ind++)
    {
      nature = cellNat[ind];
      if (K_FUNC::fEqualZero(nature-2.) == true) { mod = 2; break; }
      else if (nature < -0.2) { mod = 3; break; }
    }
  }

  E_Int nfld = FNode->getNfld();
  E_Int api = FNode->getApi();
  FldArrayF* FCenter;
    
  if (res == 1)
  {
    E_Int nil=ni, njl=nj, nkl=nk;
    if (ni == 1 && nj == 1) { nil = nk-1; nkl = 1; }
    else if (ni == 1 && nk == 1) {nil = nj-1; njl = 1;}
    else
    {
      if (ni != 1) nil = ni-1;
      if (nj != 1) njl = nj-1;
      if (nk != 1) nkl = nk-1;
    }
    tpl = K_ARRAY::buildArray3(nfld, varString, nil, njl, nkl, api);
    K_ARRAY::getFromArray3(tpl, FCenter);

    ret = K_LOC::node2centerStruct(*FNode, ni, nj, nk, cellN, mod, *FCenter);
    RELEASESHAREDS(array, FNode);
    RELEASESHAREDS(tpl, FCenter);
    if (ret == 0) return NULL;
    return tpl;
  }
  else if (res == 2)
  {
    char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
    K_ARRAY::starVarString(eltType, eltType2);
    E_Int ncells = 0;
    E_Bool compact = false;
    if (api == 1) compact = true;
    if (strcmp(eltType, "NGON") == 0)
    {
      ncells = c->getNElts();
      FCenter = new FldArrayF(ncells, nfld, compact);
      ret = K_LOC::node2centerNGon(*FNode, *c, *FCenter, sorted);
    }
    else if (strcmp(eltType, "NODE") == 0)
    {
      ncells = FNode->getSize();
      FCenter = new FldArrayF(ncells, nfld, compact);
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fn = FNode->begin(n);
        E_Float* fc = FCenter-> begin(n); 
        for (E_Int i = 0; i < ncells; i++) fc[i] = fn[i];
      }
    }
    else // autres elements basiques
    { 
      E_Int nc = c->getNConnect();
      for (E_Int ic = 0; ic < nc; ic++)
      { 
        FldArrayI& cm = *(c->getConnect(ic));
        ncells += cm.getSize();
      }

      FCenter = new FldArrayF(ncells, nfld, compact);
      ret = K_LOC::node2centerUnstruct(*FNode, *c, cellN, mod, *FCenter);
    }
        
    if (ret == 0) 
    {
      delete FCenter; RELEASESHAREDU(array, FNode, c); 
      return NULL;
    }
    tpl = K_ARRAY::buildArray3(*FCenter, varString, *c, eltType2);
    delete[] eltType2; delete FCenter; RELEASESHAREDU(array, FNode, c);
    return tpl;
  }
  else return NULL;
}
