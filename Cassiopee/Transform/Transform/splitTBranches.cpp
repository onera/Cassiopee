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

# include "transform.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Split a BAR array if branches are detected */
//=============================================================================
PyObject* K_TRANSFORM::splitTBranches(PyObject* self, PyObject* args)
{
  PyObject*  array;
  E_Float eps;
  if (!PYPARSETUPLE_(args, O_ R_,
                    &array, &eps))
  {
      return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType); 

  if (res == 1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "splitTBranches: cannot be used on a structured array.");
    return NULL;
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitTBranches: unknown type of array.");
    return NULL;
  }
  
  if (strcmp(eltType, "BAR") != 0)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitTBranches: must be used on a BAR-array.");
    return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); 
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitTBranches: cannot find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  K_CONNECT::cleanConnectivity(posx, posy, posz, eps, "BAR", *f, *cn);

  //determination des vertices de split
  E_Int npts = f->getSize(); E_Int nfld = f->getNfld();
  vector< vector<E_Int> > cVE(npts);
  K_CONNECT::connectEV2VE(*cn, cVE);
  vector<E_Int> splitVertices;
  for (E_Int nov = 0; nov < npts; nov++)
  {
    E_Int neltsVoisins = cVE[nov].size();
    if ( neltsVoisins > 2 ) // T-Branch detected
      splitVertices.push_back(nov);
  }

  if (splitVertices.size() == 0) 
  {
    PyObject* tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
    delete f; delete cn;
    return tpl;
  }

  vector<FldArrayF*> fields;
  vector<FldArrayI*> cnt;
  // construction des BAR splittes au niveau des ramifications
  E_Int nelts = cn->getSize();
  FldArrayI dejaVu(nelts); dejaVu.setAllValuesAtNull();
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Int v1, v2, etstart, et0, nop, noet;
  FldArrayF* fnew; FldArrayI* cnnew;
  // initialisation
  E_Int vnext = 0;
  E_Int vstart = 0; E_Int novstart = 0;
  E_Int nsplitPts = splitVertices.size();
  while (novstart < nsplitPts)
  {
    vstart = splitVertices[novstart];
    vector<E_Int>& eltsVoisins = cVE[vstart];
    E_Int nvoisins = eltsVoisins.size();
    etstart = -1;
    for (E_Int noet0 = 0; noet0 < nvoisins; noet0++)
    {
      et0 = eltsVoisins[noet0];
      if (dejaVu[et0] == 0) { etstart = et0; break; }
    }
    if (etstart == -1) { novstart++; goto end; }//aller au sommet multiple suivant
    fnew = new FldArrayF(2*npts, nfld);
    cnnew = new FldArrayI(nelts,2);
    nop = 0; noet = 0;
    vnext = vstart;
    next:;
    v1 = vnext;
    if (cn1[etstart]-1 == vnext) v2 = cn2[etstart]-1;
    else if (cn2[etstart]-1 == vnext) v2 = cn1[etstart]-1;
    else {printf("Error: splitTBranches: Pass.\n"); goto end;}

    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      (*fnew)(nop,eq) = (*f)(v1,eq);(*fnew)(nop+1,eq) = (*f)(v2,eq); 
      (*cnnew)(noet,1) = nop+1; (*cnnew)(noet,2) = nop+2;
    }
    dejaVu[etstart] = 1; 
    nop = nop+2; noet++;

    vnext = v2;
    if (cVE[vnext].size() == 2 && vnext != vstart) 
    {
      vector<E_Int>& eltsVoisins2 = cVE[vnext];
      E_Int nvoisins2 = eltsVoisins2.size();
      etstart = -1;
      for (E_Int noet0 = 0; noet0 < nvoisins2; noet0++)
      {
        et0 = eltsVoisins2[noet0];
        if (dejaVu[et0] == 0) etstart = et0;
      }
      goto next;
    }
    else
    {      
      fnew->reAllocMat(nop,nfld); cnnew->reAllocMat(noet,2);
      fields.push_back(fnew); cnt.push_back(cnnew);
    }
    end:;
  }
  // Sortie
  E_Int nbars = fields.size();
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  for (E_Int i = 0; i < nbars; ++i)
  {
    K_CONNECT::cleanConnectivity(posx, posy, posz, eps, "BAR",*fields[i],*cnt[i]);
    tpl = K_ARRAY::buildArray(*fields[i], varString, *cnt[i], -1, "BAR");
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  for (E_Int i = 0; i < nbars; i++)
  {delete fields[i]; delete cnt[i];}
  delete f; delete cn;
  return l;
}
