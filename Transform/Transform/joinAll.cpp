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
# include "transform.h"
# include <stdio.h>

using namespace std;
using namespace K_FLD; 
using namespace K_FUNC;
using namespace K_CONST;

// ============================================================================
/* Join all arrays */
// ============================================================================
PyObject* K_TRANSFORM::joinAll(PyObject* self, PyObject* args)
{
  PyObject* arrays; E_Float tol;
  if (!PYPARSETUPLE_(args, O_ R_,
                    &arrays, &tol))
  {
      return NULL;
  }

  // Check arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstructVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstructF;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayI*> cn; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  K_ARRAY::getFromArrays(arrays, res, structVarString,
                         unstructVarString, structF,
                         unstructF, ni, nj, nk,
                         cn, eltType, objs, obju, 
                         true, true, true, false, true);

  // Fusion des zones non-structures
  E_Int nu = unstructF.size();
  char* eltRef = NULL;
  if (nu > 0) eltRef = eltType[0];
  E_Int missed = 0;
  E_Int size = 0;
  E_Int ne = 0;
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0)
    {
      size += unstructF[i]->getSize();
      ne += cn[i]->getSize();
    }
    else missed++;
  }
  if (strcmp(eltRef, "NGON") == 0) ne -= 4*(nu-1);

  if (missed > 0) printf("Warning: joinAll: some arrays are not joined: different element types.\n");
  E_Int nfld = unstructF[0]->getNfld();
  FldArrayF field(size, nfld);

  E_Int c = 0; E_Int csav = 0;
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0)
    {
      FldArrayF& f2 = *unstructF[i];
      csav = c;
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* f2n = f2.begin(n);
        E_Float* fn = field.begin(n);
        c = csav;
        for (E_Int j = 0; j < f2.getSize(); j++)
        {
          fn[c] = f2n[j]; c++;
        }
      }
    }
  }

  // Fusion des connectivites
  E_Int nt = cn[0]->getNfld();
  FldArrayI cno(ne, nt);
  E_Int* cnop = cno.begin();

  // pour les NGON, il faut calculer les decalages avant
  E_Int nfacesTot = 0; // nbre total de faces
  E_Int sizeOfFNTot = 0; // sizeOfFN tot
  E_Int neltsTot = 0; // nbre d'elements tot
  E_Int sizeOfEFTot = 0; // sizeOfEF tot
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0 && strcmp(eltRef, "NGON") == 0)
    {
      E_Int* cnp = cn[i]->begin();
      E_Int nfaces = cnp[0];
      E_Int sizeOfFN = cnp[1];
      E_Int nelts = cnp[sizeOfFN+2];
      E_Int sizeOfEF = cnp[sizeOfFN+3];
      nfacesTot += nfaces;
      neltsTot += nelts;
      sizeOfFNTot += sizeOfFN;
      sizeOfEFTot += sizeOfEF;
    }
  }
  if (strcmp(eltRef, "NGON") == 0)
  {
    cnop[0] = nfacesTot;
    cnop[1] = sizeOfFNTot;
    cnop[sizeOfFNTot+2] = neltsTot;
    cnop[sizeOfFNTot+3] = sizeOfEFTot;
  }
  E_Int nf1 = 0; // decalage des noeuds
  c = 0; csav = 0; 
  E_Int nfa = 0; // decalage des faces
  E_Int* pt1 = cnop+2;
  E_Int* pt3 = cnop+sizeOfFNTot+4;
  
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0 && strcmp(eltRef, "NGON") == 0)
    {
      E_Int* cnp = cn[i]->begin();
      E_Int nfaces = cnp[0];
      E_Int sizeOfFN = cnp[1];
      E_Int nelts = cnp[sizeOfFN+2];
      //E_Int sizeOfEF = cnp[sizeOfFN+3];

      E_Int* pt2 = cnp+2;
      E_Int* pt4 = cnp+sizeOfFN+4;
      E_Int nv;
      for (E_Int j = 0; j < nfaces; j++)
      {
        nv = pt2[0];
        pt1[0] = nv;
        for (E_Int k = 1; k <= nv; k++) pt1[k] = pt2[k]+nf1; 
        pt1 += nv+1; pt2 += nv+1;
      }
      for (E_Int j = 0; j < nelts; j++)
      {
        nv = pt4[0];
        pt3[0] = nv;
        for (E_Int k = 1; k <= nv; k++) pt3[k] = pt4[k]+nfa; 
        pt3 += nv+1; pt4 += nv+1;
      }
      nfa += nfaces;
      nf1 += unstructF[i]->getSize();
    }
    else if (strcmp(eltRef, eltType[i]) == 0) // other elements
    {
      csav = c;
      for (E_Int n = 1; n <= nt; n++)
      {
        E_Int* cn1n = cn[i]->begin(n);
        E_Int* cnn = cno.begin(n);   
        c = csav;
        for (E_Int j = 0; j < cn[i]->getSize(); j++)
        {
          cnn[c] = cn1n[j]+nf1; c++;
        }
      }
      nf1 += unstructF[i]->getSize();
    }
  }

  for (E_Int i = 0; i < nu; i++)
    RELEASESHAREDU(obju[i], unstructF[i], cn[i]);

  E_Int posx = K_ARRAY::isCoordinateXPresent(unstructVarString[0])+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(unstructVarString[0])+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(unstructVarString[0])+1;
  K_CONNECT::cleanConnectivity(posx, posy, posz, tol, eltRef, 
                               field, cno);

  PyObject* tpl = K_ARRAY::buildArray(field, unstructVarString[0], 
                                      cno, -1, eltRef);
  return tpl;

}
// ============================================================================
/* Join all arrays and their arrays at centers */
// ============================================================================
PyObject* K_TRANSFORM::joinAllBoth(PyObject* self, PyObject* args)
{
  PyObject *arrays, *arraysc;
  E_Float tol;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &arrays, &arraysc, &tol))
  {
      return NULL;
  }

  // Check arrays for fields located at nodes 
  vector<E_Int> res;
  vector<char*> structVarString;
  vector<char*> unstructVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstructF;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayI*> cn;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  K_ARRAY::getFromArrays(arrays, res, structVarString,
			 unstructVarString, structF,
			 unstructF, ni, nj, nk,
			 cn, eltType, objs, obju, 
			 true, true, true, false, true);

  // Check arrays for fields located at centers
  vector<E_Int> resc;
  vector<char*> structVarStringc;
  vector<char*> unstructVarStringc;
  vector<FldArrayF*> structFc;
  vector<FldArrayF*> unstructFc;
  vector<E_Int> nic;
  vector<E_Int> njc; 
  vector<E_Int> nkc;
  vector<FldArrayI*> cnc;
  vector<char*> eltTypec;
  vector<PyObject*> objsc, objuc;
  K_ARRAY::getFromArrays(arraysc, resc, structVarStringc,
			 unstructVarStringc, structFc,
			 unstructFc, nic, njc, nkc,
			 cnc, eltTypec, objsc, objuc, 
			 false, false, false, false, true);

  // Fusion des zones non-structures
  E_Int nu = unstructF.size(); E_Int nuc = unstructFc.size();
  if (nu != nuc)
  {
    for (E_Int i = 0; i < nu; i++) RELEASESHAREDU(obju[i], unstructF[i], cn[i]);
    for (E_Int i = 0; i < nuc; i++) RELEASESHAREDU(objuc[i], unstructFc[i], cnc[i]);  
    PyErr_SetString(PyExc_ValueError,
                    "joinAllBoth: number of arrays at nodes and centers must be equal.");
    return NULL;
  }
  char* eltRef = NULL; char* eltRefc = NULL;
  eltRef = eltType[0]; eltRefc = eltTypec[0];

  E_Int missed = 0; E_Int missedc = 0; 
  E_Int size = 0; E_Int sizec = 0;
  E_Int ne = 0; E_Int nec = 0;
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0)
    {
      size += unstructF[i]->getSize();
      ne += cn[i]->getSize();
    }
    else missed++;

    if (strcmp(eltRefc, eltTypec[i]) == 0)
    {
      sizec += unstructFc[i]->getSize();
      nec += cnc[i]->getSize();
    }
    else missedc++;
  }
  if (missed > 0 || missedc > 0) printf("Warning: joinAllBoth: some arrays are not joined: different element types.\n");
  if (strcmp(eltRef, "NGON") == 0) ne -= 4*(nu-1);

  // Fusion des connectivites
  /*
  E_Int nt = cn[0]->getNfld();
  FldArrayI* cno = new FldArrayI(ne, nt);
  E_Int c = 0; E_Int csav = 0;
  E_Int nf1 = 0;
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0)
    {
      csav = c;
      for (E_Int n = 1; n <= nt; n++)
      {
        E_Int* cn1n = cn[i]->begin(n);
        E_Int* cnn = cno->begin(n);        
        c = csav;
        for (E_Int j = 0; j < cn[i]->getSize(); j++)
        {
          cnn[c] = cn1n[j]+nf1; c++;
        }
      }
      nf1 += unstructF[i]->getSize();
    }
  }
  */

  // Fusion des connectivites
  E_Int c = 0; E_Int csav = 0;
  E_Int nt = cn[0]->getNfld();
  FldArrayI* cno = new FldArrayI(ne, nt);
  E_Int* cnop = cno->begin();

  // pour les NGON, il faut calculer les decalages avant
  E_Int nfacesTot = 0; // nbre total de faces
  E_Int sizeOfFNTot = 0; // sizeOfFN tot
  E_Int neltsTot = 0; // nbre d'elements tot
  E_Int sizeOfEFTot = 0; // sizeOfEF tot
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0 && strcmp(eltRef, "NGON") == 0)
    {
      E_Int* cnp = cn[i]->begin();
      E_Int nfaces = cnp[0];
      E_Int sizeOfFN = cnp[1];
      E_Int nelts = cnp[sizeOfFN+2];
      E_Int sizeOfEF = cnp[sizeOfFN+3];
      nfacesTot += nfaces;
      neltsTot += nelts;
      sizeOfFNTot += sizeOfFN;
      sizeOfEFTot += sizeOfEF;
    }
  }
  if (strcmp(eltRef, "NGON") == 0)
  {
    cnop[0] = nfacesTot;
    cnop[1] = sizeOfFNTot;
    cnop[sizeOfFNTot+2] = neltsTot;
    cnop[sizeOfFNTot+3] = sizeOfEFTot;
  }
  E_Int nf1 = 0; // decalage des noeuds
  c = 0; csav = 0; 
  E_Int nfa = 0; // decalage des faces
  E_Int* pt1 = cnop+2;
  E_Int* pt3 = cnop+sizeOfFNTot+4;
  
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0 && strcmp(eltRef, "NGON") == 0)
    {
      E_Int* cnp = cn[i]->begin();
      E_Int nfaces = cnp[0];
      E_Int sizeOfFN = cnp[1];
      E_Int nelts = cnp[sizeOfFN+2];
      //E_Int sizeOfEF = cnp[sizeOfFN+3];

      E_Int* pt2 = cnp+2;
      E_Int* pt4 = cnp+sizeOfFN+4;
      E_Int nv;
      for (E_Int j = 0; j < nfaces; j++)
      {
        nv = pt2[0];
        pt1[0] = nv;
        for (E_Int k = 1; k <= nv; k++) pt1[k] = pt2[k]+nf1; 
        pt1 += nv+1; pt2 += nv+1;
      }
      for (E_Int j = 0; j < nelts; j++)
      {
        nv = pt4[0];
        pt3[0] = nv;
        for (E_Int k = 1; k <= nv; k++) pt3[k] = pt4[k]+nfa; 
        pt3 += nv+1; pt4 += nv+1;
      }
      nfa += nfaces;
      nf1 += unstructF[i]->getSize();
    }
    else if (strcmp(eltRef, eltType[i]) == 0) // other elements
    {
      csav = c;
      for (E_Int n = 1; n <= nt; n++)
      {
        E_Int* cn1n = cn[i]->begin(n);
        E_Int* cnn = cno->begin(n);   
        c = csav;
        for (E_Int j = 0; j < cn[i]->getSize(); j++)
        {
          cnn[c] = cn1n[j]+nf1; c++;
        }
      }
      nf1 += unstructF[i]->getSize();
    }
  }

  // join des champs en noeuds
  E_Int nfld = unstructF[0]->getNfld();
  FldArrayF* field = new FldArrayF(size, nfld);
  c = 0; csav = 0;
  for (E_Int i = 0; i < nu; i++)
  {
    if (strcmp(eltRef, eltType[i]) == 0)
    {
      csav = c;
      E_Int nf2 = unstructF[i]->getSize();
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* f2n = unstructF[i]->begin(n);
        E_Float* fn = field->begin(n);
        c = csav;
        for (E_Int j = 0; j < nf2; j++)
        {fn[c] = f2n[j]; c++;}
      }
    }
  }

  // join des champs en centres
  E_Int nfldc = unstructFc[0]->getNfld();
  FldArrayF* fieldc = new FldArrayF(sizec, nfldc);
  c = 0; csav = 0;
  for (E_Int i = 0; i < nuc; i++)
  {
    if (strcmp(eltRefc, eltTypec[i]) == 0)
    {
      csav = c;
      E_Int nf2 = unstructFc[i]->getSize();
      for (E_Int n = 1; n <= nfldc; n++)
      {
        E_Float* f2n = unstructFc[i]->begin(n);
        E_Float* fn = fieldc->begin(n);
        c = csav;
        for (E_Int j = 0; j < nf2; j++)
        {fn[c] = f2n[j]; c++;}
      }
    }
  }
  for (E_Int i = 0; i < nu; i++)
  {
    RELEASESHAREDU(obju[i], unstructF[i], cn[i]);
    RELEASESHAREDU(objuc[i], unstructFc[i], cnc[i]);  
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(unstructVarString[0])+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(unstructVarString[0])+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(unstructVarString[0])+1;
  K_CONNECT::cleanConnectivity(posx, posy, posz, tol, eltRef, *field, *cno);
  PyObject* l = PyList_New(0);
  PyObject* tpl1 = K_ARRAY::buildArray(*field, unstructVarString[0], 
                                       *cno, -1, eltRef);
  PyList_Append(l, tpl1); Py_DECREF(tpl1); delete field;
  PyObject* tpl2 = K_ARRAY::buildArray(*fieldc, unstructVarStringc[0], 
                                       *cno, -1, eltRefc);
  PyList_Append(l, tpl2); Py_DECREF(tpl2); delete fieldc;
  delete cno; 
  return l;
}
