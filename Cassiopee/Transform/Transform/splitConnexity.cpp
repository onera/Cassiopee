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

#include "transform.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   splitConnexity: 
   Decoupe un array non structure en partie connexes de meme type d'element.
   La connectivite doit etre propre.
*/
//=============================================================================
PyObject* K_TRANSFORM::splitConnexity(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array))
  {
      return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitConnexity: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitConnexity: cannot be used on a structured array.");
    return NULL;
  }
  
  if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitConnexity: cannot be used on a PYRA-array.");
    return NULL;
  }
  else if (K_STRING::cmp(eltType, "PENTA") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitConnexity: cannot be used on a PENTA-array.");
    return NULL;
  }
  else if (K_STRING::cmp(eltType, "MIX") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitConnexity: cannot be used on a MIX-array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitConnexity: cannot find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  PyObject* out;
  if (K_STRING::cmp(eltType, "NGON") == 0) 
    out = splitConnexityNGon(f, cn, varString, posx, posy, posz);
  else if (K_STRING::cmp(eltType, "NODE") == 0) 
    out = splitConnexityNODE(f, cn, eltType, varString, posx, posy, posz);
  else out = splitConnexityBasics(f, cn, eltType, varString, 
                                  posx, posy, posz);

  RELEASESHAREDU(array, f, cn);
  return out;
}

//=============================================================================
PyObject* K_TRANSFORM::splitConnexityBasics(
  FldArrayF* f, FldArrayI* cn, 
  char* eltType, char* varString,
  E_Int posx, E_Int posy, E_Int posz)
{
  vector< vector<E_Int> > cEEN(cn->getSize());
  K_CONNECT::connectEV2EENbrs(eltType, f->getSize(), *cn, cEEN);
  
  E_Int nt = cn->getNfld();
  E_Int ne = cn->getSize(); // nbre d'elements
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ne, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ne * sizeof(E_Int));
  E_Int mbv, p, i, ie, elt, curr;
  unsigned int iv;
  vector<FldArrayI*> components;

  mbv = 0;

  while (nev < ne)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant connexe
    FldArrayI* c = new FldArrayI(ne, nt);
 
    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    curr = 0;
    
    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      for (i = 1; i <= nt; i++) (*c)(curr,i) = (*cn)(elt,i);
      curr++;

      for (iv = 0; iv < cEEN[elt].size(); iv++)
      {
        ie = cEEN[elt][iv];
        if (isVisited[ie] == 0)
        {
          mustBeVisited[mbv] = ie;
          mbv++; nev++;
          isVisited[ie] = 1;
        }
      }
    }
    c->reAllocMat(curr, nt);
    components.push_back(c);
  }

  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);

  E_Int size = components.size();
  for (i = 0; i < size; i++)
  {
    FldArrayF* f0 = new FldArrayF(*f);
    FldArrayF& fp = *f0;
    FldArrayI& cnp = *components[i];
    K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType,
                                 fp, cnp);
    tpl = K_ARRAY::buildArray(fp, varString, cnp, -1, eltType);
    delete &fp; delete &cnp;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}

//=============================================================================
PyObject* K_TRANSFORM::splitConnexityNGon(
  FldArrayF* f, FldArrayI* cn, char* varString,
  E_Int posx, E_Int posy, E_Int posz)
{
  E_Int* ptr = cn->begin();
  E_Int sf = ptr[1];
  E_Int ne = ptr[2+sf]; // nbre d'elements

  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  
  E_Int* ptrElts = ptr + (sf+4);
  E_Int se = ptr[3+sf]; // taille connectivite elements/faces
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ne, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ne * sizeof(E_Int));
  E_Int mbv, p, i, ie, elt, curr, necurr, lt, iv;
  vector<FldArrayI*> components; // connectivite EF locale
  FldArrayI pos; K_CONNECT::getPosElts(*cn, pos);
  E_Int e1, e2, nf;

  mbv = 0;
  while (nev < ne)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant connexe
    FldArrayI* c = new FldArrayI(se+1);
    E_Int* pc = c->begin(); // current pointer
 
    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    curr = 0; necurr = 0;
    
    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      // copie elt
      ptrElts = &ptr[pos[elt]];
      lt =  ptrElts[0];
      pc[0] = lt; 
      for (i = 1; i <= pc[0]; i++) pc[i] = ptrElts[i];
      pc += lt+1;
      curr += lt+1; necurr++;

      //for (iv = 0; iv < cEEN[elt].size(); iv++)
      //{
      //  ie = cEEN[elt][iv];
      for (iv = 0; iv < lt; iv++)
      {
        nf = ptrElts[iv+1]-1;
        e1 = cFE1[nf]-1; e2 = cFE2[nf]-1;
        if (e1 == elt) ie = e2;
        else ie = e1;

        if (ie != -1 && isVisited[ie] == 0)
        {
          mustBeVisited[mbv] = ie;
          mbv++; nev++;
          isVisited[ie] = 1;
        }
      }
    }
    
    pc[0] = necurr; // sentinelle
    c->reAlloc(curr+1);
    components.push_back(c);
  }

  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  E_Int size = components.size();
  //printf("size %d\n", size);

  for (i = 0; i < size; i++)
  {
    FldArrayF* f0 = new FldArrayF(*f);
    FldArrayF& fp = *f0;
    // nouvelle connectivite
    E_Int si = components[i]->getSize()-1;
    E_Int* comp = components[i]->begin();
    E_Int size = sf+4+si;
    FldArrayI cnp(size);
    E_Int* cnpp = cnp.begin();
    for (E_Int j = 0; j < sf+2; j++) cnpp[j] = ptr[j];
    cnpp += sf+2;
    cnpp[0] = comp[si]; cnpp[1] = si; cnpp += 2;
    for (E_Int j = 0; j < si; j++) cnpp[j] = comp[j];  
    
    K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10,
                                     fp, cnp);
    tpl = K_ARRAY::buildArray(fp, varString, cnp, -1, "NGON");
    delete &fp; delete components[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}
//=============================================================================
/* splitConnexity NODE : retourne une liste de zones avec 1 elt NODE */
//=============================================================================
PyObject* K_TRANSFORM::splitConnexityNODE(FldArrayF* f, FldArrayI* cn,
                                          char* eltType, char* varString,
                                          E_Int posx, E_Int posy, E_Int posz)
{
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType,
                               *f, *cn);
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();

  for (E_Int i = 0; i < npts; i++)
  {
    FldArrayF* f0 = new FldArrayF(1,nfld);
    E_Float* fp = f0->begin();
    FldArrayI* cnp = new FldArrayI(0);
    for (E_Int v = 1; v <= nfld; v++) fp[v-1] = (*f)(i,v);

    tpl = K_ARRAY::buildArray(*f0, varString, *cnp, -1, eltType);
    delete f0; delete cnp;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l; 
}
