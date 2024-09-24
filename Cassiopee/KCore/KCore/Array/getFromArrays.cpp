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
# include "Array/Array.h"
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Extrait les donnees utiles d'une liste d'objets pythons.
   IN: o: liste d'objets pythons correspondant a des arrays
   OUT: res: type de l'array dans arrays (1:struct, 2:unstruct, 0:invalid, -1: skipped)
   OUT: structVarString: chaine des variables pour les arrays structures
   OUT: unstructVarString: chaine des variables pour les arrays 
   non-structures
   OUT: structF: vecteurs des champs structures
   OUT: unstructF: vecteurs des champs non-structures
   OUT: ni,nj,nk: vecteur des dimensions pour les champ structures
   OUT: c, eltType: vecteur des connectivites et des types d'elements
   pour les champs non-structures
   OUT: objs: liste des objets python correspondant aux arrays structures
   OUT: obju: liste des objets python correspondant aux arrays non-structures
   IN: skipNoCoord: rejette les arrays sans x,y,z
   IN: skipStructured: rejette les arrays structures
   IN: skipUnstructured: rejette les arrays non structures
   IN: skipDiffVars: selectionne les arrays ayant les memes variables 
                     que celles du premier array
                     les variables sont positionnees de la meme maniere
                     pour tous les arrays
   Retourne -1 si pb
   C'est la responsabilite de l'appelant de liberer la memoire de structF,
   unstructF et c.
*/
//=============================================================================
E_Int K_ARRAY::getFromArrays(PyObject* o,
                             vector<E_Int>& res,
                             vector<char*>& structVarString,
                             vector<char*>& unstructVarString,
                             vector<FldArrayF*>& structF,
                             vector<FldArrayF*>& unstructF,
                             vector<E_Int>& ni, 
                             vector<E_Int>& nj, 
                             vector<E_Int>& nk,
                             vector<FldArrayI*>& c,
                             vector<char*>& eltType,
                             vector<PyObject*>& objs,
                             vector<PyObject*>& obju,
                             E_Boolean skipDiffVars,
                             E_Boolean skipNoCoord,
                             E_Boolean skipStructured,
                             E_Boolean skipUnstructured,
                             E_Boolean shared)
{
  char* varString; char* eltT;
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  PyObject* tpl;
  E_Int posx, posy, posz;
  E_Int resl;

  // o doit etre une liste
  if (PyList_Check(o) == false)
  {
    PyErr_SetString( PyExc_TypeError,
                     "getFromArrays: arrays argument must be a list.");
    return -1;
  }
  
  E_Int n = PyList_Size(o);
  
  // Dim varString common
  E_Int size = 1;
  E_Int nvertex, nelt, sizeConnect;
  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    getInfoFromArray(tpl, varString, nil, njl, nkl, nvertex, nelt, 
                     sizeConnect, eltT);
    size += strlen(varString)+4;
  }
  char* varStringCommon = new char [size]; 
  varStringCommon[0] = '\0';

  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    resl = K_ARRAY::getFromArray3(tpl, varString, f,
                                  nil, njl, nkl, cn, eltT);

    if (skipNoCoord)
    {
      posx = K_ARRAY::isCoordinateXPresent(varString);
      posy = K_ARRAY::isCoordinateYPresent(varString);
      posz = K_ARRAY::isCoordinateZPresent(varString);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        printf("Warning: getFromArrays: one array has no coordinates. Skipped...\n");
        res.push_back(-1);
        if (shared == false) { delete f; if (resl == 2) delete cn; }
        else { RELEASESHAREDB(resl, tpl, f, cn); }
        goto next;
      }
    }

    if (resl == 1 && skipStructured == false)
    {
      if (skipDiffVars)
      {
        if (varStringCommon[0] == '\0') strcpy(varStringCommon, varString);
        else  
        {
          if (compareVarStrings(varString, varStringCommon) != 0)
          {
            printf("Warning: getFromArrays: one array has different variables. Skipped...\n");
            res.push_back(-1);
            if (shared == false) { delete f; if (resl == 2) delete cn; }
            else { RELEASESHAREDB(resl, tpl, f, cn); }
            goto next;
          }
        }
      }
      res.push_back(1);
      structF.push_back(f);
      structVarString.push_back(varString);
      ni.push_back(nil); nj.push_back(njl); nk.push_back(nkl);
      objs.push_back(tpl);
    }
    else if (resl == 2 && skipUnstructured == false)
    {
      if (skipDiffVars)
      {
        if (varStringCommon[0] == '\0') strcpy(varStringCommon, varString);
        else
        {
          if (compareVarStrings(varString, varStringCommon) != 0)
          {
            printf("Warning: getFromArrays: one array has different variables. Skipped...\n");
            res.push_back(-1);
            if (shared == false) { delete f; if (resl == 2) delete cn; }
            else { RELEASESHAREDB(resl, tpl, f, cn); }
            goto next;
          }
        }
      }
      res.push_back(2);
      unstructF.push_back(f);
      unstructVarString.push_back(varString);
      c.push_back(cn);
      eltType.push_back(eltT);
      obju.push_back(tpl);
    }
    else
    {
      printf("Warning: getFromArrays: one array is invalid. Skipped...\n");
      res.push_back(0);
      if (shared == false) { delete f; if (resl == 2) delete cn; }
      else { RELEASESHAREDB(resl, tpl, f, cn); }
      goto next;
    }
    next:;
  }
  delete [] varStringCommon;
  return 1;
}

//=============================================================================
/* Extrait les donnees utiles d'une liste d'objets pythons.
   IN: o: liste d'objets pythons correspondant a des arrays
   OUT: res: type de l'array dans arrays (1:struct, 2:unstruct, 0:invalid, -1: skipped). Si l'array est invalide ou skipped, il n'y pas de F, a2,...
   correspondant.
   OUT: varString: chaine des variables pour chaque array
   OUT: F: vecteurs des champs
   OUT: a2,a3,a4,: vecteur d'information pour chaque array
   Pour un array structure: a2=&ni (int*), a3=&nj (int*), a4=&nk (int*)
   Pour un array non-structure: a2=cn (FldArrayI*), a3=eltType (char*), a4=NULL
   OUT: obj: liste des objets python correspondant aux arrays
   IN: skipNoCoord: rejette les arrays sans x,y,z
   IN: skipStructured: rejette les arrays structures
   IN: skipUnstructured: rejette les arrays non structures
   IN: skipDiffVars: selectionne les arrays ayant les memes variables 
                     que celles du premier array
                     les variables sont positionnees de la meme maniere
                     pour tous les arrays
   Retourne -1 si pb
   C'est la responsabilite de l'appelant de liberer la memoire de F,
   a2,a3,a4.
*/
//=============================================================================
E_Int K_ARRAY::getFromArrays(PyObject* o,
                             vector<E_Int>& res,
                             vector<char*>& varString,
                             vector<FldArrayF*>& F,
                             vector<void*>& a2,
                             vector<void*>& a3, 
                             vector<void*>& a4,
                             vector<PyObject*>& obj,
                             E_Boolean skipDiffVars,
                             E_Boolean skipNoCoord,
                             E_Boolean skipStructured,
                             E_Boolean skipUnstructured,
                             E_Boolean shared)
{
  char* varStringl; char* eltT;
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  PyObject* tpl;
  E_Int posx, posy, posz;
  E_Int resl;

  // o doit etre une liste
  if (!PyList_Check(o))
  {
    PyErr_SetString(PyExc_TypeError, 
                    "getFromArrays: arrays argument must be a list.");
    return -1;
  }
  
  E_Int n = PyList_Size(o);
  
  // Dim varString common
  E_Int size = 1;
  E_Int nvertex, nelt, sizeConnect;
  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    getInfoFromArray(tpl, varStringl, nil, njl, nkl, nvertex, nelt, 
                     sizeConnect, eltT);
    size += strlen(varStringl)+4;
  }
  char* varStringCommon = new char [size];
  varStringCommon[0] = '\0';

  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    resl = K_ARRAY::getFromArray3(tpl, varStringl, f, 
                                  nil, njl, nkl, cn, eltT);

    if (skipNoCoord)
    {
      posx = K_ARRAY::isCoordinateXPresent(varStringl);
      posy = K_ARRAY::isCoordinateYPresent(varStringl);
      posz = K_ARRAY::isCoordinateZPresent(varStringl);
      if (posx == -1 || posy == -1 || posz == -1)
      { 
        printf("Warning: getFromArrays: one array has no coordinates. Skipped...\n");
        res.push_back(-1);
        if (!shared) { delete f; if (resl == 2) delete cn; }
        else { RELEASESHAREDB(resl, tpl, f, cn); }
        goto next;
      }
    }

    if (resl == 1 && !skipStructured)
    {
      if (skipDiffVars)
      {
        if (varStringCommon[0] == '\0') strcpy(varStringCommon, varStringl);
        else  
        {
          if (compareVarStrings(varStringl, varStringCommon) != 0)
          {
            printf("Warning: getFromArrays: one array has different variables. Skipped...\n");
            res.push_back(-1);
            if (!shared) { delete f; if (resl == 2) delete cn; }
            else { RELEASESHAREDB(resl, tpl, f, cn); }
            goto next;
          }
        }
      }
      res.push_back(1);
      F.push_back(f);
      varString.push_back(varStringl);
      E_Int* p;
      p = new E_Int; *p = nil; a2.push_back(p);
      p = new E_Int; *p = njl; a3.push_back(p);
      p = new E_Int; *p = nkl; a4.push_back(p); 
      obj.push_back(tpl);
    }
    else if (resl == 2 && !skipUnstructured)
    {
      if (skipDiffVars)
      {
        if (varStringCommon[0] == '\0')
          strcpy(varStringCommon, varStringl);
        else
        {
          if (compareVarStrings(varStringl, varStringCommon) != 0)
          {
            printf("Warning: getFromArrays: one array has different variables. Skipped...\n");
            res.push_back(-1);
            if (!shared) { delete f; if (resl == 2) delete cn; }
            else { RELEASESHAREDB(resl, tpl, f, cn); }
            goto next;
          }
        }
      }
      res.push_back(2);
      F.push_back(f);
      varString.push_back(varStringl);
      a2.push_back(cn);
      a3.push_back(eltT);
      a4.push_back(NULL);
      obj.push_back(tpl);
    }
    else
    {
      printf("Warning: getFromArrays: one array is invalid. Skipped...\n");
      res.push_back(0);
      if (!shared) { delete f; if (resl == 2) delete cn; }
      else { RELEASESHAREDB(resl, tpl, f, cn); }
      goto next;
    }
    next:;
  }
  delete [] varStringCommon;
  return 1;
}
