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
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
// Add a variable to an array: variable can be defined by a string or by
// an additional array 
//=============================================================================
PyObject* K_CONVERTER::addVar(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* additional;
  if (!PyArg_ParseTuple(args, "OO", &array, &additional)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, 
                              cn, eltType, true);
  if (res != 1 && res != 2) return NULL; // errors are alread set
  PyObject* tpl; 
  // Check additional
#if PY_VERSION_HEX >= 0x03000000
  if (PyString_Check(additional) || PyUnicode_Check(additional))
  {
    char* name;
    if (PyString_Check(additional)) name = PyString_AsString(additional);
    else name = PyBytes_AsString(PyUnicode_AsUTF8String(additional));
#else
  if (PyString_Check(additional))
  {
    // String name
    char* name = PyString_AsString(additional);
#endif
    
    // Name must be a unique var name
    E_Int i = 0;
    while (name[i] != '\0')
    {
      if (name[i] == ',')
      {
        RELEASESHAREDB(res, array, f, cn);
        PyErr_SetString(PyExc_TypeError,
                        "addVar: can only add a unique variable name.");
        return NULL;
      }
      i++;
    }

    E_Int varStringL = strlen(varString);
    E_Int nameL = strlen(name);
    char* fstring = new char [varStringL+nameL+2];

    strcpy(fstring, varString);
    E_Int pos = K_ARRAY::isNamePresent(name, varString);
    E_Int nt, sizet;
    nt = f->getNfld();

    if (pos == -1)
    { 
      sizet = nt+1; strcat(fstring, ","); strcat(fstring, name);
    }
    else sizet = nt;

    E_Int fSize = f->getSize();
    // Building array here
    if (res == 1) 
      tpl = K_ARRAY::buildArray(sizet, fstring, nil, njl, nkl);
    else
    {
      E_Int csize = cn->getSize()*cn->getNfld();
      tpl = K_ARRAY::buildArray(sizet, fstring, fSize, cn->getSize(), -1, 
                                eltType, false, csize);
      E_Int* cnRef = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnRef, cn->begin(), cn->getSize()*cn->getNfld());
    }

    E_Float* sp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF s(fSize, sizet, sp, true);

#pragma omp parallel default(shared)
    {
      for (E_Int i = 1; i <= nt; i++)
      {
        E_Float* spi = s.begin(i);
        E_Float* fp = f->begin(i);
#pragma omp for
        for (E_Int j = 0; j < fSize; j++) spi[j] = fp[j];
      }
      if (pos == -1) // on initialise que si c'est une nouvelle variable
      {
        E_Float* spi = s.begin(nt+1);
        if (K_STRING::cmp(name, "cellN") == 0 || 
            K_STRING::cmp(name, "cellNF") == 0)
        {
#pragma omp for
          for (E_Int j = 0; j < fSize; j++) spi[j] = 1.;
        }
        else
        { 
#pragma omp for
          for (E_Int j = 0; j < fSize; j++) spi[j] = 0.;
        }
      }
    }
    delete [] fstring;
  }
  else
  {
    // Additional must be an array
    E_Int res2;
    E_Int ni2, nj2, nk2;
    FldArrayF* f2; FldArrayI* cn2;
    char* varString2; char* eltType2;
    res2 = K_ARRAY::getFromArray(
      additional, varString2, f2, ni2, nj2, nk2, cn2, eltType2, true);

    if (res2 != 1 && res2 != 2) return NULL; // errors are alread set

    if (res2 == 1 && res != 1)
    {
      RELEASESHAREDB(res, array, f, cn);
      RELEASESHAREDB(res2, additional, f2, cn2);
      PyErr_SetString(PyExc_TypeError,
                      "addVar: additional must be located on the same grid as array.");
      return NULL;
    }
    if (res2 == 2 && res != 2)
    {
      RELEASESHAREDB(res, array, f, cn);
      RELEASESHAREDB(res2, additional, f2, cn2);
      PyErr_SetString(PyExc_TypeError,
                      "addVar: additional must be located on the same grid as array.");
      return NULL;
    }
    if (res2 == 1)
    {
      if (f->getSize() != f2->getSize())
      {
        RELEASESHAREDB(res, array, f, cn);
        RELEASESHAREDB(res2, additional, f2, cn2);
        PyErr_SetString(PyExc_TypeError,
                        "addVar: additional must be located on the same grid as array.");
        return NULL;
      }
    }
    if (res2 == 2)
    {
      if (cn->getSize() != cn2->getSize() || f->getSize() != f2->getSize())
      {
        RELEASESHAREDB(res, array, f, cn);
        RELEASESHAREDB(res2, additional, f2, cn2);
        PyErr_SetString(PyExc_TypeError,
                        "addVar: additional must be located on the same grid as array.");
        return NULL;
      }
    }
    
    // ExtractVars de varString2
    vector<char*> vars;
    vector<E_Int> pos1; vector<E_Int> pos2;
    K_ARRAY::extractVars(varString2, vars);
    E_Int sizet = f->getNfld();

    E_Int varStringL = strlen(varString);
    E_Int varString2L = strlen(varString2);
    char* fstring = new char [varStringL+varString2L+2];

    strcpy(fstring, varString);
    E_Int sizevars = vars.size();
    for (E_Int v = 0; v < sizevars; v++)
    {
      E_Int r = K_ARRAY::isNamePresent(vars[v], varString);
      if (r == -1)
      {
        strcat(fstring, ",");
        strcat(fstring, vars[v]);
        sizet++; pos1.push_back(sizet); pos2.push_back(v+1);
      }
      else
      {
        pos1.push_back(r+1); pos2.push_back(v+1);
      }
    }
    for (E_Int v = 0; v < sizevars; v++) delete vars[v];
    
    E_Int fSize = f->getSize();
    // Building array here
    if (res == 1) 
      tpl = K_ARRAY::buildArray(sizet, fstring, nil, njl, nkl);
    else
    {
      E_Int csize = cn->getSize()*cn->getNfld();
      tpl = K_ARRAY::buildArray(sizet, fstring, fSize, cn->getSize(), -1, eltType, false, csize);
      E_Int* cnRef = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnRef, cn->begin(), cn->getSize()*cn->getNfld());
    }
    E_Float* sp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF s(fSize, sizet, sp, true);

#pragma omp parallel default(shared)
    {
      for (E_Int i = 1; i <= sizet; i++)
      {
        E_Float* spi = s.begin(i);
        E_Float* fp = f->begin(i);
#pragma omp for
        for (E_Int j = 0; j < fSize; j++) spi[j] = fp[j];
      }
      E_Int sizepos1 = pos1.size();
      for (E_Int i = 0; i < sizepos1; i++)
      {
        E_Float* spi = s.begin(pos1[i]);
        E_Float* f2p = f2->begin(pos2[i]);
#pragma omp for
        for (E_Int j = 0; j < fSize; j++) spi[j] = f2p[j];
      }
    }
    delete [] fstring;
    RELEASESHAREDB(res2, additional, f2, cn2);
  }
 
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}

//=============================================================================
// Add vars to an array
//=============================================================================
PyObject* K_CONVERTER::addVars(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  if (!PyArg_ParseTuple(args, "O", &arrays)) return NULL;

  // Check array
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "addVars: argument must be a list.");
    return NULL;
  }
  E_Int n = PyList_Size(arrays);
  if (n == 0)
  {
    Py_INCREF(Py_None); return Py_None;
  }

  PyObject* tpl, *array;
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int sizevars;

  // dimensionnement varString
  E_Int varStringL = 0;
  for (int l = 0; l < n; l++)
  {
    array = PyList_GetItem(arrays, l);
    tpl = PyList_GetItem(array,0);
    if (PyString_Check(tpl)) varString = PyString_AsString(tpl);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl)) varString = PyBytes_AsString(PyUnicode_AsUTF8String(tpl));
#endif
    varStringL += strlen(varString)+4;
  }
  char* fstring = new char [varStringL];  // var string du array de sortie
  fstring[0] = '\0';

  vector<char*> vars; // vars du array de sortie
  vector<char*> local; // vars du array courant
  FldArrayI pos; // position des variables non communes dans le array de sortie
  E_Int* posp;
  char* localj;

  // Extraction du nombre total de variables pour chaque arrays 
  // et verification
  E_Int npts = -1;
  E_Int* cnRef = NULL;
  char eltTypeRef[256];
  E_Int sizelocal = 0;

  // Nombre de variables et taille de l'array final
  E_Int size=0, res, nvar=0, sizeRef=0;
  E_Int ni, nj, nk, nvertex, nelt, sizeConnect;
  E_Int niRef=0, njRef=0, nkRef=0, nvertexRef=0, neltRef=0, sizeConnectRef=0;
  E_Int structured=-1; // type de la sortie = type du premier array
  for (int l = 0; l < n; l++)
  {
    array = PyList_GetItem(arrays, l);
    res = K_ARRAY::getInfoFromArray(array, varString,
                                    ni, nj, nk, nvertex, nelt, 
                                    sizeConnect, eltType);
    if (res == 1) size = ni*nj*nk;
    else if (res == 2) size = nvertex;
    if (res == 1 || res == 2)
    {
      if (structured == -1)
      {
        structured = res;
        sizeRef = size;
        if (res == 1)
        {
          niRef = ni; njRef = nj; nkRef = nk;
        }
        else if (res == 2)
        {
          nvertexRef = nvertex; neltRef = nelt; sizeConnectRef = sizeConnect;
          strcpy(eltTypeRef, eltType);
        }
      }

      // Selectionne les variables non communes
      K_ARRAY::extractVars(varString, local);
      nvar = vars.size();
      sizelocal = local.size();
      for (E_Int j = 0; j < sizelocal; j++)
      {
        sizevars  = vars.size();
        localj = local[j];
        E_Boolean exist = false;
        for (E_Int i = 0; i < sizevars; i++)
        {
          if (K_STRING::cmp(vars[i], localj) == 0){exist = true; break;}
        }
        if (exist == false) // var non commune
        {
          vars.push_back(localj);
          strcat(fstring, localj);
          strcat(fstring, ",");
          nvar++;
        }
        else delete [] localj;
      }
      local.clear();
    }
  }
  E_Int nvarRef = nvar;
  // fstring final modification
  E_Int leng = strlen(fstring)-1;
  if (fstring[leng] == ',') fstring[leng] = '\0';

  // Construit le numpy de sortie
  if (structured == 1)
  {
    tpl = K_ARRAY::buildArray(nvarRef, fstring, 
                              niRef, njRef, nkRef);
  }
  else
  {  
    tpl = K_ARRAY::buildArray(nvarRef, fstring, 
                              nvertexRef, neltRef, -1, eltTypeRef,
                              false, sizeConnectRef);
    cnRef = K_ARRAY::getConnectPtr(tpl);
  }
  E_Float* fielp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF field(sizeRef, nvarRef, fielp, true);

  structured = -1; nvar = 0;
  for (int l = 0; l < n; l++) 
  { 
    array = PyList_GetItem(arrays, l);
    res = K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, 
                                cn, eltType, true);

    if (res != 1 && res != 2)
    {
      printf("Warning: addVars: array is invalid. Array %d skipped...\n",
             l+1);
      goto skip;
    }

    if (structured == 1 && res == 2)
    {
      RELEASESHAREDU(array, f, cn);
      printf("Warning: addVars: arrays must be defined on the same grid. Array %d skipped...\n", l+1);
      goto skip;
    }
    if (structured == 2 && res == 1)
    {
      RELEASESHAREDS(array, f);
      printf("Warning: addVars: arrays must be defined on the same grid. Array %d skipped...\n", l+1);
      goto skip;
    }
    if (structured != -1 && f->getSize() != sizeRef)
    {
      RELEASESHAREDB(res, array, f, cn);
      printf("Warning: addVars: arrays must be defined on the same grid. Array %d skipped...\n", l+1);
      goto skip;
    }

    if (structured == -1)
    {
      structured = res;
      if (res == 2)
      {
        E_Int* cnp = cn->begin();
        E_Int size = cn->getSize()*cn->getNfld();
        for (E_Int i = 0; i < size; i++) cnRef[i] = cnp[i];
      }
    }

    // Selectionne les variables non communes
    K_ARRAY::extractVars(varString, local);
    
    nvar = vars.size();
    sizelocal = local.size();
    pos.malloc(sizelocal); posp = pos.begin();
    for (E_Int j = 0; j < sizelocal; j++)
    {
      posp[j] = 1;
      sizevars = vars.size();
      localj = local[j];
      for (E_Int i = 0; i < sizevars; i++)
      {
        if (K_STRING::cmp(vars[i], localj) == 0) { posp[j] = i+1; break; }
      }
      delete [] localj;
    }
    local.clear();
    npts = f->getSize();

    // Init
#pragma omp parallel default(shared)
    {
      E_Int nfld = f->getNfld();
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fip = field.begin(posp[eq-1]);
        E_Float* fp = f->begin(eq);
#pragma omp for nowait
        for (E_Int i = 0; i < npts; i++) fip[i] = fp[i];
      }
    }

    RELEASESHAREDB(res, array, f, cn);
    skip: ;
    
  }

  sizevars = vars.size();
  for (E_Int i = 0; i < sizevars; i++) delete [] vars[i];
  delete [] fstring;

  return tpl;
}
