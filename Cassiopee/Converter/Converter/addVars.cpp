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
  if (!PYPARSETUPLE_(args, OO_, &array, &additional)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
                               cn, eltType);
  if (res != 1 && res != 2) return NULL; // errors are alread set
  PyObject* tpl; 
  // Check additional
#if PY_VERSION_HEX >= 0x03000000
  if (PyString_Check(additional) || PyUnicode_Check(additional))
  {
    char* name;
    if (PyString_Check(additional)) name = PyString_AsString(additional);
    else name = (char*)PyUnicode_AsUTF8(additional);
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
    E_Int api = f->getApi();
    if (res == 1) 
      tpl = K_ARRAY::buildArray3(sizet, fstring, nil, njl, nkl, api);
    else
    {
      tpl = K_ARRAY::buildArray3(sizet, fstring, fSize,
                                 *cn, eltType, false, api, true);
    }
    FldArrayF* s;
    K_ARRAY::getFromArray3(tpl, s);

    #pragma omp parallel default(shared)
    {
      for (E_Int i = 1; i <= nt; i++)
      {
        E_Float* spi = s->begin(i);
        E_Float* fp = f->begin(i);
        #pragma omp for
        for (E_Int j = 0; j < fSize; j++) spi[j] = fp[j];
      }
      if (pos == -1) // on initialise que si c'est une nouvelle variable
      {
        E_Float* spi = s->begin(nt+1);
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
    RELEASESHAREDS(tpl, s);
  }
  else
  {
    // Additional must be an array
    E_Int res2;
    E_Int ni2, nj2, nk2;
    FldArrayF* f2; FldArrayI* cn2;
    char* varString2; char* eltType2;
    res2 = K_ARRAY::getFromArray3(additional, varString2, f2, ni2, nj2, nk2, cn2, eltType2);

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
    E_Int api = f->getApi();
    // Building array here
    if (res == 1) 
      tpl = K_ARRAY::buildArray3(sizet, fstring, nil, njl, nkl, api);
    else
    {
      tpl = K_ARRAY::buildArray3(sizet, fstring, fSize,
                                 *cn, eltType, false, api, true);
    }
    FldArrayF* s;
    K_ARRAY::getFromArray3(tpl, s);

    #pragma omp parallel default(shared)
    {
      for (E_Int i = 1; i <= sizet; i++)
      {
        E_Float* spi = s->begin(i);
        E_Float* fp = f->begin(i);
        #pragma omp for
        for (E_Int j = 0; j < fSize; j++) spi[j] = fp[j];
      }
      E_Int sizepos1 = pos1.size();
      for (E_Int i = 0; i < sizepos1; i++)
      {
        E_Float* spi = s->begin(pos1[i]);
        E_Float* f2p = f2->begin(pos2[i]);
        #pragma omp for
        for (E_Int j = 0; j < fSize; j++) spi[j] = f2p[j];
      }
    }
    delete [] fstring;
    RELEASESHAREDS(tpl, s);
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
  if (!PYPARSETUPLE_(args, O_, &arrays)) return NULL;

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

  PyObject *tpl, *array;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int sizevars;

  // dimensionnement varString
  E_Int varStringL = 1;
  for (int l = 0; l < n; l++)
  {
    array = PyList_GetItem(arrays, l);
    tpl = PyList_GetItem(array,0);
    if (PyString_Check(tpl)) varString = PyString_AsString(tpl);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl)) varString = (char*)PyUnicode_AsUTF8(tpl);
#endif    
    varStringL += strlen(varString) + 4;
  }
  char* varString2 = new char [varStringL];  // var string du array de sortie
  varString2[0] = '\0';

  vector<char*> vars; // vars du array de sortie
  vector<char*> varStrings; // vars du array courant
  FldArrayI pos; // position des variables non communes dans le array de sortie
  E_Int* posp;
  char* localj;

  // Extraction du nombre total de variables pour chaque arrays 
  // et verification
  E_Int res, res0 = -1;
  E_Int npts, nfld, api, nelts, ni, nj, nk;
  E_Int sizelocal = 0, nvar = 0;
  E_Int npts0 = -1, ni0 = -1, nj0 = -1, nk0 = -1, nelts0 = -1;
  E_Int iarr0 = -1;  // index of the first valid array
  E_Bool center0 = false;
  char* eltType0 = new char[K_ARRAY::VARSTRINGLENGTH];
  eltType0[0] = '\0';
  E_Bool arrValid[n];

  for (E_Int l = 0; l < n; l++)
  {
    array = PyList_GetItem(arrays, l);
    res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

    if (res == 1)
    {
      if (res0 == -1)
      {
        arrValid[l] = true;
        iarr0 = l; res0 = res;
        ni0 = ni; nj0 = nj; nk0 = nk; npts0 = ni0*nj0*nk0;
      }
      else if (ni0 != ni || nj0 != nj || nk0 != nk)
      {
        printf("Warning: addVars: arrays must be defined on the same grid. "
             "Array " SF_D_ " skipped...\n", l+1);
        arrValid[l] = false;
        RELEASESHAREDS(array, f);
        continue;
      }
      else arrValid[l] = true;  // valid array
    }
    else if (res == 2)
    {
      npts = f->getSize();
      if (K_STRING::cmp(eltType, 4, "NGON") == 0) nelts = cn->getNElts();
      else
      {
        // Compute total number of elements
        E_Int nc = cn->getNConnect();
        nelts = 0;
        for (E_Int ic = 0; ic < nc; ic++)
        {
          K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
          nelts += cm.getSize();
        }
      }
      E_Bool center = false; //npts == nelts; TODO

      if (res0 == -1)
      {
        arrValid[l] = true;
        iarr0 = l; res0 = res;
        npts0 = npts; nelts0 = nelts;        
        strcpy(eltType0, eltType);
        if (center || strchr(eltType0, '*') != NULL) center0 = true;
      }
      else if ((center0 && nelts != nelts0) || (!center0 && npts != npts0))
      {
        printf("Warning: addVars: arrays must be defined on the same grid. "
             "Array " SF_D_ " skipped...\n", l+1);
        arrValid[l] = false;
        RELEASESHAREDU(array, f, cn);
        continue;
      }
      else arrValid[l] = true;  // valid array
    }
    else
    {
      printf("Warning: addVars: array is invalid. Array " SF_D_ " skipped...\n", l+1);
      arrValid[l] = false;
      continue;
    }

    // Selectionne les variables non communes
    K_ARRAY::extractVars(varString, varStrings);
    nvar = vars.size();
    sizelocal = varStrings.size();
    for (E_Int j = 0; j < sizelocal; j++)
    {
      sizevars = vars.size();
      localj = varStrings[j];
      E_Bool exist = false;
      for (E_Int i = 0; i < sizevars; i++)
      {
        if (K_STRING::cmp(vars[i], localj) == 0){exist = true; break;}
      }
      if (!exist) // var non commune
      {
        vars.push_back(localj);
        strcat(varString2, localj);
        strcat(varString2, ",");
        nvar++;
      }
      else delete [] localj;
    }
    varStrings.clear();
    RELEASESHAREDB(res, array, f, cn);
  }

  E_Int nvalidArrs = 0;
  for (E_Int l = 0; l < n; l++) nvalidArrs += (E_Int)arrValid[l];
  if (nvalidArrs == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addVars: none of the arrays are valid.");
    for (size_t i = 0; i < vars.size(); i++) delete [] vars[i];
    delete [] varString2;
    delete [] eltType0;
    return NULL;
  }

  // varString2 final modification
  E_Int nfld2 = nvar;
  E_Int leng = strlen(varString2) - 1;
  if (varString2[leng] == ',') varString2[leng] = '\0';

  // Get first valid array
  array = PyList_GetItem(arrays, iarr0);
  res0 = K_ARRAY::getFromArray3(array, varString, f, ni0, nj0, nk0, 
                                cn, eltType);
  api = f->getApi();
  
  // Construit le numpy de sortie
  FldArrayF* f2;
  if (res0 == 1)
  {
    tpl = K_ARRAY::buildArray3(nfld2, varString2, ni0, nj0, nk0, api);
  }
  else
  {  
    tpl = K_ARRAY::buildArray3(nfld2, varString2, f->getSize(),
                               *cn, eltType0, center0, api, true);
  }
  K_ARRAY::getFromArray3(tpl, f2);

  // Free memory from the first valid array
  RELEASESHAREDB(res0, array, f, cn);

  // Copy fields
  nvar = 0;
  for (E_Int l = 0; l < n; l++) 
  { 
    if (!arrValid[l]) continue;  // array is invalid, skip
    array = PyList_GetItem(arrays, l);
    K_ARRAY::getFromArray3(array, varString, f);
    nfld = f->getNfld();

    // Selectionne les variables non communes
    K_ARRAY::extractVars(varString, varStrings);
    nvar = vars.size();
    sizelocal = varStrings.size();
    pos.malloc(sizelocal); posp = pos.begin();
    for (E_Int j = 0; j < sizelocal; j++)
    {
      posp[j] = 1;
      sizevars = vars.size();
      localj = varStrings[j];
      for (E_Int i = 0; i < sizevars; i++)
      {
        if (K_STRING::cmp(vars[i], localj) == 0) { posp[j] = i+1; break; }
      }
      delete [] localj;
    }
    varStrings.clear();

    #pragma omp parallel default(shared)
    {
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* f2p = f2->begin(posp[eq-1]);
        E_Float* fp = f->begin(eq);
        #pragma omp for nowait
        for (E_Int i = 0; i < npts0; i++) f2p[i] = fp[i];
      }
    }

    RELEASESHAREDS(array, f);
  }

  for (size_t i = 0; i < vars.size(); i++) delete [] vars[i];
  delete [] varString2;
  delete [] eltType0;
  RELEASESHAREDS(tpl, f2);

  return tpl;
}
