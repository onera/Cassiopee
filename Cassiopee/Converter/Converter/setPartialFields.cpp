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

# include "converter.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Modify fields for points of indices in listIndices 
   fields are given in array. */
//=============================================================================
PyObject* K_CONVERTER::setPartialFields(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayF;// array 1D contenant les champs a inserer
  PyObject* listIndicesO;
  if (!PYPARSETUPLE_(args, OOO_, &array, &arrayF, &listIndicesO))
    return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: 1st array is not valid.");
    return NULL; 
  }
  // Check arrayF
  E_Int nil, njl, nkl;
  FldArrayF* fl; FldArrayI* cnl;
  char* varStringl; char* eltTypel;
  E_Int resl = K_ARRAY::getFromArray3(arrayF, varStringl, fl, nil, njl, nkl, 
                                      cnl, eltTypel);
  if (resl != 1 && resl != 2) 
  {
    RELEASESHAREDB(res,array,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: 2nd array is not valid.");
    return NULL; 
  }
  vector<E_Int> posv; vector<E_Int> posvl;
  char* varStringC; // chaine de caractere commune
  E_Int l = strlen(varString);
  varStringC = new char [l+1];
  K_ARRAY::getPosition(varString, varStringl, posv, posvl, varStringC);
  delete [] varStringC;
  
  /*--------------------------------------------*/
  /* Extraction des indices des pts a modifier */
  /*-------------------------------------------*/
  FldArrayI* listIndices;
  E_Int resi = K_NUMPY::getFromNumpyArray(listIndicesO, listIndices);
  if (resi == 0)
  {
    RELEASESHAREDB(res, array, f, cn); 
    RELEASESHAREDB(resl, arrayF, fl, cnl); 
    PyErr_SetString(PyExc_TypeError, 
                    "setPartialFields: 3rd arg must be a numpy of integers.");
    return NULL;
  }

  PyObject* tpl;
  E_Int api = f->getApi();
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray3(*f, varString, ni, nj, nk, api);
  } 
  else //unstructured 
  {    
    tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType, api);
  }
  FldArrayF* fn;
  K_ARRAY::getFromArray3(tpl, fn);

  E_Int nPts = listIndices->getSize();
  E_Int* indices = listIndices->begin();
  E_Int nfldc = posv.size();//nb de variables communes

#pragma omp parallel default(shared)
  {
    E_Int ind; E_Int pos1, posv1;
    for (E_Int eq = 0; eq < nfldc; eq++)
    {
      pos1 = posv[eq]; posv1 = posvl[eq]; 
      E_Float* foutp = fn->begin(pos1);
      E_Float* fip = fl->begin(posv1);
      #pragma omp for
      for (E_Int i = 0; i < nPts; i++)
      {
          ind = indices[i];
          foutp[ind] = fip[i];
      }
    }
  }
  RELEASESHAREDN(listIndicesO, listIndices);
  RELEASESHAREDB(res, array, f, cn); 
  RELEASESHAREDB(resl, arrayF, fl, cnl); 
  return tpl;
}
//=============================================================================
/* Modify fields for points of indices in listIndices
   IN: z: zone a modifier
   IN: arrayF: fields in array. */
//=============================================================================
PyObject* K_CONVERTER::setPartialFieldsPT(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject* zone;
  PyObject* arrayF; // array 1D contenant les champs a inserer
  PyObject* listIndicesO;
  E_Int loc;
  E_Int startFrom;
  if (!PYPARSETUPLE_(args, OOO_ I_ SSS_ I_, &zone, &arrayF, &listIndicesO, &loc, 
		     &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters, &startFrom))
    { return NULL; }

  /* zone a modifier */
  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  K_PYTREE::getFromZone(zone, 0, loc, varString,
                        fields, locs, im, jm, km,
                        cn, cnSize, cnNfld, eltType, hook,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

  E_Int nfld = fields.size();
  if (nfld == 0) 
  {
    RELEASESHAREDZ(hook, varString, eltType);
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: no field to set.");
    return NULL;
  }

  // Check arrayF: champs a inserer
  E_Int nil, njl, nkl;
  FldArrayF* fl; FldArrayI* cnl;
  char* varStringl; char* eltTypel;
  E_Int resl = K_ARRAY::getFromArray3(arrayF, varStringl, fl, nil, njl, nkl, 
                                      cnl, eltTypel);
  if (resl != 1 && resl != 2) 
  {
    RELEASESHAREDZ(hook, varString, eltType);
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: field array is not valid.");
    return NULL; 
  }
  vector<E_Int> posv; vector<E_Int> posvl;
  char* varStringC; // chaine de caractere commune
  E_Int l = strlen(varString);
  varStringC = new char [l+1];
  K_ARRAY::getPosition(varString, varStringl, posv, posvl, varStringC);
  delete [] varStringC;

  /*--------------------------------------------*/
  /* Extraction des indices des pts a modifier */
  /*-------------------------------------------*/
  FldArrayI* listIndices;
  E_Int resi = K_NUMPY::getFromNumpyArray(listIndicesO, listIndices);
  if (resi == 0)
  {
    RELEASESHAREDZ(hook, varString, eltType);
    RELEASESHAREDB(resl, arrayF, fl, cnl); 
    PyErr_SetString(PyExc_TypeError, 
                    "setPartialFields: 3rd arg must be a numpy of integers.");
    return NULL;
  }

  E_Int nPts = listIndices->getSize();
  E_Int* indices = listIndices->begin();
  E_Int nfldc = posv.size(); // nb de variables communes

  if (startFrom == 0) 
  { 
#pragma omp parallel default(shared)
    {
      E_Int ind; E_Int pos1, posv1;
      for (E_Int eq = 0; eq < nfldc; eq++)
      {
        pos1 = posv[eq]; posv1 = posvl[eq];
        E_Float* foutp = fields[pos1-1];
        E_Float* fip = fl->begin(posv1);
        #pragma omp for
        for (E_Int i = 0; i < nPts; i++)
        {
          ind = indices[i]; 
          foutp[ind] = fip[i];
        }
      }
    }
  }
  else
  {
#pragma omp parallel default(shared)
    {
      E_Int ind; E_Int pos1, posv1;
      for (E_Int eq = 0; eq < nfldc; eq++)
      {
        pos1 = posv[eq]; posv1 = posvl[eq];
        E_Float* foutp = fields[pos1-1];
        E_Float* fip = fl->begin(posv1);
        #pragma omp for
        for (E_Int i = 0; i < nPts; i++)
        {
          ind = indices[i]-startFrom; 
          foutp[ind] = fip[i];
        }
      }
    }
  }

  RELEASESHAREDN(listIndicesO, listIndices);
  RELEASESHAREDB(resl, arrayF, fl, cnl);
  RELEASESHAREDZ(hook, varString, eltType);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* SetPartialFields in a pyTree zone in place, starting from numpy fields
   (one for each field) and numpy indices  */
//=============================================================================
PyObject* K_CONVERTER::_setPartialFields(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject* zone;
  PyObject* listNumFields;// numpy contenant les champs a inserer
  PyObject* listIndicesO;
  E_Int loc;
  E_Int startFrom;
  if (!PYPARSETUPLE_(args, OOO_ I_ SSS_ I_, &zone, &listNumFields, &listIndicesO, &loc,
                     &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters, &startFrom)) return NULL;

  /*--------------------------------------------*/
  /* zone a modifier                            */
  /*--------------------------------------------*/
  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  K_PYTREE::getFromZone(zone, 0, loc, varString,
                        fields, locs, im, jm, km,
                        cn, cnSize, cnNfld, eltType, hook,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

  /*-------------------------------------------*/
  /* Extraction des champs sous forme de numpy */
  /*-------------------------------------------*/
  E_Int nfld = PyList_Size(listNumFields);
  vector<FldArrayF*> listFields;
  PyObject* tpl;
  vector<PyObject*> objs;
  vector<E_Int> isEmpty(nfld);
  for (E_Int v = 0; v < nfld; v++)
  {
    FldArrayF* oneField;
    tpl = PyList_GetItem(listNumFields, v);
    E_Int resf = K_NUMPY::getFromNumpyArray(tpl, oneField);
    if (resf == 0) { isEmpty[v] = 1; oneField = NULL; } 
    else if (oneField->getSize() == 0) isEmpty[v] = 1;
    else isEmpty[v] = 0;
    listFields.push_back(oneField);
    objs.push_back(tpl);
  }
  /*-------------------------------------------*/
  /* Extraction des indices des pts a modifier */
  /*-------------------------------------------*/
  FldArrayI* listIndices;
  K_NUMPY::getFromNumpyArray(listIndicesO, listIndices);
  E_Int nPts = listIndices->getSize();
  E_Int* indices = listIndices->begin();
  
  // no check: perfos
  if (startFrom == 0) 
  { 
#pragma omp parallel default(shared)
    {
      E_Int ind;
      for (E_Int v = 0; v < nfld; v++)
      {
        if (isEmpty[v] == 0)
        {
          E_Float* foutp = fields[v];
          E_Float* fip = listFields[v]->begin();
          #pragma omp for
          for (E_Int i = 0; i < nPts; i++)
          {
            ind = indices[i];
            foutp[ind] = fip[i];
          }
        }
      }
    }
  }
  else
  {
#pragma omp parallel default(shared)
    {
      E_Int ind;
      for (E_Int v = 0; v < nfld; v++)
      {
        if (isEmpty[v] == 0)
        {
          E_Float* foutp = fields[v];
          E_Float* fip = listFields[v]->begin();
          #pragma omp for
          for (E_Int i = 0; i < nPts; i++)
          {
            ind = indices[i]-startFrom;
            foutp[ind] = fip[i];
          }
        }
      }
    }
  }

  // sortie
  RELEASESHAREDZ(hook, varString, eltType);
  RELEASESHAREDN(listIndicesO, listIndices);
  for (E_Int v = 0; v < nfld; v++)    
  {
    if (isEmpty[v] == 0) RELEASESHAREDN(objs[v], listFields[v]);
  } 
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* SetPartialFields with accumulation  */
//=============================================================================
PyObject* K_CONVERTER::_setPartialFieldsAverage(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* indexList; PyObject* valueList;
  if (!PYPARSETUPLE_(args, OOO_,
                     &array, &indexList, &valueList)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: invalid array.");
    return NULL;
  }

  // Check index list
  if (PyList_Check(indexList) == false)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: invalid array.");
    return NULL;
  }
  E_Int size = PyList_Size(indexList);
  vector<FldArrayI*> inds(size);
  for (E_Int i = 0; i < size; i++)
  {
    FldArrayI* t;
    PyObject* index = PyList_GetItem(indexList, i);
    E_Int res2 = K_NUMPY::getFromNumpyArray(index, t);

    if (res2 == 0)
    {
      RELEASESHAREDB(res, array, f, cn);
      PyErr_SetString(PyExc_TypeError, 
                      "setPartialFields: index numpy is invalid.");
      return NULL;
    }
    inds[i] = t;
  }

  // Check value list
  if (PyList_Check(valueList) == false)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "setPartialFields: invalid value array.");
    return NULL;
  }

  vector<FldArrayF*> values(size);
  for (E_Int i = 0; i < size; i++)
  {
    PyObject* v = PyList_GetItem(valueList, i);

    E_Int ni2, nj2, nk2; 
    FldArrayF* f2; FldArrayI* cn2;
    char* varString2; char* eltType2;
    E_Int res2 = K_ARRAY::getFromArray3(v, varString2, f2, ni2, nj2, nk2, 
                                        cn2, eltType2);
    if (res2 != 1 && res2 != 2)
    {
      RELEASESHAREDB(res, array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "setPartialFields: invalid value array.");
      return NULL;
    }
    
    values[i] = f2;
  }  

  // update with accu (first pass)
  E_Int nall = f->getSize();
  E_Int nfld = f->getNfld();
  E_Int* accu = new E_Int [nall];
  for (E_Int i = 0; i < nall; i++) accu[i] = 0;
  
  for (E_Int i = 0; i < size; i++)
  {
    FldArrayI* index = inds[i];
    FldArrayF* v = values[i];
    E_Int npts = index->getSize();
    E_Int* indexp = index->begin();
    
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* vp = v->begin(n);
  
      for (E_Int j = 0; j < npts; j++)
      {
        E_Int ind = indexp[i];
        if (n == 1) accu[ind]++;
        fp[ind] += vp[j];
      }
    }
  }

  // redivise par accu
  for (E_Int n = 1; n <= nfld; n++)
  {
    E_Float* fp = f->begin(n);
    for (E_Int i = 0; i < nall; i++)
    {
      if (accu[i] == 1) fp[i] = fp[i]*0.5;
      else if (accu[i] == 2) fp[i] = fp[i]/3.;
      else fp[i] = fp[i]*0.25;
    }
  }

  // Release
  delete [] accu;
  RELEASESHAREDB(res, array, f, cn);
  for (E_Int i = 0; i < size; i++)
  {
    PyObject* index = PyList_GetItem(indexList, i);
    RELEASESHAREDN(index, inds[i]);
  }
  Py_INCREF(Py_None);
  return Py_None;
}
