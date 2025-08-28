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

# include "converter.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Modify fields for points of indices in listIndices 
   fields are given in array. */
//=============================================================================
PyObject* K_CONVERTER::updatePartialFields(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayF;// array 1D contenant les champs a inserer
  PyObject* listIndicesO;
  if (!PyArg_ParseTuple(args, "OOO", &array, &arrayF, &listIndicesO))
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
                    "updatePartialFields: 1st array is not valid.");
    return NULL; 
  }
  E_Int npts = f->getSize(); E_Int nfld = f->getNfld();
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
                    "updatePartialFields: 2nd array is not valid.");
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
                    "updatePartialFields: 3rd arg must be a numpy of integers.");
    return NULL;
  }

  PyObject* tpl; 
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varString, ni, nj, nk);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld();
    
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  fn = *f; // on copie pour ts les autres points
  E_Int nPts = listIndices->getSize();
  E_Int* indices = listIndices->begin();
  E_Int nfldc = posv.size();//nb de variables communes
#pragma omp parallel default(shared)
  {
  E_Int ind; E_Int pos1, posv1;
  for (E_Int eq = 0; eq < nfldc; eq++)
  {
    pos1 = posv[eq]; posv1 = posvl[eq]; 
    E_Float* foutp = fn.begin(pos1);
    E_Float* fip = fl->begin(posv1);
#pragma omp for
    for (E_Int i = 0; i < nPts; i++)
    {
      ind = indices[i];
      foutp[ind] += fip[i];
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
PyObject* K_CONVERTER::updatePartialFieldsPT(PyObject* self, PyObject* args)
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
                    "updatePartialFields: no field to set.");
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
                    "updatePartialFields: field array is not valid.");
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
                    "updatePartialFields: 3rd arg must be a numpy of integers.");
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
           foutp[ind] += fip[i];
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
           foutp[ind] += fip[i];
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
PyObject* K_CONVERTER::_updatePartialFields(PyObject* self, PyObject* args)
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
  if (startFrom == 0) { 
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
            foutp[ind] += fip[i];
          }
        }
      }
    }
  }
  else{
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
            foutp[ind] += fip[i];
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
