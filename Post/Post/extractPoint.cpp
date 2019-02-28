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
# include "post.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Extrait la solution en un point de coordonnees (x,y,z). Si ce point est 
   dans une zone ou des maillages se recouvrent, la solution est extraite à
   partir du maillage ayant la plus petite cellule d interpolation */
// ============================================================================
PyObject* K_POST::extractPoint(PyObject* self, PyObject* args)
{
  PyObject *arrays, *listPts; 
  E_Int interpOrder;
  E_Int extrapOrder;
  E_Float constraint;
  PyObject* hook;
  if (!PYPARSETUPLE(args,
                    "OOlldO", "OOiidO",
                    "OOllfO", "OOiifO",
                    &arrays, &listPts, &interpOrder, &extrapOrder, &constraint, &hook))
  {
      return NULL;
  }
  // Check every array in arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extractPoint: first argument must be a list.");
    return NULL;
  }
  // Check every array in arrays
  if (PyList_Check(listPts) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extractPoint: 2nd argument must be a list.");
    return NULL;
  }
  // Recuperation de la liste des points
  E_Int npts = PyList_Size(listPts);// nb d elements ds la liste
  for (E_Int i = 0; i < npts; i++)
  {
    PyObject* tpli = PyList_GetItem(listPts,i);
    if (PyTuple_Check(tpli) == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "extractPoint: each element of the list must be (x,y,z).");
      
      return NULL;
    }
    E_Int dim = PyTuple_Size(tpli);
    // verifie que chq element de la liste est un triplet (x,y,z)
    if (dim != 3)
    {        
      PyErr_SetString(PyExc_TypeError,
                      "extractPoint: 3 coordinates must be found in each element of the list.");
      return NULL;  
    }
  }

  // Interpolation type 
  E_Int ncf;
  K_INTERP::InterpData::InterpolationType interpType;
  // attention le cas purement non structure est traite apres les interpDatas 
  switch (interpOrder)
  {
    case 2:
      interpType = K_INTERP::InterpData::O2CF;
      ncf = 8;
      break;
    case 3: 
      interpType = K_INTERP::InterpData::O3ABC;
      ncf = 9;
      break;
    case 5:
      interpType = K_INTERP::InterpData::O5ABC;
      ncf = 15;
      break;
    default:
      printf("Warning: extractPoint: unknown interpolation order. Set to 2nd order.\n");
      interpType = K_INTERP::InterpData::O2CF;
      ncf = 8;
      break;
  }

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, varString, fields, a2, a3, a4, objs,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = objs.size();
  E_Int nfldTot = fields[0]->getNfld();
  if (isOk == -1 || nfldTot < 4)
  {
    if (isOk == -1)
      PyErr_SetString(PyExc_TypeError,
                      "extractPoint: invalid list of arrays.");
    if (nfldTot < 4)
      PyErr_SetString(PyExc_TypeError,
                      "extractPoint: no field found in array.");
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no], objs[no], fields[no], a2[no], a3[no], a4[no]); 
    return NULL;
  }

  FldArrayF* an = new FldArrayF(npts, nfldTot-3); // sans les coordonnees
  FldArrayF& field = *an; field.setAllValuesAtNull();
  FldArrayF coord(npts, 3);
  for (E_Int j = 0; j < npts; j++)
  {
    PyObject* tplf = PyList_GetItem(listPts,j); // on recupere les listes des elements
    for (E_Int i = 1; i <= 3; i++)
    {
      PyObject* tplij = PyTuple_GetItem(tplf, i-1);// on recupere les elements
      E_Float Pt = PyFloat_AsDouble(tplij); // conversion en double 
      coord(j,i) = Pt;      
    }  
  }
  E_Int nzonesS = 0; E_Int nzonesU = 0;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; vector<E_Int> poscs;
  vector<void*> a5;
  for (E_Int no = 0; no < nzones; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(varString[no]); posc++;
    if (a4[no] == NULL) nzonesU++;
    else nzonesS++;
    posxs.push_back(posx); 
    posys.push_back(posy); 
    poszs.push_back(posz); 
    poscs.push_back(posc); 
    a5.push_back(NULL);// PAS DE CONNECTIVITE ELTS/ELTS VOISINS
  }

  // Liste des interpDatas
  vector<E_Int> nis; vector<E_Int> njs; vector<E_Int> nks;
  vector<FldArrayI*> cnt;
  vector<K_INTERP::InterpData*> interpDatas;
  // creation des interpDatas
  if (hook == Py_None)
  {
    E_Int isBuilt;
    for (E_Int no = 0; no < nzones; no++)
    {
      K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(
        fields[no]->getSize(), 
        fields[no]->begin(posxs[no]),
        fields[no]->begin(posys[no]),
        fields[no]->begin(poszs[no]),
        a2[no], a3[no], a4[no], isBuilt);
      if ( isBuilt == 1 ) interpDatas.push_back(adt);
      else 
      {
        for (size_t noi = 0; noi < interpDatas.size(); noi++)
          delete interpDatas[noi];
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no], objs[no], fields[no], a2[no], a3[no], a4[no]); 
      }
    }
  }
  else //if (hook != Py_None) // hook fourni
  {
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    void** packet = (void**) PyCObject_AsVoidPtr(hook);
#else
    void** packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
    E_Int s1 = fields.size();
    for (E_Int i = 0; i < s1; i++) 
      interpDatas.push_back((K_INTERP::InterpAdt*)(packet[i+1])); 
  }
  
  // cas seulement non structure : ncf a 4 (minimum)
  if (nzonesU != 0)
  {
    if (interpOrder != 2) printf("Warning: extractPoint: interpolation order is 2 for tetra arrays.\n");
    if (nzonesS == 0) ncf = 4;
  }
  
  /* Interpolation de la liste des pts */
  FldArrayI indi(1); FldArrayF cf(ncf);
  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
  E_Int posx, posy, posz;
  E_Int nof = 1;
  E_Int type;
  FldArrayF fLoc(1,nfldTot);
  E_Float voli;
  E_Int noblk;
  for (E_Int i = 0; i < npts; i++)
  {
    nof = 1;    
    // Recherche de la cellule d'interpolation
    short ok = K_INTERP::getInterpolationCell(
      xt[i], yt[i], zt[i],
      interpDatas, fields,
      a2, a3, a4, a5, posxs, posys, poszs, poscs,
      voli, indi, cf, type, noblk, interpType, 0, 0);    
    if (ok < 1)
    {
      ok = K_INTERP::getExtrapolationCell(
        xt[i], yt[i], zt[i],
        interpDatas, fields,
        a2, a3, a4, a5, posxs, posys, poszs, poscs,
        voli, indi, cf, type, noblk, interpType, 0, 0, 
        constraint, extrapOrder); 
    }      
    if (ok < 1)
    {
      printf("Warning: extractPoint: cannot interpolate point %5f %5f %5f.\n", xt[i], yt[i], zt[i]);
      // field est deja a zero : on ne fait rien
    }   
    else 
    {
      noblk = noblk-1;
      E_Int ind = 0;
      K_INTERP::compInterpolatedValues(indi.begin(), cf, *fields[noblk],
                                       a2[noblk], a3[noblk], a4[noblk], 
                                       ind, type, fLoc);
      posx = posxs[noblk]; posy = posys[noblk]; posz = poszs[noblk];
      for (E_Int v = 1; v <= nfldTot; v++)
      {
        if ( v != posx && v != posy && v != posz )
        {field(i, nof) = fLoc(0,v); nof++;}
      }
    }
  }// fin parcours des pts a interpoler
  
  // nettoyage...
  for (E_Int no = 0; no < nzones; no++)
  {
    if (hook == Py_None) delete interpDatas[no];
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);     
  }
  // Construction d'un i-array de sortie
  char* tmpStr;
  if (varString.size() > 0) 
  {
    tmpStr = new char [strlen(varString[0])+1];
    strcpy(tmpStr, varString[0]);
  }
  else
  {
    tmpStr = new char [2];
    tmpStr[0] = '\0';
  }

  vector<char*> vars; K_ARRAY::extractVars(tmpStr, vars);
  E_Int varsSize = vars.size();
  char* varStringOut = new char [strlen(tmpStr)+4];
  
  for (E_Int v = 0; v < varsSize; v++)
  {
    char*& var0 = vars[v];
    if (strcmp(var0, "x") != 0 && strcmp(var0, "CoordinateX") != 0 &&
        strcmp(var0, "y") != 0 && strcmp(var0, "CoordinateY") != 0 &&
        strcmp(var0, "z") != 0 && strcmp(var0, "CoordinateZ") != 0 )
    {
      strcpy(varStringOut, var0);
    }
  }
  for (E_Int v = 0; v < varsSize; v++) delete vars[v];
  PyObject* tpl = K_ARRAY::buildArray(field, varStringOut, npts, 1, 1);
  delete an; delete [] tmpStr; delete [] varStringOut;
  return tpl;
}
