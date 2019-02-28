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

// Interpolation on zones given a solution defined on a list of zones

# include <string.h>
# include <stdio.h>
# include "post.h"
# include "kcore.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/*  
    IN: arrays of fields for all blocks
    IN: arrays of mesh on which fields are extracted
    OUT: array of interpolated values (x0,y0,z0) + other fields on 
    interpolated block 
    Rq: on interpole des domaines fournis en entree
    Si les grilles se recouvrent, on choisit la cellule la plus petite.
    On utilise uniquement de l'interpolation et de l'extrapolation.
    Les points non interpolables/extrapolables ont un champ nul.
*/
// ============================================================================
PyObject* K_POST::extractMesh(PyObject* self, PyObject* args)
{
  // listFields: liste des domaines d'interpolation (x,y,z)+sol
  // array: maillage d'extraction
  PyObject* listFields; PyObject* arrays;
  E_Int interpOrder;
  E_Int extrapOrder;
  E_Float constraint;
  PyObject* hook;
  if (!PYPARSETUPLE(args,
                    "OOlldO", "OOiidO",
                    "OOllfO", "OOiifO",
                    &listFields, &arrays, &interpOrder, &extrapOrder, &constraint, &hook))
  {
      return NULL;
  }

  // Extract infos from extraction arrays
  vector<E_Int> res0;
  vector<char*> structVarString0; vector<char*> unstrVarString0;
  vector<FldArrayF*> structF0; vector<FldArrayF*> unstrF0;
  vector<E_Int> nit0; vector<E_Int> njt0; vector<E_Int> nkt0;
  vector<FldArrayI*> cnt0; vector<char*> eltType0;
  vector<PyObject*> objst0, objut0;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res0, structVarString0, unstrVarString0,
    structF0, unstrF0, nit0, njt0, nkt0, cnt0, eltType0, objst0, objut0, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns0 = structF0.size(); E_Int nu0 = unstrF0.size();
    
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "extractMesh: invalid list of extraction arrays.");
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst0[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);
    return NULL;
  }
  
  E_Int posxi, posyi, poszi;
  vector<E_Int> posxs0; vector<E_Int> posys0; vector<E_Int> poszs0;
  vector<E_Int> posxu0; vector<E_Int> posyu0; vector<E_Int> poszu0;

  // Verification de posxi, posyi, poszi dans arrays
  for (E_Int i = 0; i < ns0; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString0[i]); posxi++;
    posyi = K_ARRAY::isCoordinateYPresent(structVarString0[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(structVarString0[i]); poszi++;
    posxs0.push_back(posxi); posys0.push_back(posyi); poszs0.push_back(poszi); 
  }
  for (E_Int i = 0; i < nu0; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString0[i]); posxi++;
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString0[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString0[i]); poszi++;
    posxu0.push_back(posxi); posyu0.push_back(posyi); poszu0.push_back(poszi); 
  }
  /* Extraction et Verification des arrays listFields (grilles a interpoler)
     - les champs x,y,z doivent etre en premier
     - le nb de champs doit etre identique */
  if (PyList_Check(listFields) == 0)
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst0[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);
    PyErr_SetString(PyExc_TypeError, 
                    "extractMesh: first argument must be a list.");
    return NULL;
  }

  // Interpolation type 
  E_Int nindi; E_Int ncf;
  K_INTERP::InterpData::InterpolationType interpType;
  // attention le cas purement non structure est traite apres les interpDatas 
  switch (interpOrder)
  {
    case 2:
      interpType = K_INTERP::InterpData::O2CF;
      ncf = 8; nindi = 1;
      break;
    case 3: 
      interpType = K_INTERP::InterpData::O3ABC;
      ncf = 9; nindi = 1;
      break;
    case 5:
      interpType = K_INTERP::InterpData::O5ABC;
      ncf = 15; nindi = 1;
      break;
    default:
      printf("Warning: extractMesh: unknown interpolation order. Set to 2nd order.\n");
      interpType = K_INTERP::InterpData::O2CF;
      ncf = 8; nindi = 1;
      break;
  }
  // Extract infos from listFields
  vector<E_Int> resl;  vector<char*> varString;
  vector<FldArrayF*> fields; 
  vector<void*> a2; // ni,nj,nk ou cnt en NS
  vector<void*> a3; // eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  skipNoCoord = true; skipStructured = false;
  skipUnstructured = false; skipDiffVars = true;
  isOk = K_ARRAY::getFromArrays(
    listFields, resl, varString, fields, a2, a3, a4, objs,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = objs.size();
  if (isOk == -1)
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst0[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "extractMesh: invalid list of arrays.");
    return NULL;
  }
  E_Int nfldTot = fields[0]->getNfld();
  
  if (nfldTot <= 3)
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst0[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]); 
    return arrays;
  }
  E_Int nvars = nfldTot;
  E_Int nzonesS = 0; E_Int nzonesU = 0;
  vector<E_Int> posxa; vector<E_Int> posya; vector<E_Int> posza; 
  vector<E_Int> posca;
  vector<void*> a5;
  for (E_Int no = 0; no < nzones; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(varString[no]); posc++;
    if (a4[no] == NULL) nzonesU++;
    else nzonesS++;
    posxa.push_back(posx); 
    posya.push_back(posy); 
    posza.push_back(posz); 
    posca.push_back(posc); 
    a5.push_back(NULL); // PAS DE CONNECTIVITE ELTS/ELTS VOISINS
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
        fields[no]->begin(posxa[no]),
        fields[no]->begin(posya[no]),
        fields[no]->begin(posza[no]),
        a2[no], a3[no], a4[no], isBuilt);
      if ( isBuilt == 1 ) interpDatas.push_back(adt);
      else 
      {
        for (size_t noi = 0; noi < interpDatas.size(); noi++)
          delete interpDatas[noi];
        for (E_Int nos = 0; nos < ns0; nos++)
          RELEASESHAREDS(objst0[nos], structF0[nos]);
        for (E_Int nos = 0; nos < nu0; nos++)
          RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]); 
        PyErr_SetString(PyExc_TypeError,
                        "extractMesh: 2D structured donor zones must be z=constant.");
        return NULL;
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
  if (interpDatas.size() == 0)
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst0[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]); 
    PyErr_SetString(PyExc_TypeError,
                    "extractMesh: no ADT built.");
    return NULL;
  }

  // cas seulement non structure: ncf a 4 (minimum)
  if (nzonesU != 0)
  {
    if (interpOrder != 2) printf("Warning: extractMesh: interpolation order is 2 for tetra arrays.\n");
    if (nzonesS == 0) ncf = 4;
  }
  // variables de sortie
  char* varStringOut;
  if (varString.size() > 0)
  {
    varStringOut = new char [strlen(varString[0])+2];
    strcpy(varStringOut, varString[0]);
  }
   else
  {
    varStringOut = new char [2]; varStringOut[0] = '\0';
  }
  E_Int posxo = K_ARRAY::isCoordinateXPresent(varStringOut); posxo++;
  E_Int posyo = K_ARRAY::isCoordinateYPresent(varStringOut); posyo++; 
  E_Int poszo = K_ARRAY::isCoordinateZPresent(varStringOut); poszo++;

  // Recherche des cellules d'interpolation
  vector<FldArrayF*> structFields;
  vector<FldArrayF*> unstructFields; vector<FldArrayI*> cnout;
  E_Int type, noblk;
  for (E_Int v = 0; v < ns0; v++)
  {
    FldArrayF& f = *structF0[v];// receptor
    E_Int nbI = f.getSize();
    E_Float vol = K_CONST::E_MAX_FLOAT;
    FldArrayF* interp = new FldArrayF(nbI, nvars); //x,y,z + interpolated field
    interp->setAllValuesAtNull();
    posxi = posxs0[v]; posyi = posys0[v]; poszi = poszs0[v];

    // Interpolation from all interpDatas
    E_Float* xt = f.begin(posxi);
    E_Float* yt = f.begin(posyi);
    E_Float* zt = f.begin(poszi);  
    E_Float* xi = interp->begin(posxo);
    E_Float* yi = interp->begin(posyo);
    E_Float* zi = interp->begin(poszo);
#pragma omp parallel default(shared) private(vol, noblk, type) if (nbI > 50)
    {
      FldArrayI indi(nindi*2); FldArrayF cf(ncf);
#pragma omp for
      for (E_Int ind = 0; ind < nbI; ind++)
      {
        E_Float x = xt[ind]; E_Float y = yt[ind]; E_Float z = zt[ind];
        short ok = K_INTERP::getInterpolationCell(
          x, y, z, interpDatas, fields,
          a2, a3, a4, a5, posxa, posya, posza, posca,
          vol, indi, cf, type, noblk, interpType, 0, 0);
        if (ok != 1)
        {
          ok = K_INTERP::getExtrapolationCell(
            x, y, z, interpDatas, fields,
            a2, a3, a4, a5, posxa, posya, posza, posca,
            vol, indi, cf, type, noblk, interpType, 0, 0, 
            constraint, extrapOrder); 
        }
        if (noblk > 0)
        {
          noblk = noblk-1;
          K_INTERP::compInterpolatedValues(indi.begin(), cf, *fields[noblk],
                                           a2[noblk], a3[noblk], a4[noblk], 
                                           ind, type, *interp);      
        }
        xi[ind] = x; yi[ind] = y; zi[ind] = z;
      }
    }
    structFields.push_back(interp);
  }
  
  for (E_Int v = 0; v < nu0; v++)
  {
    FldArrayF& f = *unstrF0[v];
    E_Int nbI = f.getSize();
    E_Float vol = K_CONST::E_MAX_FLOAT;
    FldArrayF* interp = new FldArrayF(nbI, nvars); //x,y,z + interpolated field
    FldArrayI* cninterp = new FldArrayI(*cnt0[v]);
    interp->setAllValuesAtNull();
    posxi = posxu0[v]; posyi = posyu0[v]; poszi = poszu0[v];

    // Interpolation from all interpDatas
    E_Float* xt = f.begin(posxi);
    E_Float* yt = f.begin(posyi);
    E_Float* zt = f.begin(poszi);
  
    E_Float* xi = interp->begin(posxo);
    E_Float* yi = interp->begin(posyo);
    E_Float* zi = interp->begin(poszo);
#pragma omp parallel default(shared) private(vol, noblk, type) if (nbI > 50)
    {
      FldArrayI indi(nindi*2); FldArrayF cf(ncf);
#pragma omp for
      for (E_Int ind = 0; ind < nbI; ind++)
      {
        E_Float x = xt[ind]; E_Float y = yt[ind]; E_Float z = zt[ind];
        short ok = K_INTERP::getInterpolationCell(
          x, y, z,
          interpDatas, fields,
          a2, a3, a4, a5, posxa, posya, posza, posca,
          vol, indi, cf, type, noblk, interpType, 0, 0); 
  
        if (ok != 1)
          ok = K_INTERP::getExtrapolationCell(
            x, y, z,
            interpDatas, fields,
            a2, a3, a4, a5, posxa, posya, posza, posca,
            vol, indi, cf, type, noblk, interpType,
            0, 0, constraint, extrapOrder);    
        
        if (noblk > 0)
        {
          noblk = noblk-1;
          K_INTERP::compInterpolatedValues(indi.begin(), cf, *fields[noblk],
                                           a2[noblk], a3[noblk], a4[noblk], 
                                           ind, type, *interp);
        }
        xi[ind] = x; yi[ind] = y; zi[ind] = z;
      }
    }
    unstructFields.push_back(interp); cnout.push_back(cninterp);
  }
  for (E_Int no = 0; no < nzones; no++)
  {
    if (hook == Py_None) delete interpDatas[no];
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);     
  }
 
  for (E_Int nos = 0; nos < ns0; nos++)
    RELEASESHAREDS(objst0[nos], structF0[nos]);
  for (E_Int nos = 0; nos < nu0; nos++)
    RELEASESHAREDU(objut0[nos], unstrF0[nos], cnt0[nos]);

  //------------------------------------//
  // Construction de l'arrays de sortie //
  //------------------------------------//

  // if (ns0 != 0 && nu0 != 0) 
  //   printf("Warning: extractMesh: structured then unstructured arrays are stored.\n");

  // Build arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl;
  for (E_Int nos = 0; nos < ns0; nos++)
  {
    tpl = K_ARRAY::buildArray(*structFields[nos], varStringOut,
                              nit0[nos], njt0[nos], nkt0[nos]);
    delete structFields[nos]; 
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  for (E_Int nou = 0; nou < nu0; nou++)
  {
    tpl = K_ARRAY::buildArray(*unstructFields[nou], varStringOut, *cnout[nou],
                              -1, eltType0[nou]);
    
    delete unstructFields[nou]; delete cnout[nou];
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  delete [] varStringOut;
  return l;
}
