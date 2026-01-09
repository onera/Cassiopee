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

// compute stream lines
# include "stream.h"
# include <math.h>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Cree une ligne de courant a  partir d'un point (x0,y0,z0) et d'une liste
   de grilles definies par des arrays. Les grilles d'interpolation sont celles
   qui sont structurees, contenant les infos sur la vitesse, et ayant toutes
   les m�mes variables dans le m�me ordre.*/
//=============================================================================
PyObject* K_POST::compStreamLine(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  PyObject* surfArray;
  E_Float x0, y0, z0;
  PyObject* vectorNames;
  E_Int nStreamPtsMax;
  E_Float signe;
  if (!PYPARSETUPLE_(args, OO_ TRRR_ O_ R_ I_,
                    &arrays, &surfArray, &x0, &y0, &z0, &vectorNames, &signe, &nStreamPtsMax))
  {
    return NULL;
  }
  // Check every array in arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "streamLine: first argument must be a list.");
    return NULL;
  }
  //extraction of the 3 components of the vector used in the streamline computation 
  vector<char*> vnames;
  if (PyList_Check(vectorNames) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamLine: 3 variables defined as component of the streamline vector must be defined.");
    return NULL; 
  }
  E_Int sizeOfVector = PyList_Size(vectorNames);
  if (sizeOfVector != 3 ) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamLine: vector must be defined by 3 components.");
    return NULL;
  }
  for (Py_ssize_t i = 0; i < PyList_Size(vectorNames); i++)
  {
    PyObject* tpl0 = PyList_GetItem(vectorNames, i);
    if (PyString_Check(tpl0))
    {
      char* str = PyString_AsString(tpl0);
      vnames.push_back(str);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0)) 
    {
      char* str = (char*)PyUnicode_AsUTF8(tpl0);
      vnames.push_back(str);
    }
#endif
    else  
    {
      PyErr_SetString(PyExc_TypeError,
                      "streamLine: vector component name must be a string.");
      return NULL;
    }
  }

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false; 
  E_Bool skipDiffVars = true;
  E_Int nfld = -1;
  E_Int isOk = K_ARRAY::getFromArrays(arrays, resl, 
                                      structVarString, unstrVarString,
                                      structF, unstrF, nit, njt, nkt, 
                                      cnt, eltType, objs, obju, skipDiffVars,
                                      skipNoCoord, skipStructured, 
                                      skipUnstructured, true);

  char* varStringOut;
  if (structVarString.size() > 0)
  {
    varStringOut = new char [strlen(structVarString[0])+1];
    strcpy(varStringOut, structVarString[0]);
  }
  else if (unstrVarString.size() > 0)
  {
    varStringOut = new char [strlen(unstrVarString[0])+1];
    strcpy(varStringOut, unstrVarString[0]);
  }
  else
  {
    varStringOut = new char [2];
    varStringOut[0] = '\0';
  }
  nfld = K_ARRAY::getNumberOfVariables(varStringOut);
  
  if (isOk == -1 || nfld == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamLine: invalid list of arrays.");
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }

  // get surfArray (always unstructured)
  E_Float* xSurf = NULL;
  E_Float* ySurf = NULL;
  E_Float* zSurf = NULL;
  FldArrayI* cnSurf = NULL;
  E_Int sizeSurf = 0;

  E_Int imSurf, jmSurf, kmSurf;
  FldArrayF* f = NULL;
  char* eltTypeSurf;
  char* varStringSurf;
  E_Int ress = -1000;
  if ((PyList_Check(surfArray) != 0) && (PyList_Size(surfArray) != 0))
  {
    ress = K_ARRAY::getFromArray3(surfArray, varStringSurf, f, 
                                  imSurf, jmSurf, kmSurf, cnSurf, eltTypeSurf); 
    
    E_Int nfldSurf = K_ARRAY::getNumberOfVariables(varStringSurf);
    if (ress == -1 || nfldSurf == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "streamLine: invalid list of surface arrays.");
      for (size_t nos = 0; nos < objs.size(); nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
      for (size_t nos = 0; nos < obju.size(); nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      if (ress !=-1000) RELEASESHAREDU(surfArray, f, cnSurf);

      return NULL;
    }
    if (ress != -2) // res == -2 <=> not a valid number of elts in list (not an empty list)
    {
      // keep only coordinate in surfarrays
      E_Int posx = K_ARRAY::isCoordinateXPresent(varStringSurf);
      E_Int posy = K_ARRAY::isCoordinateYPresent(varStringSurf);
      E_Int posz = K_ARRAY::isCoordinateZPresent(varStringSurf);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        for (size_t nos = 0; nos < objs.size(); nos++)
          RELEASESHAREDS(objs[nos], structF[nos]);
        for (size_t nos = 0; nos < obju.size(); nos++)
          RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
        if (ress != -1000) RELEASESHAREDU(surfArray, f, cnSurf);
        
        PyErr_SetString(PyExc_TypeError, 
                        "streamLine: coordinates not found in surface arrays.");
        return NULL;
      }
      posx++; posy++; posz++;
      xSurf = f->begin(posx);
      ySurf = f->begin(posy);
      zSurf = f->begin(posz); 
      sizeSurf = f->getSize();
    }
   }

   // Build interpData 
   E_Int nzonesS = structF.size();
   E_Int nzonesU = unstrF.size();
   // InterpData structuree
   vector<E_Int> posxs1; vector<E_Int> posys1; vector<E_Int> poszs1; vector<E_Int> poscs1;
   vector<K_INTERP::InterpData*> structInterpDatas1;
   vector<FldArrayF*> structF1;
   vector<E_Int> nis1; vector<E_Int> njs1; vector<E_Int> nks1;
   vector<char*> structVarStrings1;
   E_Int isBuilt;
   for (E_Int no = 0; no < nzonesS; no++)
   {
     E_Int posx = K_ARRAY::isCoordinateXPresent(structVarString[no]); posx++;
     E_Int posy = K_ARRAY::isCoordinateYPresent(structVarString[no]); posy++;
     E_Int posz = K_ARRAY::isCoordinateZPresent(structVarString[no]); posz++;
     E_Int posc = K_ARRAY::isCellNatureField2Present(structVarString[no]); posc++;
     posxs1.push_back(posx); posys1.push_back(posy); poszs1.push_back(posz); poscs1.push_back(posc); 
     K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(structF[no]->getSize(), 
                                                        structF[no]->begin(posx),
                                                        structF[no]->begin(posy),
                                                        structF[no]->begin(posz),
                                                        &nit[no], &njt[no], &nkt[no], isBuilt);
     nis1.push_back(nit[no]);
     njs1.push_back(njt[no]);
     nks1.push_back(nkt[no]);
     structF1.push_back(structF[no]);
     structInterpDatas1.push_back(adt);
     structVarStrings1.push_back(structVarString[no]);
   }
   // InterpData non structuree
  vector<E_Int> posxu2; vector<E_Int> posyu2; vector<E_Int> poszu2; 
  vector<E_Int> poscu2;
  vector<K_INTERP::InterpData*> unstrInterpDatas2;
  vector<FldArrayI*> cnt2;
  vector<FldArrayF*> unstrF2;
  vector<char*> unstrVarString2;
  vector<char*> eltType2;
  for (E_Int no = 0; no < nzonesU; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(unstrVarString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(unstrVarString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(unstrVarString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(unstrVarString[no]); posc++;
    posxu2.push_back(posx); posyu2.push_back(posy); poszu2.push_back(posz); poscu2.push_back(posc); 
    K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(unstrF[no]->getSize(), 
                                                       unstrF[no]->begin(posx),
                                                       unstrF[no]->begin(posy),
                                                       unstrF[no]->begin(posz),
                                                       cnt[no], NULL, NULL, isBuilt);
    unstrF2.push_back(unstrF[no]); cnt2.push_back(cnt[no]);
    unstrInterpDatas2.push_back(adt);
    unstrVarString2.push_back(unstrVarString[no]);
    eltType2.push_back(eltType[no]);
  }
  E_Int structSize = structInterpDatas1.size();
  E_Int unstrSize = unstrInterpDatas2.size();
  E_Int interpDatasSize = structSize + unstrSize;
  if (interpDatasSize == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "streamLine: no interpData built.");
    if ( ress != -1000 ) RELEASESHAREDU(surfArray, f, cnSurf);
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  } 

  // declaration de donnees
  // structure
  vector<FldArrayF*> structVector;
  vector<E_Int> nis; vector<E_Int> njs; vector<E_Int> nks;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  vector<E_Int> poscs;
  vector<char*> structVarStrings;
  vector<FldArrayF*> structFields;
  vector<K_INTERP::InterpData*> structInterpDatas;
  // non structure
  vector<FldArrayF*> unstrVector;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu; 
  vector<E_Int> poscu;
  vector<char*> unstrVarStrings;
  vector<FldArrayF*> unstrFields;
  vector<FldArrayI*> connectu;
  vector<char*> eltTypes;
  vector<K_INTERP::InterpData*> unstrInterpDatas;

  // seuls sont pris en compte les champs correspondant au vecteur
  // ts les arrays traites doivent avoir le meme nb de champs
  if (structSize > 0) 
  {
    E_Int found = extractVectorFromStructArrays(signe, nis1, njs1, nks1,
                                                posxs1, posys1, poszs1, poscs1,
                                                structVarStrings1, structF1, 
                                                structInterpDatas1,
                                                nis, njs, nks, 
                                                posxs, posys, poszs, poscs,
                                                structVarStrings, structFields, 
                                                structInterpDatas,
                                                structVector, vnames);

    if (found != 1)
    {
      if ( ress != -1000 ) RELEASESHAREDU(surfArray, f, cnSurf);
      for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
        delete structInterpDatas1[nos];
      for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
        delete unstrInterpDatas2[nos];
      for (size_t nos = 0; nos < objs.size(); nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
      for (size_t nos = 0; nos < obju.size(); nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      if (found == -1) 
        PyErr_SetString(PyExc_ValueError,
                        "streamLine: no field corresponding to vector component.");
      else // found = -2 uniquement pour l instant
        PyErr_SetString(PyExc_ValueError,
                        "streamLine: only TETRA type is valid for unstructured mesh.");
      return NULL;
    }
  }
  // extract des variables du vecteur sur les maillages non structures 
  if (unstrSize > 0) 
  {
    E_Int found = extractVectorFromUnstrArrays(signe, posxu2, posyu2, poszu2, poscu2,
                                               unstrVarString2, unstrF2, cnt2,
                                               eltType2, unstrInterpDatas2,
                                               posxu, posyu, poszu, poscu,
                                               unstrVarStrings, unstrFields, connectu,
                                               eltTypes, unstrInterpDatas,
                                               unstrVector, vnames); 
    if (found != 1)
    {
      if ( ress != -1000 ) RELEASESHAREDU(surfArray, f, cnSurf);
      for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
        delete structInterpDatas1[nos];
      for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
        delete unstrInterpDatas2[nos];
      for (size_t nos = 0; nos < objs.size(); nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
      for (size_t nos = 0; nos < obju.size(); nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      if (found == -1) 
        PyErr_SetString(PyExc_ValueError,
                        "streamLine: no field corresponding to vector component.");
      else // found = -2 uniquement pour l instant
        PyErr_SetString(PyExc_ValueError,
                        "streamLine: only TETRA type is valid for unstructured mesh.");
      return NULL;
    }
  }

  // Petit nettoyage intermediaire
  nis1.clear(); njs1.clear(); nks1.clear();
  posxs1.clear(); posys1.clear(); poszs1.clear(); poscs1.clear();
  posxu2.clear(); posyu2.clear(); poszu2.clear(); poscu2.clear();

  // Calcul de la ligne de courant
  FldArrayF* streamPts = new FldArrayF(nStreamPtsMax, nfld);
  E_Int isinterp = computeStreamLineElts(
    x0, y0, z0, 
    structInterpDatas, structFields, structVector,
    nis, njs, nks, posxs, posys, poszs, poscs, 
    unstrInterpDatas, unstrFields, unstrVector,
    connectu, posxu, posyu, poszu, poscu,
    *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
    *streamPts);
  if ( isinterp == 0 ) streamPts->malloc(0);

  E_Int api = 1;
  size_t structVectorSize = structVector.size();
  size_t unstrVectorSize = unstrVector.size();
  if (structVectorSize > 0) api = structVector[0]->getApi();
  else if (unstrVectorSize > 0) api = unstrVector[0]->getApi();

  // little cleaning
  for (size_t v = 0; v < structVectorSize; v++) delete structVector[v];
  for (size_t v = 0; v < unstrVectorSize; v++) delete unstrVector[v];

  // Build array : liste de pts, structure
  E_Int npts = streamPts->getSize(); 
  if (npts < 2)
  { 
    if ( ress != -1000 ) RELEASESHAREDU(surfArray, f, cnSurf);
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
      delete structInterpDatas1[nos];
    for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
      delete unstrInterpDatas2[nos];
    delete streamPts;
    PyErr_SetString(PyExc_ValueError,
                    "streamLine: cannot create a line.");
    return NULL;
  }

  PyObject* tpl = K_ARRAY::buildArray3(*streamPts, varStringOut, npts, 1, 1, api);
  
  delete [] varStringOut;
  
  //nettoyage...
  delete streamPts; 
  if ( ress != -1000 ) RELEASESHAREDU(surfArray, f, cnSurf);
  for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
    delete structInterpDatas1[nos];
  for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
    delete unstrInterpDatas2[nos];
  for (size_t nos = 0; nos < objs.size(); nos++)
    RELEASESHAREDS(objs[nos], structF[nos]);
  for (size_t nos = 0; nos < obju.size(); nos++)
    RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

  return tpl;
}

//=============================================================================
// Calcul de la ligne de courant 
//=============================================================================
E_Int K_POST::computeStreamLineElts(
  E_Float xp, E_Float yp, E_Float zp, 
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
  vector<E_Int>& poscu,
  FldArrayI& connectSurf, E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
  FldArrayF& streamPts)
{
  E_Int nfld = streamPts.getNfld();
  
  // nitermax  : fourni par l'utilisateur
  E_Int nitermax = streamPts.getSize();
  E_Float dt;

  // Creation of interpolation data
  K_INTERP::InterpData::InterpolationType interpType = K_INTERP::InterpData::O2CF;
  FldArrayF cf;
  FldArrayI indi;

  if (listOfStructVelocities.size() == 0) 
  {cf.malloc(4); indi.malloc(1);}
  else //ordre 2 structure
  {cf.malloc(8); indi.malloc(1);}

  // Calcul du premier point et de dt
  E_Int niter = 0;
  E_Float up, vp, wp;
  E_Int isok = initStreamLine(xp, yp, zp, listOfStructInterpData, listOfStructFields, 
                              listOfStructVelocities, nis, njs, nks, 
                              posxs, posys, poszs, poscs, 
                              listOfUnstrInterpData, listOfUnstrFields,
                              listOfUnstrVelocities, connectu,
                              posxu, posyu, poszu, poscu,
                              connectSurf, xSurf, ySurf, zSurf, sizeSurf,
                              up, vp, wp, dt, indi, cf, streamPts, 
                              interpType);
  if (isok == 0)
  {
    printf("Warning: streamLine: initialization of the ribbon failed.\n");
    return 0;
  }
  E_Float dtmin = 1.e-10;
  /* On arrete le calcul des pts de la ligne de courant des que iter=nmax
     ou que le pt x(n+1) n est plus interpolable */
  // Determination de X(1)
  E_Int type = 0;
  short next = 
    updateStreamLinePoint(xp, yp, zp, up, vp, wp, type, indi, cf, dt,
                          listOfStructInterpData, listOfStructFields, 
                          listOfStructVelocities,
                          nis, njs, nks, posxs, posys, poszs, poscs, 
                          listOfUnstrInterpData, listOfUnstrFields,
                          listOfUnstrVelocities, connectu,
                          posxu, posyu, poszu, poscu,
                          connectSurf, xSurf, ySurf, zSurf, sizeSurf,
                          interpType);
  niter++;

  if ( next == 0 ) // pt 1 non interpolable
  {
    streamPts.reAllocMat(niter, nfld);
    return 0;
  }
  // Mise a jour pour le pt X1 de son champ streamPt
  compStreamPtFields(niter, xp, yp, zp, next, type, indi, cf, nis, njs, nks,
                     posxs, posys, poszs, listOfStructFields,
                     posxu, posyu, poszu, listOfUnstrFields, connectu,
                     streamPts, interpType);// maj de streamPt            
  niter++;

  // Recherche des points X(n+1), n > 0
  while (niter < nitermax)
  {
    dt = K_FUNC::E_max(dt,dtmin);
    next = updateStreamLinePoint(xp, yp, zp, up, vp, wp, type, indi, cf, dt,
                                 listOfStructInterpData, listOfStructFields, 
                                 listOfStructVelocities,
                                 nis, njs, nks, posxs, posys, poszs, poscs, 
                                 listOfUnstrInterpData, listOfUnstrFields,
                                 listOfUnstrVelocities, connectu,
                                 posxu, posyu, poszu, poscu,
                                 connectSurf, xSurf, ySurf, zSurf, sizeSurf,
                                 interpType);
    if (next == 0) // pt 1 non interpolable
    {
      streamPts.reAllocMat(niter, nfld);
      return 1;
    }

    compStreamPtFields(niter, xp, yp, zp, next, type, indi, cf, nis, njs, nks,
                       posxs, posys, poszs, listOfStructFields,
                       posxu, posyu, poszu, listOfUnstrFields, connectu,
                       streamPts, interpType);// maj de streamPt

    niter++;
  }
  streamPts.reAllocMat(niter,nfld);
  return 1;
}

//=============================================================================
/* Calcule les coordonnees du pt X(n+1), ses cf et indi et la vitesse U(n+1)
   Si necessaire, le pas dt est modifie.
   Retourne le numero du bloc d'interpolation (demarre a 1) 
   pour le nouveau pt, et 0 si pas interpolable. */
//=============================================================================
short K_POST::updateStreamLinePoint(
  E_Float& xp, E_Float& yp, E_Float& zp,
  E_Float& up, E_Float& vp, E_Float& wp,
  E_Int& type, FldArrayI& indip, FldArrayF& cfp,
  E_Float& dt, 
  vector<K_INTERP::InterpData*>& listOfStructInterpData,
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs,
  vector<E_Int>& nks, vector<E_Int>& posxs, 
  vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs, 
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData,
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu, 
  FldArrayI& connectSurf, E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int ns = listOfStructInterpData.size();
  E_Int nu = listOfUnstrInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1;
  vector<void*> allA2;
  vector<void*> allA3;
  vector<void*> allA4;
  E_Int nzonesTot = ns+nu;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt; 
  vector<E_Int> posct;
  posxt.reserve(nzonesTot);  // preallocate memory
  posxt.insert(posxt.end(), posxs.begin(), posxs.end()); 
  posxt.insert(posxt.end(), posxu.begin(), posxu.end()); 
  posyt.reserve(nzonesTot);  // preallocate memory
  posyt.insert(posyt.end(), posys.begin(), posys.end()); 
  posyt.insert(posyt.end(), posyu.begin(), posyu.end()); 
  poszt.reserve(nzonesTot);  // preallocate memory
  poszt.insert(poszt.end(), poszs.begin(), poszs.end()); 
  poszt.insert(poszt.end(), poszu.begin(), poszu.end());
  posct.reserve(nzonesTot);  // preallocate memory
  posct.insert(posct.end(), poscs.begin(), poscs.end()); 
  posct.insert(posct.end(), poscu.begin(), poscu.end());

  allFields.reserve(nzonesTot);  // preallocate memory
  allFields.insert(allFields.end(), listOfStructFields.begin(), listOfStructFields.end()); 
  allFields.insert(allFields.end(), listOfUnstrFields.begin(), listOfUnstrFields.end());
  
  allInterpDatas.reserve(nzonesTot);  // preallocate memory
  allInterpDatas.insert(allInterpDatas.end(), listOfStructInterpData.begin(), listOfStructInterpData.end()); 
  allInterpDatas.insert(allInterpDatas.end(), listOfUnstrInterpData.begin(), listOfUnstrInterpData.end());

  for (E_Int noz = 0; noz < ns; noz++)
  {
    allA1.push_back(&nis[noz]); 
    allA2.push_back(&njs[noz]); 
    allA3.push_back(&nks[noz]); 
    allA4.push_back(NULL);
  }
  for (E_Int noz = 0; noz < nu; noz++)
  {
    allA1.push_back(connectu[noz]); 
    allA2.push_back(NULL); 
    allA3.push_back(NULL); 
    allA4.push_back(NULL);
  }

  E_Float cosmin = cos(15.);
  E_Float xn, yn, zn;
  FldArrayI indin(indip.getSize()); FldArrayF cfn(cfp.getSize());
  FldArrayI tmpIndin(indip.getSize()); FldArrayF tmpCfn(cfp.getSize());

  E_Int noblkn=0;
  short isThetaValid=0; //theta compris entre 2deg et 15deg : 1
  short isOk;
  short cnt=0;
  short cntmax=4;
  E_Int noblkn0;
  short found=0;
  FldArrayF un(3);
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];

  while (isThetaValid == 0 && cnt < cntmax)
  {
    // Si la surface n'est pas nulle, on projette le point sur cette surface
    if (sizeSurf != 0)
    {
      K_COMPGEOM::projectOrtho(xp, yp, zp,
                               xSurf, ySurf, zSurf,
                               connectSurf,
                               xp, yp, zp,
                               p0, p1, p2, p);
    }

    // 1- calcul du pt X(n+1)
    isOk = 
      compRungeKutta4(xp, yp, zp, up, vp, wp, dt, xn, yn, zn,
                      listOfStructInterpData, listOfStructFields, 
                      listOfStructVelocities, nis, njs, nks, 
                      posxs, posys, poszs, poscs,
                      listOfUnstrInterpData, listOfUnstrFields,
                      listOfUnstrVelocities, connectu,
                      posxu, posyu, poszu, poscu,
                      connectSurf, xSurf, ySurf, zSurf, sizeSurf,
                      interpType);

    if (isOk == 0) // pas de pt interpolable dans sous pas RK4
    {
      if (cnt < 3) dt = 0.5 * dt;
      else dt = dt*16;
    }
    else 
    {
      // 2- Cellule d'interpolation pour calculer U(n+1)
      type = 0; noblkn = 0; E_Float voln = 0.;
      found = K_INTERP::getInterpolationCell(
        xn, yn, zn, allInterpDatas,
        allFields, allA1, allA2, allA3, allA4,
        posxt, posyt, poszt, posct, 
        voln, indin, cfn, tmpIndin, tmpCfn, type, noblkn, interpType);

      if (found < 1) {dt = 0.5 * dt; }
      else
      {
        // Si la surface n'est pas nulle, on projette le point sur cette surface
        if (sizeSurf != 0)
        {
          K_COMPGEOM::projectOrtho(xn, yn, zn,
                                   xSurf, ySurf, zSurf,
                                   connectSurf, 
                                   xn, yn, zn,
                                   p0, p1, p2, p);
        }
        noblkn0 = noblkn-1;
        // Calcul de Un
        if (noblkn0 < ns)
        {
          K_INTERP::compInterpolatedField(
            indin.begin(), cfn, *listOfStructVelocities[noblkn0],
            &nis[noblkn0],&njs[noblkn0],&nks[noblkn0], 
            type, un);
        }
        else // non structure
        {
          E_Int noblku = noblkn0-ns;
          K_INTERP::compInterpolatedField(
            indin.begin(), cfn, *listOfUnstrVelocities[noblku],
            connectu[noblku], NULL, NULL, 
            type, un);
        }
       
        // 3- determination de cos(theta)=up.un
        E_Float unx = un[0];
        E_Float uny = un[1];
        E_Float unz = un[2];
        E_Float costheta = up*unx + vp*uny + wp*unz;
        E_Float norm = (up*up+vp*vp+wp*wp)*(unx*unx+uny*uny+unz*unz);
        costheta = costheta/sqrt(norm);
        
        if (costheta < cosmin) // theta > 15 deg
          dt = 0.5 * dt; 
        else // theta compris entre 2 et 15 degres 
        {
          isThetaValid = 1;
          xp = xn; yp = yn; zp = zn;
          indip = indin;
          cfp = cfn;
          up = unx; vp = uny; wp = unz;
          // Recalcule le pas de temps base sur la taille de la cellule
          compInitialStep(noblkn, type, indin,nis, njs, nks, 
                          posxs, posys, poszs, 
                          listOfStructFields, listOfStructVelocities,
                          posxu, posyu, poszu, listOfUnstrFields, connectu,
                          listOfUnstrVelocities, dt);          
          return noblkn;
        }
      }
    }
    cnt++;
  }
  allInterpDatas.clear(); allFields.clear(); allA1.clear(); allA2.clear(); allA3.clear(); allA4.clear();
  posxt.clear(); posyt.clear(); poszt.clear(); posct.clear();
  return 0;
}
