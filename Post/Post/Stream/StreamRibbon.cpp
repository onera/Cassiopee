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

// compute stream ribbons

# include "stream.h"
# include <math.h>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

extern "C"
{  
  void k6compmeancurlofstructcell_(
    const E_Int& ind0, const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* u, const E_Float* v, const E_Float* w, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& rotu, E_Float& rotv, E_Float& rotw);
  
  void  k6compmeancurloftetracell_(
    const E_Int& npts, const E_Int& ind1, const E_Int& ind2, 
    const E_Int& ind3, const E_Int& ind4, 
    const E_Float* u, const E_Float* v, const E_Float* w, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& rotu, E_Float& rotv, E_Float& rotw);
}
//=============================================================================
/* Cree un ruban de courant */
//=============================================================================
PyObject* K_POST::compStreamRibbon(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Float x0, y0, z0;
  E_Float n0x, n0y, n0z;
  PyObject* vectorNames;
  E_Int nStreamPtsMax;
  E_Float signe;
  if (!PYPARSETUPLE(args,
                    "O(ddd)(ddd)Odl", "O(ddd)(ddd)Odi",
                    "O(fff)(fff)Ofl", "O(fff)(fff)Ofi",
                    &arrays, &x0, &y0, &z0, &n0x, &n0y, &n0z, &vectorNames, &signe, &nStreamPtsMax))
  {
      return NULL;
  }
  // Check every array in arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "streamRibbon: first argument must be a list.");
    return NULL;
  }

  //extraction of the 3 components of the vector used in the streamline computation 
  vector<char*> vnames;
  if (PyList_Check(vectorNames) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamRibbon: 3 variables defined as component of the streamline vector must be defined.");
    return NULL; 
  }

  E_Int sizeOfVector = PyList_Size(vectorNames);
  if (sizeOfVector != 3 ) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamRibbon: vector must be defined by 3 components.");
    return NULL;
  }
  for (int i = 0; i < PyList_Size(vectorNames); i++)
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
      char* str = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0));
      vnames.push_back(str);
    }
#endif
    else  
    {
      PyErr_SetString(PyExc_TypeError,
                      "streamRibbon: vector component name must be a string.");
      return NULL;
    }
  }
  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit0; 
  vector<E_Int> njt0; 
  vector<E_Int> nkt0;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs0, obju0;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int nfld = -1;
  E_Int isOk = K_ARRAY::getFromArrays(arrays, resl, structVarString, unstrVarString,
                                      structF, unstrF, nit0, njt0, nkt0, 
                                      cnt, eltType, objs0, obju0,
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
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
                    "streamRibbon: invalid list of arrays.");
    for (unsigned int nos = 0; nos < objs0.size(); nos++)
      RELEASESHAREDS(objs0[nos], structF[nos]);
    for (unsigned int nos = 0; nos < obju0.size(); nos++)
      RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
   
    return NULL;
  }
  for (unsigned int no = 0; no < structVarString.size(); no++)
  {
    if (nit0[no] == 1 || njt0[no] == 1 || nkt0[no] == 1 ) 
    {     
      PyErr_SetString(PyExc_TypeError,
                      "streamRibbon: structured arrays must be 3D.");
      for (unsigned int nos = 0; nos < objs0.size(); nos++)
        RELEASESHAREDS(objs0[nos], structF[nos]);
      for (unsigned int nos = 0; nos < obju0.size(); nos++)
        RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);      
      return NULL;
    }
  }
  for (unsigned int no = 0; no < unstrVarString.size(); no++)
  {
    if (strcmp(eltType[no],"NODE") == 0 || strcmp(eltType[no],"BAR") == 0 ||  
        strcmp(eltType[no],"TRI") == 0 || strcmp(eltType[no],"QUAD") == 0 )
    {     
      PyErr_SetString(PyExc_TypeError,
                      "streamRibbon: unstructured arrays must be 3D.");
      for (unsigned int nos = 0; nos < objs0.size(); nos++)
        RELEASESHAREDS(objs0[nos], structF[nos]);
      for (unsigned int nos = 0; nos < obju0.size(); nos++)
        RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);      
      return NULL;
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
    K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(
      structF[no]->getSize(), 
      structF[no]->begin(posx),
      structF[no]->begin(posy),
      structF[no]->begin(posz),
      &nit0[no], &njt0[no], &nkt0[no], isBuilt);
    nis1.push_back(nit0[no]);
    njs1.push_back(njt0[no]);
    nks1.push_back(nkt0[no]);
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
    K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(
      unstrF[no]->getSize(), 
      unstrF[no]->begin(posx),
      unstrF[no]->begin(posy),
      unstrF[no]->begin(posz),
      cnt[no], NULL, NULL, isBuilt);
    unstrF2.push_back(unstrF[no]); cnt2.push_back(cnt[no]);
    unstrInterpDatas2.push_back(adt);
    unstrVarString2.push_back(unstrVarString[no]);
    eltType2.push_back(eltType[no]);
  }   

  E_Int structSize =  structInterpDatas1.size();
  E_Int unstrSize = unstrInterpDatas2.size();
  E_Int interpDatasSize = structSize + unstrSize;
  
  if (interpDatasSize == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "streamRibbon : no interpData built.");
    for (unsigned int nos = 0; nos < objs0.size(); nos++)
      RELEASESHAREDS(objs0[nos], structF[nos]);
    for (unsigned int nos = 0; nos < obju0.size(); nos++)
      RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
    return NULL;
  } 
  nit0.clear(); njt0.clear(); nkt0.clear();


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

  // seuls sont pris en compte les fields ayant les variables du vecteur
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
    if ( found != 1 ) 
    {
      for (unsigned int nos = 0; nos < objs0.size(); nos++)
        RELEASESHAREDS(objs0[nos], structF[nos]);
      for (unsigned int nos = 0; nos < obju0.size(); nos++)
        RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
      for (unsigned int nos = 0; nos < structInterpDatas1.size(); nos++)
        delete structInterpDatas1[nos];
      for (unsigned int nos = 0; nos < unstrInterpDatas2.size(); nos++)
        delete unstrInterpDatas2[nos];
      
      if ( found == -1) 
        PyErr_SetString(PyExc_ValueError,
                        "streamRibbon: no field corresponding to vector component.");
      else // found == -2
        PyErr_SetString(PyExc_ValueError,
                        "streamRibbon: only TETRA type is valid for unstructured mesh.");
      return NULL;
    }
  }
  //extract des variables du vecteur sur les maillages non structures 
  if (unstrSize > 0) 
  {
    E_Int found = extractVectorFromUnstrArrays(signe, posxu2, posyu2, poszu2, poscu2,
                                               unstrVarString2, unstrF2, cnt2,
                                               eltType2, unstrInterpDatas2,
                                               posxu, posyu, poszu, poscu,
                                               unstrVarStrings, unstrFields, connectu,
                                               eltTypes, unstrInterpDatas,
                                               unstrVector, vnames); 
    if ( found != 1 ) 
    {
      for (unsigned int nos = 0; nos < objs0.size(); nos++)
        RELEASESHAREDS(objs0[nos], structF[nos]);
      for (unsigned int nos = 0; nos < obju0.size(); nos++)
        RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
      for (unsigned int nos = 0; nos < structInterpDatas1.size(); nos++)
        delete structInterpDatas1[nos];
      for (unsigned int nos = 0; nos < unstrInterpDatas2.size(); nos++)
        delete unstrInterpDatas2[nos];
      
      if ( found == -1) 
        PyErr_SetString(PyExc_ValueError,
                        "streamRibbon: no field corresponding to vector component.");
      else // found == -2
        PyErr_SetString(PyExc_ValueError,
                        "streamRibbon: only TETRA type is valid for unstructured mesh.");
      return NULL;
    }
  }
  
  // Petit nettoyage intermediaire
  nis1.clear(); njs1.clear(); nks1.clear();
  posxs1.clear(); posys1.clear(); poszs1.clear(); poscs1.clear();
  posxu2.clear(); posyu2.clear(); poszu2.clear(); poscu2.clear();

  // Calcul des points du ruban  de courant
  FldArrayF streamPts1(nStreamPtsMax, nfld);
  FldArrayF streamPts2(nStreamPtsMax, nfld);
  E_Int isinterp = computeStreamRibbonElts(x0, y0, z0, n0x, n0y, n0z,
                                           structInterpDatas, structFields, structVector,
                                           nis, njs, nks, posxs, posys, poszs, poscs,
                                           unstrInterpDatas, unstrFields, unstrVector,
                                           connectu, posxu, posyu, poszu, poscu, 
                                           streamPts1, streamPts2);
  if ( isinterp == 0 ) {streamPts1.malloc(0); streamPts2.malloc(0);}

  // little cleaning
  E_Int structVectorSize = structVector.size();
  for (E_Int v = 0; v < structVectorSize; v++) delete structVector[v];
  E_Int unstrVectorSize = unstrVector.size();
  for (E_Int v = 0; v < unstrVectorSize; v++) delete unstrVector[v];
 
  E_Int npts1 = streamPts1.getSize();
  E_Int npts2 = streamPts2.getSize();
  E_Int npts = npts1 + npts2;
  if (npts < 3)
  { 
    for (unsigned int nos = 0; nos < objs0.size(); nos++)
      RELEASESHAREDS(objs0[nos], structF[nos]);
    for (unsigned int nos = 0; nos < obju0.size(); nos++)
      RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
    for (unsigned int nos = 0; nos < structInterpDatas1.size(); nos++)
      delete structInterpDatas1[nos];
    for (unsigned int nos = 0; nos < unstrInterpDatas2.size(); nos++)
      delete unstrInterpDatas2[nos];

    PyErr_SetString(PyExc_ValueError,
                    "streamRibbon: cannot create a ribbon.");
    return NULL;
  }

  /* Creation du tableau non structure */
  FldArrayF* streamPts = new FldArrayF(npts, nfld);
  FldArrayI* cn = new FldArrayI(npts-2,3);
  buildConnectivity(streamPts1, streamPts2, *streamPts, *cn);
  // Construction de l'array non structure de sortie 
  PyObject* tpl = K_ARRAY::buildArray(*streamPts, varStringOut, *cn, 2);
  
  //nettoyage...
  delete streamPts; delete cn;
  delete [] varStringOut;

  for (unsigned int nos = 0; nos < objs0.size(); nos++)
    RELEASESHAREDS(objs0[nos], structF[nos]);
  for (unsigned int nos = 0; nos < obju0.size(); nos++)
    RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
  
  for (unsigned int nos = 0; nos < structInterpDatas1.size(); nos++)
    delete structInterpDatas1[nos];
  for (unsigned int nos = 0; nos < unstrInterpDatas2.size(); nos++)
    delete unstrInterpDatas2[nos];
  return tpl;
}
//=============================================================================
/* Calcul des elements du ruban */
//=============================================================================
E_Int K_POST::computeStreamRibbonElts(
  E_Float xp, E_Float yp, E_Float zp, 
  E_Float nxp, E_Float nyp, E_Float nzp,
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
  FldArrayF& streamPts1, FldArrayF& streamPts2)
{
  E_Float dtmin = 1.e-10;

  E_Int nfld = streamPts1.getNfld();
  E_Int nitermax = streamPts1.getSize();
  E_Float dt;
  
  // Creation of interpolation data
  K_INTERP::InterpData::InterpolationType interpType = K_INTERP::InterpData::O2CF;
  E_Int ns = listOfStructInterpData.size();
  E_Int nu = listOfUnstrInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1; vector<void*> allA2;
  vector<void*> allA3; vector<void*> allA4;
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

  FldArrayF cf;
  FldArrayI indi;
  if (listOfStructVelocities.size() == 0) { cf.malloc(4); indi.malloc(1); }
  else //ordre 2 structure
  {cf.malloc(8); indi.malloc(1);}
  
  // Calcul des premiers pts X0 et X0' et de dt, et du vecteur vitesse de X0
  // et des normales
  E_Int niter = 0;
  E_Float up, vp, wp;
  short isok = initStreamRibbon(xp, yp, zp, listOfStructInterpData, 
                                listOfStructFields, listOfStructVelocities,
                                nis, njs, nks, posxs, posys, poszs, poscs, 
                                listOfUnstrInterpData, listOfUnstrFields,
                                listOfUnstrVelocities, connectu,
                                posxu, posyu, poszu, poscu,
                                up, vp, wp, dt, nxp, nyp, nzp, indi, cf, 
                                streamPts1, streamPts2, 
                                interpType);
  if (isok == 0)
  {
    printf("Warning: streamRibbon: initialization of the ribbon failed.\n");
    return 0;
  }

  niter++;
  /* On arrete le calcul des pts de la ligne de courant des que iter=nmax
     ou que le pt x(n+1) n est plus interpolable */
  // determination de X(1)
  E_Float thetap = 0.;
  E_Int next = 0;

  // Recalcul ici du bloc d'interpolation initial pour que noblkp soit 
  // un argument local
  E_Int type = 0; E_Int noblkp = 0; E_Float volp = 0.;
  E_Int typen = 0;
  short found = K_INTERP::getInterpolationCell( 
    xp, yp, zp, allInterpDatas,
    allFields, allA1, allA2, allA3, allA4,
    posxt, posyt, poszt, posct, 
    volp, indi, cf, type, noblkp, interpType);

  if ( found == 0 ) {printf("Warning: streamRibbon: initial point not interpolable.\n"); return 0;}
  // Recherche des points X(n+1), n > 0
  while (niter < nitermax)
  {
    dt = K_FUNC::E_max(dt, dtmin);
    found = updateStreamRibbonPoints(niter, noblkp, type, xp, yp, zp, up, vp, wp, 
                                     nxp, nyp, nzp, thetap,
                                     indi, cf, dt, streamPts1, streamPts2,
                                     listOfStructInterpData, listOfStructFields, 
                                     listOfStructVelocities,
                                     nis, njs, nks, posxs, posys, poszs, poscs,
                                     listOfUnstrInterpData, listOfUnstrFields,
                                     listOfUnstrVelocities,
                                     connectu, posxu, posyu, poszu, poscu, 
                                     interpType, next, typen);
    if (found == 0) // pt Xp+1 non interpolable
    { 
      streamPts1.reAllocMat(niter, nfld);
      streamPts2.reAllocMat(niter, nfld);
      return 1;
    }
    noblkp = next; type = typen;
    niter++;
  }
  streamPts1.reAllocMat(niter, nfld);
  streamPts2.reAllocMat(niter, nfld);
  
  return 1;
}

//=============================================================================
/* Initialisation de la normale a la ligne de courant : calcul de la normale
   Si echec : retourne -1 */
//=============================================================================
short K_POST::initNormalToStreamLine(
  const E_Float x0, const E_Float y0, const E_Float z0,
  const E_Float u0, const E_Float v0, const E_Float w0,
  E_Float& nx0, E_Float& ny0, E_Float& nz0)
{
  E_Float invu0 = sqrt(u0*u0+v0*v0+w0*w0);
  if (fEqualZero(invu0) == true) return 0;
  else invu0 = 1./invu0;
  
  E_Float t1x = u0 * invu0;
  E_Float t1y = v0 * invu0;
  E_Float t1z = w0 * invu0;
  E_Float t2x, t2y, t2z, t3x, t3y, t3z;

  // Calcul de la base (t1,t2,t3), t1 = u0
  if (fEqualZero(t1x) == false || 
      fEqualZero(t1y) == false )// rotation selon z
  {
    t2x =-t1y;
    t2y = t1x;
    t2z = t1z;
    t3x = t1y*t2z-t1z*t2x;
    t3y = t1z*t2x-t1x*t2z;
    t3z = t1x*t2y-t1y*t2x; 
  }
  else if (K_FUNC::fEqualZero(t1z) == false )// rotation selon y
  {
    t2x = t1x;
    t2y = t1z;
    t2z =-t1y;
    t3x = t1y*t2z-t1z*t2x;
    t3y = t1z*t2x-t1x*t2z;
    t3z = t1x*t2y-t1y*t2x;  
  }  
  else return 0; // normalement c est deja verifie plus haut

  // Selection de la direction de la normale
  E_Float ps1 = t1x*nx0 + t1y*ny0 + t1z*nz0;
  E_Float ps2 = t2x*nx0 + t2y*ny0 + t2z*nz0;
  E_Float ps3 = t3x*nx0 + t3y*ny0 + t3z*nz0;
  if (fEqual( E_abs(ps1 * invu0), 1) == true )
  {
    printf("Warning: streamRibbon: given normal is colinear to the initial velocity.");
    printf(" Initial normal set to a default value.\n");
    nx0 = ps2 * t2x;
    ny0 = ps2 * t2y;
    nz0 = ps2 * t2z;
    return 1;
  }
  else 
  {
    if (E_abs(ps2) >= E_abs(ps3) ) // selection de la direction
    {
      nx0 = ps2 * t2x; ny0 = ps2 * t2y; nz0 = ps2 * t2z;
    }
    else 
    {
      nx0 = ps3 * t3x; ny0 = ps3 * t3y; nz0 = ps3 * t3z;
    }
  }
  return 1;
}
//=============================================================================
/* Mise a jour du second pt X' du ruban a partir du point X et de la normale */
//=============================================================================
short K_POST::compSecondPoint(
  const E_Int nopt,
  const E_Float xp, const E_Float yp, const E_Float zp,
  const E_Float nxp, const E_Float nyp, const E_Float nzp,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  FldArrayF& streamPts2,  
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int ns = listOfStructInterpData.size();
  E_Int nu = listOfUnstrInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1; vector<void*> allA2;
  vector<void*> allA3; vector<void*> allA4;
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
  E_Float xp2 = xp + nxp;
  E_Float yp2 = yp + nyp;
  E_Float zp2 = zp + nzp;
  FldArrayF cf;
  FldArrayI indi;
  if (listOfStructFields.size() == 0) { cf.malloc(4); indi.malloc(1); }
  else //ordre 2 structure
  { cf.malloc(8); indi.malloc(1); }

  // pt 0: no du blk d'interpolation dans interpDatas: demarre a 1 
  E_Int type = 0; E_Int noblkp2 = 0; E_Float volp2 = 0.;
  short found = K_INTERP::getInterpolationCell( 
    xp2, yp2, zp2, allInterpDatas,
    allFields, allA1, allA2, allA3, allA4,
    posxt, posyt, poszt, posct, 
    volp2, indi, cf, type, noblkp2, interpType);
  
  if (found < 1) return 0;

  compStreamPtFields(nopt, xp2, yp2, zp2, noblkp2, type, 
                     indi, cf, nis, njs, nks,
                     posxs, posys, poszs, listOfStructFields,
                     posxu, posyu, poszu, listOfUnstrFields, connectu,
                     streamPts2, interpType);// maj de streamPts2
  return noblkp2;
}

//=============================================================================
/* Initialisation du ruban de courant: vitesse et dt calcules
   streamPts1: mise a jour pour le point X0
   streamPts2: mise a jour pour le point X0'
   Retourne le numero du bloc d'interpolation (demarre a 1) pour le nouveau 
   pt, et 0 si pas interpolable. */
//=============================================================================
short K_POST::initStreamRibbon(
  E_Float xp, E_Float yp, E_Float zp,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  E_Float& up, E_Float& vp, E_Float& wp, E_Float& dt,
  E_Float& nxp, E_Float& nyp, E_Float& nzp,
  FldArrayI& indi, FldArrayF& cf, 
  FldArrayF& streamPts1, FldArrayF& streamPts2, 
  K_INTERP::InterpData::InterpolationType interpType)
{
  // Donnees pour l'appel de la methode computeStreamLineElts
  FldArrayI* connectSurf = NULL;
  E_Float* xSurf = NULL;
  E_Float* ySurf = NULL;
  E_Float* zSurf = NULL;
  E_Int sizeSurf = 0;

  // Calcul du premier point et de dt
  initStreamLine(xp, yp, zp, listOfStructInterpData, listOfStructFields, 
                 listOfStructVelocities, nis, njs, nks, 
                 posxs, posys, poszs, poscs, 
                 listOfUnstrInterpData, listOfUnstrFields,
                 listOfUnstrVelocities, connectu,
                 posxu, posyu, poszu, poscu,
                 *connectSurf, xSurf, ySurf, zSurf, sizeSurf,
                 up, vp, wp, dt, indi, cf, streamPts1, 
                 interpType);
  // Calcul de N0 et X0': maj dans streamPts2
  short ok = initNormalToStreamLine(xp, yp, zp, up, vp, wp, nxp, nyp, nzp);
  if (ok == 0) 
  {
    printf("Warning: streamRibbon: initial velocity is equal to zero.\n");
    return 0;
  }
  ok = compSecondPoint(0, xp, yp, zp, nxp, nyp, nzp, 
                       listOfStructInterpData, listOfStructFields,
                       nis, njs, nks, posxs, posys, poszs, poscs,
                       listOfUnstrInterpData, listOfUnstrFields,
                       connectu, posxu, posyu, poszu, poscu,
                       streamPts2, interpType);
  if (ok == 0) // second point non interpolable
  {
    printf("Warning: streamRibbon: second starting point not interpolable: (%f, %f, %f).\n",xp+nxp,yp+nyp,zp+nzp);
    printf("Please modify the input normal magnitude.\n"); 
    return 0;
  }
  else return 1;
}
//=============================================================================
/* Calcule les coordonnees du pt X(n+1), ses cf et indi et la vitesse U(n+1)
   Si necessaire, le pas dt est modifie 
   Retourne le numero du bloc d'interpolation (demarre a 1) pour le 
   nouveau pt et le type, et 0 si pas interpolable */
//=============================================================================
short K_POST::updateStreamRibbonPoints(
  E_Int niter, E_Int noblkp, E_Int typep,
  E_Float& xp, E_Float& yp, E_Float& zp,
  E_Float& up, E_Float& vp, E_Float& wp,
  E_Float& nxp, E_Float& nyp, E_Float& nzp,
  E_Float& thetap,
  FldArrayI& indip, FldArrayF& cfp, E_Float& dt, 
  FldArrayF& streamPts1, FldArrayF& streamPts2,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int& noblkn, E_Int& type)
{
  E_Int ns = listOfStructInterpData.size();
  E_Int nu = listOfUnstrInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1; vector<void*> allA2;
  vector<void*> allA3; vector<void*> allA4;
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
  E_Float thetamax = 30.;
  E_Float xn, yn, zn, thetan;
  FldArrayI indin(indip.getSize());
  FldArrayF cfn(cfp.getSize());
  E_Float nxn, nyn, nzn;

  noblkn = 0;
  short isOk;
  short cnt = 0;
  short cntmax = 4;
  E_Int noblkn0;
  E_Float voln = 0.;
  type = 0;
  while (cnt < cntmax)
  {
    // 1- calcul du pt X(p+1): xp+1,yp+1,zp+1,thetap+1
    isOk = compRungeKutta4ForRibbon(
      noblkp, typep, xp, yp, zp, thetap, up, vp, wp, indip, cfp, dt, 
      xn, yn, zn, thetan, listOfStructInterpData, listOfStructFields, 
      listOfStructVelocities, nis, njs, nks, posxs, posys, poszs, poscs,
      listOfUnstrInterpData, listOfUnstrFields, listOfUnstrVelocities,
      connectu, posxu, posyu, poszu, poscu, interpType);
    if (isOk == 0) // pas de pt interpolable ds sous pas RK4
    {
      if (cnt < 3) dt = 0.5 * dt;
      else dt = dt*16;
    }
    else 
    {
      if (E_abs(thetan) > thetamax) dt = 0.5 * dt;
      else 
      {
        // 2- Cellule d'interpolation pour calculer U(p+1)
        short found = K_INTERP::getInterpolationCell(
          xn, yn, zn, allInterpDatas,
          allFields, allA1, allA2, allA3, allA4,
          posxt, posyt, poszt, posct, 
          voln, indin, cfn, type, noblkn, interpType);
        if (found < 1) dt = 0.5 * dt;
        else 
        {
          noblkn0 = noblkn-1;
          // calcul de U(p+1)
          FldArrayF un(3);
          if (noblkn0 < ns)
          {
            K_INTERP::compInterpolatedField(
              indin.begin(), cfn, *listOfStructVelocities[noblkn0],
              &nis[noblkn0], &njs[noblkn0], &nks[noblkn0],
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

          // calcul de N(p+1) et X(p+1)' 
          compNormalToStreamLine(xn, yn, zn, thetan,
                                 nxp, nyp, nzp, nxn, nyn, nzn);

          // Mise a jour de streamPts1
          compStreamPtFields(niter, xn, yn, zn, noblkn, type, indin, cfn, 
                             nis, njs, nks, posxs, posys, poszs,
                             listOfStructFields,
                             posxu, posyu, poszu, listOfUnstrFields, connectu,
                             streamPts1, interpType);
          // Maj de streamPts2
          short ok =
            compSecondPoint(niter, xn, yn, zn, nxn, nyn, nzn,
                            listOfStructInterpData, listOfStructFields,
                            nis, njs, nks, posxs, posys, poszs, poscs,
                            listOfUnstrInterpData, listOfUnstrFields,
                            connectu, posxu, posyu, poszu, poscu,
                            streamPts2, interpType);
          if ( ok == 0 ) // deuxieme point nonunused interpolable
            dt = 0.5 * dt;
          else 
          {
            xp = xn; yp = yn; zp = zn; thetap = thetan;
            up = un[0]; vp = un[1]; wp = un[2];
            nxp = nxn; nyp = nyn; nzp = nzn;
            indip = indin; cfp = cfn;
            // Recalcul le pas de temps base sur la taille de la cellule
            // ordre en dur ici
            compInitialStep(noblkn, type, indin,nis,njs,nks,posxs,posys,poszs,
                            listOfStructFields, listOfStructVelocities,
                            posxu, posyu, poszu, 
                            listOfUnstrFields, connectu, 
                            listOfUnstrVelocities, dt);
            return 1;
          }
        }
      }
    }
    cnt++;
  }
  return 0;
}

//=============================================================================
/* Calcul des coefficients de Runge-Kutta RK4 Xn = Xp + dt/6 * (k1+k2+k3+k4)
   avec comme X=(x,y,z,theta)
   retourne 0 si un sous-pas ne s'est pas passe correctement
   dX/dt = u(Xp) et dTheta/dt = (rot Un . Un/normUn)/2
*/
//=============================================================================
short K_POST::compRungeKutta4ForRibbon(
  const E_Int noblkp, E_Int typep,
  E_Float xp, E_Float yp, E_Float zp, E_Float thetap,
  E_Float up, E_Float vp, E_Float wp, 
  FldArrayI& indip, FldArrayF& cfp,
  E_Float& dt, E_Float& xn, E_Float& yn, E_Float& zn, E_Float& thetan,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int ns = listOfStructInterpData.size();
  E_Int nu = listOfUnstrInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1; vector<void*> allA2;
  vector<void*> allA3; vector<void*> allA4;
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
  E_Float dts6 = dt/6.;
  E_Float dts2 = dt/2.;

  /* 1- determination de k1 : k1 = f(xp)*/
  E_Float k1x = up; //coefficient "k1"
  E_Float k1y = vp;
  E_Float k1z = wp;
  E_Float k1t, k2t, k3t, k4t;
  short isok = getThetaRKCoef(noblkp, typep, up, vp, wp, indip, 
                              listOfStructVelocities, listOfStructFields, 
                              nis, njs, nks, posxs, posys, poszs, 
                              listOfUnstrVelocities, listOfUnstrFields, 
                              connectu, posxu, posyu, poszu, k1t);
  if (isok == 0) return 0;
 
  //2- determination de k2 : k2 = u(xp2), xp2 = xp + dts2 * k1
  E_Float xp2 = xp + dts2 * k1x;
  E_Float yp2 = yp + dts2 * k1y;
  E_Float zp2 = zp + dts2 * k1z;
  //Cellule d interpolation pour calculer Up2
  FldArrayI indi(indip.getSize());
  FldArrayF cf(cfp.getSize());
  E_Float voli = 0.;
  E_Int type = 0;
  E_Int noblk = 0;
  
  short found = K_INTERP::getInterpolationCell( 
    xp2, yp2, zp2, allInterpDatas,
    allFields, allA1, allA2, allA3, allA4,
    posxt, posyt, poszt, posct, 
    voli, indi, cf, type, noblk, interpType);
   
  if (found < 1) return 0; // pas de pt d'interpolation trouve
  //Calcul de Up,2
  FldArrayF up2(3);
  E_Int noblk0 = noblk-1;
  if (noblk0 < ns)
  {
    K_INTERP::compInterpolatedField(
      indi.begin(), cf, *listOfStructVelocities[noblk0],
      &nis[noblk0],&njs[noblk0],&nks[noblk0], 
      type, up2);
    
  }
  else // non structure
  {
    E_Int noblku = noblk0-ns;
    K_INTERP::compInterpolatedField(
      indi.begin(), cf, *listOfUnstrVelocities[noblku],
      connectu[noblku], NULL, NULL, 
      type, up2);
  }
  
  E_Float k2x = up2[0]; // "k2"
  E_Float k2y = up2[1];
  E_Float k2z = up2[2];
  isok = getThetaRKCoef(noblk, type, up2[0], up2[1], up2[2], indi, 
                        listOfStructVelocities, listOfStructFields, 
                        nis, njs, nks, posxs, posys, poszs, 
                        listOfUnstrVelocities, listOfUnstrFields, 
                        connectu, posxu, posyu, poszu, k2t);

  if (isok == 0) return 0;

  // 3- determination de k3: k3 = u(xp3), xp3 = xp + dts2 * k2  
  xp2 = xp + dts2 * k2x;
  yp2 = yp + dts2 * k2y;
  zp2 = zp + dts2 * k2z;

  // Cellule d interpolation pour calculer Up3
  found = K_INTERP::getInterpolationCell( xp2, yp2, zp2, allInterpDatas,
                                          allFields, allA1, allA2, allA3, allA4,
                                          posxt, posyt, poszt, posct, 
                                          voli, indi, cf, type, noblk, interpType);
  if (found < 1) return 0; // pas de pt d'interpolation trouve
  
  // Calcul de Up3
  noblk0 = noblk-1;
  if (noblk0 < ns)
  {
    K_INTERP::compInterpolatedField(
      indi.begin(), cf, *listOfStructVelocities[noblk0],
      &nis[noblk0],&njs[noblk0],&nks[noblk0], 
      type, up2);
    
  }
  else // non structure
  {
    E_Int noblku = noblk0-ns;
    K_INTERP::compInterpolatedField(
      indi.begin(), cf, *listOfUnstrVelocities[noblku],
      connectu[noblku], NULL, NULL,
      type, up2);
  }
  
  E_Float k3x = up2[0]; // "k3"
  E_Float k3y = up2[1];
  E_Float k3z = up2[2];

  isok = getThetaRKCoef(noblk, type, up2[0], up2[1], up2[2], indi, 
                        listOfStructVelocities, listOfStructFields, 
                        nis, njs, nks, posxs, posys, poszs, 
                        listOfUnstrVelocities, listOfUnstrFields, 
                        connectu, posxu, posyu, poszu, k3t);

  if (isok == 0) return 0;

  //4- determination de k4: k4 = u(xp4), xp4 = xp + dt * k3  
  xp2 = xp + dt * k3x;
  yp2 = yp + dt * k3y;
  zp2 = zp + dt * k3z;

  // Cellule d interpolation pour calculer Up4
  found = K_INTERP::getInterpolationCell( 
    xp2, yp2, zp2, allInterpDatas,
    allFields, allA1, allA2, allA3, allA4,
    posxt, posyt, poszt, posct, 
    voli, indi, cf, type, noblk, interpType);

  if (found < 1) return 0;// pas de pt d'interpolation trouve

  // Calcul de Up4
  noblk0 = noblk-1;
  if (noblk0 < ns)
  {
    K_INTERP::compInterpolatedField(
      indi.begin(), cf, *listOfStructVelocities[noblk0],
      &nis[noblk0],&njs[noblk0],&nks[noblk0], 
      type, up2);
  }
  else // non structure
  {
    E_Int noblku = noblk0-ns;
    K_INTERP::compInterpolatedField(
      indi.begin(), cf, *listOfUnstrVelocities[noblku],
      connectu[noblku], NULL, NULL, 
      type, up2);
  }
  
  E_Float k4x = up2[0]; // "k4"
  E_Float k4y = up2[1];
  E_Float k4z = up2[2];
  isok = getThetaRKCoef(noblk, type, up2[0], up2[1], up2[2], indi, 
                        listOfStructVelocities, listOfStructFields, 
                        nis, njs, nks, posxs, posys, poszs, 
                        listOfUnstrVelocities, listOfUnstrFields, 
                        connectu, posxu, posyu, poszu, k4t);
 
  if (isok == 0) return 0;
  
  // 5- Mise a jour de x(n+1)
  xn = xp + dts6 * (k1x+k2x+k3x+k4x);
  yn = yp + dts6 * (k1y+k2y+k3y+k4y);
  zn = zp + dts6 * (k1z+k2z+k3z+k4z);
  thetan = dts6 * (k1t+k2t+k3t+k4t);
  return 1;
}
//=============================================================================
/* Calcul du coefficient ki du sous pas i de Runge Kutta 4 pour l'equation 
   en theta : dtheta/dt = 0.5*(rotv . v/normv). */
//=============================================================================
short K_POST::getThetaRKCoef(
  const E_Int noblk, E_Int type,
  const E_Float up, const E_Float vp, const E_Float wp,
  const FldArrayI& indi, 
  vector<FldArrayF*>& listOfStructVelocities,
  vector<FldArrayF*>& listOfStructFields,  
  vector<E_Int>& nis, vector<E_Int>& njs,
  vector<E_Int>& nks, vector<E_Int>& posxs, 
  vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayI*>& connectu, vector<E_Int>& posxu, 
  vector<E_Int>& posyu, vector<E_Int>& poszu,       
  E_Float& kcoef)
{
  E_Int ns = listOfStructFields.size();
  //E_Int nu = listOfUnstrFields.size();
  E_Int noblk0 = noblk-1; 
  E_Float rotvx, rotvy, rotvz;
  if (noblk0 < ns) //structure
  {
    if (type != 2 && type != 3 && type != 5)
    {
      printf("Error: getThetaRKCoef: not a valid interptype: %d.\n", type);
      exit(0);
    }
    FldArrayF* velo = listOfStructVelocities[noblk0];
    FldArrayF* field = listOfStructFields[noblk0];
    E_Int ni = nis[noblk0];
    E_Int nj = njs[noblk0];
    E_Int nk = nks[noblk0];
    E_Int posx = posxs[noblk0];
    E_Int posy = posys[noblk0];
    E_Int posz = poszs[noblk0];

    // calcul de ind ordre 2 en dur
    E_Int ind = indi[0];

    /* Determination de f4(xp)=(rotv.v/normv)/2 */
    // 1-calcul de la vorticite    
    k6compmeancurlofstructcell_(
      ind, ni, nj, nk, velo->begin(1), velo->begin(2), velo->begin(3), 
      field->begin(posx), field->begin(posy), field->begin(posz),
      rotvx, rotvy, rotvz);
    
    // 2-calcul de v/normv
    E_Float normu = sqrt(up*up+vp*vp+wp*wp);
    if ( fEqualZero(normu) == true ) return 0;
    normu = 1./normu;
    E_Float sx = up*normu;
    E_Float sy = vp*normu; 
    E_Float sz = wp*normu;
    
    // 3-Calcul du coeff de runge kutta pour theta
    kcoef = rotvx*sx+rotvy*sy+rotvz*sz;
    return 1;
  }
  else // non structure
  {
    if (type != 4) 
    {
      printf("Error: getThetaRKCoef: not a valid interp type: %d.\n", type);
      exit(0);
    }
    noblk0 = noblk0-ns;// numero du bloc reel dans la liste des blocs non structures
    FldArrayI& cnEV = *(connectu[noblk0]);
    if (cnEV.getNfld() != 4)
    {
      printf("Error: getThetaRKCoef: unstructured zone must be TETRA.\n");
      exit(0);
    }
    FldArrayF* velo = listOfUnstrVelocities[noblk0];
    FldArrayF* field = listOfUnstrFields[noblk0];
    E_Int posx = posxu[noblk0];
    E_Int posy = posyu[noblk0];
    E_Int posz = poszu[noblk0];
    E_Int noet = indi[0];
    E_Int indA = cnEV(noet,1)-1;
    E_Int indB = cnEV(noet,2)-1;
    E_Int indC = cnEV(noet,3)-1;
    E_Int indD = cnEV(noet,4)-1;
    k6compmeancurloftetracell_( field->getSize(), indA, indB, indC, indD, 
                                velo->begin(1), velo->begin(2), velo->begin(3), 
                                field->begin(posx), field->begin(posy), field->begin(posz),
                                rotvx, rotvy, rotvz);
    // 2-calcul de v/normv
    E_Float normu = sqrt(up*up+vp*vp+wp*wp);
    if ( fEqualZero(normu) == true ) 
      return 0;
    normu = 1./normu;
    E_Float sx = up*normu;
    E_Float sy = vp*normu; 
    E_Float sz = wp*normu;
    
    // 3-Calcul du coeff de runge kutta pour theta
    kcoef = rotvx*sx+rotvy*sy+rotvz*sz;
    return 1;
  }
}
//-----------------------------------------------------------------------------
// Construction de la connectivite triangle a partir des deux lignes definies
// par field1 et field2 
//-----------------------------------------------------------------------------
void K_POST::buildConnectivity(FldArrayF& field1, FldArrayF& field2,
                               FldArrayF& field, FldArrayI& cn)
{
  E_Int npts1 = field1.getSize();
  E_Int npts2 = field2.getSize();
  E_Int nfld = field1.getNfld();
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    field(0,eq) = field1(0,eq); 
    field(1,eq) = field2(0,eq);
  }      

  E_Int nf = 2;
  E_Int ie1 = 1;
  E_Int ie2 = 1;
  while (ie1 < npts1 && ie2 < npts2)
  {
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      field(nf,eq) = field1(ie1,eq); 
      field(nf+1,eq) = field2(ie2,eq); 
    }
    ie1++; ie2++;
    nf = nf+2;
  }
  if (nf != field.getSize()) // si npts1 et npts2 differents
  {
    field.reAllocMat(nf,nfld);
  }

  // construction de la connectivite : demarre a 1 
  E_Int c = 0;// indice ds connect
  E_Int i = 1;//numero du pt correspondant ds field 
  
  while (c < nf-2)
  {
    cn(c,1) = i;
    cn(c,2) = i+1;
    cn(c,3) = i+2;
    c++;
    cn(c,1) = i+1;
    cn(c,2) = i+2;
    cn(c,3) = i+3;
    c++;
    i += 2;
  }
}

//=============================================================================
/* Calcul de la normale a la ligne de courant au point Xn */
//=============================================================================
void K_POST::compNormalToStreamLine(const E_Float xn, const E_Float yn, 
                                    const E_Float zn, const E_Float thetan, 
                                    const E_Float nxp, const E_Float nyp, 
                                    const E_Float nzp,
                                    E_Float& nxn, E_Float& nyn, E_Float& nzn)
{
  E_Float cost = cos(thetan);
  nxn = nxp * cost;
  nyn = nyp * cost;
  nzn = nzp * cost;
}
