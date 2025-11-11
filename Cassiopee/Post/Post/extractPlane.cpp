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
// Overset grid slice extraction 
# include <string.h>
# include <stdio.h>
# include "post.h"
# include <vector>
# include "cutPlane/cutPlane.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Extract a slice in the field:
   IN: array of fields for all blocks.
   IN: coefs: ax+by+cz+d=0 */
// ============================================================================
PyObject* K_POST::extractPlane(PyObject* self, PyObject* args)
{
  E_Float coefa, coefb, coefc, coefd;
  PyObject* listFields;
  E_Int order;
  if (!PYPARSETUPLE_(args, O_ TRRRR_ I_,
                    &listFields, &coefa, &coefb, &coefc, &coefd, &order))
  {
    return NULL;
  }

  // Check every array in listFields
  if (PyList_Check(listFields) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extractPlane: first argument must be a list.");
    return NULL;
  }

  if (fEqualZero(coefa) && fEqualZero(coefb) && fEqualZero(coefc)) 
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "extractPlane: please change equation coefficients defining plane.");
    return NULL;
  }

  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false; 
  E_Bool skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listFields, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  E_Int nfld = 0;
  if (structF.size() != 0) nfld = structF[0]->getNfld();
  else  if (unstrF.size() != 0) nfld = unstrF[0]->getNfld();

  if (isOk == -1 || nfld < 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "extractPlane: invalid list of arrays.");
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }

  // Interpolation type
  K_INTERP::InterpData::InterpolationType interpType;
  switch (order)
  {
    case 2:
      interpType = K_INTERP::InterpData::O2CF;
      break;
    case 3:
      interpType = K_INTERP::InterpData::O3ABC;
      break;
    case 5:
       interpType = K_INTERP::InterpData::O5ABC;
       break;
    default:
      printf("Warning: extractPlane: interpolation order set to 2nd order.\n");
      interpType = K_INTERP::InterpData::O2CF;
      break;
  }
  E_Int nzonesS = structF.size();
  E_Int nzonesU = unstrF.size();

  // InterpData structuree
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;  vector<E_Int> poscs;
  vector<K_INTERP::InterpData*> structInterpDatas;
  E_Int isBuilt;
  for (E_Int no = 0; no < nzonesS; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(structVarString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(structVarString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(structVarString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(structVarString[no]); posc++;
    posxs.push_back(posx); posys.push_back(posy); poszs.push_back(posz); poscs.push_back(posc); 
    K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(structF[no]->getSize(), 
                                                       structF[no]->begin(posx),
                                                       structF[no]->begin(posy),
                                                       structF[no]->begin(posz),
                                                       &nit[no], &njt[no], &nkt[no], isBuilt);
    if ( isBuilt == 1 ) structInterpDatas.push_back(adt);
    else
    {
      for (size_t noi = 0; noi < structInterpDatas.size(); noi++)
        delete structInterpDatas[noi];
      PyErr_SetString(PyExc_TypeError, "2D structured donor zones must be z=constant.");
      for (size_t nos = 0; nos < objs.size(); nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
      for (size_t nos = 0; nos < obju.size(); nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
  }


  // InterpData non structuree
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu; 
  vector<E_Int> poscu;
  vector<K_INTERP::InterpData*> unstrInterpDatas;
  for (E_Int no = 0; no < nzonesU; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(unstrVarString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(unstrVarString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(unstrVarString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(unstrVarString[no]); posc++;
    posxu.push_back(posx); posyu.push_back(posy); poszu.push_back(posz); poscu.push_back(posc); 
    K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(unstrF[no]->getSize(), 
                                                       unstrF[no]->begin(posx),
                                                       unstrF[no]->begin(posy),
                                                       unstrF[no]->begin(posz),
                                                       cnt[no], NULL, NULL, isBuilt);
    unstrInterpDatas.push_back(adt);
  }
 
  if (structInterpDatas.size() == 0 && unstrInterpDatas.size() == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "extractPlane: invalid arrays.");
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }

  /*-------------------------------------------------------------*/
  /* Search for intersection of each edge with the defined plane */
  /*-------------------------------------------------------------*/
  vector<FldArrayF*> vectOfIntersectPts; // Points intersecting the plane
  vector<FldArrayF*> vectOfVolArrays;
  vector<FldArrayI*> tagS; vector<FldArrayI*> tagU;
  E_Int nbZonesS = structF.size();
  E_Int nbZonesU = unstrF.size();
  E_Int indTab[8];
  E_Int ni, nj, nk, ninj, nic, njc, nkc, nicnjc, indcell, posx, posy, posz, posc, npts, nelts;
  E_Float sumCellN; E_Int nplus, nmoins;
  E_Int api = -1;
  // on preconditionne en mettant le tag a 1 pour les cellules intersectant le plan
  for (E_Int nob = 0; nob < nbZonesS; nob++)
  {
    ni = nit[nob]; nic = K_FUNC::E_max(ni-1, 1);
    nj = njt[nob]; njc = K_FUNC::E_max(nj-1, 1);
    nk = nkt[nob]; nkc = K_FUNC::E_max(nk-1, 1);
    ninj = ni*nj; nicnjc = nic*njc;
    posx = posxs[nob]; posy = posys[nob]; posz = poszs[nob]; posc = poscs[nob];   
    E_Float* xt = structF[nob]->begin(posx);
    E_Float* yt = structF[nob]->begin(posy);
    E_Float* zt = structF[nob]->begin(posz);
    npts = structF[nob]->getSize();
    if (api == -1) api = structF[nob]->getApi();

    FldArrayF field(npts,1); E_Float* fp = field.begin();
    for (E_Int ind = 0; ind < npts; ind++)
      fp[ind] = coefa*xt[ind]+coefb*yt[ind]+coefc*zt[ind]+coefd;
    FldArrayI* tagS1 = new FldArrayI(nic*njc*nkc); tagS1->setAllValuesAt(1); 
    E_Int* tagcp = tagS1->begin();
    if (posc == 0) 
    {
      for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
          for (E_Int i = 0; i < nic; i++)
          {
            indcell   = i+ j*nic + k*nicnjc;
            indTab[0] = i+ j*ni  + k*ninj;
            indTab[1] = indTab[0]+1;
            indTab[2] = indTab[0]+ni;
            indTab[3] = indTab[1]+ni;
            indTab[4] = indTab[0]+ninj;
            indTab[5] = indTab[1]+ninj;
            indTab[6] = indTab[2]+ninj;
            indTab[7] = indTab[3]+ninj;
            nplus = 0; nmoins = 0;
            for (E_Int noi = 0; noi < 8; noi++)
            {
              if (fp[indTab[noi]] >= 0) nplus++;
              if (fp[indTab[noi]] <= 0) nmoins++;
            }
            if (nplus == 0 || nmoins == 0) tagcp[indcell] = 0;// PAS DE CHANGEMENT DE SIGNE DS LA CELLULE                 
          }
    }
    else 
    {
      E_Float* cellN = structF[nob]->begin(posc);
      for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
          for (E_Int i = 0; i < nic; i++)
          {
            indcell   = i+ j*nic + k*nicnjc;
            indTab[0] = i+ j*ni  + k*ninj;
            indTab[1] = indTab[0]+1;
            indTab[2] = indTab[0]+ni;
            indTab[3] = indTab[1]+ni;
            indTab[4] = indTab[0]+ninj;
            indTab[5] = indTab[1]+ninj;
            indTab[6] = indTab[2]+ninj;
            indTab[7] = indTab[3]+ninj;
            sumCellN = 0.; nplus = 0; nmoins = 0;
            for (E_Int noi = 0; noi < 8; noi++)
            {
              sumCellN += cellN[indTab[noi]];
              if ( fp[indTab[noi]] >= 0 ) nplus++;
              if ( fp[indTab[noi]] <= 0 ) nmoins++;
            }
            if ( nplus == 0 || nmoins == 0 || sumCellN == 0.) tagcp[indcell] = 0;// PAS DE CHANGEMENT DE SIGNE DS LA CELLULE                 
          }
    }
    tagS.push_back(tagS1);
  }
  for ( E_Int v = 0; v < nbZonesU; v++)
  {
    posx = posxu[v]; posy = posyu[v]; posz = poszu[v];posc = poscu[v];
    E_Float* xt = unstrF[v]->begin(posx);
    E_Float* yt = unstrF[v]->begin(posy);
    E_Float* zt = unstrF[v]->begin(posz);
    npts = unstrF[v]->getSize(); nelts = cnt[v]->getSize();
    if (api == -1) api = unstrF[v]->getApi();
    FldArrayF field(npts); E_Float* fp = field.begin();
    for (E_Int ind = 0; ind < npts; ind++)
      fp[ind] = coefa*xt[ind]+coefb*yt[ind]+coefc*zt[ind]+coefd;
    FldArrayI* tagU1 = new FldArrayI(nelts); 
    tagU1->setAllValuesAt(1); E_Int* tagcp = tagU1->begin();
    E_Int nvert = cnt[v]->getNfld();
    FldArrayI& cnloc = *cnt[v];
    if ( posc > 0 )
    {
      E_Float* cellN = unstrF[v]->begin(posc);
      for (indcell = 0; indcell < nelts; indcell++)
      {
        sumCellN = 0.; nplus = 0; nmoins = 0;
        for (E_Int novert = 1; novert <= nvert; novert++)
        {
          E_Int indv = cnloc(indcell,novert)-1;
          sumCellN += cellN[indv];
          if ( fp[indv] >= 0 ) nplus++;
          if ( fp[indv] <= 0 ) nmoins++;
        }
        if ( nplus == 0 || nmoins == 0 || sumCellN == 0.) // PAS DE CHANGEMENT DE SIGNE DS LA CELLULE OU TS LES SOMMETS SONT MASQUES
          tagcp[indcell] = 0;
      }
    }
    else 
    {
      for (indcell = 0; indcell < nelts; indcell++)
      {
        nplus = 0; nmoins = 0;
        for (E_Int novert = 1; novert <= nvert; novert++)
        {
          E_Int indv = cnloc(indcell,novert)-1;
          if ( fp[indv] >= 0 ) nplus++;
          if ( fp[indv] <= 0 ) nmoins++;
        }
        if ( nplus == 0 || nmoins == 0) // PAS DE CHANGEMENT DE SIGNE DS LA CELLULE
          tagcp[indcell] = 0;
      }
    }
    tagU.push_back(tagU1);
  }
  if (api == -1) api = 1;
  compIntersectionWithPlane(coefa, coefb, coefc, coefd,
                            structInterpDatas, nit, njt, nkt,
                            posxs, posys, poszs, poscs, structF, 
                            tagS, unstrInterpDatas, cnt,
                            posxu, posyu, poszu, poscu, unstrF,
                            tagU,
                            vectOfIntersectPts, vectOfVolArrays,
                            interpType);
  for ( E_Int nob = 0; nob < nbZonesS; nob++)
    delete tagS[nob];
  for ( E_Int nob = 0; nob < nbZonesU; nob++)
    delete tagU[nob];
  //attention  ds vectOfIntersectPts: d'abord les intersectPts struct, puis ns
  //vectOfIntersectPts est un vecteur des points d'intersection pour toutes les
  //zones, y compris s'il n y a pas de pts d'intersection pour une zone donnee

  /*------------------------------------------------------*/
  /* Select the most accurate points in overlapping zones */
  /*------------------------------------------------------*/
  FldArrayF* selectedPts = new FldArrayF();
  if (vectOfIntersectPts.size() > 0)
  {
    selectPointsInOverlappingZones( 
      nit, njt, nkt, posxs, posys, poszs, poscs,
      structInterpDatas, structF,
      cnt, posxu, posyu, poszu, poscu, unstrInterpDatas, unstrF,
      vectOfIntersectPts, vectOfVolArrays, 
      *selectedPts, interpType);
  }

  E_Int vectOfIntersectPtsSize = vectOfIntersectPts.size();
  E_Int  vectOfVolArraysSize = vectOfVolArrays.size();
  for (E_Int i = 0; i < vectOfIntersectPtsSize; i++)
    delete vectOfIntersectPts[i];
  for (E_Int i = 0; i < vectOfVolArraysSize; i++)
    delete vectOfVolArrays[i];
  
  /*-----------------------------------------*/
  /* Make triangulation with selected points */
  /*-----------------------------------------*/
  npts = selectedPts->getSize();
  if (npts > 0)
  {
    FldArrayI* connect = new FldArrayI();

    makeTriangulation(coefa, coefb, coefc, coefd, 
                      nit, njt, nkt, posxs, posys, poszs, poscs, 
                      structF, structInterpDatas,
                      cnt, posxu, posyu, poszu, poscu, 
                      unstrF, unstrInterpDatas,
                      *selectedPts, *connect, interpType);
    char* tmpStr;
    if (structVarString.size() > 0) 
    {
      tmpStr = new char [strlen(structVarString[0])+1];
      strcpy(tmpStr, structVarString[0]);
    }
    else if (unstrVarString.size() > 0) 
    {
      tmpStr = new char [strlen(unstrVarString[0])+1];
      strcpy(tmpStr, unstrVarString[0]);
    }
    else
    {
      tmpStr = new char [2]; tmpStr[0] = '\0';
    }
        
    E_Int posc = K_ARRAY::isCellNatureField2Present(tmpStr);
    char* varStringOut = new char [strlen(tmpStr)+7];

    if (posc == -1) strcpy(varStringOut, tmpStr);
    else
    {
      vector<char*> vars;
      K_ARRAY::extractVars(tmpStr, vars); E_Int varsSize = vars.size();
      E_Int c = 0;
      for (E_Int v = 0; v < varsSize; v++)
      {
        if (v != posc)
        {
          if (c == 0) strcpy(varStringOut, vars[v]);
          else {strcat(varStringOut, ","); strcat(varStringOut, vars[v]);}
          c++;
        }
      }
      strcat(varStringOut, ",");
      strcat(varStringOut, vars[posc]); 
      for (size_t i = 0; i < vars.size(); i++) delete vars[i];
    }

    PyObject* t = K_ARRAY::buildArray3(*selectedPts, varStringOut,
                                       *connect, "TRI", api);
    delete selectedPts; delete connect;
    // delete des interp datas
    for ( E_Int nob = 0; nob < nbZonesS; nob++)
      delete structInterpDatas[nob];
    for ( E_Int nob = 0; nob < nbZonesU; nob++)
      delete unstrInterpDatas[nob];
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    delete [] tmpStr; delete [] varStringOut;
    return t;
  }
  else // intersection vide ou tous les pts sont masques 
  {
    PyErr_SetString(PyExc_ValueError,
                    "extractPlane: intersection is empty.");

    // Test intersection des bbox avec le plan
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    E_Float xmins, ymins, zmins, xmaxs, ymaxs, zmaxs;
    E_Float xminu, yminu, zminu, xmaxu, ymaxu, zmaxu;

    // xmins... initialises ds globalBoudingBox
    K_COMPGEOM::globalBoundingBox(
      posxs,posys,poszs,structF,xmins, ymins, zmins,xmaxs,ymaxs,zmaxs);
    K_COMPGEOM::globalBoundingBox(
      posxu,posyu,poszu,unstrF,xminu, yminu, zminu,xmaxu,ymaxu,zmaxu);
    xmin = E_min(xminu, xmins); xmax = E_max(xmaxu, xmaxs);
    ymin = E_min(yminu, ymins); ymax = E_max(ymaxu, ymaxs);
    zmin = E_min(zminu, zmins); zmax = E_max(zmaxu, zmaxs);
    
    printf("Info: Bounding box of all meshes is:\n");
    printf("Info: x is between: %f and %f\n", xmin, xmax);
    printf("Info: y is between: %f and %f\n", ymin, ymax); 
    printf("Info: z is between: %f and %f\n", zmin, zmax);

    // delete des interp datas
    for ( E_Int nob = 0; nob < nbZonesS; nob++)
      delete structInterpDatas[nob];
    for ( E_Int nob = 0; nob < nbZonesU; nob++)
      delete unstrInterpDatas[nob];
    for (size_t nos = 0; nos < objs.size(); nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (size_t nos = 0; nos < obju.size(); nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
}
