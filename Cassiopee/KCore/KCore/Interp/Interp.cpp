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

# include "Interp/Interp.h"
# include "Fld/FldArray.h"
#include "String/kstring.h"
# include "Def/DefFunction.h"
# include "Connect/connect.h"
# include "Metric/metric.h"

using namespace std;
using namespace K_FLD;
using namespace K_ARRAY;

extern "C"
{
  void k6compvoloftetracell_(const E_Int& npts, 
                             const E_Int& ind1, const E_Int& ind2, 
                             const E_Int& ind3, const E_Int& ind4,
                             const E_Float* xt, const E_Float* yt, 
                             const E_Float* zt, E_Float& vol);
}

E_Int K_INTERP::extractADTFromHooks(PyObject* allHooks, vector<K_INTERP::InterpData*>& interpDatas)
{
  E_Int nHooks = PyList_Size(allHooks);
  for (E_Int noh = 0; noh < nHooks; noh++)
  {
    PyObject* hook = PyList_GetItem(allHooks,noh);
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
      void** packet = (void**) PyCObject_AsVoidPtr(hook);
#else
      void** packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

    E_Int* type = (E_Int*)packet[0];
    if (type[0] != 1) return -1;// hook must define an ADT
    if (type[1] != 1) return -2;// one ADT per hook
    interpDatas.push_back((K_INTERP::InterpAdt*)(packet[1])); 
  }
  if ( interpDatas.size()==0 ) return 0;
  return 1;
}
//=============================================================================
/* IN: x,y,z: point a interpoler.
   IN: InterpDatas: les InterpData des domaines donneurs
   IN: fields: des champs des domaines donneurs
   IN: STRUCTURE: a1=nit, a2=njt, a3=nkt
   IN: NON STRUCTURE: a1=cEV, a4=cEEN ou NULL
   IN: posxt, posyt, poszt, posct: position de x,y,z,c dans les domaines 
   donneurs
   IN: interpType: type d'interpolation
   IN: nature: 0 (pas de cellN=0 dans la cellule d'interpolation)
               1 (pas de cellN=0 ni cellN=2 dans la cellule d'interpolation)
   IN: penalty: 0: les bords ne sont pas penalise dans la recherche
                1: les bords sont penalise (seront pris en dernier)
   OUT: voli: volume de la cellule d'interpolation
   OUT: donorIndices: indices de la cellule d'interpolation
   OUT: donorCoefs: coeff d'interpolation
   OUT: type: type d'interpolation reellement effectuee
        1: coincident, 2: O2CF, 3: O3ABC, 4: O2 avec indice de l'element pour la
        cellule d'interpolation, 5: O5ABC 
   OUT: noDonorBlk: no du bloc donneur utilise dans l'interpolation */
//=============================================================================
short K_INTERP::getInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_INTERP::InterpData*>& InterpDatas,
  vector<FldArrayF*>& fields,
  vector<void*>& a1, vector<void*>& a2, vector<void*>& a3, 
  vector<void*>& a4, 
  vector<E_Int>& posxt, vector<E_Int>& posyt, 
  vector<E_Int>& poszt, vector<E_Int>& posct,
  E_Float& voli, 
  FldArrayI& donorIndices, FldArrayF& donorCoefs,
  FldArrayI& tmpIndi, FldArrayF& tmpCf,
  E_Int& type, E_Int& noDonorBlk,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty)
{
  E_Int isExtrapolation = 0;
  E_Float sumCoefMin = K_CONST::E_MAX_FLOAT;
  E_Float best = K_CONST::E_MAX_FLOAT;
  
  // numero du bloc d'interpolation: initialise a 0 si pas trouve
  E_Int noblk = 0; 
  E_Int nzones = InterpDatas.size();
  short found = 0;
  E_Float vol;
  E_Int indcell, firstCorner, noei, ni, nj, nk, ninj;
  E_Int isBorder, tmpType, meshtype;// 1: structure , 2: non structure
  E_Int d = 0; E_Int nbNonZeroCf = 0; E_Int indNonZero=-1; E_Int indLoc; 
  E_Float coefLoc = 0.; E_Float coefNonZero=0.;
  //E_Int donorIndicesSize = donorIndices.getSize();
  //FldArrayI tmpIndi(donorIndicesSize);
  tmpIndi.setAllValuesAtNull();
  //FldArrayF tmpCf(donorCoefs.getSize());
  tmpCf.setAllValuesAtNull();
  FldArrayI* cnloc=NULL; //connectivite elts/vertex TETRA si non structure
  donorIndices.setAllValuesAtNull(); donorCoefs.setAllValuesAtNull();
  short foundSav = 0;
  E_Int posx0, posy0, posz0, posc0;

  for (E_Int noz = 0; noz < nzones; noz++)
  {
    posx0 = posxt[noz]; posy0 = posyt[noz]; 
    posz0 = poszt[noz]; posc0 = posct[noz];
    //if (a1[noz] == NULL && a2[noz] == NULL && a3[noz] == NULL)
    //{printf("Error: getInterpolationCell: not a valid donor zone.\n"); return -1;}
    
    if (a3[noz] != NULL) //structure
    {
      ni = *(E_Int*)a1[noz]; nj = *(E_Int*)a2[noz]; nk = *(E_Int*)a3[noz]; ninj = ni*nj;
      meshtype = 1;
    }
    else
    { 
      meshtype = 2;
      // a1[noz] est la connectivite elts/elts voisins si different de NULL
      cnloc = (FldArrayI*)a1[noz]; // non structure      
      ni = -1; nj = -1; nk = -1; ninj = -1;      
    }
    //nic = K_FUNC::E_max(1,ni-1);  njc = K_FUNC::E_max(1,nj-1);  
    FldArrayF& oneField = *fields[noz]; E_Float* cellNp = NULL;
    if (posc0 > 0) cellNp = oneField.begin(posc0);
    else cellNp = NULL;
    if (meshtype == 2 && interpType != K_INTERP::InterpData::O2CF) 
    {
      //printf("Warning: getInterpolationCell: interpolation order is 2 for unstructured arrays.\n");
      interpType = K_INTERP::InterpData::O2CF;
    }
    // Recherche de la cellule d'interpolation
    // on retourne isBorder = 1 si au bord
    found = getInterpolationData(
      x, y, z, InterpDatas[noz], a1[noz], a4[noz], meshtype, ni, nj, nk, 
      oneField.begin(posx0), oneField.begin(posy0), oneField.begin(posz0), 
      cellNp, isBorder, tmpType, tmpIndi, tmpCf, interpType, nature);
# include "commonTypesForExtrapAndInterp.h"

  }// fin parcours des zones
  voli = best;
  noDonorBlk = noblk;

  // check si tous les coefs sont nuls sauf 1 
  if (foundSav < 1) return foundSav;

# include "commonCheckCoefs.h"

  return foundSav;
}
//=============================================================================
/* 
   Meme routine que getInterpolationCell, mais on envoie
   les images de (x,y,z) projete sur les parois des differents
   domaines donneurs dans xt,yt,zt.
   Les pts de coordonnees (xt[noz],yt[noz],zt[noz]) sont interpoles depuis 
   l'InterpData InterpDatas[noz] uniquement. Ce cas est applique pour calculer
   les coefficients d'interpolation dans le cas double wall.
   Attention: l'ordre de (xt,yt,zt) et de InterpDatas doit etre le meme */
//=============================================================================
short K_INTERP::getInterpolationCellDW(
  E_Float* xt, E_Float* yt, E_Float* zt,
  vector<K_INTERP::InterpData*>& InterpDatas,
  vector<FldArrayF*>& fields,
  vector<void*>& a1, vector<void*>& a2, vector<void*>& a3, 
  vector<void*>& a4, 
  vector<E_Int>& posxt, vector<E_Int>& posyt, 
  vector<E_Int>& poszt, vector<E_Int>& posct,
  E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
  E_Int& type, E_Int& noDonorBlk,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty)
{
  E_Int isExtrapolation = 0;
  E_Float sumCoefMin = K_CONST::E_MAX_FLOAT;
  E_Float best = K_CONST::E_MAX_FLOAT;
  
  // numero du bloc d'interpolation: initialise a 0 si pas trouve
  E_Int noblk = 0; 
  E_Int nzones = InterpDatas.size();
  short found = 0;
  E_Float vol;
  E_Int indcell, firstCorner, noei, ni, nj, nk, ninj;
  E_Int isBorder, tmpType, meshtype;// 1: structure , 2: non structure
  E_Int d = 0; E_Int nbNonZeroCf = 0; E_Int indNonZero=-1; E_Int indLoc;
  E_Float coefLoc = 0.; E_Float coefNonZero = 0.;
  FldArrayI tmpIndi(donorIndices.getSize()); tmpIndi.setAllValuesAtNull();
  FldArrayF tmpCf(donorCoefs.getSize()); tmpCf.setAllValuesAtNull();
  FldArrayI* cnloc=NULL;//connectivite elts/vertex TETRA si non structure
  donorIndices.setAllValuesAtNull(); donorCoefs.setAllValuesAtNull();
  short foundSav = 0;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Float x = xt[noz]; E_Float y = yt[noz]; E_Float z = zt[noz];
    E_Int posx0 = posxt[noz]; E_Int posy0 = posyt[noz]; 
    E_Int posz0 = poszt[noz]; E_Int posc0 = posct[noz];
    //if (a1[noz] == NULL && a2[noz] == NULL && a3[noz] == NULL)
    //{printf("Error: getInterpolationCell: not a valid donor zone.\n"); return -1;}
    
    if (a3[noz] != NULL) //structure
    {
      ni = *(E_Int*)a1[noz]; nj = *(E_Int*)a2[noz]; nk = *(E_Int*)a3[noz]; ninj = ni*nj;
      meshtype = 1;
    }
    else
    { 
      meshtype = 2;
      // a1[noz] est la connectivite elts/elts voisins si different de NULL
      cnloc = (FldArrayI*)a1[noz]; // non structure      
      ni = -1; nj = -1; nk = -1; ninj = -1;      
    }
    //nic = K_FUNC::E_max(1,ni-1);  njc = K_FUNC::E_max(1,nj-1);  
    FldArrayF& oneField = *fields[noz]; E_Float* cellNp = NULL;
    if (posc0 > 0) cellNp = oneField.begin(posc0);
    else cellNp = NULL;
    if (meshtype == 2 && interpType != K_INTERP::InterpData::O2CF) 
    {
      //printf("Warning: getInterpolationCell: interpolation order is 2 for unstructured arrays.\n");
      interpType = K_INTERP::InterpData::O2CF;
    }
    // Recherche de la cellule d'interpolation
    // on retourne isBorder = 1 si au bord
    found = getInterpolationData(
      x, y, z, InterpDatas[noz], a1[noz], a4[noz], meshtype, ni, nj, nk, 
      oneField.begin(posx0), oneField.begin(posy0), oneField.begin(posz0), 
      cellNp, isBorder, tmpType, tmpIndi, tmpCf,
      interpType, nature);

#include "commonTypesForExtrapAndInterp.h"
  }
  voli = best;
  noDonorBlk = noblk;
  
  // check si tous les coefs sont nuls sauf 1 
  if ( foundSav < 1 ) return foundSav; 
  
# include "commonCheckCoefs.h"
  
  return foundSav;
}
//=============================================================================
/* IN: x,y,z: point a interpoler.
   IN: InterpDatas: les InterpData des domaines donneurs
   IN: fields: des champs des domaines donneurs
   IN: STRUCTURE: a1=nit, a2=njt, a3=nkt
   IN: NON STRUCTURE: a1=cEV, a4=cEEN ou NULL
   IN: posxt, posyt, poszt, posct: position de x,y,z,c dans les domaines 
   donneurs
   IN: interpType: type d'interpolation
   IN: nature: 0 (pas de cellN=0 dans la cellule d'interpolation)
               1 (pas de cellN=0 ni cellN=2 dans la cellule d'interpolation)
   IN: penalty: 0: les bords ne sont pas penalise dans la recherche
                1: les bords sont penalise (seront pris en dernier)
   IN: constraint: contrainte sur la somme des valeurs absolus des
   coeffs d'interpolation. Si cette somme est > a constraint,
   on n'extrapole pas de cette cellule.
   OUT: voli: volume de la cellule d'interpolation
   OUT: donorIndices: indices de la cellule d'interpolation
   OUT: donorCoefs: coeff d'interpolation
   OUT: type: type d'interpolation reellement effectuee
        0: echec, 1: coincident, 2: O2CF, 3: O3ABC, 
        4: O2 avec indice de l'element pour la cellule d'interpolation, 
        5: O5ABC 
   OUT: noDonorBlk: no du bloc donneur utilise dans l'interpolation  
*/
//=============================================================================
short K_INTERP::getExtrapolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_INTERP::InterpData*>& InterpDatas,
  vector<FldArrayF*>& fields,
  vector<void*>& a1, vector<void*>& a2, vector<void*>& a3, 
  vector<void*>& a4, 
  vector<E_Int>& posxt, vector<E_Int>& posyt, 
  vector<E_Int>& poszt, vector<E_Int>& posct,
  E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
  E_Int& type, E_Int& noDonorBlk,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty, E_Float constraint, E_Int extrapOrder)
{
  E_Int isExtrapolation = 1;
  E_Float sumCoefMin = K_CONST::E_MAX_FLOAT;
  E_Float best = K_CONST::E_MAX_FLOAT; 

  interpType = K_INTERP::InterpData::O2CF; // forced
  short foundSav=0;
  // numero du bloc d'interpolation: initialise a 0 si pas trouve
  E_Int noblk = 0; 
  E_Int nzones = InterpDatas.size();
  short found=0;
  E_Float vol;
  E_Int indcell, firstCorner, noei, ni, nj, nk, ninj;
  E_Int isBorder, meshtype, tmpType; // 1: structure , 2: non structure
  E_Int d = 0; E_Int nbNonZeroCf = 0; E_Int indNonZero=-1; E_Int indLoc;
  E_Float coefLoc = 0.; E_Float coefNonZero = 0.;
  FldArrayI tmpIndi(donorIndices.getSize()); tmpIndi.setAllValuesAtNull();
  FldArrayF tmpCf(donorCoefs.getSize()); tmpCf.setAllValuesAtNull();
  donorIndices.setAllValuesAtNull(); donorCoefs.setAllValuesAtNull();

  FldArrayI* cnloc=NULL; //connectivite elts/vertex TETRA si non structure
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Int posx0 = posxt[noz]; E_Int posy0 = posyt[noz]; E_Int posz0 = poszt[noz]; E_Int posc0 = posct[noz];
    if (a1[noz] == NULL && a2[noz] == NULL && a3[noz] == NULL)
    {printf("Error: getExtrapolationCell: not a valid zone.\n"); return -1;}

    if (a3[noz] != NULL) //structure
    {
      ni = *(E_Int*)a1[noz]; nj = *(E_Int*)a2[noz]; nk = *(E_Int*)a3[noz]; ninj = ni*nj;
      meshtype = 1;
    }
    else
    { 
      meshtype = 2;
      // a1[noz] est la connectivite elts/elts voisins si different de NULL
      if (a1[noz] != NULL) cnloc = (FldArrayI*)a1[noz];// non structure      
      ni = -1; nj = -1; nk = -1; ninj = -1;  
    }
    //nic = K_FUNC::E_max(1,ni-1);  njc = K_FUNC::E_max(1,nj-1); 
    FldArrayF& oneField = *fields[noz]; E_Float* cellNp = NULL;
    if (posc0 > 0) cellNp = oneField.begin(posc0);
    else cellNp = NULL;
    if (meshtype == 2 && interpType != K_INTERP::InterpData::O2CF) 
    {
      //printf("Warning: getExtrapolationCell: interpolation order is 2 for unstructured arrays.\n");
      interpType = K_INTERP::InterpData::O2CF;
    }
    // Recherche de la cellule d'interpolation
    // on retourne isBorder=1 si au bord
    found = getExtrapolationData(
      x, y, z, InterpDatas[noz], a1[noz], a4[noz], meshtype, ni, nj, nk, 
      oneField.begin(posx0), oneField.begin(posy0), oneField.begin(posz0), 
      cellNp, isBorder, tmpType, tmpIndi, tmpCf, 
      interpType, nature, constraint, extrapOrder);
#include "commonTypesForExtrapAndInterp.h"
  }
  voli = best;
  noDonorBlk = E_Int(noblk);

  // check si tous les coefs sont nuls sauf 1 
  if ( foundSav < 1 ) return foundSav; 

# include "commonCheckCoefs.h"

  return foundSav;
}
//=============================================================================
/*
  Idem routine getExtrapolationCell, mais on envoie
  les images de (x,y,z) projete sur les parois des differents
  domaines donneurs dans xt,yt,zt.
  Les pts de coordonnees (xt[noz],yt[noz],zt[noz]) sont extrapoles depuis 
  l'InterpData InterpDatas[noz] uniquement. Ce cas est applique pour calculer
  les coefficients d extrapolation dans le cas double wall 
  Attention: l'ordre de (xt,yt,zt) et de InterpDatas doit etre le meme */
//=============================================================================
short K_INTERP::getExtrapolationCellDW(
  E_Float* xt, E_Float* yt, E_Float* zt,
  vector<K_INTERP::InterpData*>& InterpDatas,
  vector<FldArrayF*>& fields,
  vector<void*>& a1, vector<void*>& a2, vector<void*>& a3, 
  vector<void*>& a4, 
  vector<E_Int>& posxt, vector<E_Int>& posyt, 
  vector<E_Int>& poszt, vector<E_Int>& posct,
  E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
  E_Int& type, E_Int& noDonorBlk,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty, E_Float constraint, E_Int extrapOrder)
{
  E_Int isExtrapolation = 1;
  E_Float sumCoefMin = K_CONST::E_MAX_FLOAT;
  E_Float best = K_CONST::E_MAX_FLOAT; 

  interpType = K_INTERP::InterpData::O2CF;
  short foundSav = 0;
  // numero du bloc d'interpolation: initialise a 0 si pas trouve
  E_Int noblk = 0; 
  E_Int nzones = InterpDatas.size();
  short found=0;
  E_Float vol;
  E_Int indcell, firstCorner, noei, ni, nj, nk, ninj;
  E_Int isBorder, meshtype, tmpType;// 1 : structure , 2 non structure
  E_Int d = 0; E_Int nbNonZeroCf = 0; E_Int indNonZero=-1; E_Int indLoc;
  E_Float coefLoc = 0.; E_Float coefNonZero = 0.;
  FldArrayI tmpIndi(donorIndices.getSize()); tmpIndi.setAllValuesAtNull();
  FldArrayF tmpCf(donorCoefs.getSize()); tmpCf.setAllValuesAtNull();
  donorIndices.setAllValuesAtNull(); donorCoefs.setAllValuesAtNull();

  FldArrayI* cnloc=NULL;//connectivite elts/vertex TETRA si non structure
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Float x = xt[noz]; E_Float y = yt[noz]; E_Float z = zt[noz];
    E_Int posx0 = posxt[noz]; E_Int posy0 = posyt[noz]; E_Int posz0 = poszt[noz]; E_Int posc0 = posct[noz];
    if (a1[noz] == NULL && a2[noz] == NULL && a3[noz] == NULL)
    {printf("Error: getExtrapolationCell: not a valid zone.\n"); return -1;}

    if (a3[noz] != NULL) //structure
    {
      ni = *(E_Int*)a1[noz]; nj = *(E_Int*)a2[noz]; nk = *(E_Int*)a3[noz]; ninj = ni*nj;
      meshtype = 1;
    }
    else
    { 
      meshtype = 2;
      // a1[noz] est la connectivite elts/elts voisins si different de NULL
      if (a1[noz] != NULL) cnloc = (FldArrayI*)a1[noz];// non structure      
      ni = -1; nj = -1; nk = -1; ninj = -1;
      
    }
    //nic = K_FUNC::E_max(1,ni-1);  njc = K_FUNC::E_max(1,nj-1); 
    FldArrayF& oneField = *fields[noz]; E_Float* cellNp = NULL;
    if (posc0 > 0) cellNp = oneField.begin(posc0);
    else cellNp = NULL;
    if (meshtype == 2 && interpType != K_INTERP::InterpData::O2CF) 
    {
      //printf("Warning: getExtrapolationCell: interpolation order is 2 for unstructured arrays.\n");
      interpType = K_INTERP::InterpData::O2CF;
    }
    // Recherche de la cellule d'interpolation
    // on retourne isBorder = 1 si au bord
    found = getExtrapolationData(
      x, y, z, InterpDatas[noz], a1[noz], a4[noz], meshtype, ni, nj, nk, 
      oneField.begin(posx0), oneField.begin(posy0), oneField.begin(posz0), 
      cellNp,
      isBorder, tmpType, tmpIndi, tmpCf, 
      interpType, nature, constraint, extrapOrder); 
#include "commonTypesForExtrapAndInterp.h"
  }
  voli = best;
  noDonorBlk = E_Int(noblk);

  // check si tous les coefs sont nuls sauf 1 
  if ( foundSav < 1 ) return foundSav; 
  
# include "commonCheckCoefs.h"

  return foundSav;
}
//=============================================================================
/* Get interpolationCell, mais uniquement a partir d'un seul maillage 
   donneur */
//=============================================================================
short K_INTERP::getInterpolationCell(
  E_Float x, E_Float y, E_Float z, K_INTERP::InterpData* InterpData, 
  FldArrayF* field, void* a1, void* a2, void* a3, void* a4, 
  E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  E_Float& voli,  
  FldArrayI& donorIndices, FldArrayF& donorCoefs,
  FldArrayI& tmpIndi, FldArrayF& tmpCf,
  E_Int& type, E_Int& noDonorBlk,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty)
{
  vector<K_INTERP::InterpData*> InterpDatas; InterpDatas.push_back(InterpData);
  vector<E_Int> posxt; posxt.push_back(posx); 
  vector<E_Int> posyt; posyt.push_back(posy); 
  vector<E_Int> poszt; poszt.push_back(posz); 
  vector<E_Int> posct; posct.push_back(posc); 
  vector<FldArrayF*> fields; fields.push_back(field);
  vector<void*> a1t; vector<void*> a2t; vector<void*> a3t; vector<void*> a4t;

  a1t.push_back(a1); a2t.push_back(a2); a3t.push_back(a3); a4t.push_back(a4);
  return K_INTERP::getInterpolationCell(
    x, y, z, InterpDatas, fields, a1t, a2t, a3t, a4t,
    posxt, posyt, poszt, posct, voli, 
    donorIndices, donorCoefs, tmpIndi, tmpCf,
    type, noDonorBlk, interpType, nature, penalty);
}
//=============================================================================
/* */
//=============================================================================
short K_INTERP::getExtrapolationCell(
  E_Float x, E_Float y, E_Float z, K_INTERP::InterpData* InterpData, 
  FldArrayF* field, void* a1, void* a2, void* a3, void* a4, 
  E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  E_Float& voli,  FldArrayI& donorIndices, FldArrayF& donorCoefs,
  E_Int& type, E_Int& noDonorBlk,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty, E_Float constraint, E_Int extrapOrder)
{
  vector<K_INTERP::InterpData*> InterpDatas; InterpDatas.push_back(InterpData);
  vector<E_Int> posxt; posxt.push_back(posx); 
  vector<E_Int> posyt; posyt.push_back(posy); 
  vector<E_Int> poszt; poszt.push_back(posz); 
  vector<E_Int> posct; posct.push_back(posc); 
  vector<FldArrayF*> fields; fields.push_back(field);
  vector<void*> a1t; vector<void*> a2t; vector<void*> a3t; vector<void*> a4t;

  a1t.push_back(a1); a2t.push_back(a2); a3t.push_back(a3); a4t.push_back(a4);
  return K_INTERP::getExtrapolationCell(
    x, y, z, InterpDatas, fields, a1t, a2t, a3t, a4t,
    posxt, posyt, poszt, posct, voli, donorIndices, donorCoefs, 
    type, noDonorBlk,
    interpType, nature, penalty, constraint, extrapOrder);
}
//=============================================================================
/* IN: a2,a3,a4: ni, nj, nk: dimensions de f0 dans le cas structure
       a2: cEV  si le donneur est non structure
   IN: f0: champ du bloc donneur pour l interpolation
   IN: indi: indice des pts de la molecule d interpolation selon le type 'type'
             peut etre defini par des sommets de la molecule donneuse
                              par l'indice de la cellule donneuse
   IN: cf: coefs d'interpolation associes 
   IN: type: permet de determiner la formule appliquee
   OUT: la valeur interpolee 
   Retourne -1 si erreur */
//=============================================================================
short K_INTERP::compOneInterpolatedValue(
  E_Int* indi, FldArrayF& cf,
  E_Float* f0, void* a2, void* a3, void* a4, 
  E_Int type, E_Float& val)
{
  E_Float* cfp = cf.begin();
  E_Int ind0; // nb de coefs par direction dans le cas OiABC
  E_Int nocf=0;
  E_Int ni, nj, ninj, i, j, k, indcell, noet, nvert;
  val = 0.;
  
  switch (type) 
  {
    case 0:
      nocf = indi[0];// nb de pts pour la formule
      for (E_Int no = 1; no <= nocf; no++)
      {
        ind0 = indi[no];
        val += cfp[no-1]*f0[ind0];
      }      
      break;

    case 1:// coincident
      ind0 = indi[0];
      val = cfp[0]*f0[ind0];
      break;
      
    case 2:
      if (a4 == NULL)// NS : indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compOneInterpolatedValue: type 2 not valid for unstructured interpolations.\n");
        return -1;
      }
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      //nk = *(E_Int*)a4;
      ninj = ni*nj;  
      ind0 = indi[0]; nocf = 0;
      k = ind0/(ninj);  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
      for (E_Int k0 = 0; k0 < 2; k0++)
        for (E_Int j0 = 0; j0 < 2; j0++)
          for (E_Int i0 = 0; i0 < 2; i0++)
          {
            E_Int ind = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
            val += cfp[nocf]*f0[ind];
            nocf+=1;
          }
      break;

    case 22:
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      ind0 = indi[0]; nocf = 0;
      j = ind0/ni; i = ind0-j*ni;
      for (E_Int j0 = 0; j0 < 2; j0++)
        for (E_Int i0 = 0; i0 < 2; i0++)
        {
          E_Int ind = (i+i0) + (j+j0)*ni;
          val += cfp[nocf]*f0[ind];
          nocf+=1;
        }
      break;

    case 3://formule directionnelle (ex O3ABC) : on stocke les indices i,j,k des sommets 
      if (a4 == NULL)// NS : indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compOneInterpolatedValue: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      //nk = *(E_Int*)a4;    
      ninj = ni*nj;  
      ind0 = indi[0]; nocf = 0;
      k = ind0/(ninj);  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
      for (E_Int k0 = 0; k0 < 3; k0++)
        for (E_Int j0 = 0; j0 < 3; j0++)
          for (E_Int i0 = 0; i0 < 3; i0++)
          {
            E_Int ind = (i+i0)+(j+j0)*ni+(k+k0)*ninj;
            val += cfp[i0]*cfp[j0+3]*cfp[k0+6]*f0[ind];  
          }
      break;

    case  4:// on stocke le centre de la cellule et on interpole des sommets    
      if (a4 != NULL)// Structure 
      {
        indcell = indi[0];
        ni = *(E_Int*)a2; E_Int nic = K_FUNC::E_max(1,ni-1);
        nj = *(E_Int*)a3; E_Int njc = K_FUNC::E_max(1,nj-1);
        //nk = *(E_Int*)a4; //E_Int nkc = K_FUNC::E_max(1,nk-1);
        E_Int nicnjc = nic*njc; ninj = ni*nj;
        k = indcell/nicnjc; j = (indcell-k*nicnjc)/nic; i = indcell-j*nic-k*nicnjc;
        E_Int ind0 = i+j*ni+k*ninj;
        E_Int ind1 = ind0+1;
        E_Int ind2 = ind0+ni;
        E_Int ind3 = ind2+1;
        E_Int ind4 = ind0+ninj;
        E_Int ind5 = ind1+ninj;
        E_Int ind6 = ind2+ninj;
        E_Int ind7 = ind3+ninj;
        val = cfp[0]*f0[ind0] + cfp[1]*f0[ind1] +
          cfp[2]*f0[ind2] + cfp[3]*f0[ind3] +
          cfp[4]*f0[ind4] + cfp[5]*f0[ind5] + 
          cfp[6]*f0[ind6] + cfp[7]*f0[ind7];              
      }
      else //cas non structure
      {
        noet = indi[0];
        FldArrayI& cn0 = *(FldArrayI*)a2;
        nvert = cn0.getNfld();
        for (E_Int nov = 1; nov <= nvert; nov++)
        {
          ind0 = cn0(noet,nov)-1;
          val += cfp[nov-1]*f0[ind0];          
        }
      }
      break;

    case  5://formule directionnelle (ex O5ABC) 
      if (a4 == NULL)// NS : indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compOneInterpolatedValue: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      //nk = *(E_Int*)a4;      
      ninj = ni*nj;
      ind0 = indi[0]; nocf = 0;
      k = ind0/(ninj); j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
      for (E_Int k0 = 0; k0 < 5; k0++)
        for (E_Int j0 = 0; j0 < 5; j0++)
          for (E_Int i0 = 0; i0 < 5; i0++)
          {
            E_Int ind = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
            val += cfp[i0]*cfp[j0+5]*cfp[k0+10]*f0[ind];  
          }
      break;
      
    default:
      printf("compOneInterpolatedValue: unknown interpolation type.");
      return -1;
  }
  return 1;
}
//=============================================================================
/* IN: ni, nj, nk: dimensions de f0 dans le cas structure
       -1, -1, -1 si le donneur est non structure
   IN: f0: champ du bloc donneur pour l'interpolation
   IN: indi: indices des pts de la molecule d interpolation
             peut etre defini par des sommets de la molecule donneuse
                              par l'indice de la cellule donneuse
             doit etre coherent avec f.
   IN: indiSize: taille de indi
   IN: cf: coefs d'interpolation associes 
   IN: ind: indice du pt a interpoler 
   IN: interpType: permet de determiner la formule appliquee
   OUT: f: champs interpoles au point ind 
   Retourne -1 si erreur 
*/
//=============================================================================
short K_INTERP::compInterpolatedValues(
  E_Int* indi, FldArrayF& cf,
  FldArrayF& f0,  void* a2, void* a3, void* a4, 
  E_Int ind, E_Int type, FldArrayF& f)
{
  E_Int nfld = f0.getNfld();
  E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7, indcell; 
  E_Int ni, nj, ninj, nvert, noet, i, j, k, nocf;
  E_Float* cfp = cf.begin();
  E_Float* fp0; E_Float* fp;

  switch (type) 
  {
    case 0: 
      nocf = indi[0]; // nb de pts pour la formule
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp = f.begin(eq);
        fp0 = f0.begin(eq);
        fp[ind] = 0.; 
        for (E_Int no = 1; no <= nocf; no++)
        {
          ind0 = indi[no];
          fp[ind] += cfp[no-1]*fp0[ind0];
        }
      }      
      break;
      
    case 1:
      ind0 = indi[0];
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp = f.begin(eq);
        fp0 = f0.begin(eq);
        fp[ind] = cfp[0] * fp0[ind0];
      }      
      break;

    case 2:// les indices et f sont localises de maniere identique, ex O2CF
      if (a4 == NULL)// NS : indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compOneInterpolatedValue: type 2 not valid for unstructured interpolations.\n");
        return -1;
      }
      else //STRUCT: indice du premier sommet en bas a gauche stockee      
      {
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        //nk = *(E_Int*)a4;      
        ninj = ni*nj;  
        ind0 = indi[0]; nocf = 0; 
        k = ind0/ninj;  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp = f.begin(eq);
          fp0 = f0.begin(eq);
          fp[ind] = 0.; nocf = 0;
          for (E_Int k0 = 0; k0 < 2; k0++)
            for (E_Int j0 = 0; j0 < 2; j0++)
              for (E_Int i0 = 0; i0 < 2; i0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                fp[ind] += cfp[nocf]*fp0[ind0];
                nocf += 1;
              }
        }
      }
      break;
  

    case 22:// les indices et f sont localises de maniere identique, ex O2CF
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      ind0 = indi[0]; nocf = 0;
      j = ind0/ni; i = ind0-j*ni;
     
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp = f.begin(eq);
        fp0 = f0.begin(eq);
        fp[ind] = 0.; nocf = 0;
        for (E_Int j0 = 0; j0 < 2; j0++)
          for (E_Int i0 = 0; i0 < 2; i0++)
          {
            ind0 = (i+i0)+(j+j0)*ni;
            fp[ind] += cfp[nocf]*fp0[ind0];
            nocf += 1;
          }
      }      
      break;

    case 3://formule directionnelle (ex O3ABC): on stocke les indices i,j,k des sommets 
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compInterpolatedValues: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      else // STRUCTURE: indice de la premiere cellule en bas a gauche stockee      
      {
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        //nk = *(E_Int*)a4;  
        ninj = ni*nj;  
        ind0 = indi[0]; nocf = 0;
        k = ind0/ninj;  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp = f.begin(eq);
          fp0 = f0.begin(eq);
          fp[ind] = 0.;
          for (E_Int i0 = 0; i0 < 3; i0++)
            for (E_Int j0 = 0; j0 < 3; j0++)
              for (E_Int k0 = 0; k0 < 3; k0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                fp[ind] += cfp[i0]*cfp[j0+3]*cfp[k0+6]* fp0[ind0];  
              }
        }
      }
      break;

   case 44://formule directionnelle (ex O4ABC): on stocke les indices i,j,k des sommets
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compInterpolatedValues: type 44 not valid for unstructured interpolations.\n");
        return -1;
      }
      else // STRUCTURE: indice de la premiere cellule en bas a gauche stockee
      {
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        ninj = ni*nj;
        ind0 = indi[0]; nocf = 0;
        k = ind0/ninj;  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp = f.begin(eq);
          fp0 = f0.begin(eq);
          fp[ind] = 0.;
          for (E_Int i0 = 0; i0 < 4; i0++)
            for (E_Int j0 = 0; j0 < 4; j0++)
              for (E_Int k0 = 0; k0 < 4; k0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                fp[ind] += cfp[i0]*cfp[j0+4]*cfp[k0+8]* fp0[ind0];
              }
        }
      }
      break;

    case  4:// on stocke le centre de la cellule et on interpole des sommets 
      if (a4 != NULL)// Structure
      {
        indcell = indi[0];
        ni = *(E_Int*)a2; E_Int nic = K_FUNC::E_max(1,ni-1);
        nj = *(E_Int*)a3; E_Int njc = K_FUNC::E_max(1,nj-1);
        //nk = *(E_Int*)a4; //E_Int nkc = K_FUNC::E_max(1,nk-1);
        E_Int nicnjc = nic*njc; ninj = ni*nj;
        k = indcell/nicnjc; j = (indcell-k*nicnjc)/nic; i = indcell-j*nic-k*nicnjc;
        ind0 = i+j*ni+k*ninj;
        ind1 = ind0+1;
        ind2 = ind0+ni;
        ind3 = ind2+1;
        ind4 = ind0+ninj;
        ind5 = ind1+ninj;
        ind6 = ind2+ninj;
        ind7 = ind3+ninj;
        for (E_Int eq = 1 ; eq <= nfld; eq++)
        {
          f(ind, eq) =  
            cfp[0]*f0(ind0,eq) + cfp[1]*f0(ind1,eq) + 
            cfp[2]*f0(ind2,eq) + cfp[3]*f0(ind3,eq) +
            cfp[4]*f0(ind4,eq) + cfp[5]*f0(ind5,eq) + 
            cfp[6]*f0(ind6,eq) + cfp[7]*f0(ind7,eq); 
        }       
      }
      else //cas non structure
      {
        noet = indi[0];
        FldArrayI& cn0 = *(FldArrayI*)a2;
        nvert = cn0.getNfld();
        for (E_Int eq = 1; eq <= nfld; eq++) 
        {
          fp = f.begin(eq);
          fp0 = f0.begin(eq);
          fp[ind] = 0.;
          for (E_Int nov = 1; nov <= nvert; nov++)
          {
            ind0 = cn0(noet,nov)-1;
            fp[ind] += cfp[nov-1]*fp0[ind0];
          }
        }
      }
      break;

    case  5://formule directionnelle (ex O5ABC) 
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d'interpolation
      {
        printf("Error: compInterpolatedValues: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      //nk = *(E_Int*)a4;      
      ninj = ni*nj;
      ind0 = indi[0]; nocf = 0;
      k = ind0/(ninj); j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);     
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp = f.begin(eq);
        fp0 = f0.begin(eq);
        fp[ind] = 0.;
        for (E_Int k0 = 0; k0 < 5; k0++)
          for (E_Int j0 = 0; j0 < 5; j0++)
            for (E_Int i0 = 0; i0 < 5; i0++)
            {
              ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
              fp[ind] += cfp[i0]*cfp[j0+5]*cfp[k0+10]* fp0[ind0];  
            }
      }      
      break;
      
    default:
      printf("compInterpolatedValues: unknown interpolation type.");
      return -1;
  }
  return 1;
}
/*---------------------------------------------------------------*/
/* Meme fonction que precedemment mais avec posvars              */
/*---------------------------------------------------------------*/
short K_INTERP::compInterpolatedValues(
  E_Int* indi, FldArrayF& cf,
  FldArrayF& f0,  void* a2, void* a3, void* a4, 
  E_Int ind, E_Int type, FldArrayF& f, vector<E_Int>& posvars)
{
  E_Int nfld = posvars.size();
  E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7, indcell; 
  E_Int ni, nj, ninj, nvert, noet, i, j, k, nocf;
  E_Float* cfp = cf.begin();
  E_Int posv;
  E_Float* fp; E_Float* fp0;

  switch (type) 
  {
    case 0: 
      nocf = indi[0];// nb de pts pour la formule
      for (E_Int noeq = 0; noeq < nfld; noeq++)
      {
        posv = posvars[noeq];
        fp = f.begin(noeq+1);
        fp0 = f0.begin(posv);
        fp[ind] = 0.; 
        for (E_Int no = 1; no <= nocf; no++)
        {
          ind0 = indi[no];
          fp[ind] += cfp[no-1]*fp0[ind0];
        }
      }      
      break;
      
    case 1:
      ind0 = indi[0];
      for (E_Int noeq = 0; noeq < nfld; noeq++)
      {
        posv = posvars[noeq];
        fp = f.begin(noeq+1);
        fp0 = f0.begin(posv);
        fp[ind] = cfp[0] * fp0[ind0];
      }      
      break;

    case 2:// les indices et f sont localises de maniere identique, ex O2CF
      if (a4 == NULL)// NS : indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compOneInterpolatedValue: type 2 not valid for unstructured interpolations.\n");
        return -1;
      }
      else //STRUCT: indice du premier sommet en bas a gauche stockee      
      {
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        //nk = *(E_Int*)a4;      
        ninj = ni*nj;  
        ind0 = indi[0]; nocf = 0; 
        k = ind0/ninj;  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int noeq = 0; noeq < nfld; noeq++)
        {
          posv = posvars[noeq];
          fp = f.begin(noeq+1);
          fp0 = f0.begin(posv);
          fp[ind] = 0.; nocf = 0;
          for (E_Int k0 = 0; k0 < 2; k0++)
            for (E_Int j0 = 0; j0 < 2; j0++)
              for (E_Int i0 = 0; i0 < 2; i0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                fp[ind] += cfp[nocf]*fp0[ind0];
                nocf += 1;
              }
        }
      }
      break;
  

    case 22:// les indices et f sont localises de maniere identique, ex O2CF
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      ind0 = indi[0]; nocf = 0;
      j = ind0/ni; i = ind0-j*ni;
     
      for (E_Int noeq = 0; noeq < nfld; noeq++)
      {
        posv = posvars[noeq];
        fp = f.begin(noeq+1);
        fp0 = f0.begin(posv);
        fp[ind] = 0.; nocf = 0;
        for (E_Int j0 = 0; j0 < 2; j0++)
          for (E_Int i0 = 0; i0 < 2; i0++)
          {
            ind0 = (i+i0)+(j+j0)*ni;
            fp[ind] += cfp[nocf]*fp0[ind0];
            nocf += 1;
          }
      }      
      break;

    case 3://formule directionnelle (ex O3ABC): on stocke les indices i,j,k des sommets 
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compInterpolatedValues: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      else // STRUCTURE: indice de la premiere cellule en bas a gauche stockee      
      {
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        //nk = *(E_Int*)a4;  
        ninj = ni*nj;  
        ind0 = indi[0]; nocf = 0;
        k = ind0/ninj;  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int noeq = 0; noeq < nfld; noeq++)
        {
          posv = posvars[noeq];
          fp = f.begin(noeq+1);
          fp0 = f0.begin(posv);       
          fp[ind] = 0.;
          for (E_Int i0 = 0; i0 < 3; i0++)
            for (E_Int j0 = 0; j0 < 3; j0++)
              for (E_Int k0 = 0; k0 < 3; k0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                fp[ind] += cfp[i0]*cfp[j0+3]*cfp[k0+6]* fp0[ind0];  
              }
        }
      }
      break;

    case  4:// on stocke le centre de la cellule et on interpole des sommets 
      if (a4 != NULL)// Structure
      {
        indcell = indi[0];
        ni = *(E_Int*)a2; E_Int nic = K_FUNC::E_max(1,ni-1);
        nj = *(E_Int*)a3; E_Int njc = K_FUNC::E_max(1,nj-1);
        //nk = *(E_Int*)a4; //E_Int nkc = K_FUNC::E_max(1,nk-1);
        E_Int nicnjc = nic*njc; ninj = ni*nj;
        k = indcell/nicnjc; j = (indcell-k*nicnjc)/nic; i = indcell-j*nic-k*nicnjc;
        ind0 = i+j*ni+k*ninj;
        ind1 = ind0+1;
        ind2 = ind0+ni;
        ind3 = ind2+1;
        ind4 = ind0+ninj;
        ind5 = ind1+ninj;
        ind6 = ind2+ninj;
        ind7 = ind3+ninj;

        for (E_Int noeq = 0; noeq < nfld; noeq++)
        {
          posv = posvars[noeq];
          fp = f.begin(noeq+1);
          fp0 = f0.begin(posv);          
          fp[ind] = 
          cfp[0]*fp0[ind0] + cfp[1]*fp0[ind1] + 
          cfp[2]*fp0[ind2] + cfp[3]*fp0[ind3] +
          cfp[4]*fp0[ind4] + cfp[5]*fp0[ind5] + 
          cfp[6]*fp0[ind6] + cfp[7]*fp0[ind7]; 
        }       
      }
      else //cas non structure
      {
        noet = indi[0];
        FldArrayI& cn0 = *(FldArrayI*)a2;
        nvert = cn0.getNfld();
        for (E_Int noeq = 0; noeq < nfld; noeq++)
        {
          posv = posvars[noeq];
          fp = f.begin(noeq+1);
          fp0 = f0.begin(posv);           
          fp[ind] = 0.;
          for (E_Int nov = 1; nov <= nvert; nov++)
          {
            ind0 = cn0(noet,nov)-1;
            fp[ind] += cfp[nov-1]*fp0[ind0];
          }
        }
      }
      break;

    case  5://formule directionnelle (ex O5ABC) 
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d'interpolation
      {
        printf("Error: compInterpolatedValues: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      //STRUCT : indice de la premiere cellule en bas a gauche stockee      
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      //nk = *(E_Int*)a4;      
      ninj = ni*nj;
      ind0 = indi[0]; nocf = 0;
      k = ind0/(ninj); j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);     
      for (E_Int noeq = 0; noeq < nfld; noeq++)
      {
        posv = posvars[noeq];
        fp = f.begin(noeq+1);
        fp0 = f0.begin(posv); 
        fp[ind] = 0.;
        for (E_Int k0 = 0; k0 < 5; k0++)
          for (E_Int j0 = 0; j0 < 5; j0++)
            for (E_Int i0 = 0; i0 < 5; i0++)
            {
              ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
              fp[ind] += cfp[i0]*cfp[j0+5]*cfp[k0+10]* fp0[ind0];  
            }
      }      
      break;
      
    default:
      printf("compInterpolatedValues: unknown interpolation type.");
      return -1;
  }
  return 1;
}
//=============================================================================
/* IN: ni, nj, nk: dimensions de f0 dans le cas structure
       -1, -1, -1 si le donneur est non structure
   IN: f0: champ du bloc donneur pour l'interpolation
   IN: indi: indices des pts de la molecule d interpolation
             peut etre defini par des sommets de la molecule donneuse
                              par l'indice de la cellule donneuse
             doit etre coherent avec f.
   IN: indiSize: taille de indi
   IN: cf: coefs d'interpolation associes 
   IN: f: tableau monodimensionnel des champs a interpoler
   IN: interpType: permet de determiner la formule appliquee
   OUT: f: champs interpoles au point ind 
   Retourne -1 si erreur 
*/
//=============================================================================
short K_INTERP::compInterpolatedField(
  E_Int* indi, FldArrayF& cf,
  FldArrayF& f0,  void* a2, void* a3, void* a4, 
  E_Int type, FldArrayF& f)
{
  E_Int nfld = f0.getNfld();
  E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7, indcell; 
  E_Int ni, nj, ninj, nvert, noet, i, j, k, nocf;
  E_Float* cfp = cf.begin();
  E_Float* fp0;

  switch (type) 
  {  
    case 0:
      nocf = indi[0];// nb de pts pour la formule
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp0 = f0.begin(eq);
        f[eq-1] = 0.; 
        for (E_Int no = 1; no <= nocf; no++)
        {
          ind0 = indi[no];
          f[eq-1] += cfp[no-1]*fp0[ind0];
        }
      }      
      break;

    case 1:
      ind0 = indi[0];
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp0 = f0.begin(eq);
        f[eq-1] = cfp[0] * fp0[ind0];
      }      
      break;

    case  2:// les indices et f sont localises de maniere identique, ex O2CF
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d'interpolation
      {
        printf("Error: compInterpolatedField: type 2 not valid for unstructured interpolations.\n");
        return -1;
      }
            
      else //STRUCT: indice de la premiere cellule en bas a gauche stockee
      {
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        //nk = *(E_Int*)a4;      
        ninj = ni*nj;  
        ind0 = indi[0]; nocf = 0;
        k = ind0/(ninj);  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp0 = f0.begin(eq);
          f[eq-1] = 0.; nocf = 0;
          for (E_Int k0 = 0; k0 < 2; k0++)
            for (E_Int j0 = 0; j0 < 2; j0++)
              for (E_Int i0 = 0; i0 < 2; i0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                f[eq-1] += cfp[nocf]*fp0[ind0];
                nocf += 1;
              }
        }
      }
      break;

      case 22:// les indices et f sont localises de maniere identique, ex O2CF
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        ind0 = indi[0]; nocf = 0;
        j = ind0/ni; i = ind0-j*ni;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp0 = f0.begin(eq);
          f[eq-1] = 0.; nocf = 0;
          for (E_Int j0 = 0; j0 < 2; j0++)
            for (E_Int i0 = 0; i0 < 2; i0++)
            {
              ind0 = (i+i0) + (j+j0)*ni;
              f[eq-1] += cfp[nocf]*fp0[ind0];
              nocf += 1;
            }
        }

      break;

    case  3://formule directionnelle (ex O3ABC): on stocke les indices i,j,k des sommets 
      if (a4 == NULL)// NS : indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compInterpolatedValues: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      ni = *(E_Int*)a2;
      nj = *(E_Int*)a3;
      //nk = *(E_Int*)a4;    
      ninj = ni*nj;  
      ind0 = indi[0]; nocf = 0;
      k = ind0/(ninj);  j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);

      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        fp0 = f0.begin(eq);
        f[eq-1] = 0.;nocf = 0;
          for (E_Int k0 = 0; k0 < 3; k0++)
            for (E_Int j0 = 0; j0 < 3; j0++)
              for (E_Int i0 = 0; i0 < 3; i0++)
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                f[eq-1] += cfp[i0]*cfp[j0+3]*cfp[k0+6]* fp0[ind0];  
                nocf+=1;
              }
      }
      break;

    case 4:// on stocke le centre de la cellule et on interpole des sommets 
      if (a4 != NULL)// Structure
      {
        indcell = indi[0];
        ni = *(E_Int*)a2; E_Int nic = K_FUNC::E_max(1,ni-1);
        nj = *(E_Int*)a3; E_Int njc = K_FUNC::E_max(1,nj-1);
        //nk = *(E_Int*)a4; //E_Int nkc = K_FUNC::E_max(1,nk-1);
        E_Int nicnjc = nic*njc; ninj = ni*nj;
        k = indcell/nicnjc; j = (indcell-k*nicnjc)/nic; i = indcell-j*nic-k*nicnjc;    
        ind0 = i+j*ni+k*ninj;
        ind1 = ind0+1;
        ind2 = ind0+ni;
        ind3 = ind2+1;
        ind4 = ind0+ninj;
        ind5 = ind1+ninj;
        ind6 = ind2+ninj;
        ind7 = ind3+ninj;
        for (E_Int eq = 1 ; eq <= nfld; eq++)
        {
          f[eq-1] =  
            cfp[0]*f0(ind0,eq) + cfp[1]*f0(ind1,eq) + 
            cfp[2]*f0(ind2,eq) + cfp[3]*f0(ind3,eq) +
            cfp[4]*f0(ind4,eq) + cfp[5]*f0(ind5,eq) + 
            cfp[6]*f0(ind6,eq) + cfp[7]*f0(ind7,eq);       
        }       
      }
      else //cas non structure
      {
        noet = indi[0];
        FldArrayI& cn0 = *(FldArrayI*)a2;
        nvert = cn0.getNfld();
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp0 = f0.begin(eq);
          f[eq-1] = 0.;
          for (E_Int nov = 1; nov <= nvert; nov++)
          {
            ind0 = cn0(noet,nov)-1;
            f[eq-1] += cfp[nov-1]*fp0[ind0];
          }
        }
      }
      break;

    case  5://formule directionnelle (ex O5ABC) 
      if (a4 == NULL)// NS: indices des nvert sommets de la molecule d interpolation
      {
        printf("Error: compInterpolatedField: type 3 not valid for unstructured interpolations.\n");
        return -1;
      }
      else
      {
        //STRUCT : indice de la premiere cellule en bas a gauche stockee      
        ni = *(E_Int*)a2;
        nj = *(E_Int*)a3;
        //nk = *(E_Int*)a4;      
        ninj = ni*nj;
        ind0 = indi[0]; nocf = 0;
        k = ind0/(ninj); j = (ind0-k*ninj)/ni; i = (ind0-j*ni-k*ninj);
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fp0 = f0.begin(eq);
          f[eq-1] = 0.;
          for (E_Int k0 = 0; k0 < 5; k0++)
            for (E_Int j0 = 0; j0 < 5; j0++)
              for (E_Int i0 = 0; i0 < 5; i0++)               
              {
                ind0 = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                f[eq-1] += cfp[i0]*cfp[j0+5]*cfp[k0+10] * fp0[ind0];  
              }
        }
      }
      break;
      
    default:
      printf("compInterpolatedField: unknown interpolation type: " SF_D_ ".\n", type);
      return -1;
  }
  return 1;
}
