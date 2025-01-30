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

// Global functions of K_KINTERP
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlkInterp.h"
# include "Fld/FldArray.h"
# include "String/kstring.h"
# include "Connect/connect.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;
using namespace K_ARRAY;

extern "C"
{
  void k6compvolofstructcell_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                              const E_Int& indcell, const E_Int& indnode,
                              const E_Float* x, const E_Float* y, const E_Float* z, 
                              E_Float& vol);
  void k6compvoloftetracell_(const E_Int& npts, 
                             const E_Int& ind1, const E_Int& ind2, 
                             const E_Int& ind3, const E_Int& ind4,
                             const E_Float* xt, const E_Float* yt, 
                             const E_Float* zt, E_Float& vol);
}

//=============================================================================
// Global bounding box of KMeshes
//=============================================================================
void K_KINTERP::boundingBox(vector<K_KINTERP::KMesh*>& KMeshList,
                            E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                            E_Float& xmin, E_Float& ymin, E_Float& zmin)
{
  E_Float xmaxl, ymaxl, zmaxl;
  E_Float xminl, yminl, zminl;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;
  xmin = +K_CONST::E_MAX_FLOAT;
  ymin = +K_CONST::E_MAX_FLOAT;
  zmin = +K_CONST::E_MAX_FLOAT;

  E_Int s = KMeshList.size();

  for (E_Int i = 0; i < s; i++)
  {
    KMeshList[i]->boundingBox(xmaxl, ymaxl, zmaxl,
                              xminl, yminl, zminl);
    xmax = E_max(xmax, xmaxl);
    ymax = E_max(ymax, ymaxl);
    zmax = E_max(zmax, zmaxl);
    xmin = E_min(xmin, xminl);
    ymin = E_min(ymin, yminl);
    zmin = E_min(zmin, zminl);
  }
}
//=============================================================================
// Interpolation cell of a list of blocks
// Retourne le no du bloc d'interpolation (commence a 1)
// 0 si rien n'est trouve.
// retourne le no du bloc ayant la plus petite cellule d interpolation
// attention si mixte (structure/non struct) indi et cf doivent etre alloues
// au plus grand
// si interp sur maillage structure: retourne le numero du bloc structure ds 
// la liste des blocs structures
// de meme pour les blocs non structures
//
// ATTENTION: cas non structure present: ne pas copier indiu,cfu dans indi,cf
// car indi et cf peuvent etre dimensionnes a partir du structure !
//
//=============================================================================
short K_KINTERP::getInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatas,
  vector<FldArrayF*>& fields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_KINTERP::BlkInterpData*>& interpDatau,
  vector<FldArrayF*>& fieldu, vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  E_Float& voli, FldArrayI& indi, FldArrayF& cf,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty)
{
  E_Int noblks = 0;
  E_Int noblku = 0;
  FldArrayI indis(indi.getSize());
  FldArrayF cfs(cf.getSize());
  FldArrayI indiu(indi.getSize());
  FldArrayF cfu(cf.getSize());
  E_Int shift = interpDatas.size();
  E_Float vols = K_CONST::E_MAX_FLOAT;
  E_Float volu = K_CONST::E_MAX_FLOAT;

  // Recherche de la meilleure cellule d interpolation structuree 
  if (interpDatas.size() > 0)
  {
    vols = selectBestStructuredInterpolationCell(
      x, y, z, interpDatas, fields, nis, njs, nks, 
      posxs, posys, poszs, poscs, 
      indis, cfs, noblks, 
      interpMeshType, interpType, nature, penalty);
  }
  // Recherche de la meilleure cellule d'interpolation non structuree
  if (interpDatau.size() > 0)
  {
    volu = selectBestUnstructuredInterpolationCell(
      x, y, z, interpDatau, fieldu, connectu, 
      posxu, posyu, poszu, poscu, 
      indiu, cfu, noblku,
      nature, penalty);
  }
  if (noblks > 0 && noblku > 0) //interpolation sur maillage mixte
  {
    if ( vols < volu ) 
    {
      for (E_Int i = 0; i < indis.getSize(); i++) indi[i] = indis[i];
      for (E_Int c = 0 ; c < cfs.getSize(); c++) cf[c] = cfs[c];
      voli = vols;
      return noblks;
    }
    else 
    {  
      for (E_Int i = 0; i < indiu.getSize(); i++) indi[i] = indiu[i];
      for (E_Int c = 0 ; c < cfu.getSize(); c++) cf[c] = cfu[c];
      voli = volu;
      return noblku+shift;
    }
  }
  else if (noblks > 0 && noblku == 0) // structure
  {   
    for (E_Int i = 0; i < indis.getSize(); i++) indi[i] = indis[i];
    for (E_Int c = 0 ; c < cfs.getSize(); c++) cf[c] = cfs[c];
    voli = vols;
    return noblks;
  }
  else if (noblks == 0 && noblku > 0) // non structure 
  {
    for (E_Int i = 0; i < indiu.getSize(); i++) indi[i] = indiu[i];
    for (E_Int c = 0 ; c < cfu.getSize(); c++) cf[c] = cfu[c];
    voli = volu;
    return noblku+shift;
  }
  else 
  {
    voli = K_CONST::E_MAX_FLOAT;
    return 0;
  }
}

//=============================================================================
// Interpolation cell of a list of blocks
// Retourne le no du bloc d'interpolation (commence a 1)
// 0 si rien n'est trouve.
// retourne le no du bloc ayant la plus petite cellule d interpolation
// attention si mixte (structure/non struct) indi et cf doivent etre alloues
// au plus grand
// si interp sur maillage structure : retourne le numero du bloc structure ds 
// la liste des blocs structures
// de meme pour les blocs non structures
//
// ATTENTION: cas non structure present : ne pas copier indiu,cfu dans indi,cf
// car indi et cf peuvent etre dimensionnes a partir du structure !
//
//=============================================================================
short K_KINTERP::getInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatas,
  vector<FldArrayF*>& fields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_KINTERP::BlkInterpData*>& interpDatau,
  vector<FldArrayF*>& fieldu, vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  FldArrayI& indi, FldArrayF& cf,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData::InterpolationType interpType,
  E_Int nature, E_Int penalty)
{
  E_Int noblks = 0;
  E_Int noblku = 0;
  FldArrayI indis(indi.getSize());
  FldArrayF cfs(cf.getSize());
  FldArrayI indiu(indi.getSize());
  FldArrayF cfu(cf.getSize());
  
  E_Float vols = K_CONST::E_MAX_FLOAT;
  E_Float volu = K_CONST::E_MAX_FLOAT;

  E_Int shift = interpDatas.size();
  // Recherche de la meilleure cellule d interpolation structuree 
  if ( interpDatas.size() > 0 )
  {
    vols = selectBestStructuredInterpolationCell(
      x, y, z, interpDatas, fields, nis, njs, nks, 
      posxs, posys, poszs, poscs, 
      indis, cfs, noblks, 
      interpMeshType, interpType, nature, penalty);
  }
  // Recherche de la meilleure cellule d interpolation non structuree
  if ( interpDatau.size() > 0 )
  {
    volu = selectBestUnstructuredInterpolationCell(
      x, y, z, interpDatau, fieldu, connectu, 
      posxu, posyu, poszu, poscu, 
      indiu, cfu, noblku, nature, penalty);
  }

  if ( noblks > 0 && noblku > 0 ) //interpolation sur maillage mixte
  {
    if ( vols < volu ) 
    {
      for (E_Int i = 0; i < indis.getSize(); i++) indi[i] = indis[i];
      for (E_Int c = 0 ; c < cfs.getSize(); c++) cf[c] = cfs[c];
      return noblks;
    }
    else 
    {
      for (E_Int i = 0; i < indiu.getSize(); i++) indi[i] = indiu[i];
      for (E_Int c = 0 ; c < cfu.getSize(); c++) cf[c] = cfu[c];
      return noblku+shift;
    }
  }
  else if ( noblks > 0 && noblku == 0 ) // structure
  {
    for (E_Int i = 0; i < indis.getSize(); i++) indi[i] = indis[i];
    for (E_Int c = 0 ; c < cfs.getSize(); c++) cf[c] = cfs[c];
    return noblks;
  }
  else if ( noblks == 0 && noblku > 0 ) // non structure 
  {
    for (E_Int i = 0; i < indiu.getSize(); i++) indi[i] = indiu[i];
    for (E_Int c = 0 ; c < cfu.getSize(); c++) cf[c] = cfu[c];
    return noblku+shift;
  }
  else return 0;
}

//=============================================================================
/* Recherche de cellule d'interpolation sur maillages structures */
//=============================================================================
short K_KINTERP::getInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatas,
  vector<FldArrayF*>& fields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  FldArrayI& indis, FldArrayF& cfs,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData::InterpolationType interpType, 
  E_Int nature, E_Int penalty)
{
  E_Int noblks;

  // Recherche de la meilleure cellule d'interpolation structuree  
  selectBestStructuredInterpolationCell(
    x, y, z, interpDatas, fields, 
    nis, njs, nks, 
    posxs, posys, poszs, poscs, 
    indis, cfs, noblks, 
    interpMeshType, interpType, nature, penalty);
  return noblks;
}
//=============================================================================
/* recherche de la cellule d'interpolation sur maillages non structures */
//=============================================================================
short K_KINTERP::getInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatau,
  vector<FldArrayF*>& fieldu, vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  FldArrayI& indiu, FldArrayF& cfu, 
  E_Int nature, E_Int penalty)
{
  E_Int noblku;

  // Recherche de la meilleure cellule d'interpolation structuree  
  selectBestUnstructuredInterpolationCell(
    x, y, z, interpDatau, fieldu, connectu,  
    posxu, posyu, poszu, poscu, indiu, cfu, noblku, 
    nature, penalty); 

  return noblku;
}

//=============================================================================
/* recherche de la cellule d'extrapolation */
//=============================================================================
short K_KINTERP::getExtrapolationCell(
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatas,
  vector<FldArrayF*>& fields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_KINTERP::BlkInterpData*>& interpDatau,
  vector<char*> eltType,
  vector<FldArrayF*>& fieldu, vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  E_Float& voli, FldArrayI& indi, FldArrayF& cf, E_Float cfMax,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Int noblks = 0;
  E_Int noblku = 0;
  FldArrayI indis(indi.getSize());
  FldArrayF cfs(cf.getSize());
  FldArrayI indiu(indi.getSize());
  FldArrayF cfu(cf.getSize());
  E_Int shift = interpDatas.size();
  E_Float vols = K_CONST::E_MAX_FLOAT;
  E_Float volu = K_CONST::E_MAX_FLOAT;
  // Recherche de la meilleure cellule d'extrapolation structuree 
  if (interpDatas.size() > 0)
  {
    vols = selectBestStructuredExtrapolationCell(
      x, y, z, interpDatas, fields, nis, njs, nks, 
      posxs, posys, poszs, poscs, 
      indis, cfs, noblks, cfMax,
      interpMeshType, interpType);
  }
  // Recherche de la meilleure cellule d'extrapolation non structuree
  if (interpDatau.size() > 0)
  {
    volu = selectBestUnstructuredExtrapolationCell(
      x, y, z, interpDatau, eltType, fieldu, connectu, posxu, posyu, poszu, poscu, 
      indiu, cfu, noblku);
  }
  if ( noblks > 0 && noblku > 0 ) //interpolation sur maillage mixte
  {
    if ( vols < volu ) 
    {
      for (E_Int i = 0; i < indis.getSize(); i++)
        indi[i] = indis[i];
      for (E_Int c = 0 ; c < cfs.getSize(); c++)
        cf[c] = cfs[c];
      voli = vols;
      return noblks;
    }
    else 
    {  
      for (E_Int i = 0; i < indiu.getSize(); i++)
        indi[i] = indiu[i];
      for (E_Int c = 0 ; c < cfu.getSize(); c++)
        cf[c] = cfu[c];
      
      voli = volu;
      return noblku+shift;
    }
  }
  else if ( noblks > 0 && noblku == 0 ) // structure
  {   
    for (E_Int i = 0; i < indis.getSize(); i++)
      indi[i] = indis[i];
    for (E_Int c = 0 ; c < cfs.getSize(); c++)
      cf[c] = cfs[c];
    voli = vols;
    return noblks;
  }
  else if ( noblks == 0 && noblku > 0 ) // non structure 
  {
    for (E_Int i = 0; i < indiu.getSize(); i++)
      indi[i] = indiu[i];
    for (E_Int c = 0 ; c < cfu.getSize(); c++)
      cf[c] = cfu[c];

    voli = volu;
    return noblku+shift;
  }
  else 
  {
    voli = K_CONST::E_MAX_FLOAT;
    return 0;
  }
}

//=============================================================================
/* recherche de la cellule d'interpolation la plus petite sur les maillages
   structures : retourne le volume de cette cellule */
//=============================================================================
E_Float K_KINTERP::selectBestStructuredInterpolationCell( 
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatas,
  vector<FldArrayF*>& fields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posx, vector<E_Int>& posy, 
  vector<E_Int>& posz, vector<E_Int>& posc,
  FldArrayI& indi, FldArrayF& cf, E_Int& noblks, 
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData::InterpolationType interpType, 
  E_Int nature, E_Int penalty)
{
  const E_Float geomCutOff = K_CONST::E_GEOM_CUTOFF;
  E_Float vol;
  // numero du bloc d'interpolation
  // initialise a 0 si pas trouve
  short noblk = 0; 
  E_Int s = interpDatas.size();
  short found;
  E_Float best = K_CONST::E_MAX_FLOAT; // for now : volume

  FldArrayI tmpIndi(indi.getSize());
  FldArrayF tmpCf(cf.getSize());
  E_Int ind, ni, nj, nk, ninj;
  E_Int posx0, posy0, posz0, posc0;
  E_Float cellN, val;
  E_Int isvalid, extrapB;
  for (E_Int i = 0; i < s; i++)
  {
    posx0 = posx[i]; posy0 = posy[i]; posz0 = posz[i]; posc0 = posc[i];
    ni = nis[i]; nj = njs[i]; nk = nks[i]; ninj = ni*nj;

    FldArrayF& oneField = *fields[i];
    E_Float* cellNp = NULL;
    if ( posc0 > 0 ) cellNp = oneField.begin(posc0);
    extrapB = 0;

    // Recherche de la cellule d'interpolation
    found = interpDatas[i]->getInterpolationCellStruct(x, y, z,
                                                       tmpIndi, 
                                                       tmpCf,
                                                       interpMeshType,
                                                       interpType);

    if (found > 0)
    {      
      // determination du cellN de la cellule d'interpolation
      if ( nature == 0 ) val = 1.;
      else val = 0.;
      E_Int order = tmpIndi[0];
      // Verif si on est au bord dans le cas penalty =1
      if ( penalty == 1 )
      {
        if ( tmpIndi[1] == 0 || tmpIndi[1] == ni-2) extrapB = 1; 
        else if  ( tmpIndi[order+1] == 0 || tmpIndi[order+1] == nj-2) extrapB = 1; 
        else if  ( tmpIndi[2*order+1] == 0 || tmpIndi[2*order+1] == nk-2) extrapB = 1; 
      }
      if ( posc0 > 0 )  // cellN existe dans oneField
      {
        for (E_Int i0 = 1; i0 <= order; i0++)
          for (E_Int j0 = order+1; j0 <= 2*order; j0++)
            for (E_Int k0 = 2*order+1; k0 <= 3*order; k0++)
            {
              ind = tmpIndi[i0] + tmpIndi[j0] * ni + tmpIndi[k0] * ninj;
              cellN = cellNp[ind];
              if ( nature == 0 )  val *= cellN;// traitement simple : pas de pt masque ds la molecule donneuse
              else val += cellN*(cellN-2.);// nature = 1 pas de pt masque ou interpole dans la cellule donneuse
                
            }
      }
      isvalid = 0;
      if ( nature == 1 && K_FUNC::fEqual(val,K_CONST::ONE,geomCutOff) == true) isvalid = 1;
      else if (nature == 0 && K_FUNC::fEqualZero(val,geomCutOff) == false) isvalid = 1;
      
      if (isvalid == 1) // pas de pt masque ou interpole dans la cellule donneuse
      {
        E_Int indnode = tmpIndi[1]+tmpIndi[order+1]*ni+tmpIndi[order*2+1]*ninj;
        // calcul de cellvol
        k6compvolofstructcell_(ni, nj, nk, -1, indnode, oneField.begin(posx0),
                               oneField.begin(posy0), oneField.begin(posz0),vol);
        if (penalty == 1 && extrapB == 1 ) vol += 1.e3;
        if ( vol < best ) 
        {
          best = vol;
          indi = tmpIndi;
          cf = tmpCf;
          noblk = i+1;
        }
      }
    }
  }
  noblks = noblk;
  return best;
}
//=============================================================================
/* recherche de la cellule d'interpolation la plus petite sur les maillages
   non structures: retourne le volume de cette cellule */
//=============================================================================
E_Float K_KINTERP::selectBestUnstructuredInterpolationCell( 
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatau,
  vector<FldArrayF*>& fieldu, vector<FldArrayI*>& connectu,
  vector<E_Int>& posx, vector<E_Int>& posy, 
  vector<E_Int>& posz, vector<E_Int>& posc,
  FldArrayI& indi, FldArrayF& cf, E_Int& noblku, 
  E_Int nature, E_Int penalty)
{
  const E_Float geomCutOff = K_CONST::E_GEOM_CUTOFF;
  E_Float vol = K_CONST::E_MAX_FLOAT;
  short noblk = 0;   // numero du bloc d interpolation
  E_Int s = interpDatau.size();
  short found;

  E_Float best = K_CONST::E_MAX_FLOAT; // for now : volume
  //E_Int nvoisinsTetra = 4; // un tetra a 4 voisins, dans le cas contraire suppose sur un bord du maillage

  FldArrayI tmpIndi(indi.getSize());
  FldArrayI indisav(indi.getSize());
  FldArrayF tmpCf(cf.getSize()); tmpCf.setAllValuesAtNull();
  FldArrayF cfsav(cf.getSize()); cfsav.setAllValuesAtNull();
  E_Int posx0, posy0, posz0, posc0;
  E_Float cellN, val;
  E_Int noet;
  E_Int isvalid, noei; 
  E_Int extrapB = 0;
  
  for (E_Int i = 0; i < s; i++)
  {
    posx0 = posx[i]; posy0 = posy[i]; posz0 = posz[i]; posc0 = posc[i];
    FldArrayF& oneField = *fieldu[i];
    FldArrayI& cnloc = *connectu[i];
    E_Int npts = oneField.getSize();
    E_Int nelts = cnloc.getNfld();
    vector< vector<E_Int> > cEEN(nelts);//elts voisins
    if (penalty == 1)
      K_CONNECT::connectEV2EENbrs("TETRA", fieldu[i]->getSize(), *connectu[i], cEEN);
    else cEEN.clear();

    // Recherche de la cellule d'interpolation
    found = interpDatau[i]->getInterpolationCellUnstruct(x, y, z, noet, tmpIndi, tmpCf);
    extrapB = 0;
    if (penalty == 1) 
    {if (cEEN[noet].size() != 4) extrapB = 1;}
    if (found > 0)
    {      
      noei = tmpIndi[1];
      // determination du cellN de la cellule d interpolation
      if (nature == 0) val = 1.;
      else val = 0.;
      if (posc0 > 0)  // cellN existe dans oneField
      {
        E_Float* cellNp = oneField.begin(posc0);
        for (E_Int nov = 1; nov <= 4; nov++)
        {
          E_Int indv = cnloc(noei,nov)-1;
          cellN = cellNp[indv];        
          if (nature == 1) val += cellN*(2.-cellN);
          else val *= cellN;
        }
      }
      
      isvalid = 0;
      if (nature == 1 && K_FUNC::fEqual(val,K_CONST::ONE,geomCutOff) == true) isvalid = 1;
      else if (nature == 0 && K_FUNC::fEqualZero(val,geomCutOff) == false) isvalid = 1;
      if (isvalid == 1)
      {        
        k6compvoloftetracell_(npts, cnloc(noei,1)-1, cnloc(noei,2)-1, cnloc(noei,3)-1, cnloc(noei,4)-1, 
                              oneField.begin(posx0), oneField.begin(posy0),oneField.begin(posz0), vol);
        if (extrapB == 1 ) vol += 1.e3;
        if ( vol < best ) 
        {
          best = vol;
          indisav = tmpIndi;
          cfsav = tmpCf;
          noblk = i+1;
        }
      }
    }
  }

  noblku = noblk;
  indi = indisav;
  cf = cfsav;
  return best;
}
//=============================================================================
/* recherche de la cellule d'extrapolation la plus petite sur les maillages
   structures: retourne le volume de cette cellule */
//=============================================================================
E_Float K_KINTERP::selectBestStructuredExtrapolationCell( 
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatas,
  vector<FldArrayF*>& fields,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posx, vector<E_Int>& posy, 
  vector<E_Int>& posz, vector<E_Int>& posc,
  FldArrayI& indi, 
  FldArrayF& cf, E_Int& noblks, E_Float cfMax,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Float vol;

  // numero du bloc d'interpolation
  // initialise a 0 si pas trouve
  short noblk = 0; 
  E_Int s = interpDatas.size();
  short found;
  E_Float best = K_CONST::E_MAX_FLOAT; // for now : volume

  FldArrayI tmpIndi(indi.getSize());
  FldArrayF tmpCf(cf.getSize());
  E_Int ni, nj, nk, ninj;
  E_Int posx0, posy0, posz0, posc0;
  E_Int order, indnode;

  for (E_Int i = 0; i < s; i++)
  {
    FldArrayF& oneField = *fields[i];
    posx0 = posx[i]; posy0 = posy[i]; posz0 = posz[i];
    posc0 = posc[i];
    ni = nis[i]; nj = njs[i]; nk = nks[i];
    ninj = ni*nj;
    
    // Recuperation du cellNatureField
    E_Int cellnsize = oneField.getSize();
    FldArrayF cellNatureField(cellnsize);
    if ( posc0 > 0 )  // cellN existe dans oneField
      for (E_Int i=0; i < cellnsize; i++)
        cellNatureField[i] = oneField(i,posc0);
    else
      cellNatureField.setAllValuesAt(1.);

    // Recherche de la cellule d'interpolation
    found = interpDatas[i]->getExtrapolationCellStruct(x, y, z,
                                                       tmpIndi, 
                                                       tmpCf,
                                                       0, cfMax,
                                                       cellNatureField,
                                                       interpType);
    if (found > 0)
    {
      order = tmpIndi[0];
      indnode = tmpIndi[1] + tmpIndi[order+1]*ni + tmpIndi[order*2+1]*ninj;
      // calcul de cellvol
      k6compvolofstructcell_(ni, nj, nk, -1, indnode, 
                             oneField.begin(posx0),
                             oneField.begin(posy0),
                             oneField.begin(posz0), vol);
      
      if (vol < best) 
      {
        best = vol;
        indi = tmpIndi;
        cf = tmpCf;
        noblk = i+1;
      }
    }
  }
  noblks = noblk;
  return best;
}
//=============================================================================
/* recherche de la cellule d'extrapolation la plus petite sur les maillages
   non structures: retourne le volume de cette cellule */
//=============================================================================
E_Float K_KINTERP::selectBestUnstructuredExtrapolationCell( 
  E_Float x, E_Float y, E_Float z,
  vector<K_KINTERP::BlkInterpData*>& interpDatau, vector<char*> eltType,
  vector<FldArrayF*>& fieldu, vector<FldArrayI*>& connectu,
  vector<E_Int>& posx, vector<E_Int>& posy, 
  vector<E_Int>& posz, vector<E_Int>& posc,
  FldArrayI& indi, FldArrayF& cf, E_Int& noblku)
{
  E_Float vol = K_CONST::E_MAX_FLOAT;
  short noblk = 0;   // numero du bloc d interpolation
  E_Int s = interpDatau.size();
  short found;
  E_Int noet, noei;

  E_Float best = K_CONST::E_MAX_FLOAT; // for now : volume
  
  FldArrayI tmpIndi(2);
  FldArrayF tmpCf(4);
  FldArrayI indisav(2);
  FldArrayF cfsav(4);
  cfsav.setAllValuesAtNull();
  cf.setAllValuesAtNull();
  E_Int posx0, posy0, posz0, posc0;

  for (E_Int i = 0; i < s; i++)
  {
    posx0 = posx[i]; posy0 = posy[i]; posz0 = posz[i]; posc0 = posc[i];
    FldArrayF& oneField = *fieldu[i];
    FldArrayI& cnloc = *connectu[i];
    E_Int npts = oneField.getSize();
    //E_Int nelts = cnloc.getNfld();

    // Recuperation du cellNatureField
    E_Int cellnsize = oneField.getSize();
    FldArrayF cellNatureField(cellnsize);
    if ( posc0 > 0 )  // cellN existe dans oneField
      for (E_Int i=0; i < cellnsize; i++)
        cellNatureField[i] = oneField(i,posc0);
    else cellNatureField.setAllValuesAt(1.);

    // Recherche de la cellule d extrapolation
    found = interpDatau[i]->getExtrapolationCellUnstr(x, y, z, noet,
                                                      tmpCf, tmpIndi, 0, cellNatureField);
    if (found > 0)
    {      
      noei = tmpIndi[1];
      k6compvoloftetracell_(npts, cnloc(noei,1)-1, cnloc(noei,2)-1, cnloc(noei,3)-1, cnloc(noei,4)-1, 
                            oneField.begin(posx0), oneField.begin(posy0),oneField.begin(posz0), vol);
      if ( vol < best ) 
      {
        best = vol;
        indisav = tmpIndi;
        cfsav = tmpCf;
        noblk = i+1;
      }
    }
  }

  noblku = noblk;
  indi = indisav;
  cf = cfsav;
  return best;
}
//===========================================================================
/* Creation de la liste des interpData pour des maillages non structures
   Skip autres que tetras. Skip les sans coords.
   IN: fieldsIn: liste des champs contenant les infos sur les grilles
   IN: varStringIn: liste des varstring correspondantes
   IN: cnIn: connectivite 
   IN: eltsIn: liste des types d elements 
   OUT: fieldsOut: liste des champs contenant les infos sur les grilles
   OUT: varStringOut: liste des varstring correspondantes
   IN: cnOut: connectivite 
   IN: eltTypesOut: liste des types d elements 
   OUT: posxt, posyt, poszt, posct: positions de x,y,z,celln (eventuellt)
   OUT: listOfUnstrMeshes: c est pour pouvoir les detruire 
   OUT: interpDatas: liste des interpDatas crees basees sur listOfMeshes
   posx, posy, posz, posc demarrent a 1 */
//===========================================================================
void K_KINTERP::
buildListOfUnstrInterpData(vector<FldArrayI*>& connectIn,
                           vector<FldArrayF*>& fieldsIn,
                           vector<char*>& varStringIn,
                           vector<char*>& eltsIn, 
                           vector<FldArrayI*>& connectOut,
                           vector<FldArrayF*>& fieldsOut,
                           vector<char*>& varStringOut,
                           vector<char*>& eltsOut, 
                           vector<E_Int>& posxt, 
                           vector<E_Int>& posyt,
                           vector<E_Int>& poszt,
                           vector<E_Int>& posct,
                           vector<K_KINTERP::KMesh*>& listOfUnstrMeshes,
                           vector<K_KINTERP::BlkInterpData*>& interpDatau)
{
  E_Int posx, posy, posz, posc;
  E_Int size = fieldsIn.size();
  if ( varStringIn.size() != fieldsIn.size() )
  {
    printf("K_KINTERP::buildListOfUnstrInterpData: size of varString and fields must be equal.\n");
    return;
  }
  for (E_Int v = 0; v < size; v++)
  {
    char* eltType = eltsIn[v];
    if ( strcmp(eltType, "TETRA") != 0 ) 
    {
      printf("K_KINTERP::buildListOfUnstrInterpData: only TETRA element type are valid.");
      printf("Array " SF_D_ " skipped.\n",v+1); 
    }
    else 
    {
      posx = isCoordinateXPresent(varStringIn[v]);
      posy = isCoordinateYPresent(varStringIn[v]);
      posz = isCoordinateZPresent(varStringIn[v]);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        printf("K_KINTERP::buildListOfUnstrInterpData: coordinate not found in array " SF_D_ ". ",v+1);
        printf("Skipped.\n");
      }
      else 
      {
        posc = isCellNatureField2Present(varStringIn[v]);
        posx++; posy++; posz++; posc++;
        
        posxt.push_back(posx);
        posyt.push_back(posy);
        poszt.push_back(posz);
        posct.push_back(posc);
  
        fieldsOut.push_back(fieldsIn[v]);
        varStringOut.push_back(varStringIn[v]);
        connectOut.push_back(connectIn[v]);
        eltsOut.push_back(eltsIn[v]);
        E_Int npts = fieldsIn[v]->getSize();
        FldArrayF coord(npts, 3);
        coord.setOneField(*fieldsIn[v], posx, 1);
        coord.setOneField(*fieldsIn[v], posy, 2);
        coord.setOneField(*fieldsIn[v], posz, 3);
     
        K_KINTERP::KMesh* msh = new K_KINTERP::KMesh(coord, *connectIn[v]);
        listOfUnstrMeshes.push_back(msh);// pour pouvoir le detruire a la fin
        
        // Creation of  InterpData
        K_KINTERP::BlkInterpData* interpData = new K_KINTERP::BlkInterpAdt(*msh);
        interpDatau.push_back(interpData);
      }
    }
  }
  return;
}

//=============================================================================
/* Construit une liste d'interpDatas pour des maillages structures: 
   les vecteurs d'entree sont modifies: seuls sont conserves ceux ayant
   une grille correcte (ni > 2 et nj > 2 et nk > 2 et possedant des coords)
   IN: fieldsIn: liste des champs contenant les infos sur les grilles structurees
   IN: varStringIn: liste des varstring correspondantes
   IN: nitin, njtin, nktin: dimensions des grilles
   OUT: fieldsOut: liste des champs contenant les infos sur les grilles
   OUT: varStringOut: liste des varstring correspondantes
   OUT: nitout, njtout, nktout: dimensions des grilles
   OUT: posxt, posyt, poszt, posct: positions de x,y,z,celln (eventuellt)
   OUT: listOfMeshes: c'est pour pouvoir les detruire 
   OUT: interpDatas: liste des interpDatas crees basees sur listOfMeshes
*/
//=============================================================================
void K_KINTERP::
buildListOfStructInterpData(vector<FldArrayF*>& fieldsIn,
                            vector<char*>& varStringIn,
                            vector<E_Int>& nitin, 
                            vector<E_Int>& njtin,
                            vector<E_Int>& nktin,
                            vector<FldArrayF*>& fieldsOut,
                            vector<char*>& varStringOut,
                            vector<E_Int>& nitout, 
                            vector<E_Int>& njtout,
                            vector<E_Int>& nktout,
                            vector<E_Int>& posxt, 
                            vector<E_Int>& posyt,
                            vector<E_Int>& poszt,
                            vector<E_Int>& posct,
                            vector<K_KINTERP::KMesh*>& listOfStructMeshes,
                            vector<K_KINTERP::BlkInterpData*>& interpDatas)
{
  E_Int posx, posy, posz, posc;
  
  E_Int size = fieldsIn.size();
  if (varStringIn.size() != fieldsIn.size())
  {
    printf("K_KINTERP::buildListOfStructInterpData: size of varString and fields must be equal.\n");
    return;
  }
  for (E_Int v = 0; v < size; v++)
  {
    E_Int ni = nitin[v]; E_Int nj = njtin[v]; E_Int nk = nktin[v];

    if (ni < 2 || nj < 2 || nk < 2)
    {
      printf("K_KINTERP::buildListOfStructInterpData: only 3D arrays are valid.");
      printf("Array " SF_D_ " skipped.\n", v+1);
    }
    else 
    {
      posx = isCoordinateXPresent(varStringIn[v]);
      posy = isCoordinateYPresent(varStringIn[v]);
      posz = isCoordinateZPresent(varStringIn[v]);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        printf("K_KINTERP::buildListOfStructInterpData: coordinate not found in array " SF_D_ ". ", v+1);
        printf("Skipped.\n");
      }
      else
      {
        posc = isCellNatureField2Present(varStringIn[v]);
        posx++; posy++; posz++; posc++;
        
        posxt.push_back(posx); posyt.push_back(posy); poszt.push_back(posz);
        posct.push_back(posc);
        
        nitout.push_back(ni); njtout.push_back(nj); nktout.push_back(nk);
        
        FldArrayF& f = *fieldsIn[v];
        E_Int npts = f.getSize();
        fieldsOut.push_back(fieldsIn[v]);
        varStringOut.push_back(varStringIn[v]);
        
        FldArrayF coord(npts, 3);
        coord.setOneField(f, posx, 1);
        coord.setOneField(f, posy, 2);
        coord.setOneField(f, posz, 3);
        K_KINTERP::KMesh* msh = new K_KINTERP::KMesh(ni, nj, nk, coord);

        listOfStructMeshes.push_back(msh);
        
        // Creation of  InterpData
        K_KINTERP::BlkInterpData* interpData = new K_KINTERP::BlkInterpAdt(*msh);  
        interpDatas.push_back(interpData);
      }
    }
  }  
  return;
}
//=============================================================================
/* Calcul d'une valeur interpolee en un pt a partir du champ nofld de f0
   IN: ni, nj, nk: dimensions de la grille d interpolation structuree
   IN: indiSize: taille de indi (indices des pts +1 pour le type)
   IN: indi: type + indices des pts d interpolation
   IN: cf: coefs d interpolation
   IN: cellNp: cellnaturefield  du bloc donneur
   OUT: valeur interpolee du champ cellN en ind
*/
//=============================================================================
E_Float K_KINTERP::compInterpolatedNature(E_Int ni, E_Int nj, E_Int nk,
                                          E_Int indiSize, E_Int* indi,
                                          FldArrayF& cf, E_Int* cellNp)
{
  E_Int ninj = ni*nj;
  E_Float val = 0.;
  E_Int ncfsize = cf.getSize();
  vector<E_Int> indt(ncfsize);
  E_Int indint, order;
  E_Int interpType = indi[0];
  switch (interpType)
  {
    case  3:
    case  5: 
      order = K_FUNC::E_abs(interpType);
      for (E_Int i0 = 1; i0 <= order; i0++)
        for (E_Int j0 = order+1; j0 <= 2*order;j0++)
          for (E_Int k0 = 2*order+1; k0 <= 3*order; k0++)
          {
            indint = indi[i0] + indi[j0] * ni + indi[k0] * ninj;
            val += cf[i0-1]*cf[j0-1]*cf[k0-1]*cellNp[indint]* (2.-cellNp[indint]);
          }
      break;
      
    case 2:
      indt[0] = indi[1] + indi[3] * ni + indi[5] * ninj;
      indt[1] = indi[2] + indi[3] * ni + indi[5] * ninj;
      indt[2] = indi[1] + indi[4] * ni + indi[5] * ninj;
      indt[3] = indi[2] + indi[4] * ni + indi[5] * ninj;
      indt[4] = indi[1] + indi[3] * ni + indi[6] * ninj;
      indt[5] = indi[2] + indi[3] * ni + indi[6] * ninj;
      indt[6] = indi[1] + indi[4] * ni + indi[6] * ninj;
      indt[7] = indi[2] + indi[4] * ni + indi[6] * ninj;
      // test interpolation pour un centre de cellule : 
      // - pas de cellule masquee, ni interpolee dans la molecule d interpolation, sauf si son cf associe est nul.
      // - somme des cf = 1.
      // -----------------
      // cellN*(2-cellN) renvoie 0 si cellN = 0 ou 2 (pt masque ou interpolee) et 1 si cellN =1 (pt calcule)
      E_Int ind0;
      for (E_Int ncf = 0; ncf < ncfsize; ncf++)
      {
        ind0 = indt[ncf]; 
        val += cf[ncf]*cellNp[ind0]*(2.-cellNp[ind0]);
      }
      break;

    default:
      printf("K_KINTERP::compInterpolatedNature : type of interpolation not implemented.\n");
      exit(0);
  }
  return val;
}
//=============================================================================
/* Calcul d'une valeur interpolee en un pt a partir du champ nofld de f0
   IN: ni, nj, nk: dimensions de la grille d interpolation structuree
   IN: indiSize: taille de indi (indices des pts +1 pour le type)
   IN: indi: type + indices des pts d interpolation
   IN: cf: coefs d interpolation
   IN: cellNp: cellnaturefield  du bloc donneur
   OUT: valeur interpolee du champ cellN en ind
*/
//=============================================================================
E_Float K_KINTERP::compInterpolatedNatureEX(E_Int ni, E_Int nj, E_Int nk,
                                            E_Int indiSize, E_Int* indi,
                                            FldArrayF& cf, E_Int* cellNp)
{
  E_Int ninj = ni*nj;
  E_Float val = 0.;
  E_Int ncfsize = cf.getSize();
  vector<E_Int> indt(ncfsize);
  E_Int indint, order;
  E_Int interpType = indi[0];
  switch (interpType)
  {
    case  3:
    case  5: 
      order = K_FUNC::E_abs(interpType);
      for (E_Int i0 = 1; i0 <= order; i0++)
        for (E_Int j0 = order+1; j0 <= 2*order;j0++)
          for (E_Int k0 = 2*order+1; k0 <= 3*order; k0++)
          {
            indint = indi[i0] + indi[j0] * ni + indi[k0] * ninj;
            val += cf[i0-1]*cf[j0-1]*cf[k0-1]*0.5*cellNp[indint]* (3.-cellNp[indint]);  
          }
      break;
      
    case 2:
      indt[0] = indi[1] + indi[3] * ni + indi[5] * ninj;
      indt[1] = indi[2] + indi[3] * ni + indi[5] * ninj;
      indt[2] = indi[1] + indi[4] * ni + indi[5] * ninj;
      indt[3] = indi[2] + indi[4] * ni + indi[5] * ninj;
      indt[4] = indi[1] + indi[3] * ni + indi[6] * ninj;
      indt[5] = indi[2] + indi[3] * ni + indi[6] * ninj;
      indt[6] = indi[1] + indi[4] * ni + indi[6] * ninj;
      indt[7] = indi[2] + indi[4] * ni + indi[6] * ninj;
      // test interpolation pour un centre de cellule : 
      // - pas de cellule masquee, ni interpolee dans la molecule d interpolation, sauf si son cf associe est nul.
      // - somme des cf = 1.
      // -----------------
      // cellN*(2-cellN) renvoie 0 si cellN = 0 ou 2 (pt masque ou interpolee) et 1 si cellN =1 (pt calcule)
      E_Int ind0;
      for (E_Int ncf = 0; ncf < ncfsize; ncf++)
      {
        ind0 = indt[ncf]; 
        val += cf[ncf]*0.5*cellNp[ind0]*(3.-cellNp[ind0]);
      }
      break;

    default:
      printf("K_KINTERP::compInterpolatedNature : type of interpolation not implemented.\n");
      exit(0);
  }
  return val;
}
//=============================================================================
/* IN: ni, nj, nk: dimensions de f0 dans le cas structure
        -1, -1, -1 si le donneur est non structure
   IN: f0: champ du bloc donneur pour l interpolation
   IN: indi: indices des pts de la molecule d interpolation
             peut etre defini par des sommets de la molecule donneuse
                              par l indice de la cellule donneuse
             doit etre coherent avec f.
   IN: indiSize: taille de indi
   IN: cf: coefs d'interpolation associes 
   IN: ind: indice du pt a interpoler 
   IN: interpType: permet de determiner la formule appliquee
   OUT: f: champs interpoles au point ind 
   Retourne -1 si erreur 
*/
//=============================================================================
short K_KINTERP::compInterpolatedValues(
  E_Int ni, E_Int nj, E_Int nk,  
  E_Int indiSize, E_Int* indi, FldArrayF& cf,
  FldArrayF& f0,  FldArrayI& cn0, 
  E_Int ind, FldArrayF& f,
  E_Int interpType)
{
  E_Int nfld = f0.getNfld();
  E_Int ind0;
  //E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7, indcell; 
  E_Int nvert, nocf, indv;
  E_Int ninj = ni*nj;
  //E_Int nic = K_FUNC::E_max(1,ni-1); //cas ou ni,nj,nk decrivent les dimensions d une grille f0 en noeuds 
  //E_Int njc = K_FUNC::E_max(1,nj-1); 

  E_Float* cfp = cf.begin();
  E_Int ncfPerDir; // nb de coefs par direction dans le cas OiABC
 
  switch (interpType) 
  {
    case  2:// les indices et f sont localises de maniere identique, ex O2CF
      if ( ni == -1 )  // NS : indices des nvert sommets de la molecule d interpolation
      {
        nvert = indiSize;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          E_Float* fp = f.begin(eq);
          E_Float* fp0 = f0.begin(eq);
          fp[ind] = 0.;
          for (E_Int nov = 0; nov < nvert; nov++)
          {
            indv = indi[nov];
            fp[ind] += cfp[nov]*fp0[indv];
          }          
        }
      }
      else //STRUCT : indices [i,i+1,j,j+1,k,k+1] stockes
      {
        ncfPerDir = 2;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          E_Float* fp = f.begin(eq);
          E_Float* fp0 = f0.begin(eq);
          fp[ind] = 0.; nocf = 0;
          for (E_Int k0 = 2*ncfPerDir; k0 < 3*ncfPerDir; k0++)
            for (E_Int j0 = ncfPerDir; j0 < 2*ncfPerDir; j0++)
              for (E_Int i0 = 0; i0 < ncfPerDir; i0++)
              {
                ind0 = indi[i0]+indi[j0]*ni+indi[k0]*ninj;
                fp[ind] += cfp[nocf]*fp0[ind0];
                nocf+=1;
              }
        }
      }
      break;
  
    case  3://formule directionnelle (ex O3ABC) : on stocke les indices i,j,k des sommets 
      if (ni == -1)
      {
        printf("BlkInterp::compInterpolatedValues: interpType 3 not valid for unstructured interpolations.");
        return -1;
      }
      ncfPerDir = 3;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f.begin(eq);
        E_Float* fp0 = f0.begin(eq);
        fp[ind] = 0.;
        for (E_Int k0 = 2*ncfPerDir; k0 < 3*ncfPerDir; k0++)
          for (E_Int j0 = ncfPerDir; j0 < 2*ncfPerDir;j0++)
            for (E_Int i0 = 0; i0 < ncfPerDir; i0++)
            {
              ind0 = indi[i0] + indi[j0] * ni + indi[k0] * ninj;
              fp[ind] += cfp[i0]*cfp[j0]*cfp[k0]* fp0[ind0];  
            }
      }
      break;

    case  5://formule directionnelle (ex O5ABC)
      if (ni== -1)
      {
        printf("BlkInterp::compInterpolatedValues: interpType 5 not valid for unstructured interpolations.");
        return -1;
      }
      ncfPerDir = 5;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f.begin(eq);
        E_Float* fp0 = f0.begin(eq);
        fp[ind] = 0.;
        for (E_Int k0 = 2*ncfPerDir; k0 <3*ncfPerDir; k0++)
          for (E_Int j0 = ncfPerDir; j0 <2*ncfPerDir;j0++)
            for (E_Int i0 = 0; i0 < ncfPerDir; i0++)
            {
              ind0 = indi[i0] + indi[j0] * ni + indi[k0] * ninj;
              fp[ind] += cfp[i0]*cfp[j0]*cfp[k0]* fp0[ind0];  
            }
      }
      break;
      
    default:
      printf("BlkInterp::compInterpolatedValues: unknown interpolation type.");
      return -1;      
  }
  return 0;
}

//========================KCore/Interp/BlkInterp.cpp===========================
