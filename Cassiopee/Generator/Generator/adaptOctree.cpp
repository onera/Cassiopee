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

#include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Adapt the octree */
//=============================================================================
PyObject* K_GENERATOR::adaptOctree(PyObject* self, PyObject* args)
{
  PyObject *octree, *indica;
  if (!PYPARSETUPLE_(args, OO_, &octree, &indica)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(octree, varString, f, ni, nj, nk, cn, 
                                     eltType);
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(octree, f); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree: array must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree: the octree must be HEXA or QUAD.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree: the octree must contain coordinates.");
    return NULL;  
  }
  posx++; posy++; posz++;
 
  // Check indicator
  E_Int nii, nji, nki;
  FldArrayF* fi; FldArrayI* cni;
  char* varStringi; char* eltTypei;
  E_Int resi = K_ARRAY::getFromArray3(indica, varStringi, fi, 
                                      nii, nji, nki, cni, eltTypei);
  if (resi != 1 && resi != 2) 
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree: indicator array must be structured.");
    return NULL;
  }
  E_Int posi = K_ARRAY::isNamePresent("indicator", varStringi);
  if (posi == -1) 
  { 
    RELEASESHAREDB(resi, indica, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: adaptOctree: no refinement indicator in the octree. Nothing done."); 
    return octree;
  }
  posi++;
  if (fi->getSize() != cn->getSize()) 
  {
    RELEASESHAREDB(resi, indica, fi, cni); RELEASESHAREDU(octree, f, cn);
    printf("Warning: adaptOctree: refinement indicator size must be equal to the number of elements. Nothing done."); 
    return octree;
  }
  /*-----------------------------*/
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Float* indict = fi->begin(posi);
  E_Int nelts = cn->getSize(); E_Int nvert = cn->getNfld();
  E_Int npts = f->getSize(); E_Int api = f->getApi();
  FldArrayF* fo = new FldArrayF(npts, 3);
  E_Float* fox = fo->begin(1); E_Float* foy = fo->begin(2); 
  E_Float* foz = fo->begin(3); 
  FldArrayI* cno = new FldArrayI(nelts, nvert);
  FldArrayF* indicout = new FldArrayF(nelts); indicout->setAllValuesAtNull();
  E_Int no = 0; E_Int eto = 0;
  E_Float indic;

  // calcul de la bbox de l'octree
  E_Float xmin = K_CONST::E_MAX_FLOAT; E_Float xmax =-K_CONST::E_MAX_FLOAT;
  E_Float ymin = K_CONST::E_MAX_FLOAT; E_Float ymax =-K_CONST::E_MAX_FLOAT;
  E_Float zmin = K_CONST::E_MAX_FLOAT; E_Float zmax =-K_CONST::E_MAX_FLOAT;
  FldArrayF bbox(nelts, 6);
  K_COMPGEOM::boundingBoxOfUnstrCells(*cn, xt, yt, zt, bbox);
  E_Float* xmint = bbox.begin(1); E_Float* xmaxt = bbox.begin(4);
  E_Float* ymint = bbox.begin(2); E_Float* ymaxt = bbox.begin(5);
  E_Float* zmint = bbox.begin(3); E_Float* zmaxt = bbox.begin(6);
  for (E_Int et = 0; et < nelts; et++)
  {
    xmin = K_FUNC::E_min(xmin,xmint[et]); xmax = K_FUNC::E_max(xmax,xmaxt[et]);
    ymin = K_FUNC::E_min(ymin,ymint[et]); ymax = K_FUNC::E_max(ymax,ymaxt[et]);
    zmin = K_FUNC::E_min(zmin,zmint[et]); zmax = K_FUNC::E_max(zmax,zmaxt[et]);
  }
  bbox.malloc(0,1);

  // calcul du pas sur chaque grille 
  FldArrayF dht(nelts);
  E_Float* dhtp = dht.begin();
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  
#pragma omp parallel default(shared) if (nelts > 100)
  {
    E_Int ind1, ind2;
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      ind1 = cn1[et]-1; ind2 = cn2[et]-1;
      dhtp[et] = xt[ind2]-xt[ind1];
    }
  }
  vector< vector<E_Int> > cEEN(nelts);
  E_Int ok = getNeighbourElts(npts, xt, yt, zt, *cn, cEEN); 
  if (ok == 0) 
  {
    RELEASESHAREDB(resi, indica, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: adaptOctree: cannot compute connectivity. Check connectivity."); 
    return octree;
  }
  // modification de l'indicateur: si une cellule de niveau l a un indicateur 
  // a 1 et une cellule adjacente est de niveau l+1, indicateur < 1: mettre 
  // indicateur(l+1) a 1 pour un saut de 2*h au pire si raffinement entre 
  // deux elements adjacents (reequilibrage)
  modifyIndicator(nelts, dht, cEEN, posi, indict);
  //raffinement / deraffinement
  FldArrayIS dejaVu(nelts); dejaVu.setAllValuesAtNull();
  short* dejaVup = dejaVu.begin();
  FldArrayI indir(npts); indir.setAllValuesAt(-1);
  E_Int* indirp = indir.begin();

  // step1 - elements a conserver tels quels
  for (E_Int et = 0; et < nelts; et++)
  {
    indic = indict[et];
    if (dejaVup[et] == 0 && K_FUNC::E_abs(indic) < 1.e-6)
    {
      for (E_Int nov = 1; nov <= nvert; nov++)
      {
        E_Int ind = (*cn)(et, nov)-1;
        if (indirp[ind] == -1)
        {
          fox[no] = xt[ind]; foy[no] = yt[ind]; foz[no] = zt[ind];
          no++; (*cno)(eto, nov) = no; indirp[ind] = no;
          if (no == fo->getSize())
          {
            fo->reAllocMat(no+npts, 3);
            fox = fo->begin(1); foy = fo->begin(2); foz = fo->begin(3);
          }
        }
        else (*cno)(eto, nov) = indirp[ind];
      }
      (*indicout)[eto] = 0.; eto++;
      if (eto == cno->getSize())
      {
        cno->reAllocMat(eto+nelts, nvert); 
        indicout->reAlloc(eto);
      }
      dejaVup[et] = 1;
    }
  }

  // raffinement/deraffinement
  for (E_Int et = 0; et < nelts; et++)
  {
    if (dejaVup[et] == 0)
    {
      indic = indict[et]; 
      if (indic > 0.5)
      {
        splitElement(et, npts, nelts, indic, indirp, *cn, xt, yt, zt, 
                     *cno, *fo, *indicout, no, eto);
        fox = fo->begin(1); foy = fo->begin(2); foz = fo->begin(3);
        dejaVup[et] = 1; 
      }
      else if (indic < -0.5)
      {                     
        // si merge: dejaVu de et et ses voisins mis a 1
        ok = mergeOctreeElement(et, npts, indic, *cn, 
                                xmin, ymin, zmin, xmax, ymax, zmax,
                                xt, yt, zt, dht.begin(1), indict,
                                *cno, *fo, *indicout, no, eto, dejaVu);
        fox = fo->begin(1); foy = fo->begin(2); foz = fo->begin(3);
        if (ok == -1)
        {
          for (E_Int nov = 1; nov <= nvert; nov++)
          {
            E_Int ind = (*cn)(et, nov)-1;
            fox[no] = xt[ind]; foy[no] = yt[ind]; foz[no] = zt[ind];
            no++; (*cno)(eto, nov) = no; 
            if (no == fo->getSize()) 
            { 
              fo->reAllocMat(fo->getSize()+npts, 3);
              fox = fo->begin(1); foy = fo->begin(2); foz = fo->begin(3);
            }
          }
          (*indicout)[eto] = 0.;
          eto++;
          if (eto == cno->getSize() )
          {
            cno->reAllocMat(cno->getSize()+nelts,nvert); 
            indicout->reAlloc(eto);
          }
          dejaVu[et] = 1;
        }
      }
    }
  }

  // Nettoyage
  for (E_Int v = 0; v < nelts; v++) cEEN[v].clear();
  cEEN.clear();

  // Realloc
  fo->reAllocMat(no,3); cno->reAllocMat(eto, nvert); indicout->reAlloc(eto);
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, eltType, *fo, *cno);

  // Sortie
  RELEASESHAREDU(octree, f, cn); RELEASESHAREDB(resi, indica, fi, cni);
  PyObject* l = PyList_New(0); 
  PyObject* tpl = K_ARRAY::buildArray3(*fo, "x,y,z", *cno, eltType, api);
  delete fo;
  PyList_Append(l, tpl); Py_DECREF(tpl);
  char eltType2[8];
  if (strcmp(eltType, "QUAD")  == 0) strcpy(eltType2, "QUAD*");
  else strcpy(eltType2, "HEXA*");
  PyObject* tpl2 = K_ARRAY::buildArray3(*indicout, "indicator", *cno, 
                                        eltType2, api);
  delete cno; delete indicout;
  PyList_Append(l, tpl2); Py_DECREF(tpl2);
  return l;
}

//=============================================================================
/* Split un element en 4/8. Le tableau indir donne la nouvelle numerotation 
   d'un vertex si deja mis */
//=============================================================================
void K_GENERATOR::splitElement(
  E_Int et, E_Int npts, E_Int nelts, E_Float indic, E_Int* indir,
  FldArrayI& cn, E_Float* xt, E_Float* yt, E_Float* zt, 
  FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout, 
  E_Int& no, E_Int& eto)
{
  E_Int neltso = cno.getSize(); E_Int nvert = cno.getNfld();
  E_Int nptso = fo.getSize(); 

  if (no >= nptso-9*nvert) 
  { nptso += npts; fo.reAllocMat(nptso,3); }
  if (eto >= neltso-17)
  { neltso += nelts; cno.reAllocMat(neltso,nvert); indicout.reAlloc(neltso); }

  E_Float* xo = fo.begin(1); E_Float* yo = fo.begin(2); E_Float* zo = fo.begin(3);
  E_Float* indico = indicout.begin(1);
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); 
  E_Int* cn3 = cn.begin(3); E_Int* cn4 = cn.begin(4);
  if (nvert == 4) // QUAD
  {
    E_Int indA = cn1[et]-1; E_Int indB = cn2[et]-1; 
    E_Int indC = cn3[et]-1; E_Int indD = cn4[et]-1;
    E_Float xmean = (xt[indB] + xt[indA])*0.5;
    E_Float ymean = (yt[indD] + yt[indA])*0.5;
    E_Float z = zt[indA];
    E_Float xA = xt[indA]; E_Float yA = yt[indA];
    E_Float xB = xt[indB]; E_Float yB = yt[indB];
    E_Float xC = xt[indC]; E_Float yC = yt[indC];
    E_Float xD = xt[indD]; E_Float yD = yt[indD];

    // 1er quad
    if (indir[indA] == -1)
    {
      xo[no] = xA; yo[no] = yA; zo[no] = z;
      no++; cno(eto,1) = no;
      indir[indA] = no;
    }
    else cno(eto,1) = indir[indA];

    xo[no] = xmean; yo[no] = yA; zo[no] = z;
    no++; cno(eto,2) = no;

    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,3) = no; 

    xo[no] = xA; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
    
    // 2eme quad
    xo[no] = xmean; yo[no] = yB; zo[no] = z;
    no++; cno(eto,1) = no; 

    if (indir[indB] == -1)
    {
      xo[no] = xB; yo[no] = yB; zo[no] = z;
      no++; cno(eto,2) = no;
      indir[indB] = no;
    }
    else cno(eto,2) = indir[indB];

    xo[no] = xB; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,3) = no; 

    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 3eme quad
    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,1) = no; 

    xo[no] = xC; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,2) = no; 

    if (indir[indC] == -1)
    {
      xo[no] = xC; yo[no] = yC; zo[no] = z;
      no++; cno(eto,3) = no;
      indir[indC] = no;
    }
    else cno(eto,3) = indir[indC];

    xo[no] = xmean; yo[no] = yC; zo[no] = z;
    no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
    
    // 4eme quad
    xo[no] = xD; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,1) = no; 

    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno(eto,2) = no; 

    xo[no] = xmean; yo[no] = yD; zo[no] = z;
    no++; cno(eto,3) = no; 

    if (indir[indD] == -1)
    {
      xo[no] = xD; yo[no] = yD; zo[no] = z;
      no++; cno(eto,4) = no;
      indir[indD] = no;
    }
    else cno(eto,4) = indir[indD];

    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
  }
  else // HEXA
  {
    E_Int* cn5 = cn.begin(5); E_Int* cn6 = cn.begin(6);  
    E_Int* cn7 = cn.begin(7);  E_Int* cn8 = cn.begin(8);
    E_Int indA = cn1[et]-1; E_Int indB = cn2[et]-1; 
    E_Int indC = cn3[et]-1; E_Int indD = cn4[et]-1;
    E_Int indE = cn5[et]-1; E_Int indF = cn6[et]-1; 
    E_Int indG = cn7[et]-1; E_Int indH = cn8[et]-1;
    E_Float xA = xt[indA]; E_Float yA = yt[indA]; E_Float zA = zt[indA];
    E_Float xB = xt[indB]; E_Float yB = yt[indB]; E_Float zB = zt[indB];
    E_Float xC = xt[indC]; E_Float yC = yt[indC]; E_Float zC = zt[indC];
    E_Float xD = xt[indD]; E_Float yD = yt[indD]; E_Float zD = zt[indD];
    E_Float xE = xt[indE]; E_Float yE = yt[indE]; E_Float zE = zt[indE];
    E_Float xF = xt[indF]; E_Float yF = yt[indF]; E_Float zF = zt[indF];
    E_Float xG = xt[indG]; E_Float yG = yt[indG]; E_Float zG = zt[indG];
    E_Float xH = xt[indH]; E_Float yH = yt[indH]; E_Float zH = zt[indH];
    E_Float xmean = (xt[indB] + xt[indA])*0.5;
    E_Float ymean = (yt[indD] + yt[indA])*0.5;
    E_Float zmean = (zt[indE] + zt[indA])*0.5;
    // 1er HEXA
    if (indir[indA] == -1)
    {
      xo[no] = xA; yo[no] = yA; zo[no] = zA;
      no++; cno(eto,1) = no;
      indir[indA] = no;
    }
    else cno(eto,1) = indir[indA];

//     xo[no] = xA; yo[no] = yA; zo[no] = zA; no++; cno(eto,1) = no; 
    xo[no] = xmean; yo[no] = yA; zo[no] = zA; no++; cno(eto,2) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zA; no++; cno(eto,3) = no; 
    xo[no] = xA; yo[no] = ymean; zo[no] = zA; no++; cno(eto,4) = no; 
    xo[no] = xA; yo[no] = yA; zo[no] = zmean; no++; cno(eto,5) = no; 
    xo[no] = xmean; yo[no] = yA; zo[no] = zmean; no++; cno(eto,6) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,7) = no; 
    xo[no] = xA; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 2eme HEXA
    xo[no] = xmean; yo[no] = yB; zo[no] = zB; no++; cno(eto,1) = no; 
    if (indir[indB] == -1)
    {
      xo[no] = xB; yo[no] = yB; zo[no] = zB;
      no++; cno(eto,2) = no;
      indir[indB] = no;
    }
    else cno(eto,2) = indir[indB];
    //xo[no] = xB; yo[no] = yB; zo[no] = zB; no++; cno(eto,2) = no; 
    xo[no] = xB; yo[no] = ymean; zo[no] = zB; no++; cno(eto,3) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zB; no++; cno(eto,4) = no; 
    xo[no] = xmean; yo[no] = yB; zo[no] = zmean; no++; cno(eto,5) = no; 
    xo[no] = xB; yo[no] = yB; zo[no] = zmean; no++; cno(eto,6) = no; 
    xo[no] = xB; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,7) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 3eme HEXA
    xo[no] = xmean; yo[no] = ymean; zo[no] = zC; no++; cno(eto,1) = no; 
    xo[no] = xC; yo[no] = ymean; zo[no] = zC; no++; cno(eto,2) = no; 
    if (indir[indC] == -1)
    {
      xo[no] = xC; yo[no] = yC; zo[no] = zC;
      no++; cno(eto,3) = no;
      indir[indC] = no;
    }
    else cno(eto,3) = indir[indC];
//     xo[no] = xC; yo[no] = yC; zo[no] = zC; no++; cno(eto,3) = no; 
    xo[no] = xmean; yo[no] = yC; zo[no] = zC; no++; cno(eto,4) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,5) = no; 
    xo[no] = xC; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,6) = no; 
    xo[no] = xC; yo[no] = yC; zo[no] = zmean; no++; cno(eto,7) = no; 
    xo[no] = xmean; yo[no] = yC; zo[no] = zmean; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 4eme HEXA
    xo[no] = xD; yo[no] = ymean; zo[no] = zD; no++; cno(eto,1) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zD; no++; cno(eto,2) = no; 
    xo[no] = xmean; yo[no] = yD; zo[no] = zD; no++; cno(eto,3) = no; 
    if (indir[indD] == -1)
    {
      xo[no] = xD; yo[no] = yD; zo[no] = zD;
      no++; cno(eto,4) = no;
      indir[indD] = no;
    }
    else cno(eto,4) = indir[indD];
//     xo[no] = xD; yo[no] = yD; zo[no] = zD; no++; cno(eto,4) = no; 
    xo[no] = xD; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,5) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,6) = no; 
    xo[no] = xmean; yo[no] = yD; zo[no] = zmean; no++; cno(eto,7) = no; 
    xo[no] = xD; yo[no] = yD; zo[no] = zmean; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 5eme HEXA
    xo[no] = xE; yo[no] = yE; zo[no] = zmean; no++; cno(eto,1) = no; 
    xo[no] = xmean; yo[no] = yE; zo[no] = zmean; no++; cno(eto,2) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,3) = no; 
    xo[no] = xE; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,4) = no; 
    if (indir[indE] == -1)
    {
      xo[no] = xE; yo[no] = yE; zo[no] = zE;
      no++; cno(eto,5) = no;
      indir[indE] = no;
    }
    else cno(eto,5) = indir[indE];
    //xo[no] = xE; yo[no] = yE; zo[no] = zE; no++; cno(eto,5) = no; 
    xo[no] = xmean; yo[no] = yE; zo[no] = zE; no++; cno(eto,6) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zE; no++; cno(eto,7) = no; 
    xo[no] = xE; yo[no] = ymean; zo[no] = zE; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 6eme HEXA
    xo[no] = xmean; yo[no] = yF; zo[no] = zmean; no++; cno(eto,1) = no; 
    xo[no] = xF; yo[no] = yF; zo[no] = zmean; no++; cno(eto,2) = no; 
    xo[no] = xF; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,3) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,4) = no; 
    xo[no] = xmean; yo[no] = yF; zo[no] = zF; no++; cno(eto,5) = no; 
    if (indir[indF] == -1)
    {
      xo[no] = xF; yo[no] = yF; zo[no] = zF;
      no++; cno(eto,6) = no;
      indir[indF] = no;
    }
    else cno(eto,6) = indir[indF];
    //xo[no] = xF; yo[no] = yF; zo[no] = zF; no++; cno(eto,6) = no; 
    xo[no] = xF; yo[no] = ymean; zo[no] = zF; no++; cno(eto,7) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zF; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 7eme HEXA
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,1) = no; 
    xo[no] = xG; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,2) = no; 
    xo[no] = xG; yo[no] = yG; zo[no] = zmean; no++; cno(eto,3) = no; 
    xo[no] = xmean; yo[no] = yG; zo[no] = zmean; no++; cno(eto,4) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zG; no++; cno(eto,5) = no; 
    xo[no] = xG; yo[no] = ymean; zo[no] = zG; no++; cno(eto,6) = no; 
    if (indir[indG] == -1)
    {
      xo[no] = xG; yo[no] = yG; zo[no] = zG;
      no++; cno(eto,7) = no;
      indir[indG] = no;
    }
    else cno(eto,7) = indir[indG];
    //xo[no] = xG; yo[no] = yG; zo[no] = zG; no++; cno(eto,7) = no; 
    xo[no] = xmean; yo[no] = yG; zo[no] = zG; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 8eme HEXA
    xo[no] = xH; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,1) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno(eto,2) = no; 
    xo[no] = xmean; yo[no] = yH; zo[no] = zmean; no++; cno(eto,3) = no; 
    xo[no] = xH; yo[no] = yH; zo[no] = zmean; no++; cno(eto,4) = no; 
    xo[no] = xH; yo[no] = ymean; zo[no] = zH; no++; cno(eto,5) = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zH; no++; cno(eto,6) = no; 
    xo[no] = xmean; yo[no] = yH; zo[no] = zH; no++; cno(eto,7) = no; 
    //xo[no] = xH; yo[no] = yH; zo[no] = zH; no++; cno(eto,8) = no; 
    if (indir[indH] == -1)
    {
      xo[no] = xH; yo[no] = yH; zo[no] = zH;
      no++; cno(eto,8) = no;
      indir[indH] = no;
    }
    else cno(eto,8) = indir[indH];
 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
  }
}

//=============================================================================
/* Split un element en 4/8 */
//=============================================================================
void K_GENERATOR::splitElement(
  E_Int et, E_Int npts, E_Int nelts, E_Float indic,
  FldArrayI& cn, E_Float* xt, E_Float* yt, E_Float* zt, 
  FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout, 
  E_Int& no, E_Int& eto)
{
  E_Int neltso = cno.getSize(); E_Int nvert = cno.getNfld();
  E_Int nptso = fo.getSize(); 

  if (no >= nptso-9*nvert) 
  { nptso += npts; fo.reAllocMat(nptso,3); }
  if (eto >= neltso-17) 
  { nelts += nelts; cno.reAllocMat(nelts, nvert); indicout.reAlloc(nelts); }

  E_Float* xo = fo.begin(1); E_Float* yo = fo.begin(2); E_Float* zo = fo.begin(3);
  E_Float* indico = indicout.begin();
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);  
  E_Int* cn3 = cn.begin(3); E_Int* cn4 = cn.begin(4);
  E_Int* cno1 = cno.begin(1); E_Int* cno2 = cno.begin(2);  
  E_Int* cno3 = cno.begin(3); E_Int* cno4 = cno.begin(4);

  if (nvert == 4) // QUAD
  {
    E_Int indA = cn1[et]-1; E_Int indB = cn2[et]-1; 
    E_Int indC = cn3[et]-1; E_Int indD = cn4[et]-1;
    E_Float xmean = (xt[indB] + xt[indA])*0.5;
    E_Float ymean = (yt[indD] + yt[indA])*0.5;
    E_Float z = zt[indA];
    E_Float xA = xt[indA]; E_Float yA = yt[indA];
    E_Float xB = xt[indB]; E_Float yB = yt[indB];
    E_Float xC = xt[indC]; E_Float yC = yt[indC];
    E_Float xD = xt[indD]; E_Float yD = yt[indD];

    // 1er quad
    xo[no] = xA; yo[no] = yA; zo[no] = z;
    no++; cno1[eto] = no;

    xo[no] = xmean; yo[no] = yA; zo[no] = z;
    no++; cno2[eto] = no; 

    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno3[eto] = no; 

    xo[no] = xA; yo[no] = ymean; zo[no] = z;
    no++; cno4[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
    
    // 2eme quad
    xo[no] = xmean; yo[no] = yB; zo[no] = z;
    no++; cno1[eto] = no; 

    xo[no] = xB; yo[no] = yB; zo[no] = z;
    no++; cno2[eto] = no; 

    xo[no] = xB; yo[no] = ymean; zo[no] = z;
    no++; cno3[eto] = no; 

    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno4[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 3eme quad
    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno1[eto] = no; 

    xo[no] = xC; yo[no] = ymean; zo[no] = z;
    no++; cno2[eto] = no; 

    xo[no] = xC; yo[no] = yC; zo[no] = z;
    no++; cno3[eto] = no; 

    xo[no] = xmean; yo[no] = yC; zo[no] = z;
    no++; cno4[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
    
    // 4eme quad
    xo[no] = xD; yo[no] = ymean; zo[no] = z;
    no++; cno1[eto] = no; 

    xo[no] = xmean; yo[no] = ymean; zo[no] = z;
    no++; cno2[eto] = no; 

    xo[no] = xmean; yo[no] = yD; zo[no] = z;
    no++; cno3[eto] = no; 

    xo[no] = xD; yo[no] = yD; zo[no] = z;
    no++; cno4[eto] = no; 
    
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
  }
  else // HEXA
  {
    E_Int* cn5 = cn.begin(5); E_Int* cn6 = cn.begin(6);  
    E_Int* cn7 = cn.begin(7);  E_Int* cn8 = cn.begin(8);
    E_Int* cno5 = cno.begin(5); E_Int* cno6 = cno.begin(6);  
    E_Int* cno7 = cno.begin(7);  E_Int* cno8 = cno.begin(8);
    E_Int indA = cn1[et]-1; E_Int indB = cn2[et]-1; 
    E_Int indC = cn3[et]-1; E_Int indD = cn4[et]-1;
    E_Int indE = cn5[et]-1; E_Int indF = cn6[et]-1; 
    E_Int indG = cn7[et]-1; E_Int indH = cn8[et]-1;
    E_Float xA = xt[indA]; E_Float yA = yt[indA]; E_Float zA = zt[indA];
    E_Float xB = xt[indB]; E_Float yB = yt[indB]; E_Float zB = zt[indB];
    E_Float xC = xt[indC]; E_Float yC = yt[indC]; E_Float zC = zt[indC];
    E_Float xD = xt[indD]; E_Float yD = yt[indD]; E_Float zD = zt[indD];
    E_Float xE = xt[indE]; E_Float yE = yt[indE]; E_Float zE = zt[indE];
    E_Float xF = xt[indF]; E_Float yF = yt[indF]; E_Float zF = zt[indF];
    E_Float xG = xt[indG]; E_Float yG = yt[indG]; E_Float zG = zt[indG];
    E_Float xH = xt[indH]; E_Float yH = yt[indH]; E_Float zH = zt[indH];
    E_Float xmean = (xt[indB] + xt[indA])*0.5;
    E_Float ymean = (yt[indD] + yt[indA])*0.5;
    E_Float zmean = (zt[indE] + zt[indA])*0.5;

    // 1er HEXA
    xo[no] = xA; yo[no] = yA; zo[no] = zA; no++; cno1[eto] = no; 
    xo[no] = xmean; yo[no] = yA; zo[no] = zA; no++; cno2[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zA; no++; cno3[eto] = no; 
    xo[no] = xA; yo[no] = ymean; zo[no] = zA; no++; cno4[eto] = no; 
    xo[no] = xA; yo[no] = yA; zo[no] = zmean; no++; cno5[eto] = no; 
    xo[no] = xmean; yo[no] = yA; zo[no] = zmean; no++; cno6[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno7[eto] = no; 
    xo[no] = xA; yo[no] = ymean; zo[no] = zmean; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 2eme HEXA
    xo[no] = xmean; yo[no] = yB; zo[no] = zB; no++; cno1[eto] = no; 
    xo[no] = xB; yo[no] = yB; zo[no] = zB; no++; cno2[eto] = no; 
    xo[no] = xB; yo[no] = ymean; zo[no] = zB; no++; cno3[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zB; no++; cno4[eto] = no; 
    xo[no] = xmean; yo[no] = yB; zo[no] = zmean; no++; cno5[eto] = no; 
    xo[no] = xB; yo[no] = yB; zo[no] = zmean; no++; cno6[eto] = no; 
    xo[no] = xB; yo[no] = ymean; zo[no] = zmean; no++; cno7[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 3eme HEXA
    xo[no] = xmean; yo[no] = ymean; zo[no] = zC; no++; cno1[eto] = no; 
    xo[no] = xC; yo[no] = ymean; zo[no] = zC; no++; cno2[eto] = no; 
    xo[no] = xC; yo[no] = yC; zo[no] = zC; no++; cno3[eto] = no; 
    xo[no] = xmean; yo[no] = yC; zo[no] = zC; no++; cno4[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno5[eto] = no; 
    xo[no] = xC; yo[no] = ymean; zo[no] = zmean; no++; cno6[eto] = no; 
    xo[no] = xC; yo[no] = yC; zo[no] = zmean; no++; cno7[eto] = no; 
    xo[no] = xmean; yo[no] = yC; zo[no] = zmean; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 4eme HEXA
    xo[no] = xD; yo[no] = ymean; zo[no] = zD; no++; cno1[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zD; no++; cno2[eto] = no; 
    xo[no] = xmean; yo[no] = yD; zo[no] = zD; no++; cno3[eto] = no; 
    xo[no] = xD; yo[no] = yD; zo[no] = zD; no++; cno4[eto] = no; 
    xo[no] = xD; yo[no] = ymean; zo[no] = zmean; no++; cno5[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno6[eto] = no; 
    xo[no] = xmean; yo[no] = yD; zo[no] = zmean; no++; cno7[eto] = no; 
    xo[no] = xD; yo[no] = yD; zo[no] = zmean; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 5eme HEXA
    xo[no] = xE; yo[no] = yE; zo[no] = zmean; no++; cno1[eto] = no; 
    xo[no] = xmean; yo[no] = yE; zo[no] = zmean; no++; cno2[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno3[eto] = no; 
    xo[no] = xE; yo[no] = ymean; zo[no] = zmean; no++; cno4[eto] = no; 
    xo[no] = xE; yo[no] = yE; zo[no] = zE; no++; cno5[eto] = no; 
    xo[no] = xmean; yo[no] = yE; zo[no] = zE; no++; cno6[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zE; no++; cno7[eto] = no; 
    xo[no] = xE; yo[no] = ymean; zo[no] = zE; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 6eme HEXA
    xo[no] = xmean; yo[no] = yF; zo[no] = zmean; no++; cno1[eto] = no; 
    xo[no] = xF; yo[no] = yF; zo[no] = zmean; no++; cno2[eto] = no; 
    xo[no] = xF; yo[no] = ymean; zo[no] = zmean; no++; cno3[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno4[eto] = no; 
    xo[no] = xmean; yo[no] = yF; zo[no] = zF; no++; cno5[eto] = no; 
    xo[no] = xF; yo[no] = yF; zo[no] = zF; no++; cno6[eto] = no; 
    xo[no] = xF; yo[no] = ymean; zo[no] = zF; no++; cno7[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zF; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 7eme HEXA
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno1[eto] = no; 
    xo[no] = xG; yo[no] = ymean; zo[no] = zmean; no++; cno2[eto] = no; 
    xo[no] = xG; yo[no] = yG; zo[no] = zmean; no++; cno3[eto] = no; 
    xo[no] = xmean; yo[no] = yG; zo[no] = zmean; no++; cno4[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zG; no++; cno5[eto] = no; 
    xo[no] = xG; yo[no] = ymean; zo[no] = zG; no++; cno6[eto] = no; 
    xo[no] = xG; yo[no] = yG; zo[no] = zG; no++; cno7[eto] = no; 
    xo[no] = xmean; yo[no] = yG; zo[no] = zG; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;

    // 8eme HEXA
    xo[no] = xH; yo[no] = ymean; zo[no] = zmean; no++; cno1[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zmean; no++; cno2[eto] = no; 
    xo[no] = xmean; yo[no] = yH; zo[no] = zmean; no++; cno3[eto] = no; 
    xo[no] = xH; yo[no] = yH; zo[no] = zmean; no++; cno4[eto] = no; 
    xo[no] = xH; yo[no] = ymean; zo[no] = zH; no++; cno5[eto] = no; 
    xo[no] = xmean; yo[no] = ymean; zo[no] = zH; no++; cno6[eto] = no; 
    xo[no] = xmean; yo[no] = yH; zo[no] = zH; no++; cno7[eto] = no; 
    xo[no] = xH; yo[no] = yH; zo[no] = zH; no++; cno8[eto] = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.);
    eto++;
  }
}

