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

#include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Balance the octree */
//=============================================================================
PyObject* K_GENERATOR::balanceOctree(PyObject* self, PyObject* args)
{
  E_Int ratio, corners;
  PyObject *octree;
  if (!PYPARSETUPLE_(args, O_ II_, &octree, &ratio, &corners)) return NULL;

  if (corners != 0 && corners != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: corners arg must be 0 or 1.");
    return NULL;
  }

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(octree, varString, f, ni, nj, nk, cn, eltType);
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(octree, f); 
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: array must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: the octree must be HEXA or QUAD.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: the octree must contain coordinates.");
    return NULL;
  }
  posx++; posy++; posz++;
  // if (ratio == 2) checkBalancing2(*cn, *f);
  if (ratio==2) balanceOctree2(*f, *cn, corners);
  else checkBalancing3(*cn, *f);


  PyObject* tpl = K_ARRAY::buildArray(*f, "x,y,z", *cn, -1, eltType);
  RELEASESHAREDU(octree, f, cn); 
  return tpl;
}
/* -----------------------------------------------------------------------
  Equilibrage de l octree en prenant en compte aussi les coins
  Deux elements partageant le meme vertex ne doivent pas etre de niveaux
  d  ecart au plus 1 
  -----------------------------------------------------------------------*/
void K_GENERATOR::balanceOctree2(FldArrayF& coords, FldArrayI& cn, E_Int corners)
{
  E_Int nvert = cn.getNfld();
  E_Int nfld = coords.getNfld();
  E_Int nptsAdd = 5;//9-4;//
  E_Int neltsAdd = 3;
  if ( nvert == 8) {nptsAdd = 19; neltsAdd = 7;}//HEXA
  E_Int splitl = 1;
  while (splitl == 1)
  {
    splitl = 0;
    E_Float* xt = coords.begin(1);
    E_Float* yt = coords.begin(2);
    E_Float* zt = coords.begin(3);
    E_Int npts = coords.getSize();
    E_Int* cn1 = cn.begin(1);
    E_Int* cn2 = cn.begin(2);
    E_Int nelts = cn.getSize();

    //determination des elts voisins
    // c est ce qui prend la majeure partie du temps CPU
    vector< vector<E_Int> > cEEN(nelts);
    getNeighbourElts(npts, xt, yt, zt, cn, cEEN, corners); 

    FldArrayI indic(nelts); indic.setAllValuesAt(-1);
    E_Int elts2Raff = 0;
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int indA = cn1[et]-1;
      E_Int indB = cn2[et]-1;
      E_Float dx = xt[indB]-xt[indA];      

      vector<E_Int>& voisins = cEEN[et];
      for (size_t noetv = 0; noetv < voisins.size(); noetv++)
      {
        E_Int etv = voisins[noetv];
        indA = cn1[etv]-1;
        indB = cn2[etv]-1;
        E_Float dxv = xt[indB]-xt[indA];     
        // au moins 2 niveaux d ecart
        if ( dx/dxv > 3. ) {indic[et] = elts2Raff; elts2Raff+=1; break;}
      }
    }

    if ( elts2Raff > 0) 
    {
      splitl = 1;
      coords.reAllocMat(npts+elts2Raff*nptsAdd,nfld);
      cn.reAllocMat(nelts+elts2Raff*neltsAdd,nvert);

      cn1 = cn.begin(1);
      cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      E_Int* cn4 = cn.begin(4);

      xt = coords.begin(1);
      yt = coords.begin(2);
      zt = coords.begin(3);

      if ( nvert == 4)//QUAD
      {
      #pragma omp parallel default(shared)
        {
        #pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          if (indic[et]>-1)
          {
            E_Int no = npts+indic[et]*nptsAdd;
            E_Int eto = nelts+indic[et]*neltsAdd;
            E_Int indTab[9];// QUAD
            E_Int indA = cn1[et]-1;
            E_Int indB = cn2[et]-1;
            E_Int indC = cn3[et]-1;
            E_Int indD = cn4[et]-1;
            E_Float xA = xt[indA]; E_Float yA = yt[indA];
            E_Float xB = xt[indB];
            E_Float yC = yt[indC];
            E_Float yD = yt[indD];
            E_Float xmean = (xB + xA)*0.5;
            E_Float ymean = (yD + yA)*0.5;
            E_Float z = zt[indA];

            indTab[0]=indA;
            indTab[1]=no;  
            indTab[2]=indB;
            indTab[3]=no+1;
            indTab[4]=no+2;
            indTab[5]=no+3;  
            indTab[6]=indD;
            indTab[7]=no+4;
            indTab[8]=indC;

            xt[indTab[0]] = xA;    yt[indTab[0]] = yA;    zt[indTab[0]] = z;
            xt[indTab[1]] = xmean; yt[indTab[1]] = yA;    zt[indTab[1]] = z;
            xt[indTab[2]] = xB;    yt[indTab[2]] = yA;    zt[indTab[2]] = z;
            xt[indTab[3]] = xA;    yt[indTab[3]] = ymean; zt[indTab[3]] = z;
            xt[indTab[4]] = xmean; yt[indTab[4]] = ymean; zt[indTab[4]] = z;
            xt[indTab[5]] = xB   ; yt[indTab[5]] = ymean; zt[indTab[5]] = z;
            xt[indTab[6]] = xA;    yt[indTab[6]] = yC;    zt[indTab[6]] = z;
            xt[indTab[7]] = xmean; yt[indTab[7]] = yC;    zt[indTab[7]] = z;
            xt[indTab[8]] = xB   ; yt[indTab[8]] = yC;    zt[indTab[8]] = z;

            // 1er quad - element et dont on remplace les sommets
            E_Int indAl = indTab[0];
            E_Int indBl = indTab[1];
            E_Int indCl = indTab[4];
            E_Int indDl = indTab[3];
            cn1[et] = indAl+1; cn2[et] = indBl+1;
            cn3[et] = indCl+1; cn4[et] = indDl+1;
        
            // 2eme quad
            indAl = indTab[1];
            indBl = indTab[2];
            indCl = indTab[5];
            indDl = indTab[4];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            eto++;

            // 3eme quad
            indAl = indTab[4];
            indBl = indTab[5];
            indCl = indTab[8];
            indDl = indTab[7];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            eto++;

            // 4eme quad
            indAl = indTab[3];
            indBl = indTab[4];
            indCl = indTab[7];
            indDl = indTab[6];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            eto++;

          }//indic=1
        }// for elts
      }
      }// QUAD
      else //HEXA
      {
        E_Int* cn5 = cn.begin(5);
        E_Int* cn6 = cn.begin(6);
        E_Int* cn7 = cn.begin(7);
        E_Int* cn8 = cn.begin(8);
      #pragma omp parallel default(shared)
        {
        #pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          if (indic[et]>-1)
          {
            E_Int no = npts+indic[et]*nptsAdd;
            E_Int eto = nelts+indic[et]*neltsAdd;
            E_Int indTab[27];
            E_Int indA = cn1[et]-1; E_Int indB = cn2[et]-1; 
            E_Int indC = cn3[et]-1; E_Int indD = cn4[et]-1;
            E_Int indE = cn5[et]-1; E_Int indF = cn6[et]-1; 
            E_Int indG = cn7[et]-1; E_Int indH = cn8[et]-1;
            E_Float xA = xt[indA]; E_Float yA = yt[indA]; E_Float zA = zt[indA];
            E_Float xB = xt[indB]; 
            E_Float yC = yt[indC];
            E_Float zE = zt[indE];
            E_Float xmean = (xt[indB] + xt[indA])*0.5;
            E_Float ymean = (yt[indD] + yt[indA])*0.5;
            E_Float zmean = (zt[indE] + zt[indA])*0.5;

            indTab[0]=indA;
            indTab[1]=no;
            indTab[2]=indB;
            indTab[3]=no+1;
            indTab[4]=no+2;
            indTab[5]=no+3;  
            indTab[6]=indD;
            indTab[7]=no+4;
            indTab[8]=indC;// fin z = zA
            xt[indTab[0]] = xA;    yt[indTab[0]] = yA;    zt[indTab[0]] = zA;
            xt[indTab[1]] = xmean; yt[indTab[1]] = yA;    zt[indTab[1]] = zA;
            xt[indTab[2]] = xB;    yt[indTab[2]] = yA;    zt[indTab[2]] = zA;
            xt[indTab[3]] = xA;    yt[indTab[3]] = ymean; zt[indTab[3]] = zA;
            xt[indTab[4]] = xmean; yt[indTab[4]] = ymean; zt[indTab[4]] = zA;
            xt[indTab[5]] = xB   ; yt[indTab[5]] = ymean; zt[indTab[5]] = zA;
            xt[indTab[6]] = xA;    yt[indTab[6]] = yC;    zt[indTab[6]] = zA;
            xt[indTab[7]] = xmean; yt[indTab[7]] = yC;    zt[indTab[7]] = zA;
            xt[indTab[8]] = xB   ; yt[indTab[8]] = yC;    zt[indTab[8]] = zA;

            indTab[9]=no+5;
            indTab[10]=no+6;
            indTab[11]=no+7;
            indTab[12]=no+8;
            indTab[13]=no+9;
            indTab[14]=no+10;  
            indTab[15]=no+11;
            indTab[16]=no+12;
            indTab[17]=no+13;//z=zmean
            xt[indTab[9]] = xA;     yt[indTab[9]] = yA;     zt[indTab[9]] = zmean;
            xt[indTab[10]] = xmean; yt[indTab[10]] = yA;    zt[indTab[10]] = zmean;
            xt[indTab[11]] = xB;    yt[indTab[11]] = yA;    zt[indTab[11]] = zmean;
            xt[indTab[12]] = xA;    yt[indTab[12]] = ymean; zt[indTab[12]] = zmean;
            xt[indTab[13]] = xmean; yt[indTab[13]] = ymean; zt[indTab[13]] = zmean;
            xt[indTab[14]] = xB   ; yt[indTab[14]] = ymean; zt[indTab[14]] = zmean;
            xt[indTab[15]] = xA;    yt[indTab[15]] = yC;    zt[indTab[15]] = zmean;
            xt[indTab[16]] = xmean; yt[indTab[16]] = yC;    zt[indTab[16]] = zmean;
            xt[indTab[17]] = xB   ; yt[indTab[17]] = yC;    zt[indTab[17]] = zmean;


            indTab[18]=indE;
            indTab[19]=no+14;
            indTab[20]=indF;
            indTab[21]=no+15;
            indTab[22]=no+16;
            indTab[23]=no+17;  
            indTab[24]=indH;
            indTab[25]=no+18;
            indTab[26]=indG;// fin z = zE
            xt[indTab[18]] = xA;    yt[indTab[18]] = yA;    zt[indTab[18]] = zE;
            xt[indTab[19]] = xmean; yt[indTab[19]] = yA;    zt[indTab[19]] = zE;
            xt[indTab[20]] = xB;    yt[indTab[20]] = yA;    zt[indTab[20]] = zE;
            xt[indTab[21]] = xA;    yt[indTab[21]] = ymean; zt[indTab[21]] = zE;
            xt[indTab[22]] = xmean; yt[indTab[22]] = ymean; zt[indTab[22]] = zE;
            xt[indTab[23]] = xB   ; yt[indTab[23]] = ymean; zt[indTab[23]] = zE;
            xt[indTab[24]] = xA;    yt[indTab[24]] = yC;    zt[indTab[24]] = zE;
            xt[indTab[25]] = xmean; yt[indTab[25]] = yC;    zt[indTab[25]] = zE;
            xt[indTab[26]] = xB   ; yt[indTab[26]] = yC;    zt[indTab[26]] = zE;

            // 1er HEXA 0143-9101312
            E_Int indAl = indTab[0]; E_Int indEl = indTab[9]; 
            E_Int indBl = indTab[1]; E_Int indFl = indTab[10];
            E_Int indCl = indTab[4]; E_Int indGl = indTab[13];
            E_Int indDl = indTab[3]; E_Int indHl = indTab[12];
            cn1[et] = indAl+1; cn2[et] = indBl+1;
            cn3[et] = indCl+1; cn4[et] = indDl+1;
            cn5[et] = indEl+1; cn6[et] = indFl+1;
            cn7[et] = indGl+1; cn8[et] = indHl+1;


            // 2eme HEXA  1254-10111413
            indAl = indTab[1]; indBl = indTab[2]; indCl = indTab[5]; indDl = indTab[4];
            indEl = indTab[10]; indFl = indTab[11]; indGl = indTab[14]; indHl = indTab[13];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;

             // 3eme HEXA 4587-13451716
            indAl = indTab[4]; indBl = indTab[5]; indCl = indTab[8]; indDl = indTab[7];
            indEl = indTab[13]; indFl = indTab[14]; indGl = indTab[17]; indHl = indTab[16];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;

    
            // 4eme HEXA 3476-12131615
            indAl = indTab[3]; indBl = indTab[4]; indCl = indTab[7]; indDl = indTab[6];
            indEl = indTab[12]; indFl = indTab[13]; indGl = indTab[16]; indHl = indTab[15];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;
       
            // 5e HEXA 9101312-18-192221
            indAl = indTab[9]; indBl = indTab[10]; indCl = indTab[13]; indDl = indTab[12];
            indEl = indTab[18]; indFl = indTab[19]; indGl = indTab[22]; indHl = indTab[21];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;

            //6e HEXA 10-11-14-13-19-20-23-22
            indAl = indTab[10]; indBl = indTab[11]; indCl = indTab[14]; indDl = indTab[13];
            indEl = indTab[19]; indFl = indTab[20]; indGl = indTab[23]; indHl = indTab[22];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;

            //7e HEXA 12-13-16-15-21-22-25-24
            indAl = indTab[12]; indBl = indTab[13]; indCl = indTab[16]; indDl = indTab[15];
            indEl = indTab[21]; indFl = indTab[22]; indGl = indTab[25]; indHl = indTab[24];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;
            //8e HEXA 13-14-17-16-22-23-26-25
            indAl = indTab[13]; indBl = indTab[14]; indCl = indTab[17]; indDl = indTab[16];
            indEl = indTab[22]; indFl = indTab[23]; indGl = indTab[26]; indHl = indTab[25];
            cn1[eto] = indAl+1; cn2[eto] = indBl+1;
            cn3[eto] = indCl+1; cn4[eto] = indDl+1;
            cn5[eto] = indEl+1; cn6[eto] = indFl+1;
            cn7[eto] = indGl+1; cn8[eto] = indHl+1;
            eto++;


            no+=nptsAdd;
          }//indic=1
        }// for elts        
      }//omp
      }//HEXA
    }//elts2Raff > 0
  }
  
  return;
}
//=============================================================================
/* Verifie l'equilibrage du maillage engendre. Cas octree 2. */
//=============================================================================
void K_GENERATOR::checkBalancing2(FldArrayI& cn, FldArrayF& coords)
{
  //E_Float eps = 1.e-10;
  E_Int ind, ind1, ind2;
  FldArrayI cno = cn; FldArrayF fo = coords; // copies
  E_Int nvert = cn.getNfld();
  E_Int npts = coords.getSize();
  E_Int nelts = cn.getSize();
  E_Int no = 0; E_Int eto = 0;
  E_Int et2, nvoisins;
  E_Float dh1, dh2, dh1s2;
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
  E_Float* xt = coords.begin(1); E_Float* yt = coords.begin(2); E_Float* zt = coords.begin(3);
  E_Int splitl = 1;
  const char* eltType;
  if (nvert == 4) eltType = "QUAD";
  else eltType = "HEXA";
  while (splitl == 1)
  {
    npts = coords.getSize(); nelts = cn.getSize();
    xt = coords.begin(1); yt = coords.begin(2); zt = coords.begin(3);
    FldArrayF indict(nelts); indict.setAllValuesAtNull();
    FldArrayI indir(npts); indir.setAllValuesAt(-1);
    E_Int* indirp = indir.begin();
    //calcul de dh
    FldArrayF dht(nelts);
    E_Float* dhtp = dht.begin();
    for (E_Int et = 0; et < nelts; et++)
    {
      ind1 = cn1[et]-1; ind2 = cn2[et]-1;
      dhtp[et] = xt[ind2]-xt[ind1];
    }  
    // determination des elements voisins
    vector< vector<E_Int> > cEEN(nelts);
    E_Int ok = getNeighbourElts(npts, xt, yt, zt, cn, cEEN); 
    if (ok == -1) return;

    eto = 0; no = 0; splitl = 0;
    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      dh1 = dhtp[et1]; dh1s2 = 0.5*dh1;
      vector<E_Int>& voisins = cEEN[et1]; nvoisins = voisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        et2 = voisins[noet2]; dh2 = dhtp[et2];
        if (dh1s2 > 1.1*dh2)
        {
          splitElement(et1, npts, nelts, indict[et1], indirp, cn, xt, yt, zt, 
                       cno, fo, indict, no, eto);
          splitl = 1; goto nextet1;
        }
      }
      // construction dans cno, fo
      for (E_Int nov = 1; nov <= nvert; nov++)
      {
        ind = cn(et1, nov)-1;
        if (indirp[ind] == -1) 
        {
          fo(no,1) = xt[ind]; fo(no,2) = yt[ind]; fo(no,3) = zt[ind];
          no++; 
          cno(eto, nov) = no; indirp[ind] = no;
          if (no+10 > fo.getSize()) fo.reAllocMat(fo.getSize()+npts, 3);
        }
        else cno(eto,nov) = indirp[ind];
      }
      eto++;
      if (eto+10 > indict.getSize()) indict.reAlloc(indict.getSize()+nelts);
      if (eto+10 > cno.getSize()) cno.reAllocMat(cno.getSize()+nelts, nvert);

      nextet1:;
    }//fin parcours et1

    // Realloc
    if (splitl == 1)
    {
      fo.reAllocMat(no,3); cno.reAllocMat(eto,nvert);
      cn = cno; coords = fo; cn1 = cn.begin(1); cn2 = cn.begin(2);
    }
  }
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, eltType, coords, cn);
  return;
}
//=============================================================================
/* Verifie l'equilibrage du maillage engendre. Cas d'un 9/27-tree */
//=============================================================================
void K_GENERATOR::checkBalancing3(FldArrayI& cn, FldArrayF& coords)
{
  E_Int ind, ind1, ind2;
  E_Float eps = 1.e-10;
  E_Int split = 1;

  while (split == 1) 
  {
    FldArrayI cno = cn; FldArrayF fo = coords; // copies
    E_Int npts = coords.getSize(); E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();
    FldArrayF indict(nelts); indict.setAllValuesAtNull();
    //calcul de dh
    FldArrayF dht(nelts);
    E_Float* dhtp = dht.begin();
    E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
    E_Float* xt = coords.begin(1); E_Float* yt = coords.begin(2); E_Float* zt = coords.begin(3);
    for (E_Int et = 0; et < nelts; et++)
    {
      ind1 = cn1[et]-1; ind2 = cn2[et]-1;
      dhtp[et] = xt[ind2]-xt[ind1];
    }  
    
    // determination des elements voisins
    vector< vector<E_Int> > cEEN(nelts);
    E_Int ok = getNeighbourElts(npts, xt, yt, zt, cn, cEEN); 
    if (ok == -1) return;
    E_Int no = 0; E_Int eto = 0;
    E_Int et2, nvoisins;
    E_Float dh1, dh2;
    split = 0;

    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      dh1 = dhtp[et1];
      vector<E_Int>& voisins = cEEN[et1]; nvoisins = voisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        et2 = voisins[noet2]; dh2 = dhtp[et2];
        if (dh1 > 3.*dh2 + eps)
        {
          splitElement27(et1, npts, indict[et1], cn, xt, yt, zt, cno, fo, 
                         indict, no, eto);
          split = 1; goto nextet1;
        }
      }
      // construction dans cno, fo
       for (E_Int nov = 1; nov <= nvert; nov++)
       {
         ind = cn(et1, nov)-1;
         fo(no,1) = xt[ind]; fo(no,2) = yt[ind]; fo(no,3) = zt[ind];
         no++; cno(eto, nov) = no; 
         if (no+10 > fo.getSize()) fo.reAllocMat(fo.getSize()+npts,3);
       }
       eto++;
       if (eto + 28 > indict.getSize()) indict.reAllocMat(indict.getSize()+nelts,1);
       if (eto + 28 > cno.getSize()) cno.reAllocMat(cno.getSize()+nelts,nvert); 
      
      nextet1:;
    }
    // Nettoyage
    for (E_Int v = 0; v < nelts; v++) cEEN[v].clear();
    cEEN.clear();
    // Realloc
    fo.reAllocMat(no,3); cno.reAllocMat(eto,nvert); indict.reAllocMat(eto,1);
    if (nvert == 4) K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, "QUAD", 
                                                 fo, cno);
    else K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, "HEXA", fo, cno);
    coords = fo; cn = cno;
  }
  return;
}
