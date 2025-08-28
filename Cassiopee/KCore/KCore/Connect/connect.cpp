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

#include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Renvoi le point image du point vert0 pour des elements QUAD/PENTA/HEXA
   point image : point "en face" de  vert0 dans la face opposee
   IN: cNG: connectivite NGON: Faces/Noeuds et Elts/Faces
   IN: et0: numero de l'element considere dans cFNEF: demarre a 0 
   IN: vert0: pt de et0, de face0 dont on cherche l'image ds et0: demarre a 1
   IN: face0: numero de la face a laquelle appartient le pt vert0: demarre a 1
   IN: vertices: liste des vertices candidats de la face 
   OUT: vert1: pt image: demarre a 1
*/
//=============================================================================
E_Int K_CONNECT::image(E_Int vert0, E_Int face0, E_Int et0, 
                       vector<E_Int>& vertices, FldArrayI& cNG)
{
  E_Int* ngon = cNG.getNGon(); E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG(); E_Int* indPH = cNG.getIndPH(); 
  return K_CONNECT::image(vert0, face0, et0, vertices, cNG,
                          ngon, nface, indPG, indPH);
}                     
                       
E_Int K_CONNECT::image(E_Int vert0, E_Int face0, E_Int et0, 
                       vector<E_Int>& vertices, FldArrayI& cNG,
                       E_Int* ngon, E_Int* nface, E_Int* indPG, E_Int* indPH)
{
  E_Int vert, vertm, vertp, nv, nfacesl;
  E_Int* elem0 = cNG.getElt(et0, nfacesl, nface, indPH); // pointeur sur l'element et0 dans cNG
  // selectionne les faces de l'element autre que face0
  vector<E_Int> faces;
  
  for (E_Int i = 0; i < nfacesl; i++)
  {
    if (elem0[i] != face0) faces.push_back(elem0[i]-1);
  }
  
  for (size_t nof = 0; nof < faces.size(); nof++)
  {
    E_Int* face = cNG.getFace(faces[nof], nv, ngon, indPG);
    // recherche du pt vert0
    for (E_Int nov = 0; nov < nv; nov++)
    {
      vert = face[nov]; // indice du point
      vertm = -1; vertp = -1; // vertm : indice du point precedent, vertp : indice du point suivant
      if (vert == vert0)
      {
        if (nov == 0) vertm = face[nv-1];
        else vertm = face[nov-1];
        if (nov == nv-1) vertp = face[0];
        else vertp = face[nov+1];
        goto fin;
      }
    }
  }
  fin:;
  // recherche du bon candidat: 
  for (size_t nov = 0; nov < vertices.size(); nov++)
  {
    if (vertices[nov] == vertm) return vertm;
    else if (vertices[nov] == vertp) return vertp;
  }
  return -1;
}

//=============================================================================
/* Detection des interfaces coincidentes
   IN: im1,jm1,km1: dimensions de l'array 1
   IN: f1: champs de 1: coordonnees incluses
   IN: posx1, posy1, posz1: positions des coordonnees ds f1
   IN: im2,jm2,km2: dimensions de l array 2
   IN: f2: champs de 2: coordonnees incluses
   IN: posx2, posy2, posz2: positions des coordonnees ds f2
   OUT: nof1: numero de la face de 1 coincidente avec une face de 2 
   OUT: nof2: numero de la face de 2 coincidente avec une face de 1  */
//=============================================================================
E_Int K_CONNECT::detectMatchInterface(E_Int im1, E_Int jm1, E_Int km1, 
                                      E_Int posx1, E_Int posy1, E_Int posz1,
                                      E_Int im2, E_Int jm2, E_Int km2,
                                      E_Int posx2, E_Int posy2, E_Int posz2,
                                      FldArrayF& f1, FldArrayF& f2,
                                      E_Int& nof1, E_Int& nof2,
                                      E_Float eps)
{
  // La fonction qui suit est moins rapide, mais fonctionne 
  // pour 2 couches
  if (im1 == 2 || jm1 == 2 || km1 == 2) 
    return detectMatchInterface2(im1, jm1, km1, 
                                 posx1, posy1, posz1,
                                 im2, jm2, km2,
                                 posx2, posy2, posz2,
                                 f1, f2,
                                 nof1, nof2,
                                 eps);

  // Face array
  E_Int face1[36]; E_Int face2[36];
  for (E_Int i = 0; i < 36; i++) face1[i] = 1;
  for (E_Int i = 0; i < 36; i++) face2[i] = 1;
  face1[0] = im1; face1[6] = im1; face1[18] = jm1; face1[30] = km1;
  face1[19] = jm1; face1[31] = km1;
  face1[8] = im1; face1[14] = jm1; face1[20] = jm1; face1[32] = km1;
  face1[9] = im1; face1[33] = km1;
  face1[10] = im1; face1[22] = jm1; face1[28] = km1; face1[34] = km1;
  face1[11] = im1; face1[23] = jm1;
  
  face2[18] = jm2; face2[30] = km2;
  face2[1] = im2; face2[7] = im2; face2[19] = jm2; face2[31] = km2;
  face2[8] = im2; face2[32] = km2;
  face2[9] = im2; face2[15] = jm2; face2[21] = jm2; face2[33] = km2;
  face2[10] = im2; face2[22] = jm2;
  face2[11] = im2; face2[23] = jm2; face2[29] = km2; face2[35] = km2;

  /* Reference (a conserver)
  FldArrayI face1(6,6); face1.setAllValuesAt(1);
  face1(0,1) = im1; face1(0,2) = im1; face1(0,4) = jm1; face1(0,6) = km1;
  face1(1,4) = jm1; face1(1,6) = km1;
  face1(2,2) = im1; face1(2,3) = jm1; face1(2,4) = jm1; face1(2,6) = km1;
  face1(3,2) = im1; face1(3,6) = km1;
  face1(4,2) = im1; face1(4,4) = jm1; face1(4,5) = km1; face1(4,6) = km1;
  face1(5,2) = im1; face1(5,4) = jm1;
  
  FldArrayI face2(6,6); face2.setAllValuesAt(1);
  face2(0,4) = jm2; face2(0,6) = km2;
  face2(1,1) = im2; face2(1,2) = im2; face2(1,4) = jm2; face2(1,6) = km2;
  face2(2,2) = im2; face2(2,6) = km2;
  face2(3,2) = im2; face2(3,3) = jm2; face2(3,4) = jm2; face2(3,6) = km2;
  face2(4,2) = im2; face2(4,4) = jm2;
  face2(5,2) = im2; face2(5,4) = jm2; face2(5,5) = km2; face2(5,6) = km2;
  */

  // Detection d'une face de 1 coincidente avec une face de 2
  E_Int a, b, ab, ind1, ind2;
  E_Int im1jm1 = im1*jm1;
  E_Int im2jm2 = im2*jm2;
  
  E_Float* xt1 = f1.begin(posx1);
  E_Float* yt1 = f1.begin(posy1);
  E_Float* zt1 = f1.begin(posz1);
  E_Float* xt2 = f2.begin(posx2);
  E_Float* yt2 = f2.begin(posy2);
  E_Float* zt2 = f2.begin(posz2);

  E_Bool match;
  nof1 = -1; nof2 = -1;
  E_Float xpt1, xpt2, ypt1, ypt2, zpt1, zpt2, d;
  E_Int ak, aj, ai;
  E_Int bk, bj, bi;
  E_Int deltaai, deltaaj, deltaak, deltabi, deltabj, deltabk;
  E_Int aip,ajp,akp;
  E_Int bip,bjp,bkp;
  E_Int orig[3*36];
  E_Int deltabai[3*36];
  E_Int deltabaj[3*36];
  E_Int deltabak[3*36];
  E_Int check[36]; // possible match
  for (E_Int i = 0; i < 36; i++) check[i] = 0;

  // Fast corner check - MT
  for (ab = 0; ab < 36; ab++)
  {
    a = ab/6; b = ab-6*a;

    match = true;
    deltaak = K_FUNC::E_max(face1[a+30]-face1[a+24], 1);
    deltaaj = K_FUNC::E_max(face1[a+18]-face1[a+12], 1);
    deltaai = K_FUNC::E_max(face1[a+6]-face1[a], 1);
    deltabk = K_FUNC::E_max(face2[b+30]-face2[b+24], 1);
    deltabj = K_FUNC::E_max(face2[b+18]-face2[b+12], 1);
    deltabi = K_FUNC::E_max(face2[b+6]-face2[b], 1);

    orig[0+3*ab] = 0; orig[1+3*ab] = 0; orig[2+3*ab] = 0; 
    deltabai[0+3*ab] = face2[b]; deltabai[1+3*ab] = face2[b+12]; deltabai[2+3*ab] = face2[b+24];
    deltabaj[0+3*ab] = face2[b]; deltabaj[1+3*ab] = face2[b+12]; deltabaj[2+3*ab] = face2[b+24];
    deltabak[0+3*ab] = face2[b]; deltabak[1+3*ab] = face2[b+12]; deltabak[2+3*ab] = face2[b+24];

    //printf("face1 %d %d | %d %d | %d %d\n", face1[a], face1[a+6],
    //       face1[a+12], face1[a+18], face1[a+24], face1[a+30]);
    //printf("face2 %d %d | %d %d | %d %d\n", face2[b], face2[b+6],
    //       face2[b+12], face2[b+18], face2[b+24], face2[b+30]);

    for (ak = face1[a+24]; ak <= face1[a+30]; ak += deltaak)
      for (aj = face1[a+12]; aj <= face1[a+18]; aj += deltaaj)
        for (ai = face1[a]; ai <= face1[a+6]; ai += deltaai)
        {
          ind1 = (ai-1) + (aj-1)*im1 + (ak-1)*im1jm1;
          xpt1 = xt1[ind1]; ypt1 = yt1[ind1]; zpt1 = zt1[ind1];

          akp = ak+1; ajp = aj+1; aip = ai+1;          
          if (ak == face1[a+24]) akp = std::min(ak+1, face1[a+30]);
          //if (km1 == 2 && ak == face1[a+24]) akp = ak;
          if (aj == face1[a+12]) ajp = std::min(aj+1, face1[a+18]);
          //if (jm1 == 2 && aj == face1[a+12]) ajp = aj;
          if (ai == face1[a]) aip = std::min(ai+1, face1[a+6]);
          //if (im1 == 2 && ai == face1[a]) aip = ai;
          if (ak == face1[a+30]) akp = std::max(ak-1, face1[a+24]);
          //if (km1 == 2 && ak == face1[a+30]) akp = ak;
          if (aj == face1[a+18]) ajp = std::max(aj-1, face1[a+12]);
          //if (jm1 == 2 && aj == face1[a+18]) ajp = aj;
          if (ai == face1[a+6]) aip = std::max(ai-1, face1[a]);
          //if (im1 == 2 && ai == face1[a]) aip = ai;
          
          ind1 = (aip-1) + (aj-1)*im1 + (ak-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          ind1 = (ai-1) + (ajp-1)*im1 + (ak-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          ind1 = (ai-1) + (aj-1)*im1 + (akp-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          ind1 = (aip-1) + (ajp-1)*im1 + (ak-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          ind1 = (aip-1) + (aj-1)*im1 + (akp-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          ind1 = (ai-1) + (ajp-1)*im1 + (akp-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          ind1 = (aip-1) + (ajp-1)*im1 + (akp-1)*im1jm1;
          xpt1 += xt1[ind1]; ypt1 += yt1[ind1]; zpt1 += zt1[ind1];
          xpt1 = xpt1/8.; ypt1 = ypt1/8.; zpt1 = zpt1/8.;

          //printf("herea %d %d %d -> %d %d %d\n",ai,aj,ak, face1[a+6],face1[a+18],face1[a+30]);

          for (bk = face2[b+24]; bk <= face2[b+30]; bk += deltabk)
            for (bj = face2[b+12]; bj <= face2[b+18]; bj += deltabj)
              for (bi = face2[b]; bi <= face2[b+6]; bi += deltabi)
              {
                ind2 = (bi-1) + (bj-1)*im2 + (bk-1)*im2jm2;
                xpt2 = xt2[ind2]; ypt2 = yt2[ind2]; zpt2 = zt2[ind2];

                bkp = bk+1; bjp = bj+1; bip = bi+1;
                if (bk == face2[b+24]) bkp = std::min(bk+1, face2[b+30]);
                //if (km2 == 2 && bk == face2[b+24]) bkp = bk;
                if (bj == face2[b+12]) bjp = std::min(bj+1, face2[b+18]);
                //if (jm2 == 2 && bj == face2[b+12]) bjp = bj;
                if (bi == face2[b]) bip = std::min(bi+1, face2[b+6]);
                //if (im2 == 2 && bi == face2[b]) bip = bi;
                if (bk == face2[b+30]) bkp = std::max(bk-1, face2[b+24]);
                //if (km2 == 2 && bk == face2[b+30]) bkp = bk;
                if (bj == face2[b+18]) bjp = std::max(bj-1, face2[b+12]);
                //if (jm2 == 2 && bj == face2[b+18]) bjp = bj;
                if (bi == face2[b+6]) bip = std::max(bi-1, face2[b]);
                //if (im2 == 2 && bi == face2[b+6]) bip = bi;
                
                ind2 = (bip-1) + (bj-1)*im2 + (bk-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                ind2 = (bi-1) + (bjp-1)*im2 + (bk-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                ind2 = (bi-1) + (bj-1)*im2 + (bkp-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                ind2 = (bip-1) + (bjp-1)*im2 + (bk-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                ind2 = (bip-1) + (bj-1)*im2 + (bkp-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                ind2 = (bi-1) + (bjp-1)*im2 + (bkp-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                ind2 = (bip-1) + (bjp-1)*im2 + (bkp-1)*im2jm2;
                xpt2 += xt2[ind2]; ypt2 += yt2[ind2]; zpt2 += zt2[ind2];
                xpt2 = xpt2/8.; ypt2 = ypt2/8.; zpt2 = zpt2/8.;

                d = (xpt1-xpt2)*(xpt1-xpt2) +
                  (ypt1-ypt2)*(ypt1-ypt2) +
                  (zpt1-zpt2)*(zpt1-zpt2);
                
                if (K_FUNC::fEqualZero(d, eps) == true)
                {
                  if (ak == face1[a+24] && aj == face1[a+12] && ai == face1[a])
                  { //printf("hereo %d %d %d\n",bi,bj,bk);
                    orig[0+3*ab]=bi; orig[1+3*ab]=bj; orig[2+3*ab]=bk; }
                  else if (ak == face1[a+24] && aj == face1[a+12] && ai == face1[a+6])
                  { //printf("here1 %d %d %d\n",bi,bj,bk);
                    deltabai[0+3*ab]=bi; deltabai[1+3*ab]=bj; deltabai[2+3*ab]=bk; }
                  else if (ak == face1[a+24] && aj == face1[a+18] && ai == face1[a])
                  { //printf("here2 %d %d %d\n",bi,bj,bk);
                    deltabaj[0+3*ab]=bi; deltabaj[1+3*ab]=bj; deltabaj[2+3*ab]=bk; }
                  else if (ak == face1[a+30] && aj == face1[a+12] && ai == face1[a])
                  { //printf("here3 %d %d %d\n",bi,bj,bk);
                    deltabak[0+3*ab]=bi; deltabak[1+3*ab]=bj; deltabak[2+3*ab]=bk; }
                  goto nextp1; //tester le point de la face 1 suivant
                }
              } //parcours des elts de face 2
          
          // un pt de face 2 non coincident avec face 1
          // alors faces (a,b) non coincidentes
          match = false;
          nextp1:;
        } // parcours des elts de face 1
    if (match == true) check[ab] = 1;
  }
  E_Int sum = 0;
  for (E_Int i = 0; i < 36; i++) sum += check[i];
  if (sum == 0) return 0; // no possible check

  // Full test
  for (ab = 0; ab < 36; ab++)
  {
    if (check[ab] == 1) // possible match
    {
      a = ab/6; b = ab-6*a;

      //printf("delta\n");
      //printf("%d %d %d\n", orig[0+3*ab], orig[1+3*ab], orig[2+3*ab]);
      //printf("%d %d %d\n", deltabai[0+3*ab], deltabai[1+3*ab], deltabai[2+3*ab]);
      //printf("%d %d %d\n", deltabaj[0+3*ab], deltabaj[1+3*ab], deltabaj[2+3*ab]);
      //printf("%d %d %d\n", deltabak[0+3*ab], deltabak[1+3*ab], deltabak[2+3*ab]);

      match = true;

      for (ak = face1[a+24]; ak <= face1[a+30]; ak++)
        for (aj = face1[a+12]; aj <= face1[a+18]; aj++)
          for (ai = face1[a]; ai <= face1[a+6]; ai++)
          {
            ind1 = (ai-1) + (aj-1)*im1 + (ak-1)*im1jm1;
            xpt1 = xt1[ind1]; ypt1 = yt1[ind1]; zpt1 = zt1[ind1];
          
            bi = orig[0+3*ab]+
            K_FUNC::E_signum(deltabak[0+3*ab]-orig[0+3*ab])*(ak-face1[a+24])+
            K_FUNC::E_signum(deltabaj[0+3*ab]-orig[0+3*ab])*(aj-face1[a+12])+
            K_FUNC::E_signum(deltabai[0+3*ab]-orig[0+3*ab])*(ai-face1[a]);
            bj = orig[1+3*ab]+
            K_FUNC::E_signum(deltabak[1+3*ab]-orig[1+3*ab])*(ak-face1[a+24])+
            K_FUNC::E_signum(deltabaj[1+3*ab]-orig[1+3*ab])*(aj-face1[a+12])+
            K_FUNC::E_signum(deltabai[1+3*ab]-orig[1+3*ab])*(ai-face1[a]);
            bk = orig[2+3*ab]+
            K_FUNC::E_signum(deltabak[2+3*ab]-orig[2+3*ab])*(ak-face1[a+24])+
            K_FUNC::E_signum(deltabaj[2+3*ab]-orig[2+3*ab])*(aj-face1[a+12])+
            K_FUNC::E_signum(deltabai[2+3*ab]-orig[2+3*ab])*(ai-face1[a]);
          
            //printf("here: %d %d %d -> %d %d %d (%d %d %d)\n",ai,aj,ak,bi,bj,
            //       bk,deltabak[2+3*ab],deltabaj[2+3*ab],deltabai[2+3*ab]);

            ind2 = (bi-1) + (bj-1)*im2 + (bk-1)*im2jm2;
            xpt2 = xt2[ind2]; ypt2 = yt2[ind2]; zpt2 = zt2[ind2];
            d = (xpt1-xpt2)*(xpt1-xpt2) +
              (ypt1-ypt2)*(ypt1-ypt2) +
              (zpt1-zpt2)*(zpt1-zpt2);
                
            if (K_FUNC::fEqualZero(d, eps) != true)
            {
                match = false;
                goto next1; //ne match pas
            }
          } //parcours des elts de face 1
      
      if (match == true) // full match
      {
        nof1 = a;
        if (nof1 == 0 || nof1 == 2 || nof1 == 4) nof1 = nof1+2;
        nof2 = b+1;
        return 1;
      }
    } // check
    
    next1: ;

  } // parcours de faces a et b
  return 0;
}
//=============================================================================
/* Detection des interfaces coincidentes
   IN: im1,jm1,km1: dimensions de l'array 1
   IN: f1: champs de 1: coordonnees incluses
   IN: posx1, posy1, posz1: positions des coordonnees ds f1
   IN: im2,jm2,km2: dimensions de l array 2
   IN: f2: champs de 2: coordonnees incluses
   IN: posx2, posy2, posz2: positions des coordonnees ds f2
   OUT: nof1: numero de la face de 1 coincidente avec une face de 2
   OUT: nof2: numero de la face de 2 coincidente avec une face de 1  */
//=============================================================================
E_Int K_CONNECT::detectMatchInterface2(E_Int im1, E_Int jm1, E_Int km1, 
                                      E_Int posx1, E_Int posy1, E_Int posz1,
                                      E_Int im2, E_Int jm2, E_Int km2,
                                      E_Int posx2, E_Int posy2, E_Int posz2,
                                      FldArrayF& f1, FldArrayF& f2,
                                      E_Int& nof1, E_Int& nof2,
                                      E_Float eps)
{
  // Face array
  E_Int face1[36]; E_Int face2[36];
  for (E_Int i = 0; i < 36; i++) face1[i] = 1;
  for (E_Int i = 0; i < 36; i++) face2[i] = 1;
  face1[0] = im1; face1[6] = im1; face1[18] = jm1; face1[30] = km1;
  face1[19] = jm1; face1[31] = km1;
  face1[8] = im1; face1[14] = jm1; face1[20] = jm1; face1[32] = km1;
  face1[9] = im1; face1[33] = km1;
  face1[10] = im1; face1[22] = jm1; face1[28] = km1; face1[34] = km1;
  face1[11] = im1; face1[23] = jm1;
 
  face2[18] = jm2; face2[30] = km2;
  face2[1] = im2; face2[7] = im2; face2[19] = jm2; face2[31] = km2;
  face2[8] = im2; face2[32] = km2;
  face2[9] = im2; face2[15] = jm2; face2[21] = jm2; face2[33] = km2;
  face2[10] = im2; face2[22] = jm2;
  face2[11] = im2; face2[23] = jm2; face2[29] = km2; face2[35] = km2;

  // Detection d'une face de 1 coincidente avec une face de 2
  E_Int a, b, ab, ind1, ind2;
  E_Int im1jm1 = im1*jm1;
  E_Int im2jm2 = im2*jm2;
 
  E_Float* xt1 = f1.begin(posx1);
  E_Float* yt1 = f1.begin(posy1);
  E_Float* zt1 = f1.begin(posz1);
  E_Float* xt2 = f2.begin(posx2);
  E_Float* yt2 = f2.begin(posy2);
  E_Float* zt2 = f2.begin(posz2);

  E_Bool match;
  nof1 = -1; nof2 = -1;
  E_Float xpt1, xpt2, ypt1, ypt2, zpt1, zpt2, d;
  E_Int ak, aj, ai;
  E_Int bk, bj, bi;
  E_Int deltaai, deltaaj, deltaak, deltabi, deltabj, deltabk;

  // Fast corner check
  for (ab = 0; ab < 36; ab++)
  {
    a = ab/6; b = ab-6*a;

    match = true;
     
    deltaak = K_FUNC::E_max(face1[a+30]-face1[a+24], 1);
    deltaaj = K_FUNC::E_max(face1[a+18]-face1[a+12], 1);
    deltaai = K_FUNC::E_max(face1[a+6]-face1[a], 1);
    deltabk = K_FUNC::E_max(face2[b+30]-face2[b+24], 1);
    deltabj = K_FUNC::E_max(face2[b+18]-face2[b+12], 1);
    deltabi = K_FUNC::E_max(face2[b+6]-face2[b], 1);

    for (ak = face1[a+24]; ak <= face1[a+30]; ak += deltaak)
      for (aj = face1[a+12]; aj <= face1[a+18]; aj += deltaaj) 
        for (ai = face1[a]; ai <= face1[a+6]; ai += deltaai)
        {
          ind1 = (ai-1) + (aj-1)*im1 + (ak-1)*im1jm1;
          xpt1 = xt1[ind1]; ypt1 = yt1[ind1]; zpt1 = zt1[ind1];
         
          for (bk = face2[b+24]; bk <= face2[b+30]; bk += deltabk)
            for (bj = face2[b+12]; bj <= face2[b+18]; bj += deltabj)
              for (bi = face2[b]; bi <= face2[b+6]; bi += deltabi)
              {
                ind2 = (bi-1) + (bj-1)*im2 + (bk-1)*im2jm2;
                xpt2 = xt2[ind2]; ypt2 = yt2[ind2]; zpt2 = zt2[ind2]; 
                d = (xpt1-xpt2)*(xpt1-xpt2) +
                  (ypt1-ypt2)*(ypt1-ypt2) +
                  (zpt1-zpt2)*(zpt1-zpt2);
               
                if (K_FUNC::fEqualZero(d, eps) == true)
                {
                  goto nextp1; //tester le point de la face 1 suivant
                }
              }//parcours des elts de face 2
         
          // un pt de face 2 non coincident avec face 1
          // alors faces (a,b) non coincidentes
          match = false;
          nextp1:;
        }// parcours des elts de face 1
    if (match == true) goto fulltest;
  }
  if (match == false) return 0;
 
  // Full test
  fulltest:;
  for (ab = 0; ab < 36; ab++)
  {
    a = ab/6; b = ab-6*a;

    match = true;

    for (ak = face1[a+24]; ak <= face1[a+30]; ak++)
      for (aj = face1[a+12]; aj <= face1[a+18]; aj++)
        for (ai = face1[a]; ai <= face1[a+6]; ai++)
        {
          ind1 = (ai-1) + (aj-1)*im1 + (ak-1)*im1jm1;
          xpt1 = xt1[ind1]; ypt1 = yt1[ind1]; zpt1 = zt1[ind1];
         
          for (bk = face2[b+24]; bk <= face2[b+30]; bk++)
            for (bj = face2[b+12]; bj <= face2[b+18]; bj++)
              for (bi = face2[b]; bi <= face2[b+6]; bi++)
              {
                ind2 = (bi-1) + (bj-1)*im2 + (bk-1)*im2jm2;
                xpt2 = xt2[ind2]; ypt2 = yt2[ind2]; zpt2 = zt2[ind2]; 
                d = (xpt1-xpt2)*(xpt1-xpt2) +
                  (ypt1-ypt2)*(ypt1-ypt2) +
                  (zpt1-zpt2)*(zpt1-zpt2);
               
                if (K_FUNC::fEqualZero(d, eps) == true)
                {
                  goto next1; //tester le point de la face 1 suivant
                }
              }//parcours des elts de face 2
         
          // un pt de face 2 non coincident avec face 1
          // alors faces (a,b) non coincidentes
          match = false;
          next1:;
        }// parcours des elts de face 1

    if (match == true)
    {
      nof1 = a;
      if (nof1 == 0 || nof1 == 2 || nof1 == 4) nof1 = nof1+2;
     
      nof2 = b+1;
      goto end1; // trouve: fin parcours
    }
  } // parcours de faces a et b
  if (match == false) return 0;
  end1:;
  //
  return 1;
}
//=============================================================================
/* Determine si un des coins de f2 coincident avec le coin (i1,j1,k1) de f1 
   attention: les indices demarrent a 1 en i j et k 
   Si pas d'ambiguite, retourne 1 et les indices i2,j2,k2 correspondant
   Si ambiguite (plusieurs points possibles), retourne le nb d'ambiguites, et 
   retourne le dernier point coincident trouve
*/
//=============================================================================
E_Int K_CONNECT::getCoincidentVertexIndices(
  E_Int i1, E_Int j1, E_Int k1, 
  E_Int im1, E_Int jm1, E_Int km1, 
  E_Int im2, E_Int jm2, E_Int km2,
  E_Int posx1, E_Int posy1, E_Int posz1, 
  E_Int posx2, E_Int posy2, E_Int posz2,
  FldArrayF& f1, FldArrayF& f2, 
  E_Int& i2, E_Int& j2, E_Int& k2,
  E_Float eps)
{
  E_Int im1jm1 = im1*jm1;
  E_Int ind1 = (i1-1) + (j1-1) * im1 + (k1-1) * im1jm1;
  E_Float* xt1 = f1.begin(posx1);
  E_Float* xt2 = f2.begin(posx2);
  E_Float* yt1 = f1.begin(posy1);
  E_Float* yt2 = f2.begin(posy2);
  E_Float* zt1 = f1.begin(posz1);
  E_Float* zt2 = f2.begin(posz2);
  E_Float dx, dy, dz;
  E_Float xt11 = xt1[ind1];
  E_Float yt11 = yt1[ind1];
  E_Float zt11 = zt1[ind1];
  E_Float eps2 = 3*eps*eps;
  
  E_Int count = 0; E_Int isav=1, jsav=1, ksav=1;
  // test (1,1,1)
  i2 = 1; j2 = 1; k2 = 1;
  dx = xt2[0] - xt11;
  dy = yt2[0] - yt11;
  dz = zt2[0] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (im2,1,1)
  i2 = im2; j2 = 1; k2 = 1;
  E_Int ind2 = im2-1;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (1,jm2,1)
  i2 = 1; j2 = jm2; k2 = 1;
  ind2 = (jm2-1) * im2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (im2,jm2,1)
  i2 = im2; j2 = jm2; k2 = 1;
  ind2 = (im2-1) + (jm2-1) * im2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (1,1,km2)
  E_Int im2jm2 = im2*jm2;
  i2 = 1; j2 = 1; k2 = km2;
  ind2 = (km2-1)*im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11; 
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (im2,1,km2)
  i2 = im2; j2 = 1; k2 = km2;
  ind2 = im2-1 + (km2-1) * im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11; 
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (1,jm2,km2)
  i2 = 1; j2 = jm2; k2 = km2;
  ind2 = (jm2-1) * im2 + (km2-1) * im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11; 
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  // test (im2,jm2,km2)
  i2 = im2; j2 = jm2; k2 = km2;
  ind2 = (im2-1) + (jm2-1) * im2 + (km2-1) * im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11; 
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2; jsav=j2; ksav=k2; count++;}

  i2 = isav; j2 = jsav; k2 = ksav;
  return count;
}