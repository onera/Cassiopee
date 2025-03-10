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

# include "CompGeom/compGeom.h"
# include "Connect/connect.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
void K_COMPGEOM::compStructCurvatureHeight1D(
  E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt, E_Float* hmaxt)
{
  E_Float eps = 1.e-10;
  E_Int ind, indp, indm;
  E_Int closed = 0;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z, sx, sy, sz, s, c;

  //test si le i-array est une boucle fermee
  if (K_FUNC::E_abs(xt[0]-xt[ni-1]) < eps && 
      K_FUNC::E_abs(yt[0]-yt[ni-1]) < eps && 
      K_FUNC::E_abs(zt[0]-zt[ni-1]) < eps ) closed = 1;
  for (ind = 1; ind < ni-1; ind++)
  {
    hmaxt[ind] = 0.;
    indm = ind-1; indp = ind+1; 
    l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
    l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
    sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
    s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
    c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
             (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
             (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));         
    if (c > eps) hmaxt[ind] = 2*s/c; //(2*aire)/base 
  }
  hmaxt[0] = 0.; hmaxt[ni-1] = 0.;
  if (closed == 1 && ni > 2) 
  {
    hmaxt[0] = 0.;
    ind = 0; indm = ni-2; indp = 1; 
    l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
    l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
    sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
    s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
    c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
             (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
             (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));         
    if (c > eps) hmaxt[ind] = 2*s/c; //(2*aire)/base     
    hmaxt[ni-1] = hmaxt[0];
  }
  else 
  {
    hmaxt[0] = hmaxt[1];
    hmaxt[ni-1] = hmaxt[ni-2];
  }
  return;
}
//=============================================================================
void K_COMPGEOM::compStructCurvatureHeight2D(
  E_Int ni, E_Int nj, 
  E_Float* xt, E_Float* yt, E_Float* zt, 
  E_Float* hmaxt)
{
  E_Float eps = 1.e-10;
  E_Int ind0, indmax, ind, indm, indp, i, ip, im, j, jp, jm;
  E_Int closed;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z, sx, sy, sz, s, c;
  //on prend le max des deux directions 
  for (ind = 0; ind < ni*nj; ind++) hmaxt[ind] = 0.;
  // parcours en j constant 
  for (j = 0; j < nj; j++)
  {
    closed = 0;
    ind0 = 0 + j*ni; indmax = ni-1 + j*ni;     
    if (K_FUNC::E_abs(xt[ind0]-xt[indmax]) < eps && 
        K_FUNC::E_abs(yt[ind0]-yt[indmax]) < eps && 
        K_FUNC::E_abs(zt[ind0]-zt[indmax]) < eps ) closed = 1;
    for (i = 1; i < ni-1; i++)
    {
      im = i-1; ip = i+1; 
      indm = im + j*ni; indp = ip + j*ni; ind = i +j*ni;
      l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
      l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
      sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
      s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
      c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
               (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
               (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));         
      if ( c > eps) hmaxt[ind] = K_FUNC::E_max(hmaxt[ind],2*s/c);//(2*aire)/base 
    }
    if ( closed == 1 && ni > 2) 
    {
      i = 0; im = ni-2; ip = 1; indm = im + j*ni; indp = ip + j*ni; ind = i +j*ni;
      l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
      l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
      sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
      s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
      c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
               (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
               (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));  
      if (c > eps) hmaxt[ind] = K_FUNC::E_max(hmaxt[ind],2*s/c);//(2*aire)/base
      hmaxt[ind+ni-1] = K_FUNC::E_max(hmaxt[ind+ni-1],hmaxt[0]);//(ni-1,j) 
    }
    else //on extrapole de l interieur
    {
      i = 0; ind = i +j*ni; 
      hmaxt[ind] = hmaxt[ind+1];
      i = ni-1; ind = i +j*ni; 
      hmaxt[ind] = hmaxt[ind-1];
    }
  }//fin parcours j constant
  for (i = 0; i < ni; i++)
  {
    closed = 0;
    ind0 = i + 0*ni; indmax = i + (nj-1)*ni;     
    if (K_FUNC::E_abs(xt[ind0]-xt[indmax]) < eps && 
        K_FUNC::E_abs(yt[ind0]-yt[indmax]) < eps && 
        K_FUNC::E_abs(zt[ind0]-zt[indmax]) < eps ) closed = 1;
    for (j = 1; j < nj-1; j++)
    {
      jm = j-1; jp = j+1; 
      indm = i + jm*ni; indp = i + jp*ni; ind = i +j*ni;
      l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
      l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
      sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
      s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
      c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
               (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
               (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));         
      if (c > eps) hmaxt[ind] = K_FUNC::E_max(hmaxt[ind],2*s/c);//(2*aire)/base
    }
    if (closed == 1 && nj > 2)
    {
      j = 0; jm = nj-2; jp = 1; indm = i + jm*ni; indp = i + jp*ni; ind = i +j*ni;
      l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
      l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
      sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
      s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
      c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
               (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
               (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));  
      if (c > eps) hmaxt[ind] = K_FUNC::E_max(hmaxt[ind],2*s/c);//(2*aire)/base
      hmaxt[ind+(nj-1)*ni] = K_FUNC::E_max(hmaxt[ind+(nj-1)*ni],hmaxt[0]);//(i,nj-1) 
    }
    else //on extrapole de l interieur
    {
      j = 0; ind = i +j*ni; 
      hmaxt[ind] = hmaxt[ind+ni];
      j = nj-1; ind = i +j*ni; 
      hmaxt[ind] = hmaxt[ind-ni];
    }
  }//fin parcours j constant
  return;
}
//=============================================================================
// la BAR ne peut pas avoir des T-Branches 
//=============================================================================
void K_COMPGEOM::compCurvatureHeightForBAR(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  FldArrayI& cn, E_Float* hmaxt)
{
  E_Float eps = 1.e-10;
  vector< vector<E_Int> > cVVN(npts); //connectivite noeuds/noeuds voisins
  K_CONNECT::connectEV2VNbrs(cn, cVVN);
  E_Int nvoisins, ind, indp, indm;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z, sx, sy, sz, s, c;

  for (ind = 0; ind < npts; ind++)
  {
    vector<E_Int>& voisins = cVVN[ind];
    nvoisins = voisins.size();
    hmaxt[ind] = 0.;
    if ( nvoisins == 2 )
    {
      indm = voisins[0]-1; indp = voisins[1]-1;
      l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
      l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
      sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
      s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
      c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
               (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
               (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));         
      if ( c > eps) hmaxt[ind] = 2*s/c;//(2*aire)/base 
    }
  }
  // deuxieme passe pour les pts n 'ayant pas 2 voisins
  for (ind = 0; ind < npts; ind++)
  {
    vector<E_Int>& voisins = cVVN[ind];
    nvoisins = voisins.size();
    if ( nvoisins != 2 ) 
    {
      indm = voisins[0]-1;
      hmaxt[ind] = hmaxt[indm];
    }
  }
  return;
}
//=============================================================================
void K_COMPGEOM::compCurvatureHeightForTRIQUAD(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  FldArrayI& cn, E_Float* hmaxt)
{
  E_Int ind, indm, indp, voisinsE1Size, voisinsE2Size;
  vector< vector<E_Int> > cVVN(npts); K_CONNECT::connectEV2VNbrs(cn, cVVN);
  vector< vector<E_Int> > cVE(npts); K_CONNECT::connectEV2VE(cn, cVE);
  E_Float l1x, l1y, l1z, l2x, l2y, l2z, sx, sy, sz, s, c, hmaxl;
  E_Float eps = 1.e-10;
  for (ind = 0; ind < npts; ind++)
  {
    hmaxl = 0.;
    vector<E_Int>& voisinsN = cVVN[ind];//indices demarrent a 1
    E_Int nvertVoisins = voisinsN.size();
    for (E_Int nov1 = 0; nov1 < nvertVoisins; nov1++)
    {
      indm = voisinsN[nov1]-1;
      vector<E_Int>& voisinsE1 = cVE[indm];
      voisinsE1Size = voisinsE1.size();
      for (E_Int nov2 = nov1+1; nov2 < nvertVoisins; nov2++)
      {
        indp = voisinsN[nov2]-1;
        vector<E_Int>& voisinsE2 = cVE[indp];
        voisinsE2Size = voisinsE2.size();
        // recherche si appartiennent a un meme element
        for (E_Int noe1 = 0; noe1 < voisinsE1Size; noe1++)
          for (E_Int noe2 = 0; noe2 < voisinsE2Size; noe2++)
          {
            if ( voisinsE1[noe1] == voisinsE2[noe2] ) goto nextv2;
          }
        //n'appartiennent pas a un meme element : on triangule        
        l1x = xt[ind]-xt[indm]; l1y = yt[ind]-yt[indm]; l1z = zt[ind]-zt[indm];
        l2x = xt[ind]-xt[indp]; l2y = yt[ind]-yt[indp]; l2z = zt[ind]-zt[indp];
        sx = (l1y*l2z-l1z*l2y); sy = (l1z*l2x-l1x*l2z); sz = (l1x*l2y-l1y*l2x);
        s = 0.5*sqrt(sx*sx+sy*sy+sz*sz);       
        c = sqrt((xt[indp]-xt[indm])*(xt[indp]-xt[indm])+
                 (yt[indp]-yt[indm])*(yt[indp]-yt[indm])+
                 (zt[indp]-zt[indm])*(zt[indp]-zt[indm]));         
        if (c > eps) hmaxl = K_FUNC::E_max(hmaxl,2*s/c);//(2*aire)/base 
        
        nextv2:;
      }
    }
    hmaxt[ind] = hmaxl;
  }
  return;
}

