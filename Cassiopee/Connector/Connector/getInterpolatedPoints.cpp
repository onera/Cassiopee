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

#include "connector.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Search for the fringe of interpolated nodes near blanked points; depth is
   the number of layers of interpolated nodes.
   IN: blankedCells: -1, point masque, 0: point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedNodesUnstr(
  E_Int depth, FldArrayI& connect,
  FldArrayI& blankedCells,
  FldArrayI& cellN)
{
  E_Int nvert = blankedCells.getSize();
  std::vector< std::vector<E_Int> > cVN(nvert);
  E_Int isNGon = connect.getNGonType();
  if (isNGon == 0)
    K_CONNECT::connectEV2VNbrs(connect, cVN);
  else
    K_CONNECT::connectNG2VNbrs(connect, cVN);

  E_Int nvoisins;

  for (E_Int ind = 0; ind < nvert; ind++)
  {
    if (blankedCells[ind] == -1)
    {
      std::vector<E_Int>& voisins = cVN[ind];
      nvoisins = voisins.size();
      for (E_Int nov = 0; nov < nvoisins; nov++)
      {
        E_Int indv = voisins[nov]-1;
        cellN[indv] = K_FUNC::E_min(0,cellN[indv]);
      }
    }
  }
  FldArrayIS tag(nvert); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int ind = 0; ind<nvert; ind++)
    {
      if (cellN[ind] == 0)// pt interpole
      {
        std::vector<E_Int>& voisins = cVN[ind];
        nvoisins = voisins.size();
        for (E_Int nov = 0; nov < nvoisins; nov++)
        {
          E_Int indv = voisins[nov]-1;
          if (cellN[indv] == 1) tag[indv] = 1;
        }
      }
    }
  }
  for (E_Int ind = 0; ind < nvert; ind++)
  { if (tag[ind] == 1) cellN[ind] = 0; }
  return;
}

//=============================================================================
/* Search for the fringe of interpolated cells near blanked points; depth is
   the number of layers of interpolated cells.
   IN: blankedCells: -1, point masque, 0 : point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsNGON(E_Int depth, FldArrayI& cNG,
                                                  FldArrayI& blankedCells,
                                                  FldArrayI& cellN)
{
  FldArrayI cFE;
  E_Int* cnp = cNG.begin();
  E_Int sizeFN = cnp[1];         // taille de la connectivite face/noeuds
  E_Int nelts = cnp[sizeFN+2];         // nombre d elements
  std::vector< std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectNG2FE(cNG, cFE);
  K_CONNECT::connectFE2EENbrs(cFE, cEEN);
  E_Int nvoisins;

  //1st layer, depth = 1
  for (E_Int et = 0; et < nelts; et++)
  {
    if (blankedCells[et] == -1)// pt masque
    {
      std::vector<E_Int>& voisins = cEEN[et];
      nvoisins = voisins.size();
      for (E_Int noev = 0; noev < nvoisins; noev++)
      {
        E_Int et2 = voisins[noev];
        cellN[et2] = K_FUNC::E_min(0,cellN[et2]);
      }
    }
  }

  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (cellN[et] == 0)// pt interpole
      {
        std::vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int noev = 0; noev < nvoisins; noev++)
        {
          E_Int et2 = voisins[noev];
          if (cellN[et2] == 1) tag[et2] = 1;
        }
      }
    }
  }
  for (E_Int et = 0; et < nelts; et++)
  { if ( tag[et] == 1) cellN[et] = 0; }
}

//=============================================================================
/* Search for the fringe of interpolated cells near blanked points; depth is
   the number of layers of interpolated cells.
   IN: blankedCells: -1, point masque, 0 : point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsUnstr(char* eltType,
                                                   E_Int depth, FldArrayI& cnEV,
                                                   FldArrayI& blankedCells,
                                                   FldArrayI& cellN)
{
  E_Int nelts = cnEV.getSize();
  E_Int nvert = nelts*cnEV.getNfld();
  std::vector< std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, nvert, cnEV, cEEN);

  E_Int nvoisins;

  //1st layer, depth = 1
  for (E_Int et = 0; et < nelts; et++)
  {
    if (blankedCells[et] == -1)// pt masque
    {
      std::vector<E_Int>& voisins = cEEN[et];
      nvoisins = voisins.size();
      for (E_Int noev = 0; noev < nvoisins; noev++)
      {
        E_Int et2 = voisins[noev];
        cellN[et2] = K_FUNC::E_min(0,cellN[et2]);
      }
    }
  }

  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (cellN[et] == 0)// pt interpole
      {
        std::vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int noev = 0; noev < nvoisins; noev++)
        {
          E_Int et2 = voisins[noev];
          if (cellN[et2] == 1) tag[et2] = 1;
        }
      }
    }
  }
  for (E_Int et = 0; et < nelts; et++)
  { if ( tag[et] == 1) cellN[et] = 0; }
}

//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsStruct(E_Int imc, E_Int jmc, E_Int kmc, E_Int depth, E_Int dir,
                                                    FldArrayI& blankedCells, FldArrayI& cellN)
{
  E_Int imjmc = imc*jmc;
  E_Int imjmkmc = imjmc*kmc;
  E_Int i, j, k, sensor, unmsensor, ind2;
  E_Int im1, ip1, jm1, jp1, km1, kp1;
  E_Int km1imjmc, kimjmc, kp1imjmc;
  E_Int nindices;

  // On n'etend que les points masques (blankedcells = -1)
  if (dir == 0) //directionnel
  {
    if (kmc == 1)
    {
      nindices = 4;
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          j = ind/imc;
          i = ind-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;

          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          indices[0] = im1 + j*imc;
          indices[1] = ip1 + j*imc;
          indices[2] = i + jm1*imc;
          indices[3] = i + jp1*imc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }
      }
    }// fin 2D
    else
    {
      nindices = 6;
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          k = ind/imjmc;
          j = ( ind-k*imjmc )/imc;
          i = ind-k*imjmc-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;

          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          km1 = K_FUNC::E_max(0,k-d); kp1 = K_FUNC::E_min(k+d,kmc-1);

          indices[0] = im1 + j*imc + k*imjmc;
          indices[1] = ip1 + j*imc + k*imjmc;
          indices[2] = i + jm1*imc + k*imjmc;
          indices[3] = i + jp1*imc + k*imjmc;
          indices[4] = i + j*imc + km1*imjmc;
          indices[5] = i + j*imc + kp1*imjmc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }
      }
    }//fin 3D dir = 0
  }//dir = 0
  else
  {
    if (kmc == 1)
    {
      nindices = 8;
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          j = ind/imc;
          i = ind-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;

          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          indices[0] = im1 + jm1*imc;
          indices[1] = i + jm1*imc;
          indices[2] = ip1 + jm1*imc;

          indices[3] = im1 + j*imc;
          indices[4] = ip1 + j*imc;

          indices[5] = im1 + jp1*imc;
          indices[6] = i  +  jp1*imc;
          indices[7] = ip1 + jp1*imc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }
      }
    }// 2D dir = 1
    else // 3D
    {
      nindices = 26;
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          k = ind/imjmc;
          j = ( ind-k*imjmc )/imc;
          i = ind-k*imjmc-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;

          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          km1 = K_FUNC::E_max(0,k-d); kp1 = K_FUNC::E_min(k+d,kmc-1);

          km1imjmc= km1*imjmc;
          kp1imjmc= kp1*imjmc;
          kimjmc= k*imjmc;

          indices[0] = im1 + jm1*imc + km1imjmc;
          indices[1] = i   + jm1*imc + km1imjmc;
          indices[2] = ip1 + jm1*imc + km1imjmc;

          indices[3] = im1 + j*imc + km1imjmc;
          indices[4] = i   + j*imc + km1imjmc;
          indices[5] = ip1 + j*imc + km1imjmc;

          indices[6] = im1 + jp1*imc + km1imjmc;
          indices[7] = i  +  jp1*imc + km1imjmc;
          indices[8] = ip1 + jp1*imc + km1imjmc;

          indices[9]  = im1 + jm1*imc + kimjmc;
          indices[10] = i   + jm1*imc + kimjmc;
          indices[11] = ip1 + jm1*imc + kimjmc;

          indices[12] = im1 + j*imc + kimjmc;
          indices[13] = ip1 + j*imc + kimjmc;

          indices[14] = im1 + jp1*imc + kimjmc;
          indices[15] = i  +  jp1*imc + kimjmc;
          indices[16] = ip1 + jp1*imc + kimjmc;

          indices[17] = im1 + jm1*imc + kp1imjmc;
          indices[18] = i   + jm1*imc + kp1imjmc;
          indices[19] = ip1 + jm1*imc + kp1imjmc;

          indices[20] = im1 + j*imc + kp1imjmc;
          indices[21] = i   + j*imc + kp1imjmc;
          indices[22] = ip1 + j*imc + kp1imjmc;

          indices[23] = im1 + jp1*imc + kp1imjmc;
          indices[24] = i  +  jp1*imc + kp1imjmc;
          indices[25] = ip1 + jp1*imc + kp1imjmc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }
      }
    }
  }//dir = 1
}

//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsStructOpt(E_Int imc, E_Int jmc, E_Int kmc, E_Int depth, E_Int dir,
                                                       E_Float* cellN, E_Float* cellN_tmp)
{
  E_Int imjmc = imc*jmc;
  E_Int imjmkmc = imjmc*kmc;
  E_Int nindices;

  if (dir == 0) //stencil par direction
  {
    if (kmc == 1) // 2D croix
    {
      nindices = 4*depth;

      #pragma omp parallel
      {
        // Def de variables privees sur les procs
        vector<E_Int> indices(nindices);
        E_Int i, j, ii, jj;
        E_Int ind2, compteur;

        #pragma omp for schedule(static)
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          if (K_FUNC::fEqual(cellN[ind],1.)) // Si cellN = 1. Changements a faire en fonction du stencil
          {
            //indices de la maille
            j = ind/imc;
            i = ind-j*imc;

            // Recherche des indices dans le stencil
            compteur = 0;
            for (E_Int d=depth; d>0; d--) // branche haute
            {
              jj = K_FUNC::E_max(0,j-d);
              indices[compteur] = i + jj*imc; compteur++;
            }
            for (E_Int d=depth; d>0; d--) // branche gauche
            {
              ii = K_FUNC::E_max(0,i-d);
              indices[compteur] = ii + j*imc; compteur++;
            }
            for (E_Int d=1; d<depth+1; d++)// branche droite
            {
              ii = K_FUNC::E_min(i+d, imc-1);
              indices[compteur] =  ii + j*imc; compteur++;
            }
            for (E_Int d=1; d<depth+1; d++) // branche basse
            {
              jj = K_FUNC::E_min(j+d, jmc-1);
              indices[compteur] = i + jj*imc; compteur++;
            }

            // Changement du cellN en fonction du stencil
            for (E_Int noi = 0; noi < nindices; noi++)
            {
              ind2 = indices[noi];
              if (K_FUNC::fEqualZero(cellN[ind2])){ cellN_tmp[ind] = 2.; break;}
            }
          }
        }
      }
    }// fin 2D
    else // 3D croix
    {
      nindices = 6*depth;

      #pragma omp parallel
      {
        // Def de variables privees sur les procs
        vector<E_Int> indices(nindices);
        E_Int i, j, k, ii, jj, kk;
        E_Int ind2, compteur;

        #pragma omp for schedule(guided)
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          if (K_FUNC::fEqual(cellN[ind],1.))
          {
            //indices de la maille
            k = ind/imjmc;
            j = ( ind-k*imjmc )/imc;
            i = ind-k*imjmc-j*imc;

            // Recherche des indices dans le stencil
            compteur = 0;
            for (E_Int d=-depth; d<0; d++) // branche arriere
            {
              kk = K_FUNC::E_max(0,k+d);
              indices[compteur] = i  + j*imc  + kk*imjmc ; compteur++;
            }
            for (E_Int d=-depth; d<0; d++) // branche haute
            {
              jj = K_FUNC::E_max(0,j+d);
              indices[compteur] = i  + jj*imc + k*imjmc  ; compteur++;
            }
            for (E_Int d=-depth; d<0; d++) // branche gauche
            {
              ii = K_FUNC::E_max(0,i+d);
              indices[compteur] = ii + j*imc  + k*imjmc  ; compteur++;
            }
            for (E_Int d=1; d<depth+1; d++)// branche droite
            {
              ii = K_FUNC::E_min(i+d, imc-1);
              indices[compteur] = ii + j*imc  + k*imjmc  ; compteur++;
            }
            for (E_Int d=1; d<depth+1; d++) // branche basse
            {
              jj = K_FUNC::E_min(j+d, jmc-1);
              indices[compteur] = i  + jj*imc + k*imjmc  ; compteur++;
            }
            for (E_Int d=1; d<depth+1; d++) // branche avant
            {
              kk = K_FUNC::E_min(k+d, kmc-1);
              indices[compteur] = i  + j*imc  + kk*imjmc ; compteur++;
            }

            // Changement du cellN en fonction du stencil
            for (E_Int noi = 0; noi < nindices; noi++)
            {
              ind2 = indices[noi];
              if (K_FUNC::fEqualZero(cellN[ind2])){ cellN_tmp[ind] = 2.; break;}
            }
          }
        }
      }
    }//fin 3D dir = 0
  }//dir = 0
  else if (dir == 1)//stencil etoile
  {
    if (kmc == 1) // 2D etoile
    {
      nindices = 8*depth;

      #pragma omp parallel
      {
        // Def de variables privees sur les procs
        vector<E_Int> indices(nindices);
        E_Int i, j, ii, jj;
        E_Int ind2, compteur;

        #pragma omp for schedule(static)
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          if (K_FUNC::fEqual(cellN[ind],1.))
          {
            //indices de la maille
            j = ind/imc;
            i = ind-j*imc;

            // Recherche des points du stencil
            compteur = 0;
            for (E_Int d=depth; d>0; d--) //stencil au dessus de la maille
            {
              // E_Int jj = K_FUNC::E_max(j-d, 0);
              jj = j-d; if (jj<0) {jj=j;}
              ii = i-d; if (ii<0) {ii=i;}
              indices[compteur] = ii + jj*imc; compteur++;
              indices[compteur] = i + jj*imc; compteur++;
              ii = i+d; if (ii>imc-1) {ii=i;}
              indices[compteur] = ii + jj*imc; compteur++;
            }
            for (E_Int d = -depth; d<depth+1; d++)
            {
              if (d!=0)
              {
                ii = i+d; if((ii<0)||(ii>imc-1)) {ii=i;}
                indices[compteur] = ii + j*imc; compteur++;
              }
            }
            for (E_Int d=1; d<depth+1; d++) //stencil au dessus de la maille
            {
              // E_Int jj = K_FUNC::E_min(j+d, jmc-1);
              jj = j+d; if (jj>jmc-1) {jj=j;}
              ii = i-d; if (ii<0) {ii=i;}
              indices[compteur] = ii + jj*imc; compteur++;
              indices[compteur] = i + jj*imc; compteur++;
              ii = i+d; if (ii>imc-1) {ii=i;}
              indices[compteur] = ii + jj*imc; compteur++;
            }

            // Changement du cellN en fonction du stencil
            for (E_Int noi = 0; noi < nindices; noi++)
            {
              ind2 = indices[noi];
              if (K_FUNC::fEqualZero(cellN[ind2])){ cellN_tmp[ind] = 2.; break;}
            }
          }
        }
      }
    }// 2D dir = 1
    else // 3D etoile
    {
      nindices = 26*depth;

      #pragma omp parallel
      {
        // Def de variables privees sur les procs
        vector<E_Int> indices(nindices);
        E_Int i, j, k, ii, jj, kk;
        E_Int ind2, compteur;

        #pragma omp for schedule(guided)
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          if (K_FUNC::fEqual(cellN[ind],1.))
          {
            //indices de la maille
            k = ind/imjmc;
            j = ( ind-k*imjmc )/imc;
            i = ind-k*imjmc-j*imc;

            // Recherche des points du stencil
            compteur = 0;
            //----------- stencil dans les Z negatifs -----------
            for (E_Int kd=-depth; kd<0; kd++) // stencil dans les z negatifs
            {
              kk = k+kd; if(kk<0) {kk=k;}
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd*K_FUNC::E_abs(kd); if((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-1; id<2; id++)
                {
                  ii = i+id*K_FUNC::E_abs(kd); if((ii<0)||(ii>imc-1)) {ii=i;}
                  indices[compteur] = ii + jj*imc + kk*imjmc; compteur++;
                }
              }
            }
            //----------- stencil dans les Z nuls -----------
            kk=k;
            for (E_Int d=depth; d>0; d--) //stencil au dessus de la maille
            {
              jj = j-d; if (jj<0) {jj=j;}
              ii = i-d; if (ii<0) {ii=i;}
              indices[compteur] = ii + jj*imc + kk*imjmc; compteur++;
              indices[compteur] = i  + jj*imc + kk*imjmc; compteur++;
              ii = i+d; if (ii>imc-1) {ii=i;}
              indices[compteur] = ii + jj*imc + kk*imjmc; compteur++;
            }
            for (E_Int d = -depth; d<depth+1; d++)
            {
              if (d!=0)
              {
                ii = i+d; if((ii<0)||(ii>imc-1)) {ii=i;}
                indices[compteur] = ii + j*imc + kk*imjmc; compteur++;
              }
            }
            for (E_Int d=1; d<depth+1; d++) //stencil au dessus de la maille
            {
              jj = j+d; if (jj>jmc-1) {jj=j;}
              ii = i-d; if (ii<0) {ii=i;}
              indices[compteur] = ii + jj*imc + kk*imjmc; compteur++;
              indices[compteur] = i  + jj*imc + kk*imjmc; compteur++;
              ii = i+d; if (ii>imc-1) {ii=i;}
              indices[compteur] = ii + jj*imc + kk*imjmc; compteur++;
            }
            //----------- stencil dans les Z positifs -----------
            for (E_Int kd=1; kd<depth+1; kd++)
            {
              kk = k+kd; if (kk>kmc-1) {kk=k;}
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd*K_FUNC::E_abs(kd); if((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-1; id<2; id++)
                {
                  ii = i+id*K_FUNC::E_abs(kd); if((ii<0)||(ii>imc-1)) {ii=i;}
                  indices[compteur] = ii + jj*imc + kk*imjmc; compteur++;
                }
              }
            }
            // Changement du cellN en fonction du stencil
            for (E_Int noi = 0; noi < nindices; noi++)
            {
              ind2 = indices[noi];
              if (K_FUNC::fEqualZero(cellN[ind2])){ cellN_tmp[ind] = 2.; break;}
            }
          }
        }
      }
    }
  }//dir = 1
  else if (dir == 2) //stencil losange (dir==2)
  {
    if (kmc == 1) // 2D losange
    {
      nindices = (depth+1)*(depth+1) + depth*depth -1;
      #pragma omp parallel
      {
        // Def de variables privees sur les procs
        vector<E_Int> indices(nindices);
        E_Int i, j, ii, jj;
        E_Int ind2, compteur, ncouche;

        #pragma omp for schedule(static)
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          if (K_FUNC::fEqual(cellN[ind],1.))
          {
            //indices de la maille
            j = ind/imc;
            i = ind-j*imc;

            // Recherche des points du stencil
            compteur = 0;
            ncouche = 0;
            for (E_Int d=depth; d>0; d--) //stencil au dessus de la maille
            {
              jj = j-d; if (jj<0) {jj=j;}
              for (E_Int id=-ncouche; id<ncouche+1; id++)
              {
                ii = i+id ; if ((ii<0)||(ii>imc-1)) {ii=i;}
                indices[compteur] = ii + jj*imc; compteur++;
              }
              ncouche += 1;
            }
            for (E_Int d = -depth; d<depth+1; d++)
            {
              if (d!=0)
              {
                ii = i+d; if((ii<0)||(ii>imc-1)) {ii=i;}
                indices[compteur] = ii + j*imc; compteur++;
              }
            }
            ncouche=depth-1;
            for (E_Int d=1; d<depth+1; d++) //stencil au dessous de la maille
            {
              jj = j+d; if (jj>jmc-1) {jj=j;}
              for (E_Int id=-ncouche; id<ncouche+1; id++)
              {
                ii = i+id ; if ((ii<0)||(ii>imc-1)) {ii=i;}
                indices[compteur] = ii + jj*imc; compteur++;
              }
              ncouche -=1;
            }

            // Changement du cellN en fonction du stencil
            for (E_Int noi = 0; noi < nindices; noi++)
            {
              ind2 = indices[noi];
              if (K_FUNC::fEqualZero(cellN[ind2])){ cellN_tmp[ind] = 2.; break;}
            }
          }
        }
      }
    }// 2D dir = 1
    else // 3D losange
    {
      // Nombre d'indices en fonction de la profondeur
      if (depth == 2) 
      {
        nindices = 33;
      }
      else if (depth == 3) 
      {
        nindices = 87;
      }
      else 
      {
        nindices = 0;
        printf("WARNING: searchMaskInterpolatedCellsStructOpt: bad choice of depth for dir=2.\n");
      }

      #pragma omp parallel
      {
        // Def de variables privees sur les procs
        vector<E_Int> indices(nindices);
        E_Int i, j, k, ii, jj, kk, l;
        E_Int ind2, im, compteur;

        #pragma omp for schedule(guided)
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          if (K_FUNC::fEqual(cellN[ind],1.))
          {
            //indices de la maille
            k = ind/imjmc;
            j = ( ind-k*imjmc )/imc;
            i = ind-k*imjmc-j*imc;

            compteur  = 0;
            if (depth==2) 
            { // Stencil de taille 2
              // premiere couche du stencil
              kk = k-2; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              ii = i ; jj = j ;
              l = ii + jj*imc + kk*imjmc;
              indices[compteur] = l; compteur++;

              // deuxieme couche du stencil
              kk = k-1; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-1; id<2; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
              }

              // troisieme couche du stencil
              kk = k; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              im = 0;
              for (E_Int jd=-depth; jd<depth+1; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-im; id<im+1; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
                if (jd<0) {im +=1;} // haut du stencil : augmente la taille de la ligne
                else      {im -=1;} // cas du stencil : diminue la taille de la ligne
              }

              // quatrieme couche du stencil
              kk = k+1; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-1; id<2; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
              }

              // derniere couche du stencil
              kk = k+2; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              ii = i ; jj = j ;
              l = ii + jj*imc + kk*imjmc;
              indices[compteur] = l; compteur++;
            } // fin depth = 2
            else if (depth==3) 
            { // stencil de taille 3
              // premiere couche du stencil
              kk = k-3; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              ii = i ; jj = j ;
              l = ii + jj*imc + kk*imjmc;
              indices[compteur] = l ; compteur++;

              // deuxieme couche du stencil
              kk = k-2; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-1; id<2; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
              }

              // troisieme couche du stencil
              kk = k-1; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              // branche haute
              jj = j-2; if ((jj<0)||(jj>jmc-1)) {jj=j;}
              for (E_Int id=-1; id<2; id++)
              {
                ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                l = ii + jj*imc + kk*imjmc;
                indices[compteur] = l; compteur++;
              }
              // milieu
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-2; id<3; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
              }
              // branche basse
              jj = j+2; if ((jj<0)||(jj>jmc-1)) {jj=j;}
              for (E_Int id=-1; id<2; id++)
              {
                ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                l = ii + jj*imc + kk*imjmc;
                indices[compteur] = l; compteur++;
              }

              // quatrieme couche du stencil
              kk = k; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              im = 0;
              for (E_Int jd=-depth; jd<depth+1; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-im; id<im+1; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
                if (jd<0) {im +=1;} // haut du stencil : augmente la taille de la ligne
                else      {im -=1;} // cas du stencil : diminue la taille de la ligne
              }

              // cinquieme couche du stencil
              kk = k+1; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              // branche haute
              jj = j-2; if ((jj<0)||(jj>jmc-1)) {jj=j;}
              for (E_Int id=-1; id<2; id++)
              {
                ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                l = ii + jj*imc + kk*imjmc;
                indices[compteur] = l; compteur++;
              }
              // milieu
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-2; id<3; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
              }
              // branche basse
              jj = j+2; if ((jj<0)||(jj>jmc-1)) {jj=j;}
              for (E_Int id=-1; id<2; id++)
              {
                ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                l = ii + jj*imc + kk*imjmc;
                indices[compteur] = l; compteur++;
              }

              // sixieme couche du stencil
              kk = k+2; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              for (E_Int jd=-1; jd<2; jd++)
              {
                jj = j+jd; if ((jj<0)||(jj>jmc-1)) {jj=j;}
                for (E_Int id=-1; id<2; id++)
                {
                  ii = i+id; if ((ii<0)||(ii>imc-1)) {ii=i;}
                  l = ii + jj*imc + kk*imjmc;
                  indices[compteur] = l; compteur++;
                }
              }

              // septieme couche du stencil
              kk = k+3; if ((kk<0)||(kk>kmc-1)) {kk=k;}
              ii = i ; jj = j ;
              l = ii + jj*imc + kk*imjmc;
              indices[compteur] = l; compteur++;
            } // fin du stencil de taille 3
            else 
            {
              printf("WARNING: searchMaskInterpolatedCellsStructOpt: bad choice of depth for dir=2.\n");
            }

            // Changement du cellN en fonction du stencil
            for (E_Int noi = 0; noi < nindices; noi++)
            {
              ind2 = indices[noi];
              if (K_FUNC::fEqualZero(cellN[ind2])){ cellN_tmp[ind] = 2.; break;}
            }
          }
        }
      }
    }
  }//dir = 2
  else if (dir == 3) //stencil octahedron (dir==3)
  {
    // Nombre d'indices en fonction de la profondeur
    if (kmc == 1)
    {
      nindices = pow(2*depth+1, 2) - 2*depth*(depth+1);
    }
    else
    {
      nindices = 0;
      for (E_Int id=0; id < depth; id++)
      { 
        nindices += pow(2*id+1, 2) - 2*id*(id+1);   
      }
      nindices = 2*nindices + pow(2*depth+1, 2) - 2*depth*(depth+1);
    }
    
    #pragma omp parallel
    {
      // Def de variables privees sur les procs
      vector<E_Int> indices(nindices);
      E_Int i, j, k;
      E_Int ind2, compteur;
      
      // Shorthands
      E_Int numKLayers = 1;
      
      // In 2D, set k to 0, otherwise 3D octahedron
      if (kmc == 1) k = 0;
      else numKLayers = depth + 1;

      #pragma omp for schedule(guided)
      for (E_Int ind = 0; ind < imjmkmc; ind++)
      {
        // Skip rest of the loop if cellN different from 1
        if (not K_FUNC::fEqual(cellN[ind], 1.)) continue;
          
        // Indices de la maille
        if (kmc != 1) k = ind/imjmc;
        j = (ind - k*imjmc)/imc;
        i = ind -k*imjmc - j*imc;
          
        compteur = 0;
        // Loop over half of all vertical layers and use the symmetry for
        // the negative half
        for (E_Int kd=0; kd<numKLayers; kd++)
        {
          E_Int dLyr = depth - abs(kd);
          E_Int kkm = K_FUNC::E_max(0, k-kd);
          E_Int kkp = K_FUNC::E_min(kmc-1, k+kd);
              
          for (E_Int id=0; id<=dLyr; id++)
          {
            E_Int iim = K_FUNC::E_max(0, i - id);
            E_Int iip = K_FUNC::E_min(imc-1, i + id);
                  
            for (E_Int jd=0; jd<=dLyr-id; jd++)
            {
              E_Int jjm = K_FUNC::E_max(0, j - jd);
              E_Int jjp = K_FUNC::E_min(jmc-1, j + jd);
                      
              indices[compteur] = iim + jjm*imc + kkm*imjmc; compteur++;
                      
              if (kd != 0)
              {
                indices[compteur] = iim + jjm*imc + kkp*imjmc; compteur++;
              }
                      
              if (jd != 0)
              {
                indices[compteur] = iim + jjp*imc + kkm*imjmc; compteur++;
                          
                if (kd != 0)
                {
                  indices[compteur] = iim + jjp*imc + kkp*imjmc; compteur++;
                }
              }

              if (id == 0) continue;
                      
              indices[compteur] = iip + jjm*imc + kkm*imjmc; compteur++;
                      
              if (kd != 0)
              {
                indices[compteur] = iip + jjm*imc + kkp*imjmc; compteur++;
              }
                      
              if (jd != 0)
              {
                indices[compteur] = iip + jjp*imc + kkm*imjmc; compteur++;
                          
                if (kd != 0)
                {
                  indices[compteur] = iip + jjp*imc + kkp*imjmc; compteur++;
                }
              }
            }
          }
        }
          
        // Changement du cellN en fonction du stencil
        for (E_Int noi = 0; noi < nindices; noi++)
        {
          ind2 = indices[noi];
          if (K_FUNC::fEqualZero(cellN[ind2]))
          {
            cellN_tmp[ind] = 2.; 
            break;
          }
        }
      }
    }
  }
}

//=============================================================================
/* Determine les noeuds interpoles a partir du cellN en noeuds
   Si le celln contient des pts masques, alors les points interpoles autour
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::getOversetHolesInterpNodes(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLE_(args, O_ II_ S_, &array, &depth, &dir, &cellNName))
  {
    return NULL;
  }
  if (dir < 0 || dir > 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: dir must be between 0 and 3.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  //E_Int res = K_ARRAY::getFromArray3(array, varString,
  //                                   field, im, jm, km, cn, eltType);
  E_Int res = K_ARRAY::getFromArray(array, varString,
                                    field, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: first argument is not recognized");
    return NULL;
  }

  if (dir == 2 && depth>3 && km>1)
  {
    printf("WARNING: getOversetHolesInterpNodes: dir=2, depth>3 and 3D are incompatibles. Force dir=1.\n");
    dir = 1;
  }
  if (dir == 2 && depth == 1)
  {
    printf("WARNING: getOversetHolesInterpNodes: dir=2, depth=1 and 3D are incompatibles. Force dir=1.\n");
    dir = 1;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);

  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: array must contain cellN variable.");
    delete field; if (res == 2) delete cn;
    //RELEASESHAREDB(res, array, field, cn);
    return NULL;
  }
  posc++;

  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int npts = field->getSize();
  E_Int api = field->getApi();

  if (res == 1)
  {
    E_Float* cellNp_tmp;
    FldArrayF cellN_tmp(npts);
    cellNp_tmp = cellN_tmp.begin();

    #pragma omp parallel
    {
      #pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        cellN_tmp[ind] = cellNp[ind];
      }
    }

    searchMaskInterpolatedCellsStructOpt(im, jm, km, depth, dir, cellNp, cellNp_tmp);

    #pragma omp parallel
    {
      #pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        cellNp[ind] = cellNp_tmp[ind];
      }
    }

    PyObject* tpl = K_ARRAY::buildArray3(*field, varString, im, jm, km);
    delete field;
    //RELEASESHAREDS(array, field);
    return tpl;
  }
  else
  {
    FldArrayI blankedCells(npts); blankedCells.setAllValuesAt(1);
    FldArrayI cellNatFld(npts); cellNatFld.setAllValuesAt(1);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (K_FUNC::fEqualZero(cellNp[ind] - 2.)) { blankedCells[ind] = 0; cellNatFld[ind] = 0;}
      else if (K_FUNC::fEqualZero(cellNp[ind])) { blankedCells[ind] = -1; cellNatFld[ind] = -1;}
    }
    // WARNING: NGON is array1 type here !!! 
    if (K_STRING::cmp(eltType, "NGON") == 0) cn->setNGonType(1);
 
    searchMaskInterpolatedNodesUnstr(depth, *cn, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.;
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.;
    }

    PyObject* tpl = K_ARRAY::buildArray3(*field, varString, *cn, eltType, api);
    //RELEASESHAREDU(array, field, cn);
    delete field; delete cn;
    return tpl;
  }
}

//=============================================================================
/* Determine les noeuds interpoles a partir du cellN en noeuds
   Si le celln contient des pts masques, alors les points interpoles autour
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::_getOversetHolesInterpNodes(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLE_(args, O_ II_ S_, &array, &depth, &dir, &cellNName))
  {
    return NULL;
  }
  if (dir != 0 && dir != 1 && dir != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_getOversetHolesInterpNodes: dir must be 0, 1 or 2.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString,
                                     field, im, jm, km, cn, eltType);
  if (res != 1)
  {
    if (res == 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpNodes: not yet implemented for unstructured zones.");
      RELEASESHAREDU(array, field, cn);
    }
    else
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpNodes: first argument is not recognized");
    return NULL;
  }

  if (dir==2 && depth>3 && km>1)
  {
    printf("WARNING: _getOversetHolesInterpNodes: dir=2, depth>3 and 3D are incompatibles. Force dir=1.\n");
    dir=1;
  }
  if (dir == 2 && depth == 1)
  {
    printf("WARNING: _getOversetHolesInterpNodes: dir=2, depth=1 and 3D are incompatibles. Force dir=1.\n");
    dir=1;
  }
  
  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);

  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_getOversetHolesInterpNodes: array must contain cellN variable.");
    RELEASESHAREDS(array, field);return NULL;
  }
  posc++;

  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int npts = field->getSize();

  E_Float* cellNp_tmp;
  FldArrayF cellN_tmp(npts);
  cellNp_tmp = cellN_tmp.begin();

  #pragma omp parallel
  {
    #pragma omp for
    for (E_Int ind=0; ind<npts; ind++)
    {
      cellN_tmp[ind] = cellNp[ind];
    }
  }

  searchMaskInterpolatedCellsStructOpt(im, jm, km, depth, dir, cellNp, cellNp_tmp);

  #pragma omp parallel
  {
    #pragma omp for
    for (E_Int ind=0; ind<npts; ind++)
    {
      cellNp[ind] = cellNp_tmp[ind];
    }
  }

  RELEASESHAREDS(array, field);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* Determine les centres interpoles a partir du cellN
   Si le celln contient des pts masques, alors les points interpoles autour
   sont construits [structure seulement] */
//=============================================================================
PyObject* K_CONNECTOR::_getOversetHolesInterpCellCenters(PyObject* self, PyObject* args)
{
  PyObject *centersArray;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLE_(args, O_ II_ S_,
                     &centersArray, &depth, &dir, &cellNName))
  {
    return NULL;
  }

  if (dir != 0 && dir != 1 && dir != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters: dir must be 0, 1 or 2.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(centersArray, varString,
                                     field, im, jm, km, cn, eltType);
  if (res != 1)
  {
    if (res == 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpCellCenters: not yet implemented for unstructured zones.");
      RELEASESHAREDU(centersArray, field, cn);
    }
    else
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpCellCenters: first argument is not recognized");
    return NULL;
  }

  if (dir == 2 && depth > 3 && km > 1)
  {
    printf("WARNING: _getOversetHolesInterpCellCenters: dir=2, depth>3 and 3D are incompatibles. Force dir=1.\n");
    dir = 1;
  }
  if (dir == 2 && depth == 1)
  {
    printf("WARNING: _getOversetHolesInterpCellCenters: dir=2, depth=1 and 3D are incompatibles. Force dir=1.\n");
    dir = 1;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_getOversetHolesInterpCellCenters: array must contain cellN variable.");
    RELEASESHAREDS(centersArray, field);return NULL;
  }
  posc++;
  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int ncells = field->getSize();

  E_Float* cellNp_tmp;
  FldArrayF cellN_tmp(ncells);
  cellNp_tmp = cellN_tmp.begin();

  #pragma omp parallel
  {
    #pragma omp for
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      cellN_tmp[ind] = cellNp[ind];
    }
  }

  searchMaskInterpolatedCellsStructOpt(im, jm, km, depth, dir, cellNp, cellNp_tmp);

  #pragma omp parallel
  {
    #pragma omp for
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      cellNp[ind] = cellNp_tmp[ind];
    }
  }

  RELEASESHAREDS(centersArray, field);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* Determine les centres interpoles a partir du cellN
   Si le celln contient des pts masques, alors les points interpoles autour
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::getOversetHolesInterpCellCenters(PyObject* self, PyObject* args)
{
  PyObject *centersArray;
  E_Int depth; E_Int dir;
  char* cellNName;
  PyObject* bindices; PyObject* bfield;
  if (!PYPARSETUPLE_(args, O_ II_ S_ OO_,
                     &centersArray, &depth, &dir, &cellNName,
                     &bindices, &bfield))
  {
    return NULL;
  }

  if (dir != 0 && dir != 1 && dir != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters: dir must be 0, 1 or 2.");
    return NULL;
  }

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  //E_Int res = K_ARRAY::getFromArray3(centersArray, varString,
  //                                   field, im, jm, km, cn, eltType);
  E_Int res = K_ARRAY::getFromArray(centersArray, varString,
                                    field, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters: first argument is not recognized");
    return NULL;
  }

  if (dir == 2 && depth>3 && km>1)
  {
    printf("WARNING: getOversetHolesInterpCellCenters: dir=2, depth>3 and 3D are incompatibles. Force dir=1.\n");
    dir = 1;
  }
  if (dir == 2 && depth == 1)
  {
    printf("WARNING: getOversetHolesInterpCellCenters: dir=2, depth=1 and 3D are incompatibles. Force dir=1.\n");
    dir = 1;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters: array must contain cellN variable.");
    //RELEASESHAREDB(res, centersArray, field, cn);
    delete field; if (res == 2) delete cn;
    return NULL;
  }
  posc++;
  E_Float* cellNp = field->begin(posc);
  
  // face fields : face index, face field=cellN
  FldArrayI* binds=NULL; FldArrayF* bf=NULL;
  if (bindices != Py_None && bfield != Py_None)
  {
    K_NUMPY::getFromNumpyArray(bindices, binds);
    K_NUMPY::getFromNumpyArray(bfield, bf);
  }

  E_Int ncells = field->getSize();
  E_Int api = field->getApi();

  if (res == 1) // structure
  {
    E_Float* cellNp_tmp;
    FldArrayF cellN_tmp(ncells);
    cellNp_tmp = cellN_tmp.begin();

    #pragma omp parallel
    {
      #pragma omp for
      for (E_Int ind=0; ind<ncells; ind++)
      {
        cellN_tmp[ind] = cellNp[ind];
      }
    }

    searchMaskInterpolatedCellsStructOpt(im, jm, km, depth, dir, cellNp, cellNp_tmp);

    #pragma omp parallel
    {
      #pragma omp for
      for (E_Int ind=0; ind<ncells; ind++)
      {
        cellNp[ind] = cellNp_tmp[ind];
      }
    }

    PyObject* tpl = K_ARRAY::buildArray3(*field, varString, im, jm, km);
    //RELEASESHAREDS(centersArray, field);
    delete field;
    return tpl;
  }
  else // non structure
  {
    FldArrayI blankedCells(ncells); blankedCells.setAllValuesAt(1);
    FldArrayI cellNatFld(ncells); cellNatFld.setAllValuesAt(1);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (K_FUNC::fEqualZero(cellNp[ind] - 2.)){ blankedCells[ind] = 0; cellNatFld[ind] = 0; }
      else if (K_FUNC::fEqualZero(cellNp[ind])){ blankedCells[ind] = -1; cellNatFld[ind] = -1; }
    }

    if (K_STRING::cmp(eltType, "NGON*") == 0)
      searchMaskInterpolatedCellsNGON(depth, *cn, blankedCells, cellNatFld);
    else
      searchMaskInterpolatedCellsUnstr(eltType, depth, *cn, blankedCells, cellNatFld);

    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.;
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.;
    }

    // modify from bc field only for ngon and depth=1
    if (binds != NULL && bf != NULL && K_STRING::cmp(eltType, "NGON*") == 0 && depth == 1)
    {
      // build cFE connectivity
      E_Int* bindsp = binds->begin();
      E_Float* bfp = bf->begin(1); 
      K_FLD::FldArrayI cFE;
      K_CONNECT::connectNG2FE(*cn, cFE);
      E_Int indF, indE;
      for (E_Int i = 0; i < binds->getSize(); i++)
      {
        indF = bindsp[i];
        indE = cFE(indF, 1)-1;
        if (K_FUNC::fEqualZero(bfp[i]) && K_FUNC::fEqualZero(cellNp[indE] - 1.)) cellNp[indE] = 2;
      }
    }

    PyObject* tpl =  K_ARRAY::buildArray3(*field, varString, *cn, eltType, api);
    //RELEASESHAREDU(centersArray, field, cn);
    delete field; delete cn;
    return tpl;
  }
}

//===============================================================================
/* Retourne le numpy des indices des pts cellN=2 et les numpys des coordonnees
  la zone en entree est le maillage des centres
  car le cellN est localise aux noeuds pour plus d efficacite */
//===============================================================================
PyObject* K_CONNECTOR::getInterpolatedPointsZ(PyObject* self, PyObject* args)
{
  PyObject* zone;
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  char* cellNName;
  if (!PYPARSETUPLE_(args, O_ SSSS_, &zone, &cellNName, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: wrong arguments.");
    return NULL;
  }
  E_Int xyz = 1; E_Int locI = 0; // tjs aux noeuds
  E_Int ni, nj, nk, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;

  E_Int zoneType = K_PYTREE::getFromZone(zone, xyz, locI, varString, fields, locs, ni, nj, nk,
                                         cn, cnSize, cnNfld, eltType, hook, GridCoordinates,
                                         FlowSolutionNodes, FlowSolutionCenters);
  if (zoneType == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: invalid zone.");
    return NULL;
  }
  if (locs.size() < 4)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: one variable missing in zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: coordinates cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: cellN cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  if ( locs[posc] != 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: cellN must be located at nodes in input zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int nptsTot;
  if (zoneType == 1) nptsTot = ni*nj*nk;
  else nptsTot = ni;
  E_Float* cellNp = fields[posc];
  
  E_Float* xp = fields[posx];
  E_Float* yp = fields[posy];
  E_Float* zp = fields[posz];

  // dimensionnement
  E_Int nthreads = __NUMTHREADS__;
  E_Int* sizeLoc = new E_Int [nthreads];

#pragma omp parallel 
  {
    E_Int ithread = __CURRENT_THREAD__;
    sizeLoc[ithread] = 0;
#pragma omp for
  for (E_Int i = 0; i < nptsTot; i++)
  {
    if ( K_FUNC::fEqualZero(cellNp[i] - 2.) )
    {
      sizeLoc[ithread] += 1;
    }
  }
 }

  E_Int size = 0; E_Int tmp;
  for (E_Int i = 0; i < nthreads; i++)
  {  tmp = sizeLoc[i]; sizeLoc[i] = size; size += tmp; }

  if (size == 0)
  {
    RELEASESHAREDZ(hook, varString, eltType);
    delete [] sizeLoc;
    Py_INCREF(Py_None);
    return Py_None;
  }

  // allocation des numpy
  PyObject* PyIndices = K_NUMPY::buildNumpyArray(size,1,1);
  PyObject* PyCoordX = K_NUMPY::buildNumpyArray(size,1,0);
  PyObject* PyCoordY = K_NUMPY::buildNumpyArray(size,1,0);
  PyObject* PyCoordZ = K_NUMPY::buildNumpyArray(size,1,0);
  E_Int* indicesInterp = K_NUMPY::getNumpyPtrI(PyIndices);
  E_Float* coordX = K_NUMPY::getNumpyPtrF(PyCoordX);
  E_Float* coordY = K_NUMPY::getNumpyPtrF(PyCoordY);
  E_Float* coordZ = K_NUMPY::getNumpyPtrF(PyCoordZ);
  
  PyObject* tpl = Py_BuildValue("[OOOO]", PyIndices, PyCoordX, PyCoordY, PyCoordZ);
  Py_DECREF(PyIndices); Py_DECREF(PyCoordX); Py_DECREF(PyCoordY); Py_DECREF(PyCoordZ);
  
  // remplissage des numpy
#pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int noi = sizeLoc[ithread];
    #pragma omp for
    for (E_Int i = 0; i < nptsTot; i++)
    {
      if ( K_FUNC::fEqualZero(cellNp[i] - 2.) )
      {
        indicesInterp[noi] = i;
        coordX[noi] = xp[i];
        coordY[noi] = yp[i];
        coordZ[noi] = zp[i];
        noi += 1;
      }
    }
  }

  delete [] sizeLoc;
  RELEASESHAREDZ(hook, varString, eltType);
  return tpl;
}

//=============================================================================
/* Recherche des pts cellN=2 */
//=============================================================================
PyObject* K_CONNECTOR::getInterpolatedPoints(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array))
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: wrong arguments.");
    return NULL;
  }
  // Check:
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  if (res != 1 && res != 2)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: invalid array.");
    return NULL;
  }
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: array must contain cellN.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posc++;

  /*fin verifs*/
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; varStringOut[0] = '\0';
  E_Int nfldOut = nfld+1;
  strcpy(varStringOut,varString); strcat(varStringOut,",indcell");
  E_Float* cellnp = f->begin(posc);

  FldArrayF* fout = new FldArrayF(npts, nfldOut);
  E_Int c=0;
  for (E_Int ind=0; ind < npts; ind++)
  {
    if (K_FUNC::fEqualZero(cellnp[ind] - 2.))
    {
      for (E_Int eq = 1; eq <= nfld; eq++) (*fout)(c,eq) = (*f)(ind,eq);
      (*fout)(c,nfldOut) = E_Float(ind);
      c++;
    }
  }
  fout->reAllocMat(c, nfldOut);

  RELEASESHAREDB(res, array, f, cn);
  FldArrayI* cnl = new FldArrayI(0);
  PyObject* tpl = K_ARRAY::buildArray3(*fout, varStringOut, *cnl, "NODE", api);
  delete fout; delete cnl;
  return tpl;
}
