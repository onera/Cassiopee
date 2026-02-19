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

# include <cstdio>
# include "loc.h"
# include <vector>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Convertit un array centres en array noeuds en structure
   Retourne 1 en cas de succes, 0 en cas d'echec.
   IN: FCenter: champ en centres
   IN: ni,nj,nk: dim de la grille en centres
   IN: cellN: position du champ cellN dans FCenter (-1 si il n'y en a pas)
   IN: mod: type de cellN
   IN: posx, posy, posz: position des coords dans FCenter
   OUT: FNode: champ aux noeuds
   OUT: nin, njn, nkn: dim de la grille en noeuds
   IN: algo: type de prise en compte du cellN.
   Si algo=0, des qu'un noeud d'une cellule a cellN=0, la cellule 
   recoit cellN=0.
   Si algo=1, une cellule recoit cellN=0 si tous ses noeuds ont cellN=0. 
   Le champ est alors pondere.
*/
//=============================================================================
E_Int K_LOC::center2nodeStruct(FldArrayF& FCenter, 
                               E_Int ni, E_Int nj, E_Int nk,
                               E_Int cellN, E_Int mod,
                               E_Int posx, E_Int posy, E_Int posz,
                               FldArrayF& FNode,
                               E_Int& nin, E_Int& njn, E_Int& nkn,
                               E_Int algo)
{
  E_Int nv = FCenter.getNfld();
  E_Int size = 0, dim = 3, im = 1;
  E_Int jm = 1; E_Int km = 1;

  if (ni == 1)
  {
    if ((nj != 1)&&(nk != 1))
    {
      size = (nj+1)*(nk+1); dim = 2;
      im = nj; jm = nk;
      nin = 1; njn = nj+1; nkn = nk+1;
    }
    else if (nj == 1)
    {
      size = (nk+1); dim = 1;
      im = nk;
      nin = 1; njn = 1; nkn = nk+1;
    }
    else if (nk == 1)
    {
      size = (nj+1); dim = 1;
      im = nj;
      nin = 1; njn = 1; nkn = nk+1;
    }
  }
  else if (nj == 1)
  {
    if ((ni != 1)&&(nk != 1))
    {
      size = (ni+1)*(nk+1); dim = 2;
      im = ni; jm = nk;
      nin = ni+1; njn = 1; nkn = nk+1;
    }
    else if (ni == 1)
    {
      size = (nk+1); dim = 1;
      im = nk;
      nin = 1; njn = 1; nkn = nk+1;
    }
    else if (nk == 1)
    {
      size = (ni+1); dim = 1;
      im = ni;
      nin = ni+1; njn = 1; nkn = 1;
    }
  }
  else if (nk == 1)
  {
    if ((ni != 1)&&(nj != 1))
    {
      size = (ni+1)*(nj+1); dim = 2;
      im = ni; jm = nj;
      nin = ni+1; njn = nj+1; nkn = 1;
    }
    else if (ni == 1)
    {
      size = (nj+1); dim = 1;
      im = nj;
      nin = 1; njn = nj+1; nkn = 1;
    }
    else if (nj == 1)
    {
      size = (ni+1); dim = 1;
      im = ni;
      nin = ni+1; njn = 1; nkn = 1;
    }
  }
  else
  {
    size = (ni+1)*(nj+1)*(nk+1); dim = 3;
    im = ni; jm = nj; km = nk;
    nin = ni+1; njn = nj+1; nkn = nk+1;
  }
 
  // On ne realloue que si FNode n'est pas alloue correctement
  if (FNode.getSize() != size || FNode.getNfld() != nv) FNode.malloc(size, nv);
    
  // Converti les champs
  if (algo == 0 || cellN == -1)
  {
    if (dim == 1)
    {
      #pragma omp parallel
      {
        E_Int alpha, ind0, ind1;
        #pragma omp for 
        for (E_Int i = 0; i <= im; i++)
        {
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          ind0 = E_max(i, 1)-1;
          ind1 = ind0 + alpha;
          for (E_Int n = 1; n <= nv; n++)
            FNode(i,n) = 0.5*(FCenter(ind0,n)+FCenter(ind1,n));
        }
      }
    }
    else if (dim == 2)
    {
      E_Int im1 = im+1; E_Int jm1 = jm+1;
      #pragma omp parallel
      {
        E_Int alpha, beta, i, j, i0, j0, ind0, ind1, ind2, ind3;
        #pragma omp for
        for (E_Int ind = 0; ind < im1*jm1; ind++)
        {
          j = ind / im1;
          i = ind - j*im1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;    
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = E_max(j, 1)-1;

          ind0 = i0 + j0 * im;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;

          for (E_Int n = 1; n <= nv; n++)
            FNode(ind,n) = 0.25*(FCenter(ind0,n)+FCenter(ind1,n)+FCenter(ind2,n)+FCenter(ind3,n));
        }
      }
    }
    else
    {
      E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
      E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
      
      #pragma omp parallel
      {
        E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
        E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
        #pragma omp for
        for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
        {
          k = ind / ijm1;
          j = (ind - k*ijm1) / im1;
          i = ind - j*im1 - k*ijm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = E_max(j, 1)-1;
          gamma = ijm;
          if (k == 0 || k == km) gamma = 0;
          k0 = E_max(k, 1)-1;
          
          ind0 = i0 + j0 * im + k0 * im*jm;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          ind4 = ind0 + gamma;
          ind5 = ind4 + alpha;
          ind6 = ind4 + beta;
          ind7 = ind6 + alpha;

          for (E_Int n = 1; n <= nv; n++)
            FNode(ind,n) = 0.125*(FCenter(ind0,n)+FCenter(ind1,n)+
              FCenter(ind2,n)+FCenter(ind3,n)+
              FCenter(ind4,n)+FCenter(ind5,n)+
              FCenter(ind6,n)+FCenter(ind7,n));
        }
      }
    }
  }
  else // algo=1 et cellN existe
  {
    if (dim == 1)
    {
      E_Int im1 = im+1;
      E_Float* cellNp = FCenter.begin(cellN);
      
      #pragma omp parallel
      {
        E_Int alpha, i, i0;
        E_Int ind0, ind1;
        E_Float cellN0, cellN1, w;
        #pragma omp for
        for (E_Int ind = 0; ind < im1; ind++)
        {
          i = ind;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;
                    
          ind0 = i0;
          ind1 = ind0 + alpha;
          
          cellN0 = E_min(cellNp[ind0], 1.);
          cellN1 = E_min(cellNp[ind1], 1.);
               
          w = cellN0 + cellN1;
          if (K_FUNC::fEqualZero(w))
          {
            w = 0.5; cellN0 = 1.; cellN1 = 1.;
          }
          else w = 1./w;
               
          for (E_Int n = 1; n <= nv; n++)
            FNode(ind,n) = w*(cellN0*FCenter(ind0,n)+cellN1*FCenter(ind1,n));
        }
      }
    }
    else if (dim == 2)
    {  
      E_Int im1 = im+1; E_Int jm1 = jm+1;
      E_Float* cellNp = FCenter.begin(cellN);
      
      #pragma omp parallel
      {
        E_Int alpha, beta, i, j, i0, j0;
        E_Int ind0, ind1, ind2, ind3;
        E_Float cellN0, cellN1, cellN2, cellN3, w;
        #pragma omp for
        for (E_Int ind = 0; ind < im1*jm1; ind++)
        {
          j = ind / im1;
          i = ind - j*im1;

          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = E_max(j, 1)-1;
                    
          ind0 = i0 + j0 * im;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          
          cellN0 = E_min(cellNp[ind0], 1.);
          cellN1 = E_min(cellNp[ind1], 1.);
          cellN2 = E_min(cellNp[ind2], 1.);
          cellN3 = E_min(cellNp[ind3], 1.);
               
          w = cellN0 + cellN1 + cellN2 + cellN3;
          if (K_FUNC::fEqualZero(w))
          {
            w = 0.5; cellN0 = 1.; cellN1 = 1.; cellN2 = 1.;
            cellN3 = 1.; 
          }
          else w = 1./w;
               
          for (E_Int n = 1; n <= nv; n++)
            FNode(ind,n) = w*(cellN0*FCenter(ind0,n)+cellN1*FCenter(ind1,n)+
              cellN2*FCenter(ind2,n)+cellN3*FCenter(ind3,n));
        }
      }
    }
    else
    {
      E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
      E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
      E_Float* cellNp = FCenter.begin(cellN);
      
      #pragma omp parallel
      {
        E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
        E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
        E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
        E_Float cellN6, cellN7, w;
        #pragma omp for
        for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
        {
          k = ind / ijm1;
          j = (ind - k*ijm1) / im1;
          i = ind - j*im1 - k*ijm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = E_max(j, 1)-1;
          gamma = ijm;
          if (k == 0 || k == km) gamma = 0;
          k0 = E_max(k, 1)-1;
          
          ind0 = i0 + j0 * im + k0 * im*jm;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          ind4 = ind0 + gamma;
          ind5 = ind4 + alpha;
          ind6 = ind4 + beta;
          ind7 = ind6 + alpha;

          cellN0 = E_min(cellNp[ind0], 1.);
          cellN1 = E_min(cellNp[ind1], 1.);
          cellN2 = E_min(cellNp[ind2], 1.);
          cellN3 = E_min(cellNp[ind3], 1.);
          cellN4 = E_min(cellNp[ind4], 1.);
          cellN5 = E_min(cellNp[ind5], 1.);
          cellN6 = E_min(cellNp[ind6], 1.);
          cellN7 = E_min(cellNp[ind7], 1.);
               
          w = cellN0 + cellN1 + cellN2 + cellN3 + cellN4 + 
              cellN5 + cellN6 + cellN7;
          if (K_FUNC::fEqualZero(w))
          {
            w = 0.125; cellN0 = 1.; cellN1 = 1.; cellN2 = 1.;
            cellN3 = 1.; cellN4 = 1.; cellN5 = 1.; cellN6 = 1.;
            cellN7 = 1.;
          }
          else w = 1./w;
               
          for (E_Int n = 1; n <= nv; n++)
            FNode(ind,n) = w*(cellN0*FCenter(ind0,n)+cellN1*FCenter(ind1,n)+
              cellN2*FCenter(ind2,n)+cellN3*FCenter(ind3,n)+
              cellN4*FCenter(ind4,n)+cellN5*FCenter(ind5,n)+
              cellN6*FCenter(ind6,n)+cellN7*FCenter(ind7,n));
        }
      }
    }
  }

  // Traitement special pour les coords
  if (posx != -1 && posy != -1 && posz != -1)
  {
    posx++; posy++; posz++;
    if (dim == 1)
    {
      E_Float* Fnx = FNode.begin(posx);
      E_Float* Fny = FNode.begin(posy);
      E_Float* Fnz = FNode.begin(posz);

      E_Float* Fcx = FCenter.begin(posx);
      E_Float* Fcy = FCenter.begin(posy);
      E_Float* Fcz = FCenter.begin(posz);
      
      E_Int im1 = im+1;
      #pragma omp parallel
      {
        E_Int alpha, ind, ind0, ind1;
        // main loop
        #pragma omp for 
        for (E_Int i = 0; i <= im; i++)
        {
          ind = i;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          ind0 = E_max(i, 1)-1;
          ind1 = ind0 + alpha;
          Fnx[ind] = 0.5*(Fcx[ind0]+Fcx[ind1]);
          Fny[ind] = 0.5*(Fcy[ind0]+Fcy[ind1]);
          Fnz[ind] = 0.5*(Fcz[ind0]+Fcz[ind1]);
        }
        // i = 0
        Fnx[0] = 2.*Fnx[1]-Fnx[2];
        Fny[0] = 2.*Fny[1]-Fny[2];
        Fnz[0] = 2.*Fnz[1]-Fnz[2];
        // i = im
        Fnx[im1-1] = 2.*Fnx[im1-2]-Fnx[im1-3];
        Fny[im1-1] = 2.*Fny[im1-2]-Fny[im1-3];
        Fnz[im1-1] = 2.*Fnz[im1-2]-Fnz[im1-3];
      }
    }
    else if (dim == 2)
    {
      E_Float* Fnx = FNode.begin(posx);
      E_Float* Fny = FNode.begin(posy);
      E_Float* Fnz = FNode.begin(posz);

      E_Float* Fcx = FCenter.begin(posx);
      E_Float* Fcy = FCenter.begin(posy);
      E_Float* Fcz = FCenter.begin(posz);
      
      E_Int im1 = im+1; E_Int jm1 = jm+1;
      #pragma omp parallel
      {
        E_Int ind0, ind1, ind2, ind3;
        E_Int i, j, i0, j0;
        E_Int alpha, beta;
        // main loop
        #pragma omp for
        for (E_Int ind = 0; ind < im1*jm1; ind++)
        {
          j = ind / im1;
          i = ind - j*im1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;    
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = E_max(j, 1)-1;

          ind0 = i0 + j0 * im;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;

          Fnx[ind] = 0.25*(Fcx[ind0]+Fcx[ind1]+Fcx[ind2]+Fcx[ind3]);
          Fny[ind] = 0.25*(Fcy[ind0]+Fcy[ind1]+Fcy[ind2]+Fcy[ind3]);
          Fnz[ind] = 0.25*(Fcz[ind0]+Fcz[ind1]+Fcz[ind2]+Fcz[ind3]);
        }
        // edges imin, imax
        #pragma omp for
        for (E_Int j = 1; j < jm; j++)
        {
          // i = 0
          ind0 = j*im1;
          ind2 = ind0+1;
          ind3 = ind2+1;
          Fnx[ind0] = 2.*Fnx[ind2]-Fnx[ind3];
          Fny[ind0] = 2.*Fny[ind2]-Fny[ind3];
          Fnz[ind0] = 2.*Fnz[ind2]-Fnz[ind3];

          // i = im
          ind0 = im + j*im1;
          ind2 = ind0-1;
          ind3 = ind2-1;
          Fnx[ind0] = 2.*Fnx[ind2]-Fnx[ind3];
          Fny[ind0] = 2.*Fny[ind2]-Fny[ind3];
          Fnz[ind0] = 2.*Fnz[ind2]-Fnz[ind3];
        }
        // edges imin, imax
        #pragma omp for
        for (E_Int i = 1; i < im; i++)
        {
          // j = 0
          ind0 = i;
          ind2 = ind0+im1;
          ind3 = ind2+im1;
          Fnx[ind0] = 2.*Fnx[ind2]-Fnx[ind3];
          Fny[ind0] = 2.*Fny[ind2]-Fny[ind3];
          Fnz[ind0] = 2.*Fnz[ind2]-Fnz[ind3];

          // j = jm
          ind0 = i + jm*im1;
          ind2 = ind0-im1;
          ind3 = ind2-im1;
          Fnx[ind0] = 2.*Fnx[ind2]-Fnx[ind3];
          Fny[ind0] = 2.*Fny[ind2]-Fny[ind3];
          Fnz[ind0] = 2.*Fnz[ind2]-Fnz[ind3];
        }
        // i = 0, j = 0
        Fnx[0] = 2.*Fnx[1]-Fnx[2];
        Fny[0] = 2.*Fny[1]-Fny[2];
        Fnz[0] = 2.*Fnz[1]-Fnz[2];
        // i = im1, j = 0
        Fnx[im1-1] = 2.*Fnx[im1-2]-Fnx[im1-3];
        Fny[im1-1] = 2.*Fny[im1-2]-Fny[im1-3];
        Fnz[im1-1] = 2.*Fnz[im1-2]-Fnz[im1-3];
        // i = 0, j = jm1
        Fnx[(jm1-1)*im1] = 2.*Fnx[1+(jm1-1)*im1]-Fnx[2+(jm1-1)*im1];
        Fny[(jm1-1)*im1] = 2.*Fny[1+(jm1-1)*im1]-Fny[2+(jm1-1)*im1];
        Fnz[(jm1-1)*im1] = 2.*Fnz[1+(jm1-1)*im1]-Fnz[2+(jm1-1)*im1];
        // i = im1, j = jm1
        Fnx[(im1-1)+(jm1-1)*im1] = 2.*Fnx[(im1-2)+(jm1-1)*im1]-Fnx[(im1-3)+(jm1-1)*im1];
        Fny[(im1-1)+(jm1-1)*im1] = 2.*Fny[(im1-2)+(jm1-1)*im1]-Fny[(im1-3)+(jm1-1)*im1];
        Fnz[(im1-1)+(jm1-1)*im1] = 2.*Fnz[(im1-2)+(jm1-1)*im1]-Fnz[(im1-3)+(jm1-1)*im1];
      }
    }
    else
    {
      E_Float* Fnx = FNode.begin(posx);
      E_Float* Fny = FNode.begin(posy);
      E_Float* Fnz = FNode.begin(posz);

      E_Float* Fcx = FCenter.begin(posx);
      E_Float* Fcy = FCenter.begin(posy);
      E_Float* Fcz = FCenter.begin(posz);

      E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
      E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
      
      #pragma omp parallel
      {
        E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
        E_Int ind, ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
        // main loop
        #pragma omp for
        for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
        {
          k = ind / ijm1;
          j = (ind - k*ijm1) / im1;
          i = ind - j*im1 - k*ijm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = E_max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = E_max(j, 1)-1;
          gamma = ijm;
          if (k == 0 || k == km) gamma = 0;
          k0 = E_max(k, 1)-1;
          
          ind0 = i0 + j0 * im + k0 * im*jm;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          ind4 = ind0 + gamma;
          ind5 = ind4 + alpha;
          ind6 = ind4 + beta;
          ind7 = ind6 + alpha;

          Fnx[ind] = 0.125*(Fcx[ind0]+Fcx[ind1]+Fcx[ind2]+Fcx[ind3]+Fcx[ind4]+Fcx[ind5]+Fcx[ind6]+Fcx[ind7]);
          Fny[ind] = 0.125*(Fcy[ind0]+Fcy[ind1]+Fcy[ind2]+Fcy[ind3]+Fcy[ind4]+Fcy[ind5]+Fcy[ind6]+Fcy[ind7]);
          Fnz[ind] = 0.125*(Fcz[ind0]+Fcz[ind1]+Fcz[ind2]+Fcz[ind3]+Fcz[ind4]+Fcz[ind5]+Fcz[ind6]+Fcz[ind7]);
        }
        //faces imin & imax
        #pragma omp for collapse(2)
        for (E_Int k = 1; k < km1-1; k++)
        {
          for (E_Int j = 1; j < jm1-1; j++)
          {
            // imin --------------------------------------------------
            i = 0;

            i0 = E_max(i, 1)-1;
            j0 = E_max(j, 1)-1;
            k0 = E_max(k, 1)-1;

            alpha = 1;
            beta  = im;
            gamma = ijm;
            
            ind = i + j*im1 + k*ijm1;

            // face i
            ind0 = i0 + j0*im + k0*ijm;
            ind2 = ind0 + beta;
            ind4 = ind0 + gamma;
            ind6 = ind4 + beta;

            // face i+1
            ind1 = ind0 + alpha;
            ind3 = ind1 + beta;
            ind5 = ind1 + gamma;
            ind7 = ind5 + beta; 

            Fnx[ind] = 0.125*(3*(Fcx[ind0]+Fcx[ind2]+Fcx[ind4]+Fcx[ind6])-(Fcx[ind1]+Fcx[ind3]+Fcx[ind5]+Fcx[ind7]));
            Fny[ind] = 0.125*(3*(Fcy[ind0]+Fcy[ind2]+Fcy[ind4]+Fcy[ind6])-(Fcy[ind1]+Fcy[ind3]+Fcy[ind5]+Fcy[ind7]));
            Fnz[ind] = 0.125*(3*(Fcz[ind0]+Fcz[ind2]+Fcz[ind4]+Fcz[ind6])-(Fcz[ind1]+Fcz[ind3]+Fcz[ind5]+Fcz[ind7]));

            // imax --------------------------------------------------
            i = im1-1;

            i0 = E_max(i, 1)-1;
            j0 = E_max(j, 1)-1;
            k0 = E_max(k, 1)-1;

            alpha = -1;
            beta  = im;
            gamma = ijm;
            
            ind = i + j*im1 + k*ijm1;

            // face i
            ind0 = i0 + j0*im + k0*ijm;
            ind2 = ind0 + beta;
            ind4 = ind0 + gamma;
            ind6 = ind4 + beta;

            // face i-1
            ind1 = ind0 + alpha;
            ind3 = ind1 + beta;
            ind5 = ind1 + gamma;
            ind7 = ind5 + beta; 

            Fnx[ind] = 0.125*(3*(Fcx[ind0]+Fcx[ind2]+Fcx[ind4]+Fcx[ind6])-(Fcx[ind1]+Fcx[ind3]+Fcx[ind5]+Fcx[ind7]));
            Fny[ind] = 0.125*(3*(Fcy[ind0]+Fcy[ind2]+Fcy[ind4]+Fcy[ind6])-(Fcy[ind1]+Fcy[ind3]+Fcy[ind5]+Fcy[ind7]));
            Fnz[ind] = 0.125*(3*(Fcz[ind0]+Fcz[ind2]+Fcz[ind4]+Fcz[ind6])-(Fcz[ind1]+Fcz[ind3]+Fcz[ind5]+Fcz[ind7]));
          }
        }
        //faces jmin & jmax
        #pragma omp for collapse(2)
        for (E_Int k = 1; k < km1-1; k++)
        {
          for (E_Int i = 1; i < im1-1; i++)
          {
            // jmin --------------------------------------------------
            j = 0;

            i0 = E_max(i, 1)-1;
            j0 = E_max(j, 1)-1;
            k0 = E_max(k, 1)-1;

            alpha = 1;
            beta  = im;
            gamma = ijm;
            
            ind = i + j*im1 + k*ijm1;

            // face j
            ind0 = i0 + j0*im + k0*ijm;
            ind2 = ind0 + alpha;
            ind4 = ind0 + gamma;
            ind6 = ind4 + alpha;

            // face j+1
            ind1 = ind0 + beta;
            ind3 = ind1 + alpha;
            ind5 = ind1 + gamma;
            ind7 = ind5 + alpha; 

            Fnx[ind] = 0.125*(3*(Fcx[ind0]+Fcx[ind2]+Fcx[ind4]+Fcx[ind6])-(Fcx[ind1]+Fcx[ind3]+Fcx[ind5]+Fcx[ind7]));
            Fny[ind] = 0.125*(3*(Fcy[ind0]+Fcy[ind2]+Fcy[ind4]+Fcy[ind6])-(Fcy[ind1]+Fcy[ind3]+Fcy[ind5]+Fcy[ind7]));
            Fnz[ind] = 0.125*(3*(Fcz[ind0]+Fcz[ind2]+Fcz[ind4]+Fcz[ind6])-(Fcz[ind1]+Fcz[ind3]+Fcz[ind5]+Fcz[ind7]));

            // jmax --------------------------------------------------
            j = jm1-1;

            i0 = E_max(i, 1)-1;
            j0 = E_max(j, 1)-1;
            k0 = E_max(k, 1)-1;

            alpha = 1;
            beta  = -im;
            gamma = ijm;
            
            ind = i + j*im1 + k*ijm1;

            // face j
            ind0 = i0 + j0*im + k0*ijm;
            ind2 = ind0 + alpha;
            ind4 = ind0 + gamma;
            ind6 = ind4 + alpha;

            // face j-1
            ind1 = ind0 + beta;
            ind3 = ind1 + alpha;
            ind5 = ind1 + gamma;
            ind7 = ind5 + alpha; 

            Fnx[ind] = 0.125*(3*(Fcx[ind0]+Fcx[ind2]+Fcx[ind4]+Fcx[ind6])-(Fcx[ind1]+Fcx[ind3]+Fcx[ind5]+Fcx[ind7]));
            Fny[ind] = 0.125*(3*(Fcy[ind0]+Fcy[ind2]+Fcy[ind4]+Fcy[ind6])-(Fcy[ind1]+Fcy[ind3]+Fcy[ind5]+Fcy[ind7]));
            Fnz[ind] = 0.125*(3*(Fcz[ind0]+Fcz[ind2]+Fcz[ind4]+Fcz[ind6])-(Fcz[ind1]+Fcz[ind3]+Fcz[ind5]+Fcz[ind7]));
          }
        }
        //faces kmin & kmax
        #pragma omp for collapse(2)
        for (E_Int j = 1; j < jm1-1; j++)
        {
          for (E_Int i = 1; i < im1-1; i++)
          {
            // kmin --------------------------------------------------
            k = 0;

            i0 = E_max(i, 1)-1;
            j0 = E_max(j, 1)-1;
            k0 = E_max(k, 1)-1;

            alpha = 1;
            beta  = im;
            gamma = ijm;
            
            ind = i + j*im1 + k*ijm1;

            // face k
            ind0 = i0 + j0*im + k0*ijm;
            ind2 = ind0 + alpha;
            ind4 = ind0 + beta;
            ind6 = ind4 + alpha;

            // face k+1
            ind1 = ind0 + gamma;
            ind3 = ind1 + alpha;
            ind5 = ind1 + beta;
            ind7 = ind5 + alpha; 

            Fnx[ind] = 0.125*(3*(Fcx[ind0]+Fcx[ind2]+Fcx[ind4]+Fcx[ind6])-(Fcx[ind1]+Fcx[ind3]+Fcx[ind5]+Fcx[ind7]));
            Fny[ind] = 0.125*(3*(Fcy[ind0]+Fcy[ind2]+Fcy[ind4]+Fcy[ind6])-(Fcy[ind1]+Fcy[ind3]+Fcy[ind5]+Fcy[ind7]));
            Fnz[ind] = 0.125*(3*(Fcz[ind0]+Fcz[ind2]+Fcz[ind4]+Fcz[ind6])-(Fcz[ind1]+Fcz[ind3]+Fcz[ind5]+Fcz[ind7]));

            // kmax --------------------------------------------------
            k = km1-1;

            i0 = E_max(i, 1)-1;
            j0 = E_max(j, 1)-1;
            k0 = E_max(k, 1)-1;

            alpha = 1;
            beta  = im;
            gamma = -ijm;
            
            ind = i + j*im1 + k*ijm1;

            // face k
            ind0 = i0 + j0*im + k0*ijm;
            ind2 = ind0 + alpha;
            ind4 = ind0 + beta;
            ind6 = ind4 + alpha;

            // face k-1
            ind1 = ind0 + gamma;
            ind3 = ind1 + alpha;
            ind5 = ind1 + beta;
            ind7 = ind5 + alpha; 

            Fnx[ind] = 0.125*(3*(Fcx[ind0]+Fcx[ind2]+Fcx[ind4]+Fcx[ind6])-(Fcx[ind1]+Fcx[ind3]+Fcx[ind5]+Fcx[ind7]));
            Fny[ind] = 0.125*(3*(Fcy[ind0]+Fcy[ind2]+Fcy[ind4]+Fcy[ind6])-(Fcy[ind1]+Fcy[ind3]+Fcy[ind5]+Fcy[ind7]));
            Fnz[ind] = 0.125*(3*(Fcz[ind0]+Fcz[ind2]+Fcz[ind4]+Fcz[ind6])-(Fcz[ind1]+Fcz[ind3]+Fcz[ind5]+Fcz[ind7]));
          }
        }
        //edges jmin/kmin, jmin/kmax, jmax/kmin, jmax/kmax 
        #pragma omp for
        for (E_Int i = 1; i < im1-1; i++)
        {
          // jmin&kmin
          j = 0;
          k = 0;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // j->k->i
          beta  = im;
          gamma = ijm;
          alpha = 1;

          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // jmin&kmax
          j = 0;
          k = km1-1;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // j->(-k)->i
          beta  = im;
          gamma = -ijm;
          alpha = 1;

          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // jmax&kmin
          j = jm1-1;
          k = 0;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // (-j)->k->i
          beta  = -im;
          gamma = ijm;
          alpha = 1;

          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // jmax&kmax
          j = jm1-1;
          k = km1-1;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // (-j)->(-k)->i
          beta  = -im;
          gamma = -ijm;
          alpha = 1;

          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

        }
        //edges imin/kmin, imin/kmax, imax/kmin, imax/kmax 
        #pragma omp for
        for (E_Int j = 1; j < jm1-1; j++)
        {
          // imin&kmin
          i = 0;
          k = 0;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // i->k->j
          beta  = 1;
          gamma = ijm;
          alpha = im;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // imin&kmax
          i = 0;
          k = km1-1;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // i->(-k)->j
          beta  = 1;
          gamma = -ijm;
          alpha = im;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // imax&kmin
          i = im1-1;
          k = 0;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // (-i)->k->j
          beta  = -1;
          gamma = ijm;
          alpha = im;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // imax&kmax
          i = im1-1;
          k = km1-1;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // (-i)->(-k)->j
          beta  = -1;
          gamma = -ijm;
          alpha = im;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));
        }
        //edges imin/jmin, imin/jmax, imax/jmin, imax/jmax 
        #pragma omp for
        for (E_Int k = 1; k < km1-1; k++)
        {
          // imin&jmin
          i = 0;
          j = 0;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // i->j->k
          beta  = 1;
          gamma = im;
          alpha = ijm;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // imin&jmax
          i = 0;
          j = jm1-1;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // i->(-j)->k
          beta  = 1;
          gamma = -im;
          alpha = ijm;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // imax&jmin
          i = im1-1;
          j = 0;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // (-i)->j->k
          beta  = -1;
          gamma = im;
          alpha = ijm;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));

          // imax&jmax
          i = im1-1;
          j = jm1-1;
          
          i0 = E_max(i, 1)-1;
          j0 = E_max(j, 1)-1;
          k0 = E_max(k, 1)-1;

          // (-i)->(-j)->k
          beta  = -1;
          gamma = -im;
          alpha = ijm;
          
          ind = i + j*im1 + k*ijm1;

          ind0 = i0 + j0*im + k0*ijm;
          ind2 = ind0 + beta;
          ind4 = ind0 + gamma;
          ind6 = ind4 + beta;

          ind1 = ind0 + alpha;
          ind3 = ind1 + beta;
          ind5 = ind1 + gamma;
          ind7 = ind5 + beta;

          Fnx[ind] = 0.125*(3*(3*Fcx[ind0]-Fcx[ind2])-(3*Fcx[ind4]-Fcx[ind6]) + 3*(3*Fcx[ind1]-Fcx[ind3])-(3*Fcx[ind5]-Fcx[ind7]));
          Fny[ind] = 0.125*(3*(3*Fcy[ind0]-Fcy[ind2])-(3*Fcy[ind4]-Fcy[ind6]) + 3*(3*Fcy[ind1]-Fcy[ind3])-(3*Fcy[ind5]-Fcy[ind7]));
          Fnz[ind] = 0.125*(3*(3*Fcz[ind0]-Fcz[ind2])-(3*Fcz[ind4]-Fcz[ind6]) + 3*(3*Fcz[ind1]-Fcz[ind3])-(3*Fcz[ind5]-Fcz[ind7]));
        }
        // i = 0, j = 0, k = 0
        Fnx[0] = 2*Fnx[im1]-Fnx[2*im1];
        Fny[0] = 2*Fny[im1]-Fny[2*im1];
        Fnz[0] = 2*Fnz[im1]-Fnz[2*im1];
        // i = 0, j = 0, k = km1-1
        Fnx[0+(km1-1)*ijm1] = 2*Fnx[im1+(km1-1)*ijm1]-Fnx[2*im1+(km1-1)*ijm1];
        Fny[0+(km1-1)*ijm1] = 2*Fny[im1+(km1-1)*ijm1]-Fny[2*im1+(km1-1)*ijm1];
        Fnz[0+(km1-1)*ijm1] = 2*Fnz[im1+(km1-1)*ijm1]-Fnz[2*im1+(km1-1)*ijm1];
        // i = 0, j = im1-1, k = 0
        Fnx[(jm1-1)*im1] = 2*Fnx[(jm1-2)*im1]-Fnx[(jm1-3)*im1];
        Fny[(jm1-1)*im1] = 2*Fny[(jm1-2)*im1]-Fny[(jm1-3)*im1];
        Fnz[(jm1-1)*im1] = 2*Fnz[(jm1-2)*im1]-Fnz[(jm1-3)*im1];
        // i = 0, j = im1-1, k = km1-1
        Fnx[(jm1-1)*im1+(km1-1)*ijm1] = 2*Fnx[(jm1-2)*im1+(km1-1)*ijm1]-Fnx[(jm1-3)*im1+(km1-1)*ijm1];
        Fny[(jm1-1)*im1+(km1-1)*ijm1] = 2*Fny[(jm1-2)*im1+(km1-1)*ijm1]-Fny[(jm1-3)*im1+(km1-1)*ijm1];
        Fnz[(jm1-1)*im1+(km1-1)*ijm1] = 2*Fnz[(jm1-2)*im1+(km1-1)*ijm1]-Fnz[(jm1-3)*im1+(km1-1)*ijm1];
        // i = im1-1, j = 0, k = 0
        Fnx[(im1-1)+0] = 2*Fnx[(im1-1)+im1]-Fnx[(im1-1)+2*im1];
        Fny[(im1-1)+0] = 2*Fny[(im1-1)+im1]-Fny[(im1-1)+2*im1];
        Fnz[(im1-1)+0] = 2*Fnz[(im1-1)+im1]-Fnz[(im1-1)+2*im1];
        // i = im1-1, j = 0, k = km1-1
        Fnx[(im1-1)+0+(km1-1)*ijm1] = 2*Fnx[(im1-1)+im1+(km1-1)*ijm1]-Fnx[(im1-1)+2*im1+(km1-1)*ijm1];
        Fny[(im1-1)+0+(km1-1)*ijm1] = 2*Fny[(im1-1)+im1+(km1-1)*ijm1]-Fny[(im1-1)+2*im1+(km1-1)*ijm1];
        Fnz[(im1-1)+0+(km1-1)*ijm1] = 2*Fnz[(im1-1)+im1+(km1-1)*ijm1]-Fnz[(im1-1)+2*im1+(km1-1)*ijm1];
        // i = im1-1, j = im1-1, k = 0
        Fnx[(im1-1)+(jm1-1)*im1] = 2*Fnx[(im1-1)+(jm1-2)*im1]-Fnx[(im1-1)+(jm1-3)*im1];
        Fny[(im1-1)+(jm1-1)*im1] = 2*Fny[(im1-1)+(jm1-2)*im1]-Fny[(im1-1)+(jm1-3)*im1];
        Fnz[(im1-1)+(jm1-1)*im1] = 2*Fnz[(im1-1)+(jm1-2)*im1]-Fnz[(im1-1)+(jm1-3)*im1];
        // i = im1-1, j = im1-1, k = km1-1
        Fnx[(im1-1)+(jm1-1)*im1+(km1-1)*ijm1] = 2*Fnx[(im1-1)+(jm1-2)*im1+(km1-1)*ijm1]-Fnx[(im1-1)+(jm1-3)*im1+(km1-1)*ijm1];
        Fny[(im1-1)+(jm1-1)*im1+(km1-1)*ijm1] = 2*Fny[(im1-1)+(jm1-2)*im1+(km1-1)*ijm1]-Fny[(im1-1)+(jm1-3)*im1+(km1-1)*ijm1];
        Fnz[(im1-1)+(jm1-1)*im1+(km1-1)*ijm1] = 2*Fnz[(im1-1)+(jm1-2)*im1+(km1-1)*ijm1]-Fnz[(im1-1)+(jm1-3)*im1+(km1-1)*ijm1];
      }
    }
  }
  
  // Traitement special pour le "cellnaturefield"
  // printf("cellN = %i, algo = %i, mod = %i\n", cellN, algo, mod);
  if (cellN != -1 && algo == 0)
  {
    switch (mod)
    {
      case 1:
        if (dim == 1)
        {
          E_Int im1 = im+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, i, i0;
            E_Int ind0, ind1;
            E_Float cellN0, cellN1, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1; ind++)
            {
              i = ind;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
                        
              ind0 = i0;
              ind1 = ind0 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.); 

              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else if (dim == 2)
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, i, j, i0, j0;
            E_Int ind0, ind1, ind2, ind3;
            E_Float cellN0, cellN1, cellN2, cellN3, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1; ind++)
            {
              j = ind / im1;
              i = ind - j*im1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
                        
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              
              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
                  
              w = cellN0*cellN1*cellN2*cellN3;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
          E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
            E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
            E_Float cellN6, cellN7, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
            {
              k = ind / ijm1;
              j = (ind - k*ijm1) / im1;
              i = ind - j*im1 - k*ijm1;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              gamma = ijm;
              if (k == 0 || k == km) gamma = 0;
              k0 = E_max(k, 1)-1;
              
              ind0 = i0 + j0 * im + k0 * im*jm;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              ind4 = ind0 + gamma;
              ind5 = ind4 + alpha;
              ind6 = ind4 + beta;
              ind7 = ind6 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
              cellN4 = E_min(cellNp[ind4], 1.);
              cellN5 = E_min(cellNp[ind5], 1.);
              cellN6 = E_min(cellNp[ind6], 1.);
              cellN7 = E_min(cellNp[ind7], 1.);
                  
              w = cellN0*cellN1*cellN2*cellN3*cellN4*cellN5*cellN6*cellN7;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        break;
        
      case 2:
        if (dim == 1)
        {
          E_Int im1 = im+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, i, i0;
            E_Int ind0, ind1;
            E_Float cellN0, cellN1, w;
            #pragma omp for
            for (E_Int ind = 1; ind < im1-1; ind++)
            {
              i = ind;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
                        
              ind0 = i0;
              ind1 = ind0 + alpha;

              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1]; 

              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(cellN0 - 2.) || K_FUNC::fEqualZero(cellN1 - 2.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            // i = 0
            cellNpn[0] = cellNp[0];
            // i = im
            cellNpn[im1-1] = cellNp[im-1];
          }
        }
        else if (dim == 2)
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          #pragma omp parallel
          {
            E_Int alpha, beta, i, j, i0, j0;
            E_Int ind, ind0, ind1, ind2, ind3;
            E_Float cellN0, cellN1, cellN2, cellN3, w;
            // main loop
            #pragma omp for
            for (E_Int ind = 1; ind < im1*jm1-1; ind++)
            {
              j = ind / im1;
              i = ind - j*im1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
                        
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
              cellN2 = cellNp[ind2];
              cellN3 = cellNp[ind3];
                  
              w = cellN0*cellN1*cellN2*cellN3;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(cellN0 - 2.) || K_FUNC::fEqualZero(cellN1 - 2.) || 
                    K_FUNC::fEqualZero(cellN2 - 2.) || K_FUNC::fEqualZero(cellN3 - 2.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            // edges imin, imax
            #pragma omp for
            for (E_Int j = 0; j < jm1; j++)
            {
              // i = 0
              i = 0;
              i0 = E_max(i, 1)-1;

              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + beta;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (w == 4.)
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }

              // i = im1-1
              i = im1-1;
              i0 = E_max(i, 1)-1;

              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + beta;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (w == 4.)
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            // edges jmin, jmax
            #pragma omp for
            for (E_Int i = 0; i < im1; i++)
            {
              // j = 0
              j = 0;
              j0 = E_max(j, 1)-1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (w == 4.)
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }

              // j = jm1-1
              j = jm1-1;
              j0 = E_max(j, 1)-1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (w == 4.)
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
          }
        }
        else
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
          E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
            E_Int ind, ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
            E_Float cellN6, cellN7, w, w2;
            // main loop
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
            {
              k = ind / ijm1;
              j = (ind - k*ijm1) / im1;
              i = ind - j*im1 - k*ijm1;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              gamma = ijm;
              if (k == 0 || k == km) gamma = 0;
              k0 = E_max(k, 1)-1;
              
              ind0 = i0 + j0 * im + k0 * im*jm;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              ind4 = ind0 + gamma;
              ind5 = ind4 + alpha;
              ind6 = ind4 + beta;
              ind7 = ind6 + alpha;

              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
              cellN2 = cellNp[ind2];
              cellN3 = cellNp[ind3];
              cellN4 = cellNp[ind4];
              cellN5 = cellNp[ind5];
              cellN6 = cellNp[ind6];
              cellN7 = cellNp[ind7];
                  
              w = cellN0*cellN1*cellN2*cellN3*cellN4*cellN5*cellN6*cellN7;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                w2 = E_max(cellN0,1.)+E_max(cellN1,1.)+E_max(cellN2,1.)+E_max(cellN3,1.)+E_max(cellN4,1.)+E_max(cellN5,1.)+E_max(cellN6,1.)+E_max(cellN7,1.);
                if (w2 > 8.5)
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            //faces imin & imax
            #pragma omp for collapse(2)
            for (E_Int k = 0; k < km1; k++)
            {
              for (E_Int j = 0; j < jm1; j++)
              {
                // imin --------------------------------------------------
                i = 0;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + beta;
                ind2 = ind0 + gamma;
                ind3 = ind2 + beta;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0*cellN1*cellN2*cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 16.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }

                // imax --------------------------------------------------
                i = im1-1;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + beta;
                ind2 = ind0 + gamma;
                ind3 = ind2 + beta;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0*cellN1*cellN2*cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 16.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }
              }
            }
            //faces jmin & jmax
            #pragma omp for collapse(2)
            for (E_Int k = 0; k < km1; k++)
            {
              for (E_Int i = 0; i < im1; i++)
              {
                // jmin --------------------------------------------------
                j = 0;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + gamma;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0*cellN1*cellN2*cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 16.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }

                // jmax --------------------------------------------------
                j = jm1-1;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + gamma;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0*cellN1*cellN2*cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 16.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }
              }
            }
            //faces kmin & kmax
            #pragma omp for collapse(2)
            for (E_Int j = 0; j < jm1; j++)
            {
              for (E_Int i = 0; i < im1; i++)
              {
                // kmin --------------------------------------------------
                k = 0;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + beta;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0*cellN1*cellN2*cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 16.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }

                // kmax --------------------------------------------------
                k = km1-1;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + beta;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0*cellN1*cellN2*cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 16.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }
              }
            }
          }
        }
        break;
        
      case 3:
        printf("Warning: center2node: this case (cellN mod 3) is not implemented yet. Using mod 1.\n");
        if (dim == 1)
        {
          E_Int im1 = im+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, i, i0;
            E_Int ind0, ind1;
            E_Float cellN0, cellN1, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1; ind++)
            {
              i = ind;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
                        
              ind0 = i0;
              ind1 = ind0 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.); 

              w = cellN0*cellN1;

              if (K_FUNC::fEqualZero(w - 0.))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else if (dim == 2)
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, i, j, i0, j0;
            E_Int ind0, ind1, ind2, ind3;
            E_Float cellN0, cellN1, cellN2, cellN3, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1; ind++)
            {
              j = ind / im1;
              i = ind - j*im1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
                        
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              
              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
                  
              w = cellN0*cellN1*cellN2*cellN3;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
          E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
            E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
            E_Float cellN6, cellN7, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
            {
              k = ind / ijm1;
              j = (ind - k*ijm1) / im1;
              i = ind - j*im1 - k*ijm1;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              gamma = ijm;
              if (k == 0 || k == km) gamma = 0;
              k0 = E_max(k, 1)-1;
              
              ind0 = i0 + j0 * im + k0 * im*jm;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              ind4 = ind0 + gamma;
              ind5 = ind4 + alpha;
              ind6 = ind4 + beta;
              ind7 = ind6 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
              cellN4 = E_min(cellNp[ind4], 1.);
              cellN5 = E_min(cellNp[ind5], 1.);
              cellN6 = E_min(cellNp[ind6], 1.);
              cellN7 = E_min(cellNp[ind7], 1.);
                  
              w = cellN0*cellN1*cellN2*cellN3*cellN4*cellN5*cellN6*cellN7;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        break;
        
      default:
        printf("Warning: center2node: unknown cellnaturefield format.\n");
        return 0;
    }
  }
  if (cellN != -1 && algo == 1)
  {
    switch (mod)
    {
      case 1:
        if (dim == 1)
        {
          E_Int im1 = im+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, i, i0;
            E_Int ind0, ind1;
            E_Float cellN0, cellN1, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1; ind++)
            {
              i = ind;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
                        
              ind0 = i0;
              ind1 = ind0 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.); 

              w = cellN0 + cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else if (dim == 2)
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, i, j, i0, j0;
            E_Int ind0, ind1, ind2, ind3;
            E_Float cellN0, cellN1, cellN2, cellN3, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1; ind++)
            {
              j = ind / im1;
              i = ind - j*im1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
                        
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              
              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
                  
              w = cellN0+cellN1+cellN2+cellN3;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
          E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
            E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
            E_Float cellN6, cellN7, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
            {
              k = ind / ijm1;
              j = (ind - k*ijm1) / im1;
              i = ind - j*im1 - k*ijm1;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              gamma = ijm;
              if (k == 0 || k == km) gamma = 0;
              k0 = E_max(k, 1)-1;
              
              ind0 = i0 + j0 * im + k0 * im*jm;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              ind4 = ind0 + gamma;
              ind5 = ind4 + alpha;
              ind6 = ind4 + beta;
              ind7 = ind6 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
              cellN4 = E_min(cellNp[ind4], 1.);
              cellN5 = E_min(cellNp[ind5], 1.);
              cellN6 = E_min(cellNp[ind6], 1.);
              cellN7 = E_min(cellNp[ind7], 1.);
                  
              w = cellN0+cellN1+cellN2+cellN3+cellN4+cellN5+cellN6+cellN7;
              
              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        break;
        
      case 2:
        if (dim == 1)
        {
          E_Int im1 = im+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, i, i0;
            E_Int ind0, ind1;
            E_Float cellN0, cellN1, w;
            #pragma omp for
            for (E_Int ind = 1; ind < im1-1; ind++)
            {
              i = ind;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
                        
              ind0 = i0;
              ind1 = ind0 + alpha;

              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1]; 

              w = cellN0+cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(cellN0 - 2.) || K_FUNC::fEqualZero(cellN1 - 2.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            // i = 0
            cellNpn[0] = cellNp[0];
            // i = im
            cellNpn[im1-1] = cellNp[im-1];
          }
        }
        else if (dim == 2)
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          #pragma omp parallel
          {
            E_Int alpha, beta, i, j, i0, j0;
            E_Int ind, ind0, ind1, ind2, ind3;
            E_Float cellN0, cellN1, cellN2, cellN3, w;
            // main loop
            #pragma omp for
            for (E_Int ind = 1; ind < im1*jm1-1; ind++)
            {
              j = ind / im1;
              i = ind - j*im1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
                        
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
              cellN2 = cellNp[ind2];
              cellN3 = cellNp[ind3];
                  
              w = cellN0+cellN1+cellN2+cellN3;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(cellN0 - 2.) || K_FUNC::fEqualZero(cellN1 - 2.) ||
                    K_FUNC::fEqualZero(cellN2 - 2.) || K_FUNC::fEqualZero(cellN3 - 2.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            // edges imin, imax
            #pragma omp for
            for (E_Int j = 0; j < jm1; j++)
            {
              // i = 0
              i = 0;
              i0 = E_max(i, 1)-1;

              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + beta;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0+cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(w - 4.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }

              // i = im1-1
              i = im1-1;
              i0 = E_max(i, 1)-1;

              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + beta;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0+cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(w - 4.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            // edges jmin, jmax
            #pragma omp for
            for (E_Int i = 0; i < im1; i++)
            {
              // j = 0
              j = 0;
              j0 = E_max(j, 1)-1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0+cellN1;

              if (K_FUNC::fEqualZero(w - 0.))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(w - 4.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }

              // j = jm1-1
              j = jm1-1;
              j0 = E_max(j, 1)-1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              
              ind = i + j *im1;
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              
              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
                  
              w = cellN0+cellN1;

              if (K_FUNC::fEqualZero(w - 0.))
                cellNpn[ind] = 0.;
              else
              {
                if (K_FUNC::fEqualZero(w - 4.))
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
          }
        }        
        else
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
          E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
            E_Int ind, ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
            E_Float cellN6, cellN7, w, w2;
            // main loop
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
            {
              k = ind / ijm1;
              j = (ind - k*ijm1) / im1;
              i = ind - j*im1 - k*ijm1;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              gamma = ijm;
              if (k == 0 || k == km) gamma = 0;
              k0 = E_max(k, 1)-1;
              
              ind0 = i0 + j0 * im + k0 * im*jm;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              ind4 = ind0 + gamma;
              ind5 = ind4 + alpha;
              ind6 = ind4 + beta;
              ind7 = ind6 + alpha;

              cellN0 = cellNp[ind0];
              cellN1 = cellNp[ind1];
              cellN2 = cellNp[ind2];
              cellN3 = cellNp[ind3];
              cellN4 = cellNp[ind4];
              cellN5 = cellNp[ind5];
              cellN6 = cellNp[ind6];
              cellN7 = cellNp[ind7];
                  
              w = cellN0+cellN1+cellN2+cellN3+cellN4+cellN5+cellN6+cellN7;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
              {
                w2 = E_max(cellN0,1.)+E_max(cellN1,1.)+E_max(cellN2,1.)+E_max(cellN3,1.)+E_max(cellN4,1.)+E_max(cellN5,1.)+E_max(cellN6,1.)+E_max(cellN7,1.);
                if (w2 > 8.5)
                  cellNpn[ind] = 2.;
                else
                  cellNpn[ind] = 1.;
              }
            }
            //faces imin & imax
            #pragma omp for
            for (E_Int k = 0; k < km1; k++)
            {
              for (E_Int j = 0; j < jm1; j++)
              {
                // imin --------------------------------------------------
                i = 0;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + beta;
                ind2 = ind0 + gamma;
                ind3 = ind2 + beta;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0+cellN1+cellN2+cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 8.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }

                // imax --------------------------------------------------
                i = im1-1;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;

                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + beta;
                ind2 = ind0 + gamma;
                ind3 = ind2 + beta;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0+cellN1+cellN2+cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 8.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }
              }
            }
            //faces jmin & jmax
            #pragma omp for
            for (E_Int k = 0; k < km1; k++)
            {
              for (E_Int i = 0; i < im1; i++)
              {
                // jmin --------------------------------------------------
                j = 0;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + gamma;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0+cellN1+cellN2+cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 8.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }

                // jmax --------------------------------------------------
                j = jm1-1;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + gamma;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0+cellN1+cellN2+cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 8.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }
              }
            }
            //faces kmin & kmax
            #pragma omp for
            for (E_Int j = 0; j < jm1; j++)
            {
              for (E_Int i = 0; i < im1; i++)
              {
                // kmin --------------------------------------------------
                k = 0;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + beta;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0+cellN1+cellN2+cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 8.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }

                // kmax --------------------------------------------------
                k = km1-1;

                alpha = 1;
                if (i == 0 || i == im) alpha = 0;
                i0 = E_max(i, 1)-1;
                beta = im;
                if (j == 0 || j == jm) beta = 0;
                j0 = E_max(j, 1)-1;
                gamma = ijm;
                if (k == 0 || k == km) gamma = 0;
                k0 = E_max(k, 1)-1;
                
                ind = i + j*im1 + k*ijm1;

                ind0 = i0 + j0*im + k0*ijm;
                ind1 = ind0 + alpha;
                ind2 = ind0 + beta;
                ind3 = ind2 + alpha;

                cellN0 = cellNp[ind0];
                cellN1 = cellNp[ind1];
                cellN2 = cellNp[ind2];
                cellN3 = cellNp[ind3];
                    
                w = cellN0+cellN1+cellN2+cellN3;

                if (K_FUNC::fEqualZero(w))
                  cellNpn[ind] = 0.;
                else
                {
                  if (K_FUNC::fEqualZero(w - 8.))
                    cellNpn[ind] = 2.;
                  else
                    cellNpn[ind] = 1.;
                }
              }
            }
          }
        }
        break;
        
      case 3:
        printf("Warning: center2node: this case (cellN mod 3) is not implemented yet. Using mod 1.\n");
        if (dim == 1)
        {
          E_Int im1 = im+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, i, i0;
            E_Int ind0, ind1;
            E_Float cellN0, cellN1, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1; ind++)
            {
              i = ind;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
                        
              ind0 = i0;
              ind1 = ind0 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.); 

              w = cellN0 + cellN1;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else if (dim == 2)
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, i, j, i0, j0;
            E_Int ind0, ind1, ind2, ind3;
            E_Float cellN0, cellN1, cellN2, cellN3, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1; ind++)
            {
              j = ind / im1;
              i = ind - j*im1;

              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
                        
              ind0 = i0 + j0 * im;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              
              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
                  
              w = cellN0+cellN1+cellN2+cellN3;

              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        else
        {
          E_Int im1 = im+1; E_Int jm1 = jm+1; E_Int km1 = km+1;
          E_Int ijm = im*jm; E_Int ijm1= im1*jm1;
          E_Float* cellNp = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          
          #pragma omp parallel
          {
            E_Int alpha, beta, gamma, i, j, k, i0, j0, k0;
            E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float cellN0, cellN1, cellN2, cellN3, cellN4, cellN5;
            E_Float cellN6, cellN7, w;
            #pragma omp for
            for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
            {
              k = ind / ijm1;
              j = (ind - k*ijm1) / im1;
              i = ind - j*im1 - k*ijm1;
              alpha = 1;
              if (i == 0 || i == im) alpha = 0;
              i0 = E_max(i, 1)-1;
              beta = im;
              if (j == 0 || j == jm) beta = 0;
              j0 = E_max(j, 1)-1;
              gamma = ijm;
              if (k == 0 || k == km) gamma = 0;
              k0 = E_max(k, 1)-1;
              
              ind0 = i0 + j0 * im + k0 * im*jm;
              ind1 = ind0 + alpha;
              ind2 = ind0 + beta;
              ind3 = ind2 + alpha;
              ind4 = ind0 + gamma;
              ind5 = ind4 + alpha;
              ind6 = ind4 + beta;
              ind7 = ind6 + alpha;

              cellN0 = E_min(cellNp[ind0], 1.);
              cellN1 = E_min(cellNp[ind1], 1.);
              cellN2 = E_min(cellNp[ind2], 1.);
              cellN3 = E_min(cellNp[ind3], 1.);
              cellN4 = E_min(cellNp[ind4], 1.);
              cellN5 = E_min(cellNp[ind5], 1.);
              cellN6 = E_min(cellNp[ind6], 1.);
              cellN7 = E_min(cellNp[ind7], 1.);
                  
              w = cellN0+cellN1+cellN2+cellN3+cellN4+cellN5+cellN6+cellN7;
              
              if (K_FUNC::fEqualZero(w))
                cellNpn[ind] = 0.;
              else
                cellNpn[ind] = 1.;
            }
          }
        }
        break;
        
      default:
        printf("Warning: center2node: unknown cellnaturefield format.\n");
        return 0;
    } 
  }

  return 1;
}

//=============================================================================
// Convertit un array centres en array noeuds (BE et ME)
// Retourne 1 en cas de succes, 0 en cas d'echec.
// IN: mod: mod du cellN (0,1 ou 0,1,2)
// IN: algo: type de traitement pour le cellN
//=============================================================================
E_Int K_LOC::center2nodeUnstruct(FldArrayF& FCenter, 
                                 FldArrayI& cn,
                                 E_Int cellN, E_Int mod,
                                 E_Int posx, E_Int posy, E_Int posz,
                                 FldArrayF& FNode,
                                 E_Int algo)
{
  // Acces universel sur BE/ME
  E_Int nc = cn.getNConnect();
  E_Int nfld = FCenter.getNfld();
  E_Int nb = FNode.getSize();
  FNode.setAllValuesAtNull();
  FldArrayI count(nb); count.setAllValuesAtNull();
  E_Int* countp = count.begin();
  
  E_Int elOffset = 0; // element offset

  // Boucle sur toutes les connectivites une premiere fois
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int ne = cm.getSize(); // nombre de centres = nombre d'elements
    E_Int nt = cm.getNfld(); // nombre de points par elements de cette connectivite

    // Boucle sur tous les champs
    for (E_Int v = 1; v <= nfld; v++)
    {
      E_Float* fnode = FNode.begin(v);
      E_Float* fcen = FCenter.begin(v);
    
      for (E_Int n = 1; n <= nt; n++)
      {
        for (E_Int e = 0; e < ne; e++)
        {
          E_Int ind = cm(e, n) - 1;
          fnode[ind] += fcen[elOffset+e];
          // Increment the number of centers connected to vertex ind
          // Must loop over all connectivities to get the right countp
          // for vertices that share several types of basic elements 
          if (v == 1) countp[ind]++;
        }
      }
    }
    elOffset += ne; // increment offset
  }

  
  // Adim fields by countp
  #pragma omp parallel
  {
    E_Float inv;
    for (E_Int v = 1; v <= nfld; v++)
    {
      E_Float* fnode = FNode.begin(v);

      #pragma omp for
      for (E_Int n = 0; n < nb; n++)
      {
        if (countp[n] > 0)
        {
          inv = 1./countp[n];
          fnode[n] *= inv;
        }
      }
    }
  }

  elOffset = 0; // reset element offset for the second loop over connectivities

  // Boucle sur toutes les connectivites une second fois pour diviser les
  // champs aux noeuds par countp et traiter le cas special du champs cellN
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int ne = cm.getSize(); // nombre de centres = nombre d'elements
    E_Int nt = cm.getNfld(); // nombre de points par elements de cette connectivite

    // Traitement special pour le champ "cellnaturefield" - reecriture
    // de fnode(cellN), ie, cellNNode
    if (cellN != -1)
    {
      // champs "cellnaturefield" aux centres et aux noeuds
      E_Float* cellNCenter = FCenter.begin(cellN);
      E_Float* cellNNode = FNode.begin(cellN);
      // tableau temporaire pour stocker les sommes ou produits de cellN
      FldArrayF temp(nb); 
      switch (mod)
      {
        case 1: // cellN=0 (blanked or interpolated) or 1 (normal) 
        case 2: // cellN=0 (blanked), 1 (normal) or 2 (interpoled)
          if (algo == 0)
          {
            /* algo=0 - produits de cellN:
              Si au moins un 0 ==> 0
              Sinon 1
            */
            temp.setAllValuesAt(1.);
            for (E_Int n = 1; n <= nt; n++)
            {
              for (E_Int e = 0; e < ne; e++)
              {
                E_Int ind = cm(e, n) - 1;
                temp[ind] *= cellNCenter[elOffset+e];
              }
            }
            #pragma omp parallel for
            for (E_Int n = 0; n < nb; n++)
              cellNNode[n] = K_FUNC::E_min(temp[n], K_CONST::ONE);
          }
          else if (algo == 1)
          {
            /* algo=1 - somme des cellN:
              Si toutes les valeurs des cellN voisins = 0 ==> 0
              Sinon 1
            */
            temp.setAllValuesAt(0.);
            for (E_Int n = 1; n <= nt; n++)
            {
              for (E_Int e = 0; e < ne; e++)
              {
                E_Int ind = cm(e, n) - 1;
                temp[ind] += cellNCenter[elOffset+e];
              }
            }
            #pragma omp parallel for
            for (E_Int n = 0; n < nb; n++)
            {
              if (temp[n] > 0.) cellNNode[n] = 1.;
              else cellNNode[n] = 0.;
            }
          }
          break;
        case 3: // cellN=0 (blanked), cellN=1 (normal), cellN=-interpolationblock (interpoled)
          printf("Warning: center2node: this case is not implemented yet.\n");
          return 0;
          
        default:
          printf("Warning: center2node: unknown cellnaturefield format.\n");
          return 0;
      }
    }
    elOffset += ne; // increment element offset
  }

  return 1;
}
//=============================================================================
// Convertit un array NGON centres en array NGON noeuds
// FNode doit etre alloue au nb de noeuds
// IN: cNG: connectivite NGON
// IN: cEV: connectivite Elts/Vertex associes aux NGON
// IN: cellN: position du champ cellN dans FCenter
// IN: mod: type de cellN
// IN: algo: type de prise en compte du cellN
// Retourne 1 en cas de succes, 0 en cas d'echec.
//=============================================================================
E_Int K_LOC::center2nodeNGon(
  FldArrayF& FCenter, FldArrayI& cNG, vector< vector<E_Int> >& cEV, 
  FldArrayF& FNode, E_Int cellN, E_Int mod, E_Int algo)
{
  E_Int nfld = FCenter.getNfld();
  E_Int nb = FNode.getSize(); FNode.setAllValuesAtNull();
  FldArrayI count(nb); count.setAllValuesAtNull();
  E_Int* countp = count.begin();
  E_Int ind, nvert; E_Float inv;
  E_Int nelts = FCenter.getSize();

  for (E_Int i = 0; i < nelts; i++)
  {
    const vector<E_Int>& vertices = cEV[i];
    nvert = vertices.size();
    for (E_Int j = 0; j < nvert; j++)
    {
      ind = vertices[j]-1;
      for (E_Int eq = 1; eq <= nfld; eq++) FNode(ind,eq) += FCenter(i,eq);
      countp[ind]++;
    }
  }

  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fnode = FNode.begin(eq);
    for (ind = 0; ind < nb; ind++)
    {
      if (countp[ind] > 0)
      {
        inv = 1./E_max(countp[ind],1);
        fnode[ind] = fnode[ind] * inv;
      }
    }
  }

  // Traitement special pour le champ "cellnaturefield"
  /* algo=0:
     Si au moins un 0 ==> 0
     Sinon 1
  */
  if (cellN != -1 && algo == 0)
  {
    // champs "cellnaturefield" aux centres et aux noeuds
    E_Float* cellNCenter = FCenter.begin(cellN);
    E_Float* cellNNode = FNode.begin(cellN);
    // tableau temporaire pour stocker les sommes ou produits de cellN
    FldArrayF temp(nb); 
    switch (mod)
    {
      case 1: // cellN=0 (blanked or interpolated) or 1 (normal) 
      case 2: // cellN=0 (blanked), 1 (normal) or 2 (interpoled)
        temp.setAllValuesAt(1.);
        for (E_Int i = 0; i < nelts; i++)
        {
          const vector<E_Int>& vertices = cEV[i];
          nvert = vertices.size();
          for (E_Int j = 0; j < nvert; j++)
          {
            ind = vertices[j]-1;
            temp[ind] *= cellNCenter[i];
          }
        }
        for (E_Int indn = 0; indn < nb; indn++)
          cellNNode[indn] = K_FUNC::E_min(temp[indn], K_CONST::ONE);
        break;
      case 3: // cellN=0 (blanked), cellN=1 (normal), cellN=-interpolationblock (interpoled)
        printf("Warning: center2node: this case is not implemented yet.\n");
        return 0;
        
      default:
        printf("Warning: center2node: unknown cellnaturefield format.\n");
        return 0;
    }
  }
  /* algo=1:
     Si toutes les valeurs des cellN voisins = 0 ==> 0
     Sinon 1
  */
  if (cellN != -1 && algo == 1)
  {
    // champs "cellnaturefield" aux centres et aux noeuds
    E_Float* cellNCenter = FCenter.begin(cellN);
    E_Float* cellNNode = FNode.begin(cellN);
    // tableau temporaire pour stocker les sommes ou produits de cellN
    FldArrayF temp(nb); 
    switch (mod)
    {
      case 1: // cellN=0 (blanked or interpolated) or 1 (normal)
      case 2: // cellN=0 (blanked), 1 (normal) or 2 (interpoled)
        temp.setAllValuesAt(0.);
        for (E_Int i = 0; i < nelts; i++)
        {
          const vector<E_Int>& vertices = cEV[i];
          nvert = vertices.size();
          for (E_Int j = 0; j < nvert; j++)
          {
            ind = vertices[j]-1;
            temp[ind] += cellNCenter[i];
          }
        }
        for (E_Int indn = 0; indn < nb; indn++)
        {
          if (temp[indn] > 0.) cellNNode[indn] = 1.;
          else cellNNode[indn] = 0.;
        }
        break;
      case 3: // cellN=0 (blanked), cellN=1 (normal), cellN=-interpolationblock (interpoled)
        printf("Warning: center2node: this case is not implemented yet.\n");
        return 0;
        
      default:
        printf("Warning: center2node: unknown cellnaturefield format.\n");
        return 0;
    }
  }

  return 1;
}
