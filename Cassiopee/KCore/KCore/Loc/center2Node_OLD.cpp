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

extern "C"
{
  void k6conv2node1_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                     const E_Int& nfld, E_Float* fieldnode, 
                     E_Float* fieldcenter);
  void k6conv2node1p_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                      const E_Int& nfld, E_Float* fieldnode, 
                      E_Float* cellN, E_Float* fieldcenter);
  void k6conv2node12d_(const E_Int& nj, const E_Int& nk, 
                       const E_Int& nfld, E_Float* fieldnode, 
                       E_Float* fieldcenter);
  void k6conv2node12dp_(const E_Int& nj, const E_Int& nk, 
                        const E_Int& nfld, E_Float* fieldnode,
                        E_Float* cellN,
                        E_Float* fieldcenter);
  void k6conv2node11d_(const E_Int& ni, 
                       const E_Int& nfld, E_Float* fieldnode, 
                       E_Float* fieldcenter);
  void k6conv2node11dp_(const E_Int& ni, 
                        const E_Int& nfld, E_Float* fieldnode, E_Float* cellN, 
                        E_Float* fieldcenter);
  void k6convcoord2node_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                         E_Float* Fcx, E_Float* Fcy, E_Float* Fcz, 
                         E_Float* X, E_Float* Y, E_Float* Z, 
                         E_Float* Xp, E_Float* Yp, E_Float* Zp, 
                         E_Float* Fnx, E_Float* Fny, E_Float* Fnz);
  void k6convcoord2node2d_(const E_Int& ni, const E_Int& nj,
                           E_Float* Fcx, E_Float* Fcy, E_Float* Fcz,
                           E_Float* Fnx, E_Float* Fny, E_Float* Fnz);
  void k6convcoord2node1d_(const E_Int& ni, 
                           E_Float* Fcx, E_Float* Fcy, E_Float* Fcz,
                           E_Float* Fnx, E_Float* Fny, E_Float* Fnz);

  void k6conv2node21_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                      E_Float* fieldcenter, 
                      E_Float* fieldnode);
  void k6conv2node21p_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                       E_Float* fieldcenter, 
                       E_Float* fieldnode);
  void k6conv2node212d_(const E_Int& nj, const E_Int& nk, 
                        E_Float* fieldcenter, 
                        E_Float* fieldnode);
  void k6conv2node212dp_(const E_Int& nj, const E_Int& nk, 
                         E_Float* fieldcenter, 
                         E_Float* fieldnode);
  void k6conv2node211d_(const E_Int& ni, 
                        E_Float* fieldcenter, 
                        E_Float* fieldnode);
  void k6conv2node211dp_(const E_Int& ni, 
                         E_Float* fieldcenter, 
                         E_Float* fieldnode);
  void k6conv2node22_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                      E_Float* fieldcenter, 
                      E_Float* fieldnode);
  void k6conv2node22p_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                       E_Float* fieldcenter, 
                       E_Float* fieldnode);
  void k6conv2node222d_(const E_Int& nj, const E_Int& nk, 
                        E_Float* fieldcenter, 
                        E_Float* fieldnode);
  void k6conv2node222dp_(const E_Int& nj, const E_Int& nk, 
                         E_Float* fieldcenter, 
                         E_Float* fieldnode);
  void k6conv2node221d_(const E_Int& ni, 
                        E_Float* fieldcenter, 
                        E_Float* fieldnode);
  void k6conv2node221dp_(const E_Int& ni, 
                         E_Float* fieldcenter, 
                         E_Float* fieldnode);
}

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
E_Int K_LOC::center2nodeStruct_OLD(FldArrayF& FCenter, 
                               E_Int ni, E_Int nj, E_Int nk,
                               E_Int cellN, E_Int mod,
                               E_Int posx, E_Int posy, E_Int posz,
                               FldArrayF& FNode,
                               E_Int& nin, E_Int& njn, E_Int& nkn,
                               E_Int algo)
{
  E_Int nv = FCenter.getNfld();
  E_Int size = 0, dim = 3, im;
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
      k6conv2node11d_(im, nv, 
                      FCenter.begin(), FNode.begin());
      /*
#pragma omp parallel
      {
        E_Int alpha, ind0, ind1;
#pragma omp for 
        for (E_Int i = 0; i <= im; i++)
        {
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          ind0 = max(i, 1)-1;
          ind1 = ind0 + alpha;
          for (E_Int n = 1; n <= nv; n++)
            FNode(i,n) = 0.5*(FCenter(ind0,n)+FCenter(ind1,n));
        }
      }
      */
    }
    else if (dim == 2)
    {
      k6conv2node12d_(im, jm, nv, 
                      FCenter.begin(), FNode.begin());
      /*
      E_Int im1 = im+1; E_Int jm1 = jm+1;
#pragma omp parallel
      {
        E_Int alpha, beta, i, j, i0, j0, ind0, ind1, ind2, ind3;
#pragma omp for
        for (E_Int ind = 0; ind < im1*jm1; ind++)
        {
          j = ind / jm1;
          i = ind - j*jm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = max(i, 1)-1;    
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = max(j, 1)-1;

          ind0 = i0 + j0 * im;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;

          for (E_Int n = 1; n <= nv; n++)
              FNode(ind,n) = 0.25*(FCenter(ind0,n)+FCenter(ind1,n)+FCenter(ind2,n)+FCenter(ind3,n));
        }
      }
      */
    }
    else
    {
      k6conv2node1_(im, jm, km, nv, 
                    FCenter.begin(), FNode.begin());
      /*
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
          j = (ind - k*ijm1) / jm1;
          i = ind - j*jm1 - k*ijm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = max(j, 1)-1;
          gamma = ijm;
          if (k == 0 || k == km) gamma = 0;
          k0 = max(k, 1)-1;
                    
          ind0 = i0 + j0 * im;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          
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
      */
    }
  }
  else // algo=1 et cellN existe
  {
    if (dim == 1)
    {
      k6conv2node11dp_(im, nv, 
                       FCenter.begin(), FCenter.begin(cellN), FNode.begin());
      /*
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
          i0 = max(i, 1)-1;
                    
          ind0 = i0;
          ind1 = ind0 + alpha;
          
          cellN0 = min(cellNp[ind0], 1.);
          cellN1 = min(cellNp[ind1], 1.);
               
          w = cellN0 + cellN1;
          if (w == 0.)
          {
            w = 0.5; cellN0 = 1.; cellN1 = 1.;
          }
          else w = 1./w;
               
          for (E_Int n = 1; n <= nv; n++)
              FNode(ind,n) = w*(cellN0*FCenter(ind0,n)+cellN1*FCenter(ind1,n));
        }
      }
      */
    }
    else if (dim == 2)
    {  
      k6conv2node12dp_(im, jm, nv, 
                       FCenter.begin(), FCenter.begin(cellN), FNode.begin());
      
      /*
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
          j = ind / jm1;
          i = ind - j*jm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = max(j, 1)-1;
                    
          ind0 = i0 + j0 * im;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          
          cellN0 = min(cellNp[ind0], 1.);
          cellN1 = min(cellNp[ind1], 1.);
          cellN2 = min(cellNp[ind2], 1.);
          cellN3 = min(cellNp[ind3], 1.);
               
          w = cellN0 + cellN1 + cellN2 + cellN3;
          if (w == 0.)
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
    */
    }
    else
    {
      k6conv2node1p_(im, jm, km, nv, 
                     FCenter.begin(), FCenter.begin(cellN), FNode.begin());
      /*
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
          j = (ind - k*ijm1) / jm1;
          i = ind - j*jm1 - k*ijm1;
          alpha = 1;
          if (i == 0 || i == im) alpha = 0;
          i0 = max(i, 1)-1;
          beta = im;
          if (j == 0 || j == jm) beta = 0;
          j0 = max(j, 1)-1;
          gamma = ijm;
          if (k == 0 || k == km) gamma = 0;
          k0 = max(k, 1)-1;
          
          ind0 = i0 + j0 * im + k0 * im*jm;
          ind1 = ind0 + alpha;
          ind2 = ind0 + beta;
          ind3 = ind2 + alpha;
          ind4 = ind0 + gamma;
          ind5 = ind4 + alpha;
          ind6 = ind4 + beta;
          ind7 = ind6 + alpha;

          cellN0 = min(cellNp[ind0], 1.);
          cellN1 = min(cellNp[ind1], 1.);
          cellN2 = min(cellNp[ind2], 1.);
          cellN3 = min(cellNp[ind3], 1.);
          cellN4 = min(cellNp[ind4], 1.);
          cellN5 = min(cellNp[ind5], 1.);
          cellN6 = min(cellNp[ind6], 1.);
          cellN7 = min(cellNp[ind7], 1.);
               
          w = cellN0 + cellN1 + cellN2 + cellN3 + cellN4 + 
              cellN5 + cellN6 + cellN7;
          if (w == 0.)
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
      */
    }
  }

  // Traitement special pour les coords
  if (posx != -1 && posy != -1 && posz != -1)
  {
    posx++; posy++; posz++;
    if (dim == 1)
    {
      k6convcoord2node1d_(im+1, 
                          FCenter.begin(posx), 
                          FCenter.begin(posy), 
                          FCenter.begin(posz),
                          FNode.begin(posx),
                          FNode.begin(posy), 
                          FNode.begin(posz));  
    }
    else if (dim == 2)
    {
      k6convcoord2node2d_(im+1, jm+1,
                          FCenter.begin(posx),
                          FCenter.begin(posy), 
                          FCenter.begin(posz),
                          FNode.begin(posx),
                          FNode.begin(posy), 
                          FNode.begin(posz)); 
    }
    else
    {
      FldArrayF X(ni*nj*nk,3); FldArrayF Xp(ni+nj+nk,3);
      k6convcoord2node_(ni, nj, nk, FCenter.begin(posx), 
                        FCenter.begin(posy), FCenter.begin(posz),
                        X.begin(1), X.begin(2), X.begin(3),
                        Xp.begin(1), Xp.begin(2), Xp.begin(3),
                        FNode.begin(posx), FNode.begin(posy), 
                        FNode.begin(posz));  
    }
  }
  
  // Traitement special pour le "cellnaturefield"
  if (cellN != -1 && algo == 0)
  {
    switch (mod)
    {
      case 1:
        if (dim == 1)
          k6conv2node211d_(im, FCenter.begin(cellN), 
                           FNode.begin(cellN));
        else if (dim == 2)
          k6conv2node212d_(im, jm, FCenter.begin(cellN), 
                           FNode.begin(cellN));
        else
          k6conv2node21_(im, jm, km, FCenter.begin(cellN), 
                         FNode.begin(cellN));
        break;
        
      case 2:
        if (dim == 1)
          k6conv2node221d_(im, FCenter.begin(cellN), 
                           FNode.begin(cellN));
        else if (dim == 2)
        {
          k6conv2node222d_(im, jm, FCenter.begin(cellN), 
                           FNode.begin(cellN));          
        }
        else
          k6conv2node22_(im, jm, km, FCenter.begin(cellN), 
                         FNode.begin(cellN));
        break;
        
      case 3:
        printf("Warning: center2node: this case (cellN mod 3) is not implemented yet. Using mod 1.\n");
        if (dim == 1)
          k6conv2node211d_(im, FCenter.begin(cellN), 
                           FNode.begin(cellN));
        else if (dim == 2)
          k6conv2node212d_(im, jm, FCenter.begin(cellN), 
                           FNode.begin(cellN));
        else
          k6conv2node21_(im, jm, km, FCenter.begin(cellN), 
                         FNode.begin(cellN));
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
          k6conv2node211dp_(im, FCenter.begin(cellN),
                            FNode.begin(cellN));
        else if (dim == 2)
          k6conv2node212dp_(im, jm, FCenter.begin(cellN), 
                            FNode.begin(cellN));
        else
          k6conv2node21p_(im, jm, km, FCenter.begin(cellN), 
                          FNode.begin(cellN));
        break;
        
      case 2:
        if (dim == 1)
          k6conv2node221dp_(im, FCenter.begin(cellN), 
                            FNode.begin(cellN));
        else if (dim == 2)
        {
          k6conv2node222dp_(im, jm, FCenter.begin(cellN), 
                            FNode.begin(cellN));          
        }
        else
          k6conv2node22p_(im, jm, km, FCenter.begin(cellN), 
                          FNode.begin(cellN));
        break;
        
      case 3:
        printf("Warning: center2node: this case (cellN mod 3) is not implemented yet. Using mod 1.\n");
        if (dim == 1)
          k6conv2node211dp_(im, FCenter.begin(cellN),
                            FNode.begin(cellN));
        else if (dim == 2)
          k6conv2node212dp_(im, jm, FCenter.begin(cellN), 
                            FNode.begin(cellN));
        else
          k6conv2node21p_(im, jm, km, FCenter.begin(cellN), 
                          FNode.begin(cellN));
        break;
        
      default:
        printf("Warning: center2node: unknown cellnaturefield format.\n");
        return 0;
    } 
  }

  return 1;
}

//=============================================================================
// Convertit un array centres en array noeuds (elements basiques)
// Retourne 1 en cas de succes, 0 en cas d'echec.
// IN: mod: mod du cellN (0,1 ou 0,1,2)
// IN: algo: type de traitement pour le cellN
//=============================================================================
E_Int K_LOC::center2nodeUnstruct_OLD(FldArrayF& FCenter, 
                                 FldArrayI& c,
                                 E_Int cellN, E_Int mod,
                                 E_Int posx, E_Int posy, E_Int posz,
                                 FldArrayF& FNode,
                                 E_Int algo)
{
  // c'est la connectivite duale
  E_Int ne = c.getSize(); // nombre de centres = nombre d'elements
  E_Int nt = c.getNfld();
  E_Int nfld = FCenter.getNfld();
  E_Int nb = FNode.getSize();
  FNode.setAllValuesAtNull();
  FldArrayI count(nb); count.setAllValuesAtNull();
  E_Int* countp = count.begin();
  E_Int ind; E_Float inv;

  for (E_Int n = 1; n <= nt; n++)
  {
    E_Int* cn = c.begin(n);
    for (E_Int e = 0; e < ne; e++)
    {
      ind = cn[e]-1;
      for (E_Int v = 1; v <= nfld; v++)
        FNode(ind, v) += FCenter(e, v);
      
      countp[ind]++;
    }
  }
  for (E_Int v = 1; v <= nfld; v++)
  {
    E_Float* fnode = FNode.begin(v);
    for (E_Int e = 0; e < nb; e++)
    {
      if (countp[e] > 0)
      {
        inv = 1./countp[e];
        fnode[e] = fnode[e] * inv;
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
        for (E_Int n = 1; n <= nt; n++)
        {
          E_Int* cn = c.begin(n); 
          for (E_Int e = 0; e < ne; e++)
          {
            ind = cn[e]-1;
            temp[ind] =  temp[ind]*cellNCenter[e];
          }
        }
        for (E_Int indn = 0; indn < nb; indn++)
          cellNNode[indn] = K_FUNC::E_min(temp[indn],K_CONST::ONE);
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
        for (E_Int n = 1; n <= nt; n++)
        {
          E_Int* cn = c.begin(n); 
          for (E_Int e = 0; e < ne; e++)
          {
            ind = cn[e]-1;
            temp[ind] = temp[ind]+cellNCenter[e];
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
E_Int K_LOC::center2nodeNGon_OLD(
  FldArrayF& FCenter, FldArrayI& cNG, vector< vector<E_Int> >& cEV, 
  FldArrayF& FNode, E_Int cellN, E_Int mod, E_Int algo)
{
  E_Int nfld = FCenter.getNfld();
  E_Int nb = FNode.getSize(); FNode.setAllValuesAtNull();
  FldArrayI count(nb); count.setAllValuesAtNull();
  E_Int* countp = count.begin();
  E_Int ind, nvert; E_Float inv;
  E_Int nelts = FCenter.getSize();

  for (E_Int et = 0; et < nelts; et++)
  {
    vector<E_Int>& vertices = cEV[et];
    nvert = vertices.size();
    for (E_Int nov = 0; nov < nvert; nov++)
    {
      ind = vertices[nov]-1;
      for (E_Int eq = 1; eq <= nfld; eq++) FNode(ind,eq) += FCenter(et,eq);
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
        inv = 1./K_FUNC::E_max(countp[ind],1);
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
        for (E_Int et = 0; et < nelts; et++)
        {
          vector<E_Int>& vertices = cEV[et];
          nvert = vertices.size();
          for (E_Int nov = 0; nov < nvert; nov++)
          {
            ind = vertices[nov]-1;
            temp[ind] =  temp[ind]*cellNCenter[et];
          }
        }
        for (E_Int indn = 0; indn < nb; indn++)
          cellNNode[indn] = K_FUNC::E_min(temp[indn],K_CONST::ONE);
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
        for (E_Int et = 0; et < nelts; et++)
        {
          vector<E_Int>& vertices = cEV[et];
          nvert = vertices.size();
          for (E_Int nov = 0; nov < nvert; nov++)
          {
            ind = vertices[nov]-1;
            temp[ind] =  temp[ind]+cellNCenter[et];
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
