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

# include "loc.h"
# include "Connect/connect.h"
# include <vector>
# include <algorithm> 

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
// Convertit un array noeuds en array centres en structure
// Retourne 1 en cas de succes, 0 en cas d'echec.
//=============================================================================
E_Int K_LOC::node2centerStruct(FldArrayF& FNode, 
                               E_Int ni, E_Int nj, E_Int nk,
                               E_Int cellN, E_Int mod, 
                               FldArrayF& FCenter)
{
  E_Int nv = FNode.getNfld();
  E_Int size, dim;  
  E_Int im, jm, km;

  if (ni != 1 && nj == 1 && nk == 1 )
  {
    size = ni-1; dim = 1; im = ni;  jm = 1; km = 1;
  }
  else if (ni == 1 && nj != 1 && nk == 1) 
  {
    size = nj-1; dim = 1; im = nj; jm = 1; km = 1;
  }
  else if (ni == 1 && nj == 1 && nk != 1)
  {
    size = nk-1; dim = 1; im = nk; jm = 1; km = 1;
  }
  else
  {
    if (ni == 1)
    {
      size = (nj-1)*(nk-1); dim = 2; im = nj; jm = nk; km = 1;
    }
    else if (nj == 1)
    {
      size = (ni-1)*(nk-1); dim = 2; im = ni; jm = nk; km = 1;
    }
    else if (nk == 1)
    {
      size = (ni-1)*(nj-1); dim = 2; im = ni; jm = nj; km = 1;
    }
    else
    {
      size = (ni-1)*(nj-1)*(nk-1); dim = 3; im = ni; jm = nj; km = nk;
    }
  }
  
  // On alloue FCenter seulement s'il n'est pas deja alloue correctement
  if (FCenter.getSize() != size || FCenter.getNfld() != nv)
    FCenter.malloc(size, nv);

  // In of each "block" i, converts field defined in nodes to field 
  // defined in centers
  if (dim == 1)
  {
    E_Int imc = im-1;
    #pragma omp parallel
    {
      E_Int ind0, ind1;
      #pragma omp for 
      for (E_Int i = 0; i < imc; i++)
      {
        ind0 = i;
        ind1 = i+1;
        for (E_Int n = 1; n <= nv; n++)
          FCenter(i,n) = 0.5*(FNode(ind0,n)+FNode(ind1,n));
      }
    }
  }
  else if (dim == 2)
  {
    E_Int imc = im-1; E_Int jmc = jm-1;
    #pragma omp parallel
    {
      E_Int i, j, ind0, ind1, ind2, ind3;
      #pragma omp for
      for (E_Int ind = 0; ind < imc*jmc; ind++)
      {
        j = ind / imc;
        i = ind - j*imc;

        ind0 = i + j*im;
        ind1 = ind0+1;
        ind2 = ind0+im;
        ind3 = ind2+1;

        for (E_Int n = 1; n <= nv; n++)
          FCenter(ind,n) = 0.25*(FNode(ind0,n)+FNode(ind1,n)+FNode(ind2,n)+FNode(ind3,n));
      }
    }
  }
  else
  {
    E_Int imc = im-1; E_Int jmc = jm-1; E_Int kmc = km-1;
    E_Int ijm = im*jm; E_Int ijmc= imc*jmc;
    #pragma omp parallel
    {
      E_Int i, j, k;
      E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
      #pragma omp for
      for (E_Int ind = 0; ind < imc*jmc*kmc; ind++)
      {
        k = ind / ijmc;
        j = (ind - k*ijmc) / imc;
        i = ind - j*imc - k*ijmc;
        
        ind0 = i + j*im + k*ijm;
        ind1 = ind0 + 1;
        ind2 = ind0 + im;
        ind3 = ind2 + 1;
        ind4 = ind0 + ijm;
        ind5 = ind4 + 1;
        ind6 = ind4 + im;
        ind7 = ind6 + 1;

        for (E_Int n = 1; n <= nv; n++)
          FCenter(ind,n) = 0.125*(FNode(ind0,n)+FNode(ind1,n)+FNode(ind2,n)+FNode(ind3,n)+
                                  FNode(ind4,n)+FNode(ind5,n)+FNode(ind6,n)+FNode(ind7,n));
      }
    }
  }
  
  // If field contains "cellnaturefield"
  if (cellN != -1)
  {
    switch (mod)
    {
      case 1:
        if (dim == 1)
        {
          E_Float* cellNpc = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          E_Int imc = im-1;
          #pragma omp parallel
          {
            E_Int ind0, ind1;
            #pragma omp for 
            for (E_Int i = 0; i < imc; i++)
            {
              ind0 = i;
              ind1 = i+1;
              cellNpc[i] = E_min(cellNpn[ind0],cellNpn[ind1]);
            }
          }
        }
        else if (dim == 2)
        {
          E_Float* cellNpc = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          E_Int imc = im-1; E_Int jmc = jm-1;
          #pragma omp parallel
          {
            E_Int i, j, ind0, ind1, ind2, ind3;
            #pragma omp for
            for (E_Int ind = 0; ind < imc*jmc; ind++)
            {
              j = ind / imc;
              i = ind - j*imc;

              ind0 = i + j*im;
              ind1 = ind0+1;
              ind2 = ind0+im;
              ind3 = ind2+1;

              cellNpc[ind] = E_min(E_min(E_min(cellNpn[ind0],cellNpn[ind1]),cellNpn[ind2]),cellNpn[ind3]);
            }
          }
        }
        else
        {
          E_Float* cellNpc = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          E_Int imc = im-1; E_Int jmc = jm-1; E_Int kmc = km-1;
          E_Int ijm = im*jm; E_Int ijmc = imc*jmc;
          #pragma omp parallel
          {
            E_Int i, j, k;
            E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            #pragma omp for
            for (E_Int ind = 0; ind < imc*jmc*kmc; ind++)
            {
              k = ind / ijmc;
              j = (ind - k*ijmc) / imc;
              i = ind - j*imc - k*ijmc;
              
              ind0 = i + j*im + k*ijm;
              ind1 = ind0 + 1;
              ind2 = ind0 + im;
              ind3 = ind2 + 1;
              ind4 = ind0 + ijm;
              ind5 = ind4 + 1;
              ind6 = ind4 + im;
              ind7 = ind6 + 1;

              cellNpc[ind] = E_min(E_min(E_min(E_min(E_min(E_min(E_min(cellNpn[ind0],cellNpn[ind1]),cellNpn[ind2]),cellNpn[ind3]),cellNpn[ind4]),cellNpn[ind5]),cellNpn[ind5]),cellNpn[ind7]);
            }
          }
        }
        break;
        
      case 2:
        if (dim == 1)
        {
          E_Float* cellNpc = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          E_Int imc = im-1;
          #pragma omp parallel
          {
            E_Int ind0, ind1;
            E_Float somme;
            #pragma omp for 
            for (E_Int i = 1; i < imc-1; i++)
            {
              ind0 = i;
              ind1 = i+1;
              somme = cellNpn[ind0]+cellNpn[ind1];
              if (K_FUNC::fEqualZero(somme))
                cellNpc[i] = 0.;
              else if (K_FUNC::fEqualZero(cellNpn[ind0] - 2.) || K_FUNC::fEqualZero(cellNpn[ind1] - 2.))
                cellNpc[i] = 2.;
              else
                cellNpc[i] = 1.;
            }
            // i = 0
            somme = cellNpn[0]+cellNpn[1];
            if (K_FUNC::fEqualZero(somme))
              cellNpc[0] = 0.;
            else if (K_FUNC::fEqualZero(somme - 4.))
              cellNpc[0] = 2.;
            else
              cellNpc[0] = 1.;
            // i = imc-1
            somme = cellNpn[imc-1]+cellNpn[imc];
            if (K_FUNC::fEqualZero(somme))
              cellNpc[imc-1] = 0.;
            else if (K_FUNC::fEqualZero(somme - 4.))
              cellNpc[imc-1] = 2.;
            else
              cellNpc[imc-1] = 1.;
          }
        }
        else if (dim == 2)
        {
          E_Float* cellNpc = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          E_Int imc = im-1; E_Int jmc = jm-1;
          #pragma omp parallel
          {
            E_Int i, j, ind, ind0, ind1, ind2, ind3;
            E_Float somme;
            // main loop
            #pragma omp for
            for (E_Int ind = 0; ind < imc*jmc; ind++)
            {
              j = ind / imc;
              i = ind - j*imc;

              ind0 = i + j*im;
              ind1 = ind0+1;
              ind2 = ind0+im;
              ind3 = ind2+1;

              somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3];

              if (K_FUNC::fEqualZero(somme))
                cellNpc[ind] = 0.;
              else if (K_FUNC::fEqualZero(cellNpn[ind0] - 2.) || K_FUNC::fEqualZero(cellNpn[ind1] - 2.) || 
                       K_FUNC::fEqualZero(cellNpn[ind2] - 2.) || K_FUNC::fEqualZero(cellNpn[ind3] - 2.))
                cellNpc[ind] = 2.;
              else
                cellNpc[ind] = 1.;
            }
            // edges imin, imax
            #pragma omp for
            for (E_Int j = 0; j < jmc; j++)
            {
              // i = 0
              i = 0;
              ind = i + j*imc;
              ind0 = i + j*im;
              ind1 = ind0+1;
              ind2 = ind0+im;
              ind3 = ind2+1;

              somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3];

              if (K_FUNC::fEqualZero(somme))
                cellNpc[ind] = 0.;
              else if (K_FUNC::fEqualZero(somme - 8.))
                cellNpc[ind] = 2.;
              else
                cellNpc[ind] = 1.;

              // i = imc-1
              i = imc-1;
              ind = i + j*imc;
              ind0 = i + j*im;
              ind1 = ind0+1;
              ind2 = ind0+im;
              ind3 = ind2+1;

              somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3];

              if (K_FUNC::fEqualZero(somme))
                cellNpc[ind] = 0.;
              else if (K_FUNC::fEqualZero(somme - 8.))
                cellNpc[ind] = 2.;
              else
                cellNpc[ind] = 1.;
            }
            // edges jmin, jmax
            #pragma omp for
            for (E_Int i = 0; i < imc; i++)
            {
              // j = 0
              j = 0;
              ind = i + j*imc;
              ind0 = i + j*im;
              ind1 = ind0+1;
              ind2 = ind0+im;
              ind3 = ind2+1;

              somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3];

              if (K_FUNC::fEqualZero(somme))
                cellNpc[ind] = 0.;
              else if (K_FUNC::fEqualZero(somme - 8.))
                cellNpc[ind] = 2.;
              else
                cellNpc[ind] = 1.;

              // j = jmc-1
              j = jmc-1;
              ind = i + j*imc;
              ind0 = i + j*im;
              ind1 = ind0+1;
              ind2 = ind0+im;
              ind3 = ind2+1;

              somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3];

              if (K_FUNC::fEqualZero(somme))
                cellNpc[ind] = 0.;
              else if (K_FUNC::fEqualZero(somme - 8.))
                cellNpc[ind] = 2.;
              else
                cellNpc[ind] = 1.;
            }
          }
        }
        else
        {
          E_Float* cellNpc = FCenter.begin(cellN);
          E_Float* cellNpn = FNode.begin(cellN);
          E_Int imc = im-1; E_Int jmc = jm-1; E_Int kmc = km-1;
          E_Int ijm = im*jm; E_Int ijmc= imc*jmc;
          #pragma omp parallel
          {
            E_Int i, j, k;
            E_Int ind, ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
            E_Float somme, somme2;
            // main loop
            #pragma omp for
            for (E_Int ind = 0; ind < imc*jmc*kmc; ind++)
            {
              k = ind / ijmc;
              j = (ind - k*ijmc) / imc;
              i = ind - j*imc - k*ijmc;
              
              ind0 = i + j*im + k*ijm;
              ind1 = ind0 + 1;
              ind2 = ind0 + im;
              ind3 = ind2 + 1;
              ind4 = ind0 + ijm;
              ind5 = ind4 + 1;
              ind6 = ind4 + im;
              ind7 = ind6 + 1;

              somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

              if (K_FUNC::fEqualZero(somme))
                cellNpc[ind] = 0.;
              else
              {
                somme2 = E_max(cellNpn[ind0],1.)+E_max(cellNpn[ind1],1.)+E_max(cellNpn[ind2],1.)+E_max(cellNpn[ind3],1.)+E_max(cellNpn[ind4],1.)+E_max(cellNpn[ind5],1.)+E_max(cellNpn[ind6],1.)+E_max(cellNpn[ind7],1.);
                if (somme2 > 8.5)
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;
              }
            }
            //faces imin & imax
            #pragma omp for
            for (E_Int k = 0; k < kmc; k++)
            {
              for (E_Int j = 0; j < jmc; j++)
              {
                // imin --------------------------------------------------
                i = 0;

                ind = i + j*imc + k*ijmc;
                ind0 = i + j*im + k*ijm;
                ind1 = ind0 + 1;
                ind2 = ind0 + im;
                ind3 = ind2 + 1;
                ind4 = ind0 + ijm;
                ind5 = ind4 + 1;
                ind6 = ind4 + im;
                ind7 = ind6 + 1;

                somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

                if (K_FUNC::fEqualZero(somme))
                  cellNpc[ind] = 0.;
                else if (K_FUNC::fEqualZero(somme - 16.))
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;

                // imax --------------------------------------------------
                i = imc-1;

                ind = i + j*imc + k*ijmc;
                ind0 = i + j*im + k*ijm;
                ind1 = ind0 + 1;
                ind2 = ind0 + im;
                ind3 = ind2 + 1;
                ind4 = ind0 + ijm;
                ind5 = ind4 + 1;
                ind6 = ind4 + im;
                ind7 = ind6 + 1;

                somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

                if (K_FUNC::fEqualZero(somme))
                  cellNpc[ind] = 0.;
                else if (K_FUNC::fEqualZero(somme - 16.))
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;
              }
            }
            //faces jmin & jmax
            #pragma omp for
            for (E_Int k = 0; k < kmc; k++)
            {
              for (E_Int i = 0; i < imc; i++)
              {
                // jmin --------------------------------------------------
                j = 0;

                ind = i + j*imc + k*ijmc;
                ind0 = i + j*im + k*ijm;
                ind1 = ind0 + 1;
                ind2 = ind0 + im;
                ind3 = ind2 + 1;
                ind4 = ind0 + ijm;
                ind5 = ind4 + 1;
                ind6 = ind4 + im;
                ind7 = ind6 + 1;

                somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

                if (K_FUNC::fEqualZero(somme))
                  cellNpc[ind] = 0.;
                else if (K_FUNC::fEqualZero(somme - 16.))
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;

                // jmax --------------------------------------------------
                j = jmc-1;

                ind = i + j*imc + k*ijmc;
                ind0 = i + j*im + k*ijm;
                ind1 = ind0 + 1;
                ind2 = ind0 + im;
                ind3 = ind2 + 1;
                ind4 = ind0 + ijm;
                ind5 = ind4 + 1;
                ind6 = ind4 + im;
                ind7 = ind6 + 1;

                somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

                if (K_FUNC::fEqualZero(somme))
                  cellNpc[ind] = 0.;
                else if (K_FUNC::fEqualZero(somme - 16.))
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;
              }
            }
            //faces kmin & kmax
            #pragma omp for
            for (E_Int j = 0; j < jmc; j++)
            {
              for (E_Int i = 0; i < imc; i++)
              {
                // kmin --------------------------------------------------
                k = 0;

                ind = i + j*imc + k*ijmc;
                ind0 = i + j*im + k*ijm;
                ind1 = ind0 + 1;
                ind2 = ind0 + im;
                ind3 = ind2 + 1;
                ind4 = ind0 + ijm;
                ind5 = ind4 + 1;
                ind6 = ind4 + im;
                ind7 = ind6 + 1;

                somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

                if (K_FUNC::fEqualZero(somme))
                  cellNpc[ind] = 0.;
                else if (K_FUNC::fEqualZero(somme - 16.))
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;

                // kmax --------------------------------------------------
                k = kmc-1;

                ind = i + j*imc + k*ijmc;
                ind0 = i + j*im + k*ijm;
                ind1 = ind0 + 1;
                ind2 = ind0 + im;
                ind3 = ind2 + 1;
                ind4 = ind0 + ijm;
                ind5 = ind4 + 1;
                ind6 = ind4 + im;
                ind7 = ind6 + 1;

                somme = cellNpn[ind0]+cellNpn[ind1]+cellNpn[ind2]+cellNpn[ind3]+cellNpn[ind4]+cellNpn[ind5]+cellNpn[ind6]+cellNpn[ind7];

                if (K_FUNC::fEqualZero(somme))
                  cellNpc[ind] = 0.;
                else if (K_FUNC::fEqualZero(somme - 16.))
                  cellNpc[ind] = 2.;
                else
                  cellNpc[ind] = 1.;
              }
            }
          }
        }
        break;
        
      case 3:
        printf("Warning: node2center: this cellN type is not implemented yet.\n");
        return 0;
        break;
        
      default:
        printf("Warning: node2center: unknown cellnaturefield format.\n");
        return 0;
    }
  }
  return 1;

}
//=============================================================================
// Convertit un array noeuds en array centres en non-structure
// Le traitement specifique cellN n'est pas implemente.
// Retourne 1 en cas de succes, 0 en cas d'echec.
// FCenter doit deja etre alloue au nb d'elements
//=============================================================================
E_Int K_LOC::node2centerUnstruct(FldArrayF& FNode, 
                                 FldArrayI& c,
                                 E_Int cellN, E_Int mod, 
                                 FldArrayF& FCenter)
{
  // Acces universel sur BE/ME
  E_Int nc = c.getNConnect();
  E_Int nfld = FNode.getNfld();
  FCenter.setAllValuesAtNull();
  E_Int elOffset = 0; // element offset

  // Boucle sur toutes les connectivites
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(c.getConnect(ic));
    E_Int ne = cm.getSize(); // nombre de centres = nombre d'elements
    E_Int nt = cm.getNfld(); // nombre de points par elements de cette connectivite
    E_Float ntinv = 1./nt;
    
    // Boucle sur tous les champs
    #pragma omp parallel
    {      
      E_Int ind;
      for (E_Int v = 1; v <= nfld; v++)
      {
        E_Float* fnode = FNode.begin(v);
        E_Float* fcen = FCenter.begin(v);

        for (E_Int n = 1; n <= nt; n++)
        {
          #pragma omp for
          for (E_Int e = 0; e < ne; e++)
          {
            ind = cm(e, n) - 1;
            fcen[elOffset+e] += fnode[ind];
          }
        }
        #pragma omp for
        for (E_Int e = 0; e < ne; e++) fcen[elOffset+e] *= ntinv;
      }
    }

    elOffset += ne; // increment offset
  }
  return 1;
}
//===============================================================================
// Convertit un champ en noeuds en champ en centres en NGON
// Retourne 1 en cas de succes, 0 en cas d'echec.
// FCenter doit deja etre alloue au nb d'elements
// sorted=1: vertices coordinates are sorted for better accuracy in summations 
//===============================================================================
E_Int K_LOC::node2centerNGon(FldArrayF& FNode, FldArrayI& cNG,
                             FldArrayF& FCenter, E_Int sorted)
{
  E_Int ncells = cNG.getNElts();
  E_Int nfld = FNode.getNfld();
  FCenter.setAllValuesAtNull();

  if (sorted == 0)
  {
    std::vector< std::vector<E_Int> > cEV(ncells);
    K_CONNECT::connectNG2EV(cNG, cEV);

    for (E_Int v = 1; v <= nfld; v++)
    {
      E_Float* fnode = FNode.begin(v);
      E_Float* fcen = FCenter.begin(v);

      #pragma omp parallel
      {
        E_Int nvert, indv;
        E_Float ntinv;
        
        #pragma omp for
        for (E_Int et = 0; et < ncells; et++)
        {
          std::vector<E_Int>& vertices = cEV[et]; // noeuds associes a l'element et
          nvert = vertices.size();
          ntinv = 1./nvert;
          for (E_Int nv = 0; nv < nvert; nv++)
          {
            indv = vertices[nv]-1;
            fcen[et] += fnode[indv];
          }
          fcen[et] *= ntinv;
        }// loop on elts
      }
    }
  }
  else
  {
    // Acces non universel sur le ptrs
    E_Int* ngon = cNG.getNGon();
    E_Int* nface = cNG.getNFace();
    E_Int* indPG = cNG.getIndPG();
    E_Int* indPH = cNG.getIndPH();

    #pragma omp parallel
    {
      std::vector<E_Int> vertices;
      E_Int indv, nfaces, nv;
      vertices.reserve(1024);
      
      #pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        // Acces universel element i
        E_Int* elt = cNG.getElt(i, nfaces, nface, indPH);
        for (E_Int n=0; n < nfaces; n++)
        {
          // Acces universel face elt[n]-1
          E_Int* face = cNG.getFace(elt[n]-1, nv, ngon, indPG);
          for (E_Int p = 0; p < nv; p++)
          {
            indv = face[p];
            vertices.push_back(indv);
          }//loop on vertices
        }
        std::sort(vertices.begin(), vertices.end());
        vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end() );

        std::vector<E_Float> fsort(vertices.size());
        for (E_Int nofld = 1; nofld <= nfld; nofld++)
        {
          E_Float* fnode = FNode.begin(nofld);
          for (size_t nov = 0; nov < vertices.size(); nov++)
          {
            E_Int indv = vertices[nov]; fsort[nov] = fnode[indv];
          }
          std::sort(fsort.begin(), fsort.end());
          E_Float* fcen = FCenter.begin(nofld);
          for (size_t nov = 0; nov < fsort.size(); nov++) {fcen[i] += fsort[nov];}
          E_Float inv = 1./E_Float(fsort.size()); fcen[i] *= inv;   
        }
      }
    }
  }
  return 1;
}
