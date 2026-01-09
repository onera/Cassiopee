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

# include "transform.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
void K_TRANSFORM::breakMixedElements(
  FldArrayF& field, FldArrayI& ce,
  vector<FldArrayI*>& cEV, vector<FldArrayF*>& fields, vector<E_Int>& eltType)
{
  E_Int npts = field.getSize(); E_Int nfld = field.getNfld();
  E_Int* cnp = ce.begin();

  E_Int netbar = 0; E_Int nptsbar = 0;
  E_Int nettri = 0; E_Int nptstri = 0;
  E_Int netquad = 0; E_Int nptsquad = 0;
  E_Int nettetra = 0; E_Int nptstetra = 0;
  E_Int nethexa = 0; E_Int nptshexa = 0;
  E_Int netpenta = 0; E_Int nptspenta = 0;
  E_Int netpyra = 0; E_Int nptspyra = 0;

  // Compte les types d'elements
  E_Int ps = 0; E_Int ntype;
  E_Int size = ce.getSize();

  while (ps < size)
  {
    ntype = cnp[0];
    if (ntype == 3) // BAR
    {
      netbar += 1; nptsbar += 2;
      ps += 3; cnp += 3;
    }
    else if (ntype == 5) // TRI
    {
      nettri += 1; nptstri += 3;
      ps += 4; cnp += 4;
    }
    else if (ntype == 7) // QUAD
    {
      netquad += 1; nptsquad += 4;
      ps += 5; cnp += 5;
    }
    else if (ntype == 10) // TETRA
    {
      nettetra += 1; nptstetra += 4;
      ps += 5; cnp += 5;
    }
    else if (ntype == 12) // PYRA
    {
      netpyra += 1; nptspyra += 5;
      ps += 6; cnp += 6;
    }
    else if (ntype == 14) // PENTA
    {
      netpenta += 1; nptspenta += 6;
      ps += 7; cnp += 7;
    }
    else if (ntype == 17) // HEXA
    {
      nethexa += 1; nptshexa += 8;
      ps += 9; cnp += 9;
    }
    else
    {
      printf("Warning: breakElements: unknow type of element.\n");
    }
  }

  // Remplit
  FldArrayI* cEVbarp = new FldArrayI(netbar, 2);
  FldArrayF* fbarp = new FldArrayF(nptsbar, nfld);
  FldArrayF& fbar = *fbarp; FldArrayI& cEVbar = *cEVbarp;
  FldArrayI indirbF(npts); indirbF.setAllValuesAt(-1);
  E_Int* indirb = indirbF.begin();

  FldArrayI* cEVtrip = new FldArrayI(nettri,3);
  FldArrayF* ftrip = new FldArrayF(nptstri,nfld);
  FldArrayF& ftri = *ftrip; FldArrayI& cEVtri = *cEVtrip;
  FldArrayI indirtF(npts); indirtF.setAllValuesAt(-1);
  E_Int* indirt = indirtF.begin();

  FldArrayI* cEVquadp = new FldArrayI(netquad,4);
  FldArrayF* fquadp = new FldArrayF(nptsquad,nfld);
  FldArrayF& fquad = *fquadp; FldArrayI& cEVquad = *cEVquadp;
  FldArrayI indirqF(npts); indirqF.setAllValuesAt(-1);
  E_Int* indirq = indirqF.begin();

  FldArrayI* cEVtetrap = new FldArrayI(nettetra,4);
  FldArrayF* ftetrap = new FldArrayF(nptstetra,nfld);
  FldArrayF& ftetra = *ftetrap; FldArrayI& cEVtetra = *cEVtetrap;
  FldArrayI indirttF(npts); indirttF.setAllValuesAt(-1);
  E_Int* indirtt = indirttF.begin();

  FldArrayI* cEVpentap = new FldArrayI(netpenta, 6);
  FldArrayF* fpentap = new FldArrayF(nptspenta, nfld);
  FldArrayF& fpenta = *fpentap; FldArrayI& cEVpenta = *cEVpentap;
  FldArrayI indirpF(npts); indirpF.setAllValuesAt(-1);
  E_Int* indirp = indirpF.begin();

  FldArrayI* cEVpyrap = new FldArrayI(netpyra, 5);
  FldArrayF* fpyrap = new FldArrayF(nptspyra,nfld);
  FldArrayF& fpyra = *fpyrap; FldArrayI& cEVpyra = *cEVpyrap;
  FldArrayI indiryF(npts); indiryF.setAllValuesAt(-1);
  E_Int* indiry = indiryF.begin();

  FldArrayI* cEVhexap = new FldArrayI(nethexa, 8);
  FldArrayF* fhexap = new FldArrayF(nptshexa,nfld);
  FldArrayF& fhexa = *fhexap; FldArrayI& cEVhexa = *cEVhexap;
  FldArrayI indirhF(npts); indirhF.setAllValuesAt(-1);
  E_Int* indirh = indirhF.begin();

  E_Int et;
  cnp = ce.begin();
  ps = 0;
  netbar = 0; nptsbar = 0;
  nettri = 0; nptstri = 0;
  netquad = 0; nptsquad = 0;
  nettetra = 0; nptstetra = 0;
  nethexa = 0; nptshexa = 0;
  netpenta = 0; nptspenta = 0;
  netpyra = 0; nptspyra = 0;
  printf("breaking MIXED\n");

  while (ps < size)
  {
    ntype = cnp[0];

    if (ntype == 3) //BAR
    {
      for (E_Int nov = 0; nov < 2; nov++)
      {
        et = cnp[nov+1]-1;
        if (indirb[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) fbar(nptsbar,eq) = field(et,eq);
          nptsbar++;
          indirb[et] = nptsbar;
        }
      }
      for (E_Int nov = 0; nov < 2; nov++)
      {
        cEVbar(netbar, nov+1) = indirb[cnp[nov+1]-1];
      }
      netbar++;
      cnp += 2+1; ps += 2+1;
    }
    else if (ntype == 5) // TRI
    {
      for (E_Int nov = 0; nov < 3; nov++)
      {
        et = cnp[nov+1]-1;
        if (indirt[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) ftri(nptstri,eq) = field(et,eq);
          nptstri++;
          indirt[et] = nptstri;
        }
      }
      for (E_Int nov = 0; nov < 3; nov++)
      {
        cEVtri(nettri,nov+1) = indirt[cnp[nov+1]-1];
      }
      nettri++;
      cnp += 3+1; ps += 3+1;
    }
    else if (ntype == 7) // QUAD
    {
      for (E_Int nov = 0; nov < 4; nov++)
      {
        et = cnp[nov+1]-1;
        if (indirq[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) fquad(nptsquad,eq) = field(et,eq);
          nptsquad++;
          indirq[et] = nptsquad;
        }
      }
      for (E_Int nov = 0; nov < 4; nov++)
      {
        cEVquad(netquad,nov+1) = indirq[cnp[nov+1]-1];
      }
      netquad++;
      cnp += 4+1; ps += 4+1;
    }
    else if (ntype == 10) // TETRA
    {
      for (E_Int nov = 0; nov < 4; nov++)
      {
        et = cnp[nov+1]-1;
        if (indirtt[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) ftetra(nptstetra,eq) = field(et,eq);
          nptstetra++;
          indirtt[et] = nptstetra;
        }
      }
      for (E_Int nov = 0; nov < 4; nov++)
      {
        cEVtetra(nettetra,nov+1) = indirtt[cnp[nov+1]-1];
      }
      nettetra++;
      cnp += 4+1; ps += 4+1;
    }
    else if (ntype == 12) // PYRA
    {
      for (E_Int nov = 0; nov < 5; nov++)
      {
        et = cnp[nov+1]-1;
        if (indiry[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) fpyra(nptspyra,eq) = field(et,eq);
          nptspyra++;
          indiry[et] = nptspyra;
        }
      }
      for (E_Int nov = 0; nov < 5; nov++)
      {
        cEVpyra(netpyra,nov+1) = indiry[cnp[nov+1]-1];
      }
      netpyra++;
      cnp += 5+1; ps += 5+1;
    }
    else if (ntype == 14) // PENTA
    {
      for (E_Int nov = 0; nov < 6; nov++)
      {
        et = cnp[nov+1]-1;
        if (indirp[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) fpenta(nptspenta,eq) = field(et,eq);
          nptspenta++;
          indirp[et] = nptspenta;
        }
      }
      for (E_Int nov = 0; nov < 6; nov++)
      {
        cEVpenta(netpenta,nov+1) = indirp[cnp[nov+1]-1];
      }
      netpenta++;
      cnp += 6+1; ps += 6+1;
    }
    else if (ntype == 17) // HEXA
    {
      for (E_Int nov = 0; nov < 8; nov++)
      {
        et = cnp[nov+1]-1;
        if (indirp[et] == -1)
        {
          for (E_Int eq = 1; eq <= nfld; eq++) fhexa(nptshexa,eq) = field(et,eq);
          nptshexa++;
          indirh[et] = nptshexa;
        }
      }
      for (E_Int nov = 0; nov < 8; nov++)
      {
        cEVhexa(nethexa,nov+1) = indirh[cnp[nov+1]-1];
      }
      nethexa++;
      cnp += 8+1; ps += 8+1;
    }
  }
  printf("found " SF_D_ " TRI\n", nettri);
  printf("found " SF_D_ " QUAD\n", netquad);
  printf("found " SF_D_ " HEXA\n", nethexa);
  printf("found " SF_D_ " PENTA\n", netpenta);
  printf("found " SF_D_ " TETRA\n", nettetra);

  // BAR
  cEVbarp->reAllocMat(netbar,2);
  cEV.push_back(cEVbarp); fields.push_back(fbarp); eltType.push_back(1);
  //TRI
  cEVtrip->reAllocMat(nettri,3);
  cEV.push_back(cEVtrip); fields.push_back(ftrip); eltType.push_back(2);
  //QUAD
  cEVquadp->reAllocMat(netquad,4);
  cEV.push_back(cEVquadp); fields.push_back(fquadp); eltType.push_back(3);
  //TETRA
  cEVtetrap->reAllocMat(nettetra,4);
  cEV.push_back(cEVtetrap); fields.push_back(ftetrap); eltType.push_back(4);
  //PYRA
  cEVpyrap->reAllocMat(netpyra,5);
  cEV.push_back(cEVpyrap); fields.push_back(fpyrap); eltType.push_back(5);
  //PENTA
  cEVpentap->reAllocMat(netpenta,6);
  cEV.push_back(cEVpentap); fields.push_back(fpentap); eltType.push_back(6);
  //HEXA
  cEVhexap->reAllocMat(nethexa,8);
  cEV.push_back(cEVhexap); fields.push_back(fhexap); eltType.push_back(7);
}
