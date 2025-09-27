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
#include <vector>
#include <stdio.h>
#include <string.h>

using namespace std;
using namespace K_CONST;
using namespace K_FUNC;
using namespace K_FLD;

extern "C"
{
  void k6stretch_(const E_Float& t1, const E_Float& t2,
                  E_Float* sn,
                  const E_Int& nbp, const E_Float& dsm,
                  const E_Float& dsp,
                  const E_Int& ityp, const E_Int& inewt);
}

void checkDistribution(FldArrayF& sn, 
                       E_Int& croissante, E_Int& monotonic, 
                       E_Float& regularity);

//=============================================================================
PyObject* K_GENERATOR::enforceMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int supp; E_Int add;
  char* name;
  E_Float eh; // enforce length
  E_Float P0;
  
  if (!PYPARSETUPLE_(args, O_ S_ RR_ II_,
                     &array, &name, &P0, &eh, &supp, &add)) return NULL;
 
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  if (res == 1)
  {
    FldArrayF* out = new FldArrayF();
    E_Int niout, njout, nkout;
    E_Int api = f->getApi();
    E_Int ret = enforceCommon(name, varString, ni, nj, nk, 
                              *f, P0, eh, supp, add, 
                              *out, niout, njout, nkout);
    RELEASESHAREDS(array, f);
    if (ret != 0) return NULL;
    PyObject* tpl = K_ARRAY::buildArray3(*out, "x,y,z", niout, njout, nkout, api);
    delete out;
    return tpl;
  }
  if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);    
    PyErr_SetString(PyExc_TypeError, 
                    "enforce: not for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "enforce: invalid array.");
    return NULL;
  }
}

//=============================================================================
// Enforce: partie commune aux enforceX,...
// IN: name: enforceX, enforceY, enforceZ, enforcePlusX, ...
// IN: auto: si autoAdd=1, ajuste automatiquement add en fonction de la
// regularite, si autoAdd=0, utilise add en entree.
// Retourne un code d'erreur: 
// 0: OK
// 1: coordonnees non presente dans array
// 2: fonction inconnue (name inconnu)
//=============================================================================
E_Int K_GENERATOR::enforceCommon(const char* name, char* varString, 
                E_Int ni, E_Int nj, E_Int nk,
                FldArrayF& coord, E_Float P0, E_Float eh, 
                E_Int supp, E_Int add, FldArrayF& out,
                E_Int& niout, E_Int& njout, E_Int& nkout,
                E_Int autoAdd)
{
  char msg[256*8];
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    sprintf(msg, "%s: can't find coordinates in array.", name);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  posx++; posy++; posz++;

  //--------------------------------------------------------
  // Parametres suivant les differents types de distribution
  //--------------------------------------------------------
  E_Int pos = 0;
  E_Int side = 0;
  E_Int bsup, binf, bsupm, binfp, nic, njc, nkc, N, P;
  if (strcmp(name, "enforceX") == 0)
  {
    pos = posx;
    side = 0; // both sides
    binf = 0;
    binfp = 1;
    bsup = ni-1;
    bsupm = ni-2;
    nic = ni;
    njc = 1;
    nkc = 1;
    N = ni;
    P = 1;
  }
  else if (strcmp(name, "enforceMoinsX") == 0)
  {
    pos = posx;
    side = -1; // minus side
    binf = 0;
    binfp = 1;
    bsup = ni-1;
    bsupm = ni-2;
    nic = ni;
    njc = 1;
    nkc = 1;
    N = ni;
    P = 1;
  }
  else if (strcmp(name, "enforcePlusX") == 0)
  {
    pos = posx;
    side = +1; // plus side
    binf = 0;
    binfp = 1;
    bsup = ni-1;
    bsupm = ni-2;
    nic = ni;
    njc = 1;
    nkc = 1;
    N = ni;
    P = 1;
  }
  else if (strcmp(name, "enforceY") == 0)
  {
    pos = posy;
    side = 0; // both sides
    binf = 0;
    binfp = ni;
    bsupm = (nj-2)*ni;
    bsup = (nj-1)*ni;
    nic = 1;
    njc = nj;
    nkc = 1;
    N = nj;
    P = ni;
  }
  else if (strcmp(name, "enforceMoinsY") == 0)
  {
    pos = posy;
    side = -1; // minus side
    binf = 0;
    binfp = ni;
    bsupm = (nj-2)*ni;
    bsup = (nj-1)*ni;
    nic = 1;
    njc = nj;
    nkc = 1;
    N = nj;
    P = ni;
  }
  else if (strcmp(name, "enforcePlusY") == 0)
  {
    pos = posy;
    side = +1; // plus side
    binf = 0;
    binfp = ni;
    bsupm = (nj-2)*ni;
    bsup = (nj-1)*ni;
    nic = 1;
    njc = nj;
    nkc = 1;
    N = nj;
    P = ni;
  }
  else if (strcmp(name, "enforceZ") == 0)
  {
    pos = posz;
    side = 0; // both sides
    binf = 0;
    binfp = ni*nj;
    bsupm = (nk-2)*ni*nj;
    bsup = (nk-1)*ni*nj;
    nic = 1;
    njc = 1;
    nkc = nk;
    N = nk;
    P = ni*nj;
  }
  else if (strcmp(name, "enforceMoinsZ") == 0)
  {
    pos = posz;
    side = -1; // minus side
    binf = 0;
    binfp = ni*nj;
    bsupm = (nk-2)*ni*nj;
    bsup = (nk-1)*ni*nj;
    nic = 1;
    njc = 1;
    nkc = nk;
    N = nk;
    P = ni*nj;
  }
  else if (strcmp(name, "enforcePlusZ") == 0)
  {
    pos = posz;
    side = +1; // plus side
    binf = 0;
    binfp = ni*nj;
    bsupm = (nk-2)*ni*nj;
    bsup = (nk-1)*ni*nj;
    nic = 1;
    njc = 1;
    nkc = nk;
    N = nk;
    P = ni*nj;
  }
  else
  {
    sprintf(msg, "%s: unknown function name.", name);
    PyErr_SetString(PyExc_TypeError, msg);
    return 2;
  }

  //----------------------------------------------------------------
  // Determination de l'indice de modification : il (seulement pour
  // les modifications side == 0)
  //----------------------------------------------------------------
  E_Int il=0;
  E_Int i=-1; E_Int j=-1; E_Int k=-1;
  E_Float hl=P0;

  E_Float* coordp = coord.begin(pos);

  if (side == 0)
  {
    // Determination de la cellule de coord contenant hl
    if ( (hl <= coordp[binf]) || (hl >= coordp[bsup]) )
    {
      sprintf(msg, 
              "%s (P0=" SF_F_ ", eh=" SF_F_ "): cannot find P0 in array, stopped.",
              name, P0, eh);
      PyErr_SetString(PyExc_TypeError, msg);
      return 3;
    }

    //printf("hl=" SF_F3_ "\n", hl, coordp[bsupm], eh);

    if ( hl <= coordp[binfp] )
    {
      side = +1; // force enforcePlus
    }
    else if ( hl >= coordp[bsupm] )
    {
      side = -1; // force enforceMoins
    }
    else if ( hl <= coordp[binfp] + eh/2. )
    {
      side = +1; // force enforcePlus
    }
    else if ( hl >= coordp[bsupm] - eh/2. )
    {
      side = -1; // force enforceMoins
    }
    else
    {
      E_Float h;
      E_Int ind = 1;
      for (k = 0; k < nkc; k++)
        for (j = 0; j < njc; j++)
          for (i = 0; i < nic; i++)
          {
            ind = i + j*ni + k*ni*nj;
            h = coordp[ind];
            if (h > hl) goto exit;
          }
      exit:

      if (pos == posx) il = i-1;
      else if (pos == posy) il = j-1;
      else if (pos == posz) il = k-1;
    }
  }

  //cout << "il : "<<il<<" side : "<<side<<" "<<ni<<endl;

  E_Int compt = 0;
  E_Int cmax = 10*N;
  if (autoAdd == 0) cmax = 1;
  E_Int suppl = supp;
  E_Int suppr = supp;
  if (side == 0)
  {
    suppl = E_min(suppl, il);
    suppr = E_min(suppr, N-1-il-1);
  }
  else if (side == 1) suppr = E_min(suppr, N-2);
  else suppl = E_min(suppl, N-2);

  //printf("supp " SF_D2_ "\n", suppl, suppr);

  E_Int addl, addr;
  if (cmax == 1) // pas de recherche
  {
    addl = add;
    addr = add;
  }
  else
  {
    addl = 3-suppl; // reglage minimum pour la recherche auto
    addr = 3-suppr;
  }
  E_Int istart = 0;
  E_Int iend = 0;
  FldArrayF snl; FldArrayF snl2;
  FldArrayF snr; FldArrayF snr2;
  FldArrayF sn1; FldArrayF sn2;
  E_Float regularity, regularity2, regularityBest;
  E_Int croissante, monotonic, croissante2, monotonic2;
  E_Int addrBest=0, addlBest=0;

  E_Float* xt = coord.begin(posx);
  E_Float* yt = coord.begin(posy);
  E_Float* zt = coord.begin(posz);
  
  regularityBest = 1.e6;

  while (compt < cmax)
  {
    if (side == 0) // =======================================
    {      
      E_Int addln = addl;
      E_Int addrn = addr;

      istart = il - suppl;
      iend = il + suppr + 1;

      E_Float pt1 = coordp[istart*P];
      E_Float pt1a = coordp[istart*P+P];
      E_Float pt2a = coordp[iend*P-P];
      E_Float pt2 = coordp[iend*P];
      E_Float pt3 = hl-eh/2.;
      E_Float pt4 = hl+eh/2.;
      E_Float deltal = pt1a - pt1;
      E_Float deltar = pt2 - pt2a;

      // Distribution a gauche
      E_Int npl = suppl + addl + 2;
      snl.malloc(npl); snl2.malloc(npl);
      k6stretch_(pt1, pt4, snl2.begin(), npl, eh, deltal, 2, 1);

      for (E_Int pp = 0; pp < npl; pp++)
        snl[pp] = -snl2[npl-1-pp] + pt4 + pt1;
      checkDistribution(snl, croissante, monotonic, regularity);

      // Distribution a droite
      E_Int npr = suppr + addr + 2;
      snr.malloc(npr);
      k6stretch_(pt3, pt2, snr.begin(), npr, eh, deltar, 2, 1);
      checkDistribution(snr, croissante2, monotonic2, regularity2);

      if (addl > addr) { addrn = addr + 1; }
      else { addln = addl + 1; }
      
      regularity = E_max(regularity, regularity2);
      if (regularity < regularityBest)
      { 
        regularityBest = regularity;
        if (monotonic == 0) regularityBest += 1000.;
        addrBest = addr; addlBest = addl;
      }

      if (compt == cmax-1) 
      {
        compt = cmax;
        //  printf("Warning: %s: the distribution steps are not monotonic.\n",
        //         name);
        addr = addrBest; addl = addlBest;
        npl = suppl + addl + 2;
        snl.malloc(npl); snl2.malloc(npl);
        k6stretch_(pt1, pt4, snl2.begin(), npl, eh, deltal, 2, 1);
        for (E_Int pp = 0; pp < npl; pp++)
          snl[pp] = -snl2[npl-1-pp] + pt4 + pt1;
        npr = suppr + addr + 2;
        snr.malloc(npr);
        k6stretch_(pt3, pt2, snr.begin(), npr, eh, deltar, 2, 1);
        printf("Info: %s: regularity " SF_F_ ".\n", name, regularityBest);

        if (pos == posx)
        {
          E_Int np = ni + addr + addl;
          out.malloc(np*nj*nk, 3);
          out.setAllValuesAtNull();
          
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < istart ; i++)
              {
                E_Int ind0 = i + j*ni + k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = istart; i < istart+npl-1; i++)
              {
                E_Int ind0 = j*ni + k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = snl[i-istart];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }  
          
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = istart+npl-1; i < istart+npl-1+npr-1; i++)
              {
                E_Int ind0 = j*ni + k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = snr[i-istart-npl+2];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
         
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = istart+npl-1+npr-1; i < ni+addr+addl; i++)
              {
                E_Int ind0 = i-addr-addl+j*ni+k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = ni+addr+addl; njout = nj; nkout = nk;
        }
        else if (pos == posy)
        {
          E_Int np = nj + addr + addl;
          out.malloc(ni*np*nk, 3);
          out.setAllValuesAtNull();
          
          for (k = 0; k < nk; k++)
            for (j = 0; j < istart; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = 0; k < nk; k++)
            for (j = istart; j < istart+npl-1 ; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = snl[j-istart];
                out(ind,3) = zt[ind0];
              }  
          
          for (k = 0; k < nk; k++)
            for (j = istart+npl-1; j < istart+npl-1+npr-1; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = snr[j-istart-npl+2];
                out(ind,3) = zt[ind0];
            }
         
          for (k = 0; k < nk; k++)
            for (j = istart+npl-1+npr-1; j < nj+addr+addl; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i+(j-addr-addl)*ni+k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = ni; njout = nj+addr+addl; nkout = nk;
        }
        else if (pos == posz)
        {
          E_Int np = nk + addr + addl;
          out.malloc(ni*nj*np, 3);
          out.setAllValuesAtNull();
          
          for (k = 0; k < istart; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = istart; k < istart+npl-1 ; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = snl[k-istart];
              }  
          
          for (k = istart+npl-1; k < istart+npl-1+npr-1; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = snr[k-istart-npl+2];
            }
         
          for (k = istart+npl-1+npr-1; k < nk+addr+addl; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i+j*ni+(k-addr-addl)*ni*nj;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = ni; njout = nj; nkout = nk+addr+addl;
        }
      }
      addr = addrn;
      addl = addln;
    }

    else if (side == +1)  //======================================
    {
      istart = 0;
      iend = (istart + suppr + 1)*P;
      E_Int addrn = addr;

      E_Float delta1 = coord(iend, pos) - coord(iend-P, pos);
      E_Int np1 = suppr+addr+2;
      sn1.malloc(np1);

      E_Float pt1 = coord(istart, pos);
      E_Float pt2 = coord(iend, pos);
      // Taille de la premiere maille = eh
      // taille de la derniere maille inchangee = x(iend) - x(iend-1)
      k6stretch_(pt1, pt2, sn1.begin(), np1, eh, delta1, 2, 1);
      checkDistribution(sn1, croissante, monotonic, regularity);
      //printf("pt1 " SF_F2_ " -> " SF_F2_ " (npts=" SF_D_ ") -> reg=" SF_F_ "," SF_D_ "\n",pt1,pt2,eh,delta1,np1,regularity,monotonic);
      addrn = addr + 1;

      if (regularity < regularityBest)
      { 
        regularityBest = regularity;
        if (monotonic == 0) regularityBest += 1000.;
        addrBest = addr;
      }

      if (compt == cmax-1) // Generation reelle
      {
        compt = cmax;
        //  printf("Warning: %s: the distribution steps are not monotonic.\n",
        //         name);
        addr = addrBest;
        np1 = suppr+addr+2;
        sn1.malloc(np1);
        k6stretch_(pt1, pt2, sn1.begin(), np1, eh, delta1, 2, 1);
        printf("Info: %s: regularity " SF_F_ ".\n", name, regularityBest);

        if (pos == posx)
        { 
          E_Int np = ni + addr;
          out.malloc(np*nj*nk, 3);
          out.setAllValuesAtNull();
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < np1; i++)
              {
                E_Int ind0 = j*ni + k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = sn1[i];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = np1; i < np; i++)
              {
                E_Int ind0 = i-addr+j*ni+k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = np; njout = nj; nkout = nk;
        }
        else if (pos == posy)
        {
          E_Int np = nj + addr;
          out.malloc(ni*np*nk, 3);
          out.setAllValuesAtNull();
          for (k = 0; k < nk; k++)
            for (j = 0; j < np1; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = sn1[j];
                out(ind,3) = zt[ind0];

              }
          
          for (k = 0; k < nk; k++)
            for (j = np1; j < np; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + (j-addr)*ni+k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = ni; njout = np; nkout = nk;
        }
        else if (pos == posz)
        {
          E_Int np = nk + addr;
          out.malloc(ni*nj*np, 3);
          out.setAllValuesAtNull();
          for (k = 0; k < np1; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = sn1[k];
              }
          
          for (k = np1; k < np; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni + (k-addr)*ni*nj;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = ni; njout = nj; nkout = np;
        }
      }
      addr = addrn;
    }
    else if (side == -1) // ================================
    {
      iend = N-1;
      istart = iend - suppl - 1;
      E_Int addln = addl;

      E_Float delta1 = coord(istart*P+P, pos) - coord(istart*P, pos);
      E_Int np1 = suppl+addl+2;
      sn1.malloc(np1); sn2.malloc(np1);

      E_Float pt1 = coord(istart*P, pos);
      E_Float pt2 = coord(iend*P, pos);
      // Taille de la premiere maille = eh ; 
      // taille de la derniere maille inchangee = x(iend) - x(iend-1)
      k6stretch_(pt1, pt2, sn2.begin(), np1, eh, delta1, 2, 1);
      for (i = 0; i < np1; i++) sn1[i] = -sn2[np1-1-i]+pt2+pt1;

      checkDistribution(sn1, croissante, monotonic, regularity);
     
      if (croissante == 0 || monotonic == 0) addln = addl + 1;

      if (regularity < regularityBest)
      { 
        regularityBest = regularity;
        if (monotonic == 0) regularityBest += 1000.;
        addlBest = addl;
      }
      
      if (compt == cmax-1) // Generation reelle
      {
        compt = cmax;
        //  printf("Warning: %s: the distribution steps are not monotonic.\n",
        //         name);
        addl = addlBest;
        np1 = suppl+addl+2;
        sn1.malloc(np1); sn2.malloc(np1);
        k6stretch_(pt1, pt2, sn2.begin(), np1, eh, delta1, 2, 1);
        for (i = 0; i < np1; i++)
          sn1[i] = -sn2[np1-1-i]+pt2+pt1;
        printf("Info: %s: regularity " SF_F_ ".\n", name, regularityBest);

        if (pos == posx)
        { 
          E_Int np = ni + addl;
          out.malloc(np*nj*nk, 3);
          out.setAllValuesAtNull();
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < istart; i++)
              {
                E_Int ind0 = i + j*ni + k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              for (i = istart; i < np; i++)
              {
                E_Int ind0 = j*ni+k*ni*nj;
                E_Int ind = i + j*np + k*np*nj;
                out(ind,1) = sn1[i-istart];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          niout = np; njout = nj; nkout = nk;
        }
        else if (pos == posy)
        {
          E_Int np = nj + addl;
          out.malloc(ni*np*nk, 3);
          out.setAllValuesAtNull();
          for (k = 0; k < nk; k++)
            for (j = 0; j < istart; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = 0; k < nk; k++)
            for (j = istart; j < np; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*np;
                out(ind,1) = xt[ind0];
                out(ind,2) = sn1[j-istart];
                out(ind,3) = zt[ind0];
              }
          niout = ni; njout = np; nkout = nk;
        }
        else if (pos == posz)
        {
          E_Int np = nk + addl;
          out.malloc(ni*nj*np, 3);
          out.setAllValuesAtNull();
          for (k = 0; k < istart; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni + k*ni*nj;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = zt[ind0];
              }
          
          for (k = istart; k < np; k++)
            for (j = 0; j < nj; j++)
              for (i = 0; i < ni; i++)
              {
                E_Int ind0 = i + j*ni;
                E_Int ind = i + j*ni + k*ni*nj;
                out(ind,1) = xt[ind0];
                out(ind,2) = yt[ind0];
                out(ind,3) = sn1[k-istart];
              }
          niout = ni; njout = nj; nkout = np;
        }
      }
      addl = addln;
    }
      
    compt++;
  }

  return 0;
}

//=============================================================================
// Verifie si la distribution est valide
// IN: np: nbre de pts de la distribution
// IN: sn: distribution
// OUT: croissante: = 1 si la distribution est croissante
// OUT: monotonic: = 1 si le pas de la distribution est croissante ou 
// decroissante tout le temps
// OUT: regularity: retourne la regularite la plus mauvaise
//=============================================================================
void checkDistribution(FldArrayF& sn, 
                       E_Int& croissante, E_Int& monotonic, 
                       E_Float& regularity)
{
  E_Int i;
  E_Int np = sn.getSize();
  // Croissance des points
  croissante = 1;
  for (i = 0; i < np-1; i++)
  {
    if (sn[i+1] <= sn[i]) { croissante = 0; break; }
  }

  // Verification de la distribution: croissance/decroissance de la taille 
  // de maille 
  monotonic = 1;
  regularity = 1.;
  E_Float eh = sn[1]-sn[0];
  E_Float delta1 = sn[np-1]-sn[np-2];
  E_Float reg, h1, h2;
  if (eh <= delta1)
  {
    for (i = 0; i < np-2; i++)
    {
      h2 = sn[i+2]-sn[i+1]; h1 = sn[i+1]-sn[i];
      if (h2 > h1) reg = h2/E_max(h1, 1.e-12);
      else reg = h1/E_max(h2, 1.e-12);
      regularity = E_max(regularity, reg);
      if (h2 < h1) { monotonic = 0; }
    }
  }
  else
  {
    for (i = 0; i < np-2; i++)
    {
      h2 = sn[i+2]-sn[i+1]; h1 = sn[i+1]-sn[i];
      if (h2 > h1) reg = h2/E_max(h1, 1.e-12);
      else reg = h1/E_max(h2, 1.e-12);
      regularity = E_max(regularity, reg);
      if (h2 > h1) { monotonic = 0; }
    }
  }
}


