/*    
    Copyright 2013-2019 Onera.

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

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Passage des données de centres etendus aux noeuds
   IN: FextCenter: champ aux centres etendus
   IN: ime, jme, kme: taille en centres etendus
   IN: im, jm, km: taille aux noeuds
   OUT: FNode: champ aux noeuds
   Retourne 1 en cas de succès, 0 en cas d'échec */
//=============================================================================
E_Int K_LOC::extCenters2NodeStruct(E_Int ime, E_Int jme, E_Int kme,
                                   FldArrayF& FextCenter,
                                   E_Int im, E_Int jm, E_Int km,
                                   FldArrayF& FNode)
{
  E_Int nfld = FextCenter.getNfld();
  E_Int ind;

  E_Int imjm = im*jm;
  E_Int im1 = im-1; E_Int jm1 = jm-1; E_Int km1 = km-1;

  E_Int imejme = ime*jme;
  E_Int ie = 1; if (ime == 1) ie = 0;
  E_Int je = 1; if (jme == 1) je = 0;
  E_Int ke = 1; if (kme == 1) ke = 0;

  E_Int ip,jp,kp,ic,jc,kc;
  E_Int ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8;

  for (E_Int v = 1; v <= nfld; v++)
  {
    E_Float* fc = FextCenter.begin(v);
    E_Float* fn = FNode.begin(v);
      
    for (E_Int k = 0; k < km; k++)
    for (E_Int j = 0; j < jm; j++)
    for (E_Int i = 0; i < im; i++)
    {
      ip = ie; jp = je; kp = ke;
      ic = 0; jc = 0; kc = 0;
      if (i == 0) ip = 0;
      else if (i == im1) ic = ie;
      if (j == 0) jp = 0;
      else if (j == jm1) jc = je;
      if (k == 0) kp = 0;
      else if (k == km1) kc = ke;

      ind = i+j*im+k*imjm;
      ind1 = i+ic+(j+jc)*ime+(k+kc)*imejme;
      ind2 = i+ip+(j+jc)*ime+(k+kc)*imejme;
      ind3 = i+ic+(j+jp)*ime+(k+kc)*imejme;
      ind4 = i+ip+(j+jp)*ime+(k+kc)*imejme;
      ind5 = i+ic+(j+jc)*ime+(k+kp)*imejme;
      ind6 = i+ip+(j+jc)*ime+(k+kp)*imejme;
      ind7 = i+ic+(j+jp)*ime+(k+kp)*imejme;
      ind8 = i+ip+(j+jp)*ime+(k+kp)*imejme;

      fn[ind] = 0.125*(fc[ind1]+fc[ind2]+fc[ind3]+fc[ind4]+
                       fc[ind5]+fc[ind6]+fc[ind7]+fc[ind8]);
    }
  }
  return 1;
}
