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

#include "Connect/connect.h"
#include "String/kstring.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/*
  Change a MIX connectivity to a set of Elts-Vertex connectivities
  IN: cMIX: MIX connectivity. For each elt, give type + vertices index.
  OUT: cEV: set of Vertex-Elts connectivities. 
  Les indices de vertex commencent a 1
  Les indices d'elements commencent a 0
*/
//=============================================================================
void K_CONNECT::connectMix2EV(FldArrayI& cMIX,
                              FldArrayI& cBAR,
                              FldArrayI& cTRI,
                              FldArrayI& cQUAD,
                              FldArrayI& cTETRA,
                              FldArrayI& cPYRA,
                              FldArrayI& cPENTA,
                              FldArrayI& cHEXA)
{
  // Calcul du nombre d'elements de chaque type
  E_Int nBAR = 0;
  E_Int nQUAD = 0;
  E_Int nTRI = 0;
  E_Int nPENTA = 0;
  E_Int nTETRA = 0;
  E_Int nHEXA = 0;
  E_Int nPYRA = 0;

  E_Int sizeTot = cMIX.getSize()*cMIX.getNfld();
  E_Int* mix = cMIX.begin();
  E_Int ntype;
  E_Int ls; E_Int* lp;

  E_Int c = 0;
  while (c < sizeTot)
  {
    ntype = mix[c];
    switch (ntype)
    {
      case 3:
        nBAR++; ls = 2; break;

      case 5:
        nTRI++; ls = 3; break;

      case 7:
        nQUAD++; ls = 4; break;

      case 10:
        nTETRA++; ls = 4; break;

      case 12:
        nPYRA++; ls = 5; break;

      case 14:
        nPENTA++; ls = 6; break;

      case 17:
        nHEXA++; ls = 8; break;

      default: // probleme
        printf("Warning: Mix2Ev: Unknow element type (" SF_D_ ").\n", ntype); 
        ls = 0; break;
    }
    c += ls+1;
  }

  printf("Info: Mix2EV: found " SF_D_ " HEXAS " SF_D_ " PENTAS\n", nHEXA, nPENTA);

  // Dimensionne
  cBAR.malloc(nBAR*2);
  cTRI.malloc(nTRI*3);
  cQUAD.malloc(nQUAD*4);
  cTETRA.malloc(nTETRA*4);
  cPYRA.malloc(nPYRA*5);
  cPENTA.malloc(nPENTA*6);
  cHEXA.malloc(nHEXA*8);

  // Remplit
  E_Int* ptBAR = cBAR.begin();
  E_Int* ptTRI = cTRI.begin();
  E_Int* ptQUAD = cQUAD.begin();
  E_Int* ptTETRA = cTETRA.begin();
  E_Int* ptPYRA = cPYRA.begin();
  E_Int* ptPENTA = cPENTA.begin();
  E_Int* ptHEXA = cHEXA.begin();

  c = 0;
  while (c < sizeTot)
  {
    ntype = mix[c];
    switch (ntype)
    {
      case 3:
        ls = 2; lp = ptBAR; ptBAR += 2; break;

      case 5:
        ls = 3; lp = ptTRI; ptTRI += 3; break;

      case 7:
        ls = 4; lp = ptQUAD; ptQUAD += 4; break;

      case 10:
        ls = 4; lp = ptTETRA; ptTETRA += 4; break;

      case 12:
        ls = 5; lp = ptPYRA; ptPYRA += 5; break;

      case 14:
        ls = 6; lp = ptPENTA; ptPENTA += 6; break;

      case 17:
        ls = 8; lp = ptHEXA; ptHEXA += 8; break;

      default:
        printf("Warning: Mix2Ev: unknown element type.\n");
        ls = 0; lp = NULL; break;

    }
    for (E_Int i = 0; i < ls; i++) lp[i] = mix[c+i+1];
    c += ls+1;
  }
}
