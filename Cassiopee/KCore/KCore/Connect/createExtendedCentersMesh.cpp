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
using namespace std;

//=============================================================================
/* Creation du maillage en centres etendus a partir du maillage en noeuds 
   le maillage est alloue ici */
//=============================================================================
void 
K_CONNECT::createExtendedCentersMesh(E_Int imo, E_Int jmo, E_Int kmo, 
                                     E_Float* xt, E_Float* yt, E_Float* zt,
                                     E_Int& ime, E_Int& jme, E_Int& kme, 
                                     K_FLD::FldArrayF& extendedCentersMesh)
{
  E_Int i, j, k;
  E_Int pos, pos2, pos2a, pos2b, pos2c, pos2d, pos2e, pos2f, pos2g;
  if (extendedCentersMesh.getSize() == 0) 
  {
    ime = imo+1; jme = jmo+1; kme = kmo+1;
    if (imo == 1) ime = 1;
    if (jmo == 1) jme = 1;
    if (kmo == 1) kme = 1;
    extendedCentersMesh.malloc(ime*jme*kme, 3);
  }
  E_Int imjme = ime*jme;
  E_Int imjmo = imo*jmo;
  E_Float* xe = extendedCentersMesh.begin(1);
  E_Float* ye = extendedCentersMesh.begin(2);
  E_Float* ze = extendedCentersMesh.begin(3);

  /* corners */
  pos  = 0; pos2 = 0;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];

  pos  = ime-1; pos2 = imo-1;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];
  
  pos  = (jme-1)*ime; pos2 = (jmo-1)*imo;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];

  pos = (kme-1)*imjme; pos2 = (kmo-1)*imjmo;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];

  pos  = ime-1 + (jme-1)*ime; pos2 = imo-1 + (jmo-1)*imo;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2]; 

  pos = ime-1 + (kme-1)*imjme; pos2 = imo-1 + (kmo-1)*imjmo;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];

  pos = (jme-1)*ime + (kme-1)*imjme; pos2 = (jmo-1)*imo + (kmo-1)*imjmo;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];
  
  pos = ime-1 + (jme-1)*ime + (kme-1)*imjme; pos2 = imo-1 + (jmo-1)*imo + (kmo-1)*imjmo;
  xe[pos] = xt[pos2];
  ye[pos] = yt[pos2];
  ze[pos] = zt[pos2];
  
  /* border lines */
  for (i = 2; i < ime; i++)
  {
    pos   = i-1; pos2  = i-1; pos2a = i-2;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);

    pos   = i-1 + (jme-1)*ime;
    pos2  = i-1 + (jmo-1)*imo;//origMesh.getPos(i,jmo,1);
    pos2a = i-2 + (jmo-1)*imo;// origMesh.getPos(i-1,jmo,1);
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);

    pos   = i-1 + (kme-1)*imjme;//getPos(i,1,kme);
    pos2  = i-1 + (kmo-1)*imjmo;//origMesh.getPos(i,1,kmo);
    pos2a = pos2 -1 ;//origMesh.getPos(i-1,1,kmo);
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);

    pos   = i-1 + (jme-1)*ime + (kme-1)*imjme;//getPos(i,jme,kme);
    pos2  = i-1 + (jmo-1)*imo + (kmo-1)*imjmo;//origMesh.getPos(i,jmo,kmo);
    pos2a = pos2-1;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);
  }
  
  for (j = 2; j < jme; j++)
  {
    pos   = (j-1)*ime;
    pos2  = (j-1)*imo; 
    pos2a = pos2 - imo;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);

    pos   = ime-1 + (j-1)*ime;
    pos2  = imo-1 + (j-1)*imo; 
    pos2a = pos2 - imo;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);
    
    pos   = (j-1)*ime + (kme-1)*imjme;//getPos(1,j,kme);
    pos2  = (j-1)*imo + (kmo-1)*imjmo;//origMesh.getPos(1,j,kmo);
    pos2a = pos2 - imo;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);

    pos   = (ime-1) + (j-1)*ime + (kme-1)*imjme;//getPos(ime,j,kme);
    pos2  = (imo-1) + (j-1)*imo + (kmo-1)*imjmo;//origMesh.getPos(imo,j,kmo);
    pos2a = pos2 - imo;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);
  }

  for (k = 2; k < kme; k++)
  {
    pos   = (k-1)*imjme; //getPos(1,1,k);
    pos2  = (k-1)*imjmo;//origMesh.getPos(1,1,k);
    pos2a = pos2 - imjmo;//origMesh.getPos(1,1,k-1);
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);
  
    pos   = (ime-1) + (k-1) * imjme;//getPos(ime,1,k);
    pos2  = (imo-1) + (k-1) * imjmo;//origMesh.getPos(imo,1,k);
    pos2a = pos2 - imjmo;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);

    pos   = (jme-1)*ime + (k-1)*imjme;//getPos(1,jme,k);
    pos2  = (jmo-1)*imo + (k-1)*imjmo;//origMesh.getPos(1,jmo,k);
    pos2a = pos2 - imjmo;//origMesh.getPos(1,jmo,k-1);
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);
    
    pos   = ime-1 + (jme-1)*ime + (k-1)*imjme;//getPos(ime,jme,k);
    pos2  = imo-1 + (jmo-1)*imo + (k-1)*imjmo;//origMesh.getPos(imo,jmo,k);
    pos2a = pos2 - imjmo;
    xe[pos] = K_CONST::ONE_HALF*(xt[pos2] + xt[pos2a]);
    ye[pos] = K_CONST::ONE_HALF*(yt[pos2] + yt[pos2a]);
    ze[pos] = K_CONST::ONE_HALF*(zt[pos2] + zt[pos2a]);
  }
  
  /* inside plan */
  for (j = 2; j < jme; j++)
    for (i = 2; i < ime; i++)
    {
      pos   = i-1+(j-1)*ime;
      pos2  = i-1+(j-1)*imo;//origMesh.getPos(i,j,1);
      pos2a = pos2-1;//origMesh.getPos(i-1,j,1);
      pos2b = pos2-imo;//sorigMesh.getPos(i,j-1,1);
      pos2c = pos2-1-imo;//origMesh.getPos(i-1,j-1,1);
      xe[pos] =
        K_CONST::ONE_FOURTH*( xt[pos2]  + xt[pos2a]
                     + xt[pos2b] + xt[pos2c]);
      ye[pos] =
        K_CONST::ONE_FOURTH*( yt[pos2]  + yt[pos2a]
                     + yt[pos2b] + yt[pos2c]);
      ze[pos] =
        K_CONST::ONE_FOURTH*( zt[pos2]  + zt[pos2a]
                     + zt[pos2b] + zt[pos2c]);

      pos   = i-1+(j-1)*ime + (kme-1)*imjme;
      pos2  = i-1+(j-1)*imo + (kmo-1)*imjmo;
      pos2a = pos2-1; pos2b = pos2-imo; pos2c = pos2-1-imo;
      xe[pos] =
        K_CONST::ONE_FOURTH*(  xt[pos2]  + xt[pos2a]
                      + xt[pos2b] + xt[pos2c]);
      ye[pos] =
        K_CONST::ONE_FOURTH*(  yt[pos2]  + yt[pos2a]
                      + yt[pos2b] + yt[pos2c]);
      ze[pos] =
        K_CONST::ONE_FOURTH*(  zt[pos2]  + zt[pos2a]
                      + zt[pos2b] + zt[pos2c]);
    }

  for (k = 2; k < kme; k++)
    for (i = 2; i < ime; i++)
    {
      pos   = i-1 + (k-1)*imjme;//getPos(i,1,k);
      pos2  = i-1 + (k-1)*imjmo;//origMesh.getPos(i,1,k);
      pos2a = pos2-1;//rigMesh.getPos(i-1,1,k);
      pos2b = pos2-imjmo;//origMesh.getPos(i,1,k-1);
      pos2c = pos2-1-imjmo;//origMesh.getPos(i-1,1,k-1);
      xe[pos] =
        K_CONST::ONE_FOURTH*(  xt[pos2]  + xt[pos2a]
                      + xt[pos2b] + xt[pos2c]);
      ye[pos] =
        K_CONST::ONE_FOURTH*(  yt[pos2]  + yt[pos2a]
                      + yt[pos2b] + yt[pos2c]);
      ze[pos] =
        K_CONST::ONE_FOURTH*(  zt[pos2]  + zt[pos2a]
                      + zt[pos2b] + zt[pos2c]);
      pos   = (i-1) + (jme-1)*ime + (k-1)*imjme;//getPos(i,jme,k);
      pos2  = (i-1) + (jmo-1)*imo + (k-1)*imjmo;//origMesh.getPos(i,jmo,k);
      pos2a = pos2-1;
      pos2b = pos2-imjmo;//origMesh.getPos(i,jmo,k-1);
      pos2c = pos2-1-imjmo;//origMesh.getPos(i-1,jmo,k-1);
      xe[pos] =
        K_CONST::ONE_FOURTH*(  xt[pos2]  + xt[pos2a]
                      + xt[pos2b] + xt[pos2c]);
      ye[pos] =
        K_CONST::ONE_FOURTH*(  yt[pos2]  + yt[pos2a]
                      + yt[pos2b] + yt[pos2c]);
      ze[pos] =
        K_CONST::ONE_FOURTH*(  zt[pos2]  + zt[pos2a]
                      + zt[pos2b] + zt[pos2c]);
    }

  for (k = 2; k < kme; k++)
    for (j = 2; j < jme; j++)
    {
      pos   = (j-1)*ime + (k-1)*imjme;//getPos(1,j,k);
      pos2  = (j-1)*imo + (k-1)*imjmo;//origMesh.getPos(1,j,k);
      pos2a = pos2-imjmo;
      pos2b = pos2-imo;
      pos2c = pos2-imo-imjmo;
      xe[pos] =
        K_CONST::ONE_FOURTH*(  xt[pos2]  + xt[pos2a]
                      + xt[pos2b] + xt[pos2c]);
      ye[pos] =
        K_CONST::ONE_FOURTH*(  yt[pos2]  + yt[pos2a]
                      + yt[pos2b] + yt[pos2c]);
      ze[pos] =
        K_CONST::ONE_FOURTH*(  zt[pos2]  + zt[pos2a]
                      + zt[pos2b] + zt[pos2c]);

      pos   = (ime-1) + (j-1)*ime + (k-1)*imjme;
      pos2  = (imo-1) + (j-1)*imo + (k-1)*imjmo;//origMesh.getPos(1,j,k);
      pos2a = pos2-imjmo;
      pos2b = pos2-imo;
      pos2c = pos2-imo-imjmo;
      xe[pos] =
        K_CONST::ONE_FOURTH*(  xt[pos2]  + xt[pos2a]
                      + xt[pos2b] + xt[pos2c]);
      ye[pos] =
        K_CONST::ONE_FOURTH*(  yt[pos2]  + yt[pos2a]
                      + yt[pos2b] + yt[pos2c]);
      ze[pos] =
        K_CONST::ONE_FOURTH*(  zt[pos2]  + zt[pos2a]
                      + zt[pos2b] + zt[pos2c]);
    }

  /* current points */
  for (k = 2; k < kme; k++)
    for (j = 2; j < jme; j++)
      for (i = 2; i < ime; i++)
      {
        pos   = i-1 + (j-1)*ime + (k-1)*imjme;
        pos2  = i-1 + (j-1)*imo + (k-1)*imjmo;
        pos2a = pos2-1;
        pos2b = pos2-imo;
        pos2c = pos2-imjmo;
        pos2d = pos2-1-imo;
        pos2e = pos2-1-imjmo;
        pos2f = pos2-1-imo-imjmo;
        pos2g = pos2-imo-imjmo;
        xe[pos] =
          K_CONST::ONE_EIGHTH*(xt[pos2] + xt[pos2a]
                       + xt[pos2b] + xt[pos2c]
                       + xt[pos2d] + xt[pos2e]
                       + xt[pos2f] + xt[pos2g]);
        ye[pos] =
          K_CONST::ONE_EIGHTH*(yt[pos2] + yt[pos2a]
                       + yt[pos2b] + yt[pos2c]
                       + yt[pos2d] + yt[pos2e]
                       + yt[pos2f] + yt[pos2g]);
        ze[pos] =
          K_CONST::ONE_EIGHTH*( zt[pos2] + zt[pos2a]
                               + zt[pos2b] + zt[pos2c]
                               + zt[pos2d] + zt[pos2e]
                               + zt[pos2f] + zt[pos2g]);      
      }
}
