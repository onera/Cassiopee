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

# include "Connect/connect.h"

//=============================================================================
// supIdPoints
// suppress identical points in coord vector
// On suppose que les coordonnees sont en position 1,2,3
// IN: posx, posy, poz: position des variables de coordonnees dans coord
// IN: tol: tolerance pour dire que 2 points sont identiques
//=============================================================================
void K_CONNECT::supIdPoints(K_FLD::FldArrayF& coord, 
                            E_Int posx, E_Int posy, E_Int posz, E_Float tol)
{
  E_Int i;
  E_Float dt, dtx, dty, dtz;
  E_Int ni = coord.getSize(); 
  E_Int nfld = coord.getNfld();
  K_FLD::FldArrayF coord2(ni, nfld);

  E_Float* cx = coord.begin(posx);
  E_Float* cy = coord.begin(posy);
  E_Float* cz = coord.begin(posz);

  E_Int np = 0;
  if (ni > 0)
  {
    for (E_Int n = 1; n <= nfld; n++)
      coord2(0, n) = coord(0, n);
    np++;
  }

  for (i = 1; i < ni; i++)
  {
    dtx = cx[i] - cx[i-1];
    dty = cy[i] - cy[i-1];
    dtz = cz[i] - cz[i-1];
    dt = dtx*dtx + dty*dty + dtz*dtz;
    if (dt > tol)
    { 
      for (E_Int n = 1; n <= nfld; n++)
        coord2(np, n) = coord(i, n);
      np++;
    }
  }

  coord2.reAllocMat(np, nfld);
  coord = coord2;
}
