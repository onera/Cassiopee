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

# include "Connect/connect.h"
# include "metric.h"
# include <vector>

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

//=============================================================================
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
//=============================================================================
E_Float K_METRIC::compVolOfStructCell2D(E_Int ni, E_Int nj, 
                                        E_Float* xt, E_Float* yt, E_Float* zt,
                                        E_Int indcell, E_Int indnode)
{
  E_Int i, j;
  E_Int nic = max(1,ni-1);
  //E_Int njc = max(1,nj-1);

  if ( indcell >  -1)
  {
    if ( indnode < 0 )
    {
      j = indcell/nic;
      i = indcell -j*nic;
      indnode = i+j*ni;
    }
  }
  else
  {
    if ( indnode < 0)
    {
      printf("INTERNAL ERROR: compVolOfStructCell2D: one of indcell or indnode must be a positive value.");
      exit(0);
    }
  }
  // E_Int j = ind/ni; E_Int i = ind-j*ni;
  // if (i==ni-1) i=i-1;
  // if (j==nj-1) j=j-1;
  // E_Int ind1 = i+j*ni;
  E_Int ind1 = indnode;
  E_Int ind2 = ind1+ni;
  E_Int ind3 = ind1+1;
  E_Int ind4 = ind2+1;
  E_Float l1x = xt[ind1]-xt[ind2];
  E_Float l1y = yt[ind1]-yt[ind2];
  E_Float l1z = zt[ind1]-zt[ind2];
  E_Float l2x = xt[ind1]-xt[ind3];
  E_Float l2y = yt[ind1]-yt[ind3];
  E_Float l2z = zt[ind1]-zt[ind3];
 
  E_Float surf1x = (l1y*l2z-l1z*l2y);
  E_Float surf1y = (l1z*l2x-l1x*l2z);
  E_Float surf1z = (l1x*l2y-l1y*l2x);
  E_Float surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z);
  E_Float l3x = xt[ind4]-xt[ind3];
  E_Float l3y = yt[ind4]-yt[ind3];
  E_Float l3z = zt[ind4]-zt[ind3];
  E_Float l4x = xt[ind4]-xt[ind2];
  E_Float l4y = yt[ind4]-yt[ind2];
  E_Float l4z = zt[ind4]-zt[ind2]; 
  E_Float surf2x = (l3y*l4z-l3z*l4y);
  E_Float surf2y = (l3z*l4x-l3x*l4z);
  E_Float surf2z = (l3x*l4y-l3y*l4x);
  E_Float surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z);
  E_Float ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z;
  return 0.5*K_FUNC::E_sign(ps)*(surface1+surface2);
}
