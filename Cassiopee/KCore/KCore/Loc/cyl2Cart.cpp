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

//=============================================================================
// Conversion repere Cylindrique -> repere cartesien
// Suivant les axes canoniques
// IN: npts: nbre de pts du maillage
// IN: rt, thetat: coord. cylindrique 
// OUT: xt, yt, zt: coord. repere cart
//=============================================================================
E_Int K_LOC::cyl2Cart(E_Int npts, E_Float* rt, E_Float* thetat, 
                      E_Float X0, E_Float Y0, E_Float Z0,
                      E_Float ex, E_Float ey, E_Float ez,
                      E_Float* xt, E_Float* yt, E_Float* zt)
{
    E_Float x0, y0;
    E_Float* xl; E_Float* yl; //E_Float* zl;
    E_Float eps = 1.e-12;

    // Choix direction suivant axe
    if (ex > eps && ey < eps && ez < eps) // axe X
    {
      xl = yt; yl = zt; //zl = xt;
      x0 = Y0; y0 = Z0;
    }
    else if (ey > eps && ex < eps && ez < eps) // axe Y
    {
      xl = zt; yl = xt; //zl = yt;
      x0 = Z0; y0 = X0;
    }
    else if (ez > eps && ey < eps && ex < eps) // axe Z
    {
      xl = xt; yl = yt; //zl = zt;
      x0 = X0; y0 = Y0;
    }
    else 
    { 
      // Not a canonical axis
      return 1; // FAILED
    }
    // Maintenant axe Z
  
#pragma omp parallel default(shared)
    {
        E_Float r;
        E_Float theta;
#pragma omp for 
      for (E_Int ind = 0; ind < npts; ind++)                   
      {
        r = rt[ind];
        theta = thetat[ind];
        xl[ind] = r*cos(theta) + x0;
        yl[ind] = r*sin(theta) + y0; 
      }
    }
    return 0; // OK
}
