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

#include "loc.h"
#include <math.h>
#include "CompGeom/compGeom.h"
//=============================================================================
// Conversion repere Cartesien -> repere cylindrique
// Suivant les axes canoniques
// IN: npts: nbre de pts du maillage
// IN: xt, yt, zt: coord. repere cart
// OUT: rt, thetat: coord. cylindrique 
//=============================================================================
E_Int K_LOC::cart2Cyl(E_Int npts, 
                      E_Float* xt, E_Float* yt, E_Float* zt,
                      E_Float X0, E_Float Y0, E_Float Z0,
                      E_Float ex, E_Float ey, E_Float ez,
                      E_Float* rt, E_Float* thetat, 
                      E_Int ni, E_Int nj, E_Int nk, E_Int depth, 
                      E_Float thetaShift)
{
    E_Float x0, y0;
    E_Float *xl, *yl;
    E_Float eps = 1.e-12;
    //E_Float eps = K_CONST::E_ZERO_MACHINE;
    //E_Float eps = K_CONST::E_GEOM_CUTOFF;
    
    // rotate of thetaShift
    if (thetaShift != 0.)
    {
      E_Float* xDR = new E_Float[npts];
      E_Float* yDR = new E_Float[npts];
      E_Float* zDR = new E_Float[npts];
      
      K_COMPGEOM::rotateMesh2(npts, thetaShift,
			      X0, Y0, Z0,
			      ex, ey, ez,
			      xt, yt, zt,
			      xDR, yDR, zDR);
      xt = xDR; yt = yDR; zt = zDR; // leak
    }
    // Choix direction suivant axe
    if (ex > eps && ey < eps && ez < eps) // axe X
    {
        xl = yt; yl = zt; //zl = xt;
        x0 = Y0; y0 = Z0; //z0 = X0;
    }
    else if (ey > eps && ex < eps && ez < eps) // axe Y
    {
        xl = zt; yl = xt; //zl = yt;
        x0 = Z0; y0 = X0; //z0 = Y0;
    }
    else if (ez > eps && ey < eps && ex < eps) // axe Z
    {
        xl = xt; yl = yt; //zl = zt;
        x0 = X0; y0 = Y0; //z0 = Z0;
    }
    else 
    { 
      // Not a canonical axis
      return 1; // FAILED
    }
    // Maintenant axe Z
    //E_Float thetaref = atan2(yl[0]-y0,xl[0]-x0);

#pragma omp parallel default(shared)
    {
      E_Float dx, dy, r;
      E_Float theta;
#pragma omp for 
      for (E_Int ind = 0; ind < npts; ind++)                   
      {
        dx = xl[ind]-x0;
        dy = yl[ind]-y0;
        r = sqrt(dx*dx+dy*dy);
        theta = atan2(dy,dx);
        rt[ind] = r; thetat[ind] = theta;
      }
    }

    // firewall
    if (depth == 0) return 0;

    // i-lines
    if (ni > 0 && nj >0 && nk > 0)
    {
        E_Float DEUXPI = 2*K_CONST::E_PI; 
        E_Float DELTATHETAMAX = 0.75*DEUXPI;
        E_Int ninj = ni*nj;

#pragma omp parallel default(shared)
    {
        E_Int ind, indm; 
        E_Float theta, thetap;
#pragma omp for     
        for (E_Int k = 0; k < nk; k++)
            for (E_Int j = 0; j < nj; j++)
            {
                //cas i = depth
                ind = depth+j*ni+k*ninj; E_Int indp = ind+1;
                theta = thetat[ind];  
                thetap = thetat[indp];  

            if (K_FUNC::E_abs(theta-thetap)>DELTATHETAMAX)
            {
                if ( theta < thetap) thetat[ind] = theta+DEUXPI;
                else thetat[ind] = theta-DEUXPI;
            }

            for (E_Int i = depth+1; i < ni; i++)       
            {
                ind = i+j*ni+k*ninj; theta = thetat[ind];  
                indm = ind-1; E_Float thetam = thetat[indm];  
                if (K_FUNC::E_abs(theta-thetam)>DELTATHETAMAX)
                {   
                   if (theta < thetam) thetat[ind] = theta+DEUXPI;
                   else thetat[ind] = theta-DEUXPI;
                }  
            }   

            // depth premieres rangees
            if (depth > 0)
            {
                for (E_Int i = depth-1; i >= 0; i--)
                {
                    ind = i+j*ni+k*ninj; E_Int indp = ind+1;
                    theta = thetat[ind];  
                    thetap = thetat[indp]; 
                    if (K_FUNC::E_abs(theta-thetap)>DELTATHETAMAX)
                    {   
                        if ( theta < thetap) thetat[ind] = theta+DEUXPI;
                        else thetat[ind] = theta-DEUXPI;
                    }                     
                }
            }
        }
#pragma omp for     
        for (E_Int k = 0; k < nk; k++)
            for (E_Int i = 0; i < ni; i++)
            {
            //cas j = depth
            ind = i+depth*ni+k*ninj; E_Int indp = ind+ni;
            theta = thetat[ind];  
            thetap = thetat[indp];  

            if (K_FUNC::E_abs(theta-thetap)>DELTATHETAMAX)
            {
                if (theta < thetap) thetat[ind] = theta+DEUXPI;
                else thetat[ind] = theta-DEUXPI;
            }

            for (E_Int j = depth+1; j < nj; j++)       
            {
                ind = i+j*ni+k*ninj; theta = thetat[ind];  
                indm = ind-ni; E_Float thetam = thetat[indm];  
                if (K_FUNC::E_abs(theta-thetam)>DELTATHETAMAX)
                {   
                    if (theta < thetam) thetat[ind] = theta+DEUXPI;
                    else thetat[ind] = theta-DEUXPI;
                }  
            }    
            //deux premieres rangees
            if (depth > 0)
            {
                for (E_Int j = depth-1; j >= 0; j--)
                {
                    ind = i+j*ni+k*ninj; E_Int indp = ind+ni;
                    theta = thetat[ind];  
                    thetap = thetat[indp]; 
                    if (K_FUNC::E_abs(theta-thetap)>DELTATHETAMAX)
                    {   
                        if (theta < thetap) thetat[ind] = theta+DEUXPI;
                        else thetat[ind] = theta-DEUXPI;
                    }                     
                }                
            }
        }

        for (E_Int j = 0; j < nj; j++)
            for (E_Int i = 0; i < ni; i++)
            {
                if (depth+1 < nk)
                {
                //cas k = depth
                ind = i+j*ni+depth*ninj; E_Int indp = ind+ninj;
                theta = thetat[ind];  
                thetap = thetat[indp];  

                if (K_FUNC::E_abs(theta-thetap)>DELTATHETAMAX)
                {
                    if (theta < thetap) thetat[ind] = theta+DEUXPI;
                    else thetat[ind] = theta-DEUXPI;
                }

                for (E_Int k = depth+1; k < nk; k++)       
                {
                    ind = i+j*ni+k*ninj; theta = thetat[ind];  
                    indm = ind-ninj; E_Float thetam = thetat[indm];  
                    if (K_FUNC::E_abs(theta-thetam)>DELTATHETAMAX)
                    {   
                        if (theta < thetam) thetat[ind] = theta+DEUXPI;
                        else thetat[ind] = theta-DEUXPI;
                    }  
                }    
            
            // depth premieres rangees
            if (depth > 0)
            {
                for (E_Int k = depth-1; k >= 0; k--)
                {
                    ind = i+j*ni+k*ninj; E_Int indp = ind+ninj;
                    theta = thetat[ind];  
                    thetap = thetat[indp]; 
                    if (K_FUNC::E_abs(theta-thetap)>DELTATHETAMAX)
                    {   
                        if (theta < thetap) thetat[ind] = theta+DEUXPI;
                        else thetat[ind] = theta-DEUXPI;
                    }                     
                }
            }    
            }               
        }    
    } 
    }
    return 0; // OK
}
