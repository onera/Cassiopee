C  
C    Copyright 2013-2025 Onera.
C
C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

      SUBROUTINE k6lamb(x0, y0, Gamma, MInf,
     &                  npts,      
     &                  xc, yc, zc,
     &                  u1, u2, u3, u4, u5)
C
      IMPLICIT NONE 
C
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      REAL_E x0, y0             ! vortex center
      REAL_E Gamma              ! vortex intensity
      REAL_E MInf               ! infinite Mach number
      INTEGER_E npts            ! number of points
      REAL_E xc(0:npts-1) ! x coord of points
      REAL_E yc(0:npts-1) ! y cooord of points
      REAL_E zc(0:npts-1) ! z cooord of points
C_OUT
      REAL_E u1(0:npts-1) ! field to be initialized
      REAL_E u2(0:npts-1) ! field to be initialized
      REAL_E u3(0:npts-1) ! field to be initialized
      REAL_E u4(0:npts-1) ! field to be initialized
      REAL_E u5(0:npts-1) ! field to be initialized

C_LOCAL  
      INTEGER_E nd
      PARAMETER (nd=50000)   !discretisation
      REAL_E pc, roc, ro0, a0
      REAL_E p0, h0, S0
      REAL_E S(2*nd)
      INTEGER_E ind
      REAL_E dr, r, aa
      REAL_E Sinf, roinf, uinf, pinf
      INTEGER_E i,m
      REAL_E rmax, vr, vrdr
      REAL_E cos_teta, sin_teta
      REAL_E va, ss, ro, pp
      REAL_E gam
      REAL_E pi

      pi = ACOS(-1.D0)

      pc = 1.D0/1.4D0           !pression critique
      roc = 1.D0
      
      ro0 = 1.D0                ! etat d'adimensionnement=etat d'arret
      p0 = 1.D0/1.4D0
      a0 = 1.D0
      h0 = 1.D0/0.4D0
      S0 = p0/(ro0**1.4D0)

      gam = Gamma*a0                ! intensite du tourbillon
      
      rmax = 0.
      DO ind = 0, npts-1
         rmax = MAX(rmax, SQRT((xc(ind)-x0)**2+(yc(ind)-y0)**2))
      ENDDO
      
      dr = rmax/(1.*nd)
      S(1) = S0
      
      DO i = 2, 2*nd
         r = (i-1)*dr

         vr = gam*(1-EXP(-r*r))/(2*pi*r)
         vrdr = gam*(1-EXP(-(r+dr)**2))/(2*pi*(r+dr))
         aa = (r+dr)*vrdr-r*vr
                              
         aa = -1.4D0*vr*aa/r
         aa = aa/(h0-0.5*vr*vr)   
         S(i) = LOG(S(i-1))+aa
         S(i) = EXP(S(i))
      ENDDO
                                
      Sinf = S(2*nd)
      roinf = 0.4D0*h0/(1.4D0*Sinf)
      roinf = roinf**(1.D0/0.4D0)
      pinf = roinf*(0.4D0/1.4D0)*(h0)
      uinf = MInf*SQRT(1.4D0*pinf/roinf)
                                
C      WRITE(*,*) 'Lamb: infinite state (conservative variables): ',
C     &     roinf, roinf*uinf, 0., 0.,
C     &     pinf/0.4+0.5*roinf*uinf*uinf

                                ! pression et densite en chaque point

      DO ind = 0, npts-1
         r = SQRT( (xc(ind)-x0)**2+(yc(ind)-y0)**2 )
         IF (r.LE.1.E-8) THEN
            cos_teta = 0.D0
            sin_teta = 0.D0
            va = 0.D0
         ELSE
            cos_teta = (yc(ind)-y0)/r
            sin_teta = -(xc(ind)-x0)/r
            va = gam*(1.D0-EXP(-r*r))/(2*pi*r)
         ENDIF
               
         m = INT(r*nd/rmax)+1
         ss = S(m)
         ro = 0.4D0*(h0-0.5D0*va*va)/(1.4D0*ss)
         ro = ro**(1.D0/0.4D0)
         pp = ro*0.4D0*(h0-0.5D0*va*va)/1.4D0
         
         u1(ind) = ro
         u2(ind) = ro*uinf+ro*cos_teta*va
         u3(ind) = ro*sin_teta*va
         u4(ind) = 0.D0
         u5(ind) = pp/0.4D0+0.5D0*(u2(ind)**2+u3(ind)**2)/ro
      ENDDO
      
      RETURN

      END




