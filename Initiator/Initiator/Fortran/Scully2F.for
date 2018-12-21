C  
C    Copyright 2013-2019 Onera.
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

C ============================================================================

      SUBROUTINE k6scully2(x0, y0, Gamma, a, MInf,
     &                     npts,
     &                     xc, yc, zc,
     &                     u)
C
      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      REAL_E x0, y0             ! vortex center
      REAL_E Gamma              ! vortex intensity
      REAL_E a                  ! vortex core radius
      REAL_E MInf               ! infinite Mach number
      INTEGER_E npts            ! number of centers
      REAL_E xc(0:npts-1)    ! x coord of centers
      REAL_E yc(0:npts-1)    ! y cooord of centers
      REAL_E zc(0:npts-1)    ! z cooord of centers
C_OUT
      REAL_E u(0:npts-1,5)     ! field to be initialized
C_LOCAL
      INTEGER_E ind
      REAL_E r, p, ptmp
      REAL_E roinf, uinf, pinf
      INTEGER_E i,j,k
      REAL_E cos_teta, sin_teta
      REAL_E va, ro, pi
C==============================================================================
      pi = 4*ATAN(1.D0)

C     Etat d'adimensionnement      
      roinf = 1.D0
      uinf = 1.D0
      pinf = 1.D0/(1.4D0*MInf*MInf)

C     Initialisation

      ptmp = (0.4D0/1.4D0)*(Gamma*Gamma)/(8*pi*pi*pinf)

      DO ind = 0, npts-1
               
         r = SQRT( (xc(ind)-x0)**2+(yc(ind)-y0)**2 )

         IF (r.LE.1.E-12) THEN
            cos_teta = 0.D0
            sin_teta = 0.D0
            va = 0.D0
         ELSE
            cos_teta = (yc(ind)-y0)/r	
            sin_teta = -(xc(ind)-x0)/r
            va = Gamma/(2.*pi)*r/(r*r+a*a)
         ENDIF

         p = pinf*((1-ptmp/(r*r+a*a))**(1.4D0/0.4D0))
         ro = ((1-ptmp/(r*r+a*a))**(1.D0/0.4D0))*roinf

         u(ind,1) = ro
         u(ind,2) = ro*uinf + ro*cos_teta*va
         u(ind,3) = ro*sin_teta*va
         u(ind,4) = 0.D0
         u(ind,5) = p/0.4D0+0.5D0*(u(ind,2)**2+u(ind,3)**2)/ro
      ENDDO
      
      RETURN

      END











