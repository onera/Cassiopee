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

C Initialization of Visbal vortex in field with density constant
C ============================================================================

      SUBROUTINE k6visbal(x0, y0, Gamma, MInf,
     &                    npts,    
     &                    xc, yc, zc,
     &                    u1, u2, u3, u4, u5)
C
      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      REAL_E x0, y0            ! vortex center
      REAL_E Gamma              ! vortex intensity
      REAL_E MInf              ! infinite Mach number
      INTEGER_E npts            ! number of centers
      REAL_E xc(0:npts-1)    ! x coord of centers
      REAL_E yc(0:npts-1)    ! y cooord of centers
      REAL_E zc(0:npts-1)    ! z cooord of centers
C_OUT
      REAL_E u1(0:npts-1) ! field to be initialized
      REAL_E u2(0:npts-1) ! field to be initialized
      REAL_E u3(0:npts-1) ! field to be initialized
      REAL_E u4(0:npts-1) ! field to be initialized
      REAL_E u5(0:npts-1) ! field to be initialized

C_LOCAL  
      
      REAL_E  ro0, ainf,hinf,ar,ardr
      REAL_E p0,  u0
      REAL_E S(12000)
      INTEGER_E ind
      REAL_E dr, r
      REAL_E Sinf, roinf, uinf, pinf
      INTEGER_E i,m
      REAL_E rmax, vr, vrdr
      REAL_E cos_teta, sin_teta
      REAL_E va, ss, ro, pp
      REAL_E gam
C==============================================================================

      roinf = 1.                ! etat d'adimensionnement : etat inf
      pinf = 1./1.4
      ainf = 1.                 
      hinf = 1./(0.4)           ! Constante
      uinf= MInf
      gam = Gamma*ainf          ! intensite du tourbillon

      rmax = 0.
      DO ind = 0, npts-1
         rmax = MAX( rmax, SQRT( (xc(ind)-x0)**2+(yc(ind)-y0)**2 ) )
      ENDDO

      dr = rmax/(1.*6000)
      r = dr 
      vr = gam/(2*3.14)*EXP(1/2.0)*r*EXP(-(r**2)/2.0)
      S(1) = EXP((1.4/0.4)*(vr**2/(r*(hinf-0.5*vr**2)))*dr/2.0)

      DO i = 2,12000
         r = (i-1)*dr
         vr = gam/(2*3.14)*EXP(1/2.0)*r*EXP(-(r**2)/2.0)
         vrdr = gam/(2*3.14)*EXP(1/2.0)*(r+dr)*EXP(-((r+dr)**2)/2.0)
         ar= vr**2/(r*(hinf-0.5*vr**2))
         ardr=vrdr**2/((r+dr)*(hinf-0.5*vrdr**2))
         S(i) = LOG(S(i-1))+((1.4/0.4)*(ardr+ar)/2)*dr
         S(i) = EXP(S(i))
      ENDDO
                                
      Sinf = S(12000)
      p0 = pinf/Sinf
      ro0 = (1.4*p0)/(0.4*hinf)
      u0 = MInf

      ! pression et densite en chaque point
      DO ind = 0, npts-1
         r = SQRT( (xc(ind)-x0)**2+(yc(ind)-y0)**2 )
         IF (r.LE.1.E-8) THEN
            cos_teta = 0.
            sin_teta = 0.
            va = 0.
         ELSE
            cos_teta = (yc(ind)-y0)/r
            sin_teta = -(xc(ind)-x0)/r
            va = gam/(2*3.14)*EXP(1/2.0)*r*EXP(-r**2/2)
         ENDIF
               
         m = INT(r*1.*6000/rmax)+1
         ss = S(m)
         pp = ss*p0
         ro = (1.4*pp)/(0.4*(hinf-0.5*va**2))
               
         u1(ind) = ro
         u2(ind) = ro*uinf+ro*cos_teta*va
         u3(ind) = ro*sin_teta*va
         u4(ind) = 0.
         u5(ind) = pp/0.4+0.5*(u2(ind)**2+u3(ind)**2)/ro
      ENDDO
      
      RETURN

      END











