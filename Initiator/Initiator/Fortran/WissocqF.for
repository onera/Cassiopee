C
C    Copyright 2013-2020 Onera.
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

      SUBROUTINE k6wissocq(x0, y0, Gamma, MInf,
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
      REAL_E MInf               ! infinite Mach number
      INTEGER_E npts            ! number of points
      REAL_E xc(0:npts-1) ! x coord of points
      REAL_E yc(0:npts-1) ! y cooord of points
      REAL_E zc(0:npts-1) ! z cooord of points
C_OUT
      REAL_E u(0:npts-1,5) ! field to be initialized
C_LOCAL
      INTEGER_E ind
      REAL_E r, rc, rc2
      REAL_E roinf, uinf
      REAL_E cos_teta, sin_teta
      REAL_E va, ro, pp
      REAL_E gam
      REAL_E c0
      REAL_E coef
      REAL_E pinf, tinf

      roinf = 1.1765
      pinf = 101320
      tinf = pinf/(roinf*287.053)

      c0 = SQRT(1.4*287.053*tinf)
      gam = Gamma*c0                ! intensite du tourbillon
      coef = 0.5D0*(gam**2)/(c0**2)
      uinf = MInf*c0
      rc = 0.1D0  !1D0
      rc2 = rc*rc

      DO ind = 0, npts-1
         r = SQRT( (xc(ind)-x0)**2+(yc(ind)-y0)**2 )
         cos_teta = -(yc(ind)-y0)/rc
         sin_teta = (xc(ind)-x0)/rc
         va = gam*EXP(-0.5D0*(r*r)/rc2)

         ! Masse volumique et pression
         ro = roinf*EXP(-coef*EXP(-(r*r)/rc2))
         pp = roinf*287.053*tinf + (ro-roinf)*(c0**2)
         !Variables conservatives
         u(ind,1) = ro
         u(ind,2) = ro*uinf+ro*cos_teta*va
         u(ind,3) = ro*sin_teta*va
         u(ind,4) = 0.D0
         u(ind,5) = pp/0.4D0 + 0.5D0*(u(ind,2)**2+u(ind,3)**2)/ro
      ENDDO

      RETURN

      END
