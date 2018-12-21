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
      SUBROUTINE k6homothety( npts, x, y, z,
     &                        xc, yc, zc,
     &                        alpha, 
     &                        xo, yo, zo)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E npts            ! mesh size
      REAL_E x(1:npts)          ! mesh coordinates
      REAL_E y(1:npts)
      REAL_E z(1:npts)
      REAL_E xc,yc,zc           ! center of homothety
      REAL_E alpha              ! factor
C_OUT
      REAL_E xo(1:npts)         ! rotated mesh
      REAL_E yo(1:npts)
      REAL_E zo(1:npts)
C_LOCAL
      INTEGER_E ind
      REAL_E rx, ry, rz
C==============================================================================
!$OMP PARALLEL PRIVATE(ind, rx, ry, rz)
!$OMP DO 
      DO ind = 1, npts
         rx = x(ind)-xc
         ry = y(ind)-yc
         rz = z(ind)-zc         
         xo(ind) = xc + alpha*rx
         yo(ind) = yc + alpha*ry
         zo(ind) = zc + alpha*rz
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END
