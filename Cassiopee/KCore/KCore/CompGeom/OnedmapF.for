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
 
C ============================================================================
C map a 1D distribution over a profile
C  ==========================================================================
      SUBROUTINE k6onedmap(ni, x, y, z, no, d, xo, yo, zo,
     &                     s, dx, dy, dz)
C=============================================================================
      IMPLICIT NONE
C_IN
      INTEGER_E ni              ! nombre de points de la ligne entrante
      REAL_E x(ni),y(ni),z(ni)  ! ligne entrante
      INTEGER_E no              ! nombre de points de la ligne sortante
      REAL_E    d(no)           ! distribution
C_OUT
      REAL_E xo(no), yo(no), zo(no) ! ligne sortante
C_LOCAL
      REAL_E small, big, stota
      REAL_E s(ni)
      REAL_E dx(ni)
      REAL_E dy(ni)
      REAL_E dz(ni)
C==============================================================================
      small = 1.e-18
      big = 1.e+15

C     Parametrisation	de la ligne entrante
      CALL k6slope(small, big, ni,
     &               x, y, z,            
     &               dx, dy, dz)

       CALL k6param(stota, small, ni, 
     &  	     x, y, z,          
     &               dx, dy, dz, s)
            
C       Projection
       CALL k6interp(ni, no, stota,       
     &                x, y, z, s,       
     &                dx, dy, dz,      
     &                xo, yo, zo, d)                

      END


C=============================================================================
      SUBROUTINE k6onedmapbar(npts, x, y, z, no, d, net, cn1, cn2, neto,
     &                        cn1o, cn2o, xo, yo, zo, s, dx, dy, dz)

      IMPLICIT NONE
C_IN
      INTEGER_E npts          ! nombre de points de la ligne entrante
      REAL_E x(0:npts-1)      ! ligne entrante
      REAL_E y(0:npts-1)      ! ligne entrante
      REAL_E z(0:npts-1)      ! ligne entrante
      INTEGER_E no            ! nombre de points de la ligne sortante
      REAL_E d(0:no-1)        ! distribution
      INTEGER_E net           ! nb d'elts BAR initiaux
      INTEGER_E cn1(0:net-1)  ! connectivite BAR : 1ers sommets
      INTEGER_E cn2(0:net-1)  ! connectivite BAR : 2nds sommets
C_OUT
      INTEGER_E neto          ! nb d'elts BAR  de sortie
      INTEGER_E cn1o(0:neto-1) ! connectivite BAR : 1ers sommets
      INTEGER_E cn2o(0:neto-1) ! connectivite BAR : 2nds sommets
      REAL_E xo(0:no-1), yo(0:no-1), zo(0:no-1) ! ligne sortante
C_LOCAL
      REAL_E small, big, stota
      REAL_E s(0:npts-1)
      REAL_E dx(0:npts-1)
      REAL_E dy(0:npts-1)
      REAL_E dz(0:npts-1)

      small = 1.e-18
      big = 1.e+15

C       Parametrisation	de la ligne entrante
      CALL k6slopebar(small, big, npts, x, y, z, net, cn1, cn2,
     &                  dx, dy, dz)

      CALL k6parambar(stota, small, npts, x, y, z, net, cn1, cn2, 
     &                  dx, dy, dz, s)
        
C       Projection
      CALL k6interpbar(npts, no, stota, net, cn1, cn2, x, y, z, s,
     &       dx, dy, dz, neto, cn1o, cn2o, xo, yo, zo, d)
               
      END
C==================  OnedmapF.for ============================================
