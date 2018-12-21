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

C  ============================================================================
C  Compute velocity  Field
C  ============================================================================

      SUBROUTINE k6computevelocity(npts, eqnb, 
     &                             posro, posrou, posrov, posrow,
     &                             w, velo)
      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E   npts ! Total number of pts 
      INTEGER_E   eqnb  ! Number of Equations
      INTEGER_E   posro ! position of rho in the w field
      INTEGER_E   posrou, posrov, posrow ! position of rou, rov, row in w
      REAL_E      w(npts, eqnb)      ! Conservative Variables

C_OUT
      REAL_E      velo(npts, 3) ! Velocity Field 
      
C_LOCAL
      INTEGER_E   l
      REAL_E      ro, roi
c------------------------------------------------------------------------------
!$OMP PARALLEL PRIVATE(l, ro, roi)
!$OMP DO
      DO l = 1, npts
         ro  = w(l, posro)
         roi = 1.D0 / ro
         velo(l, 1) = w(l, posrou) * roi
         velo(l, 2) = w(l, posrov) * roi
         velo(l, 3) = w(l, posrow) * roi
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
C ==== Post/Fortran/computeVelocityF.for ====
