C
C    Copyright 2013-2026 ONERA.
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

C Cas XRay3D + cell intersect

      SUBROUTINE k6searchblankedcellsx(ni, nj, nk,
     &     meshX, meshY, meshZ,
     &     xmin, ymin,
     &     niray, njray,
     &     hiray, hjray,
     &     indir, nz, z, isnot,
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk      ! Dimension of the domain
      REAL_E    meshX(0:ni*nj*nk-1) ! Abcisses of mesh points
      REAL_E    meshY(0:ni*nj*nk-1) ! Ordinates of mesh points
      REAL_E    meshZ(0:ni*nj*nk-1) ! Height of mesh points
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot           ! inverse mask

C_OUT
      INTEGER_E cellNatureField(0:(ni-1)*(nj-1)*(nk-1)-1) ! Give the nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l, d
      INTEGER_E ip, jp, kp
      INTEGER_E ind, indray
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xp7, yp7, zp7, xp8, yp8, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nij, nijray, nic, nicnjc
      INTEGER_E et
C==============================================================================

      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

!$OMP PARALLEL PRIVATE(d, i, j, k, l, ip, jp, kp, et,
!$OMP&   ind, indray, dx1, dy1, dz1,
!$OMP&   xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3, xp4, yp4, zp4,
!$OMP&   xp5, yp5, zp5, xp6, yp6, zp6, xp7, yp7, zp7, xp8, yp8, zp8,
!$OMP&   xmincell, ymincell, zmincell, xmaxcell, ymaxcell, zmaxcell,
!$OMP&   iray, jray, ibeg, iend, iraymin, iraymax, jraymin, jraymax,
!$OMP&   npmin, npmax, cellN)
      IF (isnot .EQ. 0) THEN
!$OMP DO REDUCTION(MAX:isMasked)
      DO d = 0, nicnjc*(nk-1)-1
          k = d/nicnjc
          j = (d - k*nicnjc)/nic
          i = d - k*nicnjc - j*nic
          cellN = 0
C
#include "../../Connector/Fortran/MaskCell3DF.for"
C
          et = i+j*nic+k*nicnjc
#include "../../Connector/Fortran/MaskCellIntersectF.for"

      ENDDO
!$OMP END DO

      ELSE
!$OMP DO REDUCTION(MAX:isMasked)
      DO d = 0, nicnjc*(nk-1)-1
          k = d/nicnjc
          j = (d - k*nicnjc)/nic
          i = d - k*nicnjc - j*nic
          cellN = 0
C
#include "../../Connector/Fortran/MaskCell3DF.for"
C
          et = i+j*nic+k*nicnjc
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"

      ENDDO
!$OMP END DO

      ENDIF
!$OMP END PARALLEL
      END
C     ===== Connector/Fortran/MaskSearchBlankedCellsXF.for =====
