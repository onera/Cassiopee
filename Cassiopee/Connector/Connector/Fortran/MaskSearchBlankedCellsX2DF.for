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

C Recherche des points masques cas : plan XRay en (x,y)
C                                    critere de masquage cell intersect

      SUBROUTINE k6searchblankedcellsx2d(
     &     ni, nj, nk, meshX, meshY,
     &     xmin, niray, njray, hiray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

      IMPLICIT NONE
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk   ! Dimension of the domain
      REAL_E    meshX(0:ni*nj*nk-1) ! Abcisses of mesh points
      REAL_E    meshY(0:ni*nj*nk-1) ! Ordinates of mesh points
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:(ni-1)*(nj-1)*nk-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l, d
      INTEGER_E ip, jp, et
      INTEGER_E ind, indray
      REAL_E    dx1, dy1
      REAL_E    xp1, xp2, xp3, xp4
      REAL_E    yp1, yp2, yp3, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E iray, ibeg, iend
      INTEGER_E iraymin, iraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nij, nijray, nic, nicnjc
C==============================================================================
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0


!$OMP PARALLEL PRIVATE(d, i, j, k, l, ip, jp, et, ind, indray, dx1, dy1,
!$OMP&                  xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4,
!$OMP&                  xmincell, ymincell, xmaxcell, ymaxcell,
!$OMP&                  iray, ibeg, iend, iraymin, iraymax,
!$OMP&                  npmin, npmax, cellN)

      IF ( isnot .EQ. 0 ) THEN
!$OMP DO REDUCTION(MAX:isMasked)
      DO d = 0, nicnjc-1
          j = d/nic
          i = d - j*nic
          cellN = 0
          et = i+j*nic
C
#include "../../Connector/Fortran/MaskCell2DF.for"
C
#include "../../Connector/Fortran/MaskCellIntersect2DF.for"
      ENDDO
!$OMP END DO

      ELSE
!$OMP DO REDUCTION(MAX:isMasked)
      DO d = 0, nicnjc-1
          j = d/nic
          i = d - j*nic
          cellN = 0
          et = i+j*nic
C
#include "../../Connector/Fortran/MaskCell2DF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectNot2DF.for"
      ENDDO
!$OMP END DO

      ENDIF
!$OMP END PARALLEL
      END
C ===== XRay/MaskSearchBlankedCellsX2DF.for =====
