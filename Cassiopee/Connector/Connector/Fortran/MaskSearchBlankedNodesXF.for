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

C XRay 3D + node in 

      SUBROUTINE k6searchblankednodesx(
     &     npts, meshX, meshY, meshZ,
     &     xmin, ymin,
     &     niray, njray,
     &     hiray, hjray,
     &     indir, nz, z, isnot,
     &     cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! Dimension of the domain
      REAL_E    meshX(0:npts-1) ! Abcisses of mesh points
      REAL_E    meshY(0:npts-1) ! Ordinates of mesh points
      REAL_E    meshZ(0:npts-1) ! Height of mesh points
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot           ! inverse mask

C_OUT
      INTEGER_E cellNatureField(0:npts-1) ! Give the nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E l       
      INTEGER_E et, indray
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xp7, yp7, zp7, xp8, yp8, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nijray
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E iray, irayc, iraycp1, iraycp2, iraycp3
      INTEGER_E jray, jrayc, jraycp1, jraycp2, jraycp3
      REAL_E ibmin, ibmax, ibc1, ibc2, ibc3, ibc4
      REAL_E    xp, yp, zp
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      INTEGER_E np, pos
      REAL_E z1, z2, z3, z4, zi

C==============================================================================

      nijray = niray * njray
      isMasked = 0
      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 

      DO et = 0, npts-1
         cellN = 0
         xmincell = meshX(et)
         ymincell = meshY(et)
         zmincell = meshZ(et)
         xmaxcell = meshX(et)
         ymaxcell = meshY(et)
         zmaxcell = meshZ(et)
         xp = xmincell
         yp = ymincell
         zp = zmincell

#include "../../Connector/Fortran/MaskNodeInF.for"
      ENDDO

      IF (isnot .NE. 0) THEN   
         DO et = 0, npts-1
            cellNatureField(et) = -cellNatureField(et)
         ENDDO
      ENDIF

      END


