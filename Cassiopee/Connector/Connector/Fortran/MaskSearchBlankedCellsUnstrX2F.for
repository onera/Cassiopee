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

      SUBROUTINE k6searchblankedcellstrix2( npts, 
     &     xt, yt, nelts, cn1, cn2, cn3,
     &     xmin, niray, njray, hiray,
     &     indir, nz, z, 
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite  tri : 1ere colonne
      INTEGER_E cn2(0:nelts-1)  ! connectivite  tri : 2eme colonne
      INTEGER_E cn3(0:nelts-1)  ! connectivite  tri : 3eme colonne
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp, et
      INTEGER_E ind, indray, indrayc, indraycp1
      REAL_E    dx1, dy1
      REAL_E    xp, yp
      REAL_E    xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1
      INTEGER_E iraymin, iraymax
      REAL_E ibmin, ibmax, ibc, ibcp1
      REAL_E xrayc, alpha, xmax
      INTEGER_E npmin, npmax, nnp
      INTEGER_E cellN, nijray, pos
      REAL_E z1, z2, zi
C=============================================================================
      nijray = niray * njray
      isMasked = 0
      xmax = xmin + (niray-1) * hiray
      pos = 0
      DO et = 0, nelts-1
         cellN = 0
c
#include "../../Connector/Fortran/MaskTriCellF.for"
C
C       Barycenter of cell
C     
         xp = xp1 + xp2 + xp3
         yp = yp1 + yp2 + yp3    
         xp = xp * ONE/THREE
         yp = yp * ONE/THREE
         
# include "../../Connector/Fortran/MaskCenterIn2DF.for"

      ENDDO
      END

C ============================================================================
C Cas QUAD
C ============================================================================
      SUBROUTINE k6searchblankedcellsquadx2( npts, 
     &     xt, yt, nelts, cn1, cn2, cn3, cn4,
     &     xmin, niray, njray, hiray,
     &     indir, nz, z, 
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite  tri : 1ere colonne
      INTEGER_E cn2(0:nelts-1)  ! connectivite  tri : 2eme colonne
      INTEGER_E cn3(0:nelts-1)  ! connectivite  tri : 3eme colonne
      INTEGER_E cn4(0:nelts-1)  ! connectivite  tri : 3eme colonne
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp, et
      INTEGER_E ind, indray, indrayc, indraycp1
      REAL_E    dx1, dy1
      REAL_E    xp, yp
      REAL_E    xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1
      INTEGER_E iraymin, iraymax
      REAL_E ibmin, ibmax, ibc, ibcp1
      REAL_E xrayc, alpha, xmax
      INTEGER_E npmin, npmax, nnp
      INTEGER_E cellN, nijray, pos
      REAL_E z1, z2, zi
C=============================================================================
      nijray = niray * njray
      isMasked = 0
      xmax = xmin + (niray-1) * hiray
      pos = 0
      DO et = 0, nelts-1
         cellN = 0
c
#include "../../Connector/Fortran/MaskQuadCellF.for"
C
C       Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4
         yp = yp1 + yp2 + yp3 + yp4    
         xp = xp * ONE_FOURTH
         yp = yp * ONE_FOURTH
         
# include "../../Connector/Fortran/MaskCenterIn2DF.for"

      ENDDO
      END

C  ============================================================================

      SUBROUTINE k6searchblankedcellstetrax2(
     &     npts, xt, yt, zt, nelts, 
     &     cn1, cn2, cn3, cn4, 
     &     xmin, ymin, niray, njray, hiray, hjray,
     &     indir, nz, z, cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      REAL_E    zt(0:npts-1)    ! Height of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite  tetra : 1ere colonne
      INTEGER_E cn2(0:nelts-1)
      INTEGER_E cn3(0:nelts-1) 
      INTEGER_E cn4(0:nelts-1) 
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l      ! For loop in the three directions
      INTEGER_E ip, jp, kp
      INTEGER_E ind, indray, et
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      REAL_E    dx1, dy1, dz1
      REAL_E    xp, yp, zp
      REAL_E    xp1, xp2, xp3, xp4, xp5, xp6, xp7, xp8
      REAL_E    yp1, yp2, yp3, yp4, yp5, yp6, yp7, yp8
      REAL_E    zp1, zp2, zp3, zp4, zp5, zp6, zp7, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1, iraycp2, iraycp3
      INTEGER_E jray, jrayc, jraycp1, jraycp2, jraycp3
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      REAL_E ibmin, ibmax, ibc1, ibc2, ibc3, ibc4
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E npmin, npmax, np
      INTEGER_E cellN, nijray, pos
      REAL_E z1, z2, z3, z4, zi
C==============================================================================
C
      nijray = niray * njray
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 
      pos = 0

      DO et = 0, nelts-1
         cellN = 0
C     
#include "../../Connector/Fortran/MaskTetraCellF.for"
C     
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4 
         yp = yp1 + yp2 + yp3 + yp4
         zp = zp1 + zp2 + zp3 + zp4
         xp = xp * ONE_FOURTH
         yp = yp * ONE_FOURTH
         zp = zp * ONE_FOURTH
         
#include "../../Connector/Fortran/MaskCenterInF.for"

      ENDDO
      END

C  ============================================================================

      SUBROUTINE k6searchblankedcellspentax2(
     &     npts, xt, yt, zt, nelts, 
     &     cn1, cn2, cn3, cn4, cn5, cn6,
     &     xmin, ymin, niray, njray, hiray, hjray,
     &     indir, nz, z, cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      REAL_E    zt(0:npts-1)    ! Height of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite  tetra : 1ere colonne
      INTEGER_E cn2(0:nelts-1)
      INTEGER_E cn3(0:nelts-1) 
      INTEGER_E cn4(0:nelts-1) 
      INTEGER_E cn5(0:nelts-1) 
      INTEGER_E cn6(0:nelts-1)  
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l      ! For loop in the three directions
      INTEGER_E ip, jp, kp
      INTEGER_E ind, indray, et
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      REAL_E    dx1, dy1, dz1
      REAL_E    xp, yp, zp
      REAL_E    xp1, xp2, xp3, xp4, xp5, xp6, xp7, xp8
      REAL_E    yp1, yp2, yp3, yp4, yp5, yp6, yp7, yp8
      REAL_E    zp1, zp2, zp3, zp4, zp5, zp6, zp7, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1, iraycp2, iraycp3
      INTEGER_E jray, jrayc, jraycp1, jraycp2, jraycp3
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      REAL_E ibmin, ibmax, ibc1, ibc2, ibc3, ibc4
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E npmin, npmax, np
      INTEGER_E cellN, nijray, pos
      REAL_E z1, z2, z3, z4, zi
      REAL_E onesix
C==============================================================================
C
      nijray = niray * njray
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 
      pos = 0
      onesix = ONE_HALF/THREE
      DO et = 0, nelts-1
         cellN = 0
C     
#include "../../Connector/Fortran/MaskPentaCellF.for"
C     
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4 + xp5 + xp6
         yp = yp1 + yp2 + yp3 + yp4 + yp5 + yp6
         zp = zp1 + zp2 + zp3 + zp4 + zp5 + zp6
         xp = xp * onesix
         yp = yp * onesix
         zp = zp * onesix
         
#include "../../Connector/Fortran/MaskCenterInF.for"

      ENDDO
      END


C  ============================================================================

      SUBROUTINE k6searchblankedcellshexax2(
     &     npts, xt, yt, zt, nelts, 
     &     cn1, cn2, cn3, cn4, cn5, cn6, cn7, cn8,
     &     xmin, ymin, niray, njray, hiray, hjray,
     &     indir, nz, z, cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      REAL_E    zt(0:npts-1)    ! Height of mesh points
      INTEGER_E nelts           ! nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite  tetra : 1ere colonne
      INTEGER_E cn2(0:nelts-1)
      INTEGER_E cn3(0:nelts-1) 
      INTEGER_E cn4(0:nelts-1) 
      INTEGER_E cn5(0:nelts-1) 
      INTEGER_E cn6(0:nelts-1) 
      INTEGER_E cn7(0:nelts-1) 
      INTEGER_E cn8(0:nelts-1) 
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E l      ! For loop in the three directions
      INTEGER_E ind, indray, et
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      REAL_E    xp, yp, zp
      REAL_E    xp1, xp2, xp3, xp4, xp5, xp6, xp7, xp8
      REAL_E    yp1, yp2, yp3, yp4, yp5, yp6, yp7, yp8
      REAL_E    zp1, zp2, zp3, zp4, zp5, zp6, zp7, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1, iraycp2, iraycp3
      INTEGER_E jray, jrayc, jraycp1, jraycp2, jraycp3
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      REAL_E ibmin, ibmax, ibc1, ibc2, ibc3, ibc4
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E npmin, npmax, np
      INTEGER_E cellN, nijray, pos
      REAL_E z1, z2, z3, z4, zi
C==============================================================================
C
      nijray = niray * njray
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 
      pos = 0

      DO et = 0, nelts-1
         cellN = 0
C     
#include "../../Connector/Fortran/MaskHexaCellF.for"
C     
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4 + xp5 + xp6 + xp7 + xp8
         yp = yp1 + yp2 + yp3 + yp4 + yp5 + yp6 + yp7 + yp8 
         zp = zp1 + zp2 + zp3 + zp4 + zp5 + zp6 + zp7 + zp8
         xp = xp * ONE_EIGHT
         yp = yp * ONE_EIGHT
         zp = zp * ONE_EIGHT
         
#include "../../Connector/Fortran/MaskCenterInF.for"

      ENDDO
      END
