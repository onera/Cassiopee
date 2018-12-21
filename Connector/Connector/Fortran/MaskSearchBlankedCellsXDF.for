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

C Search the cells which will be blanked "within" a distance of XRay mask.     
C Put the results in a cell nature field array.
C A cell is blanked if the bounding box of cell (+/- delta) intersects 
C the XRay mask (+/- delta).

      SUBROUTINE k6searchblankedcellsxd( ni, nj, nk, 
     &     meshX, meshY, meshZ,
     &     xmin, ymin,
     &     niray, njray,
     &     hiray, hjray,
     &     indir, nz, z, 
     &     listOfInterpolatedPoints, np,
     &     delta, isnot,
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
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:(ni-1)*(nj-1)*(nk-1)-1) ! cellN
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E ind, et, indray
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xp7, yp7, zp7, xp8, yp8, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E n, npmin, npmax
      INTEGER_E cellN, nij, nijray, nic, nicnjc

C==============================================================================
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN 
      DO n = 0, np-1
         cellN = 0
C     
         et = listOfInterpolatedPoints(n)
         k = et/nicnjc
         j = (et-k*nicnjc)/nic
         i = et-k*nicnjc-j*nic
C        
#include "../../Connector/Fortran/MaskCell3DF.for"
C
         xmincell = xmincell-delta
         ymincell = ymincell-delta
         xmaxcell = xmaxcell+delta
         ymaxcell = ymaxcell+delta
#include "../../Connector/Fortran/MaskCellIntersectF.for"

      ENDDO

      ELSE

      DO n = 0, np-1
         cellN = 0
C     
         et = listOfInterpolatedPoints(n)
         k = et/nicnjc
         j = (et-k*nicnjc)/nic
         i = et-k*nicnjc-j*nic
C        
#include "../../Connector/Fortran/MaskCell3DF.for"
C
         xmincell = xmincell-delta
         ymincell = ymincell-delta
         xmaxcell = xmaxcell+delta
         ymaxcell = ymaxcell+delta
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"

      ENDDO
      ENDIF
      END


