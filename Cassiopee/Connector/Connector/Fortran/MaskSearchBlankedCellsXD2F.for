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

C Search the cells which will be blanked "within" a distance of XRay mask.
C Put the results in a cell nature field array.
C      
C A cell is blanked if the bounding box of cell (+/- delta) 
C intersects the XRay mask (+/- delta)
C and the center of cell (+/- delta) is in the mask (+/- delta). 
C depth = 2 case.

      SUBROUTINE k6searchblankedcellsxd2(ni, nj, nk, 
     &     meshX, meshY, meshZ,
     &     xmin, ymin,
     &     niray, njray,
     &     hiray, hjray,
     &     indir, nz, z, 
     &     listOfInterpolatedPoints, np,
     &     delta,
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj, nk   ! Dimension of the domain
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

C_OUT
      INTEGER_E cellNatureField(0:MAX(ni-1,1)*MAX(nj-1,1)*MAX(nk-1,1)-1)
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l      ! For loop in the three directions
      INTEGER_E ip, jp, kp
      INTEGER_E ind, et, indray
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
      INTEGER_E irayminp, iraymaxp, jrayminp, jraymaxp
      REAL_E ibmin, ibmax, ibc, ibcp1, ibcp2, ibcp3, ib
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E npmin, npmax, nnp
      INTEGER_E cellN, nij, nijray, nic, nicnjc
      INTEGER_E n, tmp, ntemp1, ntemp2

C==============================================================================
      tmp = 0
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 
      
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
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4 + xp5 + xp6 + xp7 + xp8
         yp = yp1 + yp2 + yp3 + yp4 + yp5 + yp6 + yp7 + yp8 
         zp = zp1 + zp2 + zp3 + zp4 + zp5 + zp6 + zp7 + zp8
         xp = xp * ONE_EIGHT
         yp = yp * ONE_EIGHT
         zp = zp * ONE_EIGHT
               
         xmincell = xmincell-delta
         ymincell = ymincell-delta
         xmaxcell = xmaxcell+delta
         ymaxcell = ymaxcell+delta
                  
#include "../../Connector/Fortran/MaskCenterInDeltaF.for"

      ENDDO
      END
