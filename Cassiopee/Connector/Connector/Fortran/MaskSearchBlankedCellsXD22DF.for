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

C Search the cells which will be blanked "within" a distance of a XRay mask.
C Put the results in a cell nature field array.     
C A cell is blanked if the bounding box of cell (+/- delta) intersects 
C the XRay mask (+/- delta)
C and the center of cell is in the mask (+/- delta). Depth = 2 case.
C For 2D: the mesh must be nk=2 and in plane (x,y).

      SUBROUTINE k6searchblankedcellsxd22d( ni, nj, nk, 
     &     meshX, meshY, 
     &     xmin, 
     &     niray, njray,
     &     hiray,
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
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:MAX(ni-1,1)*MAX(nj-1,1)*MAX(nk-1,1)-1)
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp
      INTEGER_E ind, et, indray, indrayc, indraycp1
      INTEGER_E n
      REAL_E    dx1, dy1
      REAL_E    xp, yp
      REAL_E    xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1
      INTEGER_E iraymin, iraymax
      INTEGER_E irayminp, iraymaxp
      REAL_E ibmin, ibmax, ibc, ibcp1, ib
      REAL_E xrayc, alpha, xmax
      INTEGER_E npmin, npmax, nnp
      INTEGER_E cellN, nij, nijray, nic, nicnjc
      INTEGER_E tmp, ntemp
C
      tmp = 0
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      xmax = xmin + (niray-1) * hiray

      DO n = 0, np-1
         cellN = 0
         et = listOfInterpolatedPoints(n)
         j = et/nic
         i = et-j*nic
C     
#include "../../Connector/Fortran/MaskCell2DF.for"
C     
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4
         yp = yp1 + yp2 + yp3 + yp4    
         xp = xp * ONE_FOURTH
         yp = yp * ONE_FOURTH
         
         xmincell = xmincell-delta
         xmaxcell = xmaxcell+delta

# include "../../Connector/Fortran/MaskCenterInDelta2DF.for"
         
      ENDDO
      END




