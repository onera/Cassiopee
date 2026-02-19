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

C Search the nodes which will be blanked "within" a distance of XRay mask.     
C Put the results in a cell nature field array.
C A cell is blanked if the bounding box of cell (+/- delta) intersects 
C the XRay mask (+/- delta).

      SUBROUTINE k6searchblankednodesxd(npts, meshX, meshY, meshZ,
     &     xmin, ymin, niray, njray, hiray, hjray,
     &     indir, nz, z, 
     &     listOfInterpolatedPoints, np, diri, dirj,
     &     delta, isnot,
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts     ! Dimension of the domain
      REAL_E    meshX(0:npts-1) ! Abcisses of mesh points
      REAL_E    meshY(0:npts-1) ! Ordinates of mesh points
      REAL_E    meshZ(0:npts-1) ! Height of mesh points
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      INTEGER_E diri(0:np-1)
      INTEGER_E dirj(0:np-1)
      REAL_E    delta           ! Wall distance
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:npts-1) ! cellN
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E ind, et, pos, indray
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
      INTEGER_E cellN,  nijray
      INTEGER_E irayc, iraycp1, iraycp2, iraycp3
      INTEGER_E jrayc, jraycp1, jraycp2, jraycp3
      REAL_E ibmin, ibmax, ibc1, ibc2, ibc3, ibc4
      REAL_E    xp, yp, zp
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      REAL_E z1, z2, z3, z4, zi
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax

C==============================================================================
      nijray = niray * njray
      isMasked = 0
      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 

      DO n = 0, np-1
         cellN = 0
         et = listOfInterpolatedPoints(n)

         xmincell = meshX(et)+delta*min(diri(n),0)
         ymincell = meshY(et)+delta*min(dirj(n),0)
         zmincell = meshZ(et)
         xmaxcell = xmincell+delta
         ymaxcell = ymincell+delta  
         zmaxcell = zmincell
         xp = meshX(et)+delta*diri(n)
         yp = meshY(et)+delta*dirj(n)
         zp = zmincell

#include "../../Connector/Fortran/MaskNodeInF.for"

      ENDDO

      IF (isnot .NE. 0) THEN   
         DO n = 0, np-1
            et = listOfInterpolatedPoints(n)
            cellNatureField(et) = -cellNatureField(et)
         ENDDO
      ENDIF

      END


