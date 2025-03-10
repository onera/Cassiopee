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
C      Put the results in a cell nature field array.
C      A cell is blanked if the bounding box of cell (+/- delta) intersects 
C      the XRay mask (+/- delta).
C      For 2D: the mesh must be nk=2 and in plane (x,y).

      SUBROUTINE k6searchblankednodesxd2d(npts, meshX, meshY, 
     &     xmin, niray, njray,
     &     hiray, indir, nz, z, 
     &     listOfInterpolatedPoints, np, diri,
     &     delta, isnot,
     &     cellNatureField,
     &     isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E npts            ! Dimension of the domain
      REAL_E    meshX(0:npts-1) ! Abcisses of mesh points
      REAL_E    meshY(0:npts-1) ! Ordinates of mesh points
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      INTEGER_E diri(0:np-1)
      REAL_E    delta           ! Wall distance
      INTEGER_E isnot           ! inverse masque
C_OUT
      INTEGER_E cellNatureField(0:npts-1) ! Give the nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E ip, jp
      INTEGER_E ind, et, indray
      INTEGER_E n, l
      REAL_E    dx1, dy1
      REAL_E    xp1, xp2, xp3, xp4, yp1, yp2, yp3, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E iray, ibeg, iend
      INTEGER_E iraymin, iraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nijray, nic, nicnjc
      REAL_E xp, yp, xmax, ymax, zmax, z1, z2, xrayc,alpha, zi
      INTEGER_E ibmin, ibmax, irayc, indrayc, nnp, ibc, pos
      INTEGER_E iraycp1, indraycp1, ibcp1
C==============================================================================
      nijray = niray * njray
      isMasked = 0
      xmax = xmin + (niray-1) * hiray

      DO n = 0, np-1
         cellN = 0
         et = listOfInterpolatedPoints(n)
         xmincell = meshX(et)+delta*min(diri(n),0)
         ymincell = meshY(et)
         xmaxcell = xmincell+delta
         ymaxcell = ymincell   
         xp = meshX(et)+delta*diri(n)
         yp = ymincell

#include "../../Connector/Fortran/MaskNodeIn2DF.for"          

      ENDDO

      IF (isnot .NE. 0) THEN   
         DO n = 0, np-1
            et = listOfInterpolatedPoints(n)
            cellNatureField(et) = -cellNatureField(et)
         ENDDO
      ENDIF

      END
     




