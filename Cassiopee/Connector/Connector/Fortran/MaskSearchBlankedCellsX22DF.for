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

C center_in + depth = 2 cas 2D

      SUBROUTINE k6searchblankedcellsx22d(ni, nj, nk, 
     &     meshX, meshY,
     &     xmin, 
     &     niray, njray,
     &     hiray,
     &     indir, nz, z, 
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
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:MAX(ni-1,1)*MAX(nj-1,1)*MAX(nk-1,1)-1) 
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
      INTEGER_E cellN, nij, nijray, nic, nicnjc, pos
      REAL_E z1, z2, zi
C=============================================================================
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0
      xmax = xmin + (niray-1) * hiray
      pos = 0
      DO j = 0, nj-2
         DO i = 0, nic-1
            cellN = 0
            et = i + j * nic
c
#include "../../Connector/Fortran/MaskCell2DF.for"
C
C       Barycenter of cell
C              
            xp = xp1 + xp2 + xp3 + xp4
            yp = yp1 + yp2 + yp3 + yp4    
            xp = xp * ONE_FOURTH
            yp = yp * ONE_FOURTH
           
# include "../../Connector/Fortran/MaskCenterIn2DF.for"

         ENDDO
      ENDDO
      END




