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

C Cas XRay 2D + node in

      SUBROUTINE k6searchblankednodesx2d(      
     &     npts, meshX, meshY,
     &     xmin, niray, njray, hiray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts
      REAL_E    meshX(0:npts-1) ! Abcisses of mesh points
      REAL_E    meshY(0:npts-1) ! Ordinates of mesh points
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:npts-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E et, indray
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E iray, ibeg, iend
      INTEGER_E iraymin, iraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nijray
      REAL_E xp, yp, xmax, ymax, zmax, z1, z2, xrayc,alpha, zi
      INTEGER_E ibmin, ibmax, irayc, indrayc, nnp, ibc, pos
      INTEGER_E iraycp1, indraycp1, ibcp1
      
C==============================================================================
      nijray = niray * njray
      isMasked = 0
      xmax = xmin + (niray-1) * hiray

      DO et = 0, npts-1
         cellN = 0
         xmincell = meshX(et)
         ymincell = meshY(et)
         xmaxcell = meshX(et)
         ymaxcell = meshY(et)   
         xp = xmincell
         yp = ymincell
         
#include "../../Connector/Fortran/MaskNodeIn2DF.for"          
      ENDDO

      IF (isnot .NE. 0) THEN   
         DO et = 0, npts-1
            cellNatureField(et) = -cellNatureField(et)
         ENDDO
      ENDIF

      END
C==============================================================================
