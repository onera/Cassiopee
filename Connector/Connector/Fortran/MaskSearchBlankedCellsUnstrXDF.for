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

C  ============================================================================
C searchblankedcellsunstrxd
C Search the cells which will be blanked "within" a distance of a XRay mask.
C Put the results in a cell nature field array.     
C A cell is blanked if the bounding box of cell (+/- delta) intersects 
C the XRay mask (+/- delta)
C and the center of cell is in the mask (+/- delta). Depth = 2 case.
C For 2D: the mesh must be tri or quad and in plane (x,y).
C ============================================================================

      SUBROUTINE k6searchblankedcellstrixd( 
     &     npts, xt, yt, nelts, cn1, cn2, cn3,
     &     xmin, niray, njray, hiray, indir, nz, z, 
     &     listOfInterpolatedPoints, np, delta, isnot,
     &     cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

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
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
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
      INTEGER_E cellN, nijray
      INTEGER_E tmp, ntemp
C=============================================================================
      
      nijray = niray * njray
      isMasked = 0

      IF (isnot .EQ. 0) THEN
 
         DO n = 0, np-1
            cellN = 0
            et = listOfInterpolatedPoints(n)       
C     
#include "../../Connector/Fortran/MaskTriCellF.for"
C     
            xmincell = xmincell-delta
            xmaxcell = xmaxcell+delta

# include "../../Connector/Fortran/MaskCellIntersect2DF.for"
         
         ENDDO
      ELSE
         DO n = 0, np-1
            cellN = 0
            et = listOfInterpolatedPoints(n)       
C     
#include "../../Connector/Fortran/MaskTriCellF.for"
C     
            xmincell = xmincell-delta
            xmaxcell = xmaxcell+delta

# include "../../Connector/Fortran/MaskCellIntersectNot2DF.for"
         
         ENDDO         
      ENDIF
      END
C ============================================================================
C Cas QUAD
C ============================================================================
      SUBROUTINE k6searchblankedcellsquadxd( 
     &     npts, xt, yt, nelts, cn1, cn2, cn3, cn4,
     &     xmin, niray, njray, hiray, indir, nz, z, 
     &     listOfInterpolatedPoints, np, delta, isnot,
     &     cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           ! nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite elts/noeuds
      INTEGER_E cn2(0:nelts-1)  
      INTEGER_E cn3(0:nelts-1)  
      INTEGER_E cn4(0:nelts-1)
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
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
      INTEGER_E cellN, nijray
      INTEGER_E tmp, ntemp
C=============================================================================
      nijray = niray * njray
      isMasked = 0

      IF (isnot .EQ. 0) THEN
 
         DO n = 0, np-1
            cellN = 0
            et = listOfInterpolatedPoints(n)       
C     
#include "../../Connector/Fortran/MaskQuadCellF.for"
C     
            xmincell = xmincell-delta
            xmaxcell = xmaxcell+delta

# include "../../Connector/Fortran/MaskCellIntersect2DF.for"
         
         ENDDO
      ELSE
         DO n = 0, np-1
            cellN = 0
            et = listOfInterpolatedPoints(n)       
C     
#include "../../Connector/Fortran/MaskTriCellF.for"
C     
            xmincell = xmincell-delta
            xmaxcell = xmaxcell+delta

# include "../../Connector/Fortran/MaskCellIntersectNot2DF.for"
         
         ENDDO         
      ENDIF
      END
 
C ============================================================================
C Cas QUAD
C ============================================================================
      SUBROUTINE k6searchblankedcellstetraxd( 
     &     npts, xt, yt, zt, nelts, cn1, cn2, cn3, cn4,
     &     xmin, ymin, niray, njray, hiray, hjray, indir, nz, z, 
     &     listOfInterpolatedPoints, np, delta, isnot,
     &     cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      REAL_E    zt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite elts/noeuds
      INTEGER_E cn2(0:nelts-1)  
      INTEGER_E cn3(0:nelts-1)  
      INTEGER_E cn4(0:nelts-1)
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E ind, et, indray
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E n, npmin, npmax
      INTEGER_E cellN, nij, nijray, nic, nicnjc

C=============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN 
         DO n = 0, np-1
            cellN = 0
            et = listOfInterpolatedPoints(n)       
C     
#include "../../Connector/Fortran/MaskTetraCellF.for"
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
            et = listOfInterpolatedPoints(n)       
C     
#include "../../Connector/Fortran/MaskTetraCellF.for"
C     
            xmincell = xmincell-delta
            ymincell = ymincell-delta
            xmaxcell = xmaxcell+delta
            ymaxcell = ymaxcell+delta

#include "../../Connector/Fortran/MaskCellIntersectNotF.for"
         ENDDO         
      ENDIF
      END
C ============================================================================
C Cas QUAD
C ============================================================================
      SUBROUTINE k6searchblankedcellspentaxd( 
     &     npts, xt, yt, zt, nelts, cn1, cn2, cn3, cn4, cn5, cn6,
     &     xmin, ymin, niray, njray, hiray, hjray, indir, nz, z, 
     &     listOfInterpolatedPoints, np, delta, isnot,
     &     cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      REAL_E    zt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite elts/noeuds
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
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E ind, et, indray
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E n, npmin, npmax
      INTEGER_E cellN, nij, nijray, nic, nicnjc

C=============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN 
      DO n = 0, np-1
         cellN = 0
         et = listOfInterpolatedPoints(n)
C        
#include "../../Connector/Fortran/MaskPentaCellF.for"
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
         et = listOfInterpolatedPoints(n)
C        
#include "../../Connector/Fortran/MaskPentaCellF.for"
C
         xmincell = xmincell-delta
         ymincell = ymincell-delta
         xmaxcell = xmaxcell+delta
         ymaxcell = ymaxcell+delta
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"

      ENDDO
      ENDIF
      END
C ============================================================================
C Cas QUAD
C ============================================================================
      SUBROUTINE k6searchblankedcellshexaxd( 
     &     npts, xt, yt, zt, nelts, 
     &     cn1, cn2, cn3, cn4, cn5, cn6, cn7, cn8,
     &     xmin, ymin, niray, njray, hiray, hjray, indir, nz, z, 
     &     listOfInterpolatedPoints, np, delta, isnot,
     &     cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      REAL_E    zt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite elts/noeuds
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
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot

C_OUT
      INTEGER_E cellNatureField(0:nelts-1) 
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

C=============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN 
      DO n = 0, np-1
         cellN = 0
         et = listOfInterpolatedPoints(n)
C        
#include "../../Connector/Fortran/MaskHexaCellF.for"
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
         et = listOfInterpolatedPoints(n)
C        
#include "../../Connector/Fortran/MaskHexaCellF.for"
C
         xmincell = xmincell-delta
         ymincell = ymincell-delta
         xmaxcell = xmaxcell+delta
         ymaxcell = ymaxcell+delta
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"

      ENDDO
      ENDIF
      END
