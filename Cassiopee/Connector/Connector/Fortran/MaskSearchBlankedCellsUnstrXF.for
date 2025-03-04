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

C cell intersect + depth=1 et 2 en NS: TRI, QUAD, TETRA, HEXA, PENTA
C ============================================================================
C TRI CASE FOR CELL INTERSECT ALGORITHM 
C ============================================================================
      SUBROUTINE k6searchblankedcellstrix( 
     &     npts, xt, yt, nelts, cn1, cn2, cn3,
     &     xmin, niray, njray, hiray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

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
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:nelts-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E ind, indray, et, l
      REAL_E    dx1, dy1
      REAL_E    xp1, xp2, xp3, xp4
      REAL_E    yp1, yp2, yp3, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E iray, ibeg, iend
      INTEGER_E iraymin, iraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nij, nijray
C==============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN
      DO et = 0, nelts-1
         cellN = 0
C
#include "../../Connector/Fortran/MaskTriCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersect2DF.for"
C
      ENDDO
      ELSE
         DO et = 0, nelts-1
            cellN = 0
C     
#include "../../Connector/Fortran/MaskTriCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectNot2DF.for"
C
      ENDDO
      ENDIF
      END

C=============================================================================
C QUAD CASE FOR CELL INTERSECT ALGORITHM 
C=============================================================================
      SUBROUTINE k6searchblankedcellsquadx( 
     &     npts, xt, yt, nelts, cn1, cn2, cn3, cn4,
     &     xmin, niray, njray, hiray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb of pts in zone
      REAL_E    xt(0:npts-1)    ! Abcisses of mesh points
      REAL_E    yt(0:npts-1)    ! Ordinates of mesh points
      INTEGER_E nelts           !nb of elements in unstructured zone
      INTEGER_E cn1(0:nelts-1)  ! connectivite  quad : 1ere colonne
      INTEGER_E cn2(0:nelts-1)  ! connectivite  quad : 2eme colonne
      INTEGER_E cn3(0:nelts-1)  ! connectivite  quad: 3eme colonne
      INTEGER_E cn4(0:nelts-1)  ! connectivite  quad : 3eme colonne
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:nelts-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E ind, indray, et
      REAL_E    dx1, dy1
      REAL_E    xp1, xp2, xp3, xp4
      REAL_E    yp1, yp2, yp3, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E iray, ibeg, iend, l 
      INTEGER_E iraymin, iraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nij, nijray
C==============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN
      DO et = 0, nelts-1
         cellN = 0
C
#include "../../Connector/Fortran/MaskQuadCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersect2DF.for"
C
      ENDDO
      ELSE
         DO et = 0, nelts-1
            cellN = 0
C     
#include "../../Connector/Fortran/MaskQuadCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectNot2DF.for"
C
      ENDDO
      ENDIF
      END
C=============================================================================
C CAS TETRA + CELL INTERSECT 
C=============================================================================
      SUBROUTINE k6searchblankedcellstetrax( 
     &     npts, xt, yt, zt, nelts, cn1, cn2, cn3, cn4,
     &     xmin, ymin, niray, njray, hiray, hjray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

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
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:nelts-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E et, indray, ind
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xp7, yp7, zp7, xp8, yp8, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nijray
C==============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN
         DO et = 0, nelts-1
            cellN = 0
C
#include "../../Connector/Fortran/MaskTetraCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectF.for"
C
         ENDDO
      ELSE
         DO et = 0, nelts-1
            cellN = 0
C     
#include "../../Connector/Fortran/MaskTetraCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"
C
         ENDDO
      ENDIF
      END
C=============================================================================
C CAS PENTA + CELL INTERSECT 
C=============================================================================
      SUBROUTINE k6searchblankedcellspentax( 
     &     npts, xt, yt, zt, nelts, cn1, cn2, cn3, cn4, cn5, cn6,
     &     xmin, ymin, niray, njray, hiray, hjray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

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
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:nelts-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E et, indray, ind
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xp7, yp7, zp7, xp8, yp8, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nijray
C==============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN
         DO et = 0, nelts-1
            cellN = 0
C
#include "../../Connector/Fortran/MaskPentaCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectF.for"
C
         ENDDO
      ELSE
         DO et = 0, nelts-1
            cellN = 0
C     
#include "../../Connector/Fortran/MaskPentaCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"
C
         ENDDO
      ENDIF
      END
C=============================================================================
C CAS HEXA + CELL INTERSECT 
C=============================================================================
      SUBROUTINE k6searchblankedcellshexax( 
     &     npts, xt, yt, zt, nelts, 
     &     cn1, cn2, cn3, cn4, cn5, cn6, cn7, cn8,
     &     xmin, ymin, niray, njray, hiray, hjray, indir,
     &     nz, z, isnot, cellNatureField, isMasked )

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
      INTEGER_E cn7(0:nelts-1) 
      INTEGER_E cn8(0:nelts-1) 
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points
      INTEGER_E isnot           ! inverse mask
C_OUT
      INTEGER_E cellNatureField(0:nelts-1) ! nature of the cells ( masked or not )
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       
      INTEGER_E ip, jp, kp
      INTEGER_E et, indray, ind
      REAL_E    dx1, dy1, dz1
      REAL_E    xp1, yp1, zp1, xp2, yp2, zp2
      REAL_E    xp3, yp3, zp3, xp4, yp4, zp4
      REAL_E    xp5, yp5, zp5, xp6, yp6, zp6
      REAL_E    xp7, yp7, zp7, xp8, yp8, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E iray, jray, ibeg, iend
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E npmin, npmax
      INTEGER_E cellN, nijray
C==============================================================================
      nijray = niray * njray
      isMasked = 0

      IF ( isnot .EQ. 0 ) THEN
         DO et = 0, nelts-1
            cellN = 0
C
#include "../../Connector/Fortran/MaskHexaCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectF.for"
C
         ENDDO
      ELSE
         DO et = 0, nelts-1
            cellN = 0
C     
#include "../../Connector/Fortran/MaskHexaCellF.for"
C
#include "../../Connector/Fortran/MaskCellIntersectNotF.for"
C
         ENDDO
      ENDIF
      END
