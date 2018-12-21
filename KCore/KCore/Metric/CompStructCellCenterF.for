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
C  Compute barycenter of cells.
C  This version does not use minOfInt
C  cp de CompCenterCellF.for ou Geo/Grid/GeoCompCenterCellF.for du Kernel
C  ============================================================================
      SUBROUTINE  k6compstructcellcenter(im, jm, km, nbnode,cellt,     
     &                                   xt, yt, zt, bary )
 
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E im              ! Mesh number of nodes %i
      INTEGER_E jm              ! Mesh number of nodes %j
      INTEGER_E km              ! Mesh number of nodes %k
      INTEGER_E nbnode          ! Total number of nodes
      INTEGER_E cellt           ! Total number of cells
      REAL_E    xt(0:nbnode-1)  ! coord x of nodes
      REAL_E    yt(0:nbnode-1)  ! coord y of nodes
      REAL_E    zt(0:nbnode-1)  ! coord z of nodes
      
C_OUT
      REAL_E    bary(0:cellt-1,3) ! cell centers

C_LOCAL
      INTEGER_E lv, n1, n2, n3, n4, n5, n6, n7, n8
      INTEGER_E inci, incj, inck, incij, incjk, incik, incijk
      INTEGER_E imjm, i, j, k, im1, im1jm1

C =============================================================================
C
      imjm = im * jm
      im1 = im-1
      im1jm1 = im1*(jm-1)

      incj =  im 
      incjk = im + imjm
      inck =  imjm
      inci = 1
      incij = 1 + im
      incijk = 1+ im + imjm 
      incik = 1 + imjm

!$OMP PARALLEL PRIVATE(i,j,k,lv,n1,n2,n3,n4,n5,n6,n7,n8)
!$OMP DO
       DO k = 0, km-2
          DO j = 0, jm-2
             DO i = 0, im-2
                lv = i + j*im1 + k*im1jm1
                n1 = i + j*im + k*imjm           
                n2 = n1 + incj
                n3 = n1 + incjk
                n4 = n1 + inck
                n5 = n1 + inci
                n6 = n1 + incij
                n7 = n1 + incijk
                n8 = n1 + incik

                bary(lv,1) = ONE_EIGHT * ( 
     &               xt(n1) + xt(n2) + xt(n3) + xt(n4) +
     &               xt(n5) + xt(n6) + xt(n7) + xt(n8) )

                bary(lv,2) = ONE_EIGHT * (  
     &               yt(n1) + yt(n2) + yt(n3) + yt(n4) +
     &               yt(n5) + yt(n6) + yt(n7) + yt(n8) )  

                bary(lv,3) = ONE_EIGHT * (
     &               zt(n1) + zt(n2) + zt(n3) + zt(n4) +
     &               zt(n5) + zt(n6) + zt(n7) + zt(n8) )
             END DO
          ENDDO
       ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END
C ======================= Metric/CompStructCellCenterF.for ====================
