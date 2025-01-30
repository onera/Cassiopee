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

      SUBROUTINE  k6compmotioncenters(im, jm, km, nbnode, cellt, 
     &                                mat, r0, x0,
     &                                coord, bary )
 
      IMPLICIT NONE
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E  im             ! Mesh number of nodes %i
      INTEGER_E  jm             ! Mesh number of nodes %j
      INTEGER_E  km             ! Mesh number of nodes %k
      INTEGER_E  nbnode         ! Total number of nodes
      INTEGER_E  cellt          ! Total number of cells
      REAL_E     mat(3,3)       ! rotation matrix
      REAL_E     r0(3)          ! deplacement
      REAL_E     x0(3)          ! rotation center coordinates
      REAL_E     coord(0:nbnode-1,3) ! 3-D field of coord. of points
C_OUT
      REAL_E     bary(0:cellt-1,3) ! cell centers

C_LOCAL
      INTEGER_E lv, n1, n2, n3, n4, n5, n6, n7, n8
      INTEGER_E inci, incj, inck, incij, incjk, incik, incijk
      INTEGER_E imjm, i, j, k, im1, im1jm1
      REAL_E r11, r12, r13, r21, r22, r23, r31, r32, r33
      REAL_E dx, dy, dz
      REAL_E xp, yp, zp
C =============================================================================
C

      r11 = mat(1,1)
      r12 = mat(1,2)
      r13 = mat(1,3)
      r21 = mat(2,1)
      r22 = mat(2,2)
      r23 = mat(2,3)
      r31 = mat(3,1)
      r32 = mat(3,2)
      r33 = mat(3,3)
      dx = r0(1)
      dy = r0(2)
      dz = r0(3)

      incj = im 
      incjk = im+ im*jm
      inck = im*jm
      inci = 1
      incij = 1+im
      incijk = 1+im+im*jm
      incik = 1+im*jm
      imjm = im * jm
      im1 = im-1
      im1jm1 = im1*(jm-1)
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

                xp = ONE_EIGHT * (  coord(n1,1) + coord(n2,1)
     &               + coord(n3,1) + coord(n4,1)
     &               + coord(n5,1) + coord(n6,1)
     &               + coord(n7,1) + coord(n8,1))
                yp = ONE_EIGHT * (  coord(n1,2) + coord(n2,2)
     &               + coord(n3,2) + coord(n4,2)
     &               + coord(n5,2) + coord(n6,2)
     &               + coord(n7,2) + coord(n8,2))   
                zp =  ONE_EIGHT * (  coord(n1,3) + coord(n2,3)
     &               + coord(n3,3) + coord(n4,3)
     &               + coord(n5,3) + coord(n6,3)
     &               + coord(n7,3) + coord(n8,3)) 
             
                bary(lv,1) = dx + r11 * (xp-x0(1)) + r12 * (yp-x0(2)) 
     &               + r13 * (zp-x0(3))
                
                bary(lv,2) = dy + r21 * (xp-x0(1)) + r22 * (yp-x0(2)) 
     &               + r23 * (zp-x0(3))
                
                bary(lv,3) = dz + r31 * (xp-x0(1)) + r32 * (yp-x0(2))
     &               + r33 * (zp-x0(3))
             END DO
          ENDDO
       ENDDO
      END
C ======= Connector/Fortran/CompMotionCentersF.for
