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

      SUBROUTINE k6compmotioncentersex( nbi, im, jm, km, mat, r0, x0, 
     &                                  X, Y, Z, EXdirection, 
     &                                  EXindirectionMp, nfs, cellofint, 
     &                                  bary )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E nbi
      INTEGER_E im,jm,km
      INTEGER_E EXdirection(0:nbi-1)
      INTEGER_E EXindirectionMp(0:nbi-1)
      REAL_E     mat(3,3)       ! rotation matrix
      REAL_E     r0(3)          ! deplacement
      REAL_E     x0(3)          ! rotation center coordinates
      REAL_E X(0:im*jm*km-1)
      REAL_E Y(0:im*jm*km-1)
      REAL_E Z(0:im*jm*km-1)
      INTEGER_E nfs
      INTEGER_E cellofint(0:nfs-1,2)

C_OUT
      REAL_E bary(0:nbi-1,3)    

C_LOCAL
      INTEGER_E dir,i,dir2
      INTEGER_E ind,indc
      INTEGER_E imc,ijmc,ic,jc,kc,ijc
      INTEGER_E node1,node2,node3,node4
      INTEGER_E i1,i2,i3,i4
      INTEGER_E j1,j2,j3,j4
      INTEGER_E k1,k2,k3,k4
      INTEGER_E un, imjm

      INTEGER_E tabi1(0:5)
      INTEGER_E tabi2(0:5)
      INTEGER_E tabi3(0:5)
      INTEGER_E tabi4(0:5)
      INTEGER_E tabj1(0:5)
      INTEGER_E tabj2(0:5)
      INTEGER_E tabj3(0:5)
      INTEGER_E tabj4(0:5)
      INTEGER_E tabk1(0:5)
      INTEGER_E tabk2(0:5)
      INTEGER_E tabk3(0:5)
      INTEGER_E tabk4(0:5)

      REAL_E r11, r12, r13, r21, r22, r23, r31, r32, r33
      REAL_E dx, dy, dz
      REAL_E xp, yp, zp

      DATA tabi1 /1,0,0,0,0,0/
      DATA tabi2 /1,0,1,1,1,1/
      DATA tabi3 /1,0,0,0,0,0/
      DATA tabi4 /1,0,1,1,1,1/

      DATA tabj1 /0,0,1,0,0,0/
      DATA tabj2 /1,1,1,0,0,0/
      DATA tabj3 /0,0,1,0,1,1/
      DATA tabj4 /1,1,1,0,1,1/

      DATA tabk1 /0,0,0,0,1,0/
      DATA tabk2 /0,0,0,0,1,0/
      DATA tabk3 /1,1,1,1,1,0/
      DATA tabk4 /1,1,1,1,1,0/
C==============================================================================
      imc = im-1
      ijmc = imc*(jm-1)
      imjm = im*jm
      un = 1
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

      DO i = 0, nbi-1

C direction
         dir = EXdirection(i)
         IF (dir.GE.6) THEN
            dir = dir - 6
         ENDIF
C direction left/right
         dir2 = IAND(dir,un) + 1
C interface
         ind = EXindirectionMp(i)
C cellofint
         indc = cellofint(ind,dir2)
C indCell
         kc = indc/ijmc
         ijc = indc-kc*ijmc
         jc = ijc/imc
         ic = indc-jc*imc-kc*ijmc
         
C nodes
         i1 = ic+tabi1(dir)
         i2 = ic+tabi2(dir)
         i3 = ic+tabi3(dir)
         i4 = ic+tabi4(dir)
         j1 = jc+tabj1(dir)
         j2 = jc+tabj2(dir)
         j3 = jc+tabj3(dir)
         j4 = jc+tabj4(dir)
         k1 = kc+tabk1(dir)
         k2 = kc+tabk2(dir)
         k3 = kc+tabk3(dir)
         k4 = kc+tabk4(dir)

         node1 = i1 + j1*im + k1*imjm
         node2 = i2 + j2*im + k2*imjm
         node3 = i3 + j3*im + k3*imjm
         node4 = i4 + j4*im + k4*imjm

         xp=0.25D0*( X(node1) + X(node2) + X(node3) + X(node4) )
         yp=0.25D0*( Y(node1) + Y(node2) + Y(node3) + Y(node4) )
         zp=0.25D0*( Z(node1) + Z(node2) + Z(node3) + Z(node4) )

         bary(i,1) = dx + r11 * (xp-x0(1)) + r12 * (yp-x0(2)) 
     &        + r13 * (zp-x0(3))
                
         bary(i,2) = dy + r21 * (xp-x0(1)) + r22 * (yp-x0(2)) 
     &        + r23 * (zp-x0(3))
                
         bary(i,3) = dz + r31 * (xp-x0(1)) + r32 * (yp-x0(2))
     &        + r33 * (zp-x0(3))
      ENDDO
      
      RETURN
      
      END
C===== Connector/Fortran/CompMotionCentersEXF.for
