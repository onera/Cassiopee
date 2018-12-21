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
C  Search the minimum and maximum point of a set of cells
C  Cette version est differente de la routine dans Mask car les rotations
C  ne sont pas basées sur AgtTransfo : on rentre directement la matrice 3x3. 
C  ============================================================================
      SUBROUTINE k6compcartelembox( is1, is2, js1, js2, ks1, ks2,
     &                              i2, j2, k2,
     &                              im, jm, km, x, y, z,
     &                              rotMat,
     &                              du, dv, nbElts1, nbElts2,
     &                              dir, area, nodeMin, nodeMax,
     &                              szCartElt, cartEltMin, cartEltMax )
C
      IMPLICIT NONE
C
# include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E is1, is2, js1, js2, ks1, ks2 ! Window indices
      INTEGER_E i2, j2, k2
      INTEGER_E im, jm ,km             ! Mesh dimensions
      REAL_E    x(0:im*jm*km-1), y(0:im*jm*km-1), z(0:im*jm*km-1) ! Mesh
      REAL_E    rotMat(3,3)            ! Rotation matrix
      REAL_E    du, dv                 ! Discret. step of cart. box
      INTEGER_E dir                    ! Direction of mask
      INTEGER_E area                   ! -1/0 (both sided)/+1
      REAL_E    nodeMin(3)             ! Minimum of box.
      REAL_E    nodeMax(3)             ! Maximum of box.
      INTEGER_E nbElts1, nbElts2
      INTEGER_E szCartElt              ! Size of cartElt arrays
C_OUT
      REAL_E    cartEltMin(0:szCartElt-1) ! Minimum box of cart elem
      REAL_E    cartEltMax(0:szCartElt-1) ! Maximum box of cart elem
C_LOCAL
      INTEGER_E i,j,k, ip, jp ,kp, ind, ic, jc
      REAL_E    minpt(3), maxpt(3)
      REAL_E    xr, yr, zr, dui, dvi
      REAL_E    ptR1(3),
     &          ptR2(3),
     &          ptR3(3),
     &          ptR4(3),
     &          ptR5(3),
     &          ptR6(3),
     &          ptR7(3),
     &          ptR8(3)
      INTEGER_E icmin, icmax, jcmin, jcmax, ijc, imjm
C-----------------------------------------------------------------------------
      imjm = im*jm
      dui = ONE / du
      dvi = ONE / dv
C
      DO k = ks1, ks2
         DO j = js1, js2
            DO i = is1, is2
               minpt(1) = MAXFLOAT
               minpt(2) = MAXFLOAT
               minpt(3) = MAXFLOAT

               maxpt(1) = -MAXFLOAT
               maxpt(2) = -MAXFLOAT
               maxpt(3) = -MAXFLOAT

               ip       = MIN( i+1, i2 )
               jp       = MIN( j+1, j2 )
               kp       = MIN( k+1, k2 )

               ind      = (i-1) + (j-1)*im + (k-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR1(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR1(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR1(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (ip-1) + (j-1)*im + (k-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR2(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR2(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR2(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (ip-1) + (jp-1)*im + (k-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR3(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR3(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR3(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (i-1) + (jp-1)*im + (k-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR4(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR4(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR4(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (i-1) + (j-1)*im + (kp-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR5(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR5(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR5(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (ip-1) + (j-1)*im + (kp-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR6(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR6(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR6(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (i-1) + (jp-1)*im + (kp-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR7(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR7(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR7(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               ind      = (ip-1) + (jp-1)*im + (kp-1)*imjm
               xr       = x(ind)
               yr       = y(ind)
               zr       = z(ind)
               ptR8(1)  = rotMat(1,1)*xr + rotMat(1,2)*yr +
     &                    rotMat(1,3)*zr
               ptR8(2)  = rotMat(2,1)*xr + rotMat(2,2)*yr +
     &                    rotMat(2,3)*zr
               ptR8(3)  = rotMat(3,1)*xr + rotMat(3,2)*yr +
     &                    rotMat(3,3)*zr

               minpt(1) = MIN( minpt(1), ptR1(1) )
               minpt(1) = MIN( minpt(1), ptR2(1) )
               minpt(1) = MIN( minpt(1), ptR3(1) )
               minpt(1) = MIN( minpt(1), ptR4(1) )
               minpt(1) = MIN( minpt(1), ptR5(1) )
               minpt(1) = MIN( minpt(1), ptR6(1) )
               minpt(1) = MIN( minpt(1), ptR7(1) )
               minpt(1) = MIN( minpt(1), ptR8(1) )

               maxpt(1) = MAX( maxpt(1), ptR1(1) )
               maxpt(1) = MAX( maxpt(1), ptR2(1) )
               maxpt(1) = MAX( maxpt(1), ptR3(1) )
               maxpt(1) = MAX( maxpt(1), ptR4(1) )
               maxpt(1) = MAX( maxpt(1), ptR5(1) )
               maxpt(1) = MAX( maxpt(1), ptR6(1) )
               maxpt(1) = MAX( maxpt(1), ptR7(1) )
               maxpt(1) = MAX( maxpt(1), ptR8(1) )

               minpt(2) = MIN( minpt(2), ptR1(2) )
               minpt(2) = MIN( minpt(2), ptR2(2) )
               minpt(2) = MIN( minpt(2), ptR3(2) )
               minpt(2) = MIN( minpt(2), ptR4(2) )
               minpt(2) = MIN( minpt(2), ptR5(2) )
               minpt(2) = MIN( minpt(2), ptR6(2) )
               minpt(2) = MIN( minpt(2), ptR7(2) )
               minpt(2) = MIN( minpt(2), ptR8(2) )

               maxpt(2) = MAX( maxpt(2), ptR1(2) )
               maxpt(2) = MAX( maxpt(2), ptR2(2) )
               maxpt(2) = MAX( maxpt(2), ptR3(2) )
               maxpt(2) = MAX( maxpt(2), ptR4(2) )
               maxpt(2) = MAX( maxpt(2), ptR5(2) )
               maxpt(2) = MAX( maxpt(2), ptR6(2) )
               maxpt(2) = MAX( maxpt(2), ptR7(2) )
               maxpt(2) = MAX( maxpt(2), ptR8(2) )

               minpt(3) = MIN( minpt(3), ptR1(3) )
               minpt(3) = MIN( minpt(3), ptR2(3) )
               minpt(3) = MIN( minpt(3), ptR3(3) )
               minpt(3) = MIN( minpt(3), ptR4(3) )
               minpt(3) = MIN( minpt(3), ptR5(3) )
               minpt(3) = MIN( minpt(3), ptR6(3) )
               minpt(3) = MIN( minpt(3), ptR7(3) )
               minpt(3) = MIN( minpt(3), ptR8(3) )

               maxpt(3) = MAX( maxpt(3), ptR1(3) )
               maxpt(3) = MAX( maxpt(3), ptR2(3) )
               maxpt(3) = MAX( maxpt(3), ptR3(3) )
               maxpt(3) = MAX( maxpt(3), ptR4(3) )
               maxpt(3) = MAX( maxpt(3), ptR5(3) )
               maxpt(3) = MAX( maxpt(3), ptR6(3) )
               maxpt(3) = MAX( maxpt(3), ptR7(3) )
               maxpt(3) = MAX( maxpt(3), ptR8(3) )

               IF ( dir .EQ. 1 ) THEN
                  IF ( area .EQ. -1 ) minpt(1) = nodeMin(1)
                  IF ( area .EQ. +1 ) maxpt(1) = nodeMax(1)

                  IF ( ( minpt(2) .LE. nodeMax(2) ) .AND.
     &                 ( maxpt(2) .GE. nodeMin(2) ) .AND.
     &                 ( minpt(3) .LE. nodeMax(3) ) .AND.
     &                 ( maxpt(3) .GE. nodeMin(3) ) ) THEN
                     icmin = INT((minpt(2)-nodeMin(2))*dui) + 1
                     icmax = INT((maxpt(2)-nodeMin(2))*dui) + 1
                     icmax = MIN( icmax, nbElts1 )

                     jcmin = INT((minpt(3)-nodeMin(3))*dvi) + 1
                     jcmax = INT((maxpt(3)-nodeMin(3))*dvi) + 1
                     jcmax = MIN( jcmax, nbElts2 )

                     DO jc = jcmin, jcmax
                        DO ic = icmin, icmax
                           ijc = ic-1 + (jc-1)*nbElts1
                           cartEltMin(ijc) = MIN( cartEltMin(ijc), 
     &                                            minpt(1) )
                           cartEltMax(ijc) = MAX( cartEltMax(ijc), 
     &                                            maxpt(1) )
                        END DO
                     END DO
                  END IF
               ELSE IF  ( dir .EQ. 2 ) THEN
                  IF ( area .EQ. -1 ) minpt(2) = nodeMin(2)
                  IF ( area .EQ. +1 ) maxpt(2) = nodeMax(2)

                  IF ( ( minpt(3) .LE. nodeMax(3) ) .AND.
     &                 ( maxpt(3) .GE. nodeMin(3) ) .AND.
     &                 ( minpt(1) .LE. nodeMax(1) ) .AND.
     &                 ( maxpt(1) .GE. nodeMin(1) ) ) THEN
                     icmin = INT((minpt(3)-nodeMin(3))*dui) + 1
                     icmax = INT((maxpt(3)-nodeMin(3))*dui) + 1
                     icmax = MIN( icmax, nbElts1 )

                     jcmin = INT((minpt(1)-nodeMin(1))*dvi) + 1
                     jcmax = INT((maxpt(1)-nodeMin(1))*dvi) + 1
                     jcmax = MIN( jcmax, nbElts2 )

                     DO jc = jcmin, jcmax
                        DO ic = icmin, icmax
                           ijc = ic-1 + (jc-1)*nbElts1
                           cartEltMin(ijc) = MIN( cartEltMin(ijc), 
     &                                            minpt(2) )
                           cartEltMax(ijc) = MAX( cartEltMax(ijc), 
     &                                            maxpt(2) )
                        END DO
                     END DO
                  END IF
               ELSE IF ( dir .EQ. 3 ) THEN
                  IF ( area .EQ. -1 ) minpt(3) = nodeMin(3)
                  IF ( area .EQ. +1 ) maxpt(3) = nodeMax(3)

                  IF ( ( minpt(1) .LE. nodeMax(1) ) .AND.
     &                 ( maxpt(1) .GE. nodeMin(1) ) .AND.
     &                 ( minpt(2) .LE. nodeMax(2) ) .AND.
     &                 ( maxpt(2) .GE. nodeMin(2) ) ) THEN
                     icmin = INT((minpt(1)-nodeMin(1))*dui) + 1
                     icmax = INT((maxpt(1)-nodeMin(1))*dui) + 1
                     icmax = MIN( icmax, nbElts1 )

                     jcmin = INT((minpt(2)-nodeMin(2))*dvi) + 1
                     jcmax = INT((maxpt(2)-nodeMin(2))*dvi) + 1
                     jcmax = MIN( jcmax, nbElts2 )

                     DO jc = jcmin, jcmax
                        DO ic = icmin, icmax
                           ijc = ic-1 + (jc-1)*nbElts1
                           cartEltMin(ijc) = MIN( cartEltMin(ijc), 
     &                                            minpt(3) )
                           cartEltMax(ijc) = MAX( cartEltMax(ijc), 
     &                                            maxpt(3) )
                        END DO
                     END DO
                  END IF
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
C ========================= CompGeom/CompCEBoxF.for ===========================
