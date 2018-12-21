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

C ============================================================================
C Compute the minimum distance between 2 blocs and returns
C the corresponding cell indices
C  ===========================================================================
      SUBROUTINE k6compmindist(ni1, nj1, x1, y1, z1, ni2, nj2, 
     &                         x2, y2, z2,
     &                         ind1s, ind2s, dmin)
C
      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni1, nj1
      INTEGER_E ni2, nj2
      REAL_E  x1(0:ni1*nj1-1)
      REAL_E  y1(0:ni1*nj1-1)
      REAL_E  z1(0:ni1*nj1-1)
      REAL_E  x2(0:ni2*nj2-1)
      REAL_E  y2(0:ni2*nj2-1)
      REAL_E  z2(0:ni2*nj2-1)

C_OUT
      INTEGER_E ind1s, ind2s
      REAL_E dmin
C_LOCAL
      INTEGER_E i1, j1, i2, j2, m1, m2
      INTEGER_E ind1, ind2
      REAL_E dist
      REAL_E x10, y10, z10 
      REAL_E x20, y20, z20
      INTEGER_E ind11, ind12, ind13, ind21, ind22, ind23

      ind11 = ni1-1
      ind12 = (nj1-1)*ni1
      ind13 = ind11 + (nj1-1) * ni1
      ind21 = ni2-1
      ind22 = (nj2-1)*ni2
      ind23 = ind21 + (nj2-1) * ni2

C==============================================================================
      dmin = MAXFLOAT
      
C     distance  entre les points sauf les coins
      DO j1 = 0, nj1-1
         DO i1 = 0, ni1-1
            ind1 = i1 + j1 * ni1 

            IF ( ind1 .ne. 0 .and. ind1 .ne. ind11 .and. 
     &           ind1 .ne. ind12 .and. ind1 .ne. ind13 ) THEN 
               x10 = x1(ind1)
               y10 = y1(ind1)
               z10 = z1(ind1)
            
               DO j2 = 0, nj2-1
                  DO i2 = 0, ni2-1
                     ind2 = i2 + j2 * ni2 
                     IF (ind2 .ne. 0 .and. ind2 .ne. ind21 .and. 
     &                   ind2 .ne. ind22 .and. ind2 .ne. ind23) THEN 
                        x20 = x2(ind2)
                        y20 = y2(ind2)
                        z20 = z2(ind2)  
                        
                        dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &                       + (z20-z10)*(z20-z10)
                        
                        IF (dist .LT. dmin) THEN
                           dmin = dist
                           ind1s = ind1
                           ind2s = ind2
                           
C     Cas ou distance nulle: arret  
                           IF (dist .lt. 1.e-12) THEN
                              RETURN
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

C Distance entre les points de 1 avec les coins de 2
      DO j1 = 0, nj1-1
         DO i1 = 0, ni1-1
            ind1 = i1 + j1 * ni1 

            IF (ind1 .ne. 0 .and. ind1 .ne. ind11 .and. 
     &          ind1 .ne. ind12 .and. ind1 .ne. ind13) THEN 
               x10 = x1(ind1)
               y10 = y1(ind1)
               z10 = z1(ind1)
               
               x20 = x2(0)
               y20 = y2(0)
               z20 = z2(0)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)

               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind1
                  ind2s = 0
                  
C     Cas ou distance nulle : arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF
               
               x20 = x2(ind21)
               y20 = y2(ind21)
               z20 = z2(ind21)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)
               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind1
                  ind2s = ind21
                  
C     Cas ou distance nulle : arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF

               x20 = x2(ind22)
               y20 = y2(ind22)
               z20 = z2(ind22)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)
               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind1
                  ind2s = ind22
                  
C     Cas ou distance nulle: arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF

               x20 = x2(ind23)
               y20 = y2(ind23)
               z20 = z2(ind23)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)
               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind1
                  ind2s = ind23
                  
C     Cas ou distance nulle: arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF

            ENDIF
         ENDDO
      ENDDO


C Distance entre les points de 2 avec les coins de 1
      DO j2 = 0, nj2-1
         DO i2 = 0, ni2-1
            ind2 = i2 + j2 * ni2 

            IF (ind2 .ne. 0 .and. ind2 .ne. ind21 .and. 
     &          ind2 .ne. ind22 .and. ind2 .ne. ind23) THEN 
               x20 = x2(ind2)
               y20 = y2(ind2)
               z20 = z2(ind2)
               
               x10 = x1(0)
               y10 = y1(0)
               z10 = z1(0)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)
               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = 0
                  ind2s = ind2
                  
C     Cas ou distance nulle : arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF
               
               x10 = x1(ind11)
               y10 = y1(ind11)
               z10 = z1(ind11)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)

               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind11
                  ind2s = ind2 
                  
C     Cas ou distance nulle: arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF

               x10 = x1(ind12)
               y10 = y1(ind12)
               z10 = z1(ind12)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)

               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind12
                  ind2s = ind2 
                  
C     Cas ou distance nulle: arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF

               x10 = x1(ind13)
               y10 = y1(ind13)
               z10 = z1(ind13)  
               
               dist = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10)  
     &              + (z20-z10)*(z20-z10)
               IF (dist .LT. dmin) THEN
                  dmin = MIN(dist, dmin)
                  ind1s = ind13
                  ind2s = ind2
                  
C     Cas ou distance nulle: arret  
                  IF (dist .lt. 1.e-12) THEN
                     RETURN
                  ENDIF
               ENDIF

            ENDIF
         ENDDO
      ENDDO
            
      RETURN
      END
