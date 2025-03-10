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

C  ============================================================================
C Calcul des centres des interfaces pour une grille structuree definie 
C en noeuds
C  ============================================================================
      SUBROUTINE  k6compcenterinterface(im, jm, km, nbnode,
     &     intt, inti, intij, x, y, z, cix, ciy, ciz)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E im            ! Mesh number of nodes %i
      INTEGER_E jm            ! Mesh number of nodes %j
      INTEGER_E km            ! Mesh number of nodes %k
      INTEGER_E nbnode        ! Total number of nodes
      INTEGER_E intt          ! Total number of interfaces
      INTEGER_E inti          ! Total number of interfaces %i direction
      INTEGER_E intij         ! Total number of interfaces %i direction and j direction
      REAL_E    x(0:nbnode-1), y(0:nbnode-1), z(0:nbnode-1) ! coord of points
C_OUT
      REAL_E    cix(0:intt-1), ciy(0:intt-1), ciz(0:intt-1) ! interf centers

C_LOCAL
      REAL_E xc, yc, zc
      INTEGER_E li,lj,lk,l1, l2, l3, l4, i, j, k, imjm
      INTEGER_E im1, jm1, im1jm1
      INTEGER_E incnodeij,incnodeik,incnodeijk
      INTEGER_E incnodeji,incnodejk,incnodejik
      INTEGER_E incnodeki,incnodekj,incnodekij

C =============================================================================

      imjm = im*jm
      im1 = im-1
      jm1 = jm-1
      im1jm1 = im1*jm1

C PARAMETRES D'INCREMENTATION POUR LES INTERFACES "I"
      incnodeij  = im           ! incnode(0,1,0,im,jm,km)
      incnodeik  = imjm         ! incnode(0,0,1,im,jm,km)
      incnodeijk = im + imjm    ! incnode(0,1,1,im,jm,km)

C PARAMETRES D'INCREMENTATION POUR LES INTERFACES "J"
      incnodejk  = imjm         ! incnode(0,0,1,im,jm,km)
      incnodeji  = 1            ! incnode(1,0,0,im,jm,km)
      incnodejik = 1 + imjm     ! incnode(1,0,1,im,jm,km)

C PARAMETRES D'INCREMENTATION POUR LES INTERFACES "K"
      incnodeki  = 1            ! incnode(1,0,0,im,jm,km)
      incnodekj  = im           ! incnode(0,1,0,im,jm,km)
      incnodekij = 1 + im       ! incnode(1,1,0,im,jm,km)

!$OMP PARALLEL PRIVATE(i,j,k,li,lj,lk,l1,l2,l3,l4,xc,yc,zc)
!$OMP DO
      DO k = 0, km-2
         DO j = 0, jm-2
            DO  i = 0, im1            
               li = i + j*im + k*im*jm1
               l1 = i + j*im + k*imjm
               l2 = l1 + incnodeij
               l3 = l1 + incnodeijk
               l4 = l1 + incnodeik

               xc = ONE_FOURTH * (x(l1) + x(l2) + x(l3) + x(l4))
               yc = ONE_FOURTH * (y(l1) + y(l2) + y(l3) + y(l4))
               zc = ONE_FOURTH * (z(l1) + z(l2) + z(l3) + z(l4))

               cix(li) = xc
               ciy(li) = yc
               ciz(li) = zc
            END DO
         ENDDO
      ENDDO
!$OMP END DO NOWAIT

!$OMP DO
      DO k = 0, km-2
         DO j = 0, jm1
            DO i = 0, im-2
            
               lj = i + j*im1 + k*im1*jm + inti
               l1 = i + j*im + k*imjm
               l2 = l1 + incnodejk 
               l3 = l1 + incnodejik
               l4 = l1 + incnodeji

               xc = ONE_FOURTH * (x(l1) + x(l2) + x(l3) + x(l4))
               yc = ONE_FOURTH * (y(l1) + y(l2) + y(l3) + y(l4))
               zc = ONE_FOURTH * (z(l1) + z(l2) + z(l3) + z(l4))

               cix(lj) = xc
               ciy(lj) = yc
               ciz(lj) = zc
            END DO
         ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO k = 0, km-1
         DO j = 0, jm-2
            DO i = 0, im-2
            
               lk = i + j*im1 + k*im1jm1 + intij
               l1 = i + j*im + k*imjm
               l2 = l1  + incnodeki
               l3 = l1  + incnodekij
               l4 = l1  + incnodekj
               
               xc = ONE_FOURTH * (x(l1) + x(l2) + x(l3) + x(l4))
               yc = ONE_FOURTH * (y(l1) + y(l2) + y(l3) + y(l4))
               zc = ONE_FOURTH * (z(l1) + z(l2) + z(l3) + z(l4))

               cix(lk) = xc
               ciy(lk) = yc
               ciz(lk) = zc
            END DO
         ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      END
C ===================== KCore/Metric/CompCenterInterfaceF.for =================
