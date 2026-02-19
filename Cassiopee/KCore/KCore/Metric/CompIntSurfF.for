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

C  ============================================================================
C  Surface Vector and normal to the interfaces for structured grids
C  ============================================================================
      SUBROUTINE  k6compintsurf(im, jm, km, nbnode,
     &                          intt, inti, intij, x, y, z,
     &                          surfx, surfy, surfz, snorm)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E  im             ! Mesh number of nodes %i
      INTEGER_E  jm             ! Mesh number of nodes %j
      INTEGER_E  km             ! Mesh number of nodes %k
      INTEGER_E  nbnode         ! Total number of nodes
      INTEGER_E  intt           ! Total number of interfaces
      INTEGER_E  inti           ! Total number of interfaces %i direction
      INTEGER_E  intij          ! Total number of interfaces %i direction 
                                ! and j direction
      REAL_E x(0:nbnode-1), y(0:nbnode-1), z(0:nbnode-1) ! coordinates of points

C_OUT
      REAL_E surfx(0:intt-1), surfy(0:intt-1), surfz(0:intt-1) ! surface Vector
      REAL_E snorm(0:intt-1)  !  norm of surf

C_LOCAL
      REAL_E sgn            !  +1 if the mesh is direct, -1 if it is not
      REAL_E d13x, d13y, d13z, d24x, d24y, d24z
      REAL_E nx, ny, nz
      INTEGER_E li,lj,lk,l1, l2, l3, l4
      INTEGER_E incnodeij, incnodeik, incnodeijk
      INTEGER_E incnodeji, incnodejk, incnodejik
      INTEGER_E incnodeki, incnodekj, incnodekij
      INTEGER_E i, j, k, im1, jm1, im1jm1, imjm

      INTRINSIC SQRT
C =============================================================================
C
      sgn = 1.
C
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

!$OMP PARALLEL PRIVATE(i,j,k,li,lj,lk,l1,l2,l3,l4,d13x,d13y,d13z)
!$OMP& PRIVATE(d24x,d24y,d24z,nx,ny,nz)

!$OMP DO
      DO k = 0, km-2
         DO j = 0, jm-2
            DO  i = 0, im-1
            
               li = i + j*im + k*im*jm1
               l1 = i + j*im + k*imjm
               l2 = l1 + incnodeij
               l3 = l1 + incnodeijk
               l4 = l1 + incnodeik

               d13x = x(l3) - x(l1)
               d13y = y(l3) - y(l1)
               d13z = z(l3) - z(l1)
               
               d24x = x(l4) - x(l2)
               d24y = y(l4) - y(l2)
               d24z = z(l4) - z(l2)
               
               nx = ONE_HALF * (d13y*d24z - d13z*d24y)
               ny = ONE_HALF * (d13z*d24x - d13x*d24z)
               nz = ONE_HALF * (d13x*d24y - d13y*d24x)
         
               surfx(li) = sgn*nx
               surfy(li) = sgn*ny
               surfz(li) = sgn*nz    
               snorm(li) = SQRT(nx*nx + ny*ny + nz*nz)
            END DO
         ENDDO
      ENDDO
!$OMP END DO NOWAIT

!$OMP DO
      DO k = 0, km-2
         DO j = 0, jm-1
            DO i = 0, im-2
            
               lj = i + j*im1 + k*im1*jm + inti
               l1 = i + j*im + k*imjm
               l2 = l1 + incnodejk
               l3 = l1 + incnodejik
               l4 = l1 + incnodeji
               
               d13x = x(l3) - x(l1)
               d13y = y(l3) - y(l1)
               d13z = z(l3) - z(l1)
               
               d24x = x(l4) - x(l2)
               d24y = y(l4) - y(l2)
               d24z = z(l4) - z(l2)               
               nx = ONE_HALF * (d13y*d24z - d13z*d24y)
               ny = ONE_HALF * (d13z*d24x - d13x*d24z)
               nz = ONE_HALF * (d13x*d24y - d13y*d24x)
               
               surfx(lj) = sgn*nx
               surfy(lj) = sgn*ny
               surfz(lj) = sgn*nz
               snorm(lj) = SQRT(nx*nx + ny*ny + nz*nz)
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
               l2 = l1 + incnodeki
               l3 = l1 + incnodekij
               l4 = l1 + incnodekj

               d13x = x(l3) - x(l1)
               d13y = y(l3) - y(l1)
               d13z = z(l3) - z(l1)
               
               d24x = x(l4) - x(l2)
               d24y = y(l4) - y(l2) 
               d24z = z(l4) - z(l2)

               nx = ONE_HALF * (d13y*d24z - d13z*d24y)
               ny = ONE_HALF * (d13z*d24x - d13x*d24z)
               nz = ONE_HALF * (d13x*d24y - d13y*d24x)

               surfx(lk) = sgn*nx
               surfy(lk) = sgn*ny
               surfz(lk) = sgn*nz
               snorm(lk) = SQRT(nx*nx + ny*ny + nz*nz)
            END DO
         ENDDO
      ENDDO
!$OMP ENDDO NOWAIT

!$OMP END PARALLEL
      END
C ======= KCore/Metric/CompIntSurfF.for ------
