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
C Calcul les vecteurs normals a des triangles

      SUBROUTINE k6normunstructsurf(nt, nv, cn, xt, yt, zt, nsurf)
      IMPLICIT NONE

# include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E nt              ! nombre de triangles
      INTEGER_E nv              ! nombre de noeuds
      INTEGER_E cn(0:nt-1,3)    ! connectivite elts-noeuds
      REAL_E xt(nv)             ! coord des noeuds
      REAL_E yt(nv)
      REAL_E zt(nv)

C_OUT
      REAL_E nsurf(0:nt-1,3)              

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E l1x, l1y, l1z
      REAL_E l2x, l2y, l2z
      REAL_E l3x, l3y, l3z
      REAL_E surfx, surfy, surfz
C==============================================================================
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,l1x,l1y,l1z,l2x,l2y,l2z)
!$OMP& PRIVATE (surfx,surfy,surfz)
!$OMP DO
      DO i = 0, nt-1
         ind1 = cn(i,1)
         ind2 = cn(i,2)
         ind3 = cn(i,3)
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind1)-xt(ind3)
         l2y = yt(ind1)-yt(ind3)
         l2z = zt(ind1)-zt(ind3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         nsurf(i,1) = ONE_HALF*surfx
         nsurf(i,2) = ONE_HALF*surfy
         nsurf(i,3) = ONE_HALF*surfz
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END
C=============================================================================
      SUBROUTINE k6normunstructsurft(nt, nv, cn, xt, yt, zt, 
     &     nxt, nyt, nzt)
      IMPLICIT NONE

# include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E nt              ! nombre de triangles
      INTEGER_E nv              ! nombre de noeuds
      INTEGER_E cn(0:nt-1,3)    ! connectivite elts-noeuds
      REAL_E xt(nv)             ! coord des noeuds
      REAL_E yt(nv)
      REAL_E zt(nv)
C_OUT
      REAL_E nxt(0:nt-1)              
      REAL_E nyt(0:nt-1)              
      REAL_E nzt(0:nt-1)              

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E l1x, l1y, l1z
      REAL_E l2x, l2y, l2z
      REAL_E l3x, l3y, l3z
      REAL_E surfx, surfy, surfz
C==============================================================================
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,l1x,l1y,l1z,l2x,l2y,l2z)
!$OMP& PRIVATE (surfx,surfy,surfz)
!$OMP DO
      DO i = 0, nt-1
         ind1 = cn(i,1)
         ind2 = cn(i,2)
         ind3 = cn(i,3)
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind1)-xt(ind3)
         l2y = yt(ind1)-yt(ind3)
         l2z = zt(ind1)-zt(ind3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         nxt(i) = ONE_HALF*surfx
         nyt(i) = ONE_HALF*surfy
         nzt(i) = ONE_HALF*surfz
      ENDDO
!$OMP END DO
!$OMP END PARALLEL     
      END
C ================ KCore/Metric/CompNormUnstrSurfF.for ====================
