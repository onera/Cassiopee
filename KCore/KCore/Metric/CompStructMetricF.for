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
C  Calcul du volume de toutes les cellules  et des surface des interfaces
C  CAS STRUCTURE
C  ============================================================================
      SUBROUTINE  k6compstructmetric(ni, nj, nk, nbcell, nbInt,
     &     nbInti, nbIntj, nbIntk, 
     &     x, y, z, vol, surfx, surfy, surfz,
     &     snorm, cix, ciy, ciz)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E  ni, nj, nk
      INTEGER_E  nbInti, nbIntj, nbIntk 
      INTEGER_E  nbCell         ! number of cells
      INTEGER_E  nbInt          ! number of interfaces
      REAL_E x(0:ni*nj*nk-1), y(0:ni*nj*nk-1), z(0:ni*nj*nk-1) ! coordonnees des noeuds du maillage
C_OUT
      REAL_E surfx(0:nbInt-1), surfy(0:nbInt-1), surfz(0:nbInt-1)
      REAL_E     snorm(0:nbInt-1)  !  norm of surf
      REAL_E     vol(0:nbCell-1)
      REAL_E     cix(0:nbInt-1), ciy(0:nbInt-1), ciz(0:nbInt-1)

C_LOCAL
      INTEGER_E  i, j, k, nbIntij
      INTEGER_E  ni1, nj1, ni1nj1, nk1
      INTEGER_E  incInti, incIntj, incIntk
      INTEGER_E  lv, n1, n2, n3, n4, n5, n6
      INTEGER_E  nbnode
      REAL_E     v1, v2, v3, v4, v5, v6
      REAL_E     onethird
C==============================================================================
C
      nbnode = ni*nj*nk
      onethird = ONE/THREE
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
      ni1nj1 = ni1*nj1
      nbIntij = nbInti + nbIntj
C
      call k6compcenterinterface(
     &     ni, nj, nk, nbnode, 
     &     nbInt, nbInti, nbIntij,
     &     x, y, z, cix, ciy, ciz)
 
      call k6compintsurf( 
     &     ni, nj, nk, nbnode, 
     &     nbInt, nbInti, nbIntij, 
     &     x, y, z, surfx, surfy, surfz, snorm)
C
      incInti = 1
      incIntj = ni1
      incIntk = ni1nj1
C
!$OMP PARALLEL PRIVATE(i,j,k,lv,n1,v1,n2,v2,n3,v3,n4,v4,n5,v5,n6,v6)
!$OMP DO
      DO k = 1, nk1
      DO j = 1, nj1
      DO i = 1, ni1
         lv = i-1 + (j-1)*ni1 + (k-1)*ni1nj1
         n1 = i-1 + (j-1)*ni  + (k-1)*ni*nj1
         v1 =  cix(n1) * surfx(n1) 
     &        + ciy(n1) * surfy(n1)
     &        + ciz(n1) * surfz(n1)
         
         n2 = n1 + incInti
         v2 =  cix(n2) * surfx(n2) 
     &        + ciy(n2) * surfy(n2)
     &        + ciz(n2) * surfz(n2)
    
         n3 = i-1 + (j-1)*ni1 + (k-1)*ni1*nj + nbInti
         v3 =  cix(n3) * surfx(n3) 
     &        + ciy(n3) * surfy(n3)
     &        + ciz(n3) * surfz(n3)
    
         n4 = n3 + incIntj
         v4 =  cix(n4) * surfx(n4) 
     &        + ciy(n4) * surfy(n4)
     &        + ciz(n4) * surfz(n4)
         
         n5 =  i-1 + (j-1)*ni1 + (k-1)*ni1nj1 + nbIntij
         v5 =  cix(n5) * surfx(n5) 
     &        + ciy(n5) * surfy(n5)
     &        + ciz(n5) * surfz(n5)
         
         n6 = n5 + incIntk
         v6 =  cix(n6) * surfx(n6) 
     &        + ciy(n6) * surfy(n6)
     &        + ciz(n6) * surfz(n6)
         
         vol(lv) = (v2 - v1 + v4 - v3 + v6 - v5) * onethird
      ENDDO
      ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END
C========================= Metric/CompStructMetricF.for ======================
