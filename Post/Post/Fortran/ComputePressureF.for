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

C------------------------------------------------------------------------------
      SUBROUTINE k6computepressure(cellt, gamma, 
     &                             velo, ro, roe, p)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E   cellt ! Total Cell Number  
      REAL_E      velo(0:cellt-1,3) ! Velocity
      REAL_E      ro(0:cellt-1) ! density  
      REAL_E      roe(0:cellt-1) !  roe (5th conservative variable)
      REAL_E      gamma
C_OUT
      REAL_E      p(0:cellt-1)  ! Pressure Field 
      
C_LOCAL
      INTEGER_E   l   
      REAL_E      gam1
      REAL_E      vx2, vy2, vz2
c------------------------------------------------------------------------------

      gam1 = gamma - ONE

!$OMP PARALLEL DO PRIVATE(l, vx2,vy2,vz2)
      DO l = 0, cellt-1
         
         vx2 =  velo(l,1) * velo(l,1)
         vy2 =  velo(l,2) * velo(l,2)
         vz2 =  velo(l,3) * velo(l,3)

         p(l) = gam1 *
     &        (roe(l)-ONE_HALF*ro(l)*(vx2 + vy2 + vz2))

      END DO
!$OMP END PARALLEL DO
      END
C ========================== Post/Fortran/ComputePressureF.for =============
