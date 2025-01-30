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

C ============================================================================
C compute surface integral of the moment M (OM^ F), coordinates 
C       and field have the same size
C       I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
C       Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
C  ============================================================================
      SUBROUTINE k6integmomentunstruct(nbt, size, cn, cx, cy, cz, ratio, 
     &                                 xt, yt, zt, surf, vx, vy, vz, 
     &                                 result)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,3)
      REAL_E cx, cy, cz
      REAL_E ratio(0:size-1)
      REAL_E xt(0:size-1)
      REAL_E yt(0:size-1)      
      REAL_E zt(0:size-1)
      REAL_E surf(0:nbT-1)
      REAL_E vx(0:size-1)
      REAL_E vy(0:size-1)
      REAL_E vz(0:size-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i, ind1, ind2, ind3
      REAL_E fx, fy, fz
      REAL_E m1x, m2x, m3x
      REAL_E m1y, m2y, m3y
      REAL_E m1z, m2z, m3z
      REAL_E res1, res2, res3
      REAL_E dx, dy, dz, r1, r2, r3, si
      REAL_E ONE_THIRD
C
      ONE_THIRD = ONE/THREE

      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
C
         r1 = ratio(ind1)
         dx = xt(ind1) - cx
         dy = yt(ind1) - cy
         dz = zt(ind1) - cz
         fx = vx(ind1)
         fy = vy(ind1) 
         fz = vz(ind1) 

         m1x = dy*fz - dz*fy
         m1y = dz*fx - dx*fz
         m1z = dx*fy - dy*fx
C
         r2 = ratio(ind2)
         dx = xt(ind2) - cx
         dy = yt(ind2) - cy
         dz = zt(ind2) - cz
         fx = vx(ind2)
         fy = vy(ind2) 
         fz = vz(ind2) 

         m2x = dy*fz - dz*fy
         m2y = dz*fx - dx*fz
         m2z = dx*fy - dy*fx
C
         r3 = ratio(ind3)
         dx = xt(ind3) - cx
         dy = yt(ind3) - cy
         dz = zt(ind3) - cz
         fx = vx(ind3)
         fy = vy(ind3) 
         fz = vz(ind3) 

         m3x = dy*fz - dz*fy
         m3y = dz*fx - dx*fz
         m3z = dx*fy - dy*fx
C
         si = surf(i)
         res1 = res1 + si * (r1*m1x + r2*m2x + r3*m3x)
         res2 = res2 + si * (r1*m1y + r2*m2y + r3*m3y)
         res3 = res3 + si * (r1*m1z + r2*m2z + r3*m3z)
      ENDDO

      result(0) = res1 * ONE_THIRD
      result(1) = res2 * ONE_THIRD
      result(2) = res3 * ONE_THIRD

      END
C =============================================================================
C compute linear integral of the moment M (OM^ F), coordinates 
C       and field have the same size
C  ============================================================================
      SUBROUTINE k6integmomentunstruct1d(nbt, size, cn, cx, cy, cz,  
     &                                   ratio, xt, yt, zt, length, 
     &                                   vx, vy, vz, result)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,2)
      REAL_E cx, cy, cz
      REAL_E ratio(0:size-1)
      REAL_E xt(0:size-1)
      REAL_E yt(0:size-1)
      REAL_E zt(0:size-1)
      REAL_E length(0:nbT-1)
      REAL_E vx(0:size-1)
      REAL_E vy(0:size-1)
      REAL_E vz(0:size-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i, ind1, ind2
      REAL_E fx, fy, fz
      REAL_E m1x, m2x
      REAL_E m1y, m2y
      REAL_E m1z, m2z
      REAL_E res1, res2, res3, li
      REAL_E dx, dy, dz
C
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         li = length(i)
         
         dx = xt(ind1)-cx
         dy = yt(ind1)-cy
         dz = zt(ind1)-cz
         
         fx = ratio(ind1)*vx(ind1)
         fy = ratio(ind1)*vy(ind1)
         fz = ratio(ind1)*vz(ind1)

         m1x = dy*fz - dz*fy
         m1y = dz*fx - dx*fz
         m1z = dx*fy - dy*fx
C
         dx = xt(ind2)-cx
         dy = yt(ind2)-cy
         dz = zt(ind2)-cz

         fx = ratio(ind2)*vx(ind2)
         fy = ratio(ind2)*vy(ind2)
         fz = ratio(ind2)*vz(ind2)
         
         m2x = dy*fz - dz*fy
         m2y = dz*fx - dx*fz
         m2z = dx*fy - dy*fx
C
         res1 = res1 + li * (m1x+m2x)
         res2 = res2 + li * (m1y+m2y)
         res3 = res3 + li * (m1z+m2z)
      ENDDO

      result(0) = res1 * ONE_HALF 
      result(1) = res2 * ONE_HALF 
      result(2) = res3 * ONE_HALF 
      END
C =============================================================================
C compute surface integral of the moment M (OM^ F), coordinates 
C       are defined in nodes and F is defined in center, unstructured case
C  ============================================================================
      SUBROUTINE k6integmomentunsnodecenter(nbt, size, cn, cx, cy, cz, 
     &                                      ratio, xt, yt, zt, surf,
     &                                      vx, vy, vz, result) 

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,3)
      REAL_E cx, cy, cz
      REAL_E ratio(0:nbt-1)
      REAL_E xt(0:size-1)
      REAL_E yt(0:size-1)
      REAL_E zt(0:size-1)      
      REAL_E surf(0:nbt-1)
      REAL_E vx(0:nbt-1)
      REAL_E vy(0:nbt-1)
      REAL_E vz(0:nbt-1)
C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i, ind1, ind2, ind3
      REAL_E mx, my, mz, sri
      REAL_E centerx, centery, centerz
      REAL_E res1, res2, res3
      REAL_E ONE_THIRD
C
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      ONE_THIRD = ONE/THREE

      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1

         centerx = ONE_THIRD * ( xt(ind1)+xt(ind2)+xt(ind3) ) - cx
         centery = ONE_THIRD * ( yt(ind1)+yt(ind2)+yt(ind3) ) - cy
         centerz = ONE_THIRD * ( zt(ind1)+zt(ind2)+zt(ind3) ) - cz

         mx = centery*vz(i) - centerz*vy(i)
         my = centerz*vx(i) - centerx*vz(i)
         mz = centerx*vy(i) - centery*vx(i)        

         sri = surf(i)*ratio(i)
         res1 = res1 + sri * mx
         res2 = res2 + sri * my
         res3 = res3 + sri * mz
      ENDDO

      result(0) = res1
      result(1) = res2
      result(2) = res3

      END    
C =============================================================================
C compute surface integral of the moment M (OM^ F), coordinates 
C       are defined in nodes and F is defined in center, unstructured case
C  ============================================================================
      SUBROUTINE k6integmomentunsnodecenter1d(nbt, size, cn, cx, cy, cz, 
     &                                        ratio, xt, yt, zt, surf,
     &                                        vx, vy, vz, result) 

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,2)
      REAL_E cx, cy, cz
      REAL_E ratio(0:nbt-1)
      REAL_E xt(0:size-1)
      REAL_E yt(0:size-1)
      REAL_E zt(0:size-1)
      REAL_E surf(0:nbt-1)
      REAL_E vx(0:nbt-1)
      REAL_E vy(0:nbt-1)
      REAL_E vz(0:nbt-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i, ind1, ind2
      REAL_E mx, my, mz
      REAL_E centerx, centery, centerz
      REAL_E res1, res2, res3, sri
C
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         
         centerx = ONE_HALF * (xt(ind1) + xt(ind2) ) - cx
         centery = ONE_HALF * (yt(ind1) + yt(ind2) ) - cy
         centerz = ONE_HALF * (zt(ind1) + zt(ind2) ) - cz

         mx = centery*vz(i) - centerz*vy(i)
         my = centerz*vx(i) - centerx*vz(i)
         mz = centerx*vy(i) - centery*vx(i)

         sri = surf(i) * ratio(i)

         res1 = res1 + sri * mx
         res2 = res2 + sri * my
         res3 = res3 + sri * mz
      ENDDO

      result(0) = res1
      result(1) = res2
      result(2) = res3

      END
C ==================== Post/Fortran/IntegMomentUnstrF.for =====================
