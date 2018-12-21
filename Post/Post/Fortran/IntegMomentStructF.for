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
C  compute surface integral of the moment M (OM^ F), coordinates 
C      and field have the same size
C      I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
C      Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
C  ============================================================================
      SUBROUTINE k6integmomentstruct(ni, nj, cx, cy, cz, ratio, 
     &                               xt, yt, zt, surf, vx, vy, vz, 
     &                               result)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj
      REAL_E cx, cy, cz
      REAL_E ratio(0:ni*nj-1)
      REAL_E xt(0:ni*nj-1)
      REAL_E yt(0:ni*nj-1)
      REAL_E zt(0:ni*nj-1)
      REAL_E surf(0:(ni-1)*(nj-1)-1)
      REAL_E vx(0:ni*nj-1)
      REAL_E vy(0:ni*nj-1)
      REAL_E vz(0:ni*nj-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i, j, ni1
      INTEGER_E ind1, ind2, ind3, ind4, ind
      REAL_E f1x, f2x, f3x, f4x
      REAL_E f1y, f2y, f3y, f4y
      REAL_E f1z, f2z, f3z, f4z
      REAL_E m1x, m2x, m3x, m4x
      REAL_E m1y, m2y, m3y, m4y
      REAL_E m1z, m2z, m3z, m4z
      REAL_E dx, dy, dz
      REAL_E rind, sind
      REAL_E res1, res2, res3
C==============================================================================
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      ni1 = ni-1

      DO j = 0, nj-2
         DO i = 0, ni-2
            ind1 = i+j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind3 + ni
            ind = i + j*ni1

            rind = ratio(ind1)
            dx = xt(ind1) - cx
            dy = yt(ind1) - cy
            dz = zt(ind1) - cz

            f1x = rind * vx(ind1)
            f1y = rind * vy(ind1)
            f1z = rind * vz(ind1)
            m1x = dy*f1z - dz*f1y
            m1y = dz*f1x - dx*f1z
            m1z = dx*f1y - dy*f1x

            rind = ratio(ind2)
            dx = xt(ind2) - cx
            dy = yt(ind2) - cy 
            dz = zt(ind2) - cz
            f2x = rind * vx(ind2)
            f2y = rind * vy(ind2)
            f2z = rind * vz(ind2)
            m2x = dy*f2z - dz*f2y
            m2y = dz*f2x - dx*f2z
            m2z = dx*f2y - dy*f2x

            rind = ratio(ind3)
            dx = xt(ind3) - cx
            dy = yt(ind3) - cy
            dz = zt(ind3) - cz
            f3x = rind * vx(ind3)
            f3y = rind * vy(ind3)
            f3z = rind * vz(ind3)
            m3x = dy*f3z - dz*f3y
            m3y = dz*f3x - dx*f3z
            m3z = dx*f3y - dy*f3x

            rind = ratio(ind4)
            dx = xt(ind4) - cx
            dy = yt(ind4) - cy
            dz = zt(ind4) - cz
            f4x = rind * vx(ind4)
            f4y = rind * vy(ind4)
            f4z = rind * vz(ind4)
            m4x = dy*f4z - dz*f4y
            m4y = dz*f4x - dx*f4z
            m4z = dx*f4y - dy*f4x

            sind = surf(ind)
            res1 = res1 + sind*(m1x+m2x+m3x+m4x)
            res2 = res2 + sind*(m1y+m2y+m3y+m4y)
            res3 = res3 + sind*(m1z+m2z+m3z+m4z)
         ENDDO
      ENDDO
      result(0) = ONE_FOURTH * res1
      result(1) = ONE_FOURTH * res2 
      result(2) = ONE_FOURTH * res3 
      
      END
C =============================================================================
C  compute linear integral of the moment M (OM^ F), coordinates 
C      and field have the same size
C      I(AB) = LENGTH(ABCD)*(F(A)+F(B))/2
C  ============================================================================
      SUBROUTINE k6integmomentstruct1d(ni, cx, cy, cz, ratio, xt, yt,zt,  
     &                                 length, vx, vy, vz, result)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni
      REAL_E cx, cy, cz
      REAL_E ratio(0:ni-1)
      REAL_E xt(0:ni-1)
      REAL_E yt(0:ni-1)
      REAL_E zt(0:ni-1)      
      REAL_E length(0:ni-2)
      REAL_E vx(0:ni-1)
      REAL_E vy(0:ni-1)
      REAL_E vz(0:ni-1)
C_OUT
      REAL_E result(0:2)             
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind
      REAL_E f1x, f2x
      REAL_E f1y, f2y
      REAL_E f1z, f2z
      REAL_E m1x, m2x
      REAL_E m1y, m2y
      REAL_E m1z, m2z
      REAL_E ri, li, dx, dy, dz
      REAL_E res1, res2, res3
C
C==============================================================================
C
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      DO i = 0, ni-2
         ind1 = i
         ind2 = i+1

         ind = i
         li = length(ind)
C
         ri = ratio(ind1)
         f1x = ri * vx(ind1)
         f1y = ri * vy(ind1)
         f1z = ri * vz(ind1)

         dx = xt(ind1) - cx
         dy = yt(ind1) - cy
         dz = zt(ind1) - cz
         
         m1x = dy*f1z - dz*f1y
         m1y = dz*f1x - dx*f1z
         m1z = dx*f1y - dy*f1x
C
         ri = ratio(ind2)
         f2x = ri * vx(ind2)
         f2y = ri * vy(ind2)
         f2z = ri * vz(ind2)
         
         dx = xt(ind2) - cx
         dy = yt(ind2) - cy
         dz = zt(ind2) - cz

         m2x = dy*f2z - dz*f2y
         m2y = dz*f2x - dx*f2z 
         m2z = dx*f2y - dy*f2x
C
         res1 = res1 + li *(m1x+m2x)
         res2 = res2 + li *(m1y+m2y)
         res3 = res3 + li *(m1z+m2z)
      ENDDO

      result(0) = ONE_HALF * res1
      result(1) = ONE_HALF * res2 
      result(2) = ONE_HALF * res3 

      END

C==============================================================================
C Compute surface integral of the moment M (OM^F), coordinates 
C       are defined in nodes and F is defined in center
C  ============================================================================
      SUBROUTINE k6integmomentstructnodecenter(ni, nj, cx, cy, cz, 
     &                                         ratio, xt, yt, zt, surf,
     &                                         vx, vy, vz, result)
      IMPLICIT NONE   

C==============================================================================
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj
      REAL_E cx, cy, cz
      REAL_E ratio(0:(ni-1)*(nj-1)-1)
      REAL_E xt(0:ni*nj-1)
      REAL_E yt(0:ni*nj-1)
      REAL_E zt(0:ni*nj-1)
      REAL_E surf(0:(ni-1)*(nj-1)-1)
      REAL_E vx(0:(ni-1)*(nj-1)-1)
      REAL_E vy(0:(ni-1)*(nj-1)-1)
      REAL_E vz(0:(ni-1)*(nj-1)-1)
C_OUT
      REAL_E result(0:2)
      
C_LOCAL
      INTEGER_E i, j, ni1
      INTEGER_E ind, ind1, ind2, ind3, ind4
      REAL_E mx, my, mz
      REAL_E res1, res2, res3
      REAL_E centerx, centery, centerz
      REAL_E srind
C
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      ni1 = ni-1
      DO j = 0, nj-2
         DO i = 0, ni-2
            ind1 = i+j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind2 + 1
            ind = i + j*ni1
            srind = surf(ind)*ratio(ind)
            
            centerx = xt(ind1)+xt(ind2)+xt(ind3)+xt(ind4)
            centerx = ONE_FOURTH*centerx
            centerx = centerx - cx
            
            centery = yt(ind1)+yt(ind2)+yt(ind3)+yt(ind4)
            centery = ONE_FOURTH*centery
            centery = centery - cy
            
            centerz = zt(ind1)+zt(ind2)+zt(ind3)+zt(ind4)
            centerz = ONE_FOURTH*centerz
            centerz = centerz - cz

            mx = centery*vz(ind) - centerz*vy(ind)
            my = centerz*vx(ind) - centerx*vz(ind)
            mz = centerx*vy(ind) - centery*vx(ind)

            res1 = res1 + srind*mx
            res2 = res2 + srind*my
            res3 = res3 + srind*mz
         ENDDO
      ENDDO

      result(0) = res1
      result(1) = res2
      result(2) = res3

      END
C==============================================================================
C compute linear integral of the moment M (OM^F), coordinates 
C       are defined in nodes and F is defined in center
C  ============================================================================
      SUBROUTINE k6integmomentstructnodecenter1d(ni, cx, cy, cz, ratio, 
     &                                           xt, yt, zt, length, 
     &                                           vx, vy, vz, result)
      IMPLICIT NONE   

C
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni
      REAL_E cx, cy, cz
      REAL_E ratio(0:ni-2)
      REAL_E xt(0:ni-1)
      REAL_E yt(0:ni-1)
      REAL_E zt(0:ni-1)
      REAL_E length(0:ni-2)
      REAL_E vx(0:ni-2)
      REAL_E vy(0:ni-2)
      REAL_E vz(0:ni-2)
      
C_OUT
      REAL_E result(0:2)
      
C_LOCAL
      INTEGER_E i
      INTEGER_E ind, ind1, ind2
      REAL_E resx, resy, resz
      REAL_E fx, fy, fz
      REAL_E mx, my, mz, rlind
      REAL_E centerx, centery, centerz
C
      resx = 0.D0
      resy = 0.D0
      resz = 0.D0

      DO i = 0, ni-2
         ind1 = i
         ind2 = i+1
         ind = i
        
         rlind = length(ind) *  ratio(ind)

         centerx = ONE_HALF * (xt(ind1)+xt(ind2)) - cx
         centery = ONE_HALF * (yt(ind1)+yt(ind2)) - cy
         centerz = ONE_HALF * (zt(ind1)+zt(ind2)) - cz

         mx = centery*vz(ind) - centerz*vy(ind)
         my = centerz*vx(ind) - centerx*vz(ind)
         mz = centerx*vy(ind) - centery*vx(ind)

         resx = resx + rlind * mx
         resy = resy + rlind * my
         resz = resz + rlind * mz
      ENDDO
      result(0) = resx
      result(1) = resy
      result(2) = resz 
      
      END
C ====================Post/Fortran/IntegMomentStructF.for=====================
