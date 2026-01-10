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

C attention cette routine est 2D : 
C pas de calcul des normales sur un hexaedre!
C ============================================================================
      SUBROUTINE k6normstructsurft(ni,nj,npts,xt, yt, zt, nxt, nyt, nzt)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj, npts
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)

C_OUT
      REAL_E nxt(0:(ni-1)*(nj-1)-1)              
      REAL_E nyt(0:(ni-1)*(nj-1)-1)              
      REAL_E nzt(0:(ni-1)*(nj-1)-1)              

C_LOCAL
      INTEGER_E i, j, indcell
      INTEGER_E ind1, ind2, ind3, ind4
      REAL_E l1x, l1y, l1z
      REAL_E l2x, l2y, l2z
      REAL_E surf1x, surf1y, surf1z
      REAL_E surf2x, surf2y, surf2z
      REAL_E x1, y1, z1

C!$OMP PARALLEL
C!$OMP& DEFAULT(SHARED)
C!$OMP& PRIVATE(i,j,ind1,ind2,ind3,ind4,indcell,l1x,l1y,l1z,x1,y1,z1)
C!$OMP& PRIVATE(l2x,l2y,l2z,surf1x,surf1y,surf1z,surf2x,surf2y,surf2z)
C!$OMP DO      
      DO j = 0, nj-2
      DO i = 0, ni-2
         ind1 = i+j*ni
         ind2 = ind1+1
         ind3 = ind2+ni
         ind4 = ind1+ni
         indcell = i+j*(ni-1)

C     AB x AC
         x1 = xt(ind1)
         y1 = yt(ind1)
         z1 = zt(ind1)
         l1x = xt(ind2) - x1
         l1y = yt(ind2) - y1
         l1z = zt(ind2) - z1

         l2x = xt(ind3) - x1
         l2y = yt(ind3) - y1
         l2z = zt(ind3) - z1  

         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     AC x AD
         l1x = xt(ind3) - x1
         l1y = yt(ind3) - y1
         l1z = zt(ind3) - z1

         l2x = xt(ind4) - x1
         l2y = yt(ind4) - y1
         l2z = zt(ind4) - z1  

         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)

         nxt(indcell) = ONE_HALF*(surf1x+surf2x)
         nyt(indcell) = ONE_HALF*(surf1y+surf2y)
         nzt(indcell) = ONE_HALF*(surf1z+surf2z)
      ENDDO
      ENDDO
C!$OMP END DO
C!$OMP END PARALLEL
      END
