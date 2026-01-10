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

C ============================================================================
C Cubic hermite interpolation
C ============================================================================
      SUBROUTINE k6interp(im0, im, stota,                                
     &                    tabx0, taby0, tabz0, s0,                      
     &                    dx0, dy0, dz0,                               
     &                    tabx, taby, tabz, s)
!*************************************************************************
!
!    Cubic Hermite interpolation
!     By F.MONTIGNY-RANNOU 1994
!     
!*************************************************************************
      IMPLICIT NONE
C
      INTEGER_E im0, im, i, k
      REAL_E stota,t,coef,d1,d2,a,b,c,d,valcub
      REAL_E tabx0(im0),taby0(im0),tabz0(im0),s0(im0)
      REAL_E dx0(im0),dy0(im0),dz0(im0)
      REAL_E tabx(im),taby(im),tabz(im),s(im)
C
      valcub(a,b,c,d,t) = (a*(1.D0+2.D0*t)+c*t)*(1.D0-t)**2
     &                     +(b*(3.D0-2.D0*t)+d*(t-1.D0))*t**2

      tabx(1) = tabx0(1)
      taby(1) = taby0(1)
      tabz(1) = tabz0(1)
      tabx(im) = tabx0(im0)
      taby(im) = taby0(im0)
      tabz(im) = tabz0(im0)
      
      IF (im.EQ.2) THEN
         RETURN
      ELSE
         DO 10 i = 2, im-1
            DO 15 k = 1, im0-1
               IF (s(i).LE.s0(k+1)) GOTO 20
 15         continue
            k = im0-1
 20         continue
            t = (s(i)-s0(k))/(s0(k+1)-s0(k))
            coef = stota*(s0(k+1)-s0(k))
            d1 = dx0(k)*coef
            d2 = dx0(k+1)*coef
            tabx(i) = valcub(tabx0(k),tabx0(k+1),d1,d2,t)
            d1 = dy0(k)*coef
            d2 = dy0(k+1)*coef
            taby(i) = valcub(taby0(k),taby0(k+1),d1,d2,t)
            d1 = dz0(k)*coef
            d2 = dz0(k+1)*coef
            tabz(i) = valcub(tabz0(k),tabz0(k+1),d1,d2,t)
 10      CONTINUE
      ENDIF
      
      END

C=============================================================================
C attention: ici la BAR est isomorphe a un i-array 
      SUBROUTINE k6interpbar(im0, im, stota, 
     &     net0, cn10, cn20, tabx0, taby0, tabz0, s0,
     &     dx0, dy0, dz0, 
     &     net, cn1, cn2, tabx, taby, tabz, s)

      IMPLICIT NONE
C_IN 
      INTEGER_E im0, im
      REAL_E stota,t,coef,d1,d2,a,b,c,d,valcub
      REAL_E tabx0(0:im0-1),taby0(0:im0-1),tabz0(0:im0-1)
      REAL_E s0(0:im0-1)
      INTEGER_E net0             ! nb d'elts BAR
      INTEGER_E cn10(0:net-1)    ! connectivite BAR : 1ers sommets
      INTEGER_E cn20(0:net-1)    ! connectivite BAR : 2nds sommets
C_OUT
      REAL_E dx0(0:im0-1),dy0(0:im0-1),dz0(0:im0-1)
      REAL_E tabx(0:im-1),taby(0:im-1),tabz(0:im-1)
      INTEGER_E net 
      INTEGER_E cn1(0:net-1)    ! connectivite BAR : 1ers sommets
      INTEGER_E cn2(0:net-1)    ! connectivite BAR : 2nds sommets
      REAL_E s(0:im-1)
C_LOCAL
      INTEGER_E ind0, ind1, ind2
      INTEGER_E et, i, k
C-----------------------------------------------------------------------------
      valcub(a,b,c,d,t) = (a*(1.D0+2.D0*t)+c*t)*(1.D0-t)**2
     &                     +(b*(3.D0-2.D0*t)+d*(t-1.D0))*t**2

      tabx(0) = tabx0(0)
      taby(0) = taby0(0)
      tabz(0) = tabz0(0)
      tabx(im-1) = tabx0(im0-1)
      taby(im-1) = taby0(im0-1)
      tabz(im-1) = tabz0(im0-1)

      IF (im.EQ.2) THEN
         RETURN
      ELSE
         DO 60 i=1, im-2
            DO 65 k=0,im0-2
               IF (s(i).LE.s0(k+1)) GOTO 70
 65         continue
            k = im0-2
 70         continue
            t = (s(i)-s0(k))/(s0(k+1)-s0(k))
            coef = stota*(s0(k+1)-s0(k))
            d1 = dx0(k)*coef
            d2 = dx0(k+1)*coef
            tabx(i) = valcub(tabx0(k),tabx0(k+1),d1,d2,t)
            d1 = dy0(k)*coef
            d2 = dy0(k+1)*coef
            taby(i) = valcub(taby0(k),taby0(k+1),d1,d2,t)
            d1 = dz0(k)*coef
            d2 = dz0(k+1)*coef
            tabz(i) = valcub(tabz0(k),tabz0(k+1),d1,d2,t)
 60      CONTINUE
      ENDIF
      i = 1
      DO et = 0, net-1
         cn1(et) = i
         cn2(et) = i+1
         i = i + 1 
      ENDDO
      END
c=============================================================================
