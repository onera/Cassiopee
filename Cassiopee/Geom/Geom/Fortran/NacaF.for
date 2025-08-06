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
C Generate a naca profile (1D line)
C  ===========================================================================
C
C naca avec fermeture de Van Rouzaud
C
      SUBROUTINE k6naca1(e, npt, x, y, z)
C
      IMPLICIT NONE
C
C==============================================================================
C_IN
      REAL_E e             ! epaisseur du profil
      INTEGER_E npt        ! number of points (input and output)
C_OUT
      REAL_E x(1:npt)
      REAL_E y(1:npt)
      REAL_E z(1:npt)
C_LOCAL
      INTEGER_E np, nps2, n
      REAL_E pi, usnps2, xpn, ypn, en, epais, si, is
C==============================================================================
      epais = e/20.D0
      np = npt
      nps2 = (np-1)/2 
      usnps2 = 1.D0/nps2
      pi = 3.1415926D0
C     
C     Naca donne analytiquement (d apres Van Rouzaud)
C     
      DO n = 1, npt
         is = (n-2)/nps2 
         is = 2*is-1
         si = is
         en = n-nps2-1 
         en = SQRT(1.008930411365D0)*SIN(0.5D0*pi*en*usnps2)
         xpn = en*en
         ypn = si*epais*(0.2969D0*ABS(en)-0.126D0*xpn-0.3516D0*xpn*xpn
     &        +0.2843D0*xpn*xpn*xpn-0.1015D0*xpn*xpn*xpn*xpn) 
         x(n) = xpn/1.008930411365D0
         y(n) = ypn/1.008930411365D0
         z(n) = 0.D0
      ENDDO
      ypn = 0.5D0 * (y(1)+y(n))
      y(1) = ypn
      y(n) = ypn
  
      END
C
C-------------------------------------------------
C Idem mais avec fermeture a la leballeur
C-------------------------------------------------
C
      SUBROUTINE k6naca2( e, npt, x, y, z)
C
      IMPLICIT NONE
C
C==============================================================================
C_IN
      REAL_E e             ! epaisseur du profil
      INTEGER_E npt        ! number of points
C_OUT
      REAL_E x(1:npt)
      REAL_E y(1:npt)
      REAL_E z(1:npt)
C_LOCAL
      INTEGER_E np, nps2, n, no, nr
      REAL_E pi, usnps2, xpn, ypn, en, epais, si, is
      REAL_E hbfs2, bdabf
C==============================================================================
      epais = e/20.D0
      np = npt
      nps2 = (np-1)/2 
      usnps2 = 1./nps2
      pi = 3.1415926D0
C     
C     Naca donne analytiquement
C     
      nr = 1
      DO n = 1, npt
         is = (n-2)/nps2 
         is = 2*is-1
         si = is
         en = n-nps2-1 
         en = SQRT(1.008930411365D0)*SIN(0.5D0*pi*en*usnps2)
         xpn = en*en
         ypn = si*epais*(0.2969D0*ABS(en)-0.126D0*xpn-0.3516D0*xpn*xpn
     &        +0.2843D0*xpn*xpn*xpn-0.1015*xpn*xpn*xpn*xpn)
         IF (xpn.LE.1.) THEN
            x(nr) = xpn
            y(nr) = ypn
            z(nr) = 0.D0
            nr = nr+1
         ENDIF
      ENDDO
      npt = nr-1
C
C     Fermeture
C
      ypn = 0.5D0 * (y(1)+y(npt))
      xpn = 0.5D0 * (x(1)+x(npt))
      hbfs2 = ypn-y(1)
      bdabf = 70.D0
      IF (bdabf*ypn.GT.epais) THEN
         bdabf = epais/hbfs2
      ENDIF

      DO n = 1, npt/4
         no = npt-n+1
         y(n) = y(n)+hbfs2*EXP(MAX(bdabf*(x(n)-xpn),-30.D0))
         y(no) = y(no)-hbfs2*EXP(MAX(bdabf*(x(no)-xpn),-30.D0))
      ENDDO

      x(1) = 1.D0
      x(npt) = 1.D0

      END
