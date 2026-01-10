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

c Coefficients pour le lisseur elliptique
c
      subroutine k6coeffttm(x, y, n, m, a11, a12, a22, 
     &                      b11, b12, b22,
     &                      c11, c12, c22)
c
      IMPLICIT NONE
      INTEGER_E n, m
      REAL_E x(n,m), y(n,m)
      REAL_E a11(n,m), a12(n,m), a22(n,m)
      REAL_E b11(n,m), b12(n,m), b22(n,m)
      REAL_E c11(n,m), c12(n,m), c22(n,m)

      REAL_E g11, g12, g22, xxi, yxi, xet, yet, rg, g
      REAL_E p11, p12, p22, p11p11, p11p12, p11p22, p11prg,
     &     p12p12, p12p22, p12prg, p22p22, p22prg, prgprg
      REAL_E xm11(2,2), xm12(2,2), xm22(2,2)
      REAL_E dx, dy, tdx, tdy
      INTEGER_E i,j
C
C tdx, tdy : ?
      dx = 1.D0/m
      dy = 1.D0/n
      tdx = 2.D0*dx
      tdy = 2.D0*dy

      DO j = 2, m-1
         DO i = 2, n-1
            call k6tngts(i, j, x, y, n, m, tdx, tdy, xxi, yxi, xet, yet)
            call k6mets(xxi, yxi, xet, yet, g11, g12, g22, rg, g)
            call k6tnsr(xxi, yxi, xet, yet, xm11, xm12, xm22)
            call k6partials(g11, g12, g22, rg, g,
     &           p11, p12, p22, p11p11, p11p12, p11p22, p11prg,
     &           p12p12, p12p22, p12prg, p22p22, p22prg, prgprg)
            call k6t11(xm11, xm12, xm22,
     &           g11, g12, g22, rg, g,
     &           p11,p12,p22,p11p11,p11p12,p11p22,p11prg,
     $           p12p12,p12p22,p12prg,p22p22,p22prg,prgprg,
     &           a11(i,j), a12(i,j), a22(i,j))
            call k6t12(xm11, xm12, xm22,
     &           g11, g12, g22, rg, g,
     &           p11,p12,p22,p11p11,p11p12,p11p22,p11prg,
     $           p12p12,p12p22,p12prg,p22p22,p22prg,prgprg,
     &           b11(i,j), b12(i,j), b22(i,j))
            call k6t22(xm11, xm12, xm22,
     &           g11, g12, g22, rg, g,
     &           p11, p12, p22, p11p11, p11p12, p11p22, p11prg,
     $           p12p12, p12p22, p12prg, p22p22, p22prg, prgprg,
     &           c11(i,j), c12(i,j), c22(i,j))
C            WRITE(*,*) c11(i,j),a11(i,j),b11(i,j)
         ENDDO
      ENDDO
      RETURN
      END
