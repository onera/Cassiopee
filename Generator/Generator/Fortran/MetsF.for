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

c==============================================================================
      SUBROUTINE k6mets(xxi, yxi, xet, yet,
     &     g11, g12, g22, rg, g)
c
c     computes the local metrics g11, g12, g22, root-g, and g
c     from the local tangents
c
      REAL_E xxi, yxi, xet, yet
      REAL_E g11, g12, g22, rg, g

      g11 = xxi*xxi + yxi*yxi
      g12 = xxi*xet + yxi*yet
      g22 = xet*xet + yet*yet
      rg = xxi*yet - xet*yxi
      g = rg*rg

      RETURN
      END
c
c==============================================================================
      SUBROUTINE k6tngts(i, j, x, y, n, m, tdx, tdy,
     &     xxi, yxi, xet, yet)
c
c     computes components of local tangents using centered-
c     difference approximation
c
      INTEGER_E i, j
      INTEGER_E n, m
      REAL_E x(n,m), y(n,m)
      REAL_E tdx, tdy
      REAL_E xxi, yxi, xet, yet
c
      xxi = ( x(i+1,j)-x(i-1,j) )/tdx
      xet = ( x(i,j+1)-x(i,j-1) )/tdy
      yxi = ( y(i+1,j)-y(i-1,j) )/tdx
      yet = ( y(i,j+1)-y(i,j-1) )/tdy
c
      RETURN
      END
C
C==============================================================================
      SUBROUTINE k6tnsr(xxi, yxi, xet, yet, xm11, xm12, xm22)
c
      REAL_E xxi, yxi, xet, yet
      REAL_E xm11(2,2), xm12(2,2), xm22(2,2)
c
      xm11(1,1) = xxi*xxi
      xm11(1,2) = xxi*yxi
      xm11(2,2) = yxi*yxi
c
      xm12(1,1) = 2.*xxi*xet
      xm12(1,2) = xxi*yet+yxi*xet
      xm12(2,2) = 2.*yxi*yet
c
      xm22(1,1) = xet*xet
      xm22(1,2) = xet*yet
      xm22(2,2) = yet*yet
c
      RETURN
      END

C==============================================================================
      SUBROUTINE k6pq(i, j, p, q, n, m, dy, q0, xlam, ets)
c
      INTEGER_E i,j
      INTEGER_E n,m
      REAL_E p(n,m), q(n,m)
      REAL_E et, r, sgnr
c
      p(i,j) = 0.D0
c
      et = (j-1)*dy
      r = et - ets
      if (r.ge.0.) sgnr = 1.
      if (r.lt.0.) sgnr = -1.
      q(i,j) = -q0*sgnr*exp( -xlam*abs(r) )
c
      RETURN
      END

c==============================================================================
      SUBROUTINE k6t11(xm11, xm12, xm22,
     &     g11, g12, g22, rg, g,
     &     p11,p12,p22,p11p11,p11p12,p11p22,p11prg,
     $     p12p12,p12p22,p12prg,p22p22,p22prg,prgprg,
     &     a11, a12, a22)
c
      IMPLICIT NONE
      REAL_E xm11(2,2), xm12(2,2), xm22(2,2)
      REAL_E g11, g12, g22, rg, g
      REAL_E p11,p12,p22,p11p11,p11p12,p11p22,p11prg
      REAL_E p12p12,p12p22,p12prg,p22p22,p22prg,prgprg
      REAL_E a11, a12, a22
      REAL_E s11, s12, s22
c
      s11 = 2.*p11p11 + 2.*p11prg*g22/rg + prgprg*g22*g22/g/2.
      s12 = p11p12 - p11prg*g12/rg + p12prg*g22/rg/2. -
     $      prgprg*g12*g22/g/2.
      s22 = p12p12/2. - p12prg*g12/rg + prgprg*g12*g12/g/2.
c
      a11 = p11 + s11*xm11(1,1) + s12*xm12(1,1) + s22*xm22(1,1)
      a12 =      s11*xm11(1,2) + s12*xm12(1,2) + s22*xm22(1,2)
      a22 = p11 + s11*xm11(2,2) + s12*xm12(2,2) + s22*xm22(2,2)
c
      RETURN
      END

c==============================================================================
      SUBROUTINE k6t12(xm11, xm12, xm22,
     &     g11, g12, g22, rg, g,
     &     p11,p12,p22,p11p11,p11p12,p11p22,p11prg,
     $     p12p12,p12p22,p12prg,p22p22,p22prg,prgprg,
     &     b11, b12, b22)
c
      IMPLICIT NONE
      REAL_E xm11(2,2), xm12(2,2), xm22(2,2)
      REAL_E g11, g12, g22, rg, g
      REAL_E p11,p12,p22,p11p11,p11p12,p11p22,p11prg
      REAL_E p12p12,p12p22,p12prg,p22p22,p22prg,prgprg
      REAL_E b11, b12, b22
      REAL_E s11, s12, s22
c
      s11 = 2.*p11p12 - 2.*p11prg*g12/rg - prgprg*g12*g22/g
      s12 = 2.*p11p22 + p11prg*g11/rg + p12p12/2. + p22prg*g22/rg +
     $      prgprg*(g11*g22+g12*g12)/g/2.
      s22 = 2.*p12p22 - 2.*p22prg*g12/rg - prgprg*g11*g12/g
c
      b11 = p12 + rg*p12prg + s11*xm11(1,1) + s12*xm12(1,1) + 
     $      s22*xm22(1,1)
      b12 = s11*xm11(1,2) + s12*xm12(1,2) + s22*xm22(1,2)
      b22 = p12 + rg*p12prg + s11*xm11(2,2) + s12*xm12(2,2) + 
     $      s22*xm22(2,2)
c
      RETURN
      END

c==============================================================================
      SUBROUTINE k6t22(xm11, xm12, xm22,
     &     g11, g12, g22, rg, g,
     &     p11,p12,p22,p11p11,p11p12,p11p22,p11prg,
     $     p12p12,p12p22,p12prg,p22p22,p22prg,prgprg,
     &     c11, c12, c22)
c
      IMPLICIT NONE
      REAL_E xm11(2,2), xm12(2,2), xm22(2,2)
      REAL_E g11, g12, g22, rg, g
      REAL_E p11,p12,p22,p11p11,p11p12,p11p22,p11prg
      REAL_E p12p12,p12p22,p12prg,p22p22,p22prg,prgprg
      REAL_E c11, c12, c22
      REAL_E s11, s12, s22

      s11 = p12p12/2. - p12prg*g12/rg + prgprg*g12*g12/g/2.
      s12 = p12p22 - p22prg*g12/rg + p12prg*g11/rg/2. -
     $      prgprg*g12*g11/g/2.
      s22 = 2.*p22p22 + 2.*p22prg*g11/rg + prgprg*g11*g11/g/2.
c
      c11 = p22 + s11*xm11(1,1) + s12*xm12(1,1) + s22*xm22(1,1)
      c12 =       s11*xm11(1,2) + s12*xm12(1,2) + s22*xm22(1,2)
      c22 = p22 + s11*xm11(2,2) + s12*xm12(2,2) + s22*xm22(2,2)
c
      return
      end

C==============================================================================
      SUBROUTINE k6partials(g11, g12, g22, rg, g,
     &     p11, p12, p22, p11p11, p11p12, p11p22, p11prg,
     &     p12p12, p12p22, p12prg, p22p22, p22prg, prgprg)
c
      INTEGER_E imth
      REAL_E g11, g12, g22, rg, g
      REAL_E p11, p12, p22, p11p11, p11p12, p11p22, p11prg,
     &     p12p12, p12p22, p12prg, p22p22, p22prg, prgprg
c
      imth = 2

      p11 = 0.D0
      p12 = 0.D0
      p22 = 0.D0
c
      p11p11 = 0.D0
      p11p12 = 0.D0
      p11p22 = 0.D0
      p11prg = 0.D0
      p12p12 = 0.D0
      p12p22 = 0.D0
      p12prg = 0.D0
      p22p22 = 0.D0
      p22prg = 0.D0
      prgprg = 0.D0
c
      if (imth.eq.1) then
c
c    ao
c
         p11 = g22
         p22 = g11
         p11p22 = 1.D0
         return
      end if
c
      if (imth.eq.2) then
c
c    ttm
c
         p11 = 1.D0/rg
         p22 = 1.D0/rg
         p11prg = -1.D0/g
         p22prg = -1.D0/g
         prgprg = 2.*(g11+g22)/g/rg
         return
      end if
c
      if (imth.eq.3) then
c
c    length
c
         p11 = 1.D0
         p22 = 1.D0       
         return
      end if
c
      if (imth.eq.4) then
c
c    area
c
         p11 = g22
         p12 = -2.*g12
         p22 = g11
         p11p22 = 1.D0
         p12p12 = -2.D0
         return
      end if
c
      if (imth.eq.5) then
c
c    0.9 area + 0.1 length
c
         w = 0.1
         p11 = w*1.
         p22 = w*1.
         prgprg = (1.-w)*2.
         return
      end if
c
      if (imth.eq.6) then
c
c    g12*g12
c
         p12 = 2.*g12
         p12p12 = 2.D0
         return
      end if
c
      if (imth.eq.7) then
c
c    Liao
c
         p11 = 2.*g11
         p12 = 4.*g12
         p22 = 2.*g22
         p11p11 = 2.D0
         p12p12 = 4.D0
         p22p22 = 2.D0
         return
      end if
c
      if (imth.eq.8) then
c
c    mod Liao
c
         ca = (g11+g22)/rg
         p11 = 2.*ca/rg
         p22 = p11
         p11p11 = 2./g
         p11p22 = p11p11
         p22p22 = p11p11
         p11prg = -4.*ca/g
         p22prg = p11prg
         prgprg = 6.*ca*ca/g
         return
      end if
c
      if (imth.eq.9) then
c
c    scaled-laplace
c
         r = sqrt(g22/g11)
         p11 = r/2.D0
         p22 = 0.5D0/r
         p11p11 = -0.25D0*r/g11
         p11p22 = 0.25D0/sqrt(g11*g22)
         p22p22 = -0.25D0/r/g22
         return
      end if
c
      if (imth.eq.10) then
c
c    g12*g12/g11/g22
c
         p11 = -g12*g12/g11/g11/g22
         p12 = 2.*g12/g11/g22
         p22 = -g12*g12/g11/g22/g22
         p11p11 = 2.*g12*g12/g11/g11/g11/g22
         p11p12 = -2.*g12/g11/g11/g22
         p11p22 = g12*g12/g11/g11/g22/g22
         p12p12 = 2./g11/g22
         p12p22 = -2.*g12/g11/g22/g22
         p22p22 = 2.*g12*g12/g11/g22/g22/g22
         return
      end if
c
      if (imth.eq.11) then
c
c    g11+g22 - rg*(r + 1/r)
c
         r = sqrt(g22/g11)
         s = 1.D0/r
         u = sqrt(g/g11/g22)
         p11 = 1.D0-u + 0.5D0*u*(1.D0+r*r)
         p22 = 1.D0-u + 0.5D0*u*(1.D0+s*s)
         p11p11 = (g11-3.D0*g22)*u/4.D0/g11/g11
         p11p22 = 0.25D0*u*(1.D0/g11 + 1.D0/g22)
         p22p22 = (g22-3.*g11)*u/4.D0/g22/g22
         p11prg = (r-s)/2.D0/g11
         p22prg = (s-r)/2.D0/g22
         return
      end if
c
      if (imth.eq.12) then
c
c     AO**2
c
         p11 = 2.*g11*g22*g22
         p22 = 2.*g22*g11*g11
         p11p11 = 2.*g22*g22
         p11p22 = 4.*g11*g22
         p22p22 = 2.*g11*g11
      end if
c
      if (imth.eq.13) then
c
cc       Guest: g + g12*g12/g11/g22  (modified AO)
cc
c         p11 = g22 - g12*g12/g11/g11/g22
c         p12 = 2.*g12*(1./g11/g22 -1.)
c         p22 = g11 - g12*g12/g22/g22/g11
c         p11p11 = 2.*g12*g12/g11/g11/g11/g22
c         p11p12 = -2.*g12/g11/g11/g22
c         p11p22 = 1. + g12*g12/g11/g11/g22/g22
c         p12p12 = 2.*(1./g11/g22 - 1.)
c         p22p22 = 2.*g12*g12/g22/g22/g22/g11
c         p12p22 = -2.*g12/g22/g22/g11
c
c        Guest:  ( sqrt( g11 ) + sqrt( g22 ) )^2
c
          p11 = 1.D0 + sqrt( g22/g11 )
          p22 = 1.D0 + sqrt( g11/g22 )
          p11p11 = -0.5D0*sqrt( g22/g11 )/g11
          p11p22 = 0.5D0/sqrt( g11*g22 )
          p22p22 = -0.5D0*sqrt( g11/g22 )/g22
c
      end if
c
      return
      end
