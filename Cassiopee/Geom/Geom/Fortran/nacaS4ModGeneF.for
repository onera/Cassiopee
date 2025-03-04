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
c Naca Serie 4 modified airfoils by Dr Moens
c
      SUBROUTINE k6nacas4m(im, ip, ith, it, ii, sharpte, npt, x, y, z, 
     & xl)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E im, ip, ith, it, ii, sharpte, npt

C_OUT
      REAL_E x(2*npt),y(2*npt),z(2*npt)

C_LOCAL
      REAL_E xl(npt)
      REAL_E teta, y_camb
      REAL_E xm, xp, xth, xt
      REAL_E coef, rap, d0, d1, d2, d3
      REAL_E a0, a1, a2, a3
      REAL_E yy, xx, xkhile, y_thick, rho1
      INTEGER_E i, shift
C==============================================================================

      if (sharpte.EQ.1) then
        shift = 1
      else
        shift = 0
      endif

        xm = im*0.010
        xp = ip*0.100
        xth = ith*0.010
        xt = it*0.10 
c
c Evolution en X : loi imposee en sin**3
c
        coef = 0.92*acos(-1.)/2.
        rap = 1./sin(coef)**3
        xl(1) = 0.0
        do i = 2, npt-1
         xl(i) = rap*sin(coef*real(i-1)/real(npt-1))**3
        enddo
        xl(npt) = 1.00
c
c Serie 4 Modifiee : Uniquement la loi d'epaisseur
c (cambrure identique a la Serie 4)
c Coefficients (from Mason)
c
        d0=0.002
        if (sharpte.EQ.1) d0=0.0

        d1 = (2.24-5.42*xt+12.3*xt**2)/(10.0*(1.0-0.878*xt))
        d2 = (0.294-2.0*(1.0-xt)*d1)/(1-xt)**2
        d3 = (-0.196+(1.0-xt)*d1)/(1-xt)**3
        xkhile = float(ii)/6.
        if (ii.eq.9) xkhile = 10.3933
        a0 = 0.296904*xkhile
        rho1 = 0.2* (1-xt)**2/(0.588-2.*d1*(1.0-xt)) 
        a1 = 0.3/xt -15.0*a0/(8*sqrt(xt)) - xt/(10.0*rho1)
        a2 = -0.3/(xt**2) + 5.0*a0/(4.0*xt**1.5) +1.0/(5.0*rho1)
        a3 = 0.1/(xt**3) - 0.375*a0/(xt**2.5) -1.0/(10.0*rho1*xt)
c
c Generation profil
c
        do i=1, npt
         xx=xl(i)
         if (xx.lt.xt) then
           y_thick = a0*sqrt(xx) +a1*xx + a2*xx**2 + a3 *xx**3
         else
           y_thick = d0 + d1*(1-xx) + d2*(1-xx)**2 +d3*(1-xx)**3
         endif
         y_thick = y_thick*xth/0.20

         if (xx.lt.xp) then
           if (xp.ne.0.) then
            y_camb = (2.*xp-xx)*xm*xx/(xp**2)
            teta = (xp-xx)*2.*xm/(xp**2)
           else
            y_camb = 0.
            teta = 0.
           endif
           teta = atan(teta)
           x(2*npt-i+1-shift) = xx-y_thick*sin(teta)
           y(2*npt-i+1-shift) = y_camb+y_thick*cos(teta)
           z(2*npt-i+1-shift) = 0.
           x(i) = xx+y_thick*sin(teta)
           y(i) = y_camb-y_thick*cos(teta)
           z(i) = 0.
         else
           teta = 0.
           if (xp.ne.1.) then
            y_camb = (1.+xx-2.*xp)*xm*(1.-xx)/(1.-xp)**2
            teta = (xp-xx)*2.*xm/(1.-xp)**2
           else
            y_camb = 0.
            teta = 0.
           endif
           teta = atan(teta)
           x(2*npt-i+1-shift) = xx-y_thick*sin(teta)
           y(2*npt-i+1-shift) = y_camb+y_thick*cos(teta)
           z(2*npt-i+1-shift) = 0.
           x(i) = xx+y_thick*sin(teta)
           y(i) = y_camb-y_thick*cos(teta)
           z(i) = 0.
         endif

        enddo
c
        yy = d0*xth/0.20

        x(npt+1-shift) = 1.
        y(npt+1-shift) = yy
        z(npt+1-shift) = 0.
        x(npt) = 1.0
        y(npt) = -yy
        z(npt) = 0.
        end
