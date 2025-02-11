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
c Naca Serie 5 airfoils by Dr Moens
c
      SUBROUTINE k6nacas5g(il, ip, iq, it, sharpte, npt, x, y, z, xl)
C
      IMPLICIT NONE
C==============================================================================
C_IN
        INTEGER_E il, ip, iq, it, sharpte, npt

C_OUT
        REAL_E x(2*npt), y(2*npt), z(2*npt)
        
C_LOCAL
        REAL_E xm(5),xk1(5),xf(5)
        REAL_E xl(npt)
        REAL_E teta, y_camb, y_thick
        REAL_E a5, coef, yy, rap, xt, xk2sk1, xx, xxk1

        INTEGER_E i, shift

C==============================================================================

      if (sharpte.EQ.1) then
        shift = 1
      else
        shift = 0
      endif

c
c Evolution en X : loi imposee en sin**3
c
        xt = it*0.010
c
c Comme les digits sont des entiers, autant rentrer
c les valeurs tabulees ...
c
        xf(1) = 0.05
        xf(2) = 0.10
        xf(3) = 0.15
        xf(4) = 0.20
        xf(5) = 0.25
        if (iq.eq.0) then
         xm(1) = 0.0580
         xm(2) = 0.1260
         xm(3) = 0.2025
         xm(4) = 0.2900
         xm(5) = 0.3910
         xk1(1) = 361.4
         xk1(2) = 51.65
         xk1(3) = 15.65
         xk1(4) = 6.643
         xk1(5) = 3.230
        else
         xm(1) = 0.0580
         xm(2) = 0.1300
         xm(3) = 0.2170
         xm(4) = 0.3180
         xm(5) = 0.4410
         xk1(1) = 361.4
         xk1(2) = 51.99
         xk1(3) = 15.793
         xk1(4) = 6.520
         xk1(5) = 3.191
        endif 
c
c Calcul coefficient k1 et k2/k1 (Mason)
c
      xk2sk1 = (3.*(xm(ip)-xf(ip))**2-(xm(ip)**3))/(1.0-xm(ip))**3
c adaptation du coef k1 selon le premier digit (Cz d'adaptation)
      xxk1 = xk1(ip)*float(il)/float(2)

c ***************************
c Evolution en X : loi imposee en sin**3
c
        coef = 0.92*acos(-1.)/2.
        rap = 1./sin(coef)**3
        xl(1) = 0.0
        do i = 2, npt-1
         xl(i) = rap*sin(coef*real(i-1)/real(npt-1))**3
        enddo
        xl(npt) = 1.00
c ***************************
        a5 = -0.1015
        if (sharpte.EQ.1) a5 = -0.1036

        do i=1, npt
         xx = xl(i)
         y_thick = 0.2969*sqrt(xx)-0.1260*xx-0.3516*(xx**2)
     &            +0.2843*(xx**3)+a5*(xx**4)
         y_thick = y_thick*xt/0.20

         if (xx.lt.xm(ip)) then

           if (iq.eq.0) then
           y_camb = xx**3-3.0*xm(ip)*(xx**2)+xx*(3.0-xm(ip))*(xm(ip)**2)
           y_camb = y_camb*xxk1/6.0

           teta = 3.0*xx**2 -6.0*xm(ip)*xx+(3.0-xm(ip))*(xm(ip)**2)
           teta = teta*xxk1/6.0

         else

           y_camb = ((xx-xm(ip))**3-xk2sk1*xx*(1.0-xm(ip))**3
     &            -xx*(xm(ip))**3 + xm(ip)**3) * xxk1/6.0

           teta = xxk1/6.0 *(3.0*(xx-xm(ip))**2-xk2sk1*(1.0-xm(ip))**3
     &          -(xm(ip))**3)
           endif
           teta = atan(teta)

           x(2*npt-i+1-shift) = xx-y_thick*sin(teta)
           y(2*npt-i+1-shift) = y_camb+y_thick*cos(teta)
           z(2*npt-i+1-shift) = 0.
           x(i) = xx+y_thick*sin(teta)
           y(i) = y_camb-y_thick*cos(teta)
           z(i) = 0.
         else

           if (iq.eq.0) then
           y_camb = xxk1/6.0 *(xm(ip)**3)*(1-xx)
           teta = -xxk1/6.0 *(xm(ip)**3)
           else
           y_camb = (xk2sk1*(xx-xm(ip))**3-xk2sk1*xx*(1.0-xm(ip))**3
     &            -xx*(xm(ip))**3+ xm(ip)**3) * xxk1/6.0

           teta = xxk1/6.0 *(3.0*xk2sk1*(xx-xm(ip))**2
     &          -xk2sk1*(1.0-xm(ip))**3-(xm(ip))**3)
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
c Pour le bord de fuite, on impose x/c=1 et l'epaisseur
c
        yy = (0.1036+a5)*xt/0.20

        x(npt+1-shift) = 1.
        y(npt+1-shift) = yy
        z(npt+1-shift) = 0.
        x(npt) = 1.0
        y(npt) = -yy
        z(npt) = 0.

        end
