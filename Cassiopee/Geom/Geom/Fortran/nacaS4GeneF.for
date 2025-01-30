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
c Naca Serie 4 airfoils by Dr Moens
c
      SUBROUTINE k6nacas4g(im, ip, it, sharpte, npt, x, y, z, xl)
C
      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E im, ip, it ! digit
      INTEGER_E sharpte    ! sharp trailing edge
      INTEGER_E npt        ! number of points (input and output)
C_OUT
      REAL_E x(1:2*npt)
      REAL_E y(1:2*npt)
      REAL_E z(1:2*npt)
C_LOCAL
      INTEGER_E i, shift
      REAL_E a5, coef, rap, xx, y_thick, yy
      REAL_E teta, y_camb
      REAL_E xm            ! valeur max de la ligne de cambrure (en % corde)
      REAL_E xp            ! position du maximum de cambrure (en 10 eme de corde)
      REAL_E xth
      REAL_E xl(1:npt)

C==============================================================================

c Passage digits -> grandeurs
      xm = im*0.010
      xp = ip*0.100
      xth = it*0.010

      if (sharpte.EQ.1) then
        shift = 1
      else
        shift = 0
      endif

c ***************************
c Evolution en X : loi imposee en sin**3
c
      coef = 0.92*acos(-1.)/2.
      rap = 1./sin(coef)**3
      xl(1) = 0.0
      do i=2, npt-1
         xl(i) = rap*sin(coef*real(i-1)/real(npt-1))**3
      enddo
      xl(npt) = 1.00

c ***************************
      a5 = -0.1015
      if (sharpte.EQ.1) a5 = -0.1036
      do i=1, npt
         xx = xl(i)
         y_thick = 0.2969*sqrt(xx)-0.1260*xx-0.3516*(xx**2)
     &        +0.2843*(xx**3)+a5*(xx**4)
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

c Pour le bord de fuite, on impose x/c=1 et l'epaisseur
c
      yy = (0.1036+a5)*xth/0.20

      x(npt+1-shift) = 1.
      y(npt+1-shift) = yy
      z(npt+1-shift) = 0.
      x(npt) = 1.0
      y(npt) = -yy
      z(npt) = 0.

c        write (fic_out,'(a8,a4)') airfoil_name,'.dat'
c        open (44,file=fic_out)
c        write (44,'(2a)') 'Title=', airfoil_name
c        write (44,'(a)') 'Variables= "x/c", "z/c" '
c        write (44,'(a)') 'Zone'
c        do i=npt,1,-1
c         write (44,'(2(1x,e13.6))') x_low(i),y_low(i)
c        enddo
c        do i=2,npt
c         write (44,'(2(1x,e13.6))') x_up(i),y_up(i)
c        enddo 
c        close(44)

      end
