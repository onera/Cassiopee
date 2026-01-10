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
C Ensemble de subroutines pour generer une repartition tanhmap 
C  ==========================================================================
      SUBROUTINE k6stretch(t1,t2,sn,nbp,dsm,dsp,ityp,inewt)
c
c       ityp = 1 : one side tanh distribution
c       ityp = 2 : two sides tanh distribution
c	  nbp:	 nombre de points de la distribution
c	  t1:		 coordonnees du premier point
c	  t2:		 coordonnees du dernier point
c	  sn:		 distribution
c	  dsm:	 longeur du premier intervalle
c	  dsp:	 longueur du dernier intervalle (si ityp=2)
c	  inewt =1 : utilisation de newton pour que ??
c       
      IMPLICIT NONE
C       
C==============================================================================
C_IN
      REAL_E t1, t2, dsm, dsp
      INTEGER_E nbp, ityp, inewt
C_OUT
      REAL_E sn(nbp)
C==============================================================================
      IF (ityp.EQ.1) THEN
         CALL k6tanh1s(t1,t2,nbp,sn,dsm,inewt)
      ELSE IF (ityp.EQ.2) THEN
         CALL k6tanh2s(t1,t2,nbp,sn,dsm,dsp,inewt)
      ELSE
         WRITE(*,1000)
         STOP 007
      ENDIF
c
 1000   FORMAT(2x,'non-comprehensible stretching function type')
c
      RETURN
      end
c
c
c ======================================================================
      SUBROUTINE k6tanh1s(t1,t2,jmax,t,ds1,inewt)
c ======================================================================
c
c   subroutine : hyperbolic tangent spacing interpolation                 
c   -----------------------------------------------------                 
c                                                                             
c       t1     : starting coordinate                                          
c       t2     : ending coordinate                                            
c   jmax   : number of point in t vector ( j = 1,2,3, ... , jmax)         
c            this value is given by the user                              
c       t      : output vector containing the tanh spacing distribution       
c   ds1    : the specified lenght of initial interval                     
c            this value is given by the user                              
c ----------------------------------------------------------------------
C
      IMPLICIT NONE
C_IN
      REAL_E t1, t2, ds1
      INTEGER_E inewt, jmax
C_OUT
      REAL_E t(jmax)
C_AUX
      REAL_E b, sn1
      REAL_E dt, eps, bbar, dlt, v, w, pi, x1, f, df, dx1
      REAL_E r, s
      INTEGER_E it, j
c.....
C   common/tanh1/ b,sn1,kmax
c.....
      dt      = t2 - t1
      eps     = 1.0e-11
C   kmax    = jmax
c.....
      sn1     = ds1/dt
c.....
      b       = dt/float(jmax-1)/ds1
c.....
c.....  ** first case b>1.0001 **
c.....
      if (b.ge.1.0001) then
         if (b.lt.2.7829681) then
c.........
c.........  ** initial guess by series expansion **
c.........
            bbar  = b-1.D0
            dlt   = sqrt(6.*(bbar))*(1.-0.15*(bbar)
     &                +0.057321429*(bbar)**2
     &                -0.024907295*(bbar)**3+0.0077424461*(bbar)**4
     &                -0.0010794123*(bbar)**5 )
          else
c.........
c.........  ** initial guess by series expansion **
c.........
            v     = log(b)
            w     = 1./b - 0.02852731
            dlt   = v+(1.+1./v)*log(2.*v)-0.02041793+0.24902722*w
     &              +1.9496443*w**2-2.6294547*w**3+8.56795911*w**4
         end if
c.....
c.....  ** second case b<0.9999 **
c.....
      else if(b.lt.0.9999) then
         if(b.lt.0.26938972) then
c........
c........  ** initial guess by series expansion **
c........
            pi    = acos(0.)*2.
            dlt   = pi*(1.-b+b**2-(1.+pi**2/6.)*b**3+6.794732*b**4
     &               -13.205501*b**5+11.726095*b**6)
         else
c........
c........  ** initial guess by series expansion **
c........
            bbar  = 1.-b
            dlt   = sqrt(6.*(bbar))*(1.+0.15*(bbar)
     &                +0.057321429*(bbar)**2
     &                +0.048774238*(bbar)**3-0.0053337753*(bbar)**4
     &                +0.075845134*(bbar)**5 )
         end if
c.....
c.....  ** third case 0.9999<=b<=1.0001 **
c.....
      else 
         dlt   = b-1.0
      end if
c.....
      if (inewt.eq.1) then
c.....
c.....  *use of Newton_Raphson method for zero search ***
c.....
         x1      = dlt
         f       = 0.0
         df      = 1.0
         do it = 1,20
            call k6fnctanh1s (f, x1, b, sn1, jmax)
            call k6dfnctanh1s(df, x1, b, sn1, jmax)
            dx1    = -f/df
            x1     = x1 + dx1
            if (abs(f).lt.eps.and.abs(dx1).lt.eps) goto 1111
         end do
 1111    continue
         dlt     = x1
       end if
c..... 
c.....  ** spacing distribution **
c.....
        t(1)    = t1
        t(jmax) = t2
c..... 
c.....  ** first case b>1.0001 **
c.....
       if (b.ge.1.0001) then
         do 10 j = 2,jmax-1
            r     = float(j-1)/float(jmax-1)
            s     = 1.0 + tanh(dlt*0.5*(r-1.0))/tanh(0.5*dlt)
            t(j)  = t1 + dt*s
 10      continue 
c.....
c.....  ** second case b<0.9999 **
c.....
       else if(b.lt.0.9999) then
         do 20 j = 2,jmax-1
            r     = float(j-1)/float(jmax-1)
            s     = 1.0 + tan (dlt*0.5*(r-1.0))/tan (0.5*dlt)
            t(j)  = t1 + dt*s
 20      continue
c.....
c.....  ** third case 0.9999<=b<=1.0001 **
c.....
       else 
         do 30 j = 2,jmax-1
           r     = float(j-1)/float(jmax-1)
           s     = r*(1.-0.5*dlt*(1.-r)*(2.-r))
           t(j)  = t1 + dt*s
 30      continue 
       end if
c.....
      return
      end

c
c =====================================================================
      subroutine k6fnctanh1s(f,x1,b,sn1,jmax)
c =====================================================================
c
c   subroutine : give the function for hyperbolic tangent spacing     
c ---------------------------------------------------------------------
c
      IMPLICIT NONE
C_IN
      REAL_E x1
      REAL_E b, sn1
      INTEGER_E jmax
C_OUT
      REAL_E f
C_AUX
      REAL_E r
c.....
      r = 1.0/float(jmax-1)
c.....
      if (b.ge.1.0001) then
         f = 1.0D0 + tanh(x1*0.5D0*(r-1.0D0))/tanh(0.5D0*x1)-sn1
      else if(b.lt.0.9999) then
         f = (1.0D0 +tan(x1*0.5D0*(r-1.0D0))/tan(0.5D0*x1))-sn1
      else
         f = (r*(1.0D0-0.5D0*x1*(1.0D0-r)*(2.0D0-r)))-sn1
      end if
c.....
      return
      end
c
c ======================================================================
      subroutine k6dfnctanh1s(df,x1,b,sn1,jmax)
c ======================================================================
c                                                                            
c   function  : give the derivative for hyperbolic tangent spacing      
c
      IMPLICIT NONE
C
C_IN
      REAL_E x1
      REAL_E b, sn1
      INTEGER_E jmax
C_OUT 
      REAL_E df
C_AUX
      REAL_E r, eps, x1eps

      r = 1./float(jmax-1)
      eps = -1.0e-06
      x1eps = x1+eps
c.....
      if (b.ge.1.0001) then
         df = 0.5D0*(r-1.0D0)/(cosh(x1*0.5D0*(r-1.0D0)))**2/tanh(x1*0.5) 
     &      - 0.5*tanh(x1*0.5*(r-1.0))/(sinh(x1*0.5))**2
      else if(b.lt.0.9999) then
         df = 0.5*(r-1.0)/(cos (x1*0.5*(r-1.0)))**2/tan (x1*0.5) 
     &      - 0.5*tan (x1*0.5*(r-1.0))/((sin (x1*0.5))**2)
      else
         df = -r*0.5*(1.0-r)*(2.0-r)
      end if
c.....  
      return
      end
c
c ======================================================================
      SUBROUTINE k6tanh2s(t1,t2,jmax,t,ds1,ds2,inewt)
c ======================================================================
c                                                                             
c   subroutine  :  2 side tanh interpolation                              
c   ----------------------------------------                              
c                                                                            
c       t1     : starting coordinate                                          
c       t2     : ending coordinate                                            
c   jmax   : number of point in t vector ( j = 1,2,3, ... , jmax )        
c            this value is given by the user                              
c       t      : output vector containing the tanh spacing distribution       
c   ds1    : the specified length of initial interval                     
c            this value is given by the user                              
c   ds2    : the specified lenght of final interval                       
c            this value is given by the user                              
c ---------------------------------------------------------------------
c
      IMPLICIT NONE
C_IN
      REAL_E t1, t2, ds1, ds2
      INTEGER_E jmax, inewt
C_OUT
      REAL_E t(jmax)
      REAL_E b, s1n, s2n
C_AUX
      REAL_E dt, eps, zero2, s1, s2, a, v, w, delta
      REAL_E bbar, pi, x1, x2, f1, f2, a11, a12, a21, a22
      REAL_E omega1, omega2, rdelta, dx1, dx2, x1old, x2old
      REAL_E r, u, s
      INTEGER_E itmax, itmax1, it, j
c.....
      dt    = t2 - t1
      eps   = 1.0e-10
      zero2  = 1.0e-12
C        kmax  = jmax
c.....
      s1    = dt/ds1/float(jmax-1)
      s2    = dt/ds2/float(jmax-1)
      s1n   = ds1/dt
      s2n   = ds2/dt
c.....
      a     = sqrt(s1/s2)
      b     = sqrt(s1*s2)
c.....	WRITE(*,*) '1... a =', a, 'b =', b
c.....
c.....   ** first case b>1.0001 **
c.....
      if(b.ge.1.0001) then
        if(b.lt.2.7829681) then
c........
c........  ** initial guess by series expansion **
c........
        bbar  = b-1.D0
        delta = sqrt(6.*(bbar))*(1.-0.15*(bbar)+0.057321429*(bbar)**2
     &           -0.024907295*(bbar)**3+0.0077424461*(bbar)**4
     &           -0.0010794123*(bbar)**5 )
        else
c........
c........  ** initial guess by series expansion **
c........
         v = log(b)
         w = 1./b - 0.02852731
         delta = v+(1.+1./v)*log(2.*v)-0.02041793+0.24902722*w
     &           +1.9496443*w**2-2.6294547*w**3+8.56795911*w**4

        end if
c.....
c.....   ** second case b<0.9999 **
c.....
      else if(b.lt.0.9999) then
         if(b.lt.0.26938972) then
c........
c........   ** initial guess by series expansion **
c........
           pi    = acos(0.)*2.
           delta = pi*(1.-b+b**2-(1.+pi**2/6.)*b**3+6.794732*b**4
     &           -13.205501*b**5+11.726095*b**6)
         else
c........
c........   ** initial guess by series expansion **
c........
           bbar  = 1.-b
           delta = sqrt(6.*(bbar))*(1.+0.15*(bbar)+0.057321429*(bbar)**2
     &           +0.048774238*(bbar)**3-0.0053337753*(bbar)**4
     &           +0.075845134*(bbar)**5 )
         end if
c.....
c.....	** third case 0.9999<=b<=1.0001 **
c.....
       else 
c.......
c........   ** initial guess by series expansion **
c........
         delta = b-1.0
       end if
c.....	WRITE(*,*) '2... delta =', delta
       if (inewt.eq.1) then
c.....
c.....   ** find root delta and a to obtain exactly s1 and s2 **
c.....   ** use of Newton-Raphson method for zero search      **
c.....
        
         x1    = a
         x2    = delta
         f1    = 0.0
         f2    = 0.0
         a11   = 1.0
         a12   = 0.0
         a21   = 0.0
         a22   = 1.0
c.....	WRITE(*,*) '3... a =', a, 'delta =', delta
c.....
         itmax  = 60
         itmax1 = 20
         do it = 1, itmax
         
          call k6fnctanh(f1,f2,x1,x2,b,s1n,s2n,jmax)
          call k6dfnctanh(a11,a12,a21,a22,x1,x2,b,s1n,s2n,jmax)
c.....
          omega1= float(it)/float(itmax1-1)
          omega2= float(it)/float(itmax1-1)
          omega1= MIN(1.0,omega1)
          omega2= MIN(1.0,omega2)
          delta = (a11*a22 - a21*a12) 
          if ( abs(delta).lt.zero2 ) then
            goto 1111
          end if
c.....	  WRITE(*,*) '4... delta =', delta
          rdelta= 1.0/delta
c.....
          dx1   = -(f1*a22 - f2*a12)*rdelta
          dx2   = -(f2*a11 - f1*a21)*rdelta
          dx1   = omega1*dx1
          dx2   = omega2*dx2
c.....
          x1old = x1
          x2old = x2
          x1    = x1 + dx1
          x2    = x2 + dx2
          x1    = MAX(x1,zero2)
          x2    = MAX(x2,zero2)
c.....
          if ( (abs(dx1).lt.eps.and.abs(dx2).lt.eps)
     &    .and.(abs(f1 ).lt.eps.and.abs(f2 ).lt.eps) ) goto 1111
c.....
         end do
 1111    continue
         a      = x1
         delta  = x2
        end if

c.....
c.....   ** spacing distribution **
c.....
    
        t(1)   = t1      
        t(jmax)= t2
c.....
c.....   ** first case b>1.0001 **
c.....
        if(b.ge.1.0001) then
c.....
c.....    ** spacing distribution **
c.....
          do j = 2,jmax-1
             r    = float(j-1)/float(jmax-1)
             u    = 0.5*(1.+tanh(delta*(r-0.5))/tanh(0.5*delta))
             s    = u/(a+(1.-a)*u)
             t(j) = t1 + dt*s
          end do
c.....
c.....   ** second case b<0.9999 **
c.....
        else if(b.lt.0.9999) then
c.....
c.....    ** spacing distribution **
c.....
           do j = 2,jmax-1
             r    = float(j-1)/float(jmax-1)
             u    = 0.5*(1.+tan(delta*(r-0.5))/tan(0.5*delta))
             s    = u/(a+(1.-a)*u)
             t(j) = t1 + dt*s
           end do
c.....
c.....	** third case 0.9999<=b<=1.0001 **
c.....
         else 
c.....
c.....   ** spacing distribution **
c.....
           do j = 2,jmax-1
             r    = float(j-1)/float(jmax-1)
             u    = r*(1.0+2.0*delta*(r-0.5)*(1.0-r))
             s    = u/(a+(1.-a)*u)
             t(j) = t1 + dt*s
           end do
         end if

      return
      end


c ======================================================================
       subroutine k6fnctanh(f1,f2,x1,x2,b,s1n,s2n,jmax)
c ======================================================================
       IMPLICIT NONE
C
c
C_IN
       REAL_E x1, x2
       REAL_E b, s1n, s2n
       INTEGER_E jmax
C_OUT
       REAL_E f1, f2
C_AUX
       REAL_E r1, r2, u1, u2
c.....
C        common/tanh2/b,s1n,s2n,jmax
c.....
       r1  = float(   1  )/float(jmax-1)
       r2  = float(jmax-2)/float(jmax-1)
c.....
       if(b.ge.1.0001) then
c.....
         u1 = 0.5*(1.0+tanh(x2*(r1-0.5))/tanh(0.5*x2))
         u2 = 0.5*(1.0+tanh(x2*(r2-0.5))/tanh(0.5*x2))
         f1 = (x1 + (1.0-x1)*u1)*s1n - u1
         f2 = (x1 + (1.0-x1)*u2)*s2n - x1*(1.0-u2)
c.....
       else if(b.lt.0.9999) then
c.....
         u1 = 0.5*(1.0+tan (x2*(r1-0.5))/tan (0.5*x2))
         u2 = 0.5*(1.0+tan (x2*(r2-0.5))/tan (0.5*x2))
         f1 = (x1 + (1.0-x1)*u1)*s1n - u1
         f2 = (x1 + (1.0-x1)*u2)*s2n - x1*(1.0-u2)
c.....
       else 
c.....
         u1 = r1*(1.0+2.0*x2*(r1-0.5)*(1.0-r1))
         u2 = r2*(1.0+2.0*x2*(r2-0.5)*(1.0-r2))
         f1 = (x1 + (1.0-x1)*u1)*s1n - u1
         f2 = (x1 + (1.0-x1)*u2)*s2n - x1*(1.0-u2)
c        f1 = 0.0
c        f2 = 0.0
c.....
       end if
c.....
       return
       end

c ======================================================================
       subroutine k6dfnctanh(a11,a12,a21,a22,x1,x2,b,s1n,s2n,jmax)
c ======================================================================
       IMPLICIT NONE
C
c
C_IN
        REAL_E x1, x2, b, s1n, s2n
        INTEGER_E jmax
C_OUT
        REAL_E a11, a12, a21, a22
C_AUX
        REAL_E r1, r2, u1, u2, du1dx2, du2dx2
c.....
C        common/tanh2/b,s1n,s2n,jmax
c.....
        r1      = float(   1  )/float(jmax-1)
        r2      = float(jmax-2)/float(jmax-1)
c.....
        if (b.ge.1.0001) then
c.....
         u1     = 0.5*(1.0+tanh(x2*(r1-0.5))/tanh(0.5*x2))
         u2     = 0.5*(1.0+tanh(x2*(r2-0.5))/tanh(0.5*x2))
         du1dx2 = 0.5*(r1-0.5)/(cosh(x2*(r1-0.5)))**2/tanh(x2*0.5) 
     &          - 0.25*tanh(x2*(r1-0.5))/(sinh(x2*0.5))**2
         du2dx2 = 0.5*(r2-0.5)/(cosh(x2*(r2-0.5)))**2/tanh(x2*0.5) 
     &          - 0.25*tanh(x2*(r2-0.5))/(sinh(x2*0.5))**2
c.....
        else if(b.lt.0.9999) then
c.....
         u1     = 0.5*(1.0+tan (x2*(r1-0.5))/tan (0.5*x2))
         u2     = 0.5*(1.0+tan (x2*(r2-0.5))/tan (0.5*x2))
         du1dx2 = 0.5*(r1-0.5)/(cos (x2*(r1-0.5)))**2/tan (x2*0.5)
     &          - 0.25*tan (x2*(r1-0.5))/(sin (x2*0.5))**2
         du2dx2 = 0.5*(r2-0.5)/(cos (x2*(r2-0.5)))**2/tan (x2*0.5)
     &          - 0.25*tan (x2*(r2-0.5))/(sin (x2*0.5))**2
c.....
        else 
c.....
         u1     = r1*(1.0+2.0*x2*(r1-0.5)*(1.0-r1))
         u2     = r2*(1.0+2.0*x2*(r2-0.5)*(1.0-r2))
         du1dx2 = 2.0*r1*(r1-0.5)*(1.0-r1)
         du2dx2 = 2.0*r2*(r2-0.5)*(1.0-r2)
c        u1     = 1.0/s1n + 1.0
c        u2     = 1.0/(s2n-1.0) + 1.0
c        du1dx2 = 0.0
c        du2dx2 = 1.0/( s2n      - x1*(s2n-1.0))
c.....
        end if
c..... 
        a11     = (1.0-u1)*s1n
        a12     = du1dx2*((s1n-1.0) - x1* s1n     )
        a21     = (s2n-1.0)*(1.0-u2)
        a22     = du2dx2*( s2n      - x1*(s2n-1.0))
c.....
       return
       end
