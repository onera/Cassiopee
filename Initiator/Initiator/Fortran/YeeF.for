C ============================================================================
C THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
C ============================================================================
C File  : Initiator/Initiator/Fortran/YeeF.for
C SVN   : $Rev$ $Date$
C Cksum : Initialization of Visbal vortex in field with density constant
C ============================================================================

      SUBROUTINE k6yee(x0, y0, Gamma, Minf,
     &                    npts,    
     &                    xc, yc, zc,
     &                    u)
C
      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      REAL_E x0, y0           ! vortex center
      REAL_E Gamma            ! vortex intensity
      REAL_E Minf             ! infinite Mach number
      INTEGER_E npts          ! number of centers
      REAL_E xc(0:npts-1)    ! x coord of centers
      REAL_E yc(0:npts-1)    ! y cooord of centers
      REAL_E zc(0:npts-1)    ! z cooord of centers
C_OUT
      REAL_E u(0:npts-1,5)     ! field to be initialized
C_LOCAL  
      
      REAL_E  ro0, ainf
      REAL_E p0,  u0
      INTEGER_E ind
      REAL_E r2
      REAL_E roinf, uinf, pinf
      INTEGER_E i,j,k,m
      REAL_E ta
      REAL_E cos_teta, sin_teta
      REAL_E va, ss, ro, pp, t0
      REAL_E gam,gma,pi,rgp
C==============================================================================

      pi=3.141592653589793238462643
      gma = 1.4
      rgp = 287.53

      roinf = 1.                ! etat d'adimensionnement : etat inf
      pinf = 1./gma
      ainf = 1.                 
      uinf= Minf
      gam = Gamma*ainf          ! intensite du tourbillon

      p0 = pinf
      ro0 = roinf
      u0 = Minf
      t0 = p0 / (ro0 * rgp)

      ! pression et densite en chaque point
      DO ind = 0, npts-1
         r2 = (xc(ind)-x0)**2+(yc(ind)-y0)**2
         cos_teta = -(yc(ind)-y0)
         sin_teta = +(xc(ind)-x0)
         va = (gam/(2.*pi))*EXP((1-r2)/2)
         ta = (1.-((gma-1.)/(2.*gma*rgp*t0))*va**2)
               
         pp = ta**(gma/(gma-1.))*p0
         ro = ta**(1./(gma-1.))*ro0
               
         u(ind,1) = ro
         u(ind,2) = ro*uinf+ro*cos_teta*va
         u(ind,3) = ro*sin_teta*va
         u(ind,4) = 0.
         u(ind,5) = pp/(gma-1.)+0.5*(u(ind,2)**2+u(ind,3)**2)/ro
      ENDDO
      
      RETURN

      END











