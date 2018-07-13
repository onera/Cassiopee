! ================================================================
  module CalculateForcesAndMomentsM
! ================================================================

  implicit none

  private
  public  :: CalculateForcesAndMoments
  private :: Spalding_RGH

  contains

! ====================================================================
  subroutine CalculateForcesAndMoments(p)
! ====================================================================
! Perform the surface integration of the hydrodynamic forces and
! moments (due to pressure and viscous forces) if a solution file
! exists.  Separate coding for UNCLE and OVERFLOW is used in order
! to mimic the native solvers as closely as possible.

  use FomoArrays  ,only: cfvx,cfvy,cfvz,cmvx,cmvy,cmvz
  use FomoArrays  ,only: cfpx,cfpy,cfpz,cmpx,cmpy,cmpz
  use FomoArrays  ,only: cfmx,cfmy,cfmz,cmmx,cmmy,cmmz
  use FomoArrays  ,only: cfr,weta,ctax,ctay,ctaz
  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: num_components,num_patches,num_scalars
  use Types       ,only: panel,surf,scal
  use UserInput   ,only: cpufmt,iwfn,over_rel,overflow
  use UserInput   ,only: plotvel,rel,rough,showcpu,xcg,ycg,zcg

  implicit none

  type(panel),pointer,dimension (:) :: p
  type(panel),pointer :: q

  real(kind=rd) :: dpx,dpy,dpz
  real(kind=rd) :: dvx,dvy,dvz,dvm,dwc
  real(kind=rd) :: dmx,dmy,dmz
  real(kind=rd) :: eddy,pf
  real(kind=rd) :: rew,rgrel,ry
  real(kind=rd) :: tauwall
  real(kind=rd) :: txx,txy,txz,tyy,tyz,tzz
  real(kind=rd) :: ustar,vl,vm1,ypl
  real(kind=rd) :: u3,v3,w3
  real(kind=rd) :: ui,vi,wi
  real(kind=rd) :: xc,yc,zc
  real(kind=rd) :: ftmuj,tftmuj,zuvw,fib,rho,ffr
  real(kind=rd) :: ru,rv,rw
  real(kind=rs) :: time1,time2

  integer       :: ip,ipan,ic,ierr


  continue


  call my_cpu_time(time1)


  allocate(cfpx(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfpx"
    stop
  end if

  allocate(cmpx(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmpx"
    stop
  end if

  allocate(cfpy(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfpy"
    stop
  end if

  allocate(cmpy(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmpy"
    stop
  end if

  allocate(cfpz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfpz"
    stop
  end if

  allocate(cmpz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmpz"
    stop
  end if

  allocate(cfvx(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfvx"
    stop
  end if

  allocate(cmvx(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmvx"
    stop
  end if

  allocate(cfvy(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfvy"
    stop
  end if

  allocate(cmvy(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmvy"
    stop
  end if

  allocate(cfvz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfvz"
    stop
  end if

  allocate(cmvz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmvz"
    stop
  end if

  allocate(ctax(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on ctax"
    stop
  end if

  allocate(ctay(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on ctay"
    stop
  end if

  allocate(ctaz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on ctaz"
    stop
  end if

  allocate(weta(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on weta"
    stop
  end if

  allocate(cfr(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfr"
    stop
  end if

  allocate(cfmx(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfmx"
    stop
  end if

  allocate(cfmy(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfmy"
    stop
  end if

  allocate(cfmz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cfmz"
    stop
  end if

  allocate(cmmx(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmmx"
    stop
  end if

  allocate(cmmy(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmmy"
    stop
  end if

  allocate(cmmz(0:num_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cmmz"
    stop
  end if

  cfpx = 0.0_rd
  cfpy = 0.0_rd
  cfpz = 0.0_rd
  cfvx = 0.0_rd
  cfvy = 0.0_rd
  cfvz = 0.0_rd
  cfmx = 0.0_rd
  cfmy = 0.0_rd
  cfmz = 0.0_rd

  cmpx = 0.0_rd
  cmpy = 0.0_rd
  cmpz = 0.0_rd
  cmvx = 0.0_rd
  cmvy = 0.0_rd
  cmvz = 0.0_rd
  cmmx = 0.0_rd
  cmmy = 0.0_rd
  cmmz = 0.0_rd

  ctax = 0.0_rd
  ctay = 0.0_rd
  ctaz = 0.0_rd
  weta = 0.0_rd
  cfr  = 0.0_rd

  do ip = 1,num_patches

    do ipan = surf(ip)%first_panel,surf(ip)%last_panel

      q => p(ipan)

!     determine component number of this panel
      ic = q%icomp

!     pressure force
      pf  = scal(2,ipan)
      dpx = -pf*q%area_ref*q%xn(1)
      dpy = -pf*q%area_ref*q%xn(2)
      dpz = -pf*q%area_ref*q%xn(3)

!     viscous force
      if (rel > 0.0) then

        if (num_scalars <= 2) then
          write(0,*)
          write(0,*)'ERROR!  viscous force can not be computed'
          write(0,*)'with num_scalars <= 2'
          stop
        end if

        u3    = scal( 3,ipan)
        v3    = scal( 4,ipan)
        w3    = scal( 5,ipan)
        eddy  = scal( 6,ipan)
        vl    = scal( 7,ipan)
            
        if (over_rel) then

          rgrel = 1.0/vl/rel * q%area_ref
           
          txx = 2.0*eddy*(q%area_ref*q%xn(1)*u3)
          tyy = 2.0*eddy*(q%area_ref*q%xn(2)*v3)
          tzz = 2.0*eddy*(q%area_ref*q%xn(3)*w3)
          txy = eddy*q%area_ref*(q%xn(2)*u3+q%xn(1)*v3)
          txz = eddy*q%area_ref*(q%xn(3)*u3+q%xn(1)*w3)
          tyz = eddy*q%area_ref*(q%xn(3)*v3+q%xn(2)*w3)

          dvx = (txx*q%xn(1)+txy*q%xn(2)+txz*q%xn(3))*rgrel
          dvy = (txy*q%xn(1)+tyy*q%xn(2)+tyz*q%xn(3))*rgrel
          dvz = (txz*q%xn(1)+tyz*q%xn(2)+tzz*q%xn(3))*rgrel

        else if (overflow) then

          ftmuj  = eddy/rel/vl
          tftmuj = 2.0*ftmuj
           
          zuvw = (q%xn(1)*u3 + q%xn(2)*v3 + q%xn(3)*w3)/3.0

          txx = tftmuj * q%area_ref * (q%xn(1)*u3 -zuvw)
          tyy = tftmuj * q%area_ref * (q%xn(2)*v3 -zuvw)
          tzz = tftmuj * q%area_ref * (q%xn(3)*w3 -zuvw)
          txy = ftmuj  * q%area_ref * (q%xn(2)*u3+q%xn(1)*v3)
          txz = ftmuj  * q%area_ref * (q%xn(3)*u3+q%xn(1)*w3)
          tyz = ftmuj  * q%area_ref * (q%xn(3)*v3+q%xn(2)*w3)

          dvx = (txx*q%xn(1)+txy*q%xn(2)+txz*q%xn(3))*q%area_ref
          dvy = (txy*q%xn(1)+tyy*q%xn(2)+tyz*q%xn(3))*q%area_ref
          dvz = (txz*q%xn(1)+tyz*q%xn(2)+tzz*q%xn(3))*q%area_ref

          rho = scal( 8,ipan)
          ui  = scal( 9,ipan) 
          vi  = scal(10,ipan) 
          wi  = scal(11,ipan) 
          ru  = scal(12,ipan)
          rv  = scal(13,ipan)
          rw  = scal(14,ipan)

          fib = (ui*q%xn(1) + vi*q%xn(2) + wi*q%xn(3)) * q%area_ref

          ffr = rho * fib
          dmx = -ru * fib
          dmy = -rv * fib
          dmz = -rw * fib

        else

          write(0,*)"ERROR! force integration not available ", &
                    "in this environment."
          stop

        end if

        if (over_rel .and.(iwfn == 1)) then
          if (rough == 0.0) then
            dwc     = scal( 8,ipan)
            ui      = scal( 9,ipan)
            vi      = scal(10,ipan)
            wi      = scal(11,ipan)
            vm1     = sqrt(ui*ui + vi*vi + wi*wi)
            rew     = rel*dwc*vm1
            ry      = 0.0
            ypl     = spalding_rgh(rew,ry)
            ustar   = ypl/dwc/rel
          else
            write(unit=0,fmt=*)"ERROR! no roughness coded yet."
            stop
          end if
          tauwall = ustar**2
          dvm = sqrt(dvx**2 + dvy**2 + dvz**2)
          if (dvm /= 0) then
            dvx = dvx / dvm
            dvy = dvy / dvm
            dvz = dvz / dvm
            dvm = tauwall*q%area_ref
            dvx = dvx * dvm
            dvy = dvy * dvm
            dvz = dvz * dvm
          end if
        end if
      else ! rel <= 0.0
        dvx = 0.0
        dvy = 0.0
        dvz = 0.0
      end if

!     store the skin friction coefficient for plotting
!     (OVERWRITE scal(12,ipan) for overflow!!!)
      if (plotvel) then
        if (q%area_ref /= 0.0) then
          tauwall = sqrt(dvx**2+dvy**2+dvz**2)/q%area_ref
        else
          tauwall = 0.0
        end if
        scal(12,ipan) = 2.0*tauwall
      end if

!     center of panel
      xc  = q%xmid(1)
      yc  = q%xmid(2)
      zc  = q%xmid(3)

!     adjust forces by overlap area ratio
      dpx = dpx * q%ratio
      dpy = dpy * q%ratio
      dpz = dpz * q%ratio
      dvx = dvx * q%ratio
      dvy = dvy * q%ratio
      dvz = dvz * q%ratio
      dmx = dmx * q%ratio
      dmy = dmy * q%ratio
      dmz = dmz * q%ratio

!     adjust mass flow by overlap area ratio
      ffr = ffr * q%ratio

      cfpx(ic) = cfpx(ic) + dpx
      cfpy(ic) = cfpy(ic) + dpy
      cfpz(ic) = cfpz(ic) + dpz

      cfvx(ic) = cfvx(ic) + dvx
      cfvy(ic) = cfvy(ic) + dvy
      cfvz(ic) = cfvz(ic) + dvz

      cfmx(ic) = cfmx(ic) + dmx
      cfmy(ic) = cfmy(ic) + dmy
      cfmz(ic) = cfmz(ic) + dmz

      cmpx(ic) = cmpx(ic) - (dpy*(zc-zcg)-dpz*(yc-ycg))
      cmpy(ic) = cmpy(ic) - (dpz*(xc-xcg)-dpx*(zc-zcg))
      cmpz(ic) = cmpz(ic) - (dpx*(yc-ycg)-dpy*(xc-xcg))

      cmvx(ic) = cmvx(ic) - (dvy*(zc-zcg)-dvz*(yc-ycg))
      cmvy(ic) = cmvy(ic) - (dvz*(xc-xcg)-dvx*(zc-zcg))
      cmvz(ic) = cmvz(ic) - (dvx*(yc-ycg)-dvy*(xc-xcg))

      cmmx(ic) = cmmx(ic) - (dmy*(zc-zcg)-dmz*(yc-ycg))
      cmmy(ic) = cmmy(ic) - (dmz*(xc-xcg)-dmx*(zc-zcg))
      cmmz(ic) = cmmz(ic) - (dmx*(yc-ycg)-dmy*(xc-xcg))

      ctax(ic) = ctax(ic) + (q%area_ref * q%xn(1)) * q%ratio
      ctay(ic) = ctay(ic) + (q%area_ref * q%xn(2)) * q%ratio
      ctaz(ic) = ctaz(ic) + (q%area_ref * q%xn(3)) * q%ratio

      weta(ic) = weta(ic) + q%area_ref * q%ratio
      cfr(ic)  = cfr(ic)  + ffr
    end do
  end do


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CalculateForcesAndMoments:',time2-time1
  end if


  return
  end subroutine CalculateForcesAndMoments

! =======================================================================
  function Spalding_rgh(re,ry) result (yp)
! =======================================================================
! wall function based on spalding's profile from the
! wall to the edge of the log-layer. this version
! includes roughness effects, i.e., 
! (ry = (roughness height)/(dist. from wall))

  use IntrType,only: rd

  implicit none

  real(kind=rd),intent(in) :: re,ry
  real(kind=rd)            :: yp
  real(kind=rd)            :: yp0
  real(kind=rd)            :: f, dfdy, dyp
  real(kind=rd),parameter  :: tol = 1.0d-05
  real(kind=rd),parameter  :: xk = 0.40_rd,bb = 5.0_rd

  integer,parameter        :: maxit = 100
  integer                  :: m


  continue


! initial guess
  yp0 = (0.11451_rd*re + 1.2363_rd)**(0.8772_rd)


! newton-raphson iteration of the inverse of spalding's
! law-of-the-wall for the whole boundary layer
! (including roughness effects)

  yp = yp0
  do m = 1,maxit

    f    = yp - re/yp - (1.0_rd+0.3_rd*ry*yp)*exp(-xk*bb)* &
          (exp(xk*re/yp) - 1.0_rd - (xk*re/yp) &
          - 1.0_rd/2.0_rd*(xk*re/yp)**2  &
          - 1.0_rd/6.0_rd*(xk*re/yp)**3)

    dfdy = 1.0_rd + re/yp**2 - (1.0_rd+0.3_rd*ry*yp)*exp(-xk*bb)* &
          ((xk*re/yp**2)*(1.0_rd - exp(xk*re/yp)) &
          + (xk*re/yp)**2/yp + 1.0_rd/2.0_rd*(xk*re/yp)**3/yp) &
          -(0.3_rd*ry)*exp(-xk*bb)* &
          (exp(xk*re/yp) - 1.0_rd - (xk*re/yp) &
          - 1.0_rd/2.0_rd*(xk*re/yp)**2 &
          - 1.0_rd/6.0_rd*(xk*re/yp)**3)

    dyp  = -f/dfdy
    yp   = yp + dyp

    if (abs(dyp) <= tol) return
  end do

  write(*,*)'error - spalding_rgh iteration did not converge'
  write(*,*)'yp, dyp, yp0 =', yp,dyp,yp0
  yp = yp0

  return
  end function Spalding_rgh

  end module CalculateForcesAndMomentsM
