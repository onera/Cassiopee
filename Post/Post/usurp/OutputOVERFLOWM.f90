! ================================================================
  module OutputOVERFLOWM
! ================================================================

  implicit none

  private
  public :: OutputOVERFLOW

  contains

! ====================================================================
  subroutine OutputOVERFLOW()
! ====================================================================
! Create a text output of the integrated force and moment data.
! The format is intended to replicate the output format in OVERINT and
! is adapted from the WRSUFM and WRCOFM procedures.

  use FomoArrays  ,only: cfpx,cfpy,cfpz,cmpx,cmpy,cmpz
  use FomoArrays  ,only: cfvx,cfvy,cfvz,cmvx,cmvy,cmvz
  use FomoArrays  ,only: cfmx,cfmy,cfmz,cmmx,cmmy,cmmz
  use FomoArrays  ,only: ctax,ctay,ctaz,weta,cfr
  use GroupInfo   ,only: num_groups,group_size
  use GroupInfo   ,only: group_name,group_data,group_iref
  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use PatchInfo   ,only: num_components,refc
  use UserInput   ,only: cpufmt,showcpu,qinf,refa
  use UserInput   ,only: rhoinf,beta,alpha,vminf,refl,refq

  implicit none

  real(kind=rd),parameter :: zero = 0.0_rd
  real(kind=rd)           :: conv,cosa,cosb,sina,sinb

  real(kind=rd)           :: cxps,cyps,czps
  real(kind=rd)           :: cxvs,cyvs,czvs
  real(kind=rd)           :: cxms,cyms,czms
  real(kind=rd)           :: clps,cdps,csps
  real(kind=rd)           :: clvs,cdvs,csvs
  real(kind=rd)           :: clms,cdms,csms,cfrs
  real(kind=rd)           :: cmxps,cmyps,cmzps
  real(kind=rd)           :: cmxvs,cmyvs,cmzvs
  real(kind=rd)           :: cmxms,cmyms,cmzms
  real(kind=rd)           :: cmrs,cmps,cmys

  real(kind=rd)           :: cxpc,cypc,czpc
  real(kind=rd)           :: cxvc,cyvc,czvc
  real(kind=rd)           :: cxmc,cymc,czmc
  real(kind=rd)           :: clpc,cdpc,cspc
  real(kind=rd)           :: clvc,cdvc,csvc
  real(kind=rd)           :: clmc,cdmc,csmc,cfrc
  real(kind=rd)           :: cmxpc,cmypc,cmzpc
  real(kind=rd)           :: cmxvc,cmyvc,cmzvc
  real(kind=rd)           :: cmxmc,cmymc,cmzmc
  real(kind=rd)           :: cmrc,cmpc,cmyc

  real(kind=rd)           :: taxc,tayc,tazc,sac
  real(kind=rd)           :: rvra,qra,qralx,motype,qraly,qralz
  real(kind=rd)           :: hax,hay,haz,hm,cmhs,cmhc
  real(kind=rd)           :: xmc,ymc,zmc
  real(kind=rd)           :: tfx,tfy,tfz

  real(kind=rs)           :: time1,time2

  integer                 :: firstg,i,ic,icf,ig,isc,nr


  continue


  call my_cpu_time(time1)

  conv = acos(-1.0_rd)/180.0_rd
  cosa = cos(conv*alpha)
  sina = sin(conv*alpha)
  cosb = cos(conv*beta )
  sinb = sin(conv*beta )

  do ic = 1,num_components

!   non-dimensionalize the forces and moments
    nr     = refc(ic)
    refl   = refq(1,nr)
    refa   = refq(2,nr)
    xmc    = refq(3,nr)
    ymc    = refq(4,nr)
    zmc    = refq(5,nr)
    motype = refq(9,nr)

    if (refl == 0.0) refl = 1.0
    if (refa == 0.0) then

!     find the first group to which this surface belongs
!     (note that negative values in group_data will not register)
      firstg = 0
      group_loop: do ig = 1,num_groups
        do icf = 1,group_size(ig)
          if (group_data(icf,ig) == ic) then
            firstg = ig
            exit group_loop
          end if
        end do
      end do group_loop

!     if the component does not belong to any group then use
!     the component area for refa else use the group area.
      if (firstg == 0) then
        refa = weta(ic)
      else
        refa = 0.0
        do i = 1,group_size(firstg)
          icf = group_data(i,firstg)
          if (icf > 0) refa = refa + weta(icf)
        end do
      end if
    end if

    rvra  = rhoinf * vminf * refa
    qra   = qinf * refa
    qralx = qinf * refa * refl

    if (motype == 2.0) then
      qraly = qinf * refa * refq(6,nr)
      qralz = qinf * refa * refq(7,nr)
    else
      qraly = qralx
      qralz = qralx
    end if

    clps = -cfpx(ic)*sina + cfpz(ic)*cosa
    cdps = (cfpx(ic)*cosa + cfpz(ic)*sina)*cosb - cfpy(ic)*sinb
    csps = (cfpx(ic)*cosa + cfpz(ic)*sina)*sinb + cfpy(ic)*cosb

    clvs = -cfvx(ic)*sina + cfvz(ic)*cosa
    cdvs = (cfvx(ic)*cosa + cfvz(ic)*sina)*cosb - cfvy(ic)*sinb
    csvs = (cfvx(ic)*cosa + cfvz(ic)*sina)*sinb + cfvy(ic)*cosb

    clms = -cfmx(ic)*sina + cfmz(ic)*cosa
    cdms = (cfmx(ic)*cosa + cfmz(ic)*sina)*cosb - cfmy(ic)*sinb
    csms = (cfmx(ic)*cosa + cfmz(ic)*sina)*sinb + cfmy(ic)*cosb

    cmrs = cmpx(ic) + cmvx(ic) + cmmx(ic)
    cmps = cmpy(ic) + cmvy(ic) + cmmy(ic)
    cmys = cmpz(ic) + cmvz(ic) + cmmz(ic)

    cmxps = cmpx(ic) - ( (ymc*cfpz(ic)) - (zmc*cfpy(ic)) )
    cmyps = cmpy(ic) - ( (zmc*cfpx(ic)) - (xmc*cfpz(ic)) )
    cmzps = cmpz(ic) - ( (xmc*cfpy(ic)) - (ymc*cfpx(ic)) )
    cmxvs = cmvx(ic) - ( (ymc*cfvz(ic)) - (zmc*cfvy(ic)) )
    cmyvs = cmvy(ic) - ( (zmc*cfvx(ic)) - (xmc*cfvz(ic)) )
    cmzvs = cmvz(ic) - ( (xmc*cfvy(ic)) - (ymc*cfvx(ic)) )
    cmxms = cmmx(ic) - ( (ymc*cfmz(ic)) - (zmc*cfmy(ic)) )
    cmyms = cmmy(ic) - ( (zmc*cfmx(ic)) - (xmc*cfmz(ic)) )
    cmzms = cmmz(ic) - ( (xmc*cfmy(ic)) - (ymc*cfmx(ic)) )

    tfx = cfpx(ic) + cfvx(ic) + cfmx(ic)
    tfy = cfpy(ic) + cfvy(ic) + cfmy(ic)
    tfz = cfpz(ic) + cfvz(ic) + cfmz(ic)

    cmrs = cmrs - ( (ymc*tfz) - (zmc*tfy) )
    cmps = cmps - ( (zmc*tfx) - (xmc*tfz) )
    cmys = cmys - ( (xmc*tfy) - (ymc*tfx) )

!   non-dimensionalize the forces
    cxps = cfpx(ic) / qra
    cyps = cfpy(ic) / qra
    czps = cfpz(ic) / qra
    cxvs = cfvx(ic) / qra
    cyvs = cfvy(ic) / qra
    czvs = cfvz(ic) / qra
    cxms = cfmx(ic) / qra
    cyms = cfmy(ic) / qra
    czms = cfmz(ic) / qra
    clps = clps     / qra
    cdps = cdps     / qra
    csps = csps     / qra
    clvs = clvs     / qra
    cdvs = cdvs     / qra
    csvs = csvs     / qra
    clms = clms     / qra
    cdms = cdms     / qra
    csms = csms     / qra

!   non-dimensionalize the moments
    cmxps = cmxps / qralx
    cmyps = cmyps / qralx
    cmzps = cmzps / qralx
    cmxvs = cmxvs / qralx
    cmyvs = cmyvs / qralx
    cmzvs = cmzvs / qralx
    cmxms = cmxms / qralx
    cmyms = cmyms / qralx
    cmzms = cmzms / qralx

!   non-dimensionalize the rpy moments
    cmrs = cmrs / qralx
    cmps = cmps / qraly
    cmys = cmys / qralz

!   non-dimensionalize the mass flux
    cfrs = cfr(ic) /rvra

!=============== wrsufm
    write(*,380)
    write(*,390) ic
    write(*,400)
    write(*,401)
    write(*,402)
    write(*,403) ctax(ic),ctay(ic),ctaz(ic),weta(ic)
    write(*,410)
    write(*,412)
    write(*,414) clps,cdps,csps
    write(*,416) clvs,cdvs,csvs
    if ((clms.ne.zero).or.(cdms.ne.zero).or.(csms.ne.zero)) then
      write(*,418) clms,cdms,csms
    end if
    write(*,420) clps+clvs+clms, cdps+cdvs+cdms, csps+csvs+csms

    if ((motype == 1.0).or.(motype == 2.0)) then
      write(*,430)
      write(*,412)
      write(*,414) cmxps,cmyps,cmzps
      write(*,416) cmxvs,cmyvs,cmzvs
      if ( (cmxms.ne.zero).or.(cmyms.ne.zero).or. &
           (cmzms.ne.zero) ) then
        write(*,418) cmxms,cmyms,cmzms
      end if
      write(*,420) cmrs,cmps,cmys
    else if (motype == 3.0) then
      hax = refq(6,nr) - refq(3,nr)
      hay = refq(7,nr) - refq(4,nr)
      haz = refq(8,nr) - refq(5,nr)
      hm  = sqrt(hax*hax + hay*hay + haz*haz)
      cmhs = (hax*cmrs + hay*cmps + haz*cmys)/hm
      write(*,433) cmhs
    end if
    if (cfrs.ne.zero) then
      write(*,440) cfrs
      write(*,441) -cfrs*2
    end if

  end do ! ic = 1,num_components


  do ig = 1,num_groups

    nr = group_iref(ig)
    motype = refq(9,nr)
    refl   = refq(1,nr)
    refa   = refq(2,nr)
    xmc    = refq(3,nr)
    ymc    = refq(4,nr)
    zmc    = refq(5,nr)

    if (refa == 0.0) then
      do i = 1,group_size(ig)
        ic   = abs(group_data(i,ig))
        isc  = group_data(i,ig) / ic
        refa = refa + weta(ic)*isc
      end do
    end if

    rvra  = rhoinf * vminf * refa
    qra   = qinf * refa
    qralx = qinf * refa * refl

    if (motype == 2.0) then
      qraly = qinf * refa * refq(6,nr)
      qralz = qinf * refa * refq(7,nr)
    else
      qraly = qralx
      qralz = qralx
    end if

    taxc = 0.0_rd
    tayc = 0.0_rd
    tazc = 0.0_rd
    sac  = 0.0_rd

    cxpc = 0.0_rd
    cypc = 0.0_rd
    czpc = 0.0_rd

    cxvc = 0.0_rd
    cyvc = 0.0_rd
    czvc = 0.0_rd

    cxmc = 0.0_rd
    cymc = 0.0_rd
    czmc = 0.0_rd

    cmxpc = 0.0_rd
    cmypc = 0.0_rd
    cmzpc = 0.0_rd

    cmxvc = 0.0_rd
    cmyvc = 0.0_rd
    cmzvc = 0.0_rd

    cmxmc = 0.0_rd
    cmymc = 0.0_rd
    cmzmc = 0.0_rd

    cfrc  = 0.0_rd

    do i = 1,group_size(ig)

      ic   = abs(group_data(i,ig))
      isc  = group_data(i,ig) / ic

      tfx  = cfpx(ic) + cfvx(ic) + cfmx(ic)
      tfy  = cfpy(ic) + cfvy(ic) + cfmy(ic)
      tfz  = cfpz(ic) + cfvz(ic) + cfmz(ic)

      taxc = taxc + ctax(ic)*isc
      tayc = tayc + ctay(ic)*isc
      tazc = tazc + ctaz(ic)*isc

      sac  = sac  + weta(ic)*isc

      cxpc = cxpc + cfpx(ic)*isc
      cypc = cypc + cfpy(ic)*isc
      czpc = czpc + cfpz(ic)*isc

      cxvc = cxvc + cfvx(ic)*isc
      cyvc = cyvc + cfvy(ic)*isc
      czvc = czvc + cfvz(ic)*isc

      cxmc = cxmc + cfmx(ic)*isc
      cymc = cymc + cfmy(ic)*isc
      czmc = czmc + cfmz(ic)*isc

      cmxps = cmpx(ic) - ( (ymc*cfpz(ic)) - (zmc*cfpy(ic)) )
      cmyps = cmpy(ic) - ( (zmc*cfpx(ic)) - (xmc*cfpz(ic)) )
      cmzps = cmpz(ic) - ( (xmc*cfpy(ic)) - (ymc*cfpx(ic)) )
      cmxvs = cmvx(ic) - ( (ymc*cfvz(ic)) - (zmc*cfvy(ic)) )
      cmyvs = cmvy(ic) - ( (zmc*cfvx(ic)) - (xmc*cfvz(ic)) )
      cmzvs = cmvz(ic) - ( (xmc*cfvy(ic)) - (ymc*cfvx(ic)) )
      cmxms = cmmx(ic) - ( (ymc*cfmz(ic)) - (zmc*cfmy(ic)) )
      cmyms = cmmy(ic) - ( (zmc*cfmx(ic)) - (xmc*cfmz(ic)) )
      cmzms = cmmz(ic) - ( (xmc*cfmy(ic)) - (ymc*cfmx(ic)) )

      cmxpc = cmxpc + cmxps*isc
      cmypc = cmypc + cmyps*isc
      cmzpc = cmzpc + cmzps*isc

      cmxvc = cmxvc + cmxvs*isc
      cmyvc = cmyvc + cmyvs*isc
      cmzvc = cmzvc + cmzvs*isc

      cmxmc = cmxmc + cmxms*isc
      cmymc = cmymc + cmyms*isc
      cmzmc = cmzmc + cmzms*isc

!     mass flux
      cfrc = cfrc + cfr(ic)*isc

    end do

    clpc = -cxpc*sina + czpc*cosa
    cdpc = (cxpc*cosa + czpc*sina)*cosb - cypc*sinb
    cspc = (cxpc*cosa + czpc*sina)*sinb + cypc*cosb

    clvc = -cxvc*sina + czvc*cosa
    cdvc = (cxvc*cosa + czvc*sina)*cosb - cyvc*sinb
    csvc = (cxvc*cosa + czvc*sina)*sinb + cyvc*cosb

    clmc = -cxmc*sina + czmc*cosa
    cdmc = (cxmc*cosa + czmc*sina)*cosb - cymc*sinb
    csmc = (cxmc*cosa + czmc*sina)*sinb + cymc*cosb

    cmrc = cmxpc + cmxvc + cmxmc
    cmpc = cmypc + cmyvc + cmymc
    cmyc = cmzpc + cmzvc + cmzmc

!   non-dimensionalize the forces
    cxpc = cxpc / qra
    cypc = cypc / qra
    czpc = czpc / qra
    cxvc = cxvc / qra
    cyvc = cyvc / qra
    czvc = czvc / qra
    cxmc = cxmc / qra
    cymc = cymc / qra
    czmc = czmc / qra
    clpc = clpc / qra
    cdpc = cdpc / qra
    cspc = cspc / qra
    clvc = clvc / qra
    cdvc = cdvc / qra
    csvc = csvc / qra
    clmc = clmc / qra
    cdmc = cdmc / qra
    csmc = csmc / qra

!   non-dimensionalize the moments
    cmxpc = cmxpc / qralx
    cmypc = cmypc / qralx
    cmzpc = cmzpc / qralx
    cmxvc = cmxvc / qralx
    cmyvc = cmyvc / qralx
    cmzvc = cmzvc / qralx
    cmxmc = cmxmc / qralx
    cmymc = cmymc / qralx
    cmzmc = cmzmc / qralx

!   non-dimensionalize the rpy moments
    cmrc = cmrc / qralx
    cmpc = cmpc / qraly
    cmyc = cmyc / qralz

!   non-dimensionalize the mass flow rate
    cfrc = cfrc / rvra


! ================== wrcofm
    write(*, 80)
    write(*, 90) group_name(ig)
    write(*,100)
    write(*,101)
    write(*,102)
    write(*,103) taxc,tayc,tazc,sac

    write(*,105)
    write(*,112)
    write(*,114) cxpc,cypc,czpc
    write(*,116) cxvc,cyvc,czvc
    if ((cxmc.ne.zero).or.(cymc.ne.zero).or.(czmc.ne.zero)) then
      write(*,118) cxmc,cymc,czmc
    end if
    write(*,120) cxpc+cxvc+cxmc, cypc+cyvc+cymc, czpc+czvc+czmc

    write(*,110)
    write(*,112)
    write(*,114) clpc,cdpc,cspc
    write(*,116) clvc,cdvc,csvc
    if ((clmc.ne.zero).or.(cdmc.ne.zero).or.(csmc.ne.zero)) then
      write(*,118) clmc,cdmc,csmc
    end if
    write(*,120) clpc+clvc+clmc, cdpc+cdvc+cdmc, cspc+csvc+csmc

    if ((motype == 1.0).or.(motype == 2.0)) then
      write(*,130)
      write(*,112)
      write(*,114) cmxpc,cmypc,cmzpc
      write(*,116) cmxvc,cmyvc,cmzvc
      if ( (cmxmc.ne.zero).or.(cmymc.ne.zero).or. &
           (cmzmc.ne.zero) ) then
        write(*,118) cmxmc,cmymc,cmzmc
      end if
      write(*,120) cmrc,cmpc,cmyc
    else if (motype == 3.0) then
      hax = refq(6,nr) - refq(3,nr)
      hay = refq(7,nr) - refq(4,nr)
      haz = refq(8,nr) - refq(5,nr)
      hm  = sqrt(hax*hax + hay*hay + haz*haz)
      cmhc = (hax*cmrc + hay*cmpc + haz*cmyc)/hm
      write(*,133) cmhc
    end if
    if (cfrc.ne.zero) then
      write(*,140) cfrc
      write(*,141) -cfrc*2
    end if

  end do ! ig = 1,num_groups

! ================== wrcofm format statements ======================
  80  format('*****************************', &
             '*****************************')
  90  format('Data for component ',a,/)
 100  format(' Integrated Areas')
 101  format('       X             Y             Z            Total')
 102  format('-----------------------------', &
             '-----------------------------')
 103  format(3e14.6,e15.6,/)
 105  format(' Force coefs          X              Y              Z ')
 110  format(' Force coefs        Lift           Drag            Side')
 112  format('-----------------------------', &
             '------------------------------')
 114  format(' Pressure     ',f13.6,'   ',f13.6,'   ',f13.6)
 116  format(' Viscous      ',f13.6,'   ',f13.6,'   ',f13.6)
 118  format(' Momentum     ',f13.6,'   ',f13.6,'   ',f13.6)
 120  format(' Total        ',f13.6,'   ',f13.6,'   ',f13.6,/)
 130  format(' Moment coefs     Roll (X)       Pitch (Y)       Yaw (Z)')
 132  format('              ',f13.6,'   ',f13.6,'   ',f13.6,/)
 133  format(' Hinge moment coef.   = ',e13.6,/)
 140  format(' Mass flow rate coef. = ',e13.6,/)
 141  format(' Ram drag of this mdot= ',f13.6,/)

!=============== wrsufm format statements ======================
 380  format('*****************************', &
             '*****************************')
 390  format('Data for surface ',i2,/)
 400  format(' Integrated Areas')
 401  format('       X             Y             Z            Total')
 402  format('-----------------------------', &
             '-----------------------------')
 403  format(3e14.6,e15.6,/)
 410  format(' Force coefs        Lift           Drag            Side')
 412  format('-----------------------------', &
             '------------------------------')
 414  format(' Pressure     ',f13.6,'   ',f13.6,'   ',f13.6)
 416  format(' Viscous      ',f13.6,'   ',f13.6,'   ',f13.6)
 418  format(' Momentum     ',f13.6,'   ',f13.6,'   ',f13.6)
 420  format(' Total        ',f13.6,'   ',f13.6,'   ',f13.6,/)
 430  format(' Moment coefs     Roll (X)       Pitch (Y)       Yaw (Z)')
 432  format('              ',f13.6,'   ',f13.6,'   ',f13.6,/)
 433  format(' Hinge moment coef.   = ',e13.6,/)
 440  format(' Mass flow rate coef. = ',e13.6,/)
 441  format(' Ram drag of this mdot= ',f13.6,/)

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in OutputOVERFLOW:',time2-time1
  end if


  return
  end subroutine OutputOVERFLOW

  end module OutputOVERFLOWM
