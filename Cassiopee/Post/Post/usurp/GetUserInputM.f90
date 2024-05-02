! ================================================================
  module GetUserInputM
! ================================================================

  implicit none

  private
  public  :: GetUserInput
  private :: ReadNamelistCFDSHIP
  private :: ReadSubin
  private :: read_over_rel_nml
  private :: read_unclem

  contains

! =========================================================
  subroutine GetUserInput()
! =========================================================
! Determine what CFD environment USURP is working in and
! read any associated data.  Current supported environments
! include:
!
! 1. UNCLE / OVER-REL (if over_rel.nml is found)
! 2. OVERFLOW / FOMOCO (if grid.in is found)
! 3. CFD-SHIP (if cfd_ship.nml is found)
! 4. UNCLE-M (if uncle.inp is found)
! 5. NPHASE (if nphase_elements.dat is found)
! 6. UNCLE / UNCLE-REL (if sub.in is found)
! 7. generic (specified format for USURP; if generic.bc is found)
! 8. DPLR (if none of the first seven files is found)
!
! The files are queried in the order listed above.

  use IntrType     ,only: rd,rs
  use my_cpu_timeM ,only: my_cpu_time
  use UserInput    ,only: workdir,over_rel,overflow,cfd_ship
  use UserInput    ,only: default_disjoin,default_basis,dplr
  use UserInput    ,only: refa,refl,iwfn,rough,generic,nphase,basis
  use UserInput    ,only: overset,disjoin,unclem
  use UserInput    ,only: showcpu,cpufmt
  use UserInput    ,only: xcg,ycg,zcg,rel,dci_file_name

  implicit none

  real(kind=rs) :: time1,time2
  logical       :: uncle_rel


  continue


! initialization
  call my_cpu_time(time1)


! ascertain the CFD working environment by looking for 
! identifying files

  inquire(file=trim(workdir)//"over_rel.nml",exist=over_rel)
  inquire(file=trim(workdir)//"grid.in",exist=overflow)
  inquire(file=trim(workdir)//"cfd_ship.nml",exist=cfd_ship)
  inquire(file=trim(workdir)//"uncle.inp",exist=unclem)
  inquire(file=trim(workdir)//"nphase_elements.dat",exist=nphase)
  inquire(file=trim(workdir)//"sub.in",exist=uncle_rel)
  inquire(file=trim(workdir)//"generic.bc",exist=generic)



  if (over_rel) then

    call read_over_rel_nml()

    write(*,*)"from over_rel.nml:"
    write(*,*)"xcg,ycg,zcg    = ",xcg,ycg,zcg
    write(*,*)"iwfn,rel,rough = ",iwfn,rel,rough
    write(*,*)"overset        = ",overset

!   set default behaviors for OVER-REL (UNCLE) users
    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

  else if (overflow) then

    overset = .true.

!   set default behaviors for OVERFLOW users
!   if (default_disjoin) disjoin = .true.  !emulate mixsur
    if (default_disjoin) disjoin = .false. !DJV 3/16/06
    if (default_basis)   basis = "patch"

  else if (cfd_ship) then

    call ReadNameListCFDSHIP()

!   set default behaviors for CFD-SHIP users
    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

  else if (unclem) then

    write(*,*)"USURP found uncle.inp file for UNCLE-M:"

    call read_unclem()

!   set default behaviors for UNCLE-M users
    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

  else if (nphase) then

    write(*,*)
    write(*,*)"using unstructured nphase data in nphase_elements.dat"

!   set default behaviors for NPHASE users
    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

  else if (uncle_rel) then

    call ReadSubin()

    over_rel = .true.

    write(*,*)"from sub.in:"
    write(*,*)"xcg,ycg,zcg    = ",xcg,ycg,zcg
    write(*,*)"iwfn,rel,rough = ",iwfn,rel,rough
    write(*,*)"overset        = ",overset

!   set default behaviors for UNCLE-REL (UNCLE) users
    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

  else if (generic) then

    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

  else

!   if no other code environment has been detected, then assume
!   that this is for DPLR

    dplr = .true. 
    if (default_disjoin) disjoin = .false.
    if (default_basis)   basis = "panel"

! else

!   write(0,*)
!   write(0,*)"ERROR! could not determine cfd work environment."
!   write(0,*)
!   write(0,*)"USURP must be able to find one of the following files:"
!   write(0,*)
!   write(0,*)"over_rel.nml        (for UNCLE / OVER-REL users)"
!   write(0,*)"grid.in             (for OVERFLOW / FOMOCO users)"
!   write(0,*)"cfd_ship.nml        (for CFDSHIP users)"
!   write(0,*)"uncle.inp           (for UNCLE-M users)"
!   write(0,*)"nphase_elements.dat (for NPHASE users)"
!   write(0,*)"sub.in              (for UNCLE / UNCLE-REL users)"
!   write(0,*)"generic.bc          (for other users)"
!   write(0,*)
!   write(0,*)"The files are queried in the order listed above."
!   write(0,*)
!   write(0,*)"Specify --help on the command line for more help."
!   write(0,*)

!   stop

  end if


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in GetUserInput: ',time2-time1
  end if


  return
  end subroutine GetUserInput

! ==============================================================
  subroutine ReadNameListCFDSHIP()
! ==============================================================
! Read the namelist file for CFD-SHIP to get the
! file name prefixes for the grid and boundary condition files.

  use UserInput,only: workdir,fgrid,fnameo

  implicit none

  character(len=30)             :: fnamei
  integer                       :: ierr

  namelist /filenames/ fgrid,fnamei,fnameo


  continue


  open ( unit     = 4,                             &
         file     = trim(workdir)//'cfd_ship.nml', &
         action   = 'read',                        &
         position = 'rewind',                      &
         status   = 'old'                          )

  read(unit=4,nml=filenames,iostat=ierr)
  if (ierr /= 0) then
    write(0,*)'ERROR in namelist filenames: check for undefined variables'
    stop
  end if
  close(4)

  return
  end subroutine ReadNameListCFDSHIP

! ==============================================================
  subroutine ReadSubin()
! ==============================================================
! Read the UNCLE / UNCLE-REL input file to get parameters like
! the Reynolds number and moment origin coordinates.

  use IntrType ,only: rs
  use UserInput,only: rel,xcg,ycg,zcg,iwfn,rough,workdir

  implicit none

  real(kind=rs) :: alph1,alphadot
  real(kind=rs) :: dummy,omega
  real(kind=rs) :: pr_tot,pr_beg,pr_end
  real(kind=rs) :: xref,yref,zref
  real(kind=rs) :: rf,vinf,vplunge

  integer       :: idyn,n,nblades,ndummy
  integer       :: ictype,motion,motyp,backp_it

  namelist /maneuver/alph1,ictype,motion,motyp,rf,vinf,vplunge,&
                     xref,yref,zref,alphadot,xcg,ycg,zcg,&
                     pr_tot,backp_it,pr_beg,pr_end


  continue


! get number of blades and rotation rate

  open(unit=4,file=trim(workdir)//"sub.in",status="old",action="read")
  do n = 1,13
    read(4,*)
  end do
  read(4,*)ndummy,rel,nblades
  read(4,*)
  read(4,*)
  read(4,*)iwfn,rough
  do n = 1,5
    read(4,*)
  end do
  read(4,*)dummy,omega
  read(4,*)
  read(4,*)
  read(4,*)ndummy,idyn

! read namelist for maneuvering parameters
  motyp   = 0
  xref    = 0.
  yref    = 0.
  zref    = 0.
  xcg     = 0.
  ycg     = 0.
  zcg     = 0.

  read (4,nml=maneuver,end=1000)

 1000 continue
  close(4)

  return
  end subroutine ReadSubin

! ==============================================================
  subroutine read_unclem()
! ==============================================================
! Read the UNCLE-M namelist file to get parameters like the 
! DCI filename

  use UserInput,only: overset,dci_file_name,workdir

  implicit none

  character(len=132) :: whole_line


  continue


! find out whether this is overset and what the DCI file is called
  open(unit=3,file=trim(workdir)//"/uncle.inp")
  do
    read(unit=3,fmt="(a132)")whole_line
    if (whole_line == "end") exit
    if (index(whole_line,"OVERSET") /= 0) then
      read(3,*)overset,dci_file_name
      write(*,*)"overset = ",overset
      write(*,*)"dci_file_name = ",trim(dci_file_name)
      exit
    end if
  end do
  close(3)


  return
  end subroutine read_unclem

! =========================================================================
  subroutine read_over_rel_nml()
! =========================================================================
! Read the OVER-REL namelist file to get parameters like the Reynolds
! number and moment origin coordinates.

  use IntrType ,only: rs
  use UserInput,only: iwfn,rel,rough,xcg,ycg,zcg
  use UserInput,only: overset,dci_file_name,workdir

  implicit none

  real(kind=rs) :: xref,yref,zref
  real(kind=rs) :: alpha,chi,dxdt,dydt,dzdt,dtdt,vinf
  real(kind=rs) :: cfl,dtphys,order,time_accuracy
  real(kind=rs) :: cfl_cyc(10),cfl_cfl(10)
  real(kind=rs) :: vturb1,vturb2,order2
  real(kind=rs) :: pr_tot,pr_beg,pr_end
  real(kind=rs) :: hub_parting_lines(2),shroud_parting_lines(2)
  real(kind=rs) :: alph1,rf,vplunge,alphadot
  real(kind=rs) :: omega_wall,sptol,beta

  integer :: motyp,nblades,nblades_wake,nsteps,nsub,restart
  integer :: isgs,ifreq,istate,idyn,ictype,itm,isgs2
  integer :: backp_it,uout_freq,ns2,motion,interp_clip
  integer :: trace_ijkn(4),navg,des_blnd,uout_start,uout_end

  character(len=17) :: force_integration

  logical           :: debug_mode,analytical,menter_prdk
  logical           :: des_mode,spalding_must_converge


! define input namelists
  namelist /operating_conditions/ alpha,chi,nblades, &
    nblades_wake,dxdt,dydt,dzdt,rel,dtdt,vinf
  namelist /iteration_controls/ nsteps,nsub,restart
  namelist /solver_controls/ beta,cfl,order,dtphys,isgs,ifreq, &
                     time_accuracy,istate,idyn,sptol, &
                     cfl_cyc,cfl_cfl,ictype,spalding_must_converge
  namelist /turbulence_model/ itm,vturb1,vturb2,isgs2, &
                     order2,iwfn,rough,analytical,des_mode,des_blnd, &
                     menter_prdk
  namelist /boundary_conditions/  &
                     pr_tot,backp_it,pr_beg,pr_end,overset, &
                     hub_parting_lines,shroud_parting_lines, &
                     dci_file_name,omega_wall,interp_clip
  namelist /output_options/ debug_mode,uout_freq,ns2,navg, &
                     force_integration,trace_ijkn,uout_start, &
                     uout_end
  namelist /maneuvering/alph1,motion,motyp,rf,vplunge, &
                     xref,yref,zref,alphadot,xcg,ycg,zcg
      

  continue


! set default values for namelist parameters

! operating_conditions
  vinf    = 1.0d0
  rel     = 0.0
  nblades = 0
  dtdt    = 0.0
  alpha   = 0.0; chi = 0.0
  dxdt = 0.0; dydt = 0.0; dzdt = 0.0
  nblades_wake = 0

! iteration_controls
  nsteps  = 1
  nsub    = 3000
  restart = 1

! solver_controls
  cfl      = 5.0
  cfl_cfl  = 5.0
  cfl_cyc  = 99999
  dtphys   = 0.
  ictype   = 1
  idyn     = 0
  ifreq    = 1
  isgs     = 5
  istate   = 0
  order    = 3.0
  time_accuracy = 2

! turbulence_model
  itm        = 2
  vturb1     = 0.006
  vturb2     = 0.05
  isgs2      = 5
  order2     = 3.0
  iwfn       = 1
  rough      = 0.0
  analytical = .false.
  menter_prdk = .false.

! boundary_conditions
  overset    = .false.
  dci_file_name = 'gen_dirt.dci'
  pr_tot     = 0.
  backp_it   = 100 ; pr_beg=0.0 ; pr_end=0.0

  hub_parting_lines(1)    = -9999999.
  hub_parting_lines(2)    = +9999999.
  shroud_parting_lines(1) = -9999999.
  shroud_parting_lines(2) = +9999999.

! output_options
  debug_mode = .false.
  ns2        = 100
  uout_freq  = 999999999
  force_integration       = "all"

! maneuvering
  alph1      = 0.
  alphadot   = 0.
  motion     = 99999
  motyp      = 0
  rf         = 0.
  vplunge    = 1.0
  xref = 0.0; yref = 0.0; zref = 0.0
  xcg  = 0.0; ycg  = 0.0; zcg  = 0.0



! read the namelist file
  write(*,*)
  write(*,*)'Reading over_rel.nml ...'
  open (unit=4,file= trim(workdir)//"over_rel.nml", &
        status="old",action="read")
  read (4,nml = operating_conditions, end=10)
  read (4,nml = iteration_controls,   end=10)
  read (4,nml = solver_controls,      end=10)
  read (4,nml = turbulence_model,     end=10)
  read (4,nml = boundary_conditions,  end=10)
  read (4,nml = output_options,       end=10)
  read (4,nml = maneuvering,          end=10)
  10  continue
  close(4)

  return
  end subroutine read_over_rel_nml

  end module GetUserInputM
