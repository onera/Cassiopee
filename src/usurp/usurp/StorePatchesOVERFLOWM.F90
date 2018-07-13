! ================================================================
  module StorePatchesOVERFLOWM
! ================================================================

  implicit none

  private
  public  :: StorePatchesOVERFLOW
  private :: InterpolateIBlankFull
  private :: ReadFlow
  private :: ReadGrid
  private :: ReadHeader

  contains


! ===========================================================
  subroutine StorePatchesOVERFLOW()
! ===========================================================
! Read in the grid (from grid.in) and flow solution (from 
! q.save, if it exists) for all integration surfaces, and 
! store all necessary information on a panel-by-panel basis.

  use CalculateVolumesM,only: CalculateVolumes
  use IntrType         ,only: rd,rs
  use my_cpu_timeM     ,only: my_cpu_time
  use NetAlloc         ,only: nalloc
  use PatchInfo        ,only: i1,i2,j1,j2,k1,k2,ibdir,block
  use PatchInfo        ,only: num_panels,num_patches,num_scalars
  use Types            ,only: scal,surf
  use UserInput        ,only: alpha,beta,double,fsmach
  use UserInput        ,only: full_surface,gaminf,ibex
  use UserInput        ,only: ignore_pinf,ignore_solution
  use UserInput        ,only: multiple_grid,qinf,rel,rey,rhoinf
  use UserInput        ,only: showcpu,solution_exists,tinf
  use UserInput        ,only: cpufmt,verbosity,vminf,workdir
  use VertexData       ,only: nqv,num_nodes,qv,xv
  use VolumeArrays     ,only: ibv,qd,volv,xb,yb,zb

  implicit none

  real(kind=rd),allocatable,dimension (:)   :: rgas

  real(kind=rd) :: fsmach1,alpha1,beta1,rey1,gaminf1,tinf1
  real(kind=rd) :: htinf,ht1,ht2,time
  real(kind=rd) :: pinf,rqinf,ttinf,gigm1,gm1
  real(kind=rd) :: cinf
  real(kind=rd) :: rr,u,v,w,e0,v2,ei,qq,p,dx,dy,dz,rl
  real(kind=rd) :: ul,vl,wl,vmu,c2b,c2bp,gamma
  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: utime,ltime,vtime

  integer,allocatable,dimension (:) :: nig,njg,nkg
  integer,allocatable,dimension (:) :: niq,njq,nkq
  integer,dimension (3)             :: i1s,i1e,i2e

  integer :: ierr,ip,ipan,n,nb,nbq,nq,nqc
  integer :: i,j,k,iv,m
  integer :: ibv1,ibv2,ibv3,ibv4,igam
  integer :: iv1,iv2,iv3,iv4
  integer :: dir1,dir2,dir3
  integer :: istep,jstep,kstep
  integer :: iqtyp
  logical :: incomp


  continue


! initialization
  call my_cpu_time(time1)
  utime = 0
  ltime = 0
  vtime = 0
  incomp = .false.


  allocate(xv(3,num_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate of xv failed in StorePatchesOVERFLOW"
    write(0,*)"num_nodes = ",num_nodes
    stop
  end if


! check to see if there is flow data
  inquire(file=trim(workdir)//"q.save",exist=solution_exists)
  if (ignore_solution) solution_exists = .false.

  if (solution_exists) then
    num_scalars = 14
  else
    num_scalars = 1
  end if

  allocate(scal(num_scalars,num_panels)); nalloc = nalloc + 1


! open the input grid file
  open ( unit    = 3,                          &
         file    = trim(workdir)//"grid.in",   &
         form    = "unformatted",              &
         status  = "old",                      &
#ifdef CONVERT_BIG_ENDIAN
         convert = "big_endian" ,              &
#endif
         action  = "read"                      )


! read the grid header
  if (multiple_grid) then
     read(unit=3)nb
  else
     nb = 1
  end if
  allocate(nig(nb),njg(nb),nkg(nb)); nalloc = nalloc + 3
  read(unit=3)(nig(n),njg(n),nkg(n),n=1,nb)


  if (solution_exists) then
    if (verbosity >= 1) write(*,*)"Reading q.save ..."
    open ( unit    = 4,                          &
           file    = trim(workdir)//"q.save",    &
           form    = "unformatted",              &
           status  = "old",                      &
#ifdef CONVERT_BIG_ENDIAN
           convert = "big_endian",               &
#endif
           action  = "read"                      )

    if (multiple_grid) then  
      read(4)nbq
    else
      nbq = 1
    end if
    if (nbq < nb) then
      write(0,*)"ERROR! nbq < nb in StorePatchesOVERFLOW"
      write(0,*)"nb  = ",nb
      write(0,*)"nbq = ",nbq
      stop
    end if
    allocate(niq(nbq),njq(nbq),nkq(nbq)); nalloc = nalloc + 3
    read(4,iostat=ierr)(niq(n),njq(n),nkq(n),n=1,nbq),nq,nqc
    if (ierr /= 0) then
      rewind(4)
      if (multiple_grid) read(4)nbq
      read(4)(niq(n),njq(n),nkq(n),n=1,nbq)
      nq = 5
      nqc = 0
    end if

!   verify that flow field dimensions are equal to grid dimensions
    do n = 1,nb
      if ((niq(n) /= nig(n)).or.(njq(n) /= njg(n)).or.(nkq(n) /= nkg(n))) then
        write(0,*)"ERROR! grid point numbers do not match ", &
                  "in StorePatchesOVERFLOW"
        write(0,*)"n   = ",n
        write(0,*)"grid = ",nig(n),njg(n),nkg(n)
        write(0,*)"flow = ",niq(n),njq(n),nkq(n)
        stop
      end if
    end do
    deallocate(nkq,njq,niq); nalloc = nalloc - 3
    allocate(rgas(max(nqc,2))); nalloc = nalloc + 1
    rgas = 0.0_rd


!   read flow conditions from the q.save header

    call ReadHeader(4,iqtyp,nqc, &
                    fsmach1,alpha1,rey1,time,gaminf1,beta1, &
                    tinf1,htinf,ht1,ht2,rgas,igam)

!   initialize reference quantities

    if (fsmach <= 0.0) then
      if (fsmach1 > 0.0) then
        fsmach = fsmach1
        if (verbosity >= 1) then
          write(*,*)"fsmach (from q.save) is ",real(fsmach)
        end if
      else
        fsmach = 1.0
        incomp = .true.
        if (verbosity >= 1) then
          write(*,*)"fsmach is being set to ",real(fsmach)
        end if
      end if
    else
      if (verbosity >= 1) then
        write(*,*)"fsmach (from stdin) is ",real(fsmach)
      end if
    end if

    if ( (alpha < -360.0) .or. (alpha > 360.0) ) then
      alpha = alpha1
      if (verbosity >= 1) then
        write(*,*)"alpha (from q.save) is ",real(alpha)
      end if
    else
      if (verbosity >= 1) then
        write(*,*)"alpha (from stdin) is ",real(alpha)
      end if
    end if

    if ( (beta < -360.0) .or. (beta > 360.0) ) then
      if ( (beta1 < -360.0) .or. (beta1 > 360.0) ) then
        write(0,*)"ERROR! beta is not in q.save.  specify in input."
        stop
      else
        beta = beta1
        if (verbosity >= 1) then
          write(*,*)"beta (from q.save) is ",real(beta)
        end if
      end if
    else
      if (verbosity >= 1) then
        write(*,*)"beta (from stdin) is ",real(beta)
      end if
    end if

    if (rey < 0.0) then
      rey = rey1
      if (verbosity >=1 ) then
        write(*,*)"rey (from q.save) is ",real(rey)
      end if
    else
      if (verbosity >= 1) then
        write(*,*)"rey (from stdin) is ",real(rey)
      end if
    end if

    rey = rey / fsmach
    if (verbosity >= 1) then
      write(*,*)"scaling rey by fsmach; rey = ",real(rey)
    end if

    if (gaminf <= 0.0) then
      if (gaminf1 <= 0.0) then
        write(0,*)"ERROR! gaminf is not in q.save.  specify in input."
        stop
      else
        gaminf = gaminf1
        if (verbosity >= 1) then
          write(*,*)"gaminf (from q.save) is ",real(gaminf)
        end if
      end if
    else
       if (verbosity >= 1) then
         write(*,*)"gaminf (from stdin) is ",real(gaminf)
       end if
    end if

    if (tinf <= 0.0) then
      if (tinf1 <= 0.0) then
        write(0,*)"ERROR! tinf is not in q.save.  specify in input."
        stop
      else
        tinf = tinf1
        if (verbosity >= 1) then
          write(*,*)"tinf (from q.save) is ",real(tinf)
        end if
      end if
    else
      if (verbosity >= 1) then
        write(*,*)"tinf (from stdin) is ",real(tinf)
      end if
    end if

    if (verbosity >=1) then
      write(*,*)"Other parameters from q.save:"
      write(*,*)"time   = ",real(time)
      if (iqtyp /= 0) then
!       write(*,*)"iqtyp  = ",iqtyp
        write(*,*)"igam   = ",igam
        write(*,*)"htinf  = ",real(htinf)
        write(*,*)"ht1    = ",real(ht1)
        write(*,*)"ht2    = ",real(ht2)
        write(*,*)"rgas   = ",real(rgas)
      end if
    end if

!   dervied quantities

    cinf   = 1.0

    if (incomp.or.ignore_pinf) then
      pinf = 0.0
    else
      pinf   = rhoinf*cinf*cinf/gaminf !freestream static pressure
    end if
    vminf  = fsmach*cinf             !freestream velocity
    qinf   = 0.5*rhoinf*vminf*vminf  !freestream dynamic pressure
    rel    = rey                     !reynolds number
    rqinf  = 1.0/qinf

    if (tinf == 0.0) then
      write(0,*)"ERROR! tinf = 0"
      stop
    end if

    c2b    = 199.0/tinf
    c2bp   = c2b + 1.0

  end if !solution_exists


! allocate the vertex-based solution variables
  if (solution_exists) then
!   nqv = min(5,nq)
    if (nq < 5) then
      write(0,*)"ERROR! Not expecting nq in q.save to be less than 5"
      write(0,*)"nq = ",nq
      stop
    end if
!   (11/13/2006) hardwire nqv to 11 to match triq format from overint, 
!   in which slot 1 is Cp, slots 2-6 are Q1-Q5, slot 7 is laminar 
!   viscosity, slot 8-10 are the zeta derivatives, and slot 11 is the 
!   distance from the surface point to the next point up
    nqv = 11
    allocate(qv(nqv,num_nodes),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! allocate of qv failed in StorePatchesOVERFLOW"
      write(0,*)"num_nodes = ",num_nodes
      stop
    end if
  end if


  do n = 1,nb

    if (solution_exists) then

      read(4) !skip header

      allocate( qd(nig(n),njg(n),nkg(n),nq),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! allocate failed for qd; nq = ",nq
        stop
      end if

    end if !solution_exists


!   allocate space for the grid and i-blank variable
    allocate( xb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
    allocate( yb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
    allocate( zb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
    allocate(ibv(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1


!   read the next grid block
    call my_cpu_time(time3)

    call ReadGrid(3,nig(n),njg(n),nkg(n))
    if (solution_exists) call ReadFlow(4,nig(n),njg(n),nkg(n),nq)

    call my_cpu_time(time4)
    utime = utime + (time4-time3)

    call my_cpu_time(time3)
    if (solution_exists) then

!     calculate the cell volumes
      allocate(volv(nig(n)+1,njg(n)+1,nkg(n)+1)); nalloc = nalloc + 1
      call CalculateVolumes(nig(n),njg(n),nkg(n))

    end if
    call my_cpu_time(time4)
    vtime = vtime + (time4-time3)


!   fill in the grid for each patch
    call my_cpu_time(time3)
    do ip = 1,num_patches
      if (block(ip) == n) then


!       figure out which direction is constant
        dir1 = abs(ibdir(ip))
        dir2 = mod(dir1,3) + 1
        dir3 = mod(dir2,3) + 1


        i1s(1) = i1(ip)
        i1e(1) = i2(ip)
        i1s(2) = j1(ip)
        i1e(2) = j2(ip)
        i1s(3) = k1(ip)
        i1e(3) = k2(ip)

        if (solution_exists.and.(rel > 0.0)) then
          if (((dir1 == 1).and.(nig(n) == 1)).or. &
              ((dir1 == 2).and.(njg(n) == 1)).or. &
              ((dir1 == 3).and.(nkg(n) == 1))) then
             write(0,*)"WARNING! Only one layer of grid points was provided."
             write(0,*)"Viscous effects will be ignored."
             rel = 0.0
          end if
        end if

        if (i1s(dir2) == i1e(dir2)) then
          write(0,*)"ERROR in patch ",ip
          stop
        end if
        if (i1s(dir3) == i1e(dir3)) then
          write(0,*)"ERROR in patch ",ip
          stop
        end if
   
        i1s(dir2) = i1s(dir2) + 1
        i1s(dir3) = i1s(dir3) + 1


        istep = 0
        jstep = 0
        kstep = 0

        if (dir1 == 1) then !i-surface
          istep = -1
          if (i1(ip) == 1) istep = 1
          if (nig(n) == 1) istep = 0
        else if (dir1 == 2) then !j-surface
          jstep = -1
          if (j1(ip) == 1) jstep = 1
          if (njg(n) == 1) jstep = 0
        else if (dir1 == 3) then !k-surface
          kstep = -1
          if (k1(ip) == 1) kstep = 1
          if (nkg(n) == 1) kstep = 0
        end if


!       store the grid points for the patch
        iv = surf(ip)%first_node
        do k = k1(ip),k2(ip)
        do j = j1(ip),j2(ip)
        do i = i1(ip),i2(ip)
          xv(1,iv) = xb(i,j,k)
          xv(2,iv) = yb(i,j,k)
          xv(3,iv) = zb(i,j,k)
          iv = iv + 1
        end do
        end do
        end do


        if (solution_exists) then
!         store the vertex variables for the patch
          iv = surf(ip)%first_node
          do k = k1(ip),k2(ip)
          do j = j1(ip),j2(ip)
          do i = i1(ip),i2(ip)

!           (11/13/2006): Fill qv slots to match overint triq format

            if (qd(i,j,k,1) == 0.0) then
              write(0,*)"ERROR! zero density" 
              write(0,*)"i,j,k,ip = ",i,j,k,ip
              stop
            end if

            rr = 1./qd(i,j,k,1)              ! recip of density
            u  = qd(i,j,k,2)*rr              ! x-velocity
            v  = qd(i,j,k,3)*rr              ! y-velocity
            w  = qd(i,j,k,4)*rr              ! z-velocity
            e0 = qd(i,j,k,5)*rr              ! total energy
            v2 = u**2 + v**2 + w**2          ! velocity mag squared
            ei = e0-0.5*v2                   ! internal energy
            qq = 0.5*qd(i,j,k,1)*v2          ! dynamic pressure

            if (nq < 6) then
              gamma = gaminf
            else
              gamma = qd(i,j,k,6)
            end if

            p  = (gamma-1.)*(qd(i,j,k,5)-qq) ! pressure

            rl = qd(i+istep,j+jstep,k+kstep,1)

            if (rl == 0.0) then
              write(0,*)"ERROR! zero density" 
              write(0,*)"i,j,k,ip = ",i+istep,j+jstep,k+kstep,ip
              stop
            end if

            ul = qd(i+istep,j+jstep,k+kstep,2)/rl !x-velocity, 2nd layer
            vl = qd(i+istep,j+jstep,k+kstep,3)/rl !y-velocity, 2nd layer
            wl = qd(i+istep,j+jstep,k+kstep,4)/rl !z-velocity, 2nd layer

!           sutherlands law for viscosity
            gm1 = gamma - 1.
            gigm1 = gaminf * gm1
            ttinf = gigm1 * ei

            if (incomp) then
              vmu = 1.
            else
              if (ttinf < 0.0) then
                write(0,*)"ERROR! ttinf < 0"
                stop
              end if
              vmu = c2bp * (ttinf * sqrt(ttinf))/(c2b + ttinf)
            end if

            dx = xb(i+istep,j+jstep,k+kstep) - xb(i,j,k)
            dy = yb(i+istep,j+jstep,k+kstep) - yb(i,j,k)
            dz = zb(i+istep,j+jstep,k+kstep) - zb(i,j,k)

            qv(  1,iv)   = (p - pinf) / qinf
            qv(2:6,iv)   = qd(i,j,k,1:5)
            qv(  7,iv)   = vmu
            qv(  8,iv)   = ul - u
            qv(  9,iv)   = vl - v
            qv( 10,iv)   = wl - w
            qv( 11,iv)   = sqrt(dx*dx + dy*dy + dz*dz)

            iv = iv + 1
          end do
          end do
          end do
        end if


        if (solution_exists) then

!         loop over all panels

          ipan = surf(ip)%first_panel

          m = min(nq,6)

          do k = i1s(3),i1e(3)
          do j = i1s(2),i1e(2)
          do i = i1s(1),i1e(1)

            iv4 = surf(ip)%first_node  &
                + (i-i1(ip))  &
                + (j-j1(ip))*(i2(ip)-i1(ip)+1)  &
                + (k-k1(ip))*(j2(ip)-j1(ip)+1)*(i2(ip)-i1(ip)+1)

            if (dir1 == 1) then !i-surface
              iv3 = iv4 - (i2(ip)-i1(ip)+1)
              iv2 = iv4 - (j2(ip)-j1(ip)+1)*(i2(ip)-i1(ip)+1)
              iv1 = iv2 - (i2(ip)-i1(ip)+1)
            else if (dir1 == 2) then !j-surface
              iv3 = iv4 - (j2(ip)-j1(ip)+1)*(i2(ip)-i1(ip)+1)
              iv2 = iv4 - 1
              iv1 = iv3 - 1
            else if (dir1 == 3) then !k-surface
              iv3 = iv4 - 1
              iv2 = iv4 - (i2(ip)-i1(ip)+1)
              iv1 = iv2 - 1
            end if

            scal( 2,ipan) = 0.25_rd *(qv( 1,iv1)+qv( 1,iv2) &
                                     +qv( 1,iv3)+qv( 1,iv4)) * qinf ! p
            scal( 3,ipan) = 0.25_rd *(qv( 8,iv1)+qv( 8,iv2) &
                                     +qv( 8,iv3)+qv( 8,iv4))        ! du
            scal( 4,ipan) = 0.25_rd *(qv( 9,iv1)+qv( 9,iv2) &
                                     +qv( 9,iv3)+qv( 9,iv4))        ! dv
            scal( 5,ipan) = 0.25_rd *(qv(10,iv1)+qv(10,iv2) &
                                     +qv(10,iv3)+qv(10,iv4))        ! dw
            scal( 6,ipan) = 0.25_rd *(qv( 7,iv1)+qv( 7,iv2) &
                                     +qv( 7,iv3)+qv( 7,iv4))        ! vmu
            scal( 8,ipan) = 0.25_rd *(qv( 2,iv1)+qv( 2,iv2) &
                                     +qv( 2,iv3)+qv( 2,iv4))        ! rho
            scal( 9,ipan) = 0.25_rd *(qv( 3,iv1)/qv( 2,iv1) &
                                    + qv( 3,iv2)/qv( 2,iv2) &
                                    + qv( 3,iv3)/qv( 2,iv3) &
                                    + qv( 3,iv4)/qv( 2,iv4))        ! u
            scal(10,ipan) = 0.25_rd *(qv( 4,iv1)/qv( 2,iv1) &
                                    + qv( 4,iv2)/qv( 2,iv2) &
                                    + qv( 4,iv3)/qv( 2,iv3) &
                                    + qv( 4,iv4)/qv( 2,iv4))        ! v
            scal(11,ipan) = 0.25_rd *(qv( 5,iv1)/qv( 2,iv1) &
                                    + qv( 5,iv2)/qv( 2,iv2) &
                                    + qv( 5,iv3)/qv( 2,iv3) &
                                    + qv( 5,iv4)/qv( 2,iv4))        ! w
            scal(12,ipan) = 0.25_rd *(qv( 3,iv1)+qv( 3,iv2) &
                                     +qv( 3,iv3)+qv( 3,iv4))        ! ru
            scal(13,ipan) = 0.25_rd *(qv( 4,iv1)+qv( 4,iv2) &
                                     +qv( 4,iv3)+qv( 4,iv4))        ! rv
            scal(14,ipan) = 0.25_rd *(qv( 5,iv1)+qv( 5,iv2) &
                                     +qv( 5,iv3)+qv( 5,iv4))        ! rw

            ipan = ipan + 1

          end do
          end do
          end do
        end if !solution_exists


!       store the i-blank variable

        ipan = surf(ip)%first_panel

        do k = i1s(3),i1e(3)
        do j = i1s(2),i1e(2)
        do i = i1s(1),i1e(1)

          i2e(1) = i
          i2e(2) = j
          i2e(3) = k
          if (i2e(dir1) == 1) i2e(dir1) = i2e(dir1) + 1

          if (dir1 == 1) then
            ibv1 = ibv(i,j  ,k  )
            ibv2 = ibv(i,j-1,k  )
            ibv3 = ibv(i,j  ,k-1)
            ibv4 = ibv(i,j-1,k-1)
          else if (dir1 == 2) then
            ibv1 = ibv(i  ,j,k  )
            ibv2 = ibv(i-1,j,k  )
            ibv3 = ibv(i  ,j,k-1)
            ibv4 = ibv(i-1,j,k-1)
          else if (dir1 == 3) then
            ibv1 = ibv(i  ,j  ,k)
            ibv2 = ibv(i-1,j  ,k)
            ibv3 = ibv(i  ,j-1,k)
            ibv4 = ibv(i-1,j-1,k)
          end if

          scal(1,ipan) = InterpolateIBlankFull(ibv1,ibv2,ibv3,ibv4)
          if (solution_exists) then
            scal(7,ipan) = volv(i2e(1),i2e(2),i2e(3))
          end if

          ipan = ipan + 1
        end do
        end do
        end do

      end if
    end do
    call my_cpu_time(time4)
    ltime = ltime + (time4-time3)

    deallocate(ibv,zb,yb,xb); nalloc = nalloc - 4
    if (solution_exists) then
      deallocate(volv,qd); nalloc = nalloc - 2
    end if

  end do

  close(unit=3)
  if (solution_exists) close(unit=4)

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)"cpu time in StorePatchesOVERFLOW:",time2-time1
    write(0,cpufmt)"  > time in grid I/O:",utime
    write(0,cpufmt)"  > time in volumes:",vtime
    write(0,cpufmt)"  > time in storage:",ltime
  end if

  if (allocated(rgas)) then
    deallocate(rgas); nalloc = nalloc - 1
  end if
  deallocate(nkg,njg,nig); nalloc = nalloc - 3
  if ((allocated(qv)).and.(.not.full_surface)) then
    deallocate(qv); nalloc = nalloc - 1
  end if

  return
  end subroutine StorePatchesOVERFLOW

! =====================================================================
  subroutine ReadFlow(iu,ni,nj,nk,nq)
! =====================================================================
! Read flow data out of the q.save file (unit number iu)
! for the current block for both single- and double-precision files.

  use NetAlloc    ,only: nalloc
  use UserInput   ,only: double
  use VolumeArrays,only: qd,qs

  implicit none

  integer,intent(in) :: iu,ni,nj,nk,nq


  continue


  if (double) then
    read(iu)qd
  else
    allocate(qs(ni,nj,nk,nq)); nalloc = nalloc + 1
    read(iu)qs
    qd = qs
    deallocate(qs); nalloc = nalloc - 1
  end if


  return
  end subroutine ReadFlow

! =====================================================================
  subroutine ReadGrid(iu,ni,nj,nk)
! =====================================================================
! Read grid data out of the grid.in file (unit number iu)
! for the current block for both single- and double-precision files,
! with or without iblank data.

  use NetAlloc    ,only: nalloc
  use UserInput   ,only: double,ibex
  use VolumeArrays,only: xb,yb,zb,xs,ys,zs,ibv

  implicit none

  integer,intent(in) :: ni,nj,nk
  integer,intent(in) :: iu


  continue


  if (double) then
    if (ibex) then
      read(iu)xb,yb,zb,ibv
    else
      read(iu)xb,yb,zb
      ibv = 1
    end if
  else
    allocate(xs(ni,nj,nk)); nalloc = nalloc + 1
    allocate(ys(ni,nj,nk)); nalloc = nalloc + 1
    allocate(zs(ni,nj,nk)); nalloc = nalloc + 1
    if (ibex) then
      read(iu)xs,ys,zs,ibv
    else
      read(iu)xs,ys,zs
      ibv = 1
    end if
    xb = xs
    yb = ys
    zb = zs
    deallocate(zs,ys,xs); nalloc = nalloc - 3
  end if


  return
  end subroutine ReadGrid

! ===============================================================
  subroutine ReadHeader(iu,iqtyp,nqc, &
     fsmach,alpha,rey,time,gaminf,beta,tinf,htinf,ht1,ht2,rgas, &
     igam)
! ===============================================================
! Ascertain whether the flow conditions in q.save are in
! plot3d format or in either of two OVERFLOW solution file
! formats.  Store single-precision data into double-precision
! variables if q.save is single-precision.

  use IntrType ,only: rd,rs
  use UserInput,only: multiple_grid,double

  implicit none

  integer,intent(in) :: nqc

  real(kind=rd),dimension(nqc),intent(out) :: rgas

  real(kind=rd),intent(out) :: fsmach,alpha,rey,time,gaminf
  real(kind=rd),intent(out) :: beta,tinf,htinf,ht1,ht2

  real(kind=rs),dimension(nqc) :: rgas_s

  real(kind=rs) :: fsmach_s,alpha_s,rey_s,time_s,gaminf_s
  real(kind=rs) :: beta_s,tinf_s,htinf_s,ht1_s,ht2_s

  integer,intent(in)  :: iu
  integer,intent(out) :: iqtyp,igam
  integer             :: ierr


  continue


  if (double) then
    read(iu,iostat=ierr)fsmach,alpha,rey,time, &
                        gaminf,beta,tinf, &
                        igam,htinf,ht1,ht2,rgas
  else
    read(iu,iostat=ierr)fsmach_s,alpha_s,rey_s,time_s, &
                        gaminf_s,beta_s,tinf_s, &
                        igam,htinf_s,ht1_s,ht2_s,rgas_s
  end if

  if (ierr == 0) then
    if (.not.double) then
      fsmach = fsmach_s
      alpha  = alpha_s
      rey    = rey_s
      time   = time_s
      gaminf = gaminf_s
      beta   = beta_s
      tinf   = tinf_s
      htinf  = htinf_s
      ht1    = ht1_s
      ht2    = ht2_s
      rgas   = rgas_s
    end if
    iqtyp = 2
    backspace(iu)
    return
  end if

  rewind(iu)
  if (multiple_grid) read(iu)
  read(iu)

  if (double) then
    read(iu,iostat=ierr)fsmach,alpha,rey,time, &
                        gaminf, &
                        igam,htinf,ht1,ht2,rgas
  else
    read(iu,iostat=ierr)fsmach_s,alpha_s,rey_s,time_s, &
                        gaminf_s, &
                        igam,htinf_s,ht1_s,ht2_s,rgas_s
  end if
  if (ierr == 0) then
    if (.not.double) then
      fsmach = fsmach_s
      alpha  = alpha_s
      rey    = rey_s
      time   = time_s
      gaminf = gaminf_s
      htinf  = htinf_s
      ht1    = ht1_s
      ht2    = ht2_s
      rgas   = rgas_s
    end if
    beta = -400.0
    tinf = -1.0
    iqtyp = 1
    backspace(iu)
    return
  end if

  rewind(iu)
  if (multiple_grid) read(iu)
  read(iu)

  if (double) then
    read(iu,iostat=ierr) fsmach,alpha,rey,time
  else
    read(iu,iostat=ierr) fsmach_s,alpha_s,rey_s,time_s
  end if

  if (ierr == 0) then
    if (.not.double) then
      fsmach = fsmach_s
      alpha  = alpha_s
      rey    = rey_s
      time   = time_s
      htinf  = htinf_s
      ht1    = ht1_s
      ht2    = ht2_s
      rgas   = rgas_s
    end if
    beta = -400.0
    tinf = -1.0
    gaminf = -1.0
    iqtyp = 0
    backspace(iu)
    return
  end if

  write(0,*)"ERROR! could not resolve format for q.save."
  stop

  end subroutine ReadHeader

! =================================================================
  pure function InterpolateIBlankFull(in1,in2,in3,in4) result (ic)
! =================================================================
! Function to move the iblank variable from corners to
! the panel center for the case where iblank can include
! negative values (fringe) and flags for ophans (101),
! consistent with usage with SUGGAR.

  implicit none

  integer,intent(in) :: in1,in2,in3,in4
  integer            :: iv1,iv2,iv3,iv4
  integer            :: ic


  continue


! do not differentiate among different types of fringe nodes
  iv1 = max(in1,-1)
  iv2 = max(in2,-1)
  iv3 = max(in3,-1)
  iv4 = max(in4,-1)


! if all values are equal, return that value
  if ((iv1 == iv2).and.(iv2 == iv3).and.(iv3 == iv4)) then
    ic = iv1
    return
  end if


! if any value is 101 (orphan), return that value
  if (max(iv1,iv2,iv3,iv4) == 101) then
    ic = 101
    return
  end if

! if any value is zero (hole), return that value
  if (iv1*iv2*iv3*iv4 == 0) then
    ic = 0
    return
  end if

! we must have a mix of fringe and interior cells, return -1
  ic = -1
  return

  end function InterpolateIBlankFull

  end module StorePatchesOVERFLOWM
