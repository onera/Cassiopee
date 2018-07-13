! ================================================================
  module WritePatchesM
! ================================================================

  implicit none

  private
  public :: WritePatches

  contains

! ====================================================================
  subroutine WritePatches(p)
! ====================================================================
! This procedure writes the panel data back into a Tecplot
! file while maintaining the structured nature of the
! original boundary patches.  Note that the dummy argument p
! is a pointer, so an explicit interface is required.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: i1,i2,j1,j2,k1,k2,num_patches
  use PatchInfo   ,only: color,tag_list,icomp
  use Types       ,only: panel,scal,surf
  use UserInput   ,only: colormap,tecbinary,pathname,solution_exists
  use UserInput   ,only: showcpu,cpufmt,plotvel
  use VertexData  ,only: xv

  implicit none

  type(panel),pointer,dimension (:) :: p
  real,allocatable,dimension (:)    :: v2
  real(kind=rs)                     :: time1,time2

  integer                           :: ip,ipan,iv,idir
  integer                           :: zni,znj,iii
  integer,allocatable,dimension (:) :: varloc
#ifdef USE_TECPLOT
  integer                           :: ierr
  integer,parameter                 :: zero=0
  integer,parameter                 :: visdouble=0
  integer,parameter                 :: debugtecplot=0
  integer                           :: tecini100,tecdat100,teczne100
#endif

  character(len=132),dimension(4)   :: zonestuff
  character(len=132)                :: tecvariables
  character(len=8)                  :: colorname


  continue


  call my_cpu_time(time1)

  if (solution_exists) then
    if (plotvel) then
      tecvariables = "x,y,z,C<sub>p</sub>,iblank,ratio,ipan,"// &
                     "component,u,v,w,C<sub>f</sub>"
      zonestuff(3) = "  varlocation=([4-12]=cellcentered)"
      allocate(varloc(12)); nalloc = nalloc + 1
      varloc = (/ 1,1,1,0,0,0,0,0,0,0,0,0 /)
    else
      tecvariables = "x,y,z,C<sub>p</sub>,iblank,ratio,ipan,component"
      zonestuff(3) = "  varlocation=([4-8]=cellcentered)"
      allocate(varloc(8)); nalloc = nalloc + 1
      varloc = (/ 1,1,1,0,0,0,0,0 /)
    end if
  else
    tecvariables = "x,y,z,iblank,ratio,ipan,component"
    zonestuff(3) = "  varlocation=([4-7]=cellcentered)"
    allocate(varloc(7)); nalloc = nalloc + 1
    varloc = (/ 1,1,1,0,0,0,0 /)
  end if

  if (tecbinary) then

#ifdef USE_TECPLOT
    ierr = tecini100 ( "usurp"//char(0), &
                       trim(tecvariables)//char(0), &
                       trim(pathname)//"usurp-surfaces.plt"//char(0), &
                       trim(pathname)//char(0), &
                       debugtecplot, &
                       visdouble )
    if (colormap) then
      open(unit=43,                          &
           file=trim(pathname)//"color.mcr", &
           action="write")
      write(unit=43,fmt="(a)")"#!MC 1000"
    end if

#endif
  else
    zonestuff(1) = ",zonetype=ordered"
    zonestuff(2) = "  datapacking=block"
    open(unit=42,file=trim(pathname)//"usurp-surfaces.dat",action="write")
    write(unit=42,fmt=*)"variables = "//trim(tecvariables)
  end if


  do ip = 1,num_patches

!   alter the zone color map
    if (colormap) then
      select case (color(ip))
        case(1)
          colorname = "red"
        case(2)
          colorname = "green"
        case(3)
          colorname = "blue"
        case(4)
          colorname = "cyan"
        case(5)
          colorname = "yellow"
        case(6)
          colorname = "purple"
        case default
          colorname = "custom2"
      end select
      zonestuff(4)="  c="//trim(colorname)
    end if

!   write Tecplot zone headers
    if (i1(ip) == i2(ip)) then
      zni = j2(ip)-j1(ip)+1
      znj = k2(ip)-k1(ip)+1
    else if (j1(ip) == j2(ip)) then
      zni = i2(ip)-i1(ip)+1
      znj = k2(ip)-k1(ip)+1
    else if (k1(ip) == k2(ip)) then
      zni = i2(ip)-i1(ip)+1
      znj = j2(ip)-j1(ip)+1
    end if

    if (tecbinary) then
#ifdef USE_TECPLOT
      ierr = teczne100(trim(tag_list(icomp(ip)))//char(0), &
                       0, &
                       zni,znj,1, &
                       0,0,0, &
                       1,0,0, &
                       varloc, &
                       %val(zero), &
!                      zero, &
                       0)
       if (colormap) then
         write(43,"(a,i4,a)") &
          "$!Field [",ip,"]  Mesh{Color = "//trim(colorname)//"}"
         write(43,"(a,i4,a)") &
          "$!Field [",ip,"]  Shade{Color = "//trim(colorname)//"}"
       end if
#endif
    else
      write(unit=42,fmt=*)"zone i=",zni,",j=",znj,trim(zonestuff(1))
      write(unit=42,fmt="(a)")trim(zonestuff(2))
      write(unit=42,fmt="(a)")trim(zonestuff(3))
      if (colormap) write(unit=42,fmt="(a)")trim(zonestuff(4))
    end if


!   output the grid points
    if (tecbinary) then
      iii = surf(ip)%last_node - surf(ip)%first_node + 1
      allocate(v2(iii)); nalloc = nalloc + 1
    end if

    do idir = 1,3
      if (tecbinary) then
        do iv = surf(ip)%first_node,surf(ip)%last_node
          v2(iv - surf(ip)%first_node + 1) = xv(idir,iv)
        end do
#ifdef USE_TECPLOT
        ierr = tecdat100(iii,v2,visdouble)
#endif
      else
        write(unit=42,fmt="(5es16.7)") &
             (xv(idir,iv),iv=surf(ip)%first_node,surf(ip)%last_node)
      end if
    end do

    if (tecbinary) then
      deallocate(v2); nalloc = nalloc - 1
      iii = surf(ip)%last_panel - surf(ip)%first_panel + 1
      allocate(v2(iii)); nalloc = nalloc + 1
    end if

!   write static pressure coefficient data (if it exists)
    if (solution_exists) then
      if (tecbinary) then
        do ipan = surf(ip)%first_panel,surf(ip)%last_panel
          v2(ipan - surf(ip)%first_panel + 1) = 2.0*scal(2,ipan)
        end do
#ifdef USE_TECPLOT
        ierr = tecdat100(iii,v2,visdouble)
#endif
      else
        write(unit=42,fmt="(5es16.7)") &
          (2.0*scal(2,ipan),ipan=surf(ip)%first_panel,surf(ip)%last_panel)
      end if
    end if


!   write i-blanking data
    if (tecbinary) then
      do ipan = surf(ip)%first_panel,surf(ip)%last_panel
        v2(ipan - surf(ip)%first_panel + 1) = p(ipan)%iblank
      end do
#ifdef USE_TECPLOT
      ierr = tecdat100(iii,v2,visdouble)
#endif
    else
      write(unit=42,fmt="(5i7    )") &
        (p(ipan)%iblank,ipan=surf(ip)%first_panel,surf(ip)%last_panel)
    end if


!   write the panel weights
    if (tecbinary) then
      do ipan = surf(ip)%first_panel,surf(ip)%last_panel
        v2(ipan - surf(ip)%first_panel + 1) = p(ipan)%ratio
      end do
#ifdef USE_TECPLOT
      ierr = tecdat100(iii,v2,visdouble)
#endif
    else
      write(unit=42,fmt="(5es16.7)") &
        (p(ipan)%ratio ,ipan=surf(ip)%first_panel,surf(ip)%last_panel)
    end if


!   write the panel number (for debugging problems)
    if (tecbinary) then
      do ipan = surf(ip)%first_panel,surf(ip)%last_panel
        v2(ipan - surf(ip)%first_panel + 1) = ipan
      end do
#ifdef USE_TECPLOT
      ierr = tecdat100(iii,v2,visdouble)
#endif
    else
      write(unit=42,fmt="(5i7    )") &
        (ipan          ,ipan=surf(ip)%first_panel,surf(ip)%last_panel)
    end if


!   write the component number
    if (tecbinary) then
      do ipan = surf(ip)%first_panel,surf(ip)%last_panel
        v2(ipan - surf(ip)%first_panel + 1) = p(ipan)%icomp
      end do
#ifdef USE_TECPLOT
      ierr = tecdat100(iii,v2,visdouble)
#endif
    else
      write(unit=42,fmt="(5i7    )") &
        (p(ipan)%icomp,ipan=surf(ip)%first_panel,surf(ip)%last_panel)
    end if


!   write the velocity and skin friction data (if it exists)
    if (plotvel) then
      do idir = 1,4
      if (tecbinary) then
        do ipan = surf(ip)%first_panel,surf(ip)%last_panel
          v2(ipan - surf(ip)%first_panel + 1) = scal(8+idir,ipan)
        end do
#ifdef USE_TECPLOT
        ierr = tecdat100(iii,v2,visdouble)
#endif
      else
        write(unit=42,fmt="(5es16.7)") &
          (scal(8+idir,ipan),ipan=surf(ip)%first_panel,surf(ip)%last_panel)
      end if
      end do !idir
    end if

    if (tecbinary) then
      deallocate(v2); nalloc = nalloc - 1
    end if

  end do


  deallocate(varloc); nalloc = nalloc - 1


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in WritePatches:',time2-time1
  end if


  return
  end subroutine WritePatches

  end module WritePatchesM
