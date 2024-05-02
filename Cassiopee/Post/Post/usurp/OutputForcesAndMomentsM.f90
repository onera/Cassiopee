! ================================================================
  module OutputForcesAndMomentsM
! ================================================================

  implicit none

  private
  public  :: OutputForcesAndMoments
  private :: BodyForces

  contains

! ====================================================================
  subroutine OutputForcesAndMoments()
! ====================================================================
! Create a text output of the integrated force and moment data.
! The file is intended to facilitate import into spreadsheet software
! (e.g. Microsoft Excel) as a space delimited ASCII file and is
! consistent with component_data.txt from OVER-REL.

  use FomoArrays  ,only: cfpx,cfpy,cfpz,cmpx,cmpy,cmpz
  use FomoArrays  ,only: cfvx,cfvy,cfvz,cmvx,cmvy,cmvz
  use FomoArrays  ,only: weta
  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use PatchInfo   ,only: num_components,ncomp_surfs,tag_list
  use UserInput   ,only: cpufmt,pathname,showcpu

  implicit none

  real(kind=rd),dimension(6) :: ctot,csub
  real(kind=rd),dimension(3) :: bff,bfm
  real(kind=rs)              :: time1,time2
  integer                    :: i,ic,ic1
  character(len=20)          :: label


  continue


  call my_cpu_time(time1)

! start loops at index=0 if there are untagged surfaces
  if (ncomp_surfs(0) > 0) then
    ic1 = 0
  else
    ic1 = 1
  end if


! write an Excel file with the component data
  open(unit=3,file=trim(pathname)//"usurp_data.txt",action="write")

  label = '"PRESSURE FORCE"'
  write(3,10)label
  csub = 0.0
  do ic = ic1,num_components
    write(3,20)tag_list(ic),cfpx(ic),cfpy(ic),cfpz(ic), &
                            cmpx(ic),cmpy(ic),cmpz(ic)
    csub(1) = csub(1) + cfpx(ic)
    csub(2) = csub(2) + cfpy(ic)
    csub(3) = csub(3) + cfpz(ic)
    csub(4) = csub(4) + cmpx(ic)
    csub(5) = csub(5) + cmpy(ic)
    csub(6) = csub(6) + cmpz(ic)
  end do
  label = 'subtotal'
  write(3,20)label,(csub(i),i=1,6)


  label = '"VISCOUS FORCE"'
  write(3,10)label
  csub = 0.0
  do ic = ic1,num_components
    write(3,20)tag_list(ic),cfvx(ic),cfvy(ic),cfvz(ic), &
                            cmvx(ic),cmvy(ic),cmvz(ic)
    csub(1) = csub(1) + cfvx(ic)
    csub(2) = csub(2) + cfvy(ic)
    csub(3) = csub(3) + cfvz(ic)
    csub(4) = csub(4) + cmvx(ic)
    csub(5) = csub(5) + cmvy(ic)
    csub(6) = csub(6) + cmvz(ic)
  end do
  label = 'subtotal'
  write(3,20)label,(csub(i),i=1,6)


  label = '"COMBINED FORCE"'
  write(3,10)label
  do ic = ic1,num_components
    write(3,20)tag_list(ic), &
       cfpx(ic)+cfvx(ic),cfpy(ic)+cfvy(ic),cfpz(ic)+cfvz(ic), &
       cmpx(ic)+cmvx(ic),cmpy(ic)+cmvy(ic),cmpz(ic)+cmvz(ic)
  end do

  ctot = 0.0_rd

  do ic = ic1,num_components
    ctot(1) = cfpx(ic) + cfvx(ic) + ctot(1)
    ctot(2) = cfpy(ic) + cfvy(ic) + ctot(2)
    ctot(3) = cfpz(ic) + cfvz(ic) + ctot(3)
    ctot(4) = cmpx(ic) + cmvx(ic) + ctot(4)
    ctot(5) = cmpy(ic) + cmvy(ic) + ctot(5)
    ctot(6) = cmpz(ic) + cmvz(ic) + ctot(6)
  end do

  call BodyForces(bff,bfm)
  if ((maxval(abs(bff))/=0.0).or.(maxval(abs(bfm))/=0.0)) then
    label = 'subtotal'
    write(3,20)label,(ctot(i),i=1,6)
    write(3,*)
    label = '"BODY FORCE"'
    write(3,20)label,bff(1),bff(2),bff(3), &
                     bfm(1),bfm(2),bfm(3)
    ctot(1:3) = ctot(1:3) + bff
    ctot(4:6) = ctot(4:6) + bfm
  end if

  label = 'total'
  write(3,20)label,(ctot(i),i=1,6)


  label = '"WETTED AREA CALCS"'
  write(3,15)label
  csub = 0.0
  do ic = ic1,num_components
    write(3,20)tag_list(ic),weta(ic),cfpx(ic)/weta(ic), &
                                     cfvx(ic)/weta(ic), &
                          (cfpx(ic)+cfvx(ic))/weta(ic)
    csub(1) = csub(1) + weta(ic) 
    csub(2) = csub(2) + cfpx(ic) 
    csub(3) = csub(3) + cfvx(ic) 
  end do
  label = 'total'
  write(3,20)label,csub(1),csub(2)/csub(1), &
                           csub(3)/csub(1), &
                 (csub(2)+csub(3))/csub(1)

  close(3)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in OutputForcesAndMoments:',time2-time1
  end if


  return

 10 format(/,a20, 6x,'CFX',10x,'CFY',10x,'CFZ',10x,'CMX',10x,'CMY',10x,'CMZ')
 15 format(/,a20, 6x,'AREA',6x,'CFXP/AREA',4x,'CFVX/AREA',5x,'CFX/AREA')
 20 format(a20,6(2x,e11.4))

  end subroutine OutputForcesAndMoments

! ===================================================================
  subroutine BodyForces(bff,bfm)
! ===================================================================
! This routine integrates body forces used in UNCLE/OVER-REL cases

  use CalculateVolumesM,only: CalculateVolumes
  use IntrType         ,only: rd
  use NetAlloc         ,only: nalloc
  use UserInput        ,only: workdir,xcg,ycg,zcg
  use VolumeArrays     ,only: volv,xb,yb,zb

  implicit none

  real(kind=rd),allocatable,dimension (:,:,:) :: bfx,bfy,bfz

  real(kind=rd),dimension(3),intent(out) :: bff,bfm
  real(kind=rd) :: xc,yc,zc

  integer           :: n,nn,ni,nj,nk,i,j,k
  character(len=80) :: gname,fname
  character(len=5)  :: fileid
  logical           :: exg,exf


  continue


  n      = 0
  bff    = 0.0_rd
  bfm    = 0.0_rd
  exg    = .true.

  write(*,*)
  write(*,*)'searching for F* files ...'

  do while (exg)
    n = n + 1

    write(unit=fileid,fmt="(i5)")n*10
    gname = trim(workdir)//"IO_FILES/G"//adjustl(fileid)

    inquire(file=gname,exist=exg)
    if (exg) then

!     read the forces
      fname  = trim(workdir)//"IO_FILES/F"//adjustl(fileid)

      inquire(file=fname,exist=exf)
      if (exf) then

        open(unit=4,file=gname,form='unformatted')
        read(4)nn
        read(4)ni,nj,nk
        allocate(xb(ni,nj,nk),yb(ni,nj,nk),zb(ni,nj,nk)); nalloc = nalloc + 3
        read(4)xb,yb,zb
        close(4)

        write(*,20)'F'//adjustl(fileid),ni,nj,nk
  20    format(1x,a5,2x,'(',i4,'x',i4,'x',i4,')')

        allocate(volv(ni+1,nj+1,nk+1)); nalloc = nalloc + 1
        call CalculateVolumes(ni,nj,nk)

        allocate(bfx(ni,nj,nk),bfy(ni,nj,nk),bfz(ni,nj,nk)); nalloc = nalloc + 3
        open(unit=4,file=fname,form='unformatted')
        read(4)bfx
        read(4)bfy
        read(4)bfz
        close(4)


!       check forces
        do k = 2,nk
        do j = 2,nj
        do i = 2,ni

          xc = 0.125_rd*(xb(i,j,k)     + xb(i-1,j,k)    &
                       + xb(i,j-1,k)   + xb(i-1,j-1,k)  &
                       + xb(i,j,k-1)   + xb(i-1,j,k-1)  &
                       + xb(i,j-1,k-1) + xb(i-1,j-1,k-1))
          yc = 0.125_rd*(yb(i,j,k)     + yb(i-1,j,k)    &
                       + yb(i,j-1,k)   + yb(i-1,j-1,k)  &
                       + yb(i,j,k-1)   + yb(i-1,j,k-1)  &
                       + yb(i,j-1,k-1) + yb(i-1,j-1,k-1))
          zc = 0.125_rd*(zb(i,j,k)     + zb(i-1,j,k)    &
                       + zb(i,j-1,k)   + zb(i-1,j-1,k)  &
                       + zb(i,j,k-1)   + zb(i-1,j,k-1)  &
                       + zb(i,j-1,k-1) + zb(i-1,j-1,k-1))


!         note:  the negative sign below is because we want the
!         force of the fluid acting on the "body"

          bff(1) = bff(1) - volv(i,j,k)*bfx(i,j,k)
          bff(2) = bff(2) - volv(i,j,k)*bfy(i,j,k)
          bff(3) = bff(3) - volv(i,j,k)*bfz(i,j,k)
          bfm(1) = bfm(1) - volv(i,j,k)*(bfz(i,j,k)*(yc-ycg) &
                                        -bfy(i,j,k)*(zc-zcg))
          bfm(2) = bfm(2) - volv(i,j,k)*(bfx(i,j,k)*(zc-zcg) &
                                        -bfz(i,j,k)*(xc-xcg))
          bfm(3) = bfm(3) - volv(i,j,k)*(bfy(i,j,k)*(xc-xcg) &
                                        -bfx(i,j,k)*(yc-ycg))

        end do
        end do
        end do

        deallocate(bfz,bfy,bfx); nalloc = nalloc - 3
        deallocate(volv)       ; nalloc = nalloc - 1
        deallocate(zb,yb,xb)   ; nalloc = nalloc - 3

      end if
    end if
  end do

  return
  end subroutine BodyForces

  end module OutputForcesAndMomentsM
