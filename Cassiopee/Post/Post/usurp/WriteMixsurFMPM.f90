! ================================================================
  module WriteMixsurFMPM
! ================================================================

  implicit none

  private
  public :: WriteMixsurFMP

  contains


! ======================================================
  subroutine WriteMixsurFMP()
! ======================================================
! Write parameters file for input to OVERFLOW
! adapted from WRMXPA in chimera1.9/src/mixsur.f

  use GroupInfo   ,only: group_iref,group_name
  use GroupInfo   ,only: group_size,group_data,num_groups
  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use PatchInfo   ,only: num_components,num_patches,ncomp_surfs
  use PatchInfo   ,only: i1,i2,j1,j2,k1,k2,icomp,ibdir,block,refc
  use UserInput   ,only: cpufmt,pathname,showcpu
  use UserInput   ,only: fsmach,alpha,beta,rey,gaminf,tinf
  use UserInput   ,only: nref,refq
  use VertexData  ,only: num_nodes

  implicit none

  real(kind=rs)     :: time1,time2
  integer           :: nzip,ntri,nr,nsurf,nf,i,nc,motype
  integer           :: nsubt,ncomp,ntot,nu
  character(len=80) :: group80


  continue


  call my_cpu_time(time1)

  open (unit=7,file=trim(pathname)//'mixsur.fmp', &
        status='unknown', &
        action='write', &
        form='formatted')

! translate nomenclature from USURP to MIXSUR
  nsurf = num_components
  nsubt = num_patches
  nzip  = 0
  ntri  = 0
  ntot  = num_nodes
  ncomp = num_groups


  write(7,101) nsurf,nsubt,nzip,ntri,ntot,ncomp,nref
! if (fsmach /= 0.0) then
!   write(7,102) fsmach,alpha,beta,rey*fsmach,gaminf,tinf
! else
    write(7,102) fsmach,alpha,beta,rey,gaminf,tinf
! end if
  do nr=1,nref
    write(7,103) (refq(i,nr),i=1,5)
  end do

! write(7,104) (nsub(nf),nf=1,nsurf)
  write(7,104) (ncomp_surfs(nf),nf=1,nsurf)

! write(7,104) (irefs(nf),nf=1,nsurf)
  write(7,104) (refc(nf),nf=1,nsurf)

! do nu=1,nsubt
!   isin(nu,1) = igrd(nu)
! end do
  do nu=1,nsubt
!   write(7,105) (isin(nu,i),i=1,9)
    write(7,105) block(nu), &                                  !1
                 ibdir(nu), &                                  !2
                 i1(nu),i2(nu),j1(nu),j2(nu),k1(nu),k2(nu), &  !3-8
                 icomp(nu)                                     !9
  end do

! read and write component info.

  do nc=1,ncomp
!   write(7,106) cname(nc)
!   write(7,107) nis(nc),irefc(nc)
!   write(7,104) (icslist(nc,i),i=1,nis(nc))

!   change the name to an 80-string character to match mixsur.f
    group80 = group_name(nc)

!   write(7,106) group_name(nc)
    write(7,106) group80
    write(7,107) group_size(nc),group_iref(nc)
    write(7,104) (group_data(i,nc),i=1,group_size(nc))
  end do

! write more reference set info. do it here so that things can be
! backward compatible.

  do nr=1,nref
!   write(7,108) motype(nr)
!   if (motype(nr).eq.2) then
!     write(7,109) refq(6,nr),refq(7,nr)
!   else if (motype(nr).eq.3) then
!     write(7,110) refq(6,nr),refq(7,nr),refq(8,nr)
!   end if

    motype = nint(refq(9,nr))
    write(7,108) motype
    if (motype == 2) then
      write(7,109) refq(6,nr),refq(7,nr)
    else if (motype == 3) then
      write(7,110) refq(6,nr),refq(7,nr),refq(8,nr)
    end if

  end do

  close(7)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time to write mixsur.fmp: ',time2-time1
  end if


 101  format(7i10)
 102  format(6e16.7)
 103  format(5e16.7)
 104  format(10i5)
 105  format(9i6)
 106  format(a)
 107  format(2i5)
 108  format(i5)
 109  format(2e16.7)
 110  format(3e16.7)

  return
  end subroutine WriteMixsurFMP

  end module WriteMixsurFMPM
