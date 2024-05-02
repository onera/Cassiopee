! ================================================================
  module PreReadOVERFLOWM
! ================================================================

  implicit none

  private
  public  :: PreReadOVERFLOW
  private :: GetInput 
  private :: GridType
  private :: isanum
  private :: rmsp

  contains

! ==========================================================
  subroutine PreReadOVERFLOW()
! ==========================================================
! read the integration group data for fomoco / overflow in 
! order to acquire the number of panels, patches, edges, and 
! vertices, and to store the patch data.

  use EdgeData                 ,only: num_edges
  use GroupInfo                ,only: group_data,group_iref
  use GroupInfo                ,only: group_name,group_size
  use GroupInfo                ,only: num_groups
  use IntrType                 ,only: rd
  use NetAlloc                 ,only: nalloc
  use PatchInfo                ,only: block,ibdir,i1,i2,j1,j2,k1,k2
  use PatchInfo                ,only: max_components,num_components
  use PatchInfo                ,only: num_panels,num_patches
  use PatchInfo                ,only: icomp,ncomp_surfs
  use PatchInfo                ,only: refc,tag_list
  use PriPairs                 ,only: ipri
  use StructuredPatchProcedures,only: CountStructuredPatch
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use Types                    ,only: surf
  use UserInput                ,only: pathname,workdir,verbosity
  use UserInput                ,only: fsmach,alpha,beta,rey,gaminf,tinf
  use UserInput                ,only: nref,refq
  use UserInput                ,only: double,ibex,multiple_grid
  use UserInput                ,only: use_priority
  use VertexData               ,only: num_nodes

  implicit none

  real(kind=rd) :: xin(20)

  integer,allocatable,dimension(:) :: itmp,nig,njg,nkg

  integer :: ibdirs,ic,ierr,ip,iu
  integer :: i1s,i2s,j1s,j2s,k1s,k2s,its
  integer :: j,k,loop
  integer :: n,nb,nin,npri,nrsub,nu1,nu2
  integer :: npri_total,nsub0

  character(len=120) :: line
  character(len=5)   :: tagn


  continue


! initialization
  npri_total = 0


! need to read the top of the grid file to get the block dimensions
  open ( unit    = 3,                          &
         file    = trim(workdir)//"grid.in",   &
         form    = "unformatted",              &
! #ifdef CONVERT_BIG_ENDIAN
         ! convert = "big_endian",               &
! #endif
         status  = "old",                      &
         action  = "read"                      )

! ascertain the grid format
  call GridType(3)

  if (multiple_grid) then
    read(unit=3)nb
  else
    nb = 1
  end if
  allocate(nig(nb),njg(nb),nkg(nb)); nalloc = nalloc + 3
  read(unit=3)(nig(n),njg(n),nkg(n),n=1,nb)
  close(unit=3)

  write(*,*)"Reading from standard input."


! read operating conditions from standard input
  read(5,*)fsmach,alpha,beta,rey,gaminf,tinf


! read the reference parameter sets from standard input
  read(5,*)nref

  allocate(refq(9,nref),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(unit=0,fmt=*)"ERROR! allocate failed for refq"
    write(unit=0,fmt=*)"nref = ",nref
    stop
  end if

  do j = 1,nref

!   read(5,*)(refq(k,j),k=1,5)

    ierr = 1 !signals blank line
    do while (ierr == 1)
      read(unit=5,fmt="(a)",iostat=ierr) line
      if (ierr /= 0) then
        write(unit=0,fmt=*)"ERROR! unable to read reference conditions."
        stop
      end if
      call GetInput(line,nin,xin,ierr)
    end do

    if ((ierr /= 0).or.(nin < 5)) then
      write(unit=0,fmt=*)"ERROR! unable to read reference conditions."
      stop
    end if

    refq(1:nin,j) = xin(1:nin)
    refq(nin+1:8,j) = 0.0

!   set the "motype" (moment type)
    if (nin == 7) then
      refq(9,j) = 2.0    !scale mx,my,mz by refl,refly,reflz
    else if (nin >= 8) then
      refq(9,j) = 3.0    !hinge moment
    else
      refq(9,j) = 1.0    !scale mx,my,mz by refl
    end if

    if (refq(1,j) == 0.0) then
      write(*,*)"REFL is zero. Reset to 1.0."
      refq(1,j) = 1.0
    end if

    if (refq(2,j) == 0.0) then
      write(*,*)"REFA is zero. Reset to total surface ", &
                         "area of component or surface."
    end if

    if (refq(9,j) == 2.0) then
      if (refq(6,j) == 0.0) then
        write(*,*)"REFLY is zero. Reset to 1.0."
        refq(6,j) = 1.0
      end if
      if (refq(7,j) == 0.0) then
        write(*,*)"REFLZ is zero. Reset to 1.0."
        refq(7,j) = 1.0
      end if
    end if

  end do


! read the rest of the data stream twice so all memory can be allocated.
! the first time through, read standard input and echo to scratch.
! the second time through, read from the scratch file.

  do loop = 1,2

    if (loop == 1) then
      if (verbosity >= 1) then
        write(*,*)"Echo data to "//trim(pathname)//"usurp.scratch"
      end if
      iu = 5
      open(unit=3,file=trim(pathname)//"usurp.scratch", &
                  form="formatted", &
                  action="write")
    else
      if (verbosity >= 1) then
        write(*,*)"Reading from "//trim(pathname)//"usurp.scratch"
      end if
      iu = 3
      open(unit=3,file=trim(pathname)//"usurp.scratch", &
                  form="formatted", &
                  action="read")
    end if


!   read the number of integration subsets

    read(unit=iu,fmt=*)num_components
    if (num_components <= 0) then
      write(unit=0,fmt=*)"ERROR! nsurf not positive."
      close(unit=3,status="delete")
      stop
    end if

    if (loop == 1) then

      allocate(tag_list(0:num_components),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! failed to allocate tag_list in PreReadOVERFLOW"
        write(0,*)"num_components = ",num_components
        close(unit=3,status="delete")
        stop
      end if
      tag_list(0)    = "untagged"

      allocate(ncomp_surfs(0:num_components),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! failed to allocate ncomp_surfs in PreReadOVERFLOW"
        write(0,*)"num_components = ",num_components
        close(unit=3,status="delete")
        stop
      end if
      ncomp_surfs = 0

      allocate(refc(1:num_components),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! failed to allocate refc in PreReadOVERFLOW"
        write(0,*)"num_components = ",num_components
        close(unit=3,status="delete")
        stop
      end if
    end if

    if (loop == 1) write(unit=3,fmt=*)num_components

!   keep track of the total number of patches
    ip = 0

    do ic = 1,num_components

      read(unit=iu,fmt=*)nrsub,refc(ic)

      if (nrsub <= 0) then
        write(0,*)"ERROR! nsub must be positive in component ",ic
        write(0,*)"Program terminated."
        stop
      end if

      if ((refc(ic) <= 0).or.(refc(ic) > nref)) then
        write(0,*)"ERROR! Invalid reference set number in component ",ic
        write(0,*)"Program terminated."
        stop
      end if

      if (loop == 1) then
        write(unit=3,fmt=*)nrsub,refc(ic)
        num_patches = num_patches + nrsub
      end if

      do j = 1,nrsub

        ip = ip + 1

        read(unit=iu,fmt=*)n,ibdirs,i1s,i2s,j1s,j2s,k1s,k2s

        if ((n <= 0).or.(n > nb)) then
          write(0,*)"ERROR! Invalid grid number in component ",ic
          write(0,*)"subset ",j
          write(0,*)"grid number (ig) = ",n
          write(0,*)"max number of grids = ",nb
          write(0,*)"Program terminated."
          stop
        end if

        if ((ibdirs < -3).or.(ibdirs == 0).or.(ibdirs > 3)) then
          write(0,*)"ERROR! ibdir has invalid value in component ",ic
          write(0,*)"subset ",j
          write(0,*)"ibdir = ",ibdir
          write(0,*)"ibdir must be +/- 1, 2, or 3."
          write(0,*)"Program terminated."
          stop
        end if

        if (loop == 1) then

          write(unit=3,fmt=*)n,ibdirs,i1s,i2s,j1s,j2s,k1s,k2s

        else

          if (i1s < 0) i1s = nig(n) + 1 + i1s
          if (i2s < 0) i2s = nig(n) + 1 + i2s
          if (j1s < 0) j1s = njg(n) + 1 + j1s
          if (j2s < 0) j2s = njg(n) + 1 + j2s
          if (k1s < 0) k1s = nkg(n) + 1 + k1s
          if (k2s < 0) k2s = nkg(n) + 1 + k2s

!         check subset range
          if ((i1s < 1).or.(i1s > nig(n)).or. &
              (j1s < 1).or.(j1s > njg(n)).or. &
              (k1s < 1).or.(k1s > nkg(n)).or. &
              (i2s < 1).or.(i2s > nig(n)).or. &
              (j2s < 1).or.(j2s > njg(n)).or. &
              (k2s < 1).or.(k2s > nkg(n))) then
            write(0,*)"ERROR! Invalid subset range in component ",ic
            write(0,*)"subset ",j
            if (i1s < 1) then
              write(0,*)"1st dir1 index cannot be less than 1"
            end if
            if (j1s < 1) then
              write(0,*)"1st dir2 index cannot be less than 1"
            end if
            if (k1s < 1) then
              write(0,*)"1st dir3 index cannot be less than 1"
            end if
            if (i1s > nig(n)) then
              write(0,*)"1st dir1 index cannot be greater than ",nig(n)
            end if
            if (j1s > njg(n)) then
              write(0,*)"1st dir2 index cannot be greater than ",njg(n)
            end if
            if (k1s > nkg(n)) then
              write(0,*)"1st dir3 index cannot be greater than ",nkg(n)
            end if
            if (i2s < 1) then
              write(0,*)"2nd dir1 index cannot be less than 1"
            end if
            if (j2s < 1) then
              write(0,*)"2nd dir2 index cannot be less than 1"
            end if
            if (k2s < 1) then
              write(0,*)"2nd dir3 index cannot be less than 1"
            end if
            if (i2s > nig(n)) then
              write(0,*)"2nd dir1 index cannot be greater than ",nig(n)
            end if
            if (j2s > njg(n)) then
              write(0,*)"2nd dir2 index cannot be greater than ",njg(n)
            end if
            if (k2s > nkg(n)) then
              write(0,*)"2nd dir3 index cannot be greater than ",nkg(n)
            end if
            write(0,*)"Program terminated."
            stop
          end if

!         check for i1s > i2s, j1s > j2s, k1s > k2s
          if (i1s > i2s) then
            its = i1s
            i1s = i2s
            i2s = its
          end if
          if (j1s > j2s) then
            its = j1s
            j1s = j2s
            j2s = its
          end if
          if (k1s > k2s) then
            its = k1s
            k1s = k2s
            k2s = its
          end if

!         check consistency of ibdirs and subset ranges
          if (abs(ibdirs) == 1) then
            if (i1s /= i2s) then
              write(0,*)"ERROR! i1s /= i2s for ibdir=1."
              write(0,*)"component, subset = ",ic,j
              write(0,*)"Program terminated."
              stop
            end if
            if (nig(n) > 1) then
              if (((ibdirs == 1).and.(i1s == nig(n))).or. &
                  ((ibdirs == -1).and.(i1s == 1))) then
                write(0,*)"ERROR! Incompatible ibdir,i1s,i2s."
                write(0,*)"component, subset = ",ic,j
                write(0,*)"Program terminated."
                stop
              end if
            else
              if (i1s /= 1) then
                write(0,*)"ERROR! Only 1 slice available in dir1."
                write(0,*)"component, subset = ",ic,j
                write(0,*)"Program terminated."
                stop
              end if
            end if
          end if

          if (abs(ibdirs) == 2) then
            if (j1s /= j2s) then
              write(0,*)"ERROR! j1s /= j2s for ibdir=2."
              write(0,*)"component, subset = ",ic,j
              write(0,*)"Program terminated."
              stop
            end if
            if (njg(n) > 1) then
              if (((ibdirs == 2).and.(j1s == njg(n))).or. &
                  ((ibdirs == -2).and.(j1s == 1))) then
                write(0,*)"ERROR! Incompatible ibdir,j1s,j2s."
                write(0,*)"component, subset = ",ic,j
                write(0,*)"Program terminated."
                stop
              end if
            else
              if (j1s /= 1) then
                write(0,*)"ERROR! Only 1 slice available in dir2."
                write(0,*)"component, subset = ",ic,j
                write(0,*)"Program terminated."
                stop
              end if
            end if
          end if

          if (abs(ibdirs) == 3) then
            if (k1s /= k2s) then
              write(0,*)"ERROR! k1s /= k2s for ibdir=3."
              write(0,*)"component, subset = ",ic,j
              write(0,*)"Program terminated."
              stop
            end if
            if (nkg(n) > 1) then
              if (((ibdirs == 3).and.(k1s == nkg(n))).or. &
                  ((ibdirs == -3).and.(k1s == 1))) then
                write(0,*)"ERROR! Incompatible ibdir,k1s,k2s."
                write(0,*)"component, subset = ",ic,j
                write(0,*)"Program terminated."
                stop
              end if
            else
              if (k1s /= 1) then
                write(0,*)"ERROR! Only 1 slice available in dir3."
                write(0,*)"component, subset = ",ic,j
                write(0,*)"Program terminated."
                stop
              end if
            end if
          end if

          i1(ip) = i1s
          i2(ip) = i2s
          j1(ip) = j1s
          j2(ip) = j2s
          k1(ip) = k1s
          k2(ip) = k2s
          ibdir(ip) = ibdirs
          block(ip) = n

!         use surface index as component number
          icomp(ip) = ic
          ncomp_surfs(ic) = ncomp_surfs(ic) + 1

!         count panels
          call CountStructuredPatch(ip)

        end if
      end do
   
!     read the priority pairs

      read(unit=iu,fmt=*)npri

      if (loop == 1) then
        write(unit=3,fmt=*)npri
        npri_total = npri_total + npri
        do j = 1,npri
          read(unit=iu,fmt=*)nu1,nu2
          if ( (nu1 == 0).or.(abs(nu1) > nrsub).or. &
               (nu2 == 0).or.(abs(nu2) > nrsub) ) then
            write(0,*)"ERROR! subset number out of range."
            write(0,*)"priority pair ",j,":",nu1,nu2
            stop
          else if ( ((nu1 > 0).and.(nu2 < 0)) .or. &
                    ((nu2 > 0).and.(nu1 < 0)) ) then
            write(0,*)"ERROR! Invalid pair. Both must have same sign."
            write(0,*)"priority pair ",j,":",nu1,nu2
            stop
          else
            write(unit=3,fmt=*)nu1,nu2
          end if
        end do
      else
!       note that the patches in this subset have global numbers ranging
!       from ip-nrsub+1 to ip.  To get local number "1" to equal global
!       number ip-nrsub+1, we need to add ip-nrsub to each local number.
        nsub0 = ip - nrsub
        do j = 1,npri
          read(unit=iu,fmt=*)nu1,nu2
          if ((nu1 > 0).and.(nu2 > 0)) then
!           set up normal priority pair
            ipri(nsub0+nu1,nsub0+nu2) = 1
          else
!           set up non-communicating pair
            ipri(nsub0-nu1,nsub0-nu2) = -3
            ipri(nsub0-nu2,nsub0-nu1) = -3
          end if
        end do
      end if
   
    end do


!   read the integration groups

    read(unit=iu,fmt=*)num_groups

    if (loop == 1) then

      write(unit=3,fmt=*)num_groups

      allocate(group_name(num_groups),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(unit=0,fmt=*)"ERROR! allocate of group_name failed."
        write(unit=0,fmt=*)"intended dimensions ",num_groups
        stop
      end if

      allocate(group_size(num_groups),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(unit=0,fmt=*)"ERROR! allocate of group_size failed."
        write(unit=0,fmt=*)"intended dimensions ",num_groups
        stop
      end if

      allocate(group_iref(num_groups),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(unit=0,fmt=*)"ERROR! allocate of group_iref failed."
        write(unit=0,fmt=*)"intended dimensions ",num_groups
        stop
      end if
    end if
   
    do ic = 1,num_groups

      read(unit=iu,fmt=*)group_name(ic)
      read(unit=iu,fmt=*)group_size(ic),group_iref(ic)

      if (loop == 1) then
        write(unit=3,fmt=*)group_name(ic)
        write(unit=3,fmt=*)group_size(ic),group_iref(ic)
        allocate(itmp(group_size(ic))); nalloc = nalloc + 1
        read(unit=iu,fmt=*)(itmp(k),k=1,group_size(ic))
        write(unit=3,fmt=*)(itmp(k),k=1,group_size(ic))
        deallocate(itmp); nalloc = nalloc - 1
      else
        read(unit=iu,fmt=*)(group_data(k,ic),k=1,group_size(ic))
      end if
    end do


!   close the scratch file

    if (loop == 1) then
      close(unit=3,status="keep")
    else
      close(unit=3,status="delete")
    end if


!   allocate memory for subset and group data the first time through

    if (loop == 1) then

      call AllocateStructuredPatch()

      allocate(ibdir(num_patches),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! allocate failed for ibdir"
        stop
      end if

      allocate(group_data(maxval(group_size),num_groups),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! allocate failed for group_data"
        stop
      end if

      allocate(ipri(num_patches,num_patches),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! allocate failed for ipri"
        write(0,*)"num_patches = ",num_patches
        stop
      end if

      ipri = 0

    end if

  end do


! generate a list of component names

  do ic = 1,num_components
    write(tagn,"(i5)")ic
    tag_list(ic) = "surface_"//trim(adjustl(tagn))
  end do


! screen output to confirm results

  write(*,fmt="(a,t41,a1,i7)") &
       "Total number of surface panels","=",num_panels
  write(*,fmt="(a,t41,a1,i7)") &
       "Number of integration surfaces","=",num_components
  write(*,fmt="(a,t41,a1,i7)") &
       "Number of integration subsets","=",num_patches
  write(*,fmt="(a,t41,a1,i7)") &
       "Number of points in integration subsets","=",num_nodes
  write(*,fmt="(a,t41,a1,i7)") &
       "Number of edges in integration subsets","=",num_edges
  if (use_priority) then
    write(*,fmt="(a,t41,a1,i7)") &
       "Number of priority pairs","=",npri_total
  else
    write(*,fmt="(a,t41,a1,i7)") &
       "Number of priority pairs (ignored)","=",npri_total
  end if
  write(*,*)

  deallocate(nkg,njg,nig); nalloc = nalloc - 3

  return
  end subroutine PreReadOVERFLOW

! ====================================================================
  subroutine GridType(iu)
! ====================================================================
! Ascertain (by trial and error) whether grid.in 
! (for FOMOCO / OVERFLOW) is single- or double-precision, 
! whether it is in single-grid or multiple-grid plot3d format, 
! and whether it does or does not include iblank data.
! Rewind the file unit when finished.

  use IntrType ,only: rd,rs
  use NetAlloc ,only: nalloc
  use UserInput,only: multiple_grid,double,ibex

  implicit none

  real(kind=rd),allocatable,dimension(:,:,:) :: xd,yd,zd
  real(kind=rs),allocatable,dimension(:,:,:) :: xs,ys,zs

  integer,allocatable,dimension(:,:,:) :: ib
  integer,allocatable,dimension(:)     :: ni,nj,nk

  integer,intent(in) :: iu
  integer            :: id,jd,kd,ierr,nb,n


  continue


! see if the grid is in single grid format
  read(iu,iostat=ierr)id,jd,kd
  if (ierr /= 0) multiple_grid = .true.

  rewind(iu)

  if (multiple_grid) then
    read(iu)nb
  else
    nb = 1
  end if

  allocate(ni(nb),nj(nb),nk(nb),stat=ierr); nalloc = nalloc + 3
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for ni,nj,nk"
    write(0,*)"nb = ",nb
    stop
  end if

  read(iu)(ni(n),nj(n),nk(n),n=1,nb)

! assume double-precision with iblanks, and try to read the first block
  allocate(xd(ni(1),nj(1),nk(1))); nalloc = nalloc + 1
  allocate(yd(ni(1),nj(1),nk(1))); nalloc = nalloc + 1
  allocate(zd(ni(1),nj(1),nk(1))); nalloc = nalloc + 1
  allocate(ib(ni(1),nj(1),nk(1))); nalloc = nalloc + 1

  read(iu,iostat=ierr)xd,yd,zd,ib
  rewind(iu)

  if (ierr == 0) then
    double = .true.
    ibex   = .true.
    deallocate(nk,nj,ni)   ; nalloc = nalloc - 3
    deallocate(ib,zd,yd,xd); nalloc = nalloc - 4
    return
  end if


! assume double-precision with no iblanks

  if (multiple_grid) read(iu)nb
  read(iu) !dimensions
  read(iu,iostat=ierr)xd,yd,zd
  rewind(iu)

  if (ierr == 0) then
    double = .true.
    ibex   = .false.
    deallocate(nk,nj,ni)   ; nalloc = nalloc - 3
    deallocate(ib,zd,yd,xd); nalloc = nalloc - 4
    return
  end if


! assume single-precision with iblanks
  deallocate(xd,yd,zd)           ; nalloc = nalloc - 3
  allocate(xs(ni(1),nj(1),nk(1))); nalloc = nalloc + 1
  allocate(ys(ni(1),nj(1),nk(1))); nalloc = nalloc + 1
  allocate(zs(ni(1),nj(1),nk(1))); nalloc = nalloc + 1

  if (multiple_grid) read(iu)nb
  read(iu) !dimensions
  read(iu,iostat=ierr)xs,ys,zs,ib
  rewind(iu)

  if (ierr == 0) then
    double = .false.
    ibex   = .true.
    deallocate(nk,nj,ni)   ; nalloc = nalloc - 3
    deallocate(ib,zs,ys,xs); nalloc = nalloc - 4
    return
  end if
  

! assume single-precision with no iblanks
  if (multiple_grid) read(iu)nb
  read(iu) !dimensions
  read(iu,iostat=ierr)xs,ys,zs
  rewind(iu)

  if (ierr == 0) then
    double = .false.
    ibex   = .false.
    deallocate(nk,nj,ni)   ; nalloc = nalloc - 3
    deallocate(ib,zs,ys,xs); nalloc = nalloc - 4
    return
  end if


  write(0,*)"ERROR! failed to ascertain grid format."
  deallocate(nk,nj,ni)   ; nalloc = nalloc - 3
  deallocate(ib,zs,ys,xs); nalloc = nalloc - 4
  stop

  end subroutine GridType

! ==============================================================
  subroutine GetInput(input,nin,x,istat)
! ==============================================================
! Parses the character line input, converts words to real and
! load them into x array.  Character data to the right of valid
! real and integer words is ignored.
!
! If the line contains a '*', indicating cray run-length-encoding,
! it expands the field(s).
!
! WMC 4/15/2003 : added extra test in isanum to test for negative
!                     numbers with no digit in front of decimal, e.g., -.28
! input:
!
! input = string containing line input to be parsed
!
! output:
!
! nin   = number of valid real numbers that were read (integers
!         in input are read as reals)
! x     = array of reals containing numbers parsed from input
! istat = 0   no error (nin >=0)
!       = 1   empty input line
!       = 2   error on input
!
! ==============================================================

  use IntrType,only: rd

  implicit none

  character(len=120) :: rem,sndhlf,rest
  character(len=001) :: comma,sp,star

  character(len=*) :: input
  real(kind=rd) :: x(*)

  integer :: istat,lbeg,lend,l,mbeg,mend,nstar,nrle,lsnd
  integer :: lrest,nfld,nin,nend,ind,len,nstring

! logical :: isanum


  continue


! initialize
  istat = 0
  comma = ','
  sp = ' '
  star = '*'

! return if the line is empty
  call rmsp(input,lbeg,lend)
  if (lbeg  ==  0) then
    istat = 1
    return
  end if

! replace commas with spaces
  do l=lbeg,lend
    if ( input(l:l)  ==  comma ) input(l:l) = sp
  end do

! get rid of run-length-encoding
  do while (index(input,star) /= 0)
    nstar = index(input,star)

!   find the space on the lhs of star field
    mbeg = -1
    do l=nstar,1,-1
      if (input(l:l) == sp) then
        mbeg = l + 1
        exit
      end if
    end do
    if (mbeg == -1) mbeg = 1

!   find the space on the rhs of star field
    mend = -1
    do l=nstar,lend
      if (input(l:l) == sp) then
        mend = l - 1
        exit
      end if
    end do
    if (mend == -1) mend = lend

!   convert string to integer
    read(input(mbeg:nstar-1),*,err=50)nrle

!   now write out second half of field nrle times
    sndhlf = input(nstar+1:mend)
    lsnd = mend - nstar
    rest = input(mend+1:lend)
    lrest = lend - mend
    nend = mbeg-1
    do nfld=1,nrle
      input = input(1:nend-1)//sp//sndhlf(1:lsnd)
      nend = nend + lsnd + 1
    end do

!   append the rest of the line, and go look for more stars
    input = input(1:(nend-1))//rest(1:lrest)
    lend = lend - (mend - mbeg + 1) + nrle*(lsnd+1) - 1
  end do

   50 continue
!  50: error while reading the first half of the * field -- must be
!  something other than run-length encoding, so leave input as it is.

! process input
  nin = 0

  20 continue

  len = lend + 1 - lbeg
  rem = input(lbeg:lend)
  ind = index(rem,sp)
  if (ind  ==  1) then
    lbeg = lbeg + 1
    goto 20
  end if

  nstring = ind-1
  if ( .not. isanum(rem(1:ind-1),nstring) ) goto 100

  read(rem(1:ind-1),*,err=100) x(nin+1)
  nin = nin + 1
  if (ind < len) then
    lbeg = lbeg + ind
    goto 20
  end if

  100 continue

! check status
  if (nin >= 0) then
    istat = 0
    return
  end if

! error on input
  istat = 2

  return
  end subroutine getinput

! ==================================================================
  subroutine rmsp(string,lbeg,lend)
! ==================================================================
! Find beginning and end indices of non-space characters in string.

  implicit none

  integer,intent(out) :: lbeg,lend
  integer :: l,lenst

  character(len=*) :: string
  character(len=1) :: ch,sp,null,tab


  continue


  lenst = len(string)
  lbeg  = 0
  lend  = 0
  sp    = ' '
  null  = char(0)
  tab   = char(9)

! find begining
  do l=1,lenst
    ch = string(l:l)
    if ((ch /= sp).and.(ch /= null) .and.(ch /= tab)) then
      lbeg = l
      exit
    end if
  end do
  if (lbeg  ==  0) return

! find end
  do l=lenst,lbeg,-1
    ch = string(l:l)
    if ((ch /= sp).and.(ch /= null).and.(ch /= tab)) then
      lend = l
      exit
    end if
  end do

  return
  end subroutine rmsp

! =======================================================================
  logical function isanum(istring,nstring)
! =======================================================================
! Determines if the string of nstring characters istring represents a
! number.  If so, isnumber is set to .true..  otherwise, it returns false.
! istring should not include spaces or commas.
 
  implicit none

  integer,intent(in) :: nstring
  integer :: l

  character(len=*),intent(in) :: istring


  continue


! initialize
  isanum = .false.

! return if the line is empty
  if ( nstring  ==  0 ) return 

! first character should be numeric 0-9, a minus sign, or a decimal point
  l = 1
  if (( istring(l:l)  ==  '.' ).and.( nstring > 1 )) then
    l = 2
  else if (( istring(l:l)  ==  '-' ).and.( nstring > 1 )) then
    l = 2
    if (( istring(l:l)  ==  '.' ).and.( nstring > 2 )) l = 3
  end if

  if ( istring(l:l)  ==  '0' .or. &
       istring(l:l)  ==  '1' .or. &
       istring(l:l)  ==  '2' .or. &
       istring(l:l)  ==  '3' .or. &
       istring(l:l)  ==  '4' .or. &
       istring(l:l)  ==  '5' .or. &
       istring(l:l)  ==  '6' .or. &
       istring(l:l)  ==  '7' .or. &
       istring(l:l)  ==  '8' .or. &
       istring(l:l)  ==  '9' ) then
    isanum = .true.
  end if

  return
  end function isanum

  end module PreReadOVERFLOWM
