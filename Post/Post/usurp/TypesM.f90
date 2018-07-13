! =====================================================================
  module Types
! =====================================================================
! polygon data structure, patterned after suggestions by Larry Wigton

  use IntrType,only: rd

  implicit none
  save

  type,public :: contour
    real(kind=rd),pointer,dimension (:) :: x
    real(kind=rd),pointer,dimension (:) :: y
    integer,pointer,dimension (:)       :: node
    integer,pointer,dimension (:)       :: edge
    integer,pointer,dimension (:)       :: ie1
    integer,pointer,dimension (:)       :: ie2
    integer :: num_vert
  end type contour

  type,public :: polygon
    type(contour),pointer,dimension (:) :: cont
    integer,pointer,dimension (:)       :: hole
    integer :: num_cont
  end type  polygon

  type,public :: bcpatch
    real(kind=rd),dimension(3,2)            :: box
    integer                                 :: first_panel,last_panel
    integer                                 :: first_node,last_node
    integer                                 :: first_edge,last_edge
  end type bcpatch

  type(bcpatch),dimension(:),allocatable,public :: surf

  real(kind=rd),allocatable,dimension(:,:) :: scal

  type,public :: panel

    type(polygon)                 :: poly  ! hold subgrid after overlaps

    real(kind=rd),dimension (3,3) :: g     ! rotation matrix which maps xn
                                           ! to (/0,0,1/)
    real(kind=rd),dimension (3,2) :: bound ! reshape(
                                           ! (/xmin,ymin,zmin,xmax,ymax,zmax/),
                                           ! (/3,2/))
    real(kind=rd),dimension (3)   :: xn    ! unit normal (to panel plane)
    real(kind=rd),dimension (3)   :: xmid  ! (x(:,1)+x(:,2)+x(:,3)+x(:,4))/4

    real(kind=rd)     :: area_ref          ! area of full panel
    real(kind=rd)     :: ratio             ! final area to initial area

    integer           :: iblank            ! original iblank for this panel
    integer           :: isurf
    integer           :: icomp
    integer           :: panel_number
    integer           :: num_vertices
    integer           :: itri
    integer,dimension(4) :: node
    integer,dimension(4) :: edge

  end type panel


  end module Types
