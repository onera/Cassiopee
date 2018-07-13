! =====================================================================
  module UserInput
! =====================================================================
! Module to carry parameters associated with the different CFD work
! environments and with command line options provided during the
! execution of USURP.

  use IntrType,only: rd,rs

  implicit none
  save

  real(kind=rd),allocatable,dimension(:,:) :: refq

  real(kind=rd)     ,public :: fsmach,alpha,beta,rey,gaminf,tinf
  real(kind=rd)     ,public :: merge_tol,merge_tol2
  real(kind=rd)     ,public :: degen_tol,degen_tol2
  real(kind=rd)     ,public :: dotlimit

  real(kind=rd)     ,public :: inflate    = 1.00_rd
  real(kind=rd)     ,public :: min_factor = 0.01_rd
  real(kind=rd)     ,public :: min_size   = 0.00_rd
  real(kind=rd)     ,public :: qinf       = 0.50_rd
  real(kind=rd)     ,public :: rhoinf     = 1.00_rd
  real(kind=rd)     ,public :: vminf      = 1.00_rd

  real(kind=rs)     ,public :: refa  = 1.0_rs
  real(kind=rs)     ,public :: refl  = 1.0_rs
  real(kind=rs)     ,public :: rel   = 0.0_rs
  real(kind=rs)     ,public :: rough = 0.0_rs
  real(kind=rs)     ,public :: xcg   = 0.0_rs
  real(kind=rs)     ,public :: ycg   = 0.0_rs
  real(kind=rs)     ,public :: zcg   = 0.0_rs

  real(kind=rs),dimension(3),public :: gpcf_time

  integer,allocatable,dimension(:),public :: ltri

  integer           ,public :: iwfn      = 0
  integer           ,public :: plotpanel = 0
  integer           ,public :: verbosity = 0

  integer           ,public :: icpanel1,icpanel2
  integer           ,public :: icpatch1,icpatch2
  integer           ,public :: my_bit_size
  integer           ,public :: nref

  character(len=132),public,allocatable,dimension(:) :: arg

  character(len=80) ,public :: cpufmt        = "(1x,a,t40,f7.2)"
  character(len=80) ,public :: dci_file_name = "error"
  character(len=80) ,public :: pathname      = './'
  character(len=80) ,public :: workdir       = './'
  character(len=30) ,public :: fgrid         = "error"
  character(len=30) ,public :: fnameo        = "error"
  character(len=08) ,public :: ttype         = "Triangle"
  character(len=05) ,public :: basis         = "error"
  character(len=03) ,public :: prefix        = "RST"

  logical           ,public :: double          = .false.
  logical           ,public :: ibex            = .false.
  logical           ,public :: multiple_grid   = .false.
  logical           ,public :: overset         = .false.
  logical           ,public :: solution_exists = .false.

  logical           ,public :: allow_fringe    = .false.
  logical           ,public :: cfd_ship        = .false.
  logical           ,public :: checkpanels     = .false.
  logical           ,public :: checkpatches    = .false.
  logical           ,public :: colormap        = .false.
  logical           ,public :: debugcl         = .false.
  logical           ,public :: disjoin         = .false.
  logical           ,public :: default_disjoin = .true.
  logical           ,public :: default_basis   = .true.
  logical           ,public :: dplr            = .false.
  logical           ,public :: dumptree        = .false.
  logical           ,public :: freememory      = .false.
  logical           ,public :: full_surface    = .false.
  logical           ,public :: generic         = .false.
  logical           ,public :: ignore_solution = .false.
  logical           ,public :: ignore_pinf     = .false.
  logical           ,public :: keep_weights    = .false.
  logical           ,public :: never_skip_clip = .false.
  logical           ,public :: never_skip_tri  = .false.
  logical           ,public :: notree          = .false.
  logical           ,public :: nphase          = .false.
  logical           ,public :: over_rel        = .false.
  logical           ,public :: overflow        = .false.
  logical           ,public :: plotvel         = .false.
  logical           ,public :: showcpu         = .false.
  logical           ,public :: tecbinary       = .false.
  logical           ,public :: trap_clip       = .false.
  logical           ,public :: trap_tri        = .false.
  logical           ,public :: unclem          = .false.
  logical           ,public :: use_map         = .false.
  logical           ,public :: use_priority    = .false.
  logical           ,public :: watertight      = .false.

  end module UserInput
