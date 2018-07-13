! ================================================================
  module CommandLineM
! ================================================================

  implicit none

  private
  public :: CommandLine

  contains

! ======================================================
  subroutine CommandLine()
! ======================================================
! Deduce command line arguments for USURP.

  use IntrType ,only: rd
  use NetAlloc ,only: nalloc
  use UserInput,only: arg,allow_fringe
  use UserInput,only: basis,checkpanels,checkpatches,colormap,debugcl
  use UserInput,only: default_basis,default_disjoin,disjoin,dotlimit
  use UserInput,only: dumptree,freememory
  use UserInput,only: full_surface,icpanel1,icpanel2,icpatch1,icpatch2
  use UserInput,only: ignore_pinf,ignore_solution,inflate
  use UserInput,only: keep_weights,min_factor,min_size
  use UserInput,only: never_skip_clip,never_skip_tri,notree
  use UserInput,only: pathname,plotpanel,plotvel,prefix,showcpu
  use UserInput,only: tecbinary,trap_clip,trap_tri,ttype,use_map
  use UserInput,only: use_priority,verbosity,watertight,workdir

  implicit none

  integer            :: i,j,iargc,narg,olen
  character(len=132) :: option
  character(len= 11) :: cfmt
  logical            :: help


  continue


! initialization
  help = .false.
  olen = len(option)
  cfmt = "(000(1x,a))"


! hardwire this for now
! dotlimit = cos(90.0_rd/180.0_rd*acos(-1.0_rd))
! dotlimit = cos(75.0_rd/180.0_rd*acos(-1.0_rd))
  dotlimit = cos(78.0_rd/180.0_rd*acos(-1.0_rd))


#ifdef USE_TECPLOT
  tecbinary = .true.
#else
  tecbinary = .false.
#endif


! gather the command line arguments
  narg = iargc()
  allocate(arg(0:narg)); nalloc = nalloc + 1
  do i = 0,narg
    call getarg(i,arg(i))
  end do


! echo the command line arguments
  write(cfmt(2:4),"(i3.3)")narg+2
  write(0,cfmt)"command:",(trim(arg(i)),i=0,narg)


! parse the command line arguments
  do i = 1,narg
    option = arg(i)

    if (index(option,"--allow-fringe") == 1) then
      allow_fringe = .true.

    else if (index(option,"--basis") == 1) then
      j = index(option,"=")+1
      if (index(option,"patch") == j) then
        basis = "patch"
        default_basis = .false.
      else if (index(option,"panel") == j) then
        basis = "panel"
        default_basis = .false.
      else
        write(0,*)
        write(0,*)"ERROR! invalid command line option: ",trim(option)
        help = .true.
      end if

    else if (index(option,"--check-panels") == 1) then
      checkpanels = .true.
      j = index(option,"=")
      if (j == 0) then
        write(0,*)"ERROR! panel number arguments missing for --check-panels"
        stop
      end if
      read(option(j+1:olen),*)icpanel1,icpanel2

    else if (index(option,"--check-patches") == 1) then
      checkpatches = .true.
      j = index(option,"=")
      if (j == 0) then
        write(0,*)"ERROR! patch number arguments missing for --check-patches"
        stop
      end if
      read(option(j+1:olen),*)icpatch1,icpatch2

    else if (index(option,"--color") == 1) then
      colormap = .true.

    else if (index(option,"--debug") == 1) then
      debugcl = .true.

    else if (index(option,"--disjoin") == 1) then
      if ((index(option,"=") == 0).or.(index(option,"=yes") /= 0)) then
        disjoin = .true.
        default_disjoin = .false.
      else if (index(option,"=no") /= 0) then
        disjoin = .false.
        default_disjoin = .false.
      else
        write(0,*)
        write(0,*)"ERROR! invalid command line option: ",trim(option)
        help = .true.
      end if

    else if (index(option,"--dumptree") == 1) then
      dumptree = .true.

    else if (index(option,"--free-memory") == 1) then
      freememory = .true.

    else if (index(option,"--full-surface") == 1) then
      full_surface = .true.

    else if (index(option,"--help") == 1) then
      write(0,*)
      write(0,*)"stopping code: --help option detected on command line."
      help = .true.

    else if (index(option,"--ignore-solution") == 1) then
      ignore_solution = .true.

    else if (index(option,"--inflate") == 1) then
      j = index(option,"=")
      if (j == 0) then
        write(0,*)
        write(0,*)"ERROR!  a size must be given for the --inflate option."
        stop
      else
        read(option(j+1:olen),*)inflate
      end if

    else if (index(option,"--input") == 1) then
      if (index(option,"~") /= 0) then
        write(0,*)
        write(0,*)"ERROR!  can not interpret tildes in pathnames."
        stop
      end if
      j = index(option,"=")+1
      workdir = trim(option(j:olen))//'/'
      j = len(trim(workdir))
      if (workdir(j:j) /= '/') workdir = trim(workdir)//'/'

    else if (index(option,"--int") == 1) then
      j = index(option,"=")
      if (j == 0) then
        write(0,*)
        write(0,*)"ERROR! p or cp must be specified for the --int option."
        stop
      else if (option(j+1:j+1) == "p") then
        ignore_pinf = .true.
      else if (option(j+1:j+2) == "cp") then
        ignore_pinf = .false.
      else
        write(0,*)
        write(0,*)"ERROR! Invalid option. p or cp must be ", &
                  "specified for the --int option."
        stop
      end if

    else if (index(option,"--keep-weights") == 1) then
      keep_weights = .true.

    else if (index(option,"--min-factor") == 1) then
      j = index(option,"=")
      if (j == 0) then
        write(0,*)
        write(0,*)"ERROR!  a size must be given for the --min-factor option."
        stop
      else
        read(option(j+1:olen),*)min_factor
      end if

    else if (index(option,"--min-size") == 1) then
      j = index(option,"=")
      if (j == 0) then
        write(0,*)
        write(0,*)"ERROR!  a size must be given for the --min-size option."
        stop
      else
        read(option(j+1:olen),*)min_size
      end if

    else if (index(option,"--never-skip") == 1) then
      j = index(option,"=")
      if (j == 0) then
        never_skip_clip = .true.
        never_skip_tri  = .true.
      else
        if (index(option,"clip") > j) then
          never_skip_clip = .true.
        end if
        if (index(option,"tri") > j) then
          never_skip_tri  = .true.
        end if
      end if

    else if (index(option,"--notree") == 1) then
      notree = .true.

    else if (index(option,"--output") == 1) then
      if (index(option,"~") /= 0) then
        write(0,*)
        write(0,*)"ERROR!  can not interpret tildes in pathnames."
        stop
      end if
      j = index(option,"=")+1
      pathname = trim(option(j:olen))
      j = len(trim(pathname))
      if (pathname(j:j) /= '/') pathname = trim(pathname)//'/'

    else if (index(option,"--plot-velocity") == 1) then
      plotvel = .true.

    else if (index(option,"--prefix") == 1) then
      j = index(option,"=")+1
      if (index(option,"RMS") == j) then
        prefix = "RMS"
      else if (index(option,"RST") == j) then
        prefix = "RST"
      else if (index(option,"SAV") == j) then
        prefix = "SAV"
      else
        write(0,*)
        write(0,*)"ERROR! invalid command line option: ",trim(option)
        help = .true.
      end if

    else if (index(option,"--showcpu") == 1) then
      showcpu = .true.

    else if (index(option,"--show-panel") == 1) then
      j = index(option,"=")+1
      read(option(j:olen),*)plotpanel

    else if (index(option,"-Wl,-T") == 1) then
      if (i /= 1) then
        write(0,*)"ERROR! Lahey Fortran argument -Wl,-T must be ", &
                  "the first command line argument."
        stop
      end if

    else if (index(option,"--tecformat") == 1) then
      j = index(option,"=")+1
      if (index(option,"ascii") == j) then
        tecbinary = .false.
      else if (index(option,"binary") == j) then
#ifdef USE_TECPLOT
        tecbinary = .true.
#else
        write(0,*)
        write(0,*)"ERROR! executable was not linked to the tecplot ", &
                  "library during compilation, so the option ", &
                  "--tecformat=binary is invalid. ", &
                  "ensure /home/local/tecplot/lib/tecio.a is included ", &
                  "in the makefile and that -DUSE_TECPLOT is set."
        stop
#endif
      else
        write(0,*)
        write(0,*)"ERROR! invalid qualifier for --tecformat ", &
                  "option detected on command line."
        help = .true.
      end if

    else if (index(option,"--trap") == 1) then
      j = index(option,"=")
      if (j == 0) then
        trap_clip = .true.
        trap_tri  = .true.
      else
        if (index(option,"clip") > j) then
          trap_clip = .true.
        end if
        if (index(option,"tri") > j) then
          trap_tri  = .true.
        end if
      end if

    else if (index(option,"--ttype") == 1) then
      j = index(option,"=")+1
      if (index(option,"GPC") == j) then
        ttype = "GPC"
      else if (index(option,"Triangle") == j) then
        ttype = "Triangle"
      else if (index(option,"none") == j) then
        ttype = "none"
      else
        write(0,*)
        write(0,*)"ERROR! ttype must be GPC, Triangle, or none"
        help = .true.
      end if

    else if (index(option,"--use-priority-pairs") == 1) then
       use_priority = .true.

    else if (index(option,"--use-map") == 1) then
       use_map = .true.
       keep_weights = .true.
       full_surface = .true.

    else if (index(option,"--verbose") == 1) then
      j = index(option,"=")
      if (j == 0) then
        verbosity = 1
      else
        read(option(j+1:olen),*)verbosity
        if ((verbosity < 0).or.(verbosity > 2)) then
          write(0,*)
          write(0,*)"ERROR! index out of range:",verbosity
          write(0,*)"Valid range for verbosity is 0-2"
          stop
        end if
      end if
      if (verbosity == 0) then
        showcpu = .false.
      else
        showcpu = .true.
      end if

    else if (option == "-v") then
      verbosity = 1
      showcpu = .true.

    else if (option == "-vv") then
      verbosity = 2
      showcpu = .true.

    else if (index(option,"--watertight") == 1) then
      watertight = .true.
      full_surface = .true.

    else
      write(0,*)
      write(0,*)"ERROR! invalid command line option: ",trim(option)
      help = .true.

    end if


!   check for inconsistencies
    if ((watertight).and.(ttype /= "Triangle")) then
      write(0,*)
      write(0,*)"ERROR! --watertight is incompatible with --ttype=" &
                // trim(ttype)
      write(0,*)
      stop
    end if

    if ((full_surface).and.(ttype == "none")) then
      write(0,*)
      write(0,*)"ERROR! --full-surface is incompatible with --ttype=" &
                // trim(ttype)
      write(0,*)
      stop
    end if


    if (help) then
      write(0,*)
      write(0,*)"usage: usurp ",                 &
                "[--allow-fringe] ",             &
                "[--basis=<panel,patch>] ",      &
                "[--check-panels=<p1,p2>] ",     &
                "[--check-patches=<p1,p2>] ",    &
                "[--color] ",                    &
                "[--debug] ",                    &
                "[--disjoin=<yes,no>] ",         &
                "[--dumptree] ",                 &
                "[--free-memory] ",              &
                "[--full-surface] ",             &
                "[--help] ",                     &
                "[--ignore-solution] ",          &
                "[--inflate=<eps>] ",            &
                "[--input=<path>] ",             &
                "[--int=<cp,p>] ",               &
                "[--keep-weights] ",             &
                "[--min-factor=<eps>] ",         &
                "[--min-size=<eps>] ",           &
                "[--never-skip[=<clip,tri>]] ",  &
                "[--notree] ",                   &
                "[--output=<path>] ",            &
                "[--plot-velocity] ",            &
                "[--prefix=<RMS,RST,SAV>] ",     &
                "[--showcpu] ",                  &
                "[--show-panel=<panel>] ",       &
                "[--tecformat=<ascii,binary>] ", &
                "[--trap[=<clip,tri>]] ",        &
                "[--ttype=<GPC,Triangle,none>] ",&
                "[--use-map] ",                  &
                "[--use-priority-pairs] ",       &
                "[--verbose[=index] (or -v)] ",  &
                "[--watertight]"

      write(0,*)
      write(0,*)" --allow-fringe           do not blank surface fringe "
      write(0,*)" --basis=<panel,patch>    choose basis for ",                &
                                           "polygon ranking (panel or patch)"
      write(0,*)" --check-panels=<p1,p2>   explain how two panels interact "
      write(0,*)" --check-patches=<p1,p2>  explain how two patches interact ",&
                                           "(for --basis=patch)"
      write(0,*)" --color                  apply a 6-color graph coloring ",  &
                                           "heuristic"
      write(0,*)" --debug                  display verbose debugging ",       &
                                           "information"
      write(0,*)" --disjoin=<yes,no>       specify whether components ",      &
                                           "interact (yes or no)"
      write(0,*)" --dumptree               print the R-tree to a (possibly ", &
                                           "large) Tecplot file"
      write(0,*)" --free-memory            free up memory whenever possible"
      write(0,*)" --full-surface           triangulate entire surface ", &
                                           "upon output"
      write(0,*)" --help                   display this message and stop ",   &
                                           "the code"
      write(0,*)" --ignore-solution        do not read solution files even ", &
                                           "if they exist"
      write(0,*)" --inflate=<eps>          inflate panel bounding boxes by eps" 
      write(0,*)" --input=<path>           reap input files from working ",   &
                                           "directory path"
      write(0,*)" --int=<cp,p>             integrate pressures using cp or ", &
                                           "p (overflow only)"
      write(0,*)" --keep-weights           use weights from panel_weights.dat"
      write(0,*)" --min-factor=<eps>       set min dimension of panel ",  &
                                           "bboxes to eps*size"
      write(0,*)" --min-size=<eps>         set minimum dimension of panel ",  &
                                           "bounding boxes to eps"
      write(0,*)" --never-skip[=tri,clip]  never skip clipping and/or ",      &
                                           "triangulations"
      write(0,*)" --notree                 disable all use of the R-tree"
      write(0,*)" --output=<path>          redirect output to path"
      write(0,*)" --plot-velocity          include u,v,w in usurp-surfaces ", &
                                           "output file"
      write(0,*)" --prefix=<RMS,RST,SAV>   specify solution file type ", &
                                           "(over-rel only)"
      write(0,*)" --showcpu                display cpu time for each step ",  &
                                           "of the code"
      write(0,*)" --show-panel=<panel>     plot the evolution of a ",         &
                                           "specified panel"
      write(0,*)" --tecformat=opt          output format, opt is ascii ",     &
                                           "or binary (default)"
      write(0,*)" --trap[=tri,clip]        trap operations upon code crash"
      write(0,*)" --ttype=opt              triangle type (GPC, Triangle, ",   &
                                           "or none)"
      write(0,*)" --verbose[=index]        (or -v) provide extra ", &
                                           "diagnostic info to STDERR"
      write(0,*)" --use-map                create new grid.i.triq file ",      &
                                           "using usurp.map"
      write(0,*)" --use-priority-pairs     over-rule rankings with set ",      &
                                           "priority pairs"
      write(0,*)" --watertight             create a watertight triangulated ", &
                                           "surface "
      write(0,*)
      stop
    end if

  end do

  deallocate(arg); nalloc = nalloc - 1

  return
  end subroutine CommandLine

  end module CommandLineM
