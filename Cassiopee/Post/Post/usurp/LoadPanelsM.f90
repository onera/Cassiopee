! ================================================================
  module LoadPanelsM
! ================================================================

  implicit none

  private
  public :: LoadPanels

  contains

! =====================================================
  subroutine LoadPanels(p)
! =====================================================
! Load points from structured patches into panels.
! The order of points is chosen to ensure that the
! surface normals calculated later will point into the 
! flow.  Also, note that the dummy argument p is a 
! pointer, so an explicit interface is required.

  use EdgeData    ,only: edge,num_edges
  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use Types       ,only: panel,scal,surf
  use PatchInfo   ,only: i1,i2,j1,j2,k1,k2,num_patches,icomp
  use UserInput   ,only: cpufmt,showcpu

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rs) :: time1,time2
  integer       :: ipan,ip,i,j,k
  integer       :: horiz_edges
  integer       :: south_edge,north_edge
  integer       :: east_edge,west_edge


  continue


  call my_cpu_time(time1)

  allocate(edge(num_edges)); nalloc = nalloc + 1
  do i = 1,num_edges
    edge(i)%left_face = 0
    edge(i)%right_face = 0
  end do


  ipan = 0
  do ip = 1,num_patches

    if (i1(ip) == i2(ip)) then

      horiz_edges = (j2(ip)-j1(ip))*(k2(ip)-k1(ip)+1)

      i = i1(ip)
      do k = k1(ip)+1,k2(ip)
      do j = j1(ip)+1,j2(ip)
        ipan = ipan + 1
        p(ipan)%panel_number = ipan
        p(ipan)%iblank = scal(1,ipan)
        p(ipan)%isurf = ip
        p(ipan)%icomp = icomp(p(ipan)%isurf)
        p(ipan)%num_vertices = 4

        south_edge  = (surf(ip)%first_edge - 1)    &
                    + (ipan - surf(ip)%first_panel + 1)
        north_edge  = south_edge + (j2(ip)-j1(ip))
        west_edge   = south_edge + horiz_edges + (k - k1(ip) - 1)
        east_edge   = west_edge + 1

        p(ipan)%node(1) = (surf(ip)%first_node - 1)         &
                        + (ipan - surf(ip)%first_panel + 1) &
                        + (k - k1(ip) - 1)
        p(ipan)%node(3) = p(ipan)%node(1) + (j2(ip)-j1(ip)+1) + 1
        if (i == 1) then
          p(ipan)%node(2) = p(ipan)%node(1) + 1
          p(ipan)%node(4) = p(ipan)%node(3) - 1
          p(ipan)%edge(1) = south_edge
          p(ipan)%edge(2) = east_edge
          p(ipan)%edge(3) = north_edge
          p(ipan)%edge(4) = west_edge
          edge(east_edge)%vertex(1) = p(ipan)%node(2)
          edge(east_edge)%vertex(2) = p(ipan)%node(3)
          edge(west_edge)%vertex(1) = p(ipan)%node(1)
          edge(west_edge)%vertex(2) = p(ipan)%node(4)
          edge(south_edge)%vertex(1) = p(ipan)%node(1)
          edge(south_edge)%vertex(2) = p(ipan)%node(2)
          edge(north_edge)%vertex(1) = p(ipan)%node(4)
          edge(north_edge)%vertex(2) = p(ipan)%node(3)
        else
          p(ipan)%node(2) = p(ipan)%node(3) - 1
          p(ipan)%node(4) = p(ipan)%node(1) + 1
          p(ipan)%edge(1) = west_edge
          p(ipan)%edge(2) = north_edge
          p(ipan)%edge(3) = east_edge
          p(ipan)%edge(4) = south_edge
          edge(east_edge)%vertex(1) = p(ipan)%node(4)
          edge(east_edge)%vertex(2) = p(ipan)%node(3)
          edge(west_edge)%vertex(1) = p(ipan)%node(1)
          edge(west_edge)%vertex(2) = p(ipan)%node(2)
          edge(south_edge)%vertex(1) = p(ipan)%node(1)
          edge(south_edge)%vertex(2) = p(ipan)%node(4)
          edge(north_edge)%vertex(1) = p(ipan)%node(2)
          edge(north_edge)%vertex(2) = p(ipan)%node(3)
        end if
        edge(east_edge)%left_face = ipan
        edge(west_edge)%right_face = ipan
        edge(south_edge)%left_face = ipan
        edge(north_edge)%right_face = ipan
      end do
      end do

    else if (j1(ip) == j2(ip)) then

      horiz_edges = (i2(ip)-i1(ip))*(k2(ip)-k1(ip)+1)

      j = j1(ip)
      do k = k1(ip)+1,k2(ip)
      do i = i1(ip)+1,i2(ip)
        ipan = ipan + 1
        p(ipan)%panel_number = ipan
        p(ipan)%iblank = scal(1,ipan)
        p(ipan)%isurf = ip
        p(ipan)%icomp = icomp(p(ipan)%isurf)
        p(ipan)%num_vertices = 4

        south_edge  = (surf(ip)%first_edge - 1)    &
                    + (ipan - surf(ip)%first_panel + 1)
        north_edge  = south_edge + (i2(ip)-i1(ip))
        west_edge   = south_edge + horiz_edges + (k - k1(ip) - 1)
        east_edge   = west_edge + 1

        p(ipan)%node(1) = (surf(ip)%first_node - 1)         &
                        + (ipan - surf(ip)%first_panel + 1) &
                        + (k - k1(ip) - 1)
        p(ipan)%node(3) = p(ipan)%node(1) + (i2(ip)-i1(ip)+1) + 1
        if (j == 1) then
          p(ipan)%node(2) = p(ipan)%node(3) - 1
          p(ipan)%node(4) = p(ipan)%node(1) + 1
          p(ipan)%edge(1) = west_edge
          p(ipan)%edge(2) = north_edge
          p(ipan)%edge(3) = east_edge
          p(ipan)%edge(4) = south_edge
          edge(east_edge)%vertex(1) = p(ipan)%node(4)
          edge(east_edge)%vertex(2) = p(ipan)%node(3)
          edge(west_edge)%vertex(1) = p(ipan)%node(1)
          edge(west_edge)%vertex(2) = p(ipan)%node(2)
          edge(south_edge)%vertex(1) = p(ipan)%node(1)
          edge(south_edge)%vertex(2) = p(ipan)%node(4)
          edge(north_edge)%vertex(1) = p(ipan)%node(2)
          edge(north_edge)%vertex(2) = p(ipan)%node(3)
        else
          p(ipan)%node(2) = p(ipan)%node(1) + 1
          p(ipan)%node(4) = p(ipan)%node(3) - 1
          p(ipan)%edge(1) = south_edge
          p(ipan)%edge(2) = east_edge
          p(ipan)%edge(3) = north_edge
          p(ipan)%edge(4) = west_edge
          edge(east_edge)%vertex(1) = p(ipan)%node(2)
          edge(east_edge)%vertex(2) = p(ipan)%node(3)
          edge(west_edge)%vertex(1) = p(ipan)%node(1)
          edge(west_edge)%vertex(2) = p(ipan)%node(4)
          edge(south_edge)%vertex(1) = p(ipan)%node(1)
          edge(south_edge)%vertex(2) = p(ipan)%node(2)
          edge(north_edge)%vertex(1) = p(ipan)%node(4)
          edge(north_edge)%vertex(2) = p(ipan)%node(3)
        end if
        edge(east_edge)%left_face = ipan
        edge(west_edge)%right_face = ipan
        edge(south_edge)%left_face = ipan
        edge(north_edge)%right_face = ipan
      end do
      end do

    else if (k1(ip) == k2(ip)) then

      horiz_edges = (i2(ip)-i1(ip))*(j2(ip)-j1(ip)+1)

      k = k1(ip)
      do j = j1(ip)+1,j2(ip)
      do i = i1(ip)+1,i2(ip)
        ipan = ipan + 1
        p(ipan)%panel_number = ipan
        p(ipan)%iblank = scal(1,ipan)
        p(ipan)%isurf = ip
        p(ipan)%icomp = icomp(p(ipan)%isurf)
        p(ipan)%num_vertices = 4

        south_edge  = (surf(ip)%first_edge - 1)    &
                    + (ipan - surf(ip)%first_panel + 1)
        north_edge  = south_edge + (i2(ip)-i1(ip))
        west_edge   = south_edge + horiz_edges + (j - j1(ip) - 1)
        east_edge   = west_edge + 1

        p(ipan)%node(1) = (surf(ip)%first_node - 1)         &
                        + (ipan - surf(ip)%first_panel + 1) &
                        + (j - j1(ip) - 1)
        p(ipan)%node(3) = p(ipan)%node(1) + (i2(ip)-i1(ip)+1) + 1
        if (k == 1) then
          p(ipan)%node(2) = p(ipan)%node(1) + 1
          p(ipan)%node(4) = p(ipan)%node(3) - 1
          p(ipan)%edge(1) = south_edge
          p(ipan)%edge(2) = east_edge
          p(ipan)%edge(3) = north_edge
          p(ipan)%edge(4) = west_edge
          edge(east_edge)%vertex(1) = p(ipan)%node(2)
          edge(east_edge)%vertex(2) = p(ipan)%node(3)
          edge(west_edge)%vertex(1) = p(ipan)%node(1)
          edge(west_edge)%vertex(2) = p(ipan)%node(4)
          edge(south_edge)%vertex(1) = p(ipan)%node(1)
          edge(south_edge)%vertex(2) = p(ipan)%node(2)
          edge(north_edge)%vertex(1) = p(ipan)%node(4)
          edge(north_edge)%vertex(2) = p(ipan)%node(3)
        else
          p(ipan)%node(2) = p(ipan)%node(3) - 1
          p(ipan)%node(4) = p(ipan)%node(1) + 1
          p(ipan)%edge(1) = west_edge
          p(ipan)%edge(2) = north_edge
          p(ipan)%edge(3) = east_edge
          p(ipan)%edge(4) = south_edge
          edge(east_edge)%vertex(1) = p(ipan)%node(4)
          edge(east_edge)%vertex(2) = p(ipan)%node(3)
          edge(west_edge)%vertex(1) = p(ipan)%node(1)
          edge(west_edge)%vertex(2) = p(ipan)%node(2)
          edge(south_edge)%vertex(1) = p(ipan)%node(1)
          edge(south_edge)%vertex(2) = p(ipan)%node(4)
          edge(north_edge)%vertex(1) = p(ipan)%node(2)
          edge(north_edge)%vertex(2) = p(ipan)%node(3)
        end if
        edge(east_edge)%left_face = ipan
        edge(west_edge)%right_face = ipan
        edge(south_edge)%left_face = ipan
        edge(north_edge)%right_face = ipan

      end do
      end do
    end if
  end do

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in LoadPanels:',time2-time1
  end if

  return
  end subroutine LoadPanels

  end module LoadPanelsM
