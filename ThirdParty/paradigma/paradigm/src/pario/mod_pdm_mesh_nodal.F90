!     
! File:   mod_pdm_mesh_nodal.F90
! Author: equemera
!
! Created on July 10, 2017, 1:34 PM
!

MODULE mod_pdm_mesh_nodal
  
  use mod_pdm
  use iso_c_binding
  
  implicit none

  contains

  
!> Create a Mesh nodal structure
!!
!! @param[in]   n_part   Number of partition on the current process
!! @param[out]  idx      New nodal mesh handle
!!
  
  subroutine pdm_mesh_nodal_create (n_part, idx)
    use iso_c_binding

    implicit none
    
    integer, intent(in)  :: n_part
    integer, intent(out) :: idx
    
    interface
      function pdm_mesh_nodal_create_c (n_part) result(idx) bind(c, name='PDM_Mesh_nodal_create')
        use iso_c_binding
        implicit none
        integer(c_int), intent (in), value :: n_part
        integer(c_int)                     :: idx
      end function pdm_mesh_nodal_create_c
    end interface
 
    idx = pdm_mesh_nodal_create_c (n_part)
    
  end subroutine pdm_mesh_nodal_create

  
!> \brief Free partially a nodal mesh structure
!!
!! @param[in]  idx      Nodal mesh handle
!!

  subroutine pdm_mesh_nodal_partial_free (idx)
    use iso_c_binding

    implicit none
    
    integer, intent(in) :: idx
    
    interface
      subroutine pdm_mesh_nodal_partial_free_c (idx) &
        bind(c, name='PDM_Mesh_nodal_partial_free')

        use iso_c_binding

        implicit none

        integer(c_int), intent (in), value :: idx
      end subroutine pdm_mesh_nodal_partial_free_c
    end interface
       
    call pdm_mesh_nodal_partial_free_c (idx)
    
  end subroutine pdm_mesh_nodal_partial_free  

  
!> Free a nodal mesh structure
!!
!! @param[in]  idx      Nodal mesh handle
!!

  subroutine pdm_mesh_nodal_free (idx)
    use iso_c_binding

    implicit none
    
    integer, intent(in), value :: idx
    
    interface
      subroutine pdm_mesh_nodal_free_c (idx) & 
        bind(c, name='PDM_Mesh_nodal_free')
        
        use iso_c_binding

        implicit none

        integer(c_int), intent (in), value :: idx
      end subroutine pdm_mesh_nodal_free_c
    end interface
    
    call pdm_mesh_nodal_free_c (idx)
    
  end subroutine pdm_mesh_nodal_free  


!> Define partition vertices
!!
!! @param[in]  idx      Nodal mesh handle
!! @param[in]  id_part  Partition identifier
!! @param[in]  n_vtx    Number of vertices
!! @param[in]  coords   Interlaced coordinates (size = 3 * \ref n_vtx)
!! @param[in]  numabs   Global numbering
!!

  subroutine pdm_mesh_nodal_coord_set (idx, id_part, n_vtx, coords, numabs)
    use iso_c_binding

    implicit none
    
    integer, intent(in)                 :: idx
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_vtx
    double precision, dimension(*)      :: coords
    integer (pdm_g_num_s), dimension(*) :: numabs
    
    interface
      subroutine pdm_mesh_nodal_coord_set_c(idx, id_part, n_vtx, coords, numabs) &
        bind(c, name='PDM_Mesh_nodal_coord_set')

        use iso_c_binding
        use mod_pdm

        implicit none

        integer (c_int), intent(in), value   :: idx
        integer (c_int), intent(in), value   :: id_part
        integer (c_int), intent(in), value   :: n_vtx

        real (c_double)         :: coords(*)
        integer (kind=pdm_g_num_s)   :: numabs(*)
      end subroutine pdm_mesh_nodal_coord_set_c
    end interface
    
    call  pdm_mesh_nodal_coord_set_c(idx, id_part, n_vtx, coords, numabs)

  end subroutine  
  
END MODULE mod_pdm_mesh_nodal
