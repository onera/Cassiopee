!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2019  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_dist_cellcenter_surf

  use pdm

  implicit none


  interface

    !> \brief Create a structure to compute the distance from cell centers
    !!        to a surface mesh
    !!
    !! \param [in]    comm           MPI communicator
    !! \param [inout] id             Identifier
    !!

    subroutine pdm_dist_cellcenter_surf_create (fComm, id) &
         bind (c, name = 'PDM_dist_cellcenter_surf_create_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: fComm

      integer(c_int)        :: id

    end subroutine pdm_dist_cellcenter_surf_create

    !> \brief Set global data of a surface mesh
    !!
    !! \param [in]   id             Identifier
    !! \param [in]   n_g_face       Global number of faces
    !! \param [in]   n_g_vtx        Global number of vertices
    !! \param [in]   n_part         Number of partition
    !!

    subroutine pdm_dist_cellcenter_surf_surf_mesh_global_data_set (id, n_g_face, n_g_vtx, n_part) &
         bind (c, name = 'PDM_dist_cellcenter_surf_surf_mesh_global_data_set')

      use iso_c_binding

      implicit none

      integer(c_int), value        :: id
#ifdef PDM_LONG_G_NUM
      integer(c_long), value       :: n_g_face
      integer(c_long), value       :: n_g_vtx
#else
      integer(c_int), value        :: n_g_face
      integer(c_int), value        :: n_g_vtx
#endif
      integer(c_int), value        :: n_part
    end subroutine pdm_dist_cellcenter_surf_surf_mesh_global_data_set

    !> \brief Set a part of a surface mesh
    !!
    !! \param [in]   id            Identifier
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_face        Number of faces
    !! \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
    !! \param [in]   face_vtx      face -> vertex connectivity
    !! \param [in]   face_ln_to_gn Local face numbering to global face numbering
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!

    subroutine pdm_dist_cellcenter_surf_surf_mesh_part_set (id, i_part, n_face, face_vtx_idx, &
                                                 face_vtx, face_ln_to_gn, n_vtx, coords, &
                                                 vtx_ln_to_gn) &
      bind (c, name = 'PDM_dist_cellcenter_surf_surf_mesh_part_set')
      use iso_c_binding

      implicit none

      integer(c_int), value     :: id
      integer(c_int), value     :: i_part
      integer(c_int), value     :: n_face
      type(c_ptr), value        :: face_vtx_idx
      type(c_ptr), value        :: face_vtx
      type(c_ptr), value        :: face_ln_to_gn
      integer(c_int), value     :: n_vtx
      type(c_ptr), value        :: coords
      type(c_ptr), value        :: vtx_ln_to_gn

    end subroutine pdm_dist_cellcenter_surf_surf_mesh_part_set

    !>
    !!
    !! \brief Set global data of a surface mesh
    !!
    !! \param [in]   id             Identifier
    !! \param [in]   n_g_face       Global number of faces
    !! \param [in]   n_g_vtx        Global number of vertices
    !! \param [in]   n_part         Number of partition
    !!
    !!

    subroutine pdm_dist_cellcenter_surf_vol_mesh_global_data_set (id, n_g_cell, n_g_face, n_g_vtx, &
                                                                  n_part) &

      bind (c, name = 'PDM_dist_cellcenter_surf_vol_mesh_global_data_set')
      use iso_c_binding

      implicit none
      integer(c_int), value        :: id
#ifdef PDM_LONG_G_NUM
      integer(c_long), value       :: n_g_cell
      integer(c_long), value       :: n_g_face
      integer(c_long), value       :: n_g_vtx
#else
      integer(c_int), value        :: n_g_cell
      integer(c_int), value        :: n_g_face
      integer(c_int), value        :: n_g_vtx
#endif
      integer(c_int), value        :: n_part
    end subroutine pdm_dist_cellcenter_surf_vol_mesh_global_data_set

    !!
    !! \brief Set a part of a surface mesh
    !!
    !! \param [in]   id            Identifier
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_cell        Number of cells
    !! \param [in]   cell_face_idx Cell -> face connectivity index
    !! \param [in]   cell_face     Cell -> face connectivity
    !! \param [in]   cell_center   Cell center or NULL
    !! \param [in]   cell_ln_to_gn Local cell numbering to global cell numbering
    !! \param [in]   n_face        Number of faces
    !! \param [in]   face_vtx_idx  Face -> vtx connectivity index
    !! \param [in]   face_vtx      Face -> vtx connectivity
    !! \param [in]   face_cell     face -> cell   connectivity
    !! \param [in]   face_ln_to_gn Local face numbering to global face numbering
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!

    subroutine pdm_dist_cellcenter_surf_vol_mesh_part_set (id, &
                                                           i_part, &
                                                           n_cell, &
                                                           cell_face_idx, &
                                                           cell_face, &
                                                           cell_center, &
                                                           cell_ln_to_gn, &
                                                           n_face, &
                                                           face_vtx_idx, &
                                                           face_vtx, &
                                                           face_ln_to_gn, &
                                                           n_vtx, &
                                                           coords, &
                                                           vtx_ln_to_gn)&

      bind (c, name = 'PDM_dist_cellcenter_surf_vol_mesh_part_set')
      use iso_c_binding

      implicit none
      integer(c_int), value  :: id
      integer(c_int), value  :: i_part
      integer(c_int), value  :: n_cell
      type(c_ptr),    value  :: cell_face_idx
      type(c_ptr),    value  :: cell_face
      type(c_ptr),    value  :: cell_center
      type(c_ptr),    value  :: cell_ln_to_gn
      integer(c_int), value  :: n_face
      type(c_ptr),    value  :: face_vtx_idx
      type(c_ptr),    value  :: face_vtx
      type(c_ptr),    value  :: face_ln_to_gn
      integer(c_int), value  :: n_vtx
      type(c_ptr),    value  :: coords
      type(c_ptr),    value  :: vtx_ln_to_gn

    end subroutine pdm_dist_cellcenter_surf_vol_mesh_part_set

    !!> \brief Compute distance
    !!
    !! \param [in]   id  Identifier
    !!

    subroutine pdm_dist_cellcenter_surf_compute (id) &
         bind (c, name = 'PDM_dist_cellcenter_surf_compute')
      use iso_c_binding

      implicit none

      integer(c_int), value     :: id
    end subroutine pdm_dist_cellcenter_surf_compute

    !> \brief Get mesh distance
    !!
    !! \param [in]   id                    Identifier
    !! \param [in]   i_point_cloud         Current cloud
    !! \param [in]   i_part                Index of partition of the cloud
    !! \param [out]  closest_elt_distance  Distance
    !! \param [out]  closest_elt_projected Projected point coordinates
    !! \param [out]  closest_elt_g_num     Global number of the closest element
    !!

    subroutine pdm_dist_cellcenter_surf_get (id, i_part, &
                                            closest_elt_distance, &
                                            closest_elt_projected, &
                                            closest_elt_gnum) &
      bind (c, name = 'PDM_dist_cellcenter_surf_get')

      use iso_c_binding

      implicit none

      integer(c_int), value     :: id
      integer(c_int), value     :: i_part
      type(c_ptr)               :: closest_elt_distance
      type(c_ptr)               :: closest_elt_projected
      type(c_ptr)               :: closest_elt_gnum

    end subroutine pdm_dist_cellcenter_surf_get

    !> \brief Free a distance mesh structure
    !!
    !! \param [in]  id       Identifier
    !! \param [in]  partial  if partial is equal to 0, all data are removed.
    !!                       Otherwise, results are kept.
    !!

    subroutine pdm_dist_cellcenter_surf_free (id, partial) &
       bind (c, name = 'PDM_dist_cellcenter_surf_free')

      use iso_c_binding

      implicit none

      integer(c_int), value     :: id
      integer(c_int), value     :: partial

    end subroutine pdm_dist_cellcenter_surf_free

    !> \brief  Dump elapsed an CPU time
    !!
    !! \param [in]  id       Identifier
    !!

    subroutine pdm_dist_cellcenter_surf_dump_times (id) &
         bind (c, name = 'PDM_dist_cellcenter_surf_dump_times')

      use iso_c_binding

      implicit none

      integer(c_int), value     :: id
    end subroutine pdm_dist_cellcenter_surf_dump_times

  end interface

end module pdm_dist_cellcenter_surf
