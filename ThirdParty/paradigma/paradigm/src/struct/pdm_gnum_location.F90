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

module pdm_gnum_location

  use pdm

  implicit none

  interface

    !>
    !!
    !! \brief Build a global numbering location structure
    !!
    !! \param [in]   n_part_in      Number of local partitions for elements
    !! \param [in]   n_part_out     Number of local partitions for requested locations
    !! \param [in]   comm           PDM_MPI communicator
    !! \param [out]  id             Identifier
    !!

    subroutine pdm_gnum_location_create (n_part_in, n_part_out, fComm, id) &
         bind (c, name = 'PDM_gnum_location_create_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: n_part_in
      integer(c_int), value :: n_part_out
      integer(c_int), value :: fComm

      integer(c_int)        :: id

    end subroutine pdm_gnum_location_create


    !>
    !!
    !! \brief Set global numbering
    !!
    !! \param [in]   id          Identifier
    !! \param [in]   i_part_in   Current partition
    !! \param [in]   n_elts_in   Number of elements
    !! \param [in]   gnum_in     Global numbering
    !!

    subroutine pdm_gnum_location_elements_set (id, i_part_in, n_elts_in, gnum_in) &
         bind (c, name = 'PDM_gnum_location_elements_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_part_in
      integer(c_int), value :: n_elts_in

      type(c_ptr), value    :: gnum_in

    end subroutine pdm_gnum_location_elements_set


    !>
    !!
    !! \brief Set requested elements
    !!
    !! \param [in]   id           Identifier
    !! \param [in]   i_part_out   Current partition
    !! \param [in]   n_elts_out   Number of elements
    !! \param [in]   gnum_out     Global numbering
    !!
    !!

    subroutine pdm_gnum_location_requested_elements_set (id, i_part_out, n_elts_out, gnum_out) &
         bind (c, name = 'PDM_gnum_location_requested_elements_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_part_out
      integer(c_int), value :: n_elts_out

      type(c_ptr), value    :: gnum_out

    end subroutine pdm_gnum_location_requested_elements_set


    !>
    !!
    !! \brief Compute the location (processus, partittion, local number in the partition)
    !!
    !! \param [in]   id           Identifier
    !!
    !!

    subroutine pdm_gnum_location_compute (id) &
      bind (c, name = 'PDM_gnum_location_compute')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id

    end subroutine pdm_gnum_location_compute


    !>
    !!
    !! \brief Get localtion
    !!
    !! \param [in]    id             Identifier
    !! \param [in]    i_part_out     Current partition
    !! \param [out]   location_idx   Index in the location arrays (size = 3 !! \ref n_elts + 1)
    !! \param [out]   location       Locations of each element
    !!                                (Three informations : process, partition, element)
    !!

    subroutine pdm_gnum_location_get (id, i_part_out, location_idx, location) &
      bind (c, name = 'PDM_gnum_location_get')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_part_out

      type(c_ptr)           :: location_idx
      type(c_ptr)           :: location

    end subroutine pdm_gnum_location_get


    !>
    !!
    !! \brief Free
    !!
    !! \param [in]   id           Identifier
    !! \param [in]   partial      If partial = 1, location arrays are not freed
    !!
    !!/

    subroutine pdm_gnum_location_free (id, partial) &
      bind (c, name = 'PDM_gnum_location_free')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: partial
    end subroutine pdm_gnum_location_free

  end interface

end module pdm_gnum_location
