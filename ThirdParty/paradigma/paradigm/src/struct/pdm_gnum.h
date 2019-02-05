#ifndef __PDM_GNUM_H__
#define __PDM_GNUM_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   dim          Spatial dimension 
 * \param [in]   n_part       Number of local partitions 
 * \param [in]   merge        Merge double points or not
 * \param [in]   tolerance    Geometric tolerance (if merge double points is activated)
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier    
 */

int
PDM_gnum_create
(
 const int          dim,
 const int          n_part,
 const PDM_bool_t   merge,
 const double       tolerance,       
 const PDM_MPI_Comm comm
);

void
PROCF (pdm_gnum_create, PDM_GNUM_CREATE)
(
 const int *dim,
 const int *n_part,
 const int *merge,  
 const double *tolerance,  
 const PDM_MPI_Fint *fcomm,
       int *id
);


/**
 *
 * \brief Set from coordinates
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
 * \param [in]   char_length  Characteristic length (or NULL) 
 *                            (used if merge double points is activated)
 *
 */

void
PDM_gnum_set_from_coords
(
 const int id,
 const int i_part,
 const int n_elts,
 const double *coords,
 const double *char_length
);

void
PROCF (pdm_gnum_set_from_coords, PDM_GNUM_SET_FROM_COORDS)
(
 const int *id,
 const int *i_part,
 const int *n_elts,
 const double *coords,
 const double *char_length

);


/**
 *
 * \brief Set Parent global numbering
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   parent_gnum  Parent global numbering (size = \ref n_elts)
 *
 */

void
PDM_gnum_set_from_parents
(
 const int id,
 const int i_part,
 const int n_elts,
 const PDM_g_num_t *parent_gnum
);

void
PROCF (pdm_gnum_set_from_parents, PDM_GNUM_SET_FROM_PARENTS)
(
 const int *id,
 const int *i_part,
 const int *n_elts,
 const PDM_g_num_t *parent_gnum
);


/**
 *
 * \brief Compute
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_compute
(
 const int id
);


void
PROCF (pdm_gnum_compute, PDM_GNUM_COMPUTE)
(
 const int *id
);


/**
 *
 * \brief Set from coordinates
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
 *
 */


const PDM_g_num_t *
PDM_gnum_get
(
 const int id,
 const int i_part
);

void
PROCF (pdm_gnum_get, PDM_GNUM_GET)
(
 const int *id,
 const int *i_part,
 PDM_g_num_t *gnum
);


/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_free
(
 const int id,
 const int partial
);

void
PROCF (pdm_gnum_free, PDM_GNUM_FREE)
(
 const int *id,
 const int *partial
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GNUM_H__ */
