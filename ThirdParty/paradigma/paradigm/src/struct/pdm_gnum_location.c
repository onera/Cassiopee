/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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

/*============================================================================
 * Look for the location of a global numbering element. The location has three
 * propertie : process, partition, number of element in this partition.
 * A global numbering can be located in multiple partitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Standard headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_gnum_location.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Local structure definitions
 *============================================================================*/


/**
 * \struct _pdm_gnum_location_t
 * \brief  Define a global numbering location structure
 *
 */

typedef struct  {
  int n_part_in;                       /*!< Number of local partitions  */
  int n_part_out;                      /*!< Number of local partitions
                                            for requested locations */
  int *n_elts_in;                      /*!< Number of elements of each partition */
  const PDM_g_num_t **g_nums_in;       /*!< Global numbering  */
  int *n_elts_out;                     /*!< Number of elements requesting location */
  const PDM_g_num_t **g_nums_out;      /*!< Global numbering of elements requesting location */
  int **location_idx;                  /*!< Location index of elements requesting location */
  int **location;                      /*!< Location of elements requesting location */
  PDM_MPI_Comm comm;                   /*!< Communicator */
} _pdm_gnum_location_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_glocs   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _pdm_gnum_location_t *
_get_from_id
(
 int  id
)
{

  _pdm_gnum_location_t *gloc = (_pdm_gnum_location_t *) PDM_Handles_get (_glocs, id);

  if (gloc == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_gnum_location error : Bad identifier\n");
  }

  return gloc;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/



/**
 *
 * \brief Build a global numbering location structure
 *
 * \param [in]   n_part_in      Number of local partitions for elements
 * \param [in]   n_part_out     Number of local partitions for requested locations
 * \param [in]   comm           PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_gnum_location_create
(
 const int          n_part_in,
 const int          n_part_out,
 const PDM_MPI_Comm comm
)
{
  /*
   * Search a ppart free id
   */

  if (_glocs == NULL) {
    _glocs = PDM_Handles_create (4);
  }

  _pdm_gnum_location_t *_gloc = (_pdm_gnum_location_t *) malloc(sizeof(_pdm_gnum_location_t));
  int id = PDM_Handles_store (_glocs, _gloc);

  _gloc->n_part_in = n_part_in;
  _gloc->n_part_out = n_part_out;

  _gloc->n_elts_in = (int *) malloc (sizeof(int) * n_part_in);
  _gloc->g_nums_in = (const PDM_g_num_t **) malloc (sizeof(const PDM_g_num_t *) * n_part_in);
  for (int i = 0; i < n_part_in; i++) {
    _gloc->g_nums_in[i] = NULL;
  }
  _gloc->n_elts_out = (int *) malloc (sizeof(int) * n_part_out);
  _gloc->g_nums_out = (const PDM_g_num_t **) malloc (sizeof(const PDM_g_num_t *) * n_part_out);
  for (int i = 0; i < n_part_out; i++) {
    _gloc->g_nums_out[i] = NULL;
  }

  _gloc->location_idx = NULL;
  _gloc->location = NULL;
  _gloc->comm = comm;

  return id;

}


void
PDM_gnum_location_create_cf
(
 const int          n_part_in,
 const int          n_part_out,
 const PDM_MPI_Fint comm,
 int *id
)
{
  const PDM_MPI_Comm _comm = PDM_MPI_Comm_f2c(comm);

  *id = PDM_gnum_location_create (n_part_in, n_part_out, _comm);
}

/**
 *
 * \brief Set global numbering
 *
 * \param [in]   id          Identifier
 * \param [in]   i_part_in   Current partition
 * \param [in]   n_elts_in   Number of elements
 * \param [in]   gnum_in     Global numbering
 *
 */

void
PDM_gnum_location_elements_set
(
 const int id,
 const int i_part_in,
 const int n_elts_in,
 const PDM_g_num_t *gnum_in
)
{
  _pdm_gnum_location_t *_gloc = _get_from_id (id);
  _gloc->n_elts_in[i_part_in] = n_elts_in;
  _gloc->g_nums_in[i_part_in] = gnum_in;
}



/**
 *
 * \brief Set requested elements
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part_out   Current partition
 * \param [in]   n_elts_out   Number of elements
 * \param [in]   gnum_out     Global numbering
 *
 */

void
PDM_gnum_location_requested_elements_set
(
 const int id,
 const int i_part_out,
 const int n_elts_out,
 const PDM_g_num_t *gnum_out
)
{
  _pdm_gnum_location_t *_gloc = _get_from_id (id);
  _gloc->n_elts_out[i_part_out] = n_elts_out;
  _gloc->g_nums_out[i_part_out] = gnum_out;
}


/**
 *
 * \brief Compute the location (processus, partittion, local number in the partition)
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_location_compute
(
 const int id
)
{
  _pdm_gnum_location_t *_gloc = _get_from_id (id);

  int rank;
  PDM_MPI_Comm_rank (_gloc->comm, &rank);

  int n_rank;
  PDM_MPI_Comm_size (_gloc->comm, &n_rank);

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1,
                                                       (PDM_g_num_t **) _gloc->g_nums_in,
                                                       NULL,
                                                       _gloc->n_elts_in,
                                                       _gloc->n_part_in,
                                                       _gloc->comm);

  PDM_g_num_t *block_distrib_index = PDM_part_to_block_distrib_index_get (ptb);

  const int s_data = sizeof(int);
  const PDM_stride_t t_stride = PDM_STRIDE_VAR;
  const int cst_stride = 3;

  int  **part_stride = (int **) malloc (sizeof(int *) * _gloc->n_part_in);
  int  **part_data = (int **) malloc (sizeof(int *) * _gloc->n_part_in);

  for (int i = 0; i < _gloc->n_part_in; i++) {
    part_stride[i] = malloc (sizeof(int) * _gloc->n_elts_in[i]);
    part_data[i] = malloc (sizeof(int) * 3 * _gloc->n_elts_in[i]);
    for (int j = 0; j < _gloc->n_elts_in[i]; j++) {
      part_stride[i][j]   = 3;
      part_data[i][3*j]   = rank;
      part_data[i][3*j+1] = i;
      part_data[i][3*j+2] = j+1;
    }
  }

  int  *block_stride = NULL;
  int  *block_data = NULL;

  PDM_part_to_block_exch (ptb,
                          s_data,
                          t_stride,
                          cst_stride,
                          part_stride,
                          (void **) part_data,
                          &block_stride,
                          (void **) &block_data);

  for (int i = 0; i < _gloc->n_part_in; i++) {
    free (part_stride[i]);
    free (part_data[i]);
  }

  free (part_data);
  free (part_stride);

  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_index,
                                                       _gloc->g_nums_out,
                                                       _gloc->n_elts_out,
                                                       _gloc->n_part_out,
                                                       _gloc->comm);

  PDM_block_to_part_exch2 (btp,
                          s_data,
                          PDM_STRIDE_VAR,
                          block_stride,
                          block_data,
                          &part_stride,
                           (void ***) &_gloc->location);

  _gloc->location_idx = (int **) malloc (sizeof(int *) * _gloc->n_part_out);
  for (int i = 0; i < _gloc->n_part_out; i++) {
    _gloc->location_idx[i] = malloc (sizeof(int) * (_gloc->n_elts_out[i] + 1));
    _gloc->location_idx[i][0] = 0;
    for (int j = 0; j < _gloc->n_elts_out[i]; j++) {
      _gloc->location_idx[i][j+1] = _gloc->location_idx[i][j] + part_stride[i][j];;
    }
  }
  free (block_stride);
  free (block_data);

  for (int i = 0; i < _gloc->n_part_out; i++) {
    free (part_stride[i]);
  }

  free (part_stride);

  PDM_part_to_block_free (ptb);
  PDM_block_to_part_free (btp);

}


/**
 *
 * \brief Get location
 *
 * \param [in]    id             Identifier
 * \param [in]    i_part_out     Current partition
 * \param [out]   location_idx   Index in the location arrays (size = 3 * \ref n_elts + 1)
 * \param [out]   location       Locations of each element
 *                                (Three informations : process, partition, element)
 *
 */

void
PDM_gnum_location_get
(
 const int id,
 const int i_part_out,
       int **location_idx,
       int **location
)
{
  _pdm_gnum_location_t *_gloc = _get_from_id (id);
  *location_idx = _gloc->location_idx[i_part_out];
  *location     = _gloc->location[i_part_out];
}


/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_location_free
(
 const int id,
 const int partial
)
{
  _pdm_gnum_location_t *_gloc = _get_from_id (id);

  free (_gloc->n_elts_in);
  free (_gloc->g_nums_in);

  free (_gloc->n_elts_out);
  free (_gloc->g_nums_out);

  if (partial != 1) {
    for (int i = 0; i < _gloc->n_part_out; i++) {
      free (_gloc->location_idx[i]);
    }
    free (_gloc->location_idx);

    for (int i = 0; i < _gloc->n_part_out; i++) {
      free (_gloc->location[i]);
    }
    free (_gloc->location);

  }

  free (_gloc);

  PDM_Handles_handle_free (_glocs, id, PDM_FALSE);

  const int n_gloc = PDM_Handles_n_get (_glocs);

  if (n_gloc == 0) {
    _glocs = PDM_Handles_free (_glocs);
  }
}



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
