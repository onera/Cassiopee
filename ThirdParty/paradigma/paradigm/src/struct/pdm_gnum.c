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

/*============================================================================
 * Main structure for an I/O numbering scheme associated with mesh entities
 * (such as cells, faces, and vertices);
 *
 * In parallel mode, such a scheme is important so as to redistribute
 * locally numbered entities on n processes to files written by p
 * processes, with p <= n.
 *
 * Only the case where p = 1 is presently implemented, so the numbering
 * scheme is simply based on entity's global labels.
 *
 * For p > 1, it would probably be necessary to extend the numbering
 * schemes so as to account for the fact that a given entity may have
 * a main index on its main associated domain, but may be present
 * as a ghost entity with another index on neighboring domains.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_binary_search.h"
#include "pdm_mpi.h"
#include "pdm_points_merge.h"
#include "pdm_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_gnum.h"

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
 * \struct _pdm_gnum_t
 * \brief  Define a global numberring 
 * 
 */

typedef struct  {

  int          n_part;      /*!< Number of partitions */
  int          dim;         /*!< Spatial dimension */
  PDM_bool_t   merge;       /*!< Merge double point status */
  double       tolerance;   /*!< Geometric tolerance */
  PDM_MPI_Comm comm;        /*!< MPI communicator */
  PDM_g_num_t  n_g_elt;     /*!< Global number of elements */
  int          *n_elts;     /*!< Number of elements in partitions */
  PDM_g_num_t **g_nums;     /*!< Global numbering of elements */ 
  double      **coords;     /*!< Coordinates of elements */
  double      **char_length;/*!< Characteristic length */
  int         **index;      /*!< Index : used if merge is activated */
  PDM_g_num_t **parent;     /*!< Global n */ 

} _pdm_gnum_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_gnums   = NULL;

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

static _pdm_gnum_t *
_get_from_id
(
 int  id
)
{
  
  _pdm_gnum_t *gnum = (_pdm_gnum_t *) PDM_Handles_get (_gnums, id);
    
  if (gnum == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_gnum error : Bad identifier\n");
  }

  return gnum;
}


/**
 * \brief Use bubble sort on an expectedly short sequence of coordinates
 * to ensure lexicographical ordering.
 *
 *  \param[in]      dim        <-- spatial dimension
 *  \param[in]      start_id   <-- start id in array
 *  \param[in]      end_id     <-- past-the-end id in array
 *  \param[in]      coords     <-- pointer to entity coordinates (interlaced)
 *  \param[in, out] order      <-> ordering array base on Morton encoding, or
 *                                 lexicographical coordinate ordering for ties
 */

/*
  This function comes from "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

inline static void
_reorder_coords_lexicographic
(
int                dim,
size_t             start_id,
size_t             end_id,
const double       coords[],
PDM_l_num_t        order[]
)
{
  size_t  i;
  bool g_swap;

  do {

    g_swap = false;

    for (i = start_id + 1; i < end_id; i++) {

      size_t j_prev = order[i-1], j = order[i];
      bool l_swap = false;

      if (dim == 3) {
        if (coords[j_prev*3] < coords[j*3])
          continue;
        else if (coords[j_prev*3] > coords[j*3])
          l_swap = true;
        else if (coords[j_prev*3 + 1] < coords[j*3 + 1])
          continue;
        else if (   coords[j_prev*3 + 1] > coords[j*3 + 1]
                 || coords[j_prev*3 + 2] > coords[j*3 + 2])
          l_swap = true;
      }
      else if (dim == 2) {
        if (coords[j_prev*2] < coords[j*2 + 1])
          continue;
        else if (   coords[j_prev*2]     > coords[j*2]
                 || coords[j_prev*2 + 1] > coords[j*2 + 1])
          l_swap = true;
      }
      else { /* if (dim == 1) */
        if (coords[j_prev] > coords[j])
          l_swap = true;
      }

      if (l_swap) {
        PDM_l_num_t o_save = order[i-1];
        order[i-1] = order[i];
        order[i] = o_save;
        g_swap = true;
      }
    }

  } while (g_swap);
}


/**
 * 
 * \brief Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be different, but their order is undetermined.
 *
 *  \param [in]      dim        spatial dimension
 *  \param [in]      n_entities number of entities considered
 *  \param [in]      coords     pointer to entity coordinates (interlaced)
 *  \param [in]      m_code     Morton code associated with each entity
 *  \param [in, out] order      ordering array base on Morton encoding, or
 *                              lexicographical coordinate ordering for ties
 *
 */

/*
  This function comes from "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

static void
_check_morton_ordering
(
int                      dim,
size_t                   n_entities,
const double        coords[],
const PDM_morton_code_t  m_code[],
PDM_l_num_t               order[]
)
{
  size_t  i_prev = 0, i = 1;

  if (n_entities == 0)
    return;

  /* Check ordering; if two entities have the same Morton codes,
     use lexicographical coordinates ordering to ensure the
     final order is deterministic. */

  for (i = 1; i < n_entities; i++) {

    size_t j_prev = order[i_prev], j = order[i];

    if (   m_code[j_prev].X[0] != m_code[j].X[0]
        || m_code[j_prev].X[1] != m_code[j].X[1]
        || m_code[j_prev].X[2] != m_code[j].X[2]) {

      /* If successive values have the same Morton code,
         order them lexicographically */
      if (i_prev < i - 1)
        _reorder_coords_lexicographic(dim, i_prev, i-1, coords, order);

    }
    i_prev = i;
  }

  if (i_prev < n_entities - 1)
    _reorder_coords_lexicographic(dim, i_prev, n_entities - 1, coords, order);
}


/**
 *
 * \brief Compute from coords
 *
 * \param [in]   _gnum          Current _pdm_gnum_t structure
 *
 */

/*
  This function comes from "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

static void
_gnum_from_coords_compute
(
 _pdm_gnum_t *_gnum 
)
{
  double extents[6];
  PDM_l_num_t  *order = NULL;
  PDM_morton_code_t *m_code = NULL;

  PDM_MPI_Comm comm = _gnum->comm;

  const int level = sizeof(PDM_morton_int_t)*8 - 1;
  
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  /* Merge double points */
  
  int id_pm = 0; 
  
  int n_entities = 0;

  int iproc;
  PDM_MPI_Comm_rank (comm, &iproc);

  if (_gnum->merge) {
    
    _gnum->index = malloc (sizeof(int *) * _gnum->n_part);
  
    id_pm = PDM_points_merge_create (_gnum->n_part, _gnum->tolerance, _gnum->comm);

    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      _gnum->index[ipart] = malloc (sizeof(int) * _gnum->n_elts[ipart]);
      PDM_points_merge_cloud_set (id_pm, ipart, _gnum->n_elts[ipart], 
                                  _gnum->coords[ipart], _gnum->char_length[ipart]);
      for (int i = 0; i < _gnum->n_elts[ipart]; i++) {
        _gnum->index[ipart][i] = 0;
      }
    }
    
//    PDM_timer_t *timer = PDM_timer_create();                                                                                  
//    PDM_timer_resume(timer);                                                                                     
    PDM_points_merge_process (id_pm);
//    PDM_timer_hang_on(timer);                                                                                    
//    printf("Compute points merge %12.5es\n", PDM_timer_elapsed(timer));
//    PDM_timer_free(timer);

    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {

      int *candidates_idx;
      int *candidates_desc;
      
      PDM_points_merge_candidates_get (id_pm, ipart, &candidates_idx, &candidates_desc);
      
      for (int i = 0; i < _gnum->n_elts[ipart]; i++) {
        for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
          int idx = j;
          int distant_proc = candidates_desc[3*idx    ];
          int distant_part = candidates_desc[3*idx + 1];
//          int distant_pt = candidates_desc[3*idx + 2];
          
          if ((distant_proc < iproc) || ((distant_proc == iproc) && (distant_part < ipart))) {
            _gnum->index[ipart][i] = -1;
          }
        }
        if (_gnum->index[ipart][i] == 0) {
          _gnum->index[ipart][i] = n_entities;
          n_entities++;
        }
      }
    }
  }
  else {
    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      n_entities += _gnum->n_elts[ipart];
    }
  }
  
  double *coords = malloc (sizeof(double) * _gnum->dim * n_entities);
  
  if (_gnum->merge) {
    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      for (int i = 0; i <  _gnum->n_elts[ipart]; i++) {
        if (_gnum->index[ipart][i] != -1) {
          for (int k = 0; k < _gnum->dim; k++) {
            coords[_gnum->dim * _gnum->index[ipart][i] + k] =
            _gnum->coords[ipart][_gnum->dim * i + k];
          }
        }
      }
    }
  }
  else {
    int k = 0;
    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      for (int j = 0; j < _gnum->dim * _gnum->n_elts[ipart]; j++) {
        coords[k++] = _gnum->coords[ipart][j];
      }
    }
  }
  
  /* Build Morton encoding and order it */

  PDM_morton_get_coord_extents(_gnum->dim, n_entities, coords, extents, comm);

  m_code = malloc (n_entities * sizeof(PDM_morton_code_t));
  order = malloc (n_entities * sizeof(PDM_l_num_t));

  PDM_morton_encode_coords(_gnum->dim, level, extents, n_entities, coords, m_code);

  PDM_morton_local_order(n_entities, m_code, order);

  if (n_ranks > 1) {

    int rank_id;
//    PDM_l_num_t j; 
    PDM_l_num_t shift;

    size_t n_block_ents = 0;
    PDM_g_num_t current_global_num = 0, global_num_shift = 0;

    int *c_rank = NULL;
    int *send_count = NULL, *send_shift = NULL;
    int *recv_count = NULL, *recv_shift = NULL;
    double *send_coords = NULL, *recv_coords = NULL;
    PDM_l_num_t *weight = NULL;
    PDM_g_num_t *block_global_num = NULL, *part_global_num = NULL;
    PDM_morton_code_t *morton_index = NULL;

    weight = malloc (n_entities * sizeof(PDM_l_num_t));
    morton_index = malloc ((n_ranks + 1) * sizeof(PDM_morton_code_t));

    for (int i = 0; i < n_entities; i++) {
      weight[i] = 1;
    }

    PDM_morton_build_rank_index(_gnum->dim,
                                level,
                                n_entities,
                                m_code,
                                weight,
                                order,
                                morton_index,
                                comm);

    free(order);
    free(weight);
    c_rank = malloc (n_entities * sizeof(int));

    for (int i = 0; i < n_entities; i++) {
      size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                   m_code[i],
                                                   morton_index);
      c_rank[i] = (int) _c_rank; 
    }

    free(morton_index);
    free(m_code);

    /* Build send_buf, send_count and send_shift
       to build a rank to coords indexed list */

    send_count = malloc (n_ranks * sizeof (int));
    recv_count = malloc (n_ranks * sizeof (int));
    send_shift = malloc ((n_ranks + 1) * sizeof (int));
    recv_shift = malloc ((n_ranks + 1) * sizeof (int));

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (int i = 0; i < n_entities; i++)
      send_count[c_rank[i]] += _gnum->dim;

    /* Exchange number of coords to send to each process */

    PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);

    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
      recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
    }

    /* Build send and receive buffers */

    send_coords = malloc (send_shift[n_ranks] * sizeof(double));

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (int i = 0; i < n_entities; i++) {
      rank_id = c_rank[i];
      shift = send_shift[rank_id] + send_count[rank_id];
      for (int j = 0; j < _gnum->dim; j++)
        send_coords[shift + j] = coords[i*_gnum->dim + j];
      send_count[rank_id] += _gnum->dim;
    }

    recv_coords = malloc (recv_shift[n_ranks] * sizeof(double));

    /* Exchange coords between processes */

    PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                      recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                      comm);

    free(send_coords);

    /* Now re-build Morton codes on block distribution */

    n_block_ents = recv_shift[n_ranks] / _gnum->dim;

    m_code = malloc (n_block_ents * sizeof(PDM_morton_code_t));
    order = malloc (n_block_ents * sizeof(PDM_l_num_t));

    PDM_morton_encode_coords(_gnum->dim,
                             level,
                             extents,
                             n_block_ents,
                             recv_coords,
                             m_code);

    PDM_morton_local_order((int) n_block_ents, m_code, order);

    /* Check ordering; if two entities have the same Morton codes,
       use lexicographical coordinates ordering to ensure the
       final order is deterministic. */

    _check_morton_ordering(_gnum->dim, n_block_ents, recv_coords, m_code, order);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    free(m_code);
    free(recv_coords);
    block_global_num = malloc(n_block_ents * sizeof(PDM_g_num_t));

    for (int i = 0; i < n_block_ents; i++) {
      block_global_num[order[i]] = (PDM_g_num_t) i + 1;
    }
    
    free(order);

    current_global_num = (PDM_g_num_t) n_block_ents;

    /* At this stage, block_global_num[] is valid for this process, and
       current_global_num indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    PDM_MPI_Scan(&current_global_num, &global_num_shift, 1, PDM__PDM_MPI_G_NUM,
                 PDM_MPI_SUM, comm);
    global_num_shift -= current_global_num;

    for (int i = 0; i < n_block_ents; i++)
      block_global_num[i] += global_num_shift;

    /* Return global order to all processors */

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] /= _gnum->dim;
      recv_count[rank_id] /= _gnum->dim;
      send_shift[rank_id] /= _gnum->dim;
      recv_shift[rank_id] /= _gnum->dim;
    }

    send_shift[n_ranks] /= _gnum->dim;

    part_global_num = malloc (send_shift[n_ranks] * sizeof(PDM_g_num_t));

    PDM_MPI_Alltoallv(block_global_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                      part_global_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                      comm);

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] = 0;
    }

    PDM_g_num_t _max_loc = -1;

    if (_gnum->merge) {
      
      /*
       * Define local points
       */

      int k = 0;
      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        _gnum->g_nums[ipart] = malloc (sizeof(PDM_g_num_t) * _gnum->n_elts[ipart]);
        for (int j1 = 0; j1 < _gnum->n_elts[ipart]; j1++) {
          _gnum->g_nums[ipart][j1] = -1;
        }
        for (int j1 = 0; j1 < _gnum->n_elts[ipart]; j1++) {
          if (_gnum->index[ipart][j1] != -1) {
            rank_id = c_rank[k++];
            shift = send_shift[rank_id] + send_count[rank_id];
            _gnum->g_nums[ipart][j1] = part_global_num[shift];
            _max_loc = PDM_MAX (_max_loc, part_global_num[shift]);
            send_count[rank_id] += 1;
          }
        }  
      }

      int *send_count2 = malloc (n_ranks * sizeof (int));
      int *recv_count2 = malloc (n_ranks * sizeof (int));
      int *send_shift2 = malloc ((n_ranks + 1) * sizeof (int));
      int *recv_shift2 = malloc ((n_ranks + 1) * sizeof (int));

      /*
       * Count number of values to send 
       */ 
      
      for (rank_id = 0; rank_id < n_ranks; rank_id++) {
        send_count2[rank_id] = 0;
      }

      for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
        send_shift2[rank_id] = 0;
        recv_shift2[rank_id] = 0;
      }

      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        
        int *candidates_idx;
        int *candidates_desc;

        
        PDM_points_merge_candidates_get (id_pm, ipart, &candidates_idx, &candidates_desc);

        for (int i = 0; i < _gnum->n_elts[ipart]; i++) {
          int update_proc = iproc;
          
          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
//            int distant_part = candidates_desc[3*idx + 1];
//            int distant_pt   = candidates_desc[3*idx + 2];
            update_proc = PDM_MIN(update_proc, distant_proc);
          }
          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
            int distant_part = candidates_desc[3*idx + 1];
            int distant_pt   = candidates_desc[3*idx + 2];
            if ((update_proc == iproc) && (distant_proc != iproc)) {
              send_count2[distant_proc] += 1;
              assert (_gnum->g_nums[ipart][i] != -1);
            }              
            else if (  (update_proc == iproc)
                    && (distant_proc == iproc)
                    && (ipart <= distant_part)) {
              assert (_gnum->g_nums[ipart][i] != -1);
              _gnum->g_nums[distant_part][distant_pt] = _gnum->g_nums[ipart][i]; 
            }
          }
        }
      }

      PDM_MPI_Alltoall (send_count2, 1, PDM_MPI_INT, 
                        recv_count2, 1, PDM_MPI_INT, comm);

      for (rank_id = 0; rank_id < n_ranks; rank_id++) {
        send_shift2[rank_id + 1] = send_shift2[rank_id] + 3 * send_count2[rank_id];
        send_count2[rank_id] = 0;

        recv_count2[rank_id] *= 3;
        recv_shift2[rank_id + 1] = recv_shift2[rank_id] + recv_count2[rank_id];
      }  

      /*
       * Send values
       */ 

      PDM_g_num_t *send_buff = (PDM_g_num_t *) 
              malloc (sizeof(PDM_g_num_t)*send_shift2[n_ranks]);
      PDM_g_num_t *recv_buff = (PDM_g_num_t *) 
              malloc (sizeof(PDM_g_num_t)*recv_shift2[n_ranks]);

      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        
        int *candidates_idx;
        int *candidates_desc;

        PDM_points_merge_candidates_get (id_pm, ipart, &candidates_idx, &candidates_desc);
        
        for (int i = 0; i < _gnum->n_elts[ipart]; i++) {
          int update_proc = iproc;
          
          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
//            int distant_part = candidates_desc[3*idx + 1];
//            int distant_pt   = candidates_desc[3*idx + 2];
            update_proc = PDM_MIN (update_proc, distant_proc);
          }

          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
            int distant_part = candidates_desc[3*idx + 1];
            int distant_pt   = candidates_desc[3*idx + 2];

            if ((iproc == update_proc) && (distant_proc != iproc)){ 
              int idx2 = send_shift2[distant_proc] + send_count2[distant_proc];
              send_buff[idx2]     = distant_part;
              send_buff[idx2 + 1] = distant_pt;
              send_buff[idx2 + 2] = _gnum->g_nums[ipart][i];;

              send_count2[distant_proc] += 3;
            }
          }
        }
      }

      /* 
       * Send : distant_part, distant_pt, gnum
       */      
      
      PDM_MPI_Alltoallv(send_buff, send_count2, send_shift2, PDM__PDM_MPI_G_NUM,
                        recv_buff, recv_count2, recv_shift2, PDM__PDM_MPI_G_NUM,
                        comm);

      /* 
       * update gnum
       */      
        
      k = 0;
      while (k < recv_shift2[n_ranks]) {
        int ipart        = (int) recv_buff[k++];
        int ipt          = (int) recv_buff[k++];
        PDM_g_num_t gnum =       recv_buff[k++];
         _gnum->g_nums[ipart][ipt] = gnum;
      }

      free (send_buff);
      free (recv_buff);
      free (send_count2);
      free (recv_count2);
      free (send_shift2);
      free (recv_shift2);    
    
    }

    else {
      int k = 0;
      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        _gnum->g_nums[ipart] = malloc (sizeof(PDM_g_num_t) * _gnum->n_elts[ipart]);
        for (int j1 = 0; j1 < _gnum->n_elts[ipart]; j1++) {
          rank_id = c_rank[k++];
          shift = send_shift[rank_id] + send_count[rank_id];
          _gnum->g_nums[ipart][j1] = part_global_num[shift];
          _max_loc = PDM_MAX (_max_loc, part_global_num[shift]);
          send_count[rank_id] += 1;
        }  
      }
    }
    
    /* Free memory */

    free(c_rank);

    free(block_global_num);
    free(part_global_num);

    free(send_count);
    free(recv_count);
    free(send_shift);
    free(recv_shift);

    /* Get final maximum global number value */

    PDM_MPI_Allreduce (&_max_loc, 
                       &_gnum->n_g_elt, 
                       1, 
                       PDM__PDM_MPI_G_NUM, 
                       PDM_MPI_MAX, 
                       comm);

  }

  else if (n_ranks == 1) {

    _check_morton_ordering(_gnum->dim, n_entities, coords, m_code, order);

    free(m_code);
    
    PDM_g_num_t *tmp_gnum = malloc (sizeof(PDM_g_num_t) * n_entities);

    for (int i = 0; i < n_entities; i++) {
      tmp_gnum[order[i]] = (PDM_g_num_t) i+1;
    }

    if (_gnum->merge) {
      
      int *_entities = malloc (sizeof(int) * 2 * n_entities);
      int k = 0;
      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        _gnum->g_nums[ipart] = malloc (sizeof(PDM_g_num_t) * _gnum->n_elts[ipart]);
        for (int j1 = 0; j1 < _gnum->n_elts[ipart]; j1++) {
          _gnum->g_nums[ipart][j1] = -1;
          if (_gnum->index[ipart][j1] != -1) {
            _entities[k++] = ipart;
            _entities[k++] = j1;
          }
        }
      }
      
      k = 0;
      for (int i = 0; i < n_entities; i++) {
        _gnum->g_nums[_entities[2*k]][_entities[2*k+1]] = tmp_gnum[k];
        k += 1;
      }

      free (_entities);

      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        int *candidates_idx;
        int *candidates_desc;
        PDM_points_merge_candidates_get (id_pm, ipart, &candidates_idx, &candidates_desc);

        for (int i = 0; i < _gnum->n_elts[ipart]; i++) {
          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
            int distant_part = candidates_desc[3*idx + 1];
            int distant_pt   = candidates_desc[3*idx + 2];
            if ((iproc == distant_proc) && (ipart <= distant_part)) {              
              assert (_gnum->g_nums[ipart][i] != -1);
              _gnum->g_nums[distant_part][distant_pt] = _gnum->g_nums[ipart][i]; 
            }
          }
        }
        
      }
    }
    
    else {

      int k = 0;
      for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
        _gnum->g_nums[ipart] = malloc (sizeof(PDM_g_num_t) * _gnum->n_elts[ipart]);
        for (int j1 = 0; j1 <  _gnum->n_elts[ipart]; j1++) {
          _gnum->g_nums[ipart][j1] = tmp_gnum[k++];
        }  
      }
    }

    free(order);
    free(tmp_gnum);
    
    _gnum->n_g_elt = n_entities;

  }

  if (_gnum->merge) {
    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      free (_gnum->index[ipart]);
    }
    free(_gnum->index);
    PDM_points_merge_free (id_pm);
  }
  
  free (coords);

}


/**
 *
 * \brief Compute from coords
 *
 * \param [in]   _gnum          Current _pdm_gnum_t structure
 *
 */

static void
_gnum_from_parent_compute
(
 _pdm_gnum_t *_gnum 
)
{ 
  int n_procs = 0;
  PDM_MPI_Comm_size(_gnum->comm, &n_procs);
  
  int i_proc = 0;
  PDM_MPI_Comm_rank(_gnum->comm,
                &i_proc);

  int *sendBuffN   = (int *) malloc(sizeof(int) * n_procs);
  int *sendBuffIdx = (int *) malloc(sizeof(int) * n_procs);

  int *recvBuffN   = (int *) malloc(sizeof(int) * n_procs);
  int *recvBuffIdx = (int *) malloc(sizeof(int) * n_procs);

  /* Calcul du nombre total d'elements du bloc */

  PDM_l_num_t n_elt_loc_total = 0;
  PDM_g_num_t l_max_parent = -1;
  
  for (int j = 0; j < _gnum->n_part; j++) {
    n_elt_loc_total += _gnum->n_elts[j];
    for (int k = 0; k < _gnum->n_elts[j]; k++) {
      l_max_parent = PDM_MAX (l_max_parent, _gnum->parent[j][k]);
    }
  }

  PDM_g_num_t max_parent = 0;
  PDM_MPI_Allreduce (&l_max_parent, &max_parent, 1, 
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _gnum->comm);
  
  /* Comptage du nombre d'elements a envoyer a chaque processus */
  
  for (int j = 0; j < n_procs; j++) {
    sendBuffN[j] = 0;
    sendBuffIdx[j] = 0;
    recvBuffN[j] = 0;
    recvBuffIdx[j] = 0;
  }

  PDM_g_num_t *d_elt_proc = 
          (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_procs + 1));
  

  PDM_g_num_t div_entiere = max_parent / n_procs;
  PDM_g_num_t div_reste = max_parent % n_procs;

  d_elt_proc[0] = 1;
  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] =  div_entiere;
    if (i < div_reste) {
      d_elt_proc[i+1] += 1;
    }
  }

  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] += d_elt_proc[i]; 
  }
          
  for (int j = 0; j < _gnum->n_part; j++) {
    for (int k = 0; k < _gnum->n_elts[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long (_gnum->parent[j][k],
                                                         d_elt_proc,
                                                         n_procs + 1);
      sendBuffN[i_elt_proc] += 1;
    }
  }

  sendBuffIdx[0] = 0;
  for (int j = 1; j < n_procs; j++) {
    sendBuffIdx[j] = sendBuffIdx[j-1] + sendBuffN[j-1];
  }
   
  /* Determination du nombre d'elements recu de chaque processus */

  PDM_MPI_Alltoall(sendBuffN, 
               1, 
               PDM_MPI_INT, 
               recvBuffN, 
               1, 
               PDM_MPI_INT, 
               _gnum->comm);

  recvBuffIdx[0] = 0;
  for(int j = 1; j < n_procs; j++) {
    recvBuffIdx[j] = recvBuffIdx[j-1] + recvBuffN[j-1];
  }

  /* Transmission des numeros absolus  */
 
  PDM_g_num_t _l_numabs_tmp = d_elt_proc[i_proc+1] - d_elt_proc[i_proc];
  int l_numabs_tmp = (int) _l_numabs_tmp;
  PDM_g_num_t *numabs_tmp = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * l_numabs_tmp);
  
  PDM_g_num_t *n_elt_stocke_procs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_procs+1));
 
  PDM_g_num_t *sendBuffNumabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                                       n_elt_loc_total);
  PDM_g_num_t *recvBuffNumabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                                      (recvBuffIdx[n_procs - 1] + 
                                                       recvBuffN[n_procs - 1]));
                                                          
  for (int j = 0; j < n_procs; j++) {
    sendBuffN[j] = 0;
  }

  for (int j = 0; j < _gnum->n_part; j++) {
    for (int k = 0; k < _gnum->n_elts[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(_gnum->parent[j][k],
                                                        d_elt_proc,
                                                        n_procs+1);
      sendBuffNumabs[sendBuffIdx[i_elt_proc] + sendBuffN[i_elt_proc]] = _gnum->parent[j][k];
      sendBuffN[i_elt_proc] += 1;
    }
  }

  PDM_MPI_Alltoallv((void *) sendBuffNumabs, 
                sendBuffN, 
                sendBuffIdx, 
                PDM__PDM_MPI_G_NUM,
                (void *) recvBuffNumabs, 
                recvBuffN, 
                recvBuffIdx,
                PDM__PDM_MPI_G_NUM, 
                _gnum->comm);
  
  /* Echange du nombre d'elements stockes sur chaque processus */

  const PDM_g_num_t n_elt_stocke = 
    (PDM_g_num_t) (recvBuffIdx[n_procs - 1] + recvBuffN[n_procs - 1]);

  PDM_MPI_Allgather((void *) &n_elt_stocke, 
                1,
                PDM__PDM_MPI_G_NUM,
                (void *) (n_elt_stocke_procs + 1), 
                1,
                PDM__PDM_MPI_G_NUM,
                _gnum->comm);

  n_elt_stocke_procs[0] = 1;
  for (int j = 1; j < n_procs + 1; j++) {
    n_elt_stocke_procs[j] += n_elt_stocke_procs[j-1];
  }    

  /* Stockage du resultat et determination de la nouvelle numerotation absolue
     independante du parallelisme */
    
  for (int j = 0; j < l_numabs_tmp; j++) 
    numabs_tmp[j] = 0;

  for (int j = 0; j < n_procs; j++) {
    
    const int ideb = recvBuffIdx[j];
    const int ifin = recvBuffIdx[j] + recvBuffN[j];
    
    for (int k = ideb; k < ifin; k++) {
      
      PDM_g_num_t _idx = recvBuffNumabs[k] - d_elt_proc[i_proc];
      const int idx = (int) _idx;
      assert((idx < l_numabs_tmp) && (idx >= 0));

      numabs_tmp[idx] = 1; /* On marque les elements */
    }
  }

  int cpt_elt_proc = 0;
  for (int j = 0; j < l_numabs_tmp; j++) {
    if (numabs_tmp[j] == 1) {

      /* On fournit une numerotation independante du parallelisme */
      
      numabs_tmp[j] = n_elt_stocke_procs[i_proc] + cpt_elt_proc;
      cpt_elt_proc += 1;
    }
  }

  /* On remplit le buffer de reception qui devient le buffer d'envoi
     Le buffer d'envoi devient lui le buffer de reception */

  cpt_elt_proc = 0;
  for (int j = 0; j < n_procs; j++) {
    
    const int ideb = recvBuffIdx[j];
    const int ifin = recvBuffIdx[j] + recvBuffN[j];
    
    for (int k = ideb; k < ifin; k++) {
      
      PDM_g_num_t _idx = recvBuffNumabs[k] - d_elt_proc[i_proc];
      const int idx = (int) _idx;
      
      recvBuffNumabs[cpt_elt_proc] = numabs_tmp[idx];
      
      cpt_elt_proc += 1;
    }
  }

  PDM_MPI_Alltoallv((void *) recvBuffNumabs, 
                recvBuffN, 
                recvBuffIdx, 
                PDM__PDM_MPI_G_NUM,
                (void *) sendBuffNumabs, 
                sendBuffN, 
                sendBuffIdx,
                PDM__PDM_MPI_G_NUM, 
                _gnum->comm);

  /* On Stocke l'information recue */

  for (int j = 0; j < n_procs; j++)
    sendBuffN[j] = 0;

  for (int j = 0; j < _gnum->n_part; j++) {
    
    _gnum->g_nums[j] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _gnum->n_elts[j]);

    for (int k = 0; k < _gnum->n_elts[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(_gnum->parent[j][k],
                                                        d_elt_proc,
                                                        n_procs+1);
      _gnum->g_nums[j][k] = sendBuffNumabs[sendBuffIdx[i_elt_proc] + sendBuffN[i_elt_proc]];
      sendBuffN[i_elt_proc] += 1;
    }
  }

  /* Liberation memoire */

  free(sendBuffIdx);
  free(sendBuffN);  
  free(recvBuffIdx);
  free(recvBuffN);  
  free(sendBuffNumabs);
  free(recvBuffNumabs);
  free(d_elt_proc);
  free(numabs_tmp);
  free(n_elt_stocke_procs);
  
}
  
/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   dim          Spatial dimension 
 * \param [in]   n_part       Number of local partitions 
 * \param [in]   merge        Merge double points or not
 * \param [in]   tolerance    Geometric tolerance (used if merge double points is activated)
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
)
{
  
  /*
   * Search a ppart free id
   */

  if (_gnums == NULL) {
    _gnums = PDM_Handles_create (4);
  }

  _pdm_gnum_t *_gnum = (_pdm_gnum_t *) malloc(sizeof(_pdm_gnum_t));
  int id = PDM_Handles_store (_gnums, _gnum);

  _gnum->n_part      = n_part;
  _gnum->dim         = dim;
  _gnum->merge       = merge;
  _gnum->tolerance   = tolerance;
  _gnum->comm        = comm;
  _gnum->n_g_elt     = -1;
  _gnum->g_nums      = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * n_part); 
  _gnum->coords      = NULL;
  _gnum->char_length = NULL;
  _gnum->parent      = NULL;
  _gnum->n_elts      = (int *) malloc (sizeof(int) * n_part);  
  _gnum->index       = NULL;
  
  for (int i = 0; i < n_part; i++) {
    _gnum->g_nums[i] = NULL;
  }

  return id;

}

void
PROCF (pdm_gnum_create, PDM_GNUM_CREATE)
(
 const int          *dim,
 const int          *n_part,
 const int          *merge,
 const double       *tolerance,       
 const PDM_MPI_Fint *fcomm,
       int          *id
)
{
  const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

  *id = PDM_gnum_create (*dim, *n_part, (PDM_bool_t) *merge, *tolerance, c_comm);
}


/**
 *
 * \brief Set from coordinates
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be determined by lexicographical ordering of coordinates.
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
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);
  
  if (_gnum->coords == NULL) {
    _gnum->coords = (double **) malloc (sizeof(double * ) * _gnum->n_part);
    for (int i = 0; i < _gnum->n_part; i++) {
      _gnum->coords[i_part] = NULL;
    }
  }

  if (_gnum->merge && _gnum->char_length == NULL) {
    _gnum->char_length = (double **) malloc (sizeof(double * ) * _gnum->n_part);
    for (int i = 0; i < _gnum->n_part; i++) {
      _gnum->char_length[i_part] = NULL;
    }
  }
  _gnum->coords[i_part]      = (double *) coords;
  if (_gnum->merge) {
    _gnum->char_length[i_part] = (double *) char_length;  
  }
  _gnum->n_elts[i_part]      = n_elts;
  _gnum->g_nums[i_part]      = NULL;
 
}     

void
PROCF (pdm_gnum_set_from_coords, PDM_GNUM_SET_FROM_COORDS)
(
 const int    *id,
 const int    *i_part,
 const int    *n_elts,
 const double *coords,
 const double *char_length
)
{
  PDM_gnum_set_from_coords (*id, *i_part, *n_elts, coords, char_length);
}


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
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);
  
  if (_gnum->parent == NULL) {
    _gnum->parent = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * _gnum->n_part);
    for (int i = 0; i < _gnum->n_part; i++) {
      _gnum->parent[i_part] = NULL;
    }
  }

  _gnum->parent[i_part] = (PDM_g_num_t *) parent_gnum;
  _gnum->n_elts[i_part] = n_elts;
  _gnum->g_nums[i_part] = NULL;
  
}

void
PROCF (pdm_gnum_set_from_parents, PDM_GNUM_SET_FROM_PARENTS)
(
 const int *id,
 const int *i_part,
 const int *n_elts,
 const PDM_g_num_t *parent_gnum
)
{
  PDM_gnum_set_from_parents (*id, *i_part, *n_elts, parent_gnum);
}


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
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);
  
  if (_gnum->coords != NULL) {
    _gnum_from_coords_compute (_gnum);
  }
  
  else if (_gnum->parent != NULL) {
    _gnum_from_parent_compute (_gnum);    
  }

}

void
PROCF (pdm_gnum_compute, PDM_GNUM_COMPUTE)
(
 const int *id
)
{
  PDM_gnum_compute (*id);
}


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
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);

  return _gnum->g_nums[i_part];
}

void
PROCF (pdm_gnum_get, PDM_GNUM_GET)
(
 const int *id,
 const int *i_part,
 PDM_g_num_t *gnum
)
{
  _pdm_gnum_t *_gnum = _get_from_id (*id);
  
  const PDM_g_num_t *tmp = PDM_gnum_get (*id, *i_part);
  for (int i = 0; i < _gnum->n_elts[*i_part]; i++) {
    gnum[i] = tmp[i];    
  }
}


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
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);

  if (_gnum->coords != NULL) {
    free (_gnum->coords);
  }
  
  if (_gnum->char_length != NULL) {
    free (_gnum->char_length);
  }

  if (_gnum->parent != NULL) {
    free (_gnum->parent);
  }

  if (partial != 1) {
    for (int i = 0; i < _gnum->n_part; i++) {
      free (_gnum->g_nums[i]);
    }
  }
  
  free (_gnum->g_nums);
  free (_gnum->n_elts);
  
  free (_gnum);
  
  PDM_Handles_handle_free (_gnums, id, PDM_FALSE);

  const int n_gnum = PDM_Handles_n_get (_gnums);
  
  if (n_gnum == 0) {
    _gnums = PDM_Handles_free (_gnums);
  }

}

void
PROCF (pdm_gnum_free, PDM_GNUM_FREE)
(
 const int *id,
 const int *partial       
)
{
  PDM_gnum_free (*id, *partial);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
