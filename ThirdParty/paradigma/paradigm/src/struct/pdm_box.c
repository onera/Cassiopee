/*============================================================================
 * Handle boxes aligned with Cartesian axes.
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_box.h"
#include "pdm_box_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type and structure definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Display a histogram on leaves associated to the boxes and several
 * other pieces of information (min, max, ...)
 *
 * parameters:
 *   distrib          <-- pointer to the PDM_box_distrib_t structure
 *   n_quantiles      <-> number of quantiles requested (or NULL);
 *                        may be lower upon return
 *   quantile_start   --> start of quantile (size: n_quantiles + 1)
 *   n_quantile_boxes --> number of boxes per quantile (size: n_quantiles)
 *   imbalance        --> distribution imbalance measure (or NULL)
 *   n_ranks          --> number of ranks with boxes (or NULL)
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static void
_get_distrib_statistics(const PDM_box_distrib_t  *distrib,
                        int                *n_quantiles,
                        int                 quantile_start[],
                        int                 n_quantile_boxes[],
                        double                   *imbalance,
                        int                      *n_ranks,
                        PDM_MPI_Comm                  comm)
{
  int   i, j, k, step, delta, _n_rank_boxes;

  int  _n_ranks = 0;
  int   _min = INT_MAX, _max = 0, gmin = 0, gmax = 0;

  /* Sanity checks */

  assert(distrib != NULL);
  assert(distrib->index != NULL);

  if (n_quantiles != NULL) {

    int _n_quantiles = 1;

    /* Get global min and max number of boxes */

    for (i = 0; i < distrib->n_ranks; i++) {

      _n_rank_boxes = distrib->index[i+1] - distrib->index[i];
      _min = PDM_MIN(_min, _n_rank_boxes);
      _max = PDM_MAX(_max, _n_rank_boxes);

      if (_n_rank_boxes > 0)
        _n_ranks += 1;

    }

    gmin = _min;
    gmax = _max;

    PDM_MPI_Allreduce(&_min, &gmin, 1, PDM_MPI_INT, PDM_MPI_MIN, comm);
    PDM_MPI_Allreduce(&_max, &gmax, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

    /* Build a histogram for the distribution of boxes */

    delta = gmax - gmin;
    if (delta < _n_quantiles)
      _n_quantiles = delta;

    /* Define quantiles */

    step = delta / _n_quantiles;
    if (delta % _n_quantiles > 0)
      step++;

    for (i = 0; i < _n_quantiles; i++)
      quantile_start[i] = gmin + i*step;

    quantile_start[_n_quantiles] = gmax + 1;

    /* Count for quantiles */

    for (j = 0; j < _n_quantiles; j++)
      n_quantile_boxes[j] = 0;

    if (delta > 0) {  /* Loop on boxes */

      for (i = 0; i < distrib->n_ranks; i++) {

        _n_rank_boxes = distrib->index[i+1] - distrib->index[i];

        for (k = 1; k < _n_quantiles; k++) {
          if (_n_rank_boxes < gmin + k*step)
            break;
        }
        n_quantile_boxes[k-1] += 1;

      } /* End of loop on boxes */

    }

    *n_quantiles = _n_quantiles;
  }

  /* Set other return values */

  if (imbalance != NULL)
    *imbalance = distrib->fit;

  if (n_ranks != NULL)
    *n_ranks = _n_ranks;
}

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create a set of boxes and initialize it.
 *
 * parameters:
 *   dim              <-- spatial dimension
 *   normalize        <-- 1 if boxes are to be normalized, 0 otherwize
 *   allow_projection <-- if 1, project to lower dimension if all boxes
 *                        are cut by the median plane of the set.
 *   n_boxes          <-- number of elements to create
 *   box_gnum         <-- global numbering of boxes
 *   extents          <-- coordinate extents (size: n_boxes*dim*2, as
 *                        xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 *   origin   <--  initial location (size: n_boxes*3, as
 *                        iproc, ipart, local num, ...)
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   a new allocated pointer to a PDM_box_set_t structure.
 *---------------------------------------------------------------------------*/

PDM_box_set_t *
PDM_box_set_create(int                dim,
                   int                normalize,
                   int                allow_projection,
                   int                n_boxes,
                   const PDM_g_num_t  *box_gnum,
                   const double      *box_extents,
                   const int          n_part_orig,
                   const int         *n_boxes_orig,
                   const int         *origin,
                   PDM_MPI_Comm           comm)
{
  int j, k;
  int   i;
  PDM_g_num_t  n_g_boxes = n_boxes;
  double  g_min[3], g_max[3], g_extents[6];

  PDM_box_set_t  *boxes = NULL;

  /* Get global min/max coordinates */

  PDM_morton_get_global_extents(dim, n_boxes, box_extents, g_extents, comm);

  for (j = 0; j < 3; j++) {
    g_min[j] = g_extents[j];
    g_max[j] = g_extents[j+dim];
  }

  if (comm != PDM_MPI_COMM_NULL) {

    PDM_g_num_t  box_max = 0;

    for (i = 0; i < n_boxes; i++)
      box_max = PDM_MAX(box_max, box_gnum[i]);

    PDM_MPI_Allreduce(&box_max, &n_g_boxes, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);

  }

  /* Allocate box set structure and initialize it */

  boxes = (PDM_box_set_t  *) malloc (sizeof(PDM_box_set_t));

  boxes->dim = dim;
  boxes->n_boxes      = n_boxes;
  boxes->n_g_boxes = n_g_boxes;

  for (j = 0; j < 3; j++) {
    boxes->dimensions[j] = j;
    boxes->gmin[j] = g_min[j];
    boxes->gmax[j] = g_max[j];
  }

  boxes->g_num = NULL;
  boxes->extents = NULL;

  boxes->n_part_orig = n_part_orig;
  boxes->n_boxes_orig = NULL;
  boxes->origin = NULL;

  boxes->comm = comm;

  /* Optionally allow and detect a layout of lower
     dimension than the spatial dimension */

  if (allow_projection) {

    double g_mid[3];
    int proj[] = {1, 1, 1};

    for (j = 0; j < dim; j++)
      g_mid[j] = (g_min[j] + g_max[j]) * 0.5;

    for (i = 0; i < n_boxes; i++) {
      for (j = 0; j < dim; j++) {
        if (   box_extents[i*dim*2 + j]     > g_mid[j]
            || box_extents[i*dim*2 + j+dim] < g_mid[j])
          proj[j] = 0;
      }
    }

    if (comm != PDM_MPI_COMM_NULL) {
      int l_proj[3];
      for (j = 0; j < dim; j++)
        l_proj[j] = proj[j];
      PDM_MPI_Allreduce(l_proj, proj, dim, PDM_MPI_INT, PDM_MPI_MIN, comm);
    }

    boxes->dim = 0;
    for (j = 0; j < dim; j++) {
      if (proj[j] == 0) {
        boxes->dimensions[boxes->dim] = j;
        boxes->dim += 1;
      }
    }

  }

  for (j = boxes->dim; j < 3; j++) /* ensure all is initialized */
    boxes->dimensions[j] = -1;

  /* Now assign values */

  boxes->g_num = (PDM_g_num_t *) malloc (n_boxes * sizeof(PDM_g_num_t));
  boxes->extents = (double *) malloc(n_boxes*boxes->dim*2 * sizeof(double));
  boxes->origin = (int *) malloc(n_boxes * 3 * sizeof(int));
  boxes->n_boxes_orig = (int *) malloc(n_part_orig * sizeof(int));

  memcpy(boxes->n_boxes_orig, n_boxes_orig, n_part_orig * sizeof(int));

  for (i = 0; i < n_boxes; i++) {

    double *_min = boxes->extents + (boxes->dim*2*i);
    double *_max = _min + boxes->dim;

    boxes->g_num[i] = box_gnum[i];
    for (j = 0; j < 3; j++) {
      boxes->origin[3*i+j] = origin[3*i+j];
    }

    for (j = 0; j < boxes->dim; j++) {
      k = boxes->dimensions[j];
      _min[j] = box_extents[i*dim*2 + k];
      _max[j] = box_extents[i*dim*2 + k+dim];
      assert(_min[j] <= _max[j]);
    }
  }

  /* Define the normalized min/max coordinates of the box */

  boxes->normalized = normalize;
  if (normalize) {

    for (j = 0; j < boxes->dim; j++) {
      k = boxes->dimensions[j];
      boxes->s[j] = g_min[k];
      boxes->d[j] = g_max[k] - g_min[k];
    }

    for (i = 0; i < n_boxes; i++) {

      double *_min = boxes->extents + (boxes->dim*2*i);
      double *_max = _min + boxes->dim;

      PDM_box_set_normalize (boxes, _min, _min);
      PDM_box_set_normalize (boxes, _max, _max);

    }

  }

  /* Return pointer to structure */

  return boxes;
}


void
PDM_box_set_normalize
(
 PDM_box_set_t  *boxes,
 const double *pt_origin,
 double *pt_nomalized
 )
{
  for (int j = 0; j < boxes->dim; j++) {
    pt_nomalized[j] = (pt_origin[j] - boxes->s[j]) / boxes->d[j];
  }
}



void
PDM_box_set_normalize_inv
(
 PDM_box_set_t  *boxes,
 const double *pt_nomalized,
 double *pt_origin
 )
{
  for (int j = 0; j < boxes->dim; j++) {
    pt_origin[j] = pt_nomalized[j] * boxes->d[j] + boxes->s[j];
  }
}


/*----------------------------------------------------------------------------
 * Delete a PDM_box_set_t structure.
 *
 * parameters:
 *   boxes <-> pointer to the PDM_box_set_t structure to delete
 *---------------------------------------------------------------------------*/

void
PDM_box_set_destroy(PDM_box_set_t  **boxes)
{
  if (boxes != NULL) {

    PDM_box_set_t  *_boxes = *boxes;

    if (_boxes == NULL)
      return;

    free(_boxes->g_num);
    free(_boxes->extents);
    free(_boxes->origin);
    free(_boxes->n_boxes_orig);
    free(_boxes);
  }
}

/*----------------------------------------------------------------------------
 * Return the dimension associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   associated spatial dimension
 *---------------------------------------------------------------------------*/

int
PDM_box_set_get_dim(const PDM_box_set_t  *boxes)
{
  int retval = 0;

  if (boxes != NULL)
    retval = boxes->dim;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the local number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   local number of boxes
 *---------------------------------------------------------------------------*/

int
PDM_box_set_get_size(const PDM_box_set_t  *boxes)
{
  int retval = 0;

  if (boxes != NULL)
    retval = boxes->n_boxes;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the global number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   global number of boxes
 *---------------------------------------------------------------------------*/

PDM_g_num_t
PDM_box_set_get_global_size(const PDM_box_set_t  *boxes)
{
  PDM_g_num_t retval = 0;

  if (boxes != NULL)
    retval = boxes->n_g_boxes;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return extents associated with a set of boxes.
 *
 * The extents array is organized in the following fashion:
 * {x_min_0, y_min_0, ..., x_max_0, y_max_0, ...
 *  x_min_n, y_min_n, ..., x_max_n, y_max_n, ...}
 *
 * Its size is thus: n_boxes * dim * 2.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to extents array
 *---------------------------------------------------------------------------*/

const double *
PDM_box_set_get_extents(PDM_box_set_t  *boxes)
{
  assert(boxes != NULL);

  return boxes->extents;
}

/*----------------------------------------------------------------------------
 * Return global numbers associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to global box numbers array
 *---------------------------------------------------------------------------*/

const PDM_g_num_t *
PDM_box_set_get_g_num(PDM_box_set_t  *boxes)
{
  assert(boxes != NULL);

  return boxes->g_num;
}

/*----------------------------------------------------------------------------
 * Return initial location associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to initial location array
 *---------------------------------------------------------------------------*/

const int *
PDM_box_set_origin_get (PDM_box_set_t  *boxes)
{
  assert(boxes != NULL);

  return boxes->origin;
}

/*----------------------------------------------------------------------------
 * Build a Morton_index to get a well-balanced distribution of the boxes.
 *
 * parameters:
 *  boxes      <-- pointer to associated PDM_box_set_t structure
 *  distrib    <-> pointer to a PDM_box_distrib_t structure
 *  n_leaves   <-- number of leaves with weight > 0
 *  leaf_codes <-- Morton code for each leaf
 *  weight     <-- number of boxes related to each leaf
 *---------------------------------------------------------------------------*/

void
PDM_box_set_build_morton_index(const PDM_box_set_t  *boxes,
                               PDM_box_distrib_t    *distrib,
                               int             n_leaves,
                               PDM_morton_code_t    *leaf_codes,
                               int            *weight)
{

  int   *order = NULL;

  assert(distrib != NULL);
  assert(distrib->morton_index != NULL);

  order = (int *) malloc(n_leaves * sizeof(int));

  /* Locally order Morton encoding */

  PDM_morton_local_order(n_leaves,
                         leaf_codes,
                         order);

  /* Compute a Morton index on ranks and return the associated fit */

  if (boxes->comm != PDM_MPI_COMM_NULL)
    distrib->fit = PDM_morton_build_rank_index(boxes->dim,
                                               distrib->max_level,
                                               n_leaves,
                                               leaf_codes,
                                               weight,
                                               order,
                                               distrib->morton_index,
                                               boxes->comm);
  /* Free memory */

  free(order);

}

/*----------------------------------------------------------------------------
 * Redistribute boxes over the ranks according to the Morton index to
 * assume a better balanced distribution of the boxes.
 *
 * parameters:
 *  distrib <--  data structure on box distribution
 *  boxes   <->  pointer to the structure to redistribute
 *---------------------------------------------------------------------------*/

void
PDM_box_set_redistribute(const PDM_box_distrib_t  *distrib,
                         PDM_box_set_t            *boxes)
{

  int  rank_id;

  int   i, j;
  int  *send_count = NULL, *send_shift = NULL;
  int  *recv_count = NULL, *recv_shift = NULL;
  PDM_g_num_t  *send_g_num = NULL;
  double  *send_extents = NULL;
  int  *send_origin = NULL;

  const int stride = boxes->dim * 2;
  const int stride_origin= 3;

  /* Sanity checks */

  assert(distrib != NULL);
  assert(boxes != NULL);
  assert(distrib->n_ranks > 1);

  /* Build send_buf, send_count and send_shift
     to build a rank to boxes indexed list */

  send_count = (int *) malloc(distrib->n_ranks * sizeof(int));
  recv_count = (int *) malloc(distrib->n_ranks * sizeof(int));
  send_shift = (int *) malloc((distrib->n_ranks + 1) * sizeof(int));
  recv_shift = (int *) malloc((distrib->n_ranks + 1) * sizeof(int));

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++)
    send_count[rank_id]
      = distrib->index[rank_id+1] - distrib->index[rank_id];

  /* Exchange number of boxes to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
               recv_count, 1, PDM_MPI_INT, boxes->comm);

  for (i = 0; i < distrib->n_ranks; i++)
    send_shift[i] = distrib->index[i];

  recv_shift[0] = 0;
  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++)
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];

  /* Build send_buffers */

  send_g_num = (PDM_g_num_t *) malloc(distrib->index[distrib->n_ranks] * sizeof(PDM_g_num_t));
  send_extents = (double *) malloc(distrib->index[distrib->n_ranks] * boxes->dim * 2
                                   * sizeof(double));
  send_origin = (int *) malloc(distrib->index[distrib->n_ranks] * 3
                                   * sizeof(int));

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++)
    send_count[rank_id] = 0;

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++) {

    for (i = distrib->index[rank_id];
         i < distrib->index[rank_id+1];
         i++) {

      int   box_id = distrib->list[i];
      int   shift = distrib->index[rank_id] + send_count[rank_id];

      send_g_num[shift] = boxes->g_num[box_id];

      for (j = 0; j < stride; j++)
        send_extents[shift*stride + j] = boxes->extents[box_id*stride + j];

      for (j = 0; j < stride_origin; j++)
        send_origin[shift*stride_origin + j] =
          boxes->origin[box_id*stride_origin + j];

      send_count[rank_id] += 1;

    }

  } /* End of loop on ranks */

  /* Prepare to replace the local arrays */

  boxes->n_boxes = recv_shift[distrib->n_ranks];
  free(boxes->g_num);
  free(boxes->extents);
  free(boxes->origin);

  boxes->g_num = (PDM_g_num_t *) malloc(boxes->n_boxes * sizeof(PDM_g_num_t));
  boxes->extents = (double *) malloc(boxes->n_boxes*stride * sizeof(double));
  boxes->origin = (int *) malloc(boxes->n_boxes*stride_origin * sizeof(int));

  /* Exchange boxes between processes */

  PDM_MPI_Alltoallv(send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                boxes->g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                boxes->comm);

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++) {
    send_count[rank_id] *= stride;
    send_shift[rank_id] *= stride;
    recv_count[rank_id] *= stride;
    recv_shift[rank_id] *= stride;
  }

  PDM_MPI_Alltoallv(send_extents, send_count, send_shift, PDM_MPI_DOUBLE,
                boxes->extents, recv_count, recv_shift, PDM_MPI_DOUBLE,
                boxes->comm);

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++) {
    send_count[rank_id] = send_count[rank_id]/stride * stride_origin;
    send_shift[rank_id] = send_shift[rank_id]/stride * stride_origin;
    recv_count[rank_id] = recv_count[rank_id]/stride * stride_origin;
    recv_shift[rank_id] = recv_shift[rank_id]/stride * stride_origin;
  }

  PDM_MPI_Alltoallv(send_origin, send_count, send_shift, PDM_MPI_INT,
                boxes->origin, recv_count, recv_shift, PDM_MPI_INT,
                boxes->comm);

  /* Free buffers */

  free(send_g_num);
  free(send_extents);
  free(send_origin);
  free(send_count);
  free(send_shift);
  free(recv_count);
  free(recv_shift);

}

/*----------------------------------------------------------------------------
 * Dump a PDM_box_set_t structure.
 *
 * parameters:
 *   boxes     <-- pointer to the PDM_box_t structure
 *   verbosity <-- verbosity level (0 or 1)
 *----------------------------------------------------------------------------*/

void
PDM_box_set_dump(const PDM_box_set_t  *boxes,
                 int                   verbosity)
{
  int   i;

  const char  XYZ[3] = "XYZ";

  if (boxes == NULL)
    return;

  /* Print basic information */

  if (boxes->dim == 3)
    PDM_printf("\nBox set (3D layout):\n\n"
               "global min/max on selected faces:\n"
               "  [%7.5e %7.5e %7.5e] --> [%7.5e %7.5e %7.5e]\n",
               boxes->gmin[0], boxes->gmin[1], boxes->gmin[2],
               boxes->gmax[0], boxes->gmax[1], boxes->gmax[2]);

  else if (boxes->dim == 2) {
    PDM_printf("\nBox set (2D layout, selected axes [%c, %c]\n\n",
               XYZ[boxes->dimensions[0]],
               XYZ[boxes->dimensions[1]]);
    PDM_printf("global min/max on selected faces:\n"
               "  [%7.5e %7.5e] --> [%7.5e %7.5e]\n",
               boxes->gmin[boxes->dimensions[0]],
               boxes->gmin[boxes->dimensions[1]],
               boxes->gmax[boxes->dimensions[0]],
               boxes->gmax[boxes->dimensions[1]]);
  }

  else if (boxes->dim == 1) {
    PDM_printf("\nBox set (1D layout, selected axis [%c]\n\n",
               XYZ[boxes->dimensions[0]]);
    PDM_printf("global min/max on selected faces:\n"
               "  [%7.5e %7.5e] --> [%7.5e %7.5e]\n",
               boxes->gmin[boxes->dimensions[0]],
               boxes->gmin[boxes->dimensions[1]],
               boxes->gmax[boxes->dimensions[0]],
               boxes->gmax[boxes->dimensions[1]]);
  }
  fflush(stdout);

  /* Print detailed box information */

  if (verbosity < 1)
    return;

  if (boxes->dim == 3) {
    for (i = 0; i < boxes->n_boxes; i++) {
      const double *bmin = boxes->extents + i*6;
      const double *bmax = boxes->extents + i*6 + 3;
      PDM_printf("  id %8d, num %9llu: "
                 "[%7.5e %7.5e %7.5e] --> [%7.5e %7.5e %7.5e]\n",
                 i, (unsigned long long)(boxes->g_num[i]),
                 bmin[0], bmin[1], bmin[2],
                 bmax[0], bmax[1], bmax[2]);
    }
  }

  else if (boxes->dim == 2) {
    for (i = 0; i < boxes->n_boxes; i++) {
      const double *bmin = boxes->extents + i*4;
      const double *bmax = boxes->extents + i*4 + 2;
      PDM_printf("  id %8d, num %9llu: "
                 "[%7.5e %7.5e] --> [%7.5e %7.5e]\n",
                 i, (unsigned long long)(boxes->g_num[i]),
                 bmin[0], bmin[1], bmax[0], bmax[1]);
    }
  }

  else if (boxes->dim == 1) {
    for (i = 0; i < boxes->n_boxes; i++) {
      const double *bmin = boxes->extents + i*2;
      const double *bmax = boxes->extents + i*2 + 1;
      PDM_printf("  id %8d, num %9llu: "
                 "[%7.5e] --> [%7.5e]\n",
                 i, (unsigned long long)(boxes->g_num[i]),
                 bmin[0], bmax[0]);
    }
  }

  /* Sanity check */

  for (i = 0; i < boxes->n_boxes; i++) {
    int j;
    const double *bmin = boxes->extents + boxes->dim*2*i;
    const double *bmax = boxes->extents + boxes->dim*(2*i + 1);
    for (j = 0; j < boxes->dim; j++) {
      if (bmin[j] > bmax[j]) {
        PDM_error(__FILE__, __LINE__, 0,
                  "PDM_box_set_dump error : Inconsistent box found (min > max):\n"
                    "  global number:  %llu\n"
                    "  min       :  %10.4g\n"
                    "  max       :  %10.4g\n",
                  (unsigned long long)(boxes->g_num[i]), bmin[j], bmax[j]);
        abort();
      }
    }
  }

}


/*----------------------------------------------------------------------------
 * Receive data from origin for any box
 *
 * parameters:
 *   box_set                <-- pointer to the PDM_box_t structure
 *   t_stride               <-- Type of stride
 *   stride_cst             <-- Constant stride
 *   data_size              <-- Size of data
 *   origin_distrib_stride  <-- Origin stride distribution
 *   origin_distrib_data    <-- Origin data distribution
 *   current_distrib_stride <-> Current stride distribution (Allocate and compute if input is NULL,
 *                              otherwise nothing)
 *   current_distrib_data   --> Current data distribution
 *
 * return : size of current_distrib_data
 *----------------------------------------------------------------------------*/

void
PDM_box_set_recv_data_from_origin_distrib
(
 PDM_box_set_t  *boxes,
 PDM_stride_t   t_stride,
 int            stride_cst,
 size_t         data_size,
 int          **origin_distrib_stride,
 void         **origin_distrib_data,
 int          **current_distrib_stride,
 void         **current_distrib_data
)
{

  /*
   * If current_distrib_stride != NULL :
   * No allocation for current_distrib_stride
   */

  int alloc_distrib_stride = 1;
  if (*current_distrib_stride != NULL)
    alloc_distrib_stride = 0;

  *current_distrib_data = NULL;
  unsigned char **_origin_distrib_data = (unsigned char **) origin_distrib_data;

  /* Send origin properties to the origin process :
   *   - Compute the number element to send for any process
   *   - Exchange -> origin these numbers all_to_all
   *   - Exchange -> origin origin properties (part, num loc) all_to_all_v */

  int n_boxes = boxes->n_boxes;
  int *origin = boxes->origin;

  int s_comm = PDM_MPI_Comm_rank (boxes->comm, &s_comm);

  int *curr_count = (int *) malloc (sizeof(int) * s_comm);
  for (int i = 0; i < s_comm; i++) {
    curr_count[i] = 0;
  }

  for (int i = 0; i < n_boxes; i++) {
    curr_count[origin[3*i]] += 1;
  }

  int *orig_count = (int *) malloc (sizeof(int) * s_comm);
  PDM_MPI_Alltoall(curr_count, 1, PDM_MPI_INT,
               orig_count, 1, PDM_MPI_INT, boxes->comm);

  int *curr_shift = (int *) malloc (sizeof(int) * (s_comm + 1));
  int *orig_shift = (int *) malloc (sizeof(int) * (s_comm + 1));

  for (int i = 0; i < s_comm + 1; i++) {
    curr_shift[i] = 0;
    orig_shift[i] = 0;
  }

  for (int i = 0; i < s_comm; i++) {
    curr_shift[i+1] = curr_shift[i] + curr_count[i];
    orig_shift[i+1] = orig_shift[i] + orig_count[i];
  }

  for (int i = 0; i < s_comm + 1; i++) {
    curr_shift[i] *= 2;
    orig_shift[i] *= 2;
  }

  for (int i = 0; i < s_comm; i++) {
    orig_count[i] *= 2;
  }

  int *curr_loc      = (int *) malloc (sizeof(int) * curr_shift[s_comm]);
  int *orig_loc      = (int *) malloc (sizeof(int) * orig_shift[s_comm]);
  int *idxCurrToBuff = (int *) malloc (sizeof(int) * n_boxes);

  for (int i = 0; i < s_comm; i++) {
    curr_count[i] = 0;
  }

  for (int i = 0; i < n_boxes; i++) {
    int iProc = origin[3*i    ];
    int iPart = origin[3*i + 1];
    int iElt  = origin[3*i + 2];
    int idx   = curr_shift[iProc] + curr_count[iProc];

    idxCurrToBuff[i] = idx/2;

    curr_loc[idx++] = iPart;
    curr_loc[idx++] = iElt;

    curr_count[iProc] += 2;
  }

  PDM_MPI_Alltoallv(curr_loc, curr_count, curr_shift, PDM_MPI_INT,
                orig_loc, orig_count, orig_shift, PDM_MPI_INT,
                boxes->comm);

  free (curr_loc);

  for (int i = 0; i < s_comm+1; i++) {
    curr_shift[i] = curr_shift[i]/2;
    orig_shift[i] = orig_shift[i]/2;
  }

  for (int i = 0; i < s_comm; i++) {
    curr_count[i] = curr_count[i]/2;
    orig_count[i] = orig_count[i]/2;
  }

  /* Send variable stride to current distribution
   *   - Exchange <- origin : stride
   *   - Sort curr_stride to obtain current_distrib_stride */

  int *curr_stride = NULL;
  int *orig_stride = NULL;

  if (t_stride == PDM_STRIDE_VAR) {

    curr_stride = (int *) malloc (sizeof(int) * curr_shift[s_comm]);
    orig_stride = (int *) malloc (sizeof(int) * orig_shift[s_comm]);

    for (int i = 0; i < orig_shift[s_comm]; i++) {
      int iPart = orig_loc[2*i    ];
      int iElt  = orig_loc[2*i + 1];
      orig_stride[i] = origin_distrib_stride[iPart][iElt];
    }

    PDM_MPI_Alltoallv(orig_stride, orig_count, orig_shift, PDM_MPI_INT,
                  curr_stride, curr_count, curr_shift, PDM_MPI_INT,
                  boxes->comm);

    if (alloc_distrib_stride) {

      int *_current_distrib_stride = (int *) malloc (sizeof(int) * orig_shift[s_comm]);
      *current_distrib_stride = _current_distrib_stride;

      for (int i = 0; i < curr_shift[s_comm]; i++) {
        _current_distrib_stride[i] = curr_stride[idxCurrToBuff[i]];
      }
    }

  }

  /* Send data to current distribution
   *   - Exhange <- origin : data
   *   - Sort curr_data to obtain current_distrib_data */

  if (t_stride == PDM_STRIDE_VAR) {

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] = 0;
      orig_count[i] = 0;
    }

    for (int i = 0; i < s_comm+1; i++) {
      for (int k = curr_shift[i]; k < curr_shift[i+1]; k++) {
        curr_count[i] += curr_stride[k];
      }
      for (int k = orig_shift[i]; k < orig_shift[i+1]; k++) {
        orig_count[i] += orig_stride[k];
      }
    }

    int n_origin = orig_shift[s_comm];

    for (int i = 0; i < s_comm + 1; i++) {
      curr_shift[i] = 0;
      orig_shift[i] = 0;
    }

    for (int i = 0; i < s_comm; i++) {
      curr_shift[i+1] = curr_shift[i] + curr_count[i];
      orig_shift[i+1] = orig_shift[i] + orig_count[i];
    }

    for (int i = 0; i < s_comm + 1; i++) {
      curr_shift[i] *= (int) data_size;
      orig_shift[i] *= (int) data_size;
    }

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] *= (int) data_size;
    }

    for (int i = 0; i < s_comm; i++) {
      orig_count[i] = 0;
    }

    unsigned char *orig_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * orig_shift[s_comm]);

    unsigned char *curr_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * curr_shift[s_comm]);

    int **_origin_distrib_idx = (int **) malloc (sizeof(int) * boxes->n_part_orig);
    for (int i = 0; i < boxes->n_part_orig; i++) {
      _origin_distrib_idx[i] = (int *) malloc (sizeof(int) * (boxes->n_boxes_orig[i] + 1));
      for (int k = 0; k < boxes->n_boxes_orig[i] + 1; k++) {
        _origin_distrib_idx[i][k] = 0;
      }
      for (int k = 0; k < boxes->n_boxes_orig[i]; k++) {
        _origin_distrib_idx[i][k+1] = _origin_distrib_idx[i][k] + origin_distrib_stride[i][k];
      }
    }

    for (int i = 0; i < n_origin; i++) {
      int iPart = orig_loc[2*i    ];
      int iElt  = orig_loc[2*i + 1];

      int s_block = origin_distrib_stride[iPart][iElt] * (int) data_size;
      int idx = orig_shift[i] + orig_count[i];
      int idx1 = _origin_distrib_idx[iPart][iElt];

      for (int k = 0; k < s_block; k++) {
        orig_data[idx++] = _origin_distrib_data[iPart][idx1++];
      }

    }

    for (int i = 0; i < boxes->n_part_orig; i++) {
      free (_origin_distrib_idx[i]);
    }
    free (_origin_distrib_idx);

    PDM_MPI_Alltoallv(orig_data, orig_count, orig_shift, PDM_MPI_UNSIGNED_CHAR,
                  curr_data, curr_count, curr_shift, PDM_MPI_UNSIGNED_CHAR,
                  boxes->comm);

    unsigned char *_current_distrib_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                                     * curr_shift[s_comm]);

    int *curr_data_idx       = (int *) malloc (sizeof(int) * (boxes->n_boxes + 1));
    int *current_distrib_idx = (int *) malloc (sizeof(int) * (boxes->n_boxes + 1));

    for (int i = 0; i < boxes->n_boxes + 1; i++) {
      curr_data_idx[i]       = 0;
      current_distrib_idx[i] = 0;
    }

    for (int i = 0; i < boxes->n_boxes; i++) {
      curr_data_idx[i+1]       = curr_data_idx[i] + curr_stride[i];
      current_distrib_idx[i+1] = current_distrib_idx[i] + (*current_distrib_stride)[i];
    }

    for (int i = 0; i < boxes->n_boxes + 1; i++) {
      curr_data_idx[i]       *= (int) data_size;
      current_distrib_idx[i] *= (int) data_size;
    }

    for (int i = 0; i < boxes->n_boxes; i++) {

      int s_block = (*current_distrib_stride)[i] * (int) data_size;
      int idx     = current_distrib_idx[i];
      int idx1    = curr_data_idx[idxCurrToBuff[i]];

      for (int k = 0; k < s_block; k++) {
        _current_distrib_data[idx++] = curr_data[idx1++];
      }
    }

    /* Clean up */

    free (orig_data);
    free (curr_data);
    free (orig_stride);
    free (curr_stride);
    free (curr_data_idx);
    free (current_distrib_idx);

    *current_distrib_data = (void *) _current_distrib_data;

  }

  else {

    int s_block = stride_cst * (int) data_size;

    unsigned char *orig_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * orig_shift[s_comm]
                                                         * s_block);

    unsigned char *curr_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * curr_shift[s_comm]
                                                         * s_block);

    for (int i = 0; i < s_comm+1; i++) {
      curr_shift[i] *= s_block;
      orig_shift[i] *= s_block;
    }

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] *= s_block;
      orig_count[i] *= s_block;
    }

    for (int i = 0; i < orig_shift[s_comm]; i++) {
      int iPart = orig_loc[2*i    ];
      int iElt  = orig_loc[2*i + 1];
      for (int k = 0; k < s_block; k++) {
        orig_data[s_block * i + k] = _origin_distrib_data[iPart][s_block * iElt + k];
      }
    }

    PDM_MPI_Alltoallv(orig_data, orig_count, orig_shift, PDM_MPI_UNSIGNED_CHAR,
                  curr_data, curr_count, curr_shift, PDM_MPI_UNSIGNED_CHAR,
                  boxes->comm);

    unsigned char *_current_distrib_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                                     * curr_shift[s_comm]
                                                                     * s_block);

    for (int i = 0; i < curr_shift[s_comm]; i++) {
      for (int k = 0; k < s_block; k++) {
        _current_distrib_data[s_block * i + k] = curr_data[s_block * idxCurrToBuff[i] + k];
      }
    }

    free (orig_data);
    free (curr_data);
    *current_distrib_data = (void *) _current_distrib_data;

  }

  /* Clean up */

  if (curr_stride != NULL)
    free (curr_stride);

  if (orig_stride != NULL)
    free (orig_stride);

  free (curr_count);
  free (curr_shift);
  free (orig_count);
  free (orig_shift);
  free (orig_loc);
  free (idxCurrToBuff);

}


/*----------------------------------------------------------------------------
 * Send data to origin for any box
 *
 * parameters:
 *   box_set                <-- pointer to the PDM_box_t structure
 *   t_stride               <-- Type of stride
 *   stride_cst             <-- Constant stride
 *   data_size              <-- Size of data
 *   current_distrib_stride <-- Current stride distribution
 *   current_distrib_data   <-- Current data distribution
 *   origin_distrib         --> Origin distribution (ipart, num_loc) for any box
 *   origin_distrib_stride  --> Origin stride distribution
 *   origin_distrib_data    --> Origin data distribution
 *
 * return : size of origin_distrib_data
 *----------------------------------------------------------------------------*/

void
PDM_box_set_send_data_to_origin_distrib
(
 PDM_box_set_t  *boxes,
 PDM_stride_t   t_stride,
 int            stride_cst,
 size_t         data_size,
 int           *current_distrib_stride,
 void          *current_distrib_data,
 int          **origin_distrib_stride,
 void         **origin_distrib_data
)
{
  /*
   * Send origin properties to the origin process :
   *   - Compute the number element to send for any process
   *   - Exchange -> origin these numbers all_to_all
   *   - Exchange -> origin origin properties (part, num loc) all_to_all_v
   */

  int n_boxes = boxes->n_boxes;
  int *origin = boxes->origin;

  int s_comm = PDM_MPI_Comm_rank (boxes->comm, &s_comm);

  int *curr_count = (int *) malloc (sizeof(int) * s_comm);
  for (int i = 0; i < s_comm; i++) {
    curr_count[i] = 0;
  }

  for (int i = 0; i < n_boxes; i++) {
    curr_count[origin[3*i]] += 1;
  }

  int *orig_count = (int *) malloc (sizeof(int) * s_comm);
  PDM_MPI_Alltoall(curr_count, 1, PDM_MPI_INT,
               orig_count, 1, PDM_MPI_INT, boxes->comm);

  int *curr_shift = (int *) malloc (sizeof(int) * (s_comm + 1));
  int *orig_shift = (int *) malloc (sizeof(int) * (s_comm + 1));

  for (int i = 0; i < s_comm + 1; i++) {
    curr_shift[i] = 0;
    orig_shift[i] = 0;
  }

  for (int i = 0; i < s_comm; i++) {
    curr_shift[i+1] = curr_shift[i] + curr_count[i];
    orig_shift[i+1] = orig_shift[i] + orig_count[i];
  }

  for (int i = 0; i < s_comm + 1; i++) {
    curr_shift[i] *= 2;
    orig_shift[i] *= 2;
  }

  for (int i = 0; i < s_comm; i++) {
    orig_count[i] *= 2;
  }

  int *curr_loc      = (int *) malloc (sizeof(int) * curr_shift[s_comm]);
  int *orig_loc      = (int *) malloc (sizeof(int) * orig_shift[s_comm]);

  for (int i = 0; i < s_comm; i++) {
    curr_count[i] = 0;
  }

  for (int i = 0; i < n_boxes; i++) {
    int iProc = origin[3*i    ];
    int iPart = origin[3*i + 1];
    int iElt  = origin[3*i + 2];
    int idx   = curr_shift[iProc] + curr_count[iProc];

    curr_loc[idx++] = iPart;
    curr_loc[idx++] = iElt;

    curr_count[iProc] += 2;
  }

  PDM_MPI_Alltoallv(curr_loc, curr_count, curr_shift, PDM_MPI_INT,
                orig_loc, orig_count, orig_shift, PDM_MPI_INT,
                boxes->comm);

  free (curr_loc);

  for (int i = 0; i < s_comm+1; i++) {
    curr_shift[i] = curr_shift[i]/2;
    orig_shift[i] = orig_shift[i]/2;
  }

  for (int i = 0; i < s_comm; i++) {
    curr_count[i] = curr_count[i]/2;
    orig_count[i] = orig_count[i]/2;
  }

  /*
   * Send variable stride to current distribution
   *   - Exchange -> origin : stride
   */

  int *curr_stride = NULL;
  int *orig_stride = NULL;

  if (t_stride == PDM_STRIDE_VAR) {

    orig_stride = (int *) malloc (sizeof(int) * orig_shift[s_comm]);
    curr_stride = (int *) malloc (sizeof(int) * curr_shift[s_comm]);

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] = 0;
    }

    for (int i = 0; i < n_boxes; i++) {
      int iProc = origin[3*i    ];
      int idx   = curr_shift[iProc] + curr_count[iProc];

      curr_stride[idx++] = current_distrib_stride[i];
      curr_count[iProc] += 1;
    }

    PDM_MPI_Alltoallv(curr_stride, curr_count, curr_shift, PDM_MPI_INT,
                  orig_stride, orig_count, orig_shift, PDM_MPI_INT,
                  boxes->comm);

    for (int i = 0; i < orig_shift[s_comm]; i++) {
      int iPart = orig_loc[2*i];
      int iElt  = orig_loc[2*i+1];
      origin_distrib_stride[iPart][iElt] = orig_stride[i];
    }

    free (curr_stride);
  }

  /*
   * Send data to current distribution
   *   - Exhange -> origin : data
   */

  if (t_stride == PDM_STRIDE_VAR) {

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] = 0;
      orig_count[i] = 0;
    }

    for (int i = 0; i < s_comm+1; i++) {
      for (int k = curr_shift[i]; k < curr_shift[i+1]; k++) {
        curr_count[i] += curr_stride[k];
      }
      for (int k = orig_shift[i]; k < orig_shift[i+1]; k++) {
        orig_count[i] += orig_stride[k];
      }
    }

    for (int i = 0; i < s_comm + 1; i++) {
      curr_shift[i] = 0;
      orig_shift[i] = 0;
    }

    for (int i = 0; i < s_comm; i++) {
      curr_shift[i+1] = curr_shift[i] + curr_count[i];
      orig_shift[i+1] = orig_shift[i] + orig_count[i];
    }

    for (int i = 0; i < s_comm + 1; i++) {
      curr_shift[i] *= (int) data_size;
      orig_shift[i] *= (int) data_size;
    }

    for (int i = 0; i < s_comm; i++) {
      orig_count[i] *= (int) data_size;
    }

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] = 0;
    }

    unsigned char *curr_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * curr_shift[s_comm]);

    unsigned char *orig_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * orig_shift[s_comm]);

    int *current_distrib_idx = (int *) malloc (sizeof(int) * (boxes->n_boxes + 1));
    for (int i = 0; i < boxes->n_boxes + 1; i++) {
      current_distrib_idx[i] = 0;
    }

    for (int i = 0; i < boxes->n_boxes; i++) {
      current_distrib_idx[i+1] = current_distrib_idx[i] + current_distrib_stride[i];
    }

    for (int i = 0; i < boxes->n_boxes + 1; i++) {
      current_distrib_idx[i] *= (int) data_size;
    }

    unsigned char *_current_distrib_data = (unsigned char *) current_distrib_data;

    for (int i = 0; i < n_boxes; i++) {
      int iProc = origin[3*i    ];
      int idx   = curr_shift[iProc] + curr_count[iProc];
      int idx1  = current_distrib_idx[i];

      int s_block = current_distrib_stride[i] * (int) data_size;
      for (int k = 0; k < s_block; k++) {
        curr_data[idx++] = _current_distrib_data[idx1++];
      }
      curr_count[iProc] += s_block;
    }

    free (current_distrib_idx);

    PDM_MPI_Alltoallv(curr_data, curr_count, curr_shift, PDM_MPI_UNSIGNED_CHAR,
                  orig_data, orig_count, orig_shift, PDM_MPI_UNSIGNED_CHAR,
                  boxes->comm);

    int **_origin_distrib_idx = (int **) malloc (sizeof(int) * boxes->n_part_orig);
    for (int i = 0; i < boxes->n_part_orig; i++) {
      _origin_distrib_idx[i]   = (int *) malloc (sizeof(int) * (boxes->n_boxes_orig[i] + 1));
      for (int k = 0; k < boxes->n_boxes_orig[i] + 1; k++) {
        _origin_distrib_idx[i][k] = 0;
      }
      for (int k = 0; k < boxes->n_boxes_orig[i]; k++) {
        _origin_distrib_idx[i][k+1] = _origin_distrib_idx[i][k] + origin_distrib_stride[i][k];
      }
    }

    for (int i = 0; i < boxes->n_part_orig; i++) {
      unsigned char *_origin_distrib_data =
        (unsigned char *) malloc (sizeof(unsigned char)
                                  * _origin_distrib_idx[i][boxes->n_boxes_orig[i]]);
      origin_distrib_data[i] = (void *) _origin_distrib_data;
    }

    int idx1  = 0;
    for (int i = 0; i < orig_shift[s_comm]; i++) {
      int iPart = orig_loc[2*i];
      int iElt  = orig_loc[2*i+1];

      int idx   = _origin_distrib_idx[iPart][iElt];

      int s_block = current_distrib_stride[i] * (int) data_size;
      unsigned char *_origin_distrib_data =  origin_distrib_data[iPart];
      for (int k = 0; k < s_block; k++) {
        _origin_distrib_data[idx++] = orig_data[idx1++];
      }
    }

    /* Clean up */

    for (int i = 0; i < boxes->n_part_orig; i++) {
      free (_origin_distrib_idx[i]);
    }
    free (_origin_distrib_idx);

    free (orig_data);
    free (curr_data);
    free (orig_stride);
    free (current_distrib_idx);

  }

  else {

    int s_block = stride_cst * (int) data_size;

    for (int i = 0; i < s_comm + 1; i++) {
      curr_shift[i] *= s_block;
      orig_shift[i] *= s_block;
    }

    for (int i = 0; i < s_comm; i++) {
      orig_count[i] *= s_block;
    }

    unsigned char *orig_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * orig_shift[s_comm]);
    unsigned char *curr_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                         * curr_shift[s_comm]);

    for (int i = 0; i < s_comm; i++) {
      curr_count[i] = 0;
    }

    unsigned char *_current_distrib_data = (unsigned char *) current_distrib_data;

    int idx1 = 0;
    for (int i = 0; i < n_boxes; i++) {
      int iProc = origin[3*i    ];
      int idx   = curr_shift[iProc] + curr_count[iProc];

      for (int k = 0; k < s_block; k++) {
        curr_data[idx++] = _current_distrib_data[idx1++];
      }
      curr_count[iProc] += s_block;
    }

    PDM_MPI_Alltoallv(curr_data, curr_count, curr_shift, PDM_MPI_UNSIGNED_CHAR,
                  orig_data, orig_count, orig_shift, PDM_MPI_UNSIGNED_CHAR,
                  boxes->comm);


    for (int i = 0; i < boxes->n_part_orig; i++) {
      unsigned char *_origin_distrib_data = (unsigned char *) malloc (sizeof(unsigned char)
                                                                      * s_block
                                                                      * boxes->n_boxes_orig[i]);
      origin_distrib_data[i] = (void *) _origin_distrib_data;
    }

    idx1  = 0;
    for (int i = 0; i < orig_shift[s_comm]; i++) {
      int iPart = orig_loc[2*i];
      int iElt  = orig_loc[2*i+1];

      int idx   = iElt * s_block;

      unsigned char *_origin_distrib_data =  (unsigned char *) origin_distrib_data[iPart];
      for (int k = 0; k < s_block; k++) {
        _origin_distrib_data[idx++] = orig_data[idx1++];
      }
    }

    /*
     * Clean up
     */

    free (orig_data);
    free (curr_data);

  }

  /*
   * Clean up
   */

  free (orig_count);
  free (curr_count);
  free (orig_shift);
  free (curr_shift);
  free (curr_loc);
  free (orig_loc);

}

/*----------------------------------------------------------------------------
 * Create a PDM_box_distrib_t structure.
 *
 * parameters:
 *   n_boxes   <-- number of boxes
 *   n_g_boxes <-- global number of boxes
 *   max_level <-- max level reached locally in the related tree
 *   comm      <-- MPI communicator. on which the distribution takes place
 *
 * returns:
 *   a pointer to a new allocated PDM_box_distrib_t structure.
 *---------------------------------------------------------------------------*/

PDM_box_distrib_t *
PDM_box_distrib_create(int  n_boxes,
                       PDM_g_num_t  n_g_boxes,
                       int        max_level,
                       PDM_MPI_Comm   comm)
{
  int  i, n_ranks, gmax_level;

  PDM_box_distrib_t  *new_distrib = NULL;

  if (n_g_boxes == 0)
    return NULL;

  new_distrib = (PDM_box_distrib_t *) malloc(sizeof(PDM_box_distrib_t));

  /* Parallel parameters */

  PDM_MPI_Comm_size(comm, &n_ranks);

  new_distrib->n_ranks = n_ranks;

  new_distrib->n_boxes = n_boxes;

  assert(n_ranks > 1);

  new_distrib->morton_index = (PDM_morton_code_t *) malloc((n_ranks + 1) * sizeof(PDM_morton_code_t));

  PDM_MPI_Allreduce(&max_level, &gmax_level, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  new_distrib->max_level = gmax_level;
  new_distrib->fit = 999.0;

  new_distrib->index = (int *) malloc((n_ranks + 1) * sizeof(int));

  for (i = 0; i < n_ranks + 1; i++)
    new_distrib->index[i] = 0;

  new_distrib->list = NULL;

  return  new_distrib;
}

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_distrib_t structure.
 *
 * parameters:
 *   distrib <-> pointer to pointer to the structure to destroy
 *---------------------------------------------------------------------------*/

void
PDM_box_distrib_destroy(PDM_box_distrib_t  **distrib)
{
  if (distrib != NULL) {

    PDM_box_distrib_t  *d = *distrib;

    if (d == NULL)
      return;

    free(d->index);
    free(d->list);
    free(d->morton_index);

    free(d);
  }
}

/*----------------------------------------------------------------------------
 * Delete redundancies in box distribution
 *
 * parameters:
 *   distrib <->  pointer to the PDM_box_distrib_t structure
 *---------------------------------------------------------------------------*/

void
PDM_box_distrib_clean(PDM_box_distrib_t  *distrib)
{
  int  i, rank;

  int   *counter = NULL, *new_index = NULL;

  counter = (int *) malloc(distrib->n_boxes * sizeof(int));
  new_index = (int *) malloc((distrib->n_ranks + 1) * sizeof(int));

  for (i = 0; i < distrib->n_ranks + 1; i++)
    new_index[i] = 0;

  for (rank = 0; rank < distrib->n_ranks; rank++) {

    int   shift = new_index[rank];
    int   start = distrib->index[rank];
    int   end = distrib->index[rank + 1];

    if (end - start > 0) {

      for (i = 0; i < distrib->n_boxes; i++)
        counter[i] = 0;

      for (i = start; i < end; i++)
        counter[distrib->list[i]] += 1;

      for (i = 0; i < distrib->n_boxes; i++) {

        if (counter[i] > 0)
          distrib->list[shift++] = i;

      }

    } /* end if end - start > 0 */

    new_index[rank+1] = shift;

  } /* End of loop on ranks */

  /* Memory management */

  free(distrib->index);
  distrib->list = (int *) realloc(distrib->list, new_index[distrib->n_ranks] * sizeof(int));
  distrib->index = new_index;

  free(counter);
}

/*----------------------------------------------------------------------------
 * Display a histogram on leaves associated to the boxes and several
 * other pieces of information (min, max, ...)
 *
 * parameters:
 *   distrib <-- pointer to the PDM_box_distrib_t structure
 *   comm    <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
PDM_box_distrib_dump_statistics(const PDM_box_distrib_t  *distrib,
                                PDM_MPI_Comm                   comm)
{
  int  i;

  int  n_ranks = 0;
  int  n_quantiles = 5;
  int  quantile_start[6];
  int  n_boxes[5];

  /* Sanity checks */

  assert(distrib != NULL);
  assert(distrib->index != NULL);

  _get_distrib_statistics(distrib,
                          &n_quantiles,
                          quantile_start,
                          n_boxes,
                          NULL,
                          &n_ranks,
                          comm);

  PDM_printf("\n"
             "- Box distribution statistics -\n\n");

  PDM_printf("   Distribution imbalance:              %10.4g\n",
             distrib->fit);
  PDM_printf("   Number of ranks in distribution:     %8d\n\n",
             n_ranks);

  /* Print histogram to show the distribution of boxes */

  if (n_quantiles > 0) {

    for (i = 0; i < n_quantiles - 1; i++)
      PDM_printf("    %3d : [ %10d ; %10d [ = %10d\n",
                 i+1, quantile_start[i], quantile_start[i+1], n_boxes[i]);

    i = n_quantiles -1;
    PDM_printf("    %3d : [ %10d ; %10d ] = %10d\n",
               i+1, quantile_start[i], quantile_start[i+1] - 1, n_boxes[i]);

  }
  fflush(stdout);
}


/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
