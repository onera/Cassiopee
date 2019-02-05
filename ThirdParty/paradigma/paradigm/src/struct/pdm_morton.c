/*============================================================================
 * Morton encoding for 2D or 3D coordinates.
 *============================================================================*/
/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA
  Copyright (C) 2008-2010  EDF

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_morton.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  PDM_morton_distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int PDM_morton_distrib_n_iter_max = 5;

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

static const int _3d_children[8][3] = {{0, 0, 0},    /* child 1 */
                                       {0, 0, 1},    /* 2 */
                                       {0, 1, 0},    /* 3 */
                                       {0, 1, 1},    /* 4 */
                                       {1, 0, 0},    /* 5 */
                                       {1, 0, 1},    /* 6 */
                                       {1, 1, 0},    /* 7 */
                                       {1, 1, 1}};   /* 8 */

static const int _2d_children[4][2] = {{0, 0},   /* child 1 */
                                       {0, 1},   /* 2 */
                                       {1, 0},   /* 3 */
                                       {1, 1}};  /* 4 */

static const int _1d_children[2][1] = {{0},   /* child 1 */
                                       {1}};  /* 2 */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Transorm local extents to global extents.
 *
 * parameters:
 *   dim       <-- spatial dimension (1, 2, or 3)
 *   g_extents <-> global extents (size: dim*2)
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static void
_local_to_global_extents(int         dim,
                         double  extents[],
                         PDM_MPI_Comm    comm)
{
  int i;
  double  l_min[3], l_max[3];

  for (i = 0; i < dim; i++) {
    l_min[i] = extents[i];
    l_max[i] = extents[i + dim];
  }

  PDM_MPI_Allreduce(l_min, extents, dim, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_max, extents + dim, dim, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater or equal to Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_a_ge_b(PDM_morton_code_t  code_a,
        PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[2] << a_diff;
    code_a.X[1] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = code_b.X[0] << b_diff;
    code_b.X[1] = code_b.X[2] << b_diff;
    code_b.X[1] = code_b.X[2] << b_diff;
  }

  i = l - 1;
  while (i > 0) {
    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

  a =   ((code_a.X[0] >> i) % 2) * 4
      + ((code_a.X[1] >> i) % 2) * 2
      + ((code_a.X[2] >> i) % 2);
  b =   ((code_b.X[0] >> i) % 2) * 4
      + ((code_b.X[1] >> i) % 2) * 2
      + ((code_b.X[2] >> i) % 2);

  return (a >= b) ? true : false;
}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater than Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_a_gt_b(PDM_morton_code_t  code_a,
        PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[2] << a_diff;
    code_a.X[1] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = code_b.X[0] << b_diff;
    code_b.X[1] = code_b.X[2] << b_diff;
    code_b.X[1] = code_b.X[2] << b_diff;
  }

  i = l - 1;
  while (i > 0) {
    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

  a =   ((code_a.X[0] >> i) % 2) * 4
      + ((code_a.X[1] >> i) % 2) * 2
      + ((code_a.X[2] >> i) % 2);
  b =   ((code_b.X[0] >> i) % 2) * 4
      + ((code_b.X[1] >> i) % 2) * 2
      + ((code_b.X[2] >> i) % 2);

  return (a > b) ? true : false;
}

/*----------------------------------------------------------------------------
 * Build a heap structure or order a heap structure.
 *
 * parameters:
 *  parent       <--  parent id in the Morton code list
 *  n_codes      <--  number of codes to work with
 *  morton_codes <->  list of Morton codes to work with
 *----------------------------------------------------------------------------*/

static void
_descend_morton_heap(PDM_g_num_t        parent,
                     int                n_codes,
                     PDM_morton_code_t  morton_codes[])
{
  PDM_morton_code_t  tmp;
  PDM_g_num_t child = 2 * parent + 1;

  while (child < n_codes) {

    if (child + 1 < n_codes)
      if (_a_gt_b(morton_codes[child + 1], morton_codes[child]))
        child++;

    if (_a_ge_b(morton_codes[parent], morton_codes[child])) return;

    tmp = morton_codes[parent];
    morton_codes[parent] = morton_codes[child];
    morton_codes[child] = tmp;
    parent = child;
    child = 2 * parent + 1;

  } /* End of while */

}

/*----------------------------------------------------------------------------
 * Convert a double into a Morton code.
 *
 * parameters:
 *   dim       <-- 2D or 3D
 *   input     <-- double to convert
 *   level     <-- level of the grid on which the code has to be built
 *
 * returns:
 *  a Morton code associated to the input.
 *----------------------------------------------------------------------------*/

inline static PDM_morton_code_t
_double_to_code(int     dim,
                double  input,
                int     level)
{
  int  l, child_id;
  PDM_morton_code_t  code;
  double coords[3] = {0.0, 0.0, 0.0};
  double l_mult = 1.0;

  const int max_level = 15; /* no more than 52 bits in mantissa / 3 */

  /* Build associated Morton code */

  code.L = max_level;

  if (input <= 0.0) {
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
  }

  else if (input >= 1.0) {
    coords[0] = 1.0;
    coords[1] = 1.0;
    coords[2] = 1.0;
  }

  else if (dim == 3) {
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*8);
      if (child_id > 7) child_id = 7;
      input = input*8 - child_id;
      coords[0] += child_id/4 * l_mult;
      coords[1] += (child_id%4)/2 * l_mult;
      coords[2] += child_id%2 * l_mult;
    }
  }

  else if (dim == 2) {
    coords[2] = 0;
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*4);
      if (child_id > 3) child_id = 3;
      input = input*4 - child_id;
      coords[0] += child_id/2 * l_mult;
      coords[1] += child_id%2 * l_mult;
    }
  }

  else if (dim == 1) {
    coords[1] = 0;
    coords[2] = 0;
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*2);
      if (child_id > 1) child_id = 1;
      input = input*2 - child_id;
      coords[0] += child_id * l_mult;
    }
  }

  code = PDM_morton_encode(dim, level, coords);

  return code;
}

/*----------------------------------------------------------------------------
 * Build a heap structure or order a heap structure with a working array
 * to save the ordering.
 *
 * parameters:
 *  parent       <-- parent id in the Morton code list
 *  n_codes      <-- number of codes to work with
 *  morton_codes <-- list of Morton codes to work with
 *  order        <-> working array to save the ordering
 *----------------------------------------------------------------------------*/

static void
_descend_morton_heap_with_order(PDM_g_num_t                 parent,
                                int                 n_codes,
                                const PDM_morton_code_t   morton_codes[],
                                int                *order)
{
  int          tmp;
  PDM_g_num_t  child = 2 * parent + 1;

  while (child < n_codes) {

    if (child + 1 < n_codes) {
      if (_a_gt_b(morton_codes[order[child + 1]],
                  morton_codes[order[child]]))
        child++;
    }

    if (_a_ge_b(morton_codes[order[parent]],
                morton_codes[order[child]]))
      return;

    tmp = order[parent];
    order[parent] = order[child];
    order[child] = tmp;
    parent = child;
    child = 2 * parent + 1;

  } /* End while */
}

/*----------------------------------------------------------------------------
 * Evaluate a distribution array.
 *
 * parameters:
 *   n_ranges     <-- Number of ranges in the distribution
 *   distribution <-- Number of elements associated to each range of
 *                    the distribution
 *   optim        <-- Optimal count in each range
 *
 * returns:
 *   a fit associated to the distribution. If fit = 0,
 *   distribution is perfect.
 *----------------------------------------------------------------------------*/

static double
_evaluate_distribution(int          n_ranges,
                       PDM_g_num_t   *distribution,
                       double       optim)
{
  int  i;
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = PDM_MAX(d_up, distribution[i] - optim);
    else
      d_low = PDM_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  /* if (cs_glob_rank_id <= 0) */
  /*   PDM_printf("<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n", */
  /*              optim, fit); */
#endif

  return  fit;
}

/*----------------------------------------------------------------------------
 * Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 * parameters:
 *   dim          <-- 2D or 3D
 *   n_ranks      <-- number of ranks (= number of ranges)
 *   gmax_level   <-- level on which Morton encoding is build
 *   gsum_weight  <-- global sum of all weightings
 *   n_codes      <-- local number of Morton codes
 *   morton_codes <-- local list of Morton codes to distribute
 *   weight       <-- weighting related to each code
 *   order        <-- ordering array
 *   sampling     <-- sampling array
 *   c_freq       <-> pointer to the cumulative frequency array
 *   g_distrib    <-> pointer to a distribution array
 *   comm         <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_define_rank_distrib(int                      dim,
                     int                      n_ranks,
                     int                      gmax_level,
                     PDM_g_num_t                gsum_weight,
                     int                n_codes,
                     const PDM_morton_code_t  morton_codes[],
                     const int          weight[],
                     const int          order[],
                     const double             sampling[],
                     double                   cfreq[],
                     PDM_g_num_t                g_distrib[],
                     PDM_MPI_Comm                 comm)
{
  int  id, rank_id;
  PDM_morton_code_t  sample_code;
  int   i;

  int  bucket_id = 1;
  PDM_g_num_t   *l_distrib = NULL;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;

  /* Initialization */

  l_distrib = (PDM_g_num_t *) malloc(n_samples * sizeof(PDM_g_num_t));

  for (id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  /* morton_codes are supposed to be ordered */

  sample_code = _double_to_code(dim, sampling[bucket_id], gmax_level);

  for (i = 0; i < n_codes; i++) {

    PDM_g_num_t   o_id = order[i];

    if (_a_ge_b(sample_code, morton_codes[o_id]))
      l_distrib[bucket_id - 1] += weight[o_id];

    else {

      while (_a_gt_b(morton_codes[o_id], sample_code)) {

        bucket_id++;
        assert(bucket_id < n_samples + 1);

        sample_code = _double_to_code(dim, sampling[bucket_id], gmax_level);
      }

      l_distrib[bucket_id - 1] += weight[o_id];

    }

  } /* End of loop on elements */

  /* Define the global distribution */

  PDM_MPI_Allreduce(l_distrib, g_distrib, n_samples, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM,
                comm);

  free(l_distrib);

  /* Define the cumulative frequency related to g_distribution */

  cfreq[0] = 0.;
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  for (id = 0; id < n_samples; id++)
    cfreq[id+1] = cfreq[id] + (double)g_distrib[id]/(double)gsum_weight;
  cfreq[n_samples] = 1.0;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

#if 0 && defined(DEBUG) && !defined(DEBUG) /* For debugging purpose only */

  /* if (cs_glob_rank_id <= 0) { */

  /*   FILE  *dbg_file = NULL; */
  /*   char  *rfilename = NULL; */
  /*   int  len; */
  /*   static int  loop_id1 = 0; */

  /*   len = strlen("DistribOutput_l.dat")+1+2; */
  /*   rfilename = (char *) malloc(len * sizeof(char)); */
  /*   sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1); */

  /*   loop_id1++; */

  /*   dbg_file = fopen(rfilename, "w"); */

  /*   fprintf(dbg_file, */
  /*           "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |" */
  /*           "Global Distrib\n"); */
  /*   for (i = 0; i < n_samples; i++) */
  /*     fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n", */
  /*             i, (double)i/(double)n_samples, cfreq[i], */
  /*             sampling[i], g_distrib[i]); */
  /*   fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n", */
  /*           i, 1.0, 1.0, 1.0, 0); */

  /*   fclose(dbg_file); */
  /*   free(rfilename); */

  /* } */

#endif /* debugging output */

  /* Convert global distribution from n_samples to n_ranks */

  for (rank_id = 0; rank_id < n_ranks; rank_id++) {

    PDM_g_num_t   sum = 0;
    int   shift = rank_id * sampling_factor;

    for (id = 0; id < sampling_factor; id++)
      sum += g_distrib[shift + id];
    g_distrib[rank_id] = sum;

  } /* End of loop on ranks */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Sanity check in debug */
  {
    PDM_g_num_t   sum = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      sum += g_distrib[rank_id];

    if (sum != gsum_weight) {
      PDM_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
      abort();
    }
  }
#endif /* sanity check */

}

/*----------------------------------------------------------------------------
 * Update a distribution associated to sampling to assume a well-balanced
 * distribution of the leaves of the tree.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   n_ranks  <-- number of ranks (= number of ranges)
 *   c_freq   <-> cumulative frequency array
 *   sampling <-> pointer to pointer to a sampling array
 *   comm     <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_update_sampling(int      dim,
                 int      n_ranks,
                 double   c_freq[],
                 double  *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  double  s_low, s_high;

  double  *new_sampling = NULL, *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute new_sampling */

  new_sampling = (double *) malloc((n_samples + 1) * sizeof(double));

  new_sampling[0] = _sampling[0];
  next_id = 1;

  for (i = 0; i < n_samples; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (j = next_id; j < n_samples + 1; j++) {
      if (c_freq[j] >= target_freq) {
        next_id = j;
        break;
      }
    }

    /* Find new s such as new_s is equal to target_freq by
       a linear interpolation */

    f_low = c_freq[next_id-1];
    f_high = c_freq[next_id];

    s_low = _sampling[next_id-1];
    s_high = _sampling[next_id];

    if (f_high - f_low > 0) {
      delta = (target_freq - f_low) * (s_high - s_low) / (f_high - f_low);
      new_sampling[i+1] = s_low + delta;
    }
    else /* f_high = f_low */
      new_sampling[i+1] = s_low + 0.5 * (s_low + s_high);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    /* PDM_printf(" <_update_distrib> (rank: %d) delta: %g, target: %g," */
    /*            " next_id: %d, f_low: %g, f_high: %g, s_low: %g, s_high: %g\n" */
    /*            "\t => new_sampling: %g\n", */
    /*            cs_glob_rank_id, delta, target_freq, next_id, */
    /*            f_low, f_high, s_low, s_high, new_sampling[i+1]); */
#endif

  } /* End of loop on samples */

  new_sampling[n_samples] = 1.0;

  free(_sampling);

  /* Return pointers */

  *sampling = new_sampling;
}

/*----------------------------------------------------------------------------
 * Compute a sampling array which assumes a well-balanced distribution of
 * leaves of the tree among the ranks.
 *
 * parameters:
 *   dim          <--  2D or 3D
 *   n_ranks      <--  number of ranks
 *   gmax_level   <--  level on which Morton encoding is build
 *   n_codes      <--  local number of Morton ids
 *   morton_codes <--  local list of Morton ids to distribute
 *   weight       <--  weighting related to each code
 *   order        <--  ordering array
 *   sampling     <-> pointer to pointer to a sampling array
 *   comm         <--  mpi communicator
 *
 * returns:
 *   fit associated to the returned sampling array
 *----------------------------------------------------------------------------*/

static double
_bucket_sampling(int                      dim,
                 int                      n_ranks,
                 int                     gmax_level,
                 int                n_codes,
                 const PDM_morton_code_t  morton_codes[],
                 const int          weight[],
                 const int          order[],
                 double                  *sampling[],
                 PDM_MPI_Comm                 comm)
{
  int  i, n_iters;
  int   j;
  double  fit, best_fit, optim;

  PDM_g_num_t   lsum_weight = 0, gsum_weight = 0;
  PDM_g_num_t   *distrib = NULL;
  double  *cfreq = NULL, *best_sampling = NULL;
  double  *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute the global number of elements and the optimal number of elements
     on each rank */

  for (j = 0; j < n_codes; j++)
    lsum_weight += weight[j];

  PDM_MPI_Allreduce(&lsum_weight, &gsum_weight, 1,  PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  optim = (double)gsum_weight / (double)n_ranks;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

  /* Define a naive sampling (uniform distribution) */

  for (i = 0; i < n_samples + 1; i++)
    _sampling[i] = i*unit;

  /* Define the distribution associated to the current sampling array */

  distrib = (PDM_g_num_t *) malloc(n_samples * sizeof(PDM_g_num_t));
  cfreq = (double *) malloc((n_samples + 1) * sizeof(double));

  _define_rank_distrib(dim,
                       n_ranks,
                       gmax_level,
                       gsum_weight,
                       n_codes,
                       morton_codes,
                       weight,
                       order,
                       _sampling,
                       cfreq,
                       distrib,
                       comm);

  /* Initialize best choice */

  fit = _evaluate_distribution(n_ranks, distrib, optim);
  best_fit = fit;

  best_sampling = malloc((n_samples + 1) * sizeof(double));

  for (i = 0; i < n_samples + 1; i++)
    best_sampling[i] = _sampling[i];

  /* Loop to get a better sampling array */

  for (n_iters = 0;
       (   n_iters < PDM_morton_distrib_n_iter_max
        && fit > PDM_morton_distrib_tol);
       n_iters++)  {

    _update_sampling(dim, n_ranks, cfreq, &_sampling);

    /* Compute the new distribution associated to the new sampling */

    _define_rank_distrib(dim,
                         n_ranks,
                         gmax_level,
                         gsum_weight,
                         n_codes,
                         morton_codes,
                         weight,
                         order,
                         _sampling,
                         cfreq,
                         distrib,
                         comm);

    fit = _evaluate_distribution(n_ranks, distrib, optim);

    /* Save the best sampling array and its fit */

    if (fit < best_fit) {

      best_fit = fit;
      for (i = 0; i < n_samples + 1; i++)
        best_sampling[i] = _sampling[i];

    }

  } /* End of while */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  /* if (cs_glob_rank_id <= 0) */
  /*   PDM_printf("\n  <_bucket_sampling> n_iter: %d, opt: %g, best_fit: %g\n", */
  /*              n_iters, optim, best_fit); */
#endif

  /* Free memory */

  free(cfreq);
  free(distrib);
  free(_sampling);

  *sampling = best_sampling;

  return best_fit;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine the global extents associated with a set of coordinates
 *
 * parameters:
 *   dim       <-- spatial dimension
 *   n_coords  <-- local number of coordinates
 *   coords    <-- entity coordinates; size: n_entities*dim (interlaced)
 *   g_extents --> global extents (size: dim*2)
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
PDM_morton_get_coord_extents(int               dim,
                             size_t            n_coords,
                             const double  coords[],
                             double        g_extents[],
                             PDM_MPI_Comm          comm)
{
  size_t  i, j;

  /* Get global min/max coordinates */

  for (j = 0; j < (size_t)dim; j++) {
    g_extents[j]       = DBL_MAX;
    g_extents[j + dim] = -DBL_MAX;
  }

  for (i = 0; i < n_coords; i++) {
    for (j = 0; j < (size_t)dim; j++) {
      if (coords[i*dim + j] < g_extents[j])
        g_extents[j] = coords[i*dim + j];
      if (coords[i*dim + j] > g_extents[j + dim])
        g_extents[j + dim] = coords[i*dim + j];
    }
  }

  if (comm != PDM_MPI_COMM_NULL)
    _local_to_global_extents(dim, g_extents, comm);

}

/*----------------------------------------------------------------------------
 * Determine the global extents associated with a set of local extents
 *
 * parameters:
 *   dim       <-- spatial dimension
 *   n_extents <-- local number of coordinates
 *   extents   <-- entity coordinates; size: n_entities*dim*2 (interlaced)
 *   g_extents --> global extents (size: dim*2)
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
PDM_morton_get_global_extents(int               dim,
                              size_t            n_extents,
                              const double  extents[],
                              double        g_extents[],
                              PDM_MPI_Comm          comm)
{
  size_t  i, j;

  /* Get global min/max coordinates */

  for (i = 0; i < (size_t)dim; i++) {
    g_extents[i]       = DBL_MAX;
    g_extents[i + dim] = -DBL_MAX;
  }

  for (i = 0; i < n_extents; i++) {
    for (j = 0; j < (size_t)dim; j++) {
      g_extents[j]     = PDM_MIN(g_extents[j],
                                extents[i*dim*2 + j]);
      g_extents[j+dim] = PDM_MAX(g_extents[j + dim],
                                extents[i*dim*2 + j + dim]);
    }
  }

  if (comm != PDM_MPI_COMM_NULL)
    _local_to_global_extents(dim, g_extents, comm);

}

/*----------------------------------------------------------------------------
 * Build a Morton code according to the level in an octree grid and its
 * coordinates in the grid.
 *
 * parameters:
 *   dim    <-- 1D, 2D or 3D
 *   level  <-- level in the grid
 *   coords <-- coordinates in the grid (normalized)
 *
 * returns:
 *  a Morton code
 *----------------------------------------------------------------------------*/

PDM_morton_code_t
PDM_morton_encode(int               dim,
                  PDM_morton_int_t  level,
                  const double  coords[])
{
  int  i;
  PDM_morton_code_t  morton_code;

  PDM_morton_int_t  refinement = 1 << level;

  morton_code.L = level;

  /* Initialize last components for 1D or 2D case */

  morton_code.X[1] = 0;
  morton_code.X[2] = 0;

  for (i = 0; i < dim; i++)
    morton_code.X[i] = (PDM_morton_int_t) PDM_MIN(floor(coords[i]*refinement), refinement - 1);

  return morton_code;
}

/*----------------------------------------------------------------------------
 * Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   level    <-- level in the grid
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   m_code   --> array of corresponding Morton codes
 *----------------------------------------------------------------------------*/

void
PDM_morton_encode_coords(int                dim,
                         PDM_morton_int_t   level,
                         const double   extents[],
                         size_t             n_coords,
                         const double   coords[],
                         PDM_morton_code_t  m_code[])
{
  size_t i, j;
  double s[3], d[3], n[3];
  double d_max = 0.0;

  PDM_morton_int_t  refinement = 1 << level;

  for (i = 0; i < (size_t)dim; i++) {
    s[i] = extents[i];
    d[i] = extents[i+dim] - extents[i];
    d_max = PDM_MAX(d_max, d[i]);
  }

  for (i = 0; i < (size_t)dim; i++) { /* Reduce effective dimension */
    if (d[i] < d_max * 1e-10)
      d[i] = d_max * 1e-10;
  }

  switch(dim) {

  case 3:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 3; j++) {
        n[j] = (coords[i*dim + j] - s[j]) / d[j];
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
    }
    break;

  case 2:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 2; j++) {
        n[j] = (coords[i*dim + j] - s[j]) / d[j];
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
      m_code[i].X[2] = 0;
    }
    break;

  case 1:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      n[0] = (coords[i] - s[0]) / d[0];
      m_code[i].X[0] = (PDM_morton_int_t) PDM_MIN(floor(n[0]*refinement), refinement - 1);
      m_code[i].X[1] = 0;
      m_code[i].X[2] = 0;
    }
    break;

  default:
    assert(dim > 0 && dim < 4);
    break;
  }
}

/*----------------------------------------------------------------------------
 * Given a Morton code in the grid, compute the Morton codes of its
 * children when refining the grid by one level.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   parent   <-- Morton code associated with parent
 *   children --> array of children Morton codes
 *                (size: 8 in 3D, 4 in 2D, 2 in 1D)
 *----------------------------------------------------------------------------*/

void
PDM_morton_get_children(int                dim,
                        PDM_morton_code_t  parent,
                        PDM_morton_code_t  children[])
{
  int  i;
  PDM_morton_code_t  anchor;

  if (dim == 3) {

    for (i = 0; i < 3; i++)
      anchor.X[i] = 2 * parent.X[i];

    for (i = 0; i < 8; i++) {
      children[i].L = parent.L + 1;
      children[i].X[0] = anchor.X[0] + _3d_children[i][0];
      children[i].X[1] = anchor.X[1] + _3d_children[i][1];
      children[i].X[2] = anchor.X[2] + _3d_children[i][2];
    }

  }
  else if (dim == 2) {

    for (i = 0; i < 2; i++)
      anchor.X[i] = 2 * parent.X[i];

    for (i = 0; i < 4; i++) {
      children[i].L = parent.L + 1;
      children[i].X[0] = anchor.X[0] + _2d_children[i][0];
      children[i].X[1] = anchor.X[1] + _2d_children[i][1];
      children[i].X[2] = 0;
    }

  }

  else if (dim == 1) {

    anchor.X[0] = 2 * parent.X[0];

    for (i = 0; i < 2; i++) {
      children[i].L = parent.L + 1;
      children[i].X[0] = anchor.X[0] + _1d_children[i][0];
      children[i].X[1] = 0;
      children[i].X[2] = 0;
    }

  }
}

/*----------------------------------------------------------------------------
 * Locally order a list of Morton ids.
 *
 * parameters:
 *   n_codes      <-- number of Morton ids to order
 *   morton_codes <-- array of Morton ids to order
 *   order        --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
PDM_morton_local_order(int                n_codes,
                       const PDM_morton_code_t  morton_codes[],
                       int                order[])
{
  int   i;
  int   tmp;

  assert(n_codes == 0 || morton_codes != NULL);
  assert(n_codes == 0 || order != NULL);

  for (i = 0; i < n_codes; i++)
    order[i] = i;

  /* Build heap */

  for (i = n_codes/2 - 1; (int)i >= 0; i--)
    _descend_morton_heap_with_order(i,  n_codes, morton_codes, order);

  /* Sort array */

  for (i = n_codes - 1; (int)i >= 0; i--) {

    tmp = order[0];
    order[0] = order[i];
    order[i] = tmp;

    _descend_morton_heap_with_order(0, i, morton_codes, order);

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)   /* Check ordering */
  for (i = 1; i < n_codes; i++) {
    if (_a_gt_b(morton_codes[order[i-1]], morton_codes[order[i]])) {
      PDM_error(__FILE__, __LINE__, 0,
              "Id: %u inconsistent: bad ordering of Morton codes.",
              (unsigned)i);
      abort();
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 * Locally sort a list of Morton ids.
 *
 * parameters:
 *   n_codes      <-- number of Morton ids to order
 *   morton_codes <-> array of Morton ids to sort
 *----------------------------------------------------------------------------*/

void
PDM_morton_local_sort(int          n_codes,
                      PDM_morton_code_t  morton_codes[])
{
  int   i;
  PDM_morton_code_t  tmp;

  /* Build heap */

  for (i = n_codes/2 - 1; (int)i >= 0; i--)
    _descend_morton_heap(i, n_codes, morton_codes);

  /* Sort array */

  for (i = n_codes - 1; (int)i >= 0; i--) {

    tmp = morton_codes[0];
    morton_codes[0] = morton_codes[i];
    morton_codes[i] = tmp;

    _descend_morton_heap(0, i, morton_codes);

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)   /* Check good ordering */
  for (i = 1; i < n_codes; i++) {
    if (_a_gt_b(dim, morton_codes[i - 1], morton_codes[i])) {
      PDM_error(__FILE__, __LINE__, 0,
              "Id: %u inconsistent: bad ordering of Morton codes.",
              (unsigned)i);
      abort();
    }
  }
#endif

}

/*----------------------------------------------------------------------------
 * Compare two Morton encodings and check if these two codes are equal,
 * different or shared the same anchor.
 *
 * parameters:
 *   dim    <--  2D or 3D
 *   code_a <--  first Morton code to compare
 *   code_b <--  second Morton code to compare
 *
 * returns:
 *  a type on the kind of relation between the two Morton encodings.
 *----------------------------------------------------------------------------*/

PDM_morton_compare_t
PDM_morton_compare(int                dim,
                   PDM_morton_code_t  code_a,
                   PDM_morton_code_t  code_b)
{
  int i;

  if (code_a.L == code_b.L) {

    for (i = 0; i < dim; i++)
      if (code_a.X[i] != code_b.X[i])
        return PDM_MORTON_DIFFERENT_ID;
    return PDM_MORTON_EQUAL_ID;

  }
  else {

    if (code_a.L < code_b.L) {

      PDM_morton_int_t  delta = code_b.L - code_a.L;

      for (i = 0; i < dim; i++)
        code_a.X[i] = code_a.X[i] << delta;

    }
    else {

      PDM_morton_int_t  delta = code_a.L - code_b.L;

      for (i = 0; i < dim; i++)
        code_b.X[i] = code_b.X[i] << delta;

    }

    for (i = 0; i < dim; i++)
      if (code_a.X[i] != code_b.X[i])
        return PDM_MORTON_DIFFERENT_ID;
    return PDM_MORTON_SAME_ANCHOR;

  }

}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater than Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

_Bool
PDM_morton_a_gt_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b)
{
  return  _a_gt_b(a, b);
}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater or equal to Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

_Bool
PDM_morton_a_ge_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b)
{
  return  _a_ge_b(a, b);
}

/*----------------------------------------------------------------------------
 * Get the index associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the array.
 *
 * parameters:
 *   size  <-- size of the array
 *   code  <-- code we are searching for
 *   codes <-- array of Morton codes
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

int
PDM_morton_binary_search(int           size,
                         PDM_morton_code_t   code,
                         PDM_morton_code_t  *codes)
{
  int start = 0;
  int end = size;

  while (end - start > 1) {

    int  middle = (end - start)/2 + start;

    if (_a_gt_b(codes[middle], code))
      end = middle;
    else
      start = middle;
  }

  return start;
}

/*----------------------------------------------------------------------------
 * Get the quantile associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * parameters:
 *   n_quantiles    <-- number of quantiles
 *   code           <-- code we are searching for
 *   quantile_start <-- first Morton code in each quantile (size: n_quantiles)
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

size_t
PDM_morton_quantile_search(size_t              n_quantiles,
                           PDM_morton_code_t   code,
                           PDM_morton_code_t  *quantile_start)
{
  size_t mid_id = 0;
  size_t start_id = 0;
  size_t end_id = n_quantiles;

  /* use binary search */

  while (start_id + 1 < end_id) {
    mid_id = start_id + ((end_id -start_id) / 2);
    if (_a_gt_b(quantile_start[mid_id], code))
      end_id = mid_id;
    else
      start_id = mid_id;
  }

  /* We may have stopped short of the required value,
     or have multiple occurences of a quantile start
     (in case of empty quantiles), of which we want to
     find the find highest one */

  while (   start_id < n_quantiles - 1
         && _a_ge_b(code, quantile_start[start_id+1]))
    start_id++;

  return start_id;
}

/*----------------------------------------------------------------------------
 * Build a global Morton encoding rank index.
 *
 * The rank_index[i] contains the first Morton code assigned to rank [i].
 *
 * parameters:
 *   dim         <-- 1D, 2D or 3D
 *   gmax_level  <-- level in octree used to build the Morton encoding
 *   n_codes     <-- number of Morton codes to be indexed
 *   morton_code <-- array of Morton codes to be indexed
 *   weight      <-- weighting related to each code
 *   order       <-- ordering array
 *   rank_index  <-> pointer to the global Morton encoding rank index
 *   comm        <-- MPI communicator on which we build the global index
 *
 * returns:
 *  the fit related to the Morton encoding distribution (lower is better).
 *----------------------------------------------------------------------------*/

double
PDM_morton_build_rank_index(int                      dim,
                            int                      gmax_level,
                            PDM_l_num_t                n_codes,
                            const PDM_morton_code_t  code[],
                            const int          weight[],
                            const int          order[],
                            PDM_morton_code_t        rank_index[],
                            PDM_MPI_Comm                 comm)
{
  int  i, id, rank_id, n_ranks, n_samples;
  double  best_fit;

  double  *sampling = NULL;

  const int  sampling_factor = _sampling_factors[dim];

  /* Allocations and Initialization */

  PDM_MPI_Comm_size(comm, &n_ranks);

  n_samples = sampling_factor * n_ranks;

  sampling = (double *) malloc((n_samples + 1) * sizeof(double));

  for (i = 0; i < n_samples + 1; i++)
    sampling[i] = 0.0;

  best_fit = _bucket_sampling(dim,
                              n_ranks,
                              gmax_level,
                              n_codes,
                              code,
                              weight,
                              order,
                              &sampling,
                              comm);

  /* Define Morton index */

  for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {

    id = rank_id * sampling_factor;
    rank_index[rank_id] = _double_to_code(dim, sampling[id], gmax_level);

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  { /* Dump Morton index and associated sampling on rank 0 */
    PDM_printf("\nMorton rank index:\n\n");
    for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
      id = sampling_factor * rank_id;
      PDM_printf("rank: %5d (sampling: %7.4g)- ", rank_id, sampling[id]);
      PDM_morton_dump(dim, rank_index[rank_id]);

    }
    PDM_printf("\n");
    fflush(stdout);
  }
#endif

  /* Free memory */

  free(sampling);

  return best_fit;
}

/*----------------------------------------------------------------------------
 * Dump a Morton to standard output or to a file.
 *
 * parameters:
 *   dim  <-- 2D or 3D
 *   code <-- Morton code to dump
 *----------------------------------------------------------------------------*/

void
PDM_morton_dump(int                dim,
                PDM_morton_code_t  code)
{
  int  i;
  double  coord[3];

  const PDM_g_num_t   n = 1 << code.L;
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  const double  stride = 1/(double)n;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

  for (i = 0; i < dim; i++)
    coord[i] = stride * code.X[i];

  if (dim == 3)
    PDM_printf("Morton Code:\n"
               "L =  %3u [X, Y, Z] - [%5u %5u %5u]"
               "[%6.5lf %6.5lf %6.5lf]\n",
               code.L, code.X[0], code.X[1], code.X[2],
               coord[0], coord[1], coord[2]);

  else if (dim == 2)
    PDM_printf("Morton Code\n"
               "L =  %3u [X, Y] - [%5u %5u] [%6.5lf %6.5lf]\n",
               code.L, code.X[0], code.X[1], coord[0], coord[1]);

  fflush(stdout);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

