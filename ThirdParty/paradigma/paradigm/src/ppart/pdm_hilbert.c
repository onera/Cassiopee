/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_hilbert.h"
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

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  pdm_hilbert_distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int pdm_hilbert_distrib_n_iter_max = 5;

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

static const unsigned  _imax = ~(0U);

/* 2 dimension to nkey conversion */
static const unsigned  _cs_idata2d[]
= {0, 3, 1, 2,
   0, 1, 3, 2,
   2, 3, 1, 0,
   2, 1, 3, 0};

/* 2 dimension to nkey state transitions */
static const unsigned  _cs_istate2d[]
= {1, 2, 0, 0,
   0, 1, 3, 1,
   2, 0, 2, 3,
   3, 3, 1, 2};

/**
 * Zhang
 */

static const unsigned  _zhang_idata2d[]
= {
    0,  1,  3,  2,
    0,  2,  3,  1, 
    3,  2,  0,  1,
    3,  1,  0,  2};

static const unsigned  _zhang_istate2d[]
= {
    1,  0,  0,  3,
    0,  1,  1,  2,
    3,  2,  2,  1,
    2,  3,  3,  0};

/**
 * An Introduction with Applications in Scientific Computing
 * BADER Michael
 * 8.1 3D Hilbert Curves p.115
 */

static const unsigned _bader_idata2d[]
= {
    0,  1,  3,  2,
    0,  2,  3,  1, 
    3,  1,  0,  2,
    3,  2,  0,  1};

static const unsigned _bader_istate2d[]
= {
    1,  0,  0,  2,
    0,  1,  1,  3,
    3,  2,  2,  0,
    2,  3,  3,  1
    };

/* 3 dimension to nkey conversion */
static const unsigned  _cs_idata3d[]
= {0,  7,  3,  4,  1,  6,  2,  5,
   0,  1,  3,  2,  7,  6,  4,  5,
   0,  3,  7,  4,  1,  2,  6,  5,
   2,  3,  5,  4,  1,  0,  6,  7,
   4,  5,  3,  2,  7,  6,  0,  1,
   4,  7,  3,  0,  5,  6,  2,  1,
   6,  7,  5,  4,  1,  0,  2,  3,
   0,  1,  7,  6,  3,  2,  4,  5,
   2,  1,  5,  6,  3,  0,  4,  7,
   6,  1,  5,  2,  7,  0,  4,  3,
   0,  7,  1,  6,  3,  4,  2,  5,
   2,  1,  3,  0,  5,  6,  4,  7,
   4,  7,  5,  6,  3,  0,  2,  1,
   4,  5,  7,  6,  3,  2,  0,  1,
   6,  1,  7,  0,  5,  2,  4,  3,
   0,  3,  1,  2,  7,  4,  6,  5,
   2,  3,  1,  0,  5,  4,  6,  7,
   6,  7,  1,  0,  5,  4,  2,  3,
   2,  5,  1,  6,  3,  4,  0,  7,
   4,  3,  7,  0,  5,  2,  6,  1,
   4,  3,  5,  2,  7,  0,  6,  1,
   6,  5,  1,  2,  7,  4,  0,  3,
   2,  5,  3,  4,  1,  6,  0,  7,
   6,  5,  7,  4,  1,  2,  0,  3};

/* 3 dimension to nkey state transitions */
static const unsigned  _cs_istate3d[]
= { 1,  6,  3,  4,  2,  5,  0,  0,
    0,  7,  8,  1,  9,  4,  5,  1,
   15, 22, 23, 20,  0,  2, 19,  2,
    3, 23,  3, 15,  6, 20, 16, 22,
   11,  4, 12,  4, 20,  1, 22, 13,
   22, 12, 20, 11,  5,  0,  5, 19,
   17,  0,  6, 21,  3,  9,  6,  2,
   10,  1, 14, 13, 11,  7, 12,  7,
    8,  9,  8, 18, 14, 12, 10, 11,
   21,  8,  9,  9,  1,  6, 17,  7,
    7, 17, 15, 12, 16, 13, 10, 10,
   11, 14,  9,  5, 11, 22,  0,  8,
   18,  5, 12, 10, 19,  8, 12, 20,
    8, 13, 19,  7,  5, 13, 18,  4,
   23, 11,  7, 17, 14, 14,  6,  1,
    2, 18, 10, 15, 21, 19, 20, 15,
   16, 21, 17, 19, 16,  2,  3, 18,
    6, 10, 16, 14, 17, 23, 17, 15,
   18, 18, 21,  8, 17,  7, 13, 16,
    3,  4, 13, 16, 19, 19,  2,  5,
   16, 13, 20, 20,  4,  3, 15, 12,
    9, 21, 18, 21, 15, 14, 23, 10,
   22, 22,  6,  1, 23, 11,  4,  3,
   14, 23,  2,  9, 22, 23, 21,  0};

/**
 * Numerical Simulation in molecular dynamics
 * Author : Michael Griebel, Stephan Knapek, Gerhard Zumbush
 * 8.4 Parallel Tree Methods p. 363
 */

static const unsigned  _griebel_idata3d[]
= {0,  7,  3,  4,  1,  6,  2,  5,
   4,  3,  7,  0,  5,  2,  6,  1,
   6,  1,  5,  2,  7,  0,  4,  3,
   2,  5,  1,  6,  3,  4,  0,  7, 
   0,  1,  7,  6,  3,  2,  4,  5,
   6,  7,  1,  0,  5,  4,  2,  3,
   2,  3,  5,  4,  1,  0,  6,  7,
   4,  5,  3,  2,  7,  6,  0,  1,
   0,  3,  1,  2,  7,  4,  6,  5,
   2,  1,  3,  0,  5,  6,  4,  7,
   4,  7,  5,  6,  3,  0,  2,  1,
   6,  5,  7,  4,  1,  2,  0,  3};

static const unsigned  _griebel_istate3d[]
= { 8,  10,  3,  3,  4,  5,  4,  5,
    2,  2,  11,  9,  4,  5,  4,  5,
    7,  6,  7,  6,  8,  10,  1,  1,
    7,  6,  7,  6,  0,  0,  11,  9,
    0,  8,  1,  11,  6,  8,  6,  11,
    10,  0,  9,  1,  10,  7,  9,  7,
    10,  4,  9,  4,  10,  2,  9,  3,
    5,  8,  5,  11,  2,  8,  3,  11,
    4,  9,  0,  0,  7,  9,  2,  2,
    1,  1,  8,  5,  3,  3,  8,  6,
    11,  5,  0,  0,  11,  6,  2,  2,
    1,  1,  4,  10,  3,  3,  7,  10};
    
/**
 * Study on Pseudo-Hilbert Scan and Its Application to
 * 		HDR Tone Mapping
 * ZHANG Jian
 * 2.3 3-D Pseudo-Hilbert Scan for Arbitrarily-sized Cuboid Regions p.37
 */

static const unsigned  _zhang_idata3d[]
= {
    0,  1,  3,  2,  6,  7,  5,  4,
    6,  7,  5,  4,  0,  1,  3,  2,
    5,  7,  6,  4,  0,  2,  3,  1,
    3,  1,  0,  2,  6,  4,  5,  7,
    0,  1,  5,  4,  6,  7,  3,  2,
    6,  2,  3,  7,  5,  1,  0,  4,
    5,  1,  0,  4,  6,  2,  3,  7,
    3,  2,  6,  7,  5,  4,  0,  1,
    0,  4,  6,  2,  3,  7,  5,  1,
    6,  4,  0,  2,  3,  1,  5,  7,
    5,  7,  3,  1,  0,  2,  6,  4,
    3,  7,  5,  1,  0,  4,  6,  2};

static const unsigned  _zhang_istate3d[]
= {
    8,  4,  4,  3,  3,  5,  5,  10,
    9,  5,  5,  2,  2,  4,  4,  11,
    6,  10,  10,  1,  1,  8,  8,  7,
    7,  11,  11,  0,  0,  9,  9,  6,
    8,  0,  0,  6,  6,  1,  1,  11,
    1,  9,  9,  7,  7,  10,  10,  0,
    2,  10,  10,  4,  4,  9,  9,  3,
    11,  3,  3,  5,  5,  2,  2,  8,
    0,  4,  4,  9,  9,  7,  7,  2,
    5,  1,  1,  8,  8,  3,  3,  6,
    6,  2,  2,  11,  11,  0,  0,  5,
    3,  7,  7,  10,  10,  4,  4,  1};

/**
 * An Introduction with Applications in Scientific Computing
 * BADER Michael
 * 8.1 3D Hilbert Curves p.115
 */

static const unsigned  _bader_idata3d[]
= {
    0,  4,  6,  2,  3,  7,  5,  1,
    5,  1,  3,  7,  6,  2,  0,  4,
    0,  2,  3,  1,  5,  7,  6,  4,
    5,  7,  6,  4,  0,  2,  3,  1,
    0,  1,  5,  4,  6,  7,  3,  2,
    5,  4,  0,  1,  3,  2,  6,  7,
    3,  7,  5,  1,  0,  4,  6,  2,
    6,  2,  0,  4,  5,  1,  3,  7,
    3,  1,  0,  2,  6,  4,  5,  7,
    6,  4,  5,  7,  3,  1,  0,  2,
    3,  2,  6,  7,  5,  4,  0,  1,
    6,  7,  3,  2,  0,  1,  5,  4};

static const unsigned  _bader_istate3d[]
= {
    2,  4,  4,  7,  7,  10,  10,  3,
    3,  5,  5,  6,  6,  11,  11,  2,
    4,  0,  0,  8,  8,  1,  1,  11,
    5,  1,  1,  9,  9,  0,  0,  10,
    0,  2,  2,  5,  5,  9,  9,  6,
    1,  3,  3,  4,  4,  8,  8,  7,
    8,  10,  10,  1,  1,  4,  4,  9,
    9,  11,  11,  0,  0,  5,  5,  8,
    10,  6,  6,  2,  2,  7,  7,  5,
    11,  7,  7,  3,  3,  6,  6,  4,
    6,  8,  8,  11,  11,  3,  3,  0,
    7,  9,  9,  10,  10,  2,  2,  1};


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build a Hilbert key based on a 1-d coordinate in [0, 1].
 *
 * parameters:
 *   coord   <-- 1-d coordinate, normalized
 *
 * returns:
 *   associated Hilbert encoding
 *----------------------------------------------------------------------------*/

static PDM_hilbert_code_t
_hilbert_encode_1d(const double        coord[1])
{
  return coord[0];
}

/*----------------------------------------------------------------------------
 * Build a Hilbert key based on 2-d coordinates in [0, 1].
 *
 * parameters:
 *   coord   <-- 2-d coordinates, normalized
 *
 * returns:
 *   associated Hilbert encoding
 *----------------------------------------------------------------------------*/

static PDM_hilbert_code_t
_hilbert_encode_2d(const double  coord[2], const unsigned * _idata2d, const unsigned * _istate2d)
{
  int level;
  unsigned int c[2], temp, state;
  unsigned int key[2] = {0, 0};
  const int maxlevel = 28; /* 28 bits of significance per dimension */

  const unsigned *d[]
    = {_idata2d,  _idata2d+4, _idata2d+8, _idata2d+12};
  const unsigned *s[]
    ={_istate2d, _istate2d+4, _istate2d+8, _istate2d+12};

  assert(coord[0] >= 0.0 && coord[0] <= 1.0);
  assert(coord[1] >= 0.0 && coord[1] <= 1.0);

  /* convert x, y coordinates to integers in range [0, imax] */
  c[0] = (unsigned int) (coord[0] * (double) _imax);     /* x */
  c[1] = (unsigned int) (coord[1] * (double) _imax);     /* y */

  /* use state tables to convert nested quadrant's coordinates level by level */

  state = 0;
  for (level = 0; level < maxlevel; level++) {
    temp = (  (c[0]  >> (30-level)) & 2)   /* extract 2 bits at current level */
            | ((c[1] >> (31-level)) & 1);

    /* treat key[] as long shift register, shift in converted coordinate */
    key[0] = (key[0] << 2) | (key[1] >> 30);
    key[1] = (key[1] << 2) | *(d[state] + temp);

    state = *(s[state] + temp);
  }

  /* Convert 2 part Hilbert key to double and return;
     Note that maxlevel could be increased from 28 to 32
     by using long doubles (with a double, we have 56 significant bits,
     which allows for 28 bits per coordinate). This could be increased
     further by using 64-bit integers in intermediate calculations. */

  return ldexp ((double) key[0], -24)  +  ldexp ((double) key[1], -56);
}

/*----------------------------------------------------------------------------
 * Build a Hilbert key based on 3-d coordinates in [0, 1].
 *
 * parameters:
 *   coord   <-- 3-d coordinates, normalized
 *
 * returns:
 *   associated Hilbert encoding
 *----------------------------------------------------------------------------*/

static PDM_hilbert_code_t
_hilbert_encode_3d(const double  coord[3], const unsigned * _idata3d, const unsigned * _istate3d)
{
  int level;
  unsigned int c[3], temp, state;
  unsigned int key[3] = {0, 0, 0};
  const int maxlevel = 19; /* 32 bits of significance per dimension */

  const unsigned int *d[]
    = {_idata3d,     _idata3d+8,   _idata3d+16,  _idata3d+24,
       _idata3d+32,  _idata3d+40,  _idata3d+48,  _idata3d+56,
       _idata3d+64,  _idata3d+72,  _idata3d+80,  _idata3d+88,
       _idata3d+96,  _idata3d+104, _idata3d+112, _idata3d+120,
       _idata3d+128, _idata3d+136, _idata3d+144, _idata3d+152,
       _idata3d+160, _idata3d+168, _idata3d+176, _idata3d+184};

  const unsigned int *s[]
    = {_istate3d,     _istate3d+8,   _istate3d+16,  _istate3d+24,
       _istate3d+32,  _istate3d+40,  _istate3d+48,  _istate3d+56,
       _istate3d+64,  _istate3d+72,  _istate3d+80,  _istate3d+88,
       _istate3d+96,  _istate3d+104, _istate3d+112, _istate3d+120,
       _istate3d+128, _istate3d+136, _istate3d+144, _istate3d+152,
       _istate3d+160, _istate3d+168, _istate3d+176, _istate3d+184};

  assert(coord[0] >= 0.0 && coord[0] <= 1.0);
  assert(coord[1] >= 0.0 && coord[1] <= 1.0);
  assert(coord[2] >= 0.0 && coord[2] <= 1.0);

  /* convert x,y,z coordinates to integers in range [0, _imax] */
  c[0] = (unsigned int) (coord[0] * (double) _imax);     /* x */
  c[1] = (unsigned int) (coord[1] * (double) _imax);     /* y */
  c[2] = (unsigned int) (coord[2] * (double) _imax);     /* z */

  /* use state tables to convert nested quadrant's coordinates level by level */
  key[0] = 0; key[1] = 0; key[2] = 0;
  state = 0;
  for (level = 0; level < maxlevel; level++) {
    temp = (  (c[0]  >> (29-level)) & 4)  /* extract 3 bits at current level */
            | ((c[1] >> (30-level)) & 2)
            | ((c[2] >> (31-level)) & 1);

    /* treat key[] as long shift register, shift in converted coordinate */
    key[0] = (key[0] << 3) |  (key[1] >> 29);
    key[1] = (key[1] << 3) | *(d[state] + temp);

    state = *(s[state] + temp);
  }

  /* Convert 2 part Hilbert key to double and return;
     Note that maxlevel could be increased from 19 to 32 by using
     a 3-part key and long doubles (with a double, we have 56 significant
     bits, which allows for 19 bits per coordinate). This could be increased
     further by using 64-bit integers in intermediate calculations. */

  return ldexp ((double) key[0], -25)  +  ldexp ((double) key[1], -57);
}

static double
_hilbert_encode_3d_others(const double  coord[3], const unsigned * _idata3d, const unsigned * _istate3d)
{
  int level;
  unsigned int c[3], temp, state;
  unsigned int key[2] = {0, 0};
  const int maxlevel = 19; /* 32 bits of significance per dimension */

  const unsigned int *d[]
    = {_idata3d,     _idata3d+8,   _idata3d+16,  _idata3d+24,
       _idata3d+32,  _idata3d+40,  _idata3d+48,  _idata3d+56,
       _idata3d+64,  _idata3d+72,  _idata3d+80,  _idata3d+88};

  const unsigned int *s[]
    = {_istate3d,     _istate3d+8,   _istate3d+16,  _istate3d+24,
       _istate3d+32,  _istate3d+40,  _istate3d+48,  _istate3d+56,
       _istate3d+64,  _istate3d+72,  _istate3d+80,  _istate3d+88};

  assert(coord[0] >= 0.0 && coord[0] <= 1.0);
  assert(coord[1] >= 0.0 && coord[1] <= 1.0);
  assert(coord[2] >= 0.0 && coord[2] <= 1.0);

  /* convert x,y,z coordinates to integers in range [0, _imax] */
  c[0] = (unsigned int) (coord[0] * (double) _imax);     /* x */
  c[1] = (unsigned int) (coord[1] * (double) _imax);     /* y */
  c[2] = (unsigned int) (coord[2] * (double) _imax);     /* z */

  /* use state tables to convert nested quadrant's coordinates level by level */
  state = 0;
  for (level = 0; level < maxlevel; level++) {
    temp = (  (c[0]  >> (29-level)) & 4)  /* extract 3 bits at current level */
            | ((c[1] >> (30-level)) & 2)
            | ((c[2] >> (31-level)) & 1);

    /* treat key[] as long shift register, shift in converted coordinate */
    key[0] = (key[0] << 3) |  (key[1] >> 29);
    key[1] = (key[1] << 3) | *(d[state] + temp);
    state = *(s[state] + temp);
  }

  /* Convert 2 part Hilbert key to double and return;
 *      Note that maxlevel could be increased from 19 to 32 by using
 *           a 3-part key and long doubles (with a double, we have 56 significant
 *                bits, which allows for 19 bits per coordinate). This could be increased
 *                     further by using 64-bit integers in intermediate calculations. */

  return ldexp ((double) key[0], -25)  +  ldexp ((double) key[1], -57);
}

/*----------------------------------------------------------------------------
 * Transform local extents to global extents.
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

  PDM_MPI_Allreduce(l_min, extents, dim, PDM__PDM_MPI_REAL, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_max, extents + dim, dim, PDM__PDM_MPI_REAL, PDM_MPI_MAX, comm);
}


/*----------------------------------------------------------------------------
 * Build a heap structure or order a heap structure with a working array
 * to save the ordering.
 *
 * parameters:
 *  parent        <-- parent id in the Hilbert code list
 *  n_codes       <-- number of codes to work with
 *  hilbert_codes <-- array of Hilbert codes to work with
 *  order         <-> working array to save the ordering
 *----------------------------------------------------------------------------*/

static void
_descend_hilbert_heap(PDM_g_num_t           parent,
                      int                  n_codes,
                      const PDM_hilbert_code_t   hilbert_codes[],
                      int                 *order)
{
  int   tmp;
  PDM_g_num_t   child = 2 * parent + 1;

  while (child < n_codes) {

    if (child + 1 < n_codes) {
      if (hilbert_codes[order[child + 1]] > hilbert_codes[order[child]])
        child++;
    }

    if (hilbert_codes[order[parent]] >= hilbert_codes[order[child]])
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
      d_up = _MAX(d_up, distribution[i] - optim);
    else
      d_low = _MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    PDM_printf( "<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}

/*----------------------------------------------------------------------------
 * Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 * parameters:
 *   dim           <-- 2D or 3D
 *   n_ranks       <-- number of ranks (= number of ranges)
 *   gsum_weight   <-- global sum of all weightings
 *   n_codes       <-- local number of Hilbert codes
 *   hilbert_codes <-- local list of Hilbert codes to distribute
 *   weight        <-- weighting related to each code
 *   order         <-- ordering array
 *   sampling      <-- sampling array
 *   c_freq        <-> pointer to the cumulative frequency array
 *   g_distrib     <-> pointer to a distribution array
 *   comm          <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_define_rank_distrib(int                       dim,
                     int                       n_ranks,
                     PDM_g_num_t                 gsum_weight,
                     int                 n_codes,
                     const PDM_hilbert_code_t  hilbert_codes[],
                     const int           weight[],
                     const int           order[],
                     const PDM_hilbert_code_t  sampling[],
                     double                    cfreq[],
                     PDM_g_num_t                 g_distrib[],
                     PDM_MPI_Comm                  comm)
{
  int  id, rank_id;
  PDM_hilbert_code_t  sample_code;
  int   i;

  int  bucket_id = 1;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;

  /* Initialization */

  PDM_g_num_t   *l_distrib = (PDM_g_num_t   *) malloc (n_samples * sizeof(PDM_g_num_t));
  
  for (id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  /* hilbert_codes are supposed to be ordered */

  sample_code = sampling[bucket_id];

  for (i = 0; i < n_codes; i++) {

    PDM_g_num_t   o_id = order[i];

    if (sample_code >= hilbert_codes[o_id])
      l_distrib[bucket_id - 1] += weight[o_id];

    else {

      while (hilbert_codes[o_id] > sample_code) {
        bucket_id++;
        assert(bucket_id < n_samples + 1);
        sample_code = sampling[bucket_id];
      }

      l_distrib[bucket_id - 1] += weight[o_id];

    }

  } /* End of loop on elements */

  /* Define the global distribution */

  PDM_MPI_Allreduce(l_distrib, g_distrib, n_samples, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  free(l_distrib);

  /* Define the cumulative frequency related to g_distribution */

  cfreq[0] = 0.;
  for (id = 0; id < n_samples; id++) {
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif  
    double _g_distrib  = (double)g_distrib[id];
    double _gsum_weight = (double)gsum_weight;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
    cfreq[id+1] = cfreq[id] + _g_distrib/_gsum_weight;
  }
  cfreq[n_samples] = 1.0;

#if 0 && defined(DEBUG) && !defined(DEBUG) /* For debugging purpose only */

  if (cs_glob_rank_id <= 0) {

    FILE  *dbg_file = NULL;
    int  len;
    static int  loop_id1 = 0;

    len = strlen("DistribOutput_l.dat")+1+2;
    char  *rfilename = (char *) malloc (len * sizeof(char));
    sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1);

    loop_id1++;

    dbg_file = fopen(rfilename, "w");

    fprintf(dbg_file,
            "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |"
            "Global Distrib\n");
    for (i = 0; i < n_samples; i++)
      fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
              i, (double)i/(double)n_samples, cfreq[i],
              (double)(sampling[i]), distrib[i]);
    fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
            i, 1.0, 1.0, 1.0, 0);

    fclose(dbg_file);
    free(rfilename);

  }

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

    if (sum != gsum_weight)
      PDM_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
    exit(1);
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
_update_sampling(int                  dim,
                 int                  n_ranks,
                 double               c_freq[],
                 PDM_hilbert_code_t  *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  double  s_low, s_high;

  PDM_hilbert_code_t  *new_sampling = NULL, *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute new_sampling */

  new_sampling = ( PDM_hilbert_code_t  *) malloc (sizeof(PDM_hilbert_code_t) * (n_samples + 1));
   
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
    PDM_printf( " <_update_distrib> (rank: %d) delta: %g, target: %g,"
               " next_id: %d, f_low: %g, f_high: %g, s_low: %g, s_high: %g\n"
               "\t => new_sampling: %g\n",
               cs_glob_rank_id, delta, target_freq, next_id,
               f_low, f_high, s_low, s_high, new_sampling[i+1]);
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
 *   dim           <-- 2D or 3D
 *   n_ranks       <-- number of ranks
 *   gmax_level    <-- level on which Hilbert encoding is build
 *   n_codes       <-- local number of Hilbert ids
 *   hilbert_codes <-- local list of Hilbert ids to distribute
 *   weight        <-- weighting related to each code
 *   order         <-- ordering array
 *   sampling      <-> pointer to pointer to a sampling array
 *   comm          <-- mpi communicator
 *
 * returns:
 *   fit associated to the returned sampling array
 *----------------------------------------------------------------------------*/

static double
_bucket_sampling(int                       dim,
                 int                       n_ranks,
                 int                 n_codes,
                 const PDM_hilbert_code_t  hilbert_codes[],
                 const int           weight[],
                 const int           order[],
                 PDM_hilbert_code_t       *sampling[],
                 PDM_MPI_Comm                  comm)
{
  int  i, n_iters;
  int   j;
  double  fit, best_fit, optim;

  PDM_g_num_t   lsum_weight = 0, gsum_weight = 0;
  PDM_hilbert_code_t  *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute the global number of elements and the optimal number of elements
     on each rank */

  for (j = 0; j < n_codes; j++)
    lsum_weight += weight[j];

  PDM_MPI_Allreduce(&lsum_weight, &gsum_weight, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double _gsum_weight = (double)gsum_weight; 
  double _n_ranks = (double)n_ranks;
  optim = _gsum_weight / _n_ranks;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

  /* Define a naive sampling (uniform distribution) */

  for (i = 0; i < n_samples + 1; i++)
    _sampling[i] = i*unit;

  /* Define the distribution associated to the current sampling array */

  PDM_g_num_t   *distrib = (PDM_g_num_t   *) malloc (sizeof(PDM_g_num_t) * n_samples);
  double  *cfreq = (double *) malloc (sizeof(double) * (n_samples + 1));

  _define_rank_distrib(dim,
                       n_ranks,
                       gsum_weight,
                       n_codes,
                       hilbert_codes,
                       weight,
                       order,
                       _sampling,
                       cfreq,
                       distrib,
                       comm);

  /* Initialize best choice */

  fit = _evaluate_distribution(n_ranks, distrib, optim);
  best_fit = fit;

  PDM_hilbert_code_t  *best_sampling = (PDM_hilbert_code_t  *) malloc (sizeof(PDM_hilbert_code_t) * (n_samples + 1));

  for (i = 0; i < (n_samples + 1); i++)
    best_sampling[i] = _sampling[i];

  /* Loop to get a better sampling array */

  for (n_iters = 0;
       (   n_iters < pdm_hilbert_distrib_n_iter_max
        && fit > pdm_hilbert_distrib_tol);
       n_iters++)  {

    _update_sampling(dim, n_ranks, cfreq, &_sampling);

    /* Compute the new distribution associated to the new sampling */

    _define_rank_distrib(dim,
                         n_ranks,
                         gsum_weight,
                         n_codes,
                         hilbert_codes,
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
      for (i = 0; i < (n_samples + 1); i++)
        best_sampling[i] = _sampling[i];

    }

  } /* End of while */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0) {
    PDM_printf( "\n  <_bucket_sampling> n_iter: %d, opt: %g, best_fit: %g\n",
               n_iters, optim, best_fit);
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
PDM_hilbert_get_coord_extents_seq(int               dim,
                              size_t            n_coords,
                              const double  coords[],
                              double        g_extents[])
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
      else if (coords[i*dim + j] > g_extents[j + dim])
        g_extents[j + dim] = coords[i*dim + j];
    }
  }
}

void
PDM_hilbert_get_coord_extents_par(int               dim,
                              size_t            n_coords,
                              const double  coords[],
                              double        g_extents[],
                              PDM_MPI_Comm          comm)
{
  PDM_hilbert_get_coord_extents_seq(dim,
                                    n_coords,
                                    coords,
                                    g_extents);

  if (comm != PDM_MPI_COMM_NULL)
    _local_to_global_extents(dim, g_extents, comm);
}

/*----------------------------------------------------------------------------
 * Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   encode   <-- type of encode
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   h_code   --> array of corresponding Hilbert codes (size: n_coords)
 *----------------------------------------------------------------------------*/

void
PDM_hilbert_encode_coords(int                 dim,
			  PDM_hilbert_encode_t encode,
                          const double    extents[],
                          int           n_coords,
                          const double    coords[],
                          PDM_hilbert_code_t  h_code[])
{

  const unsigned *idata = NULL;
  const unsigned *istate = NULL;

  switch (encode) {

  case PDM_HILBERT_CS:

    if (dim == 2) {
      idata  = _cs_idata2d;
      istate = _cs_istate2d;
    }
    else if (dim == 3) {
      idata  = _cs_idata3d;
      istate = _cs_istate3d;
    }
    break;
    
  case PDM_HILBERT_GRIEBEL:

    if (dim == 2) {
      PDM_error(__FILE__, __LINE__, 0, "pdm_hilbert_encode_coords : No data for griebel 2D\n");
      exit(0);
    }
    else if (dim == 3) {
      idata  = _griebel_idata3d;
      istate = _griebel_istate3d;
    }
    break;

  case PDM_HILBERT_BADER:

    if (dim == 2) {
      idata  = _bader_idata2d;
      istate = _bader_istate2d;
    }
    else if (dim == 3) {
      idata  = _bader_idata3d;
      istate = _bader_istate3d;
    }
    break;

  case PDM_HILBERT_ZHANG:  /*!< Zhang */
    if (dim == 2) {
      idata  = _zhang_idata2d;
      istate = _zhang_istate2d;
    }
    else if (dim == 3) {
      idata  = _zhang_idata3d;
      istate = _zhang_istate3d;
    }
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "pdm_hilbert_encode_coords : Unknow encode type\n");
    exit(0);
   
  }
  
  int i, j, k;
  double s[3], d[3], n[3];

  int e_dim = 0;
  int dim_map[3] = {-1, -1, -1};

  double d_max = 0.0;
  const double epsilon = 1e-4;

  for (i = 0; i < dim; i++) {
    s[i] = extents[i];
    d[i] = extents[i+dim] - extents[i];
  }

  /* Check if box is not flat */

  for (i = 0; i < dim; i++) {
    d[i] = extents[i+dim] - extents[i];
    d_max = _MAX(d_max, d[i]);
  }

  for (i = 0; i < dim; i++) {
    if (d[i] >= d_max * epsilon) {
      dim_map[e_dim] = i;
      e_dim += 1;
    }
  }

  switch(dim) {

  case 3:
    {
      if (e_dim == 3) {
        for (i = 0; i < n_coords; i++) {
          for (j = 0; j < 3; j++)
            n[j] = (coords[i*3 + j] - s[j]) / d[j];
	  if (encode == PDM_HILBERT_CS)
	    h_code[i] = _hilbert_encode_3d(n, idata, istate);
	  else
	    h_code[i] = _hilbert_encode_3d_others(n, idata, istate);
        }
      }
      else if (e_dim == 2) {
        for (i = 0; i < n_coords; i++) {
          for (j = 0; j < 2; j++) {
            k = dim_map[j];
            n[j] = (coords[i*3 + k] - s[k]) / d[k];
          }
          h_code[i] = _hilbert_encode_2d(n, idata, istate);
        }
      }
      else if (e_dim == 1) {
        for (i = 0; i < n_coords; i++) {
          k = dim_map[0];
          n[0] = (coords[i*3 + k] - s[k]) / d[k];
          h_code[i] = _hilbert_encode_1d(n);
        }
      }
    }
    break;

  case 2:
    {
      if (e_dim == 2) {
        for (i = 0; i < n_coords; i++) {
          for (j = 0; j < 2; j++)
            n[j] = (coords[i*2 + j] - s[j]) / d[j];
          h_code[i] = _hilbert_encode_2d(n, idata, istate);
        }
      }
      else if (e_dim == 1) {
        for (i = 0; i < n_coords; i++) {
          k = dim_map[0];
          n[0] = (coords[i*3 + k] - s[k]) / d[k];
          h_code[i] = _hilbert_encode_1d(n);
        }
      }
    }
    break;

  case 1:
    {
      for (i = 0; i < n_coords; i++) {
        n[0] = (coords[i] - s[0]) / d[0];
        h_code[i] = _hilbert_encode_1d(n);
      }
    }
    break;

  default:
    assert(dim > 0 && dim < 4);
    break;
  }
}

/*----------------------------------------------------------------------------
 * Locally order a list of Hilbert ids.
 *
 * parameters:
 *   n_codes       <-- number of Hilbert ids to order
 *   hilbert_codes <-- array of Hilbert ids to order
 *   order         --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
PDM_hilbert_local_order(int                 n_codes,
                        const PDM_hilbert_code_t  hilbert_codes[],
                        int                 order[])
{
  int   i, tmp;

  assert(n_codes == 0 || hilbert_codes != NULL);
  assert(n_codes == 0 || order != NULL);

  for (i = 0; i < n_codes; i++)
    order[i] = i;

  /* Build heap */

  for (i = n_codes/2 - 1; (int)i >= 0; i--)
    _descend_hilbert_heap(i,  n_codes, hilbert_codes, order);

  /* Sort array */

  for (i = n_codes - 1; (int)i >= 0; i--) {

    tmp = order[0];
    order[0] = order[i];
    order[i] = tmp;

    _descend_hilbert_heap(0, i, hilbert_codes, order);

  }
}

/*----------------------------------------------------------------------------
 * Locally order a list of coordinates based on their Hilbert code.
 *
 * This variant may use a maximum depth of 32 levels, and switches
 * to lexicographical ordering if this is not enough.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   encode   <-- type of encode
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   order    --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
PDM_hilbert_local_order_coords(int                dim,
                               PDM_hilbert_encode_t encode,
                               const double   extents[],
                               int          n_coords,
                               const double   coords[],
                               int          order[])
{
  PDM_hilbert_code_t *h_code = (PDM_hilbert_code_t *) malloc (sizeof(PDM_hilbert_code_t) * n_coords);

  PDM_hilbert_encode_coords(dim, encode, extents, n_coords, coords, h_code);

  PDM_hilbert_local_order(n_coords, h_code, order);

  free(h_code);
}

/*----------------------------------------------------------------------------
 * Get the quantile associated to a Hilbert code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * parameters:
 *   n_quantiles    <-- number of quantiles
 *   code           <-- code we are searching for
 *   quantile_start <-- first Hilbert code in each quantile (size: n_quantiles)
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

size_t
PDM_hilbert_quantile_search(size_t              n_quantiles,
                            PDM_hilbert_code_t  code,
                            PDM_hilbert_code_t  quantile_start[])
{
  size_t mid_id = 0;
  size_t start_id = 0;
  size_t end_id = n_quantiles;

  /* use binary search */

  while (start_id + 1 < end_id) {
    mid_id = start_id + ((end_id -start_id) / 2);
    if (quantile_start[mid_id] > code)
      end_id = mid_id;
    else
      start_id = mid_id;
  }

  /* We may have stopped short of the required value,
     or have multiple occurences of a quantile start
     (in case of empty quantiles), of which we want to
     find the find highest one */

  while (   start_id < n_quantiles - 1
         && code >= quantile_start[start_id+1])
    start_id++;

  return start_id;
}

/*----------------------------------------------------------------------------
 * Build a global Hilbert encoding rank index.
 *
 * The rank_index[i] contains the first Hilbert code assigned to rank [i].
 *
 * parameters:
 *   dim          <-- 1D, 2D or 3D
 *   n_codes      <-- number of Hilbert codes to be indexed
 *   hilbert_code <-- array of Hilbert codes to be indexed
 *   weight       <-- weighting related to each code
 *   order        <-- ordering array
 *   rank_index   <-> pointer to the global Hilbert encoding rank index
 *   comm         <-- MPI communicator on which we build the global index
 *
 * returns:
 *  the fit related to the Hilbert encoding distribution (lower is better).
 *----------------------------------------------------------------------------*/

double
PDM_hilbert_build_rank_index(int                       dim,
                             int                       n_t_part,
                             int                       n_codes,
                             const PDM_hilbert_code_t  hilbert_code[],
                             const int                 weight[],
                             const int                 order[],
                             PDM_hilbert_code_t        rank_index[],
                             PDM_MPI_Comm                  comm)
{
  int  i, id, rank_id, n_samples;
  double  best_fit;

  const int  sampling_factor = _sampling_factors[dim];

  /* Allocations and Initialization */

  n_samples = sampling_factor * n_t_part;

  PDM_hilbert_code_t  *sampling = 
          (PDM_hilbert_code_t  *) malloc(sizeof(PDM_hilbert_code_t) * (n_samples + 1));

  for (i = 0; i < (n_samples + 1); i++)
    sampling[i] = 0;

  best_fit = _bucket_sampling(dim,
                              n_t_part,
                              n_codes,
                              hilbert_code,
                              weight,
                              order,
                              &sampling,
                              comm);

  /* Define Hilbert index */

  for (rank_id = 0; rank_id < n_t_part + 1; rank_id++) {
    id = rank_id * sampling_factor;
    rank_index[rank_id] = sampling[id];
  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  { /* Dump Hilbert index and associated sampling on rank 0 */
    PDM_printf( "\nHilbert rank index:\n\n");
    for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
      id = sampling_factor * rank_id;
      PDM_printf( "rank: %5d (sampling:   %f)\n"
                 "           rank_index: %f\n",
                 rank_id,
                 (double)sampling[id], (double)rank_index[rank_id]);
    }
    PDM_printf( "\n");
    fflush(stdout);
  }
#endif

  /* Free memory */

  free(sampling);

  return best_fit;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#undef _MIN
#undef _MAX
