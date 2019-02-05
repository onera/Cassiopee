#ifndef __PDM_MORTON_H__
#define __PDM_MORTON_H__

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

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef enum {

  PDM_MORTON_EQUAL_ID,
  PDM_MORTON_SAME_ANCHOR,
  PDM_MORTON_DIFFERENT_ID

} PDM_morton_compare_t;

typedef unsigned int   PDM_morton_int_t;

typedef struct {

  PDM_morton_int_t   L;     /* Level in the tree structure */
  PDM_morton_int_t   X[3];  /* X, Y, Z coordinates in Cartesian grid */

} PDM_morton_code_t;

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
                             const double      coords[],
                             double             g_extents[],
                             PDM_MPI_Comm          comm);

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
                              PDM_MPI_Comm          comm);

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
PDM_morton_encode(int                dim,
                  PDM_morton_int_t   level,
                  const double   coords[]);

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
                         PDM_morton_code_t  m_code[]);

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
                        PDM_morton_code_t  children[]);

/*----------------------------------------------------------------------------
 * Compare two Morton encoding and check if these two codes are equal,
 * different or shared the same anchor.
 *
 * parameters:
 *   dim    <-- 2D or 3D
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  a type on the kind of relation between the two Morton encodings.
 *----------------------------------------------------------------------------*/

PDM_morton_compare_t
PDM_morton_compare(int                dim,
                   PDM_morton_code_t  code_a,
                   PDM_morton_code_t  code_b);

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
                       int                order[]);

/*----------------------------------------------------------------------------
 * Locally sort a list of Morton ids.
 *
 * parameters:
 *   n_codes      <-- number of Morton ids to order
 *   morton_codes <-> array of Morton ids to sort
 *----------------------------------------------------------------------------*/

void
PDM_morton_local_sort(int          n_codes,
                      PDM_morton_code_t  morton_codes[]);

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
                  PDM_morton_code_t  b);

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
                  PDM_morton_code_t  b);

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
                         PDM_morton_code_t  *codes);

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
                           PDM_morton_code_t  *quantile_start);

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
                            PDM_MPI_Comm                 comm);

/*----------------------------------------------------------------------------
 * Dump a Morton to standard output or to a file.
 *
 * parameters:
 *   dim  <-- 2D or 3D
 *   code <-- Morton code to dump
 *----------------------------------------------------------------------------*/

void
PDM_morton_dump(int                 dim,
                PDM_morton_code_t   code);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MORTON_H__ */
