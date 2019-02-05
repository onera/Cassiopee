#ifndef __PDM_HILBERT_H__
#define __PDM_HILBERT_H__

/*============================================================================
 * Hilbert space-filling curve construction for coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
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

/* Hilbert code
   (could be switched from double to long double for extened range) */

typedef double  PDM_hilbert_code_t;


/**
 * \enum PDM_hilbert_encode_t
 * \brief Type of encode
 *
 */

typedef enum {

  PDM_HILBERT_CS      = 0,  /*!< Code Saturne */
  PDM_HILBERT_GRIEBEL = 1,  /*!< Griebel */
  PDM_HILBERT_BADER   = 2,  /*!< Bader */
  PDM_HILBERT_ZHANG   = 3,  /*!< Zhang */

} PDM_hilbert_encode_t;


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
PDM_hilbert_get_coord_extents_par(int           dim,
                                  size_t        n_coords,
                                  const double  coords[],
                                  double        g_extents[],
                                  PDM_MPI_Comm      comm);
void
PDM_hilbert_get_coord_extents_seq(int      dim,
                              size_t       n_coords,
                              const double coords[],
                              double       g_extents[]);

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
                          PDM_hilbert_code_t  h_code[]);

/*----------------------------------------------------------------------------
 * Locally order a list of Hilbert ids.
 *
 * This variant uses an encoding into floating-point numbers. In 3D, this
 * limits us to 19 levels for a double, though using a long double could
 * increase this range.
 *
 * parameters:
 *   n_codes       <-- number of Hilbert ids to order
 *   hilbert_codes <-- array of Hilbert ids to order
 *   order         --> pointer to pre-allocated ordering table
 *----------------------------------------------------------------------------*/

void
PDM_hilbert_local_order(int                 n_codes,
                        const PDM_hilbert_code_t  hilbert_codes[],
                        int                 order[]);

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
                               int          order[]);

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
                            PDM_hilbert_code_t  quantile_start[]);

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
                             int                       nPart,
                             int                       n_codes,
                             const PDM_hilbert_code_t  hilbert_code[],
                             const int                 weight[],
                             const int                 order[],
                             PDM_hilbert_code_t        rank_index[],
                             PDM_MPI_Comm                  comm);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HILBERT_H__ */
