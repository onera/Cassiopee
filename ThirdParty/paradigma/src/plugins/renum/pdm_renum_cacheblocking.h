#ifndef __PDM_RENUM_CACHE_BLOCKING_H__
#define __PDM_RENUM_CACHE_BLOCKING_H__

/*============================================================================
 * Hilbert space-filling curve construction for coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

#include "pdm_part.h"
#include "pdm_part_priv.h"

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


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Add renumbering cacheblocking method in the ppart methods
 *
 */

void 
PDM_renum_cacheblocking_ppart_add
(
void
);


/**
 * Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition 
 * \param [in]   split_method        Split method
 * \param [in]   nCellPerCacheWanted Approximate number of cells on each cache
 * \param [in]   isAsynchrone        [0 : Synchrone/1 : Asynchrone]
 * \param [in]   isVectorisation     [0 : No/1 : Yes]
 *
 */

void 
PDM_renum_cacheblocking
(
 _part_t     *part,
int          split_method, 
int          nCellPerCacheWanted, 
int          isAsynchrone,  
int          isVectorisation 
);

/**
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition 
 *
 */
void 
PDM_renum_cacheblocking_compute_loop_array
(
 _part_t     *part
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_RENUM_CACHE_BLOCKING_H__ */
