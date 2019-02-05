#ifndef __PDM_CUTHILL_H__
#define __PDM_CUTHILL_H__

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

/**
 * \enum PDM_rcm_graph_t
 * \brief Type of graph for building RCM numbering 
 *
 */

typedef enum {

  PDM_CUTHILL_GR1 = 0,  /*!< Rank 1 Graph */
  PDM_CUTHILL_GR2 = 1,  /*!< Rank 2 Graph */

} PDM_cuthill_graph_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine the reverse Cut-Hill Mac-Kee numbering
 *
 * parameters:
 *   part       --> Mesh Partition 
 *   perm       <-- Array of prmutation 
 *---------------------------------------------------------------------------*/

void 
PDM_cuthill_generate
(
 _part_t           *ppart,
 int               *perm
);

/*----------------------------------------------------------------------------
 * Compute the bandwidth of the current mesh 
 *
 * parameters:
 *   part       --> Mesh Partition 
 *---------------------------------------------------------------------------*/
int 
PDM_cuthill_checkbandwidth
(
 _part_t           *ppart
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_CUTHILL_H__ */
