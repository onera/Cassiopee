#ifndef __PDM_PART_GRAPH_H__
#define __PDM_PART_GRAPH_H__

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

/*----------------------------------------------------------------------------
 * Perform a cells renumbering with cache blocking (Synchrone)
 *
 * parameters:
 *   part       --> Mesh Partition 
 *---------------------------------------------------------------------------*/

void 
PDM_part_graph_split
(
 int         method,
 int         nPart,
 _part_t    *part_ini,
 int        *cellCellIdx,
 int        *cellCell,
 int        *cellWeight,
 int        *faceWeight,
 int       **cellPart
);

/*----------------------------------------------------------------------------
 * Perform a cells renumbering with cache blocking (Asynchrone)
 *
 * parameters:
 *   part       --> Mesh Partition 
 *---------------------------------------------------------------------------*/

void 
PDM_part_graph_compute_from_face_cell
(
  _part_t        *part_ini,
  int           **cellCellIdxCompressed,
  int           **cellCellCompressed
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_RENUM_CACHE_BLOCKING_H__ */
