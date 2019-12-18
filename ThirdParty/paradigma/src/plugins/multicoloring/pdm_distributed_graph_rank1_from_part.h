#ifndef __PDM_PDM_DISTRIBUTED_GRAPH_RANK1_FROM_PART_H__
#define __PDM_PDM_DISTRIBUTED_GRAPH_RANK1_FROM_PART_H__

/*----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"
#include "pdm_distributed_csr.h"


#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

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
_dist_csr*
PDM_compute_dist_graph_rank1_from_ppart
(
 int ppartId,
 int nPart
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTICOLORING_H__ */
