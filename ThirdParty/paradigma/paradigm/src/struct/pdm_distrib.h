#ifndef PDM_DISTRIB_H
#define PDM_DISTRIB_H

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


/*============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_distrib_compute
(
 const int           dnelt,
       PDM_g_num_t  *elt_distrib,
       int           offset,
 const PDM_MPI_Comm  comm
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // PDM_DISTRIB_H
