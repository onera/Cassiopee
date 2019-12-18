#ifndef __PDM_MULTICOLORING_H__
#define __PDM_MULTICOLORING_H__

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
void
PDM_compute_distributed_multicoloring
(
 _dist_csr* dcsr
);

int
compute_multicoloring
(
 int          sizeG,
 int*         ia,
 PDM_g_num_t* ja,
 int*         color,
 int*         saveColor
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTICOLORING_H__ */
