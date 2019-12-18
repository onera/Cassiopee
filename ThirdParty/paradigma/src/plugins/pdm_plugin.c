
/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_plugin.h"
#include "pdm_renum_cacheblocking.h"
#include "pdm_coarse_mesh_aniso_agglo.h"

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
 * \brief Load plugins
 */

void
PDM_plugin_load
(
 void
)
{
 PDM_renum_cacheblocking_ppart_add();
 PDM_renum_cacheblocking2_ppart_add();
#ifdef PDM_HAVE_ANISO_AGGLO
 PDM_coarse_mesh_aniso_agglo_add();
#endif
}

void
PROCF (pdm_plugin_load, PDM_PLUGIN_LOAD)
(
 void
)
{
 PDM_plugin_load();
}

/*----------------------------------------------------------------------------
 * Load plugins
 *---------------------------------------------------------------------------*/

/**
 * \brief Free plugins
 */

void
PDM_plugin_free
(
 void
)
{
    printf("Free plugin \n");
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
