#ifndef __PDM_COARSE_MESH_ANISO_AGGLO_H__
#define __PDM_COARSE_MESH_ANISO_AGGLO_H__

/*----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"
#include "pdm_part_coarse_mesh.h"


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

/**
 * \brief Add aniso_agglo coarse mesh method
 *
 */

void 
PDM_coarse_mesh_aniso_agglo_add
(
void
);

/**
 *
 * \brief Add isotropic array to current coarse mesh 
 * 
 * \param [in]   cmId                      Coarse mesh identifier
 * \param [in]   iPart                     Current partition
 * 
 * \param [out]  agglomerationLines
 * \param [out]  agglomerationLinesIdx
 * \param [out]  isOnFineBnd            
 *
 */

void 
PDM_part_coarse_mesh_part_set_anisotropic_info
(
 const int    cmId,
 const int    iPart,       
 const int    *agglomerationLinesInit,
 const int    *agglomerationLinesInitIdx,
 const int      agglomerationLinesInitIdx_size,
 const int    *isOnFineBndInit
);

void
PROCF (pdm_part_coarse_mesh_part_set_anisotropic_info, PDM_PART_COARSE_MESH_PART_SET_ANISOTROPIC_INFO)
(
 int          *cmId,
 int          *iPart,
 int          *agglomerationLinesInit,
 int          *agglomerationLinesInitIdx,
 int          *agglomerationLinesInitIdx_size,
 int          *isOnFineBndInit
);


/**
 *
 * \brief Add option for anisotropic mesh agglomeration
 *
 * \param [out]  cmId              Coarse mesh identifier
 * \param [in]   Option
 */
    
void 
PDM_part_coarse_mesh_add_option_anisotropic
(
 int        cmId,
 const int* anisotropicOption
);

void
PROCF (pdm_part_coarse_mesh_add_option_anisotropic, PDM_PART_COARSE_MESH_ADD_OPTION_ANISOTROPIC)
(
 int        *cmId,
 const int  *anisotropicOption
 );


/**
 *
 * \brief Return a mesh partition
 * 
 * \param [in]   cmId                      Coarse mesh identifier
 * \param [in]   iPart                     Current partition
 * 
 * \param [out]  agglomerationLines
 * \param [out]  agglomerationLinesIdx
 * \param [out]  isOnFineBnd            
 *
 */

void 
PDM_part_coarse_mesh_part_get_anisotropic_info
(
 const int    cmId,
 const int    iPart,       
 int          **agglomerationLines,
 int          **agglomerationLinesIdx,
 int           *agglomerationLinesIdx_size,
 int          **isOnFineBnd       
 );


void
PROCF (pdm_part_coarse_mesh_part_get_anisotropic_info, PDM_PART_COARSE_MESH_PART_GET_ANISOTROPIC_INFO)
(
 int          *cmId,
 int          *iPart,       
 int          *agglomerationLines,
 int          *agglomerationLinesIdx,
 int          *agglomerationLinesIdx_size,
 int          *isOnFineBnd
 );

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_COARSE_MESH_ANISO_AGGLO_H__ */
