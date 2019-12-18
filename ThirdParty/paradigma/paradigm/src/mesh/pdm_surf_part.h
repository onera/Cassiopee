#ifndef __PDM_SURF_PART_H__
#define __PDM_SURF_PART_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \struct PDM_surf_part_t
 * \brief  Surface partition
 *
 *  PDM_surf_part_t defines a surface partition
 *
 */

typedef struct _pdm_surf_part_t PDM_surf_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Return an intialized \ref PDM_surf_part_t structure
 *
 * This function returns an initialized \ref PDM_surf_part_t structure
 *
 * \param [in]  nFace       Number of faces
 * \param [in]  faceVtxIdx  Index in the face -> vertex connectivity
 * \param [in]  faceVtxIdx  face -> vertex connectivity
 * \param [in]  faceLnToGn  Local face numbering to global face numbering
 * \param [in]  nVtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtxLnToGn   Local vertex numbering to global vertex numbering
 *
 * \return      A new initialized \ref _part_t structure
 *
 */

PDM_surf_part_t *
PDM_surf_part_create
(
const int         nFace,
const int        *faceVtxIdx,
const int        *faceVtx,
const PDM_g_num_t *faceLnToGn,
const int         nVtx,
const double     *coords,
const PDM_g_num_t *vtxLnToGn
 );


/**
 * \brief Delete a \ref PDM_surf_part_t structure
 *
 * This function deletes a  PDM_surf_part_t structure
 *
 * \param [in]  part      part to delete
 *
 * \return     Null pointer
 */

PDM_surf_part_t *
PDM_surf_part_free
(
 PDM_surf_part_t * part
);


/**
 * \brief Compute partition edge entities
 *
 * This function defines edges of an initial partition and
 * computes edge connectivities
 *
 * \param [in]  part      Partition to compute
 *
 */

void
PDM_surf_part_build_edges
(
PDM_surf_part_t *part
);


/**
 * \brief Return faceLnToGn
 *
 *
 * \param [in]  part      Partition to compute
 *
 */

const PDM_g_num_t *
PDM_surf_part_faceLnToGn_get
(
PDM_surf_part_t *part
);




#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_PART_H__ */
