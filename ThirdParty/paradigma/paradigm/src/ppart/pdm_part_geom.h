#ifndef __PDM_PART_GEOM_H__
#define __PDM_PART_GEOM_H__

/*============================================================================
 * Mesh partitioning with geometric methods (SFC)
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"
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

/**
 * \enum PDM_split_geom_t
 * \brief Geometric plit method
 *
 */

typedef enum {
  PDM_PART_GEOM_HILBERT = 1,
} PDM_part_geom_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Perform geometric partitioning
 *
 * \param [in]   method         Geometric method
 * \param [in]   nPart          Number of partition to build on this process
 * \param [in]   comm           Communicator
 * \param [in]   dNCell         Number of distributed cells
 * \param [in]   dNFace         Number of distributed faces
 * \param [in]   dNVtx          Number of distributed vertices
 * \param [in]   dCellFaceIdx   Distributed cell face connectivity index or NULL
 *                              (size : dNCell + 1, numbering : 0 to n-1)
 * \param [in]   dCellFace      Distributed cell face connectivity or NULL
 *                              (size : dFaceVtxIdx[dNCell], numbering : 1 to n)
 * \param [in]   dCellWeight    Cell weight (size : nCell) or NULL
 * \param [in]   dFaceVtxIdx    Distributed face to vertex connectivity index 
 *                              (size : dNFace + 1, numbering : 0 to n-1)
 * \param [in]   dFaceVtx       Distributed face to vertex connectivity 
 *                              (size : dFaceVtxIdx[dNFace], numbering : 1 to n)
 * \param [in]   dVtxCoord      Distributed vertex coordinates 
 *                              (size : 3*dNVtx)
 * \param [inout]   dCellPart      Distributed cell partitioning 
 *                              (size = dNCell)
 *
 */

void
PDM_part_geom
(
 PDM_part_geom_t     method,
 const int           nPart,       
 const PDM_MPI_Comm      comm,       
 const int           dNCell,
 const int          *dCellFaceIdx,
 const PDM_g_num_t *dCellFace,
 const int          *dCellWeight,
 const int          *dFaceVtxIdx,
 const PDM_g_num_t *dFaceVtx,
 const PDM_g_num_t *dFaceProc,       
 const double       *dVtxCoord,
 const PDM_g_num_t *dVtxProc,       
 int                *dcellPart
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_GEOM_H__ */
