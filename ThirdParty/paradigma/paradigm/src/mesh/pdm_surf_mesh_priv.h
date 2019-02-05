#ifndef __PDM_SURF_MESH_PRIV_H__
#define __PDM_SURF_MESH_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_surf_part.h"
#include "pdm_graph_bound.h"

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
 * \struct pdm_surf_mesh_t
 * \brief  Mesh type
 * 
 * \ref pdm_surf_mesh_t defines a mesh structure
 *
 */

struct _pdm_surf_mesh_t{

  PDM_MPI_Comm    comm;   /*!< MPI communicator of mesh */
  PDM_g_num_t  nGFace; /*!< Global number of faces     */
  PDM_g_num_t  nGVtx;  /*!< Global number of vertices  */
  PDM_g_num_t  nGEdge;  /*!< Global number of edges  */
  int         nGPart; /*!< Number of global partitions */
  int         nPart;  /*!< Number of local partitions */
  double      gMinCarLgthVtx; /*!< Global min 
                                   of caracteristic length vertex */
  double      gMaxCarLgthVtx; /*!< Global max 
                                   of caracteristic length vertex */
  PDM_surf_part_t   **part;   /*!< Mesh partition             */ 

  PDM_graph_bound_t *interPartEdgeGraph; /*!< Inter partition edges graph */
  PDM_graph_bound_t *interPartVtxGraph;  /*!< Inter partition vertices graph */

  PDM_part_bound_t **vtxPartBound; /*!< pointer on bounding vertices of each part */
  PDM_part_bound_t **edgePartBound; /*!< pointer on bounding edges of each part */

} ;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_PRIV_H__ */
