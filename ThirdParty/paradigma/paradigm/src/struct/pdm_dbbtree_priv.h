#ifndef __PDM_DBBTREE_PRIV_H__
#define __PDM_DBBTREE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_writer.h"
#include "pdm_box_tree.h"
#include "pdm_box_tree.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _box_tree_stats_t
 * \brief  Statistic about bbtre
 * 
 *  _box_tree_stats_t defines statistics about bounding box tree
 *
 */

typedef struct {

  int         dim;                     /*!< Layout dimension */

  /* The following fields have 3 global values:
     mean on ranks, minimum on ranks, and maximum on ranks */

  int         depth[3];                /*!< Tree depth */
  int         n_leaves[3];             /*!< Number of leaves */
  int         n_boxes[3];              /*!< Number of associated boxes */
  int         n_threshold_leaves[3];   /*!< Number of leaves over threshold */
  int         n_leaf_boxes[3];         /*!< Number of boxes per leaf */
  size_t      mem_used[3];             /*!< Memory used */
  size_t      mem_required[3];         /*!< Memory temporarily required */

} _box_tree_stats_t;

/**
 * \struct PDM_ol_t
 * \brief  Overlay type
 * 
 * PDM_ol_t defines a overlaying structure 
 *
 */

typedef struct {

  int     maxTreeDepth;    /*!< Max tree depth for local BBTree */ 

  float   maxBoxRatio;     /*!< Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes) 
                             for local BBTree */

  int     maxBoxesLeafShared; /*!< Max number of boxes in a leaf for local BBTree */

  int     maxTreeDepthShared; /*!< Max tree depth for local BBTree */ 

  float   maxBoxRatioShared;  /*!< Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes) 
                               for local BBTree */

  int     maxBoxesLeaf;       /*!< Max number of boxes in a leaf for local BBTree */

  int     maxTreeDepthCoarse; /*!< Max tree depth for coarse shared BBTree */

  float   maxBoxRatioCoarse;  /*!< Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes) 
                                for coarse shared BBTree */

  int     maxBoxesLeafCoarse; /*!<  Max number of boxes in a leaf for coarse shared BBTree */

  PDM_box_set_t  *rankBoxes;  /*!< Rank Boxes */
  PDM_box_tree_t *btShared;   /*!< Shared Boundary box tree */
  _box_tree_stats_t btsShared;/*!< Shared Boundary box tree statistic */

  PDM_box_set_t  *boxes;      /*!< Boxes */
  PDM_box_tree_t *btLoc;      /*!< Local Boundary box tree */
  _box_tree_stats_t btsLoc;   /*!< Local Boundary box tree statistic */

  _box_tree_stats_t btsCoarse; /*!< Local Boundary box tree statistic */

  PDM_MPI_Comm   comm;             /*!< MPI communicator */
  int        dim;              /*!< Dimension */

} _PDM_dbbtree_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Static function definitions
 *============================================================================*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DBBTREE_PRIV_H__ */
