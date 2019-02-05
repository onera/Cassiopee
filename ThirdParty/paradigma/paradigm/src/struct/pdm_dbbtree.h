#ifndef __PDM_DBBTREE_H__
#define __PDM_DBBTREE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_box.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_dbbtree_t
 * \brief  Distributed boundary box tree
 * 
 *  PDM_dbbtree_t defines a distributed boundary box tree
 *
 */

typedef struct _PDM_dbbtree_t PDM_dbbtree_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_dbbtree_t structure
 *
 * This function returns an initialized \ref PDM_dbbtree_t structure
 *
 * \param [in]  comm  Associated communicator
 * \param [in]  dim   boxes dimension
 *
 * \return      A new initialized \ref PDM_dbbtree_t structure
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_create
(
 PDM_MPI_Comm          comm,
 int               dim
 );



/**
 * \brief Free a \ref PDM_dbbtree_t structure
 *
 * \param [in]  dbbt   Pointer to a distributed bounding box tree
 *
 * \return      NULL
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_free
(
PDM_dbbtree_t     *dbbt
);


/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  nPart    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to 
 * the tree location
 *
 */

PDM_box_set_t  *
PDM_dbbtree_boxes_set
(
PDM_dbbtree_t     *dbbt,
const int          nPart,
const int         *nElts,
const double     **extents,
const PDM_g_num_t **gNum
);


/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function  assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  nPart    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 * \param [out] box_index Pointer to the index array on associated tree bounding boxeq
 * \param [out] box_g_num Pointer to the list of intersecting bounding boxes
 *
 * \return associated \ref PDM_box_set_t structure distributed according 
 * to the tree intersection
 *
 */


PDM_box_set_t  *
PDM_dbbtree_intersect_boxes_set
(
PDM_dbbtree_t    *dbbt,
const int         nPart,
const int        *nElts,
const double     **extents,
const PDM_g_num_t **gNum,
int              *box_index[],
int              *box_l_num[]
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DBBTREE_H__ */
