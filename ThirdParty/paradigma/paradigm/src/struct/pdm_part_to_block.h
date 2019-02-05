#ifndef __PDM_writer_PART_TO_BLOCK_H__
#define __PDM_writer_PART_TO_BLOCK_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_writer.h"
#include "pdm.h"
#include "pdm_mpi.h"

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
 * \enum PDM_part_to_block_distrib_t
 * \brief Type of distribution
 *
 */

typedef enum {

  PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC          = 0,  /*!< Distribute block on all processes */
  PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE = 1,  /*!< Distribute block on one processe pere node */
  PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE      = 2,  /*!< Distribute block on part of nodes */

} PDM_part_to_block_distrib_t;


/**
 * \enum PDM_part_to_block_post_t
 * \brief Type of Post processing about blocks
 *
 */

typedef enum {

  PDM_PART_TO_BLOCK_POST_NOTHING       = 0,  /*!< No post processing     */
  PDM_PART_TO_BLOCk_POST_CLEANUP       = 1,  /*!< Cleanup multi-elements */
  PDM_PART_TO_BLOCK_POST_MERGE         = 2,  /*!< Merge multi-elements   */

} PDM_part_to_block_post_t;

/**
 * \enum PDM_part_to_block_distrib_t
 * \brief Type of distribution
 *
 */

//typedef enum {
//
//  PDM_PART_TO_BLOCK_STRIDE_CST = 0,  /*!< Constant stride element */
//  PDM_PART_TO_BLOCK_STRIDE_VAR = 1,  /*!< Variable stride element */               
//
//} PDM_part_to_block_stride_t;


/**
 * \struct PDM_part_to_block_t
 * \brief  Partitioning to block redistribution
 * 
 */

typedef struct _cs_part_to_block_t PDM_part_to_block_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition      
 * \param [in]   comm            MPI communicator         
 *
 * \return   Initialized cs_part_to_block         
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 float                         partActiveNode,
 PDM_g_num_t                 **gnum_elt,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);


/**
 *
 * \brief Return number of active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of active ranks
 *
 */

int 
PDM_part_to_block_n_active_ranks_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  active ranks
 *
 */

int *
PDM_part_to_block_active_ranks_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return if current rank is active
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  if current rank is active
 *
 */

int
PDM_part_to_block_is_active_rank
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return number of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of element in the current process
 *
 */

int 
PDM_part_to_block_n_elt_block_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return global numbers of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers
 *
 */

PDM_g_num_t *
PDM_part_to_block_block_gnum_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 * \param [in]   part_stride  Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   part_data    partitioned data
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */

int 
PDM_part_to_block_exch
(
 PDM_part_to_block_t       *ptb,
 size_t                     s_data,
 PDM_stride_t               t_stride,
 int                        cst_stride,
 int                      **part_stride,
 void                     **part_data,
 int                      **block_stride,
 void                     **block_data
);


/**
 *
 * \brief Free a part to block structure                       
 *
 * \param [inout] ptb         Part to block structure
 *
 * \return       NULL                             
 */

PDM_part_to_block_t *
PDM_part_to_block_free
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return block distribution                     
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Distribution (size = communicator size + 1)                            
 */

PDM_g_num_t *
PDM_part_to_block_distrib_index_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return processus destination                     
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Destination (size = sum of partition elements)                            
 */

PDM_l_num_t *
PDM_part_to_block_destination_get
(
 PDM_part_to_block_t *ptb
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_writer_PART_TO_BLOCK_H__ */
