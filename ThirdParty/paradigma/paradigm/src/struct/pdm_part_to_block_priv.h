#ifndef __PDM_writer_PART_TO_BLOCK_PRIV_H__
#define __PDM_writer_PART_TO_BLOCK_PRIV_H__

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
 * \struct _cs_part_to_block_t
 *
 * \brief  Data transfer from partitions to blocks
 * 
 */

typedef struct {

  /*
   * Block distribution properties
   */

  PDM_part_to_block_distrib_t  t_distrib;        /*!< Distribution type */
  PDM_part_to_block_post_t     t_post;           /*!< post processing type */
  int                         n_activeRanks;    /*!< Number of active ranks */
  int                        *activeRanks;      /*!< List of active ranks */
  PDM_MPI_Comm                    comm;             /*!< MSG communicator */
  int                         s_comm;           /*!< Communicator size */
  int                         myRank;           /*!< Current rank in comm */
  int                         isMyRankActive;   /*!< Is active current rank */
  float                       partActiveNode;   /*!< Part of active nodes */

  /*
   * General exchange data
   */

  int                         n_part;            /*!< Number of parts */
  int                        *n_elt;             /*!< Number of elements for any part */
  int                         n_eltProc;         /*!< Number of elements on the current processus */
  PDM_g_num_t                 **gnum_elt;          /*!< Global numbering of elements for any part */
  int                        *destProc;          /*!< Destination process for any element (size = n_eltProc) */
  int                        *idxInSendData;     /*!< Index in send data for any element (size = n_eltProc) */

  PDM_g_num_t                  *dataDistribIndex;  /*!< Data distribution on ranks 
                                                  (size = s_comm + 1) */
  int                         s_blockMin;        /*!< Minimum block size */
  int                         s_blockMax;        /*!< Maximum block size */

  int                        *i_sendData;        /*!< Data to send to other processes index
                                                   (size = s_comm) */
  int                        *i_recvData;        /*!< Received Data from other processes index
                                                   (size = s_comm) */
  int                        *n_sendData;        /*!< Number of data to send to other processes
                                                   (size = s_comm) */
  int                        *n_recvData;        /*!< Number of received Data from other processes
                                                   (size = s_comm) */

  int                         tn_sendData;       /*!< Total number of sended data */ 
  int                         tn_recvData;       /*!< Total number of received data */ 
  PDM_g_num_t                  *sorted_recvGnum;   /*!< Sorted Global number of 
                                                      reveived data (size = tn_recvData) */
  int                        *order;             /*!< Order of sorted_recvGnum
                                                   (size = tn_recvData) */
  int                         n_eltBlock ;      /*!< Number of element in current block */
  PDM_g_num_t                  *block_gnum;        /*!< Sorted Global number of 
                                                      reveived data (size = block_n_elt) */


} _cs_part_to_block_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_writer_PART_TO_BLOCK_PRIV_H__ */
