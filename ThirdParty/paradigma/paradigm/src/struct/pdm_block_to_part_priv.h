/*
 * File:   pdm_block_to_part_priv.h
 * Author: equemera
 *
 * Created on April 14, 2016, 8:16 AM
 */

#ifndef PDM_BLOCK_TO_PART_PRIV_H
#define	PDM_BLOCK_TO_PART_PRIV_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_block_to_part_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \struct _PDM_block_to_part_t
 *
 * \brief  Data transfer from blocks to partitions
 *
 */

typedef struct {
  PDM_g_num_t    *blockDistribIdx;  /*!< Block distribution
                                    * (size : \ref size of \ref comm + 1) */
  PDM_MPI_Comm       comm;             /*!< MSG communicator */
  int            s_comm;           /*!< Communicator size */
  int            myRank;           /*!< Current rank in comm */
  int            n_part;           /*!< Number of partitions */
  int           *n_elt;            /*!< Number of elements for each partition */
  int          **ind;              /*!< Ind for each element partition in distributed_data */
  int           *requested_data_n; /*!< Numer of requested data for each process index
                                    * (size : s_comm) */
  int           *requested_data_idx;/*!< Requested data for each process index
                                    * (size : s_comm) */
  int           *distributed_data_n;/*!< Numer of distributed data for each process index
                                    * (size : s_comm) */
  int           *distributed_data_idx;/*!< Distributed data for each process index
                                    * (size : s_comm) */
  int           *distributed_data;    /*!< Distributed data for each process
                                    * (size : requestd_data_idx[s_comm - 1]
                                             +requestd_data_n[s_comm - 1] ) */
  int           pttopt_comm;         /*!< Use point to point communication if pttopt_comm == 1 */

} _pdm_block_to_part_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/




#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_PART_PRIV_H */
