#ifndef __PDM_GRAPH_BOUND_PRIV_H__
#define __PDM_GRAPH_BOUND_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_bound.h"

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

static const int tag = 'P' + 'T' + 'K';

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _graph_bound_t
 * \brief  Communication graph between partitions
 *
 * _graph_bound_t defines a communication graph between partitions
 *
 */

struct _graph_bound_t {

  int              lComm;         /*!< Size of communicator */
  PDM_MPI_Comm          comm;         /*!< MPI communicator */
  int               nPart;        /*!< Number of local partitions */
  int               nExchRank;    /*!< Number of exchange ranks */
  int              *exchRank;     /*!< List of exchange ranks */
  PDM_part_bound_t **partBound;   /*!< pat_bound for each partition */

  int               nComp;          /*!< Number of exchanged data */
  int               lBuffer;        /*!< Size of buffer */
  void             *sendBuffer;     /*!< Sending Data buffer*/
  void             *recvBuffer;     /*!< Receiving Data buffer*/
  PDM_MPI_Request      *sendRequest;    /*!< Send Request (size = lComm) */
  PDM_MPI_Request      *recvRequest;    /*!< Receive Request (size = lComm) */
  void            **field;          /*!< Send Field */
  void            **ghostField;     /*!< Receive Field on ghost elements */

  int               tData;          /*!< Exchanged data type */

  int  nSendElt;    /*!< Number of sent element for ghost elements  */

  int *sendEltIdx;  /*!< Send Element index to send a field on sended point
                               (size = \ref lComm + 1 */

  int *sendElt;     /*!< Local number of sent elements for ghost elements
                               (size = \ref nSendElt */

  int *sendEltPart; /* Local number of sent elements for ghost elements
                               (size = \ref nSendElt */

  int *ghostEltIdx;     /*!< Ghost elements index
                          (size = \ref lComm + 1) */

  int nGhostElt;        /* Number of ghost element */


  int *nGhostEltPart;        /* Number of ghost element for each part
                                (size = \ref nPart */

  int **ghostEltPartIdx;  /*!< Elements connected to each part ghost elements Idx for each part
                            (size for each part = (nGhostEltPart + 1) */

  int **ghostEltPart2GhostElt;     /*!< Ghost element part to ghost element
                                     (size for each part = \ref nGhostEltPart */

  int **ghostEltPartElt;   /*!< partition number of Elements connected to each ghost element
                            (size for each part =  ghostEltPartIdx[iPart][ghostEltPartIdx[iPart][nGhostEltPart[iPart]]] */


};

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GRAPH_BOUND_PRIV_H__ */
