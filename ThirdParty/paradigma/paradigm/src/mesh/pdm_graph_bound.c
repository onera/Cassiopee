/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_graph_bound.h"
#include "pdm_graph_bound_priv.h"
#include "pdm_part_bound.h"
#include "pdm_part_bound_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

const int nDataExchCreate = 5; /*!< Number of exchanged data for creation */

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Return an intialized \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  comm         MPI communicator
 * \param [in]  nPart        Number of partitions
 * \param [in]  partBound    partition boundaries (size : \ref nPart)
 *
 * \return      A new initialized \ref PDM_graph_bound_t structure
 *
 */

PDM_graph_bound_t *
PDM_graph_bound_create
(
const PDM_MPI_Comm           comm,
const int                nPart,
      PDM_part_bound_t **partBound
)
{
  int myRank;
  int lComm;

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &lComm);

  PDM_graph_bound_t * graph_bound =
    (PDM_graph_bound_t *) malloc (sizeof(PDM_graph_bound_t));

  graph_bound->comm  = comm;
  graph_bound->nPart = nPart;

  graph_bound->exchRank = NULL;
  graph_bound->partBound = NULL;
  graph_bound->sendBuffer = NULL;
  graph_bound->sendRequest = NULL;
  graph_bound->recvRequest = NULL;
  graph_bound->field = NULL;
  graph_bound->sendEltIdx = NULL;
  graph_bound->sendElt = NULL;
  graph_bound->sendEltPart = NULL;
  graph_bound->ghostEltIdx = NULL;

  graph_bound->lComm = lComm;
  graph_bound->partBound =
    (PDM_part_bound_t **) malloc (nPart * sizeof(PDM_part_bound_t *));

  memcpy (graph_bound->partBound, partBound, nPart * sizeof(PDM_part_bound_t *));

  /*
   * Dump pdm_part_boud
   */

  if (1 == 0) {
    for (int i = 0; i < nPart; i++) {
      PDM_part_bound_dump (graph_bound->partBound[i]);
    }
  }

  /*
   * Look for connected ranks
   */

  int *offeredEltsRankIdx = (int *) malloc ((lComm + 1) * sizeof(int));
  for (int i = 0; i < lComm + 1; i++) {
    offeredEltsRankIdx[i] = 0;
  }

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *part_bound = graph_bound->partBound[i];
    int nEltPartBound = PDM_part_bound_n_elt_bound_get (part_bound);

    for (int j = 0; j < nEltPartBound; j++) {

      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (part_bound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get (part_bound,
                                        j+1,
                                        k,
                                        &iProc,
                                        &iPart,
                                        &iElt,
                                        &iDistElt);

        offeredEltsRankIdx[iProc+1]++;

      }
    }
  }

  graph_bound->nExchRank = 0;
  graph_bound->lBuffer = 0;
  for (int i = 0; i < lComm; i++) {
    if ((i != myRank) && (offeredEltsRankIdx[i+1] > 0)) {
      graph_bound->nExchRank += 1;
      graph_bound->lBuffer += offeredEltsRankIdx[i+1];
    }
    offeredEltsRankIdx[i+1] = offeredEltsRankIdx[i+1] +
                              offeredEltsRankIdx[i];
  }

  /*
   * Build list of exchanged ranks
   */

  graph_bound->exchRank= (int *) malloc (graph_bound->nExchRank * sizeof(int));
  graph_bound->nExchRank = 0;

  for (int i = 0; i < lComm; i++) {
    if ((offeredEltsRankIdx[i+1] > offeredEltsRankIdx[i]) && (i != myRank)) {
        graph_bound->exchRank[graph_bound->nExchRank++] = i;
    }
  }

  /*
   * Exchange number of offered elements for each connected element
   * to connected ranks
   */

  int *recvEltN    = (int *) malloc (lComm * sizeof(int));
  for (int i = 0; i < lComm; i++) {
    recvEltN[i] = 0;
  }

  graph_bound->sendRequest =
    (PDM_MPI_Request*) malloc (graph_bound->nExchRank * sizeof(PDM_MPI_Request));
  graph_bound->recvRequest =
    (PDM_MPI_Request*) malloc (graph_bound->nExchRank * sizeof(PDM_MPI_Request));

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    graph_bound->sendRequest[i] = PDM_MPI_REQUEST_NULL;
    graph_bound->recvRequest[i] = PDM_MPI_REQUEST_NULL;
  }

  int *sendOfferedElts =
    (int *) malloc (nDataExchCreate * offeredEltsRankIdx[lComm] * sizeof(int));
  int *recvOfferedElts =
    (int *) malloc (nDataExchCreate * offeredEltsRankIdx[lComm] * sizeof(int));

  for (int i = 0; i < nDataExchCreate * offeredEltsRankIdx[lComm]; i++) {
    sendOfferedElts[i] = -999;
    recvOfferedElts[i] = -999;
  }

  for (int i = 0; i < nPart; i++) {

    PDM_part_bound_t *_partBound = graph_bound->partBound[i];
    int nEltPartBound = PDM_part_bound_n_elt_bound_get(_partBound);

    for (int j = 0; j < nEltPartBound; j++) {

      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (_partBound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      int nOfferElt = PDM_part_bound_n_offer_elt_get (_partBound,
                                                      j+1);

      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get(_partBound,
                                       j+1,
                                       k,
                                       &iProc,
                                       &iPart,
                                       &iElt,
                                       &iDistElt);
        if (iProc != myRank) {
          int idx = nDataExchCreate *
            (offeredEltsRankIdx[iProc] + recvEltN[iProc]++);

          sendOfferedElts[idx++] = i;
          sendOfferedElts[idx++] = localElt;
          sendOfferedElts[idx++] = nOfferElt;
          sendOfferedElts[idx++] = iPart;
          sendOfferedElts[idx++] = iElt;
        }
      }
    }
  }

  free(recvEltN);

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];
    int idx =   nDataExchCreate * offeredEltsRankIdx[iProc];
    int count = nDataExchCreate * offeredEltsRankIdx[iProc+1] - idx;

    PDM_MPI_Irecv((void *) (recvOfferedElts + idx),
              count,
              PDM_MPI_INT,
              iProc,
              tag,
              comm,
              &graph_bound->recvRequest[i]);

  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];
    int idx =   nDataExchCreate * offeredEltsRankIdx[iProc];
    int count = nDataExchCreate * offeredEltsRankIdx[iProc+1] - idx;

    PDM_MPI_Issend((void *) (sendOfferedElts + idx),
               count,
               PDM_MPI_INT,
               iProc,
               tag,
               comm,
               &graph_bound->sendRequest[i]);

  }

  int idxMyrank = nDataExchCreate * offeredEltsRankIdx[myRank];
  int idxMyrank1 = nDataExchCreate * offeredEltsRankIdx[myRank];

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *_partBound = graph_bound->partBound[i];
    int nEltPartBound = PDM_part_bound_n_elt_bound_get(_partBound);
    for (int j = 0; j < nEltPartBound; j++) {
      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get(_partBound,
                                   j+1,
                                   &localElt,
                                   &nConnectedElt);

      int nOfferElt = PDM_part_bound_n_offer_elt_get (_partBound,
                                                      j+1);

      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get(_partBound,
                                       j+1,
                                       k,
                                       &iProc,
                                       &iPart,
                                       &iElt,
                                       &iDistElt);
        if (iProc == myRank) {
          sendOfferedElts[idxMyrank++] = i;
          sendOfferedElts[idxMyrank++] = localElt;
          sendOfferedElts[idxMyrank++] = nOfferElt;
          sendOfferedElts[idxMyrank++] = iPart;
          sendOfferedElts[idxMyrank++] = iElt;
          recvOfferedElts[idxMyrank1++] = i;
          recvOfferedElts[idxMyrank1++] = localElt;
          recvOfferedElts[idxMyrank1++] = nOfferElt;
          recvOfferedElts[idxMyrank1++] = iPart;
          recvOfferedElts[idxMyrank1++] = iElt;
        }
      }
    }
  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->recvRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->recvRequest[i]);
      graph_bound->recvRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->sendRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->sendRequest[i]);
      graph_bound->sendRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  free (sendOfferedElts);

  /*
   * Allocation arrays to storage offered elements
   */

  graph_bound->nGhostElt = 0;

  graph_bound->ghostEltIdx =
    (int *) malloc ((lComm + 1) * sizeof(int));

  for (int i = 0; i < lComm + 1; i++) {
    graph_bound->ghostEltIdx[i] = 0;
  }

  idxMyrank1 = 0;
  for (int i = 0; i < lComm; i++) {

    for (int k = offeredEltsRankIdx[i]; k < offeredEltsRankIdx[i+1]; k++) {

      int nOfferElt    = recvOfferedElts[nDataExchCreate*k + 2];
      graph_bound->nGhostElt += nOfferElt;
      graph_bound->ghostEltIdx[i+1] += nOfferElt;

    }
  }

  for (int i = 0; i < lComm; i++) {
    graph_bound->ghostEltIdx[i+1] += graph_bound->ghostEltIdx[i];
  }

  int *ghostEltEltIdx =
    (int *) malloc ((graph_bound->nGhostElt + 1) * sizeof(int));

  ghostEltEltIdx[0] = 0;
  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    ghostEltEltIdx[i+1] = ghostEltEltIdx[i] + 1;
  }

  int *ghostEltEltPart =
    (int *) malloc (graph_bound->nGhostElt * sizeof(int));

  int *ghostEltElt =
    (int *) malloc (graph_bound->nGhostElt * sizeof(int));

  idxMyrank1 = 0;
  graph_bound->nGhostElt = 0;

  for (int i = 0; i < lComm; i++) {

    for (int k = offeredEltsRankIdx[i]; k < offeredEltsRankIdx[i+1]; k++) {

      int nOfferElt    = recvOfferedElts[nDataExchCreate*k+2];
      int iLocalPart   = recvOfferedElts[nDataExchCreate*k+3];
      int iLocalElt    = recvOfferedElts[nDataExchCreate*k+4];

      for (int j = 0; j < nOfferElt; j++) {
        ghostEltEltPart[graph_bound->nGhostElt] = iLocalPart;
        ghostEltElt[graph_bound->nGhostElt] = iLocalElt;
        graph_bound->nGhostElt++;
      }
    }
  }

  /*
   * Exchange ghost element absolute number
   */

  for (int i = 0; i < lComm + 1; i++) {
    offeredEltsRankIdx[i] = 0;
  }

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *part_bound = graph_bound->partBound[i];
    int nEltPartBound = PDM_part_bound_n_elt_bound_get(part_bound);

    for (int j = 0; j < nEltPartBound; j++) {

      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (part_bound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get (part_bound,
                                        j+1,
                                        k,
                                        &iProc,
                                        &iPart,
                                        &iElt,
                                        &iDistElt);

        int nOfferElt = PDM_part_bound_n_offer_elt_get (part_bound, j+1);
        offeredEltsRankIdx[iProc+1] += nOfferElt;

      }
    }
  }

  for (int i = 0; i < lComm; i++) {
    offeredEltsRankIdx[i+1] = offeredEltsRankIdx[i+1] +
                              offeredEltsRankIdx[i];
  }

  graph_bound->nSendElt = offeredEltsRankIdx[lComm];

  PDM_g_num_t *gNumOfferedEltsSend =
    (PDM_g_num_t *) malloc (graph_bound->nSendElt * sizeof(PDM_g_num_t));

  graph_bound->sendEltIdx = offeredEltsRankIdx;

  graph_bound->sendElt =
    (int *) malloc (graph_bound->nSendElt * sizeof(int));

  graph_bound->sendEltPart =
    (int *) malloc (graph_bound->nSendElt * sizeof(int));

  PDM_g_num_t *gNumOfferedEltsRecv =
    (PDM_g_num_t *) malloc (graph_bound->nGhostElt * sizeof(PDM_g_num_t));

  int *nOfferedEltProc =
    (int *) malloc (lComm * sizeof(int));

  for (int i = 0; i < lComm; i++) {
    nOfferedEltProc[i] = 0;
  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    gNumOfferedEltsRecv[i] = -999;
  }

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *part_bound = graph_bound->partBound[i];
    int nEltPartBound = PDM_part_bound_n_elt_bound_get(part_bound);

    for (int j = 0; j < nEltPartBound; j++) {

      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (part_bound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get (part_bound,
                                        j+1,
                                        k,
                                        &iProc,
                                        &iPart,
                                        &iElt,
                                        &iDistElt);

        int nOfferElt = PDM_part_bound_n_offer_elt_get (part_bound, j+1);

        for (int k1 = 0; k1 < nOfferElt; k1++) {
          int        lNum;
          PDM_g_num_t gNum;
          PDM_part_bound_offer_elt_get (part_bound,
                                        j+1,
                                        k1,
                                        &lNum,
                                        &gNum);
          int idx = offeredEltsRankIdx[iProc]+nOfferedEltProc[iProc]++;
          graph_bound->sendEltPart[idx] = i;
          graph_bound->sendElt[idx]     = lNum;
          gNumOfferedEltsSend[idx]       = gNum;

        }
      }
    }
  }

  if (0 == 1) {
    PDM_printf ("gnumofferedelt : ");
    for (int i = 0; i < lComm; i++) {
      for (int j = offeredEltsRankIdx[i]; j < offeredEltsRankIdx[i+1]; j++) {
        PDM_printf (" "PDM_FMT_G_NUM" ", gNumOfferedEltsSend[j]);
      }
      PDM_printf (" §§ ");
    }
    PDM_printf ("\n");
  }

  free (nOfferedEltProc);

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];
    int idx   = graph_bound->ghostEltIdx[iProc];
    int count = graph_bound->ghostEltIdx[iProc+1] -
                graph_bound->ghostEltIdx[iProc];

    PDM_MPI_Irecv((void *) (gNumOfferedEltsRecv + idx),
              count,
              PDM__PDM_MPI_G_NUM,
              iProc,
              tag,
              comm,
              &graph_bound->recvRequest[i]);

  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];
    int idx   = offeredEltsRankIdx[iProc];
    int count = offeredEltsRankIdx[iProc+1] - offeredEltsRankIdx[iProc];

    PDM_MPI_Issend((void *) (gNumOfferedEltsSend + idx),
               count,
               PDM__PDM_MPI_G_NUM,
               iProc,
               tag,
               comm,
               &graph_bound->sendRequest[i]);

  }

  idxMyrank = offeredEltsRankIdx[myRank];
  idxMyrank1 = graph_bound->ghostEltIdx[myRank];

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *_partBound = graph_bound->partBound[i];
    int nEltPartBound = PDM_part_bound_n_elt_bound_get(_partBound);
    for (int j = 0; j < nEltPartBound; j++) {
      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (_partBound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      int nOfferElt = PDM_part_bound_n_offer_elt_get (_partBound,
                                                      j+1);

      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get (_partBound,
                                        j+1,
                                        k,
                                        &iProc,
                                        &iPart,
                                        &iElt,
                                        &iDistElt);

        if (iProc == myRank) {
          for (int k1 = 0; k1 < nOfferElt; k1++) {
            gNumOfferedEltsRecv[idxMyrank1++] =
              gNumOfferedEltsSend[idxMyrank++];
          }
        }
      }
    }
  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->recvRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->recvRequest[i]);
      graph_bound->recvRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->sendRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->sendRequest[i]);
      graph_bound->sendRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  free (gNumOfferedEltsSend);

  if (0 == 1) {
    PDM_printf ("gNumOfferedEltsRecv : ");
    int k = 0;
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      if (i == graph_bound->ghostEltIdx[k]) {
        k++;
        PDM_printf (" $$ ");
      }
      PDM_printf (" "PDM_FMT_G_NUM, gNumOfferedEltsRecv[i] -1);
    }
    PDM_printf ("\n");
  }

  /*
   * Remove double points
   *   - store ghost element in hash table
   *   - look for double element in hash table
   *   - tagGhostElt send to sender processus
   *
   */

  PDM_g_num_t _nKeys = PDM_part_bound_n_total_offer_elt_get (graph_bound->partBound[0]) / lComm;

  int nKeys = (int) (_nKeys);

  int *hashTableIdx = (int *) malloc (sizeof(int) * (nKeys+1));
  for (int i = 0; i < nKeys+1; i++) {
    hashTableIdx[i] = 0;
  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    PDM_g_num_t _key = gNumOfferedEltsRecv[i] % nKeys;
    int key = (int) (_key);
    hashTableIdx[key+1]++;
  }

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *_partBound = graph_bound->partBound[i];
    int nLocalOfferElt = PDM_part_bound_n_local_offer_elt_get (_partBound);
    const PDM_g_num_t *localOfferEltLnToGn =
      PDM_part_bound_local_offer_elt_ln_to_gn_get (_partBound);
    for (int j = 0; j < nLocalOfferElt; j++) {
      PDM_g_num_t _key = localOfferEltLnToGn[j] % nKeys;
      int key =  (int) _key;
      hashTableIdx[key+1]++;
    }
  }

  for (int i = 0; i < nKeys; i++) {
    hashTableIdx[i+1] += hashTableIdx[i];
  }

  /*
   * Storage in hash table
   */

  int *hashTableN = (int *) malloc (sizeof(int) * nKeys);
  PDM_g_num_t *hashTableGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * hashTableIdx[nKeys]);
  int *hashTableDataLoc         = (int *) malloc (sizeof(int) * hashTableIdx[nKeys]);
  int *hashTableDataNumLoc      = (int *) malloc (sizeof(int) * hashTableIdx[nKeys]);

  for (int i = 0; i < nKeys; i++) {
    hashTableN[i] = 0;
  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    PDM_g_num_t gNum = gNumOfferedEltsRecv[i];
    PDM_g_num_t _key = gNum % nKeys;
    int key = (int) _key;
    int idx = hashTableIdx[key] + hashTableN[key]++;
    hashTableGnum[idx] = gNum;
    hashTableDataLoc[idx] = -1;
    hashTableDataNumLoc[idx] = i;
  }

  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *_partBound = graph_bound->partBound[i];
    int nLocalOfferElt = PDM_part_bound_n_local_offer_elt_get (_partBound);
    const PDM_g_num_t *localOfferEltLnToGn =
      PDM_part_bound_local_offer_elt_ln_to_gn_get (_partBound);
    for (int j = 0; j < nLocalOfferElt; j++) {
      PDM_g_num_t gNum = localOfferEltLnToGn[j];
      PDM_g_num_t _key = gNum % nKeys;
      int key = (int) _key;
      int idx = hashTableIdx[key] + hashTableN[key]++;
      hashTableGnum[idx] = gNum;
      hashTableDataLoc[idx] = i;
      hashTableDataNumLoc[idx] = j;
    }
  }

  free (hashTableN);

  /*
   * TagGhost link about an other ghost element (>0) or a new local (<0)
   */

  int *tagGhostElt = malloc (graph_bound->nGhostElt * sizeof(int));
  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    tagGhostElt[i] = 0;
  }

  /*
   * Remove elements already in local offered elements
   */

  int **eltToPartBound = (int **) malloc(sizeof(int *) * nPart);
  for (int i = 0; i < nPart; i++) {
    PDM_part_bound_t *_partBound = graph_bound->partBound[i];

    int nElt = PDM_part_bound_n_elt_get (_partBound);

    eltToPartBound[i] = (int *) malloc(sizeof(int) * nElt);
    int *_eltToPartBound = eltToPartBound[i];

    for (int k = 0; k < nElt; k++) {
      _eltToPartBound[k] = -1;
    }

    int nEltPartBound = PDM_part_bound_n_elt_bound_get(_partBound);
    for (int j = 0; j < nEltPartBound; j++) {
      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (_partBound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);
      _eltToPartBound[localElt-1] = j;
    }

  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    PDM_g_num_t gNumGhost = gNumOfferedEltsRecv[i];
    for (int j = ghostEltEltIdx[i]; j < ghostEltEltIdx[i+1]; j++) {
      int iLocalPart = ghostEltEltPart[j];
      int iLocalElt = ghostEltElt[j];

      PDM_part_bound_t *_partBound = graph_bound->partBound[iLocalPart];
      int iPartBound = eltToPartBound[iLocalPart][iLocalElt-1];

      int nOfferElt = PDM_part_bound_n_offer_elt_get (_partBound, iPartBound + 1);

      for (int k1 = 0; k1 < nOfferElt; k1++) {
        int        lNumOffer;
        PDM_g_num_t gNumOffer;
        PDM_part_bound_offer_elt_get (_partBound,
                                      iPartBound +1 ,
                                      k1,
                                      &lNumOffer,
                                      &gNumOffer);

        if (gNumGhost == gNumOffer) {
          tagGhostElt[i] = -1;
          break;
        }
      }
      if (tagGhostElt[i] == -1) {
        break;
      }
    }
  }

  free (gNumOfferedEltsRecv);

  for (int i = 0; i < nPart; i++) {
    free (eltToPartBound[i]);
  }
  free (eltToPartBound);

  /*
   * Remove double elements
   */

  int sNewLocalGhost = 10;
  int nNewLocalGhost = 0;
  int *newLocalGhost = (int *) malloc (sizeof(int) * 2 * sNewLocalGhost);

  /*
   * Launch irecv about tag
   */

  int *recvTagGhostElt = (int *) malloc (graph_bound->nSendElt * sizeof(int));
  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];
    int idx   = offeredEltsRankIdx[iProc];
    int count = offeredEltsRankIdx[iProc+1] - offeredEltsRankIdx[iProc];

    PDM_MPI_Irecv ((void *) (recvTagGhostElt + idx),
               count,
               PDM_MPI_INT,
               iProc,
               tag,
               comm,
               &graph_bound->recvRequest[i]);

  }

  /*
   * Compute tag
   */

  for (int i = 0; i < nKeys; i++) {
    int iBegLocal = hashTableIdx[i+1];
    int nNewLocalGhostPrevious = nNewLocalGhost;
    for (int j = hashTableIdx[i]; j < hashTableIdx[i+1]; j++) {
      int dataLoc1 = hashTableDataLoc[j];
      if (dataLoc1 >= 0) {
        iBegLocal = j;
        break;
      }
    }

    for (int j = hashTableIdx[i]; j < iBegLocal; j++) {
      PDM_g_num_t gNum1 = hashTableGnum[j];
      int dataNumLoc1 = hashTableDataNumLoc[j];
      int iSLocalProc1 = (dataNumLoc1 >= graph_bound->ghostEltIdx[myRank]) &&
                         (dataNumLoc1 < graph_bound->ghostEltIdx[myRank+1]);

      /*
       * If not tagged : Look for in other ghost elements
       */

      if (tagGhostElt[dataNumLoc1] == 0) {

        for (int k = j+1; k < iBegLocal; k++) {

          PDM_g_num_t gNum2 = hashTableGnum[k];
          int dataNumLoc2 = hashTableDataNumLoc[k];
          int iSLocalProc2 = (dataNumLoc2 >= graph_bound->ghostEltIdx[myRank]) &&
            (dataNumLoc2 < graph_bound->ghostEltIdx[myRank+1]);

          if ((gNum1 == gNum2) && (tagGhostElt[dataNumLoc2] == 0)) {

            if (iSLocalProc1 || !iSLocalProc2) {

              tagGhostElt[dataNumLoc2] = dataNumLoc1 + 1;

            }

            else  {

              tagGhostElt[dataNumLoc1] = dataNumLoc2 + 1;

            }
          }
        }
      }

      /*
       * If not tagged : Look for in new local ghost elements
       */

      if (!iSLocalProc1 && (tagGhostElt[dataNumLoc1] == 0)) {

        for (int k = nNewLocalGhostPrevious; k < nNewLocalGhost; k++) {
          int iPartLoc = newLocalGhost[2*k];
          int iNumLoc = newLocalGhost[2*k+1];
          PDM_part_bound_t *_partBound = graph_bound->partBound[iPartLoc];
          const PDM_g_num_t gNum2 =
            PDM_part_bound_local_offer_elt_ln_to_gn_get (_partBound)[iNumLoc];
          if (gNum1 == gNum2) {
            tagGhostElt[dataNumLoc1] = -(k+1) - 1;
            break;
          }
        }

      }

      /*
       * If not tagged : Look for in other local elements
       */

      if (!iSLocalProc1 && tagGhostElt[dataNumLoc1] == 0) {

        for (int k = iBegLocal; k < hashTableIdx[i+1]; k++) {
          PDM_g_num_t gNum2 = hashTableGnum[k];
          int dataLoc2     = hashTableDataLoc[k];
          int dataNumLoc2  = hashTableDataNumLoc[k];

          if (gNum1 == gNum2) {

            if (nNewLocalGhost >= sNewLocalGhost) {
              sNewLocalGhost *= 2;
              newLocalGhost = (int *) realloc (newLocalGhost,
                                               sizeof(int) * 2 * sNewLocalGhost);
            }

            newLocalGhost[2*nNewLocalGhost]   = dataLoc2 ;
            newLocalGhost[2*nNewLocalGhost+1] = dataNumLoc2 + 1;
            nNewLocalGhost +=  1;
            tagGhostElt[dataNumLoc1] = -nNewLocalGhost - 1;
            break;

          }
        }
      }

    }
  }

  free (hashTableIdx);
  free (hashTableGnum);
  free (hashTableDataLoc);
  free (hashTableDataNumLoc);

  if (0 == 1) {

    PDM_printf ("tagghostelt : ");
    int k = 0;
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      if (i == graph_bound->ghostEltIdx[k]) {
        k++;
        PDM_printf (" $$ ");
      }
      PDM_printf (" %d", tagGhostElt[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("newLocalGhost : ");
    for (int i = 0; i < nNewLocalGhost; i++) {
      PDM_printf (" %d %d,", newLocalGhost[2*i], newLocalGhost[2*i+1]);
    }
    PDM_printf ("\n");

    PDM_printf ("localgnum : ");

    for (int i = 0; i < nPart; i++) {
      PDM_part_bound_t *_partBound = graph_bound->partBound[i];
      int nLocalOfferElt = PDM_part_bound_n_local_offer_elt_get (_partBound);
      const PDM_g_num_t *localOfferEltLnToGn =
        PDM_part_bound_local_offer_elt_ln_to_gn_get (_partBound);
      for (int j = 0; j < nLocalOfferElt; j++) {
        PDM_g_num_t gNum = localOfferEltLnToGn[j];
        PDM_printf (" %d %d " PDM_FMT_G_NUM",", i, j, gNum);
      }
    }
    PDM_printf ("\n");

    PDM_printf("ghostEltEltPart     :");
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      for (int j = ghostEltEltIdx[i]; j < ghostEltEltIdx[i+1]; j++) {
        PDM_printf(" %d", ghostEltEltPart[j]);
      }
      PDM_printf(" $$ ");
    }
    PDM_printf("\n");

    PDM_printf("ghostEltElt     :");
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      for (int j = ghostEltEltIdx[i]; j < ghostEltEltIdx[i+1]; j++) {
        PDM_printf(" %d", ghostEltElt[j]);
      }
      PDM_printf(" $$ ");
    }
    PDM_printf("\n");
  }

  /*
   * Launch issend about tag
   */

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];
    int idx   = graph_bound->ghostEltIdx[iProc];
    int count = graph_bound->ghostEltIdx[iProc+1] -
                graph_bound->ghostEltIdx[iProc];

    PDM_MPI_Issend ((void *) (tagGhostElt + idx),
                count,
                PDM_MPI_INT,
                iProc,
                tag,
                comm,
                &graph_bound->sendRequest[i]);

  }

  /*
   * Update recvTag for current Rank
   */

  int idx1MyRank   = graph_bound->ghostEltIdx[myRank];
  int count1MyRank = graph_bound->ghostEltIdx[myRank+1]
                   - graph_bound->ghostEltIdx[myRank];

  int idx2MyRank   = offeredEltsRankIdx[myRank];
  int count2MyRank = offeredEltsRankIdx[myRank+1]
                   - offeredEltsRankIdx[myRank];

  assert (count1MyRank == count2MyRank);

  for (int i = 0; i < count1MyRank; i++) {
    recvTagGhostElt[idx2MyRank + i] = tagGhostElt[idx1MyRank + i];
  }

  /*
   * Update ghost element structures :
   */

  /* - ghostEltIdx */

  int *newGhostEltIdx = (int *) malloc (sizeof(int) * (lComm+1));
  int *oldToNewGhost = (int *) malloc (sizeof(int) * graph_bound->nGhostElt);

  memcpy (oldToNewGhost, tagGhostElt, sizeof(int) * graph_bound->nGhostElt);

  int nTotalNewGhost = 0;

  for (int i = 0; i < lComm + 1; i++) {
    newGhostEltIdx[i] = 0;
  }

  for (int i = 0; i < lComm; i++) {
    for (int j = graph_bound->ghostEltIdx[i];
         j < graph_bound->ghostEltIdx[i+1];
         j++) {
      if (tagGhostElt[j] == 0) {
        newGhostEltIdx[i+1]++;
        nTotalNewGhost++;
      }
    }
  }

  int begNewLocalGhost = newGhostEltIdx[myRank+1];
  newGhostEltIdx[myRank+1] += nNewLocalGhost;
  nTotalNewGhost += nNewLocalGhost;

  for (int i = 0; i < lComm; i++) {
    newGhostEltIdx[i+1] += newGhostEltIdx[i];
  }

  int *nNewGhostElt = (int *) malloc (sizeof(int) * lComm);
  for (int i = 0; i < lComm; i++) {
    nNewGhostElt[i] = 0;
  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    if (oldToNewGhost[i] > 0) {
      while (oldToNewGhost[oldToNewGhost[i]-1] > 0) {
        oldToNewGhost[i] = oldToNewGhost[oldToNewGhost[i]-1];
      }
    }
  }

  for (int i = 0; i < lComm; i++) {
    for (int j = graph_bound->ghostEltIdx[i];
         j < graph_bound->ghostEltIdx[i+1];
         j++) {
      if (tagGhostElt[j] == 0) {
        oldToNewGhost[j] = newGhostEltIdx[i] + nNewGhostElt[i]++;
      }
      else if (tagGhostElt[j]  < -1) {
        oldToNewGhost[j] = newGhostEltIdx[myRank] + begNewLocalGhost -(tagGhostElt[j]+1+1);;
      }
    }
  }

  free (nNewGhostElt);

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    if (tagGhostElt[i] > 0) {
      oldToNewGhost[i] = oldToNewGhost[oldToNewGhost[i]-1];
    }
  }

  if (0 == 1) {

    PDM_printf ("oldToNewGhost : ");
    int k = 0;
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      if (i == graph_bound->ghostEltIdx[k]) {
        k++;
        PDM_printf (" $$ ");
      }
      PDM_printf (" %d", oldToNewGhost[i]);
    }
    PDM_printf ("\n");
  }

  int *oldGhostEltIdx = graph_bound->ghostEltIdx;
  graph_bound->ghostEltIdx = newGhostEltIdx;

  /* - ghostEltEltIdx
     - ghostEltElt
     - ghostEltEltPart */

  int *newGhostEltEltIdx =
    (int *) malloc ((nTotalNewGhost + 1) * sizeof(int));

  int *nNewGhostEltElt =
    (int *) malloc (nTotalNewGhost * sizeof(int));

  for (int i = 0; i < nTotalNewGhost + 1; i++) {
    newGhostEltEltIdx[i] = 0;
  }

  for (int i = 0; i < nTotalNewGhost; i++) {
    nNewGhostEltElt[i] = 0;
  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    int newGhost = oldToNewGhost[i];
    if (newGhost >= 0) {
      newGhostEltEltIdx[newGhost+1]++;
    }
  }

  for (int i = 0; i < nTotalNewGhost; i++) {
    newGhostEltEltIdx[i+1] += newGhostEltEltIdx[i];
  }

  int *newGhostEltEltPart =
    (int *) malloc (newGhostEltEltIdx[nTotalNewGhost] * sizeof(int));

  int *newGhostEltElt =
    (int *) malloc (newGhostEltEltIdx[nTotalNewGhost] * sizeof(int));

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    int newGhost = oldToNewGhost[i];
    if (newGhost >= 0) {
      int idx = newGhostEltEltIdx[newGhost] + nNewGhostEltElt[newGhost]++;
      newGhostEltEltPart[idx] = ghostEltEltPart[i];
      newGhostEltElt[idx] = ghostEltElt[i];
    }
  }

  /*
   * Remove double elt in ghostEltElt
   */

  int *tagEltElt = (int *) malloc (sizeof(int) * newGhostEltEltIdx[nTotalNewGhost]);
  int *newGhostEltEltIdx2 = (int *) malloc (sizeof(int) * (nTotalNewGhost + 1));

  for (int i = 0; i < nTotalNewGhost + 1; i++) {
    newGhostEltEltIdx2[i] = 0;
  }

  for (int i = 0; i < newGhostEltEltIdx[nTotalNewGhost]; i++) {
    tagEltElt[i] = 1;
  }

  for (int i = 0; i < nTotalNewGhost; i++) {
    for (int j = newGhostEltEltIdx[i]; j < newGhostEltEltIdx[i+1]; j++) {
      int part1 = newGhostEltEltPart[j];
      int elt1 = newGhostEltElt[j];
      if (tagEltElt[j] == 1) {
        for (int k = j+1; k < newGhostEltEltIdx[i+1]; k++) {
          int part2 = newGhostEltEltPart[k];
          int elt2 = newGhostEltElt[k];
          if ((part1 == part2) && (elt1 == elt2)) {
            tagEltElt[k] = 0;
            break;
          }
        }
      }
    }
  }

  for (int i = 0; i < nTotalNewGhost; i++) {
    for (int j = newGhostEltEltIdx[i]; j < newGhostEltEltIdx[i+1]; j++) {
      if (tagEltElt[j] == 1) {
        newGhostEltEltIdx2[i+1]++;
      }
    }
  }

  for (int i = 0; i < nTotalNewGhost; i++) {
    newGhostEltEltIdx2[i+1] += newGhostEltEltIdx2[i];
  }

  int idxTag = 0;
  for (int i = 0; i < newGhostEltEltIdx[nTotalNewGhost]; i++) {
    if (tagEltElt[i] == 1) {
      newGhostEltEltPart[idxTag] = newGhostEltEltPart[i];
      newGhostEltElt[idxTag++] = newGhostEltElt[i];
    }
  }

  newGhostEltEltPart = (int *) realloc (newGhostEltEltPart,
                                        sizeof (int) * newGhostEltEltIdx2[nTotalNewGhost]);

  newGhostEltElt = (int *) realloc (newGhostEltElt,
                                    sizeof (int) * newGhostEltEltIdx2[nTotalNewGhost]);

  free (tagEltElt);

  free (newGhostEltEltIdx);
  newGhostEltEltIdx = newGhostEltEltIdx2;

  free (oldToNewGhost);
  free (nNewGhostEltElt);
  free (ghostEltElt);
  free (ghostEltEltPart);
  free (ghostEltEltIdx);

  ghostEltElt = newGhostEltElt;
  ghostEltEltPart = newGhostEltEltPart;
  ghostEltEltIdx = newGhostEltEltIdx;
  graph_bound->nGhostElt = nTotalNewGhost;

  /*
   * Wait tag exchanges
   */

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->recvRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->recvRequest[i]);
        graph_bound->recvRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->sendRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->sendRequest[i]);
      graph_bound->sendRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  /*
   * Select offered element to send from received tag
   */

  /* - sendEltIdx */
  /* - sendElt */
  /* - sendEltPart */

  int *newSendEltIdx = (int *) malloc ((lComm + 1) * sizeof(int));
  int *nNewSendElt = (int *) malloc (lComm  * sizeof(int));

  newSendEltIdx[0] = 0;
  for (int i = 0; i < lComm; i++) {
    nNewSendElt[i] = 0;
    newSendEltIdx[i+1] = 0;
  }

  newSendEltIdx[myRank+1] = newGhostEltIdx[myRank+1] - newGhostEltIdx[myRank];

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];

    for (int j = offeredEltsRankIdx[iProc]; j < offeredEltsRankIdx[iProc+1]; j++) {
      if (recvTagGhostElt[j] == 0) {
        newSendEltIdx[iProc+1]++;
      }
    }
  }

  for (int i = 0; i < lComm; i++) {
    newSendEltIdx[i+1] += newSendEltIdx[i];
  }

  int *newSendEltPart = (int *) malloc (newSendEltIdx[lComm] * sizeof(int));
  int *newSendElt     = (int *) malloc (newSendEltIdx[lComm] * sizeof(int));

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];

    for (int j = offeredEltsRankIdx[iProc]; j < offeredEltsRankIdx[iProc+1]; j++) {
      if (recvTagGhostElt[j] == 0) {
        int idx = newSendEltIdx[iProc] + nNewSendElt[iProc]++;
        newSendElt[idx] = graph_bound->sendElt[j];
        newSendEltPart[idx] = graph_bound->sendEltPart[j];
      }
    }
  }

  for (int j = offeredEltsRankIdx[myRank]; j < offeredEltsRankIdx[myRank+1]; j++) {

    if (recvTagGhostElt[j] == 0) {
      int idx = newSendEltIdx[myRank] + nNewSendElt[myRank]++;
      newSendElt[idx] = graph_bound->sendElt[j];
      newSendEltPart[idx] = graph_bound->sendEltPart[j];
    }
  }

  for (int i = 0; i < nNewLocalGhost; i++) {
    int idx = newSendEltIdx[myRank] + nNewSendElt[myRank]++;
    newSendElt[idx] = newLocalGhost[2*i+1];
    newSendEltPart[idx] = newLocalGhost[2*i];
  }

  /*
   * Clean up
   */

  free (nNewSendElt);
  free (graph_bound->sendElt);
  free (graph_bound->sendEltPart);
  free (graph_bound->sendEltIdx);

  graph_bound->sendElt     = newSendElt;
  graph_bound->sendEltPart = newSendEltPart;
  graph_bound->sendEltIdx  = newSendEltIdx;
  graph_bound->nSendElt = graph_bound->sendEltIdx[lComm];

  free (oldGhostEltIdx);
  free (tagGhostElt);
  free (recvTagGhostElt);

  free (recvOfferedElts);
  free (newLocalGhost);

  /*
   * Build structure for each part
   */

  if (1 == 0) {

    PDM_printf("ghostEltEltPart     :");
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      for (int j = ghostEltEltIdx[i]; j < ghostEltEltIdx[i+1]; j++) {
        PDM_printf(" %d", ghostEltEltPart[j]);
      }
      PDM_printf(" $$ ");
    }
    PDM_printf("\n");

    PDM_printf("ghostEltElt     :");
    for (int i = 0; i < graph_bound->nGhostElt; i++) {
      for (int j = ghostEltEltIdx[i]; j < ghostEltEltIdx[i+1]; j++) {
        PDM_printf(" %d", ghostEltElt[j]);
      }
      PDM_printf(" $$ ");
    }
    PDM_printf("\n");
  }

  graph_bound->nGhostEltPart = (int *) malloc (nPart * sizeof(int));
  graph_bound->ghostEltPart2GhostElt = (int **) malloc (nPart * sizeof(int *));
  graph_bound->ghostEltPartIdx = (int **) malloc (nPart * sizeof(int *));
  graph_bound->ghostEltPartElt = (int **) malloc (nPart * sizeof(int *));

  int **tagGhostEltPart =  (int **) malloc (nPart * sizeof(int *));

  for (int i = 0; i < nPart; i++) {
    tagGhostEltPart[i] =  (int *) malloc (graph_bound->nGhostElt * sizeof(int *));
    graph_bound->nGhostEltPart[i] = 0;

    for (int j = 0; j < graph_bound->nGhostElt; j++) {
      tagGhostEltPart[i][j] = 0;
    }

  }

  for (int i = 0; i < graph_bound->nGhostElt; i++) {
    for (int j = ghostEltEltIdx[i]; j < ghostEltEltIdx[i+1]; j++) {
      tagGhostEltPart[ghostEltEltPart[j]][i]++;
    }
  }

  int lGhostEltPartElt = 0;
  for (int i = 0; i < nPart; i++) {
    for (int j = 0; j < graph_bound->nGhostElt; j++) {
      if (tagGhostEltPart[i][j] > 0) {
        lGhostEltPartElt += tagGhostEltPart[i][j];
        graph_bound->nGhostEltPart[i]++;
      }
    }

    int nGhostEltPart = graph_bound->nGhostEltPart[i];

    int *ghostEltPart2GhostElt = (int *) malloc (nGhostEltPart * sizeof(int));
    int *ghostEltPartIdx = (int *) malloc ((nGhostEltPart + 1) * sizeof(int));
    int *ghostEltPartElt = (int *) malloc (lGhostEltPartElt * sizeof(int));

    graph_bound->ghostEltPart2GhostElt[i] = ghostEltPart2GhostElt;
    graph_bound->ghostEltPartIdx[i] = ghostEltPartIdx;
    graph_bound->ghostEltPartElt[i] = ghostEltPartElt;

    int idx = 0;
    for (int j = 0; j < ghostEltEltIdx[graph_bound->nGhostElt]; j++) {
      if (ghostEltEltPart[j] == i) {
        ghostEltPartElt[idx++] = ghostEltElt[j];
      }
    }

    idx = 0;
    ghostEltPartIdx[0] = 0;
    for (int j = 0; j < graph_bound->nGhostElt; j++) {
      if (tagGhostEltPart[i][j] > 0) {
        ghostEltPart2GhostElt[idx] = j;
        ghostEltPartIdx[++idx] = tagGhostEltPart[i][j];
      }
    }

    for (int j = 0; j < nGhostEltPart; j++) {
      ghostEltPartIdx[j+1] += ghostEltPartIdx[j];
    }

    free (tagGhostEltPart[i]);

  }

  graph_bound->sendElt     = newSendElt;
  graph_bound->sendEltPart = newSendEltPart;
  graph_bound->sendEltIdx  = newSendEltIdx;

  free (ghostEltEltIdx);
  free (ghostEltEltPart);
  free (ghostEltElt);
  free (tagGhostEltPart);

  return graph_bound;

}


/**
 * \brief Return an intialized \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  graph_bound   Boundary graph
 * \param [in]  nComp         Number of composantes
 * \param [in]  tData         Data type
 * \param [in]  field         Data field
 * \param [out] ghostField    Ghost cell field
 *
 */

void
PDM_graph_bound_exch_data_init
(
      PDM_graph_bound_t *graph_bound,
const int                nComp,
const PDM_data_t         tData,
      void             **field,
      void             **ghostField
)
{
  if (graph_bound->sendBuffer != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error _graph_bound_exch_data_init :"
    "exchange already initialized\n");
    abort();
  }

  int myRank;
  int lComm;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &lComm);

  graph_bound->field      = field;
  graph_bound->ghostField = ghostField;
  graph_bound->tData      = tData;
  graph_bound->nComp      = nComp;

  size_t sizeExchType = 0;

  /*
   * Fill sending buffers, initialize communications and
   * build intra processus communication
   */

  switch (tData) {

  case PDM_INT : {
    int *_sendBuffer = (int *) malloc(sizeof(int) * graph_bound->nSendElt * nComp);
    int *_recvBuffer = (int *) malloc(sizeof(int) * graph_bound->nGhostElt * nComp);
    int **_field = (int **) field;
    sizeExchType = sizeof(int);

    graph_bound->sendBuffer = (void *) _sendBuffer;
    graph_bound->recvBuffer = (void *) _recvBuffer;

    int idx = 0;
    for (int i = 0; i < graph_bound->nSendElt; i++) {
      int iPart = graph_bound->sendEltPart[i];
      int iElt  = graph_bound->sendElt[i]-1;
      for (int k = 0; k < nComp; k++) {
        _sendBuffer[idx++] = _field[iPart][nComp*iElt+k];
      }
    }
    break;
  }

  case PDM_DOUBLE : {
    double *_sendBuffer = (double *) malloc(sizeof(double) * graph_bound->nSendElt * nComp);
    double *_recvBuffer = (double *) malloc(sizeof(double) * graph_bound->nGhostElt * nComp);
    double **_field = (double **) field;
    sizeExchType = sizeof(double);

    graph_bound->sendBuffer = (void *) _sendBuffer;
    graph_bound->recvBuffer = (void *) _recvBuffer;

    int idx = 0;
    for (int i = 0; i < graph_bound->nSendElt; i++) {
      int iPart = graph_bound->sendEltPart[i];
      int iElt  = graph_bound->sendElt[i]-1;
      for (int k = 0; k < nComp; k++) {
        _sendBuffer[idx++] = _field[iPart][nComp*iElt+k];
      }
    }
    break;
  }

  default : {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_exch_data_init :"
                    "this data type is not taking into account\n");
    abort();
  }

  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    int iProc = graph_bound->exchRank[i];

    unsigned char *sendPtr = (unsigned char *) graph_bound->sendBuffer;
    sendPtr +=  graph_bound->sendEltIdx[iProc] * sizeExchType * nComp;

    int sendCount = (int) (
      (graph_bound->sendEltIdx[iProc+1] - graph_bound->sendEltIdx[iProc])
      * (int) sizeExchType * nComp);

    unsigned char *recvPtr = (unsigned char *) graph_bound->recvBuffer;
    recvPtr += graph_bound->ghostEltIdx[iProc] * sizeExchType * nComp;

    int recvCount = (int) (
      (graph_bound->ghostEltIdx[iProc+1] - graph_bound->ghostEltIdx[iProc])
      * (int) sizeExchType * nComp);

    PDM_MPI_Irecv((void *) recvPtr,
              recvCount,
              PDM_MPI_UNSIGNED_CHAR,
              iProc,
              tag,
              graph_bound->comm,
              &graph_bound->recvRequest[i]);

    PDM_MPI_Issend((void *) sendPtr,
               sendCount,
               PDM_MPI_UNSIGNED_CHAR,
               iProc,
               tag,
               graph_bound->comm,
               &graph_bound->sendRequest[i]);
  }

  /*
   * Local copy
   */

  unsigned char *sendLocPtr = (unsigned char *) graph_bound->sendBuffer;
  sendLocPtr +=  graph_bound->sendEltIdx[myRank] * sizeExchType * nComp;

  unsigned char *recvLocPtr = (unsigned char *) graph_bound->recvBuffer;
  recvLocPtr += graph_bound->ghostEltIdx[myRank] * sizeExchType * nComp;

  int sendLocCount =
    (graph_bound->sendEltIdx[myRank+1] - graph_bound->sendEltIdx[myRank])
    * (int) sizeExchType * nComp;
  int recvLocCount =
    (graph_bound->ghostEltIdx[myRank+1] - graph_bound->ghostEltIdx[myRank])
    * (int) sizeExchType * nComp;

  assert (sendLocCount == recvLocCount);

  if (sendLocCount > 0) {
    memcpy ((void *) recvLocPtr, (void *) sendLocPtr, sendLocCount);
  }

}


/**
 * \brief Return an intialized \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  graph_bound   Boundary graph
 *
 */

void
PDM_graph_bound_exch_data_wait
(
PDM_graph_bound_t *graph_bound
)
{
  if (graph_bound->sendBuffer == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error _graph_bound_exch_data_wait :"
    "Exchange not initialized\n");
    abort();
  }

  int myRank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->recvRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->recvRequest[i]);
      graph_bound->recvRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  for (int i = 0; i < graph_bound->nExchRank; i++) {
    if (graph_bound->sendRequest[i] != PDM_MPI_REQUEST_NULL) {
      PDM_MPI_Wait(&graph_bound->sendRequest[i]);
      graph_bound->sendRequest[i] = PDM_MPI_REQUEST_NULL;
    }
  }

  switch (graph_bound->tData) {

  case PDM_INT : {
    int *_sendBuffer = (int *) graph_bound->sendBuffer;
    int *_recvBuffer = (int *) graph_bound->recvBuffer;
    int **ghostField = (int **) graph_bound->ghostField;

    free (_sendBuffer);

    for (int i = 0; i < graph_bound->nPart; i++) {
      int *ghostFieldPart        = ghostField[i];
      int *ghostEltPart2GhostElt = graph_bound->ghostEltPart2GhostElt[i];
      for (int j = 0; j < graph_bound->nGhostEltPart[i]; j++) {
        int ghostElt = ghostEltPart2GhostElt[j];
        int nComp = graph_bound->nComp;
        for (int k = 0; k < nComp; k++) {
          ghostFieldPart[j * nComp + k] = _recvBuffer[ghostElt * nComp + k];
        }
      }
    }

    free (_recvBuffer);
    break;
  }

  case PDM_DOUBLE : {
    double *_sendBuffer = (double *) graph_bound->sendBuffer;
    double *_recvBuffer = (double *) graph_bound->recvBuffer;
    double **ghostField = (double **) graph_bound->ghostField;

    free (_sendBuffer);

    for (int i = 0; i < graph_bound->nPart; i++) {
      double *ghostFieldPart        = ghostField[i];
      int *ghostEltPart2GhostElt = graph_bound->ghostEltPart2GhostElt[i];
      for (int j = 0; j < graph_bound->nGhostEltPart[i]; j++) {
        int ghostElt = ghostEltPart2GhostElt[j];
        int nComp = graph_bound->nComp;
        for (int k = 0; k < nComp; k++) {
          ghostFieldPart[j * nComp + k] = _recvBuffer[ghostElt * nComp + k];
        }
      }
    }

    free (_recvBuffer);
    break;
  }

  default : {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_exch_data_wait :"
    "this data type is not taking into account\n");
    abort();
  }

  }

  graph_bound->sendBuffer = NULL;
  graph_bound->recvBuffer = NULL;
}


/**
 * \brief free \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  graph_bound   Boundary graph
 *
 * \return      NULL
 *
 */

PDM_graph_bound_t *
PDM_graph_bound_free
(
PDM_graph_bound_t *graph_bound
)
{
  if (graph_bound != NULL) {
    if (graph_bound->exchRank != NULL)
      free (graph_bound->exchRank);
    if (graph_bound->sendBuffer != NULL)
      free (graph_bound->sendBuffer);
    if (graph_bound->sendRequest != NULL)
      free (graph_bound->sendRequest);
    if (graph_bound->recvRequest != NULL)
      free (graph_bound->recvRequest);
    if (graph_bound->sendEltIdx != NULL)
      free (graph_bound->sendEltIdx);
    if (graph_bound->sendElt != NULL)
      free (graph_bound->sendElt);
    if (graph_bound->sendEltPart != NULL)
      free (graph_bound->sendEltPart);
    if (graph_bound->ghostEltIdx != NULL)
      free (graph_bound->ghostEltIdx);
    if (graph_bound->partBound != NULL)
      free (graph_bound->partBound);
    for (int i = 0; i < graph_bound->nPart; i++) {
      free (graph_bound->ghostEltPartIdx[i]);
      free (graph_bound->ghostEltPartElt[i]);
      free (graph_bound->ghostEltPart2GhostElt[i]);
    }
    free (graph_bound->ghostEltPartIdx);
    free (graph_bound->nGhostEltPart);
    free (graph_bound->ghostEltPartElt);
    free (graph_bound->ghostEltPart2GhostElt);
    free (graph_bound);
  }
  return NULL;
}


/**
 * \brief Return the number of ghost elements for each partition
 *
 * This function returns the number of ghost elements for each partition
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  part           Partition number
 *
 * \return  Number of ghost element
 *
 */

int
PDM_graph_bound_n_ghost_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int               part
)
{

  if (part < 0 || part >= graph_bound->nPart) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_n_ghost_elt_get : "
            "part number '%d' is not available\n", part);
    abort();
  }

  return graph_bound->nGhostEltPart[part];

}


/**
 * \brief Return the number of local elements which touch the current ghost element
 *
 * This function returns the number of local elements which touch the current ghost element
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  part           Partition number
 *
 * \return  Number of local element which touch the current ghost element
 *
 */

int
PDM_graph_bound_ghost_elt_n_touch_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                part,
 int                ghostElt
)
{

  if (part < 0 || part >= graph_bound->nPart) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_ghost_elt_n_touch_elt_get : "
            "part number '%d' is not available\n", part);
    abort();
  }

  int nGhostElt = graph_bound->nGhostEltPart[part];
  int *ghostEltPartIdx = graph_bound->ghostEltPartIdx[part];

  if (ghostElt < 0 || ghostElt >= nGhostElt) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_ghost_elt_n_touch_elt_get : "
            "ghost element number '%d' is not available\n", part);
    abort();
  }

  return (ghostEltPartIdx[ghostElt+1] - ghostEltPartIdx[ghostElt]);

}


/**
 * \brief Return the list of local elements which touch the current ghost element
 *
 * This function returns the list of local elements which touch the current ghost element
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  part           Partition number
 *
 * \return  list of local elements which touch the current ghost element
 *
 */

int *
PDM_graph_bound_ghost_elt_touch_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                part,
 int                ghostElt
)
{

  if (part < 0 || part >= graph_bound->nPart) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_ghost_elt_touch_elt_get : "
            "part number '%d' is not available\n", part);
    abort();
  }

  int nGhostElt = graph_bound->nGhostEltPart[part];
  int *ghostEltPartIdx = graph_bound->ghostEltPartIdx[part];
  int *ghostEltPartElt = graph_bound->ghostEltPartElt[part];

  if (ghostElt < 0 || ghostElt >= nGhostElt) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_ghost_elt_n_touch_elt_get : "
            "ghost element number '%d' is not available\n", part);
    abort();
  }

  return ghostEltPartElt + ghostEltPartIdx[ghostElt];

}


/**
 * \brief Return the number of element to send to the specified processus
 *
 * This function returns the number of element to send to the specified processus
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  iProc          Processus to send
 *
 * \return  Number of element to send
 *
 */

int
PDM_graph_bound_n_send_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                iProc
)
{

  if (iProc < 0 || iProc >= graph_bound->lComm) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_n_send_elt_get : "
            "iProc number '%d' is not available\n", iProc);
    abort();
  }

  return graph_bound->sendEltIdx[iProc+1] - graph_bound->sendEltIdx[iProc];
}


/**
 * \brief Return the local element numbers in theirs partition to send to the specified processus
 *
 * This function returns the local element numbers in theirs partition to send to the specified processus
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  iProc          Processus to send
 *
 * \return  Local element numbers to send to this processus
 *
 */

int *
PDM_graph_bound_send_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                iProc
)
{

  if (iProc < 0 || iProc >= graph_bound->lComm) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_send_elt_get : "
            "iProc number '%d' is not available\n", iProc);
    abort();
  }

  return graph_bound->sendElt + graph_bound->sendEltIdx[iProc];
}


/**
 * \brief Return the partitions of elements to send to the specified processus
 *
 * This function returns the partitions of elements to send to the specified processus
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  iProc          Processus to send
 *
 * \return  Paritions of elements to send to this processus
 *
 */

int *
PDM_graph_bound_send_part_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                iProc
)
{

  if (iProc < 0 || iProc >= graph_bound->lComm) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_graph_bound_send_elt_get : "
            "iProc number '%d' is not available\n", iProc);
    abort();
  }

  return graph_bound->sendEltPart + graph_bound->sendEltIdx[iProc];
}


/**
 * \brief Dump the graph boundary structure
 *
 * This function dumps the graph boundary structure
 *
 * \param [in]  graph_bound    Boundary graph
 *
 */

void
PDM_graph_bound_dump
(
 PDM_graph_bound_t *graph_bound
)
{

  PDM_printf ("PDM_graph_bound :\n");
  PDM_printf ("  - Data location to send to the other processus\n");
  PDM_printf ("    - nSendElt %d\n",  graph_bound->nSendElt);

  for (int i = 0; i < graph_bound->lComm; i++) {

    PDM_printf ("    -  Element partitions for processus %d :", i);

    for (int j = graph_bound->sendEltIdx[i]; j < graph_bound->sendEltIdx[i+1]; j++) {
      PDM_printf (" %d", graph_bound->sendEltPart[j]);
    }
    PDM_printf ("\n");
    PDM_printf ("    -  Local element number in their partitions for processus %d :", i);

    for (int j = graph_bound->sendEltIdx[i]; j < graph_bound->sendEltIdx[i+1]; j++) {
      PDM_printf (" %d", graph_bound->sendElt[j]);
    }
    PDM_printf ("\n");

  }

  PDM_printf ("  - Data location received from the other processus\n");

  for (int i = 0; i < graph_bound->lComm; i++) {
    PDM_printf ("    - Number element received from processus %d : %d\n",
            i, graph_bound->ghostEltIdx[i+1] - graph_bound->ghostEltIdx[i]);
  }

  PDM_printf ("    -  Received element per part\n");

  for (int i = 0; i < graph_bound->nPart; i++) {

    PDM_printf ("      -  Part %d\n", i);

    for (int j = 0; j < graph_bound->nGhostEltPart[i]; j++) {

      PDM_printf ("        -  Local Connected elements to ghost element %d :", j);

      for (int k = graph_bound->ghostEltPartIdx[i][j]; k < graph_bound->ghostEltPartIdx[i][j+1]; k++) {
        PDM_printf (" %d", graph_bound->ghostEltPartElt[i][k]);
      }
      PDM_printf ("\n");

      PDM_printf ("        -  Indirection partition ghost element to processus ghost element : %d\n", graph_bound->ghostEltPart2GhostElt[i][j]);
    }
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
