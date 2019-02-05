
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"

#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_fortran_to_c_string.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"


/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/**
 * \def _PDM_part_MIN(a,b)
 * Computes the minimum of \a x and \a y. 
 *
 */

#define _PDM_part_MIN(a,b) ((a) > (b) ? (b) : (a))

/**
 * \def _PDM_part_MAX(a,b)
 * Computes the maximum of \a x and \a y. 
 *
 */

#define _PDM_part_MAX(a,b) ((a) < (b) ? (b) : (a))

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_pparts   = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * \param [inout]   c     Array sorted as a
 * 
 */

static void 
_quickSort_pdm_part_long_t
( 
 PDM_g_num_t a[], 
 int          l, 
 int          r, 
 int          c[]
)
{
  if (l < r) {
    int j = r+1;
    
    PDM_g_num_t t;
    int v; 
    PDM_g_num_t pivot = a[l];
    int i = l; 

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i]; 
      a[i] = a[j]; 
      a[j] = t;

      v    = c[i]; 
      c[i] = c[j];
      c[j] = v;
    }
    t    = a[l]; 
    a[l] = a[j]; 
    a[j] = t;

    v    = c[l]; 
    c[l] = c[j]; 
    c[j] = v;

    _quickSort_pdm_part_long_t(a, l  , j-1, c);
    _quickSort_pdm_part_long_t(a, j+1,   r, c);
  }
}

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * 
 */

static void 
_quickSort_int
( 
 int a[], 
 int l, 
 int r 
)
{
  if (l < r) {
    int j = r+1;
    int t; 
    int pivot = a[l];
    int i = l; 

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i]; 
      a[i] = a[j]; 
      a[j] = t;

    }
    t    = a[l]; 
    a[l] = a[j]; 
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}

/**
 *
 * \brief Quick sort 2
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * \param [inout]   c     Array sorted as a
 * 
 */

static void 
_quickSort_int2
( 
 int          a[], 
 int          l, 
 int          r, 
 int          c[]
)
{
  if (l < r) {
    int j = r+1;
    int  t, v; 
    int pivot = a[l];
    int i = l; 

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i]; 
      a[i] = a[j]; 
      a[j] = t;

      v    = c[i]; 
      c[i] = c[j];
      c[j] = v;
    }
    t    = a[l]; 
    a[l] = a[j]; 
    a[j] = t;

    v    = c[l]; 
    c[l] = c[j]; 
    c[j] = v;

    _quickSort_int2(a, l  , j-1, c);
    _quickSort_int2(a, j+1,   r, c);
  }
}

/**
 *
 * \brief Search the rank where element of distributed array is storage
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   id1          First index into array
 * \param [in]   id2          Last index into array
 *
 * \return       Rank where the element is stored 
 */


static int
 _search_rank
(
 PDM_g_num_t   elt,
 PDM_g_num_t  *array,
 int            id1,
 int            id2
)
{
  if (elt > array[id2]) {
    PDM_printf("PPART error : Element not in initial distributed array "
           PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n", 
           elt, array[id1], array[id2]);
    abort();
  } 

  if (elt < array[id1]) {
    PDM_printf("PPART error : Element not in initial distributed array "
           PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
           elt, array[id1], array[id2]);
    abort();
  } 

  if (id2 == id1 + 1) {
    return id1;
  }

  else {

    int midId = (id2 + id1) / 2;

    if (elt == array[id1])
      return id1;
    else if (elt == array[id2])
      return id2;
    else if (elt == array[midId])
      return midId;
    else if (elt < array[midId])
      return _search_rank(elt, array, id1, midId); 
    else if (elt > array[midId])
      return _search_rank(elt, array, midId, id2);
  } 
  return -1;
}


/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _PDM_part_t *
_get_from_id
(
 int  ppartId
)
{
  _PDM_part_t *part = (_PDM_part_t *) PDM_Handles_get (_pparts, ppartId);
    
  if (part == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "PPART error : Bad ppart identifier\n");
    exit(1);
  }

  return part;
}

/**
 *
 * \brief Call a couple MPI_Alltoall MPI_Alltoallv
 * 
 * \param [in]   sendBuff,            Sending buffer
 * \param [in]   sendBuffN,           Number of data to send to each process 
 *                                    (size : communicator size)
 * \param [in]   sendBuffIdx,         Index in sendBuff for each process 
 *                                    (size : communicator size)
 * \param [out]  recvBuff,            Receiving buffer 
 * \param [out]  recvBuffSize         Receiving buffer size  
 * \param [in]   exch_mpi_data_type   Data type to exchange
 * \param [in]   type_exch_size       Size of data type
 * \param [in]   comm                 Communicator
 *
 */

static void 
_alltoall
(
 void              *sendBuff,
 int               *sendBuffN,
 int               *sendBuffIdx,
 void             **recvBuff,
 int               *recvBuffN,
 int               *recvBuffIdx,
 PDM_MPI_Datatype   MPIDataType,
 size_t             MPIDataTypeSize,
 PDM_MPI_Comm       comm
)
{
  int nRank = 0;
  PDM_MPI_Comm_size(comm, &nRank);

  /* Get number data to receive from each process */

  PDM_MPI_Alltoall(sendBuffN, 
                   1, 
                   PDM_MPI_INT, 
                   recvBuffN, 
                   1, 
                   PDM_MPI_INT, 
                   comm);

  recvBuffIdx[0] = 0;
  for(int i = 0; i < nRank; i++) {
      recvBuffIdx[i+1] = recvBuffIdx[i] + recvBuffN[i];
  }

  *recvBuff = malloc(recvBuffIdx[nRank] * MPIDataTypeSize);

  /* Receive data from each process */

  PDM_MPI_Alltoallv(sendBuff, 
                    sendBuffN, 
                    sendBuffIdx, 
                    MPIDataType, 
                    *recvBuff, 
                    recvBuffN, 
                    recvBuffIdx,
                    MPIDataType, 
                    comm);

}

/**
 *
 * \brief Builds dual graph face cell connectivity
 * 
 * \param [inout] ppart       Ppart object
 *
 */

static void 
_dual_graph_from_face_cell
(
 _PDM_part_t *ppart
)
{
  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  /*
   * cellToSendN allocation
   */

  int *cellToSendN = (int *) malloc(nRank*sizeof(int));

  const int nData = 3; /* Number data to send */

  /*
   * Set cell list to send to each process
   */

  for (int i = 0; i < nRank; i++) {
    cellToSendN[i] = 0;
  }

  for (int i = 0; i < ppart->dNFace; i++) {
    PDM_g_num_t iCell1 = PDM_ABS (ppart->_dFaceCell[2*i    ]);
    PDM_g_num_t iCell2 = PDM_ABS (ppart->_dFaceCell[2*i + 1]);

    int irank1 = _search_rank(iCell1, ppart->dCellProc, 0, nRank);      
    cellToSendN[irank1] += nData;

    if (iCell2 > 0) {
      int irank2 = _search_rank(iCell2, ppart->dCellProc, 0, nRank);
      cellToSendN[irank2] += nData;
    }
  }
  
  /*
   * Create index aray
   */

  int *cellToSendIdx = (int *) malloc((nRank+1) * sizeof(int));

  cellToSendIdx[0] = 0;
  for (int i = 1; i < nRank + 1; i++) {
    cellToSendIdx[i] = cellToSendIdx[i-1] + cellToSendN[i-1];
    cellToSendN[i-1] = 0;
  }

  PDM_g_num_t *cellToSend = (PDM_g_num_t *) malloc(cellToSendIdx[nRank] * sizeof(PDM_g_num_t));

  /*
   * Stores pair of cells to send to the others processes
   */

  for (int i = 0; i < ppart->dNFace ; i++) {
    PDM_g_num_t iCell1 = PDM_ABS (ppart->_dFaceCell[2*i    ]);
    PDM_g_num_t iCell2 = PDM_ABS (ppart->_dFaceCell[2*i + 1]);
    int irank1 = _search_rank(iCell1, ppart->dCellProc, 0, nRank);

    int idx1             = cellToSendIdx[irank1] + cellToSendN[irank1];
    cellToSend[idx1  ]   = iCell1;
    cellToSend[idx1+1]   = iCell2;
    cellToSend[idx1+2]   = ppart->dFaceProc[myRank] + i;
    cellToSendN[irank1] += nData;

    if (iCell2 > 0) {
      int irank2 = _search_rank(iCell2, ppart->dCellProc, 0, nRank);
      int idx2             = cellToSendIdx[irank2] + cellToSendN[irank2];
      cellToSend[idx2  ]   = iCell2;
      cellToSend[idx2+1]   = iCell1;
      cellToSend[idx2+2]   = ppart->dFaceProc[myRank] + i;
      cellToSendN[irank2] += nData;
    }
  }
  
  /*
   * Receive pair of Cells from the others processes
   */

  int *cellToRecvN = (int *) malloc(nRank * sizeof(int));
  
  PDM_MPI_Alltoall(cellToSendN,
                   1, 
                   PDM_MPI_INT, 
                   cellToRecvN,
                   1, 
                   PDM_MPI_INT,
                   ppart->comm);

  int *cellToRecvIdx = (int *) malloc((nRank+1) * sizeof(int));

  cellToRecvIdx[0] = 0;
  for(int i = 1; i < (nRank+1); i++) {
    cellToRecvIdx[i] = cellToRecvIdx[i-1] + cellToRecvN[i-1];
  }

  PDM_g_num_t *cellToRecv = (PDM_g_num_t *) malloc(cellToRecvIdx[nRank]*sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(cellToSend,
                    cellToSendN,
                    cellToSendIdx, 
                    PDM__PDM_MPI_G_NUM, 
                    cellToRecv,
                    cellToRecvN,
                    cellToRecvIdx,
                    PDM__PDM_MPI_G_NUM, 
                    ppart->comm);

  int nRecvPair = cellToRecvIdx[nRank]/nData;

  /*
   * Free
   */

  free(cellToSendIdx);
  free(cellToSendN);
  free(cellToSend);
  free(cellToRecvIdx);
  free(cellToRecvN);

  cellToSendIdx  = NULL;
  cellToSendN    = NULL;
  cellToSend     = NULL;
  cellToRecvIdx  = NULL;
  cellToRecvN    = NULL;
  
  /*
   * Count neighbour cells for each cell
   */

  int *nNeighbour = (int *) malloc(ppart->dNCell * sizeof(int));

  int have_dCellFace = 0;
  if (ppart->_dCellFaceIdx != NULL)
    have_dCellFace = 1;

  int *dCellFaceN = NULL;
  if (!have_dCellFace) {
    dCellFaceN = (int *) malloc(ppart->dNCell * sizeof(int));
    for (int i = 0; i < ppart->dNCell; i++) {
      dCellFaceN[i] = 0;
    }
  }

  for (int i = 0; i < ppart->dNCell; i++) {
    nNeighbour[i]= 0; 
  }

  for (int i = 0; i < nRecvPair; i++) {
    PDM_g_num_t  gElt1 = cellToRecv[nData*i  ];                         // Get global numbering
    PDM_g_num_t  gElt2 = cellToRecv[nData*i+1];                         // Get global numbering
    PDM_g_num_t  _lElt1 = gElt1 - ppart->dCellProc[myRank];
    int          lElt1 = (int) _lElt1;      // Switch to local numbering

    if (gElt2 > 0) {
      nNeighbour[lElt1] += 1;
    }
    if (!have_dCellFace) {
      dCellFaceN[lElt1] += 1;
    }
  }

  /*
   * Allocate dual graph from neighbour cells for each cell
   */

  if (!have_dCellFace) {
    ppart->dCellFaceIdx  = (int *) malloc((1+ppart->dNCell) * sizeof(int));
    ppart->dCellFaceIdx[0] = 0;
    for (int i = 0; i < ppart->dNCell; i++) {
      ppart->dCellFaceIdx[i+1]     = ppart->dCellFaceIdx[i] + dCellFaceN[i];
    }
    ppart->dCellFace  = (PDM_g_num_t *) malloc(ppart->dCellFaceIdx[ppart->dNCell] * sizeof(PDM_g_num_t));
    for (int i = 0; i < ppart->dCellFaceIdx[ppart->dNCell]; i++) {
      ppart->dCellFace[i] = -1;
    }
    for (int i = 0; i < ppart->dNCell; i++) {
      dCellFaceN[i]  = 0;
    }
    ppart->_dCellFaceIdx = ppart->dCellFaceIdx;
    ppart->_dCellFace = ppart->dCellFace;
  }

  ppart->dDualGraphIdx = (PDM_g_num_t *) malloc((1+ppart->dNCell) * sizeof(PDM_g_num_t));
  ppart->dDualGraphIdx[0] = 0;
  for (int i = 0; i < ppart->dNCell; i++) {
    ppart->dDualGraphIdx[i+1] = ppart->dDualGraphIdx[i] + nNeighbour[i];
  }
  
  ppart->dDualGraph = (PDM_g_num_t *) malloc(ppart->dDualGraphIdx[ppart->dNCell] * 
                                              sizeof(PDM_g_num_t));

  for (int i = 0; i < ppart->dDualGraphIdx[ppart->dNCell]; i++) {
    ppart->dDualGraph[i] = -1;
  }

  for (int i = 0; i < ppart->dNCell; i++) {
    nNeighbour[i] = 0;
  }

  /*
   * Complete dual graph
   */

  for (int i = 0; i < nRecvPair; i++) {
    PDM_g_num_t  gCel1  = cellToRecv[nData*i];                      // global numbering
    PDM_g_num_t  _lCel1 = gCel1 - ppart->dCellProc[myRank];
    int           lCel1  = (int) _lCel1; // local numbering
    PDM_g_num_t  gCel2  = cellToRecv[nData*i+1];                    // global numbering
    PDM_g_num_t  gFace2 = cellToRecv[nData*i+2];                    // global numbering
    
    if (!have_dCellFace) {
      ppart->dCellFace[ppart->dCellFaceIdx[lCel1] + dCellFaceN[lCel1]] = gFace2;
      dCellFaceN[lCel1] += 1;
    }
        
    /*
     * Search if cel2 is already stored (To optimize for polyhedra (lot of neighbours) ?)
     */

    if (gCel2 > 0) {

      PDM_g_num_t k;
      for (k = ppart->dDualGraphIdx[lCel1]; k < ppart->dDualGraphIdx[lCel1] + nNeighbour[lCel1]; k++) {
        if (ppart->dDualGraph[k] == gCel2 - 1)
          break;
      }
      
      if (k == ppart->dDualGraphIdx[lCel1] + nNeighbour[lCel1]) {
        ppart->dDualGraph[ppart->dDualGraphIdx[lCel1] + nNeighbour[lCel1]] = gCel2 - 1;
        nNeighbour[lCel1] += 1;
      }
    }
  }

  /*
   * Compress dual graph
   */

  int k = 0;
  int k1 = 0;
  while (k < ppart->dDualGraphIdx[ppart->dNCell]) {
    if (ppart->dDualGraph[k] >= 0) {
      ppart->dDualGraph[k1] = ppart->dDualGraph[k];
      k++;
      k1++;
    }
    else
      k++;
  }
  
  /*
   * Reallocate to free unused memory
   */

  ppart->dDualGraph = realloc(ppart->dDualGraph, k1 * sizeof(PDM_g_num_t));

  ppart->dDualGraphIdx[0] = 0;
  for (int i = 1; i < ppart->dNCell + 1; i++)
    ppart->dDualGraphIdx[i] = ppart->dDualGraphIdx[i-1] + nNeighbour[i-1];

  /*
   * ppart->dCellFaceIdx is ppart->dDualGraphIdx
   */

  if (1 == 0) {
    if (!have_dCellFace) {
      PDM_printf("ppart->_dCellFace : \n");
      for (int i = 0; i < ppart->dNCell; i++) {
        for (int j = ppart->_dCellFaceIdx[i]; j < ppart->_dCellFaceIdx[i+1]; j++)
          PDM_printf(" "PDM_FMT_G_NUM, ppart->_dCellFace[j]);
        PDM_printf("\n");
      }
    }
  }

  free(cellToRecv);
  free(nNeighbour);
  if (!have_dCellFace) {
    free(dCellFaceN);
  }
}

/**
 *
 * \brief Builds dual graph from face cell connectivity
 * 
 * \param [inout] ppart       Ppart object
 *
 */

static void 
_dual_graph_from_cell_face
(
 _PDM_part_t *ppart
)
{
  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  /*
   * cellToSendN allocation
   */

  int *faceToSendN = (int *) malloc(nRank*sizeof(int));

  const int nData = 2; /* Number data to send */

  /*
   * Set cell list to send to each process
   */

  for (int i = 0; i < nRank; i++) {
    faceToSendN[i] = 0;
  }

  for (int i = 0; i < ppart->dNCell; i++) {
    for (int j = ppart->_dCellFaceIdx[i]; j < ppart->_dCellFaceIdx[i+1]; j++) {
      PDM_g_num_t iFace = PDM_ABS(ppart->_dCellFace[j]);

      int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);      
      faceToSendN[irank] += nData;
    }
  }
  
  /*
   * Create index aray
   */

  int *faceToSendIdx = (int *) malloc((nRank+1) * sizeof(int));

  faceToSendIdx[0] = 0;
  for (int i = 1; i < nRank + 1; i++) {
    faceToSendIdx[i] = faceToSendIdx[i-1] + faceToSendN[i-1];
    faceToSendN[i-1] = 0;
  }

  PDM_g_num_t *faceToSend = 
    (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));

  /*
   * Stores faces to send to the others processes
   */

  for (int i = 0; i < ppart->dNCell; i++) {
    for (int j = ppart->_dCellFaceIdx[i]; j < ppart->_dCellFaceIdx[i+1]; j++) {
      PDM_g_num_t iFace = PDM_ABS (ppart->_dCellFace[j]);

      int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);      
      int idx   = faceToSendIdx[irank] + faceToSendN[irank];
      faceToSend[idx  ]   = iFace;
      faceToSend[idx+1]   = ppart->dCellProc[myRank] + i;
      faceToSendN[irank] += nData;
    }
  }
  
  /*
   * Receive faces from the others processes
   */

  int *faceToRecvN = (int *) malloc(nRank * sizeof(int));
  
  PDM_MPI_Alltoall(faceToSendN,
                   1, 
                   PDM_MPI_INT, 
                   faceToRecvN,
                   1, 
                   PDM_MPI_INT,
                   ppart->comm);

  int *faceToRecvIdx = (int *) malloc((nRank+1) * sizeof(int));

  faceToRecvIdx[0] = 0;
  for(int i = 1; i < (nRank+1); i++) {
    faceToRecvIdx[i] = faceToRecvIdx[i-1] + faceToRecvN[i-1];
  }

  PDM_g_num_t *faceToRecv = 
    (PDM_g_num_t *) malloc(faceToRecvIdx[nRank]*sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(faceToSend,
                    faceToSendN,
                    faceToSendIdx, 
                    PDM__PDM_MPI_G_NUM, 
                    faceToRecv,
                    faceToRecvN,
                    faceToRecvIdx,
                    PDM__PDM_MPI_G_NUM, 
                    ppart->comm);

  int nRecvFace = faceToRecvIdx[nRank]/nData;

  /*
   * Rename
   */

  int *cellToSendIdx = faceToRecvIdx;
  int *cellToSendN   = faceToRecvN;

  PDM_g_num_t *cellToSend = 
    (PDM_g_num_t *) malloc(faceToRecvIdx[nRank]*sizeof(PDM_g_num_t));

  int         *cellToRecvIdx = faceToSendIdx;
  PDM_g_num_t *cellToRecv    = faceToSend;
  int         *cellToRecvN   = faceToSendN;

  /*
   * Buid ppart->dFaceCell
   */

  int have_dFaceCell = 0;

  if (ppart->dFaceCell != NULL) {
    have_dFaceCell = 1;
  }

  if (!have_dFaceCell) {
    ppart->dFaceCell = 
      (PDM_g_num_t *)  malloc((2*ppart->dNFace) * sizeof(PDM_g_num_t));
    for (int i = 0; i < 2*ppart->dNFace; i++) {
      ppart->dFaceCell[i] = 0;
    }

    for (int i = 0; i < nRecvFace; i++) {
      PDM_g_num_t  gFace = faceToRecv[nData*i  ];                    // Get global numbering
      PDM_g_num_t  gCell = faceToRecv[nData*i+1];                    // Get global numbering
      PDM_g_num_t  _lFace = gFace - ppart->dFaceProc[myRank]; 
      int          lFace = (int) _lFace; // Switch to local numbering

      if (ppart->dFaceCell[2*lFace] == 0)
        ppart->dFaceCell[2*lFace] = gCell;
      else if (ppart->dFaceCell[2*lFace + 1] == 0)
        ppart->dFaceCell[2*lFace + 1] = gCell;
      else {
        PDM_printf("PPART internal error : Face already defined in ppart->dFaceCell connectivity\n");
        exit(1);
      }
    }
    ppart->_dFaceCell = ppart->dFaceCell;
  }

  /*
   * Exchange cell neighbour
   */

  for (int i = 0; i < nRecvFace; i++) {
    PDM_g_num_t  gFace  = faceToRecv[nData*i  ];                    // Get global numbering
    PDM_g_num_t  gCell1 = faceToRecv[nData*i+1];                    // Get global numbering
    PDM_g_num_t _lFace = gFace - ppart->dFaceProc[myRank]; // Switch to local numbering
    int          lFace = (int) _lFace;
    PDM_g_num_t gCell2;

    if (ppart->dFaceCell[2*lFace] == gCell1)
      gCell2 = PDM_ABS (ppart->dFaceCell[2*lFace + 1]);
    else if (ppart->dFaceCell[2*lFace + 1] == gCell1)
      gCell2 = PDM_ABS (ppart->dFaceCell[2*lFace]);
    else {
      PDM_printf("PPART internal error : Problem in dual grah building "
              PDM_FMT_G_NUM" "
              PDM_FMT_G_NUM" "
              PDM_FMT_G_NUM" \n",
             ppart->dFaceCell[2*lFace ], ppart->dFaceCell[2*lFace + 1], gCell1);
      exit(1);
    }
    
    cellToSend[nData*i    ] = gCell1;
    cellToSend[nData*i + 1] = gCell2;
  }

  free(faceToRecv);

  PDM_MPI_Alltoallv(cellToSend,
                    cellToSendN,
                    cellToSendIdx, 
                    PDM__PDM_MPI_G_NUM, 
                    cellToRecv,
                    cellToRecvN,
                    cellToRecvIdx,
                    PDM__PDM_MPI_G_NUM, 
                    ppart->comm);

  int nRecvPair = cellToRecvIdx[nRank]/nData;

  /*
   * Allocate dual graph
   */

  ppart->dDualGraphIdx = (PDM_g_num_t *) malloc((1+ppart->dNCell) * sizeof(PDM_g_num_t));
  int *nNeighbour      = (int *) malloc(ppart->dNCell * sizeof(int));

  for (int i = 0; i < ppart->dNCell; i++) {
    nNeighbour[i] = 0;
  }

  ppart->dDualGraphIdx[0] = 0;
  
  ppart->dDualGraph = (PDM_g_num_t *) malloc(ppart->_dCellFaceIdx[ppart->dNCell] * 
                                              sizeof(PDM_g_num_t));

  for (int i = 0; i < ppart->_dCellFaceIdx[ppart->dNCell]; i++) {
    ppart->dDualGraph[i] = -1;
  }

  /*
   * Build dual graph
   */

  for (int i = 0; i < nRecvPair; i++) {
    PDM_g_num_t   gCel1  = cellToRecv[nData*i];              // global numbering
    PDM_g_num_t  _lCel1  = gCel1 - ppart->dCellProc[myRank]; // local numbering
    int           lCel1  = (int) (_lCel1);
    PDM_g_num_t   gCel2  = cellToRecv[nData*i+1];            // global numbering
    
    /*
     * Search if cel2 is already stored (To optimize for polyhedra (lot of neighbours) ?)
     */

    if (gCel2 > 0) {

      int k;
      for (k = ppart->_dCellFaceIdx[lCel1]; 
           k < ppart->_dCellFaceIdx[lCel1] + nNeighbour[lCel1]; k++) {
        if (ppart->dDualGraph[k] == gCel2 - 1)
          break;
      }
      
      if (k == ppart->_dCellFaceIdx[lCel1] + nNeighbour[lCel1]) {
        ppart->dDualGraph[ppart->_dCellFaceIdx[lCel1] + nNeighbour[lCel1]] = gCel2 - 1;
        nNeighbour[lCel1] += 1;
      }
    }
  }

  /*
   * Compress dual graph
   */

  int k = 0;
  int k1 = 0;
  while (k < ppart->_dCellFaceIdx[ppart->dNCell]) {
    if (ppart->dDualGraph[k] >= 0) {
      ppart->dDualGraph[k1] = ppart->dDualGraph[k];
      k++;
      k1++;
    }
    else
      k++;
  }
  
  /*
   * Reallocate to free unused memory
   */

  ppart->dDualGraph = realloc(ppart->dDualGraph, k1 * sizeof(PDM_g_num_t));

  ppart->dDualGraphIdx[0] = 0;
  for (int i = 1; i < ppart->dNCell + 1; i++)
    ppart->dDualGraphIdx[i] = ppart->dDualGraphIdx[i-1] + nNeighbour[i-1];

  /* Verifier tous les tableaux ..... */

  free(cellToRecv);
  free(cellToRecvIdx);
  free(cellToRecvN);
  free(cellToSend);
  free(cellToSendIdx);
  free(cellToSendN);
  free(nNeighbour);
}


/**
 *
 * \brief Splits the graph
 * 
 * \param [in]  ppart     ppart object
 * \param [out] cellPart  Cell partitioning (size : dNCell)
 *
 */

static void 
_split
(
 _PDM_part_t  *ppart,
 int          *cellPart
)
{
  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  for (int i = 0; i < ppart->dNCell; i++) {
    cellPart[i] = 0;
  }

  switch (ppart->split_method) {
  case PDM_PART_SPLIT_PARMETIS:
    {
#ifdef PDM_HAVE_PARMETIS

      /*
       * Define metis properties
       */

      int wgtflag    = 0;
      int numflag    = 0;        /* C or Fortran numbering (C = 0) */
      int edgecut;
      int ncon       = 1;
      
      double *ubvec = (double *) malloc(ncon * sizeof(double));
      for (int i = 0; i < ncon; i++) {
        ubvec[i] = 1.05;
      }

      double *tpwgts = (double *) malloc(ncon * ppart->tNPart * sizeof(double));

      for (int i = 0; i < ncon * ppart->tNPart; i++) {
        tpwgts[i] = (double) (1./ppart->tNPart);
      }

      /*
       * Call metis
       */

      PDM_g_num_t *_dCellProc = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));

      for (int i = 0; i < nRank + 1; i++) {
        _dCellProc[i] = ppart->dCellProc[i] - 1;
      }

      PDM_ParMETIS_V3_PartKway (_dCellProc,
                                ppart->dDualGraphIdx,
                                ppart->dDualGraph,
                                (int *) ppart->_dCellWeight,
                                NULL,
                                &wgtflag, 
                                &numflag, 
                                &ncon, 
                                &ppart->tNPart, 
                                tpwgts, 
                                ubvec, 
                                &edgecut, 
                                cellPart,
                                ppart->comm);
      
      free(ubvec);
      free(tpwgts);
      free(_dCellProc);

#else
      if(myRank == 0) {
        PDM_printf("PPART error : ParMETIS unavailable\n");
        exit(1);
      }
#endif
      break;
    }
  case PDM_PART_SPLIT_PTSCOTCH:
    {
#ifdef PDM_HAVE_PTSCOTCH
      int check = 0;
      int *edgeWeight = NULL;

      PDM_SCOTCH_dpart (ppart->dNCell,
                        ppart->dDualGraphIdx,
                        ppart->dDualGraph,
                        ppart->_dCellWeight,
                        edgeWeight,
                        check,        
                        ppart->comm,
                        ppart->tNPart,  
                        cellPart);

#else
      if(myRank == 0) {
        PDM_printf("PPART error : PT-Scotch unavailable\n");
        exit(1);
      }
#endif
      break;
    }
  case PDM_PART_SPLIT_HILBERT:
    {
      PDM_error(__FILE__, __LINE__, 0, "PPART error : Error in PT-Scotch graph check\n");
      exit(1);

      PDM_part_geom (PDM_PART_GEOM_HILBERT,
                     ppart->nPart,
                     ppart->comm,
                     ppart->dNCell,
                     ppart->_dCellFaceIdx,
                     ppart->_dCellFace,
                     ppart->_dCellWeight,
                     ppart->_dFaceVtxIdx,
                     ppart->dFaceProc,
                     ppart->_dFaceVtx,
                     ppart->_dVtxCoord,
                     ppart->dVtxProc,
                     cellPart);
      break;
    }
  default: 
    if(myRank == 0) {
      PDM_printf("PPART error : '%i' unknown partioning choice\n", ppart->split_method);
      exit(1);
    }
  }
}


/**
 *
 * \brief Distributes cell arrays 
 * 
 * \param [in]  ppart      ppart object
 * \param [in] cellPart  Cell partitioning (size : 3*dNCell)
 *
 */

static void
_distrib_cell
(
 _PDM_part_t  *ppart,
 int          *cellPart
)
{
  
  int myRank;
  int nRank;
  
  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);
  
  /* 
   *  For each cell faceToSend contains :
   *     - le numero de partition local
   *     - le numero global de l'element 
   *     - le nombre de faces
   *     - la liste des faces en numerotation globale
   * on envoie aussi en parallele le tableau send_nbfac qui donne 
   * le nombre de faces de chaque element envoye
   * ainsi que son numero de partition local (allant de 0 a nbMeshparProc)
   *
   */

  /* 1ere boucle pour compter le nombre d'elements qu'on envoie a chaque proc */

  int *faceToSendIdx = (int *) malloc((nRank + 1) * sizeof(int));
  for (int i = 0; i < nRank + 1; i++) {
    faceToSendIdx[i] = 0;
  }

  int nData = 3; /* Num cell, Partition locale, nbFac */
  if (ppart->_dCellTag != NULL)
    nData += 1;

  for (int i = 0; i < ppart->dNCell; i++) {
    int _part = (int) cellPart[i];
    int rankToSend = ppart->gPartTolProcPart[2*_part];
    int nbfac      = ppart->_dCellFaceIdx[i+1] - ppart->_dCellFaceIdx[i];
    faceToSendIdx[rankToSend+1]  += nData + nbfac; /* Num cell,
                                                      Partition locale
                                                      nbFac, 
                                                      liste des faces */ 
  }
  
  faceToSendIdx[0] = 0;
  for (int i = 0; i < nRank; i++) {
    faceToSendIdx[i+1] += faceToSendIdx[i] ;
  }

  int          *faceToSendN = (int *) malloc(nRank * sizeof(int));
  PDM_g_num_t *faceToSend  = 
    (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));

  for (int i = 0; i < nRank; i++)
    faceToSendN[i] = 0;

  /* 2nde boucle pour remplir le tableau a envoyer via alltoallv */

  for (int i = 0; i < ppart->dNCell; i++) {
    int _part = (int) cellPart[i];
    int rankToSend = ppart->gPartTolProcPart[2*_part    ];
    int lPart      = ppart->gPartTolProcPart[2*_part + 1];
    int nbfac      = ppart->_dCellFaceIdx[i+1] - ppart->_dCellFaceIdx[i];

    int place = faceToSendIdx[rankToSend] + faceToSendN[rankToSend];

    faceToSend[place++] = lPart;  /* Partition locale */
    faceToSendN[rankToSend] += 1;

    faceToSend[place++] = ppart->dCellProc[myRank] + i;  /* Numero global de l'elt*/
    faceToSendN[rankToSend] += 1;

    faceToSend[place++] = nbfac;  /* Nombre de faces*/
    faceToSendN[rankToSend] += 1;

    for (int j = ppart->_dCellFaceIdx[i]; j < ppart->_dCellFaceIdx[i+1]; j++) {
      faceToSend[place++] = PDM_ABS(ppart->_dCellFace[j]);   /* Numero global de ses faces */
      faceToSendN[rankToSend] += 1;
    }

    if (ppart->_dCellTag != NULL) {
      faceToSend[place++] = ppart->_dCellTag[i];  /* Tag des cellules si elles existent */
      faceToSendN[rankToSend] += 1;
    }
  }
  
  PDM_g_num_t *faceToRecv    = NULL;
  int          *faceToRecvN   = (int *) malloc(nRank * sizeof(int));
  int          *faceToRecvIdx = (int *) malloc((nRank + 1) * sizeof(int));
  
  _alltoall(faceToSend,
            faceToSendN,
            faceToSendIdx,
            (void **) &faceToRecv,
            faceToRecvN,
            faceToRecvIdx,
            PDM__PDM_MPI_G_NUM,
            sizeof(PDM_g_num_t),
            ppart->comm);

  int lFaceToRecv = faceToRecvIdx[nRank];

  free(faceToSend);
  free(faceToSendN);
  free(faceToSendIdx);
  free(faceToRecvN);
  free(faceToRecvIdx);

  /* Complete partitions */


  for (int i = 0; i < ppart->nPart; i++) {
    if (ppart->meshParts[i] == NULL)
      ppart->meshParts[i] = _part_create();
    _part_t *meshPart  = ppart->meshParts[i];
    meshPart->nVtx           = 0;
    meshPart->nFace          = 0;
    meshPart->nCell          = 0;
    meshPart->nFacePartBound = 0;
  }

  /* First loop for counting */

  int k = 0;
  while (k < lFaceToRecv) {

    _part_t *meshPart      = ppart->meshParts[faceToRecv[k++]];
    k += 1;
    int          nCellFace = (int) faceToRecv[k++];     

    k += nCellFace;
    if (ppart->_dCellTag != NULL)
      k += 1;

    meshPart->nCell += 1;
    meshPart->nFace += nCellFace;  /* Utilisation temporaire de nFace */
  }    

  /* Allocates arrays */

  for (int i = 0; i < ppart->nPart; i++) {

    _part_t *meshPart  = ppart->meshParts[i];

    meshPart->cellFaceIdx    = (int *)          malloc((meshPart->nCell + 1) * sizeof(int));
    meshPart->cellFaceIdx[0] = 0;
    meshPart->gCellFace      = (PDM_g_num_t *) malloc(meshPart->nFace * sizeof(PDM_g_num_t));
    meshPart->cellLNToGN     = (PDM_g_num_t *) malloc(meshPart->nCell * sizeof(PDM_g_num_t));
    if (ppart->_dCellTag != NULL)
      meshPart->cellTag      = (int *)          malloc(meshPart->nCell * sizeof(int));

    meshPart->nCell          = 0; /* reset temporary */

  }

  /* Second loop to complete arrays */

  k = 0;
  while (k < lFaceToRecv) {

    _part_t *meshPart  = ppart->meshParts[faceToRecv[k++]];

    PDM_g_num_t gNcell    =       faceToRecv[k++];     
    meshPart->cellLNToGN[meshPart->nCell] = gNcell;

    int          nCellFace = (int) faceToRecv[k++];     
    int idx = meshPart->cellFaceIdx[meshPart->nCell];
    meshPart->cellFaceIdx[meshPart->nCell + 1] = idx + nCellFace;
   
    for (int i = 0; i < nCellFace; i++)
      meshPart->gCellFace[idx + i] = faceToRecv[k++];

    if (ppart->_dCellTag != NULL) {
      int tag = (int) faceToRecv[k++];
      meshPart->cellTag[meshPart->nCell] = tag;
    }

    meshPart->nCell += 1;
  }    

  free(faceToRecv);

  /* Face local numbering */

  for (int i = 0; i < ppart->nPart; i++) {
    
    _part_t *meshPart  = ppart->meshParts[i];

    int *initialIdx     = (int *) malloc(meshPart->nFace * sizeof(int));
    meshPart->cellFace  = (int *) malloc(meshPart->nFace * sizeof(int));

    /* Map on gCellface */

    meshPart->faceLNToGN = meshPart->gCellFace;

    for (int k1 = 0; k1 < meshPart->nFace; k1++) {
      initialIdx[k1] = k1;
    }

    /* Sort faceLNToGN */

    if (1 == 0) {
      PDM_printf("meshPart->nFace 1 : %i\n", meshPart->nFace);
      PDM_printf("meshPart->faceLNToGN 1 : ");
      for (int i1 = 0; i1 < meshPart->nFace; i1++)
        PDM_printf(" "PDM_FMT_G_NUM, meshPart->faceLNToGN[i1]);
      PDM_printf("\n");
    }

    _quickSort_pdm_part_long_t(meshPart->faceLNToGN, /* tableau a trier */
                            0,                    /* premier elt */
                            meshPart->nFace - 1,  /* dernier elt */
                            initialIdx); 

    /* Remove duplicate faces and build local cell face connectivity*/
    
    int nDupl = 0;
    int k2 = 0;
    int kCompress = 0;

    while (k2 < meshPart->nFace) {
      PDM_g_num_t iface = meshPart->faceLNToGN[k2];
      meshPart->faceLNToGN[kCompress]   = iface;
      meshPart->cellFace[initialIdx[k2]] = kCompress + 1;
      k2 += 1;
      while (k2 < meshPart->nFace) {
        if (meshPart->faceLNToGN[k2] == iface) {
          meshPart->cellFace[initialIdx[k2]] = kCompress + 1;
          k2 += 1;
          nDupl += 1;
        }
        else
          break;
      }
      kCompress += 1;
    }

    meshPart->nFace = kCompress;

    meshPart->faceLNToGN = (PDM_g_num_t *) realloc(meshPart->faceLNToGN, 
                                                    meshPart->nFace * sizeof(PDM_g_num_t));

    if (1 == 0) {
      PDM_printf("meshPart->nCell : %i\n", meshPart->nCell);
      
      PDM_printf("meshPart->cellLNToGN : ");
      for (int i1 = 0; i1 < meshPart->nCell; i1++)
        PDM_printf(" "PDM_FMT_G_NUM, meshPart->cellLNToGN[i1]);
      PDM_printf("\n");
      
      PDM_printf("meshPart->cellFace : \n");
      for (int i1 = 0; i1 < meshPart->nCell; i1++) {
        for (int j = meshPart->cellFaceIdx[i1]; j < meshPart->cellFaceIdx[i1+1]; j++)
          PDM_printf(" %i", meshPart->cellFace[j]);
        PDM_printf("\n");
      }

      PDM_printf("meshPart->nFace : %i\n", meshPart->nFace);

      PDM_printf("meshPart->faceLNToGN : ");
      for (int i1 = 0; i1 < meshPart->nFace; i1++)
        PDM_printf(" "PDM_FMT_G_NUM, meshPart->faceLNToGN[i1]);
      PDM_printf("\n");
    }

    /* Free */

    free(initialIdx);
    meshPart->gCellFace = NULL;

  }

}


/**
 *
 * \brief Distributes face arrays 
 * 
 * \param [in]  ppart      ppart object
 *
 */

static void
_distrib_face
(
 _PDM_part_t      *ppart
)
{
  int myRank;
  int nRank;
  
  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  const int nData     = 1;
  int       nDataFace = 2; 
  if (ppart->_dFaceTag != NULL)
    nDataFace += 1;

  int          *faceToSendIdx = (int *) malloc((nRank + 1) * sizeof(int));
  int          *faceToSendN   = (int *) malloc(nRank * sizeof(int));
  PDM_g_num_t *faceToSend    = NULL;

  int          *requestedFaceN   = (int *) malloc(nRank * sizeof(int));
  int          *requestedFaceIdx = (int *) malloc((nRank + 1) * sizeof(int));

  for (int ipart = 0; ipart < ppart->mNPart; ipart++) {

    _part_t *meshPart  = NULL;
    int *allToallNToLN = NULL;

    for (int i = 0; i < nRank+1; i++)
      faceToSendIdx[i] = 0;
    
    for (int i = 0; i < nRank; i++)
      faceToSendN[i] = 0;

    faceToSend = NULL;

    if (ipart < ppart->nPart) {
  
      meshPart  = ppart->meshParts[ipart];
      allToallNToLN = (int *) malloc(meshPart->nFace * sizeof(int));

      /* 
       *  Processes exchange list of faces which they want receive information
       */

      for (int i = 0; i < meshPart->nFace; i++) {
        PDM_g_num_t iFace = meshPart->faceLNToGN[i];
        int          irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);
        faceToSendIdx[irank+1] += nData;
        
      }

      for (int i = 0; i < nRank; i++) {
        faceToSendIdx[i+1] += faceToSendIdx[i] ;
      }

      faceToSend = (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < meshPart->nFace; i++) {

        PDM_g_num_t iFace = meshPart->faceLNToGN[i];
        int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);

        int idx = faceToSendIdx[irank] + faceToSendN[irank];
        
        allToallNToLN[idx/nData] = i;

        faceToSend[idx++]   = iFace;        /* Face global numbering */
        faceToSendN[irank] += nData;
      }
    }

    PDM_g_num_t *requestedFace    = NULL;
    
    _alltoall(faceToSend,
              faceToSendN,
              faceToSendIdx,
              (void **) &requestedFace,
              requestedFaceN,
              requestedFaceIdx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    if (faceToSend != NULL)
      free(faceToSend);

    /* 
     *  Processes exchange information about requested faces
     *  For each face, information contains :
     *     - tag (if ppart->_dFaceTag != NULL)
     *     - Number of vertices
     *     - Vertices
     *
     */

    int *sFaceInfoIdx = faceToSendIdx; 
    int *sFaceInfoN   = faceToSendN;

    for (int i = 0; i < nRank+1; i++) {
      sFaceInfoIdx[i] = 0;
    }

    for (int i = 0; i < nRank; i++) {
      sFaceInfoN[i] = 0;
    }

    for (int i = 0; i < nRank; i++) {
      for (int k = requestedFaceIdx[i]; k < requestedFaceIdx[i+1]; k+=nData) {
        PDM_g_num_t gFace     = requestedFace[k];
        PDM_g_num_t _lFace    = gFace - ppart->dFaceProc[myRank];
        int          lFace     = (int) _lFace;
        int          nbVtxFace = (int) (ppart->_dFaceVtxIdx[lFace+1] 
                                      - ppart->_dFaceVtxIdx[lFace]);
        sFaceInfoIdx[i+1] += nDataFace + nbVtxFace;
      }
    } 

    for (int i = 0; i < nRank; i++) {
      sFaceInfoIdx[i+1] += sFaceInfoIdx[i] ;
    }

    PDM_g_num_t *sFaceInfo = (PDM_g_num_t *) 
      malloc(sFaceInfoIdx[nRank] * sizeof(PDM_g_num_t));

    for (int i = 0; i < nRank; i++) {
      for (int k = requestedFaceIdx[i]; k < requestedFaceIdx[i+1]; k+=nData) {
        PDM_g_num_t gFace     = requestedFace[k];
        PDM_g_num_t _lFace    = gFace - ppart->dFaceProc[myRank];
        int          lFace    = (int) _lFace;
        int          nbVtxFace = (int) (ppart->_dFaceVtxIdx[lFace+1] 
                                      - ppart->_dFaceVtxIdx[lFace]);

        int idx = sFaceInfoIdx[i] + sFaceInfoN[i]; 

        if (ppart->_dFaceTag != NULL) {
          sFaceInfo[idx++] = ppart->_dFaceTag[lFace];   /* Tag de la face */
          sFaceInfoN[i] += 1;
        }

        sFaceInfo[idx++] = nbVtxFace;                   /* Number of vertices */
        sFaceInfoN[i] += 1;
    
        for(int j = ppart->_dFaceVtxIdx[lFace]; j < ppart->_dFaceVtxIdx[lFace+1]; j++) {
          sFaceInfo[idx++] =  ppart->_dFaceVtx[j];  /*numero global du sommet qui compose la face*/
          sFaceInfoN[i] += 1;
        }
      }
    } 

    free(requestedFace);

    PDM_g_num_t *rFaceInfo    = NULL;
    int          *rFaceInfoN   = requestedFaceN;
    int          *rFaceInfoIdx = requestedFaceIdx;
  
    _alltoall(sFaceInfo,
              sFaceInfoN,
              sFaceInfoIdx,
              (void **) &rFaceInfo,
              rFaceInfoN,
              rFaceInfoIdx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    free(sFaceInfo);
    sFaceInfo = NULL;

    if (ipart < ppart->nPart) {

      /* Complete faceTag, faceVtxIdx gfaceVtx */ 

      if (ppart->_dFaceTag != NULL)
        meshPart->faceTag  = (int *) malloc(meshPart->nFace * sizeof(int)); 
      meshPart->faceVtxIdx = (int *) malloc((meshPart->nFace + 1) * sizeof(int)); 
  
      int k = 0;
      for (int i = 0; i < meshPart->nFace; i++) {
        if (ppart->_dFaceTag != NULL)
          meshPart->faceTag[allToallNToLN[i]] = (int) rFaceInfo[k++];

        int nVtx = (int) rFaceInfo[k++];
        meshPart->faceVtxIdx[allToallNToLN[i]+1] = nVtx;

        k += nVtx;
      }

      meshPart->faceVtxIdx[0] = 0;
      for (int i = 0; i < meshPart->nFace; i++)
        meshPart->faceVtxIdx[i+1] += meshPart->faceVtxIdx[i];
      
      meshPart->gFaceVtx = 
        (PDM_g_num_t *) malloc(meshPart->faceVtxIdx[meshPart->nFace] * sizeof(PDM_g_num_t)); 

      k = 0;
      for (int i = 0; i < meshPart->nFace; i++) {
        if (ppart->_dFaceTag != NULL)
          k += 1;

        int nVtx = (int) rFaceInfo[k++];
        int idx = meshPart->faceVtxIdx[allToallNToLN[i]];
        for (int j = 0; j < nVtx; j++)
          meshPart->gFaceVtx[idx + j] = rFaceInfo[k++];
      }

      if (rFaceInfo != NULL)
        free(rFaceInfo);
      rFaceInfo = NULL;

      /* Vertex local numbering vtxLNToGN */ 

      int *initialIdx   = (int *) malloc(meshPart->faceVtxIdx[meshPart->nFace] * sizeof(int));
      meshPart->faceVtx = (int *) malloc(meshPart->faceVtxIdx[meshPart->nFace] * sizeof(int));

      /* Map on gFaceVtx */

      meshPart->vtxLNToGN = meshPart->gFaceVtx;

      for (int k1 = 0; k1 <  meshPart->faceVtxIdx[meshPart->nFace]; k1++) {
        initialIdx[k1] = k1;
      }

      /* Sort faceLNToGN */

      _quickSort_pdm_part_long_t(meshPart->vtxLNToGN,                       /* Array to sort */
                                 0,                                         /* First face */
                                 meshPart->faceVtxIdx[meshPart->nFace] - 1, /* Latest face */
                                 initialIdx); 

      /* Remove duplicate Vertex and build local face vertex connectivity*/
    
      int nDupl = 0;
      int k2 = 0;
      int kCompress = 0;

      while (k2 < meshPart->faceVtxIdx[meshPart->nFace]) {
        PDM_g_num_t iVtx = meshPart->vtxLNToGN[k2];
        meshPart->vtxLNToGN[kCompress]   = iVtx;
        meshPart->faceVtx[initialIdx[k2]] = kCompress + 1;
        k2 += 1;
        while (k2 < meshPart->faceVtxIdx[meshPart->nFace]) {
          if (meshPart->vtxLNToGN[k2] == iVtx) {
            meshPart->faceVtx[initialIdx[k2]] = kCompress + 1;
            k2 += 1;
            nDupl += 1;
          }
          else
            break;
        }
        kCompress += 1;
      }


      meshPart->nVtx = kCompress;

      meshPart->vtxLNToGN = 
        (PDM_g_num_t *) realloc(meshPart->vtxLNToGN, meshPart->nVtx * sizeof(PDM_g_num_t));

      /* Free */

      free(initialIdx);
      meshPart->gFaceVtx = NULL;

      if (1 == 0) {
        PDM_printf("meshPart->nVtx 1 : %i\n", meshPart->nVtx);
        PDM_printf("meshPart->vtxLNToGN 1 : ");
        for (int i1 = 0; i1 < meshPart->nVtx; i1++)
          PDM_printf(" "PDM_FMT_G_NUM, meshPart->vtxLNToGN[i1]);
        PDM_printf("\n");
        
        PDM_printf("meshPart->faceVtx : \n");
        for (int i1 = 0; i1 < meshPart->nFace; i1++) {
          for (int j = meshPart->faceVtxIdx[i1]; j < meshPart->faceVtxIdx[i1+1]; j++)
            PDM_printf(" %i", meshPart->faceVtx[j]);
          PDM_printf("\n");
        }
      }

    }

    if (rFaceInfo != NULL)
      free(rFaceInfo);
    rFaceInfo = NULL;

    if (allToallNToLN != NULL)
      free(allToallNToLN);

  } /* For ipart */


  free(faceToSendN);
  free(faceToSendIdx);
  free(requestedFaceN);
  free(requestedFaceIdx);

}


/**
 *
 * \brief Distributes vertex arrays 
 * 
 * \param [in]  ppart      ppart object
 *
 */

static void
_distrib_vtx
(
 _PDM_part_t      *ppart
)
{

  const int nData    = 1;
  int nDataVtx = 0; 
  if (ppart->_dVtxTag != NULL)
    nDataVtx += 1;

  int myRank;
  int nRank;
  
  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  int          *vtxToSendIdx = (int *) malloc((nRank + 1) * sizeof(int));
  int          *vtxToSendN   = (int *) malloc(nRank * sizeof(int));
  PDM_g_num_t *vtxToSend    = NULL;

  int          *requestedVtxN   = (int *) malloc(nRank * sizeof(int));
  int          *requestedVtxIdx = (int *) malloc((nRank + 1) * sizeof(int));

  for (int ipart = 0; ipart < ppart->mNPart; ipart++) {

    _part_t *meshPart  = NULL;
    int *allToallNToLN = NULL;

    for (int i = 0; i < nRank+1; i++)
      vtxToSendIdx[i] = 0;
    
    for (int i = 0; i < nRank; i++)
      vtxToSendN[i] = 0;

    vtxToSend = NULL;

    if (ipart < ppart->nPart) {
  
      meshPart  = ppart->meshParts[ipart];
      allToallNToLN = (int *) malloc(meshPart->nVtx * sizeof(int));

      /* 
       *  Processes exchange list of vtxs which they want receive information
       */

      for (int i = 0; i < meshPart->nVtx; i++) {
        PDM_g_num_t iVtx = meshPart->vtxLNToGN[i];
        int irank = _search_rank(iVtx, ppart->dVtxProc, 0, nRank);
        vtxToSendIdx[irank+1] += nData;
      }

      for (int i = 0; i < nRank; i++) {
        vtxToSendIdx[i+1] += vtxToSendIdx[i] ;
      }

      vtxToSend = (PDM_g_num_t *) malloc(vtxToSendIdx[nRank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < meshPart->nVtx; i++) {

        PDM_g_num_t iVtx = meshPart->vtxLNToGN[i];
        int irank = _search_rank(iVtx, ppart->dVtxProc, 0, nRank);

        int idx = vtxToSendIdx[irank] + vtxToSendN[irank];
        
        allToallNToLN[idx/nData] = i;

        vtxToSend[idx++]   = iVtx;        /* Vtx global numbering */
        vtxToSendN[irank] += nData;
      }
    }

    PDM_g_num_t *requestedVtx    = NULL;
    
    _alltoall(vtxToSend,
              vtxToSendN,
              vtxToSendIdx,
              (void **) &requestedVtx,
              requestedVtxN,
              requestedVtxIdx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    if (vtxToSend != NULL)
      free(vtxToSend);

    /* 
     *  Processes exchange information about requested vtxs
     *  For each vtx, information contains :
     *     - Tag (if ppart->_dVtxTag != NULL)
     *     - Coordinates
     *
     */

    int *sVtxInfoIdx = vtxToSendIdx; 
    int *sVtxInfoN   = vtxToSendN;

    for (int i = 0; i < nRank+1; i++) {
      sVtxInfoIdx[i] = 0;
    }

    for (int i = 0; i < nRank; i++) {
      sVtxInfoN[i] = 0;
    }

    for (int i = 0; i < nRank; i++) {
      for (int k = requestedVtxIdx[i]; k < requestedVtxIdx[i+1]; k += nData) {
        sVtxInfoIdx[i+1] += nDataVtx * (int) sizeof(int) + 3 * (int) sizeof(double);
      }
    } 

    for (int i = 0; i < nRank; i++) {
      sVtxInfoIdx[i+1] += sVtxInfoIdx[i];
    }

    unsigned char *sVtxInfo = (unsigned char *) 
      malloc(sVtxInfoIdx[nRank] * sizeof(unsigned char));

    for (int i = 0; i < nRank; i++) {
      for (int k = requestedVtxIdx[i]; k < requestedVtxIdx[i+1]; k+=nData) {
        PDM_g_num_t gVtx     = requestedVtx[k];
        PDM_g_num_t _lVtx     = gVtx - ppart->dVtxProc[myRank];
        int          lVtx     = (int) _lVtx;

        int idx = sVtxInfoIdx[i] + sVtxInfoN[i]; 

        if (ppart->_dVtxTag != NULL) {
          int *_i_sVtxInfo = (int *) (sVtxInfo + idx);
          *_i_sVtxInfo     = ppart->_dVtxTag[lVtx];   /* Tag de la vtx */
          sVtxInfoN[i] += sizeof(int);
          idx += sizeof(int);
        }
    
        double *_d_sVtxInfo = (double *) (sVtxInfo + idx);
        _d_sVtxInfo[0] = ppart->_dVtxCoord[3*lVtx    ];  
        _d_sVtxInfo[1] = ppart->_dVtxCoord[3*lVtx + 1];  
        _d_sVtxInfo[2] = ppart->_dVtxCoord[3*lVtx + 2];  
        sVtxInfoN[i] += 3 * sizeof(double);

      }
    } 

    free(requestedVtx);

    unsigned char *rVtxInfo    = NULL;
    int           *rVtxInfoN   = requestedVtxN;
    int           *rVtxInfoIdx = requestedVtxIdx;
  
    _alltoall(sVtxInfo,
              sVtxInfoN,
              sVtxInfoIdx,
              (void **) &rVtxInfo,
              rVtxInfoN,
              rVtxInfoIdx,
              PDM_MPI_UNSIGNED_CHAR,
              sizeof(PDM_MPI_UNSIGNED_CHAR),
              ppart->comm);

    if (sVtxInfo != NULL)
      free(sVtxInfo);

    if (ipart < ppart->nPart) {

      /* Complete vtxTag, vtx */ 

      if (ppart->_dVtxTag != NULL)
        meshPart->vtxTag  = (int *) malloc(meshPart->nVtx * sizeof(int)); 
      meshPart->vtx = (double *) malloc(3 * meshPart->nVtx * sizeof(double)); 
  
      int k = 0;
      for (int i = 0; i < meshPart->nVtx; i++) {
        if (ppart->_dVtxTag != NULL) {
          int *_i_rVtxInfo = (int *) (rVtxInfo + k);
          meshPart->vtxTag[allToallNToLN[i]] = *_i_rVtxInfo;
          k += sizeof(int);
        }

        double *_d_rVtxInfo = (double *) (rVtxInfo + k);
        meshPart->vtx[3*allToallNToLN[i]    ] = _d_rVtxInfo[0];
        meshPart->vtx[3*allToallNToLN[i] + 1] = _d_rVtxInfo[1];
        meshPart->vtx[3*allToallNToLN[i] + 2] = _d_rVtxInfo[2];
        k += 3*sizeof(double);
      }
      
      if (1 == 0) {
        PDM_printf("meshPart->vtx : \n");
        for (int i1 = 0; i1 < meshPart->nVtx; i1++) {
          PDM_printf(" %12.5e %12.5e %12.5e", meshPart->vtx[3*i1 ], meshPart->vtx[3*i1+1], meshPart->vtx[3*i1+2]);
          PDM_printf("\n");
        }
      }
        
    } /* if ipart */

    if (allToallNToLN != NULL)
      free(allToallNToLN);

    if (rVtxInfo != NULL)
      free(rVtxInfo);
    rVtxInfo = NULL;

  } /* For ipart */
  
  free(vtxToSendN);
  free(vtxToSendIdx);
  free(requestedVtxN);
  free(requestedVtxIdx);
  
}


/**
 *
 * \brief Builds face-cell connectivity
 * 
 * \param [in]  ppart      ppart object
 *
 */

static void 
_build_faceCell
(
 _PDM_part_t  *ppart
)
{
  for (int ipart = 0; ipart < ppart->nPart; ipart++) {
    _part_t *meshPart  = ppart->meshParts[ipart];
  
    meshPart->faceCell = (int *) malloc(2*meshPart->nFace * sizeof(int));

    for (int i = 0; i < 2 * meshPart->nFace; i++)
      meshPart->faceCell[i] = 0;
  
    for (int i = 0; i < meshPart->nCell; i++) {
      for (int j = meshPart->cellFaceIdx[i]; j < meshPart->cellFaceIdx[i+1]; j++) {
        int idx = 2 * (PDM_ABS(meshPart->cellFace[j])-1);
        if (meshPart->faceCell[idx] == 0) 
          meshPart->faceCell[idx] = i + 1;
        else 
          meshPart->faceCell[idx + 1] = i + 1;
      }
    }
    if (1 == 0) {
      PDM_printf("meshPart->faceCell : \n");
      for (int i1 = 0; i1 < meshPart->nFace; i1++) {
        PDM_printf(" %i %i", meshPart->faceCell[2*i1],  meshPart->faceCell[2*i1+1]);
        PDM_printf("\n");
      }
    } 
  }
}


/**
 *
 * \brief Search partitioning boundary faces
 * 
 * \param [in]  ppart      ppart object
 *
 */

static void 
_search_part_bound_face
(
 _PDM_part_t *ppart
)
{
  const int nData  = 4;
  const int nData2 = 5;
 
  int myRank;
  int nRank;
  
  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  int          *faceToSendIdx = (int *) malloc((nRank + 1) * sizeof(int));
  int          *faceToSendN   = (int *) malloc(nRank * sizeof(int));
  PDM_g_num_t *faceToSend    = NULL;

  int          *requestedFaceN   = (int *) malloc(nRank * sizeof(int));
  int          *requestedFaceIdx = (int *) malloc((nRank + 1) * sizeof(int));

  int nDataPB = 6;

  ppart->dPartBound = (int *) malloc(nDataPB * ppart->dNFace * sizeof(int));
  for (int i = 0; i < nDataPB * ppart->dNFace; i++)
    ppart->dPartBound[i] = -1;

  /*
   * First loop on partitions to look for boundary faces
   */

  for (int ipart = 0; ipart < ppart->mNPart; ipart++) {

    _part_t *meshPart  = NULL;

    for (int i = 0; i < nRank+1; i++)
      faceToSendIdx[i] = 0;
    
    for (int i = 0; i < nRank; i++)
      faceToSendN[i] = 0;

    faceToSend = NULL;

    if (ipart < ppart->nPart) {
  
      meshPart  = ppart->meshParts[ipart];

      /* 
       *  Processes exchange list of faces which they want receive information
       */

      int nBoundFace = 0;

      for (int i = 0; i < meshPart->nFace; i++) {
        int icell2 = PDM_ABS (meshPart->faceCell[2*i + 1]);
        if (icell2 == 0) {
          PDM_g_num_t iFace = meshPart->faceLNToGN[i];
          int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);
          faceToSendIdx[irank+1] += nData;
          nBoundFace += 1;
        }
      }

      for (int i = 0; i < nRank; i++) {
        faceToSendIdx[i+1] += faceToSendIdx[i] ;
      }

      faceToSend = (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < meshPart->nFace; i++) {

        int icell2 = PDM_ABS (meshPart->faceCell[2*i + 1]);
        if (icell2 == 0) {
          PDM_g_num_t gFace = meshPart->faceLNToGN[i];
          int irank = _search_rank(gFace, ppart->dFaceProc, 0, nRank);

          int idx = faceToSendIdx[irank] + faceToSendN[irank];

          faceToSend[idx++]   = gFace;        /* Face global numbering */
          faceToSend[idx++]   = i+1;          /* Face local numbering  */
          faceToSend[idx++]   = myRank;       /* Rank                  */
          faceToSend[idx++]   = ipart;    /* Partition             */
          faceToSendN[irank] += nData;
        }
      }
    }

    PDM_g_num_t *requestedFace    = NULL;
    
    _alltoall(faceToSend,
              faceToSendN,
              faceToSendIdx,
              (void **) &requestedFace,
              requestedFaceN,
              requestedFaceIdx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    free(faceToSend);
    int nFace = requestedFaceIdx[nRank]/nData;


    int idx = 0;
    for(int i = 0; i < nFace; i++) {
      PDM_g_num_t  gFace     = requestedFace[idx++];
      PDM_g_num_t  _lFace    = gFace - ppart->dFaceProc[myRank];
      int          lFace     = (int) _lFace;
      int          lFaceRank = (int) requestedFace[idx++];
      int          faceRank  = (int) requestedFace[idx++];
      int          partition = (int) requestedFace[idx++];

      int idx2 = 0;
      if (ppart->dPartBound[nDataPB * lFace] != -1)
        idx2 += nDataPB/2;

      ppart->dPartBound[nDataPB * lFace + idx2++] = faceRank; 
      ppart->dPartBound[nDataPB * lFace + idx2++] = lFaceRank; 
      ppart->dPartBound[nDataPB * lFace + idx2  ] = partition; 

    }
 
    free(requestedFace);

  } /* mNPart */
  
  /* Exchange dPartBound */

  for (int i = 0; i < nRank+1; i++)
    faceToSendIdx[i] = 0;
    
  for (int i = 0; i < nRank; i++)
    faceToSendN[i] = 0;

  int idx = 0;
  for(int i = 0; i < ppart->dNFace; i++) {
    int faceRank1  = ppart->dPartBound[idx++]; 
    idx += 2;

    int faceRank2  = ppart->dPartBound[idx++]; 
    idx += 2;
      
    if ((faceRank1 > -1) && (faceRank2 > -1)) {
      faceToSendIdx[faceRank1 + 1] += nData2;
      faceToSendIdx[faceRank2 + 1] += nData2;
    }
  }

  for (int i = 0; i < nRank; i++)
    faceToSendIdx[i + 1] += faceToSendIdx[i]; 

  int *faceToSendInt = (int *) malloc(faceToSendIdx[nRank] * sizeof(int));
  
  idx = 0;
  for(int i = 0; i < ppart->dNFace; i++) {
    int faceRank1  = ppart->dPartBound[idx++]; 
    int lFaceRank1 = ppart->dPartBound[idx++]; 
    int partition1 = ppart->dPartBound[idx++]; 

    int faceRank2  = ppart->dPartBound[idx++]; 
    int lFaceRank2 = ppart->dPartBound[idx++]; 
    int partition2 = ppart->dPartBound[idx++]; 
      
    if ((faceRank1 > -1) && (faceRank2 > -1)) {
      int idx2 = faceToSendIdx[faceRank1] + faceToSendN[faceRank1]; 
      faceToSendInt[idx2++]   = lFaceRank1;
      faceToSendInt[idx2++]   = partition1;
      faceToSendInt[idx2++]   = faceRank2;
      faceToSendInt[idx2++]   = lFaceRank2;
      faceToSendInt[idx2++]   = partition2;
      faceToSendN[faceRank1] += nData2;

      int idx3 = faceToSendIdx[faceRank2] + faceToSendN[faceRank2]; 
      faceToSendInt[idx3++]   = lFaceRank2;
      faceToSendInt[idx3++]   = partition2;
      faceToSendInt[idx3++]   = faceRank1;
      faceToSendInt[idx3++]   = lFaceRank1;
      faceToSendInt[idx3++]   = partition1;
      faceToSendN[faceRank2] += nData2;
      
    }
  }

  int *requestedFaceInt = NULL;

  _alltoall(faceToSendInt,
            faceToSendN,
            faceToSendIdx,
            (void **) &requestedFaceInt,
            requestedFaceN,
            requestedFaceIdx,
            PDM_MPI_INT,
            sizeof(PDM_MPI_INT),
            ppart->comm);
  
  free(faceToSendInt);

  /* Complete facePartBound */

  int nFacePartBoundRank = requestedFaceIdx[nRank]/nData2;

  idx = 0;
  for (int i = 0; i < nFacePartBoundRank; i++) {
    idx += 1;
    int partition1 = requestedFaceInt[idx++]; 

    idx += 3;

    _part_t *meshPart  = ppart->meshParts[partition1];
    meshPart->nFacePartBound += 1;
  }

  int nDataFacePartBound = 4;

  for (int i = 0; i < ppart->nPart; i++) {
    _part_t *meshPart  = ppart->meshParts[i];
    meshPart->facePartBound = 
      (int *) malloc(nDataFacePartBound * meshPart->nFacePartBound * sizeof(int));
    meshPart->nFacePartBound = 0;

    meshPart->facePartBoundProcIdx = 
      (int *) malloc((nRank + 1) * sizeof(int));

    meshPart->facePartBoundPartIdx = 
      (int *) malloc((ppart->tNPart + 1) * sizeof(int));

    for (int j = 0; j < nRank + 1; j++) {
      meshPart->facePartBoundProcIdx[j] = 0;
    }
    for (int j = 0; j < ppart->tNPart + 1; j++) {
      meshPart->facePartBoundPartIdx[j] = 0;
    }
  }

  idx = 0;

  for (int i = 0; i < nFacePartBoundRank; i++) {
    int lFaceRank1 = requestedFaceInt[idx++]; 
    int partition1 = requestedFaceInt[idx++]; 

    int faceRank2  = requestedFaceInt[idx++]; 
    int lFaceRank2 = requestedFaceInt[idx++]; 
    int partition2 = requestedFaceInt[idx++]; 

    _part_t *meshPart  = ppart->meshParts[partition1];
    meshPart->facePartBoundProcIdx[faceRank2+1] += 1;
    meshPart->facePartBoundPartIdx[ppart->dPartProc[faceRank2] + partition2 + 1] += 1;
    meshPart->facePartBound[nDataFacePartBound * meshPart->nFacePartBound    ] = lFaceRank1;
    meshPart->facePartBound[nDataFacePartBound * meshPart->nFacePartBound + 1] = faceRank2;
    meshPart->facePartBound[nDataFacePartBound * meshPart->nFacePartBound + 2] = partition2 + 1;
    meshPart->facePartBound[nDataFacePartBound * meshPart->nFacePartBound + 3] = lFaceRank2;
    meshPart->nFacePartBound += 1;
  }

  for (int i = 0; i < ppart->nPart; i++) {
    _part_t *meshPart  = ppart->meshParts[i];
    for (int j = 1; j < nRank + 1; j++) {
      meshPart->facePartBoundProcIdx[j] = meshPart->facePartBoundProcIdx[j] + meshPart->facePartBoundProcIdx[j-1];
    }
    for (int j = 1; j <  ppart->tNPart + 1; j++) {
      meshPart->facePartBoundPartIdx[j] = meshPart->facePartBoundPartIdx[j] + meshPart->facePartBoundPartIdx[j-1];
    }
  }
   
  for (int i = 0; i < ppart->nPart; i++) {
    _part_t *meshPart  = ppart->meshParts[i];

    int *work_array  = (int *) malloc(meshPart->nFacePartBound * sizeof(int));
    PDM_g_num_t *work_array2;
    if (sizeof(PDM_g_num_t) == sizeof(int)) {
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif      
      work_array2 = (PDM_g_num_t *) work_array;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
    }
    else {
      work_array2  = (PDM_g_num_t *) malloc(meshPart->nFacePartBound * sizeof(PDM_g_num_t));
    }
    
    int *ind         = (int *) malloc(meshPart->nFacePartBound * sizeof(int));
    int *copyFacePartBound = (int *) malloc(nDataFacePartBound * meshPart->nFacePartBound * sizeof(int));

    /* Sort by procs */

    int k = 0;
    for (int j = 0; j < meshPart->nFacePartBound; j++) {
      ind[j] = k;
      k += 1;
      work_array[j] =  meshPart->facePartBound[nDataFacePartBound * j + 1];
    }

    _quickSort_int2(work_array,
                    0,
                    meshPart->nFacePartBound - 1,
                    ind);

    for (int j = 0; j < nDataFacePartBound * meshPart->nFacePartBound; j++)
      copyFacePartBound[j] =  meshPart->facePartBound[j];

    for (int j = 0; j <  meshPart->nFacePartBound; j++) {
      meshPart->facePartBound[nDataFacePartBound * j    ] = copyFacePartBound[nDataFacePartBound * ind[j]    ];
      meshPart->facePartBound[nDataFacePartBound * j + 1] = copyFacePartBound[nDataFacePartBound * ind[j] + 1];
      meshPart->facePartBound[nDataFacePartBound * j + 2] = copyFacePartBound[nDataFacePartBound * ind[j] + 2];
      meshPart->facePartBound[nDataFacePartBound * j + 3] = copyFacePartBound[nDataFacePartBound * ind[j] + 3];
    }
    
    /* Sort by part in procs */

    for (int j = 0; j < nDataFacePartBound * meshPart->nFacePartBound; j++)
      copyFacePartBound[j] =  meshPart->facePartBound[j];

    k = 0;
    for (int j = 0; j < meshPart->nFacePartBound; j++) {
      ind[j] = k;
      k += 1;
      work_array[j] =  meshPart->facePartBound[nDataFacePartBound * j + 2];
    }

    for (int j = 0; j < nRank; j++) {
      _quickSort_int2(work_array,
                      meshPart->facePartBoundProcIdx[j],
                      meshPart->facePartBoundProcIdx[j+1]-1,
                      ind);
    }

    for (int j = 0; j <  meshPart->nFacePartBound; j++) {
      meshPart->facePartBound[nDataFacePartBound * j    ] = copyFacePartBound[nDataFacePartBound * ind[j]    ];
      meshPart->facePartBound[nDataFacePartBound * j + 1] = copyFacePartBound[nDataFacePartBound * ind[j] + 1];
      meshPart->facePartBound[nDataFacePartBound * j + 2] = copyFacePartBound[nDataFacePartBound * ind[j] + 2];
      meshPart->facePartBound[nDataFacePartBound * j + 3] = copyFacePartBound[nDataFacePartBound * ind[j] + 3];
    }

    /* Sort by face absolute number in parts */

    for (int j = 0; j < nDataFacePartBound * meshPart->nFacePartBound; j++)
      copyFacePartBound[j] =  meshPart->facePartBound[j];

    k = 0;
    for (int j = 0; j < meshPart->nFacePartBound; j++) {
      ind[j] = k;
      k += 1;
      work_array2[j] =   meshPart->faceLNToGN[meshPart->facePartBound[nDataFacePartBound * j] - 1];
    }

    for (int j = 0; j < ppart->tNPart; j++) {
      _quickSort_pdm_part_long_t (work_array2,
                                  meshPart->facePartBoundPartIdx[j],
                                  meshPart->facePartBoundPartIdx[j+1]-1,
                                  ind);
    }

    for (int j = 0; j <  meshPart->nFacePartBound; j++) {
      meshPart->facePartBound[nDataFacePartBound * j    ] = copyFacePartBound[nDataFacePartBound * ind[j]    ];
      meshPart->facePartBound[nDataFacePartBound * j + 1] = copyFacePartBound[nDataFacePartBound * ind[j] + 1];
      meshPart->facePartBound[nDataFacePartBound * j + 2] = copyFacePartBound[nDataFacePartBound * ind[j] + 2];
      meshPart->facePartBound[nDataFacePartBound * j + 3] = copyFacePartBound[nDataFacePartBound * ind[j] + 3];
    }

    if (sizeof(PDM_g_num_t) != sizeof(int)) {
      free (work_array2);
    }
    free(work_array);
    free(ind);
    free(copyFacePartBound);
    
  }

  /* Sort by absolute face in parts */

  if (0 == 1) {
    for (int i = 0; i < ppart->nPart; i++) {
      _part_t *meshPart  = ppart->meshParts[i];
      PDM_printf("[%i] meshPart->nFacePartBound : %i\n",myRank, meshPart->nFacePartBound);
      PDM_printf("[%i] meshPart->facePartBound : \n", myRank);
      for (int i1 = 0; i1 < meshPart->nFacePartBound; i1++) {
        PDM_printf("[%i] %i %i %i %i", myRank, meshPart->facePartBound[4*i1    ], 
               meshPart->facePartBound[4*i1 + 1],
               meshPart->facePartBound[4*i1 + 2],
               meshPart->facePartBound[4*i1 + 3]);
        PDM_printf("\n");
      }
    }
  }

  free(requestedFaceInt);
  free(ppart->dPartBound);
  ppart->dPartBound = NULL;
  free(faceToSendIdx);
  free(faceToSendN);

  free(requestedFaceN);
  free(requestedFaceIdx);

}


/**
 *
 * \brief Distributes boundaries 
 * 
 * \param [in]  ppart      ppart object
 *
 */

static void
_distrib_face_groups
(
 _PDM_part_t      *ppart
)
{

  int myRank;
  int nRank;
    
  PDM_MPI_Comm_rank(ppart->comm, &myRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  int          *faceToSendIdx = (int *) malloc((nRank + 1) * sizeof(int));
  int          *faceToSendN   = (int *) malloc(nRank * sizeof(int));
  PDM_g_num_t  *faceToSend    = NULL;
  
  int          *requestedFaceN   = (int *) malloc(nRank * sizeof(int));
  int          *requestedFaceIdx = (int *) malloc((nRank + 1) * sizeof(int));
    
  PDM_g_num_t *dFaceGroupProc = (PDM_g_num_t *) malloc((nRank + 1) * sizeof(PDM_g_num_t)); 
  PDM_g_num_t *dFaceGroup     = (PDM_g_num_t *) malloc(ppart->dNFace * sizeof(PDM_g_num_t)); 

  for (int igroup = 0; igroup < ppart->nFaceGroup; igroup++) { 

    /* 
     *  Build dFaceGroupProc
     */

    PDM_g_num_t nFaceGroup = ppart->_dFaceGroupIdx[igroup+1] 
                           - ppart->_dFaceGroupIdx[igroup];

    PDM_MPI_Allgather(&nFaceGroup,
                      1,
                      PDM__PDM_MPI_G_NUM, 
                      (void *) (&dFaceGroupProc[1]), 
                      1, 
                      PDM__PDM_MPI_G_NUM, 
                      ppart->comm);

    dFaceGroupProc[0] = 0;
    for (int i = 1; i < nRank+1; i++) {
      dFaceGroupProc[i] = dFaceGroupProc[i] + dFaceGroupProc[i-1];
    }

    /* 
     *  Build dFaceGroup
     */

    const int nDataG = 2;
    
    for (int i = 0; i < nRank+1; i++)
      faceToSendIdx[i] = 0;
      
    for (int i = 0; i < nRank; i++)
      faceToSendN[i] = 0;

    for (int i = 0; i < ppart->dNFace; i++)
      dFaceGroup[i] = (PDM_g_num_t) -1;
    
    for (int i = ppart->_dFaceGroupIdx[igroup]; 
             i < ppart->_dFaceGroupIdx[igroup+1]; 
             i++) {
      PDM_g_num_t iFace =  ppart->_dFaceGroup[i];
      int irank = _search_rank(iFace,  ppart->dFaceProc, 0, nRank);
      faceToSendIdx[irank+1] += nDataG;
    } 
    
    for (int i = 0; i < nRank; i++) {
      faceToSendIdx[i+1] += faceToSendIdx[i] ;
    }

    faceToSend = (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));
    for (int i = ppart->_dFaceGroupIdx[igroup]; 
         i < ppart->_dFaceGroupIdx[igroup+1]; 
         i++) {
      PDM_g_num_t iFace = ppart->_dFaceGroup[i];
      int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);
      int idx = faceToSendIdx[irank] + faceToSendN[irank];
      faceToSend[idx++] = dFaceGroupProc[myRank] + (PDM_g_num_t) (i - ppart->_dFaceGroupIdx[igroup] + 1);
      faceToSend[idx++] = iFace;
      faceToSendN[irank] += nDataG;
    } 

    PDM_g_num_t *requestedFace = NULL;

    _alltoall(faceToSend,
              faceToSendN,
              faceToSendIdx,
              (void **) &requestedFace,
              requestedFaceN,
              requestedFaceIdx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    free(faceToSend);
    faceToSend = NULL;

    int idx = 0;
    for (int i = 0; i < requestedFaceIdx[nRank]/nDataG; i++) {
      PDM_g_num_t iFace      = requestedFace[idx++];
      PDM_g_num_t gFaceGroup = requestedFace[idx++];
      PDM_g_num_t _lFace = gFaceGroup - ppart->dFaceProc[myRank];
      int          lFace = (int) _lFace;
      dFaceGroup[lFace] = (PDM_g_num_t) iFace;
    }

    /* 
     *  As distributes_faces
     */

    const int nData     = 1;
    const int nDataFace = 1; 
     
    for (int ipart = 0; ipart < ppart->mNPart; ipart++) {
      
      _part_t *meshPart  = NULL;
      int *allToallNToLN = NULL;
      
      for (int i = 0; i < nRank+1; i++)
        faceToSendIdx[i] = 0;
      
      for (int i = 0; i < nRank; i++)
        faceToSendN[i] = 0;
      
      faceToSend = NULL;
      
      if (ipart < ppart->nPart) {
        
        meshPart  = ppart->meshParts[ipart];
        allToallNToLN = (int *) malloc(meshPart->nFace * sizeof(int));
        
        /* 
         *  Processes exchange list of faces which they want receive information
         */
        
        for (int i = 0; i < meshPart->nFace; i++) {
          PDM_g_num_t iFace = meshPart->faceLNToGN[i];
          int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);
          faceToSendIdx[irank+1] += nData;
        }
        
        for (int i = 0; i < nRank; i++) {
          faceToSendIdx[i+1] += faceToSendIdx[i] ;
        }
        
        faceToSend = (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));
        
        for (int i = 0; i < meshPart->nFace; i++) {
          
          PDM_g_num_t iFace = meshPart->faceLNToGN[i];
          int irank = _search_rank(iFace, ppart->dFaceProc, 0, nRank);
          
          idx = faceToSendIdx[irank] + faceToSendN[irank];
          
          allToallNToLN[idx/nData] = i;
          
          faceToSend[idx++]   = iFace;        /* Face global numbering */
          faceToSendN[irank] += nData;
        }
      }
      
      if (requestedFace != NULL)
        free(requestedFace);
      requestedFace    = NULL;
      
      PDM_g_num_t *requestedFace2 = NULL;
      
      _alltoall(faceToSend,
                faceToSendN,
                faceToSendIdx,
                (void **) &requestedFace2,
                requestedFaceN,
                requestedFaceIdx,
                PDM__PDM_MPI_G_NUM,
                sizeof(PDM_g_num_t),
                ppart->comm);
      
      if (faceToSend != NULL)
        free(faceToSend);
      
      /* 
       *  Processes exchange information about requested faces
       *  For each face, information contains :
       *     - Face global number in the current group
       *
       */
      
      int *sFaceInfoIdx = faceToSendIdx; 
      int *sFaceInfoN   = faceToSendN;
      
      for (int i = 0; i < nRank+1; i++) {
        sFaceInfoIdx[i] = 0;
      }
      
      for (int i = 0; i < nRank; i++) {
        sFaceInfoN[i] = 0;
      }
      
      for (int i = 0; i < nRank; i++) {
        for (int k = requestedFaceIdx[i]; k < requestedFaceIdx[i+1]; k+=nData) {
          sFaceInfoIdx[i+1] += nDataFace;
        }
      } 
      
      for (int i = 0; i < nRank; i++) {
        sFaceInfoIdx[i+1] += sFaceInfoIdx[i] ;
      }
      
      PDM_g_num_t *sFaceInfo = (PDM_g_num_t *) 
        malloc(sFaceInfoIdx[nRank] * sizeof(PDM_g_num_t));
      
      for (int i = 0; i < nRank; i++) {
        for (int k = requestedFaceIdx[i]; k < requestedFaceIdx[i+1]; k+=nData) {
          PDM_g_num_t gFace     = requestedFace2[k];
          PDM_g_num_t _lFace    = gFace - ppart->dFaceProc[myRank];
          int          lFace     = (int) _lFace;
          
          idx = sFaceInfoIdx[i] + sFaceInfoN[i]; 
          
          sFaceInfo[idx++] = dFaceGroup[lFace];                   /* Number of vertices */
          sFaceInfoN[i] += 1;
          
        }
      } 
      
      free(requestedFace2);
      
      PDM_g_num_t *rFaceInfo    = NULL;
      int          *rFaceInfoN   = requestedFaceN;
      int          *rFaceInfoIdx = requestedFaceIdx;
      
      _alltoall(sFaceInfo,
                sFaceInfoN,
                sFaceInfoIdx,
                (void **) &rFaceInfo,
                rFaceInfoN,
                rFaceInfoIdx,
                PDM__PDM_MPI_G_NUM,
                sizeof(PDM_g_num_t),
                ppart->comm);

      free(sFaceInfo);
      sFaceInfo = NULL;

      if (ipart < ppart->nPart) {
        
        /* Complete faceGroupIdx faceGroupeFace */ 
        
        if (igroup == 0) {
          meshPart->faceGroupIdx = (int *) malloc((ppart->nFaceGroup+1) * sizeof(int));
          for (int i = 0; i < ppart->nFaceGroup+1; i++)
            meshPart->faceGroupIdx[i] = 0;
        } 
        
        meshPart->faceGroupIdx[igroup+1] = meshPart->faceGroupIdx[igroup]; 
        for (int i = 0; i < rFaceInfoIdx[nRank]; i++)
          if (rFaceInfo[i] != -1)
            meshPart->faceGroupIdx[igroup+1] += 1;

        if (igroup == 0) {
          meshPart->faceGroup = (int *) malloc(meshPart->faceGroupIdx[igroup+1] * sizeof(int));
          meshPart->faceGroupLNToGN = 
            (PDM_g_num_t *) malloc(meshPart->faceGroupIdx[igroup+1] 
                                    * sizeof(PDM_g_num_t));
        }
        else {
          meshPart->faceGroup = 
            (int *) realloc(meshPart->faceGroup, meshPart->faceGroupIdx[igroup+1] * sizeof(int));
          meshPart->faceGroupLNToGN = 
            (PDM_g_num_t *) realloc(meshPart->faceGroupLNToGN, meshPart->faceGroupIdx[igroup+1] 
                                     * sizeof(PDM_g_num_t));
        }
          
        idx = meshPart->faceGroupIdx[igroup];
        for (int i = 0; i < rFaceInfoIdx[nRank]; i++) {
          if (rFaceInfo[i] != -1) {
            meshPart->faceGroup[idx] = allToallNToLN[i]+1;
            meshPart->faceGroupLNToGN[idx] = rFaceInfo[i];
            idx += 1;
          }

        }
      }
        
      if (rFaceInfo != NULL)
        free(rFaceInfo);

      if (allToallNToLN != NULL)
        free(allToallNToLN);

    }  /* For ipart */

  } /* For nFaceGroup */

  if (1 == 0) {
    for (int ipart = 0; ipart < ppart->nPart; ipart++) {
      
      _part_t *meshPart  = ppart->meshParts[ipart];
      
      PDM_printf("meshPart->nFaceGroup : %i\n",  ppart->nFaceGroup);
      PDM_printf("meshPart->faceGroup : \n");
      for (int i1 = 0; i1 < ppart->nFaceGroup; i1++) {
        for (int i2 = meshPart->faceGroupIdx[i1]; i2 < meshPart->faceGroupIdx[i1+1]; i2++) 
          PDM_printf(" %i", meshPart->faceGroup[i2]);
        PDM_printf(" --\n");
      }
      PDM_printf("meshPart->faceGroupLNToGN : \n");
      for (int i1 = 0; i1 < ppart->nFaceGroup; i1++) {
        for (int i2 = meshPart->faceGroupIdx[i1]; i2 < meshPart->faceGroupIdx[i1+1]; i2++) 
          PDM_printf(" "PDM_FMT_G_NUM, meshPart->faceGroupLNToGN[i2]);
        PDM_printf(" --\n");
      }
    }
  }

  free(faceToSendN);
  free(faceToSendIdx);
  free(requestedFaceN);
  free(requestedFaceIdx);
  free(dFaceGroupProc);
  free(dFaceGroup);
}

/**
 *
 * \brief Free partition
 *
 * \param [in]   part      partition
 *
 */

static void 
_part_free
(
 _part_t *part
)
{
  if (part->cellFaceIdx != NULL)
    free(part->cellFaceIdx);
  part->cellFaceIdx = NULL;    

  if (part->gCellFace != NULL)
    free(part->gCellFace);
  part->gCellFace = NULL;    

  if (part->cellFace != NULL)
    free(part->cellFace);
  part->cellFace = NULL;    

  if (part->cellLNToGN != NULL)
    free(part->cellLNToGN);
  part->cellLNToGN = NULL;    

  if (part->cellTag != NULL)
    free(part->cellTag);
  part->cellTag = NULL;    

  if (part->faceCell != NULL)
    free(part->faceCell);
  part->faceCell = NULL;    

  if (part->faceVtxIdx != NULL)
    free(part->faceVtxIdx);
  part->faceVtxIdx = NULL;    

  if (part->gFaceVtx != NULL)
    free(part->gFaceVtx);
  part->gFaceVtx = NULL;    

  if (part->faceVtx != NULL)
    free(part->faceVtx);
  part->faceVtx = NULL;    

  if (part->faceLNToGN != NULL)
    free(part->faceLNToGN);
  part->faceLNToGN = NULL;    

  if (part->faceTag != NULL)
    free(part->faceTag);
  part->faceTag = NULL;    

  if (part->facePartBoundProcIdx != NULL)
    free(part->facePartBoundProcIdx);
  part->facePartBoundProcIdx = NULL;    

  if (part->facePartBoundPartIdx != NULL)
    free(part->facePartBoundPartIdx);
  part->facePartBoundPartIdx = NULL;    

  if (part->facePartBound != NULL)
    free(part->facePartBound);
  part->facePartBound = NULL;    

  if (part->faceGroupIdx != NULL)
    free(part->faceGroupIdx);
  part->faceGroupIdx = NULL;    

  if (part->faceGroup != NULL)
    free(part->faceGroup);
  part->faceGroup = NULL;    

  if (part->faceGroupLNToGN != NULL)
    free(part->faceGroupLNToGN);
  part->faceGroupLNToGN = NULL;    

  if (part->vtx != NULL)
    free(part->vtx);
  part->vtx = NULL;    

  if (part->vtxLNToGN != NULL)
    free(part->vtxLNToGN);
  part->vtxLNToGN = NULL;    

  if (part->vtxTag != NULL)
    free(part->vtxTag);
  part->vtxTag = NULL;    

  free(part);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a initial partitioning
 *
 *  Build a initial partitioning from :
 *      - Cell block distribution with implicit global numbering 
 *         (the first cell is the first cell of the first process and 
 *          the latest cell is the latest cell of the latest process)   
 *      - Face block distribution with implicit global numbering 
 *      - Vertex block distribution with implicit global numbering 
 *  To repart an existing partition use \ref PDM_part_repart function
 * 
 * \param [out]  ppartId        ppart identifier
 * \param [in]   comm           Communicator
 * \param [in]   split_method   Split method
 * \param [in]   renum_cell_method Cell renumbering method
 * \param [in]   renum_face_method Cell renumbering method
 * \param [in]   renum_properties_cell  For cache blocking [ nCellPerCacheWanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking 
 * \param [in]   renum_face_method Cell renumbering method
 * \param [in]   renum_properties_face  NOT USE
 * \param [in]   nPart          Number of partition to build on this process
 * \param [in]   dNCell         Number of distributed cells
 * \param [in]   dNFace         Number of distributed faces
 * \param [in]   dNVtx          Number of distributed vertices
 * \param [in]   nFaceGroup     Number of face groups             
 * \param [in]   dCellFaceIdx   Distributed cell face connectivity index or NULL
 *                              (size : dNCell + 1, numbering : 0 to n-1)
 * \param [in]   dCellFace      Distributed cell face connectivity or NULL
 *                              (size : dFaceVtxIdx[dNCell], numbering : 1 to n)
 * \param [in]   dCellTag       Cell tag (size : nCell) or NULL
 * \param [in]   dCellWeight    Cell weight (size : nCell) or NULL
 * \param [in]   dCellPart      Distributed cell partitioning 
 *                              (size = dNCell) or NULL (No partitioning if != NULL)
 * \param [in]   dFaceCell      Distributed face cell connectivity or NULL
 *                              (size : 2 * dNFace, numbering : 1 to n)
 * \param [in]   dFaceVtxIdx    Distributed face to vertex connectivity index 
 *                              (size : dNFace + 1, numbering : 0 to n-1)
 * \param [in]   dFaceVtx       Distributed face to vertex connectivity 
 *                              (size : dFaceVtxIdx[dNFace], numbering : 1 to n)
 * \param [in]   dFaceTag       Distributed face tag (size : dNFace)
 *                              or NULL
 * \param [in]   dVtxCoord      Distributed vertex coordinates 
 *                              (size : 3*dNVtx)
 * \param [in]   dVtxTag        Distributed vertex tag (size : dNVtx) or NULL
 * \param [in]   dFaceGroupIdx  Index of distributed faces list of each group 
 *                              (size = nFaceGroup + 1) or NULL
 * \param [in]   dFaceGroup     distributed faces list of each group
 *                              (size = dFaceGroup[dFaceGroupIdx[nFaceGroup]], numbering : 1 to n) 
 *                              or NULL
 *
 */

void 
PDM_part_create
(
 int                         *ppartId,
 const PDM_MPI_Comm           comm,
 const PDM_part_split_t       split_method,
 const char                  *renum_cell_method,
 const char                  *renum_face_method,
 const int                    nPropertyCell,
 const int*                   renum_properties_cell,
 const int                    nPropertyFace,
 const int*                   renum_properties_face,
 const int                    nPart,
 const int                    dNCell,
 const int                    dNFace,
 const int                    dNVtx,
 const int                    nFaceGroup,
 const int                   *dCellFaceIdx,
 const PDM_g_num_t           *dCellFace,
 const int                   *dCellTag,
 const int                   *dCellWeight,
 const int                    have_dCellPart,
       int                   *dCellPart,
 const PDM_g_num_t           *dFaceCell,
 const int                   *dFaceVtxIdx,
 const PDM_g_num_t           *dFaceVtx,
 const int                   *dFaceTag,
 const double                *dVtxCoord,
 const int                   *dVtxTag,
 const int                   *dFaceGroupIdx,
 const PDM_g_num_t           *dFaceGroup
)
{
  int myRank;
  int nRank;
  
  PDM_MPI_Comm_rank(comm, &myRank);
  PDM_MPI_Comm_size(comm, &nRank);


  PDM_part_renum_method_load_local();

  /*
   * Search a ppart free id
   */

  if (_pparts == NULL) {
    _pparts = PDM_Handles_create (4); 
  }

  _PDM_part_t *ppart = (_PDM_part_t *) malloc(sizeof(_PDM_part_t));
  
  *ppartId = PDM_Handles_store (_pparts, ppart);

  /*
   * Build ppart structure
   */

  ppart->timer = PDM_timer_create();
  for (int i = 0; i < 4; i++) {
    ppart->times_elapsed[i] = 0.;
    ppart->times_cpu[i] = 0.;
    ppart->times_cpu_u[i] = 0.;
    ppart->times_cpu_s[i] = 0.;
  }
  PDM_timer_resume(ppart->timer);

  /* Local dimensions */

  ppart->dNVtx      = dNVtx;
  ppart->dNCell     = dNCell;
  ppart->dNFace     = dNFace;
  ppart->nFaceGroup = nFaceGroup;

  /* Cell definitions */

  ppart->_dCellFaceIdx = dCellFaceIdx;
  ppart->_dCellFace    = dCellFace;
  ppart->_dCellTag     = dCellTag;
  ppart->_dCellWeight  = dCellWeight;
  ppart->_dCellPart    = dCellPart;
  ppart->dCellFaceIdx  = NULL;
  ppart->dCellFace     = NULL;
  ppart->dFaceCell     = NULL;
  
  /* Set up for renumbering */
  ppart->nPropertyCell          = nPropertyCell;
  ppart->renum_properties_cell  = renum_properties_cell;
  ppart->nPropertyFace          = nPropertyFace;
  ppart->renum_properties_face  = renum_properties_face;

  ppart->dCellProc = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dNCell = (PDM_g_num_t) dNCell;
  PDM_MPI_Allgather((void *) &_dNCell,
                    1,
                    PDM__PDM_MPI_G_NUM, 
                    (void *) (&ppart->dCellProc[1]), 
                    1, 
                    PDM__PDM_MPI_G_NUM, 
                    comm);

  ppart->dCellProc[0] = 1;

  for (int i = 1; i < nRank+1; i++) {
    ppart->dCellProc[i] +=  ppart->dCellProc[i-1];
  }

  if (1 == 0) {
    PDM_printf("ppart->dCellProc : "PDM_FMT_G_NUM,  ppart->dCellProc[0]);
    for (int i = 1; i < nRank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, ppart->dCellProc[i]);
    }
    PDM_printf("\n");
  }

  /* Face definitions */

  ppart->_dFaceTag     = dFaceTag;
  ppart->_dFaceCell    = dFaceCell;
  ppart->_dFaceVtxIdx  = dFaceVtxIdx;
  ppart->_dFaceVtx     = dFaceVtx;

  ppart->dFaceProc = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  int *dNFaceProc = (int *) malloc((nRank) * sizeof(int));

  PDM_MPI_Allgather((void *) &dNFace,
                1,
                PDM_MPI_INT, 
                (void *) dNFaceProc, 
                1, 
                PDM_MPI_INT, 
                comm);
  ppart->dFaceProc[0] = 1;
  for (int i = 1; i < nRank+1; i++) {
    ppart->dFaceProc[i] = (PDM_g_num_t) dNFaceProc[i-1] + ppart->dFaceProc[i-1];
  }

  free(dNFaceProc);  

  /* Vertex definitions */

  ppart->_dVtxCoord    = dVtxCoord;
  ppart->_dVtxTag      = dVtxTag;

  ppart->dVtxProc = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  int *dNVtxProc = (int *) malloc((nRank) * sizeof(int));

  PDM_MPI_Allgather((void *) &dNVtx,
                    1,
                    PDM_MPI_INT, 
                    (void *) dNVtxProc, 
                    1, 
                    PDM_MPI_INT, 
                    comm);
  ppart->dVtxProc[0] = 1;
  for (int i = 1; i < nRank+1; i++) {
    ppart->dVtxProc[i] = dNVtxProc[i-1] + ppart->dVtxProc[i-1];
  }
  
  free(dNVtxProc);

  /* Boundaries definitions */

  ppart->_dFaceGroupIdx  = dFaceGroupIdx;
  ppart->_dFaceGroup = dFaceGroup;

  /* Dual graph */

  ppart->dDualGraphIdx = NULL;
  ppart->dDualGraph    = NULL;

  /* Partitions */

  ppart->nPart = nPart;
  
  ppart->mNPart = -1;

  ppart->dPartProc = (int *) malloc((nRank + 1) * sizeof(int));
  PDM_MPI_Allgather((void *) &nPart,
                    1,
                    PDM_MPI_INT, 
                    (void *) (&ppart->dPartProc[1]), 
                    1, 
                    PDM_MPI_INT, 
                    comm);

  ppart->dPartProc[0] = 0;
  for (int i = 1; i < nRank+1; i++) {
    ppart->mNPart = _PDM_part_MAX(ppart->mNPart, ppart->dPartProc[i]);
    ppart->dPartProc[i] = ppart->dPartProc[i] + ppart->dPartProc[i-1];
  }

  ppart->tNPart =  ppart->dPartProc[nRank];

  ppart->gPartTolProcPart = (int *) malloc(2*ppart->tNPart * sizeof(int));

  for (int i = 0; i < nRank; i++) {
    for (int j = ppart->dPartProc[i]; j < ppart->dPartProc[i+1]; j++) {
      ppart->gPartTolProcPart[2*j    ] = i;
      ppart->gPartTolProcPart[2*j + 1] = j - ppart->dPartProc[i];
    }
  }

  ppart->meshParts = (_part_t **) malloc(ppart->nPart * sizeof(_part_t *));

  for (int i = 0; i < ppart->nPart; i++) 
    ppart->meshParts[i] = NULL;

  /* Communicator */

  ppart->comm = comm;

  /* Method */

  ppart->split_method = split_method;

  int _method = PDM_part_renum_method_cell_idx_get(renum_cell_method);
  
  if (_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering cell method\n", renum_cell_method);
  }
  
  ppart->renum_cell_method = _method;
  
  _method = PDM_part_renum_method_face_idx_get(renum_face_method);
  
  if (_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering face method\n", renum_face_method);
  }
  ppart->renum_face_method = _method;

  ppart->dPartBound = NULL;

  /*
   * Build dual graph
   */

  if (dCellFace != NULL)
    _dual_graph_from_cell_face(ppart);
  else if (dFaceCell != NULL)
    _dual_graph_from_face_cell(ppart);
  else {
    PDM_printf("PDM_part_part_create error : dCellFace and dFaceCell are undefined, define one of two\n");
    exit(1);
  }

  int itime = 1;
  PDM_timer_hang_on(ppart->timer);
  ppart->times_elapsed[itime] = PDM_timer_elapsed(ppart->timer);
  ppart->times_cpu[itime]     = PDM_timer_cpu(ppart->timer);
  ppart->times_cpu_u[itime]   = PDM_timer_cpu_user(ppart->timer);
  ppart->times_cpu_s[itime]   = PDM_timer_cpu_sys(ppart->timer);
  itime += 1;

  /*
   * Graph partitioning
   */

  PDM_timer_resume(ppart->timer);

  int *cellPart;

  if (have_dCellPart == 0) {
    cellPart = (int *) malloc(dNCell * sizeof(int));
    _split(ppart,
           cellPart);
    for (int i = 0; i < dNCell; i++) {
      dCellPart[i] = cellPart[i];
    }
  }
  
  else {
    cellPart = (int *) ppart->_dCellPart;
  }
    
  if (1 == 0) {
    PDM_printf("cellPart : ");
    for (int i = 0; i <dNCell; i++)
      PDM_printf(" %d", cellPart[i]);
    PDM_printf("\n");
  }
  
  PDM_timer_hang_on(ppart->timer);
  ppart->times_elapsed[itime] = PDM_timer_elapsed(ppart->timer);
  ppart->times_cpu[itime]     = PDM_timer_cpu(ppart->timer);
  ppart->times_cpu_u[itime]   = PDM_timer_cpu_user(ppart->timer);
  ppart->times_cpu_s[itime]   = PDM_timer_cpu_sys(ppart->timer);
  itime += 1;
  
  /*
   * Cell distribution to build local connectivities
   *     - cellFaceIdx : ok
   *     - gCellFace   : ok
   *     - cellFace    : ok
   *     - cellLNToGN  : ok
   *     - cellTag     : ok
   *     - faceLNToGN  : ok
   */

  PDM_timer_resume(ppart->timer);
  
  _distrib_cell(ppart,
                cellPart);
  
  if (have_dCellPart == 0) {
    free(cellPart);
    cellPart = NULL;
  }

  /*
   * Face distribution to build local connectivities
   *     - faceVtxIdx  : ok
   *     - gFaceVtx    : ok 
   *     - faceVtx     : ok
   *     - faceTag     : ok
   */

  _distrib_face(ppart);
  _build_faceCell(ppart);

  /*
   * Vertex distribution to build local connectivities
   */

  _distrib_vtx(ppart); 
    
  /*
   * Cell renumbering
   */

  PDM_part_renum_cell (ppart); 
    
  /*
   * Face renumbering
   */

  PDM_part_renum_face (ppart); 

  /*
   * Look for partitioning boundary faces
   */

  _search_part_bound_face(ppart);

  /*
   * Face group distribution to build local connectivities
   */

  _distrib_face_groups(ppart);

  PDM_timer_hang_on(ppart->timer);
  ppart->times_elapsed[itime] = PDM_timer_elapsed(ppart->timer);
  ppart->times_cpu[itime]     = PDM_timer_cpu(ppart->timer);
  ppart->times_cpu_u[itime]   = PDM_timer_cpu_user(ppart->timer);
  ppart->times_cpu_s[itime]   = PDM_timer_cpu_sys(ppart->timer);

  ppart->times_elapsed[0]     = ppart->times_elapsed[itime];
  ppart->times_cpu[0]         = ppart->times_cpu[itime];
  ppart->times_cpu_u[0]       = ppart->times_cpu_u[itime];
  ppart->times_cpu_s[0]       = ppart->times_cpu_s[itime];

  for (int i = itime; i > 1; i--) {
    ppart->times_elapsed[i] -= ppart->times_elapsed[i-1];
    ppart->times_cpu[i]     -= ppart->times_cpu[i-1];
    ppart->times_cpu_u[i]   -= ppart->times_cpu_u[i-1];
    ppart->times_cpu_s[i]   -= ppart->times_cpu_s[i-1];
  }
}


void
PROCF (pdm_part_create_cf, PDM_PART_CREATE_CF)
(
 int                *ppartId,
 const PDM_MPI_Fint *fcomm,
 const int          *split_method,
 const char         *renum_cell_method,
 const int          *l_renum_cell_method,
 const char         *renum_face_method,
 const int          *l_renum_face_method,
 const int          *nPropertyCell,
 const int          *renum_properties_cell,
 const int          *nPropertyFace,
 const int          *renum_properties_face,
 const int          *nPart,
 const int          *dNCell,
 const int          *dNFace,
 const int          *dNVtx,
 const int          *nFaceGroup,
 const int          *have_dCellFace,
 const int          *dCellFaceIdx,
 const PDM_g_num_t  *dCellFace,
 const int          *have_dCellTag,
 const int          *dCellTag,
 const int          *have_dCellWeight,
 const int          *dCellWeight,
 const int          *have_dCellPart,
       int          *dCellPart,
 const int          *have_dFaceCell,
 const PDM_g_num_t  *dFaceCell,
 const int          *dFaceVtxIdx,
 const PDM_g_num_t  *dFaceVtx,
 const int          *have_dFaceTag,
 const int          *dFaceTag,
 const double       *dVtxCoord,
 const int          *have_dVtxTag,
 const int          *dVtxTag,
 const int          *dFaceGroupIdx,
 const PDM_g_num_t  *dFaceGroup
)
{
  
  const PDM_MPI_Comm c_comm    = PDM_MPI_Comm_f2c (*fcomm);
  const int *_dCellFaceIdx = dCellFaceIdx;
  const PDM_g_num_t *_dCellFace = dCellFace;
  const int *_dCellTag           = dCellTag;
  const int *_dCellWeight        = dCellWeight;
        int *_dCellPart          = dCellPart;
  const int *_dFaceTag           = dFaceTag;
  const int *_dVtxTag            = dVtxTag;
  const PDM_g_num_t *_dFaceCell = dFaceCell;
  
  if (*have_dCellFace == 0) {
    _dCellFaceIdx = NULL;
    _dCellFace    = NULL;
  }

  if (*have_dCellTag == 0)
    _dCellTag = NULL; 
  if (*have_dCellWeight == 0)
    _dCellWeight = NULL;

  if (*have_dFaceTag == 0)
    _dFaceTag = NULL;
  if (*have_dVtxTag == 0)
    _dVtxTag = NULL;

  if (*have_dFaceCell == 0)
    _dFaceCell = NULL;

  char *_renum_cell_method = 
      PDM_fortran_to_c_string(renum_cell_method, *l_renum_cell_method); 

  char *_renum_face_method = 
      PDM_fortran_to_c_string(renum_face_method, *l_renum_face_method); 
  
  PDM_part_create(ppartId,
                  c_comm,
                  (PDM_part_split_t) *split_method,
                  _renum_cell_method,
                  _renum_face_method,
                  *nPropertyCell,
                   renum_properties_cell,
                  *nPropertyFace,
                   renum_properties_face,
                  *nPart,
                  *dNCell,
                  *dNFace,
                  *dNVtx,
                  *nFaceGroup,
                  _dCellFaceIdx,
                  _dCellFace,
                  _dCellTag,
                  _dCellWeight,
                  *have_dCellPart,
                  _dCellPart,
                  _dFaceCell,
                  dFaceVtxIdx,
                  dFaceVtx,
                  _dFaceTag,
                  dVtxCoord,
                  _dVtxTag,
                  dFaceGroupIdx,
                  dFaceGroup);

  free (_renum_cell_method); 
  free (_renum_face_method);
  
}

/**
 *
 * \brief Return a mesh partition dimensions
 * 
 * \param [in]   ppartId            ppart identifier
 * \param [in]   ipart              Current partition
 * \param [out]  nCell              Number of cells
 * \param [out]  nFace              Number of faces
 * \param [out]  nFacePartBound     Number of partitioning boundary faces
 * \param [out]  nVtx               Number of vertices 
 * \param [out]  nProc              Number of processus
 * \param [out]  nTPart             Number of partitions
 * \param [out]  sCellFace          Size of cell-face connectivity 
 * \param [out]  sFaceVtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound     Size of facePartBound array
 * \param [out]  sFaceGroup         Size of faceGroup array 
 * \param [out]  nFaceGroup         Number of boundary
 *
 */

void 
PDM_part_part_dim_get
(
const   int  ppartId,
const   int  ipart,
 int        *nCell,
 int        *nFace,
 int        *nFacePartBound,
 int        *nVtx,
 int        *nProc,
 int        *nTPart,
 int        *sCellFace,
 int        *sFaceVtx,
 int        *sFaceGroup,
 int        *nFaceGroup
)
{
  _PDM_part_t *ppart = _get_from_id(ppartId);
  int numProcs;
  PDM_MPI_Comm_size(ppart->comm, &numProcs);

  _part_t *meshPart = NULL;
  if (ipart < ppart->nPart)
    meshPart  = ppart->meshParts[ipart];
  
  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_get error : unknown partition\n");
    exit(1);
  }

  *nCell           = meshPart->nCell;
  *nFace           = meshPart->nFace;
  *nFacePartBound  = meshPart->nFacePartBound;
  *nProc           = numProcs;
  *nTPart          = ppart->tNPart;  
  *nVtx            = meshPart->nVtx;
  *sCellFace       = meshPart->cellFaceIdx[*nCell];
  *sFaceVtx        = meshPart->faceVtxIdx[*nFace];
  *sFaceGroup      = 0;
  if (ppart->nFaceGroup > 0)
    *sFaceGroup    = meshPart->faceGroupIdx[ppart->nFaceGroup];
    *nFaceGroup    = ppart->nFaceGroup;
}

void 
PROCF (pdm_part_part_dim_get, PDM_PART_PART_DIM_GET)
(
 int           *ppartId,
 int           *ipart,
 int           *nCell,
 int           *nFace,
 int           *nFacePartBound,
 int           *nVtx,
 int           *nProc,
 int           *nTPart,
 int           *sCellFace,
 int           *sFaceVtx,
 int           *sFaceGroup,
 int           *nFaceGroup
)
{
  PDM_part_part_dim_get(*ppartId,
                        *ipart,
                        nCell,
                        nFace,
                        nFacePartBound,
                        nVtx,
                        nProc,
                        nTPart,
                        sCellFace,
                        sFaceVtx,
                        sFaceGroup,
                        nFaceGroup
                        );
}

/**
 *
 * \brief Return a mesh partition
 * 
 * \param [in]   ppartId            ppart identifier
 * \param [in]   ipart              Current partition
 * \param [out]  cellTag            Cell tag (size = nCell)
 * \param [out]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1)
 * \param [out]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace)
 * \param [out]  cellLNToGN         Cell local numbering to global numbering (size = nCell)
 * \param [out]  faceTag            Face tag (size = nFace)
 * \param [out]  faceCell           Face to cell connectivity  (size = 2 * nFace)
 * \param [out]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1) 
 * \param [out]  faceVtx            Face to Vertex connectivity (size = faceVtxIdx[nFace])
 * \param [out]  faceLNToGN         Face local numbering to global numbering (size = nFace)
 * \param [out]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
 * \param [out]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
 * \param [out]  facePartBound      Partitioning boundary faces (size = 4 * nFacePartBound)
                                         For each face :
                                          - Face local number
                                          - Connected process
                                          - Connected Partition 
                                            on the connected process 
                                          - Connected face local number 
                                            in the connected partition
 * \param [out]  vtxTag             Vertex tag (size = nVtx)
 * \param [out]  vtx                Vertex coordinates (size = 3 * nVtx)
 * \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx)
 * \param [out]  faceGroupIdx       face group index (size = nFaceGroup + 1)
 * \param [out]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup)
 * \param [out]  faceGroupLNToGN    faces global numbering for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup)
 *
 */

void PDM_part_part_val_get
(
const  int      ppartId,
const  int      ipart,
 int          **cellTag,
 int          **cellFaceIdx,
 int          **cellFace,
 PDM_g_num_t  **cellLNToGN,
 int          **faceTag,
 int          **faceCell,
 int          **faceVtxIdx,
 int          **faceVtx,
 PDM_g_num_t  **faceLNToGN,
 int          **facePartBoundProcIdx,
 int          **facePartBoundPartIdx,
 int          **facePartBound,
 int          **vtxTag,
 double       **vtx,
 PDM_g_num_t  **vtxLNToGN,
 int          **faceGroupIdx,
 int          **faceGroup,
 PDM_g_num_t  **faceGroupLNToGN
)
{
  _PDM_part_t *ppart = _get_from_id(ppartId);

  _part_t *meshPart = NULL;
  if (ipart < ppart->nPart)
    meshPart  = ppart->meshParts[ipart];
  
  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_val_get error : unknown partition\n");
    exit(1);
  }

  *cellTag              = meshPart->cellTag;
  *cellFaceIdx          = meshPart->cellFaceIdx;
  *cellFace             = meshPart->cellFace;
  *cellLNToGN           = meshPart->cellLNToGN;
  *faceTag              = meshPart->faceTag;
  *faceCell             = meshPart->faceCell;
  *faceVtxIdx           = meshPart->faceVtxIdx;
  *faceVtx              = meshPart->faceVtx;
  *faceLNToGN           = meshPart->faceLNToGN;
  *facePartBoundProcIdx = meshPart->facePartBoundProcIdx;
  *facePartBoundPartIdx = meshPart->facePartBoundPartIdx;
  *facePartBound        = meshPart->facePartBound;
  *vtxTag               = meshPart->vtxTag;
  *vtx                  = meshPart->vtx;
  *vtxLNToGN            = meshPart->vtxLNToGN;
  *faceGroupIdx         = meshPart->faceGroupIdx;
  *faceGroup            = meshPart->faceGroup;
  *faceGroupLNToGN      = meshPart->faceGroupLNToGN;
}


void 
PROCF (pdm_part_part_val_get, PDM_PART_PART_VAL_GET)
(
 int           *ppartId,
 int           *ipart,
 int           *cellTag,
 int           *cellFaceIdx,
 int           *cellFace,
 PDM_g_num_t   *cellLNToGN,
 int           *faceTag,
 int           *faceCell,
 int           *faceVtxIdx,
 int           *faceVtx,
 PDM_g_num_t   *faceLNToGN,
 int           *facePartBoundProcIdx,
 int           *facePartBoundPartIdx,
 int           *facePartBound,
 int           *vtxTag,
 double        *vtx,
 PDM_g_num_t   *vtxLNToGN,
 int           *faceGroupIdx,
 int           *faceGroup,
 PDM_g_num_t   *faceGroupLNToGN
)
{
  _PDM_part_t *ppart = _get_from_id(*ppartId);
  int numProcs;
  PDM_MPI_Comm_size(ppart->comm, &numProcs);

  _part_t *meshPart = NULL;
  if (*ipart < ppart->nPart)
    meshPart  = ppart->meshParts[*ipart];
  
  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_val_get error : unknown partition\n");
    exit(1);
  }

  for (int i = 0; i < meshPart->nCell; i++){
    if (meshPart->cellTag != NULL)
      cellTag[i]    = meshPart->cellTag[i];
    cellLNToGN[i] = meshPart->cellLNToGN[i];
  }

  for (int i = 0; i < meshPart->nCell + 1; i++)
    cellFaceIdx[i] = meshPart->cellFaceIdx[i];

  for (int i = 0; i < meshPart->cellFaceIdx[meshPart->nCell]; i++)
    cellFace[i] = meshPart->cellFace[i];

  for (int i = 0; i < meshPart->nFace; i++){
    if (meshPart->faceTag != NULL) 
      faceTag[i]    = meshPart->faceTag[i];
    faceLNToGN[i] = meshPart->faceLNToGN[i];
  }

  for (int i = 0; i < 2 * meshPart->nFace; i++){
    faceCell[i]    = meshPart->faceCell[i];
  }

  for (int i = 0; i < meshPart->nFace + 1; i++)
    faceVtxIdx[i] = meshPart->faceVtxIdx[i];

  for (int i = 0; i < faceVtxIdx[meshPart->nFace]; i++)
    faceVtx[i] = meshPart->faceVtx[i];

  for (int i = 0; i < 4 * meshPart->nFacePartBound; i++)
    facePartBound[i] = meshPart->facePartBound[i];

  for (int i = 0; i < numProcs + 1; i++)
    facePartBoundProcIdx[i] = meshPart->facePartBoundProcIdx[i];

  for (int i = 0; i < ppart->tNPart + 1; i++)
    facePartBoundPartIdx[i] = meshPart->facePartBoundPartIdx[i];

  for (int i = 0; i < meshPart->nVtx; i++){
    if (meshPart->vtxTag != NULL) 
      vtxTag[i]    = meshPart->vtxTag[i];
    vtxLNToGN[i] = meshPart->vtxLNToGN[i];
  }

  for (int i = 0; i < 3 * meshPart->nVtx; i++){
    vtx[i] = meshPart->vtx[i];
  }

  for (int i = 0; i < ppart->nFaceGroup + 1; i++)
    faceGroupIdx[i] = meshPart->faceGroupIdx[i];

  for (int i = 0; i < meshPart->faceGroupIdx[ppart->nFaceGroup]; i++) {
    faceGroup[i]       = meshPart->faceGroup[i];
    faceGroupLNToGN[i] = meshPart->faceGroupLNToGN[i];
  }
}

/**
 *
 * \brief Return a mesh partition
 * 
 * \param [in]   ppartId            ppart identifier
 * \param [in]   ipart              Current partition
 * \param [out]  cellColor          Cell tag (size = nCell)
 * \param [out]  faceColor          Face tag (size = nFace)

 */

void PDM_part_part_color_get
(
const  int      ppartId,
const  int      ipart,
 int          **cellColor,
 int          **faceColor
)
{
  _PDM_part_t *ppart = _get_from_id(ppartId);

  _part_t *meshPart = NULL;
  if (ipart < ppart->nPart)
    meshPart  = ppart->meshParts[ipart];
  
  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_val_get error : unknown partition\n");
    exit(1);
  }

  *cellColor = meshPart->cellColor;
  *faceColor = meshPart->faceColor;
}

void 
PROCF (pdm_part_part_color_get, PDM_PART_PART_COLOR_GET)
(
 int           *ppartId,
 int           *ipart,
 int           *cellColor,
 int           *faceColor
)
{
  _PDM_part_t *ppart = _get_from_id(*ppartId);
  int numProcs;
  PDM_MPI_Comm_size(ppart->comm, &numProcs);

  _part_t *meshPart = NULL;
  if (*ipart < ppart->nPart)
    meshPart  = ppart->meshParts[*ipart];
  
  if (meshPart == NULL) {
    PDM_printf("PDM_part_part_val_get error : unknown partition\n");
    exit(1);
  }

  if (meshPart->cellColor != NULL){
    for (int i = 0; i < meshPart->nCell; i++){
        cellColor[i]    = meshPart->cellColor[i];
    }
  }
  
  if (meshPart->faceColor != NULL){
    for (int i = 0; i < meshPart->nFace; i++){
        faceColor[i]    = meshPart->faceColor[i];
    }
  }
  
}

/**
 *
 * \brief Free ppart
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

void 
PDM_part_free
(
 int  ppartId
)
{
  _PDM_part_t *ppart = _get_from_id(ppartId);


  if (ppart->dCellFaceIdx != NULL)
    free(ppart->dCellFaceIdx);
  ppart->dCellFaceIdx = NULL;

  if (ppart->dCellFace != NULL)
    free(ppart->dCellFace);
  ppart->dCellFace = NULL;

  if (ppart->dCellProc != NULL)
    free(ppart->dCellProc);
  ppart->dCellProc = NULL;

  if (ppart->dFaceProc != NULL)
    free(ppart->dFaceProc);
  ppart->dFaceProc = NULL;

  if (ppart->dFaceCell != NULL)
    free(ppart->dFaceCell);
  ppart->dFaceCell = NULL;

  if (ppart->dVtxProc != NULL)
    free(ppart->dVtxProc);
  ppart->dVtxProc = NULL;

  if (ppart->dPartProc != NULL)
    free(ppart->dPartProc);
  ppart->dPartProc = NULL;
    
  if (ppart->gPartTolProcPart != NULL)
    free(ppart->gPartTolProcPart);
  ppart->gPartTolProcPart = NULL;

  if (ppart->dPartBound != NULL)
    free(ppart->dPartBound);
  ppart->dPartBound = NULL;

  if (ppart->dDualGraphIdx != NULL)
    free(ppart->dDualGraphIdx);
  ppart->dDualGraphIdx = NULL;

  if (ppart->dDualGraph != NULL)
    free(ppart->dDualGraph);
  ppart->dDualGraph = NULL;

  for (int i = 0; i < ppart->nPart; i++) { 
    _part_free(ppart->meshParts[i]);
    ppart->meshParts[i] = NULL;
  }

  PDM_timer_free(ppart->timer);
  ppart->timer = NULL;

  if (ppart->meshParts != NULL)
    free(ppart->meshParts);
  ppart->meshParts = NULL;
  
  free (ppart);

  PDM_Handles_handle_free (_pparts, ppartId, PDM_FALSE);

  const int n_part = PDM_Handles_n_get (_pparts);
  
  if (n_part == 0) {
    _pparts = PDM_Handles_free (_pparts);
  }
  
}

void 
PROCF (pdm_part_free, PDM_PART_FREE)
(
 int                *ppartId
 )
{
  PDM_part_free(*ppartId);
}

/**
 *
 * \brief Return times
 * 
 * \param [in]   ppartId     ppart identifier
 * \param [out]  elapsed     elapsed times (size = 4)
 * \param [out]  cpu         cpu times (size = 4)
 * \param [out]  cpu_user    user cpu times (size = 4)
 * \param [out]  cpu_sys     system cpu times (size = 4)
 *
 */

void PDM_part_time_get
(
 int       ppartId,
 double  **elapsed,
 double  **cpu,
 double  **cpu_user,
 double  **cpu_sys
)
{
  _PDM_part_t *ppart = _get_from_id(ppartId);

  *elapsed  = ppart->times_elapsed;
  *cpu      = ppart->times_cpu;
  *cpu_user = ppart->times_cpu_u;
  *cpu_sys  = ppart->times_cpu_s;
}

void 
PROCF (pdm_part_time_get, PDM_PART_TIME_GET)
(
 int      *ppartId,
 double   *elapsed,
 double   *cpu,
 double   *cpu_user,
 double   *cpu_sys
 )
{
  _PDM_part_t *ppart = _get_from_id(*ppartId);

  for (int i = 0; i < 4; i++) {
    elapsed[i]  = ppart->times_elapsed[i];
    cpu[i]      = ppart->times_cpu[i];
    cpu_user[i] = ppart->times_cpu_u[i];
    cpu_sys[i]  = ppart->times_cpu_s[i];
  }
}

/**
 *
 * \brief Return statistic
 * 
 * \param [in]   ppartId                        ppart identifier
 * \param [out]  cells_average                  average of cells number 
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells nummber
 * \param [out]  cells_max                      maximum of cells nummber
 * \param [out]  bound_part_faces_average       average of partitioning boundary faces 
 * \param [out]  bound_part_faces_median        median of partitioning boundary faces
 * \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 * \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 * \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 *
 */

void
PDM_part_stat_get
(
const int       ppartId,
      int      *cells_average, 
      int      *cells_median, 
      double   *cells_std_deviation, 
      int      *cells_min, 
      int      *cells_max,
      int      *bound_part_faces_average, 
      int      *bound_part_faces_median, 
      double   *bound_part_faces_std_deviation, 
      int      *bound_part_faces_min, 
      int      *bound_part_faces_max,
      int      *bound_part_faces_sum
)
{
  _PDM_part_t *ppart = _get_from_id(ppartId);
  int numProcs;
  PDM_MPI_Comm_size(ppart->comm, &numProcs);

  int *n_loc = (int *) malloc(ppart->nPart * sizeof(int));
  int *n_tot = (int *) malloc(ppart->dPartProc[numProcs] * sizeof(int));

  int *s_loc = (int *) malloc(ppart->nPart * sizeof(int));
  int *s_tot = (int *) malloc(ppart->dPartProc[numProcs] * sizeof(int));

  for (int i = 0; i < ppart->nPart; i++) {
    n_loc[i] = 0;
    s_loc[i] = 0;  
  }

  for (int i = 0; i < ppart->dPartProc[numProcs]; i++) {
    n_tot[i] = 0;
    s_tot[i] = 0;  
  }

  for (int i = 0; i < ppart->nPart; i++) {
    n_loc[i] = ppart->meshParts[i]->nCell;
    s_loc[i] = ppart->meshParts[i]->nFacePartBound; 
  } 

  int *nPartProc = (int *) malloc((numProcs) * sizeof(int));

  for (int i = 0; i < numProcs; i++) {
    nPartProc[i] = ppart->dPartProc[i+1] - ppart->dPartProc[i];
  }

  PDM_MPI_Allgatherv((void *) n_loc, 
                     ppart->nPart, 
                     PDM_MPI_INT,
                     (void *) n_tot, 
                     nPartProc, 
                     ppart->dPartProc, 
                     PDM_MPI_INT,
                     ppart->comm);

  PDM_MPI_Allgatherv((void *) s_loc, 
                     ppart->nPart, 
                     PDM_MPI_INT,
                     (void *) s_tot, 
                     nPartProc, 
                     ppart->dPartProc, 
                     PDM_MPI_INT,
                     ppart->comm);

  _quickSort_int(s_tot, 0, ppart->dPartProc[numProcs]-1);
  _quickSort_int(n_tot, 0, ppart->dPartProc[numProcs]-1);

  double   _cells_average; 

  double   _bound_part_faces_average; 

  *bound_part_faces_min = -1;
  *bound_part_faces_max = -1;
  *cells_min = -1;
  *cells_max = -1;
  _cells_average = 0;
  _bound_part_faces_average = 0;

  for (int i = 0; i < ppart->dPartProc[numProcs]; i++) {
    if (*bound_part_faces_min < 0)
      *bound_part_faces_min = s_tot[i];
    else
      *bound_part_faces_min = _PDM_part_MIN(*bound_part_faces_min, s_tot[i]);
    if (*bound_part_faces_max < 0)
      *bound_part_faces_max = s_tot[i];
    else
      *bound_part_faces_max = _PDM_part_MAX(*bound_part_faces_max, s_tot[i]);
    if (*cells_min < 0)
      *cells_min = n_tot[i];
    else
      *cells_min = _PDM_part_MIN(*cells_min, n_tot[i]);
    if (*cells_max < 0)
      *cells_max = n_tot[i];
    else
      *cells_max = _PDM_part_MAX(*cells_max, n_tot[i]);

    _cells_average += n_tot[i];
    _bound_part_faces_average += s_tot[i];
  }

  _cells_average = (_cells_average/((double) ppart->dPartProc[numProcs]));
  *bound_part_faces_sum = (int) _bound_part_faces_average;
  _bound_part_faces_average = 
    _bound_part_faces_average/((double) ppart->dPartProc[numProcs]);

  *cells_average = (int) round(_cells_average);
  *bound_part_faces_average = (int) round(_bound_part_faces_average);
  
  *cells_std_deviation = 0.;
  *bound_part_faces_std_deviation = 0.;
  for (int i = 0; i < ppart->dPartProc[numProcs]; i++) {
    *cells_std_deviation += (n_tot[i] - _cells_average) * (n_tot[i] - _cells_average);
    *bound_part_faces_std_deviation += (s_tot[i] - _bound_part_faces_average) * 
                                      (s_tot[i] - _bound_part_faces_average);
  }

  *cells_std_deviation = sqrt(*cells_std_deviation/ppart->dPartProc[numProcs]);
  *bound_part_faces_std_deviation = 
    sqrt(*bound_part_faces_std_deviation/ppart->dPartProc[numProcs]);

  int mid = ppart->dPartProc[numProcs]/2;
  if (ppart->dPartProc[numProcs] % 2 == 1) {
    *cells_median = n_tot[mid];
    *bound_part_faces_median = s_tot[mid];
  }

  else {
    *cells_median =(int) round((n_tot[mid-1] + n_tot[mid])/2.);
    *bound_part_faces_median = (int) ((s_tot[mid-1] + s_tot[mid])/2.);
  }

  free(n_tot);
  free(s_tot);
  free(n_loc);
  free(s_loc);
  free(nPartProc);
}

void 
PROCF (pdm_part_stat_get, PDM_PART_STAT_GET)    
(
const int      *ppartId,
      int      *cells_average, 
      int      *cells_median, 
      double   *cells_std_deviation, 
      int      *cells_min, 
      int      *cells_max,
      int      *bound_part_faces_average, 
      int      *bound_part_faces_median, 
      double   *bound_part_faces_std_deviation, 
      int      *bound_part_faces_min, 
      int      *bound_part_faces_max,
      int      *bound_part_faces_sum
)
{
  const int _ppartId = *ppartId;
  
  PDM_part_stat_get(_ppartId,
                 cells_average, 
                 cells_median, 
                 cells_std_deviation, 
                 cells_min, 
                 cells_max,
                 bound_part_faces_average, 
                 bound_part_faces_median, 
                 bound_part_faces_std_deviation, 
                 bound_part_faces_min, 
                 bound_part_faces_max,
                 bound_part_faces_sum);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
