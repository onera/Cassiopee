
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_geom_elem.h"
#include "pdm_order.h"
#include "pdm_part_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_distributed_graph_rank1_from_part.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 *
 * \brief
 *
 */
static
void debug_part
(
 int ppartId,
 int nPart
)
{
  int iRank;
  int nRank;

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &iRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &nRank);

  for (int ipart = 0; ipart < nPart; ipart++) {

    int nCell;
    int nFace;
    int nFacePartBound;
    int nVtx;
    int nProc;
    int nTPart;
    int sCellFace;
    int sFaceVtx;
    int sFaceGroup;
    int nFaceGroup2;

    PDM_part_part_dim_get(ppartId,
                          ipart,
                          &nCell,
                          &nFace,
                          &nFacePartBound,
                          &nVtx,
                          &nProc,
                          &nTPart,
                          &sCellFace,
                          &sFaceVtx,
                          &sFaceGroup,
                          &nFaceGroup2);

    int          *cellTag;
    int          *cellFaceIdx;
    int          *cellFace;
    PDM_g_num_t  *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t  *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t  *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t  *faceGroupLNToGN;

    PDM_part_part_val_get(ppartId,
                          ipart,
                          &cellTag,
                          &cellFaceIdx,
                          &cellFace,
                          &cellLNToGN,
                          &faceTag,
                          &faceCell,
                          &faceVtxIdx,
                          &faceVtx,
                          &faceLNToGN,
                          &facePartBoundProcIdx,
                          &facePartBoundPartIdx,
                          &facePartBound,
                          &vtxTag,
                          &vtx,
                          &vtxLNToGN,
                          &faceGroupIdx,
                          &faceGroup,
                          &faceGroupLNToGN);


    PDM_printf("[%i] nFaceGroup2    : %i\n", iRank, nFaceGroup2);
    PDM_printf("[%i] nCell          : %i\n", iRank, nCell);
    PDM_printf("[%i] nFace          : %i\n", iRank, nFace);
    PDM_printf("[%i] nVtx           : %i\n", iRank, nVtx);
    PDM_printf("[%i] nFacePartBound : %i\n", iRank, nFacePartBound);

    PDM_printf("[%i] cellFace     : ", iRank);
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        PDM_printf(" %i", cellFace[j]);
      }
      PDM_printf("\n");
    }

    PDM_printf("\n");

    PDM_printf("[%i]  cellLNToGN    : ", iRank);
    for (int i = 0; i < nCell; i++)
      PDM_printf(" "PDM_FMT_G_NUM, cellLNToGN[i]);
    PDM_printf("\n");

    PDM_printf("[%i] faceCell     : ", iRank);
    for (int i = 0; i < 2 * nFace; i++)
      PDM_printf(" %i", faceCell[i]);
    PDM_printf("\n");

    PDM_printf("[%i] faceVtx      : ", iRank);
    for (int i = 0; i < nFace; i++) {
      for (int j = faceVtxIdx[i]; j < faceVtxIdx[i+1]; j++) {
        PDM_printf(" %i", faceVtx[j]);
      }
      PDM_printf("\n");
    }

    PDM_printf("[%i]  faceLNToGN    : ", iRank);
    for (int i = 0; i < nFace; i++)
      PDM_printf(" "PDM_FMT_G_NUM, faceLNToGN[i]);
    PDM_printf("\n");

    PDM_printf("[%i] vtx           : ", iRank);
    for (int i = 0; i < 3 * nVtx; i++)
      PDM_printf(" %12.5e", vtx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] vtxLNToGN     : ", iRank);
    for (int i = 0; i <  nVtx; i++)
      PDM_printf(" "PDM_FMT_G_NUM, vtxLNToGN[i]);
    PDM_printf("\n");

    PDM_printf("[%i] faceGroupIdx : ", iRank);
    for (int i = 0; i < nFaceGroup2 + 1; i++)
      PDM_printf(" %i", faceGroupIdx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] faceGroup    : ", iRank);
    for (int i = 0; i < nFaceGroup2; i++) {
      for (int j = faceGroupIdx[i]; j < faceGroupIdx[i+1]; j++) {
        PDM_printf(" %i", faceGroup[j]);
      }
      PDM_printf("\n");
    }

    PDM_printf("[%i] nFacePartBound    : / nTPart = %i \n ", nFacePartBound, nTPart);
    for (int i = 0; i < nTPart; i++) {
      PDM_printf(" facePartBoundPartIdx[%i] = %i  -> %i \n ", i, facePartBoundPartIdx[i], facePartBoundPartIdx[i+1]);
      for (int j = facePartBoundPartIdx[i]; j < facePartBoundPartIdx[i+1]; j++) {
        PDM_printf(" %i  %i  %i  %i \n ", facePartBound[4*j], facePartBound[4*j+1], facePartBound[4*j+2], facePartBound[4*j+3]);
      }
      PDM_printf("\n");
    }

    PDM_printf("[%i] faceGroupLNToGN   : ", iRank);
    for (int i = 0; i < nFaceGroup2; i++) {
      for (int j = faceGroupIdx[i]; j < faceGroupIdx[i+1]; j++) {
        PDM_printf(" "PDM_FMT_G_NUM, faceGroupLNToGN[j]);
      }
      PDM_printf("\n");
    }
  }
}

/**
 *
 * \brief
 *
 */
static int
_unique_inplace
(
 PDM_g_num_t* a,
 int          n
)
{
  PDM_g_num_t tmp[n];
  for(int i = 0; i < n; i++){
    tmp[i] = a[i];
  }

  PDM_sort_long(tmp, 0, n);

  /* Setup first element */
  PDM_g_num_t oldE = tmp[0];

  for(int i = 0; i < n; i++){
    a[i] = -1;
  }

  // printf(" oooooooooooooooooooooooooooooooooooo \n");

  a[0] = tmp[0];
  int nNew = 1;
  for(int i = 1; i < n; i++ ){
    PDM_g_num_t newE = tmp[i];
    // printf(" newE = %i | oldE = %i | nNew = %i \n", newE, oldE, nNew);
    if(newE != oldE){
      a[nNew++] = tmp[i];
      oldE = tmp[i];
    }
  }
  // printf(" oooooooooooooooooooooooooooooooooooo \n");
  return nNew;
}

/**
 *
 * \brief
 *
 */
static void
compute_new_lntogn
(
 PDM_part_t*  ppart,
 PDM_g_num_t* shiftToGN,
 int*         shiftPart
)
{

  int iRank;
  int nRank;

  PDM_MPI_Comm_rank(ppart->comm, &iRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  printf("compute_new_lntogn\n");

  // Init
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    shiftPart[ipart] = 0;
  }

  // ------------------------------------
  // I - Step 1 - Collect nCell for current proc
  PDM_g_num_t nCellLoc = 0;
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {

    _part_t* meshPart = ppart->meshParts[ipart];

    printf(" [%i] nCell = %i  \n", ipart, meshPart->nCell);

    shiftPart[ipart] += nCellLoc;
    nCellLoc         += meshPart->nCell;
  }


  PDM_MPI_Allgather((void *) &nCellLoc,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) &shiftToGN[1],
                    1,
                    PDM__PDM_MPI_G_NUM,
                    ppart->comm);

  shiftToGN[0] = 1;
  for (int i = 1; i < nRank+1; i++) {
    shiftToGN[i] +=  shiftToGN[i-1];
  }


  if (1 == 1) {
    PDM_printf("shiftToGN : "PDM_FMT_G_NUM,  shiftToGN[0]);
    for (int i = 1; i < nRank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, shiftToGN[i]);
    }
    PDM_printf("\n");
  }

  printf("compute_new_lntogn end \n");
}

/**
 *
 * \brief
 *
 */
static _dist_csr*
compute_distributed_first_rank
(
 PDM_part_t*  ppart,
 PDM_g_num_t* shiftToGN,
 int*         shiftPart,
 PDM_g_num_t* cellPartBound,
 int          nCellBound
)
{
  printf("compute_distributed_first_rank \n");

  int iRank, nRank;
  PDM_MPI_Comm_rank(ppart->comm, &iRank);
  PDM_MPI_Comm_size(ppart->comm, &nRank);

  // Compute size for all arrays (shared by part )
  int nCellTot     = 0;
  int nCellFaceTot = 0;
  for (int ipart = 0; ipart < ppart->nPart; ipart++) {
    _part_t* meshPart = ppart->meshParts[ipart];
    nCellTot     += meshPart->nCell;
    nCellFaceTot += meshPart->cellFaceIdx[meshPart->nCell];
  }
  nCellFaceTot += nCellBound;

  printf(" nCellTot     = %i \n", nCellTot);
  printf(" nCellFaceTot = %i \n", nCellFaceTot);


  PDM_g_num_t cellCellN[nCellTot];

  /* -------------------------------------------------------------- */
  /* First step --> Prepare and count */
  for (int i = 0; i < nCellTot; i++) {
    cellCellN[i] = 0;
  }

  int shiftCellPartBound = 0;
  int nCellCellMax       = 0;
  for (int ipart = 0; ipart < ppart->nPart; ipart++) {
    _part_t* meshPart = ppart->meshParts[ipart];

    int sPart = shiftPart[ipart];
    // printf("[%i] sPart = %i \n", ipart, sPart);

    /* Internal to current part */
    for (int i = 0; i < meshPart->nFace; i++) {
      int iCell1 = PDM_ABS (meshPart->faceCell[2*i    ]);
      int iCell2 = PDM_ABS (meshPart->faceCell[2*i + 1]);
      if (iCell2 > 0) {
        cellCellN[iCell1+sPart-1] += 1;
        cellCellN[iCell2+sPart-1] += 1;
        nCellCellMax += 2;
      }
    }

    // Parse all
    for (int i = 0; i < ppart->tNPart; i++) {
      for (int j = meshPart->facePartBoundPartIdx[i]; j < meshPart->facePartBoundPartIdx[i+1]; j++) {

        int iFaceLoc  = meshPart->facePartBound[4*j];
        int iCellLoc  = PDM_ABS (meshPart->faceCell[2*(iFaceLoc-1)]); // [1, n]
        // int iCellOppG = cellPartBound[shiftCellPartBound++];
        // assert(iCellLoc+sPart-1 >= 0);
        // assert(iCellLoc+sPart-1 < nCellTot);
        cellCellN[iCellLoc+sPart-1] += 1;
        nCellCellMax += 1;
        // PDM_printf(" iFaceLoc = %i | iCellLoc = %i | iCellOppG = %i\n ", iFaceLoc, iCellLoc, iCellOppG);
      }
    }

  } /* End part loop */

  /* Verbose */
  if(0 == 1){
    for(int i = 0; i < nCellTot; i++){
      printf(" cellCellN[%d]   = "PDM_FMT_G_NUM" \n ", i, cellCellN[i]);
    }
  }
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  int cellCellIdx[nCellTot+1];
  PDM_g_num_t cellCell[nCellCellMax];
  cellCellIdx[0] = 0;
  for(int i = 0; i < nCellTot; i++) {
    cellCellIdx[i+1] = cellCellIdx[i] + cellCellN[i];
  }

  /* Verbose */
  if(0 == 1){
    printf(" nCellCellMax = %i \n", nCellCellMax);
    for(int i = 0; i < nCellTot+1; i++){
      printf(" cellCellIdx[%i]   = %i \n ", i, cellCellIdx[i]);
    }
  }
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  for (int i = 0; i < nCellTot; i++) {
    cellCellN[i] = 0;
  }
  for (int i = 0; i < nCellTot; i++) {
    cellCell[i] = -1000000;
  }

  shiftCellPartBound = 0;
  for (int ipart = 0; ipart < ppart->nPart; ipart++) {
    _part_t* meshPart = ppart->meshParts[ipart];

    int sPart = shiftPart[ipart];

    /* Internal to current part */
    for (int i = 0; i < meshPart->nFace; i++) {
      int iCell1 = PDM_ABS (meshPart->faceCell[2*i    ]);
      int iCell2 = PDM_ABS (meshPart->faceCell[2*i + 1]);
      if (iCell2 > 0) {
        int iCell1Loc = iCell1+sPart-1; // En "local"
        int iCell2Loc = iCell2+sPart-1;
        int idxCell1 = cellCellIdx[iCell1Loc] + cellCellN[iCell1Loc];
        int idxCell2 = cellCellIdx[iCell2Loc] + cellCellN[iCell2Loc];
        // Save
        // cellCell[idxCell1] = iCell2Loc+1; // [1,N]
        // cellCell[idxCell2] = iCell1Loc+1; // [1,N]
        cellCell[idxCell1] = iCell2Loc+shiftToGN[iRank]; // [1,N]
        cellCell[idxCell2] = iCell1Loc+shiftToGN[iRank]; // [1,N]
        // Go to next
        cellCellN[iCell1Loc] += 1;
        cellCellN[iCell2Loc] += 1;
      }
    }

    // Parse all
    for (int i = 0; i < ppart->tNPart; i++) {
      for (int j = meshPart->facePartBoundPartIdx[i]; j < meshPart->facePartBoundPartIdx[i+1]; j++) {

        int iFaceLoc  = meshPart->facePartBound[4*j];
        int iCellLoc  = PDM_ABS (meshPart->faceCell[2*(iFaceLoc-1)])+sPart-1; // [1, n]
        int iCellOppG = cellPartBound[shiftCellPartBound++];
        assert(iCellLoc >= 0);
        assert(iCellLoc < nCellTot);
        int idxCell1 = cellCellIdx[iCellLoc] + cellCellN[iCellLoc];

        // Save
        // printf(" idxCell1/iCellOppG = %i/%i \n", idxCell1, iCellOppG);
        cellCell[idxCell1] = iCellOppG;

        // Go to next
        cellCellN[iCellLoc] += 1;
      }
    }

  } /* End part loop */

  /* -------------------------------------------------------------- */


  /* -------------------------------------------------------------- */
  /* Verbose */
  if(0 == 1){
    printf(" -------------------------------\n");
    printf(" nCellCellMax = %i \n", nCellCellMax);
    for(int i = 0; i < nCellTot; i++){
      printf(" cellCell[%i] --> ", i);
      for(int j = cellCellIdx[i]; j < cellCellIdx[i+1]; j++){
        printf(PDM_FMT_G_NUM" ", cellCell[j]);
      }
      printf("\n");
    }
    printf(" -------------------------------\n");
  }
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  // Now we need to sort the CellCell and add diag
  //
  int* cellCellTmpIdx      = (int *) malloc( sizeof(int) * (nCellTot+1) );
  PDM_g_num_t* cellCellTmp = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (nCellCellMax+nCellTot) );

  // for(int i = 0; i < nCellCellMax+nCellTot; i++){
  //   cellCellTmp[i] = -1;
  // }

  int idxN = 0;
  cellCellTmpIdx[0] = 0;
  for(int i = 0; i < nCellTot; i++){
    // First copy diag
    int nCellConnect = 0;
    cellCellTmp[idxN+nCellConnect++] = i+shiftToGN[iRank]; // [1,N]
    for(int j = cellCellIdx[i]; j < cellCellIdx[i+1]; j++){
      cellCellTmp[idxN+nCellConnect++] = cellCell[j];
    }

    // Panic verbose
    // printf(" --------------------------- \n");
    // printf(" cellCell = ");
    // for(int j = idxN; j < idxN+nCellConnect; j++){
    //   printf("%i ", cellCellTmp[j]);
    // }
    // printf("\n");

    /* Sort and unique */
    nCellConnect = _unique_inplace(&cellCellTmp[idxN], nCellConnect);

    // Panic verbose
    // printf(" --------------------------- \n");
    // printf(" cellCell Sorted = ");
    // for(int j = idxN; j < idxN+nCellConnect; j++){
    //   printf("%i ", cellCellTmp[j]);
    // }
    // printf("\n");

    idxN += nCellConnect;
    cellCellTmpIdx[i+1] = cellCellTmpIdx[i] + nCellConnect;
  }

  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* Verbose */
  if(0 == 1){
    printf(" -------------------------------\n");
    printf(" cellCellTmpIdx[nCellTot] = %i \n", cellCellTmpIdx[nCellTot]);
    for(int i = 0; i < nCellTot; i++){
      printf(" cellCellS[%i] --> ", i);
      for(int j = cellCellTmpIdx[i]; j < cellCellTmpIdx[i+1]; j++){
        printf(PDM_FMT_G_NUM" ", cellCellTmp[j]);
      }
      printf("\n");
    }
    printf(" -------------------------------\n");
  }
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  // Setup struct to return
  cellCellTmp = (PDM_g_num_t *) realloc (cellCellTmp, sizeof(PDM_g_num_t) * cellCellTmpIdx[nCellTot]);
  // _dist_csr* dCsr = (_dist_csr *) malloc( sizeof(_dist_csr) );
  _dist_csr* dCsr = PDM_dist_csr_create(ppart->comm);

  dCsr->lSize  = nCellTot;
  dCsr->gSize  = shiftToGN[nRank]-1;
  dCsr->ia     = cellCellTmpIdx;
  dCsr->ja     = cellCellTmp;
  dCsr->shiftG = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (nRank+1) );
  for(int ii = 0; ii < nRank+1; ii++){
    dCsr->shiftG[ii] = shiftToGN[ii];
  }
  dCsr->dnnzIdx = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (nRank+1) );

  PDM_MPI_Allgather((void *) &dCsr->ia[dCsr->lSize],
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) &dCsr->dnnzIdx[1],
                    1,
                    PDM__PDM_MPI_G_NUM,
                    ppart->comm);

  dCsr->dnnzIdx[0] = 1;
  for (int i = 1; i < nRank+1; i++) {
    dCsr->dnnzIdx[i] +=  dCsr->dnnzIdx[i-1];
  }

  if (0 == 1) {
    PDM_printf("dnnzIdx : "PDM_FMT_G_NUM,  dCsr->dnnzIdx[0]);
    for (int i = 1; i < nRank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, dCsr->dnnzIdx[i]);
    }
    PDM_printf("\n");
  }

  // dist_csr_free(dCsr);
  // free(dCsr);

  /* -------------------------------------------------------------- */

  printf("compute_distributed_first_rank end \n");
  return dCsr;

}

/*============================================================================
 * Public function definitions
 *============================================================================*/


_dist_csr*
PDM_compute_dist_graph_rank1_from_ppart
(
 int ppartId,
 int nPart
)
{
  int iRank;
  int nRank;

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &iRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &nRank);

  if(0 == 1){
    debug_part(ppartId, nPart);
  }

  // Direct access to partition !!!
  PDM_part_t *ppart = PDM_get_meshpart_from_id(ppartId);

  PDM_g_num_t shiftToGN[nRank+1];
  int         shiftPart[ppart->nPart+1];
  compute_new_lntogn(ppart, shiftToGN, shiftPart);

  // Compute Buffer Size
  int nSendBuff = 0;
  for(int ipart = 0; ipart < nPart; ++ipart) {
    _part_t *meshPart = ppart->meshParts[ipart];
    for (int i = 0; i < ppart->tNPart; i++) {
      nSendBuff += (meshPart->facePartBoundPartIdx[i+1] - meshPart->facePartBoundPartIdx[i]);
    }
  }

  printf(" nSendBuff = %i  \n", nSendBuff);
  PDM_g_num_t sendBuff[nSendBuff];

  int sendBuffN[nRank];
  int recvBuffN[nRank];
  for(int i = 0; i < nRank; i++){
    sendBuffN[i] = 0;
  }

  // Second part ---> Exchnage the first rank
  for(int ipart = 0; ipart < nPart; ++ipart) {
    /* Get current part id */
    _part_t *meshPart = ppart->meshParts[ipart];

    printf(" [%i] nCell = %i  \n", ipart, meshPart->nCell);
    printf(" [%i] facePartBoundProcIdx = %i  \n", ipart, meshPart->facePartBoundProcIdx[nRank]);

    // Panic verbose
    // for (int i = 0; i < ppart->tNPart; i++) {
    //   for (int j = meshPart->facePartBoundPartIdx[i]; j < meshPart->facePartBoundPartIdx[i+1]; j++) {
    //     PDM_printf(" %i  %i  %i  %i \n ", meshPart->facePartBound[4*j], meshPart->facePartBound[4*j+1], meshPart->facePartBound[4*j+2], meshPart->facePartBound[4*j+3]);
    //   }
    // }

    // Prepare Send
    // On envoie le numero de cellules Globals dans la nouvelle numerotation ---> à faire
    for (int i = 0; i < ppart->tNPart; i++) {
      for (int j = meshPart->facePartBoundPartIdx[i]; j < meshPart->facePartBoundPartIdx[i+1]; j++) {
        int iFaceLoc  = meshPart->facePartBound[4*j];
        assert(meshPart->faceCell[2*(iFaceLoc-1)+1] == 0); // is a boudary
        int iProcCon  = meshPart->facePartBound[4*j+1];
        // int iPartOpp  = meshPart->facePartBound[4*j+2];
        // int iCellGlob = iCellLoc + shiftPart[ipart] + shiftToGN[iRank];

        // PDM_printf(" iProcCon = %i ", iProcCon);
        // PDM_printf(" iFaceLoc = %i / iCellGlob = %i | sPart = %i | sProc = %i \n ", iFaceLoc, iCellGlob, shiftPart[ipart], shiftToGN[iRank]);
        sendBuffN[iProcCon]   += 1;
      }
    }
  }

  /* Get number data to receive from each process */
  PDM_MPI_Alltoall(sendBuffN,
                   1,
                   PDM_MPI_INT,
                   recvBuffN,
                   1,
                   PDM_MPI_INT,
                   PDM_MPI_COMM_WORLD);

  int recvBuffIdx[nRank + 1];
  int sendBuffIdx[nRank + 1];
  recvBuffIdx[0] = 0;
  sendBuffIdx[0] = 0;
  for(int i = 0; i < nRank; i++) {
    recvBuffIdx[i+1] = recvBuffIdx[i] + recvBuffN[i];
    sendBuffIdx[i+1] = sendBuffIdx[i] + sendBuffN[i];
  }

  for(int i = 0; i < nRank; i++){
    sendBuffN[i] = 0;
  }

  // Second part ---> Exchnage the first rank
  for(int ipart = 0; ipart < nPart; ++ipart) {
    /* Get current part id */
    _part_t *meshPart = ppart->meshParts[ipart];

    printf(" [%i] nCell = %i  \n", ipart, meshPart->nCell);
    printf(" [%i] facePartBoundProcIdx = %i  \n", ipart, meshPart->facePartBoundProcIdx[nRank]);

    // Prepare Send
    // On envoie le numero de cellules Globals dans la nouvelle numerotation ---> à faire
    for (int i = 0; i < ppart->tNPart; i++) {
      for (int j = meshPart->facePartBoundPartIdx[i]; j < meshPart->facePartBoundPartIdx[i+1]; j++) {
        int iFaceLoc  = meshPart->facePartBound[4*j];
        int iCellLoc  = PDM_ABS (meshPart->faceCell[2*(iFaceLoc-1)])-1;
        assert(meshPart->faceCell[2*(iFaceLoc-1)+1] == 0); // is a boudary
        int iProcCon  = meshPart->facePartBound[4*j+1];
        // int iPartOpp  = meshPart->facePartBound[4*j+2];
        int iCellGlob = iCellLoc + shiftPart[ipart] + shiftToGN[iRank];

        // PDM_printf(" iProcCon = %i ", iProcCon);
        // PDM_printf(" iFaceLoc = %i / iCellGlob = %i | sPart = %i | sProc = %i \n ", iFaceLoc, iCellGlob, shiftPart[ipart], shiftToGN[iRank]);
        int idxBuff = sendBuffIdx[iProcCon] + sendBuffN[iProcCon];
        sendBuff [idxBuff]  = iCellGlob;
        sendBuffN[iProcCon] += 1;
      }
    }
  }

  if(0 == 1){
    for(int i = 0; i < nRank; i++){
      printf(" recvBuffN[%d]   = %d \n ", i, recvBuffN[i]);
      printf(" sendBuffN[%d]   = %d \n ", i, sendBuffN[i]);
      printf(" recvBuffIdx[%i] = %i | %i \n", i, recvBuffIdx[i], recvBuffIdx[i+1]);
    }
  }

  PDM_g_num_t recvBuff[recvBuffIdx[nRank]];

  /* Receive data from each process */
  PDM_MPI_Alltoallv(sendBuff,
                    sendBuffN,
                    sendBuffIdx,
                    PDM__PDM_MPI_G_NUM,
                    recvBuff,
                    recvBuffN,
                    recvBuffIdx,
                    PDM__PDM_MPI_G_NUM,
                    PDM_MPI_COMM_WORLD);

  /* Verbose */
  if(0 == 1){
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
    for(int jj = 0; jj < nRank; jj++){
      PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
      if(jj == iRank){
        for(int i = 0; i < recvBuffIdx[nRank]; i++){
          printf("[%i] recvBuff[%i]   = "PDM_FMT_G_NUM" \n ", iRank, i, recvBuff[i]);
        }
      }
    }
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
  }

  /* Allocate and return */
  // Now we have in recvBuffer the connected cells in other proc/part in global numbering
  _dist_csr* dCsrO1 = compute_distributed_first_rank(ppart, shiftToGN, shiftPart, recvBuff, recvBuffIdx[nRank]);

  return dCsrO1;

};


/*============================================================================
 * Private function definitions
 *============================================================================*/


#ifdef  __cplusplus
}
#endif
