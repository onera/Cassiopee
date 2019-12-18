/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef PDM_HAVE_OPENMP
#include <omp.h>
#endif
/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_order.h"
#include "pdm_cuthill.h"
#include "pdm_renum_cacheblocking.h"
#include "pdm_sort.h"
#include "pdm_part_graph.h"
#include "pdm_part_renum.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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

/**
 *
 * \brief Perform a cells renumbering with cache blocking
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_renum_cells_cacheblocking
(
 _part_t **meshParts,
 int       nPart,
 void     *specific_data
)
{
  // Voir avec erie mais il faut recuperer le method face ici ...
  int methodface = 0;

  int *renum_properties_cell = (int *) specific_data;

  int nCellPerCacheWanted = renum_properties_cell[0];
  int isAsynchrone        = renum_properties_cell[1];
  int isVectorisation     = renum_properties_cell[2];
  int nVectFace           = renum_properties_cell[3];
  int split_method        = renum_properties_cell[4];

  const char *_name = PDM_part_renum_method_face_name_get (methodface);

  if (strcmp (_name, "PDM_PART_RENUM_FACE_NONE")) {
   PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : face numbering for cacheblocking need to be set to PDM_PART_RENUM_FACE_NONE \n");
  }

  /* Loop over all part of the current process */
  for(int ipart = 0; ipart < nPart; ++ipart) {
    /* Get current part id */
    _part_t *part = meshParts[ipart];

    if(part->newToOldOrderCell == NULL){
      // printf(" SetUp newToOldOrderCell \n");
      part->newToOldOrderCell = (int *) malloc (sizeof(int) * part->nCell);
      for (int i = 0; i < part->nCell; i++){
        part->newToOldOrderCell[i] = i;
      }
    }

    if(part->newToOldOrderFace == NULL){
      part->newToOldOrderFace = (int *) malloc (sizeof(int) * part->nFace);
      for (int i = 0; i < part->nFace; i++){
        part->newToOldOrderFace[i] = i;
      }
    }


    PDM_renum_cacheblocking(part,
                            split_method,
                            nCellPerCacheWanted,
                            isAsynchrone,
                            isVectorisation,
                            nVectFace);

  }
}


/**
 *
 * \brief Perform a cells renumbering with cache blocking and OpenMP and Face Vectorisation
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_renum_cells_cacheblocking2
(
 _part_t **meshParts,
 int       nPart,
 void     *specific_data
)
{
  // Voir avec erie mais il faut recuperer le method face ici ...
  int methodface = 0;

  int *renum_properties_cell = (int *) specific_data;

  int nCellPerCacheWanted = renum_properties_cell[0];
  int isAsynchrone        = 0; //renum_properties_cell[1];
  int isVectorisation     = renum_properties_cell[2];
  int nVectFace           = renum_properties_cell[3];
  int split_method        = renum_properties_cell[4];

  const char *_name = PDM_part_renum_method_face_name_get (methodface);

  if (strcmp (_name, "PDM_PART_RENUM_FACE_NONE")) {
   PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking2 Error : face numbering for cacheblocking need to be set to PDM_PART_RENUM_FACE_NONE \n");
  }

  // We need to SetUp properly newToOldArray to store the final reordering order
  // For complex renumbering we use a reorder multiple timesand we want the ordering 1--->N
  for(int ipart = 0; ipart < nPart; ++ipart) {
    _part_t *part = meshParts[ipart];
    if(part->newToOldOrderCell == NULL){
      part->newToOldOrderCell = (int *) malloc (sizeof(int) * part->nCell);
      for (int i = 0; i < part->nCell; i++){
        part->newToOldOrderCell[i] = i;
      }
    }

    if(part->newToOldOrderFace == NULL){
      part->newToOldOrderFace = (int *) malloc (sizeof(int) * part->nFace);
      for (int i = 0; i < part->nFace; i++){
        part->newToOldOrderFace[i] = i;
      }
    }
  }
  // This pre-renumbering is necessary to the other method
  PDM_part_renum_cell (        meshParts,
                               nPart,
                               3,
                       (void*) renum_properties_cell);

  PDM_part_renum_face (        meshParts,
                               nPart,
                               2,
                       (void*) renum_properties_cell);

  /* Loop over all part of the current process */
  for(int ipart = 0; ipart < nPart; ++ipart) {
    /* Get current part id */
    _part_t *part = meshParts[ipart];

    PDM_renum_cacheblocking2(part,
                             split_method,
                             nCellPerCacheWanted,
                             isAsynchrone,
                             isVectorisation,
                             nVectFace);

  }
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \brief Add renumbering cacheblocking method in the ppart methods
 *
 */

static void
_computeIdx
(
int* arrayIdx,
int  size
)
{
  int tmp1 = arrayIdx[0];
  int tmp2 = -1;
  arrayIdx[0] = 0;
  for (int i = 0; i < size; i++){
    tmp2 = arrayIdx[i+1];
    arrayIdx[i+1] = arrayIdx[i] + tmp1;
    tmp1 = tmp2;
  }

}


/**
 *
 * \brief Prepare subdomain
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout
 *
 */

static void
_prepare_subdomain
(
 _part_t          *part,
 _subpartlayout_t *subpartlayout
)
{
  int tmp1     = -1;
  int tmp2     = -1;

  /* Init to zero */
  subpartlayout->nFaceInt =  0;
  subpartlayout->nFaceExt =  0;

  /* Identify max number of subdomaie for faces */
  for (int iface = 0; iface < part->nFace; iface++){
    tmp1 = PDM_MAX(tmp1, part->faceColor[iface]);

    int iCell2 = part->faceCell[2*iface + 1];
    if(iCell2 < 1){
      subpartlayout->nFaceExt = subpartlayout->nFaceExt + 1;
    }
    else{
      subpartlayout->nFaceInt = subpartlayout->nFaceInt + 1;
    }
  }

  /* Identify max number of subdomaie for faces */
  for (int iface = 0; iface < part->nFace; iface++){
    tmp1 = PDM_MAX(tmp1, part->faceColor[iface]);
  }

  /* Identify max number of subdomain for cells  */
  for (int icell = 0; icell < part->nCell; icell++){
    tmp2 = PDM_MAX(tmp2, part->cellColor[icell]);
  }

  /* Assert that subdomain are the same size for face and cell */
  if(tmp1 != tmp2)
  {
    PDM_error(__FILE__, __LINE__, 0, "Inconsistency betwenn cell and face tag \n");
  }

  subpartlayout->nSdom = tmp1 + 1;

  /* Verbose */
  if(0 == 1)
  {
    PDM_printf("nSdom     : " PDM_FMT_L_NUM" \n ", subpartlayout->nSdom);
    PDM_printf("nFaceExt  : " PDM_FMT_L_NUM" \n ", subpartlayout->nFaceExt);
    PDM_printf("nFaceInt  : " PDM_FMT_L_NUM" \n ", subpartlayout->nFaceInt);
  }

}


/**
 *
 * \brief Compute subdomain
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout
 *
 */

static void
_compute_subdomain
(
 _part_t          *part,
 _subpartlayout_t *subpartlayout
)
{
  int nSdom = subpartlayout->nSdom;

  if( subpartlayout->cellTileIdx != NULL){
    free(subpartlayout->cellTileIdx);
    subpartlayout->cellTileIdx = NULL;
  }
  if( subpartlayout->faceTileIdx != NULL){
    free(subpartlayout->faceTileIdx);
    subpartlayout->faceTileIdx = NULL;
  }
  if( subpartlayout->faceBndTileIdx != NULL){
    free(subpartlayout->faceBndTileIdx);
    subpartlayout->faceBndTileIdx = NULL;
  }
  /* Allocate */
  if( subpartlayout->cellTileIdx == NULL)
    subpartlayout->cellTileIdx    = (int *) malloc (sizeof(int) * (nSdom + 1) );
  if( subpartlayout->faceTileIdx == NULL)
    subpartlayout->faceTileIdx    = (int *) malloc (sizeof(int) * (nSdom + 1) );
  if( subpartlayout->faceBndTileIdx == NULL)
    subpartlayout->faceBndTileIdx = (int *) malloc (sizeof(int) * (nSdom + 1) );

  /* Initialise the Idx array  */
  for (int iSub = 0; iSub < nSdom+1; iSub++){
    subpartlayout->cellTileIdx   [iSub] = 0;
    subpartlayout->faceTileIdx   [iSub] = 0;
    subpartlayout->faceBndTileIdx[iSub] = 0;
  }

  /* Loop on internal faces */
  for (int iface = 0; iface < subpartlayout->nFaceInt; iface++){
    int color = part->faceColor[iface];
    subpartlayout->faceTileIdx[color]++;
  }

  /* Loop on external faces */
  for (int iface = subpartlayout->nFaceInt; iface < part->nFace; iface++){
    int color = part->faceColor[iface];
    subpartlayout->faceBndTileIdx[color]++;
  }

  /* Loop on cell */
  for (int icell = 0; icell < part->nCell; icell++){
    int color = part->cellColor[icell];
    subpartlayout->cellTileIdx[color]++;
  }

  /* Switch to Idx array */
  _computeIdx(subpartlayout->faceTileIdx   , nSdom);
  _computeIdx(subpartlayout->cellTileIdx   , nSdom);

  // Shift FaceBndTile
  // _computeIdx(subpartlayout->faceBndTileIdx, nSdom);
  int tmp1 = subpartlayout->faceBndTileIdx[0];
  subpartlayout->faceBndTileIdx[0] = subpartlayout->faceTileIdx[nSdom];
  for (int iSub = 0; iSub < nSdom; iSub++){
    int tmp2 = subpartlayout->faceBndTileIdx[iSub+1];
    subpartlayout->faceBndTileIdx[iSub+1] = subpartlayout->faceBndTileIdx[iSub] + tmp1;
    tmp1 = tmp2;
  }

  /* Panic verbose */
  if(0 == 1)
  {
    /* Debug faceTileIdx */
    PDM_printf(" --------------------- faceTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int iface = subpartlayout->faceTileIdx[iSub]; iface < subpartlayout->faceTileIdx[iSub+1]; iface++){
        printf("iface : %i \n", iface);
      }
    }

    /* Debug faceBndTileIdx */
    PDM_printf(" --------------------- faceBndTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int iface = subpartlayout->faceBndTileIdx[iSub]; iface < subpartlayout->faceBndTileIdx[iSub+1]; iface++){
        printf("iface : %i \n", iface);
      }
    }

    /* Debug faceBndTileIdx */
    PDM_printf(" --------------------- cellTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int icell = subpartlayout->cellTileIdx[iSub]; icell < subpartlayout->cellTileIdx[iSub+1]; icell++){
        printf("icell : %i \n", icell);
      }
    }

  }

}

/**
 *
 * \brief Compute maskage
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout
 *
 */

static void
_compute_mask
(
 _part_t          *part,
 _subpartlayout_t *subpartlayout
)
{
  int nSdom = subpartlayout->nSdom;

  if( subpartlayout->maskTileIdx != NULL){
    free(subpartlayout->maskTileIdx);
    subpartlayout->maskTileIdx = NULL;
  }
  if( subpartlayout->cellVectTileIdx != NULL){
    free(subpartlayout->cellVectTileIdx);
    subpartlayout->cellVectTileIdx = NULL;
  }
  if( subpartlayout->maskTileN != NULL){
    free(subpartlayout->maskTileN);
    subpartlayout->maskTileN = NULL;
  }
  if( subpartlayout->cellVectTileN != NULL){
    free(subpartlayout->cellVectTileN);
    subpartlayout->cellVectTileN = NULL;
  }

  if( subpartlayout->maskTile != NULL){
    free(subpartlayout->maskTile);
    subpartlayout->maskTile = NULL;
  }

  /* Allocate */
  if( subpartlayout->maskTileIdx == NULL)
    subpartlayout->maskTileIdx     = (int *) malloc (sizeof(int) * (nSdom + 1   ) );

  if( subpartlayout->cellVectTileIdx == NULL)
    subpartlayout->cellVectTileIdx = (int *) malloc (sizeof(int) * (nSdom + 1   ) );

  if( subpartlayout->maskTileN == NULL)
    subpartlayout->maskTileN       = (int *) malloc (sizeof(int) * (nSdom       ) );

  if( subpartlayout->cellVectTileN == NULL)
    subpartlayout->cellVectTileN   = (int *) malloc (sizeof(int) * (nSdom       ) );

  if( subpartlayout->maskTile == NULL)
    subpartlayout->maskTile        = (int *) malloc (sizeof(int) * (part->nCell ) );

  int* alreadyInit               = (int *) malloc (sizeof(int) * (part->nCell ) );

  /* Initialize */
  for (int iSub = 0; iSub < nSdom; iSub++){
    subpartlayout->maskTileN[iSub]     = 0;
    subpartlayout->cellVectTileN[iSub] = 0;
  }

  for (int iSub = 0; iSub < nSdom+1; iSub++){
    subpartlayout->maskTileIdx[iSub]     = 0;
  }

  for (int icell = 0; icell < part->nCell; icell++){
    alreadyInit[icell]             = -1;
    subpartlayout->maskTile[icell] = -100;
  }

  /* Begin computation of mask and vectTile */
  int offset = 0;
  int imask  = 0;

  for (int iSub = 0; iSub < nSdom; iSub++){
    int flagbeg = 0;
    int flag2   = 0;
    subpartlayout->maskTileIdx[iSub] = offset;

    /* Because ordering is design like this */
    subpartlayout->cellVectTileIdx[iSub] = subpartlayout->cellTileIdx[iSub];

    /* Step 1 */
    for (int icell = subpartlayout->cellTileIdx[iSub]; icell < subpartlayout->cellTileIdx[iSub+1]; icell++){

      flagbeg = 0;

      /* Face loop to flag border cell... */
      for (int idxface = part->cellFaceIdx[icell]; idxface < part->cellFaceIdx[icell+1]; idxface++){

        int iface = part->cellFace[idxface]-1;

        int iCell1 = part->faceCell[2*iface  ];
        int iCell2 = part->faceCell[2*iface+1];

        if(iCell2 == 0)
          continue;

        int color1 = part->cellColor[iCell1-1];
        int color2 = part->cellColor[iCell2-1];
        // PDM_printf("iSub : %i /C1 = %i / C2 = %i \n", iSub, color1, color2 );
        // if( (color1 != iSub) || (color2 != iSub) ){
        if( (color1 < iSub) || (color2 < iSub) ){
          // PDM_printf(" Cas 1 \n ");
          flagbeg = 1;
        }
        else{
          // PDM_printf(" Cas 2 \n ");
          flagbeg = PDM_MAX(0, flagbeg);
        }

      } /* End face loop */

      // PDM_printf("icell : %i /flagbeg = %i / flag2 = %i \n", icell, flagbeg, flag2 );

      /* Is the first element out ? */
      if( (flagbeg == 1) && (flag2 == 0) ){
        flag2 = 2;
        imask = offset;
        subpartlayout->maskTileIdx[iSub] = offset;
      }
      else if( flagbeg == 0){
        subpartlayout->cellVectTileN[iSub] += 1;
        alreadyInit[icell] = 1;
      }
      else if( (flagbeg == 0) && (flag2 == 2) ){
        PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error in compute mask \n");}
      }

      /* Interior face loop */
      for (int iface = subpartlayout->faceTileIdx[iSub]; iface < subpartlayout->faceTileIdx[iSub+1]; iface++){

        int iCell1 = part->faceCell[2*iface  ]-1;
        int iCell2 = part->faceCell[2*iface+1]-1;

        /* Left */
        if( alreadyInit[iCell1] == -1 ){
          alreadyInit[iCell1] = 1;
          subpartlayout->maskTileN[iSub] += 1;
          subpartlayout->maskTile[imask]  = iCell1;
          imask++;
        }
        /* Right */
        if( alreadyInit[iCell2] == -1 ){
          alreadyInit[iCell2] = 1;
          subpartlayout->maskTileN[iSub] += 1;
          subpartlayout->maskTile[imask]  = iCell2;
          imask++;
        }
      } /* End interior face loop */

      /* Exterior face loop */
      for (int iface = subpartlayout->faceBndTileIdx[iSub]; iface < subpartlayout->faceBndTileIdx[iSub+1]; iface++){

        int iCell1 = part->faceCell[2*iface  ]-1;
        /* Left */
        if( alreadyInit[iCell1] == -1 ){
          alreadyInit[iCell1] = 1;
          subpartlayout->maskTileN[iSub] += 1;
          subpartlayout->maskTile[imask]  = iCell1;
          imask++;
        }
      } /* End Exterior face loop */


      offset += subpartlayout->maskTileN[iSub];

  } /* End Sub dom */

  free(alreadyInit);

  for (int iSub = 0; iSub < nSdom; iSub++){
    int beg = part->subpartlayout->maskTileIdx[iSub];
    PDM_sort_int(&subpartlayout->maskTile[beg], NULL, subpartlayout->maskTileN[iSub]);
  }

  /* Panic verbose */
  if(0 == 1)
  {
    /* Debug faceBndTileIdx */
    PDM_printf(" --------------------- cellVectTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int icell = subpartlayout->cellVectTileIdx[iSub];
               icell < subpartlayout->cellVectTileIdx[iSub]+subpartlayout->cellVectTileN[iSub]; icell++){
        printf("icell : %i \n", icell);
      }
    }

  }
}

/**
 *
 * \brief Seach a candidate
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout
 *
 */
// static int
// _search_candidatesubdom
// (
//  int begS,
//  int endS,
//  int nThreadThreated,
//  int *isAssigneSdom,
//  int *ColorColorCellRank2Idx,
//  int *ColorColorCellRank2Arr,
//  int *isLock,
//  int *flagHaveBnd,
//  int  isBlocking,
//  int *cellTileIdx,
//  int *isLockCell,
//  int *CellCellIdx,
//  int *CellCell
// )
// {
//   int isOkToContinue = -1;
//   int iSubCandidate  = -1;

//   for (int iSubT = begS; iSubT < endS; iSubT++)
//   {
//     if(isAssigneSdom[iSubT] == -1 && flagHaveBnd[iSubT] == isBlocking )
//     {
//       isOkToContinue = 1;
//       for(int ptColor = ColorColorCellRank2Idx[iSubT  ];
//               ptColor < ColorColorCellRank2Idx[iSubT+1]; ptColor++)
//       {
//         // printf("isLock[%i] = %i \n", ColorColorCellRank2Arr[ptColor], isLock[ColorColorCellRank2Arr[ptColor]]);
//         if(isLock[ColorColorCellRank2Arr[ptColor]] != nThreadThreated)
//         {
//           // if(isLock[ColorColorCellRank2Arr[ptColor]] == 1)
//           if(isLock[ColorColorCellRank2Arr[ptColor]] >= 0)
//           {
//             isOkToContinue = -1;
//             break;
//           }
//         }
//       }

//       if(isOkToContinue == 1)
//       {
//         iSubCandidate = iSubT;
//         break;
//       }
//     }
//     if(isOkToContinue == 1){break;}
//   }
//   if(isOkToContinue == -1){iSubCandidate = -1;}
//   return iSubCandidate;
// }

static int
_search_candidatesubdom
(
 int begS,
 int endS,
 int nThreadThreated,
 int *isAssigneSdom,
 int *ColorColorCellRank2Idx,
 int *ColorColorCellRank2Arr,
 int *isLock,
 int *flagHaveBnd,
 int  isBlocking,
 int *cellTileIdx,
 int *isLockCell,
 int *CellCellIdx,
 int *CellCell
)
{
  int isOkToContinue = -1;
  int iSubCandidate  = -1;

  for (int iSubT = begS; iSubT < endS; iSubT++)
  {
    if(isAssigneSdom[iSubT] == -1 && flagHaveBnd[iSubT] == isBlocking )
    {

      isOkToContinue = 1;

      for(int icell = cellTileIdx[iSubT  ]; icell < cellTileIdx[iSubT+1]; icell++){
        if(isLockCell[icell] != nThreadThreated)
        {
          if(isLockCell[icell] >= 0)
          {
            isOkToContinue = -1;
            break;
          }
        }
        if(isOkToContinue == -1){break;}
      }

      if(isOkToContinue == 1)
      {
        for(int icell = cellTileIdx[iSubT  ]; icell < cellTileIdx[iSubT+1]; icell++){
          for(int idxCell = CellCellIdx[icell  ]; idxCell < CellCellIdx[icell+1]; idxCell++ ){

            if(isLockCell[CellCell[idxCell]] != nThreadThreated)
            {
              if(isLockCell[CellCell[idxCell]] >= 0)
              {
                isOkToContinue = -1;
                break;
              }
            }
            if(isOkToContinue == -1){break;}
          }
          if(isOkToContinue == -1){break;}
        }
      }

      if(isOkToContinue == 1)
      {
        iSubCandidate = iSubT;
        break;
      }
    }
    if(isOkToContinue == 1){break;}
  }
  if(isOkToContinue == -1){iSubCandidate = -1;}
  return iSubCandidate;
}

static int
_search_candidatesubdom_and_exclude
(
 int begS,
 int endS,
 int nThreadThreated,
 int *isAssigneSdom,
 int *ColorColorCellRank2Idx,
 int *ColorColorCellRank2Arr,
 int *isLock,
 int *flagHaveBnd,
 int  isBlocking,
 int *cellTileIdx,
 int *isLockCell,
 int *CellCellIdx,
 int *CellCell,
 int idxExcludeBeg
)
{
  int isOkToContinue = -1;
  int iSubCandidate  = -1;

  for (int iSubT = begS; iSubT < endS; iSubT++)
  {
    if(isAssigneSdom[iSubT] == -1 && flagHaveBnd[iSubT] == isBlocking )
    {

      isOkToContinue = 1;

      for(int ptColor = ColorColorCellRank2Idx[iSubT  ];
              ptColor < ColorColorCellRank2Idx[iSubT+1]; ptColor++)
      {
        if(ColorColorCellRank2Arr[ptColor] != iSubT)
        {
          // printf(" iSubCandidate : %i \n", iSubT);
          // printf(" ColorColorCellRank2Arr[ptColor] : %i \n", ColorColorCellRank2Arr[ptColor]);
          // printf(" idxExcludeBeg : %i \n", idxExcludeBeg);
          if(ColorColorCellRank2Arr[ptColor] < idxExcludeBeg)
          {
            isOkToContinue = -1;
          }
        }
      }

      if(isOkToContinue == 1)
      {
        for(int icell = cellTileIdx[iSubT  ]; icell < cellTileIdx[iSubT+1]; icell++){
          if(isLockCell[icell] != nThreadThreated)
          {
            if(isLockCell[icell] >= 0)
            {
              isOkToContinue = -1;
              break;
            }
          }
          if(isOkToContinue == -1){break;}
        }
      }

      if(isOkToContinue == 1)
      {
        for(int icell = cellTileIdx[iSubT  ]; icell < cellTileIdx[iSubT+1]; icell++){
          for(int idxCell = CellCellIdx[icell  ]; idxCell < CellCellIdx[icell+1]; idxCell++ ){

            if(isLockCell[CellCell[idxCell]] != nThreadThreated)
            {
              if(isLockCell[CellCell[idxCell]] >= 0)
              {
                isOkToContinue = -1;
                break;
              }
            }
            if(isOkToContinue == -1){break;}
          }
          if(isOkToContinue == -1){break;}
        }
      }

      if(isOkToContinue == 1)
      {
        iSubCandidate = iSubT;
        break;
      }
    }
    if(isOkToContinue == 1){break;}
  }
  if(isOkToContinue == -1){iSubCandidate = -1;}
  return iSubCandidate;
}

static void
_assign_candidatesubdom
(
 int  iSubCandidate,
 int  nThreadThreated,
 int *isAssigneSdom,
 int *ColorColorCellRank2Idx,
 int *ColorColorCellRank2Arr,
 int *isTreatedSdom,
 int *isLock
)
{
  isAssigneSdom[iSubCandidate] = 1;
  isLock[iSubCandidate]        = nThreadThreated;
  // isLock[iSubCandidate]        = 1;

  for(int ptColor = ColorColorCellRank2Idx[iSubCandidate  ];
          ptColor < ColorColorCellRank2Idx[iSubCandidate+1]; ptColor++)
  {
    if(isTreatedSdom[ColorColorCellRank2Arr[ptColor]] == -1)
    {
      // isLock[ColorColorCellRank2Arr[ptColor]] = 1;
      isLock[ColorColorCellRank2Arr[ptColor]] = nThreadThreated;
    }
  }
}

/**
 *
 * \brief Compute flag have bnd for asynchronous
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout
 *
 */

static void
_compute_flaghavebnd
(
 _part_t   *part,
 int       *flagHaveBnd,
 int        nSdom
)
{
  for(int i = 0; i < nSdom; i++) {
    flagHaveBnd[i] = 0;
  }

  /*
   * Loop on faces : if at least one is border, flagHaveBnd is set to 1
   */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];
    int color1 = part->cellColor[iCell1-1];

    /* Solution 2 */
    if(iCell2 == 0){flagHaveBnd[color1] = 1;}

  }

}

/* static int */
/* _get_color */
/* ( */
/*  int  *saveColor, */
/*  int  nColorFound, */
/*  int  firstColor, */
/*  int *lastColor */
/* ) */
/* { */
/*   int iColor = -1; */
/*   PDM_sort_int(saveColor, NULL, nColorFound); */


/*   // printf(" firstColor : %i \n", firstColor); */
/*   // printf(" lastColor  : %i \n", *lastColor); */
/*   saveColor[nColorFound] = saveColor[nColorFound-1]; */

/*   // Faire une recherche simple */
/*   int iColorMax = saveColor[nColorFound-1]; */
/*   int iColorMin = saveColor[0]; */

/*   if(nColorFound == 0){ */
/*     iColor = firstColor; */
/*     if(  (*lastColor) == -1 ) */
/*       (*lastColor)++; */
/*   } */
/*   else if(iColorMin > firstColor) */
/*   { */
/*     iColor = firstColor; */
/*     // printf(" Cas 1 : %i/%i\n", iColor, firstColor); */
/*   } */
/*   else */
/*   { */
/*     // Search gap */
/*     for(int isc = 0; isc < nColorFound; isc++){ */
/*       if( (saveColor[isc+1] - saveColor[isc] ) > 1 ) */
/*       { */
/*         iColor = saveColor[isc+1]-1; // Because Gap */
/*         // printf(" Cas 3 %i \n", iColor); */
/*         break; */
/*       } */
/*     } */

/*   } */


/*   if( (iColor == -1) && (iColorMax < *lastColor)) */
/*   { */
/*     iColor = *lastColor; */
/*     // printf(" Cas 2 : %i/%i\n", iColor, *lastColor); */
/*   } */

/*   if(iColor == -1 ) */
/*   { */
/*     // printf(" Add color T1 : %i \n", *lastColor); */
/*     (*lastColor)++; */
/*     iColor = *lastColor; */
/*   } */
/*   else if( iColor == *lastColor + 1) */
/*   { */
/*     // printf(" Add color T2 : %i \n", *lastColor); */
/*     abort(); */
/*     (*lastColor)++; */
/*     iColor = *lastColor; */
/*   } */

/*   if(1 == 1) */
/*   { */
/*     // printf(" iColor = %i -----> saveColor : ", iColor); */
/*     // printf(" -----> saveColor : "); */
/*     // for(int isc = 0; isc < nColorFound; isc++){ */
/*     //   printf(" %i ", saveColor[isc]); */
/*     // } */
/*     // printf("\n"); */

/*     for(int isc = 0; isc < nColorFound; isc++){ */
/*       assert(saveColor[isc] != iColor); */
/*     } */
/*   } */

/*   return iColor; */
/* } */


/* static void */
/* _compute_index_color */
/* ( */
/*  int  *begColor, */
/*  int  *flagInt, */
/*  int   colorMax, */
/*  int   BegV, */
/*  int   EndV */
/* ) */
/* { */

/*   for (int iColor = 0; iColor < colorMax+1; iColor++){ */
/*     begColor[iColor] = 0; */
/*   } */
/*   for (int icell = BegV; icell < EndV; icell++){ */
/*     begColor[flagInt[icell]]++; */
/*   } */

/*   if(0 == 1) */
/*   { */
/*     for (int iColor = 0; iColor < colorMax+1; iColor++){ */
/*       printf("nColor[%i] = %i \n", iColor, begColor[iColor]); */
/*     } */
/*   } */

/*   int nC = begColor[0]; */
/*   begColor[0] = BegV; */
/*   for (int iColor = 1; iColor < colorMax+1; iColor++){ */
/*     int tmp = begColor[iColor]; */
/*     begColor[iColor] = begColor[iColor-1]+nC; */
/*     nC = tmp; */
/*   } */

/*   if(0 == 1) */
/*   { */
/*     for (int iColor = 0; iColor < colorMax+1; iColor++){ */
/*       printf("[%i/%i] ----> begColor[%i] = %i \n", BegV, EndV, iColor, begColor[iColor]); */
/*     } */
/*   } */
/* } */

static void
_prepare_color
(
 _part_t  *part,
 int      *cellFaceLowerIdx,
 int      *cellFaceLowerArr,
 int      *flagInt,
 int      *flagBnd,
 int      nBlkCacheWanted
)
{

  int *multiColorCell = (int *) malloc(           (part->nCell) * sizeof(int));
  for(int icell = 0; icell < part->nCell; icell++){ flagInt[icell] = -1; };
  for(int icell = 0; icell < part->nCell; icell++){ flagBnd[icell] = -1; };
  for(int icell = 0; icell < part->nCell; icell++){ multiColorCell[icell] = -1; };

  for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){

    int BegV = part->subpartlayout->cellVectTileIdx[iSub];
    int EndV = part->subpartlayout->cellVectTileIdx[iSub]+part->subpartlayout->cellVectTileN[iSub];
    int EndS = part->subpartlayout->cellTileIdx[iSub+1];
    int nextMultiColor = 0;

    /* Interior */
    for (int icell = BegV; icell < EndV; icell++){

      for (int idxface = cellFaceLowerIdx[icell  ];
               idxface < cellFaceLowerIdx[icell+1]; idxface++){

        int iface     = cellFaceLowerArr[idxface];
        int iCell1    = part->faceCell[2*iface  ]-1;
        int iCell2    = part->faceCell[2*iface+1]-1;
        int lowerCell = PDM_MIN(iCell1, iCell2);
        if(lowerCell > -1){
          if(multiColorCell[lowerCell] == nextMultiColor){
            nextMultiColor++;
          }
        }
      }

      multiColorCell[icell] = nextMultiColor;
      flagInt[icell] = nextMultiColor;
    }

    for (int icell = BegV; icell < EndV; icell++){
      multiColorCell[icell] = -1;
    }

    nextMultiColor = 0;
    /* Exterior */
    for (int icell = EndV; icell < EndS; icell++){

      for (int idxface = cellFaceLowerIdx[icell  ];
               idxface < cellFaceLowerIdx[icell+1]; idxface++){

        int iface     = cellFaceLowerArr[idxface];
        int iCell1    = part->faceCell[2*iface  ]-1;
        int iCell2    = part->faceCell[2*iface+1]-1;
        int lowerCell = PDM_MIN(iCell1, iCell2);
        if(lowerCell > -1){
          if(multiColorCell[lowerCell] == nextMultiColor){
            nextMultiColor++;
          }
        }
      }

      multiColorCell[icell] = nextMultiColor;
      flagBnd[icell] = nextMultiColor;
    }

    for (int icell = EndV; icell < EndS; icell++){
      multiColorCell[icell] = -1;
    }

  }

  // for(int icell = 0; icell < part->nCell; icell++){
  //   printf("multiColorCell[%i] = %i | flagInt = %i | flagBnd = %i \n", icell, multiColorCell[icell], flagInt[icell], flagBnd[icell]);
  // };
  free(multiColorCell);
}


static void
_compute_lower_part
(
 _part_t  *part,
 int      *cellFaceLowerIdx,
 int      *cellFaceLowerArr
)
{
  int* tmpLower = (int *) malloc( part->nCell * sizeof(int));

  for (int icell = 0; icell < part->nCell; icell++){
    tmpLower[icell] = 0;
  }

  for (int icell = 0; icell < part->nCell; icell++){
    for (int idxface = part->cellFaceIdx[icell];
             idxface < part->cellFaceIdx[icell+1]; idxface++){

      int iface  = PDM_ABS(part->cellFace[idxface])-1;
      int iCell1 = part->faceCell[2*iface  ]-1;
      int iCell2 = part->faceCell[2*iface+1]-1;

      // printf("face Cell[%i] = %i/%i\n", iface, iCell1, iCell2);

      if(iCell1 < icell){tmpLower[icell]++;}
      if(iCell2 > -1 && iCell2 < icell){
        tmpLower[icell]++;
      }

    }
  }

  int maxNeightSize = -1;
  cellFaceLowerIdx[0] = 0;
  for (int icell = 0; icell < part->nCell; icell++){
    cellFaceLowerIdx[icell+1] = cellFaceLowerIdx[icell] + tmpLower[icell];
    maxNeightSize = PDM_MAX(maxNeightSize, tmpLower[icell]);
  }
  // int nLower = cellFaceLowerIdx[part->nCell];

  int* cellCellLowerArr = (int *) malloc( cellFaceLowerIdx[part->nCell] * sizeof(int));

  for (int icell = 0; icell < part->nCell; icell++){

    int begLower = cellFaceLowerIdx[icell  ];
    // int endLower = cellFaceLowerIdx[icell+1];

    for (int idxface = part->cellFaceIdx[icell];
             idxface < part->cellFaceIdx[icell+1]; idxface++){

      int iface  = PDM_ABS(part->cellFace[idxface])-1;
      int iCell1 = part->faceCell[2*iface  ]-1;
      int iCell2 = part->faceCell[2*iface+1]-1;

      if(iCell1 < icell){
        cellFaceLowerArr[begLower] = iface;
        cellCellLowerArr[begLower] = iCell1;
        begLower++;
      }
      if(iCell2 > -1 && iCell2 < icell){
        cellFaceLowerArr[begLower] = iface;
        cellCellLowerArr[begLower] = iCell2;
        begLower++;
      }
    }
  }

  // printf("maxNeightSize : %i \n", maxNeightSize);
  int locOrder[maxNeightSize];
  int tmp[maxNeightSize];

  for (int icell = 0; icell < part->nCell; icell++){

    int beg = cellFaceLowerIdx[icell];
    int nLo = cellFaceLowerIdx[icell+1]-beg;

    for(int i=0; i < nLo; i++ ){
      locOrder[i] = i;
    }

    if(nLo > 0)
      PDM_sort_int(&cellCellLowerArr[beg], locOrder, nLo);

    // for(int i=0; i < nLo; i++ ){
    //   printf("locOrder[%i] = %i \n", i, locOrder[i]);
    // }

    /* Reorder face */
    int idxL = 0;
    for (int idxface = cellFaceLowerIdx[icell];
             idxface < cellFaceLowerIdx[icell+1]; idxface++){

      int idxN = locOrder[idxL++];
      tmp[idxN] = cellFaceLowerArr[idxface];
    }

    idxL = 0;
    for (int idxface = cellFaceLowerIdx[icell];
             idxface < cellFaceLowerIdx[icell+1]; idxface++){
      cellFaceLowerArr[idxface] = tmp[idxL++];
    }

  }

  if(0 == 1)
  {
    printf(" ----------------------------- \n");
    for (int icell = 0; icell < part->nCell; icell++){

      printf(" cellFaceLowerIdx[%i] = %i is link to ", icell, cellFaceLowerIdx[icell] );

      for (int idxface = cellFaceLowerIdx[icell];
               idxface < cellFaceLowerIdx[icell+1]; idxface++){
        printf(" %i ", cellFaceLowerArr[idxface]);
        // printf(" %i ", cellCellLowerArr[idxface]);
      }
      printf(" \n");
    }
  }
  free(tmpLower);
  free(cellCellLowerArr);
}


/* static int */
/* _compute_multi_coloring_sdom */
/* ( */
/*  int  BegV, */
/*  int  EndV, */
/*  int *CellCell, */
/*  int *CellCellIdx, */
/*  int *flagInt, */
/*  int *saveColor, */
/*  int *begColor, */
/*  int *CellOrder */
/* ) */
/* { */
/*   int firstColor = 0; */
/*   int lastColor  = -1; */
/*   int colorMax   = -1; */

/*   for (int icell = BegV; icell < EndV; icell++){ */
/*     flagInt[icell] = -1; */
/*   } */

/*   for (int icell = BegV; icell < EndV; icell++){ */

/*     int nColorFound = 0; */
/*     for(int ptNeight = CellCellIdx[icell]; ptNeight < CellCellIdx[icell+1]; ptNeight++){ */
/*       int iCellLink = CellCell[ptNeight]; */

/*       if( (iCellLink < EndV ) && (iCellLink >= BegV ) && flagInt[iCellLink] != -1 ) */
/*       { */
/*         saveColor[nColorFound++] = flagInt[iCellLink]; */
/*       } */
/*     } */

/*     int iColor = -1; */
/*     iColor = _get_color(saveColor, nColorFound, firstColor, &lastColor); */

/*     // printf(" flagInt[%i] = %i \n", icell, iColor); */
/*     flagInt[icell] = iColor; */
/*     colorMax = PDM_MAX(iColor, colorMax); */
/*   } */

/*   _compute_index_color(begColor, flagInt, colorMax, BegV, EndV); */

/*   for (int icell = BegV; icell < EndV; icell++){ */
/*     flagInt[icell] = -1; */
/*   } */

/*   lastColor  = -1; */
/*   for (int icell = BegV; icell < EndV; icell++){ */

/*     int nColorFound = 0; */
/*     for(int ptNeight = CellCellIdx[icell]; ptNeight < CellCellIdx[icell+1]; ptNeight++){ */
/*       int iCellLink = CellCell[ptNeight]; */

/*       if( (iCellLink < EndV ) && (iCellLink >= BegV ) && flagInt[iCellLink] != -1 ) */
/*       { */
/*         saveColor[nColorFound++] = flagInt[iCellLink]; */
/*       } */
/*     } */

/*     int iColor = -1; */
/*     iColor = _get_color(saveColor, nColorFound, firstColor, &lastColor); */

/*     // printf(" flagInt[%i] = %i \n", icell, iColor); */
/*     // printf(" CellOrder[%i] = %i \n", begColor[iColor], icell); */
/*     CellOrder[begColor[iColor]++] = icell ; */
/*     flagInt[icell] = iColor; */
/*   } */

/*   return colorMax+1; */
/* } */


static void
_compute_sdom_graph
(
 int  BegV,
 int  EndV,
 int *CellCell,
 int *CellCellIdx,
 int *CellCellSub,
 int *CellCellSubIdx
)
{
  int nCel = EndV - BegV; // +1 ??
  // int nCel = EndV - BegV - 1; // +1 ??

  // printf(" _compute_sdom_graph ---> [%i/%i] \n", BegV, EndV);

  for(int icell = 0; icell < nCel +1; icell++){CellCellSubIdx[icell] = 0;}

  /* First loop to count */
  int idx = 0;
  CellCellSubIdx[0] = idx;
  for (int icell = BegV; icell < EndV; icell++){
    int icellSub = icell - BegV;
    CellCellSubIdx[icellSub+1] = idx;
    for(int ptNeight = CellCellIdx[icell]; ptNeight < CellCellIdx[icell+1]; ptNeight++){
      int iCellLink = CellCell[ptNeight];

      if( (iCellLink < EndV ) && (iCellLink >= BegV )){
        CellCellSubIdx[icellSub+1]++;
        CellCellSub[idx++] = iCellLink - BegV;
        // printf("adddd \n");
      }
    }
  }

  if(0 == 1){
    for(int icell = 0; icell < nCel; icell++){
      printf(" icell : %i ---> ", icell);
      // printf(" CellCellSubIdx[%i] : %i ---> ", icell  , CellCellSubIdx[icell]);
      // printf(" CellCellSubIdx[%i] : %i ---> ", icell+1, CellCellSubIdx[icell+1]);
      for(int ptNeight = CellCellSubIdx[icell]; ptNeight < CellCellSubIdx[icell+1]; ptNeight++){
        int iCellLink = CellCellSub[ptNeight];
        printf(" %i ", iCellLink);
      }
      printf("\n");
    }
  }

  for (int i = 0; i < nCel+1; i++){
    CellCellSubIdx[i] = CellCellSubIdx[i]+1;
  }
  for (int i = 0; i < CellCellSubIdx[nCel]; i++){
    CellCellSub[i] = CellCellSub[i]+1;
  }


}


/* static void */
/* _compute_sdom_graph_rank2 */
/* ( */
/*  int  BegV, */
/*  int  EndV, */
/*  int *CellCell, */
/*  int *CellCellIdx, */
/*  int *CellCellSub, */
/*  int *CellCellSubIdx */
/* ) */
/* { */
/*   int nCel = EndV - BegV; // +1 ?? */
/*   // int nCel = EndV - BegV - 1; // +1 ?? */

/*   // printf(" _compute_sdom_graph ---> [%i/%i] \n", BegV, EndV); */

/*   for(int icell = 0; icell < nCel +1; icell++){CellCellSubIdx[icell] = 0;} */

/*   /\* First loop to count *\/ */
/*   int idx = 0; */
/*   CellCellSubIdx[0] = idx; */
/*   for (int icell = BegV; icell < EndV; icell++){ */
/*     int icellSub = icell - BegV; */
/*     CellCellSubIdx[icellSub+1] = idx; */
/*     for(int ptNeight = CellCellIdx[icell]; ptNeight < CellCellIdx[icell+1]; ptNeight++){ */
/*       int iCellLink = CellCell[ptNeight]; */

/*       if( (iCellLink < EndV ) && (iCellLink >= BegV )){ */
/*         CellCellSubIdx[icellSub+1]++; */
/*         CellCellSub[idx++] = iCellLink - BegV; */
/*         // printf("adddd \n"); */
/*       } */
/*       for(int ptNeight2 = CellCellIdx[iCellLink]; ptNeight2 < CellCellIdx[iCellLink+1]; ptNeight2++){ */
/*         int iCellLink2 = CellCell[ptNeight2]; */
/*         if( (iCellLink2 < EndV ) && (iCellLink2 >= BegV )){ */
/*           CellCellSubIdx[icellSub+1]++; */
/*           CellCellSub[idx++] = iCellLink2 - BegV; */
/*           // printf("adddd \n"); */
/*         } */
/*       } */
/*     } */
/*   } */

/*   if(0 == 1){ */
/*     for(int icell = 0; icell < nCel; icell++){ */
/*       printf(" icell : %i ---> ", icell); */
/*       // printf(" CellCellSubIdx[%i] : %i ---> ", icell  , CellCellSubIdx[icell]); */
/*       // printf(" CellCellSubIdx[%i] : %i ---> ", icell+1, CellCellSubIdx[icell+1]); */
/*       for(int ptNeight = CellCellSubIdx[icell]; ptNeight < CellCellSubIdx[icell+1]; ptNeight++){ */
/*         int iCellLink = CellCellSub[ptNeight]; */
/*         printf(" %i ", iCellLink); */
/*       } */
/*       printf("\n"); */
/*     } */
/*   } */

/*   for (int i = 0; i < nCel+1; i++){ */
/*     CellCellSubIdx[i] = CellCellSubIdx[i]+1; */
/*   } */
/*   for (int i = 0; i < CellCellSubIdx[nCel]; i++){ */
/*     CellCellSub[i] = CellCellSub[i]+1; */
/*   } */


/* } */


static void
_compute_multi_coloring_sdom_rcm
(
 int  BegV,
 int  EndV,
 int *CellCell,
 int *CellCellIdx,
 int *flagInt,
 int *saveColor,
 int *begColor,
 int *CellOrder
)
{
  int nCel = EndV - BegV; // +1 ??

  if(nCel > 0){
    int *CellOrderSub = (int * ) malloc( nCel    * sizeof(int *));
    int *ColorSub     = (int * ) malloc( nCel    * sizeof(int *));
    // int *level_row    = (int * ) malloc( (nCel+1) * sizeof(int *));

    for (int icell = BegV; icell < EndV; icell++){
      flagInt[icell] = -1;
    }

    for (int icell = 0; icell < nCel; icell++  ){ColorSub[icell] = -10000000;}
    for (int icell = 0; icell < nCel; icell++  ){CellOrderSub[icell] = -10000000;}
    // for (int icell = 0; icell < nCel+1; icell++){level_row[icell] = -10000000;}

    // printf(" nCel : %i\n", nCel);

    PDM_genrcm(nCel, CellCellIdx, CellCell, CellOrderSub);

    // int nColor = 1;
    // PDM_compute_cyclic_rcm(nColor, CellOrderSub, ColorSub, nCel, level_row, level_num);

    for( int k = 0; k < nCel; k++){
      CellOrderSub[k] = CellOrderSub[k]-1;
    }

    for (int icell = BegV; icell < EndV; icell++){
      // flagInt[icell] = ColorSub[icell-BegV];
      flagInt[icell] = -1;
      // flagInt[icell] = ColorSub[CellOrderSub[icell-BegV]];
      // CellOrder[CellOrderSub[icell-BegV]+BegV] = icell;
      CellOrder[icell] = CellOrderSub[icell-BegV]+BegV;
      // CellOrder[icell] = icell;
    }

    // for (int icell = 0; icell < nCel; icell++){
    //   assert(ColorSub[icell] != -10000000);
    // }

    free(CellOrderSub);
    free(ColorSub);
  }


}

static void
_reorder_independant_faces
(
 int  begF,
 int  endF,
 int *flagFace,
 int *flagCell,
 int *faceCell,
 int *faceOrder,
 int *facePack,
 int *iPack
)
{
  int nFacSub = endF - begF;
  for (int iface = begF; iface < endF; iface++){
    flagFace[iface] = -1;
  }

  int nFacTreated = 0;
  int nVectFace   = nFacSub+10;
  // int nVectFace   = 64;
  while(nFacTreated != nFacSub){

    /* Reset for one pacquet all flag */
    for (int iface = begF; iface < endF; iface++){
      int iCell1 = faceCell[2*iface  ];
      int iCell2 = faceCell[2*iface+1];
      flagCell[iCell1-1] = -1;
      if(iCell2 > 0 )
        flagCell[iCell2-1] = -1;
    }

    (*iPack)++;

    /* Dramatic verbose */
    if(0 == 1){
      printf("nFacTreated/nFacSub : %i -> %i \n", nFacTreated, nFacSub);
    }

    int nPackFaceLoc = 0;

    /* Loop on face */
    for (int iface = begF; iface < endF; iface++){

      if(flagFace[iface] == -1){
        int iCell1 = faceCell[2*iface  ];
        int iCell2 = faceCell[2*iface+1];

        int t1 = flagCell[iCell1-1];
        int t2 = -1; // Si frontière pas de pb
        if(iCell2 > 0 )
          t2 = flagCell[iCell2-1];

        // printf("t1/t2 : %i/%i \n", t1, t2);
        if( (t1 == -1) && (t2 == -1) && nPackFaceLoc < nVectFace ){
          flagCell[iCell1-1] = 1;
          if(iCell2 > 0 )
            flagCell[iCell2-1] = 1;

          int bFac = begF;

          faceOrder[bFac+nFacTreated] = iface;

          flagFace[iface] = 1;
          facePack[bFac+nFacTreated] = (*iPack);

          nFacTreated++;
          nPackFaceLoc++;
        }
      } /* End Face loop */

    } /* End While */

    if(0 == 1){
      printf("begF/endF : %i -> %i --> %i \n", begF, endF, nPackFaceLoc);
    }
  } /* End While */


  for (int iface = begF; iface < endF; iface++){
    int iCell1 = faceCell[2*iface  ];
    int iCell2 = faceCell[2*iface+1];
    flagCell[iCell1-1] = -1;
    if(iCell2 > 0 )
      flagCell[iCell2-1] = -1;
  }

}

static void
_compute_multi_coloring
(
 _part_t   *part,
 int       *CellCellIdx,
 int       *CellCell,
 int       *CellOrder,
 int       *partCellIdx,
 int        nBlkCacheWanted
)
{
  PDM_renum_cacheblocking_compute_loop_array(part);

  int *flagInt   = (int *) malloc(  part->nCell * sizeof(int));
  int *flagBnd   = (int *) malloc(  part->nCell * sizeof(int));
  int *saveColor = (int *) malloc(  part->nCell * sizeof(int));
  int *begColor  = (int *) malloc(  part->nCell * sizeof(int));

  int *isTreatedF = (int *) malloc(  part->nFace * sizeof(int));
  int *faceOrder  = (int *) malloc(  part->nFace * sizeof(int));
  int *faceType   = (int *) malloc(  part->nFace * sizeof(int));

  int* CouplingIdx    = (int *) malloc( (nBlkCacheWanted+1 ) * sizeof(int)); // SurAlloc
  int* CouplingArr    = (int *) malloc( nBlkCacheWanted      * sizeof(int)); // SurAlloc
  int* Coupling1to1   = (int *) malloc( nBlkCacheWanted      * sizeof(int)); // SurAlloc


  for (int icell = 0; icell < part->nCell; icell++){
    flagInt[icell] = -1;
    flagBnd[icell] = -1;
  }

  for (int i = 0; i < part->nFace; i++){
     faceOrder[i] = -1;
     isTreatedF[i] = -1;
  }


  // int verifCell = 0;
  // int nFaceType = 0;
  int multiColorAlgo = 1;
  if(multiColorAlgo == 0){
    /* int colorMaxInt[nBlkCacheWanted]; */
    /* int colorMaxExt[nBlkCacheWanted]; */
    for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){

      /* int BegV = part->subpartlayout->cellVectTileIdx[iSub]; */
      /* int EndV = part->subpartlayout->cellVectTileIdx[iSub]+part->subpartlayout->cellVectTileN[iSub]; */
      /* int EndS = part->subpartlayout->cellTileIdx[iSub+1]; */
      // nFaceType = 0;

      /* colorMaxInt[iSub] = _compute_multi_coloring_sdom(BegV, EndV, */
      /*                                                  CellCell, CellCellIdx, flagInt, saveColor, begColor, */
      /*                                                  CellOrder); */

      // printf("[%i] ----> colorMax Int = %i \n", iSub, colorMaxInt[iSub]);


      /* colorMaxExt[iSub] = _compute_multi_coloring_sdom(EndV, EndS, */
      /*                                                  CellCell, CellCellIdx, flagBnd, saveColor, begColor, */
      /*                                                  CellOrder); */
      // printf("[%i] ----> colorMax Ext  = %i \n", iSub, colorMaxExt[iSub]);
    }
  }

  else{

    int* CellCellSub    = (int *) malloc( CellCellIdx[part->nCell] * sizeof(int *));
    int* CellCellSubIdx = (int *) malloc( (part->nCell + 1)        * sizeof(int *));
    // Si rank 2 is used ---> Need to allocate more memory !!!
    for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){

      // printf(" iSub ---> [%i/%i] \n", iSub, nBlkCacheWanted);
      int BegV = part->subpartlayout->cellVectTileIdx[iSub];
      int EndV = part->subpartlayout->cellVectTileIdx[iSub]+part->subpartlayout->cellVectTileN[iSub];
      int EndS = part->subpartlayout->cellTileIdx[iSub+1];
      // nFaceType = 0;
      _compute_sdom_graph(BegV, EndV, CellCell, CellCellIdx, CellCellSub, CellCellSubIdx);
      // _compute_sdom_graph_rank2(BegV, EndV, CellCell, CellCellIdx, CellCellSub, CellCellSubIdx);
      _compute_multi_coloring_sdom_rcm(BegV, EndV,
                                       CellCellSub, CellCellSubIdx, flagInt, saveColor, begColor,
                                       CellOrder);

      // printf("[%i] ----> colorMax Int = %i \n", iSub, colorMaxInt[iSub]);

      _compute_sdom_graph(EndV, EndS, CellCell, CellCellIdx, CellCellSub, CellCellSubIdx);
      // _compute_sdom_graph_rank2(EndV, EndS, CellCell, CellCellIdx, CellCellSub, CellCellSubIdx);
      _compute_multi_coloring_sdom_rcm(EndV, EndS,
                                       CellCellSub, CellCellSubIdx, flagBnd, saveColor, begColor,
                                       CellOrder);
      // colorMaxExt[iSub] = _compute_multi_coloring_sdom(EndV, EndS,
      //                                                  CellCell, CellCellIdx, flagBnd, saveColor, begColor,
      //                                                  CellOrder);
      // printf("[%i] ----> colorMax Ext  = %i \n", iSub, colorMaxExt[iSub]);
    }
    free(CellCellSub);
    free(CellCellSubIdx);


  }


  if(0 == 1)
  {
    int* verif2          = (int * ) malloc( part->nCell * sizeof(int) );
    for (int i = 0; i < part->nCell; i++){
      verif2[i] = 0;
    }
    for (int i = 0; i < part->nCell; i++){
      // printf("CellOrder[%i] = %i \n", i,  CellOrder[i]);
      // assert(CellOrder[i] != -1);
      verif2[CellOrder[i]]++;
    }
    for (int i = 0; i < part->nCell; i++){
      assert(verif2[i] == 1);
    }
    free(verif2);
  }
  // abort();

  /* Reorder all array to make easy the other part */
  PDM_part_reorder_cell(part, CellOrder);

  if(multiColorAlgo == 0){
    PDM_order_array(part->nCell,
                    sizeof(int),
                    CellOrder,
                    flagInt);

    PDM_order_array(part->nCell,
                    sizeof(int),
                    CellOrder,
                    flagBnd);
  }


  // for (int i = 0; i < part->nCell; i++){
  //   printf("flagInt[%i] = %i | flagExt[%i] : %i \n", i,  flagInt[i], i, flagBnd[i]);
  // }


  PDM_renum_cacheblocking_compute_loop_array(part);

  // Compute Lower part
  int *cellFaceLowerIdx = (int *) malloc(                 (part->nCell+1) * sizeof(int));
  int *cellFaceLowerArr = (int *) malloc(part->cellFaceIdx[part->nCell ]  * sizeof(int));
  _compute_lower_part(part, cellFaceLowerIdx, cellFaceLowerArr);

  _prepare_color(part,
                 cellFaceLowerIdx,
                 cellFaceLowerArr,
                 flagInt,
                 flagBnd, nBlkCacheWanted);

  /* Organize face with LowerPart */
  int iNewFace = 0;
  for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){

    int BegV = part->subpartlayout->cellVectTileIdx[iSub];
    int EndV = part->subpartlayout->cellVectTileIdx[iSub]+part->subpartlayout->cellVectTileN[iSub];
    int EndS = part->subpartlayout->cellTileIdx[iSub+1];

    // printf(" BegV : %i \n", BegV);
    // printf(" EndV : %i \n", EndV);
    // printf(" EndS : %i \n", EndS);
    // printf(" part->nCell : %i \n", part->nCell);
    // printf(" part->nFace : %i \n", part->nFace);
    // printf(" part->subpartlayout->faceTileIdx[%i] : %i \n", iSub, part->subpartlayout->faceTileIdx[iSub]);
    // printf(" part->subpartlayout->faceTileIdx[%i] : %i \n", iSub+1, part->subpartlayout->faceTileIdx[iSub+1]);
    // printf(" part->subpartlayout->faceBndTileIdx[%i] : %i \n", iSub, part->subpartlayout->faceBndTileIdx[iSub]);
    // printf(" part->subpartlayout->faceBndTileIdx[%i] : %i \n", iSub+1, part->subpartlayout->faceBndTileIdx[iSub+1]);

    int tag     = flagInt[0];
    int tagFace = 0;
    for (int icell = BegV; icell < EndS; icell++){

      // int begLower = cellFaceLowerIdx[icell  ];
      // int endLower = cellFaceLowerIdx[icell+1];

      if(icell >= EndV)
      {
        if(icell == EndV){
          tag     = flagBnd[icell];
          tagFace++;
        }
        else if(flagBnd[icell] != tag)
        {
          // printf(" Changement de tag Ext : \n");
          tag = flagBnd[icell];
          tagFace++;
        }
      }
      else if(flagInt[icell] != tag)
      {
        // printf(" Changement de tag : \n");
        tag = flagInt[icell];
        tagFace++;
      }


      // printf(" flagInt[%i] = %i || flagExt[%i] = %i || \n", icell, flagInt[icell], icell, flagBnd[icell]);

      for (int idxface = cellFaceLowerIdx[icell  ];
               idxface < cellFaceLowerIdx[icell+1]; idxface++){

        int iface  = cellFaceLowerArr[idxface];
        // int iCell1 = part->faceCell[2*iface  ]-1;
        // int iCell2 = part->faceCell[2*iface+1]-1;

        if(isTreatedF[iface] == -1)
        {
          // assert(part->faceColor[iface] == iSub);
          // printf(" -------------------------------------- \n");
          // printf(" Interior : il = %i || ir = %i \n", iCell1, iCell2);
          // printf(" Interior faceOrder[%i] = %i \n", iNewFace, iface);
          // if(part->faceColor[iface] != iSub)
          // {
          //   printf(" part->faceColor[%i] = %i -> %i/%i \n", iSub, part->faceColor[iface], iCell1, iCell2);
          //   abort();
          // }
          faceOrder[iNewFace] = iface;
          isTreatedF[iface] = 1;
          faceType[iNewFace] = tagFace;
          iNewFace++;
        }

      }
    }

    int firstCell = part->subpartlayout->maskTile[part->subpartlayout->maskTileIdx[iSub]];
    tag = flagBnd[firstCell];
    tagFace++;
    for (int mskcell = part->subpartlayout->maskTileIdx[iSub];
             mskcell < part->subpartlayout->maskTileIdx[iSub]+part->subpartlayout->maskTileN[iSub];
             mskcell++){

      int icell = part->subpartlayout->maskTile[mskcell];
      // printf(" maskcell : %i --> icell : %i  \n", mskcell, icell);

      if(flagBnd[icell] != tag)
      {
        // printf(" Changement de tag Boundary : %i ----> %i \n", tag, flagBnd[icell]);
        tag = flagBnd[icell];
        // tagFace++;
      }

      // printf(" flagInt[%i] = %i || flagExt[%i] = %i || \n", icell, flagInt[icell], icell, flagBnd[icell]);

      for (int idxface = cellFaceLowerIdx[icell  ];
               idxface < cellFaceLowerIdx[icell+1]; idxface++){

        int iface  = cellFaceLowerArr[idxface];
        // int iCell1 = part->faceCell[2*iface  ]-1;
        // int iCell2 = part->faceCell[2*iface+1]-1;

        if(isTreatedF[iface] == -1 && part->faceColor[iface] == iSub)
        {
          // printf(" -------------------------------------- \n");
          // printf(" Exterior : il = %i || ir = %i\n", iCell1, iCell2);
          // printf(" Exterior faceOrder[%i] = %i \n", iNewFace, iface);
          faceOrder[iNewFace] = iface;
          isTreatedF[iface] = 1;
          faceType[iNewFace] = tagFace;
          iNewFace++;
        }
      }
    }

    tagFace++;
    for (int iface = part->subpartlayout->faceTileIdx[iSub];
             iface < part->subpartlayout->faceTileIdx[iSub+1]; iface++){

      if(isTreatedF[iface] == -1 && part->faceColor[iface] == iSub)
      {
        // printf(" -------------------------------------- \n");
        // printf(" Exterior faceOrder[%i] = %i \n", iNewFace, iface);
        faceOrder[iNewFace] = iface;
        isTreatedF[iface] = 1;
        faceType[iNewFace] = tagFace;
        iNewFace++;
      }
    }

    tagFace++;
    for (int iface = part->subpartlayout->faceBndTileIdx[iSub];
             iface < part->subpartlayout->faceBndTileIdx[iSub+1]; iface++){

      // printf(" Boundary faceOrder[%i] = %i \n", iface, iface);
      faceOrder[iface] = iface;
      isTreatedF[iface] = 1;
      faceType[iface] = tagFace+1;
      // Tmp for elsA we shift by one to know the sign
      // part->faceColor[iface] = (part->faceColor[iface]+1);
    }

    // abort();


  }

  for (int i = 0; i < part->nFace; i++){
    // printf("faceOrder[%i] = %i \n", i,  faceOrder[i]);
    if(faceOrder[i] == -1)
    {
      abort();
    }
  }


  PDM_part_reorder_face(part, faceOrder);

  /* Organize independant pacquet */
  int *isTreatedC = (int *) malloc(  part->nCell * sizeof(int));
  int *facePack = (int *) malloc(  part->nFace * sizeof(int));

  for (int icell = 0; icell < part->nCell; icell++){
    isTreatedC[icell] = -1;
  }
  for (int iface = 0; iface < part->nFace; iface++){
    isTreatedF[iface] = -1;
    faceOrder [iface] = iface;
  }

  int begFaceType[part->nCell];
  // int faceColorTagSign = 1;
  int maxPack = -1;
  for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){

    int nFT = 1;
    int iPack = -1;
    begFaceType[0] = part->subpartlayout->faceTileIdx[iSub];
    int tmp1 = faceType[begFaceType[0]];
    for (int iface = part->subpartlayout->faceTileIdx[iSub];
             iface < part->subpartlayout->faceTileIdx[iSub+1]; iface++){
      if(tmp1 != faceType[iface]){
        begFaceType[nFT++] = iface;
      }
      tmp1 = faceType[iface];
    }
    begFaceType[nFT++] = part->subpartlayout->faceTileIdx[iSub+1];

    for (int iT = 0; iT < nFT-1; iT++){
      // printf("[%i] - begFaceType[%i] = %i\n", iSub, iT, begFaceType[iT]);

      _reorder_independant_faces(begFaceType[iT], begFaceType[iT+1],
                                 isTreatedF,
                                 isTreatedC,
                                 part->faceCell,
                                 faceOrder,
                                 facePack,
                                 &iPack);

      // Tmp for elsA we shift by one to know the sign
      // for (int iface = begFaceType[iT]; iface < begFaceType[iT+1]; iface++){
      //   part->faceColor[iface] = (part->faceColor[iface]+1)*faceColorTagSign;
      // }
      // faceColorTagSign = -1*faceColorTagSign;
    }
    maxPack = PDM_MAX(maxPack, iPack);

    int begFaceBnd = part->subpartlayout->faceBndTileIdx[iSub  ];
    int endFaceBnd = part->subpartlayout->faceBndTileIdx[iSub+1];
    _reorder_independant_faces(begFaceBnd, endFaceBnd,
                               isTreatedF,
                               isTreatedC,
                               part->faceCell,
                               faceOrder,
                               facePack,
                               &iPack);

  }

  // for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){
  //   for (int iface = part->subpartlayout->faceTileIdx[iSub];
  //            iface < part->subpartlayout->faceTileIdx[iSub+1]; iface++){
  //     printf("[%i] ---> faceType[%i] = %i \n", iSub, iface, faceType[iface]);
  //   }
  // }

  PDM_part_reorder_face(part, faceOrder);

  for (int iface = 0; iface < part->nFace; iface++){
    faceOrder [iface] = -1;
  }

  int *faceCellTmp = (int *) malloc( 2*part->nFace * sizeof(int));
  int *orderLoc    = (int *) malloc(   part->nFace * sizeof(int));
  int begPack[maxPack+2];

  // _order_lnum_s (faceCellTmp, 2, order, nFace);
  for (int iSub = 0; iSub < nBlkCacheWanted; iSub++){

    // printf(" [%i] ---> %i / %i \n", iSub, part->subpartlayout->faceTileIdx[iSub], part->subpartlayout->faceTileIdx[iSub+1]);
    for (int iPack = 0; iPack < maxPack+1; iPack++){
      begPack[iPack] = 0;
    }

    int maxPackLoc = -1;
    for (int iface = part->subpartlayout->faceTileIdx[iSub];
             iface < part->subpartlayout->faceTileIdx[iSub+1]; iface++){
      // printf("[%i] ---> facePack[%i] = %i \n", iSub, iface, facePack[iface]);
      begPack[facePack[iface]]++;
      maxPackLoc = PDM_MAX(maxPackLoc, facePack[iface]+1);
    }

    _computeIdx(begPack, maxPackLoc);


    for (int iPack = 0; iPack < maxPackLoc; iPack++){
      // printf("begPack[%i] = %i \n", iPack, begPack[iPack]);

      int begSub   = part->subpartlayout->faceTileIdx[iSub];
      int beg      = part->subpartlayout->faceTileIdx[iSub]+begPack[iPack];
      int nFaceLoc = begPack[iPack+1]-begPack[iPack];
      PDM_order_lnum_s (&part->faceCell[2*beg], 2, &faceOrder[beg], nFaceLoc);

      for (int iface = begPack[iPack];
               iface < begPack[iPack+1]; iface++){
        faceOrder[begSub+iface] = beg+faceOrder[begSub+iface];
      }

    }

    // for (int iface = part->subpartlayout->faceTileIdx[iSub];
    //          iface < part->subpartlayout->faceTileIdx[iSub+1]; iface++){
    //   printf("[%i] ---> faceOrder[%i] = %i \n", iSub, iface, faceOrder[iface]);
    // }

    for (int iface = part->subpartlayout->faceBndTileIdx[iSub];
             iface < part->subpartlayout->faceBndTileIdx[iSub+1]; iface++){
      faceOrder[iface] = iface;
    }

  }

  PDM_part_reorder_face(part, faceOrder);

  // for (int i = 0; i < part->nFace; i++){
  //   printf("faceCell[%i] = %i / %i ----> %i \n",i, part->faceCell[2*i]-1, part->faceCell[2*i+1]-1, part->faceColor[i]);
  // }

  // /* Connectivity change -> Recompute -- Caution CellCell and CellCellIdx is not valide*/
  // for (int i = 0; i < part->nFace; i++){
  //   if(faceOrder[i] == -1)
  //   {
  //     printf("faceOrder[%i] = %i \n", i,  faceOrder[i]);
  //     abort();
  //   }
  // }

  // if(1 == 1)
  // {
  //   int* veriff          = (int * ) malloc( part->nFace * sizeof(int) );
  //   for (int i = 0; i < part->nFace; i++){
  //     veriff[i] = 0;
  //   }
  //   for (int i = 0; i < part->nFace; i++){
  //     assert(faceOrder[i] != -1);
  //     veriff[faceOrder[i]]++;
  //     // printf("CellOrder[%i] = %i \n", i,  CellOrder[i]);
  //   }
  //   for (int i = 0; i < part->nFace; i++){
  //     if(veriff[i] != 1)
  //     {
  //       printf("veriff[%i] = %i \n", i,  veriff[i]);
  //       printf("faceOrder[%i] = %i \n", i,  faceOrder[i]);
  //     }
  //     assert(veriff[i] == 1);
  //   }
  //   free(veriff);
  // }

  free(begColor);
  free(flagInt);
  free(flagBnd);
  free(saveColor);
  free(faceOrder);
  free(faceType);
  free(isTreatedF);
  free(isTreatedC);
  free(CouplingIdx);
  free(CouplingArr);
  free(Coupling1to1);
  free(cellFaceLowerIdx);
  free(cellFaceLowerArr);
  free(facePack);
  free(faceCellTmp);
  free(orderLoc);
}

static void
_compute_distrib
(
 int  nCell,
 int  nPart,
 int* distribIdx
)
{
  // int ptCol = 0;
  int nCeT = nCell/nPart;
  int rCeT = nCell%nPart;

  // printf("nCell  = %i \n", nCell);
  // printf("nCeT   = %i \n", nCeT);
  // printf("rCeT   = %i \n", rCeT);

  for(int iPart = 0; iPart < nPart; iPart++){
    distribIdx[iPart] = nCeT;
    // printf("1) nCeTD[%i]   = %i \n", iPart, distribIdx[iPart]);
    if(rCeT > 0){
      distribIdx[iPart] += 1;
      rCeT--;
      // printf("2) nCeTD[%i]   = %i \n", iPart, distribIdx[iPart]);
    }
  }
  distribIdx[nPart] = 0;

  _computeIdx(distribIdx, nPart);

  // printf("nCell  = %i \n", nCell);
  // for(int iPart = 0; iPart < nPart+1; iPart++){
  //   printf("nCeTD[%i] = %i \n", iPart, distribIdx[iPart]);
  // }
}

static int
_compute_GPU_renumbering
(
 _part_t   *part,
 int        split_method,
 int        nCellPerCacheWanted,
 int       *CellOrder,
 int       *CellCellIdx,
 int       *CellCell
)
{
  /* Get nFac and nCel */
  const int nCell = part->nCell;

  int nThread;
  /* Split the graph */
#ifdef PDM_HAVE_OPENMP
#pragma omp parallel default(shared)
  {
    #pragma omp master
    {
     nThread = omp_get_num_threads();
    }
  }
#else
  nThread = 1;
#endif
  /* Compute graph associate to mesh */
  PDM_part_graph_compute_from_face_cell(          part,
                                        (int **) &CellCellIdx,
                                        (int **) &CellCell);


  int nBlkCacheWanted;
  if(nCellPerCacheWanted == 0){nBlkCacheWanted = 1;}
  else                        {nBlkCacheWanted = PDM_MAX(nCell/nCellPerCacheWanted,1);}

  nBlkCacheWanted += (nThread - nBlkCacheWanted%nThread );

  // int* CellWeight    = (int *) malloc( nCell              * sizeof(int));

  // for (int i = 0; i < nCell; i++){
  //   CellWeight[i] = CellCellIdx[i+1]-CellCellIdx[i];
  // }


  /* Split the graph */
  // PDM_part_graph_split(          split_method,
  //                                nThread,
  //                                part,
  //                                CellCellIdx,
  //                                CellCell,
  //                      (int *)   NULL,
  //                      (int *)   NULL,
  //                      (int **) &part->threadColor);

  if(part->threadColor == NULL){
    part->threadColor = (int *) malloc (sizeof(int) * nCell);
  }
  int nCeTD[nThread+1];
  _compute_distrib(nCell, nThread, nCeTD);

  int ptCol = 0;
  for (int i = 0; i < nCell; i++){
    if( i >= nCeTD[ptCol+1]){
      ptCol++;
    }
    // printf(" part->threadColor[%i] = %i \n", i, ptCol);
    part->threadColor[i] = ptCol;
    assert(part->threadColor[i] <= nThread );
  }
  // free(CellWeight);


  /*
   *  Identify the differents types of Cell :
   *              O -> Interior
   *              1 -> Boundary/join
   *              2 -> Thread Coupling
   */
  int* cellType    = (int *) malloc( nCell              * sizeof(int));

  for (int icell = 0; icell < nCell; icell++){
    cellType[icell] = -1;
  }


  // +-----------------------------+-----------------+------------------------------+
  // |         cellType            |    size         |           begin              |
  // +-----------------------------+-----------------+------------------------------+
  // | Full Coupling.............. |  1              |  0                           |
  // | Coupling 1to1.............. |  nCoupling1to1  |  1                           |
  // | Full Coupling Couronne..... |  1              |  1+nCoupling1to1             |
  // | Coupling 1to1 Couronne..... |  nCoupling1to1  |  2+nCoupling1to1             |
  // | Synchrone.................. |  nThread        |  2+2*nCoupling1to1           |
  // | Asynchrone................. |  nThread        |  2+2*nCoupling1to1+nThread   |
  // +-----------------------------+-----------------+------------------------------+
  const int nCoupling1to1 = nThread*(nThread-1)/2;
  const int nSdomType     = 2+nCoupling1to1+nCoupling1to1+nThread+nThread;
  int nCellByType[nSdomType];

  const int BeginFullCoupling   = 0;
  const int BeginCoupling1to1   = 1;
  const int BeginFullCouplingC  = 1+nCoupling1to1;
  const int BeginCoupling1to1C  = 2+nCoupling1to1;
  const int BeginSynchrone      = 2+2*nCoupling1to1;
  const int BeginAsynchrone     = 2+2*nCoupling1to1+nThread;

  const int EndFullCoupling   = BeginCoupling1to1;
  const int EndCoupling1to1   = BeginFullCouplingC;
  // const int EndFullCouplingC  = BeginCoupling1to1C;
  // const int EndCoupling1to1C  = BeginSynchrone;
  // const int EndSynchrone      = BeginAsynchrone;
  // const int EndAsynchrone     = nSdomType;

  if(0 == 1)
  {
    printf(" nThread            = %i \n", nThread);
    printf(" nSdomType          = %i \n", nSdomType);
    printf(" nCoupling1to1      = %i \n", nCoupling1to1);
    printf(" BeginFullCoupling  = %i \n", BeginFullCoupling );
    printf(" BeginFullCouplingC = %i \n", BeginFullCouplingC);
    printf(" BeginCoupling1to1  = %i \n", BeginCoupling1to1 );
    printf(" BeginCoupling1to1C = %i \n", BeginCoupling1to1C);
    printf(" BeginSynchrone     = %i \n", BeginSynchrone    );
    printf(" BeginAsynchrone    = %i \n", BeginAsynchrone   );
  }
  // int CouplingIdx[nThread+1];
  // int CouplingArr[nCoupling1to1];
  // int Coupling[nCoupling1to1];
  int* CouplingIdx    = (int *) malloc( (nThread+1 )    * sizeof(int));
  int* CouplingArr    = (int *) malloc( nCoupling1to1   * sizeof(int));
  int* Coupling       = (int *) malloc( nCoupling1to1   * sizeof(int));

  for (int iThr = 0; iThr < nThread+1; iThr++){
    CouplingIdx[iThr] = 0;
  }
  for (int iThr = 0; iThr < nThread; iThr++){
    CouplingIdx[iThr] = nThread - iThr - 1;
  }
  _computeIdx(CouplingIdx, nThread);

  int ptArrTh = 0;
  for (int iThr = 0; iThr < nThread; iThr++){
    int lThread = iThr+1;
    for(int iThrC = CouplingIdx[iThr]; iThrC < CouplingIdx[iThr+1]; iThrC++){
      CouplingArr[ptArrTh++] = lThread++;
    }
  }

  if(0 == 1)
  {
    printf("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n");
    for (int iThr = 0; iThr < nThread; iThr++){
      printf("CouplingIdx[%i] = %i / %i \n", iThr, CouplingIdx[iThr], CouplingIdx[iThr+1]);
      for(int iThrC=CouplingIdx[iThr]; iThrC < CouplingIdx[iThr+1]; iThrC++){
        printf("\t \t CouplingArr[%i] = %i \n", iThrC, CouplingArr[iThrC]);
      }
    }
    printf("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n");
  }


  for (int icell = 0; icell < nCell; icell++){
    int Tcolor = part->threadColor[icell];

    int isInterior       = 1;
    int nConnectedThread = 1;
    int ColorLink        = -1;
    for(int idxCell = CellCellIdx[icell  ];
            idxCell < CellCellIdx[icell+1];
            idxCell++ )
    {
      int ConnectedCell  = CellCell[idxCell];
      int ConnectedColor = part->threadColor[ConnectedCell];
      if(ConnectedColor != Tcolor){
        isInterior = -1;
        if(ColorLink == - 1){
          ColorLink = ConnectedColor;
          nConnectedThread++;
        }
        else if(ColorLink != ConnectedColor)
        {
          nConnectedThread++;
        }
        // break;
      }
    }
    int isBnd = -1;
    for(int idxCell = part->cellFaceIdx[icell  ];
            idxCell < part->cellFaceIdx[icell+1];
            idxCell++ )
    {
      int ConnectedFace  = PDM_ABS(part->cellFace[idxCell])-1;
      // int iCell1 = part->faceCell[2*ConnectedFace  ];
      int iCell2 = part->faceCell[2*ConnectedFace+1];
      if(iCell2 == 0){
        isBnd = 1;
        break;
      }
    }

    // if(nConnectedThread > 1)
    //   printf(" nConnectedThread : %i \n", nConnectedThread);

    if(isInterior ==  1 && isBnd == -1)
    {
      // -> Asynchrone (interior and not bnd)
      cellType[icell] = BeginAsynchrone + Tcolor;
    }
    else if(isInterior ==  1 && isBnd ==  1)
    {
      // -> Synchrone (interior and bnd)
      cellType[icell] =  BeginSynchrone + Tcolor;
    }
    else if(isInterior == -1 && isBnd == 1 && nConnectedThread == 2)
    {
      // -> Coupling 1to1 Couronne (not interior and bnd and connect to another thread)
      int iLin = PDM_MIN(Tcolor, ColorLink);
      int iCol = PDM_MAX(Tcolor, ColorLink);

      int idx = -1;
      for(int iThrC=CouplingIdx[iLin]; iThrC < CouplingIdx[iLin+1]; iThrC++){
        if(CouplingArr[iThrC] == iCol){
          idx = iThrC;
          break;
        }
      }
      // printf(" --> %i / %i --> %i \n", Tcolor, ColorLink, idx);

      cellType[icell] = BeginCoupling1to1C + idx;
    }
    else if(isInterior == -1 && isBnd ==  1)
    {
      // -> Full Coupling Couronne (not interior and bnd)
      cellType[icell] = BeginFullCouplingC;
    }
    else if(isInterior == -1 && isBnd ==  -1 && nConnectedThread == 2)
    {
      // -> Coupling 1to1 (not interior and not bnd and connect to another thread)
      int iLin = PDM_MIN(Tcolor, ColorLink);
      int iCol = PDM_MAX(Tcolor, ColorLink);

      int idx = -1;
      for(int iThrC=CouplingIdx[iLin]; iThrC < CouplingIdx[iLin+1]; iThrC++){
        if(CouplingArr[iThrC] == iCol){
          idx = iThrC;
          break;
        }
      }
      // printf(" --> %i / %i --> %i \n", Tcolor, ColorLink, idx);

      cellType[icell] = BeginCoupling1to1 + idx;
    }
    else
    {
      cellType[icell] = BeginFullCoupling;
    }
  }

  for (int i = 0; i < nSdomType; i++){
    nCellByType[i] = 0;
  }

  for (int icell = 0; icell < nCell; icell++){
    // printf("cellType[%i] = %i \n", icell, cellType[icell]);
    nCellByType[cellType[icell]]++;
  }

  // printf(" nCell = %i \n", nCell);

  int nCellVerif = 0;
  for (int i = 0; i < nSdomType; i++){
    nCellVerif += nCellByType[i];
  }
  assert(nCellVerif == nCell);


  /*
   * Reorganisz Cell to have [Couplage, Interior, Bnd]
   */
  int begCellByType[nSdomType];
  begCellByType[0] = 0;
  int offsetInt = nCellByType[0];
  for (int i = 1; i < nSdomType; i++){
    begCellByType[i] = offsetInt;
    offsetInt += nCellByType[i];
    // printf("offsetInt : %i \n", offsetInt );
  }

  if(0 == 1)
  {
    for (int i = 0; i < nSdomType; i++){
      printf(" begCellByType[%i] = %i --> %i \n", i, begCellByType[i], nCellByType[i]);
    }
  }

  for (int icell = 0; icell < nCell; icell++){
    cellType[icell] = -1;
    CellOrder[icell] = -1;
  }


  for (int icell = 0; icell < nCell; icell++){
    int Tcolor = part->threadColor[icell];

    if(cellType[icell] == -1)
    {

      int isInterior       = 1;
      int nConnectedThread = 1;
      int ColorLink        = -1;
      for(int idxCell = CellCellIdx[icell  ];
              idxCell < CellCellIdx[icell+1];
              idxCell++ )
      {
        int ConnectedCell  = CellCell[idxCell];
        int ConnectedColor = part->threadColor[ConnectedCell];
        if(ConnectedColor != Tcolor){
          isInterior = -1;
          if(ColorLink == - 1){
            ColorLink = ConnectedColor;
            nConnectedThread++;
          }
          else if(ColorLink != ConnectedColor)
          {
            nConnectedThread++;
          }
          // break;
        }
      }


      int isBnd = -1;
      for(int idxCell = part->cellFaceIdx[icell  ];
              idxCell < part->cellFaceIdx[icell+1];
              idxCell++ )
      {
        int ConnectedFace  = PDM_ABS(part->cellFace[idxCell])-1;
        // int iCell1 = part->faceCell[2*ConnectedFace  ];
        int iCell2 = part->faceCell[2*ConnectedFace+1];
        if(iCell2 == 0){
          isBnd = 1;
          break;
        }
      }

      // if(nConnectedThread > 0)
      //   printf(" nConnectedThread : %i \n", nConnectedThread);
      if(isInterior ==  1 && isBnd == -1)
      {
        // -> Asynchrone (interior and not bnd)
        cellType[icell] = BeginAsynchrone + Tcolor;
        CellOrder[begCellByType[BeginAsynchrone + Tcolor]++] = icell;
      }
      else if(isInterior ==  1 && isBnd ==  1)
      {
        // -> Synchrone (interior and bnd)
        cellType[icell] =  BeginSynchrone + Tcolor;
        CellOrder[begCellByType[BeginSynchrone + Tcolor]++] = icell;
      }
      else if(isInterior == -1 && isBnd == 1 && nConnectedThread == 2)
      {
        // -> Coupling 1to1 Couronne (not interior and bnd and connect to another thread)
        int iLin = PDM_MIN(Tcolor, ColorLink);
        int iCol = PDM_MAX(Tcolor, ColorLink);

        int idx = -1;
        for(int iThrC=CouplingIdx[iLin]; iThrC < CouplingIdx[iLin+1]; iThrC++){
          if(CouplingArr[iThrC] == iCol){
            idx = iThrC;
            break;
          }
        }
        // printf(" --> %i / %i --> %i \n", Tcolor, ColorLink, idx);

        cellType[icell] = BeginCoupling1to1C + idx;
        CellOrder[begCellByType[BeginCoupling1to1C + idx]++] = icell;
      }
      else if(isInterior ==  -1 && isBnd ==  1)
      {
        // -> Full Coupling Couronne (not interior and bnd)
        cellType[icell] = BeginFullCouplingC;
        CellOrder[begCellByType[BeginFullCouplingC]++] = icell;

      }
      else if(isInterior == -1 && isBnd == -1 && nConnectedThread == 2)
      {
        // -> Coupling 1to1 (not interior and not bnd and connect to another thread)
        int iLin = PDM_MIN(Tcolor, ColorLink);
        int iCol = PDM_MAX(Tcolor, ColorLink);

        int idx = -1;
        for(int iThrC=CouplingIdx[iLin]; iThrC < CouplingIdx[iLin+1]; iThrC++){
          if(CouplingArr[iThrC] == iCol){
            idx = iThrC;
            break;
          }
        }
        // printf(" --> %i / %i --> %i \n", Tcolor, ColorLink, idx);

        cellType[icell] = BeginCoupling1to1 + idx;
        CellOrder[begCellByType[BeginCoupling1to1 + idx]++] = icell;
      }
      else
      {
        cellType[icell] = BeginFullCoupling;
        CellOrder[begCellByType[BeginFullCoupling]++] = icell;
      }

    }
  }

  if(0 == 1)
  {
    int* verif          = (int * ) malloc( part->nCell * sizeof(int) );
    for (int i = 0; i < part->nCell; i++){
      verif[i] = 0;
    }
    // for (int i = 0; i < part->nCell; i++){
    //   printf("CellOrder[%i] = %i \n", i,  CellOrder[i]);
    // }
    for (int i = 0; i < part->nCell; i++){
      assert(CellOrder[i] != -1);
      verif[CellOrder[i]]++;
      // printf("CellOrder[%i] = %i \n", i,  CellOrder[i]);
    }
    for (int i = 0; i < part->nCell; i++){
      assert(verif[i] == 1);
    }
    free(verif);
  }


  PDM_part_reorder_cell(part, CellOrder);

  /* The graph changed so we rebuild the graph */
  free(CellCellIdx);
  free(CellCell);
  CellCellIdx = NULL;
  CellCell = NULL;
  /* Compute graph associate to mesh */
  PDM_part_graph_compute_from_face_cell(          part,
                                        (int **) &CellCellIdx,
                                        (int **) &CellCell);

  free(cellType);

  /*
   * Split the tree type of group to fit in Cache
   */
  if(part->cellColor == NULL){
    part->cellColor = (int * ) malloc( nCell              * sizeof(int) );
  }
  begCellByType[0] = 0;
  offsetInt = nCellByType[0];
  for (int i = 1; i < nSdomType; i++){
    begCellByType[i] = offsetInt;
    offsetInt += nCellByType[i];
    // printf("offsetInt : %i \n", offsetInt );
  }

  int* CellCellSubIdx = (int *) malloc( (nCell+1)          * sizeof(int) );
  int* CellCellSub    = (int *) malloc( CellCellIdx[nCell] * sizeof(int) );

  int offsetColor = 0;
  int colorByTypeIdx[nSdomType+1];

  for (int iSdom = 0; iSdom < nSdomType+1; iSdom++){colorByTypeIdx[iSdom] = 0;}

  for (int iSdom = 0; iSdom < nSdomType; iSdom++){
    int nC   = nCellByType[iSdom];
    int BegC = begCellByType[iSdom];
    int EndC = begCellByType[iSdom]+nC;
    colorByTypeIdx[iSdom] = 0;
    // printf(" iSdom -> %i ----> [%i,%i] \n", iSdom, BegC, EndC);

    /* Init */
    for (int icell = 0; icell < nC+1; icell++){CellCellSubIdx[icell] = 0;}

    /* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
    for (int icell = BegC; icell < EndC; icell++){
      for(int idxCell = CellCellIdx[icell  ];
              idxCell < CellCellIdx[icell+1]; idxCell++ ){

        int ConnectedCell   = CellCell[idxCell];

        if(ConnectedCell >= BegC && ConnectedCell < EndC )
        {
          CellCellSubIdx[icell-BegC]++;
        }
      }
    }
    /* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

    /* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
    _computeIdx(CellCellSubIdx, nC);
    /* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

    /* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
    int ptArr = 0;
    for (int icell = BegC; icell < EndC; icell++){
      for(int idxCell = CellCellIdx[icell  ];
              idxCell < CellCellIdx[icell+1]; idxCell++ ){

        int ConnectedCell = CellCell[idxCell];

        if(ConnectedCell >= BegC && ConnectedCell < EndC )
        {
          CellCellSub[ptArr] = CellCell[idxCell]-BegC;
          // printf(" [%i/%i] - CellCellSub[%i] = %i \n", iThr, iThrC, ptArr, CellCell[idxCell]-BegC);
          ptArr++;
        }
      }
    }
    /* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

    int nBlkCacheWantedSub;
    if(nCellPerCacheWanted == 0){nBlkCacheWantedSub = 1;}
    else                        {nBlkCacheWantedSub = PDM_MAX(nC/nCellPerCacheWanted,1);}

    int* ColorSub = NULL;
    if(nC != 0)
    {
      // PDM_part_graph_split_bis(          split_method,
      //                                    nBlkCacheWantedSub,
      //                                    nC,
      //                                    CellCellSubIdx,
      //                                    CellCellSub,
      //                          (int *)   NULL,
      //                          (int *)   NULL,
      //                                   &ColorSub);
       ColorSub = (int * ) malloc( sizeof(int) * nC);

       int distribIdx[nBlkCacheWantedSub+1];
       _compute_distrib(nC, nBlkCacheWantedSub, distribIdx);

       int ptCol2 = 0;
       for (int i = BegC; i < EndC; i++){
         // if( i/nCeT > ptCol){
         if( (i-BegC) >= distribIdx[ptCol2+1]){
           ptCol2++;
         }
         ColorSub[i-BegC] = ptCol2;
       }
     }

    if(0 == 1)
    {
      printf(" ----------- ColorSub : %i \n", nC);
      for (int i = 0; i < nC; i++){
        printf(" =====> ColorSub[%i]  = %i\n", i, ColorSub[i]);
      }
    }

    int offsetColorTmp = offsetColor;
    int nColorLoc = 0;
    for (int icell = 0; icell < nC; icell++){
      int newColor =  offsetColorTmp+ColorSub[icell];
      part->cellColor[BegC+icell] = newColor;
      // printf(" =====> cellColor[%i]  = %i\n", BegC+icell, newColor);
      offsetColor = PDM_MAX(offsetColor, newColor+1);
      nColorLoc   = PDM_MAX(nColorLoc, ColorSub[icell]+1);
    }
    colorByTypeIdx[iSdom] = nColorLoc;
    free(ColorSub);

  }

  // for (int iSdom = 0; iSdom < nSdomType+1; iSdom++){
  //   printf("colorByTypeIdx[%i] = %i\n", iSdom, colorByTypeIdx[iSdom]);
  // }

  _computeIdx(colorByTypeIdx, nSdomType);
  // for (int iSdom = 0; iSdom < nSdomType+1; iSdom++){
  //   printf("colorByTypeIdx[%i] = %i\n", iSdom, colorByTypeIdx[iSdom]);
  // }

  // int souDomType[nNewBlkCacheWanted];
  // for (int iSdomType = 0; iSdomType < nSdomType+1; iSdomType++){
  //   for (int iSdom = colorByTypeIdx[iSdomType]; iSdom < colorByTypeIdx[iSdomType+1]; iSdom++){
  //     souDomType[iSdom] = iSdomType;
  //     printf("souDomType[%i] = %i \n", iSdom, souDomType[iSdom] );
  //   }
  // }



  free(CellCellSubIdx);
  free(CellCellSub);


  int nNewBlkCacheWanted = 0;
  for (int icell = 0; icell < nCell; icell++){
    nNewBlkCacheWanted = PDM_MAX(nNewBlkCacheWanted, part->cellColor[icell]);
  }
  nNewBlkCacheWanted = nNewBlkCacheWanted + 1;

  // Renumber properly
  int* nCellPerCache    = (int *) malloc( nNewBlkCacheWanted      * sizeof(int));
  int* cellTileIdx      = (int *) malloc((nNewBlkCacheWanted + 1) * sizeof(int));

  /* Init to Zero */
  for (int i = 0; i < nNewBlkCacheWanted; i++){
    nCellPerCache[i] = 0;
  }

  for (int iSub = 0; iSub < nNewBlkCacheWanted+1; iSub++){
    cellTileIdx[iSub]=0;
  }

  /* First loop to identify the size of each color */
  for (int i = 0; i < nCell; i++){
    int color = part->cellColor[i];
    nCellPerCache[color]++;
    CellOrder[i] = -1;
  }

  /* Compute the index of each Sub-Domain **/
  cellTileIdx[0] = 0;
  for (int i = 0; i < nNewBlkCacheWanted; i++){
    cellTileIdx[i + 1] = cellTileIdx[i] + nCellPerCache[i];
  }

  /* Reset array */
  for (int i = 0; i < nNewBlkCacheWanted; i++){
    nCellPerCache[i] = 0;
  }

  for (int i = 0; i < nCell; i++){
    int color  = part->cellColor[i];
    int idx    = cellTileIdx[color] + nCellPerCache[color];

    CellOrder[idx] = i;
    nCellPerCache[color]++;
  }

  PDM_part_reorder_cell(part, CellOrder);

  /* The graph changed so we rebuild the graph */
  free(CellCellIdx);
  free(CellCell);
  CellCellIdx = NULL;
  CellCell = NULL;

  // printf("oooooooooo\n");
  // return nNewBlkCacheWanted;

  /* Compute graph associate to mesh */
  PDM_part_graph_compute_from_face_cell(          part,
                                        (int **) &CellCellIdx,
                                        (int **) &CellCell);



  int *ColorColorIdx           = (int *) malloc( (nNewBlkCacheWanted + 1) * sizeof(int));
  int *ColorColorRank1Idx      = (int *) malloc( (nNewBlkCacheWanted + 1) * sizeof(int));
  int *ColorColorCellRank2Idx  = (int *) malloc( (nNewBlkCacheWanted + 1) * sizeof(int));
  int *isAlreadyAssignedRank1  = (int *) malloc( (nNewBlkCacheWanted    ) * sizeof(int));
  int *isAlreadyAssignedRank2  = (int *) malloc( (nNewBlkCacheWanted    ) * sizeof(int));

  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
    isAlreadyAssignedRank1[iSub] = 0;
    isAlreadyAssignedRank2[iSub] = 0;
  }

  for (int iSub = 0; iSub < nNewBlkCacheWanted+1; iSub++){
    ColorColorIdx         [iSub] = 0;
    ColorColorRank1Idx    [iSub] = 0;
    ColorColorCellRank2Idx[iSub] = 0;
  }

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  /* First Loop to count */
  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){

    for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){
      isAlreadyAssignedRank1 [iSub2] = -1;
      isAlreadyAssignedRank2 [iSub2] = -1;
    }

    /* Loop on cell of current SubDomain */
    for(int icell = cellTileIdx[iSub  ];
            icell < cellTileIdx[iSub+1];
            icell++){
      for(int idxCell = CellCellIdx[icell  ];
              idxCell < CellCellIdx[icell+1];
              idxCell++ ){

        int ConnectedCell  = CellCell[idxCell]; //-1;
        int ConnectedColor = part->cellColor[ConnectedCell];

        /* Check if already assigned for Rank 2 */
        if(isAlreadyAssignedRank1[ConnectedColor] == -1)
        {
          isAlreadyAssignedRank1[ConnectedColor] = 1;
          ColorColorRank1Idx[iSub]++;
        }

        /* Check if already assigned for Rank 2 */
        if(isAlreadyAssignedRank2[ConnectedColor] == -1)
        {
          isAlreadyAssignedRank2[ConnectedColor] = 1;
          ColorColorCellRank2Idx[iSub]++;
        } /* End if already Treated */

        if(ConnectedColor != iSub)
        {
          /* Check if a cell neighbour can be shared by multiple subdmain */
          for(int idxCell2 = CellCellIdx[ConnectedCell  ];
                  idxCell2 < CellCellIdx[ConnectedCell+1];
                  idxCell2++ )
          {
            int ConnectedCell2  = CellCell[idxCell2];//-1;
            int ConnectedColor2 = part->cellColor[ConnectedCell2];

            // if(isAlreadyAssignedRank2[ConnectedColor2] == -1 && ConnectedColor2 < iSub)
            if(isAlreadyAssignedRank2[ConnectedColor2] == -1)
            {
              isAlreadyAssignedRank2[ConnectedColor2] = 1;
              ColorColorCellRank2Idx[iSub]++;
            }
          } /* End for second rank */
        }
      }
    }
    // printf("iSub : %i // R1 = %i // R2 = %i end  \n", iSub, ColorColorRank1Idx[iSub], ColorColorCellRank2Idx[iSub]);
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  _computeIdx(ColorColorRank1Idx    , nNewBlkCacheWanted);
  _computeIdx(ColorColorCellRank2Idx, nNewBlkCacheWanted);

  int *ColorColorRank1Arr     = (int *) malloc( ColorColorRank1Idx    [nNewBlkCacheWanted] * sizeof(int));
  int *ColorColorCellRank2Arr = (int *) malloc( ColorColorCellRank2Idx[nNewBlkCacheWanted] * sizeof(int));

  int ptRank1 = 0;
  int ptRank2 = 0;
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  /* Second Loop to fill */
  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){

    for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){
      isAlreadyAssignedRank1 [iSub2] = -1;
      isAlreadyAssignedRank2 [iSub2] = -1;
    }

    /* Loop on cell of current SubDomain */
    for(int icell = cellTileIdx[iSub  ];
            icell < cellTileIdx[iSub+1];
            icell++){
      for(int idxCell = CellCellIdx[icell  ];
              idxCell < CellCellIdx[icell+1];
              idxCell++ ){

        int ConnectedCell  = CellCell[idxCell];//-1;
        int ConnectedColor = part->cellColor[ConnectedCell];

        /* Check if already assigned for Rank 2 */
        if(isAlreadyAssignedRank1[ConnectedColor] == -1)
        {
          isAlreadyAssignedRank1[ConnectedColor] = 1;
          ColorColorRank1Arr[ptRank1] = ConnectedColor;
          ptRank1++;
        }

        /* Check if already assigned for Rank 2 */
        if(isAlreadyAssignedRank2[ConnectedColor] == -1)
        {
          isAlreadyAssignedRank2[ConnectedColor] = 1;
          ColorColorCellRank2Arr[ptRank2] = ConnectedColor;
          ptRank2++;
        } /* End if already Treated */

        if(ConnectedColor != iSub)
        {
          /* Check if a cell neighbour can be shared by multiple subdmain */
          for(int idxCell2 = CellCellIdx[ConnectedCell  ];
                  idxCell2 < CellCellIdx[ConnectedCell+1];
                  idxCell2++ )
          {
            int ConnectedCell2  = CellCell[idxCell2];//-1;
            int ConnectedColor2 = part->cellColor[ConnectedCell2];

            // if(isAlreadyAssignedRank2[ConnectedColor2] == -1 && ConnectedColor2 < iSub)
            if(isAlreadyAssignedRank2[ConnectedColor2] == -1)
            {
              isAlreadyAssignedRank2[ConnectedColor2] = 1;
              ColorColorCellRank2Arr[ptRank2] = ConnectedColor2;
              ptRank2++;
            }
          } /* End for second rank */
        }
      }
    }
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  if(0 == 1)
  {
    for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
      printf(" [%i] link to ", iSub);
      for(int ptColor = ColorColorCellRank2Idx[iSub];
              ptColor < ColorColorCellRank2Idx[iSub+1]; ptColor++)
      {
        printf("%i ", ColorColorCellRank2Arr[ptColor]);
      }
      printf("\n");
    }
  }

  /*
   * For debugging :
   * Le return reste ici pour debogger les ss-dom facilement dans un numérotaton :)
   */
  // return nNewBlkCacheWanted;

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  int* isLock        = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) );
  int* isTreatedSdom = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) );
  int* isAssigneSdom = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) );

  int* OldToNewColor = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) );
  int* saveSdom      = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) );

  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){isTreatedSdom[iSub]     = -1;}
  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){isAssigneSdom[iSub]     = -1;}

  if(part->threadColor == NULL){
    part->threadColor = (int *) malloc (sizeof(int) * nCell);
  }
  if(part->hyperPlaneColor == NULL){
    part->hyperPlaneColor = (int *) malloc (sizeof(int) * nCell);
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  int *flagHaveBnd     = (int *) malloc(nNewBlkCacheWanted * sizeof(int));

  _compute_flaghavebnd(part, flagHaveBnd, nNewBlkCacheWanted);

  int nSdomB = 0;                     /* Number of   Blocking domain */
  int nSdomU = 0;                     /* Number of Unblocking domain */

  for(int i = 0; i < nNewBlkCacheWanted; i++) {
    if(flagHaveBnd[i] == 0){ nSdomU++; }
    else                   { nSdomB++; }
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  int isLockCell[part->nCell];
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  // Resolution des couplages en priorité !!
  // nSubCouplingAsync = nb of Full Coupling + nb of Coupling 1to1
  int nSubCouplingAsync = (colorByTypeIdx[EndFullCoupling]-colorByTypeIdx[BeginFullCoupling])
                        + (colorByTypeIdx[EndCoupling1to1]-colorByTypeIdx[BeginCoupling1to1]);
  // printf("nSubCouplingAsync : %i \n", nSubCouplingAsync);
  int iSubTreated     = 0;
  int nThreadThreated = 0;
  int needReset       = 1;
  int lastCouplingThread = -1;
  int reminderC       = 0;
  int nPack           = 0;
  int nHp             = 0;
  int idxSave         =  0;


  int nSdomAsynchByThread[nThread]; // Number of asynchronous sub domain per thread
  int nSdomSynchByThread[nThread];  // Number of synchronous sub domain per thread

  int cptSdomAsynchByThread[nThread];
  int cptSdomSynchByThread[nThread];
  int cptSdom[nThread];

  for (int iThr = 0; iThr < nThread; iThr++){
    // nb of Asynchrone sub domain per thread
    nSdomAsynchByThread  [iThr] = colorByTypeIdx[BeginAsynchrone+iThr+1]-colorByTypeIdx[BeginAsynchrone+iThr];
    cptSdomAsynchByThread[iThr] = 0;
  }
  for (int iThr = 0; iThr < nThread; iThr++){
    nSdomSynchByThread  [iThr] = colorByTypeIdx[BeginSynchrone+iThr+1]-colorByTypeIdx[BeginSynchrone+iThr];
    cptSdomSynchByThread[iThr] = 0;
  }
  for (int iThr = 0; iThr < nThread; iThr++){
    cptSdom[iThr] = 0;
  }

  while(iSubTreated != nSubCouplingAsync+reminderC)
  {

    // printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n");
    // printf(" Loop on : %i \n", iSubTreated);
    if(needReset == 1)
    {
      for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){isLock[iSub2] = -1;}
      for (int icell = 0; icell < part->nCell;        icell++){isLockCell[icell] = -1;}

      for(int ptSdom = 0; ptSdom < nPack; ptSdom++)
      {
        int iSub = saveSdom[ptSdom];
        isTreatedSdom[iSub] = 1;
      }

      needReset = -1;
      lastCouplingThread = -1;

      nPack     = 0;
      idxSave   = 0;
    }

    /** Seek on coupling uniquely **/
    nThreadThreated = 0;
    while(nThreadThreated != nThread)
    {
      // printf(" Treatment of Thread : %i --> %i \n", nThreadThreated, lastCouplingThread);

      /* ---------------------------------------------------------------- */
      // Search candidate in Full Coupling or Coupling 1to1
      int iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginFullCoupling  ],
                                                  colorByTypeIdx[BeginFullCoupling+1],
                                                  nThreadThreated,
                                                  isAssigneSdom,
                                                  ColorColorCellRank2Idx,
                                                  ColorColorCellRank2Arr,
                                                  isLock,
                                                  flagHaveBnd,
                                                  0,
                                                  cellTileIdx,
                                                  isLockCell,
                                                  CellCellIdx,
                                                  CellCell);
      if(iSubCandidate == -1)
      {
        // printf("CouplingIdx[%i] = %i     ---> %i \n", nThreadThreated, CouplingIdx[nThreadThreated]);
        // printf("CouplingIdx[%i] = %i     ---> %i \n", nThreadThreated+1, CouplingIdx[nThreadThreated+1]);
        for(int ptThrC=CouplingIdx[nThreadThreated]; ptThrC < CouplingIdx[nThreadThreated+1]; ptThrC++){

          int idxSearchCoupling1to1 = BeginCoupling1to1 + ptThrC;
          // printf("ptThrC                ---> %i \n", ptThrC);
          // printf("BeginCoupling1to1     ---> %i \n", BeginCoupling1to1);
          // printf("idxSearchCoupling1to1 ---> %i \n", idxSearchCoupling1to1);
          // printf("Search in [%i/%i ] \n", colorByTypeIdx[idxSearchCoupling1to1], colorByTypeIdx[idxSearchCoupling1to1+1]);
          // printf("Search in Old [%i/%i ] \n", colorByTypeIdx[BeginFullCoupling], colorByTypeIdx[EndCoupling1to1]);
          iSubCandidate = _search_candidatesubdom(colorByTypeIdx[idxSearchCoupling1to1],
                                                  colorByTypeIdx[idxSearchCoupling1to1+1],
                                                  nThreadThreated,
                                                  isAssigneSdom,
                                                  ColorColorCellRank2Idx,
                                                  ColorColorCellRank2Arr,
                                                  isLock,
                                                  flagHaveBnd,
                                                  0,
                                                  cellTileIdx,
                                                  isLockCell,
                                                  CellCellIdx,
                                                  CellCell);
          // abort();
          if(iSubCandidate != -1){break;}
        }
      }
      /* ---------------------------------------------------------------- */

      /* ---------------------------------------------------------------- */
      if(iSubCandidate                != -1 &&
         isTreatedSdom[iSubCandidate] == -1)
      {
        OldToNewColor[iSubCandidate] = iSubTreated;
        _assign_candidatesubdom(iSubCandidate,
                                nThreadThreated,
                                isAssigneSdom,
                                ColorColorCellRank2Idx,
                                ColorColorCellRank2Arr,
                                isTreatedSdom,
                                isLock);

        for(int icell = cellTileIdx[iSubCandidate  ];
                icell < cellTileIdx[iSubCandidate+1];
                icell++)
        {
          part->threadColor[icell]     = nThreadThreated;
          part->hyperPlaneColor[icell] = nHp;
          isLockCell[icell] = nThreadThreated;

          for(int idxCell = CellCellIdx[icell  ]; idxCell < CellCellIdx[icell+1]; idxCell++ ){
            isLockCell[CellCell[idxCell]] = nThreadThreated;
          }
        }
        // printf("iSubTreated : %i ----> Th : %i : Candidate = %i \n", iSubTreated, nThreadThreated, iSubCandidate);
        saveSdom[idxSave++] = iSubCandidate;
        cptSdom[nThreadThreated]++;
        iSubTreated++;
        nPack++;
        lastCouplingThread = nThreadThreated;
      }
      else
      {
        // printf("Cannot find a good color max pack : %i \n", nPack);
        // abort();
        // Si il nous reste des Sdom couplage  on fait rien !!!
        // Sion on cherche des sous-dom  dans les inerieur (en priorite dans lee ss-dom qui traitent le couplage)
        /* ---------------------------------------------------------------- */
        iSubCandidate = -1;
        // Attention ici il faut boucler sur le dernier Thread qui a trouve un couplage !!!
        // printf("lastCouplingThread : %i \n", lastCouplingThread);

        if(lastCouplingThread != -1){
          for (int iThr = lastCouplingThread; iThr <= lastCouplingThread; iThr++){
             // printf(" iThr : %i ---> %i/%i \n", iThr, colorByTypeIdx[BeginAsynchrone+iThr], colorByTypeIdx[BeginAsynchrone+iThr+1]);
             // Steal Asynchrone sub domain to ensure parallelism
             iSubCandidate = _search_candidatesubdom_and_exclude(colorByTypeIdx[BeginAsynchrone+iThr],
                                                                 colorByTypeIdx[BeginAsynchrone+iThr+1],
                                                                 nThreadThreated,
                                                                 isAssigneSdom,
                                                                 ColorColorCellRank2Idx,
                                                                 ColorColorCellRank2Arr,
                                                                 isLock,
                                                                 flagHaveBnd,
                                                                 0,
                                                                 cellTileIdx,
                                                                 isLockCell,
                                                                 CellCellIdx,
                                                                 CellCell,
                                                                 colorByTypeIdx[EndCoupling1to1]
                                                                 );

             if(iSubCandidate != -1){break;}
           }
         }
         else
         {
           for (int iThr = nThreadThreated+1; iThr < nThread-1 ; iThr++){

             // printf(" iThr : %i ---> %i/%i \n", iThr, colorByTypeIdx[BeginAsynchrone+iThr], colorByTypeIdx[BeginAsynchrone+iThr+1]);
             // Steal Asynchrone sub domain to ensure parallelism
             iSubCandidate = _search_candidatesubdom_and_exclude(colorByTypeIdx[BeginAsynchrone+iThr],
                                                                 colorByTypeIdx[BeginAsynchrone+iThr+1],
                                                                 nThreadThreated,
                                                                 isAssigneSdom,
                                                                 ColorColorCellRank2Idx,
                                                                 ColorColorCellRank2Arr,
                                                                 isLock,
                                                                 flagHaveBnd,
                                                                 0,
                                                                 cellTileIdx,
                                                                 isLockCell,
                                                                 CellCellIdx,
                                                                 CellCell,
                                                                 colorByTypeIdx[EndCoupling1to1]
                                                                 );

             if(iSubCandidate != -1){lastCouplingThread = iThr; break;}
           }
         }
        /* ---------------------------------------------------------------- */

        // /* ---------------------------------------------------------------- */
        if(iSubCandidate                != -1 &&
           isTreatedSdom[iSubCandidate] == -1)
        {
          OldToNewColor[iSubCandidate] = iSubTreated;
          _assign_candidatesubdom(iSubCandidate,
                                  nThreadThreated,
                                  isAssigneSdom,
                                  ColorColorCellRank2Idx,
                                  ColorColorCellRank2Arr,
                                  isTreatedSdom,
                                  isLock);

          cptSdomAsynchByThread[lastCouplingThread]++;

          for(int icell = cellTileIdx[iSubCandidate  ];
                  icell < cellTileIdx[iSubCandidate+1];
                  icell++)
          {
            part->threadColor[icell]     = nThreadThreated;
            part->hyperPlaneColor[icell] = nHp;
            isLockCell[icell] = nThreadThreated;

            for(int idxCell = CellCellIdx[icell  ]; idxCell < CellCellIdx[icell+1]; idxCell++ ){
              isLockCell[CellCell[idxCell]] = nThreadThreated;
            }
          }
          /* ---------------------------------------------------------------- */
          // printf("iSubTreated : %i ----> Th : %i : Candidate = %i \n", iSubTreated, nThreadThreated, iSubCandidate);

          saveSdom[idxSave++] = iSubCandidate;
          iSubTreated++;
          reminderC++;
          nPack++;
          cptSdom[nThreadThreated]++;
          lastCouplingThread = nThreadThreated;
        }

        // printf("Find a possible interior SubDom : %i \n", iSubCandidate);
        // printf("Cannot find a good color max pack for thread : %i \n", nThreadThreated);
        if(iSubCandidate == -1)
          needReset = 1;

      }
      /* ---------------------------------------------------------------- */

      nThreadThreated++;
    }
    needReset = 1;
  }

  // printf("reminderC   : %i \n", reminderC);
  // printf("iSubTreated : %i \n", iSubTreated);
  // abort();

  if(0 == 1)
  {
    for (int iThr = 0; iThr < nThread; iThr++){
      printf(" cptSdomAsynchByThread[%i] = %i / %i \n", iThr, cptSdomAsynchByThread[iThr], nSdomAsynchByThread[iThr]);
    }
    for (int iThr = 0; iThr < nThread; iThr++){
      printf(" cptSdomSynchByThread[%i] = %i / %i \n", iThr, cptSdomSynchByThread[iThr], nSdomSynchByThread[iThr]);
    }
    for (int iThr = 0; iThr < nThread; iThr++){
      printf(" cptSdom[%i] = %i \n", iThr, cptSdom[iThr]);
    }
  }
  // abort();

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  // nThreadThreated = 0;
  needReset       = 1;
  nPack           = 0;
  // nHp             = -1;
  idxSave         =  0;
  int firstTry    =  1;
  int doBlocking  =  0;

  // int begT = 0;
  int endT = nThread;
  int steT = 1;
  // int nSubCouplingBndAsync = FAUX !!!colorByTypeIdx[nCoupling1to1+2]-colorByTypeIdx[nCoupling1to1+1];
  // printf("nSubCouplingBndAsync : %i \n", nSubCouplingBndAsync);
  // printf("nSdomU           : %i \n", nSdomU);
  // printf("nSdomU-reminderC : %i \n", nSdomU-reminderC);
  lastCouplingThread = -1;
  // printf("nNewBlkCacheWanted : %i // %i / %i \n", nNewBlkCacheWanted, nSdomU, nSdomB);

  // printf(" Begin Interior : %i \n", );

  while(iSubTreated != nNewBlkCacheWanted)
  {
    if(needReset == 1)
    {
      for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){isLock[iSub2] = -1;}
      for (int icell = 0; icell < part->nCell;        icell++){isLockCell[icell] = -1;}

      int rPack     = nPack%nThread;
      int nSdomPack = nPack/nThread;

      for(int ptSdom = 0; ptSdom < nPack; ptSdom++)
      {
        int iSub = saveSdom[ptSdom];
        isTreatedSdom[iSub] = 1;
      }

      if( (nSdomPack > 0 && rPack > 0) || ( firstTry == 1 && rPack != 0 )  )
      {
        firstTry = -1;
      }
      else
      {
        firstTry = 1;
      }

      needReset = -1;
      nPack     = 0;
      idxSave   = 0;
      nHp++;

      // printf("iSubTreated / iHp : %i // %i  \n", iSubTreated, nHp);
    }

    nThreadThreated = 0;

    if(iSubTreated >= nSdomU-reminderC && doBlocking == 0)
    {
      // printf(" Begin couronne !!!! \n");
      doBlocking = 2;
      // for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){flagHaveBnd[iSub2] = doBlocking;}
      for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){
        flagHaveBnd[iSub2] = 1;
      }
      // On reset TOUT
      needReset = 1;
      nThreadThreated = nThread;
    }
    else if(doBlocking == 2)
    {
      nThreadThreated = nThread-1;
      endT = -1;
      steT = -1;
    }
    // else if(iSubTreated == nSdomU+nSubCouplingBndAsync )
    // {
    //   printf(" Begin couronne purely Async !!!! \n");
    //   doBlocking = 1;
    // }

    while(nThreadThreated != endT)
    {
      // printf("iSubTreated : %i - %i \n", iSubTreated, firstTry);
      // On cherche un sous-domaine pas dependant
      // int isOkToContinue = -1;
      int iSubCandidate  = -1;

      // printf(" cptSdomAsynchByThread[%i] = %i / %i \n", nThreadThreated, cptSdomAsynchByThread[nThreadThreated], nSdomAsynchByThread[nThreadThreated]);
      // printf(" cptSdomAsynchByThread[%i] = %i / %i \n", nThreadThreated, cptSdomSynchByThread[nThreadThreated], nSdomSynchByThread[nThreadThreated]);
      /* ---------------------------------------------------------------- */
      // printf(" \t Search in [%i/%i] \n", colorByTypeIdx[nThreadThreated+nThread+nCoupling1to1+1], colorByTypeIdx[nThreadThreated+nThread+nCoupling1to1+2]);
      if(doBlocking == 0)
      {
        // Asynchrone part
        iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginAsynchrone+nThreadThreated],
                                                colorByTypeIdx[BeginAsynchrone+nThreadThreated+1],
                                                nThreadThreated,
                                                isAssigneSdom,
                                                ColorColorCellRank2Idx,
                                                ColorColorCellRank2Arr,
                                                isLock,
                                                flagHaveBnd,
                                                doBlocking,
                                                cellTileIdx,
                                                isLockCell,
                                                CellCellIdx,
                                                CellCell);

        if(iSubCandidate != -1){cptSdomAsynchByThread[nThreadThreated]++;}
      }
      else if(doBlocking == 2)
      {
        if(cptSdomAsynchByThread[nThreadThreated] != nSdomAsynchByThread[nThreadThreated])
        {
          // printf(" Finish Asynchrone Sub as Synchrone : %i \n", nThreadThreated);
          iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginAsynchrone+nThreadThreated],
                                                  colorByTypeIdx[BeginAsynchrone+nThreadThreated+1],
                                                  nThreadThreated,
                                                  isAssigneSdom,
                                                  ColorColorCellRank2Idx,
                                                  ColorColorCellRank2Arr,
                                                  isLock,
                                                  flagHaveBnd,
                                                  1,
                                                  cellTileIdx,
                                                  isLockCell,
                                                  CellCellIdx,
                                                  CellCell);

          if(iSubCandidate != -1){cptSdomAsynchByThread[nThreadThreated]++;}
        }
        else
        {
          // printf(" Couronne : %i \n", nThreadThreated);
          iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginFullCouplingC],
                                                  colorByTypeIdx[BeginFullCouplingC+1],
                                                  nThreadThreated,
                                                  isAssigneSdom,
                                                  ColorColorCellRank2Idx,
                                                  ColorColorCellRank2Arr,
                                                  isLock,
                                                  flagHaveBnd,
                                                  1,
                                                  cellTileIdx,
                                                  isLockCell,
                                                  CellCellIdx,
                                                  CellCell);

          if(iSubCandidate == -1)
          {
            // printf("CouplingIdx[%i] = %i     ---> %i \n", nThreadThreated, CouplingIdx[nThreadThreated]);
            // printf("CouplingIdx[%i] = %i     ---> %i \n", nThreadThreated+1, CouplingIdx[nThreadThreated+1]);
            for(int ptThrC=CouplingIdx[nThreadThreated]; ptThrC < CouplingIdx[nThreadThreated+1]; ptThrC++){

              int idxSearchCoupling1to1C = BeginCoupling1to1C + ptThrC;
              // printf("ptThrC                ---> %i \n", ptThrC);
              // printf("BeginCoupling1to1     ---> %i \n", BeginCoupling1to1);
              // printf("idxSearchCoupling1to1 ---> %i \n", idxSearchCoupling1to1);
              // printf("Search in [%i/%i ] \n", colorByTypeIdx[idxSearchCoupling1to1], colorByTypeIdx[idxSearchCoupling1to1+1]);
              // printf("Search in Old [%i/%i ] \n", colorByTypeIdx[BeginFullCoupling], colorByTypeIdx[EndCoupling1to1]);
              iSubCandidate = _search_candidatesubdom(colorByTypeIdx[idxSearchCoupling1to1C],
                                                      colorByTypeIdx[idxSearchCoupling1to1C+1],
                                                      nThreadThreated,
                                                      isAssigneSdom,
                                                      ColorColorCellRank2Idx,
                                                      ColorColorCellRank2Arr,
                                                      isLock,
                                                      flagHaveBnd,
                                                      1,
                                                      cellTileIdx,
                                                      isLockCell,
                                                      CellCellIdx,
                                                      CellCell);
              // abort();
              if(iSubCandidate != -1){break;}
            }
          }


          if(iSubCandidate == -1)
          {
            // printf(" Couronne fail so try asynchrone part : %i --> %i \n", nThreadThreated, lastCouplingThread);
            iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginSynchrone+lastCouplingThread],
                                                    colorByTypeIdx[BeginSynchrone+lastCouplingThread+1],
                                                    // colorByTypeIdx[lastCouplingThread+nCoupling1to1+2],
                                                    // colorByTypeIdx[lastCouplingThread+nCoupling1to1+3],
                                                    nThreadThreated,
                                                    isAssigneSdom,
                                                    ColorColorCellRank2Idx,
                                                    ColorColorCellRank2Arr,
                                                    isLock,
                                                    flagHaveBnd,
                                                    1,
                                                    cellTileIdx,
                                                    isLockCell,
                                                    CellCellIdx,
                                                    CellCell);
            if(iSubCandidate == -1)
            {
              // printf(" Couronne fail so try asynchrone part (last synchrone): %i --> %i \n", nThreadThreated, lastCouplingThread);
              iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginSynchrone+nThreadThreated],
                                                      colorByTypeIdx[BeginSynchrone+nThreadThreated+1],
                                                      nThreadThreated,
                                                      isAssigneSdom,
                                                      ColorColorCellRank2Idx,
                                                      ColorColorCellRank2Arr,
                                                      isLock,
                                                      flagHaveBnd,
                                                      1,
                                                      cellTileIdx,
                                                      isLockCell,
                                                      CellCellIdx,
                                                      CellCell);
            }
          }
          if(iSubCandidate != -1){cptSdomSynchByThread[nThreadThreated]++;}
        }
      }
      else
      {
        iSubCandidate = _search_candidatesubdom(colorByTypeIdx[BeginSynchrone+nThreadThreated],
                                                colorByTypeIdx[BeginSynchrone+nThreadThreated+1],
                                                nThreadThreated,
                                                isAssigneSdom,
                                                ColorColorCellRank2Idx,
                                                ColorColorCellRank2Arr,
                                                isLock,
                                                flagHaveBnd,
                                                doBlocking,
                                                cellTileIdx,
                                                isLockCell,
                                                CellCellIdx,
                                                CellCell);
      }
      /* ---------------------------------------------------------------- */

      if(iSubCandidate                != -1 &&
         isTreatedSdom[iSubCandidate] == -1)
      {
        OldToNewColor[iSubCandidate] = iSubTreated;

        // if(doBlocking == 0)
        // {
        //   cptSdomAsynchByThread[nThreadThreated]++;
        // }
        // else if(doBlocking == 1)
        // {
        //   cptSdomSynchByThread[nThreadThreated]++;
        // }

        _assign_candidatesubdom(iSubCandidate,
                                nThreadThreated,
                                isAssigneSdom,
                                ColorColorCellRank2Idx,
                                ColorColorCellRank2Arr,
                                isTreatedSdom,
                                isLock);

        for(int icell = cellTileIdx[iSubCandidate  ];
                icell < cellTileIdx[iSubCandidate+1];
                icell++)
        {
          part->threadColor[icell]     = nThreadThreated;
          part->hyperPlaneColor[icell] = nHp;
          isLockCell[icell] = nThreadThreated;

          for(int idxCell = CellCellIdx[icell  ]; idxCell < CellCellIdx[icell+1]; idxCell++ ){
            isLockCell[CellCell[idxCell]] = nThreadThreated;
          }
        }
        // printf("iSubTreated : %i ----> Th : %i : Candidate = %i \n", iSubTreated, nThreadThreated, iSubCandidate);

        saveSdom[idxSave++] = iSubCandidate;
        cptSdom[nThreadThreated]++;
        iSubTreated++;
        nPack++;
        lastCouplingThread = nThreadThreated;
      }
      else
      {
        // printf("Cannot find a good color max pack : %i \n", nPack);
        needReset = 1;
      }
      // nThreadThreated++;
      nThreadThreated += steT;

    }

    // if(nPack >= nThread){needReset = 1;}

  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  // for (int iThr = 0; iThr < nThread; iThr++){
  //   printf(" cptSdomAsynchByThread[%i] = %i / %i \n", iThr, cptSdomAsynchByThread[iThr], nSdomAsynchByThread[iThr]);
  // }
  // for (int iThr = 0; iThr < nThread; iThr++){
  //   printf(" cptSdomSynchByThread[%i] = %i / %i \n", iThr, cptSdomSynchByThread[iThr], nSdomSynchByThread[iThr]);
  // }

  if(0 == 1)
  {
    for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){
      printf("OldToNewColor[%i] = %i \n", iSub2,  OldToNewColor[iSub2]);
    }
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */


  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  int *newToOldOrder    = (int *) malloc (nNewBlkCacheWanted * sizeof(int));
  int *newToOldOrderTmp = (int *) malloc (nNewBlkCacheWanted * sizeof(int));
  for(int i = 0; i < nNewBlkCacheWanted; i++) {
    newToOldOrder[OldToNewColor[i]] = i;
    newToOldOrderTmp[OldToNewColor[i]] = i;
  }

  nHp++;
  // printf(" nHp : %i \n", nHp);
  int sdomByHpIdx[nHp+1];
  for (int iHp = 0; iHp < nHp+1; iHp++){
    sdomByHpIdx[iHp] = 0;
  }

  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
    // printf("OldToNewColor[%i] = %i \n", iSub,  OldToNewColor[iSub]);
    // printf("newToOldOrder[%i] = %i \n", iSub,  newToOldOrder[iSub]);
    sdomByHpIdx[part->hyperPlaneColor[cellTileIdx[newToOldOrder[iSub]]]]++;
    // printf("part->hyperPlaneColor[%i] = %i \n", iSub, part->hyperPlaneColor[cellTileIdx[iSub]]);
    // printf("part->hyperPlaneColor[%i] = %i \n", iSub, part->hyperPlaneColor[cellTileIdx[newToOldOrder[iSub]]]);
    // printf("part->threadColor[%i] = %i \n", iSub, part->threadColor[cellTileIdx[iSub]]);
  }

  _computeIdx(sdomByHpIdx, nHp);

  if(0 == 1)
  {
    for (int iHp = 0; iHp < nHp+1; iHp++){
      printf("sdomByHpIdx[%i] = %i \n", iHp, sdomByHpIdx[iHp]);
    }
  }

  int nSdomByThreadInHpIdx[nThread+1];
  int cptSdomByThreadInHpIdx[nThread+1];
  for (int iHp = 0; iHp < nHp; iHp++){

    // printf(" [%i] Hp is compose with : \n", iHp );
    // for (int iSub = sdomByHpIdx[iHp]; iSub < sdomByHpIdx[iHp+1]; iSub++){
    //   printf(" %i ", iSub);
    // }
    // printf("\n");

    for (int iThr = 0; iThr < nThread+1; iThr++){
      nSdomByThreadInHpIdx[iThr]   = 0;
      cptSdomByThreadInHpIdx[iThr] = 0;
    }

    // printf(" [%i] Hp is compose with : \n", iHp );
    for (int iSub = sdomByHpIdx[iHp]; iSub < sdomByHpIdx[iHp+1]; iSub++){
      int iSubO = newToOldOrder[iSub];
      // printf("iSub -> %i | iThr -> %i \n ", iSub, part->threadColor[cellTileIdx[iSubO]]);
      nSdomByThreadInHpIdx[ part->threadColor[cellTileIdx[iSubO]] ]++;
    }
    _computeIdx(nSdomByThreadInHpIdx, nThread);

    if(0 == 1)
    {
      for (int iThr = 0; iThr < nThread+1; iThr++){
        printf("nSdomByThreadInHpIdx[%i] = %i \n ", iThr, nSdomByThreadInHpIdx[iThr]);
      }
    }

    for (int iSub = sdomByHpIdx[iHp]; iSub < sdomByHpIdx[iHp+1]; iSub++){
      int iSubO = newToOldOrder[iSub];
      int iThr  = part->threadColor[cellTileIdx[iSubO]];

      int idx = nSdomByThreadInHpIdx[iThr]+cptSdomByThreadInHpIdx[iThr];
      // newToOldOrderTmp[sdomByHpIdx[iHp]+idx] = iSub;
      newToOldOrderTmp[sdomByHpIdx[iHp]+idx] = iSubO;
      cptSdomByThreadInHpIdx[iThr]++;
    }



  }

  for(int i = 0; i < nNewBlkCacheWanted; i++) {
    OldToNewColor[newToOldOrderTmp[i]] = i;
  }

  if(0 == 1)
  {
    for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
      printf("OldToNewColor[%i] = %i \n", iSub,  OldToNewColor[iSub]);
      printf("newToOldOrder[%i] = %i \n", iSub,  newToOldOrder[iSub]);
      printf("newToOldOrderTmp[%i] = %i \n", iSub,  newToOldOrderTmp[iSub]);
      printf("part->hyperPlaneColor[%i] = %i \n", iSub, part->hyperPlaneColor[cellTileIdx[iSub]]);
      printf("part->hyperPlaneColor[%i] = %i \n", iSub, part->hyperPlaneColor[cellTileIdx[newToOldOrder[iSub]]]);
      printf("part->threadColor[%i] = %i \n", iSub, part->threadColor[cellTileIdx[iSub]]);
    }
  }

  free(newToOldOrder);
  free(newToOldOrderTmp);
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */


  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  int *partCellNewIdx      = (int *) malloc((nNewBlkCacheWanted + 1) * sizeof(int));

  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
    nCellPerCache[iSub] = 0;
  }
  for (int i = 0; i < nCell; i++){
    int color = OldToNewColor[part->cellColor[i]];
    // printf(" cellColor[%i] = %i -> %i \n", i, part->cellColor[i], color );
    nCellPerCache[color]++;
    CellOrder[i] = -1;
  }
    /* Compute the index of each Sub-Domain **/
  partCellNewIdx[0] = 0;
  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
    partCellNewIdx[iSub + 1] = partCellNewIdx[iSub] + nCellPerCache[iSub];
  }

  for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){
    int subNew = OldToNewColor[iSub];
    int begNew = partCellNewIdx[subNew];
    int nNew   = partCellNewIdx[subNew+1] - partCellNewIdx[subNew];
    int nOld   = cellTileIdx[iSub+1] - cellTileIdx[iSub];
    // printf(" nNew = %i \n", nNew);
    // printf(" nOld = %i \n", nOld);
    assert(nNew == nOld);
    for(int icell=cellTileIdx[iSub]; icell < cellTileIdx[iSub+1]; icell++)
    {
      CellOrder[begNew] = icell;
      begNew++;
    }
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */



  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  if(0 == 1)
  {
    int* verif          = (int * ) malloc( part->nCell * sizeof(int) );
    for (int i = 0; i < part->nCell; i++){
      verif[i] = 0;
    }
    for (int i = 0; i < part->nCell; i++){
      assert(CellOrder[i] != -1);
      verif[CellOrder[i]]++;
      // printf("CellOrder[%i] = %i \n", i,  CellOrder[i]);
    }
    for (int i = 0; i < part->nCell; i++){
      assert(verif[i] == 1);
    }
    free(verif);
  }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  PDM_part_reorder_cell(part, CellOrder);

  for(int i = 0; i < part->nCell; i++) {
    part->cellColor[i] = OldToNewColor[part->cellColor[i]];
  }

  if(0 == 1)
  {
    for (int i = 0; i < nCell; i++){
      printf(" cellColor[%i] = %i \n", i, part->cellColor[i]);
    }
  }
  // for (int i = 0; i < nCell; i++){
  //   printf("part->hyperPlaneColor[%i] = %i \n", i, part->hyperPlaneColor[i]);
  // }
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
  free(partCellNewIdx         );
  free(ColorColorIdx          );
  free(ColorColorRank1Idx     );
  free(ColorColorCellRank2Idx );
  free(isAlreadyAssignedRank1 );
  free(isAlreadyAssignedRank2 );
  free(ColorColorRank1Arr );
  free(ColorColorCellRank2Arr );
  free(cellTileIdx );
  free(OldToNewColor );
  free(isAssigneSdom );
  free(isTreatedSdom );
  free(isLock );
  free(saveSdom );
  free(flagHaveBnd );
  free(CouplingIdx);
  free(CouplingArr);
  free(Coupling   );
  free(nCellPerCache);
  free(CellCellIdx);
  free(CellCell);
  // free(ColorSub);
  /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

  return nNewBlkCacheWanted;

}


/**
 *
 * \brief Compute OpenMP renumbering constraint
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      CellOrder      new Cell Order
 *
 */
/* static int */
/* _compute_openMP_renumbering */
/* ( */
/*  _part_t   *part, */
/*  int        split_method, */
/*  int        nCellPerCacheWanted, */
/*  int       *CellOrder, */
/*  int       *CellCellIdx, */
/*  int       *CellCell */
/* ) */
/* { */
/*   /\* Get nFac and nCel *\/ */
/*   const int nCell = part->nCell; */

/*   int nThread; */
/*   /\* Split the graph *\/ */
/* #ifdef PDM_HAVE_OPENMP */
/* #pragma omp parallel default(shared) */
/*   { */
/*     #pragma omp master */
/*     { */
/*      nThread = omp_get_num_threads(); */
/*     } */
/*   } */
/* #else */
/*   nThread = 1; */
/* #endif */
/*   /\* Compute graph associate to mesh *\/ */
/*   PDM_part_graph_compute_from_face_cell(          part, */
/*                                         (int **) &CellCellIdx, */
/*                                         (int **) &CellCell); */

/*   /\* */
/*    * Determine the optimal size for cache blocking */
/*    *   -> The user specify the desired number of cells on each subdomain */
/*    *   -> An other idea is to specify directly the number of block he want */
/*    *\/ */
/*   int nBlkCacheWanted; */
/*   if(nCellPerCacheWanted == 0){nBlkCacheWanted = 1;} */
/*   else                        {nBlkCacheWanted = PDM_MAX(nCell/nCellPerCacheWanted,1);} */

/*   /\* Split the graph *\/ */
/*   PDM_part_graph_split(          split_method, */
/*                                  nThread, */
/*                                  part, */
/*                                  CellCellIdx, */
/*                                  CellCell, */
/*                        (int *)   NULL, */
/*                        (int *)   NULL, */
/*                        (int **) &part->threadColor); */

/*   /\* Allocate *\/ */
/*   int *nCellPerThread    = (int *) malloc( nThread      * sizeof(int)); */
/*   int *nCellPerThreadBeg = (int *) malloc( nThread      * sizeof(int)); */
/*   int *threadCellIdx     = (int *) malloc((nThread + 1) * sizeof(int)); */


/*   if(0 == 1) */
/*   { */
/*     printf(" ----------- nThread : %i \n", nThread); */
/*     for (int i = 0; i < nCell; i++){ */
/*       printf(" =====> part->threadColor[%i]  = %i\n", i, part->threadColor[i]); */
/*     } */
/*   } */

/*   /\* Init to Zero *\/ */
/*   for (int i = 0; i < nThread; i++){ */
/*     nCellPerThread[i] = 0; */
/*   } */

/*   /\* First loop to identify the size of each color *\/ */
/*   for (int i = 0; i < nCell; i++){ */
/*     int color = part->threadColor[i]; */
/*     nCellPerThread[color]++; */
/*   } */

/*   /\* Compute the index of each Sub-Domain **\/ */
/*   threadCellIdx[0] = 0; */
/*   for (int i = 0; i < nThread; i++){ */
/*     threadCellIdx[i + 1] = threadCellIdx[i] + nCellPerThread[i]; */
/*   } */


/*   /\* Reset array *\/ */
/*   for (int i = 0; i < nThread; i++){ */
/*     nCellPerThread[i] = 0; */
/*   } */

/*   for (int i = 0; i < nCell; i++){ */
/*     int color  = part->threadColor[i]; */
/*     int idx    = threadCellIdx[color] + nCellPerThread[color]; */

/*     CellOrder[idx] = i; */
/*     nCellPerThread[color]++; */
/*   } */

/*   PDM_part_reorder_cell(part, CellOrder); */


/*   if(0 == 1) */
/*   { */
/*     printf(" ----------- nThread : %i \n", nThread); */
/*     for (int i = 0; i < nCell; i++){ */
/*       printf(" =====> part->threadColor[%i]  = %i\n", i, part->threadColor[i]); */
/*     } */
/*   } */

/*   /\* */
/*    * Contruction of Thread coupling subDomain */
/*    *\/ */
/*   free(CellCellIdx); */
/*   free(CellCell); */
/*   CellCellIdx = NULL; */
/*   CellCell = NULL; */
/*   PDM_part_graph_compute_from_face_cell(          part, */
/*                                         (int **) &CellCellIdx, */
/*                                         (int **) &CellCell); */

/*   int *isTreatedCell  = (int *) malloc( nCell  * sizeof(int)); */

/*   for (int icell = 0; icell < nCell; icell++){isTreatedCell[icell] = -1;} */


  /* // int **nCellByPoolByThread = (int *) malloc( nThread      * sizeof(int)); */
  /* int nCellByPoolByThread  [nThread][nThread]; */
  /* int nCellByPoolByThreadR2[nThread][nThread]; */

  /* /\* ----------------------------------------------------------- *\/ */
  /* for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*   for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*     nCellByPoolByThread[iThr][iThr2] = 0; */
  /*     nCellByPoolByThreadR2[iThr][iThr2] = 0; */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* /\* ----------------------------------------------------------- *\/ */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = threadCellIdx[iThr  ]; */
  /*           icell < threadCellIdx[iThr+1]; */
  /*           icell++) */
  /*   { */
  /*     int connectedToCurrentThread =  1; */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; */
  /*             idxCell++ ) */
  /*     { */
  /*       int ConnectedCell   = CellCell[idxCell]; */
  /*       int ConnectedThread = part->threadColor[ConnectedCell]; */

  /*       if( ConnectedThread           != iThr){ */
  /*         if(connectedToCurrentThread != -1){ */

  /*           nCellByPoolByThread[iThr][ConnectedThread]++; */
  /*           connectedToCurrentThread = -1; */
  /*           isTreatedCell[icell] = 1; */
  /*         } */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */


  /* /\* ----------------------------------------------------------- *\/ */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = threadCellIdx[iThr  ]; */
  /*           icell < threadCellIdx[iThr+1]; */
  /*           icell++) */
  /*   { */
  /*     if(isTreatedCell[icell] == -1) */
  /*     { */
  /*       int connectedToCurrentThread =  1; */
  /*       for(int idxCell = CellCellIdx[icell  ]; */
  /*               idxCell < CellCellIdx[icell+1]; */
  /*               idxCell++ ) */
  /*       { */
  /*         int ConnectedCell   = CellCell[idxCell]; */
  /*         int ConnectedThread = part->threadColor[ConnectedCell]; */

  /*         if(isTreatedCell[ConnectedCell] == 1 && connectedToCurrentThread != -1){ */
  /*           nCellByPoolByThreadR2[iThr][ConnectedThread]++; */
  /*           connectedToCurrentThread = -1; */
  /*         } */
  /*       } */
  /*       if(connectedToCurrentThread == 1) */
  /*       { */
  /*         nCellByPoolByThread[iThr][iThr]++; */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* int nCellCheck = 0; */
  /* for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*   for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*     nCellCheck +=  nCellByPoolByThread[iThr][iThr2]; */
  /*     nCellCheck +=  nCellByPoolByThreadR2[iThr][iThr2]; */
  /*     printf(" nCellByPoolByThread    [%i][%i] = %i \n", iThr, iThr2, nCellByPoolByThread[iThr][iThr2]); */
  /*     printf(" nCellByPoolByThreadR2  [%i][%i] = %i \n", iThr, iThr2, nCellByPoolByThreadR2[iThr][iThr2]); */
  /*   } */
  /* } */

  /* // int nInt = */

  /* assert(nCellCheck == nCell ); */

  /* printf("nCell : %i \n", nCell); */
  /* printf("nCellCheck : %i \n", nCellCheck); */

  /* /\* ----------------------------------------------------------- *\/ */
  /* int begCellByPoolByThread  [nThread][nThread]; */
  /* int begCellByPoolByThreadR2[nThread][nThread]; */

  /* for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*   for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*     begCellByPoolByThread[iThr][iThr2] = threadCellIdx[iThr]; */
  /*   } */
  /* } */

  /* for (int iThr  = 0; iThr < nThread; iThr++ ){ */
  /*   begCellByPoolByThread[iThr][nThread-1] = threadCellIdx[iThr]; */
  /*   for (int iThr2 = nThread-2; iThr2 > -1; iThr2--){ */
  /*     begCellByPoolByThread[iThr][iThr2] = begCellByPoolByThread[iThr][iThr2+1] + nCellByPoolByThread[iThr][iThr2+1]; */
  /*   } */
  /* } */


  /* for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*   for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*     if (iThr2<=iThr){ */
  /*       begCellByPoolByThreadR2[iThr][iThr2] = threadCellIdx[iThr+1] - nCellByPoolByThreadR2[iThr][iThr]; */
  /*     } */
  /*     else{ */
  /*       begCellByPoolByThreadR2[iThr][iThr2] = threadCellIdx[iThr+1]; */
  /*     } */
  /*   } */
  /* } */

  /* /\* ----------------------------------------------------------- *\/ */


  /* if(1 == 1) */
  /* { */
  /*   /\* ----------------------------------------------------------- *\/ */
  /*   printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /*   for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*     for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*       printf(" nCellByPoolByThread  [%i][%i] = %i \n", iThr, iThr2, nCellByPoolByThread[iThr][iThr2]); */
  /*     } */
  /*   } */
  /*   printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /*   /\* ----------------------------------------------------------- *\/ */

  /*   /\* ----------------------------------------------------------- *\/ */
  /*   printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /*   for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*     for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*       printf(" begCellByPoolByThread[%i][%i] = %i \n", iThr, iThr2, begCellByPoolByThread[iThr][iThr2]); */
  /*     } */
  /*   } */
  /*   printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /*   /\* ----------------------------------------------------------- *\/ */

  /*   /\* ----------------------------------------------------------- *\/ */
  /*   printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /*   for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*     for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*       printf(" begCellByPoolByThreadR2[%i][%i] = %i \n", iThr, iThr2, begCellByPoolByThreadR2[iThr][iThr2]); */
  /*     } */
  /*   } */
  /*   printf(" oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /*   /\* ----------------------------------------------------------- *\/ */
  /* } */

  /* /\* ----------------------------------------------------------- *\/ */
  /* /\* Reset *\/ */
  /* for (int iThr  = 0; iThr  < nThread; iThr++ ){ */
  /*   for (int iThr2 = 0; iThr2 < nThread; iThr2++){ */
  /*     nCellByPoolByThread  [iThr][iThr2] = 0; */
  /*     nCellByPoolByThreadR2[iThr][iThr2] = 0; */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* for (int icell = 0; icell < nCell; icell++){isTreatedCell[icell] = -1;} */

  /* /\* ----------------------------------------------------------- *\/ */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = threadCellIdx[iThr  ]; */
  /*           icell < threadCellIdx[iThr+1]; */
  /*           icell++) */
  /*   { */
  /*     int connectedToCurrentThread =  1; */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; */
  /*             idxCell++ ) */
  /*     { */
  /*       int ConnectedCell   = CellCell[idxCell]; */
  /*       int ConnectedThread = part->threadColor[ConnectedCell]; */

  /*       if( ConnectedThread           != iThr){ */
  /*         if(connectedToCurrentThread != -1){ */

  /*           int idxBeg = begCellByPoolByThread[iThr][ConnectedThread]+nCellByPoolByThread[iThr][ConnectedThread]; */
  /*           CellOrder[idxBeg] = icell; */

  /*           nCellByPoolByThread[iThr][ConnectedThread]++; */
  /*           connectedToCurrentThread = -1; */
  /*           isTreatedCell[icell] = 1; */
  /*         } */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* /\* ----------------------------------------------------------- *\/ */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = threadCellIdx[iThr  ]; */
  /*           icell < threadCellIdx[iThr+1]; */
  /*           icell++) */
  /*   { */
  /*     if(isTreatedCell[icell] == -1) */
  /*     { */
  /*       int connectedToCurrentThread =  1; */
  /*       for(int idxCell = CellCellIdx[icell  ]; */
  /*               idxCell < CellCellIdx[icell+1]; */
  /*               idxCell++ ) */
  /*       { */
  /*         int ConnectedCell   = CellCell[idxCell]; */
  /*         int ConnectedThread = part->threadColor[ConnectedCell]; */

  /*         if(isTreatedCell[ConnectedCell] == 1 && connectedToCurrentThread != -1){ */

  /*           int idxBeg = begCellByPoolByThreadR2[iThr][ConnectedThread]+nCellByPoolByThreadR2[iThr][ConnectedThread]; */
  /*           CellOrder[idxBeg] = icell; */

  /*           nCellByPoolByThreadR2[iThr][ConnectedThread]++; */
  /*           connectedToCurrentThread = -1; */
  /*         } */
  /*       } */
  /*       if(connectedToCurrentThread == 1) */
  /*       { */
  /*         int idxBeg = begCellByPoolByThread[iThr][iThr]+nCellByPoolByThread[iThr][iThr]; */
  /*         CellOrder[idxBeg] = icell; */

  /*         nCellByPoolByThread[iThr][iThr]++; */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* /\* */
  /*  * Verbose */
  /*  *\/ */
  /* if(1 == 1) */
  /* { */
  /*   printf(" ----------- nBlkCacheWanted : %i \n", nBlkCacheWanted); */
  /*   int check[nCell]; */
  /*   for (int i = 0; i < nCell; i++){ */
  /*     check[i] = 0; */
  /*   } */

  /*   for (int i = 0; i < nCell; i++){ */
  /*     // printf("~> %i\n", i); */
  /*     // printf(" =====> CellOrder[%i]  = %i\n", i, CellOrder[i]); */
  /*     check[CellOrder[i]]++; */
  /*   } */

  /*   for (int i = 0; i < nCell; i++){ */
  /*     // printf("Check[%i] = %i \n", i, check[i]); */
  /*     assert(check[i] == 1); */
  /*   } */
  /* } */

  /* /\* ----------------------------------------------------------- *\/ */
  /* PDM_part_reorder_cell(part, CellOrder); */
  /* /\* ----------------------------------------------------------- *\/ */

  /* if(0 == 1) */
  /* { */
  /*   printf(" vvvvvvvvv  ----------- nThread : %i \n", nThread); */
  /*   for (int i = 0; i < nCell; i++){ */
  /*     printf(" =====> part->threadColor[%i]  = %i\n", i, part->threadColor[i]); */
  /*   } */
  /* } */


  /* /\* ----------------------------------------------------------- *\/ */
  /* free(CellCellIdx); */
  /* free(CellCell); */
  /* CellCellIdx = NULL; */
  /* CellCell = NULL; */
  /* PDM_part_graph_compute_from_face_cell(          part, */
  /*                                       (int **) &CellCellIdx, */
  /*                                       (int **) &CellCell); */
  /* /\* ----------------------------------------------------------- *\/ */

  /* int* CellCellSubIdx = (int *) malloc( nCell              * sizeof(int) ); */
  /* int* CellCellSub    = (int *) malloc( CellCellIdx[nCell] * sizeof(int) ); */

  /* part->cellColor = (int * ) malloc( nCell              * sizeof(int) ); */

  /* for (int i = 0; i < nCell; i++){part->cellColor[i] = -1;}; */

  /* /\* ----------------------------------------------------------- *\/ */
  /* int offsetColor = 0; */
  /* int nCellVerif = 0; */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */

  /* for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /*   int nC   = nCellByPoolByThreadR2[iThr][iThrC]; */
  /*   int BegC = begCellByPoolByThreadR2[iThr][iThrC]; */
  /*   int EndC = begCellByPoolByThreadR2[iThr][iThrC]+nC; */

  /*   nCellVerif += nC; */

  /*   /\* Init *\/ */
  /*   for (int icell = 0; icell < nC+1; icell++){CellCellSubIdx[icell] = 0;} */

  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; idxCell++ ){ */

  /*       int ConnectedCell   = CellCell[idxCell]; */

  /*       if(ConnectedCell >= BegC && ConnectedCell < EndC ) */
  /*       { */
  /*         CellCellSubIdx[icell-BegC]++; */
  /*       } */
  /*     } */
  /*   } */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   _computeIdx(CellCellSubIdx, nC); */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */


  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   int ptArr = 0; */
  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; idxCell++ ){ */

  /*       int ConnectedCell = CellCell[idxCell]; */

  /*       if(ConnectedCell >= BegC && ConnectedCell < EndC ) */
  /*       { */
  /*         CellCellSub[ptArr] = CellCell[idxCell]-BegC; */
  /*         // printf(" [%i/%i] - CellCellSub[%i] = %i \n", iThr, iThrC, ptArr, CellCell[idxCell]-BegC); */
  /*         ptArr++; */
  /*       } */
  /*     } */
  /*   } */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   int nBlkCacheWantedSub; */
  /*   if(nCellPerCacheWanted == 0){nBlkCacheWantedSub = 1;} */
  /*   else                        {nBlkCacheWantedSub = PDM_MAX(nC/nCellPerCacheWanted,1);} */

  /*   // nBlkCacheWantedSub = 1; */
  /*   /\* Split the graph *\/ */
  /*   int* ColorSub = NULL; */
  /*   printf("[%i/%i] PDM_part_graph_split_bis : %i \n", iThr, iThrC, nC ); */
  /*   if(nC != 0) */
  /*   { */
  /*     PDM_part_graph_split_bis(          split_method, */
  /*                                        nBlkCacheWantedSub, */
  /*                                        nC, */
  /*                                        CellCellSubIdx, */
  /*                                        CellCellSub, */
  /*                              (int *)   NULL, */
  /*                              (int *)   NULL, */
  /*                                       &ColorSub); */
  /*   } */
  /*   // printf("[%i/%i] PDM_part_graph_split_bis end \n", iThr, iThrC); */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /*   if(0 == 1) */
  /*   { */
  /*     printf(" ----------- ColorSub : %i \n", nC); */
  /*     for (int i = 0; i < nC; i++){ */
  /*       printf(" =====> ColorSub[%i]  = %i\n", i, ColorSub[i]); */
  /*     } */
  /*   } */

  /*   int offsetColorTmp = offsetColor; */
  /*   for (int icell = 0; icell < nC; icell++){ */
  /*     int newColor =  offsetColorTmp+ColorSub[icell]; */
  /*     part->cellColor[BegC+icell] = newColor; */
  /*     // printf(" =====> cellColor[%i]  = %i\n", BegC+icell, newColor); */
  /*     offsetColor = PDM_MAX(offsetColor, newColor+1); */
  /*   } */

  /*   // for (int icell = 0; icell < nC; icell++){ */
  /*   //   offsetColor = PDM_MAX(offsetColor, newColor+1); */
  /*   // } */

  /*   free(ColorSub); */

  /* } */

  /* for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /*   int nC   = nCellByPoolByThread[iThr][iThrC]; */
  /*   int BegC = begCellByPoolByThread[iThr][iThrC]; */
  /*   int EndC = begCellByPoolByThread[iThr][iThrC]+nC; */

  /*   nCellVerif += nC; */

  /*   /\* Init *\/ */
  /*   for (int icell = 0; icell < nC+1; icell++){CellCellSubIdx[icell] = 0;} */

  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; idxCell++ ){ */

  /*       int ConnectedCell   = CellCell[idxCell]; */

  /*       if(ConnectedCell >= BegC && ConnectedCell < EndC ) */
  /*       { */
  /*         CellCellSubIdx[icell-BegC]++; */
  /*       } */
  /*     } */
  /*   } */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   _computeIdx(CellCellSubIdx, nC); */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */


  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   int ptArr = 0; */
  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; idxCell++ ){ */

  /*       int ConnectedCell = CellCell[idxCell]; */

  /*       if(ConnectedCell >= BegC && ConnectedCell < EndC ) */
  /*       { */
  /*         CellCellSub[ptArr] = CellCell[idxCell]-BegC; */
  /*         // printf(" [%i/%i] - CellCellSub[%i] = %i \n", iThr, iThrC, ptArr, CellCell[idxCell]-BegC); */
  /*         ptArr++; */
  /*       } */
  /*     } */
  /*   } */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   int nBlkCacheWantedSub; */
  /*   if(nCellPerCacheWanted == 0){nBlkCacheWantedSub = 1;} */
  /*   else                        {nBlkCacheWantedSub = PDM_MAX(nC/nCellPerCacheWanted,1);} */

  /*   // if(iThr != iThrC) */
  /*   // { */
  /*   //   nBlkCacheWantedSub = 1; */
  /*   // } */

  /*   /\* Split the graph *\/ */
  /*   int* ColorSub = NULL; */
  /*   printf("[%i/%i] PDM_part_graph_split_bis R2 : %i \n", iThr, iThrC, nC ); */
  /*   if(nC != 0) */
  /*   { */
  /*     PDM_part_graph_split_bis(          split_method, */
  /*                                        nBlkCacheWantedSub, */
  /*                                        nC, */
  /*                                        CellCellSubIdx, */
  /*                                        CellCellSub, */
  /*                              (int *)   NULL, */
  /*                              (int *)   NULL, */
  /*                                       &ColorSub); */
  /*   } */
  /*   // printf("[%i/%i] PDM_part_graph_split_bis end \n", iThr, iThrC); */
  /*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /*   if(0 == 1) */
  /*   { */
  /*     printf(" ----------- ColorSub : %i \n", nC); */
  /*     for (int i = 0; i < nC; i++){ */
  /*       printf(" =====> ColorSub[%i]  = %i\n", i, ColorSub[i]); */
  /*     } */
  /*   } */

  /*   int offsetColorTmp = offsetColor; */
  /*   for (int icell = 0; icell < nC; icell++){ */
  /*     int newColor =  offsetColorTmp+ColorSub[icell]; */
  /*     part->cellColor[BegC+icell] = newColor; */
  /*     // printf(" =====> cellColor[%i]  = %i\n", BegC+icell, newColor); */
  /*     offsetColor = PDM_MAX(offsetColor, newColor+1); */
  /*   } */

  /*   // for (int icell = 0; icell < nC; icell++){ */
  /*   //   offsetColor = PDM_MAX(offsetColor, newColor+1); */
  /*   // } */

  /*   free(ColorSub); */

  /* } */


  /* // for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /* //   int nC   = nCellByPoolByThreadR2[iThr][iThrC]; */
  /* //   int BegC = begCellByPoolByThreadR2[iThr][iThrC]; */
  /* //   int EndC = begCellByPoolByThreadR2[iThr][iThrC]+nC; */

  /* //   int offsetColorTmp = offsetColor; */
  /* //   for (int icell = 0; icell < nC; icell++){ */
  /* //     int newColor =  offsetColorTmp+ColorSub[icell]; */
  /* //     part->cellColor[BegC+icell] = part->cellColor[BegC+icell]+offsetColorTmp; */
  /* //     // printf(" =====> cellColor[%i]  = %i\n", BegC+icell, newColor); */
  /* //     offsetColor = PDM_MAX(offsetColor, newColor+1); */
  /* //   } */

  /* // } */




  /* } */
  /* /\* ----------------------------------------------------------- *\/ */
  /* printf("nCellVerif : %i \n",nCellVerif ); */
  /* printf("nCell : %i \n", nCell); */
  /* assert(nCellVerif == nCell); */

  /* if(0 == 1) */
  /* { */
  /*   printf(" ----------- nCell : %i \n", nCell); */
  /*   for (int i = 0; i < nCell; i++){ */
  /*     printf(" =====> part->cellColor[%i]  = %i\n", i, part->cellColor[i]); */
  /*     // assert(part->cellColor[i] != -1); */
  /*   } */
  /*   for (int i = 0; i < nCell; i++){ */
  /*     assert(part->cellColor[i] != -1); */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* /\* ----------------------------------------------------------- *\/ */
  /* int *maxColorByThread = (int *) malloc( nThread * sizeof(int) ); */
  /* int *minColorByThread = (int *) malloc( nThread * sizeof(int) ); */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   maxColorByThread[iThr] = -1; */
  /*   minColorByThread[iThr] = 900000000; */
  /* } */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */

  /*   // int maxColorByThread = -1; */
  /*   // int minColorByThread = 900000000; */

  /*   for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /*     int nC   = nCellByPoolByThread[iThr][iThrC]; */
  /*     int BegC = begCellByPoolByThread[iThr][iThrC]; */
  /*     int EndC = begCellByPoolByThread[iThr][iThrC]+nC; */

  /*     for (int icell = BegC; icell < EndC; icell++){ */
  /*       maxColorByThread[iThr] = PDM_MAX(maxColorByThread[iThr], part->cellColor[icell]); */
  /*       minColorByThread[iThr] = PDM_MIN(minColorByThread[iThr], part->cellColor[icell]); */
  /*     } */
  /*   } */

  /*   for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /*     int nC   = nCellByPoolByThreadR2[iThr][iThrC]; */
  /*     int BegC = begCellByPoolByThreadR2[iThr][iThrC]; */
  /*     int EndC = begCellByPoolByThreadR2[iThr][iThrC]+nC; */

  /*     for (int icell = BegC; icell < EndC; icell++){ */
  /*       maxColorByThread[iThr] = PDM_MAX(maxColorByThread[iThr], part->cellColor[icell]); */
  /*       minColorByThread[iThr] = PDM_MIN(minColorByThread[iThr], part->cellColor[icell]); */
  /*     } */
  /*   } */
  /*   printf("[%i] - maxColorByThread[iThr] = %i // minColorByThread[iThr] = %i \n", iThr, maxColorByThread[iThr], minColorByThread[iThr]); */


  /*   for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /*     int nC   = nCellByPoolByThread[iThr][iThrC]; */
  /*     int BegC = begCellByPoolByThread[iThr][iThrC]; */
  /*     int EndC = begCellByPoolByThread[iThr][iThrC]+nC; */


  /*     for (int icell = BegC; icell < EndC; icell++){ */
  /*       part->cellColor[icell] = maxColorByThread[iThr] + minColorByThread[iThr] - part->cellColor[icell]; */
  /*     } */
  /*   } */

  /*   for (int iThrC = 0; iThrC < nThread; iThrC++){ */
  /*     int nC   = nCellByPoolByThreadR2[iThr][iThrC]; */
  /*     int BegC = begCellByPoolByThreadR2[iThr][iThrC]; */
  /*     int EndC = begCellByPoolByThreadR2[iThr][iThrC]+nC; */


  /*     for (int icell = BegC; icell < EndC; icell++){ */
  /*       part->cellColor[icell] = maxColorByThread[iThr] + minColorByThread[iThr] - part->cellColor[icell]; */
  /*     } */
  /*   } */
  /* } */
  /* /\* ----------------------------------------------------------- *\/ */

  /* // printf("offsetColor : %i \n", offsetColor); */

  /* if(0 == 1) */
  /* { */
  /*   printf(" ----------- offsetColor : %i \n", offsetColor); */
  /*   for (int i = 0; i < nCell; i++){ */
  /*     printf(" =====> part->cellColor[%i]  = %i\n", i, part->cellColor[i]); */
  /*     assert(part->cellColor[i] != -1); */
  /*   } */
  /* } */


  /* int nNewBlkCacheWanted = 0; */

  /* for (int icell = 0; icell < nCell; icell++){ */
  /*   nNewBlkCacheWanted = PDM_MAX(nNewBlkCacheWanted, part->cellColor[icell]); */
  /* } */
  /* nNewBlkCacheWanted = nNewBlkCacheWanted + 1; */


  /* int* ColorByThread    = (int *) malloc( nNewBlkCacheWanted * sizeof(int)); */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   for (int iSub = minColorByThread[iThr]; iSub < maxColorByThread[iThr]+1; iSub++){ */
  /*     ColorByThread[iSub] = iThr; */
  /*     // printf("ColorByThread[%i] = %i \n", iSub, iThr ); */
  /*   } */
  /* } */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* int *nCellPerCache    = (int *) malloc( nNewBlkCacheWanted      * sizeof(int)); */
  /* int* cellTileIdx      = (int *) malloc((nNewBlkCacheWanted + 1) * sizeof(int)); */

  /* /\* Init to Zero *\/ */
  /* for (int i = 0; i < nNewBlkCacheWanted; i++){ */
  /*   nCellPerCache[i] = 0; */
  /* } */

  /* for (int iSub = 0; iSub < nNewBlkCacheWanted+1; iSub++){ */
  /*   cellTileIdx[iSub]=0; */
  /* } */

  /* /\* First loop to identify the size of each color *\/ */
  /* for (int i = 0; i < nCell; i++){ */
  /*   int color = part->cellColor[i]; */
  /*   nCellPerCache[color]++; */
  /*   CellOrder[i] = -1; */
  /* } */

  /* /\* Compute the index of each Sub-Domain **\/ */
  /* cellTileIdx[0] = 0; */
  /* for (int i = 0; i < nNewBlkCacheWanted; i++){ */
  /*   cellTileIdx[i + 1] = cellTileIdx[i] + nCellPerCache[i]; */
  /* } */

  /* /\* Reset array *\/ */
  /* for (int i = 0; i < nNewBlkCacheWanted; i++){ */
  /*   nCellPerCache[i] = 0; */
  /* } */

  /* for (int i = 0; i < nCell; i++){ */
  /*   int color  = part->cellColor[i]; */
  /*   int idx    = cellTileIdx[color] + nCellPerCache[color]; */

  /*   CellOrder[idx] = i; */
  /*   nCellPerCache[color]++; */
  /* } */

  /* PDM_part_reorder_cell(part, CellOrder); */

  /* free(nCellPerCache); */

  /* /\* ----------------------------------------------------------- *\/ */
  /* free(CellCellIdx); */
  /* free(CellCell); */
  /* CellCellIdx = NULL; */
  /* CellCell = NULL; */
  /* PDM_part_graph_compute_from_face_cell(          part, */
  /*                                       (int **) &CellCellIdx, */
  /*                                       (int **) &CellCell); */
  /* /\* ----------------------------------------------------------- *\/ */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* /\* */
  /*  * Prefer renumber completely interior first for each thread */
  /*  *\/ */
  /* int *ColorColorIdx           = (int *) malloc( (nNewBlkCacheWanted + 1) * sizeof(int)); */
  /* int *ColorColorRank1Idx      = (int *) malloc( (nNewBlkCacheWanted + 1) * sizeof(int)); */
  /* int *ColorColorCellRank2Idx  = (int *) malloc( (nNewBlkCacheWanted + 1) * sizeof(int)); */
  /* int *isAlreadyAssignedRank1  = (int *) malloc( (nNewBlkCacheWanted    ) * sizeof(int)); */
  /* int *isAlreadyAssignedRank2  = (int *) malloc( (nNewBlkCacheWanted    ) * sizeof(int)); */

  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*   isAlreadyAssignedRank1[iSub] = 0; */
  /*   isAlreadyAssignedRank2[iSub] = 0; */
  /* } */


  /* for (int iSub = 0; iSub < nNewBlkCacheWanted+1; iSub++){ */
  /*   ColorColorIdx         [iSub] = 0; */
  /*   ColorColorRank1Idx    [iSub] = 0; */
  /*   ColorColorCellRank2Idx[iSub] = 0; */
  /* } */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* /\* First Loop to count *\/ */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */

  /*   for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){ */
  /*     isAlreadyAssignedRank1 [iSub2] = -1; */
  /*     isAlreadyAssignedRank2 [iSub2] = -1; */
  /*   } */

  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = cellTileIdx[iSub  ]; */
  /*           icell < cellTileIdx[iSub+1]; */
  /*           icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; */
  /*             idxCell++ ){ */

  /*       int ConnectedCell  = CellCell[idxCell]; //-1; */
  /*       int ConnectedColor = part->cellColor[ConnectedCell]; */

  /*       /\* Check if already assigned for Rank 2 *\/ */
  /*       if(isAlreadyAssignedRank1[ConnectedColor] == -1) */
  /*       { */
  /*         isAlreadyAssignedRank1[ConnectedColor] = 1; */
  /*         ColorColorRank1Idx[iSub]++; */
  /*       } */

  /*       /\* Check if already assigned for Rank 2 *\/ */
  /*       if(isAlreadyAssignedRank2[ConnectedColor] == -1) */
  /*       { */
  /*         isAlreadyAssignedRank2[ConnectedColor] = 1; */
  /*         ColorColorCellRank2Idx[iSub]++; */
  /*       } /\* End if already Treated *\/ */

  /*       if(ConnectedColor != iSub) */
  /*       { */
  /*         /\* Check if a cell neighbour can be shared by multiple subdmain *\/ */
  /*         for(int idxCell2 = CellCellIdx[ConnectedCell  ]; */
  /*                 idxCell2 < CellCellIdx[ConnectedCell+1]; */
  /*                 idxCell2++ ) */
  /*         { */
  /*           int ConnectedCell2  = CellCell[idxCell2];//-1; */
  /*           int ConnectedColor2 = part->cellColor[ConnectedCell2]; */

  /*           // if(isAlreadyAssignedRank2[ConnectedColor2] == -1 && ConnectedColor2 < iSub) */
  /*           if(isAlreadyAssignedRank2[ConnectedColor2] == -1) */
  /*           { */
  /*             isAlreadyAssignedRank2[ConnectedColor2] = 1; */
  /*             ColorColorCellRank2Idx[iSub]++; */
  /*           } */
  /*         } /\* End for second rank *\/ */
  /*       } */
  /*     } */
  /*   } */
  /*   // printf("iSub : %i // R1 = %i // R2 = %i end  \n", iSub, ColorColorRank1Idx[iSub], ColorColorCellRank2Idx[iSub]); */
  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */


  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* _computeIdx(ColorColorRank1Idx    , nNewBlkCacheWanted); */
  /* _computeIdx(ColorColorCellRank2Idx, nNewBlkCacheWanted); */

  /* int *ColorColorRank1Arr     = (int *) malloc( ColorColorRank1Idx    [nNewBlkCacheWanted] * sizeof(int)); */
  /* int *ColorColorCellRank2Arr = (int *) malloc( ColorColorCellRank2Idx[nNewBlkCacheWanted] * sizeof(int)); */

  /* int ptRank1 = 0; */
  /* int ptRank2 = 0; */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* /\* Second Loop to fill *\/ */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */

  /*   for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){ */
  /*     isAlreadyAssignedRank1 [iSub2] = -1; */
  /*     isAlreadyAssignedRank2 [iSub2] = -1; */
  /*   } */

  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = cellTileIdx[iSub  ]; */
  /*           icell < cellTileIdx[iSub+1]; */
  /*           icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; */
  /*             idxCell++ ){ */

  /*       int ConnectedCell  = CellCell[idxCell];//-1; */
  /*       int ConnectedColor = part->cellColor[ConnectedCell]; */

  /*       /\* Check if already assigned for Rank 2 *\/ */
  /*       if(isAlreadyAssignedRank1[ConnectedColor] == -1) */
  /*       { */
  /*         isAlreadyAssignedRank1[ConnectedColor] = 1; */
  /*         ColorColorRank1Arr[ptRank1] = ConnectedColor; */
  /*         ptRank1++; */
  /*       } */

  /*       /\* Check if already assigned for Rank 2 *\/ */
  /*       if(isAlreadyAssignedRank2[ConnectedColor] == -1) */
  /*       { */
  /*         isAlreadyAssignedRank2[ConnectedColor] = 1; */
  /*         ColorColorCellRank2Arr[ptRank2] = ConnectedColor; */
  /*         ptRank2++; */
  /*       } /\* End if already Treated *\/ */

  /*       if(ConnectedColor != iSub) */
  /*       { */
  /*         /\* Check if a cell neighbour can be shared by multiple subdmain *\/ */
  /*         for(int idxCell2 = CellCellIdx[ConnectedCell  ]; */
  /*                 idxCell2 < CellCellIdx[ConnectedCell+1]; */
  /*                 idxCell2++ ) */
  /*         { */
  /*           int ConnectedCell2  = CellCell[idxCell2];//-1; */
  /*           int ConnectedColor2 = part->cellColor[ConnectedCell2]; */

  /*           // if(isAlreadyAssignedRank2[ConnectedColor2] == -1 && ConnectedColor2 < iSub) */
  /*           if(isAlreadyAssignedRank2[ConnectedColor2] == -1) */
  /*           { */
  /*             isAlreadyAssignedRank2[ConnectedColor2] = 1; */
  /*             ColorColorCellRank2Arr[ptRank2] = ConnectedColor2; */
  /*             ptRank2++; */
  /*           } */
  /*         } /\* End for second rank *\/ */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* int* isInterior    = (int *) malloc( nNewBlkCacheWanted * sizeof(int)); */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++) */
  /* { */
  /*   int colorLocThread = ColorByThread[iSub]; */
  /*   isInterior[iSub] = 1; */
  /*   for(int ptColor = ColorColorCellRank2Idx[iSub]; ptColor < ColorColorCellRank2Idx[iSub+1]; ptColor++) */
  /*   { */
  /*     if(ColorByThread[ColorColorCellRank2Arr[ptColor]] != colorLocThread) */
  /*     { */
  /*       isInterior[iSub] = -1; */
  /*     } */
  /*   } */
  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* /\* */
  /*  * Prepare */
  /*  *\/ */
  /* int* nSubIntByThread  = (int *) malloc(   nThread     * sizeof(int)); */
  /* int* nSubExtByThread  = (int *) malloc(   nThread     * sizeof(int)); */
  /* int* nSubByThread     = (int *) malloc(   nThread     * sizeof(int)); */
  /* int* ThreadTileSubIdx = (int *) malloc( ( nThread+1 ) * sizeof(int)); */
  /* int nInt = 0; */
  /* int nExt = 0; */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   // nSubPerThread[iThr]   = 0; */
  /*   nSubIntByThread[iThr] = 0; */
  /*   nSubExtByThread[iThr] = 0; */
  /* } */

  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++) */
  /* { */
  /*   int color = ColorByThread[iSub]; */
  /*   if(isInterior[iSub] == 1){nInt++;nSubIntByThread[color]++;} */
  /*   else                     {nExt++;nSubExtByThread[color]++;} */
  /* } */

  /* printf("nInt = [%i] // nExt = %i \n", nInt, nExt); */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   printf("Thread[%i] : nInt = [%i] // nExt = %i \n", iThr, nSubIntByThread[iThr], nSubExtByThread[iThr]); */
  /* } */

  /* ThreadTileSubIdx[nThread] = 0; */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   nSubByThread[iThr]     = nSubIntByThread[iThr] + nSubExtByThread[iThr]; */
  /*   ThreadTileSubIdx[iThr] = nSubByThread[iThr]; */
  /* } */
  /* _computeIdx(ThreadTileSubIdx, nThread); */







  /* // for (int iSub = 0; iSub < nNewBlkCacheWanted+1; iSub++){ */
  /* //   printf("ColorColorRank1ExtIdx[%i] = %i \n", iSub, ColorColorRank1ExtIdx[iSub]); */
  /* //   // printf("ColorColorRank1Idx[%i] = %i \n", iSub, ColorColorRank1Idx[iSub]); */
  /* //   // printf("ColorColorRank2Idx[%i] = %i \n", iSub, ColorColorRank2Idx[iSub]); */
  /* // } */



  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */


  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* int* OldToNewColor = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) ); */
  /* int* isTreatedSdom = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) ); */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */

  /*   int nC   = nCellByPoolByThread  [iThr][iThr]; */
  /*   int BegC = begCellByPoolByThread[iThr][iThr]; */
  /*   int EndC = begCellByPoolByThread[iThr][iThr]+nC; */

  /*   int minSub = 900000000; */
  /*   int maxSub = -1; */

  /*   int firstCouplingColor = 0; */
  /*   int lastCouplingColor  = 0; */

  /*   for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*     isTreatedSdom[iSub] = -1; */
  /*     OldToNewColor[iSub] = -1; */
  /*   } */

  /*   if(0 == 1) */
  /*   { */
  /*     for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*       printf(" [%i] link to ", iSub); */
  /*       for(int ptColor = ColorColorCellRank2Idx[iSub]; */
  /*               ptColor < ColorColorCellRank2Idx[iSub+1]; ptColor++) */
  /*       { */
  /*         printf("%i ", ColorColorCellRank2Arr[ptColor]); */
  /*       } */
  /*       printf("\n"); */
  /*     } */
  /*   } */

  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     int iSubColor = part->cellColor[icell]; */
  /*     maxSub = PDM_MAX(maxSub, iSubColor); */
  /*     minSub = PDM_MIN(minSub, iSubColor); */

  /*     if(isTreatedSdom[iSubColor] == -1) */
  /*     { */
  /*       if  (isInterior[iSubColor]  == -1) */
  /*       { */
  /*         for(int ptColor = ColorColorCellRank2Idx[iSubColor]; */
  /*                 ptColor < ColorColorCellRank2Idx[iSubColor+1]; ptColor++) */
  /*         { */
  /*           if(ColorByThread[ColorColorCellRank2Arr[ptColor]] > iThr && isTreatedSdom[iSubColor] == -1) */
  /*           { */
  /*             firstCouplingColor++; */
  /*             isTreatedSdom[iSubColor] = 1; */
  /*           } */
  /*           else if(ColorByThread[ColorColorCellRank2Arr[ptColor]] < iThr  && isTreatedSdom[iSubColor] == -1) */
  /*           { */
  /*             lastCouplingColor++; */
  /*             isTreatedSdom[iSubColor] = 1; */
  /*           } */
  /*         } */
  /*       } */
  /*     } */
  /*   } */


  /*   printf("[%i] : firstCouplingColor : %i \n", iThr, firstCouplingColor); */
  /*   printf("[%i] : lastCouplingColor      : %i \n", iThr, lastCouplingColor); */
  /*   printf("[%i] : minSub  : %i \n", iThr, minSub); */
  /*   printf("[%i] : maxSub  : %i \n", iThr, maxSub); */

  /*   /\* Deduced *\/ */
  /*   int begFirstCouplingColor = minSub; */
  /*   int begInteriorColor      = minSub+firstCouplingColor; */
  /*   int begLastCouplingColor  = maxSub-lastCouplingColor+1; */

  /*   int nFirstCouplingColor = 0; */
  /*   int nInteriorColor      = 0; */
  /*   int nLastCouplingColor  = 0; */

  /*   printf("[%i] : begFirstCouplingColor : %i \n", iThr, begFirstCouplingColor); */
  /*   printf("[%i] : begInteriorColor      : %i \n", iThr, begInteriorColor); */
  /*   printf("[%i] : begLastCouplingColor  : %i \n", iThr, begLastCouplingColor); */

  /*   for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*     isTreatedSdom[iSub] = -1; */
  /*   } */


  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     int iSubColor = part->cellColor[icell]; */

  /*     if(isTreatedSdom[iSubColor] == -1) */
  /*     { */
  /*       if  (isInterior[iSubColor]  == 1){ */
  /*         OldToNewColor[iSubColor] = begInteriorColor+nInteriorColor; */
  /*         // OldToNewColor[begInteriorColor+nInteriorColor] = iSubColor; */
  /*         nInteriorColor++; */
  /*         isTreatedSdom[iSubColor] = 1; */
  /*       } */
  /*       else */
  /*       { */
  /*         for(int ptColor = ColorColorCellRank2Idx[iSubColor]; */
  /*                 ptColor < ColorColorCellRank2Idx[iSubColor+1]; ptColor++) */
  /*         { */
  /*           if(ColorByThread[ColorColorCellRank2Arr[ptColor]] > iThr&& isTreatedSdom[iSubColor] == -1) */
  /*           { */
  /*             OldToNewColor[iSubColor] = begFirstCouplingColor+nFirstCouplingColor; */
  /*             // OldToNewColor[begFirstCouplingColor+nFirstCouplingColor] = iSubColor; */
  /*             nFirstCouplingColor++; */
  /*             isTreatedSdom[iSubColor] = 1; */
  /*           } */
  /*           else if(ColorByThread[ColorColorCellRank2Arr[ptColor]] < iThr && isTreatedSdom[iSubColor] == -1) */
  /*           { */
  /*             OldToNewColor[iSubColor] = begLastCouplingColor+nLastCouplingColor; */
  /*             // OldToNewColor[begLastCouplingColor+nLastCouplingColor] = iSubColor; */
  /*             nLastCouplingColor++; */
  /*             isTreatedSdom[iSubColor] = 1; */
  /*           } */
  /*         } */
  /*       } */
  /*     } */
  /*   } */

  /*   int isInteriorTmp[nNewBlkCacheWanted]; */
  /*   for (int iSub = minSub; iSub < maxSub+1; iSub++){ */
  /*     printf("OldToNewColor[%i] = %i // isInterior[%i] = %i \n", iSub, OldToNewColor[iSub], iSub, isInterior[iSub]); */
  /*     // isInterior[OldToNewColor[iSub]] = isInterior[iSub]; */
  /*     isInteriorTmp[OldToNewColor[iSub]] = isInterior[iSub]; */
  /*     // isInteriorTmp[iSub] = isInterior[OldToNewColor[iSub]]; */
  /*   } */

  /*   for (int iSub = minSub; iSub < maxSub+1; iSub++){ */
  /*     printf(" Check isInterior[%i] = %i  \n", iSub, isInteriorTmp[iSub]); */
  /*     isInterior[iSub] = isInteriorTmp[iSub]; */
  /*   } */

  /*   for (int icell = BegC; icell < EndC; icell++){ */
  /*     part->cellColor[icell] = OldToNewColor[part->cellColor[icell]]; */
  /*   } */

  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */


  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* /\* First Loop to count *\/ */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted+1; iSub++){ */
  /*   ColorColorCellRank2Idx[iSub] = 0; */
  /* } */

  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */

  /*   for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){ */
  /*     isAlreadyAssignedRank2 [iSub2] = -1; */
  /*   } */

  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = cellTileIdx[iSub  ]; */
  /*           icell < cellTileIdx[iSub+1]; */
  /*           icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; */
  /*             idxCell++ ){ */

  /*       int ConnectedCell  = CellCell[idxCell]; //-1; */
  /*       int ConnectedColor = part->cellColor[ConnectedCell]; */

  /*       /\* Check if already assigned for Rank 2 *\/ */
  /*       if(isAlreadyAssignedRank2[ConnectedColor] == -1) */
  /*       { */
  /*         isAlreadyAssignedRank2[ConnectedColor] = 1; */
  /*         ColorColorCellRank2Idx[iSub]++; */
  /*       } /\* End if already Treated *\/ */

  /*       if(ConnectedColor != iSub) */
  /*       { */
  /*         /\* Check if a cell neighbour can be shared by multiple subdmain *\/ */
  /*         for(int idxCell2 = CellCellIdx[ConnectedCell  ]; */
  /*                 idxCell2 < CellCellIdx[ConnectedCell+1]; */
  /*                 idxCell2++ ) */
  /*         { */
  /*           int ConnectedCell2  = CellCell[idxCell2];//-1; */
  /*           int ConnectedColor2 = part->cellColor[ConnectedCell2]; */

  /*           // if(isAlreadyAssignedRank2[ConnectedColor2] == -1 && ConnectedColor2 < iSub) */
  /*           if(isAlreadyAssignedRank2[ConnectedColor2] == -1) */
  /*           { */
  /*             isAlreadyAssignedRank2[ConnectedColor2] = 1; */
  /*             ColorColorCellRank2Idx[iSub]++; */
  /*           } */
  /*         } /\* End for second rank *\/ */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* _computeIdx(ColorColorCellRank2Idx, nNewBlkCacheWanted); */
  /* ptRank2 = 0; */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* /\* Second Loop to fill *\/ */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */

  /*   for (int iSub2 = 0; iSub2 < nNewBlkCacheWanted; iSub2++){ */
  /*     isAlreadyAssignedRank2 [iSub2] = -1; */
  /*   } */

  /*   /\* Loop on cell of current SubDomain *\/ */
  /*   for(int icell = cellTileIdx[iSub  ]; */
  /*           icell < cellTileIdx[iSub+1]; */
  /*           icell++){ */
  /*     for(int idxCell = CellCellIdx[icell  ]; */
  /*             idxCell < CellCellIdx[icell+1]; */
  /*             idxCell++ ){ */

  /*       int ConnectedCell  = CellCell[idxCell];//-1; */
  /*       int ConnectedColor = part->cellColor[ConnectedCell]; */

  /*       /\* Check if already assigned for Rank 2 *\/ */
  /*       if(isAlreadyAssignedRank2[ConnectedColor] == -1) */
  /*       { */
  /*         isAlreadyAssignedRank2[ConnectedColor] = 1; */
  /*         ColorColorCellRank2Arr[ptRank2] = ConnectedColor; */
  /*         ptRank2++; */
  /*       } /\* End if already Treated *\/ */

  /*       if(ConnectedColor != iSub) */
  /*       { */
  /*         /\* Check if a cell neighbour can be shared by multiple subdmain *\/ */
  /*         for(int idxCell2 = CellCellIdx[ConnectedCell  ]; */
  /*                 idxCell2 < CellCellIdx[ConnectedCell+1]; */
  /*                 idxCell2++ ) */
  /*         { */
  /*           int ConnectedCell2  = CellCell[idxCell2];//-1; */
  /*           int ConnectedColor2 = part->cellColor[ConnectedCell2]; */

  /*           // if(isAlreadyAssignedRank2[ConnectedColor2] == -1 && ConnectedColor2 < iSub) */
  /*           if(isAlreadyAssignedRank2[ConnectedColor2] == -1) */
  /*           { */
  /*             isAlreadyAssignedRank2[ConnectedColor2] = 1; */
  /*             ColorColorCellRank2Arr[ptRank2] = ConnectedColor2; */
  /*             ptRank2++; */
  /*           } */
  /*         } /\* End for second rank *\/ */
  /*       } */
  /*     } */
  /*   } */
  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* printf(" ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   for (int iSubI = ThreadTileSubIdx[iThr]; iSubI < ThreadTileSubIdx[iThr+1]; iSubI++){ */
  /*     printf("[%i] ---> isInterior[%i] = %i \n", iThr, iSubI, isInterior[iSubI]); */
  /*   } */
  /* } */
  /* printf(" ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */

  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++) */
  /* { */
  /*   int colorLocThread = ColorByThread[iSub]; */
  /*   isInterior[iSub] = 1; */
  /*   for(int ptColor = ColorColorCellRank2Idx[iSub]; ptColor < ColorColorCellRank2Idx[iSub+1]; ptColor++) */
  /*   { */
  /*     if(ColorByThread[ColorColorCellRank2Arr[ptColor]] != colorLocThread) */
  /*     { */
  /*       isInterior[iSub] = -1; */
  /*     } */
  /*   } */
  /* } */
  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* printf(" ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   for (int iSubI = ThreadTileSubIdx[iThr]; iSubI < ThreadTileSubIdx[iThr+1]; iSubI++){ */
  /*     printf("[%i] ---> isInterior[%i] = %i \n", iThr, iSubI, isInterior[iSubI]); */
  /*   } */
  /* } */
  /* printf(" ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */



  /* /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* int* currentSub = (int * ) malloc( nThread * sizeof(int) ); */
  /* int* waitThr    = (int * ) malloc( nThread * sizeof(int) ); */
  /* int* waitThrMax = (int * ) malloc( nThread * sizeof(int) ); */

  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*   isTreatedSdom[iSub] = -1; */
  /*   OldToNewColor [iSub] = iSub; */
  /* } */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   waitThr[iThr]    = -1; */
  /*   waitThrMax[iThr] = -1; */
  /*   currentSub[iThr] = ThreadTileSubIdx[iThr]; */
  /* } */

  /* for (int iThr = 0; iThr < nThread+1; iThr++){ */
  /*   printf("ThreadTileSubIdx[%i] = %i \n", iThr, ThreadTileSubIdx[iThr]); */
  /* } */


  /* int  nSdomTreated = 0; */
  /* while(nSdomTreated != nNewBlkCacheWanted) */
  /* { */
  /*   // printf("nSdomTreated : %i // %i \n", nSdomTreated, nNewBlkCacheWanted); */
  /*   for (int iThr = 0; iThr < nThread; iThr++){ */

  /*     int iSub           = currentSub[iThr]; */
  /*     int isOkToContinue = 1; */

  /*     // printf("[%i] : iSub = %i \n", iThr, iSub); */

  /*     /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*     // if(iSub > -1 && isTreatedSdom[iSub] == -1) */
  /*     if(iSub > -1 ) */
  /*     { */
  /*       /\* -------------------------------------------------------------- *\/ */

  /*       // printf("\t\t Link to "); */
  /*       for(int ptColor = ColorColorCellRank2Idx[iSub  ]; */
  /*               ptColor < ColorColorCellRank2Idx[iSub+1]; ptColor++) */
  /*       { */
  /*         int iSubLink = ColorColorCellRank2Arr[ptColor]; */

  /*         // printf("%i ", iSubLink); */
  /*         if(iSubLink < ThreadTileSubIdx[iThr]) */
  /*         { */
  /*           if(isTreatedSdom[iSubLink] == -1) */
  /*           { */
  /*             isOkToContinue = -1; */
  /*             break; */
  /*           } */
  /*         } */
  /*       } */

  /*       // printf("\n "); */
  /*       /\* -------------------------------------------------------------- *\/ */

  /*       // printf("\t\t isOkToContinue : %i\n", isOkToContinue); */

  /*       /\* -------------------------------------------------------------- *\/ */
  /*       if(isOkToContinue == 1) */
  /*       { */
  /*         isTreatedSdom[iSub] = 1; */
  /*         currentSub[iThr] = currentSub[iThr] + 1; */
  /*         nSdomTreated++; */

  /*         waitThrMax[iThr] = PDM_MAX(waitThr[iThr], waitThrMax[iThr]); */
  /*         waitThr[iThr]    = 0; */

  /*         if(currentSub[iThr] == ThreadTileSubIdx[iThr+1]) */
  /*         { */
  /*           currentSub[iThr] = -100000; */
  /*         } */
  /*       } */
  /*       else */
  /*       { */
  /*         waitThr[iThr]++; */
  /*       } */
  /*       /\* -------------------------------------------------------------- *\/ */
  /*     } */
  /*     /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */
  /*   } */
  /* } */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   printf("waitThrMax[%i] = %i \n", iThr, waitThrMax[iThr]); */
  /* } */



  /* /\* */
  /*  * */
  /*  *\/ */
  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*   isTreatedSdom[iSub] = -1; */
  /*   OldToNewColor[iSub] = iSub; */
  /* } */

  /* for (int iThr = 0; iThr < nThread; iThr++){ */
  /*   waitThr[iThr]    = -1; */
  /*   waitThrMax[iThr] = -1; */
  /*   currentSub[iThr] = ThreadTileSubIdx[iThr]; */
  /* } */

  /* nSdomTreated = 0; */

  /* /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */

  /* int* NewToOldColor = (int * ) malloc( nNewBlkCacheWanted * sizeof(int) ); */

  /* for (int iSub = 0; iSub < nNewBlkCacheWanted; iSub++){ */
  /*   isTreatedSdom[iSub] = -1; */
  /*   NewToOldColor[iSub] = iSub; */
  /*   OldToNewColor[iSub] = iSub; */
  /* } */


/*   for (int iThr = 0; iThr < nThread; iThr++){ */

/*     for (int iSubI = ThreadTileSubIdx[iThr]; iSubI < ThreadTileSubIdx[iThr+1]; iSubI++){ */


/*       // int iSub = OldToNewColor[iSubI]; */
/*       int iSub = NewToOldColor[iSubI]; */
/*       int rSub = iSubI - ThreadTileSubIdx[iThr]; */

/*       printf("[%i] ---------------------------------------------- : [%i/%i] \n", iThr, iSubI, iSub); */
/*       // printf(" - iSub : %i \n", ); */
/*       printf("[%i] - rSub : %i \n", iThr, rSub); */

/*       /\* ----------------------------------------------------------------------------------- *\/ */
/*       int rSubLink = -1; */
/*       for(int ptColor = ColorColorCellRank2Idx[iSub  ]; */
/*               ptColor < ColorColorCellRank2Idx[iSub+1]; ptColor++) */
/*       { */
/*         // int iSubLink = ColorColorCellRank2Arr[ptColor]; */
/*         int iSubLink = OldToNewColor[ColorColorCellRank2Arr[ptColor]]; */
/*         int iThrLink = ColorByThread[iSubLink]; */
/*         if(iThrLink < iThr) */
/*         { */
/*           // int rSubLink = (iSubLink - ThreadTileSubIdx[iThrLink]) - rSub; */
/*           rSubLink = PDM_MAX( rSubLink, (iSubLink - ThreadTileSubIdx[iThrLink]) - rSub); */
/*           printf("[%i] - rSubLink : %i // iSubLink : %i \n", iThr, rSubLink, iSubLink); */
/*         } */

/*       } */
/*       /\* ----------------------------------------------------------------------------------- *\/ */

/*       /\* ----------------------------------------------------------------------------------- *\/ */
/*       if(rSubLink > 0) */
/*       { */
/*         printf("[%i] - \t Maybe a blocking %i \n", iThr, rSubLink); */

/*         int iSubP = -1; */

/*         /\* On recherche d'abord parmi les bloquant : iThrLink > iThr *\/ */
/*         int rSubPLinkSave = nNewBlkCacheWanted+1; */
/*         for (int iSubO = iSubI+1; iSubO < ThreadTileSubIdx[iThr+1]; iSubO++){ */

/*           printf(" ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n"); */
/*           int iSubCandidate = NewToOldColor[iSubO]; */
/*           int isBlockant = -1; */
/*           // int isBlocker  = -1; */
/*           for(int ptColor = ColorColorCellRank2Idx[iSubCandidate  ]; */
/*                   ptColor < ColorColorCellRank2Idx[iSubCandidate+1]; ptColor++) */
/*           { */
/*             int iSubLink = OldToNewColor[ColorColorCellRank2Arr[ptColor]]; */
/*             int iThrLink = ColorByThread[iSubLink]; */
/*             if(iThrLink > iThr) */
/*             { */
/*               // printf("Find a bloquant : %i \n", iSubCandidate); */
/*               // iSubP = iSubCandidate; */
/*               isBlockant = 1; */
/*               break; */
/*             } */
/*           } */


/*           int rSubPLink = -1; */
/*           for(int ptColor = ColorColorCellRank2Idx[iSubCandidate  ]; */
/*                   ptColor < ColorColorCellRank2Idx[iSubCandidate+1]; ptColor++) */
/*           { */
/*             int iSubLink = OldToNewColor[ColorColorCellRank2Arr[ptColor]]; */
/*             int iThrLink = ColorByThread[iSubLink]; */
/*             if(iThrLink < iThr) */
/*             { */
/*               // printf("Find a blocker : %i \n", iSubCandidate); */
/*               // iSubP = iSubCandidate; */
/*               rSubPLink = PDM_MAX( rSubPLink, (iSubLink - ThreadTileSubIdx[iThrLink]) - rSub); */
/*               // printf("\t [%i] - rSubPLink : %i // iSubLink : %i \n", iThr, rSubPLink, iSubLink); */
/*               // isBlocker = 1; */
/*               // break; */
/*             } */
/*           } */

/*           if( rSubPLink < 0 && isBlockant == 1) */
/*           { */
/*             iSubP = iSubCandidate; */
/*             break; */
/*           } */
/*           else if( rSubPLink < 0) */
/*           { */
/*             iSubP = iSubCandidate; */
/*             break; */
/*           } */
/*           else if( rSubPLink < rSubPLinkSave) */
/*           { */
/*             rSubPLinkSave = rSubPLink; */
/*             iSubP = iSubCandidate; */
/*           } */
/*           // if(iSubP != -1){break;} */
/*         } */

/*         printf(" \t \t Found the funcking best permutation with %i  \n", iSubP); */
/*         // assert(iSubP > iSub); */

/*         if(iSubP != -1){ */

/*             printf(" \t \t Found a permutation with %i  \n", iSubP); */

/*             // OldToNewColor[iSubP] = iSub; */
/*             OldToNewColor[iSubP] = iSubI; */
/*             printf("\t Assign OldToNewColor : OldToNewColor[%i] = %i \n", iSubP, OldToNewColor[iSubP]); */

/*             int dSub = 0; */
/*             for (int iSubA = iSub; iSubA < iSubP; iSubA++){ */

/*               // OldToNewColor[iSubA] = OldToNewColor[iSubA] + 1; */

/*               if(isTreatedSdom[iSubA] == -1) */
/*               { */
/*                 // OldToNewColor[iSubA] = iSubI + iSubA - iSub + 1; */
/*                 OldToNewColor[iSubA] = iSubI + dSub + 1; */
/*                 dSub++; */
/*               } */
/*               printf("\t Assign OldToNewColor : OldToNewColor[%i] = %i \n", iSubA, OldToNewColor[iSubA]); */
/*             } */

/*             for (int iSubA = ThreadTileSubIdx[iThr]; iSubA < ThreadTileSubIdx[iThr+1]; iSubA++){ */
/*               NewToOldColor[OldToNewColor[iSubA]] = iSubA; */
/*             } */

/*             isTreatedSdom[iSubP] = 1; */
/*             // isTreatedSdom[iSubI] = 1; */

/*             if(1 == 1) */
/*             { */
/*               printf("OldToNewColor Apres : [%i/%i]\n", ThreadTileSubIdx[iThr], ThreadTileSubIdx[iThr+1]); */
/*               for (int iSubA = ThreadTileSubIdx[iThr]; iSubA < ThreadTileSubIdx[iThr+1]; iSubA++){ */
/*                 printf(" %i ", OldToNewColor[iSubA]); */
/*               } */
/*               printf("\n"); */

/*               printf("NewToOldColor Apres : [%i/%i]\n", ThreadTileSubIdx[iThr], ThreadTileSubIdx[iThr+1]); */
/*               for (int iSubA = ThreadTileSubIdx[iThr]; iSubA < ThreadTileSubIdx[iThr+1]; iSubA++){ */
/*                 printf(" %i ", NewToOldColor[iSubA]); */
/*               } */
/*               printf("\n"); */
/*             } */
/*             // break; */

/*           // } */
/*         } */


/*       } */
/*       /\* ----------------------------------------------------------------------------------- *\/ */


/*     } */



/*   } */


/*   for (int icell = 0; icell < nCell; icell++){ */
/*     part->cellColor[icell] = OldToNewColor[part->cellColor[icell]]; */
/*   } */


/*   /\* oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */





/*   /\* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo *\/ */


/*   /\* */
/*    * Free */
/*    *\/ */
/*   free(nCellPerThread); */
/*   free(nCellPerThreadBeg); */
/*   free(threadCellIdx); */
/*   free(CellCellSubIdx); */
/*   free(CellCellSub); */
/*   free(CellCellIdx); */
/*   free(CellCell); */
/*   free(cellTileIdx         ); */
/*   free(ColorColorIdx         ); */
/*   free(ColorColorRank1Idx    ); */
/*   free(ColorColorCellRank2Idx    ); */
/*   free(ColorColorRank1Arr    ); */
/*   free(ColorColorCellRank2Arr    ); */
/*   free(isAlreadyAssignedRank1); */
/*   free(isAlreadyAssignedRank2); */
/*   free(isInterior); */
/*   free(ColorByThread); */
/*   free(maxColorByThread); */
/*   free(minColorByThread); */
/*   free(nSubIntByThread); */
/*   free(nSubExtByThread); */
/*   free(OldToNewColor); */
/*   free(isTreatedSdom); */
/*   free(currentSub); */
/*   free(waitThr); */
/*   free(waitThrMax); */
/*   free(nSubByThread    ); */
/*   free(ThreadTileSubIdx); */
/*   free(NewToOldColor); */





/*   return nNewBlkCacheWanted; */
/* } */


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Add renumbering cacheblocking method in the ppart methods
 *
 */

void
PROCF(pdm_renum_cacheblocking_ppart_add, PDM_RENUM_CACHEBLOCKING_PPART_ADD)
(
void
 )
{
  PDM_renum_cacheblocking_ppart_add() ;
}

void
PDM_renum_cacheblocking_ppart_add
(
void
)
{
  PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_CACHEBLOCKING", _renum_cells_cacheblocking);
}


/**
 * \brief Add renumbering cacheblocking method in the ppart methods
 *
 */

void
PROCF(pdm_renum_cacheblocking2_ppart_add, PDM_RENUM_CACHEBLOCKING2_PPART_ADD)
(
void
 )
{
  PDM_renum_cacheblocking2_ppart_add() ;
}


void
PDM_renum_cacheblocking2_ppart_add
(
void
)
{
PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_CACHEBLOCKING2", _renum_cells_cacheblocking2);
}




/**
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition
 * \param [in]   split_method        Split method
 * \param [in]   nCellPerCacheWanted Approximate number of cells on each cache
 * \param [in]   isAsynchrone        [0 : Synchrone/1 : Asynchrone]
 * \param [in]   isVectorisation     [0 : No/1 : Yes]
 *
 */

void
PDM_renum_cacheblocking
(
 _part_t     *part,
int           split_method,
int           nCellPerCacheWanted,
int           isAsynchrone,
int           isVectorisation,
int           nVectFace
)
{
  /* Get nFac and nCel */
  const int nCell = part->nCell;
  const int nFace = part->nFace;

  /** Allocate reoerdering/permutation array **/
  int *CellOrder = (int *) malloc (sizeof(int) * nCell);
  int *FaceOrder = (int *) malloc (sizeof(int) * nFace);

  /*
   * I/ Coloring the cells with classic graph library
   */
  int *CellCellIdx = NULL;
  int *CellCell    = NULL;

  assert(part->faceColor == NULL);
  assert(part->cellColor == NULL);

  /* Compute graph associate to mesh */
  PDM_part_graph_compute_from_face_cell(          part,
                                        (int **) &CellCellIdx,
                                        (int **) &CellCell);

  /*
   * Determine the optimal size for cache blocking
   *   -> The user specify the desired number of cells on each subdomain
   *   -> An other idea is to specify directly the number of block he want
   */
  int nBlkCacheWanted;
  if(nCellPerCacheWanted == 0){nBlkCacheWanted = 1;}
  else                        {nBlkCacheWanted = PDM_MAX(nCell/nCellPerCacheWanted,1);}

  /* Split the graph */
  PDM_part_graph_split(          split_method,
                                 nBlkCacheWanted,
                                 part,
                                 CellCellIdx,
                                 CellCell,
                       (int *)   NULL,
                       (int *)   NULL,
                       (int **) &part->cellColor);

  if(0 == 1)
  {
    printf(" ----------- nBlkCacheWanted : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> part->cellColor[%i]  = %i\n", i, part->cellColor[i]);
    }
  }
  /*
   * II/ Create a proper cells order :
   *       -> [ [SubDom1] , [SubDom2], ..., [SubDomN] ]
   */

  /* Allocate */
  int *nCellPerCache    = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int *nCellPerCacheBeg = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int *partCellIdx      = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));

  /* II-Bis - Asynchrone : We change only number of each sub-domain
   * in order to have interior sub-domain first and exterior after
   */
  if(isAsynchrone == 1)
  {
    /* Allocate */
    int *flagHaveBnd     = (int *) malloc(nBlkCacheWanted * sizeof(int));
    int *syncToAsynColor = (int *) malloc(nBlkCacheWanted * sizeof(int));

    /* All sub-domain is initialize */
    for(int i = 0; i < nBlkCacheWanted; i++) {
      flagHaveBnd[i] = 0;
    }

    /*
     * Loop on faces : if at least one is border, flagHaveBnd is set to 1
     */
    for (int i = 0; i < part->nFace; i++) {
      int iCell1 = part->faceCell[2*i    ];
      int iCell2 = part->faceCell[2*i + 1];
      int color1 = part->cellColor[iCell1-1];

      /* Solution 1 */
      // if(iCell2 > 0 ){
      //   int color2 = part->cellColor[iCell2-1];
      //   int color  = _PDM_part_MIN(color1, color2);
      // }
      // else
      // {
      //   flagHaveBnd[color1] = 1;
      // }

      /* Solution 2 */
      if(iCell2 == 0){flagHaveBnd[color1] = 1;}

    }

    /*
     * Built the syncToAsynColor array
     */
    int nSdomB = 0;                     /* Number of   Blocking domain */
    int nSdomU = 0;                     /* Number of Unblocking domain */
    int nextColorB = nBlkCacheWanted-1; /* Next Blocking   color       */
    int nextColorU = 0;                 /* Next UnBlocking color       */

    for(int i = 0; i < nBlkCacheWanted; i++) {
      if(flagHaveBnd[i] == 0)
      {
        syncToAsynColor[i] = nextColorU;
        nextColorU++;
        nSdomU++;
      }
      else
      {
        syncToAsynColor[i] = nextColorB;
        nextColorB--;
        nSdomB++;
      }
    }

    /* Dramatic verbose */
    if(0 == 1)
    {
      printf(" ----------- : %i \n", nBlkCacheWanted);
      for (int i = 0; i < nBlkCacheWanted; i++){
        printf("SyncToAsynColor[%i] = %i\n", i, syncToAsynColor[i]);
      }
    }

    /*
     * Apply the computed array on cellPart
     */
    for(int i = 0; i < part->nCell; i++) {
      part->cellColor[i] = syncToAsynColor[part->cellColor[i]];
    }


    /* Free memory */
    free(flagHaveBnd);
    free(syncToAsynColor);

  } /* End asynchrone */


  /* Init to Zero */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nCellPerCache[i] = 0;
  }

  /* First loop to identify the size of each color */
  for (int i = 0; i < nCell; i++){
    int color = part->cellColor[i]; // A color is a number of partition (output of Metis or Scotch)
    nCellPerCache[color]++;
  }

  /* Compute the index of each Sub-Domain **/
  partCellIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partCellIdx[i + 1] = partCellIdx[i] + nCellPerCache[i];
  }


  if(isVectorisation == 0){
    /** First method - Standard **/

    /* Reset array */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCache[i] = 0;
    }

    for (int i = 0; i < nCell; i++){
      int color  = part->cellColor[i];
      int idx    = partCellIdx[color] + nCellPerCache[color];

      CellOrder[idx] = i;
      nCellPerCache[color]++;
    }
  }
  else
  {
    /** Second method - Vectorisation on cell on current domain **/

    /* Store the current Index to add a vectorisable cell */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCacheBeg[i] = 0;
    }

    /* Loop on cells :
     *     -> If current cell is adjacent to another color : add to end
     *     -> Else add to begin
     * In fact this function sort the interior cell and exterior cell leads to the following numbering :
     *         [ [ Interior Cells] , [ Exterior Cells]]
     */
    int verifCell = 0;
    for (int i = 0; i < nCell; i++){

      /* Get color and prepare flag */
      int color = part->cellColor[i];
      int flag  = -1;
      for(int j = CellCellIdx[i]; j < CellCellIdx[i+1]; j++){
        int iCell = CellCell[j];
        // if(part->cellColor[iCell] != color){flag = 1;} // Alors on est au bord d'un sous domaine !
        if(part->cellColor[iCell] < color){flag = 1;} // Alors on est au bord d'un sous domaine !
      }

      if(flag == -1){  // Cell is interior : add to begin

        int idx = partCellIdx[color] + nCellPerCacheBeg[color];

        CellOrder[idx] = i;
        nCellPerCacheBeg[color]++;
        verifCell++;

      }
      else{ // Cell is exterior : add to end

        int idx = partCellIdx[color] + nCellPerCache[color]-1;

        CellOrder[idx] = i;
        nCellPerCache[color]--;
        verifCell++;

      }

      /* Panic verbose */
      // if(0 == 1){
      //   PDM_printf("Begin : %i %i %i \n", color, idx, nCellPerCacheBeg[color]);
      // }

    }

    assert(verifCell == part->nCell);

  }

  /*
   * Verbose
   */
  if(0 == 1)
  {
    printf(" ----------- nBlkCacheWanted : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> CellOrder[%i]  = %i\n", i, CellOrder[i]);
    }
  }

  /*
   * We apply renumbering here because we need it for faces renumbering
   */
  PDM_part_reorder_cell(part, CellOrder);

  /*
   * III/ Create a proper faces order :
   *       -> [ [SubDom1Int/SubDom1Ext] , [SubDom2Int/SubDom2Ext], ..., [SubDomNInt/SubDomNExt],  ]
   *      Renumbering the faces according to the sub block new ordering
   */

  /* Allocate */
  int *nFacePerCache    = (int *) malloc(nBlkCacheWanted * sizeof(int));
  int *nFaceBndPerCache = (int *) malloc(nBlkCacheWanted * sizeof(int));

  /* Init array */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }

  /*
   * First pass :
   *      -> the face is associate with the subdomain that have the lowest color number
   */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];

    int color1 = part->cellColor[iCell1-1];
    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      nFacePerCache[color]++;
    }
    else
    {
      nFaceBndPerCache[color1]++;
    }

    /* Dramatic test */
    if(iCell1 < 0 ){
      printf("PPART internal error \n");
      exit(1);
    }
  }

  /* Second pass :
   *      -> the face is associate with the subdomain that have the lowest color number
   */

  /* Allocate */
  int* partFaceIdx    = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  int* partFaceBndIdx = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));

  /* Initialise */
  partFaceIdx[0]    = 0;
  partFaceBndIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partFaceIdx[i + 1]    = partFaceIdx[i]    + nFacePerCache[i];
    partFaceBndIdx[i + 1] = partFaceBndIdx[i] + nFaceBndPerCache[i];
  }

  for (int i = 0; i < nBlkCacheWanted+1; i++){
    partFaceBndIdx[i] = partFaceBndIdx[i] + partFaceIdx[nBlkCacheWanted];
  }

  /*
   * Verbose
   */
  if(0 == 1)
  {
    printf(" ----------- : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nBlkCacheWanted+1; i++){
      printf("~> %i\n", i);
      printf(" =====> partFaceIdx    %i\n", partFaceIdx[i]);
      printf(" =====> partFaceBndIdx %i\n", partFaceBndIdx[i]);
    }

    for (int i = 0; i < nBlkCacheWanted; i++){
      printf(" =====> nFacePerCache    %i\n", nFacePerCache[i]);
      printf(" =====> nFaceBndPerCache %i\n", nFaceBndPerCache[i]);
    }
  }

  /* Reset Idx */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }

  /* Determine for each subDomain the associate faces */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];

    int color1 = part->cellColor[iCell1-1];

    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      int idx    = partFaceIdx[color]    + nFacePerCache[color];

      FaceOrder[idx] = i;
      nFacePerCache[color]++;
    }
    else
    {
      int color  = color1;
      int idxBnd = partFaceBndIdx[color] + nFaceBndPerCache[color];
      FaceOrder[idxBnd] = i;
      nFaceBndPerCache[color]++;
    }

    if(iCell1 < 0 ){
      printf("PPART internal error \n");
      exit(1);
    }

  } /* End file FaceOrder */

  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrder[%i] = %i \n", i, FaceOrder[i]);
    }
  }

  /*
   * Apply renumbering
   */
  PDM_part_reorder_face(part, FaceOrder);

  // if(part->newToOldOrderFace == NULL){
  //   printf(" SetUp newToOldOrderFace \n");

  //   part->newToOldOrderFace = (int *) malloc (sizeof(int) * nFace);
  //   for (int i = 0; i < nFace; i++){
  //     part->newToOldOrderFace[i] = FaceOrder[i];
  //   }
  // }

  /*
   * Save in array the color of each faces
   */
  if(part->faceColor == NULL)
    part->faceColor = (int *) malloc (sizeof(int) * nFace);

  for (int i = 0; i < nBlkCacheWanted; i++){
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
    for (int iface = partFaceBndIdx[i]; iface < partFaceBndIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
  }


  /*
   * IV/ Create a proper faces order for vectorisation :
   *         -> The idea is to found a pool of faces that not pointed same cells : "independant faces"
   *         Rmk : Boundary faces are not vecotrised
   */
  for (int i = 0; i < part->nFace; i++) {
    FaceOrder[i] = i;
  }

  /* Allocate */
  int *flagCell = (int *) malloc (sizeof(int) * nCell);
  int *flagFace = (int *) malloc (sizeof(int) * nFace);

  /* Core loop of reordenencing */
  for (int i = 0; i < nBlkCacheWanted; i++){

    /* Get local size of SubDomain */
    int nFacSub = nFacePerCache[i];
    // int nCelSub = nCellPerCache[i];

    /* Reset for one pacquet all flag */
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      flagFace[iface] = -1;
    }

    int nFacTreated = 0;
    // printf("nFacSub : %i \n", nFacSub);
    while(nFacTreated != nFacSub){

      /* Reset for one pacquet all flag */
      // for (int iCel = 0; iCel < nCell; iCel++){
      //   flagCell[iCel] = -1;
      // }
      // TEST : A remplacer par : (Beacoup moins couteux)
      for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
        int iCell1 = part->faceCell[2*iface  ];
        int iCell2 = part->faceCell[2*iface+1];
        flagCell[iCell1-1] = -1;
        flagCell[iCell2-1] = -1;
      }

      /* Dramatic verbose */
      if(0 == 1){
        printf("nFacTreated/nFacSub : %i -> %i \n", nFacTreated, nFacSub);
      }

      int nPackFaceLoc = 0;

      /* Loop on face */
      for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){

        if(flagFace[iface] == -1){
          int iCell1 = part->faceCell[2*iface  ];
          int iCell2 = part->faceCell[2*iface+1];

          int t1 = flagCell[iCell1-1];
          int t2 = flagCell[iCell2-1];

          // printf("t1/t2 : %i/%i \n", t1, t2);
          if( (t1 == -1) && (t2 == -1) && nPackFaceLoc < nVectFace ){
            flagCell[iCell1-1] = 1;
            flagCell[iCell2-1] = 1;

            int bFac = partFaceIdx[i];

            FaceOrder[bFac+nFacTreated] = iface;

            flagFace[iface] = 1;

            nFacTreated++;
            nPackFaceLoc++;
          }
        } /* End Face loop */
      } /* End While */
    } /* End SubDom loop */

  }

  /*
   * Apply renumbering -> Attention au tableau en plus à reordencer !
   */
  PDM_part_reorder_face(part, FaceOrder);

  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrderVect[%i] = %i \n", i, FaceOrder[i]);
    }
  }

  /* Copy in partition */
  // if(part->newToOldOrderCell == NULL)
  //   part->newToOldOrderCell = (int *) malloc (sizeof(int) * nCell);

  // for (int i = 0; i < nCell; i++){
  //   part->newToOldOrderCell[i] = CellOrder[i];
  // }

  // if(part->newToOldOrderFace == NULL)
  //   part->newToOldOrderFace = (int *) malloc (sizeof(int) * nFace);

  // for (int i = 0; i < nFace; i++){
  //   part->newToOldOrderFace[i] = FaceOrder[i];
  // }

  if(0 == 1){
    int iMax = -1;
    int iMin = 100000000;

    int* tmpc = (int *) malloc (sizeof(int) * nCell);
    int* tmpf = (int *) malloc (sizeof(int) * nFace);

    for (int i = 0; i < nCell; i++){
      tmpc[i] = -1;
    }
    for (int i = 0; i < nFace; i++){
      tmpf[i] = -1;
    }

    for (int i = 0; i < nCell; i++){
      iMax = PDM_MAX(iMax, CellOrder[i]);
      iMin = PDM_MIN(iMin, CellOrder[i]);
      tmpc[CellOrder[i]] = i;
    }
    printf("iMax = %i / iMin = %i / nCell : %i \n", iMin, iMax, nCell);

    iMax = -1;
    iMin = 100000000;
    for (int i = 0; i < nFace; i++){
      iMax = PDM_MAX(iMax, FaceOrder[i]);
      iMin = PDM_MIN(iMin, FaceOrder[i]);
      tmpf[FaceOrder[i]] = i;
    }
    printf("iMax = %i / iMin = %i / nCell : %i \n", iMin, iMax, nFace);

    for (int i = 0; i < nCell; i++){
      if(tmpc[i] == -1){
        printf("Bizarre ... Cell %i %i \n", i, tmpc[i] );
      }
    }
    for (int i = 0; i < nFace; i++){
      if(tmpf[i] == -1){
        printf("Bizarre ... Face %i %i \n", i, tmpf[i] );
      }
    }
    free(tmpc);
    free(tmpf);

  }

  /* Free memory */
  free(CellOrder);
  free(FaceOrder);
  free(CellCellIdx);
  free(CellCell);
  free(nCellPerCacheBeg);
  free(nCellPerCache);
  free(nFacePerCache);
  free(nFaceBndPerCache);
  free(partFaceIdx);
  free(partFaceBndIdx);
  free(flagCell);
  free(flagFace);

}


/**
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition
 * \param [in]   split_method        Split method
 * \param [in]   nCellPerCacheWanted Approximate number of cells on each cache
 * \param [in]   isAsynchrone        [0 : Synchrone/1 : Asynchrone]
 * \param [in]   isVectorisation     [0 : No/1 : Yes]
 *
 */

void
PDM_renum_cacheblocking2
(
 _part_t     *part,
int           split_method,
int           nCellPerCacheWanted,
int           isAsynchrone,
int           isVectorisation,
int           nVectFace
)
{
  /* Get nFac and nCel */
  const int nCell = part->nCell;
  const int nFace = part->nFace;

  /** Allocate reoerdering/permutation array **/
  int *CellOrder = (int *) malloc (sizeof(int) * nCell);
  int *FaceOrder = (int *) malloc (sizeof(int) * nFace);

  /*
   * I/ Coloring the cells with classic graph library
   */
  int *CellCellIdx = NULL;
  int *CellCell    = NULL;

  assert(part->faceColor == NULL);
  assert(part->cellColor == NULL);



  // int nBlkCacheWanted = _compute_openMP_renumbering(part,
  //                                                   split_method,
  //                                                   nCellPerCacheWanted,
  //                                                   CellOrder,
  //                                                   CellCellIdx,
  //                                                   CellCell);
  int nBlkCacheWanted = _compute_GPU_renumbering(part,
                                                 split_method,
                                                 nCellPerCacheWanted,
                                                 CellOrder,
                                                 CellCellIdx,
                                                 CellCell);

  PDM_part_renum_face (&part, 1, 2, NULL);
  /* Compute graph associate to mesh */
  // free(CellCellIdx);
  // free(CellCell);
  // CellCellIdx = NULL;
  // CellCell = NULL;
  PDM_part_graph_compute_from_face_cell(          part,
                                        (int **) &CellCellIdx,
                                        (int **) &CellCell);

  /*
   * Determine the optimal size for cache blocking
   *   -> The user specify the desired number of cells on each subdomain
   *   -> An other idea is to specify directly the number of block he want
   */
  // int nBlkCacheWanted;
  // if(nCellPerCacheWanted == 0){nBlkCacheWanted = 1;}
  // else                        {nBlkCacheWanted = PDM_MAX(nCell/nCellPerCacheWanted,1);}

  // /* Split the graph */
  // PDM_part_graph_split(          split_method,
  //                                nBlkCacheWanted,
  //                                part,
  //                                CellCellIdx,
  //                                CellCell,
  //                      (int *)   NULL,
  //                      (int *)   NULL,
  //                      (int **) &part->cellColor);

  if(0 == 1)
  {
    printf(" ----------- nBlkCacheWanted : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> part->cellColor[%i]  = %i\n", i, part->cellColor[i]);
    }
  }


  if(0 == 1)
  {
    // printf(" ----------- nThread : %i \n", nThread);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> part->threadColor[%i]  = %i\n", i, part->threadColor[i]);
    }
  }

  /*
   * II/ Create a proper cells order :
   *       -> [ [SubDom1] , [SubDom2], ..., [SubDomN] ]
   */

  /* Allocate */
  int *nCellPerCache    = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int *nCellPerCacheBeg = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int *partCellIdx      = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));

  /* II-Bis - Asynchrone : We change only number of each sub-domain
   * in order to have interior sub-domain first and exterior after
   */
  if(isAsynchrone == 1)
  {
    /* Allocate */
    int *flagHaveBnd     = (int *) malloc(nBlkCacheWanted * sizeof(int));
    int *syncToAsynColor = (int *) malloc(nBlkCacheWanted * sizeof(int));

    _compute_flaghavebnd(part, flagHaveBnd, nBlkCacheWanted);

    /*
     * Built the syncToAsynColor array
     */
    int nSdomB = 0;                     /* Number of   Blocking domain */
    int nSdomU = 0;                     /* Number of Unblocking domain */
    int nextColorB = nBlkCacheWanted-1; /* Next Blocking   color       */
    int nextColorU = 0;                 /* Next UnBlocking color       */

    for(int i = 0; i < nBlkCacheWanted; i++) {
      if(flagHaveBnd[i] == 0)
      {
        syncToAsynColor[i] = nextColorU;
        nextColorU++;
        nSdomU++;
      }
      else
      {
        syncToAsynColor[i] = nextColorB;
        nextColorB--;
        nSdomB++;
      }
    }

    /* Dramatic verbose */
    if(0 == 1)
    {
      printf(" ----------- : %i \n", nBlkCacheWanted);
      for (int i = 0; i < nBlkCacheWanted; i++){
        printf("SyncToAsynColor[%i] = %i\n", i, syncToAsynColor[i]);
      }
    }

    /*
     * Apply the computed array on cellPart
     */
    for(int i = 0; i < part->nCell; i++) {
      part->cellColor[i] = syncToAsynColor[part->cellColor[i]];
    }


    /* Free memory */
    free(flagHaveBnd);
    free(syncToAsynColor);

  } /* End asynchrone */
  /* Init to Zero */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nCellPerCache[i] = 0;
  }

  /* First loop to identify the size of each color */
  for (int i = 0; i < nCell; i++){
    int color = part->cellColor[i]; // A color is a number of partition (output of Metis or Scotch)
    nCellPerCache[color]++;
  }

  /* Compute the index of each Sub-Domain **/
  partCellIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partCellIdx[i + 1] = partCellIdx[i] + nCellPerCache[i];
  }

  if(isVectorisation == 0){
    /** First method - Standard **/

    /* Reset array */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCache[i] = 0;
    }

    for (int i = 0; i < nCell; i++){
      int color  = part->cellColor[i];
      int idx    = partCellIdx[color] + nCellPerCache[color];

      CellOrder[idx] = i;
      nCellPerCache[color]++;
    }
  }
  else
  {
    /** Second method - Vectorisation on cell on current domain **/

    /* Store the current Index to add a vectorisable cell */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCacheBeg[i] = 0;
    }

    /* Loop on cells :
     *     -> If current cell is adjacent to another color : add to end
     *     -> Else add to begin
     * In fact this function sort the interior cell and exterior cell leads to the following numbering :
     *         [ [ Interior Cells] , [ Exterior Cells]]
     */
    int verifCell = 0;
    for (int i = 0; i < nCell; i++){

      /* Get color and prepare flag */
      int color = part->cellColor[i];
      int flag  = -1;
      for(int j = CellCellIdx[i]; j < CellCellIdx[i+1]; j++){
        int iCell = CellCell[j];
        // if(part->cellColor[iCell] != color){flag = 1;} // Alors on est au bord d'un sous domaine !
        if(part->cellColor[iCell] < color){flag = 1;} // Alors on est au bord d'un sous domaine !
      }

      if(flag == -1){  // Cell is interior : add to begin

        int idx = partCellIdx[color] + nCellPerCacheBeg[color];

        CellOrder[idx] = i;
        nCellPerCacheBeg[color]++;
        verifCell++;

      }
      else{ // Cell is exterior : add to end

        int idx = partCellIdx[color] + nCellPerCache[color]-1;

        CellOrder[idx] = i;
        nCellPerCache[color]--;
        verifCell++;
      }
    }
  }


  /*
   * Verbose
   */
  if(0 == 1)
  {
    printf(" ----------- nBlkCacheWanted : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> CellOrder[%i]  = %i\n", i, CellOrder[i]);
    }
  }

  /*
   * We apply renumbering here because we need it for faces renumbering
   */
  PDM_part_reorder_cell(part, CellOrder);
  free(CellCellIdx);
  free(CellCell);
  CellCellIdx = NULL;
  CellCell = NULL;
  PDM_part_graph_compute_from_face_cell(          part,
                                        (int **) &CellCellIdx,
                                        (int **) &CellCell);

  for (int i = 0; i < nCell; i++){
    CellOrder[i] = -1;
  }
  if(part->threadColor == NULL){
    part->threadColor = (int *) malloc (sizeof(int) * nCell);
  }

  // printf("Shut threadColor \n");
  // for (int icell = 0; icell < nCell; icell++)
  // {
  //   part->threadColor[icell] = 0;
  // }

  /*
   * III/ Create a proper faces order :
   *       -> [ [SubDom1Int/SubDom1Ext] , [SubDom2Int/SubDom2Ext], ..., [SubDomNInt/SubDomNExt],  ]
   *      Renumbering the faces according to the sub block new ordering
   */

  /* Allocate */
  int *nFacePerCache    = (int *) malloc(nBlkCacheWanted * sizeof(int));
  int *nFaceBndPerCache = (int *) malloc(nBlkCacheWanted * sizeof(int));

  /* Init array */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }

  /*
   * First pass :
   *      -> the face is associate with the subdomain that have the lowest color number
   */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];

    int color1 = part->cellColor[iCell1-1];
    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      nFacePerCache[color]++;
    }
    else
    {
      nFaceBndPerCache[color1]++;
    }

    /* Dramatic test */
    if(iCell1 < 0 ){
      printf("PPART internal error \n");
      exit(1);
    }
  }

  /* Second pass :
   *      -> the face is associate with the subdomain that have the lowest color number
   */

  /* Allocate */
  int* partFaceIdx    = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  int* partFaceBndIdx = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));

  /* Initialise */
  partFaceIdx[0]    = 0;
  partFaceBndIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partFaceIdx[i + 1]    = partFaceIdx[i]    + nFacePerCache[i];
    partFaceBndIdx[i + 1] = partFaceBndIdx[i] + nFaceBndPerCache[i];
  }

  for (int i = 0; i < nBlkCacheWanted+1; i++){
    partFaceBndIdx[i] = partFaceBndIdx[i] + partFaceIdx[nBlkCacheWanted];
  }

  /*
   * Verbose
   */
  if(0 == 1)
  {
    printf(" ----------- : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nBlkCacheWanted+1; i++){
      printf("~> %i\n", i);
      printf(" =====> partFaceIdx    %i\n", partFaceIdx[i]);
      printf(" =====> partFaceBndIdx %i\n", partFaceBndIdx[i]);
    }

    for (int i = 0; i < nBlkCacheWanted; i++){
      printf(" =====> nFacePerCache    %i\n", nFacePerCache[i]);
      printf(" =====> nFaceBndPerCache %i\n", nFaceBndPerCache[i]);
    }
  }

  /* Reset Idx */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }

  /* Determine for each subDomain the associate faces */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];

    int color1 = part->cellColor[iCell1-1];

    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      int idx    = partFaceIdx[color]    + nFacePerCache[color];

      FaceOrder[idx] = i;
      nFacePerCache[color]++;
    }
    else
    {
      int color  = color1;
      int idxBnd = partFaceBndIdx[color] + nFaceBndPerCache[color];
      FaceOrder[idxBnd] = i;
      nFaceBndPerCache[color]++;
    }

    if(iCell1 < 0 ){
      printf("PPART internal error \n");
      exit(1);
    }

  } /* End file FaceOrder */

  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrder[%i] = %i \n", i, FaceOrder[i]);
    }
  }

  /*
   * Apply renumbering
   */
  PDM_part_reorder_face(part, FaceOrder);

  /*
   * Save in array the color of each faces
   */
  if(part->faceColor == NULL)
    part->faceColor = (int *) malloc (sizeof(int) * nFace);

  for (int i = 0; i < nBlkCacheWanted; i++){
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
    for (int iface = partFaceBndIdx[i]; iface < partFaceBndIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
  }

  _compute_multi_coloring(part, CellCellIdx, CellCell, CellOrder, partCellIdx, nBlkCacheWanted);

  // PDM_part_reorder_cell(part, CellOrder);

  /*
   * IV/ Create a proper faces order for vectorisation :
   *         -> The idea is to found a pool of faces that not pointed same cells : "independant faces"
   *         Rmk : Boundary faces are not vecotrised
   */
  for (int i = 0; i < part->nFace; i++) {
    FaceOrder[i] = i;
  }

  /* Allocate */
  int *flagCell = (int *) malloc (sizeof(int) * nCell);
  int *flagFace = (int *) malloc (sizeof(int) * nFace);

  /* Core loop of reordenencing */
  // for (int i = 0; i < nBlkCacheWanted; i++){

  //   /* Get local size of SubDomain */
  //   int nFacSub = nFacePerCache[i];
  //   // int nCelSub = nCellPerCache[i];

  //   /* Reset for one pacquet all flag */
  //   for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
  //     flagFace[iface] = -1;
  //   }

  //   int nFacTreated = 0;
  //   // printf("nFacSub : %i \n", nFacSub);
  //   while(nFacTreated != nFacSub){

  //     /* Reset for one pacquet all flag */
  //     for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
  //       int iCell1 = part->faceCell[2*iface  ];
  //       int iCell2 = part->faceCell[2*iface+1];
  //       flagCell[iCell1-1] = -1;
  //       flagCell[iCell2-1] = -1;
  //     }

  //     /* Dramatic verbose */
  //     if(0 == 1){
  //       printf("nFacTreated/nFacSub : %i -> %i \n", nFacTreated, nFacSub);
  //     }

  //     int nPackFaceLoc = 0;

  //     /* Loop on face */
  //     for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){

  //       if(flagFace[iface] == -1){
  //         int iCell1 = part->faceCell[2*iface  ];
  //         int iCell2 = part->faceCell[2*iface+1];

  //         int t1 = flagCell[iCell1-1];
  //         int t2 = flagCell[iCell2-1];

  //         // printf("t1/t2 : %i/%i \n", t1, t2);
  //         if( (t1 == -1) && (t2 == -1) && nPackFaceLoc < nVectFace ){
  //           flagCell[iCell1-1] = 1;
  //           flagCell[iCell2-1] = 1;

  //           int bFac = partFaceIdx[i];

  //           FaceOrder[bFac+nFacTreated] = iface;

  //           flagFace[iface] = 1;

  //           nFacTreated++;
  //           nPackFaceLoc++;
  //         }
  //       } /* End Face loop */
  //     } /* End While */
  //   } /* End SubDom loop */

  // }

  /*
   * Apply renumbering -> Attention au tableau en plus à reordencer !
   */
  // PDM_part_reorder_face(part, FaceOrder);

  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrderVect[%i] = %i \n", i, FaceOrder[i]);
    }
  }

  /*
   * Free
   */

  /* Copy in partition */
  /* A reflechir ici avec le multigrille on a aplliquer plusieurs numérotation mais on prendre la derniere permutation ?? */
  // if(part->newToOldOrderCell == NULL)
  //   part->newToOldOrderCell = (int *) malloc (sizeof(int) * nCell);

  // for (int i = 0; i < nCell; i++){
  //   part->newToOldOrderCell[i] = CellOrder[i];
  // }

  // if(part->newToOldOrderFace == NULL)
  //   part->newToOldOrderFace = (int *) malloc (sizeof(int) * nFace);

  // for (int i = 0; i < nFace; i++){
  //   part->newToOldOrderFace[i] = FaceOrder[i];
  // }

  if(0 == 1){
    int iMax = -1;
    int iMin = 100000000;

    int* tmpc = (int *) malloc (sizeof(int) * nCell);
    int* tmpf = (int *) malloc (sizeof(int) * nFace);

    for (int i = 0; i < nCell; i++){
      tmpc[i] = -1;
    }
    for (int i = 0; i < nFace; i++){
      tmpf[i] = -1;
    }

    for (int i = 0; i < nCell; i++){
      iMax = PDM_MAX(iMax, CellOrder[i]);
      iMin = PDM_MIN(iMin, CellOrder[i]);
      tmpc[CellOrder[i]] = i;
    }
    printf("iMax = %i / iMin = %i / nCell : %i \n", iMin, iMax, nCell);

    iMax = -1;
    iMin = 100000000;
    for (int i = 0; i < nFace; i++){
      iMax = PDM_MAX(iMax, FaceOrder[i]);
      iMin = PDM_MIN(iMin, FaceOrder[i]);
      tmpf[FaceOrder[i]] = i;
    }
    printf("iMax = %i / iMin = %i / nCell : %i \n", iMin, iMax, nFace);

    for (int i = 0; i < nCell; i++){
      if(tmpc[i] == -1){
        printf("Bizarre ... Cell %i %i \n", i, tmpc[i] );
      }
    }
    for (int i = 0; i < nFace; i++){
      if(tmpf[i] == -1){
        printf("Bizarre ... Face %i %i \n", i, tmpf[i] );
      }
    }
    free(tmpc);
    free(tmpf);

  }

  /* Free memory */
  free(CellOrder);
  free(FaceOrder);
  free(CellCellIdx);
  free(CellCell);
  free(nCellPerCacheBeg);
  free(nCellPerCache);
  free(nFacePerCache);
  free(nFaceBndPerCache);
  free(partFaceIdx);
  free(partFaceBndIdx);
  free(flagCell);
  free(flagFace);
  free(partCellIdx);
}

/*
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition
 *
 */
void
PDM_renum_cacheblocking_compute_loop_array
(
 _part_t     *part
)
{
  if(part->faceColor == NULL){PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : You need to faceColor \n");}
  if(part->cellColor == NULL){PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : You need to cellColor \n");}

  if(part->subpartlayout == NULL)
  {
    part->subpartlayout = (_subpartlayout_t *) malloc( sizeof(_subpartlayout_t));

    part->subpartlayout->cellTileIdx = NULL;
    part->subpartlayout->faceTileIdx = NULL;
    part->subpartlayout->faceBndTileIdx = NULL;
    part->subpartlayout->maskTileIdx = NULL;
    part->subpartlayout->cellVectTileIdx = NULL;
    part->subpartlayout->maskTileN = NULL;
    part->subpartlayout->cellVectTileN = NULL;
    part->subpartlayout->maskTile = NULL;
  }

  /* Prepare sub-domain */
  _prepare_subdomain(part, part->subpartlayout);

  /* Compute sub-domain */
  _compute_subdomain(part, part->subpartlayout);

  /* Compute mask */
  _compute_mask(part, part->subpartlayout);


}

#ifdef  __cplusplus
}
#endif
