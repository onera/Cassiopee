
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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multicoloring.h"

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
static int
_get_color
(
 int  *saveColor,
 int  nColorFound,
 int  firstColor,
 int *lastColor
)
{
  int iColor = -1;
  PDM_sort_int(saveColor, NULL, nColorFound);


  // printf(" firstColor : %i \n", firstColor);
  // printf(" lastColor  : %i \n", *lastColor);
  saveColor[nColorFound] = saveColor[nColorFound-1];

  // Faire une recherche simple
  int iColorMax = saveColor[nColorFound-1];
  int iColorMin = saveColor[0];

  if(nColorFound == 0){
    iColor = firstColor;
    if(  (*lastColor) == -1 )
      (*lastColor)++;
  }
  else if(iColorMin > firstColor)
  {
    iColor = firstColor;
    // printf(" Cas 1 : %i/%i\n", iColor, firstColor);
  }
  else
  {
    // Search gap
    for(int isc = 0; isc < nColorFound; isc++){
      if( (saveColor[isc+1] - saveColor[isc] ) > 1 )
      {
        iColor = saveColor[isc+1]-1; // Because Gap
        // printf(" Cas 3 %i \n", iColor);
        break;
      }
    }

  }


  if( (iColor == -1) && (iColorMax < *lastColor))
  {
    iColor = *lastColor;
    // printf(" Cas 2 : %i/%i\n", iColor, *lastColor);
  }

  if(iColor == -1 )
  {
    // printf(" Add color T1 : %i \n", *lastColor);
    (*lastColor)++;
    iColor = *lastColor;
  }
  else if( iColor == *lastColor + 1)
  {
    // printf(" Add color T2 : %i \n", *lastColor);
    abort();
    (*lastColor)++;
    iColor = *lastColor;
  }

  if(1 == 1)
  {
    // printf(" iColor = %i -----> saveColor : ", iColor);
    // printf(" -----> saveColor : ");
    // for(int isc = 0; isc < nColorFound; isc++){
    //   printf(" %i ", saveColor[isc]);
    // }
    // printf("\n");

    for(int isc = 0; isc < nColorFound; isc++){
      assert(saveColor[isc] != iColor);
    }
  }

  return iColor;
}

/**
 *
 * \brief  Usage
 *
 */
/* static */
/* int */
/* compute_local_multicoloring */
/* ( */
/*  PDM_MPI_Comm comm, */
/*  int          sizeG, */
/*  int*         ia, */
/*  PDM_g_num_t* ja, */
/*  PDM_g_num_t* shiftG, */
/*  int*         color, */
/*  int*         saveColor, */
/*  int*         cpl_strid_idx, */
/*  PDM_g_num_t* cpl_data */
/* ) */
/* { */
/*   int nColor = 0; */
/*   int firstColor = 0; */
/*   printf("compute_local_multicoloring\n"); */

/*   /\* MPI Stuff *\/ */
/*   int iRank, nRank; */
/*   PDM_MPI_Comm_rank(comm, &iRank); */
/*   PDM_MPI_Comm_size(comm, &nRank); */

/*   /\* for all in the graph *\/ */
/*   for(int icol1 = 0; icol1 < sizeG; icol1++){ */

/*     int nColorFound = 0; */

/*     int beg1 = ia[icol1  ]; */
/*     int end1 = ia[icol1+1]; */

/*     /\* First search inside domain dependancy *\/ */
/*     for(int icol2 = beg1; icol2 < end1; icol2++){ */
/*       PDM_g_num_t iCellG = ja[icol2]; */
/*       // printf(" iCellG : %i \n", iCellG); */
/*       int oppRank  = PDM_binary_search_gap_long(iCellG, shiftG, nRank+1); */
/*       if(oppRank == iRank){ */
/*         // Tt est local */
/*         /\* Search the rank 2 *\/ */
/*         int beg2 = ia[ja[icol2]  -1]; */
/*         int end2 = ia[ja[icol2]+1-1]; */
/*         for(int icol3 = beg2; icol3 < end2; icol3++){ */
/*           // printf(" color[%i] = %i \n", ja[icol3]-1, color[ja[icol3]-1]); */
/*           //PDM_g_num_t iCellG2 = ja[icol3]; */
/*           // printf(" iCellG2 : %i | icol3 = %i \n", iCellG2, icol3); */
/*           // int oppRank2 = PDM_binary_search_gap_long(iCellG2, shiftG, nRank+1); */

/*           // if(oppRank2 == iRank){ */
/*           //   if(color[ja[icol3]-1] != -1){ */
/*           //     saveColor[nColorFound++] = color[ja[icol3]-1]; */
/*           //   } */
/*           // } */
/*           // } else { */
/*           //   printf("Manage coupling \n"); */
/*           // } */
/*         } */
/*       } */
/*     } */


/*     /\* Here saveColor contains all possible color aroung icol1 */
/*      *    We need to found out the proper color in order to be independant of the other around */
/*      *\/ */
/*     int lColor = _get_color(saveColor, nColorFound, firstColor, &nColor); */

/*     color[icol1] = lColor; */
/*     nColor = PDM_MAX(nColor, lColor); */


/*   } */

/*   return nColor+1; */
/* }; */

/**
 *
 * \brief  Usage
 *
 */
static
void
compute_local_multicoloringT
(
 PDM_MPI_Comm comm,
 int          sizeG,
 int*         ia,
 PDM_g_num_t* ja,
 PDM_g_num_t* shiftG,
 PDM_g_num_t* LNToGNCoupling,
 int          nCoupling,
 int*         color,
 int*         saveColor,
 int*         toTreat,
 int*         nToTreat,
 int*         nColor
)
{
  // int nColor = 0;
  // int nColor = -1;
  int firstColor = 0;
  printf("compute_local_multicoloring\n");

  /* MPI Stuff */
  int iRank, nRank;
  PDM_MPI_Comm_rank(comm, &iRank);
  PDM_MPI_Comm_size(comm, &nRank);

  // int positionInLNToGN[sizeG];

  /* for all in the graph */
  int nToTreatCur = (*nToTreat);
  (*nToTreat) = 0;
  //int idxToTreat = 0;
  for(int ii = 0; ii < nToTreatCur; ii++){

    int icol1 = toTreat[ii];
    int nColorFound = 0;

    int beg1 = ia[icol1  ];
    int end1 = ia[icol1+1];

    /* First search inside domain dependancy */
    int allInSameRank = 1;
    int allInUpperOrSameRank = 1;
    for(int icol2 = beg1; icol2 < end1; icol2++){
      PDM_g_num_t iCellG = ja[icol2];
      // printf(" iCellG : %i \n", iCellG);
      int oppRank  = PDM_binary_search_gap_long(iCellG, shiftG, nRank+1);
      if(iRank != oppRank){
        allInSameRank = 0;
        // break;
      }
      if(oppRank < iRank){
        allInUpperOrSameRank = 0;
      }
    }

    // Quid du cas ou   3 est connecté a 0 et 5 ????
    if(allInSameRank == 1){
      for(int icol2 = beg1; icol2 < end1; icol2++){
        int col3 = ja[icol2] - shiftG[iRank];
        if(color[col3] != -1){
           saveColor[nColorFound++] = color[col3];
        }
      }
    }
    else if ( allInUpperOrSameRank == 1){ // All is upper or egal to current Rank --> Priority for coloring

      for(int icol2 = beg1; icol2 < end1; icol2++){
        PDM_g_num_t iCellG = ja[icol2];
        int oppRank  = PDM_binary_search_gap_long(iCellG, shiftG, nRank+1);
        if( iRank == oppRank ){
          int col3 = ja[icol2] - shiftG[iRank];
          if(color[col3] != -1){
            saveColor[nColorFound++] = color[col3];
          }
        }
        else {
          int idx = PDM_binary_search_long(iCellG, LNToGNCoupling, nCoupling);
          int col3 = sizeG + idx;
          if(color[col3] != -1){
            saveColor[nColorFound++] = color[col3];
          }
        }
      }
    }

    int isColorisable = 1;
    if( allInSameRank == 0 && allInUpperOrSameRank != 1){

      for(int icol2 = beg1; icol2 < end1; icol2++){
        PDM_g_num_t iCellG = ja[icol2];
        int oppRank  = PDM_binary_search_gap_long(iCellG, shiftG, nRank+1);
        if( iRank == oppRank ){
          int col3 = ja[icol2] - shiftG[iRank];
          if(color[col3] != -1){
            saveColor[nColorFound++] = color[col3];
          }
        }
        else if( iRank > oppRank ){
          int idx = PDM_binary_search_long(iCellG, LNToGNCoupling, nCoupling);
          int col3 = sizeG + idx;
          if(color[col3] != -1){
            saveColor[nColorFound++] = color[col3];
          }
          else{
            isColorisable = 0;
          }
        }
      }
    }

    // Pour pouvoir colorier il faut que :
    //   --> Le rang est inférieur et tout les couleur de couplage sont == -1 !
    //   --> Tt les rang infériuer sont traités ---> ie tt les couplages sont != -1
    //   --> Tt est local et indépendant et color == -1
    //
    //   Reflechir a traiter uniquement des cellules non traité !!!
    //        for(int icol1 = 0; icol1 < sizeG; icol1++){
    // -->    for(int i = 0; i < nUncolorLoc; i++){
    //             idx = toColor[i]


    /* Here saveColor contains all possible color aroung icol1
     *    We need to found out the proper color in order to be independant of the other around
     */
    printf("[%i] [sR = %i / aU = %i / iC = %i\n", iRank, allInSameRank, allInUpperOrSameRank, isColorisable);
    if( allInSameRank == 1 || allInUpperOrSameRank == 1 || isColorisable == 1){
      // int lColor = _get_color(saveColor, nColorFound, firstColor, &nColor);
      int lColor = _get_color(saveColor, nColorFound, firstColor, nColor);

      color[icol1] = lColor;
      (*nColor) = PDM_MAX((*nColor), lColor);
    }
    else
    {
      toTreat[(*nToTreat)++] = icol1;
    }
  }
};


/**
 *
 * \brief  Usage
 *
 */
static
void
compute_distribute_multi_coloring
(
 _dist_csr*        dcsr,
 _dist_csr*        dcsrNext
)
{
  printf(" compute_distribute_multi_coloring \n");

  /* MPI Stuff */
  int iRank, nRank;
  PDM_MPI_Comm_rank(dcsr->comm, &iRank);
  PDM_MPI_Comm_size(dcsr->comm, &nRank);

  /* Compute coupling */
  PDM_g_num_t* LNToGN;
  // Pour l'instant avec Ghost + Interior ---> To Optim after
  //  Je pense que pour le part_to_block on a besoin de tt !!
  int nElmtWithGhost = PDM_generate_part_LNToGN(dcsrNext, &LNToGN, 1);
  // int nCoupling = PDM_generate_part_LNToGN(dcsrNext, &LNToGN, 0);

  /* Now we have LNToGN for the block_to_part in order to update color */
  PDM_g_num_t dCellCell[nRank+1];
  for(int i = 0; i < nRank+1; i++){
    dCellCell[i] = dcsr->shiftG[i]-1;
    // printf(" dCellCell[%i] = %i\n", i, dCellCell[i]);
  }
  PDM_g_num_t* ptLNToGNBnd = NULL;
  int nCoupling = (int) (nElmtWithGhost - dcsr->lSize);
  printf(" nCoupling = %d \n", nCoupling);
  if(nCoupling > 0){
    ptLNToGNBnd = &LNToGN[dcsr->lSize];
  }
  // PDM_block_to_part_t *btp = PDM_block_to_part_create( dCellCell,
  //                                                     &LNToGN,
  //                                                     &nElmtWithGhost,
  //                                                      1,
  //                                                      dcsr->comm);
  PDM_block_to_part_t *btp = PDM_block_to_part_create( dCellCell,
                                                     (const PDM_g_num_t **) &ptLNToGNBnd,
                                                      &nCoupling,
                                                       1,
                                                       dcsr->comm);


  PDM_part_to_block_t* ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       0.,
                                                       &LNToGN,
                                                       dCellCell,
                                                       &nElmtWithGhost,
                                                       1,
                                                       dcsr->comm);

  /* First exchange for color */
  // int color[nElmtWithGhost*8];
  // int color_strid[nElmtWithGhost];
  int* color       = (int *) malloc( nElmtWithGhost * sizeof(int));
  int* color_strid = (int *) malloc( nElmtWithGhost * sizeof(int));
  int saveColor[dcsr->lSize];
  int toTreat[dcsr->lSize];

  printf("nElmtWithGhost = %i \n", nElmtWithGhost);
  /* Init phase */
  for(int i = 0; i < dcsr->lSize; i++){
    color[i] = -1;
    color_strid[i] = 1;
    toTreat[i] = i;
  }
  for(int i = dcsr->lSize; i < nElmtWithGhost; i++){
    color[i] = -1;
    color_strid[i] = 1;
  }

  int nColor = -1;
  int nToTreat = dcsr->lSize;

  int* ptColorBnd = NULL;
  if(nCoupling > 0){
    ptColorBnd = &color[dcsr->lSize];
  }

  int* part_data_color = (int *) malloc( nCoupling * sizeof(int));

  // for(int iStep = 0; iStep < 2; iStep++){
  int nToTreatG = dcsr->gSize;
  int iStep = 0;
  while(nToTreatG != 0){

    iStep++;
    // if(iStep > 4){
    //   abort();
    // }
    /* Local multicoloring with respect to Coupling */
    compute_local_multicoloringT(dcsrNext->comm,
                                 dcsrNext->lSize,
                                 dcsrNext->ia,
                                 dcsrNext->ja,
                                 dcsrNext->shiftG,
                                 ptLNToGNBnd,
                                 nCoupling,
                                 color,
                                 saveColor,
                                 toTreat,
                                 &nToTreat,
                                 &nColor);
    printf("[%i] nColor+1 = %i \n", iRank, nColor+1);
    printf("[%i] nToTreat = %i \n", iRank, nToTreat);

    int stride_one = 1;
    PDM_block_to_part_exch(btp,
                           sizeof(int),
                           PDM_STRIDE_CST,
                           &stride_one,
                           (void*) color,
                           NULL,
                           (void **) &ptColorBnd);

    if( 0 == 1){
      fflush(stdout);
      fflush(stdout);
      PDM_MPI_Barrier(dcsr->comm);
      for(int jj = 0; jj < nRank; jj++){
        fflush(stdout);
        fflush(stdout);
        PDM_MPI_Barrier(dcsr->comm);
        if(jj == iRank){
          for(int i = 0; i < nCoupling; i++){
            printf(" part_data[%i] = %i | color = %i \n", i, ptColorBnd[i], color[i]);
          }
        }
        fflush(stdout);
        PDM_MPI_Barrier(dcsr->comm);
      }
    }

    int nColorMax = -1;
    PDM_MPI_Allreduce (&nColor, &nColorMax, 1, PDM_MPI_INT, PDM_MPI_MAX, dcsr->comm);
    nColor = nColorMax;

    int nToTreatLoc = nToTreat;
    PDM_MPI_Allreduce (&nToTreatLoc, &nToTreatG, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dcsr->comm);

    fflush(stdout);
    fflush(stdout);
    PDM_MPI_Barrier(dcsr->comm);
    printf("[%i] nToTreatG = %i n", iRank, nToTreatG);
    fflush(stdout);
    fflush(stdout);
    PDM_MPI_Barrier(dcsr->comm);

  }

  if( 1 == 1){
    printf("[%i] Check coloring ------------------------- iStep = %i \n", iRank, iStep);
    fflush(stdout);
    fflush(stdout);
    PDM_MPI_Barrier(dcsr->comm);
    for(int jj = 0; jj < nRank; jj++){
      fflush(stdout);
      fflush(stdout);
      PDM_MPI_Barrier(dcsr->comm);
      if(jj == iRank){
        for(int i = 0; i < dcsr->lSize; i++){
          printf("[%i] color[%i] = %i \n", iRank, i, color[i]);
        }
        fflush(stdout);
      }
      fflush(stdout);
      PDM_MPI_Barrier(dcsr->comm);
    }
    printf("[%i] Check coloring ------------------------- end \n", iRank);
  }


  /*
   * Free
   */
  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  free(color);
  free(color_strid);
  free(part_data_color);
  free(LNToGN);

}



  // int* blk_color_strid;
  // int* blk_color_data;
  // int nMaxSize = PDM_part_to_block_exch(         ptb,
  //                                                sizeof(int),
  //                                                PDM_STRIDE_VAR,
  //                                                -1,
  //                                                &color_strid,
  //                                       (void**) &color,
  //                                                &blk_color_strid,
  //                                       (void**) &blk_color_data);


  // int blk_color_size = PDM_part_to_block_n_elt_block_get(ptb);

  // /*
  //  * Verbose
  //  */
  // if(1 == 1){
  //   fflush(stdout);
  //   PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
  //   for(int jj = 0; jj < nRank; jj++){
  //     printf("blk_color_strid : %d\n", blk_color_size);
  //     printf("nMaxSize : %d\n", nMaxSize);
  //     PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
  //     fflush(stdout);
  //     if(jj == iRank){
  //       int idx = 0;
  //       for(int i = 0; i < blk_color_size; i++) {
  //         printf("blk_color_strid[%d] : %d | Data -> ", i, blk_color_strid[i]);
  //         for(int ii = 0; ii < blk_color_strid[i]; ii++){
  //           printf("%i ", (int)blk_color_data[ii+idx]);
  //           // printf("%i ", jj+idx);
  //         }
  //         printf("\n");
  //         idx += blk_color_strid[i];
  //       }
  //     }
  //   }
  //   fflush(stdout);
  //   PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
  // }

  // // Refill color
  // assert(blk_color_size == dcsr->lSize);
  // int idxB = 0;
  // for(int i = 0; i < blk_color_size; i++) {
  //   /* 3 choses à verifier
  //    *     --> Tt le monde à -1
  //    *     --> Si il y a un non égale à -1 -> tt le monde doit avoir la meme valuer
  //    */
  //   int lastColor = -1;
  //   int newColor  = -1;
  //   for(int ii = 0; ii < blk_color_strid[i]; ii++){
  //     newColor = blk_color_data[ii+idxB];
  //     if( lastColor == -1){
  //       if( newColor != -1){
  //         lastColor = newColor;
  //       }
  //     }
  //     else{
  //       if(newColor != -1){
  //         assert(newColor == lastColor); // Avoid coloring of two cells comming from two different porc
  //       }
  //     }
  //   }
  //   color[i] = newColor;
  //   idxB += blk_color_strid[i];
  // }
  // free(blk_color_strid);
  // free(blk_color_data);

  // // Second block_to_part to exchange stencil rank 1
  // printf(" Second block_to_part to exchange stencil rank 1 \n ");

  // // block_to_part works with stride not index ...
  // int* block_strid = (int *) malloc( dcsr->lSize * sizeof(int));
  // for(int i = 0; i < dcsr->lSize; i++){
  //   block_strid[i] = dcsr->ia[i+1] - dcsr->ia[i];
  // }
  // int**         cpl_stride;
  // PDM_g_num_t** cpl_data;
  // PDM_block_to_part_exch2 (btp,
  //                          sizeof(PDM_g_num_t),
  //                          PDM_STRIDE_VAR,
  //                          block_strid,
  //                          (void*) dcsr->ja,
  //                          &cpl_stride,
  //                          (void ***) &cpl_data);
  // free(block_strid);
  // printf(" Second block_to_part to exchange stencil rank 1 end \n ");

  // int cpl_strid_idx[nElmtWithGhost+1];
  // cpl_strid_idx[0] = 0;
  // for (int i = 0; i < nElmtWithGhost; i++) {
  //   cpl_strid_idx[i+1] = cpl_stride[0][i] + cpl_strid_idx[i];
  // }

  // if( 1 == 1){
  //   fflush(stdout);
  //   PDM_MPI_Barrier(dcsr->comm);
  //   for(int jj = 0; jj < nRank; jj++){
  //     fflush(stdout);
  //     PDM_MPI_Barrier(dcsr->comm);
  //     if(jj == iRank){
  //       // for(int i = 0; i < nElmtWithGhost+1; i++){
  //       for(int i = 0; i < nCoupling+1; i++){
  //         printf(" cpl_strid[%i] = %i \n ", i, cpl_strid_idx[i]);
  //         // printf(" dcsr->ia[%i]    = %i \n ", i, dcsr->ia[i]);
  //       }
  //     }
  //     fflush(stdout);
  //     PDM_MPI_Barrier(dcsr->comm);
  //   }
  // }
/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */
void
PDM_compute_distributed_multicoloring
(
 _dist_csr* dcsr
)
{
  printf("compute_distribute_multicoloring\n");

  // int nColor = compute_multicoloring(dcsr->lSize,
  //                                    dcsr->ia,
  //                                    dcsr->ja,
  //                                    color,
  //                                    saveColor);
  // printf("compute_distribute_multicoloring : nColor = %i \n", nColor);

  // PDM_g_num_t* LNToGN;
  // int nElmtWithGhost = prepare_distribute_multi_coloring(dcsr, &LNToGN);
  // _dist_csr* dcsrNext = generate_superior_rank_csr(dcsr, LNToGN, nElmtWithGhost);
  // LNToGNNext to comute also
  _dist_csr* dcsrNext = PDM_dist_csr_gen_superior_rank(dcsr);

  PDM_dist_csr_print(dcsrNext);
  PDM_dist_csr_dump(dcsrNext, "dcsr_O2");
  // Beautiful this one works !!
  // _dist_csr* dcsrNext2 = PDM_dist_csr_gen_superior_rank(dcsrNext);
  // PDM_dist_csr_dump(dcsrNext2, "dcsr_O3");
  compute_distribute_multi_coloring(dcsr, dcsrNext);

  PDM_dist_csr_free(dcsrNext);
  // free(LNToGN);
};

/**
 *
 * \brief  Usage
 *
 */
int compute_multicoloring
(
 int          sizeG,
 int*         ia,
 PDM_g_num_t* ja,
 int*         color,
 int*         saveColor
)
{
  int nColor = 0;
  int firstColor = 0;
  printf("compute_multicoloring\n");

  for(int i = 0; i < sizeG; i++){
    color[i] = -1;
  }

  /* for all in the graph */
  for(int icol1 = 0; icol1 < sizeG; icol1++){

    int nColorFound = 0;

    int beg1 = ia[icol1  ];
    int end1 = ia[icol1+1];

    /* Search the rank 2 */
    for(int icol2 = beg1; icol2 < end1; icol2++){
      int beg2 = ia[ja[icol2]  -1];
      int end2 = ia[ja[icol2]+1-1];

      for(int icol3 = beg2; icol3 < end2; icol3++){
        // printf(" color[%i] = %i \n", ja[icol3]-1, color[ja[icol3]-1]);
        if(color[ja[icol3]-1] != -1){
           saveColor[nColorFound++] = color[ja[icol3]-1];
        }
      }
    }

    /* Here saveColor contains all possible color aroung icol1
     *    We need to found out the proper color in order to be independant of the other around
     */
    int lColor = _get_color(saveColor, nColorFound, firstColor, &nColor);

    color[icol1] = lColor;
    nColor = PDM_MAX(nColor, lColor);


  }

  return nColor+1;
};



/*============================================================================
 * Private function definitions
 *============================================================================*/


#ifdef  __cplusplus
}
#endif
