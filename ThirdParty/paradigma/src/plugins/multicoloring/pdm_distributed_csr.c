
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
#include "pdm_part_renum.h"
#include "pdm_order.h"
#include "pdm_part_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_distributed_csr.h"

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


/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief
 *
 */
_dist_csr*
PDM_dist_csr_create
(
  const PDM_MPI_Comm           comm
)
{
  _dist_csr* dCsr = (_dist_csr *) malloc( sizeof(_dist_csr) );

  dCsr->comm  = comm;

  dCsr->lSize  = -1;
  dCsr->gSize  = -1;
  dCsr->ia     = NULL;
  dCsr->ja     = NULL;
  dCsr->shiftG = NULL;
  dCsr->dnnzIdx= NULL;

  return dCsr;
}

/**
 *
 * \brief
 *
 */
void
PDM_dist_csr_dump
(
 _dist_csr* dcsr,
 char*      filename
)
{
  printf("PDM_dist_csr_dump \n");
  int iRank, nRank;
  PDM_MPI_Comm_rank(dcsr->comm, &iRank);
  PDM_MPI_Comm_size(dcsr->comm, &nRank);

  const char* extension_ia = "_ia.bin";
  const char* extension_ja = "_ja.bin";

  char* name_with_extension_ia;
  char* name_with_extension_ja;
  name_with_extension_ia = malloc(strlen(filename)+1+strlen(extension_ia)); /* make space for the new string (should check the return value ...) */
  name_with_extension_ja = malloc(strlen(filename)+1+strlen(extension_ja)); /* make space for the new string (should check the return value ...) */
  strcpy(name_with_extension_ia, filename);                                 /* copy name into the new var */
  strcpy(name_with_extension_ja, filename);                                 /* copy name into the new var */
  strcat(name_with_extension_ia, extension_ia);                             /* add the extension */
  strcat(name_with_extension_ja, extension_ja);                             /* add the extension */

  /* ----------------------------------------------- */
  PDM_MPI_File ia_file;
  PDM_MPI_File_open(dcsr->comm,
                    name_with_extension_ia,
                    PDM_MPI_MODE_WRONLY_CREATE,
                    &ia_file);


  int sizeWriteIA = dcsr->shiftG[iRank+1] - dcsr->shiftG[iRank];
  if( iRank == nRank-1){
    sizeWriteIA += 1;
  }
  // int nOctetIA    = sizeWriteIA * sizeof(PDM__PDM_MPI_G_NUM);
  PDM_MPI_Offset offsetIA = (dcsr->shiftG[iRank] - 1 ) * sizeof(PDM__PDM_MPI_G_NUM); // = dCellProcIdx[iRank]*sizeof(E_Int);
  // printf(" sizeWriteIA = %i \n", sizeWriteIA);
  // printf(" nOctetIA    = %i \n", nOctetIA);

  PDM_g_num_t iaG[sizeWriteIA];
  for(int i = 0; i < sizeWriteIA; i++){
    iaG[i] = dcsr->ia[i] + dcsr->dnnzIdx[iRank] - 1;
  }

  int n_octet_ecrits;
  PDM_MPI_File_write_at_all(ia_file,
                            offsetIA,
                            iaG,
                            sizeWriteIA,
                            PDM__PDM_MPI_G_NUM,
                            &n_octet_ecrits);

  PDM_MPI_File_close(&ia_file);
  /* ----------------------------------------------- */


  /* ----------------------------------------------- */
  PDM_MPI_File ja_file;
  PDM_MPI_File_open(dcsr->comm,
                    name_with_extension_ja,
                    PDM_MPI_MODE_WRONLY_CREATE,
                    &ja_file);


  int sizeWriteJA = (dcsr->dnnzIdx[iRank+1] - dcsr->dnnzIdx[iRank]);
  // int nOctetJA    = sizeWriteJA * sizeof(PDM__PDM_MPI_G_NUM);
  PDM_MPI_Offset offsetJA = (dcsr->dnnzIdx[iRank] - 1 ) * sizeof(PDM__PDM_MPI_G_NUM);
  // printf(" offsetJA = %i \n",  (dcsr->dnnzIdx[iRank] - 1 ));
  // printf(" sizeWriteIA = %i \n", sizeWriteIA);

  PDM_MPI_File_write_at_all(ja_file,
                            offsetJA,
                            dcsr->ja,
                            sizeWriteJA,
                            PDM__PDM_MPI_G_NUM,
                            &n_octet_ecrits);

  PDM_MPI_File_close(&ja_file);
  /* ----------------------------------------------- */


  printf("PDM_dist_csr_dump end \n");
  fflush(stdout);
  PDM_MPI_Barrier(dcsr->comm);

  free(name_with_extension_ia);
  free(name_with_extension_ja);
}

/**
 *
 * \brief
 *
 */
void
PDM_dist_csr_print
(
 _dist_csr* dcsr
)
{
  int iRank, nRank;
  PDM_MPI_Comm_rank(dcsr->comm, &iRank);
  PDM_MPI_Comm_size(dcsr->comm, &nRank);

  PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);

  for(int jj = 0; jj < nRank; jj++){
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
    if(jj == iRank){
      printf("PDM_dist_csr_print \n");
      printf(" lSize = %i \n", dcsr->lSize);
      printf(" gSize = %i \n", dcsr->gSize);
      for(int i = 0; i < dcsr->lSize; i++){
        printf(" [%i] -> ", i);
        for(int j = dcsr->ia[i]; j < dcsr->ia[i+1]; j++){
          printf(PDM_FMT_G_NUM" ", dcsr->ja[j]);
        }
        printf("\n");
      }

      PDM_printf("dnnzIdx : "PDM_FMT_G_NUM,  dcsr->dnnzIdx[0]);
      for (int i = 1; i < nRank+1; i++) {
        PDM_printf(" "PDM_FMT_G_NUM, dcsr->dnnzIdx[i]);
      }
      PDM_printf("\n");
      printf("PDM_dist_csr_print end \n");
    }
  }
}


/**
 *
 * \brief  Usage
 *
 */
int
PDM_generate_coupling_data
(
 _dist_csr*         dcsr,
 int                nCoupling,
 PDM_g_num_t*       cpl_LNToGN,
 int**              cpl_strid,
 PDM_g_num_t**      cpl_data
)
{
  /* MPI Stuff */
  int iRank, nRank;
  PDM_MPI_Comm_rank(dcsr->comm, &iRank);
  PDM_MPI_Comm_size(dcsr->comm, &nRank);

  *cpl_strid = (int * ) malloc( nCoupling * sizeof(int) );

  for(int i = 0; i < nCoupling; i++){
    (*cpl_strid)[i] = 0;
  }

  /* Count the peripheric size */
  // int nCoupling = 0;
  for(int icol1 = 0; icol1 < dcsr->lSize; icol1++){
    int beg1 = dcsr->ia[icol1  ];
    int end1 = dcsr->ia[icol1+1];

    for(int icol2 = beg1; icol2 < end1; icol2++){
      int iCellG = dcsr->ja[icol2];
      int oppRank = PDM_binary_search_gap_long(iCellG, dcsr->shiftG, nRank+1);
      if(oppRank != iRank){
        int idx = PDM_binary_search_long(iCellG, cpl_LNToGN, nCoupling);
        (*cpl_strid)[idx] += (end1 - beg1);
      }
    }
  }

  int cpl_strid_idx[nCoupling+1];
  cpl_strid_idx[0] = 0;
  for(int i = 0; i < nCoupling; i++){
    cpl_strid_idx[i+1] = cpl_strid_idx[i] + (*cpl_strid)[i];
    (*cpl_strid)[i] = 0;
  }

  printf(" cpl_strid_idx[nCoupling] = %i \n", cpl_strid_idx[nCoupling]);

  *cpl_data = (PDM_g_num_t * ) malloc( cpl_strid_idx[nCoupling] * sizeof(PDM_g_num_t) );
  /* Count the peripheric size */
  for(int icol1 = 0; icol1 < dcsr->lSize; icol1++){
    int beg1 = dcsr->ia[icol1  ];
    int end1 = dcsr->ia[icol1+1];

    for(int icol2 = beg1; icol2 < end1; icol2++){
      int iCellG = dcsr->ja[icol2];
      int oppRank = PDM_binary_search_gap_long(iCellG, dcsr->shiftG, nRank+1);
      if(oppRank != iRank){
        int idx = PDM_binary_search_long(iCellG, cpl_LNToGN, nCoupling);
        // int idx_data = cpl_strid_idx[idx] + (*cpl_strid)[idx];
        // (*cpl_data)[idx_data] = icol1+dcsr->shiftG[iRank];
        for(int icol3 = beg1; icol3 < end1; icol3++){
          int idx_data = cpl_strid_idx[idx] + (*cpl_strid)[idx];
          (*cpl_data)[idx_data] = dcsr->ja[icol3];
          (*cpl_strid)[idx]++;
        }
      }
    }
  }

  return cpl_strid_idx[nCoupling];
}
/**
 *
 * \brief  Usage
 *
 */
int
PDM_generate_part_LNToGN
(
 _dist_csr*         dcsr,
 PDM_g_num_t**      LNToGN,
 int                withInterior
)
{
  /* MPI Stuff */
  int iRank, nRank;
  PDM_MPI_Comm_rank(dcsr->comm, &iRank);
  PDM_MPI_Comm_size(dcsr->comm, &nRank);

  int nInterior = 0;
  if(withInterior == 1){
    nInterior = dcsr->lSize;
  }

  /* Count the peripheric size */
  int nCoupling = 0;
  for(int icol1 = 0; icol1 < dcsr->lSize; icol1++){
    int beg1 = dcsr->ia[icol1  ];
    int end1 = dcsr->ia[icol1+1];

    for(int icol2 = beg1; icol2 < end1; icol2++){
      int iCellG = dcsr->ja[icol2];
      int oppRank = PDM_binary_search_gap_long(iCellG, dcsr->shiftG, nRank+1);
      if(oppRank != iRank){
        nCoupling++;
      }
    }
  }

  /* Tentative exchange MPI Color */
  int nElmtWithGhost = nInterior+nCoupling;
  *LNToGN = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * nElmtWithGhost );
  for(int i = 0; i < nInterior; i++){
    PDM_g_num_t iG = i;
    (*LNToGN)[i] = iG + dcsr->shiftG[iRank];
  }

  /* Fill border */
  int nextCell = nInterior;
  for(int icol1 = 0; icol1 < dcsr->lSize; icol1++){

    int beg1 = dcsr->ia[icol1  ];
    int end1 = dcsr->ia[icol1+1];

    for(int icol2 = beg1; icol2 < end1; icol2++){
      int iCellG = dcsr->ja[icol2];
      int oppRank = PDM_binary_search_gap_long(iCellG, dcsr->shiftG, nRank+1);
      if(oppRank != iRank){
        (*LNToGN)[nextCell++] = iCellG;
      }
    }
  }

  // Quick sort and unique
  printf("nCoupling = %i \n", nCoupling);
  int nCouplingUnique = 0;
  if( nCoupling > 0 ){
    PDM_sort_long(&(*LNToGN)[nInterior], NULL, nCoupling);

    PDM_g_num_t lastVal = -1;
    int idxWrite = nInterior;
    for(int i = nInterior; i < nElmtWithGhost; i++ ){
      if( lastVal != (*LNToGN)[i]){
        lastVal = (*LNToGN)[i];
        (*LNToGN)[idxWrite] = (*LNToGN)[i];
        idxWrite++;
        nCouplingUnique++;
      }
    }
  }
  nElmtWithGhost = nInterior+nCouplingUnique;

  /* Verbose */
  fflush(stdout);
  fflush(stdout);
  fflush(stdout);
  PDM_MPI_Barrier(dcsr->comm);
  if(1 == 1){
    PDM_MPI_Barrier(dcsr->comm);
    for(int jj = 0; jj < nRank; jj++){
      fflush(stdout);
      PDM_MPI_Barrier(dcsr->comm);
      if(jj == iRank){
        printf("[%d] LNToGN : ", iRank);
        for(int i = 0; i < nInterior+nCouplingUnique; i++){
          printf(PDM_FMT_G_NUM" ", (*LNToGN)[i]);
        }
        printf("\n");
      }
    }
    PDM_MPI_Barrier(dcsr->comm);
  }
  return nElmtWithGhost;
}



/**
 *
 * \brief
 *
 */
_dist_csr*
PDM_dist_csr_gen_superior_rank
(
 _dist_csr* dcsr
)
{
  _dist_csr* dcsrNext = PDM_dist_csr_create(dcsr->comm);
  printf(" PDM_dist_csr_gen_superior_rank \n");

  /* MPI Stuff */
  int iRank, nRank;
  PDM_MPI_Comm_rank(dcsr->comm, &iRank);
  PDM_MPI_Comm_size(dcsr->comm, &nRank);

  PDM_g_num_t* LNToGN;
  int nElmtWithGhost = PDM_generate_part_LNToGN(dcsr, &LNToGN, 1);
  int nCoupling = nElmtWithGhost - dcsr->lSize;
  int*         cpl_strid = NULL;
  PDM_g_num_t* cpl_data  = NULL;
  PDM_g_num_t* pt_cpl_LNToGN = NULL;
  if( nCoupling > 0){
    pt_cpl_LNToGN = &LNToGN[dcsr->lSize];
  }
  int nCouplingRankSup = PDM_generate_coupling_data(dcsr, nCoupling, pt_cpl_LNToGN, &cpl_strid, &cpl_data);

  dcsrNext->lSize  = dcsr->lSize;
  dcsrNext->gSize  = dcsr->gSize;
  dcsrNext->shiftG = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (nRank+1) );
  for(int i = 0; i < nRank+1; i++){
    dcsrNext->shiftG[i] = dcsr->shiftG[i];
  }

  /* ------------------------------------------------------------------------------- */
  /* Generate the superior rank interior */
  int*         iars = (int *         ) malloc( (dcsr->lSize+1) * sizeof(int));
  int*         flag = (int *         ) malloc( dcsr->lSize     * sizeof(int));
  int*  saveIdxflag = (int *         ) malloc( dcsr->lSize     * sizeof(int)); // Pb possible for small test case

  for(int i = 0; i < dcsr->lSize; i++){
    flag[i] = -1;
  }

  int nConnect = 0;
  iars[0] = 0;
  for(int icol1 = 0; icol1 < dcsr->lSize; icol1++){
    int idxSave = 0;
    for(int icol2 = dcsr->ia[icol1  ]; icol2 < dcsr->ia[icol1+1]; icol2++){
      PDM_g_num_t iCellG = dcsr->ja[icol2];
      int oppRank  = PDM_binary_search_gap_long(iCellG, dcsr->shiftG, nRank+1);
      if(oppRank == iRank){
        int iColLoc1 = dcsr->ja[icol2]     - dcsr->shiftG[iRank];
        int iColLoc2 = dcsr->ja[icol2] + 1 - dcsr->shiftG[iRank];
        int beg2     = dcsr->ia[iColLoc1];
        int end2     = dcsr->ia[iColLoc2];

        /* Panic Verbose */
        // printf("[%i] beg2 : %i | end2 = %i \n", iRank, beg2, end2);
        // printf("[%i] idx1 : %i | idx2 = %i \n", iRank, iColLoc1, iColLoc2);
        for(int icol3 = beg2; icol3 < end2; icol3++){
          PDM_g_num_t iCellG2 = dcsr->ja[icol3];
          int oppRank2 = PDM_binary_search_gap_long(iCellG2, dcsr->shiftG, nRank+1);
          if(oppRank2 == iRank){
            int iCellLoc = iCellG2 - dcsr->shiftG[iRank];
            if(flag[iCellLoc] == -1){
              flag[iCellLoc] = 1;
              saveIdxflag[idxSave++] = iCellLoc;
              nConnect++;
            }
            // int iColLoc3 = dcsr->ja[icol3]     - dcsr->shiftG[iRank];
            // int iColLoc4 = dcsr->ja[icol3] + 1 - dcsr->shiftG[iRank];
            // int beg3     = dcsr->ia[iColLoc3];
            // int end3     = dcsr->ia[iColLoc4];
          }
        }
      }
    }
    /* Reset */
    for(int ii = 0; ii < idxSave; ii++){
      flag[saveIdxflag[ii]] = -1;
    }
    iars[icol1+1] = iars[icol1] + idxSave;
  }
  /* ------------------------------------------------------------------------------- */

  /* ------------------------------------------------------------------------------- */
  printf("[%i] nConnect : %i | nConnectR1 = %i \n", iRank, nConnect, dcsr->ia[dcsr->lSize]);
  // PDM_g_num_t* jars = (PDM_g_num_t * ) malloc( dcsr->lSize * sizeof(PDM_g_num_t));
  PDM_g_num_t* jars = (PDM_g_num_t * ) malloc( nConnect*4 * sizeof(PDM_g_num_t));
  /* ------------------------------------------------------------------------------- */

  /* ------------------------------------------------------------------------------- */
  // Compute
  int idxWriteJA = 0;
  for(int icol1 = 0; icol1 < dcsr->lSize; icol1++){
    int idxSave = 0;
    int begWrite = idxWriteJA;
    for(int icol2 = dcsr->ia[icol1  ]; icol2 < dcsr->ia[icol1+1]; icol2++){
      PDM_g_num_t iCellG = dcsr->ja[icol2];
      int oppRank  = PDM_binary_search_gap_long(iCellG, dcsr->shiftG, nRank+1);
      if(oppRank == iRank){
        int iColLoc1 = dcsr->ja[icol2]     - dcsr->shiftG[iRank];
        int iColLoc2 = dcsr->ja[icol2] + 1 - dcsr->shiftG[iRank];
        int beg2     = dcsr->ia[iColLoc1];
        int end2     = dcsr->ia[iColLoc2];

        /* Panic Verbose */
        // printf("[%i] beg2 : %i | end2 = %i \n", iRank, beg2, end2);
        // printf("[%i] idx1 : %i | idx2 = %i \n", iRank, iColLoc1, iColLoc2);
        for(int icol3 = beg2; icol3 < end2; icol3++){
          PDM_g_num_t iCellG2 = dcsr->ja[icol3];
          int oppRank2 = PDM_binary_search_gap_long(iCellG2, dcsr->shiftG, nRank+1);
          if(oppRank2 == iRank){
            int iCellLoc = iCellG2 - dcsr->shiftG[iRank];
            if(flag[iCellLoc] == -1){
              flag[iCellLoc] = 1;
              saveIdxflag[idxSave++] = iCellLoc;
              jars[idxWriteJA++] = iCellG2;
            }
          }
        }
      }
    }
    /* Reset */
    for(int ii = 0; ii < idxSave; ii++){
      flag[saveIdxflag[ii]] = -1;
    }

    // Sort and normally unique is ensured
    PDM_sort_long(&jars[begWrite], NULL, idxSave);

  }
  free(flag);
  free(saveIdxflag);
  /* ------------------------------------------------------------------------------- */

  /* ------------------------------------------------------------------------------- */
  if(1 == 1){
    for(int jj = 0; jj < nRank; jj++){
      PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
      if(jj == iRank){
        printf(" DEBUG  \n");
        printf(" lSize = %i \n", dcsr->lSize);
        printf(" gSize = %i \n", dcsr->gSize);
        for(int i = 0; i < dcsr->lSize; i++){
          printf(" [%i] -> ", i);
          for(int j = iars[i]; j < iars[i+1]; j++){
            printf(PDM_FMT_G_NUM" ", jars[j]);
          }
          printf("\n");
        }
      }
    }
  }
  /* ------------------------------------------------------------------------------- */


  /* Put coupling data after main data */
  // int nSizeTot = dcsr->ia[dcsr->lSize] + nCouplingRankSup;
  int nSizeTot = iars[dcsr->lSize] + nCouplingRankSup;
  int*         part_strid = (int         *) malloc( (dcsr->lSize + nCoupling) * sizeof(int        ) );
  PDM_g_num_t* part_data  = (PDM_g_num_t *) malloc(  nSizeTot                 * sizeof(PDM_g_num_t) );
  for(int i = 0; i < dcsr->lSize; i++){
    // part_strid[i] = dcsr->ia[i+1] - dcsr->ia[i];
    part_strid[i] = iars[i+1] - iars[i];
  }
  for(int i = 0; i < nCoupling; i++){
    part_strid[i+dcsr->lSize] = cpl_strid[i];
  }
  // for(int i = 0; i < dcsr->ia[dcsr->lSize]; i++){
  //   part_data[i] = dcsr->ja[i];
  // }
  for(int i = 0; i < iars[dcsr->lSize]; i++){
    part_data[i] =jars[i];
  }
  if(nCoupling > 0){
    int offSet = iars[dcsr->lSize];
    for(int i = 0; i < nCouplingRankSup; i++){
      part_data[i+offSet] = cpl_data[i];
    }
  }

  free(iars);
  free(jars);

  /* Now we have LNToGN for the block_to_part in order to update color */
  PDM_g_num_t dCellCell[nRank+1];
  for(int i = 0; i < nRank+1; i++){
    dCellCell[i] = dcsr->shiftG[i]-1;
  }

  // PDM_g_num_t nCoupling = nElmtWithGhost - dcsr->lSize;
  // PDM_block_to_part_t *btp = PDM_block_to_part_create( dCellCell,
  //                                                     &LNToGN,
  //                                                     &nElmtWithGhost,
  //                                                      1,
  //                                                      dcsr->comm);

  PDM_part_to_block_t* ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       0.,
                                                       &LNToGN,
                                                       dCellCell,
                                                       &nElmtWithGhost,
                                                       1,
                                                       dcsr->comm);
  // Second block_to_part to exchange stencil rank 1
  printf(" Second block_to_part to exchange stencil rank 1 \n ");

  int*         blk_strid;
  PDM_g_num_t* blk_data;
  int nMaxSize = PDM_part_to_block_exch(         ptb,
                                                 sizeof(PDM_g_num_t),
                                                 PDM_STRIDE_VAR,
                                                 -1,
                                                 &part_strid,
                                        (void**) &part_data,
                                                 &blk_strid,
                                        (void**) &blk_data);


  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  /*
   * Verbose
   */
  if(1 == 1){
    fflush(stdout);
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
    for(int jj = 0; jj < nRank; jj++){
      printf("BlkSize : %d\n", blk_size);
      printf("nMaxSize : %d\n", nMaxSize);
      PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
      fflush(stdout);
      if(jj == iRank){
        int idx = 0;
        for(int i = 0; i < blk_size; i++) {
          printf("BlkStri[%d] : %d | Data -> ", i, blk_strid[i]);
          for(int ii = 0; ii < blk_strid[i]; ii++){
            printf("%i ", (int)blk_data[ii+idx]);
            // printf("%i ", jj+idx);
          }
          printf("\n");
          idx += blk_strid[i];
        }
      }
    }
    fflush(stdout);
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
  }

  dcsrNext->ia = (int * ) malloc( ( blk_size + 1 ) * sizeof(int));

  int nMaxSizeT = 0;
  for(int i = 0; i < blk_size; i++){
    nMaxSizeT += blk_strid[i];
  }
  printf("nMaxSizeT : %d\n", nMaxSizeT);
  dcsrNext->ja = (PDM_g_num_t * ) malloc( ( nMaxSizeT ) * sizeof(PDM_g_num_t) );

  // Setup and sort
  for(int i = 0; i < blk_size; i++){
    dcsrNext->ia[i] = 0;
  }

  assert(blk_size == dcsrNext->lSize);
  int idxRead = 0;
  for(int i = 0; i < blk_size; i++){
    int nDataRead = blk_strid[i];

    /* Sort */
    PDM_sort_long(&blk_data[idxRead], NULL, nDataRead);

    /* Unique and copy in the proper array */
    PDM_g_num_t lastVal = -1;
    int idxWrite = dcsrNext->ia[i];
    int nDataWrite = 0;
    for(int jj = 0; jj < nDataRead; jj++){
      // printf("[%i] --- blk_data[%i] = %i \n", iRank, idxRead+jj,blk_data[idxRead+jj]);
      if(lastVal != blk_data[idxRead+jj]){
        lastVal = blk_data[idxRead+jj];
        dcsrNext->ja[idxWrite++] = lastVal;
        nDataWrite++;
      }
    }
    /* Update ia */
    dcsrNext->ia[i+1] = dcsrNext->ia[i] + nDataWrite;

    idxRead += nDataRead;
  }

  dcsrNext->ja = realloc( dcsrNext->ja, dcsrNext->ia[dcsrNext->lSize] * sizeof(PDM_g_num_t));

  dcsrNext->dnnzIdx = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (nRank+1) );

  PDM_MPI_Allgather((void *) &dcsrNext->ia[dcsrNext->lSize],
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) &dcsrNext->dnnzIdx[1],
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dcsrNext->comm);

  dcsrNext->dnnzIdx[0] = 1;
  for (int i = 1; i < nRank+1; i++) {
    dcsrNext->dnnzIdx[i] +=  dcsrNext->dnnzIdx[i-1];
  }

  // block_to_part works with stride not index ...
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

  /*
   * Free
   */
  // PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  free(blk_strid);
  free(blk_data);
  free(part_strid);
  free(part_data);
  if(cpl_strid != NULL){
    free(cpl_strid);
  }
  if(cpl_data != NULL){
    free(cpl_data);
  }
  free(LNToGN);
  printf(" PDM_dist_csr_gen_superior_rank end \n");
  return dcsrNext;
}

/**
 *
 * \brief
 *
 */
void
PDM_dist_csr_free
(
 _dist_csr* dcsr
)
{
  if(dcsr->ia != NULL){
    free(dcsr->ia);
  }
  if(dcsr->ja != NULL){
    free(dcsr->ja);
  }
  if(dcsr->shiftG != NULL){
    free(dcsr->shiftG);
  }
  if(dcsr->dnnzIdx != NULL){
    free(dcsr->dnnzIdx);
  }
  free(dcsr);
}


/*============================================================================
 * Private function definitions
 *============================================================================*/


#ifdef  __cplusplus
}
#endif
