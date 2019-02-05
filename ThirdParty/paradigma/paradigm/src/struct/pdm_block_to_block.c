/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_block_to_block.h"
#include "pdm_block_to_block_priv.h"
#include "pdm_binary_search.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   blockDistribIdx Block distribution (size : \ref size of \ref comm + 1)
 *                               C numbering (blockDistribIdx[0] = 0)
 * \param [in]   gnum_elt        Element global number (size : \ref n_part)
 * \param [in]   n_elt           Local number of elements (size : \ref n_part)
 * \param [in]   n_part          Number of partition      
 * \param [in]   comm            MPI communicator         
 *
 * \return   Initialized \ref PDM_block_to_block instance         
 *
 */

PDM_block_to_block_t *
PDM_block_to_block_create
(
 PDM_g_num_t   *blockDistribIniIdx,
 PDM_g_num_t   *blockDistribEndIdx,
 PDM_MPI_Comm   comm
)
{

  _pdm_block_to_block_t *btb = 
    (_pdm_block_to_block_t *) malloc (sizeof(_pdm_block_to_block_t));

  btb->comm = comm;
  PDM_MPI_Comm_size (comm, &btb->nRank);
  PDM_MPI_Comm_rank (comm, &btb->iRank);

  /*
   * Define requested data for each process
   */
  btb->blockDistribIniIdx = malloc (sizeof(PDM_g_num_t) * (btb->nRank + 1));
  for (int i = 0; i < btb->nRank + 1; i++) {
    btb->blockDistribIniIdx[i] = blockDistribIniIdx[i];
  }
  
  btb->blockDistribEndIdx = malloc (sizeof(PDM_g_num_t) * (btb->nRank + 1));
  for (int i = 0; i < btb->nRank + 1; i++) {
    btb->blockDistribEndIdx[i] = blockDistribEndIdx[i];
  }
  
  /*
   * Setup count 
   */
  btb->blockDistribIniN = malloc (sizeof(PDM_g_num_t) * (btb->nRank ));
  btb->blockDistribEndN = malloc (sizeof(PDM_g_num_t) * (btb->nRank ));
  
  /* Attention si distribution commence a 1 ou 0 */
  for (int i = 0; i < btb->nRank; i++) {
    btb->blockDistribIniN[i] = blockDistribIniIdx[i+1]-blockDistribIniIdx[i];
  }
  for (int i = 0; i < btb->nRank; i++) {
    btb->blockDistribEndN[i] = blockDistribEndIdx[i+1]-blockDistribEndIdx[i];
  }
  
  /*
   * Verbose 
   */
  if(1 == 1){
    
    PDM_printf("blockDistribIniIdx : ");
    for(int i = 0; i < btb->nRank+1; i++){
      PDM_printf("%i ", btb->blockDistribIniIdx[i]);
    }
    PDM_printf("\n");
    
    PDM_printf("blockDistribIniN : ");
    for(int i = 0; i < btb->nRank; i++){
      PDM_printf("%i ", btb->blockDistribIniN[i]);
    }
    PDM_printf("\n");
    
    PDM_printf("blockDistribEndN : ");
    for(int i = 0; i < btb->nRank; i++){
      PDM_printf("%i ", btb->blockDistribEndN[i]);
    }
    PDM_printf("\n");
    
    PDM_printf("blockDistribEndIdx : ");
    for(int i = 0; i < btb->nRank+1; i++){
      PDM_printf("%i ", btb->blockDistribEndIdx[i]);
    }
    PDM_printf("\n");
  }
  
  return (PDM_block_to_block_t *) btb; 
  
}


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   btb          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR  
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data 
 *
 */

int 
PDM_block_to_block_exch
(
 PDM_block_to_block_t *btb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride_ini,
 void                *block_data_ini,
 int                 *block_stride_end,
 void                *block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;
 
  // unsigned char *_block_data = (unsigned char *) block_data; 
  // unsigned char **_part_data = (unsigned char **) part_data;
  
  int *i_sendBuffer = (int *) malloc (sizeof(int) * _btb->nRank);
  int *i_recvBuffer = (int *) malloc (sizeof(int) * _btb->nRank);
  int *n_sendBuffer = (int *) malloc (sizeof(int) * _btb->nRank);
  int *n_recvBuffer = (int *) malloc (sizeof(int) * _btb->nRank);

  for (int i = 0; i < _btb->nRank; i++) {
    n_sendBuffer[i] = 0;
    n_recvBuffer[i] = 0;
    i_sendBuffer[i] = 0;
    i_recvBuffer[i] = 0;
  }

  unsigned char *sendBuffer = NULL;
  unsigned char *recvBuffer = NULL;
  
  int s_sendBuffer = 0;
  int s_recvBuffer = 0;
  int s_sendBufferdebug = 0;
  int s_recvBufferdebug = 0;
  
  int nRank1 = _btb->nRank - 1;
  
//  int s_distributed_data = _btb->blockDistribIniIdx[_btb->nRank];

  /*
   * Exchange Stride and build buffer properties
   */

  // int *recvStride = NULL;
  if (t_stride == PDM_STRIDE_VAR)
  {
    exit(2);
  }    
  else if (t_stride == PDM_STRIDE_CST) {
  
    int cst_stride = *block_stride_ini;
    int s_block_unit = cst_stride * (int) s_data;
    
    for (int i = 0; i < _btb->nRank; i++) {
      
      i_sendBuffer[i] = _btb->blockDistribIniIdx[i] * cst_stride * (int) s_data;
      i_recvBuffer[i] = _btb->blockDistribEndIdx[i] * cst_stride * (int) s_data;

      n_sendBuffer[i] = _btb->blockDistribIniN[i]   * cst_stride * (int) s_data;
      n_recvBuffer[i] = _btb->blockDistribEndN[i]   * cst_stride * (int) s_data;
 
    }

    s_sendBuffer = i_sendBuffer[nRank1] + n_sendBuffer[nRank1];
    s_recvBuffer = i_recvBuffer[nRank1] + n_recvBuffer[nRank1];
    
    s_sendBufferdebug = s_sendBuffer/(cst_stride * (int) s_data);
    s_recvBufferdebug = s_recvBuffer/(cst_stride * (int) s_data);

    sendBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_sendBuffer);
    recvBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recvBuffer);
    
    // int idx1 = 0;
    // for (int i = 0; i < s_distributed_data; i++) {
    //   int ind = _btb->blockDistribIniIdx[i];
    //   printf("ind : %d \n", ind);
    //   unsigned char *_block_data_deb = block_data_ini + ind * cst_stride * (int) s_data;  
    //   for (int k = 0; k < s_block_unit; k++) {
    //     sendBuffer[idx1++] = _block_data_deb[k];
    //   }
    // }  
    
    int idx1 = 0;
    for (int i = _btb->blockDistribIniIdx[_btb->iRank]; i < _btb->blockDistribIniIdx[_btb->iRank+1]; i++) {
      int ind = i - _btb->blockDistribIniIdx[_btb->iRank];
      printf("ind : %d \n", ind);
      unsigned char *_block_data_ini = (unsigned char *) block_data_ini;
      unsigned char *_block_data_deb = _block_data_ini + ind * cst_stride * (int) s_data;
      for (int k = 0; k < s_block_unit; k++) {
        sendBuffer[idx1++] = _block_data_deb[k];
      }
    }  
    
    
  }
  
  int* tmp = (int *) sendBuffer;

  PDM_printf("tmp : %d ", s_sendBufferdebug);
  for(int i = 0; i < s_sendBufferdebug; i++){
    PDM_printf("%i ", tmp[i]);
  }
  PDM_printf("\n");

  /*
   * Data exchange
   */
  PDM_printf("PDM_MPI_Alltoallv \n ");
  PDM_MPI_Alltoallv(sendBuffer,
                    n_sendBuffer,
                    i_sendBuffer,
                    PDM_MPI_BYTE, 
                    recvBuffer,
                    n_recvBuffer,
                    i_recvBuffer,
                    PDM_MPI_BYTE, 
                    _btb->comm);
  PDM_printf("PDM_MPI_Alltoallv end \n ");
    
  free(sendBuffer);
  free(n_sendBuffer);
  free(i_sendBuffer);
  free(n_recvBuffer);
  free(i_recvBuffer);
  
  
    
  PDM_printf("s_sendBufferdebug : %d ", s_sendBufferdebug);
  for(int i = 0; i < s_sendBufferdebug; i++){
    PDM_printf("%i ", recvBuffer[i]);
  }
  PDM_printf("\n");
  
  
  /*
   * Partitions filling
   */

  if (t_stride == PDM_STRIDE_VAR) {
  }
  else if (t_stride == PDM_STRIDE_CST) {

    // const int cst_stride = *block_stride;
    // const int s_block_unit = cst_stride * (int) s_data;

  //   for (int i = 0; i < _btb->n_part; i++) {

  //     for (int j = 0; j < _btb->n_elt[i]; j++) {

  //       int idx1  = j * s_block_unit;
  //       int idx2 = _btb->ind[i][j] * s_block_unit;

  //       for (int k = 0; k < s_block_unit; k++) {
  //          _part_data[i][idx1+k] = recvBuffer[idx2+k];
  //       }
  //     }
  //   }    
  }

  free(recvBuffer);
  
  // Attention il faut diviser par s_size_t non ????
  // Attention il faut revoir l'implementation en stride variable !!! ????
  return s_recvBuffer/ (int) s_data;

}


int 
PDM_block_to_block_exch_int
(
 PDM_block_to_block_t *btb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride_ini,
 int                 *block_data_ini,
 int                 *block_stride_end,
 int                 **block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;
 
  // unsigned char *_block_data = (unsigned char *) block_data; 
  // unsigned char **_part_data = (unsigned char **) part_data;
  
  int *i_sendBuffer = (int *) malloc (sizeof(int) * _btb->nRank);
  int *i_recvBuffer = (int *) malloc (sizeof(int) * _btb->nRank);
  int *n_sendBuffer = (int *) malloc (sizeof(int) * _btb->nRank);
  int *n_recvBuffer = (int *) malloc (sizeof(int) * _btb->nRank);

  for (int i = 0; i < _btb->nRank; i++) {
    n_sendBuffer[i] = 0;
    n_recvBuffer[i] = 0;
    i_sendBuffer[i] = 0;
    i_recvBuffer[i] = 0;
  }

  int *sendBuffer = NULL;
  int *recvBuffer = NULL;
  
  int s_sendBuffer = 0;
  int s_recvBuffer = 0;
//  int s_sendBufferdebug = 0;
//  int s_recvBufferdebug = 0;
  
//  int nRank1 = _btb->nRank - 1;
  
//  int s_distributed_data = _btb->blockDistribIniIdx[_btb->nRank];

  /*
   * Exchange Stride and build buffer properties
   */

  // int *recvStride = NULL;
  if (t_stride == PDM_STRIDE_VAR)
  {
    exit(2);
  }    
  else if (t_stride == PDM_STRIDE_CST) {
  
//    int cst_stride = *block_stride_ini;
//    int s_block_unit = cst_stride * (int) s_data;
    
    // for (int i = 0; i < _btb->nRank; i++) {
      
    //   i_sendBuffer[i] = _btb->blockDistribIniIdx[i] * cst_stride;
    //   i_recvBuffer[i] = _btb->blockDistribEndIdx[i] * cst_stride;

    //   n_sendBuffer[i] = _btb->blockDistribIniN[i]   * cst_stride;
    //   n_recvBuffer[i] = _btb->blockDistribEndN[i]   * cst_stride;
 
    // }

    // s_sendBuffer = i_sendBuffer[nRank1] + n_sendBuffer[nRank1];
    // s_recvBuffer = i_recvBuffer[nRank1] + n_recvBuffer[nRank1];

    // sendBuffer = (int *) malloc(sizeof(int ) * s_sendBuffer);
    // recvBuffer = (int *) malloc(sizeof(int ) * s_recvBuffer);
    
    // for (int i = 0; i < s_sendBuffer; i++) {
    //   sendBuffer[i] = block_data_ini[i];
    // }
    
    /* First loop to identifie which is used for exchange */
    // for (int i = 0; i < _btb->blockDistribIniIdx[_btb->nRank]; i++) {
    for (int i = _btb->blockDistribIniIdx[_btb->iRank]; i < _btb->blockDistribIniIdx[_btb->iRank+1]; i++) {
      
      /* ---------------------- */
      
      int SendRank = PDM_binary_search_gap_long (i,
                                                 _btb->blockDistribEndIdx,
                                                 _btb->nRank + 1);
      // int RecvRank = PDM_binary_search_gap_long (i,
      //                                            _btb->blockDistribEndIdx,
      //                                            _btb->nRank + 1);
      
      n_sendBuffer[SendRank] += 1;
      // n_recvBuffer[RecvRank] += 1;
      // printf("idx/ind = %d/%d \n", i, ind);
       
    }    
    
    
    // for (int i = 0; i < _btb->blockDistribEndIdx[_btb->nRank]; i++) {
    for (int i = _btb->blockDistribEndIdx[_btb->iRank]; i < _btb->blockDistribEndIdx[_btb->iRank+1]; i++) {
      
      /* ---------------------- */
      int RecvRank = PDM_binary_search_gap_long (i,
                                                 _btb->blockDistribIniIdx,
                                                 _btb->nRank + 1);
      
      n_recvBuffer[RecvRank] += 1;
      // printf("idx/ind = %d/%d \n", i, ind);
       
    }  
    
    
    
    
    // if(_btb->iRank == 0)
    // {
    //   for(int idxRank = 0; idxRank < _btb->nRank; idxRank++){
    //     printf("n_sendBuffer[%d] = %d \n", idxRank, n_sendBuffer[idxRank]);
    //     printf("n_recvBuffer[%d] = %d \n", idxRank, n_recvBuffer[idxRank]);
    //   }
    // }
    // PDM_MPI_Barrier(_btb->comm);
    
    
    for(int i = 1; i < _btb->nRank; i++){
      i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1];
      i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1];
    }
    
    
    // if(_btb->iRank == 0)
    // {
    //   for(int idxRank = 0; idxRank < _btb->nRank; idxRank++){
    //     printf("i_sendBuffer[%d] = %d \n", idxRank, i_sendBuffer[idxRank]);
    //     printf("i_recvBuffer[%d] = %d \n", idxRank, i_recvBuffer[idxRank]);
    //   }
    // }
    // PDM_MPI_Barrier(_btb->comm);
    
    
    
    // for (int i = 0; i < _btb->nRank; i++) {
      
    //   if(i == _btb->iRank)
    //   {
    //     // n_sendBuffer[i] = 0;
    //     i_sendBuffer[i] = 0;
    //     i_recvBuffer[i] = 0;
    //   }
    //   else
    //   {
    //     n_sendBuffer[i] = 0;
    //     i_sendBuffer[i] = 0;
    //     i_recvBuffer[i] = 0;
        
    //   }
    // }
    
  
    // for(int idxRank = 0; idxRank < _btb->nRank; idxRank++){
    //   PDM_MPI_Barrier(_btb->comm);
    //   if(idxRank == _btb->iRank){
    //     for(int idxRank = 0; idxRank < _btb->nRank; idxRank++){
    //       // printf("OO - i_sendBuffer[%d] = %d \n", idxRank, i_sendBuffer[idxRank]);
    //       // printf("OO - i_recvBuffer[%d] = %d \n", idxRank, i_recvBuffer[idxRank]);
    //       printf("OO - n_sendBuffer[%d] = %d \n", idxRank, n_sendBuffer[idxRank]);
    //       printf("OO - n_recvBuffer[%d] = %d \n", idxRank, n_recvBuffer[idxRank]);
    //     }
    //   }
    // }
    
    PDM_MPI_Barrier(_btb->comm);
    
    
    s_sendBuffer = i_sendBuffer[_btb->nRank-1] + n_sendBuffer[_btb->nRank-1];
    s_recvBuffer = i_recvBuffer[_btb->nRank-1] + n_recvBuffer[_btb->nRank-1];
    /* Alloc field */
    // sendBuffer = (int *) malloc(sizeof(int ) * n_sendBuffer[_btb->iRank]);
    // recvBuffer = (int *) malloc(sizeof(int ) * n_recvBuffer[_btb->iRank]);
    sendBuffer = (int *) malloc(sizeof(int ) * s_sendBuffer);
    recvBuffer = (int *) malloc(sizeof(int ) * s_recvBuffer);
    
    
    /* Fill up buffer send */
    // s_sendBuffer = i_sendBuffer[_btb->iRank] + n_sendBuffer[nRank1];
    // for (int i = 0; i < n_sendBuffer[_btb->iRank]; i++) {
    for (int i = 0; i < s_sendBuffer; i++) {
      sendBuffer[i] = block_data_ini[i];
    }
    
    // for (int i = 0; i < n_recvBuffer[_btb->iRank]; i++) {
    for (int i = 0; i < s_recvBuffer; i++) {
      recvBuffer[i] = -1;
    }
    
    // int idx1 = 0;
    // for (int i = 0; i < s_distributed_data; i++) {
    //   int ind = _btb->blockDistribIniIdx[i];
    //   printf("ind : %d \n", ind);
    //   unsigned char *_block_data_deb = block_data_ini + ind * cst_stride * (int) s_data;  
    //   for (int k = 0; k < s_block_unit; k++) {
    //     sendBuffer[idx1++] = _block_data_deb[k];
    //   }
    // }      
    
  }
  
  // for(int idxRank = 0; idxRank < _btb->nRank; idxRank++){
  //   PDM_MPI_Barrier(_btb->comm);
  //   if(idxRank == _btb->iRank){
  //   PDM_printf("tmp : %d ", n_sendBuffer[_btb->iRank]);
  //   // for(int i = 0; i < n_sendBuffer[_btb->iRank]; i++){
  //   for(int i = 0; i < s_sendBuffer; i++){
  //     PDM_printf("%i ", sendBuffer[i]);
  //   }
  //   PDM_printf("\n");
  //   }
  // }

  /*
   * Data exchange
   */
  // PDM_printf("[%d/%d] - PDM_MPI_Alltoallv %d \n ", _btb->iRank, _btb->nRank, n_sendBuffer[_btb->iRank]);
  PDM_MPI_Barrier(_btb->comm);
  
  PDM_MPI_Alltoallv(sendBuffer,
                    n_sendBuffer,
                    i_sendBuffer,
                    PDM_MPI_INT, 
                    recvBuffer,
                    n_recvBuffer,
                    i_recvBuffer,
                    PDM_MPI_INT, 
                    _btb->comm);
  // PDM_printf("[%d/%d] -  PDM_MPI_Alltoallv end \n ",_btb->iRank, _btb->nRank);
    

  // PDM_printf("s_recvBuffer : %d ", n_recvBuffer[_btb->iRank]);
  
  // for(int idxRank = 0; idxRank < _btb->nRank; idxRank++){
  //   PDM_MPI_Barrier(_btb->comm);
  //   if(idxRank == _btb->iRank){
  //   PDM_printf("s_recvBuffer : %d ", s_recvBuffer);
  //   // for(int i = 0; i < n_recvBuffer[_btb->iRank]; i++){
  //   for(int i = 0; i < s_recvBuffer; i++){
  //     PDM_printf("%i ", recvBuffer[i]);
  //   }
  //   PDM_printf("\n");    
  //   PDM_MPI_Barrier(_btb->comm);
  //   }
  // }
    

  free(sendBuffer);
  free(n_sendBuffer);
  free(i_sendBuffer);
  free(n_recvBuffer);
  free(i_recvBuffer);
  
  
  
  
  /*
   * Partitions filling
   */

  if (t_stride == PDM_STRIDE_VAR) {
  }
  else if (t_stride == PDM_STRIDE_CST) {
    
    // block_data_end = (int *) malloc(sizeof(int) * s_recvBuffer);
    int* _block_data_end = (int *) malloc(sizeof(int) * s_recvBuffer);
    for(int i = 0; i < s_recvBuffer; i++){
      // *block_data_end[i] = recvBuffer[i];
      _block_data_end[i] = recvBuffer[i];
      // PDM_printf("oo %i ", _block_data_end[i]);
    }
    *block_data_end = _block_data_end;
    // PDM_printf("\n");    
  }

  free(recvBuffer);
  // PDM_printf("s_recvBuffer :  %i \n ", s_recvBuffer);
  return s_recvBuffer;
}


/**
 *
 * \brief Free a block to part structure                       
 *
 * \param [inout] btb  Block to part structure
 *
 * \return       NULL                             
 */

PDM_block_to_block_t *
PDM_block_to_block_free
(
 PDM_block_to_block_t *btb
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;
  
  // for (int i = 0; i < _btb->n_part; i++) {
  //   free (_btb->ind[i]);
  // }
  
  free (_btb->blockDistribIniIdx);
  free (_btb->blockDistribEndIdx);
  free (_btb->blockDistribIniN  );
  free (_btb->blockDistribEndN  );
  // free (_btb->ind);
  // free (_btb->n_elt);
  // free (_btb->distributed_data);
  // free (_btb->distributed_data_idx);
  // free (_btb->distributed_data_n);
  // free (_btb->requested_data_idx);
  // free (_btb->requested_data_n);
  
  free (_btb);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
