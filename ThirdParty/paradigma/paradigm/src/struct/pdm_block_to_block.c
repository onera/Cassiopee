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
#include "pdm_error.h"

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

  //btb->blockDistribIniN = malloc (sizeof(PDM_g_num_t) * (btb->nRank ));
  //btb->blockDistribEndN = malloc (sizeof(PDM_g_num_t) * (btb->nRank ));

  /* Attention si distribution commence a 1 ou 0 */
  //for (int i = 0; i < btb->nRank; i++) {
  //  btb->blockDistribIniN[i] = blockDistribIniIdx[i+1]-blockDistribIniIdx[i];
  //}
  //for (int i = 0; i < btb->nRank; i++) {
  //  btb->blockDistribEndN[i] = blockDistribEndIdx[i+1]-blockDistribEndIdx[i];
  //}

  /*
   * Verbose
   */

  if(0 == 1){

    PDM_printf("blockDistribIniIdx : ");
    for(int i = 0; i < btb->nRank+1; i++){
      PDM_printf("%i ", btb->blockDistribIniIdx[i]);
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
 * \brief Initialize an exchange (Variable stride is not yet available)
 *
 * \param [in]   btb          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Constant stride
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
 int                  cst_stride,
 int                 *block_stride_ini,
 void                *block_data_ini,
 int                 *block_stride_end,
 void                **block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

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

  unsigned char *sendBuffer = (unsigned char *) block_data_ini;

  int s_recvBuffer = 0;

  /*
   * Exchange Stride and build buffer properties
   */

  for (PDM_g_num_t i = _btb->blockDistribIniIdx[_btb->iRank]; i < _btb->blockDistribIniIdx[_btb->iRank+1]; i++) {

    int SendRank = PDM_binary_search_gap_long (i,
                                               _btb->blockDistribEndIdx,
                                               _btb->nRank + 1);
    n_sendBuffer[SendRank] += 1 * s_data;
  }

  for (PDM_g_num_t i = _btb->blockDistribEndIdx[_btb->iRank]; i < _btb->blockDistribEndIdx[_btb->iRank+1]; i++) {

    int RecvRank = PDM_binary_search_gap_long (i,
                                               _btb->blockDistribIniIdx,
                                               _btb->nRank + 1);
    n_recvBuffer[RecvRank] += 1 * s_data;
  }

  if (t_stride == PDM_STRIDE_VAR) {

    PDM_error(__FILE__, __LINE__, 0, "Error : PDM_STRIDE_VAR is not yet available \n");

    for(int i = 1; i < _btb->nRank; i++){
      i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1];
      i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1];
    }

    PDM_MPI_Alltoallv(block_stride_ini,
                      n_sendBuffer,
                      i_sendBuffer,
                      PDM_MPI_INT,
                      block_stride_end,
                      n_recvBuffer,
                      i_recvBuffer,
                      PDM_MPI_INT,
                      _btb->comm);

    for (int i = 0; i < _btb->nRank; i++) {
      n_sendBuffer[i] = 0;
      n_recvBuffer[i] = 0;
      i_sendBuffer[i] = 0;
      i_recvBuffer[i] = 0;
    }

    int k = 0;
    for (PDM_g_num_t i = _btb->blockDistribIniIdx[_btb->iRank]; i < _btb->blockDistribIniIdx[_btb->iRank+1]; i++) {

      int SendRank = PDM_binary_search_gap_long (i,
                                                 _btb->blockDistribEndIdx,
                                                 _btb->nRank + 1);
      n_sendBuffer[SendRank] += block_stride_ini[k++] * s_data;
    }

    k = 0;
    for (PDM_g_num_t i = _btb->blockDistribEndIdx[_btb->iRank]; i < _btb->blockDistribEndIdx[_btb->iRank+1]; i++) {

      int RecvRank = PDM_binary_search_gap_long (i,
                                                 _btb->blockDistribIniIdx,
                                                 _btb->nRank + 1);
      n_recvBuffer[RecvRank] += block_stride_end[k++] * s_data;
    }

    for(int i = 1; i < _btb->nRank; i++){
      i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1];
      i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1];
    }

  }

  else if (t_stride == PDM_STRIDE_CST) {

    for(int i = 1; i < _btb->nRank; i++){
      // i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1] * s_data * cst_stride;
      // i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1] * s_data * cst_stride;
      i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1] * cst_stride;
      i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1] * cst_stride;
    }

  }

  s_recvBuffer = i_recvBuffer[_btb->nRank-1] + n_recvBuffer[_btb->nRank-1];

  unsigned char *recvBuffer =
    (unsigned char *) malloc(sizeof(unsigned char) * s_recvBuffer);

  /*
   * Data exchange
   */

  PDM_MPI_Alltoallv(sendBuffer,
                    n_sendBuffer,
                    i_sendBuffer,
                    PDM_MPI_BYTE,
                    recvBuffer,
                    n_recvBuffer,
                    i_recvBuffer,
                    PDM_MPI_BYTE,
                    _btb->comm);

  free(n_sendBuffer);
  free(i_sendBuffer);
  free(n_recvBuffer);
  free(i_recvBuffer);

  *block_data_end = recvBuffer;

  return s_recvBuffer/( (int)s_data );

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

  free (_btb->blockDistribIniIdx);
  free (_btb->blockDistribEndIdx);

  free (_btb);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
