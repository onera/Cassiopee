/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_block_to_part.h"
#include "pdm_block_to_part_priv.h"
#include "pdm_binary_search.h"

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
 * \return   Initialized \ref PDM_block_to_part instance         
 *
 */

PDM_block_to_part_t *
PDM_block_to_part_create
(
 PDM_g_num_t     *blockDistribIdx,
 PDM_g_num_t    **gnum_elt,
 int            *n_elt,
 int             n_part,
 PDM_MPI_Comm        comm
)
{

  _pdm_block_to_part_t *btp = 
    (_pdm_block_to_part_t *) malloc (sizeof(_pdm_block_to_part_t));

  btp->comm = comm;
  PDM_MPI_Comm_size (comm, &btp->s_comm);
  PDM_MPI_Comm_rank (comm, &btp->myRank);

  /*
   * Define requested data for each process
   */
    
  btp->blockDistribIdx = malloc (sizeof(PDM_g_num_t) * (btp->s_comm + 1));
  for (int i = 0; i < btp->s_comm + 1; i++) {
    btp->blockDistribIdx[i] = blockDistribIdx[i];
  }
  
  btp->n_part = n_part;
  
  btp->requested_data_idx = malloc (sizeof(int) * (btp->s_comm + 1));
  btp->requested_data_n = malloc (sizeof(int) * btp->s_comm);
  for (int i = 0; i < btp->s_comm; i++) {
    btp->requested_data_idx[i] = 0;
    btp->requested_data_n[i] = 0;
  }

  btp->n_elt = malloc (sizeof(int) * n_part);
  btp->ind = malloc (sizeof(int *) * n_part);
 
  for (int i = 0; i < n_part; i++) {

    btp->n_elt[i] = n_elt[i];
    btp->ind[i] = malloc (sizeof(int) * n_elt[i]);
    
    PDM_g_num_t *_gnum_elt = gnum_elt[i];
    
    for (int j = 0; j < n_elt[i]; j++) {
    
      int ind = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                            blockDistribIdx,
                                            btp->s_comm + 1);
      btp->requested_data_n[ind]++;

    }
  }

  for (int i = 0; i < btp->s_comm; i++) {
    btp->requested_data_idx[i+1] = btp->requested_data_idx[i] + 
                                   btp->requested_data_n[i];
  }
  
  int *requested_data = malloc (sizeof(int) * (
                        btp->requested_data_idx[btp->s_comm - 1] + 
                        btp->requested_data_n[btp->s_comm - 1]));
  
  for (int i = 0; i < btp->s_comm; i++) {
    btp->requested_data_n[i] = 0;
  }
  
  for (int i = 0; i < n_part; i++) {

    PDM_g_num_t *_gnum_elt = gnum_elt[i];
    
    for (int j = 0; j < n_elt[i]; j++) {
    
      int ind = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                            blockDistribIdx,
                                            btp->s_comm + 1);
      int idx = btp->requested_data_idx[ind] + btp->requested_data_n[ind]++;

      btp->ind[i][j] = idx;
      
      PDM_g_num_t _requested_data = _gnum_elt[j] - 1 - blockDistribIdx[ind];
      requested_data[idx] = (int) _requested_data;

    }
  }

  btp->distributed_data_n = malloc (sizeof(int) * btp->s_comm);
  
  PDM_MPI_Alltoall (btp->requested_data_n,   1, PDM_MPI_INT, 
                btp->distributed_data_n, 1, PDM_MPI_INT, 
                comm);
  
  btp->distributed_data_idx = malloc (sizeof(int) * (btp->s_comm + 1));
  btp->distributed_data_idx[0] = 0;

  for (int i = 0; i < btp->s_comm; i++) {
    btp->distributed_data_idx[i+1] = btp->distributed_data_n[i] + 
                                     btp->distributed_data_idx[i];
  }

  btp->distributed_data = malloc (sizeof(int) * 
                                  btp->distributed_data_idx[btp->s_comm]);

  PDM_MPI_Alltoallv (requested_data,
                 btp->requested_data_n,
                 btp->requested_data_idx,
                 PDM_MPI_INT,
                 btp->distributed_data,
                 btp->distributed_data_n,
                 btp->distributed_data_idx,
                 PDM_MPI_INT,
                 comm);

  free (requested_data);
  
  return (PDM_block_to_part_t *) btp; 
  
}


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR  
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data 
 *
 */

void 
PDM_block_to_part_exch
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int                **part_stride,
 void               **part_data
)
{
  _pdm_block_to_part_t *_btp = (_pdm_block_to_part_t *) btp;
 
  unsigned char *_block_data = (unsigned char *) block_data; 
  unsigned char **_part_data = (unsigned char **) part_data;
  
  int *i_sendBuffer = (int *) malloc (sizeof(int) * _btp->s_comm);
  int *i_recvBuffer = (int *) malloc (sizeof(int) * _btp->s_comm);
  int *n_sendBuffer = (int *) malloc (sizeof(int) * _btp->s_comm);
  int *n_recvBuffer = (int *) malloc (sizeof(int) * _btp->s_comm);

  for (int i = 0; i < _btp->s_comm; i++) {
    n_sendBuffer[i] = 0;
    n_recvBuffer[i] = 0;
    i_sendBuffer[i] = 0;
    i_recvBuffer[i] = 0;
  }

  unsigned char *sendBuffer = NULL;
  unsigned char *recvBuffer = NULL;
  
  int s_sendBuffer = 0;
  int s_recvBuffer = 0;
  
  int s_comm1 = _btp->s_comm - 1;
  
  int s_distributed_data = _btp->distributed_data_idx[_btp->s_comm];

  /*
   * Exchange Stride and build buffer properties
   */

  int *recvStride = NULL;
  if (t_stride == PDM_STRIDE_VAR) {
    
    int s_sendStride = _btp->distributed_data_idx[_btp->s_comm];
    
    int s_recvStride = _btp->requested_data_idx[_btp->s_comm];

    int *sendStride = (int *) malloc (sizeof(int) * s_sendStride);
    recvStride = (int *) malloc (sizeof(int) * s_recvStride);
    
    for (int i = 0; i < s_sendStride; i++) {
      sendStride[i] = block_stride[_btp->distributed_data[i]];
    }

    PDM_MPI_Alltoallv (sendStride,
                   _btp->distributed_data_idx,
                   _btp->distributed_data_n,
                   PDM_MPI_INT, 
                   recvStride,
                   _btp->requested_data_idx,
                   _btp->requested_data_n,
                   PDM_MPI_INT, 
                   _btp->comm);

    for (int i = 0; i < _btp->n_part; i++) {

      for (int j = 0; j < _btp->n_elt[i]; j++) {

        int ielt = _btp->ind[i][j];
        part_stride[i][j] = recvStride[ielt];

      }
    }
   
    /*
     * Build buffers
     */
    
    for (int i = 0; i < _btp->s_comm; i++) {
      int iBeg = _btp->distributed_data_idx[i];
      int iEnd = _btp->distributed_data_idx[i] +
                 _btp->distributed_data_n[i];
      
      n_sendBuffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++)  {
        n_sendBuffer[i] += sendStride[k];
      }
            
      n_sendBuffer[i] *= (int) s_data;
            
      if (i > 0) {
        i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1];
      }
      else {
        i_sendBuffer[i] = 0;
      }         

      iBeg = _btp->requested_data_idx[i];
      iEnd = _btp->requested_data_idx[i] + 
             _btp->requested_data_n[i];
            
      n_recvBuffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++) {
        n_recvBuffer[i] += recvStride[k];
      }
            
      n_recvBuffer[i] *= (int) s_data;
            
      if (i > 0) {
        i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1];
      }
      else {
        i_recvBuffer[i] = 0;
      }

    }

    free(sendStride);

    s_sendBuffer = i_sendBuffer[s_comm1] + n_sendBuffer[s_comm1];
    s_recvBuffer = i_recvBuffer[s_comm1] + n_recvBuffer[s_comm1];

    sendBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_sendBuffer);
    recvBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recvBuffer);

    for (int i = 0; i < _btp->s_comm; i++) {
      n_sendBuffer[i] = 0;
    }

    int *sendStride_idx = (int *) malloc(sizeof(int) * (s_distributed_data+1));
    sendStride_idx[0] = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      sendStride_idx[i+1] = sendStride_idx[i] + sendStride[i]; 
    }
    
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int ind = sendStride_idx[_btp->distributed_data[i]] * (int) s_data;
      int s_block_unit = (sendStride_idx[_btp->distributed_data[i+1]] 
                        - sendStride_idx[_btp->distributed_data[i]]) * (int) s_data;
      unsigned char *_block_data_deb = _block_data + ind;  
      for (int k = 0; k < s_block_unit; k++) {
        sendBuffer[idx1++] = _block_data_deb[k];
      }
    }
  }
      
  else if (t_stride == PDM_STRIDE_CST) {
  
    int cst_stride = *block_stride;
    int s_block_unit = cst_stride * (int) s_data;
    
    for (int i = 0; i < _btp->s_comm; i++) {
      
      i_sendBuffer[i] = _btp->distributed_data_idx[i] * cst_stride * (int) s_data;
      i_recvBuffer[i] = _btp->requested_data_idx[i] * cst_stride * (int) s_data;

      n_sendBuffer[i] = _btp->distributed_data_n[i] * cst_stride * (int) s_data;
      n_recvBuffer[i] = _btp->requested_data_n[i] * cst_stride * (int) s_data;
 
    }

    s_sendBuffer = i_sendBuffer[s_comm1] + n_sendBuffer[s_comm1];
    s_recvBuffer = i_recvBuffer[s_comm1] + n_recvBuffer[s_comm1];

    sendBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_sendBuffer);
    recvBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recvBuffer);
    
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int ind = _btp->distributed_data[i];
      unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;  
      for (int k = 0; k < s_block_unit; k++) {
        sendBuffer[idx1++] = _block_data_deb[k];
      }
    }  
  }

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
                _btp->comm);
    
  free(sendBuffer);
  free(n_sendBuffer);
  free(i_sendBuffer);
  free(n_recvBuffer);
  free(i_recvBuffer);
  
  /*
   * Partitions filling
   */

  if (t_stride == PDM_STRIDE_VAR) {
 
    int s_recvElt = _btp->requested_data_idx[s_comm1] +
                    _btp->requested_data_n[s_comm1];

    int **part_idx = malloc (sizeof(int *) * _btp->n_part);
    int *recv_idx = malloc (sizeof(int) * (s_recvElt + 1));
    
    recv_idx[0] = 0;
    for (int i = 0; i < s_recvElt; i++) {
      recv_idx[i+1] = recv_idx[i] + recvStride[i];
    }

    for (int i = 0; i < _btp->n_part; i++) {

      part_idx[i] = malloc (sizeof(int) * (_btp->n_elt[i]+ 1));
      part_idx[i][0] = 0;

      for (int j = 0; j < _btp->n_elt[i]; j++) {
        part_idx[i][j+1] = part_idx[i][j] + part_stride[i][j]; 
      }

    }

    for (int i = 0; i < _btp->n_part; i++) {

      for (int j = 0; j < _btp->n_elt[i]; j++) {

        int idx1  = part_idx[i][j] * (int) s_data;
        int n_elt = part_stride[i][j] * (int) s_data;

        int idx2 = recv_idx[_btp->ind[i][j]] * (int) s_data;

        for (int k = 0; k < n_elt; k++) {
           _part_data[i][idx1+k] = recvBuffer[idx2+k];
        }
      }
    }

    for (int i = 0; i < _btp->n_part; i++) {
      free (part_idx[i]);
    }

    free(recv_idx);
    free(part_idx);
    free (recvStride);
  }
  
  else if (t_stride == PDM_STRIDE_CST) {

    const int cst_stride = *block_stride;
    const int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < _btp->n_part; i++) {

      for (int j = 0; j < _btp->n_elt[i]; j++) {

        int idx1  = j * s_block_unit;
        int idx2 = _btp->ind[i][j] * s_block_unit;

        for (int k = 0; k < s_block_unit; k++) {
           _part_data[i][idx1+k] = recvBuffer[idx2+k];
        }
      }
    }    
  }

  free(recvBuffer);

}


/**
 *
 * \brief Free a block to part structure                       
 *
 * \param [inout] btp  Block to part structure
 *
 * \return       NULL                             
 */

PDM_block_to_part_t *
PDM_block_to_part_free
(
 PDM_block_to_part_t *btp
)
{
  _pdm_block_to_part_t *_btp = (_pdm_block_to_part_t *) btp;
  
  for (int i = 0; i < _btp->n_part; i++) {
    free (_btp->ind[i]);
  }
  
  free (_btp->ind);
  free (_btp->n_elt);
  free (_btp->blockDistribIdx);
  free (_btp->distributed_data);
  free (_btp->distributed_data_idx);
  free (_btp->distributed_data_n);
  free (_btp->requested_data_idx);
  free (_btp->requested_data_n);
  
  free (_btp);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
