#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);

  int myRank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  int n_elt_proc = 10;

  PDM_g_num_t *numabs = malloc(sizeof(PDM_g_num_t) * n_elt_proc);
  double *weights = malloc(sizeof(double) * n_elt_proc);
  int *stride = malloc(sizeof(stride) * n_elt_proc);

  for (int i = 0; i < n_elt_proc; i++) {
    numabs[i] = myRank * n_elt_proc + i + 1;
    weights[i] = myRank+1;
    stride[i] = 1;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &numabs,
                                                       &weights,
                                                       &n_elt_proc,
                                                       1,
                                                       PDM_MPI_COMM_WORLD);


  PDM_g_num_t *distrib_index = PDM_part_to_block_distrib_index_get (ptb);

  double *block_weights;
  int *block_stride;

  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &stride,
                          (void **) &weights,
                           &block_stride,
                          (void **) &block_weights);

  int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);
  double weight_sum = 0.;

  for (int i = 0; i < n_elt_block; i++) {
    weight_sum += block_weights[i];
  }

  double *weights_sum_procs = malloc(sizeof(double) * numProcs);

  printf("distrib_index : ");
  for (int i = 0; i < numProcs + 1; i++) {
    printf(PDM_FMT_G_NUM" ", distrib_index[i]);
  }
  printf("\n");

  PDM_MPI_Allgather (&weight_sum, 1, PDM_MPI_DOUBLE,
                     weights_sum_procs, 1, PDM_MPI_DOUBLE,
                     PDM_MPI_COMM_WORLD);

  printf("weights procs :");
  for (int i = 0; i < numProcs; i++) {
    printf(" %12.5e", weights_sum_procs[i]);
  }
  printf("\n");

  PDM_part_to_block_free (ptb);

  free (numabs);
  free (weights);
  free (weights_sum_procs);
  free (stride);
  free (block_stride);

  PDM_MPI_Finalize ();

  PDM_printf ("\nfin Test\n");

  return 0;

}
