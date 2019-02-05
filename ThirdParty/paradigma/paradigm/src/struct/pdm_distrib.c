/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#include "pdm_distrib.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Compute distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */

void
PDM_distrib_compute
(
 const int           dnelt,
       PDM_g_num_t  *elt_distrib,
       int           offset,
 const PDM_MPI_Comm  comm
)
{
  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(comm, &myRank);
  PDM_MPI_Comm_size(comm, &nRank);

  /* Compute distribution for element */

  // PDM_g_num_t* elt_distrib = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t  _dnelt      = (PDM_g_num_t) dnelt;

  PDM_MPI_Allgather((void *) &_dnelt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&elt_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  elt_distrib[0] = 1+offset;
  for (int i = 1; i < nRank+1; i++) {
    elt_distrib[i] +=  elt_distrib[i-1];
  }

  /* Verbose */
  if (1 == 0) {
    PDM_printf("elt_distrib : "PDM_FMT_G_NUM,  elt_distrib[0]);
    for (int i = 1; i < nRank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, elt_distrib[i]);
    }
    PDM_printf("\n");
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
