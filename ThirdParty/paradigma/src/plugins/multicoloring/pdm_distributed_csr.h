#ifndef __PDM_DISTRIBUTED_CSR_H__
#define __PDM_DISTRIBUTED_CSR_H__

/*----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"


#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef struct _dist_csr {
  PDM_MPI_Comm      comm;               /*!< Communicator                  */
  int lSize;                            /*!< Local size of current graph   */
  int gSize;                            /*!< Global size of current graph  */
  // PDM_g_num_t* ia;
  int*         ia;
  PDM_g_num_t* ja;
  PDM_g_num_t* shiftG;
  PDM_g_num_t* dnnzIdx;
} _dist_csr;

// Distributed DomainDecomposition Compressed Sparse Row
typedef struct _dist_ddcsr {
  PDM_MPI_Comm      comm;               /*!< Communicator                  */
  int lSize;                            /*!< Local size of current graph   */
  int gSize;                            /*!< Global size of current graph  */
  // PDM_g_num_t* ia;
  int*         ddb_ia;
  int*         ddb_ja;
  int*         ia;
  PDM_g_num_t* ja;
  PDM_g_num_t* shiftG;
  PDM_g_num_t* dnnzIdx;
} _dist_ddcsr;

/*============================================================================
 * Public function definitions
 *============================================================================*/

_dist_csr*
PDM_dist_csr_create
(
 const PDM_MPI_Comm           comm
);


_dist_csr*
PDM_dist_csr_gen_superior_rank
(
 _dist_csr* dcsrlow
);


void
PDM_dist_csr_print
(
 _dist_csr* dcsr
);

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
 );

void
PDM_dist_csr_dump
(
 _dist_csr* dcsr,
 char*      filename
);

void
PDM_dist_csr_free
(
 _dist_csr* dcsr
);

int
PDM_generate_part_LNToGN
(
 _dist_csr*         dcsr,
 PDM_g_num_t**      LNToGN,
 int                withInterior
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTICOLORING_H__ */
