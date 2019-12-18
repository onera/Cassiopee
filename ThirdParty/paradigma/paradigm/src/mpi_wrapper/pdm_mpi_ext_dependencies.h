/*
 * File:   pdm_mpi_ext_dependencies.h
 * Author: equemera
 *
 * Created on November 16, 2016, 1:49 PM
 */

#ifndef PDM_MPI_EXT_DEPENDENCIES_H
#define	PDM_MPI_EXT_DEPENDENCIES_H

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef PDM_HAVE_PARMETIS

int
PDM_ParMETIS_V3_PartKway
(
const PDM_g_num_t *vtxdist,
const PDM_g_num_t *xadj,
const PDM_g_num_t *adjncy,
const int *vwgt,
const int *adjwgt,
const int *wgtflag,
const int *numflag,
const int *ncon,
const int *nparts,
const double *tpwgts,
const double *ubvec,
const int *edgecut,
      int *part,
const PDM_MPI_Comm comm
);

#endif

#ifdef PDM_HAVE_PTSCOTCH

void
PDM_SCOTCH_dpart
(
const PDM_g_num_t dNCell,
const PDM_g_num_t *dDualGraphIdx,
const PDM_g_num_t *dDualGraph,
const int *cellWeight,
const int *edgeWeight,
const int check,
const PDM_MPI_Comm comm,
const int nPart,
int *part
);

#endif

void 
PDM_kaffpaE
(
int* n, 
int* vwgt,
int* xadj, 
int* adjcwgt, 
int* adjncy,
int* nparts, 
double* inbalance,  
int time_limit, 
int seed, 
int mode, 
PDM_MPI_Comm communicator, 
int* edgecut, 
double* balance, 
int* part
); 




#ifdef	__cplusplus
}
#endif

#endif	/* PDM_MPI_EXT_DEPENDENCIES_H */

