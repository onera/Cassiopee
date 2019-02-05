/* 
 * File:   pdm_ext_wrapper.h
 * Author: equemera
 *
 * Created on December 6, 2016, 9:50 AM
 */

#ifndef PDM_EXT_WRAPPER_H
#define	PDM_EXT_WRAPPER_H

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


#ifdef	__cplusplus
extern "C" {
#endif


#ifdef PDM_HAVE_PARMETIS

int 
PDM_METIS_PartGraphRecursive
(
int *nvtxs, 
int *ncon, 
int *xadj, 
int *adjncy, 
int *vwgt, 
int *adjwgt, 
int *nparts, 
double *tpwgts, 
double *ubvec, 
int *edgecut, 
int *part
);

int 
PDM_METIS_PartGraphKway
(
int *nvtxs, 
int *ncon, 
int *xadj, 
int *adjncy, 
int *vwgt, 
int *adjwgt, 
int *nparts, 
double *tpwgts, 
double *ubvec, 
int *edgecut, 
int *part
);

#endif

#ifdef PDM_HAVE_PTSCOTCH

void  
PDM_SCOTCH_part   
(
const int nCell,
int *dualGraphIdx,
int *dualGraph,        
int *cellWeight,
int *edgeWeight,
int check,        
const int nPart,        
int *part
);

#endif

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_EXT_WRAPPER_H */

