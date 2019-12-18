/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_ext_wrapper.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef PDM_HAVE_PARMETIS
#include <metis.h>
#endif
#ifdef PDM_HAVE_PTSCOTCH
#include <scotch.h>
#endif

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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


#ifdef PDM_HAVE_PARMETIS

int
PDM_METIS_PartGraphRecursive
(
int    *nvtxs,
int    *ncon,
int    *xadj,
int    *adjncy,
int    *vwgt,
int    *adjwgt,
int    *nparts,
double *tpwgts,
double *ubvec,
int    *edgecut,
int    *part
)
{
  idx_t options[METIS_NOPTIONS]; /* Options */
  METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NUMBERING] = 0; //C numbering = 0 (Fortran = 1)
  options[METIS_OPTION_MINCONN]   = 1; //Minimize the maximum connectivity
  options[METIS_OPTION_CONTIG]    = 1; //Force contiguous partitions
  // The graph should be compressed by combining together vertices that have identical adjacency lists.
  options[METIS_OPTION_COMPRESS]  = 1;


  //METIS provide the METIS SetDefaultOptions routine to set the options to their default values.
  //After that, the application can just modify the options that is interested in modifying.
  //options[METIS_OPTION_NSEPS] = 10;
  //options[METIS_OPTION_UFACTOR] = 100;

  idx_t _nvtxs = *nvtxs;

  idx_t _ncon = *ncon;

  int nEdge = xadj[_nvtxs];

  idx_t _nparts = *nparts;

  idx_t _edgecut = (idx_t) (*edgecut);

  real_t *_tpwgts, *__tpwgts;
  real_t *_ubvec, *__ubvec;

  __tpwgts = NULL;
  __ubvec  = NULL;

  _tpwgts  = NULL;
  _ubvec   = NULL;

  if (sizeof (double) == sizeof(real_t)) {
    _tpwgts = (real_t *) tpwgts;
    _ubvec = (real_t *) ubvec;
  }

  else {

    if(tpwgts != NULL){

      __tpwgts = malloc (sizeof(real_t) * _ncon * _nparts);
      _tpwgts  = __tpwgts;

      for (int i = 0; i < _ncon * _nparts; i++) {
        __tpwgts[i] = (real_t) tpwgts[i];
      }
    } /* End if tpwgts */

    if(ubvec != NULL){
      __ubvec = malloc (sizeof(real_t) * _ncon);
      _ubvec  = __ubvec;

      for (int i = 0; i < _ncon; i++) {
        __ubvec[i] = (real_t) ubvec[i];
      }
    } /* End if ubvec */

  }

  int *_vsize = NULL;

  idx_t *__xadj, *_xadj;
  idx_t *__adjncy, *_adjncy;
  idx_t *__vwgt, *_vwgt;
  idx_t *__adjwgt, *_adjwgt;
  idx_t *__part, *_part;

  if (sizeof(int) == sizeof(idx_t)) {
    _vwgt     = (idx_t *) vwgt;
    _adjwgt   = (idx_t *) adjwgt;
    _part     = (idx_t *) part;
    _xadj     = (idx_t *) xadj;
    _adjncy   = (idx_t *) adjncy;

    __vwgt    = NULL;
    __adjwgt  = NULL;
    __part    = NULL;
    __xadj    = NULL;
    __adjncy  = NULL;
  }

  else {
    __xadj    = (idx_t *) malloc (sizeof(idx_t) * (_nvtxs + 1));
    __adjncy  = (idx_t *) malloc (sizeof(idx_t) * nEdge);
    _xadj    = __xadj;
    _adjncy  = __adjncy;

    for (int i = 0; i < _nvtxs + 1; i++) {
      __xadj[i] = xadj[i];
    }

    for (int i = 0; i < nEdge; i++) {
      __adjncy[i] =  adjncy[i];
    }

    if (vwgt != NULL) {
      __vwgt = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);
      for (int i = 0; i < _nvtxs; i++) {
        __vwgt[i] = vwgt[i];
      }
    }
    else {
      __vwgt = NULL;
    }

    if (adjwgt != NULL) {
      __adjwgt = (idx_t *) malloc (sizeof(idx_t) * nEdge);
      for (int i = 0; i < nEdge; i++) {
        __adjwgt[i] = adjwgt[i];
      }
    }
    else {
      __adjwgt = NULL;
    }

    __part = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);

    _vwgt   = __vwgt;
    _adjwgt = __adjwgt;
    _part   = __part;

  }

  int rval = (int) METIS_PartGraphRecursive (&_nvtxs,
                                             &_ncon,
                                              _xadj,
                                              _adjncy,
                                              _vwgt,
                                              _vsize,
                                              _adjwgt,
                                              &_nparts,
                                              _tpwgts,
                                              _ubvec,
                                              options,
                                              &_edgecut,
                                              _part);

    if (sizeof(int) != sizeof(idx_t)) {
    for (int i = 0; i < _nvtxs; i++) {
      part[i] = _part[i];
    }
  }

  if (__xadj != NULL) {
    free (__xadj);
  }

  if (__adjncy != NULL) {
    free (__adjncy);
  }

  if (__part != NULL) {
    free (__part);
  }

  if (__vwgt != NULL) {
    free (__vwgt);
  }

  if (__adjwgt != NULL) {
    free (__adjwgt);
  }

  if (__tpwgts != NULL) {
    free (__tpwgts);
  }

  if (__ubvec != NULL) {
    free (__ubvec);
  }

  return rval;
}


int
PDM_METIS_PartGraphKway
(
int    *nvtxs,
int    *ncon,
int    *xadj,
int    *adjncy,
int    *vwgt,
int    *adjwgt,
int    *nparts,
double *tpwgts,
double *ubvec,
int    *edgecut,
int    *part
)
{
  idx_t options[METIS_NOPTIONS]; /* Options */
  METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NUMBERING] = 0; // C numbering = 0 (Fortran = 1)
  options[METIS_OPTION_MINCONN]   = 1; // Minimize the maximum connectivity
  options[METIS_OPTION_CONTIG]    = 1; // Force contiguous partitions
  // The graph should be compressed by combining together vertices that have identical adjacency lists.
  options[METIS_OPTION_COMPRESS]  = 1;


  //METIS provide the METIS SetDefaultOptions routine to set the options to their default values.
  //After that, the application can just modify the options that is interested in modifying.
  //options[METIS_OPTION_NSEPS] = 10;
  //options[METIS_OPTION_UFACTOR] = 100;

  idx_t _nvtxs = *nvtxs;

  idx_t _ncon = *ncon;

  int nEdge = xadj[_nvtxs];

  idx_t _nparts = *nparts;

  idx_t _edgecut = (idx_t) *edgecut;


  real_t *_tpwgts, *__tpwgts;
  real_t *_ubvec, *__ubvec;

  __tpwgts = NULL;
  __ubvec  = NULL;

  _tpwgts  = NULL;
  _ubvec   = NULL;

  if (sizeof (double) == sizeof(real_t)) {
    _tpwgts = (real_t *) tpwgts;
    _ubvec = (real_t *) ubvec;
  }

  else {

    if(tpwgts != NULL){

      __tpwgts = malloc (sizeof(real_t) * _ncon * _nparts);
      _tpwgts  = __tpwgts;

      for (int i = 0; i < _ncon * _nparts; i++) {
        __tpwgts[i] = (real_t) tpwgts[i];
      }
    } /* End if tpwgts */

    if(ubvec != NULL){
      __ubvec = malloc (sizeof(real_t) * _ncon);
      _ubvec  = __ubvec;

      for (int i = 0; i < _ncon; i++) {
        __ubvec[i] = (real_t) ubvec[i];
      }
    } /* End if ubvec */

  }

  int *_vsize = NULL;

  idx_t *__xadj, *_xadj;
  idx_t *__adjncy, *_adjncy;
  idx_t *__vwgt, *_vwgt;
  idx_t *__adjwgt, *_adjwgt;
  idx_t *__part, *_part;

  if (sizeof(int) == sizeof(idx_t)) {
    _vwgt     = (idx_t *) vwgt;
    _adjwgt   = (idx_t *) adjwgt;
    _part     = (idx_t *) part;
    _xadj     = (idx_t *) xadj;
    _adjncy   = (idx_t *) adjncy;

    __vwgt    = NULL;
    __adjwgt  = NULL;
    __part    = NULL;
    __xadj    = NULL;
    __adjncy  = NULL;
  }

  else {
    __xadj    = (idx_t *) malloc (sizeof(idx_t) * (_nvtxs + 1));
    __adjncy  = (idx_t *) malloc (sizeof(idx_t) * nEdge);
    _xadj    = __xadj;
    _adjncy  = __adjncy;

    for (int i = 0; i < _nvtxs + 1; i++) {
      __xadj[i] = xadj[i];
    }

    for (int i = 0; i < nEdge; i++) {
      __adjncy[i] =  adjncy[i];
    }

    if (vwgt != NULL) {
      __vwgt = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);
      for (int i = 0; i < _nvtxs; i++) {
        __vwgt[i] = vwgt[i];
      }
    }
    else {
      __vwgt = NULL;
    }

    if (adjwgt != NULL) {
      __adjwgt = (idx_t *) malloc (sizeof(idx_t) * nEdge);
      for (int i = 0; i < nEdge; i++) {
        __adjwgt[i] = adjwgt[i];
      }
    }
    else {
      __adjwgt = NULL;
    }

    __part = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);

    _vwgt   = __vwgt;
    _adjwgt = __adjwgt;
    _part   = __part;

  }

  int rval = (int) METIS_PartGraphKway (&_nvtxs,
                                        &_ncon,
                                         _xadj,
                                         _adjncy,
                                         _vwgt,
                                         _vsize,
                                         _adjwgt,
                                        &_nparts,
                                         _tpwgts,
                                         _ubvec,
                                         options,
                                        &_edgecut,
                                         _part);

  if (__part != NULL) {
    for (int i = 0; i < _nvtxs; i++) {
      part[i] = __part[i];
    }
  }

  if (__xadj != NULL) {
    free (__xadj);
  }

  if (__adjncy != NULL) {
    free (__adjncy);
  }

  if (__part != NULL) {
    free (__part);
  }

  if (__vwgt != NULL) {
    free (__vwgt);
  }

  if (__adjwgt != NULL) {
    free (__adjwgt);
  }

  if (__tpwgts != NULL) {
    free (__tpwgts);
  }

  if (__ubvec != NULL) {
    free (__ubvec);
  }

  return rval;
}

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
)
{

  //            Define Scotch properties
  SCOTCH_Graph grafptr;
  SCOTCH_Strat straptr;   //partitioning strategy
  int ierr = 0;

  ierr = SCOTCH_graphInit (&grafptr);
  if(ierr){
    PDM_printf("PART error : Error in PT-Scotch graph initialization\n");
    exit(1);
  }

  SCOTCH_Num _baseval = 0;
  SCOTCH_Num _vertnbr = (SCOTCH_Num) nCell;
  SCOTCH_Num *_verttab, *__verttab;
  SCOTCH_Num *_vendtab, *__vendtab;
  SCOTCH_Num *_velotab, *__velotab;
  SCOTCH_Num *_vlbltab = NULL;
  SCOTCH_Num _edgenbr = (SCOTCH_Num) dualGraphIdx[nCell];
  SCOTCH_Num _edgesiz = (SCOTCH_Num) dualGraphIdx[nCell];
  SCOTCH_Num *_edgetab, *__edgetab;
  SCOTCH_Num *_edlotab, *__edlotab;
  SCOTCH_Num *_part, *__part;

  if (sizeof(int) == sizeof(SCOTCH_Num)) {
    _verttab = dualGraphIdx;
    _vendtab = dualGraphIdx + 1;
    _edgetab = dualGraph;
    _velotab = cellWeight;
    _edlotab = edgeWeight;
    _part = part;

    __verttab = NULL;
    __vendtab = NULL;
    __edgetab = NULL;
    __velotab = NULL;
    __edlotab = NULL;
    __part = NULL;
  }

  else {
    __verttab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * (_vertnbr + 1));
    __vendtab = __verttab + 1;
    __edgetab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _edgesiz);
    __part    = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _vertnbr);

    for (int i = 0; i < _vertnbr + 1; i++) {
      __verttab[i] = dualGraphIdx[i];
    }

    for (int i = 0; i < _edgesiz; i++) {
      __edgetab[i] = dualGraph[i];
    }

    __velotab = NULL;
    __edlotab = NULL;

    if (cellWeight != NULL) {
      __velotab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _vertnbr);
      for (int i = 0; i < _vertnbr; i++) {
        __velotab[i] = cellWeight[i];
      }
    }

    if (edgeWeight != NULL) {
      __edlotab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _edgesiz);
      for (int i = 0; i < _edgesiz; i++) {
        __edlotab[i] = edgeWeight[i];
      }
    }

    _velotab = __velotab;
    _edlotab = __edlotab;
    _part    = __part;
    _verttab = __verttab;
    _vendtab = __vendtab;
    _edgetab = __edgetab;

  }

  ierr = SCOTCH_graphBuild(&grafptr,
                           _baseval,
                           _vertnbr,
                           _verttab,
                           _vendtab,
                           _velotab,
                           _vlbltab,
                           _edgenbr,
                           _edgetab,
                           _edlotab);

  if (ierr) {
    PDM_error(__FILE__, __LINE__, 0, "PART error : Error in SCOTCH_graphBuild\n");
    exit(1);
  }

  if (check) {
    ierr = SCOTCH_graphCheck (&grafptr);
    if (ierr) {
      PDM_error(__FILE__, __LINE__, 0, "PART error : Error in Scotch graph check\n");
      exit(1);
    }
  }

  //Partitioning strategy
  SCOTCH_stratInit (&straptr);

  SCOTCH_Num _nPart = (SCOTCH_Num) nPart;

  ierr = SCOTCH_graphPart (&grafptr,
                           _nPart, //nombre de partitions
                           &straptr,
                           _part);

  if (ierr) {
    PDM_error(__FILE__, __LINE__, 0, "PART error : Error in SCOTCH_graphPart\n");
    exit(1);
  }

  if (sizeof(int) != sizeof(SCOTCH_Num)) {
    for (int i = 0; i < _vertnbr; i++) {
      part[i] = _part[i];
    }
  }

  SCOTCH_stratExit (&straptr);
  SCOTCH_graphExit (&grafptr);

  if (__verttab != NULL) {
    free (__verttab);
  }

  if (__edgetab != NULL) {
    free (__edgetab);
  }

  if (__velotab != NULL) {
    free (__velotab);
  }

  if (__part != NULL) {
    free (__part);
  }

  if (__edlotab != NULL) {
    free (__edlotab);
  }

}

#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */
