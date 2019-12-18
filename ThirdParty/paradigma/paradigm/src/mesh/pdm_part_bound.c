/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_bound.h"
#include "pdm_part_bound_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
 * \brief Return an intialized \ref PDM_part_bound_t structure
 *
 * This function returns an initialized \ref PDM_part_bound_t structure
 *
 * \param [in]  lComm           Size of MPI communicator
 * \param [in]  nElt            Number of elements
 * \param [in]  nEltPartBound   Number of elements in the structure
 * \param [in]  cplx            Complexity
 * \param [in]  nConnectedElt   Number of connected elements of each element
 *                              (size : 1 or nEltPartBound according complexity)
 * \param [in]  nOfferElt       Number of connected elements of each element
 *                              (size : 1 or nEltPartBound according complexity)
 * \param [in]  nTotalOfferElt  Total number of offered elements
 * \param [in]  nLocalOfferElt  Local number of offered eleemnts
 * \param [in]  localOfferLnToGn Global num of local offer element
 *                              (size : nLocalOfferElt)
 *
 * \return      A new initialized \ref PDM_part_bound_t structure
 *
 */

PDM_part_bound_t *
PDM_part_bound_create
(
const int                    lComm,
const int                    nElt,
const int                    nEltPartBound,
const PDM_part_bound_cplx_t  cplx,
const int                   *nConnectedElt,
const int                   *nOfferElt,
const PDM_g_num_t            nTotalOfferElt,
const int                   nLocalOfferElt,
const PDM_g_num_t            *localOfferLnToGn
)

{
  PDM_part_bound_t * part_bound =
    (PDM_part_bound_t *) malloc(sizeof(_part_bound_t));
  _part_bound_t *_part_bound = (_part_bound_t *)  part_bound;

  _part_bound->nElt = nElt;
  _part_bound->nEltPartBound = nEltPartBound;
  _part_bound->lComm = lComm;
  _part_bound->eltPartBoundIdx = (int *) malloc(sizeof(int) * (nEltPartBound + 1));
  _part_bound->eltPartBound = NULL;
  _part_bound->cplx = cplx;
  _part_bound->eltPartBoundIdx[0] = 0;
  _part_bound->nTotalOfferElt = nTotalOfferElt;
  _part_bound->nLocalOfferElt = nLocalOfferElt;
  _part_bound->localOfferLnToGn = localOfferLnToGn;

  if (cplx == PDM_PART_BOUND_SIMPLE) {
    _part_bound->nConnectedElt = (int *) malloc(sizeof(int));
    _part_bound->connectedEltIdx = (int *) malloc(sizeof(int) * (nEltPartBound + 1));
    _part_bound->nConnectedElt[0] = *nConnectedElt;
    _part_bound->eltPartBound = (int *) malloc(sizeof(int) *
                                               (nDataEltPartBoundIni +
                                                nDataEltPartBoundElt *
                                                (*nConnectedElt)) *
                                               _part_bound->nEltPartBound);
    _part_bound->connectedEltIdx[0] = 0;

    _part_bound->nOfferElt = (int *) malloc(sizeof(int));
    _part_bound->offerEltIdx = (int *) malloc(sizeof(int) * (nEltPartBound + 1));
    _part_bound->nOfferElt[0] = *nOfferElt;
    _part_bound->offerEltIdx[0] = 0;
    _part_bound->offerElt = (int *) malloc(sizeof(int) * nEltPartBound *
					   (*nOfferElt));
    _part_bound->offerLnToGn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) *
						      nEltPartBound * (*nOfferElt));

    for (int i = 0; i < nEltPartBound; i++) {
      _part_bound->connectedEltIdx[i+1] = _part_bound->connectedEltIdx[i] +
	*nConnectedElt;
      _part_bound->offerEltIdx[i+1] =  _part_bound->offerEltIdx[i] + *nOfferElt;
      _part_bound->eltPartBoundIdx[i+1] = _part_bound->eltPartBoundIdx[i] +
                                          nDataEltPartBoundIni +
                                          nDataEltPartBoundElt * (*nConnectedElt);
    }
  }
  else {
    _part_bound->nConnectedElt = (int *) malloc(sizeof(int) * nEltPartBound);
    _part_bound->connectedEltIdx = (int *) malloc(sizeof(int) * (nEltPartBound + 1));
    _part_bound->connectedEltIdx[0] = 0;
    memcpy(_part_bound->nConnectedElt, nConnectedElt, sizeof(int)*nEltPartBound);

    _part_bound->nOfferElt = (int *) malloc(sizeof(int) * nEltPartBound);
    _part_bound->offerEltIdx = (int *) malloc(sizeof(int) * (nEltPartBound + 1));
    _part_bound->offerEltIdx[0] = 0;
    memcpy(_part_bound->nOfferElt, nOfferElt, sizeof(int)*nEltPartBound);

    int tConnectedElt = 0;
    for (int i = 0; i < nEltPartBound; i++) {
      _part_bound->connectedEltIdx[i+1] = _part_bound->connectedEltIdx[i] +
	nConnectedElt[i];
      _part_bound->offerEltIdx[i+1] = _part_bound->offerEltIdx[i] + nOfferElt[i];
      tConnectedElt += nConnectedElt[i];
    }

    _part_bound->offerElt = (int *) malloc(sizeof(int) *
                                           _part_bound->offerEltIdx[nEltPartBound]);
    _part_bound->offerLnToGn =
      (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) *
			    _part_bound->offerEltIdx[nEltPartBound]);

    _part_bound->eltPartBound =
      (int *) malloc(sizeof(int) *
		     (nDataEltPartBoundIni * _part_bound->nEltPartBound +
		      nDataEltPartBoundElt * tConnectedElt));

    for (int i = 0; i < nEltPartBound; i++) {
      _part_bound->eltPartBoundIdx[i+1] = _part_bound->eltPartBoundIdx[i] +
                                           nDataEltPartBoundIni +
                                           nDataEltPartBoundElt * nConnectedElt[i];
    }
  }

  _part_bound->localElt2BoundElt = malloc (sizeof(int) * nElt);

  return part_bound;
}


/**
 * \brief Set local connected element
 *
 * This function set local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt             Element in part_bound
 * \param [in]  iElt                  Vertex number
 *
 */

void
PDM_part_bound_local_elt_set
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         localElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->eltPartBoundIdx[iBoundElt];
  _part_bound->eltPartBound[idx] = localElt;
  _part_bound->localElt2BoundElt[localElt - 1] = boundElt;
}


/**
 * \brief Return number of offered elements
 *
 * This function returns the numner of offered elements
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt        Element in part_bound
 *
 * \return Number of offered element
 *
 */

int
PDM_part_bound_n_offer_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         boundElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int nOffer;
  if (_part_bound->cplx == PDM_PART_BOUND_CPLX) {
    nOffer = _part_bound->nOfferElt[iBoundElt];
  }
  else {
    nOffer = _part_bound->nOfferElt[0];
  }
  return nOffer;
}


/**
 * \brief Return total number of offered elements
 *
 * This function returns the numner of offered elements
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return Total number of offered element
 *
 */

PDM_g_num_t
PDM_part_bound_n_total_offer_elt_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->nTotalOfferElt;
}


/**
 * \brief Return number of local offered elements
 *
 * This function returns the numner of offered elements
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return number of local offered element
 *
 */

int
PDM_part_bound_n_local_offer_elt_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->nLocalOfferElt;
}


/**
 * \brief Return global numbers of local offered elements
 *
 * This function returns the global numbers of local offered elements
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return number of local offered element
 *
 */

const PDM_g_num_t *
PDM_part_bound_local_offer_elt_ln_to_gn_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->localOfferLnToGn;
}


/**
 * \brief Set local connected element
 *
 * This function set local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt        Element in part_bound
 * \param [in]  iOfferElt       index of offered element
 * \param [in]  iNum            local number
 * \param [in]  gNum            global number
 *
 */

void
PDM_part_bound_offer_elt_set
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iOfferElt,
 const int         lNum,
 const PDM_g_num_t  gNum
 )
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->offerEltIdx[iBoundElt] + iOfferElt;
  _part_bound->offerElt[idx] = lNum;
  _part_bound->offerLnToGn[idx] = gNum;
}


/**
 * \brief Get offered element
 *
 * This function gets an offered element
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt        Element in part_bound
 * \param [in]  iOfferElt       index of offered element
 * \param [out]  iNum           local number
 * \param [out]  gNum           global number
 *
 */

void
PDM_part_bound_offer_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iOfferElt,
 int              *lNum,
 PDM_g_num_t       *gNum
 )
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->offerEltIdx[iBoundElt] + iOfferElt;
  *lNum = _part_bound->offerElt[idx];
  *gNum = _part_bound->offerLnToGn[idx];
}


/**
 * \brief Get the number of elements in partition boundary
 *
 * This function returns the number of elements in partition boundary
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return Number of elements in partition boundary
 *
 */

int
PDM_part_bound_n_elt_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->nElt;
}


/**
 * \brief Get the number of elements in partition boundary
 *
 * This function returns the number of elements in partition boundary
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return Number of elements in partition boundary
 *
 */

int
PDM_part_bound_n_elt_bound_get
(
 PDM_part_bound_t *part_bound
)
{
 _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
 return _part_bound->nEltPartBound;
}


/**
 * \brief Get the complexity of the partition boundary
 *
 * This function returns the complexity of the partition boundary
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return Complexity
 *
 */

PDM_part_bound_cplx_t
PDM_part_bound_cplx_get
(
 PDM_part_bound_t *part_bound
)
{
 _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
 return _part_bound->cplx;
}


/**
 * \brief Get local connected element
 *
 * This function returns local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt       Element in part_bound
 * \param [out] iElt            Local element number
 * \param [out] nConnectedElt   Number of connected elements to local element
 *
 */

void
PDM_part_bound_bound_elt_get
(
 PDM_part_bound_t *part_bound,
 const int      boundElt,
       int     *localElt,
       int     *nConnectedElt
 )
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->eltPartBoundIdx[iBoundElt];
  *localElt = _part_bound->eltPartBound[idx];
  if (_part_bound->cplx == PDM_PART_BOUND_CPLX) {
    *nConnectedElt = _part_bound->nConnectedElt[iBoundElt];
  }
  else {
    *nConnectedElt = _part_bound->nConnectedElt[0];
  }
}


/**
 * \brief Get local connected element
 *
 * This function returns local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  localElt        Local element number
 * \param [out] boundElt        Element in part_bound
 * \param [out] nConnectedElt   Number of connected elements to local element
 *
 */

void
PDM_part_bound_local_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         localElt,
       int        *boundElt,
       int        *nConnectedElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iLocalElt = localElt - 1;
  *boundElt = _part_bound->localElt2BoundElt[iLocalElt];

  int localElt2;
  PDM_part_bound_bound_elt_get (part_bound,
                                *boundElt,
                                &localElt2,
                                nConnectedElt);

}


/**
 * \brief Set distant connected vertex
 *
 * This function set local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt             Element in part_bound
 * \param [in]  iConnectedElt         Connected element
 * \param [in]  iProc                 Connected processus
 * \param [in]  iProcPart             Connected partition in the connected processus
 * \param [in]  iProcPartElt          Connected vertex in the connected partition
 *
 */

void
PDM_part_bound_distant_elt_set
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iConnectedElt,
 const int         iProc,
 const int         iProcPart,
 const int         iProcPartElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;

  int nConnectedElt = _part_bound->nConnectedElt[0];

  if (_part_bound->cplx == PDM_PART_BOUND_CPLX)
    nConnectedElt = _part_bound->nConnectedElt[iBoundElt];

  if (iConnectedElt >= nConnectedElt) {
    PDM_error(__FILE__, __LINE__, 0, "Error part_bound_distant_elt_set :"
	    "Error in edgeFace computing\n");
    abort();
  }

  int idx = _part_bound->eltPartBoundIdx[iBoundElt] + nDataEltPartBoundIni +
    iConnectedElt * nDataEltPartBoundElt;

  _part_bound->eltPartBound[idx++] = iProc;
  _part_bound->eltPartBound[idx++] = iProcPart;
  _part_bound->eltPartBound[idx++] = iProcPartElt;
  _part_bound->eltPartBound[idx++] = _part_bound->connectedEltIdx[iBoundElt] +
    iConnectedElt;

}


/**
 * \brief Set distant connected vertex
 *
 * This function set local connected vertex
 *
 * \param [in]  part_bound            Inter partition boundary structure
 * \param [in]  boundElt             Element in part_bound
 * \param [in]  iConnectedElt         Connected element
 * \param [out] iProc                 Connected processus
 * \param [out] iProcPart             Connected partition in the connected processus
 * \param [out] iProcPartElt          Connected vertex in the connected partition
 * \param [out] iDistElt              Global Index of distant connected element
 *
 */

void
PDM_part_bound_distant_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iConnectedElt,
       int        *iProc,
       int        *iProcPart,
       int        *iProcPartElt,
       int        *iDistElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;

  int nConnectedElt = _part_bound->nConnectedElt[0];

  if (_part_bound->cplx == PDM_PART_BOUND_CPLX)
    nConnectedElt = _part_bound->nConnectedElt[iBoundElt];

  if (iConnectedElt >= nConnectedElt) {
    PDM_error(__FILE__, __LINE__, 0, "Error part_bound_distant_elt_get :"
	    "iConnectedElt > nConnectedElt\n");
    abort();
  }

  int idx = _part_bound->eltPartBoundIdx[iBoundElt] +
    nDataEltPartBoundIni +
    iConnectedElt * nDataEltPartBoundElt;

  *iProc        = _part_bound->eltPartBound[idx++];
  *iProcPart    = _part_bound->eltPartBound[idx++];
  *iProcPartElt = _part_bound->eltPartBound[idx++];
  *iDistElt     = _part_bound->eltPartBound[idx++];
}


/**
 * \brief Set distant connected vertex
 *
 * This function set local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  nEltPartBound Number of elements in the structure
 *
 */

void
PDM_part_bound_adjust_size
(
 PDM_part_bound_t *part_bound,
 const int         nEltPartBound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  if (nEltPartBound > _part_bound->nEltPartBound) {
    PDM_error(__FILE__, __LINE__, 0, "Error _part_bound_adjust_size : Error this function"
                    "can't increase the size of part_bound structure\n");
    abort();
  }

  _part_bound->nEltPartBound = nEltPartBound;
  _part_bound->eltPartBoundIdx = (int *) realloc(_part_bound->eltPartBoundIdx,
                                                 sizeof(int) * (nEltPartBound + 1));

  if (_part_bound->cplx == PDM_PART_BOUND_CPLX)
    _part_bound->nConnectedElt = (int *) realloc(_part_bound->nConnectedElt,
                                                 sizeof(int) * nEltPartBound);

  _part_bound->eltPartBound =
    (int *) realloc(_part_bound->eltPartBound,
		    sizeof(int) *
		    _part_bound->eltPartBoundIdx[nEltPartBound]);

}


/**
 * \brief Free a part_bound_t object
 *
 * This function frees an part_bound_t object
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 * \return      NULL
 *
 */

PDM_part_bound_t *
PDM_part_bound_free
(
PDM_part_bound_t *part_bound
)
{

  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
 if (_part_bound != NULL) {
    if (_part_bound->eltPartBoundIdx != NULL)
      free(_part_bound->eltPartBoundIdx);
    if (_part_bound->eltPartBound != NULL)
      free(_part_bound->eltPartBound);
    if (_part_bound->connectedEltIdx != NULL)
      free(_part_bound->connectedEltIdx);
    if (_part_bound->nConnectedElt != NULL)
      free(_part_bound->nConnectedElt);
    if (_part_bound->localElt2BoundElt != NULL)
      free(_part_bound->localElt2BoundElt);
    if (_part_bound->nOfferElt != NULL)
      free(_part_bound->nOfferElt);
    if (_part_bound->offerEltIdx != NULL)
      free(_part_bound->offerEltIdx);
    if (_part_bound->offerElt != NULL)
      free(_part_bound->offerElt);
    if (_part_bound->offerLnToGn != NULL)
      free(_part_bound->offerLnToGn);

    free(_part_bound);
  }
  return NULL;
}


/**
 * \brief Dump a part_bound_t object
 *
 * This function dumps a part bound structure
 *
 * \param [in]  part_bound      Inter partition boundary to free
 *
 */

void
PDM_part_bound_dump
(
 PDM_part_bound_t *part_bound
 )
{
  if (part_bound != NULL) {

    PDM_part_bound_cplx_t cplx = PDM_part_bound_cplx_get (part_bound);
    int nEltPartBound = PDM_part_bound_n_elt_bound_get (part_bound);
    PDM_printf("PDM_part_bound :\n");
    PDM_printf("  - cplx %d\n", cplx);
    PDM_printf("  - nEltPartBound %d\n", nEltPartBound);
    for (int j = 0; j < nEltPartBound; j++) {

      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (part_bound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      PDM_printf("    - localElt %d\n", localElt);
      PDM_printf("      - nConnectedElt %d\n", nConnectedElt);
      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int iPart;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get (part_bound,
                                        j+1,
                                        k,
                                        &iProc,
                                        &iPart,
                                        &iElt,
                                        &iDistElt);
        PDM_printf("        - %d %d %d %d (iproc iPart, iElt, iDistElt)\n", iProc, iPart, iElt, iDistElt);

      }
      int nOfferElt = PDM_part_bound_n_offer_elt_get (part_bound, j+1);
      PDM_printf("      - nOfferElt %d\n", nOfferElt);

      for (int k = 0; k < nOfferElt; k++) {
        int        lNum;
        PDM_g_num_t gNum;
        PDM_part_bound_offer_elt_get (part_bound,
                                      j+1,
                                      k,
                                      &lNum,
                                      &gNum);
        PDM_printf("        - %d "PDM_FMT_G_NUM" (lNum gNum)\n", lNum, gNum);

      }
    }

    fflush(stdout);
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
