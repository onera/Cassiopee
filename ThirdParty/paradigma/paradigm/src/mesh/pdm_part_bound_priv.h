#ifndef __PDM_PART_BOUND_PRIV_H__
#define __PDM_PART_BOUND_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_bound.h"
#include "pdm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

static const int nDataEltPartBoundIni = 1; /*!< Number of Initial data
                                                in \ref eltPartBound */

static const int nDataEltPartBoundElt = 4; /*!< Number of data for each connected
                                                element in \ref eltPartBound */

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _part_bound_t
 * \brief  inter partition boundary
 *
 * _part_bound_t defines a mesh partition structure
 *
 */

typedef struct {

  int  lComm;             /*!< Size of MPI communicator */
  int  nEltPartBound;     /*!< Number of partitioning boundary elts */
  int  nElt;              /*!< Number of elts */
  PDM_part_bound_cplx_t cplx;   /*!< Complexity */
  int *eltPartBoundIdx; /*!< Partitioning boundary bloc distribution from processus
                          (size = \ref nEltpartBound + 1) */
  int *eltPartBound;  /*!< Partitioning boundary elts sorted by
                           proc, sorted by part in the proc, and
                           sorted by absolute face number in the part
                           For each face :
                             - Elt local number
                             - for each connected element :
                                - Connected process
                                - Connected Partition
                                  on the connected process
                                - Connected elt local number
                                  in the connected partition
                                  (size = \ref nDataEltPartBound * nEltPartBound)
                                - local ghost number */
  int *nConnectedElt;     /*!< Number of connected elements
                            (size = \ref nEltpartBound)*/
  int *connectedEltIdx;    /*!< Number of connected elements
                            (size = \ref nEltpartBound)*/
  int *localElt2BoundElt;  /*!< Indrection from localElt to element part boundary
                            (size = \ref nLocalElt)*/
  int nTotalConnectedElt;  /*!< Total number of connected elements */

  int *nOfferElt;        /*!< offer element for each connected element index
                            (size = \ref nEltPartBound + 1)*/
  PDM_g_num_t nTotalOfferElt;  /*!< Total number of offered element */

  int        nLocalOfferElt; /*!< Number of local offered element */

  int *offerEltIdx;        /*!< Offer element for each connected element index
                             (size = \ref nEltPartBound + 1)*/
  int *offerElt;           /*!< oOffer element for each connected element
                             (size = \ref offerEltIdx[nEltPartBound] */
  PDM_g_num_t *offerLnToGn; /*!< Global num of offer element for each connected element
                             (size = \ref offerEltIdx[nEltPartBound] */

  const PDM_g_num_t *localOfferLnToGn; /*!< Global num of local offer element
                                        (size = \ref nLocalOfferElt */

} _part_bound_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_BOUND_PRIV_H__ */
