#ifndef __PDM_PART_BOUND_H__
#define __PDM_PART_BOUND_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum PDM_part_bound_cplx_t
 * \brief 2 complexities
 *
 */

typedef enum {

  PDM_PART_BOUND_SIMPLE  = 0,  /*!< Simple */
  PDM_PART_BOUND_CPLX = 1,  /*!< Complex */

} PDM_part_bound_cplx_t;


/**
 * \struct PDM_part_bound_t
 * \brief Partition boundary
 *
 *  PDM_part_bound_t define a partition boundary structure
 *
 */

typedef struct _part_bound_t PDM_part_bound_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Return an initialized \ref _part_bound_t structure
 *
 * This function returns an initialized \ref _part_bound_t structure
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
 * \param [in]  gNLocalOfferElt Global number of local offered elements
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
const PDM_g_num_t            *gNLocalOfferElt
);


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
);


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
);


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
);


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
);


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
 );


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
 );


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
 );


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
 * \param [in]  procPartElt          Connected vertex in the connected partition
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
 const int         procPartElt
);


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
);


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
);


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
);


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
 );

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
 );


/**
 * \brief Set distant connected vertex
 *
 * This function set local connected vertex
 *
 * \param [in]  part_bound      Inter partition boundary to free
 * \param [in]  boundElt             Element in part_bound
 * \param [in]  iConnectedElt         Connected element
 * \param [out] iProc                 Connected processus
 * \param [out] iProcPart             Connected partition in the connected processus
 * \param [out] procPartElt          Connected vertex in the connected partition
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
       int        *procPartElt,
       int        *iDistElt

);


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
);


/**
 * \brief Free a \ref PDM_part_bound_t object
 *
 * This function frees an \ref PDM_part_bound_t object
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
);


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
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_BOUND_H__ */
