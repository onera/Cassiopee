#ifndef __PDM_GRAPH_BOUND_H__
#define __PDM_GRAPH_BOUND_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_bound.h"

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
 * \struct PDM_graph_bound_t
 * \brief Graph boundary 
 * 
 *  PDM_graph_bound_t define the inter boundary exchange graph
 *
 */

typedef struct _graph_bound_t PDM_graph_bound_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes 
 *============================================================================*/


/**
 * \brief Return an intialized \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  comm         MPI communicator
 * \param [in]  nPart        Number of partitions
 * \param [in]  partBound    partition boundaries (size : \ref nPart)
 *
 * \return      A new initialized \ref PDM_graph_bound_t structure
 *
 */

PDM_graph_bound_t *
PDM_graph_bound_create
(
const PDM_MPI_Comm        comm,
const int             nPart,
      PDM_part_bound_t **partBound
);


/**
 * \brief Return an intialized \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  graph_bound   Boundary graph
 * \param [in]  nComp         Number of composantes
 * \param [in]  tData         Data type
 * \param [inout]  field         Data field (with ghost cell space) 
 *
 */

void
PDM_graph_bound_exch_data_init
(
      PDM_graph_bound_t *graph_bound,
const int                nComp,
const PDM_data_t         tData,
      void             **field,
      void             **ghostField
);


/**
 * \brief Return an intialized \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [out]  recvField     Data field to receive
 *
 */

void
PDM_graph_bound_exch_data_wait
(
PDM_graph_bound_t *graph_bound
);


/**
 * \brief Return the number of ghost elements for each partition
 *
 * This function returns the number of ghost elements for each partition
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  part           Partition number
 *
 * \return  Number of ghost element
 *
 */

int
PDM_graph_bound_n_ghost_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int               part
);


/**
 * \brief Return the number of local elements which touch the current ghost element
 *
 * This function returns the number of local elements which touch the current ghost element
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  part           Partition number
 *
 * \return  Number of local element which touch the current ghost element
 *
 */

int
PDM_graph_bound_ghost_elt_n_touch_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                part,
 int                ghostElt
);


/**
 * \brief Return the list of local elements which touch the current ghost element
 *
 * This function returns the list of local elements which touch the current ghost element
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  part           Partition number
 *
 * \return  list of local elements which touch the current ghost element
 *
 */

int *
PDM_graph_bound_ghost_elt_touch_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                part,
 int                ghostElt
);


/**
 * \brief free \ref PDM_graph_bound_t structure
 *
 * This function returns an initialized \ref PDM_graph_bound_t structure
 *
 * \param [in]  graph_bound   Boundary graph
 *
 * \return      NULL
 *
 */

PDM_graph_bound_t *
PDM_graph_bound_free
(
PDM_graph_bound_t *graphBound
);


/**
 * \brief Return the number of element to send to the specified processus
 *
 * This function returns the number of element to send to the specified processus
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  iProc          Processus to send
 *
 * \return  Number of element to send
 *
 */

int
PDM_graph_bound_n_send_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                iProc
);


/**
 * \brief Return the local element numbers in theirs partition to send to the specified processus
 *
 * This function returns the local element numbers in theirs partition to send to the specified processus
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  iProc          Processus to send
 *
 * \return  Local element numbers to send to this processus
 *
 */

int *
PDM_graph_bound_send_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                iProc
 );


/**
 * \brief Return the partitions of elements to send to the specified processus
 *
 * This function returns the partitions of elements to send to the specified processus
 *
 * \param [in]  graph_bound    Boundary graph
 * \param [in]  iProc          Processus to send
 *
 * \return  Paritions of elements to send to this processus
 *
 */

int *
PDM_graph_bound_send_part_elt_get
(
 PDM_graph_bound_t *graph_bound,
 int                iProc
 );


/**
 * \brief Dump the graph boundary structure
 *
 * This function dumps the graph boundary structure
 *
 * \param [in]  graph_bound    Boundary graph
 *
 */

void
PDM_graph_bound_dump
(
 PDM_graph_bound_t *graph_bound
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GRAPH_BOUND_H__ */
