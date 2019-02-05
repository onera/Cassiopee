#ifndef __PDM_ELT_PARENT_FIND_H__
#define __PDM_ELT_PARENT_FIND_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

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
 * Public function prototypes 
 *============================================================================*/

/**
 * \brief Find parent in a set of elements 
 *
 * \param [in]     elt_distrib          Distribution of elements on processes 
 * \param [in]     elt_def_idx          Element definition index 
 * \param [in]     elt_def              Element definition
 * \param [in]     n_elt_to_find        Number of elements to find
 * \param [in]     elt_to_find_def_idx  Element to find definition index 
 * \param [in]     elt_to_find_def      Element to find definition
 * \param [in]     comm                 MPI Communicator
 * \param [inout]  parent               Parent element of found element, 0 otherwise

 */

void
PDM_elt_parent_find_from_distrib
(
 const PDM_g_num_t  *elt_distrib,
 const int          *elt_def_idx,
 const PDM_g_num_t  *elt_def,
 const PDM_g_num_t  *elt_to_find_distrib,
 const int          *elt_to_find_def_idx,
 const PDM_g_num_t  *elt_to_find_def,
 const PDM_MPI_Comm  comm,     
       PDM_g_num_t  *parent
);

/**
 * \brief Find parent in a set of elements 
 *
 * \param [in]     dnelt                Number of elt on current process
 * \param [in]     elt_def_idx          Element definition index 
 * \param [in]     elt_def              Element definition
 * \param [in]     n_elt_to_find        Number of elements to find
 * \param [in]     elt_to_find_def_idx  Element to find definition index 
 * \param [in]     elt_to_find_def      Element to find definition
 * \param [in]     comm                 MPI Communicator
 * \param [inout]  parent               Parent element of found element, 0 otherwise

 */

void
PDM_elt_parent_find
(
 const int           dnelt,
 const int          *elt_def_idx,
 const PDM_g_num_t  *elt_def,
 const int           dnelt_to_find,
 const int          *elt_to_find_def_idx,
 const PDM_g_num_t  *elt_to_find_def,
 const PDM_MPI_Comm  comm,     
       PDM_g_num_t  *parent
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ELT_PARENT_FIND_H__ */
