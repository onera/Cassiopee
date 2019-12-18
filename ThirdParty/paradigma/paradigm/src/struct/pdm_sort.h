#ifndef __PDM_SORT_H__
#define __PDM_SORT_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
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
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout] array        Array to sort
 * \param [inout] order        new indice to old indice  (or NULL)
 * \param [in]    lArray       Array length
 *
 */

void
PDM_sort_long
(
 PDM_g_num_t  *array,
 int         *order,
 int          lArray
);

/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout] array        Array to sort
 * \param [inout] order        new indice to old indice (or NULL)
 * \param [in]    lArray       Array length
 *
 */

void
PDM_sort_int
(
 int         *array,
 int         *order,
 int          lArray
);

/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout] array        Array to sort
 * \param [inout] order        new indice to old indice (or NULL)
 * \param [in]    lArray       Array length
 *
 */

void
PDM_sort_double
(
 double     *array,
 int        *order,
 int         lArray
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_SORT_H__ */
