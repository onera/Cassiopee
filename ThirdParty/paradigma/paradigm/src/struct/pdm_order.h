#ifndef __PDM_ORDER_H__
#define __PDM_ORDER_H__

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
 * \brief Order an array
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      newToOldOrder   New order (size = \ref nElt
 * \param [in, out] Array           Array to renumber
 *
 */

void
PDM_order_array
(
const int     sizeArray,
const size_t  elt_size,
const int    *newToOldOrder,
void         *array
);

/**
 * This function is part of Code_Saturne, a general-purpose CFD tool.
 *  Copyright (C) 1998-2014 EDF S.A.
 *
 * \brief Order a strided array of global numbers lexicographically.
 *
 * \param [in]     number array of entity numbers (if NULL, a default 1 to n numbering is considered)
 * \param [in]     stride stride of array (number of values to compare)
 * \param [in,out] order  pre-allocated ordering table
 * \param [in]     nb_ent number of entities considered
 */

void
PDM_order_lnum_s
(
const int    number[],
size_t       stride,
int          order[],
const size_t nb_ent
);

/**
 * This function is part of Code_Saturne, a general-purpose CFD tool.
 *  Copyright (C) 1998-2014 EDF S.A.
 *
 * \brief Descend binary tree for the lexicographical ordering of a strided array
 *
 * \param [in]     number pointer to numbers of entities that should be ordered.
 * \param [in]            (if NULL, a default 1 to n numbering is considered)
 * \param [in]     stride stride of array (number of values to compare)
 * \param [in]     level  level of the binary tree to descend
 * \param [in]     nb_ent number of entities in the binary tree to descend
 * \param [in,out] order  ordering array
 */
void
PDM_order_lnum_descend_tree_s
(
const int    number[],
size_t       stride,
size_t       level,
const size_t nb_ent,
int          order[]
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_BINARY_SEARCH_H__ */
