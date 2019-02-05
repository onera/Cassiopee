/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_order.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
)
{
  unsigned char *oldArray = (unsigned char *) malloc (sizeArray * elt_size);
  unsigned char *_array = (unsigned char *) array;
  
  for (int i = 0; i < sizeArray; ++i) {
    for (int j = 0; j < elt_size; ++j) {
      oldArray[elt_size * i + j] = _array[elt_size * i + j];
    }
  }
  
  for (int i = 0; i < sizeArray; ++i) {
    for (int j = 0; j < elt_size; ++j) {
      _array[elt_size * i + j] = oldArray[elt_size * newToOldOrder[i] +j];
    }
  }
  
  free(oldArray);
}

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
)
{
  size_t i;
  int o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = (int) i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    PDM_order_lnum_descend_tree_s(number, stride, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    PDM_order_lnum_descend_tree_s(number, stride, 0, i, order);
  }
}


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
)
{
  size_t i_save, i1, i2, j, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      for (j = 0; j < stride; j++) {
        if (number[i1*stride + j] != number[i2*stride + j])
          break;
      }

      if (j < stride) {
        if (number[i1*stride + j] > number[i2*stride + j])
          lv_cur++;
      }

    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    for (j = 0; j < stride; j++) {
      if (number[i1*stride + j] != number[i2*stride + j])
        break;
    }

    if (j == stride) break;
    if (number[i1*stride + j] >= number[i2*stride + j]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = (int) i_save;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */


