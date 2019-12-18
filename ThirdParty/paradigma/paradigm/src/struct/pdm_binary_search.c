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
#include "pdm_binary_search.h"

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
 *
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_long
(
 const PDM_g_num_t   elt,
 const PDM_g_num_t  *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_long
(
 const PDM_g_num_t   elt,
 const PDM_g_num_t  *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if (elt == array[ind]) {
    return ind;
  }
  if (elt == array[right]) {
    return right;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_int
(
 const int   elt,
 const int  *array,
 const int   lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_int
(
 const int          elt,
 const int         *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if (elt == array[ind]) {
    return ind;
  }
  if (elt == array[right]) {
    return right;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_double
(
 const double   elt,
 const double  *array,
 const int      lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}



#ifdef __cplusplus
}
#endif /* __cplusplus */


