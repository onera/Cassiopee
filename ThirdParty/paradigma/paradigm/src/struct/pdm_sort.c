/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_sort.h"

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

/**
 *
 * \brief Swap element values
 *
 * \param [inout]   a       Pointer on first value
 * \param [inout]   b       Pointer on second value
 *
 */

static inline void
_swap_long
(
PDM_g_num_t *a,
PDM_g_num_t *b
)
{
  PDM_g_num_t tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

/**
 *
 * \brief Swap element values
 *
 * \param [inout]   a       Pointer on first value
 * \param [inout]   b       Pointer on second value
 *
 */

static inline void
_swap_int
(
int *a,
int *b
)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

/**
 *
 * \brief Swap element values
 *
 * \param [inout]   a       Pointer on first value
 * \param [inout]   b       Pointer on second value
 *
 */

static inline void
_swap_double
(
double *a,
double *b
)
{
  double tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/



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
PDM_sort_long
(
 PDM_g_num_t  *array,
 int         *order,
 int          lArray
)
{
  /* size of subarray sorted by straight insertion */
  const int M = 7;
  /* default size of the stack */
  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  PDM_g_num_t a;
  int        b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  if (order != NULL) {

    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          b = order[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
            order[i+1] = order[i];
          }
          array[i+1] = a;
          order[i+1] = b;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else{
        int k = (l+ir) / 2;
        _swap_long (&(array[k]), &(array[l+1]));
        _swap_int (&(order[k]), &(order[l+1]));
        if (array[l] > array[ir]){
          _swap_long (&(array[l]), &(array[ir]));
          _swap_int (&(order[l]), &(order[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_long (&(array[l+1]), &(array[ir]));
          _swap_int (&(order[l+1]), &(order[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_long (&(array[l]), &(array[l+1]));
          _swap_int (&(order[l]), &(order[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        b = order[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_long (&(array[i]), &(array[j]));
          _swap_int (&(order[i]), &(order[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        order[l+1] = order[j];
        order[j] = b;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  else {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else{
        int k = (l+ir) / 2;
        _swap_long (&(array[k]), &(array[l+1]));
        if (array[l] > array[ir]){
          _swap_long (&(array[l]), &(array[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_long (&(array[l+1]), &(array[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_long (&(array[l]), &(array[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_long (&(array[i]), &(array[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  free (istack);
  return;

}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout]   array        Array to sort
 * \param [inout]   order        new indice to old indice (or NULL)
 * \param [in]      lArray       Array length
 *
 */

void
PDM_sort_int
(
 int         *array,
 int         *order,
 int          lArray
)
{
  /* size of subarray sorted by straight insertion */
  const int M = 7;
  /* default size of the stack */
  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  int  a;
  int  b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  if (order != NULL) {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          b = order[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
            order[i+1] = order[i];
          }
          array[i+1] = a;
          order[i+1] = b;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_int (&(array[k]), &(array[l+1]));
        _swap_int (&(order[k]), &(order[l+1]));
        if (array[l] > array[ir]){
          _swap_int (&(array[l]), &(array[ir]));
          _swap_int (&(order[l]), &(order[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_int (&(array[l+1]), &(array[ir]));
          _swap_int (&(order[l+1]), &(order[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_int (&(array[l]), &(array[l+1]));
          _swap_int (&(order[l]), &(order[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        b = order[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_int (&(array[i]), &(array[j]));
          _swap_int (&(order[i]), &(order[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        order[l+1] = order[j];
        order[j] = b;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  else {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_int (&(array[k]), &(array[l+1]));
        if (array[l] > array[ir]){
          _swap_int (&(array[l]), &(array[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_int (&(array[l+1]), &(array[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_int (&(array[l]), &(array[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_int (&(array[i]), &(array[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  free (istack);
  return;

}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout]   array        Array to sort
 * \param [inout]   order        new indice to old indice (or NULL)
 * \param [in]      lArray       Array length
 *
 */

void
PDM_sort_double
(
 double    *array,
 int       *order,
 int        lArray
)
{
  /* size of subarray sorted by straight insertion */

  const int M = 7;
  /* default size of the stack */

  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  double  a;
  int  b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  if (order != NULL) {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          b = order[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
            order[i+1] = order[i];
          }
          array[i+1] = a;
          order[i+1] = b;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_double (&(array[k]), &(array[l+1]));
        _swap_int (&(order[k]), &(order[l+1]));
        if (array[l] > array[ir]){
          _swap_double (&(array[l]), &(array[ir]));
          _swap_int (&(order[l]), &(order[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_double (&(array[l+1]), &(array[ir]));
          _swap_int (&(order[l+1]), &(order[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_double (&(array[l]), &(array[l+1]));
          _swap_int (&(order[l]), &(order[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        b = order[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_double (&(array[i]), &(array[j]));
          _swap_int (&(order[i]), &(order[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        order[l+1] = order[j];
        order[j] = b;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  else {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_double (&(array[k]), &(array[l+1]));
        if (array[l] > array[ir]){
          _swap_double (&(array[l]), &(array[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_double (&(array[l+1]), &(array[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_double (&(array[l]), &(array[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_double (&(array[i]), &(array[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  free (istack);
  return;

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
