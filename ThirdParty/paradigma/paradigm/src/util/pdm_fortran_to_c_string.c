/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_fortran_to_c_string.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 *
 * Conversion d'une chaine fortran en chaine C en retirant les blancs
 *
 * parameters:
 *   chaine_f    <-- Chaine Fortran
 *   l_chaine_f  <-- Longueur de la chaine fortran
 *
 * return:
 *   C string
 *
 *----------------------------------------------------------------------------*/

char *
PDM_fortran_to_c_string
(
 const char *chaine_f,
 const PDM_l_num_t l_chaine_f
)
{
  char *chaine_c = NULL;
  int imin = 0;
  int imax = 0;

  while ((imin < l_chaine_f) && (chaine_f[imin] == ' '))
    imin++;

  while ((imax < l_chaine_f) && (chaine_f[l_chaine_f-imax-1] == ' '))
    imax++;

  imax = l_chaine_f-imax-1;

  assert(imax >= imin);

  if ((imax == l_chaine_f) || (imin == l_chaine_f)) {
    chaine_c = (char *) malloc(sizeof(char));
    chaine_c[0] = '\0';
  }
  else {
    int size = imax - imin + 2;
    chaine_c = (char *) malloc(sizeof(char) * size);;
    int index = 0;
    for (int k = imin; k <= imax; k++)
      chaine_c[index++] = chaine_f[k];
    chaine_c[index] = '\0';
  }

  return chaine_c;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
