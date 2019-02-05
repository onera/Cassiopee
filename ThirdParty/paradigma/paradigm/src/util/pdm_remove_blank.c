/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_remove_blank.h"

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
PDM_remove_blank
(
const char *str1
)
{
  char *str_without_blank = NULL;

  if (str1 != NULL) {
    int l_str = strlen(str1);
    
    int imin = 0;
    int imax = 0;
    
    while ((imin < l_str) && (str1[imin] == ' '))
      imin++;
    
    while ((imax < l_str) && (str1[l_str-imax-1] == ' '))
      imax++;
    
    imax = l_str-imax-1;
    
    assert(imax >= imin);
    
    if ((imax == l_str) || (imin == l_str)) {
      str_without_blank = (char *) malloc(sizeof(char));
      str_without_blank[0] = '\0';
    }
    else {
      int size = imax - imin + 2;
      str_without_blank = (char *) malloc(sizeof(char) * size);
      int index = 0;
      for (int k = imin; k <= imax; k++)
        str_without_blank[index++] = str1[k];
      str_without_blank[index] = '\0';
    }
  }

  return str_without_blank;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
