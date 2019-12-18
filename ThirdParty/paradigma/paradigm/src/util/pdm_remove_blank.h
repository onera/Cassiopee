#ifndef __PDM_REMOVE_BLANK_H__
#define __PDM_REMOVE_BLANK_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* REMOVE_BLANK */
