#ifndef __PDM_FORTRAN_TO_C_STRING_H__
#define __PDM_FORTRAN_TO_C_STRING_H__

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
PDM_fortran_to_c_string
(
 const char *chaine_f,
 const PDM_l_num_t l_chaine_f
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_FORTRAN_TO_C_STRING_H__ */
