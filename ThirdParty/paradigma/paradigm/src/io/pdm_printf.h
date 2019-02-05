#ifndef __PDM_PRINTF_H__
#define __PDM_PRINTF_H__

/*============================================================================
 * Base user-definable printf() wrapper or replacement
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */

#include <stdarg.h>

/* BFT library headers */

#include "pdm_config.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/* Function pointers for printf() and fflush(stdout) type functions */

typedef int (PDM_printf_proxy_t) (const char     *const format,
                                  va_list               arg_ptr);

typedef int (PDM_printf_flush_proxy_t) (void);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Replacement for printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by PDM_printf_proxy_set().
 *
 * parameters:
 *   format: <-- format string, as printf() and family.
 *   ... :   <-- variable arguments based on format string.
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0' used
 *   to end output strings
 */

int
PDM_printf(const char  *const format,
           ...);

/*
 * Flush for output of PDM_printf() with modifiable behavior.
 *
 * This function calls fflush(stdout) if PDM_printf()'s default behavior is
 * used. If PDM_printf's behavior is modified with PDM_printf_proxy_set(),
 * PDM_printf_flush()'s behavior may have to be also adjusted with
 * PDM_printf_flush_proxy_set().
 *
 * returns:
 *   using the default behavior, the return value is that of
 *   fflush(stdout): O upon successful completion, EOF otherwise
 *   (with errno set to indicate the error).
 */

int
PDM_printf_flush(void);

/*
 * Returns function associated with the PDM_printf() function.
 *
 * returns:
 *   pointer to the vprintf() or replacement function.
 */

PDM_printf_proxy_t *
PDM_printf_proxy_get(void);

/*
 * Associates a vprintf() type function with the PDM_printf() function.
 *
 * parameters:
 *   fct: <-- pointer to a vprintf() type function.
 */

void
PDM_printf_proxy_set(PDM_printf_proxy_t  *const fct);

/*
 * Returns function associated with PDM_printf_flush().
 *
 * returns:
 *   pointer to the PDM_printf_flush() proxy.
 */

PDM_printf_flush_proxy_t *
PDM_printf_flush_proxy_get(void);

/*
 * Associates a proxy function with PDM_printf_flush().
 *
 * warning:
 *   PDM_printf() is called by the default PDM_error() error handler
 *   (so as to ensure that the error text appears at the end of the
 *   program output), so a PDM_print_flush replacement must not itself
 *   call (directly or indirectly) PDM_error() if the default error
 *   handler is used.
 *
 * parameter:
 *   fct <-- pointer to a function similar to {return fflush(stdout)}.
 */

void
PDM_printf_flush_proxy_set(PDM_printf_flush_proxy_t  *const fct);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PRINTF_H__ */
