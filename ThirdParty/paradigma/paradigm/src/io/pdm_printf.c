/*============================================================================
 * Base user-definable PDM_printf() wrapper or replacement.
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2017 ONERA

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

#include "pdm_config.h"

/*
 * Standard C library headers
 */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Optional library and BFT headers
 */

#include "pdm_printf.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/* Associated typedef documentation (for PDM_printf.h) */

/*!
 * \typedef PDM_printf_proxy_t
 *
 * \brief Function pointer for PDM_printf() type functions.
 *
 * \param [in] format       format string, as PDM_printf() and family.
 * \param [in, out] arg_ptr pointer to variable argument list based on format
 *                          string.
 */

/*!
 * \typedef PDM_printf_flush_proxy_t
 *
 * \brief Function pointer for fflush(stdout) type functions.
 */

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default PDM_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_PDM_printf_flush_proxy_default(void);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static PDM_printf_proxy_t        *_PDM_printf_proxy = vprintf;
static PDM_printf_flush_proxy_t  *_PDM_printf_flush_proxy
                                    = _PDM_printf_flush_proxy_default;

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Default PDM_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_PDM_printf_flush_proxy_default(void)
{
  return fflush(stdout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Replacement for PDM_printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by PDM_printf_proxy_set().
 *
 * \param [in] format format string, as PDM_printf() and family.
 * \param [in] ...    variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\\0'
 *         used to end output strings
 */

int
PDM_printf(const char *const format,
           ...)
{
  int  retval;
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  retval = _PDM_printf_proxy(format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

/*!
 * \brief Flush for output of PDM_printf() with modifiable behavior.
 *
 * This function calls fflush(stdout) if PDM_printf()'s default behavior is
 * used. If PDM_printf's behavior is modified with PDM_printf_proxy_set(),
 * PDM_printf_flush()'s behavior may have to be also adjusted with
 * PDM_printf_flush_proxy_set().
 *
 * \return using the default behavior, the return value is that of
 *         fflush(stdout): O upon successful completion, EOF otherwise
 *         (with errno set to indicate the error).
 */

int
PDM_printf_flush(void)
{
  return _PDM_printf_flush_proxy();
}

/*!
 * \brief Returns function associated with the PDM_printf() function.
 *
 * \return pointer to the vprintf() or replacement function.
 */

PDM_printf_proxy_t *
PDM_printf_proxy_get(void)
{
  return _PDM_printf_proxy;
}

/*!
 * \brief Associates a vprintf() type function with the PDM_printf() function.
 *
 * \param [in] fct pointer to a vprintf() type function.
 */

void
PDM_printf_proxy_set(PDM_printf_proxy_t *const fct)
{
  _PDM_printf_proxy = fct;
}

/*!
 * \brief Returns function associated with PDM_printf_flush().
 *
 * \return pointer to the PDM_printf_flush() proxy.
 */

PDM_printf_flush_proxy_t *
PDM_printf_flush_proxy_get(void)
{
  return _PDM_printf_flush_proxy;
}

/*!
 * \brief Associates a proxy function with PDM_printf_flush().
 *
 * \warning
 *   PDM_printf() is called by the default PDM_error() error handler
 *   (so as to ensure that the error text appears at the end of the
 *   program output), so a PDM_print_flush replacement must not itself
 *   call (directly or indirectly) PDM_error() if the default error
 *   handler is used.
 *
 * \param [in] fct pointer to a function similar to {return fflush(stdout)}.
 */

void
PDM_printf_flush_proxy_set(PDM_printf_flush_proxy_t *const fct)
{
  _PDM_printf_flush_proxy = fct;
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
