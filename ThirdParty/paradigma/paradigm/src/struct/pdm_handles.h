#ifndef __PDM_HANDLES_H__
#define __PDM_HANDLES_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

/**
 * \brief Handles storage
 */

typedef struct _PDM_Handles_t PDM_Handles_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/


/**
 * \brief Create handles storage
 *
 * \param [in] init_size      Initial size of storage array
 *
 * \return  New handles storage
 */

PDM_Handles_t *
PDM_Handles_create 
(
 const int init_size 
);


/**
 * \brief Free handles storage
 *
 * \param [in] handles  Current handles storage
 *
 * \return               NULL
 */

PDM_Handles_t *
PDM_Handles_free 
(
 PDM_Handles_t *handles
);

/**
 * \brief Store a new handle pointer
 *
 * \param [in] handles  Current handles storage
 * \param [in] handle_ptr   Handle pointer
 *
 * \return  Handle index
 */

int
PDM_Handles_store 
(
 PDM_Handles_t *handles, 
 const void *handle_ptr
);


/**
 * \brief Store a new pointer
 *
 * \param [in] handles  Current handles storage
 * \param [in] handle_idx   Handle index
 *
 * \return  Handle pointer
 */

const void *
PDM_Handles_get 
(
  PDM_Handles_t *handles,
  const int handle_idx
);


/**
 * \brief Get handles index
 *
 * \param [in] handles  Current handles storage
 *
 * \return  Handles index
 */

const int *
PDM_Handles_idx_get 
(
 PDM_Handles_t *handles
 );


/**
 * \brief Free a handle
 *
 * \param [in] handles      Current handles storage
 * \param [in] handle_idx   Handle index
 * \param [in] st_free_data Free data or not
 *
 * \return  Handle pointer
 */

void
PDM_Handles_handle_free 
(
 PDM_Handles_t *handles,
 const int handle_idx,
 const PDM_bool_t st_free_data
 );


/**
 * \brief Get number of stored handles
 *
 * \param [in] handles  Current handles storage
 *
 * \return  Number of stored handles
 */

int
PDM_Handles_n_get 
(
 PDM_Handles_t *handles
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_HANDLES_H__ */
