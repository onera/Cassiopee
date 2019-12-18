#ifndef __PDM_WRITER_ENSIGHT_CASE_H__
#define __PDM_WRITER_ENSIGHT_CASE_H__

/*============================================================================
 * Manage case files associated with the EnSight Gold writer
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_writer.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Opaque structure to manage case file */

typedef struct _PDM_writer_ensight_case_t PDM_writer_ensight_case_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a new case file structure.
 *
 * parameters:
 *   name            <-- case name
 *   restart         <-- if restart == 1, case file is read
 *   dir_prefix      <-- associated local or absolute directory name
 *   time_dependency <-- indicates if and how meshes will change with time
 *
 * returns:
 *   pointer to new case file structure
 *----------------------------------------------------------------------------*/

PDM_writer_ensight_case_t *
PDM_writer_ensight_case_cree
(
const char                   *const name,
const int                           restart,
const char                   *const dir_prefix,
const PDM_writer_topologie_t                time_dependency
);

/*----------------------------------------------------------------------------
 * Destroy a case file structure.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

PDM_writer_ensight_case_t *
PDM_writer_ensight_case_lib
(
PDM_writer_ensight_case_t  *this_case
);

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

PDM_writer_topologie_t
PDM_writer_ensight_case_geo_time_dep_get
(
PDM_writer_ensight_case_t  *this_case
);

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight var
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

PDM_writer_statut_t
PDM_writer_ensight_case_var_time_dep_get
(
PDM_writer_ensight_case_t  *this_case,
const char         *name
);

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

char *
PDM_writer_ensight_case_var_file_name_get
(
PDM_writer_ensight_case_t  *this_case,
const char* name
);

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

char *
PDM_writer_ensight_case_geo_file_name_get
(
PDM_writer_ensight_case_t  *this_case
);

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_case_var_cree
(
PDM_writer_ensight_case_t  *this_case,
const char         *const name,
const PDM_writer_var_dim_t  dimension,
const PDM_writer_statut_t   time_dep,
const PDM_writer_var_loc_t  location
);


/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_case_time_step_add
(
PDM_writer_ensight_case_t  *this_case,
const double time_value
);

/*----------------------------------------------------------------------------
 * Write an EnSight Gold case file.
 *
 * This function should only be called by one process in parallel mode.
 *
 * parameters:
 *   this_case  <-- case structure
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_case_write
(
PDM_writer_ensight_case_t  *const this_case,
int                       rank
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_ENSIGHT_CASE_H__ */
