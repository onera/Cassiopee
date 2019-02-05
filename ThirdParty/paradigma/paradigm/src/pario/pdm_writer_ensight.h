#ifndef __PDM_WRITER_ENSIGHT_H__
#define __PDM_WRITER_ENSIGHT_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_writer_priv.h"

/*=============================================================================
 * Definitions des macro
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*=============================================================================
 * Variables globales
 *============================================================================*/

/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Particularise la creation d'un objet CS (Cedre Sortie)
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_create 
(
PDM_writer_t *cs
);

/*----------------------------------------------------------------------------
 * Particularise la liberation d'un objet CS (Cedre Sortie)
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_free
(
PDM_writer_t *cs
);

/*----------------------------------------------------------------------------
 * Particularise le debut d'increment
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_step_beg
(
PDM_writer_t *cs
);

/*----------------------------------------------------------------------------
 * Particularise la fin d'increment
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_step_end
(
PDM_writer_t *cs
);


/*----------------------------------------------------------------------------
 * Particularise la creation de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_geom_create
(
PDM_writer_geom_t *geom
);


/*----------------------------------------------------------------------------
 * Particularise la creation d'une variable
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_var_create
(
PDM_writer_var_t *var
);


/*----------------------------------------------------------------------------
 * Particularise l'ecriture de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_geom_write
(
 PDM_writer_geom_t    *geom
);


/*----------------------------------------------------------------------------
 * Particularise la libération de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_geom_free
(
 PDM_writer_geom_t    *geom
);
 

/*----------------------------------------------------------------------------
 * Particularise l'ecriture d'une variable 
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_var_write
(
 PDM_writer_var_t        *var
);
 

/*----------------------------------------------------------------------------
 * Particularise la libération d'une variable
 *
 * parameters :
 *      var          <-> Objet variable à libérer
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_var_free
(
 PDM_writer_var_t     *var
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_ENSIGHT_H__ */
