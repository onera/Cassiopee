#ifndef __PDM_IO_TAB_H__
#define __PDM_IO_TAB_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_io.h"

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

/*----------------------------------------------------------------------------
 * Types decrivant un tableau
 *----------------------------------------------------------------------------*/

typedef struct _PDM_io_tab_t PDM_io_tab_t;

/*=============================================================================
 * Variables globales
 *============================================================================*/


/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Initialise une phase d'écriture parallèle de tableaux de données associées
 * aux numéros de variables PDM
 * Chaque tableau a ses propres caractéristiques :
 *         - taille de données
 *         - nombre de donnée
 *         - indirection (numérotation absolue)
 *
 * arguments :
 *   unite             <-- Unite du fichier
 *   t_rangement       <-- Type de rangement
 *   num_var_cedre_max <-- Numéro max de variable PDM
 *   n_partition_local <-- Nombre de partitions locales
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_ecr_debut, PDM_IO_TAB_ECR_DEBUT)
(const PDM_l_num_t *unite,
 const int            *t_rangement,
 const PDM_l_num_t *num_var_cedre_max,
 const PDM_l_num_t *n_partition_local
);

void PDM_io_tab_ecr_debut
(const PDM_l_num_t       unite,
 const PDM_io_rangement_t t_rangement,
 const PDM_l_num_t       num_var_cedre_max,
 const PDM_l_num_t       n_partition_local
);


/*----------------------------------------------------------------------------
 * Ajoute une partie des donnees dans un tableau associés à une variable
 * PDM
 *
 * arguments :
 *   num_var_cedre         <-- Numéro de variable PDM
 *   i_part                <-- indice de partition
 *   n_composantes         <-- Nombre de composantes pour chaque donnee
 *   n_donnees             <-- Nombre de donnees a lire
 *   indirection           <-- Indirection de redistribition des donnees
 *   donnees               <-- Donnees a écrire
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_ecr_ajout_donnees, PDM_IO_TAB_ECR_AJOUT_DONNEES)
(const PDM_l_num_t  *num_var_cedre,
 const PDM_l_num_t  *i_part,
 const PDM_l_num_t  *n_composantes,
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *indirection,
 void                  *donnees
 );

void PDM_io_tab_ecr_ajout_donnees
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            i_part,
 const PDM_l_num_t           *n_composantes,
 const PDM_l_num_t            n_donnees,
 const PDM_g_num_t          *indirection,
 void                           *donnees
 );

/*----------------------------------------------------------------------------
 * Definition d'une variable en ecriture
 *
 *  arguments :
 *   num_var_cedre         <-- Numéro de variable PDM
 *   num_indirection_cedre <-- Numéro d'indirection PDM
 *   t_n_composantes       <-- Type de tailles composantes
 *                             (PDM_IO_N_COMPOSANTE_CONSTANT
 *                           ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes         <-- Nombre de composantes pour chaque donnee
 *   taille_donnee         <-- Taille unitaire de la donnnee
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_ecr_def_var, PDM_IO_TAB_ECR_DEF_VAR)
(const PDM_l_num_t  *num_var_cedre,
 const PDM_l_num_t  *num_indirection_cedre,
 const int             *t_n_composantes,
 const PDM_l_num_t  *n_composantes,
 const PDM_l_num_t  *taille_donnee
 );

void PDM_io_tab_ecr_def_var
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            num_indirection_cedre,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            n_composantes,
 const PDM_l_num_t            taille_donnee
 );

/*----------------------------------------------------------------------------
 * Finalise une phase d'écriture parallèle de tableaux de données associées
 * aux numéros de variables PDM. Cette fonction déclenche réellement
 * les écritures
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_ecr_fin, PDM_IO_TAB_ECR_FIN)
(void);

void PDM_io_tab_ecr_fin
(void);


/*----------------------------------------------------------------------------
 * Initialise une phase de lecture parallèle de tableaux de données associées
 * aux numéros de variables PDM
 * Chaque tableau a ses propres caractéristiques :
 *         - taille de données
 *         - nombre de donnée
 *         - indirection (numérotation absolue)
 *
 * arguments :
 *   unite             <-- Unite du fichier
 *   t_rangement       <-- Type de rangement
 *   num_var_cedre_max <-- Numéro max de variable PDM
 *   n_partition_local <-- Nombre de partitions locales
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_lec_debut, PDM_IO_TAB_LEC_DEBUT)
(const PDM_l_num_t *unite,
 const int            *t_rangement,
 const PDM_l_num_t *num_var_cedre_max,
 const PDM_l_num_t *n_partition_local
);

void PDM_io_tab_lec_debut
(const PDM_l_num_t       unite,
 const PDM_io_rangement_t t_rangement,
 const PDM_l_num_t       num_var_cedre_max,
 const PDM_l_num_t       n_partition_local
);


/*----------------------------------------------------------------------------
 * Ajoute une partie des donnees dans un tableau associés à une variable
 * PDM
 *
 * arguments :
 *   num_var_cedre         <-- Numéro de variable PDM
 *   num_indirection_cedre <-- Numéro d'indirection PDM
 *   i_part                <-- indice de partition
 *   t_n_composantes       <-- Type de tailles composantes
 *                             (PDM_IO_N_COMPOSANTE_CONSTANT
 *                           ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes         <-- Nombre de composantes pour chaque donnee
 *   taille_donnee         <-- Taille unitaire de la donnnee
 *   n_donnees             <-- Nombre de donnees a lire
 *   indirection           <-- Indirection de redistribition des donnees
 *   donnees               <-- Donnees a écrire
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_lec_ajout_donnees, PDM_IO_TAB_LEC_AJOUT_DONNEES)
(const PDM_l_num_t  *num_var_cedre,
 const PDM_l_num_t  *i_part,
 const PDM_l_num_t  *n_composantes,
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *indirection,
 void                  *donnees
 );

void PDM_io_tab_lec_ajout_donnees
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            i_part,
 const PDM_l_num_t           *n_composantes,
 const PDM_l_num_t            n_donnees,
 const PDM_g_num_t          *indirection,
 void                           *donnees
 );

/*----------------------------------------------------------------------------
 * Definition d'une variable en ecriture
 *
 *  arguments :
 *   num_var_cedre         <-- Numéro de variable PDM
 *   num_indirection_cedre <-- Numéro d'indirection PDM
 *   t_n_composantes       <-- Type de tailles composantes
 *                             (PDM_IO_N_COMPOSANTE_CONSTANT
 *                           ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes         <-- Nombre de composantes pour chaque donnee
 *   taille_donnee         <-- Taille unitaire de la donnnee
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_lec_def_var, PDM_IO_TAB_LEC_DEF_VAR)
(const PDM_l_num_t  *num_var_cedre,
 const PDM_l_num_t  *num_indirection_cedre,
 const int             *t_n_composantes,
 const PDM_l_num_t  *n_composantes,
 const PDM_l_num_t  *taille_donnee
 );

void PDM_io_tab_lec_def_var
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            num_indirection_cedre,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            n_composantes,
 const PDM_l_num_t            taille_donnee
 );

/*----------------------------------------------------------------------------
 * Finalise une phase de lecture parallèle de tableaux de données associées
 * aux numéros de variables PDM. Cette fonction déclenche réellement
 * les écritures
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_lec_fin, PDM_IO_TAB_LEC_FIN)
(void);

void PDM_io_tab_lec_fin
(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_IO_TAB_H__ */
