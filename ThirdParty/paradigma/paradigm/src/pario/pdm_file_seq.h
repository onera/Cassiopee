#ifndef __PDM_FILE_SEQ_H__
#define __PDM_FILE_SEQ_H__

#include "pdm.h"

/*============================================================================
 * Description d'un fichier sequentiel
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
 * Mode d'acces lecture, ecriture, lecture/ecriture
 *----------------------------------------------------------------------------*/

typedef enum {

  FICHIER_SEQ_MODE_LECTURE,
  FICHIER_SEQ_MODE_ECRITURE,
  FICHIER_SEQ_MODE_AJOUT

} PDM_file_seq_mode_t;

/*----------------------------------------------------------------------------
 * Mode de positionnement dans le fichier
 *----------------------------------------------------------------------------*/

typedef enum {

  FICHIER_SET_SEEK_SET,   /* Position a partir du debut du fichier */
  FICHIER_SET_SEEK_CUR,   /* Position a partir de la position courante */
  FICHIER_SET_SEEK_END    /* Position a partir de la fin du fichier */

} PDM_file_seq_seek_t;

/*----------------------------------------------------------------------------
 * Type decrivant un fichier sequentiel
 *----------------------------------------------------------------------------*/

typedef struct _PDM_file_seq_t PDM_file_seq_t;

/*============================================================================
 * Prototype des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui ouvre un fichier binaire a acces direct
 *
 *  parameters :
 *    nom          <-- Nom du fichier
 *    mode         <-- Mode d'acces
 *
 *----------------------------------------------------------------------------*/

PDM_file_seq_t *PDM_file_seq_open(const char *nom,
                                const PDM_file_seq_mode_t mode);

/*----------------------------------------------------------------------------
 *  Fonction d'ecriture
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *    taille_donnee      <-- Taille des donnees
 *    n_donnees          <-- Nombre de donnees
 *  Return
 *    n_donnees_ecrites      Nombre de donnees reellements ecrites
 *                           Erreur si n_donnees != n_donnees_ecrites
 *
 *----------------------------------------------------------------------------*/

PDM_g_num_t PDM_file_seq_write(PDM_file_seq_t *fichier,
                      const size_t      taille_donnee,
                      const PDM_g_num_t n_donnees,
                      void              *donnees);

/*----------------------------------------------------------------------------
 *  Fonction de lecture
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *    taille_donnee      <-- Taille des donnees
 *    n_donnees          <-- Nombre de donnees
 *  Return
 *    n_donnees_lues         Nombre de donnees reellements lues
 *                           Erreur si n_donnees != n_donnees_lues
 *
 *----------------------------------------------------------------------------*/

PDM_g_num_t PDM_file_seq_read(PDM_file_seq_t *fichier,
                     const size_t           taille_donnee,
                     const PDM_g_num_t      n_donnees,
                     void                   *donnees);

/*----------------------------------------------------------------------------
 *  Defini la position courante du fichier
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *    offset             <-- Position
 *    whence             <-- A partir :
 *                              - du debut du fichier : FICHIER_SEQ_SEEK_SET
 *                              - de la position courante : FICHIER_SEQ_SEEK_CUR
 *                              - de la fin du fchier : FICHIER_SEQ_SEEK_END
 *
 *----------------------------------------------------------------------------*/

void PDM_file_seq_seek
(PDM_file_seq_t     *fichier,
 long               offset,
 PDM_file_seq_seek_t whence);

/*----------------------------------------------------------------------------
 *  Retourne a la position courante du fichier
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *  Return
 *    offset                 Position courante du fichier
 *
 *----------------------------------------------------------------------------*/

long PDM_file_seq_tell(PDM_file_seq_t *fichier);

/*----------------------------------------------------------------------------
 *  Fonction qui ferme le fichier de donnees
 *
 *----------------------------------------------------------------------------*/

void PDM_file_seq_close(PDM_file_seq_t *fichier);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_FILE_SEQ_H__ */
