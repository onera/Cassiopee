#ifndef __FILE_PAR_H__
#define __FILE_PAR_H__

/*============================================================================
 * Description d'un fichier parallele
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

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
 * Mode d'acces lecture, ecriture, ajout
 *----------------------------------------------------------------------------*/

typedef enum {

  FICHIER_PAR_MODE_LECTURE,
  FICHIER_PAR_MODE_ECRITURE,
  FICHIER_PAR_MODE_AJOUT

} PDM_file_par_mode_t;

/*----------------------------------------------------------------------------
 * Type d'acces MPIIO
 *----------------------------------------------------------------------------*/

typedef enum {

  FICHIER_PAR_ACCES_EO,
  FICHIER_PAR_ACCES_IP

} PDM_file_par_acces_t;

/*----------------------------------------------------------------------------
 * Mode de positionnement dans le fichier
 *----------------------------------------------------------------------------*/

typedef enum {

  FICHIER_PAR_SEEK_SET,   /* Position a partir du debut du fichier */
  FICHIER_PAR_SEEK_CUR,   /* Position a partir de la position courante */
  FICHIER_PAR_SEEK_END    /* Position a partir de la fin du fichier */

} PDM_file_par_seek_t;

/*----------------------------------------------------------------------------
 * Type decrivant un fichier sequentiel
 *----------------------------------------------------------------------------*/

typedef struct _PDM_file_par_t PDM_file_par_t;

/*============================================================================
 * Prototype des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ouverture d'un fichier parallele
 *
 * parameters :
 *   nom             <-- Nom du fichier
 *   acces           <-- Access (implicit pointer, explicit offset)
 *   mode            <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
 *   pdm_mpi_comm        <-- Communicateur lie au fichier
 * return :
 *   PDM_file_par         Pointeur sur le fichier
 *
 *----------------------------------------------------------------------------*/

PDM_file_par_t *
PDM_file_par_open
(const char                *nom,
 const PDM_file_par_acces_t  acces,
 const PDM_file_par_mode_t   mode,
 PDM_MPI_Comm                   comm
 );

/*----------------------------------------------------------------------------
 * Lecture globale : Chaque processus lit la meme zone du fichier
 *
 * parameters :
 *   PDM_file_par     <-- Pointeur sur le fichier
 *   taille_donnee   <-- Taille unitaire de la donnnee
 *   n_donnees       <-- Nombre de donnees a lire
 *   donnees         --> Donnees lues
 * return
 *   n_donnees_lues      Nombre de donnees lues
 *                       Erreur de lecture si n_donnees != n_donnees_lues
 *
 *----------------------------------------------------------------------------*/

int
PDM_file_par_lecture_globale
(PDM_file_par_t *PDM_file_par,
 const size_t   taille_donnee,
 const int      n_donnees,
 void          *donnees
 );

/*----------------------------------------------------------------------------
 * Ecriture globale : Le processus maitre accede seul au fichier
 *
 * parameters :
 *   PDM_file_par       <-- Pointeur sur le fichier
 *   taille_donnee     <-- Taille unitaire de la donnnee
 *   n_donnees         <-- Nombre de donnees a ecrire
 *   donnees            --> Donnees lues
 * return
 *   n_donnees_ecrites      Nombre de donnees ecrites
 *                          Erreur d'ecriture si n_donnees != n_donnees_ecrites
 *
 *----------------------------------------------------------------------------*/

int
PDM_file_par_ecriture_globale
(PDM_file_par_t *PDM_file_par,
 const size_t   taille_donnee,
 const int      n_donnees,
 void          *donnees
 );

/*----------------------------------------------------------------------------
 * Lecture parallele de blocs de donnees suivie d'une redistribution des
 * des donnees suivant l'indirection
 *
 * parameters :
 *   PDM_file_par     <-- Pointeur sur le fichier
 *   n_composantes   <-- Nombre de composantes de la donnee
 *                       (Taille du sous bloc dans l'enregistrement)
 *   taille_donnee   <-- Taille unitaire de la donnnee en octet
 *   n_donnees_bloc  <-- Nombre de donnees a lire dans le bloc
 *                       (Nombre d'enregistrement du bloc)
 *   donnees         --> Donnees lues
 *   debut_bloc      <-- Debut du bloc dans l'ensemble des donnees
 *                       (Adresse en nombre d'enregistrements pour ce bloc)
 * return
 *   n_donnees_lues      Nombre de donnees lues
 *                       Erreur de lecture si n_donnees != n_donnees_lues
 *
 *----------------------------------------------------------------------------*/

int
PDM_file_par_lecture_parallele
(PDM_file_par_t *PDM_file_par,
 const size_t   taille_donnee,
 const int      n_donnees_bloc,
 void          *donnees,
 const long     debut_bloc
 );

/*----------------------------------------------------------------------------
 * Ecriture parallele de blocs de donnees tries
 *
 * parameters :
 *   PDM_file_par       <-- Pointeur sur le fichier
 *   n_composantes     <-- Nombre de composantes de la donnee
 *   taille_donnee     <-- Taille unitaire de la donnnee
 *   n_donnees         <-- Nombre de donnees a lire
 *   donnees           <-- Donnees a ecrire
 *   debut_bloc        <-- Debut du bloc dans l'ensemble des donnees
 * return
 *   n_donnees_ecrites      Nombre de donnees ecrites
 *                          Erreur d'ecriture si n_donnees != n_donnees_ecrites
 *
 *----------------------------------------------------------------------------*/

int
PDM_file_par_ecriture_parallele
(PDM_file_par_t *PDM_file_par,
 const size_t   taille_donnee,
 const int      n_donnees_bloc,
 void          *donnees,
 const long     debut_bloc
 );

/*----------------------------------------------------------------------------
 *  Defini la position courante du fichier
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *    offset             <-- Position
 *    whence             <-- A partir :
 *                              - du debut du fichier : FICHIER_PAR_SEEK_SET
 *                              - de la position courante : FICHIER_PAR_SEEK_CUR
 *                              - de la fin du fchier : FICHIER_PAR_SEEK_END
 *
 *----------------------------------------------------------------------------*/

void PDM_file_par_seek
(PDM_file_par_t     *fichier,
 PDM_MPI_Offset         offset,
 PDM_file_par_seek_t whence);

/*----------------------------------------------------------------------------
 *  Retourne a la position courante du fichier
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *  Return
 *    offset                 Position courante du fichier
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Offset PDM_file_par_tell(PDM_file_par_t *fichier);

/*----------------------------------------------------------------------------
 * Fermeture du fchier
 *
 * parameters :
 *   fichier            <-- Fichier courant
 *
 *----------------------------------------------------------------------------*/

void PDM_file_par_close (PDM_file_par_t *PDM_file_par);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_PAR_H__ */
