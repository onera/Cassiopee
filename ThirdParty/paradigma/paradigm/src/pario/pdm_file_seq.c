/*============================================================================
 * Description d'un fichier scalaire
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_file_seq.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Type decrivant un fichier de type MSG (MPI_IO)
 *----------------------------------------------------------------------------*/

struct _PDM_file_seq_t {

  FILE               *fichier;     /* Pointeur sur le fichier C */
  char               *nom;     /* Nom du fichier */
  PDM_file_seq_mode_t  mode;    /* Mode */

};

/*============================================================================
 * Variables globales locales
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  PDM_file_par_sek -> pdm_mpi_file_seek_t
 *----------------------------------------------------------------------------*/

static int _std_file_seek[3] = {SEEK_SET,
                                      SEEK_CUR,
                                      SEEK_END};

/*============================================================================
 * Definitions des fonctions locales
 *============================================================================*/

/*============================================================================
 * Definitions des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui ouvre un fichier binaire a acces direct
 *
 *  parameters :
 *    nom            <-- Nom du fichier
 *    mode           <-- Mode d'acces
 *
 *----------------------------------------------------------------------------*/

PDM_file_seq_t *PDM_file_seq_open(const char *nom,
                                const PDM_file_seq_mode_t mode)
{

  PDM_file_seq_t *fichier = (PDM_file_seq_t *) malloc(sizeof(PDM_file_seq_t));

  fichier->nom = (char *) malloc(strlen(nom) + 1);
  strcpy(fichier->nom, nom);
  fichier->mode = mode;

  switch (mode) {
  case FICHIER_SEQ_MODE_LECTURE:
    fichier->fichier = fopen(nom, "r");
    break;
  case FICHIER_SEQ_MODE_ECRITURE:
    fichier->fichier = fopen(nom, "w");
    break;
  case FICHIER_SEQ_MODE_AJOUT:
    fichier->fichier = fopen(nom, "a");
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_file_seq_open :\n"
                    "Mode de fichier inconnu\n");
    exit(EXIT_FAILURE);
  }

  if (fichier->fichier == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_file_seq_open :\n"
            "Erreur à l'ouverture du fichier %s\n", nom);
    exit(EXIT_FAILURE);
  }

  return fichier;

}

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
                      const size_t   taille_donnee,
                      const PDM_g_num_t n_donnees,
                      void          *donnees)
{

  if (fichier->mode == FICHIER_SEQ_MODE_LECTURE) {
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_file_seq_open :\n"
                    "Ecriture interdite pour le fichier %s "
                    "ouvert en mode lecture\n",
                    fichier->nom);
    exit(EXIT_FAILURE);
  }

  size_t _n_donnees = (size_t) n_donnees ;
  size_t _n_donnees_ecrites = fwrite(donnees, taille_donnee,
                                       _n_donnees, fichier->fichier);
  PDM_g_num_t n_donnees_ecrites = (PDM_g_num_t) _n_donnees_ecrites;

  return n_donnees_ecrites;

}

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
                     const size_t   taille_donnee,
                     const PDM_g_num_t n_donnees,
                     void          *donnees)
{

  if (fichier->mode == FICHIER_SEQ_MODE_ECRITURE ||
      fichier->mode == FICHIER_SEQ_MODE_AJOUT) {
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_file_seq_open :\n"
                   "Lecture interdite pour le fichier %s "
                   "ouvert en mode ecriture/ajout\n",
                   fichier->nom);
    exit(EXIT_FAILURE);
  }

  size_t _n_donnees = (size_t) n_donnees ;
  size_t _n_donnees_lues = fread(donnees, taille_donnee,
                                   _n_donnees, fichier->fichier);
  PDM_g_num_t n_donnees_lues = (PDM_g_num_t) _n_donnees_lues;

  return n_donnees_lues;

}

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
 PDM_file_seq_seek_t whence)
{
  fseek(fichier->fichier, offset, _std_file_seek[whence]);
}

/*----------------------------------------------------------------------------
 *  Retourne a la position courante du fichier
 *
 *  parameters :
 *    fichier            <-- Fichier courant
 *  Return
 *    offset                 Position courante du fichier
 *
 *----------------------------------------------------------------------------*/

long PDM_file_seq_tell(PDM_file_seq_t *fichier)
{
  return ftell(fichier->fichier);
}

/*----------------------------------------------------------------------------
 *  Fonction qui ferme le fichier de donnees
 *
 *----------------------------------------------------------------------------*/

void PDM_file_seq_close(PDM_file_seq_t *fichier)
{
  free(fichier->nom);
  fclose(fichier->fichier);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
