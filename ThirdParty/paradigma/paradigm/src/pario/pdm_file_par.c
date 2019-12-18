/*============================================================================
 * Description d'un fichier parallele
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_file_par.h"
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

struct _PDM_file_par_t {
  char                *nom;          /* Nom */
  PDM_file_par_mode_t   mode;         /* Mode */
  PDM_file_par_acces_t  acces;        /* Access */
  int                  rang;         /* Rang MSG */
  int                  n_rangs;      /* Nombre de rangs MSG  */
  PDM_MPI_Comm             comm;         /* Communicateur MSG associe a pa_io_fichier_t */
  PDM_MPI_File             fichier;      /* MSG file*/
  PDM_MPI_Offset           offset;       /* Offset du fichier */
};

/*============================================================================
 * Variables globales locales
 *============================================================================*/


/*============================================================================
 * Definitions des fonctions locales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Message d'erreur MPI IO.
 *
 * parameters:
 *   nom_fichier  <-- Nom du fichier
 *   code_erreur  <-- Code d'erreur
 *
 *----------------------------------------------------------------------------*/

static void
_pdm_mpi_io_error_message
(
 const char  *nom_fichier,
 int          code_erreur
)
{
  char buffer[PDM_MPI_get_max_error_string()];
  int  buffer_len;

  PDM_MPI_Error_string(code_erreur, buffer, &buffer_len);

  PDM_error(__FILE__, __LINE__, 0, "Erreur MPI IO pour le fichier: %s\n"
          "Error type: %s", nom_fichier, buffer);

  exit(EXIT_FAILURE);
}

/*============================================================================
 * Definitions des fonctions publiques
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

PDM_file_par_t
*PDM_file_par_open
(const char                *nom,
 const PDM_file_par_acces_t  acces,
 const PDM_file_par_mode_t   mode,
 PDM_MPI_Comm                   comm
 )
{
  int _mode = PDM_MPI_MODE_RDWR;

  PDM_file_par_t *PDM_file_par = (PDM_file_par_t *) malloc(sizeof(PDM_file_par_t));

  PDM_file_par->nom = (char *) malloc(strlen(nom) + 1);
  strcpy(PDM_file_par->nom, nom);

  PDM_file_par->mode = mode;

  switch (PDM_file_par->mode) {
  case FICHIER_PAR_MODE_AJOUT:
    _mode = PDM_MPI_MODE_WRONLY_APPEND;
    break;
  case FICHIER_PAR_MODE_ECRITURE:
    _mode = PDM_MPI_MODE_WRONLY_CREATE;
    break;
  case FICHIER_PAR_MODE_LECTURE:
    _mode = PDM_MPI_MODE_RDONLY;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0,"Erreur Fichier_par_open : mode inconnu\n");
    exit(EXIT_FAILURE);
  }

  PDM_file_par->acces = acces;

  PDM_MPI_Comm_rank(comm, &(PDM_file_par->rang));

  PDM_MPI_Comm_size(comm, &(PDM_file_par->n_rangs));

  PDM_file_par->comm = comm;

  PDM_MPI_File_open(PDM_file_par->comm,
                PDM_file_par->nom,
                _mode,
                &(PDM_file_par->fichier));

  PDM_file_par->offset = 0;

  return PDM_file_par;

}

/*----------------------------------------------------------------------------
 * Lecture globale : Chaque processus lit la meme zone du fichier
 *
 * parameters :
 *   PDM_file_par     <-- Pointeur sur le fichier
 *   taille_donnee   <-- Taille unitaire du type de la donnee
 *   n_donnees       <-- Nombre de donnees a lire
 *   donnees         --> Donnees lues
 * return
 *   n_donnees_lues      Nombre de donnees lues
 *                       Erreur de lecture si n_donnees != n_donnees_lues
 *
 *----------------------------------------------------------------------------*/

int
PDM_file_par_lecture_globale
(PDM_file_par_t   *PDM_file_par,
 const size_t     taille_donnee,
 const int        n_donnees,
 void            *donnees
)
{

  size_t n_donnees_lues;
  int errcode = PDM_MPI_SUCCESS;
  int n_octet_lus = 0;

  /* Variables pour acces ip */

  PDM_MPI_Datatype fichier_type;
  PDM_MPI_Aint disps[1];
  int lengths[1];

  int _taille_donnee = (int) taille_donnee;

  /* Lecture */

  switch(PDM_file_par->acces) {

  case FICHIER_PAR_ACCES_EO :
    errcode = PDM_MPI_File_read_at_all(PDM_file_par->fichier,
                                   PDM_file_par->offset,
                                   donnees,
                                   _taille_donnee * n_donnees,
                                   PDM_MPI_BYTE,
                                   &n_octet_lus);
/*     PDM_MPI_Get_count(&status, PDM_MPI_BYTE, &n_octet_lus); */

    break;

  case FICHIER_PAR_ACCES_IP :
    lengths[0] = (int) (_taille_donnee * n_donnees);
    disps[0] = 0;
    PDM_MPI_Type_create_hindexed(1, lengths, disps, PDM_MPI_BYTE, &fichier_type);
    PDM_MPI_Type_commit(&fichier_type);
    PDM_MPI_File_set_view(PDM_file_par->fichier,
                      PDM_file_par->offset,
                      PDM_MPI_BYTE,
                      fichier_type,
                      "native");
    errcode = PDM_MPI_File_read_all(PDM_file_par->fichier,
                                donnees,
                                _taille_donnee * n_donnees,
                                PDM_MPI_BYTE,
                                &n_octet_lus);

/*     PDM_MPI_Get_count(&status, PDM_MPI_BYTE, &n_octet_lus); */

    PDM_MPI_Type_free(&fichier_type);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0,"Erreur Fichier_par_lecture : type d'acces inconnu\n");
    exit(EXIT_FAILURE);

  }

  /* Traitement des erreurs de lecture */

  if (errcode != PDM_MPI_SUCCESS)
    _pdm_mpi_io_error_message(PDM_file_par->nom, errcode);


  n_donnees_lues = n_octet_lus / taille_donnee;

  /* Mise à jour de l'offset */

  PDM_file_par->offset += n_octet_lus;
  int _n_donnees_lues = (int) n_donnees_lues;

  return _n_donnees_lues;
}

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
(PDM_file_par_t      *PDM_file_par,
 const size_t        taille_donnee,
 const int           n_donnees,
 void               *donnees
)
{
  int n_donnees_ecrites;
  int sorties [2] = {PDM_MPI_SUCCESS, 0};

  /* Variables pour acces ip */

  PDM_MPI_Datatype fichier_type;
  PDM_MPI_Aint disps[1];
  int lengths[1];
  int _taille_donnee = (int) taille_donnee;

  switch(PDM_file_par->acces) {

  case FICHIER_PAR_ACCES_EO :
    if (PDM_file_par->rang == 0) {
      sorties[0] = PDM_MPI_File_write_at(PDM_file_par->fichier,
                                     PDM_file_par->offset,
                                     donnees,
                                     _taille_donnee * n_donnees,
                                     PDM_MPI_BYTE,
                                     &(sorties[1]) );
/*       PDM_MPI_Get_count(&status, PDM_MPI_BYTE, &(sorties[1])); */
    }

    break;

  case FICHIER_PAR_ACCES_IP :
    lengths[0] = _taille_donnee * n_donnees;
    disps[0] = 0;
    PDM_MPI_Type_create_hindexed(1, lengths, disps, PDM_MPI_BYTE, &fichier_type);
    PDM_MPI_Type_commit(&fichier_type);
    PDM_MPI_File_set_view(PDM_file_par->fichier,
                      PDM_file_par->offset,
                      PDM_MPI_BYTE,
                      fichier_type,
                      "native");
    if (PDM_file_par->rang == 0) {
      sorties[0] = PDM_MPI_File_write(PDM_file_par->fichier,
                                  donnees,
                                  _taille_donnee * n_donnees,
                                  PDM_MPI_BYTE,
                                   &(sorties[1]));
/*       PDM_MPI_Get_count(&status, PDM_MPI_BYTE, &(sorties[1])); */
    }
    PDM_MPI_Type_free(&fichier_type);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0,"Erreur Fichier_par_ecriture_global : type d'acces inconnu\n");
    exit(EXIT_FAILURE);

  }

  PDM_MPI_Bcast(sorties, 2, PDM_MPI_INT, 0, PDM_file_par->comm);

  /* Traitement des erreurs d'ecriture */

  if (sorties[0] != PDM_MPI_SUCCESS)
    _pdm_mpi_io_error_message(PDM_file_par->nom, sorties[0]);

  n_donnees_ecrites = sorties[1] / _taille_donnee;

  /* Mise à jour de l'offset */

  PDM_file_par->offset += sorties[1];

  return n_donnees_ecrites;

}

/*----------------------------------------------------------------------------
 * Lecture parallele de blocs de donnees
 *
 * parameters :
 *   PDM_file_par     <-- Pointeur sur le fichier
 *   taille_donnee   <-- Taille unitaire de la donnnee
 *   n_donnees_bloc  <-- Nombre de donnees a lire dans le bloc
 *   donnees         --> Donnees lues
 *   debut_bloc      <-- Debut du bloc dans l'ensemble des donnees
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
 )
{
  int n_donnees_lues = 0;
  int errcode = PDM_MPI_SUCCESS;
  int _taille_donnee = (int) taille_donnee;

  int n_octet = (int) (n_donnees_bloc * _taille_donnee); /* Passage en int
                                                           pour appel MPI */

  long debut_bloc_octet = debut_bloc * _taille_donnee;

  long n_octet_total = debut_bloc_octet + n_octet;

  /* Variables pour acces eo */

  PDM_MPI_Offset position;

  /* Variables pour acces ip */

  PDM_MPI_Datatype fichier_type;
  PDM_MPI_Aint positions[1];
  int lengths[1];
  int n_octet_lus;

  /* Lecture */

  switch(PDM_file_par->acces) {

  case FICHIER_PAR_ACCES_EO :
    position = PDM_file_par->offset + debut_bloc_octet;
    errcode = PDM_MPI_File_read_at_all(PDM_file_par->fichier,
                                   position,
                                   donnees,
                                   n_octet,
                                   PDM_MPI_BYTE,
                                   &n_octet_lus);
    break;

  case FICHIER_PAR_ACCES_IP :
    lengths[0] = n_octet;
    positions[0] = debut_bloc_octet;
    PDM_MPI_Type_create_hindexed(1, lengths, positions, PDM_MPI_BYTE, &fichier_type);
    PDM_MPI_Type_commit(&fichier_type);
    PDM_MPI_File_set_view(PDM_file_par->fichier,
                      PDM_file_par->offset,
                      PDM_MPI_BYTE,
                      fichier_type,
                      "native");
    errcode = PDM_MPI_File_read_all(PDM_file_par->fichier,
                                donnees,
                                n_octet,
                                PDM_MPI_BYTE,
                                &n_octet_lus);
    PDM_MPI_Type_free(&fichier_type);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0,"Erreur Fichier_par_lecture_parallele : type d'acces inconnu\n");
    exit(EXIT_FAILURE);

  }

  /* Traitement des erreurs de lecture */

  if (errcode != PDM_MPI_SUCCESS)
    _pdm_mpi_io_error_message(PDM_file_par->nom, errcode);

  /* Calcul du nombre de donnees lues */

/*   PDM_MPI_Get_count(&status, PDM_MPI_BYTE, &n_octet_lus); */

  n_donnees_lues = n_octet_lus / _taille_donnee;

  /* Mise à jour de l'offset */

  PDM_MPI_Bcast(&n_octet_total,
            1,
            PDM_MPI_LONG,
            PDM_file_par->n_rangs - 1,
            PDM_file_par->comm);

  PDM_file_par->offset += n_octet_total;

  return n_donnees_lues;
}

/*----------------------------------------------------------------------------
 * Ecriture parallele de blocs de donnees tries
 *
 * parameters :
 *   PDM_file_par       <-- Pointeur sur le fichier
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
 )
{
  int  n_donnees_ecrites = 0;
  int errcode = PDM_MPI_SUCCESS;
  int _taille_donnee = (int) taille_donnee;

  int n_octet = (int) (n_donnees_bloc * _taille_donnee); /* Passage en int
                                                           pour appel MPI */
  long debut_bloc_octet = debut_bloc * _taille_donnee;

  long n_octet_total = debut_bloc_octet + n_octet;

  /* Variables pour acces eo */

  PDM_MPI_Offset position;

  /* Variables pour acces ip */

  PDM_MPI_Datatype fichier_type;
  PDM_MPI_Aint positions[1];
  int lengths[1];
  int n_octet_ecrits;

  /* Ecriture */

  switch(PDM_file_par->acces) {

  case FICHIER_PAR_ACCES_EO :
    position = PDM_file_par->offset + debut_bloc_octet;
    errcode = PDM_MPI_File_write_at_all(PDM_file_par->fichier,
                                    position,
                                    donnees,
                                    n_octet,
                                    PDM_MPI_BYTE,
                                    &n_octet_ecrits);
    break;

  case FICHIER_PAR_ACCES_IP :
    lengths[0] = n_octet;
    positions[0] = debut_bloc_octet;
    PDM_MPI_Type_create_hindexed(1, lengths, positions, PDM_MPI_BYTE, &fichier_type);
    PDM_MPI_Type_commit(&fichier_type);
    PDM_MPI_File_set_view(PDM_file_par->fichier,
                      PDM_file_par->offset,
                      PDM_MPI_BYTE,
                      fichier_type,
                      "native");
    errcode = PDM_MPI_File_write_all(PDM_file_par->fichier,
                                 donnees,
                                 n_octet,
                                 PDM_MPI_BYTE,
                                 &n_octet_ecrits);
    PDM_MPI_Type_free(&fichier_type);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0,"Erreur Fichier_par_ecriture_parallele : type d'acces inconnu\n");
    exit(EXIT_FAILURE);

  }

  /* Traitement des erreurs  */

  if (errcode != PDM_MPI_SUCCESS)
    _pdm_mpi_io_error_message(PDM_file_par->nom, errcode);

  /* Calcul du nombre de donnees lues */

/*   PDM_MPI_Get_count(&status, PDM_MPI_BYTE, &n_octet_ecrits); */

  n_donnees_ecrites = n_octet_ecrits / _taille_donnee;

  /* Mise à jour de l'offset */

  PDM_MPI_Bcast(&n_octet_total,
            1,
            PDM_MPI_LONG,
            PDM_file_par->n_rangs - 1,
            PDM_file_par->comm);

  PDM_file_par->offset += n_octet_total;

  return n_donnees_ecrites;
}

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
(PDM_file_par_t     *PDM_file_par,
 PDM_MPI_Offset         offset,
 PDM_file_par_seek_t whence)
{
  int errcode = PDM_MPI_SUCCESS;

  switch(whence) {
  case FICHIER_PAR_SEEK_SET:
    PDM_file_par->offset = offset;
    break;
  case FICHIER_PAR_SEEK_CUR:
    PDM_file_par->offset += offset;
    break;
  case FICHIER_PAR_SEEK_END:
    {
      PDM_MPI_Offset f_taille = 0;
      errcode = PDM_MPI_File_get_size(PDM_file_par->fichier, &f_taille);
      PDM_file_par->offset = f_taille + offset;
    }
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_file_par_seek : whence non reconnu\n");
    exit(EXIT_FAILURE);
  }

  /* Definition de l'offset uniquement dans le cas IP
     en EO, l'offset est toujours donne explicitement
     a chaque acces par fichier->offset */

  if (PDM_file_par->acces == FICHIER_PAR_ACCES_IP) {
    int _whence;
    switch(whence) {
    case FICHIER_PAR_SEEK_SET:
      _whence = PDM_MPI_SEEK_SET;
      break;
    case FICHIER_PAR_SEEK_CUR:
      _whence = PDM_MPI_SEEK_CUR;
      break;
    case FICHIER_PAR_SEEK_END:
      _whence = PDM_MPI_SEEK_END;
      break;
    default :
      PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_file_par_seek : whence non reconnu\n");
      exit(EXIT_FAILURE);
    }

    errcode = PDM_MPI_File_seek(PDM_file_par->fichier,
                            offset,
                            _whence);
  }

  if (errcode != PDM_MPI_SUCCESS)
    _pdm_mpi_io_error_message(PDM_file_par->nom, errcode);
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

PDM_MPI_Offset PDM_file_par_tell(PDM_file_par_t *PDM_file_par)
{

  /* Theoriquement en mode IP, il faut recuperer la position par
     PDM_MPI_File_get_position + MPI_File_get_byte_offset + MPI_Allreduce max
     Mais la valeur est deja par PDM_file_par->offset qui est mise a jour
     a chaque acces en lecture ou ecriture  */

  return PDM_file_par->offset;

}

/*----------------------------------------------------------------------------
 * Fermeture du fchier
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void
PDM_file_par_close
(PDM_file_par_t       *PDM_file_par
)
{

  if (PDM_file_par != NULL) {
    free(PDM_file_par->nom);

    PDM_MPI_File_close(&(PDM_file_par->fichier));

  }

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
