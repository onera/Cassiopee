/*============================================================================
 * Traitement des entrees/sorties binaires sequentielles et paralleles
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_io_tab.h"
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
 * Definition des macros locales
 *============================================================================*/

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure partition locale
 *----------------------------------------------------------------------------*/

typedef struct {

  PDM_l_num_t         n_donnees;        /* Nombre de données de la partition */
  const PDM_l_num_t   *n_composantes;   /* Nombre variable de composantes */
  const PDM_g_num_t  *indirection;     /* Numérotation absolue pour
                                              données entrelacées */
  PDM_g_num_t  debut_bloc;             /* Adresse de début de bloc pour
                                              données rangées par blocs */
  const void             *donnees;         /* Données */

} PDM_io_partition_locale_t ;

/*----------------------------------------------------------------------------
 * Structure tableau
 *----------------------------------------------------------------------------*/

struct _PDM_io_tab_t {

  PDM_l_num_t               taille_donnee;         /* Taille unitaire de
                                                         la donnée */
  PDM_io_n_composantes_t     t_n_composantes;       /* Nombre de composante
                                                         variable ou cst */
  PDM_l_num_t               n_composantes_cst;     /* Nombre constant
                                                         de composantes */
  PDM_l_num_t               num_indirection_cedre; /* Numéro d'indirection
                                                         CEDRE */
  PDM_io_partition_locale_t **partitions_locales;   /* Description des
                                                         partitions locales */
};

/*============================================================================
 * Variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Numero max de variable CEDRE
 *----------------------------------------------------------------------------*/

static int _num_var_cedre_max = -1;

/*----------------------------------------------------------------------------
 * Unité du fichier
 *----------------------------------------------------------------------------*/

static PDM_l_num_t  _unite = -1;

/*----------------------------------------------------------------------------
 * Type de rangement
 *----------------------------------------------------------------------------*/

static PDM_io_rangement_t _t_rangement;

/*----------------------------------------------------------------------------
 * Nombre de partitions locales
 *----------------------------------------------------------------------------*/

static PDM_l_num_t _n_partition_local;

/*----------------------------------------------------------------------------
 * Stockage des structures indirections
 *----------------------------------------------------------------------------*/

static PDM_io_tab_t **PDM_io_tabs = NULL;


/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ajoute une partie des donnees dans un tableau associés à une variable
 * CEDRE
 *
 * arguments :
 *   num_var_cedre         <-- Numéro de variable CEDRE
 *   i_part                <-- indice de partition
 *   n_composantes         <-- Nombre de composantes pour chaque donnee
 *   n_donnees             <-- Nombre de donnees a lire
 *   indirection           <-- Indirection de redistribition des donnees
 *   donnees               <-- Donnees a écrire
 *
 *----------------------------------------------------------------------------*/

static void _ajout_donnees
(const PDM_l_num_t           num_var_cedre,
 const PDM_l_num_t           i_part,
 const PDM_l_num_t          *n_composantes,
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t         *indirection,
 void                          *donnees
 )
{
  if (PDM_io_tabs == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur _ajout_donnees :"
            "La phase d'acces n'a pas été initialisée\n");
    exit(EXIT_FAILURE);
  }

  int _num_var_cedre = num_var_cedre - 1;
  PDM_io_tab_t *tab = PDM_io_tabs[_num_var_cedre];

  /* Au premier appel : création du tableau lié à la variable CEDRE num_var_cedre */

  if (tab == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur _ajout_donnees :"
            "La variable n'a pas ete initialisee\n");
    exit(EXIT_FAILURE);
  }

  int _i_part = i_part - 1;

  if (tab->partitions_locales[_i_part] != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur _ajout_donnees : "
            "partition '%i' déjà traitée\n", i_part);
    exit(EXIT_FAILURE);
  }

  /* Mise à jour de la structure par prise en compte de la partition */

  tab->partitions_locales[_i_part] = (PDM_io_partition_locale_t *)
    malloc(sizeof(PDM_io_partition_locale_t));

  PDM_io_partition_locale_t *partition = tab->partitions_locales[_i_part];

  partition->n_donnees     = n_donnees;
  if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT)
    partition->n_composantes = NULL;
  else
    partition->n_composantes = n_composantes;
  if (_t_rangement == PDM_IO_RANGEMENT_BLOC) {
    partition->indirection = NULL;
    partition->debut_bloc  = indirection[0];
  }
  else
    partition->indirection   = indirection;
  partition->donnees       = donnees;
}

/*----------------------------------------------------------------------------
 * Definition d'une variable
 *
 *  arguments :
 *   num_var_cedre         <-- Numéro de variable CEDRE
 *   num_indirection_cedre <-- Numéro d'indirection CEDRE
 *   t_n_composantes       <-- Type de tailles composantes
 *                             (PDM_IO_N_COMPOSANTE_CONSTANT
 *                           ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes         <-- Nombre de composantes pour chaque donnee
 *   taille_donnee         <-- Taille unitaire de la donnnee
 *
 *----------------------------------------------------------------------------*/

static void _def_var
(const PDM_l_num_t           num_var_cedre,
 const PDM_l_num_t           num_indirection_cedre,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t           n_composantes,
 const PDM_l_num_t           taille_donnee
 )
{
  if (PDM_io_tabs == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur _def_var : "
            "La phase d'acces n'a pas été initialisée\n");
    exit(EXIT_FAILURE);
  }

  int _num_var_cedre = num_var_cedre - 1;
  PDM_io_tab_t *tab = PDM_io_tabs[_num_var_cedre];

  /* Au premier appel : création du tableau lié à la variable CEDRE num_var_cedre */

  if (tab == NULL) {
    PDM_io_tabs[_num_var_cedre] = (PDM_io_tab_t *) malloc(sizeof(PDM_io_tab_t));
    tab = PDM_io_tabs[_num_var_cedre];
    tab->taille_donnee = taille_donnee;
    tab->t_n_composantes = t_n_composantes;
    if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT)
      tab->n_composantes_cst = n_composantes;
    tab->t_n_composantes = t_n_composantes;
    tab->num_indirection_cedre = num_indirection_cedre;
    tab->partitions_locales = (PDM_io_partition_locale_t **)
      malloc(_n_partition_local * sizeof(PDM_io_partition_locale_t *));
    for (int i = 0; i < _n_partition_local; i++)
      tab->partitions_locales[i] = NULL;
  }
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialise une phase d'écriture parallèle de tableaux de données associées
 * aux numéros de variables CEDRE
 * Chaque tableau a ses propres caractéristiques :
 *         - taille de données
 *         - nombre de donnée
 *         - indirection (numérotation absolue)
 *
 * arguments :
 *   unite             <-- Unite du fichier
 *   rangement         <-- Type de rangement
 *   num_var_cedre_max <-- Numéro max de variable CEDRE
 *   n_partition_local <-- Nombre de partitions locales
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_ecr_debut, PDM_IO_TAB_ECR_DEBUT)
(const PDM_l_num_t *unite,
 const int            *t_rangement,
 const PDM_l_num_t *num_var_cedre_max,
 const PDM_l_num_t *n_partition_local
)
{
  PDM_io_rangement_t __t_rangement = PDM_IO_RANGEMENT_BLOC;

  if (*t_rangement == 0)
    __t_rangement = PDM_IO_RANGEMENT_BLOC;
  else if (*t_rangement == 1)
    __t_rangement = PDM_IO_RANGEMENT_ENTRELACE;

  PDM_io_tab_ecr_debut(*unite,
                         __t_rangement,
                         *num_var_cedre_max,
                         *n_partition_local);
}

void  PDM_io_tab_ecr_debut
(const PDM_l_num_t       unite,
 const PDM_io_rangement_t t_rangement,
 const PDM_l_num_t       num_var_cedre_max,
 const PDM_l_num_t       n_partition_local
)
{
  if (PDM_io_tabs == NULL) {
    _num_var_cedre_max = num_var_cedre_max;
    PDM_io_tabs = (PDM_io_tab_t **) malloc(_num_var_cedre_max *
                                              sizeof(PDM_io_tab_t *));

    for (int i = 0; i < _num_var_cedre_max; i++)
      PDM_io_tabs[i] = NULL;

    _unite = unite;
    _t_rangement = t_rangement;
    _n_partition_local = n_partition_local;
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_tab_ecr_debut :"
            "Une phase de lecture ou d'écriture est déjà en cours\n");
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------
 * Ajoute une partie des donnees dans un tableau associés à une variable
 * CEDRE
 *
 * arguments :
 *   num_var_cedre         <-- Numéro de variable CEDRE
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
 )
{
  PDM_io_tab_ecr_ajout_donnees(*num_var_cedre,
                                 *i_part,
                                 n_composantes,
                                 *n_donnees,
                                 indirection,
                                 donnees);
}

void PDM_io_tab_ecr_ajout_donnees
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            i_part,
 const PDM_l_num_t           *n_composantes,
 const PDM_l_num_t            n_donnees,
 const PDM_g_num_t          *indirection,
 void                           *donnees
 )
{
  _ajout_donnees(num_var_cedre,
                 i_part,
                 n_composantes,
                 n_donnees,
                 indirection,
                 donnees);
}

/*----------------------------------------------------------------------------
 * Definition d'une variable en ecriture
 *
 *  arguments :
 *   num_var_cedre         <-- Numéro de variable CEDRE
 *   num_indirection_cedre <-- Numéro d'indirection CEDRE
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
 )
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0)
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1)
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;

  PDM_io_tab_ecr_def_var(*num_var_cedre,
                           *num_indirection_cedre,
                           _t_n_composantes,
                           *n_composantes,
                           *taille_donnee);
}

void PDM_io_tab_ecr_def_var
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            num_indirection_cedre,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            n_composantes,
 const PDM_l_num_t            taille_donnee
 )
{
  _def_var(num_var_cedre,
           num_indirection_cedre,
           t_n_composantes,
           n_composantes,
           taille_donnee);
}


/*----------------------------------------------------------------------------
 * Finalise une phase d'écriture parallèle de tableaux de données associées
 * aux numéros de variables CEDRE. Cette fonction déclenche réellement
 * les écritures
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_ecr_fin, PDM_IO_TAB_ECR_FIN)
(void)
{
  PDM_io_tab_ecr_fin();
}

void PDM_io_tab_ecr_fin
(void)
{
  if (PDM_io_tabs == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_tab_ecr_fin :"
            "La phase d'écriture n'a pas été initialisée\n");
    exit(EXIT_FAILURE);
  }

  /* Definition des tableaux temporaires de concatenation */

  PDM_l_num_t  *n_composantes_concatene     = NULL;
  int t_n_composantes_concatene = 0;

  PDM_g_num_t *indirection_concatene       = NULL;
  int t_indirection_concatene = 0;

  unsigned char   *donnees_concatene           = NULL;
  int t_donnees_concatene = 0;

  /* Definition des pointeurs des tableaux utiles à l'écriture */

  const PDM_l_num_t  *n_composantes_ecrit   = NULL;
  const PDM_g_num_t *indirection_ecrit     = NULL;
  const unsigned char   *donnees_ecrit         = NULL;

  /* Boucle sur les variables */

  for (int i = 0; i < _num_var_cedre_max; i++) {

    if (PDM_io_tabs[i] != NULL) {

      PDM_io_tab_t *tab = PDM_io_tabs[i];
      PDM_io_partition_locale_t **partitions = tab->partitions_locales;

      /* Concaténation */

      int n_donnees_total = 0;
      int n_composantes_total = 0;

      if (_n_partition_local > 1) {

        /* Calcul des grandeurs concaténées */

        for (int j = 0; j < _n_partition_local; j++) {
          PDM_io_partition_locale_t *partition = partitions[j];
          if (partition != NULL) {

            n_donnees_total += partition->n_donnees;

            if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
              n_composantes_total += tab->n_composantes_cst * partition->n_donnees;
            }
            else {
              for (int k = 0; k < partition->n_donnees; k++) {
                n_composantes_total += partition->n_composantes[k];
              }
            }
          }
        }

        /* Allocation ou reallocation du tableau n_composantes_concatene */

        if (n_composantes_concatene == NULL &&
            tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

          t_n_composantes_concatene = n_donnees_total;
          n_composantes_concatene = (PDM_l_num_t *)
            malloc(t_n_composantes_concatene * sizeof(PDM_l_num_t));

        }

        else if (n_composantes_concatene != NULL &&
                 tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

          if (n_donnees_total > t_n_composantes_concatene) {
            t_n_composantes_concatene = n_donnees_total;
            n_composantes_concatene = (PDM_l_num_t *)
              realloc((void *) n_composantes_concatene,
                      t_n_composantes_concatene * sizeof(PDM_l_num_t));
          }
        }

        /* Allocation ou reallocation du tabeau indirection_concatene */

        if (indirection_concatene == NULL &&
            _t_rangement == PDM_IO_RANGEMENT_ENTRELACE) {

          t_indirection_concatene = n_donnees_total;
          indirection_concatene = (PDM_g_num_t *)
            malloc(t_indirection_concatene * sizeof(PDM_g_num_t));

        }

        else if (indirection_concatene != NULL &&
                 _t_rangement == PDM_IO_RANGEMENT_ENTRELACE) {

          if (n_donnees_total > t_indirection_concatene) {
            t_indirection_concatene = n_donnees_total;
            indirection_concatene = (PDM_g_num_t *)
              realloc((void *)indirection_concatene,
                      t_indirection_concatene * sizeof(PDM_g_num_t));
          }

        }

        /* Allocation ou reallocation du tabeau donnees_concatene */

        if (donnees_concatene == NULL) {

          t_donnees_concatene = n_composantes_total * tab->taille_donnee;
          donnees_concatene = (unsigned char*)
            malloc(t_donnees_concatene * sizeof(unsigned char));
        }
        else {
          if (t_donnees_concatene < n_composantes_total * tab->taille_donnee) {
            t_donnees_concatene = n_composantes_total * tab->taille_donnee;
            donnees_concatene = (unsigned char*)
              realloc((void *) donnees_concatene,
                      t_donnees_concatene * sizeof(unsigned char));
          }
        }

        /* Traitement du cas PDM_IO_RANGEMENT_BLOC */

        if (_t_rangement == PDM_IO_RANGEMENT_BLOC) {
          for (int j = 0; j < _n_partition_local; j++) {
            if (partitions[j] != NULL) {
              indirection_ecrit = &(partitions[j]->debut_bloc);
              break;
            }
          }
        }
        else
          indirection_ecrit = indirection_concatene;

        /* Traitement du cas PDM_IO_N_COMPOSANTE_CONSTANT */

        if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
          n_composantes_ecrit = &(tab->n_composantes_cst);
        }
        else
          n_composantes_ecrit = n_composantes_concatene;

        donnees_ecrit = donnees_concatene;

        /* Copie */

        int k1 = 0;
        int k2 = 0;

        for (int j = 0; j < _n_partition_local; j++) {

          PDM_io_partition_locale_t *partition = partitions[j];

          if (partition != NULL) {

            int l1 = 0;
            unsigned char *_donnees = (unsigned char *) partition->donnees;

            for (int k = 0; k < partition->n_donnees; k++) {
              if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE)
                n_composantes_concatene[k1] = partition->n_composantes[k];

              if (_t_rangement == PDM_IO_RANGEMENT_ENTRELACE)
                indirection_concatene[k1] = partition->indirection[k];

              k1 += 1;

              if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
                for (int l = 0; l < (tab->taille_donnee * partition->n_composantes[k]); l++) {
                  donnees_concatene[k2] = _donnees[l1];
                  l1 += 1;
                  k2 += 1;
                }
              }
              else {
                for (int l = 0; l < (tab->taille_donnee * tab->n_composantes_cst); l++) {
                  donnees_concatene[k2] = _donnees[l1];
                  l1 += 1;
                  k2 += 1;
                }
              }
            }

            free(partition);

          }

        }

      }

      else if (_n_partition_local == 1) {

        if (partitions[0] != NULL) {
          n_donnees_total     = partitions[0]->n_donnees;
          if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
            n_composantes_ecrit = &(tab->n_composantes_cst);
          }
          else
            n_composantes_ecrit = partitions[0]->n_composantes;
          indirection_ecrit   = partitions[0]->indirection;
          donnees_ecrit       = (unsigned char*) partitions[0]->donnees;
        }

        else {
          n_donnees_total     = 0;
          if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
            n_composantes_ecrit = &(tab->n_composantes_cst);
          }
          else
            n_composantes_ecrit = NULL;
          indirection_ecrit   = NULL;
          donnees_ecrit       = NULL;
        }

      }

      else {

        n_donnees_total     = 0;
        if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
          n_composantes_ecrit = &(tab->n_composantes_cst);
        }
        else
          n_composantes_ecrit = NULL;
        indirection_ecrit   = NULL;
        donnees_ecrit       = NULL;
      }

     /* Ecriture */

      if (_t_rangement == PDM_IO_RANGEMENT_BLOC) {

        PDM_g_num_t pt_bloc;
        if (indirection_ecrit == NULL)
          pt_bloc = 0;
        else
          pt_bloc = indirection_ecrit[0];

        PDM_io_ecr_par_bloc(_unite,
                              tab->t_n_composantes,
                              n_composantes_ecrit,
                              tab->taille_donnee,
                              n_donnees_total,
                              pt_bloc,
                              (const void *) donnees_ecrit);
      }
      else if (_t_rangement == PDM_IO_RANGEMENT_ENTRELACE) {
         PDM_io_ecr_par_entrelacee(_unite,
                                    tab->t_n_composantes,
                                    n_composantes_ecrit,
                                    tab->taille_donnee,
                                    n_donnees_total,
                                    indirection_ecrit,
                                    (const void *) donnees_ecrit);
      }
      if (_n_partition_local == 1)
        free(partitions[0]);

      free(tab->partitions_locales);
      free(tab);
    }
  }

  /* Libération mémoire de la structure  */

  free(PDM_io_tabs);
  PDM_io_tabs = NULL;

  /* Libération mémoire des tableaux de concatenation  */

  if (n_composantes_concatene != NULL)
    free(n_composantes_concatene);
  if (indirection_concatene != NULL)
    free(indirection_concatene);
  if (donnees_concatene != NULL)
    free(donnees_concatene);

}


/*----------------------------------------------------------------------------
 * Initialise une phase de lecture parallèle de tableaux de données associées
 * aux numéros de variables CEDRE
 * Chaque tableau a ses propres caractéristiques :
 *         - taille de données
 *         - nombre de donnée
 *         - indirection (numérotation absolue)
 *
 * arguments :
 *   unite             <-- Unite du fichier
 *   t_rangement       <-- Type de rangement
 *   num_var_cedre_max <-- Numéro max de variable CEDRE
 *   n_partition_local <-- Nombre de partitions locales
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_lec_debut, PDM_IO_TAB_LEC_DEBUT)
(const PDM_l_num_t *unite,
 const int            *t_rangement,
 const PDM_l_num_t *num_var_cedre_max,
 const PDM_l_num_t *n_partition_local
)
{
  PDM_io_rangement_t __t_rangement = PDM_IO_RANGEMENT_BLOC;

  if (*t_rangement == 0)
    __t_rangement = PDM_IO_RANGEMENT_BLOC;
  else if (*t_rangement == 1)
    __t_rangement = PDM_IO_RANGEMENT_ENTRELACE;

  PDM_io_tab_lec_debut(*unite,
                         __t_rangement,
                         *num_var_cedre_max,
                         *n_partition_local);
}

void PDM_io_tab_lec_debut
(const PDM_l_num_t       unite,
 const PDM_io_rangement_t t_rangement,
 const PDM_l_num_t       num_var_cedre_max,
 const PDM_l_num_t       n_partition_local
)
{
  if (PDM_io_tabs == NULL) {
    _num_var_cedre_max = num_var_cedre_max;
    PDM_io_tabs = (PDM_io_tab_t **) malloc(_num_var_cedre_max *
                                               sizeof(PDM_io_tab_t *));

    for (int i = 0; i < _num_var_cedre_max; i++)
      PDM_io_tabs[i] = NULL;

    _unite = unite;
    _t_rangement = t_rangement;
    _n_partition_local = n_partition_local;
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_tab_lec_debut :"
            "Une phase de lecture ou d'écriture est déjà en cours\n");
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------
 * Ajoute une partie des donnees dans un tableau associés à une variable
 * CEDRE
 *
 * arguments :
 *   num_var_cedre         <-- Numéro de variable CEDRE
 *   num_indirection_cedre <-- Numéro d'indirection CEDRE
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
 )
{

  PDM_io_tab_lec_ajout_donnees(*num_var_cedre,
                                 *i_part,
                                 n_composantes,
                                 *n_donnees,
                                 indirection,
                                 donnees);
}

void PDM_io_tab_lec_ajout_donnees
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            i_part,
 const PDM_l_num_t           *n_composantes,
 const PDM_l_num_t            n_donnees,
 const PDM_g_num_t          *indirection,
 void                           *donnees
 )
{
  _ajout_donnees(num_var_cedre,
                 i_part,
                 n_composantes,
                 n_donnees,
                 indirection,
                 donnees);
}

/*----------------------------------------------------------------------------
 * Definition d'une variable en ecriture
 *
 *  arguments :
 *   num_var_cedre         <-- Numéro de variable CEDRE
 *   num_indirection_cedre <-- Numéro d'indirection CEDRE
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
 )
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0)
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1)
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;

  PDM_io_tab_lec_def_var(*num_var_cedre,
                           *num_indirection_cedre,
                           _t_n_composantes,
                           *n_composantes,
                           *taille_donnee);
}

void PDM_io_tab_lec_def_var
(const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            num_indirection_cedre,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            n_composantes,
 const PDM_l_num_t            taille_donnee
 )
{
  _def_var(num_var_cedre,
           num_indirection_cedre,
           t_n_composantes,
           n_composantes,
           taille_donnee);
}

/*----------------------------------------------------------------------------
 * Finalise une phase de lecture parallèle de tableaux de données associées
 * aux numéros de variables CEDRE. Cette fonction déclenche réellement
 * les écritures
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tab_lec_fin, PDM_IO_TAB_LEC_FIN)
(void)
{
  PDM_io_tab_lec_fin();
}

void PDM_io_tab_lec_fin
(void)
{
  if (PDM_io_tabs == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_tab_lec_fin :"
            "La phase de lecture n'a pas été initialisée\n");
    exit(EXIT_FAILURE);
  }

  /* Definition des tableaux temporaires de concatenation */

  PDM_l_num_t  *n_composantes_concatene     = NULL;
  int t_n_composantes_concatene = 0;

  PDM_g_num_t *indirection_concatene       = NULL;
  int t_indirection_concatene = 0;

  unsigned char   *donnees_concatene           = NULL;
  int t_donnees_concatene = 0;

  /* Definition des pointeurs des tableaux utiles à l'écriture */

  const PDM_l_num_t  *n_composantes_lu   = NULL;
  const PDM_g_num_t *indirection_lu     = NULL;
  const unsigned char   *donnees_lu         = NULL;

  /* Boucle sur les variables */

  for (int i = 0; i < _num_var_cedre_max; i++) {

    if (PDM_io_tabs[i] != NULL) {

      PDM_io_tab_t *tab = PDM_io_tabs[i];
      PDM_io_partition_locale_t **partitions = tab->partitions_locales;

      /* Concaténation */

      int n_donnees_total = 0;
      int n_composantes_total = 0;

      if (_n_partition_local > 1) {

        /* Calcul des grandeurs concaténées */

        for (int j = 0; j < _n_partition_local; j++) {
          PDM_io_partition_locale_t *partition = partitions[j];
          if (partition != NULL) {

            n_donnees_total += partition->n_donnees;

            if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
              n_composantes_total += tab->n_composantes_cst * partition->n_donnees;
            }
            else {
              for (int k = 0; k < partition->n_donnees; k++) {
                n_composantes_total += partition->n_composantes[k];
              }
            }
          }
        }

        /* Allocation ou reallocation du tableau n_composantes_concatene */

        if (n_composantes_concatene == NULL &&
            tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

          t_n_composantes_concatene = n_donnees_total;
          n_composantes_concatene = (PDM_l_num_t *)
            malloc(t_n_composantes_concatene * sizeof(PDM_l_num_t));

        }

        else if (n_composantes_concatene != NULL &&
                 tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

          if (n_donnees_total > t_n_composantes_concatene) {
            t_n_composantes_concatene = n_donnees_total;
            n_composantes_concatene = (PDM_l_num_t *)
              realloc((void *) n_composantes_concatene,
                      t_n_composantes_concatene * sizeof(PDM_l_num_t));
          }
        }

        /* Allocation ou reallocation du tabeau indirection_concatene */

        if (indirection_concatene == NULL &&
            _t_rangement == PDM_IO_RANGEMENT_ENTRELACE) {

          t_indirection_concatene = n_donnees_total;
          indirection_concatene = (PDM_g_num_t *)
            malloc(t_indirection_concatene * sizeof(PDM_g_num_t));

        }

        else if (indirection_concatene != NULL &&
                 _t_rangement == PDM_IO_RANGEMENT_ENTRELACE) {

          if (n_donnees_total > t_indirection_concatene) {
            t_indirection_concatene = n_donnees_total;
            indirection_concatene = (PDM_g_num_t *)
              realloc((void *)indirection_concatene,
                      t_indirection_concatene * sizeof(PDM_g_num_t));
          }

        }

        /* Allocation ou reallocation du tabeau donnees_concatene */

        if (donnees_concatene == NULL) {

          t_donnees_concatene = n_composantes_total * tab->taille_donnee;
          donnees_concatene = (unsigned char*)
            malloc(t_donnees_concatene * sizeof(unsigned char));
        }
        else {
          if (t_donnees_concatene < n_composantes_total * tab->taille_donnee) {
            t_donnees_concatene = n_composantes_total * tab->taille_donnee;
            donnees_concatene = (unsigned char*)
              realloc((void *) donnees_concatene,
                      t_donnees_concatene * sizeof(unsigned char));
          }
        }

        /* Traitement du cas PDM_IO_RANGEMENT_BLOC */

        if (_t_rangement == PDM_IO_RANGEMENT_BLOC) {
          for (int j = 0; j < _n_partition_local; j++) {
            if (partitions[j] != NULL) {
              indirection_lu = &(partitions[j]->debut_bloc);
              break;
            }
          }
        }
        else
          indirection_lu = indirection_concatene;

        /* Traitement du cas PDM_IO_N_COMPOSANTE_CONSTANT */

        if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
          n_composantes_lu = &(tab->n_composantes_cst);
        }
        else
          n_composantes_lu = n_composantes_concatene;

        donnees_lu = donnees_concatene;

        /* Copie */

        int k1 = 0;

        for (int j = 0; j < _n_partition_local; j++) {

          PDM_io_partition_locale_t *partition = partitions[j];

          if (partition != NULL) {

             for (int k = 0; k < partition->n_donnees; k++) {
              if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE)
                n_composantes_concatene[k1] = partition->n_composantes[k];

              if (_t_rangement == PDM_IO_RANGEMENT_ENTRELACE)
                indirection_concatene[k1] = partition->indirection[k];

              k1 += 1;
            }
          }
        }
      }

      else if (_n_partition_local == 1) {
        if (partitions[0] != NULL) {
          n_donnees_total     = partitions[0]->n_donnees;
          if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
            n_composantes_lu = &(tab->n_composantes_cst);
          }
          else
            n_composantes_lu = partitions[0]->n_composantes;
          indirection_lu   = partitions[0]->indirection;
          donnees_lu       = (unsigned char*) partitions[0]->donnees;
        }
        else {

          n_donnees_total     = 0;
          if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
            n_composantes_lu = &(tab->n_composantes_cst);
          }
          else
            n_composantes_lu  = NULL;
          indirection_lu      = NULL;
          donnees_lu          = NULL;
        }
      }

      else {

        n_donnees_total     = 0;
        if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
          n_composantes_lu = &(tab->n_composantes_cst);
        }
        else
          n_composantes_lu  = NULL;
        indirection_lu      = NULL;
        donnees_lu          = NULL;
      }

      /* Lecture */

      if (_t_rangement == PDM_IO_RANGEMENT_BLOC) {

        PDM_g_num_t pt_bloc;
        if (indirection_lu == NULL)
          pt_bloc = 0;
        else
          pt_bloc = indirection_lu[0];

        PDM_io_lec_par_bloc(_unite,
                              tab->t_n_composantes,
                              n_composantes_lu,
                              tab->taille_donnee,
                              n_donnees_total,
                              pt_bloc,
                              (void *) donnees_lu);
      }
      else if (_t_rangement == PDM_IO_RANGEMENT_ENTRELACE)
        PDM_io_lec_par_entrelacee(_unite,
                                    tab->t_n_composantes,
                                    n_composantes_lu,
                                    tab->taille_donnee,
                                    n_donnees_total,
                                    indirection_lu,
                                    (void *) donnees_lu);

     /* Copie des données si nécessaire  */

      if (_n_partition_local > 1) {

        int k1 = 0;

        for (int j = 0; j < _n_partition_local; j++) {

          PDM_io_partition_locale_t *partition = partitions[j];

          if (partition != NULL) {

            int l1 = 0;
            unsigned char *_donnees = (unsigned char *) partition->donnees;

            for (int k = 0; k < partition->n_donnees; k++) {
              if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
                for (int l = 0; l < (tab->taille_donnee * partition->n_composantes[k]); l++) {
                  _donnees[l1] = donnees_concatene[k1];
                  l1 += 1;
                  k1 += 1;
                }
              }
              else if (tab->t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
                for (int l = 0; l < (tab->taille_donnee * tab->n_composantes_cst); l++) {
                  _donnees[l1] = donnees_concatene[k1];
                  l1 += 1;
                  k1 += 1;
                }
              }
            }

            free(partition);

          }
        }
      }

      /* Libération mémoire du tableau  */

      if (_n_partition_local == 1)
        free(partitions[0]);

      free(tab->partitions_locales);
      free(tab);

    }
  }

  /* Libération mémoire de la structure  */

  free(PDM_io_tabs);
  PDM_io_tabs = NULL;

  /* Libération mémoire des tableaux de concatenation  */

  if (n_composantes_concatene != NULL)
    free(n_composantes_concatene);
  if (indirection_concatene != NULL)
    free(indirection_concatene);
  if (donnees_concatene != NULL)
    free(donnees_concatene);

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
