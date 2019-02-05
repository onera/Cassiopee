/*============================================================================
 * Traitement des entrees/sorties binaires sequentielles et paralleles
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_io.h"
#include "pdm_file_seq.h"
#include "pdm_file_par.h"
#include "pdm_timer.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi_node_first_rank.h" 
#include "pdm_fortran_to_c_string.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"

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

#define PDM_IO_MAX(a,b) ((a) > (b) ? (a) : (b))
#define PDM_IO_MIN(a,b) ((a) < (b) ? (a) : (b))

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Type decrivant un fichier de type parallele io
 *----------------------------------------------------------------------------*/

struct _PDM_io_fichier_t {

  char               *nom;                /* Nom du fichier */
  PDM_io_mode_t     mode;               /* Mode du fichier */
  PDM_io_acces_t    acces;              /* Type d'entrees/sorties */

  int                 swap_endian;        /* Active le swap little/big endian */

  PDM_MPI_Comm            comm;               /* Communicateur lie 
                                             a cette structure */
  PDM_MPI_Comm            scomm;           /* Sub-communicator reduced to active ranks */
  int                 rang;               /* Rang MSG */
  int                 n_rangs;            /* Nombre de rangs MSG  */

  PDM_timer_t      *timer_fichier;      /* Mesure des temps d'acces 
                                             aux fichiers  */
  PDM_timer_t      *timer_swap_endian;  /* Mesure des temps de swap */
  PDM_timer_t      *timer_total;        /* Mesure des temps de swap */
  PDM_timer_t      *timer_distribution; /* Mesure des temps de distribution 
                                             des donnees */
  PDM_file_seq_t      *PDM_file_seq;        /* Fichier sequentiel */
  PDM_file_par_t      *PDM_file_par;        /* Fichier parallele */

  PDM_io_fmt_t      fmt_t;              /* Type de format */
  char               *fmt;                /* Format */
  PDM_l_num_t      n_char_fmt;         /* Nb de caractères du format */
  PDM_io_type_t     data_type;          /* Type de données  du format */
  PDM_io_backup_t   backup;             /* Backup du fichier en cas de reecriture */
  
  double      prop_noeuds_actifs;         /* Proportion de noeuds actifs */
  int                 n_rangs_actifs;     /* Nombre de rangs actifs */
  int                 n_rangs_inactifs;     /* Number of inactive ranks */
  int                *rangs_actifs;       /* Active ranks */
  int                *rangs_inactifs;       /* Inactive ranks */
  int                *tag_rangs_actifs;   /* Tag des rangs actifs */
  int                 rang_actif;         /* Indique si rang courant est actif */

};

/*============================================================================
 * Variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Stockage des objets PDM_io_fichiers
 *----------------------------------------------------------------------------*/

static PDM_Handles_t *PDM_io_fichiers = NULL; 

/*----------------------------------------------------------------------------
 * tag pour Echanges MPI
 *----------------------------------------------------------------------------*/

static const int PDM_io_tag = 'p'+'a'+'r'+'i'+'o'+'t'+'a'+'g';

/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Partitionnement pour le tri quick sort
 *
 * parameters :
 *   tableau          <-> tableau a trier
 *   p                <-- Indice de debut
 *   r                <-- Indice de fin
 * return
 *   PDM_io_version      Description version CEDRE      
 *----------------------------------------------------------------------------*/

static int
_mkdir
(
const char* path
 )
{
  if (mkdir(path, S_IRWXU|S_IRWXG|S_IRWXO) != 0) {

    if (errno == EEXIST) {
      struct stat buf;

      if (stat(path, &buf) != 0) {
        PDM_printf("  Fichier ou repertoire existant "
               "et son statut n'est pas valable\n");
        abort();
      }
      else if (S_ISDIR(buf.st_mode) != 1) {
        PDM_printf("  Fichier existant et ce n'est pas un repertoire\n");
        abort();
      }
      else
        return 0;

      errno = EEXIST; /* In case modified by stat() */

    }
    else {
      PDM_printf("  Fichier existant et ce n'est pas un repertoire\n");
      abort();
    }

    return -1;

  }

  return 0;

}

/*----------------------------------------------------------------------------
 * Detemination de liste des rangs actifs qui accedent reellement aux fichiers
 * et le nombre de donnees traitees par chaque rang (0 pour les rangs inactifs)
 *
 * parameters :
 *   fichier         <-- fichier traite
 *   n_donnnees      <-- nombre de donnees a traiter pour ce rang
 *   rang_actif      --> 1 si le rang courant est actif, 0 sinon
 *   n_rang_actif    --> Nombre de rangs actifs
 *   rangs_actifs    --> Liste des rangs actifs (Alloue dans la fonction)
 *
 *----------------------------------------------------------------------------*/

static void _rangs_actifs
(
PDM_io_fichier_t  *fichier
)
{

  int tout_rang_actif = 0;
  
  if (fichier->n_rangs_actifs == 0) {
  
    if (fichier->n_rangs == 1) {
      fichier->rang_actif  = 1; /* Par defaut le processus est inactif pour les acces 
                                 paralleles */
      fichier->n_rangs_actifs = 1; /* Nombre de processus actifs 
                                    pour les acces paralleles */
      fichier->n_rangs_inactifs = 0;
  
      fichier->rangs_actifs = (int *) malloc(sizeof(int) * fichier->n_rangs_actifs);
      fichier->rangs_inactifs = NULL;
      fichier->rangs_actifs[0] = 0;
    }    
 
    else {
      /* Intialisation */

      fichier->rang_actif  = 0; /* Par defaut le processus est inactif pour les acces 
				   paralleles */
      fichier->n_rangs_actifs = 0; /* Nombre de processus actifs 
				      pour les acces paralleles */
  
      fichier->n_rangs_inactifs = fichier->n_rangs;
      fichier->rangs_actifs = NULL;
      fichier->rangs_inactifs = NULL;

      /* Repartition sur les procs 0 de chaque noeud */
      /* On peut ameliorer le principe en ne prenant qu'un noeud sur 2, 4, 8, 16 */

      /* Determination du rang dans le noeud (si prop_noeud negatifs tous les coeurs sont actifs)*/

      if (fichier->prop_noeuds_actifs <= 0) {
      	fichier->rang_actif = 1;
      	tout_rang_actif = 1;
      }
      else {

      /* Determination des rangs 0 des noeuds */

        int rang_noeud = PDM_io_mpi_node_rank(fichier->comm);
        if (rang_noeud == 0)
          fichier->rang_actif = 1;
      }

      fichier->tag_rangs_actifs = (int *) malloc(sizeof(int) * fichier->n_rangs);
    
      PDM_MPI_Allgather((void *) &fichier->rang_actif, 1, PDM_MPI_INT, 
		    (void *) fichier->tag_rangs_actifs, 1, PDM_MPI_INT, 
		    fichier->comm);
    
      fichier->n_rangs_actifs = 0;
      for (int i = 0; i < fichier->n_rangs; i++) {
        if (fichier->tag_rangs_actifs[i] == 1) {
          fichier->n_rangs_actifs += 1;
        }
      }

      /* Prise en compte de la proportion de noeuds actifs */

      if ((tout_rang_actif == 0) && (fichier->prop_noeuds_actifs < 1.)) {
        int n_rang_actif_pondere =  (int) (fichier->n_rangs_actifs * fichier->prop_noeuds_actifs);
        n_rang_actif_pondere = PDM_IO_MAX(n_rang_actif_pondere, 1);
        int pas = fichier->n_rangs_actifs / n_rang_actif_pondere;
        pas = PDM_IO_MAX(pas, 1); 
      
        int cpt = 0;
        fichier->n_rangs_actifs = 0;
        for (int i = 0; i < fichier->n_rangs; i++) {
          if (fichier->tag_rangs_actifs[i] == 1) {
            if (fichier->n_rangs_actifs < n_rang_actif_pondere) {
              if (cpt != 0) {
                fichier->tag_rangs_actifs[i] = 0;
              }
              else {
                fichier->n_rangs_actifs += 1;
              }
              cpt = (cpt + 1) % pas;
            }
            else {
              fichier->tag_rangs_actifs[i] = 0;
            }
          }
        }
        if (fichier->tag_rangs_actifs[fichier->rang] == 0) {
          fichier->rang_actif = 0;
        }
      }
  
      /* Determination des rangs actifs */

      fichier->rangs_actifs = (int *) malloc(sizeof(int) * fichier->n_rangs_actifs);
      fichier->n_rangs_inactifs = fichier->n_rangs - fichier->n_rangs_actifs;
      fichier->rangs_inactifs = (int *) malloc(sizeof(int) * fichier->n_rangs_inactifs);

      fichier->n_rangs_actifs = 0;
      fichier->n_rangs_inactifs = 0;

      for (int i = 0; i < fichier->n_rangs; i++) {
        if (fichier->tag_rangs_actifs[i] == 1) {
          fichier->rangs_actifs[fichier->n_rangs_actifs++] = i;
        }
        else {
          fichier->rangs_inactifs[fichier->n_rangs_inactifs++] = i;
        }
      }

      /* Affichage */

      if (0 == 1) {
        if (fichier->rang == 0) {
          PDM_printf("rangs actifs : ");
          for(int i = 0; i < fichier->n_rangs_actifs; i++)
            PDM_printf(" %d", fichier->rangs_actifs[i]);
          PDM_printf("\n");
        }
      }
      
      PDM_MPI_Comm_split (fichier->comm,
                         fichier->rang_actif,
                         fichier->rang,
                         &(fichier->scomm));

    }
  }
}


/*----------------------------------------------------------------------------
 * Detemination de liste des rangs actifs qui accedent reellement aux fichiers
 * et le nombre de donnees traitees par chaque rang (0 pour les rangs inactifs)
 *
 * parameters :
 *   fichier         <-- fichier traite
 *   n_donnnees_total<-- nombre de donnees a traiter pour ce rang
 *   n_donnees_rangs --> Nombre de donnees traitees par chaque rang
 *                       (Taille en n_rangs, valeur nulle pour les rangs 
 *                        inactifs)
 *   n_donnees_rang_min --> nombre de donnees min pour l'ensemble des rangs
 *   n_donnees_rang_max --> nombre de donnees max pour l'ensemble des rangs
 *
 *----------------------------------------------------------------------------*/

static void _n_donnees_rang
(PDM_io_fichier_t  *fichier,
 const PDM_g_num_t n_donnees_total,
 PDM_g_num_t     *n_donnees_rangs,
 int                 *n_donnees_rang_min,
 int                 *n_donnees_rang_max
)
{

  /* Determination des données traitées par chaque rang */
  
  PDM_g_num_t _n_donnees_rang = n_donnees_total / fichier->n_rangs_actifs;
  int n_donnees_rang = (int) _n_donnees_rang;

  PDM_g_num_t _reste = n_donnees_total % fichier->n_rangs_actifs;
  int reste = (int) _reste;                    
  
  *n_donnees_rang_max = n_donnees_rang;
  *n_donnees_rang_min = n_donnees_rang;
  
  if (reste != 0)
    *n_donnees_rang_max = n_donnees_rang + 1;

  for (int i = 0; i < fichier->n_rangs + 1; i++) 
    n_donnees_rangs[i] = 0;

  int k = 0;
  for (int i = 0; i < fichier->n_rangs; i++) {
    n_donnees_rangs[i+1] += n_donnees_rangs[i];
    if (fichier->tag_rangs_actifs[i] == 1) {
      n_donnees_rangs[i+1] += n_donnees_rang;
      if (k < reste)
        n_donnees_rangs[i+1] += 1;
      k += 1;
    }
  }

  /* Affichage */

  if (0 == 1) {
    if (fichier->rang == 0) {
      PDM_printf("n_donnees_rangs : ");
      for(int i = 0; i < fichier->n_rangs + 1; i++)
        PDM_printf(PDM_FMT_G_NUM" ",n_donnees_rangs[i]);
      PDM_printf("\n");
    }
  }
}


/*----------------------------------------------------------------------------
 * Dertermine les parametres de distributions des donnees pour les ecritures
 * et lectures par blocs
 *
 * parameters :
 *   fichier             <-- fichier traite
 *   t_n_composantes     <-- Type de tailles composantes 
 *                             (PDM_IO_N_COMPOSANTE_CONSTANT
 *                           ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes       <-- Nombre de composantes pour chaque donnee
 *   debut_bloc          <-- Adresse de debut de bloc dans la numerotation 
 *                           absolue
 *   n_donnnees          <-- nombre de donnees a traiter pour ce rang
 *   rang_actif          <-- 1 si le rang courant est actif, 0 sinon
 *   n_donnees_rangs     <-- Nombre de donnees traitees par chaque rang
 *                          (Taille en n_rangs, valeur nulle pour les rangs 
 *                           inactifs)
 *   n_donnees_a_envoyer  --> Nombre de donnees a envoyer a chaque processus 
 *   i_donnees_a_envoyer  --> Index correspondant
 *   n_donnees_a_recevoir --> Nombre de donnees recues de chaque processus
 *   i_donnees_a_recevoir --> Index correspondant
 *
 *----------------------------------------------------------------------------*/

static void _calcul_parametres_distribution_bloc
(PDM_io_fichier_t             *fichier,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t           *n_composantes,         
 const PDM_g_num_t           debut_bloc,
 const PDM_l_num_t            n_donnees,
 const int                       rang_actif,
 const PDM_g_num_t          *n_donnees_traitees_rangs,
 int                            *n_donnees_a_envoyer,
 int                            *i_donnees_a_envoyer,
 int                            *n_donnees_a_recevoir,
 int                            *i_donnees_a_recevoir)

{
  /*------------------------------------------------------------
   * Repartition des donnees sur les processus actifs
   *------------------------------------------------------------ */
  
  for (int i = 0; i < fichier->n_rangs; i++) {
    n_donnees_a_envoyer[i]  = 0;
    n_donnees_a_recevoir[i] = 0;
  }
  
  /* Partage des debut de blocs et du nombre de donnees par bloc */
  
  PDM_g_num_t  *n_absolue_bloc_rangs  = 
    (PDM_g_num_t *) malloc((fichier->n_rangs + 1) * sizeof(PDM_g_num_t));
  n_absolue_bloc_rangs[0] = 1;
  
  PDM_g_num_t n_absolue_bloc = debut_bloc + n_donnees;
  
  PDM_MPI_Allgather(&n_absolue_bloc, 1, PDM__PDM_MPI_G_NUM, 
                n_absolue_bloc_rangs + 1, 1, PDM__PDM_MPI_G_NUM, 
                fichier->comm);
  
  /* Determination du nombre de donnees a envoyer a chaque processus 
     irang_min - irang_max Ã©tant la plage de repartition du bloc du 
     processus courant */
  
  PDM_g_num_t n_absolue_min = debut_bloc;
  PDM_g_num_t n_absolue_max = debut_bloc + n_donnees - 1;
  int irang = 0;
  
  while (n_absolue_min > n_donnees_traitees_rangs[irang + 1])
    irang++;
  
  int irang_min = irang;
  
  while (n_absolue_max > n_donnees_traitees_rangs[irang + 1])
    irang++;
  
  int irang_max = irang;

  PDM_g_num_t _n_donnees_a_envoyer_proc = PDM_IO_MIN(n_donnees_traitees_rangs[irang_min + 1] + 1 - n_absolue_min, n_donnees);   
  n_donnees_a_envoyer[irang_min] = (int) _n_donnees_a_envoyer_proc; 
  
  int n_donnees_non_traitees = n_donnees - n_donnees_a_envoyer[irang_min];
  for (int i = irang_min + 1; i < irang_max + 1; i++) {
    _n_donnees_a_envoyer_proc = PDM_IO_MIN(n_donnees_traitees_rangs[i + 1] - n_donnees_traitees_rangs[i], 
                   n_donnees_non_traitees);

    n_donnees_a_envoyer[i] = (int) _n_donnees_a_envoyer_proc; 
    n_donnees_non_traitees += - n_donnees_a_envoyer[i];
  }
  
  /* Determination du nombre de donnees a recevoir de chaque processus */
  
  if (rang_actif == 1) {
    
    n_absolue_min = n_donnees_traitees_rangs[fichier->rang] + 1;
    n_absolue_max = n_donnees_traitees_rangs[fichier->rang + 1];
    irang = 0;
    
    while (n_absolue_min > (n_absolue_bloc_rangs[irang+1] - 1))
      irang++;
    
    irang_min = irang;
    
    while (n_absolue_max > (n_absolue_bloc_rangs[irang+1] - 1))
      irang++;
    
    irang_max = irang;
    
    PDM_g_num_t _n_donnees_a_recevoir = PDM_IO_MIN(n_absolue_bloc_rangs[irang_min+1] - n_absolue_min, 
                   n_absolue_bloc_rangs[irang_min+1] - 
                   n_absolue_bloc_rangs[irang_min]);
 
    n_donnees_a_recevoir[irang_min] = (int) _n_donnees_a_recevoir; 
    
     PDM_g_num_t _n_donnees_non_traitees = n_absolue_max - n_absolue_min - 
      n_donnees_a_recevoir[irang_min] + 1;
     
     n_donnees_non_traitees = (int) _n_donnees_non_traitees;
    
    for (int i = irang_min + 1; i < irang_max + 1; i++) {
      _n_donnees_a_recevoir = PDM_IO_MIN(n_absolue_bloc_rangs[i + 1] - n_absolue_bloc_rangs[i], 
                     n_donnees_non_traitees);
      n_donnees_a_recevoir[i] = (int) _n_donnees_a_recevoir;  
      n_donnees_non_traitees += - n_donnees_a_recevoir[i];
    }
    
  }
  
  free(n_absolue_bloc_rangs);
  
  /* Calcul des index  */
  
  for (int i = 0; i < fichier->n_rangs; i++) {
    i_donnees_a_envoyer[i]  = 0;
    i_donnees_a_recevoir[i] = 0;
  }
  
  for (int i = 1; i < fichier->n_rangs; i++) {
    i_donnees_a_envoyer[i] += i_donnees_a_envoyer[i-1] + 
      n_donnees_a_envoyer[i-1];
    i_donnees_a_recevoir[i] += i_donnees_a_recevoir[i-1] + 
      n_donnees_a_recevoir[i-1];
  }
  
  /*------------------------------------------------------------
   * Prise en compte du nombre de composantes
   *------------------------------------------------------------ */
          
  if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
    
    int *n_composantes_recues = 
      (int *) malloc(sizeof(int) * n_absolue_max - n_absolue_min + 1);
    
    PDM_MPI_Alltoallv((void *) n_composantes,
                  n_donnees_a_envoyer,
                  i_donnees_a_envoyer,
                  PDM_MPI_INT, 
                  (void *) n_composantes_recues,
                  n_donnees_a_recevoir,
                  i_donnees_a_recevoir,
                  PDM_MPI_INT, 
                  fichier->comm);
    
    int k = 0;
    for(int i = 0; i < fichier->n_rangs; i++) {
      int n_donnees_a_envoyer_save = n_donnees_a_envoyer[i];
      i_donnees_a_envoyer[i] = 0;
      n_donnees_a_envoyer[i] = 0;
      for(int j = 0; j < n_donnees_a_envoyer_save; j++) { 
        n_donnees_a_envoyer[i] += n_composantes[k];
        k++;
      }
    }
    
    k = 0;
    for(int i = 0; i < fichier->n_rangs; i++) {
      int n_donnees_a_recevoir_save = n_donnees_a_recevoir[i];
      i_donnees_a_recevoir[i] = 0;
      n_donnees_a_recevoir[i] = 0;
      for(int j = 0; j < n_donnees_a_recevoir_save; j++) { 
        n_donnees_a_recevoir[i] += n_composantes_recues[k];
        k++;
      }
    }
    
    for(int i = 1; i < fichier->n_rangs; i++) {
      i_donnees_a_envoyer[i] = i_donnees_a_envoyer[i-1] +
        n_donnees_a_envoyer[i-1];
      i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i-1] +
        n_donnees_a_recevoir[i-1];
    }
    
    free(n_composantes_recues);
    
  }
  
  else {
    for(int i = 0; i < fichier->n_rangs; i++) {
      n_donnees_a_envoyer[i]  = *n_composantes * n_donnees_a_envoyer[i];
      i_donnees_a_envoyer[i]  = *n_composantes * i_donnees_a_envoyer[i];
      n_donnees_a_recevoir[i] = *n_composantes * n_donnees_a_recevoir[i];
      i_donnees_a_recevoir[i] = *n_composantes * i_donnees_a_recevoir[i];
    }  
  }
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Retourne un pointeur sur un fichier a partir de son unite
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *
 * return :
 *   fichier         --> fichier
 *
 *----------------------------------------------------------------------------*/

PDM_io_fichier_t *PDM_io_get_fichier
(const PDM_l_num_t  unite)
{
  return (PDM_io_fichier_t *) PDM_Handles_get (PDM_io_fichiers, unite);
}


/*----------------------------------------------------------------------------
 * Retourne le nomn du fichier ou NULL si pas de fichier
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *
 * return :
 *   fichier         --> fichier 
 *
 *----------------------------------------------------------------------------*/

const char* PDM_io_get_nom_fichier
(const PDM_l_num_t unite)
{
  char *nom = NULL;
  PDM_io_fichier_t *fichier = 
          (PDM_io_fichier_t *) PDM_Handles_get (PDM_io_fichiers, unite);  
  if (fichier != NULL)
    nom = fichier->nom;
  return nom;
}

/*----------------------------------------------------------------------------
 * Ouverture d'un fichier pour acces parallele
 *
 * parameters :
 *   nom             <-- Nom du fichier
 *   fmt             <-- Fichier text ou binaire
 *   suff_t          <-- Type de suffixe (manuel ou automatique)
 *   suff_u          <-- Suffixe (si suffixe manuel)
 *   s_backup        <-- Active le backup d'un fichier preexistant en mode ecriture
 *   accesio         <-- Type (parallele avec mpiio, parallele sans mpiio,
 *                             sequentiel)
 *   mode            <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
 *   pdm_mpi_comm        <-- Communicateur lie au fichier
 *   unite           --> Unite du fichier
 *   ierr            --> Indique si le fichier est de type PDM_io ou non     
 *                       Utiliser uniquement pour une ouverture en lecture
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_open_cf, PDM_IO_OPEN_CF)
(const char            *nom,
 const PDM_l_num_t  *l_nom,
 const int             *fmt,
 const int             *suff_t,
 const char            *suff_u,
 const PDM_l_num_t  *l_suff_u,
 const int             *s_backup,
 const int             *acces,
 const int             *mode,
 const int          *endian,
 PDM_MPI_Fint              *comm,
 double                *prop_noeuds_actifs,
 PDM_l_num_t        *unite,
 PDM_l_num_t        *ierr  
 ARGF_SUPP_CHAINE
 )
{
  char *nom_c    = PDM_fortran_to_c_string(nom, *l_nom);

  char *suff_u_c = NULL;
  if (*suff_t == PDM_IO_SUFF_MAN)
    suff_u_c  = PDM_fortran_to_c_string(suff_u, *l_suff_u);

  const PDM_io_acces_t _acces     = (PDM_io_acces_t) *acces;
  const PDM_io_mode_t _mode       = (PDM_io_mode_t) *mode;
  const PDM_io_suff_t _suff_t     = (PDM_io_suff_t) *suff_t;
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(*comm);
  const PDM_io_fmt_t _fmt         = (PDM_io_fmt_t) *fmt;
  const PDM_io_backup_t _s_backup = (PDM_io_backup_t) *s_backup;
  const PDM_io_endian_t _endian   = (PDM_io_endian_t) *endian;

  PDM_io_open(nom_c, 
              _fmt, 
              _suff_t, 
              suff_u_c, 
              _s_backup,
              _acces, 
              _mode,
              _endian,
              _comm, 
              *prop_noeuds_actifs, 
              unite, 
              ierr);

  free(nom_c);
  if (suff_u_c != NULL)
    free(suff_u_c);
}

void PDM_io_open
(const char             *nom,
 const PDM_io_fmt_t    fmt,
 const PDM_io_suff_t   suff_t,
 const char             *suff_u,
 const PDM_io_backup_t s_backup,
 const PDM_io_acces_t  acces,
 const PDM_io_mode_t   mode,
 const PDM_io_endian_t  endian,
 PDM_MPI_Comm                comm,
 double                  prop_noeuds_actifs,
 PDM_l_num_t         *unite,
 PDM_l_num_t         *ierr  
)
{

  /* Mise a jour du tableau de stockage des fichiers */
  
  *ierr = 0;

  if (PDM_io_fichiers == NULL) {
    PDM_io_fichiers = PDM_Handles_create (4);
  } 

  /* Initialisation de la structure PDM_io_fichier_t */

  PDM_io_fichier_t *nouveau_fichier = 
    (PDM_io_fichier_t*) malloc(sizeof(PDM_io_fichier_t));

  /* Initialisation des timer */

  nouveau_fichier->timer_fichier = PDM_timer_create();
  nouveau_fichier->timer_distribution = PDM_timer_create();
  nouveau_fichier->timer_swap_endian = PDM_timer_create();
  nouveau_fichier->timer_total = PDM_timer_create();

  PDM_timer_resume(nouveau_fichier->timer_total);

  nouveau_fichier->mode       = mode;                    
  nouveau_fichier->acces      = acces;                  
  nouveau_fichier->fmt_t      = fmt;
  nouveau_fichier->fmt        = NULL;
  nouveau_fichier->n_char_fmt = -1;
  nouveau_fichier->n_rangs_actifs = 0;     /* Nombre de rangs actifs */
  nouveau_fichier->n_rangs_inactifs = 0;     /* Nombre de rangs actifs */
  nouveau_fichier->rangs_actifs = NULL;       /* Rangs actifs */
  nouveau_fichier->rangs_inactifs = NULL;       /* Rangs actifs */
  nouveau_fichier->tag_rangs_actifs = NULL;   /* Tag des rangs actifs */
  nouveau_fichier->rang_actif = 1;         /* Indique si rang courant est actif */  

  /* Definition des attributs lies au communicateur MSG */

  PDM_MPI_Comm_rank(comm, &(nouveau_fichier->rang));
  PDM_MPI_Comm_size(comm, &(nouveau_fichier->n_rangs));

  /* Definition des attributs acces, mode */

  /*  if ((acces == PDM_IO_ACCES_SEQ) && (nouveau_fichier->n_rangs > 1)) {*/
  if (suff_t == PDM_IO_SUFF_AUTO) {
    if (acces == PDM_IO_ACCES_SEQ) {
      char format[8];
      int ncharint = 0;
      double _n_rangs = nouveau_fichier->n_rangs;
      while (_n_rangs >= 1.) { 
        _n_rangs = _n_rangs / 10.;
        ncharint++;
      }
      if (ncharint > 9) {
        PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_open :"
                " en mode sequentiel le format d'ecriture limite a 1 milliard de fichier\n");
        abort();
      }

      size_t l_nom = strlen(nom);
      int taille = (int) l_nom + ncharint + 1;
      nouveau_fichier->nom = (char *) malloc(taille + 1);
      sprintf(format,"%cs.%c%1.1d.%1.1dd", '%','%', ncharint, ncharint);
      sprintf(nouveau_fichier->nom, format, nom, nouveau_fichier->rang+1);
    }
    else {
      size_t l_nom = strlen(nom);
      nouveau_fichier->nom = (char *) malloc((int) l_nom + 1 + 2);
      strcpy(nouveau_fichier->nom, nom);
      strcpy(nouveau_fichier->nom + strlen(nom), ".0");
    }
  }
  
  else {
    nouveau_fichier->nom = (char *) malloc(strlen(nom) + 1 + strlen(suff_u));
    strcpy(nouveau_fichier->nom, nom);
    strcpy(nouveau_fichier->nom + strlen(nom), suff_u);
  }

  nouveau_fichier->comm = comm;
  nouveau_fichier->prop_noeuds_actifs = prop_noeuds_actifs;

  /* Test d'existence du fichier en lecture */

  if (mode == PDM_IO_MODE_LECTURE) {
    FILE *testf2 = fopen (nouveau_fichier->nom, "r");
    if (testf2 == NULL) {
      *ierr = 1;
      free (nouveau_fichier->nom);
      free (nouveau_fichier);
      return;
    }
    else {
      fclose (testf2);
    }
  }

  _rangs_actifs(nouveau_fichier);
 
  /* Ouverture du fichier en parallele ou sequentiel suivant la situation.
     En mode Ajout, ouverture dans un premier temps du fichier en mode lecture 
     pour lecture de l'entete */

  PDM_file_seq_mode_t _mode_seq;
  if (mode == PDM_IO_MODE_AJOUT)
    _mode_seq = FICHIER_SEQ_MODE_AJOUT;
  else if (mode == PDM_IO_MODE_LECTURE)
    _mode_seq = FICHIER_SEQ_MODE_LECTURE;
  else if (mode == PDM_IO_MODE_ECRITURE)
    _mode_seq = FICHIER_SEQ_MODE_ECRITURE;

  PDM_file_par_mode_t _mode_par;
  if (mode == PDM_IO_MODE_AJOUT)
    _mode_par = FICHIER_PAR_MODE_AJOUT;
  else if (mode == PDM_IO_MODE_LECTURE)
    _mode_par = FICHIER_PAR_MODE_LECTURE;
  else if (mode == PDM_IO_MODE_ECRITURE)
    _mode_par = FICHIER_PAR_MODE_ECRITURE;

  PDM_file_par_acces_t _acces_par;
  if (acces == PDM_IO_ACCES_MPIIO_EO)
    _acces_par = FICHIER_PAR_ACCES_EO;
  else if (acces == PDM_IO_ACCES_MPIIO_IP)
    _acces_par = FICHIER_PAR_ACCES_IP;

  /* Backup d'un fichier preexistant en mode ecriture */

  nouveau_fichier->backup = PDM_IO_BACKUP_OFF;
  if ((s_backup == PDM_IO_BACKUP_ON) && (mode == PDM_IO_MODE_ECRITURE)) {
    
    int proc_actif = 0;
    if (acces != PDM_IO_ACCES_SEQ) {
      if (nouveau_fichier->rang == 0) {
        proc_actif = 1;
      }
    }
    else {
      proc_actif = 1;
    } 
      
    if (proc_actif == 1) {

      FILE *testf = fopen(nouveau_fichier->nom, "r");
      
      /* Test si un fichier de même nom existe */
      
      if (testf != NULL) {
        nouveau_fichier->backup = PDM_IO_BACKUP_ON;
        fclose(testf);
        char *fichier_backup = malloc(sizeof(char) * (strlen(nouveau_fichier->nom) + 2)); /* caratère ~ + \0 */
        strcpy(fichier_backup, nouveau_fichier->nom);
        strcat(fichier_backup, "~");
        
        /* Si un fichier backup existe, il est supprimé */
        
        FILE *testbackup = fopen(fichier_backup, "r");
        if (testbackup != NULL) {
          fclose(testbackup);
          remove(fichier_backup);
        }
          
        /* renomme le fichier existant en fichier~ */
        
        int s_rename = rename(nouveau_fichier->nom, fichier_backup);
        if (s_rename != 0) {
          PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_open : Impossible de renommer le fichier %s en %s\n",nouveau_fichier->nom, fichier_backup);
          abort();
        }
        else {
          PDM_printf("PDM_io_open : backup du fichier %s avant reecriture\n", nouveau_fichier->nom);
        }
        free(fichier_backup);
      }
    }
  }

  switch(acces) {
  case PDM_IO_ACCES_MPIIO_EO:
  case PDM_IO_ACCES_MPIIO_IP:
    if (nouveau_fichier->n_rangs == 1) {
      nouveau_fichier->PDM_file_seq = PDM_file_seq_open(nouveau_fichier->nom, 
                                                        _mode_seq);
      nouveau_fichier->PDM_file_par = NULL;
    } 
    else {
      nouveau_fichier->PDM_file_seq = NULL;
      nouveau_fichier->PDM_file_par = NULL;
      if (nouveau_fichier->rang_actif) {
        nouveau_fichier->PDM_file_par = PDM_file_par_open(nouveau_fichier->nom,
                                                          _acces_par,
                                                          _mode_par,
                                                          nouveau_fichier->scomm);
      }
    } 
    break;

  case PDM_IO_ACCES_SEQ:
    nouveau_fichier->PDM_file_seq = PDM_file_seq_open(nouveau_fichier->nom, 
                                                    _mode_seq);
    nouveau_fichier->PDM_file_par = NULL;
    break;

  case PDM_IO_ACCES_MPI_SIMPLE:
    if (nouveau_fichier->rang == 0)
      nouveau_fichier->PDM_file_seq = PDM_file_seq_open(nouveau_fichier->nom, 
                                                      _mode_seq);
    else
      nouveau_fichier->PDM_file_seq = NULL;
    nouveau_fichier->PDM_file_par = NULL;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_open : Acces non valide");
    abort();
  }

  /* Test endian */

  int int_endian = 0;
  *((char *)(&int_endian)) = '\1'; /* Ecriture little-endian */

  nouveau_fichier->swap_endian = 0;
  if ((int_endian == 1) && (endian == PDM_IO_BIGENDIAN)) {
    nouveau_fichier->swap_endian = 1;
  }
  else if ((int_endian != 1) && (endian == PDM_IO_LITTLEENDIAN)) {
    nouveau_fichier->swap_endian = 1;
  }
    
  /* Stockage du fichier cree */


  *unite = PDM_Handles_store (PDM_io_fichiers, nouveau_fichier);

  PDM_timer_hang_on(nouveau_fichier->timer_total);

}


/*----------------------------------------------------------------------------
 * pdm_io_seek sets the file position indicator
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   offset          <-- Adresse
 *   seek            <-- Type d'origine
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_seek, PDM_IO_SEEK)
(
const PDM_l_num_t   *unite,
const PDM_g_num_t   *offset,
const PDM_io_seek_t *seek 
)
{
  PDM_io_seek (*unite, *offset, *seek);
}

void PDM_io_seek
(
const PDM_l_num_t   unite,
const PDM_g_num_t   offset,
const PDM_io_seek_t seek 
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
    if (fichier->PDM_file_seq != NULL) {
      long _offset = (long) offset;
      PDM_file_seq_seek (fichier->PDM_file_seq,
                         _offset,
                         (PDM_file_seq_seek_t) seek);
    }
    else if (fichier->PDM_file_par != NULL) {
      PDM_MPI_Offset _offset = (PDM_MPI_Offset) offset;
      PDM_file_par_seek (fichier->PDM_file_par,
                         _offset,
                        (PDM_file_par_seek_t) seek);
    }
    
  }

  else {
    err_code = 1;
  }
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_tell :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}


/*----------------------------------------------------------------------------
 * pdm_io_tell returns the current file position
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   offset          --> Adresse
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_tell, PDM_IO_TELL)
(
const PDM_l_num_t    *unite,
      PDM_g_num_t    *offset
)
{
  PDM_g_num_t  _offset = PDM_io_tell (*unite);
  *offset = _offset;
}

PDM_g_num_t
PDM_io_tell
(
const PDM_l_num_t     unite
)
{
  PDM_g_num_t offset;
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
    if (fichier->PDM_file_seq != NULL) {
      offset = (PDM_g_num_t) PDM_file_seq_tell (fichier->PDM_file_seq);
    }
    else if (fichier->PDM_file_par != NULL) {
      offset = (PDM_g_num_t) PDM_file_par_tell (fichier->PDM_file_par);      
    }
    
  }

  else {
    err_code = 1;
  }
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_tell :"
            " unite '%d' non valide\n", unite);
    abort();
  }
    
  return offset;

}

/*----------------------------------------------------------------------------
 * Lecture globale : Le processus maitre accede seul au fichier et redistribue
 * l'information a l'ensemble des processus du communicateur
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   taille_donnee   <-- Taille unitaire de la donnnee
 *   n_donnees       <-- Nombre de donnees a lire
 *   donnees         --> Donnees lues
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_lecture_globale, PDM_IO_LECTURE_GLOBALE)
(const PDM_l_num_t *unite,
 const PDM_l_num_t *taille_donnee,
 const PDM_g_num_t *n_donnees,
 void                 *donnees
)
{
  PDM_io_lecture_globale(*unite,
                           *taille_donnee,
                           *n_donnees,
                           donnees);
}

void PDM_io_lecture_globale
(const PDM_l_num_t  unite,
 const PDM_l_num_t  taille_donnee,
 const PDM_g_num_t  n_donnees,
 void                 *donnees
 )
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
    
    if (fichier->fmt_t == PDM_IO_FMT_TXT) {
      PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_lecture_globale :\n"
              "Format text non traite\n");
      abort();
    }

    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_swap_endian = fichier->timer_fichier;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;
    
    PDM_timer_resume(timer_total);
    
    PDM_timer_resume(timer_fichier);
    
    /* Lecture */
    
    int n_donnees_lues = 0;
    if (fichier->PDM_file_seq != NULL) {
      PDM_g_num_t n_donnees_lues_gnum = PDM_file_seq_read(fichier->PDM_file_seq,
                                        taille_donnee,
                                          n_donnees,
                                        (void *) donnees);

    
      /* Traitement de l'erreur de lecture */
      
      if (n_donnees_lues_gnum != n_donnees) {
	PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lecture_globale :"
		  " Erreur de lecture dans le fichier '%s' \n", fichier->nom);
	abort();
      }

      if (((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) &&
         fichier->n_rangs > 1) || (fichier->n_rangs_inactifs > 0)
	  || fichier->swap_endian) {
	n_donnees_lues = (int) n_donnees_lues_gnum ;
	if (n_donnees_lues_gnum > 2147483647) {
	  PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lecture_globale :"
		    " Erreur : n_donnees dépasse la taille maximale autorisée (2147483647) dans le fichier '%s' \n", fichier->nom);
	  abort() ;
	}	
      }	
    }
    else if (fichier->PDM_file_par != NULL) {
      if (fichier->rang_actif) {
	/* Vérification de non dépassement de la taille maximale pour n_donnees */
	if (n_donnees > 2147483647) {
	  PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lecture_globale :"
		  " Erreur : n_donnees dépasse la taille maximale autorisée en parallèle (2147483647) dans le fichier '%s' \n", fichier->nom);
	  abort() ;
	}
	int n_donnees_shortint = (int) n_donnees ;
	n_donnees_lues = PDM_file_par_lecture_globale(fichier->PDM_file_par,
                                                      taille_donnee,
                                                      n_donnees_shortint,
                                                      (void *) donnees);
	/* Traitement de l'erreur de lecture */
	
	if (n_donnees_lues != n_donnees_shortint) {
	  PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lecture_globale :"
		    " Erreur de lecture dans le fichier '%s' \n", fichier->nom);
	  abort();
	}
      }
    }
    
    /* Communication des valeurs autres processus si necessaire */

    if (((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) &&
         fichier->n_rangs > 1) || (fichier->n_rangs_inactifs > 0)) {
      PDM_MPI_Bcast(&n_donnees_lues, 1, PDM_MPI_INT, 0, fichier->comm);
      PDM_MPI_Bcast(donnees, n_donnees_lues * taille_donnee,
                PDM_MPI_BYTE, 0, fichier->comm);
    }
    
    PDM_timer_hang_on(timer_fichier);
    
    /* Swap endian */
    
    if (fichier->swap_endian) {
      PDM_timer_resume(timer_swap_endian);
      PDM_io_swap_endian(taille_donnee,
                         n_donnees,
                         donnees,
                   donnees);
      PDM_timer_hang_on(timer_swap_endian);
    }

    PDM_timer_hang_on(timer_total);
  }

  else 
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lecture_globale :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Ecriture globale : Le processus maitre accede seul au fichier
 *
 * parameters :
 *   unite             <-- Unite du fichier
 *   taille_donnee     <-- Taille unitaire de la donnnee
 *   n_donnees         <-- Nombre de donnees a ecrire
 *   donnees            --> Donnees lues
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_ecriture_globale, PDM_IO_ECRITURE_GLOBALE)
(const PDM_l_num_t *unite,
 const PDM_l_num_t *taille_donnee,
 const PDM_g_num_t *n_donnees,
 const void           *donnees
)
{
  PDM_io_ecriture_globale(*unite,
                            *taille_donnee,
                            *n_donnees,
                            donnees);
}

void PDM_io_ecriture_globale
(const PDM_l_num_t  unite,
 const PDM_l_num_t  taille_donnee,
 const PDM_g_num_t  n_donnees,
 const void           *donnees
)
{
  int n_donnees_ecrites = 0;
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {

    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;
    
    PDM_timer_resume(timer_total);
    PDM_timer_resume(timer_fichier);
    
    /* Ecriture globale - ecriture native */

    
    if (fichier->fmt_t == PDM_IO_FMT_TXT) {
        
      unsigned char* _donnees = (unsigned char*) donnees;
      const PDM_g_num_t l_string_donnee = n_donnees * fichier->n_char_fmt + 2;
      char *string_donnee = (char *) malloc(sizeof(char) * l_string_donnee);
      for (PDM_g_num_t i = 0; i < l_string_donnee; i++) {
        string_donnee[i] = '\0';
      }

      char *s_tmp = string_donnee;
      unsigned char* t_buffer = _donnees;
      for (PDM_g_num_t i = 0; i < n_donnees; i++) {
        switch (fichier->data_type) {
        case PDM_IO_T_INT :
          sprintf(s_tmp, fichier->fmt, *((int *) t_buffer));
          break;
        case PDM_IO_T_LONG :
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
          sprintf(s_tmp, fichier->fmt, *((long *) t_buffer));
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
          break;
        case PDM_IO_T_DOUBLE :
          sprintf(s_tmp, fichier->fmt, *((double *) t_buffer));
          break;
        case PDM_IO_T_FLOAT :
          sprintf(s_tmp, fichier->fmt, *((float *) t_buffer));
          break;
        case PDM_IO_T_CHAR :
          sprintf(s_tmp, fichier->fmt, *((char* *) t_buffer));
          break;
        }
        t_buffer += taille_donnee;
        s_tmp += fichier->n_char_fmt;
      }
      string_donnee[strlen(string_donnee)] = '\n';

      if (fichier->PDM_file_seq != NULL) {
        PDM_g_num_t n_donnees_ecrites_gnum = PDM_file_seq_write(fichier->PDM_file_seq,
                                              sizeof(char),
                                              l_string_donnee - 1,
                                              (void *) string_donnee);
	/* Traitement de l'erreur de lecture */
    
	if (n_donnees_ecrites_gnum !=  l_string_donnee - 1) {
	  PDM_error(__FILE__, __LINE__, 0,"[%d] Erreur PDM_io_ecriture_globale :"
		    " Erreur d'ecriture dans le fichier '%s'\n", fichier->rang, fichier->nom);
	  abort();
	  PDM_file_seq_close(fichier->PDM_file_seq);
	}
	if (fichier->acces != PDM_IO_ACCES_SEQ) {
	  if ((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) || (fichier->n_rangs_inactifs > 0)) { 
	    n_donnees_ecrites = (int) n_donnees_ecrites_gnum ;
	    if (n_donnees_ecrites_gnum > 2147483647) {
	      PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
			" Erreur : l_string_donnee dépasse la taille maximale autorisée en parallèle (2147483647) dans le fichier '%s' \n", fichier->nom);
	      abort() ;
	    }
	  }
	}
      }
      else if (fichier->PDM_file_par != NULL) {
        if (fichier->rang_actif) {
	  /* Vérification de non dépassement de la taille maximale pour l_string_donnee */
	  if (l_string_donnee > 2147483647) {
	    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
		      " Erreur : l_string_donnee dépasse la taille maximale autorisée en parallèle (2147483647) dans le fichier '%s' \n", fichier->nom);
	    abort() ;
	  }
	  int l_string_donnee_shortint = (int) l_string_donnee ;

          n_donnees_ecrites = PDM_file_par_ecriture_globale(fichier->PDM_file_par,
                                                         sizeof(char),
                                                          l_string_donnee_shortint - 1,
                                                         (void *) string_donnee);
	  /* Traitement de l'erreur de lecture */
    
	  if (n_donnees_ecrites !=  l_string_donnee_shortint - 1) {
	    PDM_error(__FILE__, __LINE__, 0,"[%d] Erreur PDM_io_ecriture_globale :"
		      " Erreur d'ecriture dans le fichier '%s'\n", fichier->rang, fichier->nom);
	    abort();
	    PDM_file_seq_close(fichier->PDM_file_seq);
	  }
        }
      }
      
      if (fichier->acces != PDM_IO_ACCES_SEQ) {
        if ((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) || (fichier->n_rangs_inactifs > 0)) { 
          PDM_MPI_Bcast(&n_donnees_ecrites, 1, PDM_MPI_INT, 0, fichier->comm);
        }
      }
      
      free(string_donnee);
    }
    else {
      if (fichier->PDM_file_seq != NULL) {
        PDM_g_num_t n_donnees_ecrites_gnum = PDM_file_seq_write(fichier->PDM_file_seq,
                                              taille_donnee,
                                              n_donnees,
                                              (void *) donnees);
	/* Traitement de l'erreur de lecture */
    
	if (n_donnees_ecrites_gnum != n_donnees) {
	  PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
		    " Erreur d'ecriture dans le fichier '%s' \n", fichier->nom);
	  abort();
	}
	if (fichier->acces != PDM_IO_ACCES_SEQ) {
	  if ((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) || (fichier->n_rangs_inactifs > 0)) { 
	    n_donnees_ecrites = (int) n_donnees_ecrites_gnum ;
	    if (n_donnees_ecrites_gnum > 2147483647) {
	      PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
			" Erreur : l_string_donnee dépasse la taille maximale autorisée en parallèle (2147483647) dans le fichier '%s' \n", fichier->nom);
	      abort() ;
	    }
	  }
	}
      }
      else if (fichier->PDM_file_par != NULL) {
        if (fichier->rang_actif) {
	  /* Vérification de non dépassement de la taille maximale pour n_donnees */
	  if (n_donnees > 2147483647) {
	    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lecture_globale :"
		      " Erreur : n_donnees dépasse la taille maximale autorisée en parallèle (2147483647) dans le fichier '%s' \n", fichier->nom);
	    abort() ;
	  }
	  int n_donnees_shortint = (int) n_donnees ;
          n_donnees_ecrites = PDM_file_par_ecriture_globale(fichier->PDM_file_par,
                                                            taille_donnee,
                                                            n_donnees_shortint,
                                                            (void *) donnees);

	  /* Traitement de l'erreur de lecture */
    
	  if (n_donnees_ecrites != n_donnees_shortint) {
	    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
		      " Erreur d'ecriture dans le fichier '%s' \n", fichier->nom);
	    abort();
	  }
        }
      }
      
      if (fichier->acces != PDM_IO_ACCES_SEQ) {
        if ((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) || (fichier->n_rangs_inactifs > 0)) {
          PDM_MPI_Bcast(&n_donnees_ecrites, 1, PDM_MPI_INT, 0, fichier->comm);
        }
      }        
    }
    PDM_timer_hang_on(timer_fichier);
    PDM_timer_hang_on(timer_total);
  }
  else 
    err_code = 1;

  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Lecture parallele de blocs de donnees suivie d'une redistribution des 
 * des donnees suivant l'indirection
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   t_n_composantes <-- Type de tailles composantes 
 *                       (PDM_IO_N_COMPOSANTE_CONSTANT
 *                     ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes   <-- Nombre de composantes pour chaque donnee
 *   taille_donnee   <-- Taille unitaire de la donnnee
 *   n_donnees       <-- Nombre de donnees a lire
 *   indirection     <-- Indirection de redistribition des donnees
 *   donnees         --> Donnees lues
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_lec_par_entrelacee, PDM_IO_LEC_PAR_ENTRELACEE)
(const PDM_l_num_t  *unite,
 const int             *t_n_composantes,         
 const PDM_l_num_t  *n_composantes,
 const PDM_l_num_t  *taille_donnee,
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *indirection,
 void                  *donnees
)
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;
  
  PDM_io_lec_par_entrelacee(*unite,
                              _t_n_composantes,
                              n_composantes,
                              *taille_donnee,
                              *n_donnees,
                              indirection,
                              donnees);
}

void PDM_io_lec_par_entrelacee
(const PDM_l_num_t           unite,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t          *n_composantes,         
 const PDM_l_num_t           taille_donnee,
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t         *indirection,
 void                          *donnees
 )
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  unsigned char* buffer = NULL;
  PDM_g_num_t *index = NULL; 

  if (fichier != NULL) {
    
    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_swap_endian = fichier->timer_fichier;
    PDM_timer_t *timer_distribution = fichier->timer_distribution;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;
    
    if (fichier->fmt_t == PDM_IO_FMT_TXT) {
      PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_lec_par_entrelacee :\n"
              "Format text non traite\n");
      abort();
    }

    PDM_timer_resume(timer_total);
    
    /* Acces sequentiel : sortie en erreur */
    
    /* if (fichier->acces == PDM_IO_ACCES_SEQ) { */
      
    /*   PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_entrelacee :" */
    /*           " Fonction indisponible en mode sÃ©quentiel\n"); */
    /*   abort(); */
      
    /* } */

    /* Processus unique : tri local et appel a une ecriture globale */
    
    if (fichier->n_rangs == 1) {

      PDM_timer_resume(timer_distribution);
 
      int           _n_donnees_buff;
      int            n_octet;
      unsigned char *_donnees = (unsigned char*) donnees;

      /* Calcul de l'indice max */
      PDM_g_num_t _id_max = 0;
      
      for (int i = 0; i < n_donnees; i++) {
        _id_max = PDM_IO_MAX(_id_max, indirection[i]);
      }

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

        index = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t)*
                                           (n_donnees + 1));

        for (int i = 0; i < n_donnees + 1; i++) {
          index[i] = 0;
        }
        for (int i = 0; i < n_donnees; i++) {
          index[indirection[i] - 1 + 1] = n_composantes[i] * taille_donnee;
        }

        for (int i = 1; i < n_donnees + 1; i++) {
          index[i] = index[i] + index[i-1];
        }

        PDM_g_num_t __n_donnees_buff = index[n_donnees] / taille_donnee;
        _n_donnees_buff = (int) __n_donnees_buff;

        buffer = (unsigned char*) malloc(sizeof(unsigned char) * 
                                         index[n_donnees]);

      }

      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        PDM_g_num_t __n_donnees_buff =  _id_max * n_composantes[0];
        _n_donnees_buff = (int) __n_donnees_buff;
        n_octet = taille_donnee * n_composantes[0];
        buffer = (unsigned char*) malloc(sizeof(unsigned char) * 
                                         taille_donnee *
                                         _n_donnees_buff);
      }
			
			else {
				PDM_error(__FILE__, __LINE__, 0,"PDM_io_lec_par_entrelacee Error : unknown PDM_io_n_composantes_t \n");
				abort();
			}
     
      PDM_timer_hang_on(timer_distribution);
      PDM_timer_hang_on(timer_total);

      PDM_io_lecture_globale(unite,
                               taille_donnee,
                               _n_donnees_buff,
                               buffer);

      PDM_timer_resume(timer_total);
      PDM_timer_resume(timer_distribution);

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        int k = 0;
        for (int i = 0; i < n_donnees; i++){
          for (int j = 0; j < n_composantes[i] * taille_donnee; j++){
            _donnees[k] = buffer[index[indirection[i] - 1] + j];
            k++;
          }
        }

        free(index);
      }

      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        for (int i = 0; i < n_donnees; i++){
          for (int j = 0; j < n_octet; j++){
            _donnees[i * n_octet +j] = buffer[(indirection[i] - 1) * n_octet + j];
          }
        }
      }

      free(buffer);

      PDM_timer_hang_on(timer_distribution);
    }

    /* Cas general : Echanges MPI pour prise en compte de 
       l'indirection puis ecriture dans un fichier
       unique */

    else {

      PDM_timer_resume(timer_distribution);
        
      /*---------------------------------------------------------- 
       *  Determination des rangs actifs qui accedent reellement
       *  aux fichiers et du nombre de donnee traitee par
       *  chaque rang 
       *----------------------------------------------------------*/
        
      PDM_g_num_t *n_donnees_rangs = (PDM_g_num_t *) 
        malloc(sizeof(PDM_g_num_t) * (fichier->n_rangs + 1));
        
      /* int *rangs_actifs = NULL; */
      /* int rang_actif = 0; */
      /* int n_rang_actif = 0; */
      int n_donnees_rang_min = 0;
      int n_donnees_rang_max = 0;
        
      PDM_g_num_t _id_max = 0;
      PDM_g_num_t _id_max_max = 0;

      for (int i = 0; i < n_donnees; i++) {
        _id_max = PDM_IO_MAX(_id_max, indirection[i]);
      }

      PDM_MPI_Allreduce(&_id_max,
                    &_id_max_max, 
                    1, 
                    PDM__PDM_MPI_G_NUM, 
                    PDM_MPI_MAX, 
                    fichier->comm);

      _n_donnees_rang(fichier,
                      _id_max_max,
                      n_donnees_rangs,
                      &n_donnees_rang_min,
                      &n_donnees_rang_max);

      /* Allocation du buffer si le rang est actif */
        
      PDM_g_num_t __n_donnees_rang = (n_donnees_rangs[fichier->rang + 1] - 
                                      n_donnees_rangs[fichier->rang]); 
      int _n_donnees_rang = (int) __n_donnees_rang;
      
      /*---------------------------------------
       *  Envoi/Reception des numeros absolus 
       *  decrits par l'indirection 
       *---------------------------------------*/
      
      int *n_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                                fichier->n_rangs);
      int *n_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                 fichier->n_rangs);
        
      /* Pour chaque donnee le proc ou elle va etre envoyee */

      int *donnees_proc = (int *) malloc(sizeof(int) * n_donnees);
      int *indirection_locale = (int *) malloc(sizeof(int) * n_donnees);
        
      /* Calcul du nombre de donnees a envoyer a chaque procesus */
        
      for (int i = 0; i < fichier->n_rangs; i++)
        n_donnees_a_envoyer[i] = 0;
        
      for (int i = 0; i < n_donnees; i++) {
          
        /* Recherche du processus le plus probable */
          
        PDM_g_num_t n_absolu = indirection[i] - 1;
        int irang_actif = fichier->n_rangs_actifs - 1;  

      	if (n_donnees_rang_min > 0) {
          PDM_g_num_t _irang_actif = n_absolu / (PDM_g_num_t) n_donnees_rang_min;
          int __irang_actif = (int) _irang_actif;
          irang_actif = PDM_IO_MIN(__irang_actif, 
			                             fichier->n_rangs_actifs - 1);
	      }
          
        /* Ajustement suivant la numerotation absolue */
          
        while (n_absolu < n_donnees_rangs[fichier->rangs_actifs[irang_actif]]) 
          irang_actif -= 1;
        
        assert(n_absolu < (n_donnees_rangs[fichier->rangs_actifs[irang_actif] + 1]));
        
        n_donnees_a_envoyer[fichier->rangs_actifs[irang_actif]] += 1;
        donnees_proc[i] = fichier->rangs_actifs[irang_actif];
      }

      PDM_MPI_Alltoall(n_donnees_a_envoyer,  1, PDM_MPI_INT, 
                   n_donnees_a_recevoir, 1, PDM_MPI_INT, 
                   fichier->comm);

      int *i_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                                fichier->n_rangs);
      int *i_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                 fichier->n_rangs);
        
      i_donnees_a_envoyer[0] = 0;
      for (int i = 1; i < fichier->n_rangs; i++)
        i_donnees_a_envoyer[i] = i_donnees_a_envoyer[i-1] + 
          n_donnees_a_envoyer[i-1];
      
      i_donnees_a_recevoir[0] = 0;
      for (int i = 1; i < fichier->n_rangs; i++)
        i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i-1] +
          n_donnees_a_recevoir[i-1];
      
      
      PDM_g_num_t *num_absolue_envoyee = 
        (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                   (i_donnees_a_envoyer[fichier->n_rangs - 1] + 
                                    n_donnees_a_envoyer[fichier->n_rangs - 1]));
      
      for (int i = 0; i < fichier->n_rangs; i++)
        n_donnees_a_envoyer[i] = 0;
      
      for (int i = 0; i < n_donnees; i++) {
        int iproc = donnees_proc[i];
        num_absolue_envoyee[i_donnees_a_envoyer[iproc]+
                            n_donnees_a_envoyer[iproc]] =  indirection[i];
        indirection_locale[i_donnees_a_envoyer[iproc]+
                           n_donnees_a_envoyer[iproc]] = i;
        n_donnees_a_envoyer[iproc] += 1;
      }      

      int l_num_absolue_recues = i_donnees_a_recevoir[fichier->n_rangs - 1] +
        n_donnees_a_recevoir[fichier->n_rangs - 1];

      PDM_g_num_t *num_absolue_recues = 
        (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * l_num_absolue_recues);
      
      PDM_MPI_Alltoallv(num_absolue_envoyee, 
                    n_donnees_a_envoyer, 
                    i_donnees_a_envoyer, 
                    PDM__PDM_MPI_G_NUM,
                    num_absolue_recues, 
                    n_donnees_a_recevoir, 
                    i_donnees_a_recevoir, 
                    PDM__PDM_MPI_G_NUM,
                    fichier->comm);

      free(num_absolue_envoyee); 

      /*------------------------------------------------------------
       * Determination de la taille du bloc de chaque rang 
       * (prise en compte du nombre de composantes de chaque donnees)
       *------------------------------------------------------------ */

      int *n_donnees_blocs = (int *) malloc(sizeof(int) * fichier->n_rangs);
      int  n_donnees_bloc = 0;
      int *n_composantes_recues = NULL;
        
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
          
        /* Envoi/Reception des composantes */
          
        int *n_composantes_envoyee = 
          (int *) malloc(sizeof(int) * 
                         (i_donnees_a_envoyer[fichier->n_rangs - 1] + 
                          n_donnees_a_envoyer[fichier->n_rangs - 1]));
          
        for (int i = 0; i < fichier->n_rangs; i++)
          n_donnees_a_envoyer[i] = 0;
          
        for (int i = 0; i < n_donnees; i++) {
          int iproc = donnees_proc[i];
          n_composantes_envoyee[i_donnees_a_envoyer[iproc] +
                                n_donnees_a_envoyer[iproc]]= n_composantes[i];
          n_donnees_a_envoyer[iproc] += 1;
        }      
          
        int l_n_composantes_recues = 
          i_donnees_a_recevoir[fichier->n_rangs - 1] + 
          n_donnees_a_recevoir[fichier->n_rangs - 1];

        n_composantes_recues = (int *) malloc(sizeof(int) * 
                                              l_n_composantes_recues);

        PDM_MPI_Alltoallv(n_composantes_envoyee, 
                      n_donnees_a_envoyer, 
                      i_donnees_a_envoyer, 
                      PDM_MPI_INT, 
                      n_composantes_recues, 
                      n_donnees_a_recevoir, 
                      i_donnees_a_recevoir, 
                      PDM_MPI_INT, 
                      fichier->comm);

        free(n_composantes_envoyee);
          
        /* Determination du debut et de la taille du bloc en octet 
           en fonction de la taille du type et du nombre de composantes 
           de chaque donnees */
        
        int *tag = malloc(sizeof(int) * _n_donnees_rang);
        for (int i = 0; i < _n_donnees_rang; i++) 
          tag[i] = 0;

        for (int i = 0; i < l_n_composantes_recues; i++) {
          PDM_g_num_t _num_abs = num_absolue_recues[i] - 1 - 
                                (PDM_g_num_t) n_donnees_rangs[fichier->rang];
          const int num_abs = (int) _num_abs;
          if (tag[num_abs] == 0) {
            n_donnees_bloc += n_composantes_recues[i];
            tag[num_abs] = 1;
          }
        }          

        free(tag);

        PDM_MPI_Allgather(&n_donnees_bloc, 1, PDM_MPI_INT, 
                           n_donnees_blocs, 1, PDM_MPI_INT, 
                           fichier->comm);

      }

      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        int _n_composantes = *n_composantes;
        for (int i = 0; i < fichier->n_rangs; i++) {
          PDM_g_num_t _n_donnees_rangs = n_donnees_rangs[i + 1] - n_donnees_rangs[i];
          n_donnees_blocs[i] = _n_composantes * (int) _n_donnees_rangs;
        }
        
        PDM_g_num_t _n_donnees_rangs1 = n_donnees_rangs[fichier->rang + 1]-
                                        n_donnees_rangs[fichier->rang];
        
        n_donnees_bloc = _n_composantes * (int) (_n_donnees_rangs1);
        
      }
      
      PDM_timer_hang_on(timer_distribution);

      /*---------------------------------------------------------
       * Lecture parallele des blocs 
       * ou Lecture bloc par bloc via le processus maitre 
       * si acces sequentiel
       *---------------------------------------------------------*/
        
      PDM_timer_resume(timer_fichier);

      int max_n_donnees_bloc = n_donnees_bloc;

      if (fichier->n_rangs > 1) { 
          PDM_MPI_Allreduce(&n_donnees_bloc, &max_n_donnees_bloc, 1, 
                        PDM_MPI_INT, PDM_MPI_MAX, fichier->comm);
      } 
      
      buffer = (unsigned char*) malloc(taille_donnee * max_n_donnees_bloc);

      switch (fichier->acces) {
          
        /* Lecture en parallele des blocs */
          
      case PDM_IO_ACCES_MPIIO_EO:
      case PDM_IO_ACCES_MPIIO_IP:
        {
          PDM_g_num_t debut_bloc = 0;
          for (int i = 0; i < fichier->rang; i++)
            debut_bloc += n_donnees_blocs[i];
            
          if (fichier->rang_actif) {
            int n_donnees_lues = 
              PDM_file_par_lecture_parallele(fichier->PDM_file_par,
                                             taille_donnee,
                                             n_donnees_bloc,
                                             buffer,
                                             debut_bloc);
            
            if (n_donnees_lues != n_donnees_bloc) {
              PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_entrelacee :"
                      " Erreur de lecture du fichier '%s' \n", fichier->nom);
              abort();
            }
          }
          break;
        }
          
        /* Lecture sequentielle des blocs puis envoie 
           au rangs actifs cibles */
          
      case PDM_IO_ACCES_MPI_SIMPLE:       
        {
          int etat_lecture = 1; /* Indicateur permettant de determiner
                                   une erreur de lecture */

          if (fichier->rang == 0) {
            assert(fichier->rang_actif == 1);
              
	    PDM_g_num_t n_donnees_bloc_tmp = (PDM_g_num_t) n_donnees_blocs[0] ;
            PDM_g_num_t _n_donnees_lues = PDM_file_seq_read(fichier->PDM_file_seq,
                                                   taille_donnee,
                                                   n_donnees_bloc_tmp,
                                                   buffer);

            if (_n_donnees_lues != n_donnees_bloc_tmp)
              etat_lecture = 0;

            unsigned char *buffer_tmp = 
              (unsigned char*) malloc(taille_donnee * max_n_donnees_bloc);
               
            for (int i = 1; i < fichier->n_rangs_actifs; i++) {
              int l_buffer = n_donnees_blocs[fichier->rangs_actifs[i]] * taille_donnee;

	      n_donnees_bloc_tmp = (PDM_g_num_t) n_donnees_blocs[fichier->rangs_actifs[i]] ;
              _n_donnees_lues = 
                PDM_file_seq_read(fichier->PDM_file_seq,
                                 taille_donnee,
                                 n_donnees_bloc_tmp,
                                 buffer_tmp);
                
              PDM_MPI_Send(buffer_tmp, l_buffer, PDM_MPI_BYTE, fichier->rangs_actifs[i], 
                       PDM_io_tag, fichier->comm);
                
              if (_n_donnees_lues != n_donnees_bloc_tmp)
                etat_lecture = 0;
                
            }

            free(buffer_tmp);
          }
            
          else if (fichier->rang_actif == 1) {
            int l_buffer = n_donnees_bloc * taille_donnee;
            PDM_MPI_Recv(buffer, l_buffer, PDM_MPI_BYTE, 0, 
                     PDM_io_tag, fichier->comm);
          }
            
          PDM_MPI_Bcast(&etat_lecture, 1, PDM_MPI_INT, 0, fichier->comm);
            
          if (etat_lecture == 0) {
            PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_entrelacee :"
                    " Erreur de lecture du fichier '%s' \n", fichier->nom);
            abort();
          }
            
          break;
        }
      default :
	break;
      }
        
      PDM_timer_hang_on(timer_fichier);

      free(n_donnees_blocs);
        
      /*------------------------------------
       * Distribution suivant l'indirection
       *------------------------------------*/
        
      PDM_timer_resume(timer_distribution);
        
      /* Ordonnancement du buffer pour envoi alltoall */
        
      unsigned char *buffer_ordonne = NULL;
      
      PDM_g_num_t _n_donnees_rang1 = n_donnees_rangs[fichier->rang + 1] - 
                                     n_donnees_rangs[fichier->rang];

      int n_donnees_rang = (int) _n_donnees_rang1;
     
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
          
        int *n_composantes_ordonnees = NULL;
          
        /* Tri du tableau decrivant le nombre de composantes */
          
        n_composantes_ordonnees = (int *) malloc(sizeof(int) * 
                                                 (n_donnees_rang + 1));
          
        for (int i = 0; i < l_num_absolue_recues; i++) {
          PDM_g_num_t _idx = num_absolue_recues[i] - 1 - n_donnees_rangs[fichier->rang];
          int idx = (int) _idx;
          n_composantes_ordonnees[idx+1] = n_composantes_recues[i] * taille_donnee;
        }
          
        n_composantes_ordonnees[0] = 0;
          
        for (int i = 1; i < n_donnees_rang + 1; i++) 
          n_composantes_ordonnees[i] = n_composantes_ordonnees[i] + 
            n_composantes_ordonnees[i-1];

        
        buffer_ordonne = (unsigned char*) malloc(sizeof(unsigned char) *
                                                 n_composantes_ordonnees[n_donnees_rang]);
          
        /* Ordonnancement du buffer pour echange alltoall */
          
        int k1 = 0;
        for (int i = 0; i < l_num_absolue_recues; i++) {
          const PDM_g_num_t _num_loc = num_absolue_recues[i] - 1 - n_donnees_rangs[fichier->rang]; 
          const int num_loc = (int) _num_loc;
            
          const int idx = n_composantes_ordonnees[num_loc];
          const int _n_composantes = n_composantes_ordonnees[num_loc+1] - 
                                     n_composantes_ordonnees[num_loc];
            
          for (int k = 0; k < _n_composantes; k++) {
            buffer_ordonne[k1] = buffer[idx + k];
            k1 += 1;
          }
        }

        free(buffer);
        free(n_composantes_ordonnees);

        /* n_donnees_a_envoyer et n_donnees_a_recevoir sont inverses 
           par rapports au premier alltoallv */
          
        for (int i = 0; i < fichier->n_rangs; i++) {
          int _n_donnees_a_recevoir = n_donnees_a_recevoir[i];
          int _idx_a_recevoir = i_donnees_a_recevoir[i];
            
          int _n_donnees_a_envoyer = n_donnees_a_envoyer[i];
          int _idx_a_envoyer = i_donnees_a_envoyer[i];
            
          n_donnees_a_recevoir[i] = 0;
          for (int k = 0; k < _n_donnees_a_recevoir; k++)
            n_donnees_a_recevoir[i] += 
              n_composantes_recues[_idx_a_recevoir + k] * taille_donnee;
            
          n_donnees_a_envoyer[i] = 0;
          for (int k = 0; k < _n_donnees_a_envoyer; k++) {
            n_donnees_a_envoyer[i] += 
              n_composantes[indirection_locale[_idx_a_envoyer + k]] * 
              taille_donnee;
          }
        }

        i_donnees_a_envoyer[0]  = 0;
        i_donnees_a_recevoir[0] = 0;
        for (int i = 1; i < fichier->n_rangs; i++) {
          i_donnees_a_envoyer[i]  = i_donnees_a_envoyer[i-1] + 
            n_donnees_a_envoyer[i-1];
          i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i-1] + 
            n_donnees_a_recevoir[i-1];
        }
          
        /* n_donnees_a_envoyer et n_donnees_a_recevoir sont inverses 
           par rapports au premier alltoallv */
          
        unsigned char* donnees_tmp = 
          (unsigned char*) malloc(sizeof(unsigned char) *  
                                  (n_donnees_a_envoyer[fichier->n_rangs - 1] +
                                   i_donnees_a_envoyer[fichier->n_rangs -1]));

        PDM_MPI_Alltoallv(buffer_ordonne, 
                      n_donnees_a_recevoir, 
                      i_donnees_a_recevoir, 
                      PDM_MPI_BYTE, 
                      donnees_tmp, 
                      n_donnees_a_envoyer, 
                      i_donnees_a_envoyer, 
                      PDM_MPI_BYTE, 
                      fichier->comm);

        free(buffer_ordonne);
          
        /* Tri des donnees recues */
          
        int *i_composantes = (int*) malloc(sizeof(int) * (n_donnees + 1));
        i_composantes[0] = 0;
        for (int i = 1; i < n_donnees + 1; i++) {
          i_composantes[i] = i_composantes[i-1] + 
            n_composantes[i-1] * taille_donnee;
        }
          
        unsigned char* _donnees = (unsigned char*) donnees;
          
        k1 = 0;
        for (int i = 0; i < n_donnees; i++) {
          const int pas = n_composantes[indirection_locale[i]] * 
            taille_donnee;
          for (int k = 0; k < pas; k++) {
            _donnees[i_composantes[indirection_locale[i]] + k] = 
              donnees_tmp[k1++];
          }
        }
          
        free(i_composantes);
        free(donnees_tmp);
          
      }
        
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        const int _n_octet_composantes = *n_composantes * taille_donnee;
          
        buffer_ordonne = (unsigned char*) malloc(sizeof(unsigned char) *
                                                 l_num_absolue_recues  * _n_octet_composantes);

        for (int i = 0; i < l_num_absolue_recues; i++) {
          const PDM_g_num_t _idx = 
            num_absolue_recues[i] - 1 - n_donnees_rangs[fichier->rang];
          const int idx = (int) _idx;
          for (int k = 0; k < _n_octet_composantes; k++)
            buffer_ordonne[i * _n_octet_composantes + k] = 
              buffer[idx * _n_octet_composantes + k];
        }
          
        free(buffer);
          
        for (int i = 0; i < fichier->n_rangs; i++) {
          n_donnees_a_envoyer[i]  = n_donnees_a_envoyer[i]  * 
            _n_octet_composantes;
          n_donnees_a_recevoir[i] = n_donnees_a_recevoir[i] * 
            _n_octet_composantes;
          i_donnees_a_envoyer[i]  = i_donnees_a_envoyer[i]  * 
            _n_octet_composantes;
          i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i] * 
            _n_octet_composantes;
        }
          
        /* n_donnees_a_envoyer et n_donnees_a_recevoir sont inverses 
           par rapports au premier alltoallv */
          
        unsigned char* donnees_tmp = 
          (unsigned char*) malloc(sizeof(unsigned char) *  
                                  (n_donnees_a_envoyer[fichier->n_rangs - 1] +
                                   i_donnees_a_envoyer[fichier->n_rangs - 1]));
          
        PDM_MPI_Alltoallv(buffer_ordonne, 
                      n_donnees_a_recevoir,
                      i_donnees_a_recevoir,
                      PDM_MPI_BYTE, 
                      donnees_tmp,
                      n_donnees_a_envoyer,
                      i_donnees_a_envoyer,
                      PDM_MPI_BYTE, 
                      fichier->comm);
          
        free(buffer_ordonne);
          
        /* Tri des donnees recues */
          
        unsigned char* _donnees = (unsigned char*) donnees;
          
        const int pas = _n_octet_composantes;
        for (int i = 0; i < n_donnees; i++) {
          for (int k = 0; k < pas; k++) {
            _donnees[indirection_locale[i] * pas + k] = 
              donnees_tmp[i * pas + k];
          }
        }
          
        free(donnees_tmp);
      }
        
      PDM_timer_hang_on(timer_distribution);
        
      /* Endianness */ 
        
      if (fichier->swap_endian) {
          
        PDM_timer_resume(timer_swap_endian);
        PDM_g_num_t l_donnees = 0;
        if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
          l_donnees = 0;
          for (int i = 0; i < n_donnees; i++)
            l_donnees += n_composantes[i];
        }
        else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
          const int _n_composantes = *n_composantes;
          l_donnees = _n_composantes * n_donnees;
        }
          
        PDM_io_swap_endian(taille_donnee,
                     l_donnees,
                     donnees,
                     donnees);
        PDM_timer_hang_on(timer_swap_endian);
      }
        
      PDM_timer_resume(timer_distribution);
        
      free(n_donnees_rangs);       /* n_rangs */
        
      free(n_donnees_a_envoyer);   /* n_rangs */
      free(n_donnees_a_recevoir);  /* n_rangs */
      free(i_donnees_a_envoyer);   /* n_rangs */
      free(i_donnees_a_recevoir);  /* n_rangs */
        
      free(indirection_locale);    /* n_donnees */
        
      free(donnees_proc);          /* n_donnees */
        
      free(num_absolue_recues);    /* n_donnees_buffer */
        
      /* free(rangs_actifs);          /\* n_rangs_actifs *\/ */
        
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        free(n_composantes_recues);
      }
        
      PDM_timer_hang_on(timer_distribution);
    }
    PDM_timer_hang_on(timer_total);
  }

  else 
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_entrelacee :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Lecture parallele de blocs de donnees
 * Les blocs doivent etre rangÃ©s par ordre croissant suivant la numÃ©rotation
 * des processus
 *
 * parameters :
 *   unite             <-- Unite du fichier
 *   t_n_composantes   <-- Type de tailles composantes 
 *                        (PDM_IO_N_COMPOSANTE_CONSTANT
 *                     ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes     <-- Nombre de composantes pour chaque donnee         
 *   taille_donnee     <-- Taille unitaire de la donnnee
 *   debut_bloc        <-- Adresse relative du debut de bloc
 *   n_donnees         <-- Nombre de donnees a lire
 *   donnees           --> Donnees lues
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_lec_par_bloc, PDM_IO_LEC_PAR_BLOC)
(const PDM_l_num_t  *unite,
 const int             *t_n_composantes,         
 const PDM_l_num_t  *n_composantes,         
 const PDM_l_num_t  *taille_donnee,
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *debut_bloc,
 void                  *donnees
)
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;
  
  PDM_io_lec_par_bloc(*unite,
                        _t_n_composantes,
                        n_composantes,
                        *taille_donnee,
                        *n_donnees,
                        *debut_bloc,
                        donnees);
}

void PDM_io_lec_par_bloc
(const PDM_l_num_t           unite,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t          *n_composantes,         
 const PDM_l_num_t           taille_donnee,
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t          debut_bloc,
 void                          *donnees
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  unsigned char* buffer = NULL;

  if (fichier != NULL) {
      
    if (fichier->fmt_t == PDM_IO_FMT_TXT) {
      PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_lec_par_bloc :\n"
              "Format text non traite\n");
      abort();
    }

    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_distribution = fichier->timer_distribution;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;
    PDM_timer_t *timer_swap_endian = fichier->timer_swap_endian;
      
    PDM_timer_resume(timer_total);
      
    /* En acces purement sequentiel sortie en erreur */

    /* if (fichier->acces == PDM_IO_ACCES_SEQ) { */

    /*   PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_bloc :" */
    /*           " Fonction indisponible en acces sequentiel (PDM_IO_ACCES_SEQ) \n"); */
    /*   abort(); */

    /* } */

    if (fichier->n_rangs == 1) {

      int l_donnees = 0;
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        l_donnees = 0;
        for (int i = 0; i < n_donnees; i++)
          l_donnees += n_composantes[i];
      }
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        const int _n_composantes = *n_composantes;
        l_donnees = _n_composantes * n_donnees;
      }

      PDM_timer_hang_on(timer_total);
        
      PDM_io_lecture_globale(unite,
                               taille_donnee,
                               l_donnees,
                               donnees);

      PDM_timer_resume(timer_total);

    }

    else {

      PDM_timer_resume(timer_distribution);
        
      PDM_l_num_t n_donnees_bloc_actif = 0;
      PDM_g_num_t debut_bloc_actif    = 0;

      /*---------------------------------------------------------- 
       *  Determination des rangs actifs qui accedent reellement
       *  aux fichiers et du nombre de donnees traitees par
       *  chaque rang 
       *----------------------------------------------------------*/
        
      PDM_g_num_t *n_donnees_traitees_rangs = (PDM_g_num_t *) 
        malloc(sizeof(PDM_g_num_t) * (fichier->n_rangs + 1));
        
      /* int *rangs_actifs = NULL; */
      /* int  rang_actif = 0; */
      /* int  n_rang_actif = 0; */
      int  n_donnees_traitees_rang_min = 0;
      int  n_donnees_traitees_rang_max = 0;
        
      PDM_g_num_t _id_max = n_donnees;
      PDM_g_num_t _id_max_max = 0;

      PDM_MPI_Allreduce(&_id_max,
                    &_id_max_max, 
                    1, 
                    PDM__PDM_MPI_G_NUM, 
                    PDM_MPI_SUM, 
                    fichier->comm);

      _n_donnees_rang(fichier,
                      _id_max_max,
                      n_donnees_traitees_rangs,
                      &n_donnees_traitees_rang_min,
                      &n_donnees_traitees_rang_max);

      int *n_donnees_a_envoyer = NULL;
      int *n_donnees_a_recevoir = NULL;        
      int *i_donnees_a_envoyer = NULL;
      int *i_donnees_a_recevoir = NULL;

      if (fichier->n_rangs_actifs != fichier->n_rangs) {
 
        /*------------------------------------------------------------
         * Repartition des donnees sur les processus actifs
         *------------------------------------------------------------ */
          
        n_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                             fichier->n_rangs);
        n_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                              fichier->n_rangs);
          
        i_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                             fichier->n_rangs);
        i_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                              fichier->n_rangs);
          
        _calcul_parametres_distribution_bloc(fichier,
                                             t_n_composantes,
                                             n_composantes,
                                             debut_bloc,
                                             n_donnees,
                                             fichier->rang_actif,
                                             n_donnees_traitees_rangs,
                                             n_donnees_a_envoyer,
                                             i_donnees_a_envoyer,
                                             n_donnees_a_recevoir,
                                             i_donnees_a_recevoir);
 
        /*------------------------------------------------------------
         * Calcul de la longueur du bloc et l'adresse de celui-ci
         * dans le fichier pour le processus courant
         *------------------------------------------------------------ */

        int *n_donnees_blocs_actifs = 
          (int *) malloc(sizeof(int) * fichier->n_rangs);
          
        n_donnees_bloc_actif = i_donnees_a_recevoir[fichier->n_rangs - 1] + 
          n_donnees_a_recevoir[fichier->n_rangs - 1];
          
        PDM_MPI_Allgather(&n_donnees_bloc_actif, 1, PDM_MPI_INT, 
                      n_donnees_blocs_actifs, 1, PDM_MPI_INT, 
                      fichier->comm);
            
        debut_bloc_actif = 0;
          
        for (int i = 0; i < fichier->rang; i++) {
          debut_bloc_actif += n_donnees_blocs_actifs[i];
        }
          
        free(n_donnees_blocs_actifs);

        if (0 == 1) {
          PDM_printf("distribution lec: %i ", fichier->rang);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   n_donnees_a_envoyer[i]);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   i_donnees_a_envoyer[i]);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   n_donnees_a_recevoir[i]);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   i_donnees_a_recevoir[i]);
          PDM_printf("\n");
        }

        /*------------------------------------------------------------
         * Prise en compte de la taille de la donnee
         *------------------------------------------------------------ */
          
        for(int i = 0; i < fichier->n_rangs; i++) {
          n_donnees_a_envoyer[i]  = n_donnees_a_envoyer[i] * taille_donnee;
          i_donnees_a_envoyer[i]  = i_donnees_a_envoyer[i] * taille_donnee;
          n_donnees_a_recevoir[i] = n_donnees_a_recevoir[i] * taille_donnee;
          i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i] * taille_donnee;
        }  
          
        /*------------------------------------------------------------
         * Allocation du buffer
         *------------------------------------------------------------ */
          
        buffer = malloc(sizeof(unsigned char) *  
                        (n_donnees_a_recevoir[fichier->n_rangs-1] + 
                         i_donnees_a_recevoir[fichier->n_rangs-1])); 

      }

      else {

        /*------------------------------------------------------------
         * Tous les blocs sont actifs. Pas de deplacement de donnees
         * Il faut juste mettre a jour l'index de debut de bloc et le
         * nombre de donnees en fonction du nombre de composantes
         *------------------------------------------------------------ */

        buffer = (unsigned char *) donnees;

        if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
          n_donnees_bloc_actif = 0;
          for (int i = 0; i < n_donnees; i++) {
            n_donnees_bloc_actif += n_composantes[i];
          }

          int *n_donnees_blocs_actifs = 
            (int *) malloc(sizeof(int) * fichier->n_rangs);
        
          PDM_MPI_Allgather(&n_donnees_bloc_actif, 1, PDM_MPI_INT, 
                        n_donnees_blocs_actifs, 1, PDM_MPI_INT, 
                        fichier->comm);

          debut_bloc_actif = 0;

          for (int i = 0; i < fichier->rang; i++) {
            debut_bloc_actif += n_donnees_blocs_actifs[i];
          }

          free(n_donnees_blocs_actifs);
        }

        else {

          debut_bloc_actif = (debut_bloc - 1) * n_composantes[0];
          n_donnees_bloc_actif = n_donnees * n_composantes[0];

        }
      }

      PDM_timer_hang_on(timer_distribution);

      /*---------------------------------------------------------- 
       *  Lecture du buffer
       *----------------------------------------------------------*/
        
      PDM_timer_resume(timer_fichier);
        
      switch (fichier->acces) {
          
        /* Ecriture parallele des blocs */
          
      case PDM_IO_ACCES_MPIIO_EO:
      case PDM_IO_ACCES_MPIIO_IP:
        if (fichier->rang_actif) {
          PDM_file_par_lecture_parallele(fichier->PDM_file_par,
                                         taille_donnee,
                                         n_donnees_bloc_actif,
                                         buffer,
                                         debut_bloc_actif);
        }
        break;
          
        /* Lecture sequentielle des blocs puis envoi 
           au rangs actifs cibles */
          
      case PDM_IO_ACCES_MPI_SIMPLE:       
        {
            
          int *n_donnees_blocs_actifs = 
            (int *) malloc(sizeof(int) * fichier->n_rangs);
        
          PDM_MPI_Allgather(&n_donnees_bloc_actif, 1, PDM_MPI_INT, 
                        n_donnees_blocs_actifs, 1, PDM_MPI_INT, 
                        fichier->comm);

          int max_n_donnees_blocs_actif = 0;
          for (int i = 0; i < fichier->n_rangs; i++)
            max_n_donnees_blocs_actif = PDM_IO_MAX(max_n_donnees_blocs_actif,
                                                     n_donnees_blocs_actifs[i]);

          int etat_lecture = 1;

          if (fichier->rang == 0) {
            assert(fichier->rang_actif == 1);

            /* Lecture du buffer du proc maitre */
	    PDM_g_num_t n_donnees_blocs_tmp = (PDM_g_num_t) n_donnees_blocs_actifs[0];
              
            PDM_g_num_t donnees_lues = PDM_file_seq_read(fichier->PDM_file_seq,
                                                taille_donnee,
                                                n_donnees_blocs_tmp,
                                                buffer);
              
            unsigned char *buffer_tmp = 
              (unsigned char*) malloc(sizeof(unsigned char) * 
                                      taille_donnee * 
                                      max_n_donnees_blocs_actif);

            if (donnees_lues != n_donnees_blocs_tmp)
              etat_lecture = 0;

            /* Lecture et envoi des buffer */
              
            for (int i = 1; i < fichier->n_rangs_actifs; i++) {
              
	      n_donnees_blocs_tmp = (PDM_g_num_t) n_donnees_blocs_actifs[fichier->rangs_actifs[i]];
              donnees_lues = 
                PDM_file_seq_read(fichier->PDM_file_seq,
                                 taille_donnee,
                                 n_donnees_blocs_tmp,
                                 buffer_tmp);
                
              int l_buffer = n_donnees_blocs_actifs[fichier->rangs_actifs[i]] * taille_donnee;
              PDM_MPI_Send(buffer_tmp, l_buffer, PDM_MPI_BYTE, fichier->rangs_actifs[i], 
                       PDM_io_tag, fichier->comm);
                
              if (donnees_lues != n_donnees_blocs_tmp)
                etat_lecture = 0;
                
            }

            free(buffer_tmp);

          }
          else if (fichier->rang_actif == 1) {

            /* Envoi du buffer au processu maitre si actif */

            int l_buffer = n_donnees_bloc_actif * taille_donnee;
            PDM_MPI_Recv(buffer, l_buffer, PDM_MPI_BYTE, 0, 
                     PDM_io_tag, fichier->comm);
          }

          PDM_MPI_Bcast(&etat_lecture, 1, PDM_MPI_INT, 0, fichier->comm);
            
          if (etat_lecture == 0) {
            PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_bloc :"
                    " Erreur de lecture du fichier '%s' \n", fichier->nom);
            abort();
          }

          free(n_donnees_blocs_actifs);

          break;
        }
      default :
	break;
      }

      /*------------------------------------------------------------
       * Distribution des donnees
       *------------------------------------------------------------ */
        
      if (fichier->n_rangs_actifs != fichier->n_rangs) {
         
        PDM_MPI_Alltoallv(buffer,
                      n_donnees_a_recevoir,
                      i_donnees_a_recevoir,
                      PDM_MPI_BYTE, 
                      donnees,
                      n_donnees_a_envoyer,
                      i_donnees_a_envoyer,
                      PDM_MPI_BYTE, 
                      fichier->comm);
      
        free(n_donnees_a_envoyer);
        free(i_donnees_a_envoyer);
          
        free(n_donnees_a_recevoir);
        free(i_donnees_a_recevoir);

      }

      if (donnees != buffer)
        free(buffer);

      PDM_timer_hang_on(timer_fichier);

      /* Liberation memoire */
        
      PDM_timer_resume(timer_distribution);
        
      PDM_timer_hang_on(timer_distribution);
    }

    /* Endianness */ 
        
    if (fichier->swap_endian) {
        
      PDM_timer_resume(timer_swap_endian);
      PDM_g_num_t l_donnees = 0;
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        l_donnees = 0;
        for (int i = 0; i < n_donnees; i++)
          l_donnees += n_composantes[i];
      }
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        const int _n_composantes = *n_composantes;
        l_donnees = _n_composantes * n_donnees;
      }
        
      PDM_io_swap_endian(taille_donnee,
                         l_donnees,
                         donnees,
                         donnees);

      PDM_timer_hang_on(timer_swap_endian);
    }

    PDM_timer_hang_on(timer_total);
  }

  else 
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_lec_par_bloc :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Tri des donnees suivant l'indirection puis ecriture parallele des blocs de 
 * donnees
 *
 * parameters :
 *   unite             <-- Unite du fichier
 *   t_n_composantes   <-- Type de tailles composantes 
 *                        (PDM_IO_N_COMPOSANTE_CONSTANT
 *                     ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes     <-- Nombre de composantes pour chaque donnee         
 *   taille_donnee     <-- Taille unitaire de la donnnee
 *   n_donnees         <-- Nombre de donnees a lire
 *   indirection       <-- Indirection de redistribition des donnees
 *                       Attention cet argument est un int64
 *   donnees           <-- Donnees a ecrire
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_ecr_par_entrelacee, PDM_IO_ECR_PAR_ENTRELACEE)
(const PDM_l_num_t  *unite,
 const int             *t_n_composantes,         
 const PDM_l_num_t  *n_composantes,         
 const PDM_l_num_t  *taille_donnee,
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *indirection,
 const void            *donnees
)
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;

  PDM_io_ecr_par_entrelacee(*unite,
                              _t_n_composantes,
                              n_composantes,
                              *taille_donnee,
                              *n_donnees,
                              indirection,
                              donnees);
}

void PDM_io_ecr_par_entrelacee
(const PDM_l_num_t           unite,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t          *n_composantes,         
 const PDM_l_num_t           taille_donnee,
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t         *indirection,
 const void                    *donnees
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  unsigned char* buffer = NULL;
  char* s_buffer = NULL;
  int *n_composante_trie = NULL;

  if (fichier != NULL) {


    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_distribution = fichier->timer_distribution;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;
     
    PDM_l_num_t _taille_donnee = taille_donnee;
    if (fichier->fmt_t == PDM_IO_FMT_TXT) {
      _taille_donnee = sizeof(char);
    }

    PDM_timer_resume(timer_total);
      
    /* Processus unique : tri local et appel a une ecriture globale */

    if (fichier->n_rangs == 1) {

      PDM_timer_resume(timer_distribution);
      
      int            _n_donnees;
      unsigned char* _donnees = (unsigned char*) donnees;
      
      /* Calcul de l'indice max */
      PDM_g_num_t _id_max = 0;
      
      for (int i = 0; i < n_donnees; i++) {
        _id_max = PDM_IO_MAX(_id_max, indirection[i]);
      }
      
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        
        PDM_g_num_t *index = 
          (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t)*
                                 (n_donnees + 1));
        
        for (int i = 0; i < n_donnees + 1; i++) {
          index[i] = 0;
        }
        for (int i = 0; i < n_donnees; i++) {
          index[indirection[i] - 1 + 1] = n_composantes[i] * taille_donnee;
        }
        
        for (int i = 1; i < n_donnees + 1; i++) {
          index[i] = index[i] + index[i-1];
        }
        
        PDM_g_num_t __n_donnees = index[n_donnees] / taille_donnee;
        _n_donnees = (int) __n_donnees;
        
        buffer = (unsigned char*) malloc(sizeof(unsigned char) * 
                                         index[n_donnees]);
        
        if (fichier->fmt_t == PDM_IO_FMT_TXT) {
          n_composante_trie =  (int*) malloc(sizeof(int) * _id_max);

	  for (int i = 0; i < _id_max; i++){
	    n_composante_trie[i] = 0;
	  }
	}

        int k = 0;
        for (int i = 0; i < n_donnees; i++){
          if (fichier->fmt_t == PDM_IO_FMT_TXT) {
            n_composante_trie[indirection[i] - 1] = n_composantes[i];
          }
          for (int j = 0; j < n_composantes[i] * taille_donnee; j++){
            buffer[index[indirection[i] - 1] + j] = _donnees[k];
            k++;
          }
        }
        
        free(index);
      }
      
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        _n_donnees = (int) _id_max * n_composantes[0];
        int n_octet = taille_donnee * n_composantes[0];
        buffer = (unsigned char*) malloc(sizeof(unsigned char) * 
                                         taille_donnee *
                                         _n_donnees);
        for (int i = 0; i < n_donnees; i++){
          for (int j = 0; j < n_octet; j++){
            buffer[(indirection[i] - 1) * n_octet + j] = _donnees[i * n_octet +j];
          }
        }
      }
      
      PDM_timer_hang_on(timer_distribution);
      PDM_timer_hang_on(timer_total);

      if (fichier->fmt_t == PDM_IO_FMT_TXT) {
        
        int l_string_donnee = 0;
        if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
          for (int i = 0; i < _id_max; i++) {
            l_string_donnee += n_composante_trie[i] * fichier->n_char_fmt + 1; /* +1 pour '\n' */
          }
        }
        else {
          l_string_donnee = (n_composantes[0] * fichier->n_char_fmt + 1) * (int) _id_max; /* +1 pour '\n' */
        }
        l_string_donnee += 1; /* +1 pour '\0' en fin de chaine */
        
        char *string_donnee = (char *) malloc(sizeof(char) * l_string_donnee);
        for (int i = 0; i < l_string_donnee; i++) {
          string_donnee[i] = '\0';
        }
        
        char *s_tmp = string_donnee;
        unsigned char *t_buffer = buffer;
        for (int i = 0; i < _id_max; i++) {
          if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
            for (int j = 0; j < n_composante_trie[i]; j++) {
              switch (fichier->data_type) {
              case PDM_IO_T_INT :
                sprintf(s_tmp, fichier->fmt, *((int *) t_buffer));
                break;
              case PDM_IO_T_LONG :
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
                sprintf(s_tmp, fichier->fmt, *((long *) t_buffer));
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
                break;
              case PDM_IO_T_DOUBLE :
                sprintf(s_tmp, fichier->fmt, *((double *) t_buffer));
                break;
              case PDM_IO_T_FLOAT :
                sprintf(s_tmp, fichier->fmt, *((float *) t_buffer));
                break;
              case PDM_IO_T_CHAR :
                sprintf(s_tmp, fichier->fmt, *((char* *) t_buffer));
                break;
              }
              t_buffer += taille_donnee;
              s_tmp+=fichier->n_char_fmt;
            }
            *s_tmp = '\n';
            s_tmp++;
          }
          else {
            for (int j = 0; j < n_composantes[0]; j++) {
              switch (fichier->data_type) {
              case PDM_IO_T_INT :
                sprintf(s_tmp, fichier->fmt, *((int *) t_buffer));
                break;
              case PDM_IO_T_LONG :
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
                sprintf(s_tmp, fichier->fmt, *((long *) t_buffer));
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
                break;
              case PDM_IO_T_DOUBLE :
                sprintf(s_tmp, fichier->fmt, *((double *) t_buffer));
                break;
              case PDM_IO_T_FLOAT :
                sprintf(s_tmp, fichier->fmt, *((float *) t_buffer));
                break;
              case PDM_IO_T_CHAR :
                sprintf(s_tmp, fichier->fmt, *((char* *) t_buffer));
                break;
              }
              t_buffer += taille_donnee;
              s_tmp+=fichier->n_char_fmt;
            }
            *s_tmp = '\n';
            s_tmp++;
            assert((s_tmp - string_donnee) < l_string_donnee);
          }
        }
        
        int n_donnees_ecrites;

        if (fichier->PDM_file_seq != NULL) {
	  PDM_g_num_t l_string_donnee_gnum = l_string_donnee - 1 ;
	  PDM_g_num_t n_donnees_ecrites_gnum = PDM_file_seq_write(fichier->PDM_file_seq,
                                                 sizeof(char),
                                                 l_string_donnee_gnum,
                                                 (void *) string_donnee);
	  // Fonction intrinsèquement parallèle <=> pas de vérification
	  // du dépassement de la taille max d'un int
	  n_donnees_ecrites = (int) n_donnees_ecrites_gnum ;
        }
        else if (fichier->PDM_file_par != NULL) {
          if (fichier->rang_actif) {
            n_donnees_ecrites = PDM_file_par_ecriture_globale(fichier->PDM_file_par,
                                                              sizeof(char),
                                                              l_string_donnee - 1,
                                                             (void *) string_donnee);
          }
        }

        free(string_donnee);
        
        if (fichier->fmt_t == PDM_IO_FMT_TXT) {
          free(n_composante_trie);
        }

        if (fichier->acces != PDM_IO_ACCES_SEQ) {
          if ((fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) || (fichier->n_rangs_inactifs > 0)) {
            PDM_MPI_Bcast(&n_donnees_ecrites, 1, PDM_MPI_INT, 0, fichier->comm);
          }
        }
        
        if (n_donnees_ecrites != l_string_donnee - 1) {
          PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecriture_globale :"
                  " Erreur d'ecriture dans le fichier '%s' \n", fichier->nom);
          abort();
        }
      }

      /* Ecriture binaire */

      else {
        PDM_io_ecriture_globale(unite,
                                taille_donnee,
                                _n_donnees,
                                buffer);
      }
     
      PDM_timer_resume(timer_total);
      free(buffer);

    }

    else {
     
      PDM_timer_resume(timer_distribution);
        
      /*---------------------------------------------------------- 
       *  Determination des rangs actifs qui accedent reellement
       *  aux fichiers et du nombre de donnee traitee par
       *  chaque rang 
       *----------------------------------------------------------*/
        
      PDM_g_num_t *n_donnees_rangs = (PDM_g_num_t *) 
        malloc(sizeof(PDM_g_num_t) * (fichier->n_rangs + 1));
        
      /* int *rangs_actifs = NULL; */
      /* int  rang_actif = 0; */
      /* int  n_rang_actif = 0; */
      int  n_donnees_rang_min = 0;
      int  n_donnees_rang_max = 0;
        
      PDM_g_num_t _id_max = 0;
      PDM_g_num_t _id_max_max = 0;

      for (int i = 0; i < n_donnees; i++) {
        _id_max = PDM_IO_MAX(_id_max, indirection[i]);
      }

      PDM_MPI_Allreduce(&_id_max,
                        &_id_max_max, 
                        1, 
                        PDM__PDM_MPI_G_NUM, 
                        PDM_MPI_MAX, 
                        fichier->comm);

      _n_donnees_rang(fichier,
                      _id_max_max,
                      n_donnees_rangs,
                      &n_donnees_rang_min,
                      &n_donnees_rang_max);

      /*---------------------------------------
       *  Envoi/Reception des numeros absolus 
       *  decrits par l'indirection 
       *---------------------------------------*/
      
      int *n_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                                fichier->n_rangs);
      int *n_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                 fichier->n_rangs);
        
      /* Pour chaque donnee le proc ou elle va etre envoyee */
   
      int *donnees_proc = (int *) malloc(sizeof(int) * n_donnees); 
      int *indirection_locale = (int *) malloc(sizeof(int) * n_donnees);
        
      /* Calcul du nombre de donnees a envoyer a chaque procesus */
        
      for (int i = 0; i < fichier->n_rangs; i++)
        n_donnees_a_envoyer[i] = 0;
        
      for (int i = 0; i < n_donnees; i++) {
      
        /* Recherche du processus le plus probable */

        PDM_g_num_t n_absolu = indirection[i] - 1;
        int irang_actif = fichier->n_rangs_actifs - 1;

        if (n_donnees_rang_min > 0) {
          PDM_g_num_t step = n_absolu / n_donnees_rang_min;
          irang_actif = PDM_IO_MIN((int) step, 
                                   fichier->n_rangs_actifs - 1) ;
        }

        /* Ajustement */

        while (n_absolu < n_donnees_rangs[fichier->rangs_actifs[irang_actif]])
          irang_actif -= 1;

        assert(n_absolu < (n_donnees_rangs[fichier->rangs_actifs[irang_actif] + 1])); 

        n_donnees_a_envoyer[fichier->rangs_actifs[irang_actif]] += 1;
        donnees_proc[i] = fichier->rangs_actifs[irang_actif];
      }
        
      PDM_MPI_Alltoall(n_donnees_a_envoyer,  1, PDM_MPI_INT, 
                       n_donnees_a_recevoir, 1, PDM_MPI_INT, 
                       fichier->comm);

      int *i_donnees_a_envoyer  = (int *) malloc(sizeof(int) * 
                                                 fichier->n_rangs);
      int *i_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                 fichier->n_rangs);
        
      i_donnees_a_envoyer[0] = 0;
      for (int i = 1; i < fichier->n_rangs; i++)
        i_donnees_a_envoyer[i] = i_donnees_a_envoyer[i-1] + 
          n_donnees_a_envoyer[i-1];
       
      i_donnees_a_recevoir[0] = 0;
      for (int i = 1; i < fichier->n_rangs; i++)
        i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i-1] + 
          n_donnees_a_recevoir[i-1];
      
      PDM_g_num_t *num_absolue_envoyee = 
        (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                               (i_donnees_a_envoyer[fichier->n_rangs - 1] + 
                                n_donnees_a_envoyer[fichier->n_rangs - 1]));
 
      for (int i = 0; i < fichier->n_rangs; i++)
        n_donnees_a_envoyer[i] = 0;
      
      for (int i = 0; i < n_donnees; i++) {
        int iproc = donnees_proc[i];
        num_absolue_envoyee[i_donnees_a_envoyer[iproc]+
                            n_donnees_a_envoyer[iproc]] = indirection[i];
        indirection_locale[i_donnees_a_envoyer[iproc] + 
                           n_donnees_a_envoyer[iproc]]  = i;
        n_donnees_a_envoyer[iproc] += 1;
      }      

      int _n_quantites = i_donnees_a_recevoir[fichier->n_rangs - 1] +
        n_donnees_a_recevoir[fichier->n_rangs - 1];

      PDM_g_num_t *num_absolue_recues = 
        (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_quantites);

      PDM_MPI_Alltoallv(num_absolue_envoyee, 
                        n_donnees_a_envoyer,
                        i_donnees_a_envoyer,
                        PDM__PDM_MPI_G_NUM,
                        num_absolue_recues,
                        n_donnees_a_recevoir,
                        i_donnees_a_recevoir,
                        PDM__PDM_MPI_G_NUM,
                        fichier->comm);

      free(num_absolue_envoyee);
        
      int *n_donnees_blocs = (int *) malloc(sizeof(int) * fichier->n_rangs);
      int  n_donnees_bloc = 0;
        
      int *n_composantes_recues  = NULL;
      int *n_composantes_envoyee = NULL;
        
      unsigned char* donnees_alltoall = NULL;
      unsigned char* blocs_alltoall = NULL;
        
      int l_blocs_alltoall;

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

        /*------------------------------------------------------------
         * Envoi/reception du nombre de composantes
         *------------------------------------------------------------ */

        n_composantes_envoyee = 
          (int *) malloc(sizeof(int) * 
                         (i_donnees_a_envoyer[fichier->n_rangs - 1] + 
                          n_donnees_a_envoyer[fichier->n_rangs - 1]));
        
        for (int i = 0; i < fichier->n_rangs; i++)
          n_donnees_a_envoyer[i] = 0;
          
        for (int i = 0; i < n_donnees; i++) {
          int iproc = donnees_proc[i];
          n_composantes_envoyee[i_donnees_a_envoyer[iproc] + 
                                n_donnees_a_envoyer[iproc]] = 
            n_composantes[i];
          n_donnees_a_envoyer[iproc] += 1;
        }      

        int l_n_composantes_recues = 
          i_donnees_a_recevoir[fichier->n_rangs - 1] + 
          n_donnees_a_recevoir[fichier->n_rangs - 1];
          
        n_composantes_recues = (int *) malloc(sizeof(int) * 
                                              l_n_composantes_recues );

        PDM_MPI_Alltoallv(n_composantes_envoyee,
                          n_donnees_a_envoyer,
                          i_donnees_a_envoyer,
                          PDM_MPI_INT, 
                          n_composantes_recues,
                          n_donnees_a_recevoir,
                          i_donnees_a_recevoir,
                          PDM_MPI_INT, 
                          fichier->comm);

        /*---------------------------------------------------------- 
         *  preparation des arguments pour l'echange croise des 
         *  donnees
         *----------------------------------------------------------*/

        for (int i = 0; i < fichier->n_rangs; i++) {
          int ideb = i_donnees_a_envoyer[i];
          int ifin = i_donnees_a_envoyer[i] + n_donnees_a_envoyer[i];
            
          n_donnees_a_envoyer[i] = 0;
          for (int k = ideb;  k < ifin; k++)
            n_donnees_a_envoyer[i] += n_composantes_envoyee[k];
            
          n_donnees_a_envoyer[i] *= taille_donnee;
            
          if (i > 0)
            i_donnees_a_envoyer[i] = i_donnees_a_envoyer[i-1] + 
              n_donnees_a_envoyer[i-1];
          else
            i_donnees_a_envoyer[i] = 0;
         
          ideb = i_donnees_a_recevoir[i];
          ifin = i_donnees_a_recevoir[i] + n_donnees_a_recevoir[i];
            
          n_donnees_a_recevoir[i] = 0;
          for (int k = ideb;  k < ifin; k++)
            n_donnees_a_recevoir[i] += n_composantes_recues[k];
            
          n_donnees_a_recevoir[i] *= taille_donnee;
            
          if (i > 0)
            i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i-1] + 
              n_donnees_a_recevoir[i-1];
          else
            i_donnees_a_recevoir[i] = 0;

        }

        free(n_composantes_envoyee);
          
        /*---------------------------------------------------------- 
         *  Rangement des donnees pour l'echange croise
         *  (regroupement des donnees par rang)
         *----------------------------------------------------------*/
          
        const int l_donnees_alltoall = 
          i_donnees_a_envoyer[fichier->n_rangs - 1] +
          n_donnees_a_envoyer[fichier->n_rangs -1];
        
        l_blocs_alltoall = i_donnees_a_recevoir[fichier->n_rangs - 1] + 
          n_donnees_a_recevoir[fichier->n_rangs -1];

        donnees_alltoall = (unsigned char *) 
          malloc(sizeof(unsigned char) * 
                 l_donnees_alltoall);

        int *idx_init = (int*) malloc(sizeof(int) * (n_donnees+1));

        idx_init[0] = 0;
        for (int i = 1; i < n_donnees + 1; i++)
          idx_init[i] = idx_init[i-1] + (n_composantes[i-1] * taille_donnee);

        const unsigned char* _donnees = (const unsigned char*) donnees;
          
        int k = 0;
        for (int i = 0; i < n_donnees; i++) {
          int _n_octets = n_composantes[indirection_locale[i]] * 
            taille_donnee;
          for (int j = 0; j < _n_octets; j++) {
            donnees_alltoall[k] = _donnees[idx_init[indirection_locale[i]]+j];
            k = k + 1;
          }
        }
        
        free(idx_init);
 
        blocs_alltoall = (unsigned char *) 
          malloc(sizeof(unsigned char) * l_blocs_alltoall);

      }
      
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
  
        /*---------------------------------------------------------- 
         *  Preparation des arguments pour l'echange croise des 
         *  donnees
         *----------------------------------------------------------*/

        const int _n_composantes = n_composantes[0];
          
        for (int i = 0; i < fichier->n_rangs; i++) {
          
          i_donnees_a_envoyer[i] = i_donnees_a_envoyer[i] * 
            _n_composantes * 
            taille_donnee;
          i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i] *
            _n_composantes * 
            taille_donnee;

          n_donnees_a_envoyer[i] = n_donnees_a_envoyer[i] * 
            _n_composantes * 
            taille_donnee;
          n_donnees_a_recevoir[i] = n_donnees_a_recevoir[i] * 
            _n_composantes * 
            taille_donnee;
        }

        /*---------------------------------------------------------- 
         *  Rangement des donnees pour l'echange croise
         *  (regroupement des donnees par rang)
         *----------------------------------------------------------*/

        const int l_donnees_alltoall = 
          i_donnees_a_envoyer[fichier->n_rangs - 1] +
          n_donnees_a_envoyer[fichier->n_rangs -1];
       
        l_blocs_alltoall = i_donnees_a_recevoir[fichier->n_rangs - 1] + 
          n_donnees_a_recevoir[fichier->n_rangs -1];
          
        donnees_alltoall = (unsigned char *) 
          malloc(sizeof(unsigned char) * l_donnees_alltoall * taille_donnee);
          
        const unsigned char* _donnees = (const unsigned char*) donnees;
          
        int k = 0;
        for (int i = 0; i < n_donnees; i++) {
          const int idx = _n_composantes * taille_donnee * 
            indirection_locale[i];
          for (int j = 0; j < (_n_composantes * taille_donnee); j++)
            donnees_alltoall[k++] = _donnees[idx + j];
        }
        
        blocs_alltoall = (unsigned char *) 
          malloc(sizeof(unsigned char) * l_blocs_alltoall);
          
      }

      free(donnees_proc);
      free(indirection_locale);
        
      PDM_MPI_Alltoallv(donnees_alltoall,
                        n_donnees_a_envoyer,
                        i_donnees_a_envoyer,
                        PDM_MPI_BYTE, 
                        blocs_alltoall,
                        n_donnees_a_recevoir,
                        i_donnees_a_recevoir,
                        PDM_MPI_BYTE, 
                        fichier->comm);

      free(donnees_alltoall);
      free(n_donnees_a_envoyer);
      free(i_donnees_a_envoyer);
      
      free(n_donnees_a_recevoir);
      free(i_donnees_a_recevoir);
        

      /*---------------------------------------------------------- 
       *  Tri du buffer suivant la numerotation absolue 
       *----------------------------------------------------------*/

      size_t l_s_buffer = 0;

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE && 
          fichier->acces == PDM_IO_ACCES_MPI_SIMPLE) {

        int max_l_blocs_alltoall;
        PDM_MPI_Allreduce(&l_blocs_alltoall, &max_l_blocs_alltoall, 1, 
                          PDM_MPI_INT, PDM_MPI_MAX, fichier->comm);

        if (fichier->rang == 0) {
          buffer = (unsigned char *) 
            malloc(sizeof(unsigned char) * max_l_blocs_alltoall);
          if (fichier->fmt_t == PDM_IO_FMT_TXT) {
            l_s_buffer = (max_l_blocs_alltoall / taille_donnee) * (fichier->n_char_fmt + 1) + 1;
            s_buffer = (char *) malloc(sizeof(char) * l_s_buffer);
            for (int i = 0; i < l_s_buffer; i++)
              s_buffer[i] = '\0';
          }
        }
        else {
          buffer = (unsigned char *) 
            malloc(sizeof(unsigned char) * l_blocs_alltoall);
          if (fichier->fmt_t == PDM_IO_FMT_TXT) {
            l_s_buffer = (l_blocs_alltoall / taille_donnee) * (fichier->n_char_fmt + 1) + 1;
            s_buffer = (char *) malloc(sizeof(char) * l_s_buffer);
            for (int i = 0; i < l_s_buffer; i++)
              s_buffer[i] = '\0';
          }
        }
      }
      else  {
        buffer = (unsigned char *)  
          malloc(sizeof(unsigned char) * l_blocs_alltoall);
        
        if (fichier->fmt_t == PDM_IO_FMT_TXT) {
          s_buffer = (char *) malloc(sizeof(char) * ((l_blocs_alltoall / taille_donnee) 
                                                     * (fichier->n_char_fmt + 1) + 1));
        }
      }

      n_donnees_bloc = 0;

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {

        const PDM_g_num_t __n_donnees_rang = 
          (n_donnees_rangs[fichier->rang + 1] - 
           n_donnees_rangs[fichier->rang]);

        const int _n_donnees_rang = (int) __n_donnees_rang;
        int *tag = malloc(sizeof(int) * _n_donnees_rang);
        for (int i = 0; i < _n_donnees_rang; i++) 
          tag[i] = 0;

        int *index_buffer = (int *) malloc(sizeof(int) * 
                                           (_n_quantites + 1));

        if (fichier->fmt_t == PDM_IO_FMT_TXT) {
          n_composante_trie =  (int*) malloc(sizeof(int) * _n_donnees_rang);
        }

        for (int i = 0; i < _n_quantites + 1; i++)
          index_buffer[i] = 0;

        for (int i = 0; i < _n_quantites; i++) {
          const PDM_g_num_t _num_abs = num_absolue_recues[i] - 1 - 
            n_donnees_rangs[fichier->rang];
          
          const int num_abs = (int) _num_abs;
          index_buffer[num_abs + 1] = n_composantes_recues[i] * taille_donnee;
        }
        
        for (int i = 1; i < _n_quantites + 1; i++) {
          index_buffer[i] = index_buffer[i-1] + index_buffer[i];
        }

        int k1 = 0;

        for (int i = 0; i < _n_quantites; i++) {
          const PDM_g_num_t _num_abs = num_absolue_recues[i] - 1 - 
            n_donnees_rangs[fichier->rang];

          const int num_abs = (int) _num_abs;
          if (tag[num_abs] == 0) {
            n_donnees_bloc += n_composantes_recues[i];
            tag[num_abs] = 1;
          }
          if (fichier->fmt_t == PDM_IO_FMT_TXT) {
            n_composante_trie[num_abs] = n_composantes_recues[i];
          }
          for (int k = 0; k < n_composantes_recues[i] * taille_donnee; k++) {
            buffer[index_buffer[num_abs] + k] = blocs_alltoall[k1++];
          }
        }

        /* s_buffer pour le cas formate */

        if (fichier->fmt_t == PDM_IO_FMT_TXT) {
          n_donnees_bloc = 0;
          char *s_tmp = s_buffer;
          unsigned char *t_buffer = buffer;
          for (int i = 0; i < _n_donnees_rang; i++) {
            for (int j = 0; j < n_composante_trie[i]; j++) {
              n_donnees_bloc += fichier->n_char_fmt;
              switch (fichier->data_type) {
              case PDM_IO_T_INT :
                sprintf(s_tmp, fichier->fmt, *((int *) t_buffer));
                break;
              case PDM_IO_T_LONG :
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
                sprintf(s_tmp, fichier->fmt, *((long *) t_buffer));
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
                break;
              case PDM_IO_T_DOUBLE :
                sprintf(s_tmp, fichier->fmt, *((double *) t_buffer));
                break;
              case PDM_IO_T_FLOAT :
                sprintf(s_tmp, fichier->fmt, *((float *) t_buffer));
                break;
              case PDM_IO_T_CHAR :
                sprintf(s_tmp, fichier->fmt, *((char* *) t_buffer));
                break;
              }
              t_buffer += taille_donnee;
              s_tmp += fichier->n_char_fmt;
            }
            *s_tmp = '\n';
            s_tmp++;
            n_donnees_bloc += 1;
          }
          free(buffer);
          buffer = (unsigned char *) s_buffer;
        }

        PDM_MPI_Allgather(&n_donnees_bloc, 1, PDM_MPI_INT, 
                          n_donnees_blocs, 1, PDM_MPI_INT, 
                          fichier->comm);
 
        free(index_buffer);
        free(n_composantes_recues);
        free(tag);
      }

      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        const int _n_octet_composantes = n_composantes[0] * taille_donnee;
        const PDM_g_num_t __n_donnees_rang = 
          (n_donnees_rangs[fichier->rang + 1] - 
           n_donnees_rangs[fichier->rang]);

        const int _n_donnees_rang = (int) __n_donnees_rang;
        int *tag = malloc(sizeof(int) * _n_donnees_rang);
        for (int i = 0; i < _n_donnees_rang; i++) 
          tag[i] = 0;

        n_donnees_bloc = 0;

        for (int i = 0; i < _n_quantites; i++) {
          const PDM_g_num_t _num_abs = num_absolue_recues[i] - 1 - 
            n_donnees_rangs[fichier->rang];
          const int num_abs = (int) _num_abs;
          if (tag[num_abs] == 0) {
            n_donnees_bloc += n_composantes[0];
            tag[num_abs] = 1;
          }
          for (int k = 0; k < _n_octet_composantes; k++) {
            buffer[num_abs * _n_octet_composantes + k] = 
              blocs_alltoall[i * _n_octet_composantes + k];
          }
        }

        /* s_buffer pour le cas formate */

        if (fichier->fmt_t == PDM_IO_FMT_TXT) {
          n_donnees_bloc = 0;
          char *s_tmp = s_buffer;
          unsigned char* t_buffer = buffer;
          for (int i = 0; i < _n_donnees_rang ; i++) {
            for (int j = 0; j < n_composantes[0]; j++) {
              switch (fichier->data_type) {
              case PDM_IO_T_INT :
                sprintf(s_tmp, fichier->fmt, *((int *) t_buffer));
                break;
              case PDM_IO_T_LONG :
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
                sprintf(s_tmp, fichier->fmt, *((long *) t_buffer));
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
                break;
              case PDM_IO_T_DOUBLE :
                sprintf(s_tmp, fichier->fmt, *((double *) t_buffer));
                break;
              case PDM_IO_T_FLOAT :
                sprintf(s_tmp, fichier->fmt, *((float *) t_buffer));
                break;
              case PDM_IO_T_CHAR :
                sprintf(s_tmp, fichier->fmt, *((char* *) t_buffer));
                break;
              }
              t_buffer += taille_donnee;
              s_tmp += fichier->n_char_fmt;
              n_donnees_bloc += fichier->n_char_fmt;
            }
            *s_tmp = '\n';
            s_tmp++;
            n_donnees_bloc += 1;
          }
          free(buffer);
          buffer = (unsigned char *) s_buffer;
        }

        PDM_MPI_Allgather(&n_donnees_bloc, 1, PDM_MPI_INT, 
                          n_donnees_blocs, 1, PDM_MPI_INT, 
                          fichier->comm);
        free(tag);
      }

      if (n_composante_trie != NULL)
        free(n_composante_trie);

      free(blocs_alltoall);
      free(n_donnees_rangs);

      PDM_timer_hang_on(timer_distribution);
      
      /*---------------------------------------------------------- 
       *  Ecriture du buffer
       *----------------------------------------------------------*/
        
      PDM_timer_resume(timer_fichier);
        
      switch (fichier->acces) {
        
        /* Ecriture parallele des blocs */
        
      case PDM_IO_ACCES_MPIIO_EO:
      case PDM_IO_ACCES_MPIIO_IP:
        {
          PDM_g_num_t debut_bloc = 0;
          for (int i = 0; i < fichier->rang; i++) 
            debut_bloc += n_donnees_blocs[i]; 

          if (fichier->rang_actif) {
            PDM_file_par_ecriture_parallele(fichier->PDM_file_par,
                                            _taille_donnee,
                                            n_donnees_bloc,
                                            buffer,
                                            debut_bloc);
          }
          break;
        }
          
        /* Ecriture sequentielle des blocs puis envoi
           au rangs actifs cibles */

      case PDM_IO_ACCES_MPI_SIMPLE:       
        {
          int etat_ecriture = 1;

          if (fichier->rang == 0) {
            assert((fichier->rang_actif == 1) ||
                   ((fichier->rang_actif == 0) && (n_donnees_blocs[0] == 0)));

            /* Ecriture du buffer du proc maitre */
            PDM_g_num_t n_donnees_blocs_tmp = n_donnees_blocs[0] ;  
            PDM_g_num_t donnees_ecrites = PDM_file_seq_write(fichier->PDM_file_seq,
                                                     _taille_donnee,
                                                     n_donnees_blocs_tmp,
                                                     buffer);
              
            if (donnees_ecrites != n_donnees_blocs_tmp)
              etat_ecriture = 0;

            /* Reception des buffers des autres processus actifs 
               et ecriture */
              
            for (int i = 1; i < fichier->n_rangs_actifs; i++) {
              int l_buffer = n_donnees_blocs[fichier->rangs_actifs[i]] * _taille_donnee;
              PDM_MPI_Recv(buffer, l_buffer, PDM_MPI_BYTE, fichier->rangs_actifs[i], 
                           PDM_io_tag, fichier->comm);
                
	      n_donnees_blocs_tmp = n_donnees_blocs[fichier->rangs_actifs[i]];
              donnees_ecrites = 
                PDM_file_seq_write(fichier->PDM_file_seq,
                                   _taille_donnee,
                                   n_donnees_blocs_tmp,
                                   buffer);
              if (donnees_ecrites != n_donnees_blocs_tmp)
                etat_ecriture = 0;
                
            }
          }
          else if (fichier->rang_actif == 1) {

            /* Envoi du buffer au processus maitre si actif */

            int l_buffer = n_donnees_bloc * _taille_donnee;
            PDM_MPI_Send(buffer, l_buffer, PDM_MPI_BYTE, 0, 
                         PDM_io_tag, fichier->comm);
          }

          PDM_MPI_Bcast(&etat_ecriture, 1, PDM_MPI_INT, 0, fichier->comm);
            
          if (etat_ecriture == 0) {
            PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecr_par_entrelacee :"
                    " Erreur d'ecriture du fichier '%s' \n", fichier->nom);
            abort();
          }

          break;
        }
      default:
        break;
      }

      PDM_timer_hang_on(timer_fichier);

      /* Liberation memoire */
        
      PDM_timer_resume(timer_distribution);
      
      if (buffer != NULL)
        free(buffer);
        
      free(n_donnees_blocs);       /* n_rangs */
      free(num_absolue_recues);    /* n_donnees_buffer */
        
      PDM_timer_hang_on(timer_distribution);
    }

    PDM_timer_hang_on(timer_total);
  }

  else 
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecr_par_entrelacee :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Ecriture parallele de blocs de donnees 
 * Les blocs doivent etre rangés par ordre croissant suivant la numérotation
 * des processus
 *
 * parameters :
 *   unite             <-- Unite du fichier
 *   t_n_composantes   <-- Type de tailles composantes 
 *                        (PDM_IO_N_COMPOSANTE_CONSTANT
 *                     ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes     <-- Nombre de composantes pour chaque donnee         
 *   taille_donnee     <-- Taille unitaire de la donnnee
 *   debut_bloc        <-- Adresse relative du debut de bloc
 *   n_donnees         <-- Nombre de donnees a lire
 *   donnees           <-- Donnees a ecrire
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_ecr_par_bloc, PDM_IO_ECR_PAR_BLOC)
(const PDM_l_num_t  *unite,
 const int             *t_n_composantes,         
 const PDM_l_num_t  *n_composantes,         
 const PDM_l_num_t  *taille_donnee,
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *debut_bloc,
 const void            *donnees
)
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;

  PDM_io_ecr_par_bloc(*unite,
                        _t_n_composantes,
                        n_composantes,
                        *taille_donnee,
                        *n_donnees,
                        *debut_bloc,
                        donnees);
}

void PDM_io_ecr_par_bloc
(const PDM_l_num_t           unite,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t          *n_composantes,         
 const PDM_l_num_t           taille_donnee,
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t          debut_bloc,
 const void                    *donnees
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  unsigned char* buffer = NULL;

  if (fichier != NULL) {
      
    if (fichier->fmt_t == PDM_IO_FMT_TXT) {
      PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_io_ecr_par_bloc :\n"
              "Format text non traite\n");
      abort();
    }

    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_distribution = fichier->timer_distribution;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;
      
    PDM_timer_resume(timer_total);

    if (fichier->n_rangs == 1) {

      int l_donnees = 0;
      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        l_donnees = 0;
        for (int i = 0; i < n_donnees; i++)
          l_donnees += n_composantes[i];
      }
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        const int _n_composantes = *n_composantes;
        l_donnees = _n_composantes * n_donnees;
      }

      PDM_timer_hang_on(timer_total);
        
      PDM_io_ecriture_globale(unite,
                                taille_donnee,
                                l_donnees,
                                donnees);

      PDM_timer_resume(timer_total);

    }

    else {

      PDM_timer_resume(timer_distribution);
        
      PDM_l_num_t n_donnees_bloc_actif = 0;
      PDM_g_num_t debut_bloc_actif    = 0;

      /*---------------------------------------------------------- 
       *  Determination des rangs actifs qui accedent reellement
       *  aux fichiers et du nombre de donnees traitees par
       *  chaque rang 
       *----------------------------------------------------------*/
        
      PDM_g_num_t *n_donnees_traitees_rangs = (PDM_g_num_t *) 
        malloc(sizeof(PDM_g_num_t) * (fichier->n_rangs + 1));
        
      /* int *rangs_actifs = NULL; */
      /* int  rang_actif = 0; */
      /* int  n_rang_actif = 0; */
      int  n_donnees_traitees_rang_min = 0;
      int  n_donnees_traitees_rang_max = 0;
        
      PDM_g_num_t _id_max = n_donnees;
      PDM_g_num_t _id_max_max = 0;

      PDM_MPI_Allreduce(&_id_max,
                    &_id_max_max, 
                    1, 
                    PDM__PDM_MPI_G_NUM, 
                    PDM_MPI_SUM, 
                    fichier->comm);

      _n_donnees_rang(fichier,
                      _id_max_max,
                      n_donnees_traitees_rangs,
                      &n_donnees_traitees_rang_min,
                      &n_donnees_traitees_rang_max);

      if (fichier->n_rangs_actifs != fichier->n_rangs) {
          
        /*------------------------------------------------------------
         * Repartition des donnees sur les processus actifs
         *------------------------------------------------------------ */

        int *n_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                                  fichier->n_rangs);
        int *n_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                   fichier->n_rangs);
        int *i_donnees_a_envoyer = (int *) malloc(sizeof(int) * 
                                                  fichier->n_rangs);
        int *i_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                   fichier->n_rangs);
          
        _calcul_parametres_distribution_bloc(fichier,
                                             t_n_composantes,
                                             n_composantes,
                                             debut_bloc,
                                             n_donnees,
                                             fichier->rang_actif,
                                             n_donnees_traitees_rangs,
                                             n_donnees_a_envoyer,
                                             i_donnees_a_envoyer,
                                             n_donnees_a_recevoir,
                                             i_donnees_a_recevoir);
  
         
        /*------------------------------------------------------------
         * Calcul de la longueur du bloc et l'adresse de celui-ci
         * dans le fichier pour le processus courant
         *------------------------------------------------------------ */

        n_donnees_bloc_actif = i_donnees_a_recevoir[fichier->n_rangs-1] +
          n_donnees_a_recevoir[fichier->n_rangs-1];

        int *n_donnees_blocs_actifs = 
          (int *) malloc(sizeof(int) * fichier->n_rangs);
          
        PDM_MPI_Allgather(&n_donnees_bloc_actif, 1, PDM_MPI_INT, 
                      n_donnees_blocs_actifs, 1, PDM_MPI_INT, 
                      fichier->comm);
          
        debut_bloc_actif = 0;
          
        for (int i = 0; i < fichier->rang; i++) {
          debut_bloc_actif += n_donnees_blocs_actifs[i];
        }

        free(n_donnees_blocs_actifs);

        /*------------------------------------------------------------
         * Prise en compte de la taille de la donnee
         *------------------------------------------------------------ */
 
        for(int i = 0; i < fichier->n_rangs; i++) {
          n_donnees_a_envoyer[i]  = n_donnees_a_envoyer[i] * taille_donnee;
          i_donnees_a_envoyer[i]  = i_donnees_a_envoyer[i] * taille_donnee;
          n_donnees_a_recevoir[i] = n_donnees_a_recevoir[i] * taille_donnee;
          i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i] * taille_donnee;
        }  

        /*------------------------------------------------------------
         * Envoi des donnees
         *------------------------------------------------------------ */
          
        buffer = malloc(sizeof(unsigned char) *  
                        (n_donnees_a_recevoir[fichier->n_rangs-1] + 
                         i_donnees_a_recevoir[fichier->n_rangs-1])); 

        if (0 == 1) {
          PDM_printf("distribution ecr: %i ", fichier->rang);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   n_donnees_a_envoyer[i]);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   i_donnees_a_envoyer[i]);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   n_donnees_a_recevoir[i]);
          PDM_printf("/");
          for (int i = 0; i <  fichier->n_rangs; i++)
            PDM_printf(" %i ",
                   i_donnees_a_recevoir[i]);
          PDM_printf("\n");
          abort();
        }

        PDM_MPI_Alltoallv((void *) donnees,
                      n_donnees_a_envoyer,
                      i_donnees_a_envoyer,
                      PDM_MPI_BYTE, 
                      buffer,
                      n_donnees_a_recevoir,
                      i_donnees_a_recevoir,
                      PDM_MPI_BYTE, 
                      fichier->comm);
      
        free(n_donnees_a_envoyer);
        free(i_donnees_a_envoyer);
          
        free(n_donnees_a_recevoir);
        free(i_donnees_a_recevoir);

      }
        
      else {

        /*------------------------------------------------------------
         * Tous les blocs sont actifs. Pas de deplacement de donnees
         * Il faut juste mettre a jour l'index de debut de bloc et le
         * nombre de donnees en fonction du nombre de composantes
         *------------------------------------------------------------ */
        
        buffer = (unsigned char *) donnees;

        if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
          n_donnees_bloc_actif = 0;
          for (int i = 0; i < n_donnees; i++) {
            n_donnees_bloc_actif += n_composantes[i];
          }

          int *n_donnees_blocs_actifs = 
            (int *) malloc(sizeof(int) * fichier->n_rangs);
        
          PDM_MPI_Allgather(&n_donnees_bloc_actif, 1, PDM_MPI_INT, 
                        n_donnees_blocs_actifs, 1, PDM_MPI_INT, 
                        fichier->comm);

          debut_bloc_actif = 0;

          for (int i = 0; i < fichier->rang; i++) {
            debut_bloc_actif += n_donnees_blocs_actifs[i];
          }
            
          free(n_donnees_blocs_actifs);
        }

        else {

          debut_bloc_actif = (debut_bloc - 1) * n_composantes[0];
          n_donnees_bloc_actif = n_donnees * n_composantes[0];

        }
      }

      free(n_donnees_traitees_rangs);
 
      PDM_timer_hang_on(timer_distribution);

      /*---------------------------------------------------------- 
       *  Ecriture du buffer
       *----------------------------------------------------------*/
        
      PDM_timer_resume(timer_fichier);

      switch (fichier->acces) {
          
        /* Ecriture parallele des blocs */
          
      case PDM_IO_ACCES_MPIIO_EO:
      case PDM_IO_ACCES_MPIIO_IP:
        if (fichier->rang_actif) {
          
          PDM_file_par_ecriture_parallele(fichier->PDM_file_par,
                                          taille_donnee,
                                          n_donnees_bloc_actif,
                                          buffer,
                                          debut_bloc_actif);
        }
        break;
          
        /* Ecriture sequentielle des blocs provenant des rangs actifs */
          
      case PDM_IO_ACCES_MPI_SIMPLE:       
        {
            
          int *n_donnees_blocs_actifs = 
            (int *) malloc(sizeof(int) * fichier->n_rangs);
        
          PDM_MPI_Allgather(&n_donnees_bloc_actif, 1, PDM_MPI_INT, 
                        n_donnees_blocs_actifs, 1, PDM_MPI_INT, 
                        fichier->comm);

          int max_n_donnees_blocs_actif = 0;
          for (int i = 0; i < fichier->n_rangs; i++)
            max_n_donnees_blocs_actif = PDM_IO_MAX(max_n_donnees_blocs_actif,
                                                     n_donnees_blocs_actifs[i]);

          int etat_ecriture = 1;

          if (fichier->rang == 0) {
            assert(fichier->rang_actif == 1);

            /* Ecriture du buffer du proc maitre */
              
	    PDM_g_num_t n_donnees_blocs_tmp = n_donnees_blocs_actifs[0];
            PDM_g_num_t donnees_ecrites = PDM_file_seq_write(fichier->PDM_file_seq,
                                                    taille_donnee,
                                                    n_donnees_blocs_tmp,
                                                    buffer);
            if (buffer == donnees)
              buffer = (unsigned char *) malloc (sizeof(unsigned char) * 
                                                 max_n_donnees_blocs_actif * taille_donnee);
            else
              buffer = (unsigned char *) realloc (buffer, 
                                                  sizeof(unsigned char) * 
                                                  max_n_donnees_blocs_actif * taille_donnee);
              
            if (donnees_ecrites != n_donnees_blocs_tmp)
              etat_ecriture = 0;

            /* Reception des buffers des autres processus actifs 
               et ecriture */
              
            for (int i = 1; i < fichier->n_rangs_actifs; i++) {
              int l_buffer = n_donnees_blocs_actifs[fichier->rangs_actifs[i]] * taille_donnee;
              PDM_MPI_Recv(buffer, l_buffer, PDM_MPI_BYTE, fichier->rangs_actifs[i], 
                       PDM_io_tag, fichier->comm);
                
	      n_donnees_blocs_tmp = n_donnees_blocs_actifs[fichier->rangs_actifs[i]];
              donnees_ecrites = 
                PDM_file_seq_write(fichier->PDM_file_seq,
                                  taille_donnee,
                                  n_donnees_blocs_tmp,
                                  buffer);
              if (donnees_ecrites != n_donnees_blocs_tmp)
                etat_ecriture = 0;
                
            }

          }
          else if (fichier->rang_actif == 1) {

            /* Envoi du buffer au processu maitre si actif */

            int l_buffer = n_donnees_bloc_actif * taille_donnee;
            PDM_MPI_Send(buffer, l_buffer, PDM_MPI_BYTE, 0, 
                     PDM_io_tag, fichier->comm);
          }

          PDM_MPI_Bcast(&etat_ecriture, 1, PDM_MPI_INT, 0, fichier->comm);
            
          if (etat_ecriture == 0) {
            PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecr_par_bloc :"
                    " Erreur d'ecriture du fichier '%s' \n", fichier->nom);
            abort();
          }

          free(n_donnees_blocs_actifs);

          break;
        }
      default:
	break;
      }
      
      PDM_timer_hang_on(timer_fichier);

      /* Liberation memoire */
        
      PDM_timer_resume(timer_distribution);
        
      PDM_timer_hang_on(timer_distribution);
    }

    if (buffer != donnees)
      free(buffer);
      
    PDM_timer_hang_on(timer_total);
  }

  else 
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_ecr_par_bloc :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Fermeture du fichier sans destruction de la structure PDM_io associee a
 * l'unite
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_close, PDM_IO_CLOSE)
(const PDM_l_num_t *unite
)
{
  PDM_io_close(*unite);
}

void PDM_io_close
(const PDM_l_num_t unite
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {

    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_fichier = fichier->timer_fichier;

    PDM_timer_resume(timer_total);
    PDM_timer_resume(timer_fichier);
      
    if (fichier->PDM_file_seq != NULL)  {
      PDM_file_seq_close(fichier->PDM_file_seq);
      free(fichier->PDM_file_seq);
      fichier->PDM_file_seq = NULL;
    }
      
    if (fichier->PDM_file_par != NULL) {
      PDM_file_par_close(fichier->PDM_file_par);
      free(fichier->PDM_file_par);
      fichier->PDM_file_par = NULL;
    }
     
    if (fichier->fmt != NULL) {
      free(fichier->fmt);
      fichier->fmt = NULL;
    }

    /* suppression du fichier backup */

    if (fichier->backup == PDM_IO_BACKUP_ON) {

      int proc_actif = 0;
      if (fichier->acces != PDM_IO_ACCES_SEQ) {
        if (fichier->rang == 0) {
          proc_actif = 1;
        }
      }
      else {
        proc_actif = 1;
      } 
      
      if (proc_actif == 1) {

        char *fichier_backup = malloc(sizeof(char) * (strlen(fichier->nom) + 2)); /* caratère ~ + \0 */
        strcpy(fichier_backup, fichier->nom);
        strcat(fichier_backup, "~");
        
        FILE *testbackup = fopen(fichier_backup, "r");
        if (testbackup != NULL) {
          fclose(testbackup);
          remove(fichier_backup);
        }

        free(fichier_backup);
      }
    }

    PDM_timer_hang_on(timer_fichier);
    PDM_timer_hang_on(timer_total);

  }

  else {
    err_code = 1;
  }
    
  PDM_MPI_Barrier (fichier->comm);
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_close : unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Destruction de la structure PDM_io associee a l'unite
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_detruit, PDM_IO_DETRUIT)
(const PDM_l_num_t *unite
)
{
  PDM_io_detruit(*unite);
}

void PDM_io_detruit
(const PDM_l_num_t unite
)
{
  /* Fermeture du fichier */

  PDM_io_close(unite);

  /* Liberation de la structure */

  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
    
    free(fichier->nom);
    
    PDM_timer_free(fichier->timer_fichier);
    
    PDM_timer_free(fichier->timer_distribution);
    
    PDM_timer_free(fichier->timer_swap_endian);

    PDM_timer_free(fichier->timer_total);
    
    if (fichier->rangs_actifs != NULL) {
      free(fichier->rangs_actifs);
      fichier->rangs_actifs = NULL;
    }

    if (fichier->rangs_inactifs != NULL) {
      free(fichier->rangs_inactifs);
      fichier->rangs_inactifs = NULL;
    }

    if (fichier->tag_rangs_actifs != NULL) {
      free(fichier->tag_rangs_actifs);
      fichier->tag_rangs_actifs = NULL;
    }

    free(fichier);

    PDM_Handles_handle_free (PDM_io_fichiers, unite, PDM_FALSE);
    
    int n_file = PDM_Handles_n_get (PDM_io_fichiers);  
    
    if (n_file == 0) {
      PDM_io_fichiers = PDM_Handles_free (PDM_io_fichiers); 
    }
  }

  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_detruit : unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Retourne le temps cumule d'acces aux fichiers
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   t_cpu           --> Temps CPU
 *   t_elapsed       --> Temps elapsed
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_get_timer_fichier, PDM_IO_GET_TIMER_FICHIER)
(const PDM_l_num_t *unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  PDM_io_get_timer_fichier(*unite, t_cpu, t_elapsed);
}

void PDM_io_get_timer_fichier
(const PDM_l_num_t  unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
    PDM_timer_t *timer = fichier->timer_fichier;
    *t_cpu = PDM_timer_cpu(timer);
    *t_elapsed = PDM_timer_elapsed(timer);
  }

  else 
    err_code = 1;

  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_get_timer_fichier"
            " : unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Retourne le temps cumule pour la distribution des donnees
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   t_cpu           --> Temps CPU
 *   t_elapsed       --> Temps elapsed
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_get_timer_distrib, PDM_IO_GET_TIMER_DISTRIB)
(const PDM_l_num_t *unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  PDM_io_get_timer_distrib(*unite, t_cpu, t_elapsed);
}

void PDM_io_get_timer_distrib
(const PDM_l_num_t unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
      
    PDM_timer_t *timer = fichier->timer_distribution;
    *t_cpu = PDM_timer_cpu(timer);
    *t_elapsed = PDM_timer_elapsed(timer);
  }

  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_get_timer_distribution"
            " : unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Retourne le temps cumule pour le swap des donnees
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   t_cpu           --> Temps CPU
 *   t_elapsed       --> Temps elapsed
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_get_timer_swap_endian, PDM_IO_GET_TIMER_SWAP_ENDIAN)
(const PDM_l_num_t *unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  PDM_io_get_timer_swap_endian(*unite, t_cpu, t_elapsed);
}

void PDM_io_get_timer_swap_endian
(const PDM_l_num_t unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
      
    PDM_timer_t *timer = fichier->timer_swap_endian;
    *t_cpu = PDM_timer_cpu(timer);
    *t_elapsed = PDM_timer_elapsed(timer);
  }

  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_get_timer_distribution"
            " : unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Retourne le temps cumule total
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *   t_cpu           --> Temps CPU
 *   t_elapsed       --> Temps elapsed
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_get_timer_total, PDM_IO_GET_TIMER_TOTAL)
(const PDM_l_num_t *unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  PDM_io_get_timer_total(*unite, t_cpu, t_elapsed);
}

void PDM_io_get_timer_total
(const PDM_l_num_t unite,
 double               *t_cpu,
 double               *t_elapsed
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
      
    PDM_timer_t *timer = fichier->timer_total;
    *t_cpu = PDM_timer_cpu(timer);
    *t_elapsed = PDM_timer_elapsed(timer);
    
  }
  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_get_timer_total"
            " : unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Retourne les informations sur le fichire
 *
 * parameters :
 *   unite           <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_dump, PDM_IO_DUMP)
(const PDM_l_num_t *unite
)
{
  PDM_io_dump(*unite);
}

void PDM_io_dump
(const PDM_l_num_t unite
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {
    PDM_printf("Propriete du fichier d'unite '%i'\n", unite);
    PDM_printf("   - nom                           : %s\n", fichier->nom);
    PDM_printf("   - mode                          : ");
    if (fichier->mode == PDM_IO_MODE_LECTURE)
      PDM_printf("PDM_io_mode_lecture\n");
    else if (fichier->mode == PDM_IO_MODE_ECRITURE)
      PDM_printf("PDM_io_mode_ecriture\n");
    else if (fichier->mode == PDM_IO_MODE_AJOUT)
      PDM_printf("PDM_io_mode_ajout\n");
    PDM_printf("   - acces                         : ");
    if (fichier->acces == PDM_IO_ACCES_MPIIO_EO)
      PDM_printf("PDM_io_acces_mpiio_eo\n");
    else if (fichier->acces == PDM_IO_ACCES_MPIIO_IP)
      PDM_printf("PDM_io_acces_mpiio_ip\n");
    else if (fichier->acces == PDM_IO_ACCES_MPI_SIMPLE)
      PDM_printf("PDM_io_acces_mpi_simple\n");
    else if (fichier->acces == PDM_IO_ACCES_SEQ)
      PDM_printf("PDM_io_acces_seq\n");
    
    PDM_printf("   - swap_endian                   : ");
    if (fichier->swap_endian == 1)
      PDM_printf("actif\n");
    else if (fichier->swap_endian == 0)
      PDM_printf("inactif\n");
  }
  else 
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_dump :"
            " unite '%d' non valide\n", unite);
    abort();
  }
}

/*----------------------------------------------------------------------------
 * Retourne le communicateur du fichier               
 *
 * parameters :
 *   unite   <-- Unite du fichier
 *   pdm_mpi_comm--> Communicateur mpi       
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_get_comm, PDM_IO_GET_COMM)
(PDM_l_num_t *unite,
 PDM_MPI_Fint       *pdm_mpi_comm   
)
{
  PDM_MPI_Comm              comm;

  PDM_io_get_comm(*unite, &comm);
  *pdm_mpi_comm = PDM_MPI_Comm_c2f(comm);
}

void PDM_io_get_comm
(const PDM_l_num_t  unite,
 PDM_MPI_Comm             *pdm_mpi_comm
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL)  
    *pdm_mpi_comm   = fichier->comm;
  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_get_comm"   
            " : unite '%d' non valide\n", unite);
    abort();
  }
}


/*----------------------------------------------------------------------------
 * Active le swap endian
 *
 * parameters :
 *   unite   <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_swap_endian_on, PDM_IO_SWAP_ENDIAN_ON)
(
PDM_l_num_t *unite
)
{
  PDM_io_swap_endian_on(*unite);
}

void PDM_io_swap_endian_on
(
const PDM_l_num_t unite
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL)  
    fichier->swap_endian = 1;
  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_get_comm"   
            " : unite '%d' non valide\n", unite);
    abort();
  }
}


/*----------------------------------------------------------------------------
 * Désactive le swap endian
 *
 * parameters :
 *   unite   <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_swap_endian_off, PDM_IO_SWAP_ENDIAN_OFF)
(
PDM_l_num_t *unite
)
{
  PDM_io_swap_endian_off(*unite);
}

void PDM_io_swap_endian_off
(
const PDM_l_num_t unite
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL)  
    fichier->swap_endian = 0;
  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_swap_endian_off"   
            " : unite '%d' non valide\n", unite);
    abort();
  }
}


/*----------------------------------------------------------------------------
 * swap endian pour conversion little endian <-> big endian
 *
 * parameters :
 *   nom             <-- Nom du fichier
 *   longueur_nom    <-- Longueur du nom de fichier
 *   type_io         <-- Type (parallele avec mpiio, parallele sans mpiio,
 *                             sequentiel)
 *   mode            <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
 *   pdm_mpi_comm        <-- Communicateur lie au fichier
 *   unite           --> Unite du fichier
 *
 *----------------------------------------------------------------------------*/

 void PROCF (pdm_io_swap_endian, PDM_IO_SWAP_ENDIAN) 
 (
  const int          *taille_donnee,
  const PDM_g_num_t  *n_donnees, 
  const void         *donnees, 
  void               *resultats 
)
{
  size_t _taille_donnee = (size_t) *taille_donnee; 
  size_t _n_donnees = (size_t) *n_donnees; 
   
  PDM_io_swap_endian (_taille_donnee,
                      _n_donnees,
                      donnees,
                      resultats);   
}

void PDM_io_swap_endian(const size_t   taille_donnee,
                        const size_t   n_donnees,
                        const void    *donnees,
                        void          *resultats)
                         
{

  unsigned char  *presultats = (unsigned char *) resultats;
  const unsigned char  *pdonnees = (const unsigned char *) donnees;

  switch(taille_donnee) {

  case 1:
    if (resultats != donnees)
      memcpy(resultats, (void *) donnees, n_donnees);
    break;

  case 2:
    {
      uint16_t *_donnees = (uint16_t *) donnees;
      uint16_t *_resultats = (uint16_t *) resultats;
      for (size_t i = 0 ; i < n_donnees ; i++) 
        _resultats[i] = (uint16_t) (((_donnees[i] & 0x00FF) << 8) | 
          ((_donnees[i] & 0xFF00) >> 8));
      break;
    }
  case 4:
    {
      uint32_t *_donnees = (uint32_t *) donnees;
      uint32_t *_resultats = (uint32_t *) resultats;
      for (size_t i = 0 ; i < n_donnees ; i++) 
        _resultats[i] = ((_donnees[i] & 0x000000FF) << 24) | 
          ((_donnees[i] & 0x0000FF00) << 8)  | 
          ((_donnees[i] & 0x00FF0000) >> 8)  | 
          ((_donnees[i] & 0xFF000000) >> 24) ;
      break;
    }

  case 8:
    {
      uint64_t *_donnees = (uint64_t *) donnees;
      uint64_t *_resultats = (uint64_t *) resultats;
      for (size_t i = 0 ; i < n_donnees ; i++) 
        
      _resultats[i] = (0x00000000000000FFULL & (_donnees[i] >> 56)) | 
        (0x000000000000FF00ULL & (_donnees[i] >> 40)) | 
        (0x0000000000FF0000ULL & (_donnees[i] >> 24)) | 
        (0x00000000FF000000ULL & (_donnees[i] >> 8))  | 
        (0x000000FF00000000ULL & (_donnees[i] << 8))  | 
        (0x0000FF0000000000ULL & (_donnees[i] << 24)) | 
        (0x00FF000000000000ULL & (_donnees[i] << 40)) | 
        (0xFF00000000000000ULL & (_donnees[i] << 56));

      break;
    }

  default :
    {
      size_t  shift;
      unsigned char  tmpswap;
      
      for (size_t i = 0 ; i < n_donnees ; i++) {
        
        shift = i * taille_donnee;
        
        for (size_t ib = 0 ; ib < (taille_donnee / 2) ; ib++) {
          
          tmpswap = *(pdonnees + shift + ib);
          *(presultats + shift + ib) = *(pdonnees + shift + 
                                         (taille_donnee - 1) - ib);
          *(presultats + shift + (taille_donnee - 1) - ib) = tmpswap;
          
        }
      }
    }
  } /* switch */
}


/*----------------------------------------------------------------------------
 * Définit le format de la donnée indviduelle pour la sortie text
 *
 * parameters :
 *   unite   <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_fmt_donnee_set_cf, PDM_IO_FMT_DONNEE_SET_CF)
(
 const PDM_l_num_t *unite,
 const PDM_l_num_t *n_char_fmt,
 const PDM_l_num_t *data_type,
 const char           *fmt,
 const PDM_l_num_t *l_fmt
 ARGF_SUPP_CHAINE
)
{
  char *fmt_c    = PDM_fortran_to_c_string(fmt, *l_fmt);

  PDM_io_fmt_donnee_set(*unite, *n_char_fmt, (PDM_io_type_t) *data_type, fmt_c);
  
  free(fmt_c);
}

void PDM_io_fmt_donnee_set
(
 const PDM_l_num_t unite,
 const PDM_l_num_t n_char_fmt,
 const PDM_io_type_t data_type,
 const char           *fmt
)
{
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) { 
    fichier->swap_endian = 0;
    if (fichier->fmt != NULL)
      free(fichier->fmt);
    fichier->fmt        = (char *) malloc(sizeof(char) * (strlen(fmt) + 1));
    strcpy(fichier->fmt, fmt);
    fichier->n_char_fmt = n_char_fmt;
    fichier->data_type  = data_type;
  }
  else
    err_code = 1;
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_fmt_donnee_set"   
            " : unite '%d' non valide\n", unite);
    abort();
  }
}


/*----------------------------------------------------------------------------
 * Creation d'un directory
 *
 * parameters :
 *   unite   <-- Unite du fichier
 *
 *----------------------------------------------------------------------------*/

int PDM_io_mkdir
(
const char* path
)
{
  char *tmp_path = (char *) malloc((strlen(path) + 1)*sizeof(char));
  strcpy(tmp_path, path);
  int idx = 0;
  size_t _l_path = strlen(path);
  int l_path = (int) _l_path;
  int err = 0;

  if ((l_path > 1) && tmp_path[l_path-1] == '/')
    tmp_path[l_path - 1] = '\0';
    
  while(1) {
    while((tmp_path[idx] != '/') && (tmp_path[idx] != '\0')) {
      idx += 1;
    }
    if (tmp_path[idx] == '/') {
      tmp_path[idx] = '\0';
      err = _mkdir(tmp_path);
      if (err != 0) {
        free(tmp_path);
        return 0;
      }
      tmp_path[idx] = '/';
    }
    else {
      err = _mkdir(tmp_path);
      break;
    }
    idx += 1;
  }
  free(tmp_path);
  return err;
}

/*----------------------------------------------------------------------------
 * Calcul de la taille totale d'un champ de donnees
 *
 * parameters :
 *   unite             <-- Unite du fichier
 *   t_n_composantes   <-- Type de tailles composantes 
 *                        (PDM_IO_N_COMPOSANTE_CONSTANT
 *                     ou PDM_IO_N_COMPOSANTE_VARIABLE)
 *   n_composantes     <-- Nombre de composantes pour chaque donnee         
 *   n_donnees         <-- Nombre de donnees a lire
 *   indirection       <-- Indirection de redistribition des donnees
 *                       Attention cet argument est un int64
 *   t_n_donnee        --> Nombre total de donnees (Elimination des doublons)
 *  
 *----------------------------------------------------------------------------*/

void PROCF (pdm_io_n_donnees_get, PDM_IO_N_DONNEES_GET)

(const PDM_l_num_t  *unite,
 const int             *t_n_composantes,         
 const PDM_l_num_t  *n_composantes,         
 const PDM_l_num_t  *n_donnees,
 const PDM_g_num_t *indirection,
       PDM_g_num_t *t_n_donnees
)
{
  PDM_io_n_composantes_t _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;

  if (*t_n_composantes == 0) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_CONSTANT;
  else if (*t_n_composantes == 1) 
    _t_n_composantes = PDM_IO_N_COMPOSANTE_VARIABLE;

  *t_n_donnees = PDM_io_n_donnees_get (*unite,
                                         _t_n_composantes,
                                         n_composantes,
                                         *n_donnees,
                                         indirection);
}

PDM_g_num_t
PDM_io_n_donnees_get
(const PDM_l_num_t           unite,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t          *n_composantes,         
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t         *indirection
)
{
 
  PDM_g_num_t t_n_donnees = 0;
  
  int err_code = 0;
  PDM_io_fichier_t *fichier = PDM_io_get_fichier(unite);

  if (fichier != NULL) {

    
    PDM_timer_t *timer_total = fichier->timer_total;
    PDM_timer_t *timer_distribution = fichier->timer_distribution;
     
    PDM_timer_resume(timer_total);
      
    if (fichier->n_rangs == 1) {

      PDM_timer_resume(timer_distribution);
      
      int _n_donnees = 0;

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        
        for (int i = 0; i < n_donnees; i++) {
          _n_donnees += n_composantes[i];
        }
        
      }
        
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
        _n_donnees = n_donnees * n_composantes[0];
      }
      
      PDM_timer_hang_on(timer_distribution);

      t_n_donnees = (PDM_g_num_t) _n_donnees;

    }

    else {
        
      PDM_timer_resume(timer_distribution);
        
      /*---------------------------------------------------------- 
       *  Determination des rangs actifs qui accedent reellement
       *  aux fichiers et du nombre de donnee traitee par
       *  chaque rang 
       *----------------------------------------------------------*/
        
      PDM_g_num_t *n_donnees_rangs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (fichier->n_rangs + 1));
        
      int  n_donnees_rang_min = 0;
      int  n_donnees_rang_max = 0;
        
      PDM_g_num_t _id_max = 0;
      PDM_g_num_t _id_max_max = 0;

      for (int i = 0; i < n_donnees; i++) {
        _id_max = PDM_IO_MAX(_id_max, indirection[i]);
      }

      PDM_MPI_Allreduce(&_id_max,
                    &_id_max_max, 
                    1, 
                    PDM_MPI_LONG, 
                    PDM_MPI_MAX, 
                    fichier->comm);

      if (t_n_composantes == PDM_IO_N_COMPOSANTE_VARIABLE) {
        
        _n_donnees_rang(fichier,
                        _id_max_max,
                        n_donnees_rangs,
                        &n_donnees_rang_min,
                        &n_donnees_rang_max);

        /*---------------------------------------
         *  Envoi/Reception des numeros absolus 
         *  decrits par l'indirection 
         *---------------------------------------*/
      
        int *n_donnees_a_envoyer = (int *) malloc(sizeof(int) * fichier->n_rangs);
        int *n_donnees_a_recevoir = (int *) malloc(sizeof(int) * fichier->n_rangs);
        
        /* Pour chaque donnee le proc ou elle va etre envoyee */
   
        int *donnees_proc = (int *) malloc(sizeof(int) * n_donnees); 
        int *indirection_locale = (int *) malloc(sizeof(int) * n_donnees);
        
        /* Calcul du nombre de donnees a envoyer a chaque procesus */
        
        for (int i = 0; i < fichier->n_rangs; i++)
          n_donnees_a_envoyer[i] = 0;
        
        for (int i = 0; i < n_donnees; i++) {
          
          /* Recherche du processus le plus probable */
        
          PDM_g_num_t n_absolu = indirection[i] - 1;
          int irang_actif = fichier->n_rangs_actifs - 1;

          if (n_donnees_rang_min > 0) { 
            irang_actif = PDM_IO_MIN((int) (n_absolu / 
                                              (PDM_g_num_t) n_donnees_rang_min), 
                                       fichier->n_rangs_actifs - 1) ;
          }
        
          /* Ajustement */
   
          while (n_absolu < n_donnees_rangs[fichier->rangs_actifs[irang_actif]])
            irang_actif -= 1;

          assert(n_absolu < (n_donnees_rangs[fichier->rangs_actifs[irang_actif] + 1])); 
          
          n_donnees_a_envoyer[fichier->rangs_actifs[irang_actif]] += 1;
          donnees_proc[i] = fichier->rangs_actifs[irang_actif];
        }
        
        PDM_MPI_Alltoall(n_donnees_a_envoyer,  1, PDM_MPI_INT, 
                     n_donnees_a_recevoir, 1, PDM_MPI_INT, 
                     fichier->comm);

        int *i_donnees_a_envoyer  = (int *) malloc(sizeof(int) * 
                                                   fichier->n_rangs);
        int *i_donnees_a_recevoir = (int *) malloc(sizeof(int) * 
                                                   fichier->n_rangs);
        
        i_donnees_a_envoyer[0] = 0;
        for (int i = 1; i < fichier->n_rangs; i++)
          i_donnees_a_envoyer[i] = i_donnees_a_envoyer[i-1] + 
            n_donnees_a_envoyer[i-1];
       
        i_donnees_a_recevoir[0] = 0;
        for (int i = 1; i < fichier->n_rangs; i++)
          i_donnees_a_recevoir[i] = i_donnees_a_recevoir[i-1] + 
            n_donnees_a_recevoir[i-1];
      
        PDM_g_num_t *num_absolue_envoyee = 
          (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                     (i_donnees_a_envoyer[fichier->n_rangs - 1] + 
                                      n_donnees_a_envoyer[fichier->n_rangs - 1]));
 
        for (int i = 0; i < fichier->n_rangs; i++)
          n_donnees_a_envoyer[i] = 0;
      
        for (int i = 0; i < n_donnees; i++) {
          int iproc = donnees_proc[i];
          num_absolue_envoyee[i_donnees_a_envoyer[iproc]+
                              n_donnees_a_envoyer[iproc]] = indirection[i];
          indirection_locale[i_donnees_a_envoyer[iproc] + 
                             n_donnees_a_envoyer[iproc]]  = i;
          n_donnees_a_envoyer[iproc] += 1;
        }      

        int _n_quantites = i_donnees_a_recevoir[fichier->n_rangs - 1] +
          n_donnees_a_recevoir[fichier->n_rangs - 1];

        PDM_g_num_t *num_absolue_recues = 
          (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_quantites);

        PDM_MPI_Alltoallv(num_absolue_envoyee, 
                      n_donnees_a_envoyer,
                      i_donnees_a_envoyer,
                      PDM_MPI_LONG,
                      num_absolue_recues,
                      n_donnees_a_recevoir,
                      i_donnees_a_recevoir,
                      PDM_MPI_LONG,
                      fichier->comm);

        free (num_absolue_envoyee);
        
        int *n_composantes_recues  = NULL;
        int *n_composantes_envoyee = NULL;
        
        unsigned char* donnees_alltoall = NULL;
        

        int l_data_recv = n_donnees_a_recevoir[fichier->n_rangs-1] + i_donnees_a_recevoir[fichier->n_rangs-1];
        int l_bloc = n_donnees_rangs[fichier->rang+1] - n_donnees_rangs[fichier->rang];
      
        int *n_comp_bloc = (int *) malloc(sizeof(int) * l_bloc);

        for (int i = 0; i < l_bloc; i++) {
          n_comp_bloc[i] = 0;
        }
        
        /*------------------------------------------------------------
         * Envoi/reception du nombre de composantes
         *------------------------------------------------------------ */

        n_composantes_envoyee = 
          (int *) malloc(sizeof(int) * 
                         (i_donnees_a_envoyer[fichier->n_rangs - 1] + 
                          n_donnees_a_envoyer[fichier->n_rangs - 1]));
        
        for (int i = 0; i < fichier->n_rangs; i++)
          n_donnees_a_envoyer[i] = 0;
          
        for (int i = 0; i < n_donnees; i++) {
          int iproc = donnees_proc[i];
          n_composantes_envoyee[i_donnees_a_envoyer[iproc] + 
                                n_donnees_a_envoyer[iproc]] = 
            n_composantes[i];
          n_donnees_a_envoyer[iproc] += 1;
        }      

        n_composantes_recues = (int *) malloc(sizeof(int) * l_data_recv);

        PDM_MPI_Alltoallv(n_composantes_envoyee,
                      n_donnees_a_envoyer,
                      i_donnees_a_envoyer,
                      PDM_MPI_INT, 
                      n_composantes_recues,
                      n_donnees_a_recevoir,
                      i_donnees_a_recevoir,
                      PDM_MPI_INT, 
                      fichier->comm);

        free(n_composantes_envoyee);
          
        /*---------------------------------------------------------- 
         *  Rangement des donnees pour l'echange croise
         *  (regroupement des donnees par rang)
         *----------------------------------------------------------*/
          
        for (int i = 0; i < l_data_recv; i++) {
          int j = (int) (num_absolue_recues[i] - n_donnees_rangs[fichier->rang]) - 1;
          n_comp_bloc[j] = n_composantes_recues[i];
        }

        PDM_g_num_t sum_n_comp = 0;
        for (int i = 0; i < l_bloc; i++) {
          sum_n_comp += n_comp_bloc[i];
        }

        PDM_MPI_Allreduce(&sum_n_comp, &t_n_donnees, 1,
                      PDM_MPI_LONG, PDM_MPI_SUM, fichier->comm); 

        
        free(donnees_proc);
        free(indirection_locale);
        
        free(donnees_alltoall);
        free(n_donnees_a_envoyer);
        free(i_donnees_a_envoyer);
        
        free(n_donnees_a_recevoir);
        free(i_donnees_a_recevoir);

        free (n_comp_bloc);
        free (num_absolue_recues);
        free (n_composantes_recues);
        
      }
      
      else if (t_n_composantes == PDM_IO_N_COMPOSANTE_CONSTANT) {
  
        const int _n_composantes = n_composantes[0];
          
        t_n_donnees = _n_composantes * _id_max_max;

          
      }

      PDM_timer_hang_on(timer_distribution);
  
      free (n_donnees_rangs);

    }

    PDM_timer_hang_on(timer_total);

  }

  else { 
    err_code = 1;
  }
  
  if (err_code){
    PDM_error(__FILE__, __LINE__, 0,"Erreur PDM_io_n_donnees_get :"
            " unite '%d' non valide\n", unite);
    abort();
  }

  return t_n_donnees;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
