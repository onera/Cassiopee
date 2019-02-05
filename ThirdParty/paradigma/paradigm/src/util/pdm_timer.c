/*============================================================================
 * Mesure des temps CPU et elapsed
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include "pdm_config.h"

#if defined (PDM_HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#elif defined(_POSIX_SOURCE)
#include <sys/times.h>
#include <unistd.h>
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_timer.h"
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
 * Structure de mesure des temps d'execution
 *----------------------------------------------------------------------------*/

struct _pdm_timer_t {

  double  t_cpu;                   /* Temps CPU cumule */
  double  t_elapsed;               /* Temps elapsed cumule */
#if defined (PDM_HAVE_GETRUSAGE)
  double  t_cpu_u;                 /* Temps CPU utilisateur */
  double  t_cpu_s;                 /* Temps CPU system */
  double  t_cpu_debut;
  double  t_cpu_u_debut;
  double  t_cpu_s_debut;
#else
  clock_t t_cpu_debut;             /* Marque de debut de mesure 
                                      du temps CPU */
#endif
  struct timeval t_elaps_debut;    /* Marque de debut de mesure 
                                      du temps elapsed */
  int     indic;                   /* Indique si une mesure d'une tranche
                                      est en cours */ 
};

/*============================================================================
 * Definition des fonctions locales
 *============================================================================*/

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation d'un objet timer
 *
 * return
 *   timer
 *
 *----------------------------------------------------------------------------*/

PDM_timer_t *PDM_timer_create(void)
{
  PDM_timer_t *timer = (PDM_timer_t *) malloc(sizeof(PDM_timer_t));
  PDM_timer_init(timer);

  return timer;

}

/*----------------------------------------------------------------------------
 * Debut la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_init(PDM_timer_t *timer)
{
#if defined (PDM_HAVE_GETRUSAGE)
  timer->t_cpu_u = 0.;
  timer->t_cpu_s = 0.;
#endif
  timer->t_cpu = 0.;
  timer->t_elapsed = 0.;
  timer->indic = 0;
}

/*----------------------------------------------------------------------------
 * Reprend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_resume(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_reprise : \n"
            "La mesure d'une tranche est deja en cours\n");
    exit(EXIT_FAILURE);
  } 
#if defined (PDM_HAVE_GETRUSAGE)
 {
   struct rusage  usage;

   if (getrusage(RUSAGE_SELF, &usage) == 0) {
     timer->t_cpu_u_debut = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-6;
     timer->t_cpu_s_debut = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1.e-6;
     timer->t_cpu_debut   = timer->t_cpu_u_debut + timer->t_cpu_s_debut;
   }
 }
#else
  timer->t_cpu_debut = clock();
#endif
  gettimeofday(&(timer->t_elaps_debut), NULL);
  timer->indic = 1;
}

/*----------------------------------------------------------------------------
 * Suspend la mesure du temps ecoule et incremente le temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_hang_on(PDM_timer_t *timer)
{

  if (!timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_suspend : \n"
            "La mesure de temps n'a pas ete declenchee par PDM_timer_reprise\n");
    exit(EXIT_FAILURE);
  } 

  /* Recuperation du temps CPU et elaps courant */

  struct timeval t_elaps_fin; 
  gettimeofday(&t_elaps_fin, NULL);
  
  /* Ajout de la tranche mesuree au temps cumule */
  
  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) - 
                         (timer->t_elaps_debut.tv_usec + 1000000 * 
                          timer->t_elaps_debut.tv_sec);

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double tranche_elapsed_max = (double) tranche_elapsed; 
  timer->t_elapsed += tranche_elapsed_max/1000000.;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
  
#if defined (PDM_HAVE_GETRUSAGE)
 {
   struct rusage  usage;

   if (getrusage(RUSAGE_SELF, &usage) == 0) {
     timer->t_cpu_u += usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-6 - timer->t_cpu_u_debut;
     timer->t_cpu_s += usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1.e-6 - timer->t_cpu_s_debut;
     timer->t_cpu    = timer->t_cpu_u + timer->t_cpu_s;
   }
 }
#else
  clock_t t_cpu_fin = clock(); 
  double tranche_cpu = (double) (t_cpu_fin - timer->t_cpu_debut);
  double tranche_cpu_max = tranche_cpu; 
  timer->t_cpu += tranche_cpu_max/CLOCKS_PER_SEC;
#endif
  timer->indic = 0;
}

/*----------------------------------------------------------------------------
 * Retourne le temps CPU en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_cpu : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_cpu\n");
    exit(EXIT_FAILURE);
  } 
  return timer->t_cpu;
}

/*----------------------------------------------------------------------------
 * Retourne le temps CPU user en secondes (-1 si indisponible)
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu_user(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_cpu_user : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_cpu\n");
    exit(EXIT_FAILURE);
  } 
#if defined (PDM_HAVE_GETRUSAGE)
  return timer->t_cpu_u;
#else
  return -1.;
#endif
}

/*----------------------------------------------------------------------------
 * Retourne le temps CPU systeme en secondes (-1 si indisponible)
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu_sys(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_cpu_user : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_cpu\n");
    exit(EXIT_FAILURE);
  } 
#if defined (PDM_HAVE_GETRUSAGE)
  return timer->t_cpu_s;
#else
  return -1.;
#endif
}

/*----------------------------------------------------------------------------
 * Retourne le temps elaps en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_elapsed(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_elapsed : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_elapsed\n");
    exit(EXIT_FAILURE);
  } 
  return timer->t_elapsed;
}

/*----------------------------------------------------------------------------
 * Destruction d'un objet timer
 *
 * parameters :
 *   timer            <-- Timer
 *
 *----------------------------------------------------------------------------*/

void PDM_timer_free(PDM_timer_t *timer)
{
  free(timer);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
