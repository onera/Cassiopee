#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_para_octree.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *nPts,
 double        *radius,
 int           *local,
 int           *rand
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nPts = atol(argv[i]);
        *nPts = (PDM_g_num_t) _nPts;
      }
    }

    else if (strcmp(argv[i], "-radius") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *radius = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


/**
 *
 * \brief  Random value
 *
 *  \return a random double in [-1, 1]
 */

static double
_random01
(
)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / PDM_ABS (rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);

  int myRank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  PDM_g_num_t nPts   = 10;
  double radius = 10.;
  int local = 0;
  int rand = 0;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(myRank);
  }

  /* Define the number of points */

  int _nPts_l = 0;
  if (local) {
    _nPts_l = (int) nPts;
  }
  else {
    _nPts_l = (int) (nPts/numProcs);
    if (myRank < nPts%numProcs) {
      _nPts_l += 1;
    }
  }

  /* Points definition : coordinates + gnum */

  double *coords = malloc(sizeof(double) * 3 * _nPts_l);

  for (int i = 0; i < _nPts_l; i++) {
    for (int j = 0; j < 3; j++) {
      coords[3*i+j] = _random01() * radius;
    }
  }

  int id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *char_length = malloc(sizeof(double) * _nPts_l);

  for (int i = 0; i < _nPts_l; i++) {
    char_length[i] = radius * 1.e-6;
  }

  PDM_gnum_set_from_coords (id, 0, _nPts_l, coords, char_length);

  PDM_gnum_compute (id);

  PDM_g_num_t *gnum = PDM_gnum_get(id, 0);

  PDM_gnum_free (id, 1);

  /* Parallel octree */

  const int n_point_cloud = 1;
  const int depth_max = 31;
  const int points_in_leaf_max = 1;

  int id2 = PDM_para_octree_create (n_point_cloud,
                                    depth_max,
                                    points_in_leaf_max,
                                    PDM_MPI_COMM_WORLD);

  PDM_para_octree_point_cloud_set (id2, 0, _nPts_l, coords, gnum);

  PDM_para_octree_build (id2);

  //  PDM_para_octree_dump (id2);

  PDM_para_octree_dump_times (id2);

  PDM_para_octree_free (id2);

  /* Free */

  free (coords);
  free (char_length);
  free (gnum);

  if (myRank == 0) {

    PDM_printf ("\nfin Test\n");

  }

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
