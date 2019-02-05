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
#include "pdm_part.h"
#include "pdm_part_coarse_mesh.h"
#include "pdm_gnum.h"

#include "pdm_writer.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"


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
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");
     

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   nVtxSeg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t  *nVtxSeg,
 double        *length
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp (argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nVtxSeg = atol (argv[i]);
        *nVtxSeg = (PDM_g_num_t) _nVtxSeg;
      }
    }
    else if (strcmp (argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *length = atof (argv[i]);
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      nVtxSeg  Number of arguments
 * \param [in]      length   Lenght of square
 * \param [in]      nPart    Number to obtain on this processus
 * \param [in]      post     mesh export status
 * \param [in]      method   Split method
 *
 * \return ppartId  ppart identifier
 *
 */

static int
_create_split_mesh
(
 int           imesh,
 PDM_MPI_Comm  pdm_mpi_comm,
 PDM_g_num_t  nVtxSeg,
 double        length,
 int           nPart,
PDM_part_split_t           method,
 int           haveRandom,
 PDM_g_num_t   *nGFace,
 PDM_g_num_t   *nGVtx,
 PDM_g_num_t   *nGEdge,
 int           *nTPart,
 int           *nEdgeGroup
)
{
  struct timeval t_elaps_debut;

  int myRank;
  int numProcs;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &myRank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

  double        xmin = 0.;
  double        xmax = length;
  double        ymin = 0.;
  double        ymax = length;
  PDM_g_num_t  nx = nVtxSeg;
  PDM_g_num_t  ny = nVtxSeg;
  int           dNFace;
  int           dNVtx;
  int           dNEdge;
  int          *dFaceVtxIdx;
  PDM_g_num_t *dFaceVtx;
  double       *dVtxCoord;
  PDM_g_num_t *dFaceEdge;
  PDM_g_num_t *dEdgeVtx;
  PDM_g_num_t *dEdgeFace;
  int          *dEdgeGroupIdx;
  PDM_g_num_t   *dEdgeGroup;

  int           initRandom = 0;

  /*
   *  Create mesh i
   */

  if (imesh == 1) {
    nx *= 2;
    ny *= 2;
  }

  ++initRandom;

  gettimeofday(&t_elaps_debut, NULL);

  PDM_poly_surf_gen (pdm_mpi_comm,
                     xmin,
                     xmax,
                     ymin,
                     ymax,
                     haveRandom,
                     initRandom,
                     nx,       
                     ny,
                     nGFace,
                     nGVtx,  
                     nGEdge,
                     &dNVtx,
                     &dVtxCoord,
                     &dNFace,
                     &dFaceVtxIdx,
                     &dFaceVtx,
                     &dFaceEdge,
                     &dNEdge,
                     &dEdgeVtx,
                     &dEdgeFace,
                     nEdgeGroup,
                     &dEdgeGroupIdx,
                     &dEdgeGroup);
  
  // validation
  int id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, pdm_mpi_comm);
// fin validation
       
  PDM_gnum_set_from_coords (id, 0, dNVtx, dVtxCoord, NULL);
  
  PDM_gnum_compute (id);

  const PDM_g_num_t *_numabs = PDM_gnum_get (id, 0);
    
  for (int j = 0; j < dNVtx; j++) {
    PDM_printf (PDM_FMT_G_NUM" %12.5e %12.5e %12.5e\n", _numabs[j], dVtxCoord[3*j], 
                                                     dVtxCoord[3*j+1], 
                                                     dVtxCoord[3*j+2]);
  }

  PDM_gnum_free (id, 0);

  struct timeval t_elaps_fin; 

  gettimeofday (&t_elaps_fin, NULL);

  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
    (t_elaps_debut.tv_usec + 1000000 *
     t_elaps_debut.tv_sec);
  long tranche_elapsed_max = tranche_elapsed;
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double t_elapsed = (double) tranche_elapsed_max/1000000.;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
  if (myRank == 0)
    PDM_printf("[%d] Temps dans creeMaillagePolygone2D %d : %12.5e\n",
           myRank, imesh, t_elapsed);

  if (0 == 1) {

    PDM_printf ("edgegroup : ");
    for (int i = 0; i < *nEdgeGroup; i++) {
      for (int j = dEdgeGroupIdx[i]; j <  dEdgeGroupIdx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dEdgeGroup[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dfacevtx : ");
    for (int i = 0; i < dNFace; i++) {
      for (int j = dFaceVtxIdx[i]; j <  dFaceVtxIdx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dFaceVtx[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dfaceedge : ");
    for (int i = 0; i < dNFace; i++) {
      for (int j = dFaceVtxIdx[i]; j <  dFaceVtxIdx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dFaceEdge[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedgevtx : ");
    for (int i = 0; i < dNEdge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i+1]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedgeface : ");
    for (int i = 0; i < dNEdge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeFace[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeFace[2*i+1]);
      PDM_printf ("\n");
    }
  }

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc (dNFace*sizeof(int));
  int *dEdgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

  dEdgeVtxIdx[0] = 0;
  for (int i = 0; i < dNEdge; i++) {
    dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
  }

  /*
   *  Split mesh i
   */

  int ppartId;

  int nPropertyCell = 0;
  int *renum_properties_cell = NULL;
  int nPropertyFace = 0;
  int *renum_properties_face = NULL;

  PDM_part_create (&ppartId,
                   pdm_mpi_comm,
                   method,
                   "PDM_PART_RENUM_CELL_NONE",
                   "PDM_PART_RENUM_FACE_NONE",
                   nPropertyCell,
                   renum_properties_cell,
                   nPropertyFace,
                   renum_properties_face,
                   nPart,
                   dNFace,
                   dNEdge,
                   dNVtx,
                   *nEdgeGroup,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   have_dCellPart,
                   dCellPart,
                   dEdgeFace,
                   dEdgeVtxIdx,
                   dEdgeVtx,
                   NULL,
                   dVtxCoord,
                   NULL,
                   dEdgeGroupIdx,
                   dEdgeGroup);

  free (dCellPart);

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_part_time_get (ppartId,
                  &elapsed,
                  &cpu,
                  &cpu_user,
                  &cpu_sys);

  if (myRank == 0)
    PDM_printf("[%d] Temps dans ppart %d : %12.5e\n",
           myRank, imesh, elapsed[0]);

  /* Statistiques */

  int    cells_average;
  int    cells_median;
  double cells_std_deviation;
  int    cells_min;
  int    cells_max;
  int    bound_part_faces_average;
  int    bound_part_faces_median;
  double bound_part_faces_std_deviation;
  int    bound_part_faces_min;
  int    bound_part_faces_max;
  int    bound_part_faces_sum;

  PDM_part_stat_get (ppartId,
                  &cells_average,
                  &cells_median,
                  &cells_std_deviation,
                  &cells_min,
                  &cells_max,
                  &bound_part_faces_average,
                  &bound_part_faces_median,
                  &bound_part_faces_std_deviation,
                  &bound_part_faces_min,
                  &bound_part_faces_max,
                  &bound_part_faces_sum);

  /* if (myRank == 0) { */
  /*   PDM_printf ("Statistics :\n"); */
  /*   PDM_printf ("  - Number of cells :\n"); */
  /*   PDM_printf ("       * average            : %i\n", cells_average);    */
  /*   PDM_printf ("       * median             : %i\n", cells_median);    */
  /*   PDM_printf ("       * standard deviation : %12.5e\n", cells_std_deviation);    */
  /*   PDM_printf ("       * min                : %i\n", cells_min);    */
  /*   PDM_printf ("       * max                : %i\n", cells_max);    */
  /*   PDM_printf ("  - Number of faces exchanging with another partition :\n"); */
  /*   PDM_printf ("       * average            : %i\n", bound_part_faces_average);    */
  /*   PDM_printf ("       * median             : %i\n", bound_part_faces_median);    */
  /*   PDM_printf ("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);    */
  /*   PDM_printf ("       * min                : %i\n", bound_part_faces_min);    */
  /*   PDM_printf ("       * max                : %i\n", bound_part_faces_max);    */
  /*   PDM_printf ("       * total              : %i\n", bound_part_faces_sum);    */
  /* } */

  free (dVtxCoord);
  free (dFaceVtxIdx);
  free (dFaceVtx);
  free (dFaceEdge);
  free (dEdgeVtxIdx);
  free (dEdgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);

  for (int ipart = 0; ipart < nPart; ipart++) {

    int nFace;
    int nEdge;
    int nEdgePartBound;
    int nVtx;
    int nProc;
    int sFaceEdge;
    int sEdgeVtx;
    int sEdgeGroup;
    int nEdgeGroup2;
    
    PDM_part_part_dim_get (ppartId,
                           ipart,
                           &nFace,
                           &nEdge,
                           &nEdgePartBound,
                           &nVtx,
                           &nProc,
                           nTPart,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);
    
  }
  
  return ppartId;
}

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
  
  /*
   *  Set default values
   */
  
  PDM_g_num_t   nVtxSeg = 4;
  double        length  = 1.;
  int           nPart   = 1;
#ifdef PDM_HAVE_PARMETIS  
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH  
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif  
  int           haveRandom = 0;

  int           myRank;
  int           numProcs;
  
  /*
   *  Read args
   */
  
  _read_args (argc,
              argv,
              &nVtxSeg,
              &length);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Create a partitioned mesh
   */
  
  PDM_g_num_t nGFace;
  PDM_g_num_t nGVtx;
  PDM_g_num_t nGEdge;
  int imesh = 0;
  int nTPart;
  int nEdgeGroup;
  
 _create_split_mesh (imesh,
                     PDM_MPI_COMM_WORLD,
                     nVtxSeg,
                     length,
                     nPart,
                     method,
                     haveRandom,
                     &nGFace,
                     &nGVtx,
                     &nGEdge,
                     &nTPart,
                     &nEdgeGroup);
  
 PDM_MPI_Finalize ();
 
 PDM_printf ("\nfin Test\n");
 
 return 0;

}


