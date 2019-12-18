#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_multicoloring.h"
#include "pdm_distributed_graph_rank1_from_part.h"

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
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *nVtxSeg,
           double        *length,
           int           *n_part,
     int           *post,
     int           *method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nVtxSeg = atol(argv[i]);
        *nVtxSeg = (PDM_g_num_t) _nVtxSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        nVtxSeg = 4;
  double             length  = 1.;
  int                nPart   = 2;
  int                post    = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nVtxSeg,
             &length,
             &nPart,
             &post,
             (int *) &method);

  /*
   *  Init
   */

  int myRank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  int          dNCell;
  int          dNFace;
  int          dNVtx;
  int          nFaceGroup;
  PDM_g_num_t *dFaceCell = NULL;
  int         *dFaceVtxIdx = NULL;
  PDM_g_num_t *dFaceVtx = NULL;
  double      *dVtxCoord = NULL;
  int         *dFaceGroupIdx = NULL;
  PDM_g_num_t *dFaceGroup = NULL;
  int          dFaceVtxL;
  int          dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int          id;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_gen_init(&id,
                      comm,
                      nVtxSeg,
                      length,
                      0.,
                      0.,
                      0.);

  PDM_dcube_gen_dim_get(id,
                        &nFaceGroup,
                        &dNCell,
                        &dNFace,
                        &dNVtx,
                        &dFaceVtxL,
                        &dFaceGroupL);

  PDM_dcube_gen_data_get(id,
                         &dFaceCell,
                         &dFaceVtxIdx,
                         &dFaceVtx,
                         &dVtxCoord,
                         &dFaceGroupIdx,
                         &dFaceGroup);
  int ppartId = 0;

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc(dNCell*sizeof(int));
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int nPropertyCell = 0;
  int nPropertyFace = 0;

  PDM_part_create(&ppartId,
                  comm,
                  method,
                  "PDM_PART_RENUM_CELL_CUTHILL",
                  "PDM_PART_RENUM_FACE_NONE",
                  nPropertyCell,
                  renum_properties_cell,
                  nPropertyFace,
                  renum_properties_face,
                  nPart,
                  dNCell,
                  dNFace,
                  dNVtx,
                  nFaceGroup,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  have_dCellPart,
                  dCellPart,
                  dFaceCell,
                  dFaceVtxIdx,
                  dFaceVtx,
                  NULL,
                  dVtxCoord,
                  NULL,
                  dFaceGroupIdx,
                  dFaceGroup);

  // I/ Compute first rank
  _dist_csr* dCsrO1 = PDM_compute_dist_graph_rank1_from_ppart(ppartId, nPart);
  PDM_dist_csr_print(dCsrO1);
  PDM_dist_csr_dump(dCsrO1, "dcsr_O1");

  // ------------------------------------------------------------------------------------
  // Multicoloring
  PDM_compute_distributed_multicoloring(dCsrO1);

  // ------------------------------------------------------------------------------------
  // Free
  PDM_dist_csr_free(dCsrO1);
  PDM_part_free(ppartId);
  PDM_dcube_gen_free(id);


  PDM_MPI_Finalize();

  return 0;
}
