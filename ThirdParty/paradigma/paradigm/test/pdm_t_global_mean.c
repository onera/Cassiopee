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
#include "pdm_global_mean.h"
#include "pdm_priv.h"

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

  PDM_g_num_t        nVtxSeg = 10;
  double             length  = 1.;
  int                nPart   = 1;
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

  struct timeval t_elaps_debut;

  int myRank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  int           dNCell;
  int           dNFace;
  int           dNVtx;
  int           nFaceGroup;
  PDM_g_num_t *dFaceCell = NULL;
  int          *dFaceVtxIdx = NULL;
  PDM_g_num_t *dFaceVtx = NULL;
  double       *dVtxCoord = NULL;
  int          *dFaceGroupIdx = NULL;
  PDM_g_num_t *dFaceGroup = NULL;
  int           dFaceVtxL;
  int           dFaceGroupL;

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

  if (0 == 1) {

    PDM_printf("[%i] nFaceGroup    : %i\n", myRank, nFaceGroup);
    PDM_printf("[%i] dNCell        : %i\n", myRank, dNCell);
    PDM_printf("[%i] dNFace        : %i\n", myRank, dNFace);
    PDM_printf("[%i] dNVtx         : %i\n", myRank, dNVtx);

    PDM_printf("[%i] dFaceCell     : ", myRank);
    for (int i = 0; i < 2 * dNFace; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dFaceCell[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dFaceVtxIdx   : ", myRank);
    for (int i = 0; i < dNFace + 1; i++)
      PDM_printf(" %i", dFaceVtxIdx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dFaceVtx      : ", myRank);
    for (int i = 0; i < dFaceVtxIdx[dNFace]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dFaceVtx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dVtxCoord     : ", myRank);
    for (int i = 0; i < 3*dNVtx; i++)
      PDM_printf(" %12.5e", dVtxCoord[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dFaceGroupIdx : ", myRank);
    for (int i = 0; i < nFaceGroup + 1; i++)
      PDM_printf(" %i", dFaceGroupIdx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dFaceGroup    : ", myRank);
    for (int i = 0; i < dFaceGroupIdx[nFaceGroup]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dFaceGroup[i]);
    PDM_printf("\n");

  }
  int ppartId = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc(dNCell*sizeof(int));
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int nPropertyCell = 0;
  int nPropertyFace = 0;

  int gmid = PDM_global_mean_create (nPart, comm);

  PDM_part_create(&ppartId,
                  comm,
                  method,
                  "PDM_PART_RENUM_CELL_NONE",
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

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_part_time_get(ppartId,
                 &elapsed,
                 &cpu,
                 &cpu_user,
                 &cpu_sys);

  PDM_printf("[%i]   - elapsed total                    : %12.5e\n", myRank, elapsed[0]);
  PDM_printf("[%i]   - elapsed building graph           : %12.5e\n", myRank, elapsed[1]);
  PDM_printf("[%i]   - elapsed splitting graph          : %12.5e\n", myRank, elapsed[2]);
  PDM_printf("[%i]   - elapsed building mesh partitions : %12.5e\n", myRank, elapsed[3]);

  PDM_printf("[%i]   - cpu total                        : %12.5e\n", myRank, cpu[0]);
  PDM_printf("[%i]   - cpu building graph               : %12.5e\n", myRank, cpu[1]);
  PDM_printf("[%i]   - cpu splitting graph              : %12.5e\n", myRank, cpu[2]);
  PDM_printf("[%i]   - cpu building mesh partitions     : %12.5e\n", myRank, cpu[3]);

  PDM_printf("[%i]   - cpu_user total                   : %12.5e\n", myRank, cpu_user[0]);
  PDM_printf("[%i]   - cpu_user building graph          : %12.5e\n", myRank, cpu_user[1]);
  PDM_printf("[%i]   - cpu_user splitting graph         : %12.5e\n", myRank, cpu_user[2]);
  PDM_printf("[%i]   - cpu_user building mesh partitions: %12.5e\n", myRank, cpu_user[3]);

  PDM_printf("[%i]   - cpu_sys total                    : %12.5e\n", myRank, cpu_sys[0]);
  PDM_printf("[%i]   - cpu_sys building graph           : %12.5e\n", myRank, cpu_sys[1]);
  PDM_printf("[%i]   - cpu_sys splitting graph          : %12.5e\n", myRank, cpu_sys[2]);
  PDM_printf("[%i]   - cpu_sys building mesh partitions : %12.5e\n", myRank, cpu_sys[3]);

  struct timeval t_elaps_fin;
  gettimeofday(&t_elaps_fin, NULL);

  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
                         (t_elaps_debut.tv_usec + 1000000 *
                          t_elaps_debut.tv_sec);
  long tranche_elapsed_max = tranche_elapsed;
#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double t_elapsed = (double) tranche_elapsed_max/1000000.;
#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif

  PDM_printf("[%i]   - TEMPS DANS PART_CUBE  : %12.5e\n", myRank,  t_elapsed);

  PDM_g_num_t **cellVtxGN = malloc (sizeof(PDM_g_num_t *) * nPart);

  for (int ipart = 0; ipart < nPart; ipart++) {

    int nCell;
    int nFace;
    int nFacePartBound;
    int nVtx;
    int nProc;
    int nTPart;
    int sCellFace;
    int sFaceVtx;
    int sFaceGroup;
    int nFaceGroup2;

    PDM_part_part_dim_get(ppartId,
                       ipart,
                       &nCell,
                       &nFace,
                       &nFacePartBound,
                       &nVtx,
                       &nProc,
                       &nTPart,
                       &sCellFace,
                       &sFaceVtx,
                       &sFaceGroup,
                       &nFaceGroup2);

    int          *cellTag;
    int          *cellFaceIdx;
    int          *cellFace;
    PDM_g_num_t *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;

    PDM_part_part_val_get(ppartId,
                       ipart,
                       &cellTag,
                       &cellFaceIdx,
                       &cellFace,
                       &cellLNToGN,
                       &faceTag,
                       &faceCell,
                       &faceVtxIdx,
                       &faceVtx,
                       &faceLNToGN,
                       &facePartBoundProcIdx,
                       &facePartBoundPartIdx,
                       &facePartBound,
                       &vtxTag,
                       &vtx,
                       &vtxLNToGN,
                       &faceGroupIdx,
                       &faceGroup,
                       &faceGroupLNToGN);

    int nCellVtx = 0;
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        int face = cellFace[j] - 1;
        nCellVtx += faceVtxIdx[face+1] - faceVtxIdx[face];
      }
    }

    cellVtxGN[ipart] = malloc (sizeof(PDM_g_num_t) * nCellVtx);

  }

  for (int ipart = 0; ipart < nPart; ipart++) {

    int nCell;
    int nFace;
    int nFacePartBound;
    int nVtx;
    int nProc;
    int nTPart;
    int sCellFace;
    int sFaceVtx;
    int sFaceGroup;
    int nFaceGroup2;

    PDM_part_part_dim_get(ppartId,
                       ipart,
                       &nCell,
                       &nFace,
                       &nFacePartBound,
                       &nVtx,
                       &nProc,
                       &nTPart,
                       &sCellFace,
                       &sFaceVtx,
                       &sFaceGroup,
                       &nFaceGroup2);

    int          *cellTag;
    int          *cellFaceIdx;
    int          *cellFace;
    PDM_g_num_t *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;

    PDM_part_part_val_get(ppartId,
                       ipart,
                       &cellTag,
                       &cellFaceIdx,
                       &cellFace,
                       &cellLNToGN,
                       &faceTag,
                       &faceCell,
                       &faceVtxIdx,
                       &faceVtx,
                       &faceLNToGN,
                       &facePartBoundProcIdx,
                       &facePartBoundPartIdx,
                       &facePartBound,
                       &vtxTag,
                       &vtx,
                       &vtxLNToGN,
                       &faceGroupIdx,
                       &faceGroup,
                       &faceGroupLNToGN);

    int nCellVtx = 0;
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        int face = cellFace[j] - 1;
        for (int k = faceVtxIdx[face]; k < faceVtxIdx[face+1]; k++) {
          int iVtx = faceVtx[k] - 1;
          cellVtxGN[ipart][nCellVtx++] = vtxLNToGN[iVtx];
        }
      }
    }

    PDM_global_mean_set (gmid, ipart, nCellVtx, cellVtxGN[ipart]);

  }

  double **local_field = malloc (sizeof(double *) * nPart);
  double **local_weight = malloc (sizeof(double *) * nPart);
  double **global_mean_field_ptr = malloc (sizeof(double *) * nPart);

  for (int ipart = 0; ipart < nPart; ipart++) {

    int nCell;
    int nFace;
    int nFacePartBound;
    int nVtx;
    int nProc;
    int nTPart;
    int sCellFace;
    int sFaceVtx;
    int sFaceGroup;
    int nFaceGroup2;

    PDM_part_part_dim_get(ppartId,
                       ipart,
                       &nCell,
                       &nFace,
                       &nFacePartBound,
                       &nVtx,
                       &nProc,
                       &nTPart,
                       &sCellFace,
                       &sFaceVtx,
                       &sFaceGroup,
                       &nFaceGroup2);

    int          *cellTag;
    int          *cellFaceIdx;
    int          *cellFace;
    PDM_g_num_t *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;

    PDM_part_part_val_get(ppartId,
                       ipart,
                       &cellTag,
                       &cellFaceIdx,
                       &cellFace,
                       &cellLNToGN,
                       &faceTag,
                       &faceCell,
                       &faceVtxIdx,
                       &faceVtx,
                       &faceLNToGN,
                       &facePartBoundProcIdx,
                       &facePartBoundPartIdx,
                       &facePartBound,
                       &vtxTag,
                       &vtx,
                       &vtxLNToGN,
                       &faceGroupIdx,
                       &faceGroup,
                       &faceGroupLNToGN);

    int nCellVtx = 0;
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        int face = cellFace[j] - 1;
        nCellVtx += faceVtxIdx[face+1] - faceVtxIdx[face];
      }
    }

    local_field[ipart] = malloc (sizeof(double) * nCellVtx * 3);
    local_weight[ipart] = malloc (sizeof(double) * nCellVtx);
    global_mean_field_ptr[ipart] = malloc (sizeof(double) * nCellVtx * 3);

    nCellVtx = 0;
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        int face = cellFace[j] - 1;
        for (int k = faceVtxIdx[face]; k < faceVtxIdx[face+1]; k++) {
          int iVtx = faceVtx[k] - 1;
          for (int l = 0; l < 3; l++) {
            local_field[ipart][3*nCellVtx + l] = vtx[3*iVtx+l];
            global_mean_field_ptr[ipart][3*nCellVtx + l] = 0;
          }
          local_weight[ipart][nCellVtx] = 1.;
          nCellVtx += 1;
        }
      }
    }

    PDM_global_mean_field_set (gmid, ipart, 3,
                               local_field[ipart],
                               local_weight[ipart],
                               global_mean_field_ptr[ipart]);

  }

  PDM_global_mean_field_compute (gmid);

  for (int ipart = 0; ipart < nPart; ipart++) {

    int nCell;
    int nFace;
    int nFacePartBound;
    int nVtx;
    int nProc;
    int nTPart;
    int sCellFace;
    int sFaceVtx;
    int sFaceGroup;
    int nFaceGroup2;

    PDM_part_part_dim_get(ppartId,
                       ipart,
                       &nCell,
                       &nFace,
                       &nFacePartBound,
                       &nVtx,
                       &nProc,
                       &nTPart,
                       &sCellFace,
                       &sFaceVtx,
                       &sFaceGroup,
                       &nFaceGroup2);

    int          *cellTag;
    int          *cellFaceIdx;
    int          *cellFace;
    PDM_g_num_t *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;

    PDM_part_part_val_get(ppartId,
                       ipart,
                       &cellTag,
                       &cellFaceIdx,
                       &cellFace,
                       &cellLNToGN,
                       &faceTag,
                       &faceCell,
                       &faceVtxIdx,
                       &faceVtx,
                       &faceLNToGN,
                       &facePartBoundProcIdx,
                       &facePartBoundPartIdx,
                       &facePartBound,
                       &vtxTag,
                       &vtx,
                       &vtxLNToGN,
                       &faceGroupIdx,
                       &faceGroup,
                       &faceGroupLNToGN);

    int nCellVtx = 0;
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        int face = cellFace[j] - 1;
        nCellVtx += faceVtxIdx[face+1] - faceVtxIdx[face];
      }
    }

    for (int i = 0; i < 3 * nCellVtx; i++) {
      if (PDM_ABS (local_field[ipart][i] - global_mean_field_ptr[ipart][i]) > 1e-5) {
        PDM_error (__FILE__, __LINE__, 0, "Error in global mean\n");
        abort();
      }
    }

    free (local_field[ipart]);
    free (local_weight[ipart]);
    free (global_mean_field_ptr[ipart]);
    free (cellVtxGN[ipart]);
  }

  free (local_field);
  free (local_weight);
  free (global_mean_field_ptr);
  free (cellVtxGN);

  PDM_part_free(ppartId);

  PDM_global_mean_free (gmid);

  PDM_dcube_gen_free(id);

  PDM_MPI_Finalize();

  return 0;
}
