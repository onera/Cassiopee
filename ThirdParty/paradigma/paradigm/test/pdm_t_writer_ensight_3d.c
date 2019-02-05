#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"

#include "pdm_writer.h"
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
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
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
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
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

  PDM_g_num_t  nVtxSeg = 10;
  double        length  = 1.;
  int           nPart   = 1;
  int           post    = 0;
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

  PDM_dcube_gen_init(&id,
                      PDM_MPI_COMM_WORLD,
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

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc(dNCell*sizeof(int));
//                  "PDM_PART_RENUM_CELL_CUTHILL",
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int nPropertyCell = 0;
  int nPropertyFace = 0;

  PDM_part_create(&ppartId,
                  PDM_MPI_COMM_WORLD,
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

  free(dCellPart);

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

  int id_cs = PDM_writer_create("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_3d_ens",
                                "chrd3d",
                                PDM_MPI_COMM_WORLD,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  /* Creation de la geometrie */

  int id_geom = PDM_writer_geom_create(id_cs,
                             "test3d_geom",
                             PDM_WRITER_OFF,
                             PDM_WRITER_OFF,
                             nPart);

  /* Creation des variables */

  int id_var_num_part = PDM_writer_var_create(id_cs,
                                    PDM_WRITER_OFF,
                                    PDM_WRITER_VAR_SCALAIRE,
                                    PDM_WRITER_VAR_ELEMENTS,
                                    "num_part");

  int id_var_coo_x = PDM_writer_var_create(id_cs,
                                 PDM_WRITER_ON,
                                 PDM_WRITER_VAR_SCALAIRE,
                                 PDM_WRITER_VAR_SOMMETS,
                                 "coo_x");

  int id_var_coo_xyz = PDM_writer_var_create(id_cs,
                                   PDM_WRITER_ON,
                                   PDM_WRITER_VAR_VECTEUR,
                                   PDM_WRITER_VAR_SOMMETS,
                                   "coo_xyz");

  /* Debut d'ecritures */

  PDM_writer_step_beg(id_cs, 0.);

  int **faceVtxNb = (int **) malloc(sizeof(int *) * nPart);
  int **cellFaceNb = (int **) malloc(sizeof(int *) * nPart);
  
  PDM_real_t **val_num_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * nPart);
  PDM_real_t **val_coo_x    = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * nPart);
  PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * nPart);
  int *nsom_part  = (int *) malloc(sizeof(int) * nPart);

  int *nPartProcs = (int *) malloc(sizeof(int) * numProcs);

  PDM_MPI_Allgather((void *) &nPart,      1, PDM_MPI_INT,
                (void *) nPartProcs, 1, PDM_MPI_INT,
                PDM_MPI_COMM_WORLD);

  int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

  debPartProcs[0] = 0;
  for (int i = 0; i < numProcs; i++) {
    debPartProcs[i+1] = debPartProcs[i] + nPartProcs[i];
  }

  free(nPartProcs);

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
    int nEdgeGroup2;

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
                       &nEdgeGroup2);

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

    assert(sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

    faceVtxNb[ipart] = (int *) malloc(sizeof(int) * nFace);
    cellFaceNb[ipart] = (int *) malloc(sizeof(int) * nCell);
    
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
    
    for (int i = 0; i < nCell; i++) {
      cellFaceNb[ipart][i] = cellFaceIdx[i+1] - cellFaceIdx[i]; 
    }

    for (int i = 0; i < nFace; i++) {
      faceVtxNb[ipart][i] = faceVtxIdx[i+1] - faceVtxIdx[i]; 
    }

    val_num_part[ipart] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * nCell);
    val_coo_x[ipart]    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * nVtx);
    val_coo_xyz[ipart]  = (PDM_real_t *) malloc(sizeof(PDM_real_t) * 3 * nVtx);
    nsom_part[ipart]    = nVtx;

    for (int i = 0; i < nCell; i++) {
      val_num_part[ipart][i] = ipart + 1 + debPartProcs[myRank];
    }

    for (int i = 0; i < nVtx; i++) {
      val_coo_x[ipart][i]       = vtx[3*i];
      val_coo_xyz[ipart][3*i  ] = vtx[3*i  ];
      val_coo_xyz[ipart][3*i+1] = vtx[3*i+1];
      val_coo_xyz[ipart][3*i+2] = vtx[3*i+2];
    }

    PDM_writer_geom_coord_set(id_cs,
                      id_geom,
                      ipart,
                      nVtx,
                      vtx,
                      vtxLNToGN);



    /* Construction de la connectivite pour sortie graphique */

    PDM_writer_geom_cell3d_cellface_add (id_cs,
                                         id_geom,
                                         ipart, 
                                         nCell,
                                         nFace,
                                         faceVtxIdx,
                                         faceVtxNb[ipart],
                                         faceVtx,
                                         cellFaceIdx,
                                         cellFaceNb[ipart],
                                         cellFace,
                                         cellLNToGN);  

  }

  free(debPartProcs);

  PDM_writer_geom_write(id_cs,
                        id_geom);

  /* Creation des variables :
      - numero de partition
      - scalaire
      - vecteur
      - tenseur
   */

  for (int ipart = 0; ipart < nPart; ipart++) {

    PDM_writer_var_set(id_cs,
               id_var_num_part,
               id_geom,
               ipart,
               val_num_part[ipart]);
  }

  PDM_writer_var_write(id_cs,
             id_var_num_part);

  PDM_writer_var_free(id_cs,
             id_var_num_part);

  for (int ipart = 0; ipart < nPart; ipart++) {
    free(val_num_part[ipart]);
  }
  free(val_num_part);

  for (int nstep = 0; nstep < 10; nstep++) {

    double tstep = nstep * 0.01;

    if (nstep > 0)
      PDM_writer_step_beg(id_cs, tstep);

    for (int ipart = 0; ipart < nPart; ipart++) {

      for (int i = 0; i < nsom_part[ipart]; i++) {
        val_coo_x[ipart][i]       = val_coo_x[ipart][i]       + val_coo_x[ipart][i]/length;
        val_coo_xyz[ipart][3*i]   = val_coo_xyz[ipart][3*i]   + val_coo_xyz[ipart][3*i]/length;
        val_coo_xyz[ipart][3*i+1] = val_coo_xyz[ipart][3*i+1] + val_coo_xyz[ipart][3*i+1]/length;
        val_coo_xyz[ipart][3*i+2] = val_coo_xyz[ipart][3*i+2] + val_coo_xyz[ipart][3*i+2]/length;
      }

      PDM_writer_var_set(id_cs,
                 id_var_coo_x,
                 id_geom,
                 ipart,
                 val_coo_x[ipart]);

      PDM_writer_var_set(id_cs,
                 id_var_coo_xyz,
                 id_geom,
                 ipart,
                 val_coo_xyz[ipart]);

    }

    PDM_writer_var_write(id_cs,
               id_var_coo_x);

    PDM_writer_var_write(id_cs,
               id_var_coo_xyz);

    PDM_writer_var_data_free(id_cs,
                    id_var_coo_x);

    PDM_writer_var_data_free(id_cs,
                    id_var_coo_xyz);

    PDM_writer_step_end(id_cs);
  }

  for (int ipart = 0; ipart < nPart; ipart++) {
    free(cellFaceNb[ipart]);
    free(faceVtxNb[ipart]);
    free(val_coo_x[ipart]);
    free(val_coo_xyz[ipart]);
  }
  free(val_coo_x);
  free(val_coo_xyz);
  free(cellFaceNb);
  free(faceVtxNb);
  free(nsom_part);

  PDM_writer_var_free(id_cs,
             id_var_coo_x);

  PDM_writer_var_free(id_cs,
             id_var_coo_xyz);

  /* Liberation memoire */

  PDM_writer_geom_data_free(id_cs,
                   id_geom);

  PDM_writer_geom_free(id_cs,
              id_geom);

  PDM_writer_free(id_cs);


  /* Calculs statistiques */

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

  PDM_part_stat_get(ppartId,
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

  if (myRank == 0) {
    PDM_printf("Statistics :\n");
    PDM_printf("  - Number of cells :\n");
    PDM_printf("       * average            : %i\n", cells_average);
    PDM_printf("       * median             : %i\n", cells_median);
    PDM_printf("       * standard deviation : %12.5e\n", cells_std_deviation);
    PDM_printf("       * min                : %i\n", cells_min);
    PDM_printf("       * max                : %i\n", cells_max);
    PDM_printf("  - Number of faces exchanging with another partition :\n");
    PDM_printf("       * average            : %i\n", bound_part_faces_average);
    PDM_printf("       * median             : %i\n", bound_part_faces_median);
    PDM_printf("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
    PDM_printf("       * min                : %i\n", bound_part_faces_min);
    PDM_printf("       * max                : %i\n", bound_part_faces_max);
    PDM_printf("       * total              : %i\n", bound_part_faces_sum);
  }

  PDM_part_free(ppartId);

  PDM_dcube_gen_free(id);

  PDM_MPI_Finalize();

  return 0;
}

