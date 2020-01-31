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
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part.h"

#include "pdm_writer.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_overlay.h"
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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
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
 double        *length,
 int           *n_part,
 int           *post,
 int           *method,
 int           *haveRandom
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
    else if (strcmp (argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_part = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-no_random") == 0) {
      *haveRandom = 0;
    }
    else if (strcmp (argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp (argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp (argv[i], "-parmetis") == 0) {
      *method = 1;
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
 PDM_g_num_t   *nGEdge
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
  int           nEdgeGroup;
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
                     &nEdgeGroup,
                     &dEdgeGroupIdx,
                     &dEdgeGroup);

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
    for (int i = 0; i < nEdgeGroup; i++) {
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
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i+1]);
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
                   nEdgeGroup,
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
  return ppartId;
}


/**
 *
 * \brief  Compute faceVtx connectivity
 *
 * \param [in]      ppartId  ppart identifier
 * \param [in]      nPart    Number of partitions
 *
 * \return          faceVtx connectivity for each partition of each mesh
 */

static int ***
_compute_faceVtx
(
 int *ppartId,
 int nPart
)
{

  int ***faceVtx = (int ***) malloc(sizeof(int **) * 2);

  for (int imesh = 0; imesh < 2; imesh++) {

    int id_ppart = ppartId[imesh];

    faceVtx[imesh] = (int **) malloc (sizeof(int *) * nPart);

    for (int ipart = 0; ipart < nPart; ipart++) {

      int nFace;
      int nEdge;
      int nEdgePartBound;
      int nVtx;
      int nProc;
      int nTPart;
      int sFaceEdge;
      int sEdgeVtx;
      int sEdgeGroup;
      int nEdgeGroup2;

      PDM_part_part_dim_get (id_ppart,
                          ipart,
                          &nFace,
                          &nEdge,
                          &nEdgePartBound,
                          &nVtx,
                          &nProc,
                          &nTPart,
                          &sFaceEdge,
                          &sEdgeVtx,
                          &sEdgeGroup,
                          &nEdgeGroup2);

      int          *faceTag;
      int          *faceEdgeIdx;
      int          *faceEdge;
      PDM_g_num_t *faceLNToGN;
      int          *edgeTag;
      int          *edgeFace;
      int          *edgeVtxIdx;
      int          *edgeVtx;
      PDM_g_num_t *edgeLNToGN;
      int          *edgePartBoundProcIdx;
      int          *edgePartBoundPartIdx;
      int          *edgePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t *vtxLNToGN;
      int          *edgeGroupIdx;
      int          *edgeGroup;
      PDM_g_num_t *edgeGroupLNToGN;

      assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

      PDM_part_part_val_get (id_ppart,
                          ipart,
                          &faceTag,
                          &faceEdgeIdx,
                          &faceEdge,
                          &faceLNToGN,
                          &edgeTag,
                          &edgeFace,
                          &edgeVtxIdx,
                          &edgeVtx,
                          &edgeLNToGN,
                          &edgePartBoundProcIdx,
                          &edgePartBoundPartIdx,
                          &edgePartBound,
                          &vtxTag,
                          &vtx,
                          &vtxLNToGN,
                          &edgeGroupIdx,
                          &edgeGroup,
                          &edgeGroupLNToGN);

      faceVtx[imesh][ipart] = (int *) malloc(sizeof(int) * sFaceEdge);
      int *_faceVtx = faceVtx[imesh][ipart];

      int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (nVtx + 1));

      for (int i = 0; i < nVtx + 1; i++) {
        vtxEdgeIdx[i] = 0;
      }

      for (int i = 0; i < nEdge; i++) {
        int ivtx1 = edgeVtx[2*i];
        int ivtx2 = edgeVtx[2*i + 1];

        vtxEdgeIdx[ivtx1] += 1;
        vtxEdgeIdx[ivtx2] += 1;
      }

      for (int i = 1; i < nVtx + 1; i++) {
        vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i-1];
      }

      int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[nVtx]);
      int *vtxEdgeN = (int *) malloc(sizeof(int) * nVtx);
      for (int i = 0; i < nVtx; i++) {
        vtxEdgeN[i] = 0;
      }

      for (int i = 0; i < nEdge; i++) {
        int ivtx1 = edgeVtx[2*i] - 1;
        int ivtx2 = edgeVtx[2*i + 1] - 1;
        int iedge = i + 1;

        vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
        vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
        vtxEdgeN[ivtx1] += 1;
        vtxEdgeN[ivtx2] += 1;
      }

      free(vtxEdgeN);

      for (int i = 0; i < nFace; i++) {
        int idx = faceEdgeIdx[i];
        int _nEdge = faceEdgeIdx[i+1] - idx;
        int *_edges = faceEdge + idx;
        int *_vertices = _faceVtx + idx;

        int edge_cur = _edges[0];
        int vtx_deb =  edgeVtx[2*(edge_cur - 1)];
        _vertices[0] = vtx_deb;
        int vtx_cur =  edgeVtx[2*(edge_cur - 1) + 1];
        int idxVtx = 0;

        while (vtx_deb != vtx_cur) {
          _vertices[++idxVtx] = vtx_cur;
          int find_vtx = 0;

          for (int j = vtxEdgeIdx[vtx_cur - 1]; j <  vtxEdgeIdx[vtx_cur]; j++) {
            for (int k = 0; k < _nEdge; k++) {
              if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
                edge_cur = _edges[k];
                if (edgeVtx[2*(_edges[k]-1)] == vtx_cur) {
                  vtx_cur = edgeVtx[2*(_edges[k]-1) + 1];
                }
                else {
                  vtx_cur = edgeVtx[2*(_edges[k]-1)];
                }
                find_vtx = 1;
                break;
              }
            }
            if (find_vtx)
              break;
          }
          if (!find_vtx) {
            PDM_error(__FILE__, __LINE__, 0,"Error to compute vtxedge !!!!\n");
            abort();
          }
        }
      }

      free (vtxEdge);
      free (vtxEdgeIdx);

    }
  }
  return faceVtx;
}


/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_export_ini_mesh
(
 const PDM_MPI_Comm pdm_mpi_comm,
 int*           ppartId,
 const int      nPart
)
{

  int myRank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Export Mesh to Ensight
   */

  int id_cs[2];

  id_cs[0] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "mesh1",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  id_cs[1] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "mesh2",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1,
                                NULL);

  /*
   * Creation des variables
   */

  int id_var_num_part[2];
  int id_var_coo_x[2];
  int id_var_coo_xyz[2];
  int id_geom[2];

  for (int imesh = 0; imesh < 2; imesh++) {

    id_var_num_part[imesh] = PDM_writer_var_create (id_cs[imesh],
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAIRE,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "num_part");

    id_var_coo_x[imesh] = PDM_writer_var_create (id_cs[imesh],
                                       PDM_WRITER_ON,
                                       PDM_WRITER_VAR_SCALAIRE,
                                       PDM_WRITER_VAR_SOMMETS,
                                       "coo_x");

    id_var_coo_xyz[imesh] = PDM_writer_var_create (id_cs[imesh],
                                         PDM_WRITER_ON,
                                         PDM_WRITER_VAR_VECTEUR,
                                         PDM_WRITER_VAR_SOMMETS,
                                         "coo_xyz");

    /*
     * Creation de la geometrie
     */

    char nom_geom[6];
    if (imesh == 0)
      strcpy (nom_geom,"mesh1");
    else
      strcpy (nom_geom,"mesh2");

    id_geom[imesh] = PDM_writer_geom_create (id_cs[imesh],
                                   nom_geom,
                                   PDM_WRITER_OFF,
                                   PDM_WRITER_OFF,
                                   nPart);
    /*
     * Debut des ecritures
     */

    int       **edgeVtxIdx1  = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * nPart);
    int       **edgeVtxNB1   = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * nPart);
    int       **faceEdgeIdx1 = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * nPart);
    int       **faceEdgeNB1  = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * nPart);

    int *nsom_part  = (int *) malloc(sizeof(int) * nPart);

    int *nPartProcs = (int *) malloc(sizeof(int) * numProcs);

    PDM_MPI_Allgather ((void *) &nPart,      1, PDM_MPI_INT,
                   (void *) nPartProcs, 1, PDM_MPI_INT,
                   PDM_MPI_COMM_WORLD);

    int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + nPartProcs[i];
    }

    free(nPartProcs);

    PDM_writer_step_beg (id_cs[imesh], 0.);

    for (int ipart = 0; ipart < nPart; ipart++) {

      int nFace;
      int nEdge;
      int nEdgePartBound;
      int nVtx;
      int nProc;
      int nTPart;
      int sFaceEdge;
      int sEdgeVtx;
      int sEdgeGroup;
      int nEdgeGroup2;

      PDM_part_part_dim_get (ppartId[imesh],
                          ipart,
                          &nFace,
                          &nEdge,
                          &nEdgePartBound,
                          &nVtx,
                          &nProc,
                          &nTPart,
                          &sFaceEdge,
                          &sEdgeVtx,
                          &sEdgeGroup,
                          &nEdgeGroup2);

      int          *faceTag;
      int          *faceEdgeIdx;
      int          *faceEdge;
      PDM_g_num_t *faceLNToGN;
      int          *edgeTag;
      int          *edgeFace;
      int          *edgeVtxIdx;
      int          *edgeVtx;
      PDM_g_num_t *edgeLNToGN;
      int          *edgePartBoundProcIdx;
      int          *edgePartBoundPartIdx;
      int          *edgePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t *vtxLNToGN;
      int          *edgeGroupIdx;
      int          *edgeGroup;
      PDM_g_num_t *edgeGroupLNToGN;

      assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

      PDM_part_part_val_get (ppartId[imesh],
                          ipart,
                          &faceTag,
                          &faceEdgeIdx,
                          &faceEdge,
                          &faceLNToGN,
                          &edgeTag,
                          &edgeFace,
                          &edgeVtxIdx,
                          &edgeVtx,
                          &edgeLNToGN,
                          &edgePartBoundProcIdx,
                          &edgePartBoundPartIdx,
                          &edgePartBound,
                          &vtxTag,
                          &vtx,
                          &vtxLNToGN,
                          &edgeGroupIdx,
                          &edgeGroup,
                          &edgeGroupLNToGN);

      edgeVtxIdx1[ipart]  = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nEdge);
      edgeVtxNB1[ipart]   = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nEdge);
      faceEdgeIdx1[ipart] = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nFace);
      faceEdgeNB1[ipart]  = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nFace);

      for (int i = 0; i < nFace; i++) {
        faceEdgeNB1[ipart][i] = faceEdgeIdx[i+1] - faceEdgeIdx[i];
        faceEdgeIdx1[ipart][i] = faceEdgeIdx[i] + 1;
      }

      for (int i = 0; i < nEdge; i++) {
        edgeVtxNB1[ipart][i] = edgeVtxIdx[i+1] - edgeVtxIdx[i];
        edgeVtxIdx1[ipart][i] = edgeVtxIdx[i] + 1;
      }

      PDM_writer_geom_coord_set (id_cs[imesh],
                         id_geom[imesh],
                         ipart,
                         nVtx,
                         vtx,
                         vtxLNToGN);

      PDM_writer_geom_cell2d_cellface_add (id_cs[imesh],
                                   id_geom[imesh],
                                   ipart,
                                   nFace,
                                   nEdge,
                                   edgeVtxIdx1[ipart],
                                   edgeVtxNB1[ipart],
                                   edgeVtx,
                                   faceEdgeIdx1[ipart],
                                   faceEdgeNB1[ipart],
                                   faceEdge,
                                   faceLNToGN);
    }

    PDM_writer_geom_write(id_cs[imesh],
                id_geom[imesh]);

    for (int ipart = 0; ipart < nPart; ipart++) {
      free (edgeVtxIdx1[ipart]);
      free (edgeVtxNB1[ipart]);
      free (faceEdgeIdx1[ipart]);
      free (faceEdgeNB1[ipart]);
    }

    free (edgeVtxIdx1);
    free (edgeVtxNB1);
    free (faceEdgeIdx1);
    free (faceEdgeNB1);

    /* Creation des variables :
       - numero de partition
       - scalaire
       - vecteur
       - tenseur
    */

    PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * nPart);
    PDM_real_t **val_coo_x    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * nPart);
    PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * nPart);

    for (int ipart = 0; ipart < nPart; ipart++) {

      int nFace;
      int nEdge;
      int nEdgePartBound;
      int nVtx;
      int nProc;
      int nTPart;
      int sFaceEdge;
      int sEdgeVtx;
      int sEdgeGroup;
      int nEdgeGroup2;

      PDM_part_part_dim_get (ppartId[imesh],
                          ipart,
                          &nFace,
                          &nEdge,
                          &nEdgePartBound,
                          &nVtx,
                          &nProc,
                          &nTPart,
                          &sFaceEdge,
                          &sEdgeVtx,
                          &sEdgeGroup,
                          &nEdgeGroup2);

      int          *faceTag;
      int          *faceEdgeIdx;
      int          *faceEdge;
      PDM_g_num_t *faceLNToGN;
      int          *edgeTag;
      int          *edgeFace;
      int          *edgeVtxIdx;
      int          *edgeVtx;
      PDM_g_num_t *edgeLNToGN;
      int          *edgePartBoundProcIdx;
      int          *edgePartBoundPartIdx;
      int          *edgePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t *vtxLNToGN;
      int          *edgeGroupIdx;
      int          *edgeGroup;
      PDM_g_num_t *edgeGroupLNToGN;

      assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

      PDM_part_part_val_get (ppartId[imesh],
                          ipart,
                          &faceTag,
                          &faceEdgeIdx,
                          &faceEdge,
                          &faceLNToGN,
                          &edgeTag,
                          &edgeFace,
                          &edgeVtxIdx,
                          &edgeVtx,
                          &edgeLNToGN,
                          &edgePartBoundProcIdx,
                          &edgePartBoundPartIdx,
                          &edgePartBound,
                          &vtxTag,
                          &vtx,
                          &vtxLNToGN,
                          &edgeGroupIdx,
                          &edgeGroup,
                          &edgeGroupLNToGN);

      val_num_part[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nFace);
      val_coo_x[ipart]    = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nVtx);
      val_coo_xyz[ipart]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * 3 * nVtx);
      nsom_part[ipart]    = nVtx;

      for (int i = 0; i < nFace; i++) {
        val_num_part[ipart][i] = ipart + 1 + debPartProcs[myRank];
      }

      for (int i = 0; i < nVtx; i++) {
        val_coo_x[ipart][i]       = vtx[3*i];
        val_coo_xyz[ipart][3*i  ] = vtx[3*i  ];
        val_coo_xyz[ipart][3*i+1] = vtx[3*i+1];
        val_coo_xyz[ipart][3*i+2] = vtx[3*i+2];
      }

      PDM_writer_var_set (id_cs[imesh],
                  id_var_num_part[imesh],
                  id_geom[imesh],
                  ipart,
                  val_num_part[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                  id_var_coo_x[imesh],
                  id_geom[imesh],
                  ipart,
                  val_coo_x[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                  id_var_coo_xyz[imesh],
                  id_geom[imesh],
                  ipart,
                  val_coo_xyz[ipart]);

    }

    PDM_writer_var_write (id_cs[imesh],
                id_var_num_part[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                id_var_num_part[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                id_var_coo_x[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                id_var_coo_x[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                id_var_coo_xyz[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                id_var_coo_xyz[imesh]);

    for (int ipart = 0; ipart < nPart; ipart++) {
      free (val_num_part[ipart]);
      free (val_coo_x[ipart]);
      free (val_coo_xyz[ipart]);
    }

    free (val_num_part);
    free (val_coo_x);
    free (val_coo_xyz);
    free (nsom_part);

    PDM_writer_step_end (id_cs[imesh]);
    PDM_writer_geom_data_free (id_cs[imesh],
                      id_geom[imesh]);

    PDM_writer_geom_free (id_cs[imesh],
                 id_geom[imesh]);
    PDM_writer_free (id_cs[imesh]);

    free (debPartProcs);

  }
}

/**
 *
 * \brief  Export overlay mesh
 *
 * \param [in]    pdm_id    PDM identifier
 *
 */

static void
_export_ol_mesh
(
const int pdm_id
)
{
  pdm_id;
  PDM_error(__FILE__, __LINE__, 0, "Error _export_ol_mesh : Not implemented yet\n");
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
  int           post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
  int           haveRandom = 1;
  int           myRank;
  int           numProcs;
  int           ppartId[2];

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &nVtxSeg,
              &length,
              &nPart,
              &post,
              (int *) &method,
              &haveRandom);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);


  /*
   *  Create a partitioned mesh
   */

  PDM_g_num_t nGFace[2];
  PDM_g_num_t nGVtx[2];
  PDM_g_num_t nGEdge[2];

  for (int imesh = 0; imesh < 2; imesh++) {

    ppartId[imesh] = _create_split_mesh (imesh,
                                         PDM_MPI_COMM_WORLD,
                                         nVtxSeg,
                                         length,
                                         nPart,
                                         method,
                                         haveRandom,
                                         nGFace + imesh,
                                         nGVtx + imesh,
                                         nGEdge + imesh);

  }

  _export_ini_mesh (PDM_MPI_COMM_WORLD,
                    ppartId,
                    nPart);

  /*
   *  Appel des fonctions d'intersection
   */

  double projectCoeff = 0.;

  /*
   *  Creation de l'objet PDM
   */

  int pdm_id = PDM_ol_create (nPart,
                              nGFace[0],
                              nGVtx[0],
                              nPart,
                              nGFace[1],
                              nGVtx[1],
                              projectCoeff,
                              PDM_MPI_COMM_WORLD);

  PDM_ol_parameter_set (pdm_id,
                        PDM_OL_CAR_LENGTH_TOL,
                        1e-4);

  PDM_ol_parameter_set (pdm_id,
                        PDM_OL_EXTENTS_TOL,
                        1e-4);
  /*
   *  faceVtx connectivity computation
   */

  int ***faceVtx = _compute_faceVtx (ppartId, nPart);

  /*
   *  Initial meshes definition
   */

  for (int imesh = 0; imesh < 2; imesh++) {

    int id_ppart = ppartId[imesh];
    PDM_ol_mesh_t mesh;
    if (imesh == 0) {
      mesh = PDM_OL_MESH_A;
    }
    else {
      mesh = PDM_OL_MESH_B;
    }

    for (int ipart = 0; ipart < nPart; ipart++) {

      int nFace;
      int nEdge;
      int nEdgePartBound;
      int nVtx;
      int nProc;
      int nTPart;
      int sFaceEdge;
      int sEdgeVtx;
      int sEdgeGroup;
      int nEdgeGroup2;

      PDM_part_part_dim_get (id_ppart,
                          ipart,
                          &nFace,
                          &nEdge,
                          &nEdgePartBound,
                          &nVtx,
                          &nProc,
                          &nTPart,
                          &sFaceEdge,
                          &sEdgeVtx,
                          &sEdgeGroup,
                          &nEdgeGroup2);

      int          *faceTag;
      int          *faceEdgeIdx;
      int          *faceEdge;
      PDM_g_num_t *faceLNToGN;
      int          *edgeTag;
      int          *edgeFace;
      int          *edgeVtxIdx;
      int          *edgeVtx;
      PDM_g_num_t *edgeLNToGN;
      int          *edgePartBoundProcIdx;
      int          *edgePartBoundPartIdx;
      int          *edgePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t *vtxLNToGN;
      int          *edgeGroupIdx;
      int          *edgeGroup;
      PDM_g_num_t *edgeGroupLNToGN;

      PDM_printf ("nEdge 1 %d\n", nEdge);

      assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

      PDM_part_part_val_get (id_ppart,
                          ipart,
                          &faceTag,
                          &faceEdgeIdx,
                          &faceEdge,
                          &faceLNToGN,
                          &edgeTag,
                          &edgeFace,
                          &edgeVtxIdx,
                          &edgeVtx,
                          &edgeLNToGN,
                          &edgePartBoundProcIdx,
                          &edgePartBoundPartIdx,
                          &edgePartBound,
                          &vtxTag,
                          &vtx,
                          &vtxLNToGN,
                          &edgeGroupIdx,
                          &edgeGroup,
                          &edgeGroupLNToGN);

      int *_faceVtx = faceVtx[imesh][ipart];
      const int *faceVtxIdx = faceEdgeIdx;

      PDM_ol_input_mesh_set (pdm_id,
                             mesh,
                             ipart,
                             nFace,
                             faceVtxIdx,
                             _faceVtx,
                             faceLNToGN,
                             nVtx,
                             vtx,
                             vtxLNToGN);

    }

  }

  /*
   *  Calcul
   */

  PDM_ol_compute (pdm_id);

  /*
   *  Export overlay mesh
   */

  _export_ol_mesh (pdm_id);

  /*
   *  Free Pdm
   */

  PDM_ol_del (pdm_id);

  /*
   *  Free meshes
   */

  for (int imesh = 0; imesh < 2; imesh++) {
    for (int ipart = 0; ipart < nPart; ipart++) {
      free (faceVtx[imesh][ipart]);
    }
    free (faceVtx[imesh]);

    PDM_part_free (ppartId[imesh]);
  }

  free (faceVtx);
  PDM_MPI_Finalize ();

  PDM_printf ("\nfin test\n");

  return 0;

}
