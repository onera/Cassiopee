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
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"

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

  PDM_g_num_t  nVtxSeg  = 3;
  double        length  = 1.;
  int           nPart   = 1;
  int           post    = 0;
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#else
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
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

  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  const double xmax = xmin + length;
  const double ymax = ymin + length;
  const double zmax = zmin + length;

  if (myRank == 0) {
    printf("-- Build cube\n");
    fflush(stdout);
  }

  PDM_dcube_gen_init(&id,
                     PDM_MPI_COMM_WORLD,
                     nVtxSeg,
                     length,
                     xmin,
                     ymin,
                     zmin);

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

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int nPropertyCell = 0;
  int nPropertyFace = 0;

  if (myRank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  PDM_part_create(&ppartId,
                  PDM_MPI_COMM_WORLD,
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

  free(dCellPart);

  /* free (dFaceCell); */
  /* free (dFaceVtxIdx); */
  /* free (dFaceVtx); */
  /* free (dVtxCoord); */
  /* free (dFaceGroupIdx); */
  /* free (dFaceGroup); */

  int n_point_cloud = 1;
  int id_dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_SURFACE_MESH,
                                      n_point_cloud,
                                      PDM_MPI_COMM_WORLD);

  int **select_face = malloc (sizeof(int *) * nPart);
  int *n_select_face = malloc (sizeof(int) * nPart);
  int **select_vtx = malloc (sizeof(int *) * nPart);
  int *n_select_vtx = malloc (sizeof(int) * nPart);

  int **surface_face_vtx_idx =  malloc (sizeof(int *) * nPart);
  int **surface_face_vtx =  malloc (sizeof(int *) * nPart);
  double **surface_coords = malloc (sizeof(double *) * nPart);

  PDM_g_num_t **surface_face_parent_gnum = malloc (sizeof(PDM_g_num_t *) * nPart);
  PDM_g_num_t **surface_vtx_parent_gnum = malloc (sizeof(PDM_g_num_t *) * nPart);

  const PDM_g_num_t **surface_face_gnum = malloc (sizeof(PDM_g_num_t *) * nPart);
  const PDM_g_num_t **surface_vtx_gnum = malloc (sizeof(PDM_g_num_t *) * nPart);

  int id_gnum_face = PDM_gnum_create (3, nPart, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);
  int id_gnum_vtx = PDM_gnum_create (3, nPart, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double **cell_volume = malloc (sizeof(double *) * nPart);
  double **cell_center = malloc (sizeof(double *) * nPart);

  if (myRank == 0) {
    printf("-- mesh dist set\n");
    fflush(stdout);
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
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppartId,
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

    n_select_face[ipart] = 0;
    n_select_vtx[ipart] = 0;

    select_face[ipart] = malloc (sizeof(int) * nFace);

    for (int i = 0; i < nFace; i++) {
      select_face[ipart][i] = 0;
    }

    select_vtx[ipart] = malloc (sizeof(int) * nVtx);

    for (int i = 0; i < nVtx; i++) {
      select_vtx[ipart][i] = 0;
    }

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

    PDM_part_part_val_get (ppartId,
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

    int iii = 0;
    for (int i = 0; i < nFace; i++) {
      int icel2 = faceCell[2*i+1];
      if (icel2 == 0) {
        iii++;
        select_face[ipart][i] = 1;
      }
    }

    for (int i = 0; i < facePartBoundProcIdx[numProcs]; i++) {
      select_face[ipart][facePartBound[4*i]-1] = 0;
    }

    int idx = 1;
    int s_face_vtx = 0;
    for (int i = 0; i < nFace; i++) {
      if (select_face[ipart][i] == 1) {
        select_face[ipart][i] = idx;
        s_face_vtx += (faceVtxIdx[i+1] - faceVtxIdx[i]);
        idx += 1;
      }
    }
    n_select_face[ipart] = idx - 1;

    for (int i = 0; i < nFace; i++) {
      if (select_face[ipart][i] != 0) {
        for (int j = faceVtxIdx[i]; j < faceVtxIdx[i+1]; j++) {
          select_vtx[ipart][faceVtx[j]-1] = 1;
        }
      }
    }

    idx = 1;
    for (int i = 0; i < nVtx; i++) {
      if (select_vtx[ipart][i] == 1) {
        select_vtx[ipart][i] = idx;
        idx += 1;
      }
    }
    n_select_vtx[ipart] = idx - 1;

    surface_face_vtx_idx[ipart] = malloc (sizeof(int) * (n_select_face[ipart] + 1));
    surface_face_vtx_idx[ipart][0] = 0;
    surface_face_vtx[ipart] = malloc (sizeof(int) * s_face_vtx);

    surface_coords[ipart] = malloc (sizeof(double) * 3 * n_select_vtx[ipart]);

    surface_face_parent_gnum[ipart] =
      malloc (sizeof(PDM_g_num_t) * n_select_face[ipart]);
    surface_vtx_parent_gnum[ipart] =
      malloc (sizeof(PDM_g_num_t) * n_select_vtx[ipart]);

    surface_face_gnum[ipart] = NULL;
    surface_vtx_gnum[ipart] = NULL;

    idx = 0;
    int idx1 = 0;
    for (int i = 0; i < nFace; i++) {
      if (select_face[ipart][i] > 0) {
        surface_face_vtx_idx[ipart][idx+1] =
          surface_face_vtx_idx[ipart][idx] + (faceVtxIdx[i+1] - faceVtxIdx[i]);

        surface_face_parent_gnum[ipart][idx] = faceLNToGN[i];

        idx += 1;

        for (int j = faceVtxIdx[i]; j < faceVtxIdx[i+1]; j++) {
          surface_face_vtx[ipart][idx1++] = select_vtx[ipart][faceVtx[j]-1];
        }
      }
    }

    idx = 0;
    for (int i = 0; i < nVtx; i++) {
      if (select_vtx[ipart][i] > 0) {

        surface_vtx_parent_gnum[ipart][idx] = vtxLNToGN[i];
        surface_coords[ipart][3*idx  ] = vtx[3*i];
        surface_coords[ipart][3*idx+1] = vtx[3*i+1];
        surface_coords[ipart][3*idx+2] = vtx[3*i+2];

        idx += 1;

      }
    }

    PDM_gnum_set_from_parents (id_gnum_face,
                               ipart,
                               n_select_face[ipart],
                               surface_face_parent_gnum[ipart]);


    PDM_gnum_set_from_parents (id_gnum_vtx,
                               ipart,
                               n_select_vtx[ipart],
                               surface_vtx_parent_gnum[ipart]);

    for (int i = 0; i <  n_select_vtx[ipart]; i++) {

    }

  }

  PDM_gnum_compute (id_gnum_face);

  PDM_gnum_compute (id_gnum_vtx);

  PDM_g_num_t n_g_face_loc = 0;
  PDM_g_num_t n_g_vtx_loc = 0;

  PDM_g_num_t n_g_face = 0;
  PDM_g_num_t n_g_vtx = 0;

  for (int ipart = 0; ipart < nPart; ipart++) {
    surface_face_gnum[ipart] = PDM_gnum_get (id_gnum_face, ipart);
    surface_vtx_gnum[ipart] = PDM_gnum_get (id_gnum_vtx, ipart);

    for (int i = 0; i <  n_select_face[ipart]; i++) {
      n_g_face_loc = PDM_MAX(n_g_face_loc, surface_face_gnum[ipart][i]);
    }

    for (int i = 0; i <  n_select_vtx[ipart]; i++) {
      n_g_vtx_loc = PDM_MAX(n_g_vtx_loc, surface_vtx_gnum[ipart][i]);
    }

  }

  PDM_MPI_Allreduce (&n_g_face_loc, &n_g_face, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_MPI_Allreduce (&n_g_vtx_loc, &n_g_vtx, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (id_dist,
                                           n_g_face,
                                           n_g_vtx,
                                           nPart);

  PDM_dist_cloud_surf_n_part_cloud_set (id_dist, 0, nPart);

  for (int ipart = 0; ipart < nPart; ipart++) {

    PDM_dist_cloud_surf_surf_mesh_part_set (id_dist,
                                      ipart,
                                      n_select_face[ipart],
                                      surface_face_vtx_idx[ipart],
                                      surface_face_vtx[ipart],
                                      surface_face_gnum[ipart],
                                      n_select_vtx[ipart],
                                      surface_coords[ipart],
                                      surface_vtx_gnum[ipart]);

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

    PDM_part_part_dim_get (ppartId,
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

    PDM_part_part_val_get (ppartId,
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

    const int     isOriented = 0;
    cell_volume[ipart] = malloc(sizeof(double) * nCell);
    cell_center[ipart] = malloc(sizeof(double) * 3 * nCell);

    PDM_geom_elem_polyhedra_properties (isOriented,
                                        nCell,
                                        nFace,
                                        faceVtxIdx,
                                        faceVtx,
                                        cellFaceIdx,
                                        cellFace,
                                        nVtx,
                                        vtx,
                                        cell_volume[ipart],
                                        cell_center[ipart],
                                        NULL,
                                        NULL);

    PDM_dist_cloud_surf_cloud_set (id_dist,
                             0,
                             ipart,
                             nCell,
                             cell_center[ipart],
                             cellLNToGN);

  }

  if (myRank == 0) {
    printf("-- Dist compute\n");
    fflush(stdout);
  }

  PDM_dist_cloud_surf_compute (id_dist);

  if (myRank == 0) {
    printf("-- Dist check\n");
    fflush(stdout);
  }

  for (int ipart = 0; ipart < nPart; ipart++) {
    double      *distance;
    double      *projected;
    PDM_g_num_t *closest_elt_gnum;

    PDM_dist_cloud_surf_get (id_dist,
                       0,
                       ipart,
                       &distance,
                       &projected,
                       &closest_elt_gnum);

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

    PDM_part_part_dim_get (ppartId,
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

    PDM_part_part_val_get (ppartId,
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

    int ierr = 0;
    for (int i = 0; i < nCell; i++) {
      double d1 = PDM_MIN (PDM_ABS (cell_center[ipart][3*i] - xmin),
                           PDM_ABS (cell_center[ipart][3*i] - xmax));
      double d2 = PDM_MIN (PDM_ABS (cell_center[ipart][3*i+1] - ymin),
                           PDM_ABS (cell_center[ipart][3*i+1] - ymax));
      double d3 = PDM_MIN (PDM_ABS (cell_center[ipart][3*i+2] - zmin),
                           PDM_ABS (cell_center[ipart][3*i+2] - zmax));
      double d = PDM_MIN (PDM_MIN (d1,d2), d3);
      d = d * d;
      if (PDM_ABS(distance[i] - d) > 1e-6) {
        ierr += 1;
        printf ("Erreur distance %d (%12.5e %12.5e %12.5e) : %12.5e %12.5e "PDM_FMT_G_NUM"\n", i,
                cell_center[ipart][3*i],
                cell_center[ipart][3*i+1],
                cell_center[ipart][3*i+2],
                distance[i],
                d,
                closest_elt_gnum[i]);
      }
      /* else { */
      /*   if ((i+1 == 874) || */
      /*       (i+1 == 1266) || */
      /*       (i+1 == 1069) || */
      /*       (i+1 == 1071) || */
      /*       (i+1 == 1056) || */
      /*       (i+1 == 1084) || */
      /*       (i+1 == 1070)){ */

      /*     printf ("Affiche distance %d (%12.5e %12.5e %12.5e) : %12.5e %12.5e %ld\n", i+1, */
      /*             cell_center[ipart][3*i], */
      /*             cell_center[ipart][3*i+1], */
      /*             cell_center[ipart][3*i+2], */
      /*             distance[i], */
      /*             d, */
      /*             closest_elt_gnum[i]); */
      /*   } */
      /* } */
    }

    if (ierr > 0) {
      printf ("Erreur distance pour %d points\n", ierr);
      abort();
    }

    if (myRank == 0) {
      printf ("elements surfaciques : "PDM_FMT_G_NUM"\n", 6*(nVtxSeg-1)*(nVtxSeg-1));
      printf ("nombre de points     : "PDM_FMT_G_NUM"\n", nVtxSeg*nVtxSeg*nVtxSeg);
      fflush(stdout);
    }
  }

  PDM_part_free(ppartId);

  PDM_dcube_gen_free(id);
  PDM_dist_cloud_surf_dump_times(id_dist);
  int partial = 0;
  PDM_dist_cloud_surf_free (id_dist, partial);

  for (int ipart = 0; ipart < nPart; ipart++) {
    free (select_face[ipart]);
    free (select_vtx[ipart]);

    free (surface_face_vtx_idx[ipart]);
    free (surface_face_vtx[ipart]);
    free (surface_coords[ipart]);

    free (surface_face_parent_gnum[ipart]);
    free (surface_vtx_parent_gnum[ipart]);

  }

  free (select_face);
  free (select_vtx);

  free (n_select_face);
  free (n_select_vtx);

  free (surface_face_vtx_idx);
  free (surface_face_vtx);
  free (surface_coords);

  free (surface_face_parent_gnum);
  free (surface_vtx_parent_gnum);

  free (surface_face_gnum);
  free (surface_vtx_gnum);

  PDM_gnum_free(id_gnum_face, 0);
  PDM_gnum_free(id_gnum_vtx, 0);

  PDM_MPI_Finalize();

   if (myRank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
