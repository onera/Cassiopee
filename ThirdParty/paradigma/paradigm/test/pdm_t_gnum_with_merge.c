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
#include "pdm_timer.h"


/*============================================================================
 * Type definition
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
 int          *nPart,
 double       *length
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

    else if (strcmp (argv[i], "-nPart") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        int _nPart = atoi (argv[i]);
        *nPart = (PDM_g_num_t) _nPart;
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

  int id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, pdm_mpi_comm);

  PDM_gnum_set_from_coords (id, 0, dNVtx, dVtxCoord, NULL);

  PDM_gnum_compute (id);

  const PDM_g_num_t *_numabs2 = PDM_gnum_get (id, 0);

  /* const PDM_g_num_t *_numabs2 = NULL; */
  /* int id; */

  /* int nn = 10000; */

  /* for (int i = 0; i < nn; i++) { */

  /*   if (i < nn - 1) { */

  /*     id = PDM_gnum_create (3, 1, PDM_TRUE, 1e-3, pdm_mpi_comm); */
  /*   } */

  /*   else { */
  /*     id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, pdm_mpi_comm); */

  /*   } */

  /*   double *dd = malloc (sizeof(double) * dNVtx); */

  /*   for (int j = 0; j < dNVtx; j++) { */
  /*     dd[j] = 1e-5; */
  /*   } */

  /*   if (i < nn - 1) { */
  /*     PDM_gnum_set_from_coords (id, 0, dNVtx, dVtxCoord, dd); */
  /*   } */
  /*   else { */
  /*     PDM_gnum_set_from_coords (id, 0, dNVtx, dVtxCoord, NULL); */
  /*   } */

  /*   PDM_gnum_compute (id); */

  /*   _numabs2 = PDM_gnum_get (id, 0); */

  /*   free(dd); */

  /*   if (i < nn - 1) { */
  /*     PDM_gnum_free (id, 0); */
  /*   } */

  /*   FILE *f = fopen("/proc/self/statm", "r"); */

  /*   long int mvirt, mres, mshared, val1, val2, val3; */
  /*   fscanf(f, "%ld %ld %ld %ld %ld %ld", &mvirt, &mres, &mshared, &val1, &val2, &val3); */

  /*   long int m_mvirt, m_mres, m_mshared; */

  /*   PDM_MPI_Allreduce (&mvirt, &m_mvirt, 1, PDM_MPI_LONG, PDM_MPI_MAX, pdm_mpi_comm); */
  /*   PDM_MPI_Allreduce (&mres, &m_mres, 1, PDM_MPI_LONG, PDM_MPI_MAX, pdm_mpi_comm); */
  /*   PDM_MPI_Allreduce (&mshared, &m_mshared, 1, PDM_MPI_LONG, PDM_MPI_MAX, pdm_mpi_comm); */

  /*   if (myRank == 0) { */
  /*     printf("mem %d %d : %ld Ko %ld Ko %ld Ko\n", i, dNVtx, 4*m_mvirt , 4*m_mres, 4*m_mshared); */
  /*     //      printf("mem %d %d : %ld Mo %ld Mo %ld Mo\n", i, dNVtx, 4*m_mvirt/1024 , 4*m_mres/1024, 4*m_mshared/1024); */
  /*   } */
  /*   fclose(f); */

  /* } */

//  for (int j = 0; j < dNVtx; j++) {
//    PDM_printf (PDM_FMT_G_NUM" %12.5e %12.5e %12.5e\n", _numabs2[j], dVtxCoord[3*j],
//                                                     dVtxCoord[3*j+1],
//                                                     dVtxCoord[3*j+2]);
//  }

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

  PDM_g_num_t *distrib = (PDM_g_num_t *) malloc((numProcs+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dNVtx = (PDM_g_num_t) dNVtx;

  PDM_MPI_Allgather((void *) &_dNVtx,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) &(distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    pdm_mpi_comm);

  // mesh->face_distrib[0] = 1;
  distrib[0] = 0;

  for (int i = 1; i < numProcs; i++) {
    distrib[i] +=  distrib[i-1];
  }


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



  free (dVtxCoord);
  free (dFaceVtxIdx);
  free (dFaceVtx);
  free (dFaceEdge);
  free (dEdgeVtxIdx);
  free (dEdgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);

  int id2 = PDM_gnum_create (3, nPart, PDM_TRUE, 1e-3, pdm_mpi_comm);

  double **char_length = malloc(sizeof(double *) * nPart);

  int *nVtxs = malloc (sizeof(int) * nPart);
  PDM_g_num_t **vtxLNToGNs = malloc (sizeof(PDM_g_num_t *) * nPart);

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

    nVtxs[ipart] = nVtx;

    char_length[ipart] = malloc (sizeof(double) * nVtx);

    for (int j = 0; j < nVtx; j++) {
      char_length[ipart][j] = HUGE_VAL;
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

    vtxLNToGNs[ipart] = vtxLNToGN;

    for (int j = 0; j < nEdge; j++) {
      int i1 = faceVtxIdx[j];
      int n = faceVtxIdx[j+1] - i1;
      for (int k = 0; k < n; k++) {
        int vtx1 = faceVtx[i1+k] - 1;
        int vtx2 = faceVtx[i1+(k+1)%n] - 1;
        double _lEdge = (vtx[3*vtx2  ] - vtx[3*vtx1  ]) * (vtx[3*vtx2  ] - vtx[3*vtx1  ]) +
                        (vtx[3*vtx2+1] - vtx[3*vtx1+1]) * (vtx[3*vtx2+1] - vtx[3*vtx1+1]) +
                        (vtx[3*vtx2+2] - vtx[3*vtx1+2]) * (vtx[3*vtx2+2] - vtx[3*vtx1+2]);
        char_length[ipart][vtx1] = PDM_MIN (char_length[ipart][vtx1], _lEdge);
        char_length[ipart][vtx2] = PDM_MIN (char_length[ipart][vtx2], _lEdge);
      }
    }

    PDM_gnum_set_from_coords (id2, ipart, nVtx, vtx, char_length[ipart]);

  }

  PDM_timer_t *timer = PDM_timer_create();
  PDM_timer_resume(timer);
  PDM_gnum_compute (id2);
  PDM_timer_hang_on(timer);
  printf("Compute gnum end %12.5es\n", PDM_timer_elapsed(timer));
  fflush(stdout);
  PDM_timer_free(timer);

  const PDM_g_num_t **_numabs = malloc (sizeof(PDM_g_num_t *) * nPart);

  // Check

  PDM_g_num_t *numabs_init = malloc(sizeof(PDM_g_num_t) * dNVtx);

  for (int i = 0; i < dNVtx; i++) {
    numabs_init[i] = distrib[myRank] + 1 + i;
  }

  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP, 1.,
                                                        &numabs_init,
                                                        NULL,
                                                        &dNVtx,
                                                        1,
                                                        pdm_mpi_comm);


  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP, 1.,
                                                        vtxLNToGNs,
                                                        NULL,
                                                        nVtxs,
                                                        nPart,
                                                        pdm_mpi_comm);


  for (int ipart = 0; ipart < nPart; ipart++) {

    _numabs[ipart] = PDM_gnum_get (id2, ipart);

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

  }

  int nElb1 =  PDM_part_to_block_n_elt_block_get (ptb1);

  int nElb2 =  PDM_part_to_block_n_elt_block_get (ptb2);

  assert (nElb1 == nElb2);

  PDM_g_num_t *block_numabs2;
  PDM_g_num_t *block_numabs;
  int *block_stride;

  PDM_part_to_block_exch (ptb1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                          (void **) &_numabs2,
                          &block_stride,
                          (void **) &block_numabs2);

  PDM_part_to_block_exch (ptb2,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                          (void **) _numabs,
                          &block_stride,
                          (void **) &block_numabs);

  for (int i = 0; i < nElb1; i++) {
    if (block_numabs[i] != block_numabs2[i]) {
      PDM_printf("-- diff %d : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" \n",
             i, block_numabs2[i], block_numabs[i]);
      PDM_error (__FILE__, __LINE__, 0, "Error in the generated numbering\n");
    }

  }

//  PDM_g_num_t *n1 = PDM_part_to_block_block_gnum_get (ptb1);
//  PDM_g_num_t *n2 = PDM_part_to_block_block_gnum_get (ptb2);

  for (int ipart = 0; ipart < nPart; ipart++) {
    free (char_length[ipart]);
  }
  free (char_length);

  free(_numabs);
  free(block_numabs);
  free(block_numabs2);
  free(nVtxs);
  free(vtxLNToGNs);
  free(numabs_init);
  free(distrib);

  PDM_gnum_free (id2, 0);

  PDM_gnum_free (id, 0);

  PDM_part_to_block_free (ptb1);
  PDM_part_to_block_free (ptb2);

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
              &nPart,
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
