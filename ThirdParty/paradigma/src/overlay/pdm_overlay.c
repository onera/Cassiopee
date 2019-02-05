
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_overlay.h"
#include "pdm_overlay_priv.h"
#include "pdm_dbbtree.h"
#include "pdm_sort.h"
#include "pdm_box_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_hash_tab.h"
#include "pdm_dhash_tab.h"
#include "pdm_poly_clipp.h"
#include "pdm_surf_mesh_priv.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _DOT_PRODUCT(v1, v2) \
  (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2) ( \
 cross_v1_v2[0] = v1[1]*v2[2] - v1[2]*v2[1],   \
 cross_v1_v2[1] = v1[2]*v2[0] - v1[0]*v2[2],   \
 cross_v1_v2[2] = v1[0]*v2[1] - v1[1]*v2[0]  )

#define _MODULE(v) \
  sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_data_t
 * \brief Type of data
 *
 */

typedef enum {

  XMIN    = 0,  /*!< XMIN */
  YMIN    = 1,  /*!< YMIN */
  ZMIN    = 2,  /*!< ZMIN */
  XMAX    = 3,  /*!< XMAX */
  YMAX    = 4,  /*!< YMAX */
  ZMAX    = 5,  /*!< ZMAX */

} _extents_t;

/**
 * \struct _sub_edge_t
 * \brief  sub edge definition to store in a hash table structure
 *
 */

typedef struct {

  PDM_g_num_t vtx1;       /*!< First Vertex */
  double     coords1[3];/*!< Vertex1 Coordinates */

  PDM_g_num_t vtx2;      /*!< second Vertex */
  double     coords2[3];/*!< Vertex 2 Coordinates */

} _sub_edge_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 * \brief Compute edge entities
 *
 * This function defines edges of initial meshes and computes their conectivities
 *
 * \param [in]  ol      overlay objet
 *
 */

static void
_build_edges
(
 PDM_ol_t *ol
)
{
  
  /*
   * Build edges on each part
   */ 

  PDM_surf_mesh_build_edges (ol->meshA);

  PDM_surf_mesh_build_edges (ol->meshB);

  /*
   * Define global numbering
   */ 

  PDM_surf_mesh_build_edges_gn_and_edge_part_bound (ol->meshA);

  PDM_surf_mesh_build_edges_gn_and_edge_part_bound (ol->meshB);

}


/**
 * \brief Build communication graph between internal partitions of each initial mesh 
 *
 * This function builds the communication graph between internal partitions 
 * of each initial mesh 
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_build_exchange_graph
(
 PDM_ol_t *ol
)
{

  PDM_surf_mesh_build_exchange_graph (ol->meshA);

  PDM_surf_mesh_build_exchange_graph (ol->meshB);

}


/**
 * \brief Build ghost faces and ghost edges
 *
 * This function builds ghost faces and ghost edges
 * of each initial mesh 
 *
 * \param [in]  ol       overlay object                     
 *
 */

static void
_compute_carLgthVtx
(
 PDM_ol_t *ol
)
{
  PDM_surf_mesh_compute_carLgthVtx (ol->meshA);
  
  PDM_surf_mesh_compute_carLgthVtx (ol->meshB);
}


/**
 * \brief Compute face extents
 *
 * This function computes face extents
 *
 * \param [in]  ol       overlay object                     
 *
 */

static void
_compute_faceExtents
(
 PDM_ol_t *ol
)
{
  PDM_surf_mesh_compute_faceExtentsMesh(ol->meshA,
                                        ol->extentsTol);

  PDM_surf_mesh_compute_faceExtentsMesh(ol->meshB,
                                        ol->extentsTol);
}


/**
 * \brief Chek if the two meshes are plane surfaces
 *
 * This function cheks the two meshes are plane surfaces
 * and returns the distance between the two meshes
 *
 * \param [in]   ol     overlay object                     
 * \param [out]  dist   Distance between the two planes
 *
 */

static int
_is_same_plane
(
 PDM_ol_t *ol,
 double   *dist
)
{
  double planeEquationA[4] = {0, 0, 0, 0};
  double planeEquationB[4] = {0, 0, 0, 0};
  double barycenterA[3] = {0, 0, 0};
  double barycenterB[3] = {0, 0, 0};
  
  int isPlaneMeshA = PDM_surf_mesh_is_plane_surface (ol->meshA,
                                                     ol->samePlaneTol,
                                                     planeEquationA,
                                                     barycenterA);
  
  int isPlaneMeshB = PDM_surf_mesh_is_plane_surface (ol->meshB,
                                                     ol->samePlaneTol,
                                                     planeEquationB,
                                                     barycenterB);
  
  int isSamePlane = isPlaneMeshA && isPlaneMeshB;

  if (isSamePlane) {
    if (   (1 - fabs (_DOT_PRODUCT (planeEquationA, planeEquationB))
         > ol->samePlaneTol)) {
      isSamePlane = 0;
    }
  }

  if (isSamePlane) {
    double dist1 = fabs (_DOT_PRODUCT (planeEquationA, barycenterB));
    double dist2 = fabs (_DOT_PRODUCT (planeEquationB, barycenterA));
    *dist = _MAX (dist1, dist2);
  }
    
  return isSamePlane;
}

/**
 * \brief Compute overlay planes
 *
 * This function computes overlay mesh of two plane meshes
 *
 * \param [in]  ol       overlay object                     
 *
 */

static void
_compute_overlay_planes
(
 PDM_ol_t *ol
)
{
  
  //TODO: Cette fonction est très grosse. Elle sera decoupee en petites fonctions
  //      lors de l'optimisation de l'algorithme.
  
  int myRank;
  PDM_MPI_Comm_rank (ol->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (ol->comm, &lComm);

  /*****************************************************************************
   *                                                                           *
   * Create a DBBtree structure (dbbtreeA) to store mesh A elements boxes      *
   * (boxes B)                                                                 *
   *                                                                           *
   ****************************************************************************/ 

  const int dim = 3;
  ol->dbbtreeA = PDM_dbbtree_create (ol->comm, dim);

  PDM_surf_mesh_t *meshA = ol->meshA;
  const int nPartA = PDM_surf_mesh_n_part_get (meshA);

  int *nEltsA = (int *) malloc (sizeof(int) * nPartA);
  const PDM_g_num_t **gNumA = 
                  (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * nPartA);
  const double **extentsA = (const double **) malloc (sizeof(double *) * nPartA);

  for (int i = 0; i < nPartA; i++) {
    nEltsA[i] = PDM_surf_mesh_part_n_face_get (meshA, i); 
    extentsA[i] = PDM_surf_mesh_part_extents_get (meshA, i); 
    gNumA[i] = PDM_surf_mesh_part_face_g_num_get (meshA, i);
  }

  PDM_box_set_t  *boxesA = PDM_dbbtree_boxes_set (ol->dbbtreeA,
                                                  nPartA,
                                                  nEltsA,
                                                  extentsA,
                                                  gNumA);
  
  free (nEltsA);
  free (gNumA);
  free (extentsA);

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_BUILD_DBBTREE] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_BUILD_DBBTREE]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_BUILD_DBBTREE]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_BUILD_DBBTREE]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*****************************************************************************
   *                                                                           *
   * Look for intersections between each mesh A boxes and mesh B boxes         *
   * the result is stored inside boxesB_intersection_index (size = n_eltA + 1);*
   *                             boxesB_intersection_index_l_num               *
   *                             (local number of B boxes)                     *
   *                                                                           *
   ****************************************************************************/ 

  PDM_surf_mesh_t *meshB = ol->meshB;
  const int nPartB = PDM_surf_mesh_n_part_get (meshB);

  int *nEltsB = (int *) malloc (sizeof(int) * nPartB);
  const PDM_g_num_t **gNumB = (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * nPartB);
  const double **extentsB = (const double **) malloc (sizeof(double *) * nPartB);

  for (int i = 0; i < nPartB; i++) {
    nEltsB[i] = PDM_surf_mesh_part_n_face_get (meshB, i);
    extentsB[i] = PDM_surf_mesh_part_extents_get (meshB, i);
    gNumB[i] = PDM_surf_mesh_part_face_g_num_get (meshB, i);
  }

  int *boxesB_intersection_index;
  int *boxesB_intersection_l_num;

  PDM_box_set_t  *boxesB = PDM_dbbtree_intersect_boxes_set (ol->dbbtreeA,
                                                            nPartB,
                                                            nEltsB,
                                                            extentsB,
                                                            gNumB,
                                                            &boxesB_intersection_index,
                                                            &boxesB_intersection_l_num);
  
  free (nEltsB);
  free (gNumB);
  free (extentsB);

  /*
   * Check boxesA and boxesB entities 
   */ 
  
  const int *originA = PDM_box_set_origin_get (boxesA);

  const int *originB = PDM_box_set_origin_get (boxesB);

  const PDM_g_num_t *gnum_eltA = (const PDM_g_num_t *) PDM_box_set_get_g_num (boxesA);
  int              n_eltA    = PDM_box_set_get_size (boxesA);

  const PDM_g_num_t *gnum_eltB = (const PDM_g_num_t *) PDM_box_set_get_g_num (boxesB);
  int              n_eltB    = PDM_box_set_get_size (boxesB);

  if (1 == 0) {
    
    PDM_printf ("gnum_eltA :");
    for (int i = 0; i < n_eltA; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, gnum_eltA[i]);
    }
    PDM_printf ("\n");
    
    PDM_printf ("gnum_eltB :");
    for (int i = 0; i < n_eltB; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, gnum_eltB[i]);
    }
    PDM_printf ("\n");
    
  }

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_BOXES_INTERSECT] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_BOXES_INTERSECT]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_BOXES_INTERSECT]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_BOXES_INTERSECT]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*****************************************************************************
   *                                                                           *
   * Redistribute boxesA intersections to ensure a good load balacing          *
   * in ol->comm MPI communicator                                              *
   * This step removes intersections found many times on different ranks       *
   *                                                                           *
   * After this step, data are stored in a block with n_elt_blockA, block_gnumA, 
   * part_strideA                                                              *
   *                                                                           *
   ****************************************************************************/ 

  PDM_part_to_block_t *ptb_boxesA = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_MERGE,
                                                            1.,
                                                            (PDM_g_num_t **) &gnum_eltA,
                                                            &n_eltA,
                                                            1,
                                                            ol->comm);

  int n_elt_blockA = PDM_part_to_block_n_elt_block_get (ptb_boxesA);

  PDM_g_num_t *block_gnumA = PDM_part_to_block_block_gnum_get (ptb_boxesA);

  int *part_strideA = (int *) malloc (sizeof(int) * n_eltA);

  for (int i = 0; i < n_eltA; i++) {
    part_strideA[i] = 
            boxesB_intersection_index[i+1] - boxesB_intersection_index[i];
  }

  /*
   * Redistribute data boxesA from blockB distribution with a PDM_box_distrib_t
   * structure
   * TODO: - Build a new PDM_box_distrib_t structure more simple
   *       - Hide PDM_box_distrib_t attributes
   */

  PDM_l_num_t *destination = PDM_part_to_block_destination_get (ptb_boxesA);

  PDM_g_num_t n_g_eltA = PDM_box_set_get_global_size(boxesA);

  PDM_box_distrib_t *distribA = 
          PDM_box_distrib_create (n_eltA,
                                  n_g_eltA,
                                  1, // Don't use in this case
                                  ol->comm);

  PDM_g_num_t n_g_eltB = PDM_box_set_get_global_size(boxesB);

  PDM_box_distrib_t *distribB = PDM_box_distrib_create (n_eltB,
                                                        n_g_eltB,
                                                        1, // Don't use in this case
                                                        ol->comm);

  int *countEltsA = (int *) malloc (sizeof(int) * n_eltA);
  int *countEltsB = (int *) malloc (sizeof(int) * n_eltB);
  
  for (int i = 0; i < lComm + 1; i++) {
    distribA->index[i] = 0;
    countEltsA[i] = 0;
    distribB->index[i] = 0;
    countEltsB[i] = 0;
  }
  
  for (int i = 0; i < n_eltA; i++) {
    int iProc = destination[i] + 1;
    distribA->index[iProc]++;
    distribB->index[iProc] += part_strideA[i];
  }

  for (int i = 0; i < lComm; i++) {
    distribA->index[i+1] += distribA->index[i];
    distribB->index[i+1] += distribB->index[i];
  }
  
  distribA->list = (int *) malloc (sizeof(int) * distribA->index[lComm]);
  distribB->list = (int *) malloc (sizeof(int) * distribB->index[lComm]);
  
  for (int i = 0; i < n_eltA; i++) {
    int iProc = destination[i] + 1;
    int idxA = distribA->index[iProc] + (countEltsA[iProc]++);
    distribA->list[idxA] = i;
    int idxB = distribB->index[iProc] + countEltsB[iProc];
    countEltsB[iProc] += part_strideA[i];
    for (int j = boxesB_intersection_index[i]; j < boxesB_intersection_index[i+1]; j++) {
      distribB->list[idxB] = boxesB_intersection_l_num[j];
    }
  }
  
  free (countEltsA);
  free (countEltsB);

  PDM_box_set_redistribute (distribA, boxesA);
  PDM_box_set_redistribute (distribB, boxesB);
  
  PDM_box_distrib_destroy (&distribA);
  PDM_box_distrib_destroy (&distribB);

  /*****************************************************************************
   *                                                                           *
   *  Transfer intersection information from partitions to blocks              *
   * with PDM_part_to_block_exch function                                       *
   *                                                                           *
   *  Results :                                                                *
   *      - blockA_boxesB_idx                                                  *
   *      - blockA_boxesB_gnum_data                                            *
   *                                                                           *
   ****************************************************************************/

  if (1 == 0) {
    PDM_printf ("box_l_num :\n");
    for (int i = 0; i < n_eltA; i++) {
      PDM_printf ("box_l_numA "PDM_FMT_G_NUM" : ", gnum_eltA[i]);
      for (int k = boxesB_intersection_index[i]; k < boxesB_intersection_index[i+1]; k++) {
        PDM_printf (" %d", boxesB_intersection_l_num[k]);
      }
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }
  
  PDM_g_num_t *boxesB_intersection_g_num = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *
                                                               boxesB_intersection_index[n_eltA]);

  for (int k = 0; k < boxesB_intersection_index[n_eltA]; k++) {
    boxesB_intersection_g_num[k] =  gnum_eltB[boxesB_intersection_l_num[k]];
  }

  if (1 == 0) {
    PDM_printf ("box_g_num :");
    for (int k = 0; k < boxesB_intersection_index[n_eltA]; k++) {
      PDM_printf (" "PDM_FMT_G_NUM, boxesB_intersection_g_num[k]);
    }
    PDM_printf ("\n");
  }

  int       *blockA_boxesB_stride;
  PDM_g_num_t *blockA_boxesB_gnum_data;

  PDM_part_to_block_exch (ptb_boxesA, 
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         0,
                         &part_strideA,
                         (void **) &boxesB_intersection_g_num,
                         &blockA_boxesB_stride,
                         (void **) &blockA_boxesB_gnum_data);
  
  if (1 == 0) {
    PDM_printf ("box_g_num :\n");
    for (int i = 0; i < n_eltA; i++) {
      PDM_printf ("[%d] box_g_numA "PDM_FMT_G_NUM" : ", myRank, gnum_eltA[i]);
      for (int k = boxesB_intersection_index[i]; k < boxesB_intersection_index[i+1]; k++) {
        PDM_printf (" "PDM_FMT_G_NUM, boxesB_intersection_g_num[k]);
      }
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  int *blockA_boxesB_idx = (int *) malloc (sizeof(int) * (n_elt_blockA + 1));
  for (int i = 0; i < n_elt_blockA + 1; i++) {
    blockA_boxesB_idx[i] = 0;
  }
  for (int i = 0; i < n_elt_blockA; i++) {
    blockA_boxesB_idx[i+1] = blockA_boxesB_idx[i] + blockA_boxesB_stride[i];
  }

  int idx = 0;
  int *blockA_boxesB_idx_new = (int *) malloc (sizeof(int) * (n_elt_blockA + 1));
  for (int i = 0; i < n_elt_blockA + 1; i++) {
    blockA_boxesB_idx_new[i] = 0;
  }

  for (int i = 0; i < n_elt_blockA; i++) {

    PDM_g_num_t *ideb = blockA_boxesB_gnum_data + blockA_boxesB_idx[i];
    int length = blockA_boxesB_idx[i+1] - blockA_boxesB_idx[i];

    PDM_sort_long (ideb, NULL, length);

    PDM_g_num_t pre = -1;
    for (int k = 0; k < length; k++) {
      if (pre != ideb[k]) {
        blockA_boxesB_gnum_data[idx++] = ideb[k];
        blockA_boxesB_idx_new[i+1] += 1;
        pre = ideb[k];
      }
    }
  }

  for (int i = 0; i < n_elt_blockA; i++) {
    blockA_boxesB_idx_new[i+1] += blockA_boxesB_idx_new[i];
  }

  for (int i = 0; i < n_elt_blockA + 1; i++) {
    blockA_boxesB_idx[i] = blockA_boxesB_idx_new[i];
  }

  if (1 == 1) {
    for (int i = 0; i < n_elt_blockA; i++) {
      PDM_printf ("[%d] gnum_block "PDM_FMT_G_NUM" : ", myRank, block_gnumA[i]);
      for (int k = blockA_boxesB_idx[i]; k < blockA_boxesB_idx[i+1]; k++) {
        PDM_printf (" "PDM_FMT_G_NUM, blockA_boxesB_gnum_data[k]);
      }
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  free (part_strideA);
  free (blockA_boxesB_idx_new);
  free (boxesB_intersection_g_num);
  free (boxesB_intersection_l_num);
  free (boxesB_intersection_index);

  PDM_part_to_block_free (ptb_boxesA);
  
  /***************************************************************************** 
   *                                                                           *
   *  Transfer from origin data to compute polygon intersections               *
   * (function PDM_box_set_recv_data_from_origin_distrib)                      *
   *      - faceA->edgesA connectivities (absolute number) (mapping)           *
   *      - faceB->edgesB connectivities (absolute number) (mapping)           *
   *      - faceA->vtxA connectivities (absolute number) (mapping)             *
   *      - faceB->vtxB connectivities (absolute number) (mapping)             *   
   *      - faceA->coordsVtxA (To be allocate)                                 *
   *      - faceB->coordsVtxB (To be allocate)                                 *
   *      - faceA->epsVtxA (To be allocate)                                    *
   *      - faceB->epsVtxB (To be allocate)                                    *
   *                                                                           *
   ****************************************************************************/  
  
  int        *faceStrideCurrent[2];
  int        *faceStrideCurrent3[2];
  PDM_g_num_t *faceToEdgeCurrent[2]; 
  PDM_g_num_t *faceToVtxCurrent[2];
  double     *faceVtxCooCurrent[2];
  double     *faceVtxEpsCurrent[2];
  
  PDM_surf_mesh_t *meshes[2] = {meshA, meshB};
  
  for (int imesh = 0; imesh < 2; imesh++) {
  
    size_t data_size_l = sizeof(PDM_g_num_t);
    size_t data_size_r = sizeof(double);
    PDM_surf_mesh_t *mesh = meshes[imesh];

    int nPart = PDM_surf_mesh_n_part_get (mesh);
    
    int         **faceStrideOrigin = (int **) malloc (sizeof(int *) * nPart);
    PDM_g_num_t  **faceToEdgeOrigin = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * nPart);
    PDM_g_num_t  **faceToVtxOrigin  = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * nPart);
    double      **faceVtxCooOrigin = (double **) malloc (sizeof(double *) * nPart);
    double      **faceVtxEpsOrigin = (double **) malloc (sizeof(double *) * nPart);

    for (int i = 0; i < nPart; i++) {
      int               nEltPart        = PDM_surf_mesh_part_n_face_get (mesh, i); 
      const int        *partFaceEdgeIdx = PDM_surf_mesh_part_face_edge_idx_get (mesh, i);
      const int        *partFaceEdge    = PDM_surf_mesh_part_face_edge_get (mesh, i);
      const int        *partFaceVtx     = PDM_surf_mesh_part_face_vtx_get (mesh, i);
      const PDM_g_num_t *partEdgeGNum    = PDM_surf_mesh_part_edge_g_num_get (mesh, i);
      const PDM_g_num_t *partVtxGNum     = PDM_surf_mesh_part_vtx_g_num_get (mesh, i);
      const double     *partVtxCoord    = PDM_surf_mesh_part_vtx_get (mesh, i);
      const double     *partVtxEps      = PDM_surf_mesh_part_carLgthVtx_get (mesh, i);

      faceStrideOrigin[i] = (int *) malloc (sizeof(int) * nEltPart);
      faceToEdgeOrigin[i] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * partFaceEdgeIdx[nEltPart]);
      faceToVtxOrigin[i]  = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *partFaceEdgeIdx[nEltPart]);
      faceVtxCooOrigin[i] = (double *) malloc (sizeof(double) * 3 * partFaceEdgeIdx[nEltPart]);
      faceVtxEpsOrigin[i] = (double *) malloc (sizeof(double) * 3 * partFaceEdgeIdx[nEltPart]);
      
      int        *_faceStrideOrigin = faceStrideOrigin[i];
      PDM_g_num_t *_faceToEdgeOrigin = faceToEdgeOrigin[i];
      PDM_g_num_t *_faceToVtxOrigin  = faceToVtxOrigin[i];
      double     *_faceVtxCooOrigin = faceVtxCooOrigin[i];
      double     *_faceVtxEpsOrigin = faceVtxEpsOrigin[i];
    
      for (int j = 0; j < nEltPart; j++) {
        _faceStrideOrigin[j] = partFaceEdgeIdx[j+1] - partFaceEdgeIdx[j];
      }
    
      for (int j = 0; j < partFaceEdgeIdx[nEltPart]; j++) {
        int edge1 = partFaceEdge[j] - 1;
        int vtx1 = partFaceVtx[j] - 1;
        _faceToEdgeOrigin[j] = partEdgeGNum[edge1];
        _faceToVtxOrigin[j] = partVtxGNum[vtx1];
        _faceVtxEpsOrigin[j] = partVtxEps[vtx1];
        for (int k = 0; k < 3; k++) {
          _faceVtxCooOrigin[3*j+k] = partVtxCoord[3*vtx1+k];
        }
      }
    }
  
    faceStrideCurrent[imesh]  = NULL;
    faceStrideCurrent3[imesh] = NULL;
    faceToEdgeCurrent[imesh]  = NULL; 
    faceToVtxCurrent[imesh]   = NULL;
    faceVtxCooCurrent[imesh]  = NULL;
    faceVtxEpsCurrent[imesh]  = NULL;
   
    PDM_box_set_recv_data_from_origin_distrib (boxesA,
                                               PDM_STRIDE_VAR,
                                               1,
                                               data_size_l,
                                               faceStrideOrigin,
                                               (void **) faceToEdgeOrigin,
                                               &(faceStrideCurrent[imesh]),
                                               (void **) &faceToEdgeCurrent[imesh]);
    
    PDM_box_set_recv_data_from_origin_distrib (boxesA,
                                               PDM_STRIDE_VAR,
                                               1,
                                               data_size_l,
                                               faceStrideOrigin,
                                               (void **) faceToEdgeOrigin,
                                               &(faceStrideCurrent[imesh]),
                                               (void **) &faceToVtxCurrent[imesh]);

    PDM_box_set_recv_data_from_origin_distrib (boxesA,
                                               PDM_STRIDE_VAR,
                                               1,
                                               data_size_r,
                                               faceStrideOrigin,
                                               (void **) faceVtxEpsOrigin,
                                               &(faceStrideCurrent[imesh]),
                                               (void **) &faceVtxEpsCurrent[imesh]);
  
    for (int i = 0; i < nPart; i++) {
      int nEltPart        = PDM_surf_mesh_part_n_face_get (mesh, i); 
      for (int j = 0; j < nEltPart; j++) {
        faceStrideOrigin[i][j] = 3 * faceStrideOrigin[i][j];
      }
    }
    
    PDM_box_set_recv_data_from_origin_distrib (boxesA,
                                               PDM_STRIDE_VAR,
                                               1,
                                               data_size_r,
                                               faceStrideOrigin,
                                               (void **) faceVtxCooOrigin,
                                               &(faceStrideCurrent3[imesh]),
                                               (void **) &faceVtxEpsCurrent[imesh]);
   
    for (int i = 0; i < nPart; i++) {
      free (faceStrideOrigin[i]);
      free (faceToEdgeOrigin[i]);
      free (faceToVtxOrigin[i]);
      free (faceVtxCooOrigin[i]);
      free (faceVtxEpsOrigin[i]);
    }

    free (faceStrideOrigin);
    free (faceToEdgeOrigin);
    free (faceToVtxOrigin);
    free (faceVtxCooOrigin);
    free (faceVtxEpsOrigin);
    
  }

  /***************************************************************************** 
   *                                                                           *
   *   Find BoxesB local number to build : blockA_boxesB_lnum_data             *
   *                                                                           *
   *  Steps :                                                                  *
   *      - Sort Boxes B gnum                                                  *
   *      - Binary search in sorted boxes B gnum array                         *
   *                                                                           *
   ****************************************************************************/ 

  int *blockA_boxesB_lnum_data =
    (int *) malloc (sizeof(int) * blockA_boxesB_idx[n_elt_blockA]);
  
  gnum_eltB = (const PDM_g_num_t *) PDM_box_set_get_g_num (boxesB);
  n_eltB    = PDM_box_set_get_size (boxesB);

  int *lnum = (int *) malloc (sizeof(int) * n_eltB);
  PDM_g_num_t *gnum_eltB_cp = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_eltB);

  for (int i = 0; i < n_eltB; i++) {
    lnum[i] = i + 1;
    gnum_eltB_cp[i] = gnum_eltB[i];
  }

  PDM_sort_long (gnum_eltB_cp, lnum, n_eltB);

  for (int i = 0; i <  blockA_boxesB_idx[n_elt_blockA]; i++) {
    PDM_g_num_t box_gnum = blockA_boxesB_gnum_data[i];
    int idx1 = PDM_binary_search_long (box_gnum, gnum_eltB_cp, n_eltB);

    blockA_boxesB_lnum_data[i] = lnum[idx1];
  }
  
  free (lnum);
  free (gnum_eltB_cp);

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_DISTRIB_BOXESA_BLOCK] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_DISTRIB_BOXESA_BLOCK]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_DISTRIB_BOXESA_BLOCK]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_DISTRIB_BOXESA_BLOCK]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /***************************************************************************** 
   *                                                                           *
   *   Perform polygon clippings                                               *
   *      - Use :                                                              *
   *              * faceStrideCurrent                                          *
   *              * faceStrideCurrent3                                         *
   *              * faceToEdgeCurrent                                          *
   *              * faceToVtxCurrent                                           *
   *              * faceVtxCooCurrent                                          *
   *              * faceVtxEpsCurrent                                          *
   *                                                                           *
   ****************************************************************************/            
    
  /* Stocker les intersection dans une structure */
  
  int *faceIdxCurrent[2];
  faceIdxCurrent[0] = (int *) malloc (sizeof(int) * (n_elt_blockA + 1));
  faceIdxCurrent[1] = (int *) malloc (sizeof(int) * (n_eltB + 1));

  for (int i = 0; i < 2; i++) {
    faceIdxCurrent[i][0] = 0;
  }
  for (int j = 0; j < n_elt_blockA; j++) {
    faceIdxCurrent[0][j+1] += faceIdxCurrent[0][j] + faceStrideCurrent[0][j]; 
  }
  for (int j = 0; j < n_eltB; j++) {
    faceIdxCurrent[1][j+1] += faceIdxCurrent[1][j] + faceStrideCurrent[1][j]; 
  }

  /*****************************************************************************
   *                                                                           *
   *                  Pre-clipping : Perform intersections                     *
   *                                                                           *
   ****************************************************************************/
                                                                   
  PDM_g_num_t maxGNEdgeA = PDM_surf_mesh_n_g_edge_get (meshA);
  PDM_g_num_t maxGNEdgeB = PDM_surf_mesh_n_g_edge_get (meshB);

  PDM_edges_intersect_t *intersect = PDM_edges_intersect_create (maxGNEdgeA,
                                                                 maxGNEdgeB,
                                                                 ol->vtxCarLengthTol,
                                                                 lComm);

  for (int i = 0; i < n_elt_blockA; i++) {
    int nVtxA = faceIdxCurrent[0][i + 1] - faceIdxCurrent[0][i];
    PDM_g_num_t *_faceToEdgeA = faceToEdgeCurrent[0] + faceIdxCurrent[0][i];
    PDM_g_num_t *_faceToVtxA  = faceToVtxCurrent[0] + faceIdxCurrent[0][i];
    double     *_faceVtxCooA = faceVtxCooCurrent[0] + 3 * faceIdxCurrent[0][i];
    double     *_faceVtxEpsA = faceVtxEpsCurrent[0] + faceIdxCurrent[0][i];
    
    /*  
     * Perform each edge-edge intersection :
     *      - Each elementA edge intersects each elementB edge
     *      - Storage in a hash table the result of each intersection 
     *        (key = sum of global number of vertices)
     *      - Build sub face element
     */
    
    for (int j = blockA_boxesB_idx[i]; j < blockA_boxesB_idx[i+1]; j++) {
      int lnum_boxB = blockA_boxesB_lnum_data[j];
      int nVtxB = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToEdgeB = faceToEdgeCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToVtxB  = faceToVtxCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      double     *_faceVtxCooB = faceVtxCooCurrent[1] + 3 * faceIdxCurrent[1][lnum_boxB];
      double     *_faceVtxEpsB = faceVtxEpsCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      
      /*
       * Build polygon structure according to clipping polygon structure 
       */
      
      PDM_edges_intersect_poly_add (intersect,
                                    nVtxA,
                                    _faceToEdgeA,
                                    _faceToVtxA,
                                    _faceVtxCooA,
                                    _faceVtxEpsA,
                                    nVtxB,
                                    _faceToEdgeB,
                                    _faceToVtxB,
                                    _faceVtxCooB,
                                    _faceVtxEpsB);
      
      
    }

  }

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_EDGES_INTERSECTION] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_EDGES_INTERSECTION]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_EDGES_INTERSECTION]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_EDGES_INTERSECTION]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  
  /*****************************************************************************
   *                                                                           *
   *  Synchronize intersections                                                *
   *                                                                           *
   *      Synchronize function Defines global number for :                     *
   *      - new A vertices in A overlay mesh                                   *
   *      - new B vertices in B overlay mesh                                   *
   *                                                                           *
   ****************************************************************************/            

  PDM_g_num_t n_g_vtxA = PDM_surf_mesh_n_g_vtx_get (meshA);
  PDM_g_num_t n_g_vtxB = PDM_surf_mesh_n_g_vtx_get (meshB);
  
  PDM_g_num_t n_g_newVtxA = 0;
  PDM_g_num_t n_g_newVtxB = 0;

  PDM_edges_intersect_synchronize (intersect, n_g_vtxA, n_g_vtxB, 
                                   &n_g_newVtxA, &n_g_newVtxB);  

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_EDGES_SYNCHRO] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_EDGES_SYNCHRO]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_EDGES_SYNCHRO]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_EDGES_SYNCHRO]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*****************************************************************************
   *                                                                           *
   *                             Clipping                                      *
   *                                                                           *
   ****************************************************************************/

  int iclipp = 0;

  int s_subFacesToFaces = 4 * n_elt_blockA;
  PDM_g_num_t *subFacesToFaces = malloc (sizeof(PDM_g_num_t) * s_subFacesToFaces);
    
  int s_subFacesConnecIdx = (1 + n_elt_blockA);
  int *subFacesConnecIdx = malloc (sizeof(int) * s_subFacesConnecIdx);
  subFacesConnecIdx[0] = 0;
    
  int s_subFacesConnecA = 4 * n_elt_blockA;
  int s_subFacesConnecB = 4 * n_elt_blockA;
  PDM_g_num_t *subFacesConnecA = malloc (sizeof(PDM_g_num_t) * s_subFacesConnecA);
  int s_gNumSubFacesA = n_elt_blockA;
  PDM_g_num_t *gNumSubFacesA   = malloc (sizeof(PDM_g_num_t) * s_gNumSubFacesA);
  PDM_g_num_t *subFacesConnecB = malloc (sizeof(PDM_g_num_t) * s_subFacesConnecB);

  int s_subFacesCoordsA = 3 * s_subFacesConnecA;
  double *subFacesCoordsA = malloc (sizeof(double) * s_subFacesCoordsA); 
  
  int *facesToSubFacesAIdx = malloc (sizeof(int) * (1 + n_elt_blockA));
  int *facesToSubFacesBIdx = malloc (sizeof(int) * (1 + blockA_boxesB_idx[n_elt_blockA]));
  
  facesToSubFacesAIdx[0] = 0;
  facesToSubFacesBIdx[0] = 0;
  
  for (int i = 0; i < n_elt_blockA; i++) {

    PDM_g_num_t gnum_boxA = block_gnumA[i];
    int nVtxA = faceIdxCurrent[0][i + 1] - faceIdxCurrent[0][i];
    PDM_g_num_t *_faceToEdgeA = faceToEdgeCurrent[0] + faceIdxCurrent[0][i];
    PDM_g_num_t *_faceToVtxA  = faceToVtxCurrent[0] + faceIdxCurrent[0][i];
    double     *_faceVtxCooA = faceVtxCooCurrent[0] + 3 * faceIdxCurrent[0][i];

    facesToSubFacesAIdx[i+1] = facesToSubFacesAIdx[i];
    
    for (int j = blockA_boxesB_idx[i]; j < blockA_boxesB_idx[i+1]; j++) {

      PDM_g_num_t gnum_boxB = blockA_boxesB_gnum_data[j];
      int lnum_boxB = blockA_boxesB_lnum_data[j];
      int nVtxB = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToEdgeB = faceToEdgeCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToVtxB  = faceToVtxCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      double     *_faceVtxCooB = faceVtxCooCurrent[1] + 3 * faceIdxCurrent[1][lnum_boxB];

      int         nPolyClippA      = 0;
      int        *polyClippIdxA    = NULL;
      PDM_g_num_t *polyClippConnecA = NULL;
      double     *polyClippCoordsA = NULL;

      int         nPolyClippB      = 0;
      int        *polyClippIdxB    = NULL;
      PDM_g_num_t *polyClippConnecB = NULL;
      double     *polyClippCoordsB = NULL;
            
      PDM_poly_clipp (intersect,
                      nVtxA,
                      _faceToEdgeA,
                      _faceToVtxA,
                      _faceVtxCooA,
                      nVtxB,
                      _faceToEdgeB,
                      _faceToVtxB,
                      _faceVtxCooB,
                      PDM_POLY_CLIPP_CLIP,                     
                      &nPolyClippA,        
                      &polyClippIdxA,
                      &polyClippConnecA,
                      &polyClippCoordsA,
                      &nPolyClippB,        
                      &polyClippIdxB,
                      &polyClippConnecB,
                      &polyClippCoordsB);
 
      assert (nPolyClippA == nPolyClippB);
      
      facesToSubFacesAIdx[i+1] += nPolyClippA;

      facesToSubFacesBIdx[j+1] = nPolyClippB;

      int newSize = iclipp + nPolyClippA + 1;
      if (newSize > s_subFacesConnecIdx) {
        while (newSize > s_subFacesConnecIdx) {
          s_subFacesConnecIdx *= 2;
        }
        subFacesConnecIdx = realloc (subFacesConnecIdx, 
                                  sizeof(int) * s_subFacesConnecIdx);
      }      

      newSize = iclipp + nPolyClippA;
      if (newSize > s_gNumSubFacesA) {
        while (newSize > s_gNumSubFacesA) {
          s_gNumSubFacesA *= 2;
        }
        gNumSubFacesA  = realloc (gNumSubFacesA, 
                                  sizeof(PDM_g_num_t) * s_gNumSubFacesA);
      }
      
      int ibeg = subFacesConnecIdx[iclipp];
      int ibeg2 = 3 *ibeg;
      newSize = ibeg + polyClippIdxA[nPolyClippA];      
      if (newSize > s_subFacesConnecA) {
        while (newSize > s_subFacesConnecA) {
          s_subFacesConnecA *= 2;
        }
        subFacesConnecA = realloc (subFacesConnecA, 
                                  sizeof(PDM_g_num_t) * s_subFacesConnecA);
      }
 
      if (newSize > s_subFacesConnecB) {
        while (newSize > s_subFacesConnecB) {
          s_subFacesConnecB *= 2;
        }
        subFacesConnecB = realloc (subFacesConnecB, 
                                  sizeof(PDM_g_num_t) * s_subFacesConnecB);        
      }
     
      newSize = 3 * (ibeg + polyClippIdxA[nPolyClippA]);      
      if (newSize > s_subFacesCoordsA) {
        while (newSize > s_subFacesCoordsA) {
          s_subFacesCoordsA *= 2;
        }
        subFacesCoordsA = realloc (subFacesCoordsA, 
                                  sizeof(double) * s_subFacesCoordsA);
      }

      newSize = 4*(iclipp + nPolyClippA);
      if (newSize > s_subFacesToFaces) {
        while (newSize > s_subFacesToFaces) {
          s_subFacesToFaces *= 2;
        }
        subFacesToFaces = realloc (subFacesToFaces, 
                                  sizeof(PDM_g_num_t) * s_subFacesToFaces);
      }      
      
      for (int k = 0; k < nPolyClippA; k++) {
        subFacesConnecIdx[iclipp+1] =  subFacesConnecIdx[iclipp] 
                                     + (polyClippIdxA[k+1] - polyClippIdxA[k]);

        for (int k1 = polyClippIdxA[k]; k1 < polyClippIdxA[k+1]; k1++) {
          subFacesConnecA[ibeg] = polyClippConnecA[k1]; 
          subFacesConnecB[ibeg] = polyClippConnecB[k1];
          ibeg += 1;
        }

        for (int k1 = 3*polyClippIdxA[k]; k1 < 3*polyClippIdxA[k+1]; k1++) {
          subFacesCoordsA[ibeg2++] = polyClippCoordsA[k1]; 
        }
        
        subFacesToFaces[4*iclipp ]  = i;         // A local number
        subFacesToFaces[4*iclipp+1] = gnum_boxA; // A global number
        subFacesToFaces[4*iclipp+2] = j;         // B local number
        subFacesToFaces[4*iclipp+3] = gnum_boxB; // B global number
        
        iclipp += 1;
      }     

      if (polyClippIdxB != NULL) {
        free (polyClippIdxB);
      }

      if (polyClippConnecB != NULL) {
        free (polyClippConnecB);
      }

      if (polyClippIdxA != NULL) {
        free (polyClippIdxA);
      }

      if (polyClippConnecA != NULL) {
        free (polyClippConnecA);
      }
    }
  }
     
  const int nSharedSubFaces = iclipp;

  /*
   * Global numbering Definition of sub-faces 
   * This numbering is shared between A and B
   * Partial split faces are stored behind sub-faces
   * Not split faces are stored behind partial split faces
   * 
   */
   
  PDM_g_num_t nSharedSubFaces_l = nSharedSubFaces;
  PDM_g_num_t beg_nSharedSubFaces;

  PDM_MPI_Scan(&nSharedSubFaces_l, &beg_nSharedSubFaces, 
           1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);
  
  beg_nSharedSubFaces += -nSharedSubFaces;
  
  for (int i = 0; i < nSharedSubFaces; i++) {
    gNumSubFacesA[i] = beg_nSharedSubFaces + i + 1; 
  }
  
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_CLIPPING] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_CLIPPING]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_CLIPPING]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_CLIPPING]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  //
  // TODO: Optimisation : Commencer les communications du part_to_block
  //                      non-bloquantes ici
  //
  
  /*****************************************************************************
   *                                                                           *
   * Compute sub-faces coming from A faces with partial B cover.               *
   *                                                                           *
   * Check that each sub-edge is referenced twice                              *
   *    - Take into sub-edges of sub-faces                                     *
   *    - Take into sub-edges of origin elements                               *
   *                                                                           *
   * If sub-edges are not referenced twice, the origin face is partial covered *
   * by B. Additional sub-faces have to compute.                               *
   *                                                                           *
   ****************************************************************************/
    
  int *facesToAddSubFacesAIdx = malloc (sizeof(int) * (1 + n_elt_blockA));
  facesToAddSubFacesAIdx[0] = facesToSubFacesAIdx[n_elt_blockA];  

  int *faceIniVtxIdxA = malloc (sizeof(int) * (1 + n_elt_blockA));
  faceIniVtxIdxA[0] = 0;

  int s_faceIniVtxA = 4 * n_elt_blockA;
  PDM_g_num_t *faceIniVtxA = malloc (sizeof(PDM_g_num_t) * s_faceIniVtxA);  

  idx = 0;  
  for (int i = 0; i < n_elt_blockA; i++) {
  
    PDM_g_num_t gnum_boxA = block_gnumA[i];

    PDM_g_num_t _max_key = (2*n_g_newVtxA) % lComm;
    int max_key = (int) _max_key;
    PDM_hash_tab_t *htSubEdgeA = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                                      &max_key);
    
    /* Take into account sub-faces */
    
    for (int j = facesToSubFacesAIdx[i]; j < facesToSubFacesAIdx[i+1]; j++) {
      int nElt = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      int iBeg = subFacesConnecIdx[j];
      int iEnd = subFacesConnecIdx[j+1];
      for (int k = iBeg; k < iEnd; k++) {
        _sub_edge_t *se = malloc(sizeof(_sub_edge_t));
        int next = iBeg + ((k + 1 - iBeg) % nElt);

        int imin = (subFacesConnecA[k] < subFacesConnecA[next]) ?  k : next;
        int imax = (k == imin) ? next : k; 
        
        se->vtx1 = subFacesConnecA[imin];
        se->vtx2 = subFacesConnecA[imax];
        
        for (int k1 = 0; k1 < 3; k1++) {
          se->coords1[k1] = subFacesCoordsA[3*imin+k1]; 
          se->coords2[k1] = subFacesCoordsA[3*imax+k1]; 
        }

//FIXME : key = se->vtx1 + se->vtx2 or key = (se->vtx1 + se->vtx2) % lComm
        
        PDM_g_num_t _key = se->vtx1 + se->vtx2;
        int key = (int) _key;
				PDM_hash_tab_data_add (htSubEdgeA, &key, (void *) se);
      }
    }

    int nVtxA = faceIdxCurrent[0][i + 1] - faceIdxCurrent[0][i];
    PDM_g_num_t *_faceToEdgeA = faceToEdgeCurrent[0] + faceIdxCurrent[0][i];
    PDM_g_num_t *_faceToVtxA  = faceToVtxCurrent[0] + faceIdxCurrent[0][i];
    
    /* Take into account origin edges */

    faceIniVtxIdxA[i+1] = faceIniVtxIdxA[i];
    for (int j = 0; j < nVtxA; j++) {
    
      int n_intersect = -1;        
      PDM_edges_intersect_res_t **eirA = PDM_edges_intersect_get (intersect,
                                                                  PDM_EDGES_GET_FROM_A, 
                                                                  _faceToEdgeA[j],
                                                                  0,
                                                                  &n_intersect);

      int vtx_intersect = 0;
      for (int k = 0; k < n_intersect; k++) {
        
        PDM_line_intersect_t         _tIntersect;

        PDM_g_num_t                   _nGEdgeA;
        PDM_g_num_t                   _originEdgeA;
        int                          _nNewPointsA;
        PDM_edges_intersect_point_t *_oNewPointsA;                               
        PDM_g_num_t                  *_linkA;
        PDM_g_num_t                  *_gNumVtxA;
        double                      *_coordsA; 
        double                      *_uA; 

        /*
         * Get intersection properties
         */
        
        PDM_edges_intersect_res_data_get (eirA[k],        
                                          PDM_EDGES_INTERSECT_MESHA,
                                          &_nGEdgeA,
                                          &_originEdgeA,
                                          &_tIntersect,
                                          &_nNewPointsA,
                                          &_oNewPointsA,                               
                                          &_linkA,
                                          &_gNumVtxA,
                                          &_coordsA,
                                          &_uA);

        vtx_intersect += _nNewPointsA;
      }
      
      int s_tab = vtx_intersect + 2;
      double *u_inter = malloc(sizeof(double) * s_tab);
      double *coords_inter = malloc(sizeof(double) * 3 * s_tab);
      PDM_g_num_t *nG_inter = malloc(sizeof(PDM_g_num_t) * s_tab);

      faceIniVtxIdxA[i+1] += s_tab - 1;
      
      int newSize = faceIniVtxIdxA[i+1];
      if (newSize > s_faceIniVtxA) {
        while (newSize > s_faceIniVtxA) {
          s_faceIniVtxA *= 2;
        }
        faceIniVtxA = realloc (faceIniVtxA, 
                               sizeof(PDM_g_num_t) * s_faceIniVtxA);
      }      
      
      int next = (j+1) % nVtxA;
      
      u_inter[0]       = 0;
      u_inter[s_tab-1] = 1.;
      nG_inter[0]      = _faceToVtxA[j];
      nG_inter[s_tab-1]  = _faceToVtxA[next];
      
      vtx_intersect = 1;
      for (int k = 0; k < n_intersect; k++) {
        
        PDM_line_intersect_t         _tIntersect;

        PDM_g_num_t                   _nGEdgeA;
        PDM_g_num_t                   _originEdgeA;
        int                          _nNewPointsA;
        PDM_edges_intersect_point_t *_oNewPointsA;                               
        PDM_g_num_t                  *_linkA;
        PDM_g_num_t                  *_gNumVtxA;
        double                      *_coordsA; 
        double                      *_uA; 

        /*
         * Get intersection properties
         */

        PDM_edges_intersect_res_data_get (eirA[k],        
                                          PDM_EDGES_INTERSECT_MESHA,
                                          &_nGEdgeA,
                                          &_originEdgeA,
                                          &_tIntersect,
                                          &_nNewPointsA,
                                          &_oNewPointsA,                               
                                          &_linkA,
                                          &_gNumVtxA,
                                          &_coordsA,
                                          &_uA);

        bool reverse = false;
        if (_originEdgeA != _faceToVtxA[j]) {
          assert (_faceToVtxA[next] == _originEdgeA);
          reverse = true; 
        }

        for (int k1 = 0; k1 < _nNewPointsA; k1++) {
        
          u_inter[vtx_intersect]  = reverse ? (1 - _uA[k1]) : _uA[k1];
          nG_inter[vtx_intersect] = _gNumVtxA[k1];
          
          for (int k2 = 0; k2 < 3; k2++) {
            coords_inter[3*vtx_intersect + k2] = _coordsA[3*k1+k2];
          }
                    
          vtx_intersect += 1;
          
        }

      }
            
      int *order = malloc(sizeof(double) * s_tab);
      
      for (int k = 0; k < s_tab; k++) {
        order[k] = k;
      }
           
      PDM_sort_double (u_inter, order, s_tab); 
      
      PDM_g_num_t *nG_inter_tmp = malloc(sizeof(PDM_g_num_t) * s_tab);
      double *coords_inter_tmp = malloc(sizeof(double) * s_tab);
      
      for (int k = 0; k < s_tab; k++) {
        nG_inter_tmp[k] = nG_inter[order[k]];
        for (int k2 = 0; k2 < 3; k2++) {
          coords_inter_tmp[3*k+k2] = coords_inter[3*order[k] + k2];
        }
      }
      
      free (nG_inter);
      free (coords_inter);
      
      nG_inter = nG_inter_tmp;
      coords_inter = coords_inter_tmp;
      
      int k = 0;
      int idx1 = 0;
      while (true) {
        PDM_g_num_t val = nG_inter[k];
        nG_inter[idx1] = nG_inter[k];
        u_inter[idx1] = u_inter[k];
        for (int k2 = 0; k2 < 3; k2++) {
          coords_inter[3*idx1+k2] = coords_inter[3*k + k2];
        }
        idx1 += 1;
        while (nG_inter[k] == val) {
          k += 1;
          if (k >= s_tab) {
            break;
          }
        }
        if (k >= s_tab) {
          break;
        }
      }

      s_tab = idx1;
      
      for (int k1 = 0; k1 < s_tab; k1++) {
        int next1 = (k1 + 1) % s_tab;

        _sub_edge_t *se = malloc(sizeof(_sub_edge_t));

        int imin = (nG_inter[k1] < nG_inter[next1]) ?  k1 : next1;
        int imax = (k1 == imin) ? next1 : k1; 

        se->vtx1 = nG_inter[imin];
        se->vtx2 = nG_inter[imax];
        
        if (next1 != 0) {

//FIXME: verifier si indice est idx ou idx1  et aussi next et next1!!!!

          faceIniVtxA[idx1++] = nG_inter[k1];
        }
                  
        for (int k2 = 0; k2 < 3; k2++) {
          se->coords1[k2] = coords_inter[3*imin+k2];
          se->coords2[k2] = coords_inter[3*imax+k2];
        }

        PDM_g_num_t _key = (se->vtx1 + se->vtx2) % lComm;
        int key = (int) _key;
				PDM_hash_tab_data_add (htSubEdgeA, &key, (void *) se);
      }
      
      free (order);
      free (eirA);
      free (u_inter);
      free (coords_inter);
      free (nG_inter);
      
    }

    /* Count the number of reference of each sub-edge */
    
    int keyMax = *((int *) PDM_hash_tab_keyMax_get (htSubEdgeA));

    int t_n_data = 0;
    for (int j = 0; j < keyMax; j++) {
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeA, &j);
      t_n_data += n_data;
    }
    
    PDM_g_num_t *oneRef = malloc(sizeof(PDM_g_num_t) * 2 * t_n_data);
    double     *coordsOneRef = malloc(sizeof(double) * 6 * t_n_data);
    int nOneRef = 0;
    
    for (int j = 0; j < keyMax; j++) {
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeA, &j);
      int tag[n_data];
      for (int k1 = 0; k1 < n_data; k1++) {
        tag[k1] = 0;
      }
      _sub_edge_t **se = (_sub_edge_t **) PDM_hash_tab_data_get (htSubEdgeA, &j);      
      for (int k1 = 0; k1 < n_data; k1++) {
        if (tag[k1] == 1) continue;
        _sub_edge_t *_se1 = se[k1];
        for (int k2 = k1+1; k2 < n_data; k2++) {
          if (tag[k2] == 1) continue;
          _sub_edge_t *_se2 = se[k2];
          if ((_se1->vtx1 == _se2->vtx1) &&
              (_se1->vtx2 == _se2->vtx2)) {
            tag[k1] = 1;
            tag[k2] = 1;
          }
        }
      }
      for (int k1 = 0; k1 < n_data; k1++) {
        if (tag[k1] == 0) {
          _sub_edge_t *_se1 = se[k1];
          oneRef[2 * nOneRef    ] = _se1->vtx1; 
          oneRef[2 * nOneRef + 1] = _se1->vtx2;
          for (int k2 = 0; k2 < 3; k2++) {
            coordsOneRef[6 * nOneRef + k2] = _se1->coords1[k2];
            coordsOneRef[6 * nOneRef + 3 + k2] = _se1->coords2[k2];
          }
          nOneRef += 1;
        }
      }
    }
    
    PDM_g_num_t *tag = malloc(sizeof(PDM_g_num_t) * nOneRef);
    int nAddSubFace = 0;
    int *addSubFaceIdx = malloc(sizeof(int) * (nOneRef + 1));
    PDM_g_num_t *addSubFace = malloc(sizeof(PDM_g_num_t) * nOneRef);
    double *addSubFaceCoords = malloc(sizeof(double) * 3 * nOneRef);
    
    for (int k1 = 0; k1 < nOneRef; k1++) {
      tag[k1] = 0;
      addSubFaceIdx[k1] = 0;
    }

    /* Add additional sub-faces */

    PDM_g_num_t iniVal = -1;
    PDM_g_num_t nextVal = -1;
    idx = 0;
    int idx2 = 0;
    for (int k1 = 0; k1 < nOneRef; k1++) {
      if (tag[k1] == 1) continue;
      if (iniVal == 0){
        nAddSubFace += 1;
        iniVal = oneRef[2*k1];
        nextVal = oneRef[2*k1 + 1];        
        addSubFace[idx++] = iniVal;
        addSubFace[idx++] = oneRef[2*k1 + 1];
        
        for (int k2 = 0; k2 < 3; k2++) {
          addSubFaceCoords[idx2++] = coordsOneRef[6*k1+k2];
        }
        for (int k2 = 0; k2 < 3; k2++) {
          addSubFaceCoords[idx2++] = coordsOneRef[6*k1+3+k2];
        }
        
        addSubFaceIdx[nAddSubFace] += 2;
      }
      for (int k2 = 0; k2 < nOneRef; k2++) {
        if ((tag[k2] == 1) || (k1 == k2)) continue;
        if (nextVal == oneRef[2*k2]) {
          nextVal = oneRef[2*k2 + 1];
          if (nextVal != iniVal) {
            addSubFace[idx++] = nextVal;
            for (int k3 = 0; k3 < 3; k3++) {
              addSubFaceCoords[idx2++] = coordsOneRef[6*k2+3+k3];
            }
          }
          tag[k2] = 1;
          break;
        }
        else if (nextVal == oneRef[2*k2 + 1]) {
          nextVal = oneRef[2*k2];
          if (nextVal != iniVal) {
            addSubFace[idx++] = nextVal;
            for (int k3 = 0; k3 < 3; k3++) {
              addSubFaceCoords[idx2++] = coordsOneRef[6*k2+3+k3];
            }
          }
          tag[k2] = 1;
          break;
        }
      }
      if (nextVal == iniVal) {
        iniVal = -1;
      }
    }
    
    free (oneRef);
    free (coordsOneRef);
    
    for (int k1 = 0; k1 < nAddSubFace; k1++) {
      addSubFaceIdx[k1+1] += addSubFaceIdx[k1];   
    }
    
    int newSize = iclipp + nAddSubFace + 1;
    if (newSize > s_subFacesConnecIdx) {
      while (newSize > s_subFacesConnecIdx) {
        s_subFacesConnecIdx *= 2;
      }
      subFacesConnecIdx = realloc(subFacesConnecIdx, 
                                sizeof(int) * s_subFacesConnecIdx);
    }      
    
    int ibeg = subFacesConnecIdx[iclipp];
    newSize = ibeg + addSubFaceIdx[nAddSubFace];      
    if (newSize > s_subFacesConnecA) {
      while (newSize > s_subFacesConnecA) {
        s_subFacesConnecA *= 2;
      }
      subFacesConnecA = realloc(subFacesConnecA, 
                                sizeof(PDM_g_num_t) * s_subFacesConnecA);
    }
    
    int ibeg2 = 3 * ibeg;
    newSize = 3 * (ibeg + addSubFaceIdx[nAddSubFace]);      
    if (newSize > s_subFacesCoordsA) {
      while (newSize > s_subFacesCoordsA) {
        s_subFacesCoordsA *= 2;
      }
      subFacesCoordsA = realloc (subFacesCoordsA, 
                                sizeof(double) * s_subFacesCoordsA);
    }
    
    newSize = 4*(iclipp + nAddSubFace);
    if (newSize > s_subFacesToFaces) {
      while (newSize > s_subFacesToFaces) {
        s_subFacesToFaces *= 2;
      }
      subFacesToFaces = realloc(subFacesToFaces, 
                                sizeof(PDM_g_num_t) * s_subFacesToFaces);
    }      
    
    for (int k1 = 0; k1 < addSubFaceIdx[nAddSubFace]; k1++) {
      subFacesConnecA[ibeg++] = addSubFace[k1];   
    }
    
    for (int k1 = 0; k1 < 3*addSubFaceIdx[nAddSubFace]; k1++) {
      subFacesCoordsA[ibeg2++] = addSubFaceCoords[k1]; 
    }

    for (int k1 = 0; k1 < nAddSubFace; k1++) {
      subFacesConnecIdx[iclipp+1] = subFacesConnecIdx[iclipp] +
                                    (addSubFaceIdx[k1+1] - addSubFaceIdx[k1]);
      subFacesToFaces[4*iclipp  ] = i;          // A local number
      subFacesToFaces[4*iclipp+1] = gnum_boxA;  // A global number
      subFacesToFaces[4*iclipp+2] = -1;         // B local number
      subFacesToFaces[4*iclipp+3] = -1;         // B global number
              
    }

    iclipp += nAddSubFace;
    
    facesToAddSubFacesAIdx[i+1] = facesToAddSubFacesAIdx[i] + nAddSubFace;     

    /* Free */

    PDM_hash_tab_free (htSubEdgeA);
    free (tag);
    free (addSubFaceIdx);
    free (addSubFace);
    
  }
  
  const int nSubFaces = iclipp;

  /*
   * Update Memory
   */
    
  faceIniVtxA = realloc (faceIniVtxA, 
                          sizeof(PDM_g_num_t) * faceIniVtxIdxA[n_elt_blockA]);
  
  subFacesConnecA = realloc (subFacesConnecA, 
                             sizeof(PDM_g_num_t) * subFacesConnecIdx[nSubFaces]);
  
  subFacesCoordsA = realloc (subFacesCoordsA, 
                             3 * sizeof(double) * subFacesConnecIdx[nSubFaces]);

  subFacesConnecB = realloc (subFacesConnecB, 
                             sizeof(PDM_g_num_t) * subFacesConnecIdx[nSharedSubFaces]);        
  
  subFacesConnecIdx = realloc (subFacesConnecIdx,
                               sizeof(int) * (nSubFaces + 1));
  
  subFacesToFaces = realloc (subFacesToFaces,
                             sizeof(PDM_g_num_t) *  4 * nSubFaces); 
    
  /*
   * Global numbering definition of sub-faces 
   * This numbering is shared between A and B
   * Partial split faces are stored behind sub-faces
   * Not split faces are stored behind partial split faces
   * 
   */
   
  PDM_g_num_t nAddSubFaces_l = nSubFaces - nSharedSubFaces;
  PDM_g_num_t beg_nAddSubFaces;
  PDM_g_num_t end_subFaces;
  
  PDM_MPI_Gather(&(gNumSubFacesA[nSharedSubFaces-1]), 1, PDM__PDM_MPI_G_NUM, 
             &end_subFaces,                       1, PDM__PDM_MPI_G_NUM, 
             lComm - 1, ol->comm);

  PDM_MPI_Scan (&nAddSubFaces_l, &beg_nAddSubFaces,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);
  
  beg_nAddSubFaces += -nAddSubFaces_l;
  
  for (int i = nSharedSubFaces; i < nSubFaces; i++) {
    gNumSubFacesA[i] = end_subFaces + beg_nAddSubFaces + i + 1; 
  }
  
  PDM_g_num_t nTSubFacesA;
  PDM_MPI_Gather(&(gNumSubFacesA[nSubFaces-1]), 1, PDM__PDM_MPI_G_NUM, 
             &nTSubFacesA,                       1, PDM__PDM_MPI_G_NUM, 
             lComm - 1, ol->comm);
    
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_ADD_SUB_FACESA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_ADD_SUB_FACESA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  
  /*****************************************************************************
   *                                                                           *
   * Compute sub-faces coming from B faces with partial A cover.               *
   *                                                                           *
   * Check that each sub-edge is referenced twice                              *
   *    - Take into sub-edges of sub-faces                                     *
   *    - Take into sub-edges of origin elements                               * 
   *                                                                           *
   * If sub-edges are not referenced twice, the origin face is partial covered *
   * by A. Additional sub-faces have to compute.                               *
   *                                                                           *
   * Two steps :                                                               *
   *    - Redistribute B cells (cs_part_to_block) to merge results             *
   *    - Merge Data                                                           *  
   *                                                                           *
   * TODO: Deplacer cette partie dans une fonction permettant d'enchainer les  * 
   *  communications en mode non-bloquant en plusieurs temps pendant le calcul *
   *  des faces complementaires de A                                           *
   *                                                                           *
   * TODO: Commencer a faire le transfert du resultat du maillage A pendant    *
   * le calcul des sous-facettes complementaires de B                          *
   *                                                                           *
   ****************************************************************************/  

  PDM_part_to_block_t *ptb_boxesB = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_MERGE,
                                                            1.,
                                                            (PDM_g_num_t **) &blockA_boxesB_gnum_data,
                                                            &blockA_boxesB_idx[n_elt_blockA],
                                                            1,
                                                            ol->comm);
  
  /* Send Number of vtx and global Number for each sub-face */
  
  // TODO: optimiser car on transforme tout en PDM_g_num_t alors que 3 valeurs sur
  // 4 sont des entiers
  
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    facesToSubFacesBIdx[i+1] += facesToSubFacesBIdx[i];
  }
  
  int *_tmp_stride = malloc (sizeof(int) * blockA_boxesB_idx[n_elt_blockA]);
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    _tmp_stride[i] = 4 *(facesToSubFacesBIdx[i+1] - facesToSubFacesBIdx[i]);
  }

  int *recv_stride1 = NULL;

  PDM_g_num_t *subFacesConnecN = 
          (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 4* nSharedSubFaces);
       
  originA = PDM_box_set_origin_get (boxesA);

  originB = PDM_box_set_origin_get (boxesB);

  for (int i = 0; i < nSharedSubFaces; i++) {
    int lParentFace  = (int) subFacesToFaces[4*i]; 
    PDM_g_num_t iProc = originA[3*lParentFace]; 
    PDM_g_num_t iPart = originA[3*lParentFace+1];     
    subFacesConnecN[2*i    ] = subFacesConnecIdx[i+1] - subFacesConnecIdx[i];     
    subFacesConnecN[2*i + 1] = iProc;     
    subFacesConnecN[2*i + 2] = iPart;     
    subFacesConnecN[2*i + 3] = gNumSubFacesA[i];     
  }

  int *recvSubFacesConnecBN;  
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(int),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &subFacesConnecN,
                         &recv_stride1,
                         (void **) &recvSubFacesConnecBN);

  /* Send Connectivity for each sub-face */

  int *recv_stride2 = NULL;

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    _tmp_stride[i] = 0;
  }
  
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    for (int j = facesToSubFacesBIdx[i]; j < facesToSubFacesBIdx[i+1]; j++) {
      _tmp_stride[i] += subFacesConnecIdx[j+1] - subFacesConnecIdx[j]; 
    }
  }
  
  PDM_g_num_t *recvSubFacesConnecB;
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &subFacesConnecB,
                         &recv_stride2,
                         (void **) &recvSubFacesConnecB);

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    _tmp_stride[i] *= 3;
  }
  
  int *recv_stride21;
  double *recvSubFacesCoordsB;
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &subFacesCoordsA,
                         &recv_stride21,
                         (void **) &recvSubFacesCoordsB);

  /* Send result of intersection */

  /* Send edge and vertex connectivity */
  
  //
  // FIXME: Verifier que recvFaceToEdgeAndVtx a une utilite (a priori aucune)
  //
  
  
  int *recv_stride3 = NULL;

  int idx1 = 0;
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int lnum_boxB = blockA_boxesB_lnum_data[i];
    _tmp_stride[i] = (faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB]);
    idx1 += _tmp_stride[i];
  }
  
  PDM_g_num_t *sendFaceToEdgeAndVtx = malloc (sizeof(PDM_g_num_t) * idx1);
  
  idx1 = 0;
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int lnum_boxB = blockA_boxesB_lnum_data[i];
    for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {
      sendFaceToEdgeAndVtx[idx1] = faceToEdgeCurrent[1][j];
      idx1 += 1;
    }
  }
  
  PDM_g_num_t *recvFaceToEdgeAndVtx;
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeAndVtx,
                         &recv_stride3,
                         (void **) &recvFaceToEdgeAndVtx);

  free (sendFaceToEdgeAndVtx);

  /* 
   * Send : 
   *   - number of intersection points per edge per element 
   *   - edge vertices
   *   - edge vertices coordinates (in an other exchange)
   */
  	
  int *recv_stride4 = NULL;

  idx1 = idx1 * 3;
  PDM_g_num_t *sendFaceToEdgeNPtInt = malloc (sizeof(PDM_g_num_t) * idx1);

  for (int i = 0; i < idx1; i++) {
    sendFaceToEdgeNPtInt[i] = 0;
  }
            
  idx1 = 0;
  int s_properties = 0;
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int lnum_boxB = blockA_boxesB_lnum_data[i];
    _tmp_stride[i] = (_tmp_stride[i]/2) * 3;
    int iBeg = faceIdxCurrent[1][lnum_boxB];
    int nEdge = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
    for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {
      int n_intersect = -1;        

      int next = iBeg + ((j + 1 - iBeg) % nEdge);

      PDM_g_num_t vtx1 = faceToVtxCurrent[1][j];
      PDM_g_num_t vtx2 = faceToVtxCurrent[1][next];

      PDM_edges_intersect_res_t **eirB = PDM_edges_intersect_get (intersect,
                                                                  PDM_EDGES_GET_FROM_B, 
                                                                  faceToEdgeCurrent[1][j],
                                                                  0,
                                                                  &n_intersect);

      for (int k = 0; k < n_intersect; k++) {

        PDM_line_intersect_t         _tIntersect;

        PDM_g_num_t                   _nGEdgeB;
        PDM_g_num_t                   _originEdgeB;
        int                          _nNewPointsB;
        PDM_edges_intersect_point_t *_oNewPointsB;                               
        PDM_g_num_t                  *_linkB;
        PDM_g_num_t                  *_gNumVtxB;
        double                      *_coordsB; 
        double                      *_uB; 

        PDM_edges_intersect_res_data_get (eirB[k],        
                                          PDM_EDGES_INTERSECT_MESHB,
                                          &_nGEdgeB,
                                          &_originEdgeB,
                                          &_tIntersect,
                                          &_nNewPointsB,
                                          &_oNewPointsB,                               
                                          &_linkB,
                                          &_gNumVtxB,
                                          &_coordsB,
                                          &_uB);
        
        sendFaceToEdgeNPtInt[idx1] += _nNewPointsB;
      }

      s_properties += (int) sendFaceToEdgeNPtInt[idx1++];
      sendFaceToEdgeNPtInt[idx1++] = vtx1;
      sendFaceToEdgeNPtInt[idx1++] = vtx2;
    }
  }
  	
  PDM_g_num_t *recvFaceToEdgeNPtInt;
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeNPtInt,
                         &recv_stride4,
                         (void **) &recvFaceToEdgeNPtInt);

  int *recv_stride41 = NULL;

  idx1 = idx1 * 2;
  double *sendFaceToEdgeCoordsvtx = malloc (sizeof(double) * idx1);
  double *recvFaceToEdgeCoordsvtx = malloc (sizeof(double) * idx1);

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int lnum_boxB = blockA_boxesB_lnum_data[i];
    _tmp_stride[i] = (_tmp_stride[i]/2) * 3;
    int iBeg = faceIdxCurrent[1][lnum_boxB];
    int nEdge = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
    for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {

      int next = iBeg + ((j + 1 - iBeg) % nEdge);

      for (int k1 = 0; k1 < 3; k1++) {        
        sendFaceToEdgeCoordsvtx[idx1++] = faceVtxCooCurrent[1][3*j+k1];
      }

      for (int k1 = 0; k1 < 3; k1++) {        
        sendFaceToEdgeCoordsvtx[idx1++] = faceVtxCooCurrent[1][3*next+k1];
      }
    }
  }
  
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeCoordsvtx,
                         &recv_stride41,
                         (void **) &recvFaceToEdgeCoordsvtx);
  
  /* Send origin, global number and U for each intersection point */

  int *recv_stride5 = NULL;
  int *recv_stride6 = NULL;
  
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    _tmp_stride[i] = 0;
  }
	
  PDM_g_num_t *sendFaceToEdgeOrAndGnumPtInt = 
                                 malloc (sizeof(PDM_g_num_t) * 2 * s_properties);

  double *sendFaceToEdgeUPtInt = malloc (sizeof(double) * 4 * s_properties);
            
  idx1 = 0;
  int idx2 = 0;
  int idx3 = 0;
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    
    int lnum_boxB = blockA_boxesB_lnum_data[i];
    for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {

      _tmp_stride[i] += 2 * (int) sendFaceToEdgeNPtInt[idx1];

      int n_intersect = -1;        
      PDM_edges_intersect_res_t **eirB = PDM_edges_intersect_get (intersect,
                                                                  PDM_EDGES_GET_FROM_B, 
                                                                  faceToEdgeCurrent[1][j],
                                                                  0,
                                                                  &n_intersect);

      for (int k = 0; k < n_intersect; k++) {

        PDM_line_intersect_t         _tIntersect;

        PDM_g_num_t                   _nGEdgeB;
        PDM_g_num_t                   _originEdgeB;
        int                          _nNewPointsB;
        PDM_edges_intersect_point_t *_oNewPointsB;                               
        PDM_g_num_t                  *_linkB;
        PDM_g_num_t                  *_gNumVtxB;
        double                      *_coordsB; 
        double                      *_uB; 

        PDM_edges_intersect_res_data_get (eirB[k],        
                                          PDM_EDGES_INTERSECT_MESHB,
                                          &_nGEdgeB,
                                          &_originEdgeB,
                                          &_tIntersect,
                                          &_nNewPointsB,
                                          &_oNewPointsB,                               
                                          &_linkB,
                                          &_gNumVtxB,
                                          &_coordsB,
                                          &_uB);
        
        for (int k1 = 0; k1 < _nNewPointsB; k1++) {
          sendFaceToEdgeOrAndGnumPtInt[idx2++] = _originEdgeB;
          sendFaceToEdgeOrAndGnumPtInt[idx2++] = _gNumVtxB[k1];
          sendFaceToEdgeUPtInt[idx3++]         = _uB[k1];
          sendFaceToEdgeUPtInt[idx3++]         = _coordsB[3*k1  ];
          sendFaceToEdgeUPtInt[idx3++]         = _coordsB[3*k1+1];
          sendFaceToEdgeUPtInt[idx3++]         = _coordsB[3*k1+2];
        }
      }
      idx1 += 1;
    }
  }

  PDM_g_num_t *recvFaceToEdgeOrAndGnumPtInt;
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeOrAndGnumPtInt,
                         &recv_stride5,
                         (void **) &recvFaceToEdgeOrAndGnumPtInt);
  
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    _tmp_stride[i] = _tmp_stride[i]*2;
  }

  double *recvFaceToEdgeUPtInt;
  PDM_part_to_block_exch (ptb_boxesB, 
                         sizeof(double),
                         PDM_STRIDE_VAR,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeUPtInt,
                         &recv_stride6,
                         (void **) &recvFaceToEdgeUPtInt);
  
  /* Cleanup */
  
  free (_tmp_stride);
  free (subFacesConnecB);
      
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_DISTRIB_BOXESB_BLOCK] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_DISTRIB_BOXESB_BLOCK]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_DISTRIB_BOXESB_BLOCK]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_DISTRIB_BOXESB_BLOCK]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*
   * Loop on faces
   */ 

  int n_elt_blockB = PDM_part_to_block_n_elt_block_get (ptb_boxesB);
  
  PDM_g_num_t *block_gnumB = PDM_part_to_block_block_gnum_get (ptb_boxesB);
 
  int nSubFacesB = 0;

  s_subFacesConnecB = 4 * n_elt_blockB;
  subFacesConnecB = malloc (sizeof(PDM_g_num_t) * s_subFacesConnecB);

  int s_subFacesCoordsB = 3 * 4 * n_elt_blockB;
  double *subFacesCoordsB = malloc (sizeof(double) * s_subFacesCoordsB);

  int s_subFacesConnecIdxB = 4 * n_elt_blockB;
  int *subFacesConnecIdxB = malloc (sizeof(int) * s_subFacesConnecIdxB);
  subFacesConnecIdxB[0] = 0;

  int s_subFacesToFaceB = 4 * n_elt_blockB;
  //FIXME : Le tableau subFacesToFaceB est probablement inutile (A verifier)
  PDM_g_num_t *subFacesToFaceB = malloc (sizeof(PDM_g_num_t) * s_subFacesToFaceB);
  PDM_g_num_t *gNumSubFacesB   = malloc (sizeof(PDM_g_num_t) * s_subFacesToFaceB);
  
  int s_subFacesToLinkA = 4 * n_elt_blockB;  
  PDM_g_num_t *subFacesToLinkA = malloc (sizeof(PDM_g_num_t) * s_subFacesToLinkA);
  
      idx1 = 0;
      idx2 = 0;
      idx3 = 0;
  int idx4 = 0;
  int idx41 = 0;
  int idx5 = 0;
  int idx6 = 0;

  int sIntEdge = 4;
  PDM_g_num_t *vtxIntEdge       = malloc (sizeof(PDM_g_num_t) * sIntEdge);
  PDM_g_num_t *vtxIntEdgeSorted = malloc (sizeof(PDM_g_num_t) * sIntEdge);
  double     *uIntEdge         = malloc (sizeof(double) * sIntEdge);
  double     *coordsIntEdge    = malloc (sizeof(double) * 3 * sIntEdge);
  double     *coordsIntEdgeSorted = malloc (sizeof(double) * 3 * sIntEdge);
  int        *order            = malloc (sizeof(double) * sIntEdge);
  
  facesToSubFacesBIdx = realloc (facesToSubFacesBIdx, sizeof(int) * (n_elt_blockB + 1));
  facesToSubFacesBIdx[0] = 0;
  
  int n_t_nAddSubFace = 0;
  
  int *faceIniVtxIdxB = malloc (sizeof(int) * (1 + n_elt_blockB));
  faceIniVtxIdxB[0] = 0;

  int s_faceIniVtxB = 4 * n_elt_blockB;
  PDM_g_num_t *faceIniVtxB = malloc (sizeof(PDM_g_num_t) * s_faceIniVtxB);  

  idx = 0;  
  for (int i = 0; i < n_elt_blockB; i++) {

    faceIniVtxIdxB[i+1] = faceIniVtxIdxB[i];

    PDM_g_num_t _max_key = (2*n_g_newVtxB) % lComm;
    int max_key = (int) _max_key;
    PDM_hash_tab_t *htSubEdgeB = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                                      &max_key);

    PDM_g_num_t gnum_boxB = block_gnumB[i];
    
    int n_SubFaceFace = recv_stride1[i]/4;
    
    int newSize =  nSubFacesB + n_SubFaceFace;
    if (newSize > s_subFacesToFaceB) {
      while (newSize > s_subFacesToFaceB) {
        s_subFacesToFaceB *= 2;
      }
      subFacesToFaceB = realloc(subFacesToFaceB, 
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
      gNumSubFacesB = realloc(gNumSubFacesB, 
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
    }      

    newSize =  3 * (nSubFacesB + n_SubFaceFace);
    if (newSize > s_subFacesToLinkA) {
      while (newSize > s_subFacesToLinkA) {
        s_subFacesToLinkA *= 2;
      }
      subFacesToLinkA = realloc(subFacesToLinkA, 
                                sizeof(PDM_g_num_t) * s_subFacesToLinkA);
    }      

    newSize =  nSubFacesB + n_SubFaceFace + 1;
    if (newSize > s_subFacesConnecIdxB) {
      while (newSize > s_subFacesConnecIdxB) {
        s_subFacesConnecIdxB *= 2;
      }
      subFacesConnecIdxB = realloc(subFacesConnecIdxB, 
                                   sizeof(int) * s_subFacesConnecIdxB);
    }      

    int ideb  = nSubFacesB;
    int ideb1 = 3 *nSubFacesB;
    int ideb2 = nSubFacesB + 1;
    
    /* Loop on sub-faces */ 
            
    for (int j = 0; j < n_SubFaceFace; j++) {
      
      int        nVtx   = (int) recvSubFacesConnecBN[idx1++];
      int        iProcA = (int) recvSubFacesConnecBN[idx1++];
      int        iPartA = (int) recvSubFacesConnecBN[idx1++];
      PDM_g_num_t _gNumA  =       recvSubFacesConnecBN[idx1++];
     
      subFacesToFaceB[ideb]    = gnum_boxB;
      gNumSubFacesB[ideb++]    = _gNumA;
      
      subFacesToLinkA[ideb1++] = iProcA;
      subFacesToLinkA[ideb1++] = iPartA;
      subFacesToLinkA[ideb1++] = _gNumA;
      
      subFacesConnecIdxB[ideb2] = subFacesConnecIdxB[ideb2-1] + nVtx;
      ideb2 += 1;
      
      /* Loop on sub-edges */ 

      int iBeg = idx2;
      
      newSize =  subFacesConnecIdxB[ideb2];
      if (newSize > s_subFacesConnecB) {
        while (newSize > s_subFacesConnecB) {
          s_subFacesConnecB *= 2;
        }
        subFacesConnecB = realloc(subFacesConnecB, 
                                     sizeof(PDM_g_num_t) * s_subFacesConnecB);
      }      

      newSize = 3 * subFacesConnecIdxB[ideb2];
      if (newSize > s_subFacesCoordsB) {
        while (newSize > s_subFacesCoordsB) {
          s_subFacesCoordsB *= 2;
        }
        subFacesCoordsB = realloc(subFacesCoordsB, 
                                     sizeof(double) * s_subFacesCoordsB);
      }      

      int iBeg1 = subFacesConnecIdxB[ideb2-1];
      int iBeg2 = 3 * subFacesConnecIdxB[ideb2-1];
      for (int k = 0; k < nVtx; k++) {

        int next = iBeg + (k + 1) % nVtx;
        _sub_edge_t *se = malloc(sizeof(_sub_edge_t));
        
        int imin = (recvSubFacesConnecB[idx2] < recvSubFacesConnecB[next]) ?  idx2 : next;
        int imax = (idx2 == imin) ? next : idx2; 

        se->vtx1 = recvSubFacesConnecB[imin];
        se->vtx2 = recvSubFacesConnecB[imax];
        
        for (int k2 = 0; k2 < 3; k2++) {
          se->coords1[k2] = recvSubFacesCoordsB[3*imin+k2];
          se->coords2[k2] = recvSubFacesCoordsB[3*imax+k2];
        }

        PDM_g_num_t _key = (se->vtx1 + se->vtx2) % lComm;
        int key = (int) _key;

				PDM_hash_tab_data_add (htSubEdgeB, &key, (void *) se);

        subFacesConnecB[iBeg1++] = recvSubFacesConnecB[idx2];
        
//FIXME : Verifier cette boucle        
        
        for (int k1 = 0; k1 < 3; k1++) {
          subFacesCoordsB[iBeg2++] = recvSubFacesCoordsB[3*idx2+k1];
        }
        
        idx2 += 1;        
      }
    }

    nSubFacesB += n_SubFaceFace;
    
    /* 
     * Add additional B sub-faces
     * Take into account initial edges 
     */    
    
    int nEdge = recv_stride3[i];
    
    for (int j = 0; j < nEdge; j++) {
      
      int nIntEdge = 0;
      
//      PDM_g_num_t gEdge = recvFaceToEdgeAndVtx[idx3++]; // Ca sert a qqch ?
      
      int        nInt  = (int) recvFaceToEdgeNPtInt[idx4++];      
      PDM_g_num_t nVtx1 = recvFaceToEdgeNPtInt[idx4++];
      PDM_g_num_t nVtx2 = recvFaceToEdgeNPtInt[idx4++];
      
      vtxIntEdge[nIntEdge] = nVtx1;
      
      for (int k = 0; k < 3; k++) {
        coordsIntEdge[3 * nIntEdge + k] = recvFaceToEdgeCoordsvtx[idx41++]; 
      }
      
      uIntEdge[nIntEdge++] = 0.;
      
      vtxIntEdge[nIntEdge] = nVtx2;

      for (int k = 0; k < 3; k++) {
        coordsIntEdge[3 * nIntEdge + k] = recvFaceToEdgeCoordsvtx[idx41++]; 
      }

      uIntEdge[nIntEdge++] = 1.;

      for (int k = 0; k < nInt; k++) {
        
        PDM_g_num_t _originEdgeB = recvFaceToEdgeOrAndGnumPtInt[idx5++];
        PDM_g_num_t _gnumB = recvFaceToEdgeOrAndGnumPtInt[idx5++];
        double u          = recvFaceToEdgeUPtInt[idx6++];
        double coords[3]  = {recvFaceToEdgeUPtInt[idx6],
                             recvFaceToEdgeUPtInt[idx6+1],
                             recvFaceToEdgeUPtInt[idx6+2]};
        idx6 += 3;
        
        if (nVtx1 != _originEdgeB) {
          assert (nVtx2 != _originEdgeB);
          u = 1 - u;
        }
        
        if (nIntEdge >= sIntEdge) {
          while (nIntEdge >= sIntEdge) {
            sIntEdge *= 2;
          }
          vtxIntEdge = realloc (vtxIntEdge, sizeof(PDM_g_num_t) * sIntEdge);
          vtxIntEdgeSorted = realloc (vtxIntEdgeSorted, sizeof(PDM_g_num_t) * sIntEdge);
          uIntEdge   = realloc (uIntEdge,   sizeof(double) * sIntEdge);
          coordsIntEdge   = realloc (uIntEdge,   sizeof(double) * 3 * sIntEdge);
          coordsIntEdgeSorted   = realloc (uIntEdge,   sizeof(double) * 3 * sIntEdge);
          order      = realloc (order,      sizeof(int) * sIntEdge);
        }
        
        vtxIntEdge[nIntEdge] = _gnumB;
        uIntEdge[nIntEdge]   = u;
       
        for (int k1 = 0; k1 < 3; k1++) {
          coordsIntEdge[3*nIntEdge+k] = coords[k1];
        }
        
        nIntEdge += 1;
      }
      
      faceIniVtxIdxB[i+1] += nIntEdge - 1; 
      
      int _newSize = faceIniVtxIdxB[i+1];
      if (_newSize > s_faceIniVtxB) {
        while (_newSize > s_faceIniVtxB) {
          s_faceIniVtxB *= 2;
        }
        faceIniVtxB = realloc (faceIniVtxB, 
                               sizeof(PDM_g_num_t) * s_faceIniVtxB);
      }      
      
      for (int k = 0; k < nIntEdge; k++) {
        order[k] = k;
      }
           
      PDM_sort_double (uIntEdge, order, nIntEdge); 

      for (int k = 0; k < nIntEdge; k++) {
        vtxIntEdgeSorted[k] = vtxIntEdge[order[k]];
        for (int k1 = 0; k1 < 3; k1++) {
          coordsIntEdgeSorted[3*k+k1] = coordsIntEdge[3*order[k]+k1];
        }
      }

      int k1 = 0;
      int K2 = 0;
      
      while (true) {
        PDM_g_num_t val = vtxIntEdgeSorted[k1];
        vtxIntEdgeSorted[K2] = vtxIntEdgeSorted[k1];
        uIntEdge[K2] = uIntEdge[k1];
        for (int k3 = 0; k3 < 3; k3++) {
          coordsIntEdgeSorted[3*K2+k3] = coordsIntEdgeSorted[3*k1+k3]; 
        }
        K2 += 1;
        while (vtxIntEdgeSorted[k1] == val) {
          k1 += 1;
          if (k1 >= nIntEdge) {
            break;
          }
        }
        if (k1 >= nIntEdge) {
          break;
        }
      }
      
      nIntEdge = K2;

      for (int k3 = 0; k3 < nIntEdge; k3++) {
        int next = (k3 + 1) % nIntEdge;

        _sub_edge_t *se = malloc(sizeof(_sub_edge_t));

        int imin = (vtxIntEdgeSorted[k3] < vtxIntEdgeSorted[next]) ?  k3 : next;
        int imax = (k3 == imin) ? next : k3; 

        if (next != 0) {
          faceIniVtxB[idx++] = vtxIntEdgeSorted[k3];
        }
        
        se->vtx1 = vtxIntEdgeSorted[imin];
        se->vtx2 = vtxIntEdgeSorted[imax];
        
        for (int k2 = 0; k2 < 3; k2++) {
          se->coords1[k2] = coordsIntEdgeSorted[3*imin+k2];
          se->coords2[k2] = coordsIntEdgeSorted[3*imax+k2];
        }
        
        PDM_g_num_t _key = se->vtx1 + se->vtx2;
        int key = (int) _key;
                
				PDM_hash_tab_data_add (htSubEdgeB, &key, (void *) se);
      }
      
    }
    
    /* Count the number of reference of each sub-edge */
    
    int keyMax = *((int *) PDM_hash_tab_keyMax_get (htSubEdgeB));

    int t_n_data = 0;
    for (int j = 0; j < keyMax; j++) {
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeB, &j);
      t_n_data += n_data;
    }
    
    PDM_g_num_t *oneRef = malloc(sizeof(PDM_g_num_t) * 2 * t_n_data);
    double     *coordsOneRef = malloc(sizeof(double) * 6 * t_n_data);
    
    int nOneRef = 0;
    
    for (int j = 0; j < keyMax; j++) {
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeB, &j);
      int tag[n_data];
      for (int k1 = 0; k1 < n_data; k1++) {
        tag[k1] = 0;
      }
      _sub_edge_t **se = (_sub_edge_t **) PDM_hash_tab_data_get (htSubEdgeB, &j);      
      for (int k1 = 0; k1 < n_data; k1++) {
        if (tag[k1] == 1) continue;
        _sub_edge_t *_se1 = se[k1];
        for (int k2 = k1+1; k2 < n_data; k2++) {
          if (tag[k2] == 1) continue;
          _sub_edge_t *_se2 = se[k2];
          if ((_se1->vtx1 == _se2->vtx1) &&
              (_se1->vtx2 == _se2->vtx2)) {
            tag[k1] = 1;
            tag[k2] = 1;
          }
        }
      }
      for (int k1 = 0; k1 < n_data; k1++) {
        if (tag[k1] == 0) {
          _sub_edge_t *_se1 = se[k1];
          oneRef[2 * nOneRef    ] = _se1->vtx1; 
          oneRef[2 * nOneRef + 1] = _se1->vtx2; 

          for (int k2 = 0; k2 < 3; k2++) {
            coordsOneRef[6 * nOneRef + k2] = _se1->coords1[k2];
            coordsOneRef[6 * nOneRef + 3 + k2] = _se1->coords2[k2];
          }

          nOneRef += 1;
        }
      }
    }
    
    PDM_g_num_t *tag = malloc(sizeof(PDM_g_num_t) * nOneRef);
    int nAddSubFace = 0;
    int *addSubFaceIdx = malloc(sizeof(int) * (nOneRef + 1));
    PDM_g_num_t *addSubFace = malloc(sizeof(PDM_g_num_t) * nOneRef);
    double *addSubFaceCoords = malloc(sizeof(double) * 3 * nOneRef);
    for (int k1 = 0; k1 < nOneRef; k1++) {
      tag[k1] = 0;
      addSubFaceIdx[k1] = 0;
    }

    /* Add additional sub-faces */

    PDM_g_num_t iniVal = -1;
    PDM_g_num_t nextVal = -1;
    idx = 0;
    idx2 = 0;
    for (int k1 = 0; k1 < nOneRef; k1++) {
      if (tag[k1] == 1) continue;
      if (iniVal == 0){
        nAddSubFace += 1;
        iniVal = oneRef[2*k1];
        nextVal = oneRef[2*k1 + 1];        
        
        addSubFace[idx++] = iniVal;
        addSubFace[idx++] = oneRef[2*k1 + 1];
        for (int k3 = 0; k3 < 3; k3++) {
          addSubFaceCoords[idx2++] = coordsOneRef[6*k1+k3];
        }
        for (int k3 = 0; k3 < 3; k3++) {
          addSubFaceCoords[idx2++] = coordsOneRef[6*k1+3+k3];
        }
        
        addSubFaceIdx[nAddSubFace] += 2;
      }
      for (int k2 = 0; k2 < nOneRef; k2++) {
        if ((tag[k2] == 1) || (k1 == k2)) continue;
        if (nextVal == oneRef[2*k2]) {
          nextVal = oneRef[2*k2 + 1];
          if (nextVal != iniVal) {
            addSubFace[idx++] = nextVal;
            for (int k3 = 0; k3 < 3; k3++) {
              addSubFaceCoords[idx2++] = coordsOneRef[6*k2+3+k3];
            }
          }
          tag[k2] = 1;
          break;
        }
        else if (nextVal == oneRef[2*k2 + 1]) {
          nextVal = oneRef[2*k2];
          if (nextVal != iniVal) {
            addSubFace[idx++] = nextVal;
            for (int k3 = 0; k3 < 3; k3++) {
              addSubFaceCoords[idx2++] = coordsOneRef[6*k2+3+k3];
            }
          }
          tag[k2] = 1;
          break;
        }
      }
      if (nextVal == iniVal) {
        iniVal = -1;
      }
    }
    
    for (int k1 = 0; k1 < nAddSubFace; k1++) {
      addSubFaceIdx[k1+1] += addSubFaceIdx[k1];   
    }
    
    newSize = nSubFacesB + nAddSubFace + 1;
    if (newSize > s_subFacesConnecIdxB) {
      while (newSize > s_subFacesConnecIdxB) {
        s_subFacesConnecIdxB *= 2;
      }
      subFacesConnecIdxB = realloc(subFacesConnecIdxB, 
                                   sizeof(int) * s_subFacesConnecIdxB);
    }

    newSize =  nSubFacesB + nAddSubFace;
    if (newSize > s_subFacesToFaceB) {
      while (newSize > s_subFacesToFaceB) {
        s_subFacesToFaceB *= 2;
      }
      subFacesToFaceB = realloc(subFacesToFaceB, 
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
      gNumSubFacesB = realloc(gNumSubFacesB, 
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
    }

    newSize =  3 * (nSubFacesB + nAddSubFace);
    if (newSize > s_subFacesToLinkA) {
      while (newSize > s_subFacesToLinkA) {
        s_subFacesToLinkA *= 2;
      }
      subFacesToLinkA = realloc(subFacesToLinkA, 
                                sizeof(PDM_g_num_t) * s_subFacesToLinkA);
    }      

    n_t_nAddSubFace += nAddSubFace;
    for (int j = 0; j < nAddSubFace; j++) {

      int nVtx = addSubFaceIdx[j+1] - addSubFaceIdx[j];
      subFacesToFaceB[nSubFacesB] = gnum_boxB;
      gNumSubFacesB[nSubFacesB] = -1;
      
      subFacesToLinkA[3*nSubFacesB] = -1;
      subFacesToLinkA[3*nSubFacesB+1] = -1;
      subFacesToLinkA[3*nSubFacesB+2] = -1;

      subFacesConnecIdxB[nSubFacesB++] += nVtx;

      newSize = subFacesConnecIdxB[nSubFacesB];
      if (newSize > s_subFacesConnecB) {
        while (newSize > s_subFacesConnecB) {
          s_subFacesConnecB *= 2;
        }
        subFacesConnecB = realloc(subFacesConnecB, 
                                  sizeof(PDM_g_num_t) * s_subFacesConnecB);
      }      
      
      newSize = 3 * subFacesConnecIdxB[nSubFacesB];
      if (newSize > s_subFacesCoordsB) {
        while (newSize > s_subFacesCoordsB) {
          s_subFacesCoordsB *= 2;
        }
        subFacesCoordsB = realloc(subFacesCoordsB, 
                                     sizeof(double) * s_subFacesCoordsB);
      }      

      int _ideb2 =     subFacesConnecIdxB[nSubFacesB-1];
      int ideb3 = 3 * subFacesConnecIdxB[nSubFacesB-1];
      
      for (int k = addSubFaceIdx[j]; k < addSubFaceIdx[j+1]; k++) {
        subFacesConnecB[_ideb2++] = addSubFace[j]; 
        for (int k1 = 0; k1 < 3; k1++) {
          subFacesCoordsB[ideb3++] = addSubFaceCoords[3*j+k1];
        }
      }
    }

    facesToSubFacesBIdx[i+1] = facesToSubFacesBIdx[i] 
                             + n_SubFaceFace
                             + nAddSubFace;    
    /* Cleanup */
    
    free (oneRef);
    free (coordsOneRef);
    free (htSubEdgeB);
    free (addSubFaceIdx);
    free (addSubFace);
    free (addSubFaceCoords);
    
  }

  /*
   * Update memory
   */
  
  subFacesConnecB = realloc (subFacesConnecB, sizeof(PDM_g_num_t) * subFacesConnecIdxB[nSubFacesB]);
  subFacesCoordsB = realloc (subFacesCoordsB, sizeof(double) * 3 * subFacesConnecIdxB[nSubFacesB]);
  subFacesConnecIdxB = realloc (subFacesConnecIdxB, sizeof(int) * (nSubFacesB + 1));
  subFacesToFaceB = realloc (subFacesToFaceB, sizeof(PDM_g_num_t) * nSubFacesB);
  gNumSubFacesB = realloc (gNumSubFacesB, sizeof(PDM_g_num_t) * nSubFacesB);
  subFacesToLinkA = realloc (subFacesToLinkA, 3 * nSubFacesB); 
  
  PDM_MPI_Scan (&n_t_nAddSubFace, &beg_nAddSubFaces,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  beg_nAddSubFaces += -n_t_nAddSubFace + 1;
  int j1 = 0;
  for (int i = 0; i < nSubFacesB; i++) {
    if (gNumSubFacesB[i] == -1) {
      gNumSubFacesB[i] = beg_nAddSubFaces + j1++;
    }
  }
  
  PDM_g_num_t nTSubFacesB;
  PDM_MPI_Gather(&(gNumSubFacesB[nSubFaces-1]), 1, PDM__PDM_MPI_G_NUM, 
             &nTSubFacesB,                       1, PDM__PDM_MPI_G_NUM, 
             lComm - 1, ol->comm);

  /*
   * Cleanup
   */

  free (vtxIntEdge);
  free (coordsIntEdge);
  free (coordsIntEdgeSorted);
  free (vtxIntEdgeSorted);
  free (uIntEdge);
  free (order);

  // sous-facettes calculees
  
  free (recv_stride1);
  free (subFacesConnecN);
  free (recvSubFacesConnecBN);

  free (recv_stride2);
  free (recvSubFacesConnecB);

  free (recv_stride21);
  free (recvSubFacesCoordsB);

  // arete initiales du contour

  free (recv_stride3);
  free (sendFaceToEdgeAndVtx);
  free (recvFaceToEdgeAndVtx);

  free (recv_stride4);
  free (sendFaceToEdgeNPtInt);
  free (recvFaceToEdgeNPtInt);
  
  free (recv_stride41);
  free (sendFaceToEdgeCoordsvtx);
  free (recvFaceToEdgeCoordsvtx);
  
  free (recv_stride5);
  free (sendFaceToEdgeOrAndGnumPtInt);
  free (recvFaceToEdgeOrAndGnumPtInt);

  free (recv_stride6);
  free (sendFaceToEdgeUPtInt);
  free (recvFaceToEdgeUPtInt);
  
  free (blockA_boxesB_idx);
  free (blockA_boxesB_gnum_data);
  free (blockA_boxesB_lnum_data);
  
  for (int i = 0; i < 2; i++) {
    free (faceStrideCurrent[i]);
    free (faceStrideCurrent3[i]);
    free (faceToEdgeCurrent[i]); 
    free (faceToVtxCurrent[i]);
    free (faceVtxCooCurrent[i]);
    free (faceVtxEpsCurrent[i]);
  }
        
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESB] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESB]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_ADD_SUB_FACESB]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_ADD_SUB_FACESB]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  
  /*****************************************************************************
   *                                                                           *
   * Transfer results to origin distribution for mesh A :                      *
   * For each boxe :                                                           *
   *     - Size of connectivity of each subface + numabs + grpah comm (iproc ipart) *
   *     - Connectivity of each subface (numabs)                               *
   *     - Coordinates of vertices                                             *
   *                                                                           *
   * For each face, additional faces and sub-faces have to be merge            *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  PDM_g_num_t *firstSend = malloc (sizeof(PDM_g_num_t) * 
                          facesToAddSubFacesAIdx[n_elt_blockA] * 5);
  
  int        *firstSendStride = malloc (sizeof(int) * n_elt_blockA);
  
  int **firstRecvStrideA = malloc (sizeof(int *) * nPartA);
  
  for (int i = 0; i < nPartA; i++) {
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    firstRecvStrideA[i] = malloc(sizeof(int) * nFace);
  }
  
  PDM_g_num_t **firstRecvA = malloc (sizeof(PDM_g_num_t *) * nPartA);
  
  idx = 0;
  int n_T_vertex = 0;
  
  //  PDM_g_num_t *block_gnumA = PDM_part_to_block_block_gnum_get (ptb_boxesA);
  
  for (int i = 0; i < n_elt_blockA; i++) {
    
    int nAddSubFaces = facesToAddSubFacesAIdx[i+1] - facesToAddSubFacesAIdx[i];
    int _nSubFaces = facesToSubFacesAIdx[i+1] - facesToSubFacesAIdx[i];
    
    firstSendStride[i] =  5 * (_nSubFaces + nAddSubFaces);
    
    for (int j = facesToSubFacesAIdx[i]; j < facesToSubFacesAIdx[i+1]; j++) {
      int lParentFaceB  = (int) subFacesToFaces[4*j+3];
 
      int nVtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      n_T_vertex += nVtx;
      
      PDM_g_num_t iProcB = originB[3*lParentFaceB]; 
      PDM_g_num_t iPartB = originB[3*lParentFaceB+1];     
      PDM_g_num_t numAbs = gNumSubFacesA[j];
      PDM_g_num_t oNumAbs = block_gnumA[i];
      
      firstSend[idx++] = nVtx;
      firstSend[idx++] = numAbs;
      firstSend[idx++] = oNumAbs;
      firstSend[idx++] = iProcB;
      firstSend[idx++] = iPartB;
  
    }
    
    for (int j = facesToAddSubFacesAIdx[i]; j < facesToAddSubFacesAIdx[i+1]; j++) {

      int nVtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j]; 
      PDM_g_num_t iProcB = -1; 
      PDM_g_num_t iPartB = -1;     
      PDM_g_num_t numAbs = gNumSubFacesA[j];
      PDM_g_num_t oNumAbs = block_gnumA[i];
      
      firstSend[idx++] = nVtx;
      firstSend[idx++] = numAbs;
      firstSend[idx++] = oNumAbs;
      firstSend[idx++] = iProcB;
      firstSend[idx++] = iPartB;
  
    }
  }
  
  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride, 
                                           firstSend, 
                                           firstRecvStrideA, 
                                           (void ** )firstRecvA); 
  
  free (firstSend);
  PDM_g_num_t *secondSend = malloc (sizeof(PDM_g_num_t) * n_T_vertex);
  double     *thirdSend = malloc (3 * sizeof(double) * n_T_vertex);

  idx = 0;
  idx2 = 0;
  for (int i = 0; i < n_elt_blockA; i++) {
    
    firstSendStride[i] = 0;
    
    for (int j = facesToSubFacesAIdx[i]; j < facesToSubFacesAIdx[i+1]; j++) {
      
      int nVtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j]; 
      firstSendStride[i] += nVtx;
      
      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        secondSend[idx++] = subFacesConnecA[k];
        for (int k1 = 0; k1 < 3; k1++) {
          thirdSend[idx2++] = subFacesCoordsA[3*k+k1];
        }
      }

    }
    
    for (int j = facesToAddSubFacesAIdx[i]; j < facesToAddSubFacesAIdx[i+1]; j++) {

      int nVtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j]; 
      firstSendStride[i] += nVtx;
      
      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        secondSend[idx++] = subFacesConnecA[k];
        for (int k1 = 0; k1 < 3; k1++) {
          thirdSend[idx2++] = subFacesCoordsA[3*k+k1];
        }
      }
    }
  }
  
  int **secondRecvStrideA = malloc (sizeof(int *) * nPartA);
  
  for (int i = 0; i < nPartA; i++) {
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    secondRecvStrideA[i] = malloc(sizeof(int) * nFace);
  }
    
  PDM_g_num_t **secondRecvA = malloc (sizeof(PDM_g_num_t *) * nPartA);

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride, 
                                           secondSend, 
                                           secondRecvStrideA, 
                                           (void ** ) secondRecvA); 
  
  for (int i = 0; i < n_elt_blockA; i++) {
    firstSendStride[i] *= 3;
  }
  
  int **thirdRecvStrideA = malloc (sizeof(int *) * nPartA);
  
  for (int i = 0; i < nPartA; i++) {
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    thirdRecvStrideA[i] = malloc(sizeof(int) * nFace);
  }
    
  double **thirdRecvA = malloc (sizeof(double *) * nPartA);

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(double),
                                           firstSendStride, 
                                           thirdSend, 
                                           thirdRecvStrideA, 
                                           (void ** )thirdRecvA); 
  
  PDM_g_num_t *fourthSend = faceIniVtxA;

  int **fourthRecvStrideA = malloc (sizeof(int *) * nPartA);
  PDM_g_num_t **fourthRecvA = malloc (sizeof(PDM_g_num_t *) * nPartA);
  
  for (int i = 0; i < n_elt_blockA; i++) {
    
    firstSendStride[i] = faceIniVtxIdxA[i+1] - faceIniVtxIdxA[i];
    
  }

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride, 
                                           fourthSend, 
                                           fourthRecvStrideA, 
                                           (void ** )fourthRecvA); 
  
  free (firstSend);
  free (firstSendStride);
  free (secondSend);
  free (thirdSend);
  free (fourthSend);

  free (subFacesConnecIdx);  
  free (subFacesCoordsA);
  free (facesToSubFacesAIdx);
  free (facesToAddSubFacesAIdx);
  free (gNumSubFacesA);
  free (subFacesToFaces);
  free (subFacesConnecA);
  
  PDM_box_set_destroy (&boxesA);
          
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_SEND_RESULTS_TO_INIT_PARTA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_SEND_RESULTS_TO_INIT_PARTA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*****************************************************************************
   *                                                                           *
   * Redistribute B Boxes to send results to B                                 *
   *                                                                           *
   ****************************************************************************/

  PDM_l_num_t *destinationB = PDM_part_to_block_destination_get (ptb_boxesB);
  
  
  n_g_eltB = PDM_box_set_get_global_size(boxesB);

  distribB = PDM_box_distrib_create (n_eltB,
                                     n_g_eltB,
                                     1, // Don't use in this case
                                     ol->comm);

  countEltsB = (int *) malloc (sizeof(int) * n_eltB);
  
  for (int i = 0; i < lComm + 1; i++) {
    distribB->index[i] = 0;
    countEltsB[i] = 0;
  }
  
  for (int i = 0; i < n_elt_blockB; i++) {
    int iProc = destinationB[i] + 1;
    distribB->index[iProc]++;
  }

  for (int i = 0; i < lComm; i++) {
    distribB->index[i+1] += distribB->index[i];
  }
  
  distribB->list = (int *) malloc (sizeof(int) * distribB->index[lComm]);
  
  for (int i = 0; i < n_elt_blockB; i++) {
    int iProc = destinationB[i] + 1;
    int idxB = distribB->index[iProc] + (countEltsB[iProc]++);
    distribB->list[idxB] = i;
  }
  
  free (countEltsB);

  PDM_box_set_redistribute (distribB, boxesB);
    
  /*
   * Cleanup
   */

  PDM_box_distrib_destroy (&distribB);
  
  PDM_part_to_block_free (ptb_boxesB);

  /*****************************************************************************
   *                                                                           *
   * Transfer results to origin distribution for mesh B                        *
   *                                                                           *
   ****************************************************************************/

  firstSend = malloc (sizeof(PDM_g_num_t) * 
                      facesToSubFacesBIdx[n_elt_blockB] * 5);
  
  firstSendStride = malloc (sizeof(int) * n_elt_blockA);
  
  int **firstRecvStrideB = malloc (sizeof(int *) * nPartB);
  
  for (int i = 0; i < nPartB; i++) {
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    firstRecvStrideB[i] = malloc(sizeof(int) * nFace);
  }
  
  PDM_g_num_t **firstRecvB = malloc (sizeof(PDM_g_num_t *) * nPartB);
  
  idx = 0;
  n_T_vertex = 0;
    
  for (int i = 0; i < n_elt_blockB; i++) {
    
    int _nSubFaces = facesToSubFacesBIdx[i+1] - facesToSubFacesBIdx[i];

    firstSendStride[i] = 5 * _nSubFaces;
    
    for (int j = facesToSubFacesBIdx[i]; j < facesToSubFacesBIdx[i+1]; j++) {

      int nVtx = subFacesConnecIdxB[j+1] - subFacesConnecIdxB[j];
      n_T_vertex += nVtx;
      int iProcA = (int) subFacesToLinkA[3*j];
      int iPartA = (int) subFacesToLinkA[3*j + 1];
      PDM_g_num_t numAbs = gNumSubFacesB[j];
      PDM_g_num_t oNumAbs = block_gnumB[i];

      firstSend[idx++] = nVtx;
      firstSend[idx++] = numAbs;
      firstSend[idx++] = oNumAbs;
      firstSend[idx++] = iProcA;
      firstSend[idx++] = iPartA;

    }
  }
  
  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride, 
                                           firstSend, 
                                           firstRecvStrideB, 
                                           (void ** )firstRecvB); 
 
  
  free (firstSend);

  secondSend = malloc (sizeof(PDM_g_num_t) * n_T_vertex);
  thirdSend = malloc (3 * sizeof(double) * n_T_vertex);

  idx = 0;
  idx2 = 0;

  for (int i = 0; i < n_elt_blockB; i++) {

    firstSendStride[i] = 0;
    
    for (int j = facesToSubFacesBIdx[i]; j < facesToSubFacesBIdx[i+1]; j++) {

      int nVtx = subFacesConnecIdxB[j+1] - subFacesConnecIdxB[j];

      firstSendStride[i] += nVtx;
      
      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        secondSend[idx++] = subFacesConnecB[k];
        for (int k1 = 0; k1 < 3; k1++) {
          thirdSend[idx2++] = subFacesCoordsB[3*k+k1];
        }
      }
    }
  }
  
  int **secondRecvStrideB = malloc (sizeof(int *) * nPartB);
  
  for (int i = 0; i < nPartB; i++) {
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    secondRecvStrideB[i] = malloc(sizeof(int) * nFace);
  }
    
  PDM_g_num_t **secondRecvB = malloc (sizeof(PDM_g_num_t *) * nPartB);

  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride, 
                                           secondSend, 
                                           secondRecvStrideB, 
                                           (void ** )secondRecvB); 
  
  for (int i = 0; i < n_elt_blockA; i++) {
    firstSendStride[i] *= 3;
  }
  
  int **thirdRecvStrideB = malloc (sizeof(int *) * nPartB);
  
  for (int i = 0; i < nPartB; i++) {
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    thirdRecvStrideB[i] = malloc(sizeof(int) * nFace);
  }
    
  double **thirdRecvB = malloc (sizeof(double *) * nPartB);

  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(double),
                                           firstSendStride, 
                                           thirdSend, 
                                           thirdRecvStrideB, 
                                           (void ** )thirdRecvB);
    
  fourthSend = faceIniVtxB;

  int **fourthRecvStrideB = malloc (sizeof(int *) * nPartB);
  PDM_g_num_t **fourthRecvB = malloc (sizeof(PDM_g_num_t *) * nPartB);
  
  for (int i = 0; i < n_elt_blockB; i++) {
    
    firstSendStride[i] = faceIniVtxIdxB[i+1] - faceIniVtxIdxB[i];
    
  }

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride, 
                                           fourthSend, 
                                           fourthRecvStrideB, 
                                           (void ** )fourthRecvB); 
  
  free (firstSendStride);
  free (secondSend);
  free (thirdSend);  

  free (subFacesConnecB);
  free (subFacesConnecIdxB);
  
  free (subFacesToFaceB);  
  free (gNumSubFacesB);  
  free (subFacesToLinkA);
  
  free (facesToSubFacesBIdx);
    
  PDM_box_set_destroy (&boxesB);
  
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTB] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTB]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_SEND_RESULTS_TO_INIT_PARTB]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_SEND_RESULTS_TO_INIT_PARTB]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  
  /*****************************************************************************
   *                                                                           *
   * Compute A local numbering for vertices and faces                          * 
   *                                                                           *
   ****************************************************************************/

  PDM_g_num_t nUnChangedFaceA = 0;
  PDM_g_num_t nSubFaceA = 0;
  int *nUnChangedFacePartA = malloc(sizeof(int) * nPartA);
  int *nSubFacePartA = malloc(sizeof(int) * nPartA);
  int *s_olFaceVtxA = malloc(sizeof(int) * nPartA);

  /* 
   * Compute dimensions
   */

  for (int i = 0; i < nPartA; i++) {
    nUnChangedFacePartA[i] = 0;
    nSubFacePartA[i] = 0;
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    const PDM_g_num_t *gNumFaceA = PDM_surf_mesh_part_face_g_num_get (ol->meshA, i);
    const int *iniFaceVtxIdxA =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshA, i);

    idx = 0;
    for (int j = 0; j < nFace; j++) {
      if (firstRecvStrideA[i][j] == 0) {
        nUnChangedFacePartA[i] += 1;
        s_olFaceVtxA[i] += iniFaceVtxIdxA[j+1] - iniFaceVtxIdxA[j];
      }
      else {
        nSubFacePartA[i] += firstRecvStrideA[i][j]/5;
        for (int k = 0; k < firstRecvStrideA[i][j]; k++) {
          int nVtx    = (int) firstRecvA[i][idx++];
          idx++;
          PDM_g_num_t oNumAbs = firstRecvA[i][idx++];
          idx += 2;

          s_olFaceVtxA[i] += nVtx;

          assert(oNumAbs == gNumFaceA[j]);
        }
      }
    }
    
    nUnChangedFaceA += nUnChangedFacePartA[i];
    nSubFaceA += nSubFacePartA[i];

  }
  
  PDM_g_num_t begUnChangedFaceA = 0;
  PDM_MPI_Scan (&nUnChangedFaceA, &begUnChangedFaceA,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  begUnChangedFaceA += -nUnChangedFaceA + 1;

  PDM_g_num_t nFaceA = nUnChangedFaceA + nSubFaceA;
  PDM_g_num_t nGFaceA;
  PDM_MPI_Allreduce(&nFaceA, &nGFaceA, 1, 
                PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm); 

  PDM_g_num_t nGVtxA = n_g_newVtxA;
  ol->olMeshA = _ol_mesh_create (nGFaceA, nGVtxA, nPartA);

  /* 
   * Build local connectivity 
   */
  
  int idxIpart = 0;
  
  PDM_g_num_t iniNGVtx = PDM_surf_mesh_n_g_vtx_get (ol->meshA);
  
  for (int i = 0; i < nPartA; i++) {
    
    _ol_part_t *olp = _ol_part_create (meshA->part[i]);
    ol->olMeshA->part[i] = olp;
     
    olp->sInitToOlFace   = 0;
    olp->initToOlFaceIdx = NULL;       
    olp->initToOlFace    = NULL;
    olp->nLinkedFace     = 0;
    olp->linkedFaces     = NULL;
    
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    const int *iniFaceVtxIdxA =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshA, i);
    const int *iniFaceVtxA =  PDM_surf_mesh_part_face_vtx_get (ol->meshA, i);

    int nFacePart = nUnChangedFacePartA[i] + nSubFacePartA[i];
    olp->sInitToOlFace = nFacePart;
    olp->nLinkedFace = nSubFacePartA[i];

    olp->initToOlFaceIdx = malloc (sizeof(int) * (nFace + 1));
    olp->initToOlFaceIdx[0] = 0;
    
    olp->initToOlFace = malloc (sizeof(int) * nFacePart);
    olp->linkedFaces =  malloc (sizeof(int) * 4 * nSubFacePartA[i]); 
    
    int *faceVtxIdxPart = malloc (sizeof(int) * (nFacePart + 1));
    faceVtxIdxPart[0] = 0;
    
    PDM_g_num_t *gFaceVtxPart = malloc (sizeof(PDM_g_num_t) * s_olFaceVtxA[i]);
    int *faceVtxPart = malloc(sizeof(int) * s_olFaceVtxA[i]);
    PDM_g_num_t *faceLnToGn = malloc (sizeof(PDM_g_num_t) * nFacePart);

    int nFacePartStored = 0;
    
    idx1 = 0;
    int *initToOlTmp1 = malloc (sizeof(int) * nFacePart);
    int *initToOlTmp2 = malloc (sizeof(int) * nFacePart);
    int *initToOlTmpN = malloc (sizeof(int) * nFace);

    int itmp1 = 0;
    for (int j = 0; j < nFace; j++) {
      initToOlTmpN[j] = 0;
      if (firstRecvStrideA[i][j] == 0) {
        initToOlTmpN[j] = 1;
        faceVtxIdxPart[nFacePartStored+1] = 
                faceVtxIdxPart[nFacePartStored] + 
                (iniFaceVtxIdxA[j+1] - iniFaceVtxIdxA[j]); 
        for (int k = iniFaceVtxIdxA[j]; k < iniFaceVtxIdxA[j+1]; k++) {
          faceVtxPart[idx1++] = iniFaceVtxA[k]; 
        }
        initToOlTmp1[nFacePartStored] = nFacePartStored; 
        faceLnToGn[nFacePartStored] = nTSubFacesA 
                                    + begUnChangedFaceA 
                                    + idxIpart
                                    + 1;
        idxIpart += 1;
        nFacePartStored += 1;  
      }
    }
    
    const PDM_g_num_t *iniGNumVtxA = PDM_surf_mesh_part_vtx_g_num_get (ol->meshA, i);
    const int ini_n_vtx = PDM_surf_mesh_part_n_vtx_get (ol->meshA, i);
    
    PDM_g_num_t *cpIniGNumVtxA = malloc (sizeof(PDM_g_num_t) * ini_n_vtx);
    int *orderVtxA = malloc (sizeof(int) * ini_n_vtx);
    for (int j = 0; j < ini_n_vtx; j++) {
      cpIniGNumVtxA[j] = iniGNumVtxA[j];
      orderVtxA[j] = j;
    }
    
    PDM_sort_long (cpIniGNumVtxA, orderVtxA, ini_n_vtx);
    
    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    int itmp2 = 0;
    for (int j = 0; j < nFace; j++) {
      if (firstRecvStrideA[i][j] > 0) {
        initToOlTmpN[j] = firstRecvStrideA[i][j];
        for (int k = 0; k < firstRecvStrideA[i][j]; k++) {
          int nVtx    = (int) firstRecvA[i][idx++];
          PDM_g_num_t numAbs  = firstRecvA[i][idx++];
          idx++;
          int iProcB  = (int) firstRecvA[i][idx++];
          int iPartB  = (int) firstRecvA[i][idx++];
          faceVtxIdxPart[nFacePartStored+1] = 
                  faceVtxIdxPart[nFacePartStored] + nVtx;
          
          faceLnToGn[nFacePartStored] = numAbs;
          
          olp->linkedFaces[4*itmp2]   = nFacePartStored + 1; 
          olp->linkedFaces[4*itmp2+1] = iProcB;
          olp->linkedFaces[4*itmp2+2] = iPartB;
          olp->linkedFaces[4*itmp2+3] = -1; // Come from B in the next step
          
          initToOlTmp2[itmp2++] = nFacePartStored;  
          nFacePartStored += 1;  
          
          for (int k1 = 0; k1 < nVtx; k1++) {
            PDM_g_num_t currentVtx = secondRecvA[i][idx2++];
            if (currentVtx <= iniNGVtx) {
              int idxVtx = PDM_binary_search_long (currentVtx, 
                                                   cpIniGNumVtxA, 
                                                   ini_n_vtx);
              assert (idxVtx != -1);
              faceVtxPart[idx1++] = orderVtxA[idxVtx];
            }
            else {
              gFaceVtxPart[idx3++] = currentVtx;
              faceVtxPart[idx1++] = -idx3;
            }

          }
        }
      }
    }
    
    for (int j = 0; j < nFace; j++) {
      olp->initToOlFaceIdx[j+1] = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      initToOlTmpN[0]; //FIXME: Problem with initToOlTmpN[0] initialization 
    }
    
    itmp1 = 0;
    itmp2 = 0;
    for (int j = 0; j < nFace; j++) {
      int ii = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      if (firstRecvStrideA[i][j] == 0) {
        olp->initToOlFace[ii] = initToOlTmp1[itmp1++];
      }
      else {
        for (int k = 0; k < firstRecvStrideA[i][j]; k++) {
          olp->initToOlFace[ii++] = initToOlTmp2[itmp2++];
        }
      }
      initToOlTmpN[j]++;
    }
    
    free (initToOlTmp1);
    free (initToOlTmp2);
    free (initToOlTmpN);
    
    PDM_g_num_t *cpgFaceVtxPart = malloc (sizeof(PDM_g_num_t) * idx3);
    int *ordergFaceVtxPart = malloc (sizeof(int) * idx3);
    for (int j = 0; j < idx3; j++) {
      cpgFaceVtxPart[j] = gFaceVtxPart[j];
      ordergFaceVtxPart [j] = j;
    }
    
    PDM_sort_long (cpgFaceVtxPart, ordergFaceVtxPart, idx3);
    
    int k1 = 0;
    int k2 = 0;
    while (true) {
      PDM_g_num_t val = cpgFaceVtxPart[k1];
      cpgFaceVtxPart[k2] = cpgFaceVtxPart[k1];
      k2 += 1;
      while (cpgFaceVtxPart[k1] == val) {
        k1 += 1;
        if (k1 >= idx3) {
          break;
        }
      }
      if (k1 >= idx3) {
        break;
      }
    }
    
    int iniNVtx = PDM_surf_mesh_part_n_vtx_get (ol->meshA, i);
    for (int j = 0; j < s_olFaceVtxA[i]; j++) {
      int iniIdx = faceVtxPart[j];
      if (iniIdx < 0) {
        PDM_g_num_t iniVal = gFaceVtxPart[-(iniIdx+1)];
        int idxVtx = PDM_binary_search_long (iniVal, cpgFaceVtxPart, k2);
        assert(idxVtx != -1);
        faceVtxPart[j] = iniNVtx + idxVtx + 1;
      }
    }
    
    free (ordergFaceVtxPart);
    free (cpgFaceVtxPart);

    int nVtxPart = iniNVtx + k2;
    
    double *coordsPart = malloc (sizeof(double) * 3 * nVtxPart);
    PDM_g_num_t *vtxLnToGnPart = malloc (sizeof(PDM_g_num_t) * nVtxPart);

    for (int j = 0; j < iniNVtx; j++) {
      vtxLnToGnPart[j] = iniGNumVtxA[j];
    }

    for (int j = 0; j < k2; j++) {
      vtxLnToGnPart[iniNVtx + j] = cpgFaceVtxPart[j];
    }

    PDM_g_num_t *cpVtxLnToGnPart = malloc (sizeof(PDM_g_num_t) * nVtxPart);
    int *orderVtxLnToGnPart = malloc (sizeof(int) * nVtxPart);
    for (int j = 0; j < nVtxPart; j++) {
      cpVtxLnToGnPart[j] = vtxLnToGnPart[j];
      orderVtxLnToGnPart[j] = j;
    }
    
    PDM_sort_long (cpVtxLnToGnPart, orderVtxLnToGnPart, nVtxPart);
    
    const double *iniVtx = PDM_surf_mesh_part_vtx_get (ol->meshA, i);
    
    for (int j = 0; j < iniNVtx; j++) {
      for (int k = 0; k < 3; k++) {
        coordsPart[3*j+k] = iniVtx[3*j+k];
      }
    }
    
    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    for (int j = 0; j < nFace; j++) {
      if (firstRecvStrideA[i][j] > 0) {
        for (int k = 0; k < firstRecvStrideA[i][j]; k++) {
          PDM_g_num_t nVtx    = firstRecvA[i][idx++];
          idx += 4;
          
          for (int k4 = 0; k4 < nVtx; k4++) {
            PDM_g_num_t currentVtx = secondRecvA[i][idx2++];
            int idxVtx = PDM_binary_search_long (currentVtx, 
                                                 cpVtxLnToGnPart, 
                                                 nVtxPart);

            coordsPart[3*orderVtxLnToGnPart[idxVtx]    ] = thirdRecvA[i][idx3++];
            coordsPart[3*orderVtxLnToGnPart[idxVtx] + 1] = thirdRecvA[i][idx3++];
            coordsPart[3*orderVtxLnToGnPart[idxVtx] + 2] = thirdRecvA[i][idx3++];
            
          }
        }
      }
    }

    olp->part = PDM_surf_part_create (nFacePart,
                                      faceVtxIdxPart,
                                      faceVtxPart,
                                      faceLnToGn,
                                      nVtxPart,
                                      coordsPart,
                                      vtxLnToGnPart);
    
    olp->faceIniVtxIdx = malloc (sizeof(int) * (nFace + 1));
    
    olp->faceIniVtxIdx[0] = 0;
    for (int j = 0; j < nFace; j++) {
      if (fourthRecvStrideA[i][j] > 0) {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] + fourthRecvStrideA[i][j];
      }
      else {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] + 
                iniFaceVtxIdxA[j+1] - iniFaceVtxIdxA[j];
                
      }
    }
    
    olp->faceIniVtx    = malloc (sizeof(int) * olp->faceIniVtxIdx[nFace]);

    idx = 0;
    idx2 = 0;
    for (int j = 0; j < nFace; j++) {
      if (fourthRecvStrideA[i][j] > 0) {
        for (int k = 0; k < fourthRecvStrideA[i][j]; k++) {
          PDM_g_num_t currentVtx = fourthRecvStrideA[i][idx2++];
          int idxVtx = PDM_binary_search_long (currentVtx, 
                                               cpVtxLnToGnPart, 
                                               nVtxPart);
          
          olp->faceIniVtx[idx++] = orderVtxLnToGnPart[idxVtx];
        }
      }
      else {
        for (int k = iniFaceVtxIdxA[j]; k < iniFaceVtxIdxA[j+1]; k++) {
          olp->faceIniVtx[idx++] = iniFaceVtxA[k]; 

        }
      }
    }
    
    free (orderVtxLnToGnPart);
    free (cpVtxLnToGnPart);
  }
 
  free (firstRecvA);
  free (firstRecvStrideA);

  free (secondRecvA);
  free (secondRecvStrideA);

  free (thirdRecvA);
  free (thirdRecvStrideA);
  
  free (fourthRecvA);
  free (fourthRecvStrideA);

  free (nUnChangedFacePartA);
  free (nSubFacePartA);
  free (s_olFaceVtxA);
    
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_LOCAL_CONNECTA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_LOCAL_CONNECTA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*****************************************************************************
   *                                                                           *
   * Compute B local numbering for vertices and faces                          * 
   *                                                                           *
   ****************************************************************************/

  PDM_g_num_t nUnChangedFaceB = 0;
  PDM_g_num_t nSubFaceB = 0;
  int *nUnChangedFacePartB = malloc(sizeof(int) * nPartB);
  int *nSubFacePartB = malloc(sizeof(int) * nPartB);
  int *s_olFaceVtxB = malloc(sizeof(int) * nPartB);

  /* 
   * Compute dimensions
   */

  idx = 0;
  for (int i = 0; i < nPartB; i++) {
    nUnChangedFacePartB[i] = 0;
    nSubFacePartB[i] = 0;
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    const PDM_g_num_t *gNumFaceB = PDM_surf_mesh_part_face_g_num_get (ol->meshB, i);
    const int *iniFaceVtxIdxB =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshB, i);

    idx = 0;
    for (int j = 0; j < nFace; j++) {
      if (firstRecvStrideB[i][j] == 0) {
        nUnChangedFacePartB[i] += 1;
        s_olFaceVtxB[i] += iniFaceVtxIdxB[j+1] - iniFaceVtxIdxB[j];
      }
      else {
        nSubFacePartB[i] += firstRecvStrideB[i][j]/5;
        for (int k = 0; k < firstRecvStrideB[i][j]; k++) {
          int nVtx    = (int) firstRecvB[i][idx++];
          idx++;
          PDM_g_num_t oNumAbs = firstRecvB[i][idx++];
          idx++;
          idx++;

          s_olFaceVtxB[i] += nVtx;

          assert(oNumAbs == gNumFaceB[j]);
        }
      }
    }
    
    nUnChangedFaceB += nUnChangedFacePartB[i];
    nSubFaceB += nSubFacePartB[i];

  }
  
  PDM_g_num_t begUnChangedFaceB = 0;
  PDM_MPI_Scan (&nUnChangedFaceB, &begUnChangedFaceB,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  begUnChangedFaceB += -nUnChangedFaceB + 1;

  PDM_g_num_t nFaceB = nUnChangedFaceB + nSubFaceB;
  PDM_g_num_t nGFaceB;
  PDM_MPI_Allreduce(&nFaceB, &nGFaceB, 1, 
                PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm); 

  PDM_g_num_t nGVtxB = n_g_newVtxB;
  ol->olMeshB = _ol_mesh_create (nGFaceB, nGVtxB, nPartB);

  /* 
   * Build local connectivity 
   */
  
  idxIpart = 0;
  
  iniNGVtx = PDM_surf_mesh_n_g_vtx_get (ol->meshB);
  
  for (int i = 0; i < nPartB; i++) {
    
    _ol_part_t *olp = _ol_part_create (meshB->part[i]);
    ol->olMeshB->part[i] = olp;
     
    olp->sInitToOlFace   = 0;
    olp->initToOlFaceIdx = NULL;       
    olp->initToOlFace    = NULL;
    olp->nLinkedFace     = 0;
    olp->linkedFaces     = NULL;
    
    int nFace = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    const int *iniFaceVtxIdxB =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshB, i);
    const int *iniFaceVtxB =  PDM_surf_mesh_part_face_vtx_get (ol->meshB, i);

    int nFacePart = nUnChangedFacePartB[i] + nSubFacePartB[i];
    olp->sInitToOlFace = nFacePart;
    olp->nLinkedFace = nSubFacePartB[i];

    olp->initToOlFaceIdx = malloc (sizeof(int) * (nFace + 1));
    olp->initToOlFaceIdx[0] = 0;
    
    olp->initToOlFace = malloc (sizeof(int) * nFacePart);
    olp->linkedFaces =  malloc (sizeof(int) * 4 * nSubFacePartB[i]); 
    
    int *faceVtxIdxPart = malloc (sizeof(int) * (nFacePart + 1));
    faceVtxIdxPart[0] = 0;
    
    PDM_g_num_t *gFaceVtxPart = malloc (sizeof(PDM_g_num_t) * s_olFaceVtxB[i]);
    int *faceVtxPart = malloc(sizeof(int) * s_olFaceVtxB[i]);
    PDM_g_num_t *faceLnToGn = malloc (sizeof(PDM_g_num_t) * nFacePart);

    int nFacePartStored = 0;
    
    idx1 = 0;
    int *initToOlTmp1 = malloc (sizeof(int) * nFacePart);
    int *initToOlTmp2 = malloc (sizeof(int) * nFacePart);
    int *initToOlTmpN = malloc (sizeof(int) * nFace);

    int itmp1 = 0;
    for (int j = 0; j < nFace; j++) {
      initToOlTmpN[j] = 0;
      if (firstRecvStrideB[i][j] == 0) {
        initToOlTmpN[j] = 1;
        faceVtxIdxPart[nFacePartStored+1] = 
                faceVtxIdxPart[nFacePartStored] + 
                (iniFaceVtxIdxB[j+1] - iniFaceVtxIdxB[j]); 
        for (int k = iniFaceVtxIdxB[j]; k < iniFaceVtxIdxB[j+1]; k++) {
          faceVtxPart[idx1++] = iniFaceVtxB[k]; 
        }
        initToOlTmp1[nFacePartStored] = nFacePartStored; 
        faceLnToGn[nFacePartStored] = nTSubFacesB 
                                    + begUnChangedFaceB 
                                    + idxIpart
                                    + 1;
        idxIpart += 1;
        nFacePartStored += 1;  
      }
    }
    
    const PDM_g_num_t *iniGNumVtxB = PDM_surf_mesh_part_vtx_g_num_get (ol->meshB, i);
    const int ini_n_vtx = PDM_surf_mesh_part_n_vtx_get (ol->meshB, i);
    
    PDM_g_num_t *cpIniGNumVtxB = malloc (sizeof(PDM_g_num_t) * ini_n_vtx);
    int *orderVtxB = malloc (sizeof(int) * ini_n_vtx);
    for (int j = 0; j < ini_n_vtx; j++) {
      cpIniGNumVtxB[j] = iniGNumVtxB[j];
      orderVtxB[j] = j;
    }
    
    PDM_sort_long (cpIniGNumVtxB, orderVtxB, ini_n_vtx);
    
    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    int itmp2 = 0;
    for (int j = 0; j < nFace; j++) {
      if (firstRecvStrideB[i][j] > 0) {
        initToOlTmpN[j] = firstRecvStrideB[i][j];
        for (int k = 0; k < firstRecvStrideB[i][j]; k++) {
          int nVtx    = (int) firstRecvB[i][idx++];
          PDM_g_num_t numAbs  = firstRecvB[i][idx++];
          idx++;
          int iProcB  = (int) firstRecvB[i][idx++];
          int iPartB  = (int) firstRecvB[i][idx++];
          faceVtxIdxPart[nFacePartStored+1] = 
                  faceVtxIdxPart[nFacePartStored] + nVtx;
          
          faceLnToGn[nFacePartStored] = numAbs;
          
          olp->linkedFaces[4*itmp2]   = nFacePartStored + 1; 
          olp->linkedFaces[4*itmp2+1] = iProcB;
          olp->linkedFaces[4*itmp2+2] = iPartB;
          olp->linkedFaces[4*itmp2+3] = -1; // Come from A in the next step
          
          initToOlTmp2[itmp2++] = nFacePartStored;  
          nFacePartStored += 1;  
          
          for (int k1 = 0; k1 < nVtx; k1++) {
            PDM_g_num_t currentVtx = secondRecvB[i][idx2++];
            if (currentVtx <= iniNGVtx) {
              int idxVtx = PDM_binary_search_long (currentVtx, 
                                                   cpIniGNumVtxB, 
                                                   ini_n_vtx);
              assert (idxVtx != -1);
              faceVtxPart[idx1++] = orderVtxB[idxVtx];
            }
            else {
              gFaceVtxPart[idx3++] = currentVtx;
              faceVtxPart[idx1++] = -idx3;
            }

          }
        }
      }
    }
    
    for (int j = 0; j < nFace; j++) {
      olp->initToOlFaceIdx[j+1] = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      initToOlTmpN[0];
    }
    
    itmp1 = 0;
    itmp2 = 0;
    for (int j = 0; j < nFace; j++) {
      int ii = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      if (firstRecvStrideB[i][j] == 0) {
        olp->initToOlFace[ii] = initToOlTmp1[itmp1++];
      }
      else {
        for (int k = 0; k < firstRecvStrideB[i][j]; k++) {
          olp->initToOlFace[ii++] = initToOlTmp2[itmp2++];
        }
      }
      initToOlTmpN[j]++;
    }
    
    free (initToOlTmp1);
    free (initToOlTmp2);
    free (initToOlTmpN);
    
    PDM_g_num_t *cpgFaceVtxPart = malloc (sizeof(PDM_g_num_t) * idx3);
    int *ordergFaceVtxPart = malloc (sizeof(int) * idx3);
    for (int j = 0; j < idx3; j++) {
      cpgFaceVtxPart[j] = gFaceVtxPart[j];
      ordergFaceVtxPart [j] = j;
    }
    
    PDM_sort_long (cpgFaceVtxPart, ordergFaceVtxPart, idx3);
    
    int k1 = 0;
    int k2 = 0;
    while (true) {
      PDM_g_num_t val = cpgFaceVtxPart[k1];
      cpgFaceVtxPart[k2] = cpgFaceVtxPart[k1];
      k2 += 1;
      while (cpgFaceVtxPart[k1] == val) {
        k1 += 1;
        if (k1 >= idx3) {
          break;
        }
      }
      if (k1 >= idx3) {
        break;
      }
    }
    
    int iniNVtx = PDM_surf_mesh_part_n_vtx_get (ol->meshB, i);
    for (int j = 0; j < s_olFaceVtxB[i]; j++) {
      int iniIdx = faceVtxPart[j];
      if (iniIdx < 0) {
        PDM_g_num_t iniVal = gFaceVtxPart[-(iniIdx+1)];
        int idxVtx = PDM_binary_search_long (iniVal, cpgFaceVtxPart, k2);
        assert(idxVtx != -1);
        faceVtxPart[j] = iniNVtx + idxVtx + 1;
      }
    }
    
    free (ordergFaceVtxPart);
    free (cpgFaceVtxPart);

    int nVtxPart = iniNVtx + k2;
    
    double *coordsPart = malloc (sizeof(double) * 3 * nVtxPart);
    PDM_g_num_t *vtxLnToGnPart = malloc (sizeof(PDM_g_num_t) * nVtxPart);

    for (int j = 0; j < iniNVtx; j++) {
      vtxLnToGnPart[j] = iniGNumVtxB[j];
    }

    for (int j = 0; j < k2; j++) {
      vtxLnToGnPart[iniNVtx + j] = cpgFaceVtxPart[j];
    }

    PDM_g_num_t *cpVtxLnToGnPart = malloc (sizeof(PDM_g_num_t) * nVtxPart);
    int *orderVtxLnToGnPart = malloc (sizeof(int) * nVtxPart);
    for (int j = 0; j < nVtxPart; j++) {
      cpVtxLnToGnPart[j] = vtxLnToGnPart[j];
      orderVtxLnToGnPart[j] = j;
    }
    
    PDM_sort_long (cpVtxLnToGnPart, orderVtxLnToGnPart, nVtxPart);
    
    const double *iniVtx = PDM_surf_mesh_part_vtx_get (ol->meshB, i);

    for (int j = 0; j < iniNVtx; j++) {
      for (int k = 0; k < 3; k++) {
        coordsPart[3*j+k] = iniVtx[3*j+k];
      }
    }
    
    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    for (int j = 0; j < nFace; j++) {
      if (firstRecvStrideB[i][j] > 0) {
        for (int k = 0; k < firstRecvStrideB[i][j]; k++) {
          PDM_g_num_t nVtx    = firstRecvB[i][idx++];
          idx += 4;
          
          for (int k4 = 0; k4 < nVtx; k4++) {
            PDM_g_num_t currentVtx = secondRecvB[i][idx2++];
            int idxVtx = PDM_binary_search_long (currentVtx, 
                                                 cpVtxLnToGnPart, 
                                                 nVtxPart);

            coordsPart[3*orderVtxLnToGnPart[idxVtx]    ] = thirdRecvB[i][idx3++];
            coordsPart[3*orderVtxLnToGnPart[idxVtx] + 1] = thirdRecvB[i][idx3++];
            coordsPart[3*orderVtxLnToGnPart[idxVtx] + 2] = thirdRecvB[i][idx3++];
            
          }
        }
      }
    }
    
    olp->part = PDM_surf_part_create (nFacePart,
                                      faceVtxIdxPart,
                                      faceVtxPart,
                                      faceLnToGn,
                                      nVtxPart,
                                      coordsPart,
                                      vtxLnToGnPart);
    
        
    olp->faceIniVtxIdx = malloc (sizeof(int) * (nFace + 1));
    
    olp->faceIniVtxIdx[0] = 0;
    for (int j = 0; j < nFace; j++) {
      if (fourthRecvStrideB[i][j] > 0) {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] + fourthRecvStrideB[i][j];
      }
      else {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] + 
                iniFaceVtxIdxB[j+1] - iniFaceVtxIdxB[j];
                
      }
    }
    
    olp->faceIniVtx    = malloc (sizeof(int) * olp->faceIniVtxIdx[nFace]);

    idx = 0;
    idx2 = 0;
    for (int j = 0; j < nFace; j++) {
      if (fourthRecvStrideB[i][j] > 0) {
        for (int k = 0; k < fourthRecvStrideB[i][j]; k++) {
          PDM_g_num_t currentVtx = fourthRecvStrideB[i][idx2++];
          int idxVtx = PDM_binary_search_long (currentVtx, 
                                               cpVtxLnToGnPart, 
                                               nVtxPart);
          
          olp->faceIniVtx[idx++] = orderVtxLnToGnPart[idxVtx];
        }
      }
      else {
        for (int k = iniFaceVtxIdxB[j]; k < iniFaceVtxIdxB[j+1]; k++) {
          olp->faceIniVtx[idx++] = iniFaceVtxB[k]; 

        }
      }
    }
    
    free (orderVtxLnToGnPart);
    free (cpVtxLnToGnPart);

  }
 
  free (firstRecvB);
  free (firstRecvStrideB);

  free (secondRecvB);
  free (secondRecvStrideB);

  free (thirdRecvB);
  free (thirdRecvStrideB);
  
  free (nUnChangedFacePartB);
  free (nSubFacePartB);
  free (s_olFaceVtxB);

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTB] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTB]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_LOCAL_CONNECTB]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_LOCAL_CONNECTB]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  
  /*****************************************************************************
   *                                                                           *
   * Update boundary graph with local face number (link)                       * 
   *                                                                           *
   ****************************************************************************/ 
  
  int *sendIdx = malloc (sizeof(int)*lComm);
  int *sendN   = malloc (sizeof(int)*lComm);
  PDM_g_num_t *sendBuff = NULL;
  
  int *recvIdx = malloc (sizeof(int)*lComm);
  int *recvN   = malloc (sizeof(int)*lComm);
  PDM_g_num_t *recvBuff = NULL;
  
  for (int i = 0; i < lComm; i++) {
    sendN[i] = 0;
  }
  
  int n_t_send = 0;
  for (int i = 0; i < nPartA; i++) {
    
    _ol_part_t *olp = ol->olMeshA->part[i];
    n_t_send += olp->nLinkedFace;
    
    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iProc    = olp->linkedFaces[4*j + 1];
      sendN[iProc] += 1;
    }  

  }

  int n_t_recv = 0;
  for (int i = 0; i < nPartB; i++) {
    
    _ol_part_t *olp = ol->olMeshB->part[i];
    n_t_recv += olp->nLinkedFace;
  }

  sendIdx[0] = 0;
  for (int i = 1; i < lComm; i++) {
    sendIdx[i] = sendIdx[i-1] + sendN[i-1];
    sendN[i-1] = 0;    
  }
  sendN[lComm-1] = 0;

  sendBuff = malloc(sizeof(PDM_g_num_t) * n_t_send * 3);
  recvBuff = malloc(sizeof(PDM_g_num_t) * n_t_recv * 3);

  for (int i = 0; i < nPartA; i++) {
    
    _ol_part_t *olp = ol->olMeshA->part[i];
    
    const PDM_g_num_t *numAbs = PDM_surf_part_faceLnToGn_get (olp->part);
    
    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iFacLoc  = olp->linkedFaces[4*j    ];
      int iProc    = olp->linkedFaces[4*j + 1];
      int iPart    = olp->linkedFaces[4*j + 2];

      int id = sendIdx[iProc] + sendN[iProc]; 
      ++sendN[iProc];

      sendBuff[3*id    ] = iPart; 
      sendBuff[3*id + 1] = numAbs[iFacLoc-1];
      sendBuff[3*id + 2] = iFacLoc;
      
    }  

  }
  
  //TODO: Corriger bug Intialiser les recvN et recvIdx
  
  
  PDM_MPI_Alltoallv (sendBuff, sendN, sendIdx, PDM__PDM_MPI_G_NUM,
                 recvBuff,recvN, recvIdx, PDM__PDM_MPI_G_NUM, ol->comm);  
  
  PDM_g_num_t **faceLnToGnBSorted = malloc (sizeof(PDM_g_num_t *) * nPartB);
  int **faceLnToGnBOrder = malloc (sizeof(int *) * nPartB);
  int **faceToLinked = malloc (sizeof(int *) * nPartB);

  for (int i = 0; i < nPartB; i++) {

    _ol_part_t *olp = ol->olMeshB->part[i];
    PDM_surf_part_t *_part = olp->part; 

    faceLnToGnBSorted[i] = malloc (sizeof(PDM_g_num_t) * _part->nFace);
    faceLnToGnBOrder[i] = malloc (sizeof(int) * _part->nFace);

    for (int j = 0; j < _part->nFace; j++) {
      faceLnToGnBSorted[i][j] = _part->faceLnToGn[j];
      faceLnToGnBOrder[i][j]  = j;
      faceToLinked[i][j] = -1;
    }
    
    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iFacLoc  = olp->linkedFaces[4*j    ];

      faceToLinked[i][iFacLoc-1] = j;
    }  

    PDM_sort_long (faceLnToGnBSorted[i],
                   faceLnToGnBOrder[i],
                   _part->nFace);
    
  }
  
  int k = 0;
  for (int i = 0; i < n_t_recv; i++) {

    int iPart    = (int) recvBuff[3*i    ];
    PDM_g_num_t numAbs   = recvBuff[3*i + 1];
    int iFacDist = (int) recvBuff[3*i + 2];
    
    _ol_part_t *olp = ol->olMeshB->part[i];
    PDM_surf_part_t *_part = olp->part; 
    int id = PDM_binary_search_long (numAbs, 
                                     faceLnToGnBSorted[iPart], 
                                     _part->nFace);
   
    int iFacIdx = faceLnToGnBOrder[i][id];

    int iLinked = faceToLinked[i][iFacIdx];
    
    olp->linkedFaces[4*iLinked + 3] = iFacDist;
    
    recvBuff[k++] = iFacIdx + 1;
    
  }

  for (int i = 0; i < nPartB; i++) {
    free (faceLnToGnBSorted[i]);
    free (faceLnToGnBOrder[i]);    
    free (faceToLinked[i]);
  }

  free (faceLnToGnBSorted);
  free (faceLnToGnBOrder);
  free (faceToLinked);

  for (int i = 0; i < lComm; i++) {
    sendN[i]   = sendN[i]/3;
    sendIdx[i] = sendIdx[i]/3;
    recvN[i]   = recvN[i]/3;
    recvIdx[i] = recvIdx[i]/3;
  }
  
  PDM_MPI_Alltoallv (recvBuff, recvN, recvIdx, PDM__PDM_MPI_G_NUM,
                 sendBuff, sendN, sendIdx, PDM__PDM_MPI_G_NUM,
                 ol->comm);  

  k = 0;
  for (int i = 0; i < nPartA; i++) {
    
    _ol_part_t *olp = ol->olMeshA->part[i];
    
    for (int j = 0; j < olp->nLinkedFace; j++) {
//FIXME: Verifier l'indice 3*k+2 qui etait a k++
      olp->linkedFaces[4*j + 3] = (int) sendBuff[3*k+2];
      
    }  

  }
  
  free (sendN);
  free (sendIdx);
  free (sendBuff);

  free (recvN);
  free (recvIdx);
  free (recvBuff);
  
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_UPDATE_A_B_CONNECT_GRAPH] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_UPDATE_A_B_CONNECT_GRAPH]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_UPDATE_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_UPDATE_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*****************************************************************************
   *                                                                           *
   *                     Sort boundary graph                                   * 
   *                                                                           *
   ****************************************************************************/

  for (int i = 0; i < nPartA; i++) {

    _ol_part_t *olp = ol->olMeshA->part[i];
    
    int *sortGraph = malloc (sizeof(int) * olp->nLinkedFace);
    int *orderGraph = malloc (sizeof(int) * olp->nLinkedFace);
    
    for (int j = 0; j < olp->nLinkedFace; j++) {
      orderGraph[j] = j;
      sortGraph[j]  = olp->linkedFaces[4*j + 1];
    }  

    PDM_sort_int (sortGraph, orderGraph, olp->nLinkedFace);
    
    int *tmpLinkedFaces = malloc (sizeof(int) * 4 * olp->nLinkedFace);
        
    for (int j = 0; j < olp->nLinkedFace; j++) {
      orderGraph[j] = j;
      sortGraph[j]  = olp->linkedFaces[4*j + 1];
      
      int newId = orderGraph[j];
      tmpLinkedFaces[4*j    ] = olp->linkedFaces[4*newId    ];
      tmpLinkedFaces[4*j + 1] = olp->linkedFaces[4*newId + 1];
      tmpLinkedFaces[4*j + 2] = olp->linkedFaces[4*newId + 2];
      tmpLinkedFaces[4*j + 3] = olp->linkedFaces[4*newId + 3];
      
    }  
    
    free (olp->linkedFaces);
    olp->linkedFaces = tmpLinkedFaces; 

    free (sortGraph);
    free (orderGraph);

  }

  for (int i = 0; i < nPartB; i++) {

    _ol_part_t *olp = ol->olMeshB->part[i];
    
    int *sortGraph = malloc (sizeof(int) * olp->nLinkedFace);
    int *orderGraph = malloc (sizeof(int) * olp->nLinkedFace);
    
    for (int j = 0; j < olp->nLinkedFace; j++) {
      orderGraph[j] = j;
      sortGraph[j]  = olp->linkedFaces[4*j + 1];
    }  

    PDM_sort_int (sortGraph, orderGraph, olp->nLinkedFace);
    
    int *tmpLinkedFaces = malloc (sizeof(int) * 4 * olp->nLinkedFace);
        
    for (int j = 0; j < olp->nLinkedFace; j++) {
      orderGraph[j] = j;
      sortGraph[j]  = olp->linkedFaces[4*j + 1];
      
      int newId = orderGraph[j];
      tmpLinkedFaces[4*j    ] = olp->linkedFaces[4*newId    ];
      tmpLinkedFaces[4*j + 1] = olp->linkedFaces[4*newId + 1];
      tmpLinkedFaces[4*j + 2] = olp->linkedFaces[4*newId + 2];
      tmpLinkedFaces[4*j + 3] = olp->linkedFaces[4*newId + 3];
      
    }  
    
    free (olp->linkedFaces);
    olp->linkedFaces = tmpLinkedFaces; 

    free (sortGraph);
    free (orderGraph);

  }

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_SORT_A_B_CONNECT_GRAPH] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_SORT_A_B_CONNECT_GRAPH]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_SORT_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_SORT_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  
  /*****************************************************************************
   *                                                                           *
   *                           Finish                                          * 
   *                                                                           *
   ****************************************************************************/

}


/**
 * \brief Compute overlay general surfaces
 *
 * This function computes overlay mesh of two surface meshes
 *
 * \param [in]  ol       overlay object                     
 *
 */

static void
_compute_overlay_surfaces
(
 PDM_ol_t *ol
)
{
  ol;
  PDM_error(__FILE__, __LINE__, 0, "Error _compute_overlay_surfaces : Not yet implemented\n");
  exit(0);
}


/**
 * \brief Compute overlay mesh
 *
 * This function computes overlayh mesh after plane surface checking
 *
 * \param [in]  ol       overlay object                     
 *
 */

static void
_compute_overlay
(
 PDM_ol_t *ol
)
{
  /* 
   * Compute faces extents
   */
  
  _compute_faceExtents(ol);
  
  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_FACE_EXTENTS] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_FACE_EXTENTS]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_FACE_EXTENTS]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_FACE_EXTENTS]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /* 
   * Check if initial meshes are plane sufaces ou general surfaces 
   */

  double dist;
  int isPlaneSurfaces = _is_same_plane (ol, &dist);

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[OL_IS_PLANE] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_IS_PLANE]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_IS_PLANE]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_IS_PLANE]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  double gMinCarLgthVtxA = PDM_surf_mesh_gMinCarLgthVtx_get (ol->meshA);
  double gMinCarLgthVtxB = PDM_surf_mesh_gMinCarLgthVtx_get (ol->meshB);

  double carLgth = _MIN (gMinCarLgthVtxA, gMinCarLgthVtxB);

  if (dist > carLgth) {
    PDM_error(__FILE__, __LINE__, 0, "Warning _compute_overlay : distance (%12.5e) between the two planes may"
                    " be too far to overlay its : \n"
                    " translate one of them and recompute\n", dist);
  }
  
  if (isPlaneSurfaces)
    _compute_overlay_planes(ol);
  else
    _compute_overlay_surfaces(ol);

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Build and initialize an overlaying object
 *
 * This function builds an initializes an overlaying surface meshes object 
 *
 * \param [in]  nPartMeshA   Number of local partitions of the meshA input
 * \param [in]  nGFaceA      Number of global faces of the meshA input
 * \param [in]  nGVtxA       Number of global vertices of the meshA input
 * \param [in]  nPartMeshB   Number of local partitions of the meshB input
 * \param [in]  nGFaceB      Number of global faces of the meshB input
 * \param [in]  nGVtxB       Number of global vertices of the meshB input
 * \param [in]  projectCoeff Projection coefficient to define the overlay surface projection
 *                           If value == 0, the surface projection is MeshA 
 *                           If value == 1, the surface projection is MeshB
 *                           If 0 < value < 1 , the projection surface is an intermediate surface
 * \param [in]  tolerance    Absolute tolerance used to define 
 *                           local geometric tolerance (0 < tolerance).
 * \param [in]  comm         MPI communicator.
 *
 * \return      id           Overlay object identifier.
 *
 */

int
PDM_ol_create
(
 const int         nPartMeshA,       
 const PDM_g_num_t  nGFaceMeshA,
 const PDM_g_num_t  nGVtxMeshA,
 const int         nPartMeshB,
 const PDM_g_num_t  nGFaceMeshB,
 const PDM_g_num_t  nGVtxMeshB,
 const double      projectCoeff,
 const PDM_MPI_Comm    comm
)
{
  /*
   * Update olArray size
   */

  if (olArray == NULL) {
    olArray = PDM_Handles_create (4);
  }

  /* 
   * Allocate structure
   */

  PDM_ol_t *ol = (PDM_ol_t *) malloc(sizeof(PDM_ol_t));

  int id = PDM_Handles_store (olArray, ol);
  
  /* 
   * Initialization
   */

  ol->timer                = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    ol->times_elapsed[i] = 0.;
    ol->times_cpu[i] = 0.;
    ol->times_cpu_u[i] = 0.;
    ol->times_cpu_s[i] = 0.;
  }

  PDM_timer_resume(ol->timer);
  
  ol->projectCoeff         = projectCoeff;  
  ol->vtxCarLengthTol      = 1e-3;
  ol->extentsTol           = 1e-3;
  ol->samePlaneTol         = 1e-6;  
  ol->comm                 = comm;
  ol->meshA                = PDM_surf_mesh_create(nGFaceMeshA, nGVtxMeshA, nPartMeshA, comm);
  ol->meshB                = PDM_surf_mesh_create(nGFaceMeshB, nGVtxMeshB, nPartMeshB, comm);
  ol->olMeshA              = NULL;
  ol->olMeshB              = NULL;
  ol->dbbtreeA             = NULL;

  PDM_timer_hang_on(ol->timer);

  return id;
}  


void
PROCF (pdm_ol_create, PDM_OL_CREATE) 
(
 int        *nPartMeshA,       
 PDM_g_num_t *nGFaceMeshA,
 PDM_g_num_t *nGVtxMeshA,
 int        *nPartMeshB,
 PDM_g_num_t *nGFaceMeshB,
 PDM_g_num_t *nGVtxMeshB,
 double     *projectCoeff,
 PDM_MPI_Fint   *comm,
 int        *id
)
{
  PDM_MPI_Comm comm_c = PDM_MPI_Comm_f2c(*comm);
  
  *id = PDM_ol_create (*nPartMeshA,       
                       *nGFaceMeshA,
                       *nGVtxMeshA,
                       *nPartMeshB,
                       *nGFaceMeshB,
                       *nGVtxMeshB,
                       *projectCoeff,
                       comm_c);
}


/**
 * \brief Set an overlay parameter
 *
 * This function sets en overlay parameter 
 *
 * \param [in]  id          PDM_ol identifier          
 * \param [in]  parameter   Parameter to define 
 * \param [in]  value       Parameter value 
 *
 */

void
PDM_ol_parameter_set
(
 const int                id,
 const PDM_ol_parameter_t parameter,
 const double             value
)
{
  PDM_ol_t *ol = _ol_get(id);
  PDM_timer_resume(ol->timer);

  switch (parameter) {

  case PDM_OL_CAR_LENGTH_TOL :
    ol->vtxCarLengthTol = value;
    break;

  case PDM_OL_EXTENTS_TOL :
    ol->extentsTol = value;
    break;

  case PDM_OL_SAME_PLANE_TOL :
    ol->samePlaneTol = value;
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_parameter_set :"
            " Unavailable parameter to set\n");
    abort();

  }
  PDM_timer_hang_on(ol->timer);
 
}


void
PROCF (pdm_ol_parameter_set, PDM_OL_PARAMETER_SET) 
(
 int                *id,
 PDM_ol_parameter_t *parameter,
 double             *value
)
{
  PDM_ol_parameter_set (*id,
                        *parameter,
                        *value);
}


/**
 * \brief Define input meshes properties
 *
 * This function defines the input meshes properties 
 *
 * \param [in]  id          ol identifier          
 * \param [in]  mesh        Input mesh to define 
 *                          (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  iPart       Partition to define  
 * \param [in]  nFace       Number of faces                     
 * \param [in]  faceVtxIdx  Index in the face -> vertex connectivity
 * \param [in]  faceVtxIdx  face -> vertex connectivity
 * \param [in]  faceLnToGn  Local face numbering to global face numbering 
 * \param [in]  nVtx        Number of vertices              
 * \param [in]  coords      Coordinates       
 * \param [in]  vtxLnToGn   Local vertex numbering to global vertex numbering 
 *
 */

void
PDM_ol_input_mesh_set
(
 const int            id,
 const PDM_ol_mesh_t  mesh,
 const int            iPart,
 const int            nFace,
 const int           *faceVtxIdx,
 const int           *faceVtx,
 const PDM_g_num_t    *faceLnToGn,
 const int            nVtx, 
 const double        *coords,
 const PDM_g_num_t    *vtxLnToGn
)
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 
  PDM_timer_resume(ol->timer);

  PDM_surf_mesh_t *_mesh = NULL;

  switch(mesh) {

  case PDM_OL_MESH_A:
    _mesh = ol->meshA;
    break;

  case PDM_OL_MESH_B:
    _mesh = ol->meshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_input_mesh_set :"
            " unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_surf_mesh_part_input (_mesh,
                            iPart,
                            nFace,
                            faceVtxIdx,
                            faceVtx,
                            faceLnToGn,
                            nVtx, 
                            coords,
                            vtxLnToGn);
  PDM_timer_hang_on(ol->timer);
}  

void
PROCF (pdm_ol_input_mesh_set, PDM_OL_INPUT_MESH_SET)
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 int           *ipart,
 int           *nFace,
 int           *faceVtxIdx,
 int           *faceVtx,
 PDM_g_num_t    *faceLnToGn,
 int           *nVtx, 
 double        *coords,
 PDM_g_num_t    *vtxLnToGn
)
{
  PDM_ol_input_mesh_set (*id,
                         *mesh,
                         *ipart,
                         *nFace,
                         faceVtxIdx,
                         faceVtx,
                         faceLnToGn,
                         *nVtx, 
                         coords,
                         vtxLnToGn);
}

    
/**
 * \brief Overlaying the input surface meshes
 *
 * This function overlays the input surface meshes 
 *
 * \param [in]  id       ol identifier          
 *
 */

void
PDM_ol_compute
(
 const int id
)
{

  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 

  PDM_timer_resume(ol->timer);

  PDM_timer_hang_on(ol->timer);  
  ol->times_elapsed[INIT_DEF_DATA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[INIT_DEF_DATA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[INIT_DEF_DATA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[INIT_DEF_DATA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /* 
   * First computation
   */

  if (ol->olMeshA == NULL || ol->olMeshB== NULL) {

    /* 
     * Computation of edges and inter partition exchange graph
     */
    
    _build_edges(ol);

    PDM_timer_hang_on(ol->timer);  
    ol->times_elapsed[INIT_BUILD_EDGES] = PDM_timer_elapsed(ol->timer);
    ol->times_cpu[INIT_BUILD_EDGES]     = PDM_timer_cpu(ol->timer);
    ol->times_cpu_u[INIT_BUILD_EDGES]   = PDM_timer_cpu_user(ol->timer);
    ol->times_cpu_s[INIT_BUILD_EDGES]   = PDM_timer_cpu_sys(ol->timer);
    PDM_timer_resume(ol->timer);

    /* 
     * Inter partition exchange graph (exchange to vertices)
     */

    _build_exchange_graph(ol);

    PDM_timer_hang_on(ol->timer);  
    ol->times_elapsed[INIT_BUILD_EXCH_GRAPH] = PDM_timer_elapsed(ol->timer);
    ol->times_cpu[INIT_BUILD_EXCH_GRAPH]     = PDM_timer_cpu(ol->timer);
    ol->times_cpu_u[INIT_BUILD_EXCH_GRAPH]   = PDM_timer_cpu_user(ol->timer);
    ol->times_cpu_s[INIT_BUILD_EXCH_GRAPH]   = PDM_timer_cpu_sys(ol->timer);
    PDM_timer_resume(ol->timer);

    /* 
     * Build ghost faces and edges (To compute caracteristic length 
     * and normal to vertices)
     */

    _compute_carLgthVtx(ol);

    PDM_timer_hang_on(ol->timer);  
    ol->times_elapsed[INIT_COMPUTE_CAR_LGTH] = PDM_timer_elapsed(ol->timer);
    ol->times_cpu[INIT_COMPUTE_CAR_LGTH]     = PDM_timer_cpu(ol->timer);
    ol->times_cpu_u[INIT_COMPUTE_CAR_LGTH]   = PDM_timer_cpu_user(ol->timer);
    ol->times_cpu_s[INIT_COMPUTE_CAR_LGTH]   = PDM_timer_cpu_sys(ol->timer);
    PDM_timer_resume(ol->timer);

  }
  
  else {
    
    ol->times_elapsed[INIT_BUILD_EDGES] = ol->times_elapsed[INIT_DEF_DATA];
    ol->times_cpu[INIT_BUILD_EDGES]     = ol->times_cpu[INIT_DEF_DATA];
    ol->times_cpu_u[INIT_BUILD_EDGES]   = ol->times_cpu_u[INIT_DEF_DATA];
    ol->times_cpu_s[INIT_BUILD_EDGES]   = ol->times_cpu_s[INIT_DEF_DATA];

    ol->times_elapsed[INIT_BUILD_EXCH_GRAPH] = ol->times_elapsed[INIT_DEF_DATA];
    ol->times_cpu[INIT_BUILD_EXCH_GRAPH]     = ol->times_cpu[INIT_DEF_DATA];
    ol->times_cpu_u[INIT_BUILD_EXCH_GRAPH]   = ol->times_cpu_u[INIT_DEF_DATA];
    ol->times_cpu_s[INIT_BUILD_EXCH_GRAPH]   = ol->times_cpu_s[INIT_DEF_DATA];

    ol->times_elapsed[INIT_COMPUTE_CAR_LGTH] = ol->times_elapsed[INIT_DEF_DATA];
    ol->times_cpu[INIT_COMPUTE_CAR_LGTH]     = ol->times_cpu[INIT_DEF_DATA];
    ol->times_cpu_u[INIT_COMPUTE_CAR_LGTH]   = ol->times_cpu_u[INIT_DEF_DATA];
    ol->times_cpu_s[INIT_COMPUTE_CAR_LGTH]   = ol->times_cpu_s[INIT_DEF_DATA];

  }

  /* 
   * Compute overlay
   */
  
  _compute_overlay(ol);

  PDM_timer_hang_on (ol->timer);

}  


void
PROCF (pdm_ol_compute, PDM_OL_COMPUTE) 
(
 int         *id
)
{
  PDM_ol_compute (*id);
}

 
/**
 * \brief Define the type of a mesh moving     
 *
 * This function defines the type of a mesh moving.
 * Only a mesh can move 
 *
 * \param [in]  id       PDM_ol identifier          
 * \param [in]  mesh     Moving mesh          
 * \param [in]  mv       Type of moving          
 *
 */

void
PDM_ol_moving_type_set
(
 const int           id,
 const PDM_ol_mesh_t mesh,
 const PDM_ol_mv_t   mv
)
{
  /* 
   * Get overlay object
   */
  
  PDM_ol_t *ol = _ol_get(id); 
  PDM_timer_resume(ol->timer);
  mesh;
  mv;
  PDM_error(__FILE__, __LINE__, 0, "PDM ERROR : PDM_ol_moving_type_set not implemented yet");
  exit(1);
  PDM_timer_hang_on(ol->timer);
}


void
PROCF (pdm_ol_moving_type_set, PDM_OL_MOVING_TYPE_SET) 
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 PDM_ol_mv_t   *mv
)
{
  PDM_ol_moving_type_set (*id,
                          *mesh,
                          *mv);
}

/**
 * \brief Define a translation
 *
 * This function defines a translation for the moving mesh 
 *
 * \param [in]  id       PDM_ol identifier          
 * \param [in]  vect     Translation vector          
 * \param [in]  center   Translation center
 *
 */

void
PDM_ol_translation_set
(
 const int       id,
 const double   *vect,
 const double   *center
 )
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 
  PDM_timer_resume(ol->timer);

  vect;
  center;
  PDM_error(__FILE__, __LINE__, 0, "PDM ERROR : PDM_ol_translation_set not implemented yet");
  exit(1);
  
  PDM_timer_hang_on(ol->timer);
}


void
PROCF (pdm_ol_translation_set, PDM_OL_TRANSLATION_SET)
(
 int          *id,
 double       *vect,
 double       *center
)
{
  PDM_ol_translation_set (*id,
                          vect,
                          center);
}

 

/**
 * \brief Define a rotation
 *
 * This function defines a rotation for the moving mesh 
 *
 * \param [in]  id        PDM_ol identifier          
 * \param [in]  direction Rotation direction         
 * \param [in]  center    Rotation center
 * \param [in]  angle      Rotation center (degrees)
 *
 */

void
PDM_ol_rotation_set
(
 const int      id,
 const double  *direction,
 const double  *center,
 const double   angle
)
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 

  PDM_timer_resume(ol->timer);

  direction;
  center;
  angle;
  PDM_error(__FILE__, __LINE__, 0, "PDM ERROR : PDM_ol_rotation_set not implemented yet");
  exit(EXIT_FAILURE);
  PDM_timer_hang_on(ol->timer);
}


void
PROCF (pdm_ol_rotation_set, PDM_OL_ROTATION_SET)
(
 int     *id,
 double  *direction,
 double  *center,
 double  *angle
)
{
  PDM_ol_rotation_set (*id,
                       direction,
                       center,
                       *angle);
}

 
/**
 * \brief Return the entitie sizes of the overlay mesh 
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id        PDM_ol identifier          
 * \param [in]  mesh      Input mesh
 *                        (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [out] nGOlFace  Global number of faces of the overlay mesh
 * \param [out] nGOlVtx   Global number of vertices of the overlay mesh
 *
 */

void
PDM_ol_mesh_dim_get
(
 const int            id,
 const PDM_ol_mesh_t  mesh,
       PDM_g_num_t    *nGOlFace,
       PDM_g_num_t    *nGOlVtx
)
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id);
  PDM_timer_resume(ol->timer);

  switch(mesh) {

  case PDM_OL_MESH_A:
    *nGOlFace = ol->olMeshA->nGFace;
    *nGOlVtx  = ol->olMeshA->nGVtx;  
    break;

  case PDM_OL_MESH_B:
    *nGOlFace = ol->olMeshB->nGFace;
    *nGOlVtx  = ol->olMeshB->nGVtx;  
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_timer_hang_on(ol->timer);
}  


void
PROCF (pdm_ol_mesh_dim_get, PDM_OL_MESH_DIM_GET) 
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 PDM_g_num_t    *nGOlFace,
 PDM_g_num_t    *nGOlVtx  
)
{
  PDM_ol_mesh_dim_get (*id,
                       *mesh,
                       nGOlFace,
                       nGOlVtx);  
}

 


/**
 * \brief Return the entitie sizes of the overlay mesh 
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id            PDM_ol identifier          
 * \param [in]  mesh          Input mesh
 *                            (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  ipart         Partition to define  
 * \param [out] nOlFace       Number of faces of the overlay mesh
 * \param [out] nOlLinkedFace Number of linked faces
 * \param [out] nOlVtx        Number of vertices of the overlay mesh
 * \param [out] sOlFaceVtx    Size of olFaceVtx for each partition 
 * \param [out] sInitToOlFace Size of initToOlFace for each partition
 *
 */

void
PDM_ol_part_mesh_dim_get
(
 const int            id,
 const PDM_ol_mesh_t  mesh,
 const int            ipart,
       int           *nOlFace,
       int           *nOlLinkedFace,
       int           *nOlVtx,  
       int           *sOlFaceVtx,
       int           *sInitToOlFace
)
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id);

  PDM_timer_resume(ol->timer);

  _ol_mesh_t *ol_mesh;
  
  switch(mesh) {

  case PDM_OL_MESH_A:
    ol_mesh = ol->olMeshA;
    break;

  case PDM_OL_MESH_B:
    ol_mesh = ol->olMeshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_surf_part_t *sp = ol_mesh->part[ipart]->part;
  
  *nOlFace       = sp->nFace;
  *nOlLinkedFace = ol_mesh->part[ipart]->nLinkedFace;
  *nOlVtx        = sp->nVtx;
  *sOlFaceVtx    = sp->sFaceVtx;
  *sInitToOlFace = ol_mesh->part[ipart]->sInitToOlFace;

  PDM_timer_hang_on(ol->timer);

}


void
PROCF (pdm_ol_part_mesh_dim_get, PDM_OL_PART_MESH_DIM_GET)
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 int           *ipart,
 int           *nOlFace,
 int           *nOlLinkedFace,       
 int           *nOlVtx,  
 int           *sOlFaceVtx,
 int           *sInitToOlFace
)
{
  PDM_ol_part_mesh_dim_get (*id,
                            *mesh,
                            *ipart,
                            nOlFace,
                            nOlLinkedFace,
                            nOlVtx,  
                            sOlFaceVtx,
                            sInitToOlFace);
}


/**
 * \brief Return the entitie of the overlay mesh 
 *
 * This function returns the entities of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id                   PDM_overlay identifier          
 * \param [in]  mesh                 Input mesh
 *                                   (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  iPart                Mesh partition identifier          
 * \param [out] olFaceVtxIdx    Array adress of \ref olFaceVtx index 
 *                                   (size : \ref nOlFace + 1) 
 * \param [out] olFaceVtx       Array adress of face vertex connectivity 
 *                                   (size : \ref sOlFaceVtx[\ref iPart])
 * \param [out] olLinkedFaceProcIdx olLinkedFace Index (size = nProc + 1)  
 * \param [out] olLinkedFace  Array adress of linked face in other mesh
 *                              For each face, 4 link properties :
 *                                    - local face number
 *                                    - connected process,
 *                                    - connected part number, 
 *                                    - connected local face number
 *                              (size : \ref 4 * nOlFace)
 * \param [out] olFaceLnToGn    Array adress of local to global face numbering 
 *                                   (size : \ref nOlFace)
 * \param [out] olCoords        Array adress of vertex coodinates 
 *                                   (size : 3 * \ref nOlVtx)
 * \param [out] olVtxLnToGn     Array adress of local to global vertex numbering array 
 *                                   (size : \ref nOlVtx)
 * \param [out] initToOlFaceIdx Array adress of \ref initToOlFace index 
 *                                   (size : \ref nOlVtx + 1)
 * \param [out] initToOlFace    Array adress of initial to ol faces
 * \param [out] initToOlVtx     Array adress of initial to ol vertices
 *
 */

void
PDM_ol_mesh_entities_get
(
 const int              id,
 const PDM_ol_mesh_t    mesh,
 const int              iPart,
 const int            **olFaceIniVtxIdx,
 const int            **olFaceIniVtx,
 const int            **olFaceVtxIdx,
 const int            **olFaceVtx,
 const int            **olLinkedFaceProcIdx,
 const int            **olLinkedFace,
 const PDM_g_num_t     **olFaceLnToGn,
 const double         **olCoords,
 const PDM_g_num_t     **olVtxLnToGn,
 const int            **initToOlFaceIdx,
 const int            **initToOlFace
)
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 

  PDM_timer_resume(ol->timer);
  
  _ol_mesh_t *ol_mesh;
  
  switch(mesh) {

  case PDM_OL_MESH_A:
    ol_mesh = ol->olMeshA;
    break;

  case PDM_OL_MESH_B:
    ol_mesh = ol->olMeshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_surf_part_t *sp = ol_mesh->part[iPart]->part;
  
  *olFaceIniVtxIdx      = ol_mesh->part[iPart]->faceIniVtxIdx;
  *olFaceIniVtx         = ol_mesh->part[iPart]->faceIniVtx;
  *olFaceVtxIdx         = sp->faceVtxIdx; 
  *olFaceVtx            = sp->faceVtx;
  *olLinkedFaceProcIdx  = ol_mesh->part[iPart]->linkedFacesProcIdx;
  *olLinkedFace         = ol_mesh->part[iPart]->linkedFaces;
  *olFaceLnToGn         = sp->faceLnToGn;
  *olCoords             = sp->coords;
  *olVtxLnToGn          = sp->vtxLnToGn;
  *initToOlFaceIdx      = ol_mesh->part[iPart]->initToOlFaceIdx;
  *initToOlFace         = ol_mesh->part[iPart]->initToOlFace;

  PDM_timer_hang_on(ol->timer);
}  

 
void
PROCF (pdm_ol_mesh_entities_get, PDM_OL_MESH_ENTITIES_GET)
(
 int            *id,
 PDM_ol_mesh_t  *mesh,
 int            *ipart,
 int            *olFaceIniVtxIdx,
 int            *olFaceIniVtx,
 int            *olFaceVtxIdx,
 int            *olFaceVtx,
 int            *olLinkedFaceProcIdx,
 int            *olLinkedFace,
 PDM_g_num_t     *olFaceLnToGn,
 double         *olCoords,
 PDM_g_num_t     *olVtxLnToGn,
 int            *initToOlFaceIdx,
 int            *initToOlFace
)
{
  const int          *_olFaceIniVtxIdx;
  const int          *_olFaceIniVtx;
  const int          *_olFaceVtxIdx;
  const int          *_olFaceVtx;
  const int          *_olFaceLinkFaceProcIdx;
  const int          *_olFaceLinkFace;
  const PDM_g_num_t   *_olFaceLnToGn;
  const double       *_olCoords;
  const PDM_g_num_t   *_olVtxLnToGn;
  const int          *_initToOlFaceIdx;
  const int          *_initToOlFace;

  int           nOlFace;
  int           nOlLinkedFace;
  int           nOlVtx;
  int           sOlFaceVtx;
  int           sInitToOlFace;

  PDM_ol_t *ol = _ol_get(*id);
  
  int lComm;
  PDM_MPI_Comm_size (ol->comm, &lComm);
  
  PDM_surf_mesh_t *_mesh;
 
  switch (*mesh) {

  case PDM_OL_MESH_A:
    _mesh = ol->meshA;
    break;

  case PDM_OL_MESH_B:
    _mesh = ol->meshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }
  
  int nFace_ini = PDM_surf_mesh_part_n_face_get (_mesh, *ipart);

  PDM_ol_part_mesh_dim_get (*id,
                            *mesh,
                            *ipart,
                            &nOlFace,
                            &nOlLinkedFace,
                            &nOlVtx,  
                            &sOlFaceVtx,
                            &sInitToOlFace);

  
  PDM_ol_mesh_entities_get (*id,
                            *mesh,
                            *ipart,
                            &_olFaceIniVtxIdx,
                            &_olFaceIniVtx,
                            &_olFaceVtxIdx,
                            &_olFaceVtx,
                            &_olFaceLinkFaceProcIdx,
                            &_olFaceLinkFace,
                            &_olFaceLnToGn,
                            &_olCoords,
                            &_olVtxLnToGn,
                            &_initToOlFaceIdx,
                            &_initToOlFace);

  for (int i = 0; i < nOlFace + 1; i++) {
    olFaceVtxIdx[i] = _olFaceVtxIdx[i];
  }

  for (int i = 0; i < _olFaceVtxIdx[nOlFace]; i++) {
    olFaceVtx[i] = _olFaceVtx[i];
  }
  
  for (int i = 0; i < nFace_ini + 1; i++) {
    olFaceIniVtxIdx[i] = _olFaceIniVtxIdx[i];
  }
  
  for (int i = 0; i < _olFaceIniVtxIdx[nFace_ini]; i++) {
    olFaceIniVtx[i] = _olFaceIniVtx[i];
  }

  for (int i = 0; i < lComm + 1; i++) {
    olLinkedFaceProcIdx[i] = _olFaceLinkFaceProcIdx[i];
  }

  for (int i = 0; i < 4 * olLinkedFaceProcIdx[lComm]; i++) {
    olLinkedFace[i] = _olFaceLinkFace[i];
  }
 
  for (int i = 0; i < nOlFace; i++) {
    olFaceLnToGn[i] = _olFaceLnToGn[i];
  }

  for (int i = 0; i < 3 * nOlVtx; i++) {
    olCoords[i] = _olCoords[i];
  }

  for (int i = 0; i < nOlVtx; i++) {
    olVtxLnToGn[i] = _olVtxLnToGn[i];
  }

  for (int i = 0; i < nFace_ini + 1; i++) {
    initToOlFaceIdx[i] = _initToOlFaceIdx[i];
  }
  
  for (int i = 0; i < _initToOlFaceIdx[nFace_ini]; i++) {
    initToOlFace[i] = _initToOlFace[i];
  }
  
}

 
/**
 * \brief Delete an overlay object
 *
 * This function deletes an overlay object 
 *
 * \param [in]  id                PDM_ol identifier.          
 *
 */

void
PDM_ol_del
(
 const int     id
)  
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 
  
  /* 
   * Delete
   */
  
  if (ol != NULL) {
    
    ol->meshA = PDM_surf_mesh_free(ol->meshA);
    ol->meshB = PDM_surf_mesh_free(ol->meshB);
    
    ol->olMeshA = _ol_mesh_free(ol->olMeshA);
    ol->olMeshB = _ol_mesh_free(ol->olMeshB);
    
    PDM_dbbtree_free (ol->dbbtreeA);
    
    /* 
     * Update storaged array
     */
    
    free(ol);
    
    PDM_Handles_handle_free (olArray, id, PDM_FALSE);
    
    const int n_ol = PDM_Handles_n_get (olArray);
    
    if (n_ol == 0) {
      olArray = PDM_Handles_free (olArray);
    }
  }
}  

void
PROCF (pdm_ol_del, PDM_OL_DEL)
(
 int     *id
)
{
  PDM_ol_del (*id);
}


/**
 * \brief Dump elapsed an CPU time
 *
 *
 * \param [in]  id                PDM_ol identifier.          
 *
 */

void
PDM_ol_dump_times 
(
 const int     id
)
{
  /* 
   * Get overlay object
   */

  PDM_ol_t *ol = _ol_get(id); 

  ol->times_cpu[END]     = ol->times_cpu[INIT_DEF_DATA];
  ol->times_cpu_u[END]   = ol->times_cpu_u[INIT_DEF_DATA];
  ol->times_cpu_s[END]   = ol->times_cpu_s[INIT_DEF_DATA];

  double t1;
  double t2;

  t1 = ol->times_elapsed[END] - ol->times_elapsed[BEGIN];
  t2 = ol->times_cpu[END] - ol->times_cpu[BEGIN];
  
  PDM_printf( "ol times ALL : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[INIT_DEF_DATA] - ol->times_elapsed[BEGIN];
  t2 = ol->times_cpu[INIT_DEF_DATA] - ol->times_cpu[BEGIN];
  
  PDM_printf( "ol times INIT_DEF_DATA : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[INIT_BUILD_EDGES] - ol->times_elapsed[INIT_DEF_DATA];
  t2 = ol->times_cpu[INIT_BUILD_EDGES] - ol->times_cpu[INIT_DEF_DATA];
  
  PDM_printf( "ol times INIT_BUILD_EDGES : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[INIT_BUILD_EXCH_GRAPH] - ol->times_elapsed[INIT_BUILD_EDGES];
  t2 = ol->times_cpu[INIT_BUILD_EXCH_GRAPH] - ol->times_cpu[INIT_BUILD_EDGES];
  
  PDM_printf( "ol times INIT_BUILD_EXCH_GRAPH : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[INIT_COMPUTE_CAR_LGTH] - ol->times_elapsed[INIT_BUILD_EXCH_GRAPH];
  t2 = ol->times_cpu[INIT_COMPUTE_CAR_LGTH] - ol->times_cpu[INIT_BUILD_EXCH_GRAPH];
  
  PDM_printf( "ol times INIT_COMPUTE_CAR_LGTH : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_FACE_EXTENTS] - ol->times_elapsed[INIT_COMPUTE_CAR_LGTH];
  t2 = ol->times_cpu[OL_FACE_EXTENTS] - ol->times_cpu[INIT_COMPUTE_CAR_LGTH];
  
  PDM_printf( "ol times OL_FACE_EXTENTS : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_IS_PLANE] - ol->times_elapsed[OL_FACE_EXTENTS];
  t2 = ol->times_cpu[OL_IS_PLANE] - ol->times_cpu[OL_FACE_EXTENTS];
  
  PDM_printf( "ol times OL_IS_PLANE : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_BUILD_DBBTREE] - ol->times_elapsed[OL_IS_PLANE];
  t2 = ol->times_cpu[OL_BUILD_DBBTREE] - ol->times_cpu[OL_IS_PLANE];
  
  PDM_printf( "ol times OL_BUILD_DBBTREE : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_BOXES_INTERSECT] - ol->times_elapsed[OL_BUILD_DBBTREE];
  t2 = ol->times_cpu[OL_BOXES_INTERSECT] - ol->times_cpu[OL_BUILD_DBBTREE];
  
  PDM_printf( "ol times OL_BOXES_INTERSECT : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_DISTRIB_BOXESA_BLOCK] - ol->times_elapsed[OL_BOXES_INTERSECT];
  t2 = ol->times_cpu[OL_DISTRIB_BOXESA_BLOCK] - ol->times_cpu[OL_BOXES_INTERSECT];
  
  PDM_printf( "ol times OL_DISTRIB_BOXESA_BLOCK : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_EDGES_INTERSECTION] - ol->times_elapsed[OL_DISTRIB_BOXESA_BLOCK];
  t2 = ol->times_cpu[OL_EDGES_INTERSECTION] - ol->times_cpu[OL_DISTRIB_BOXESA_BLOCK];
  
  PDM_printf( "ol times OL_EDGES_INTERSECTION : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_EDGES_SYNCHRO] - ol->times_elapsed[OL_EDGES_INTERSECTION];
  t2 = ol->times_cpu[OL_EDGES_SYNCHRO] - ol->times_cpu[OL_EDGES_INTERSECTION];
  
  PDM_printf( "ol times OL_EDGES_SYNCHRO : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_CLIPPING] - ol->times_elapsed[OL_EDGES_SYNCHRO];
  t2 = ol->times_cpu[OL_CLIPPING] - ol->times_cpu[OL_EDGES_SYNCHRO];
  
  PDM_printf( "ol times OL_CLIPPING : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESA] - ol->times_elapsed[OL_CLIPPING];
  t2 = ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESA] - ol->times_cpu[OL_CLIPPING];
  
  PDM_printf( "ol times OL_COMPUTE_ADD_SUB_FACESA : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_DISTRIB_BOXESB_BLOCK] - ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESA];
  t2 = ol->times_cpu[OL_DISTRIB_BOXESB_BLOCK] - ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESA];
  
  PDM_printf( "ol times OL_DISTRIB_BOXESB_BLOCK : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESB] - ol->times_elapsed[OL_DISTRIB_BOXESB_BLOCK];
  t2 = ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESB] - ol->times_cpu[OL_DISTRIB_BOXESB_BLOCK];
  
  PDM_printf( "ol times OL_COMPUTE_ADD_SUB_FACESB : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTA] - ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESB];
  t2 = ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTA] - ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESB];
  
  PDM_printf( "ol times OL_SEND_RESULTS_TO_INIT_PARTA : %12.5e %12.5e\n", t1, t2);
  
  t1 = ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTB] - ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTA];
  t2 = ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTB] - ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTA];
  
  PDM_printf( "ol times OL_SEND_RESULTS_TO_INIT_PARTB : %12.5e %12.5e\n", t1, t2);
  
  t1 = ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTA] - ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTB];
  t2 = ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTA] - ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTB];
  
  PDM_printf( "ol times OL_COMPUTE_LOCAL_CONNECTA : %12.5e %12.5e\n", t1, t2);
  
  t1 = ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTB] - ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTA];
  t2 = ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTB] - ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTA];
  
  PDM_printf( "ol times OL_COMPUTE_LOCAL_CONNECTB : %12.5e %12.5e\n", t1, t2);
  
  t1 = ol->times_elapsed[OL_UPDATE_A_B_CONNECT_GRAPH] - ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTB];
  t2 = ol->times_cpu[OL_UPDATE_A_B_CONNECT_GRAPH] - ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTB];
  
  PDM_printf( "ol times OL_UPDATE_A_B_CONNECT_GRAPH : %12.5e %12.5e\n", t1, t2);
  
  t1 = ol->times_elapsed[OL_SORT_A_B_CONNECT_GRAPH] - ol->times_elapsed[OL_UPDATE_A_B_CONNECT_GRAPH];
  t2 = ol->times_cpu[OL_SORT_A_B_CONNECT_GRAPH] - ol->times_cpu[OL_UPDATE_A_B_CONNECT_GRAPH];
  
  PDM_printf( "ol times OL_SORT_A_B_CONNECT_GRAPH : %12.5e %12.5e\n", t1, t2);
}



void
PROCF (pdm_ol_dump_times, PDM_OL_DUMP_TIMES)
(
 int     *id
)
{
  PDM_ol_dump_times (*id);
}


#undef _DOT_PRODUCT
#undef _MODULE
#undef _MIN
#undef _MAX
#undef _CROSS_PRODUCT_3D
 
#ifdef __cplusplus
}
#endif /* __cplusplus */


