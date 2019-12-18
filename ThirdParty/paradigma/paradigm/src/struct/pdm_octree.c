

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_mpi.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_box_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_octree.h"
#include "pdm_octree_seq.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _box_tree_stats_t
 * \brief  Statistic about bbtre
 *
 *  _box_tree_stats_t defines statistics about bounding box tree
 *
 */

typedef struct {

  int         dim;                     /*!< Layout dimension */

  /* The following fields have 3 global values:
     mean on ranks, minimum on ranks, and maximum on ranks */

  int         depth[3];                /*!< Tree depth */
  int         n_leaves[3];             /*!< Number of leaves */
  int         n_boxes[3];              /*!< Number of associated boxes */
  int         n_threshold_leaves[3];   /*!< Number of leaves over threshold */
  int         n_leaf_boxes[3];         /*!< Number of boxes per leaf */
  size_t      mem_used[3];             /*!< Memory used */
  size_t      mem_required[3];         /*!< Memory temporarily required */

} _box_tree_stats_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 *
 */

typedef struct  {
  int    octree_seq_id;             /*!< Identifier of the associated octree seq */

  PDM_MPI_Comm comm;             /*!< MPI communicator */
  PDM_MPI_Comm rankComm;             /*!< MPI communicator */

  int     maxBoxesLeafShared; /*!<  Max number of boxes in a leaf for coarse shared BBTree */

  int     maxTreeDepthShared; /*!< Max tree depth for coarse shared BBTree */

  float   maxBoxRatioShared;  /*!< Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
                                for coarse shared BBTree */

  PDM_box_set_t  *rankBoxes;  /*!< Rank Boxes */
  int             nUsedRank;  /*!< Number of used ranks */
  int            *usedRank;   /*!< used ranks */
  double         *usedRankExtents;   /*!< Extents of processes */

  PDM_box_tree_t *btShared;   /*!< Shared Boundary box tree */
  _box_tree_stats_t btsShared;/*!< Shared Boundary box tree statistic */

  int          n_point_cloud; /*!< Number of point cloud */
  int            *n_points;   /*!< Number of points */
  PDM_g_num_t   **g_num;      /*!< Point global number */

} _octree_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _octree_t *
_get_from_id
(
 int  id
)
{
  _octree_t *octree = (_octree_t *) PDM_Handles_get (_octrees, id);

  if (octree == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return octree;
}


/**
 * \brief  Initialize box_tree statistics
 *
 * \param [inout]  bts  pointer to box tree statistics structure
 *
 */

static void
_init_bt_statistics
(
_box_tree_stats_t  *bts
)
{
  size_t i;

  assert(bts != NULL);

  bts->dim = 0;

  for (i = 0; i < 3; i++) {
    bts->depth[i] = 0;
    bts->n_leaves[i] = 0;
    bts->n_boxes[i] = 0;
    bts->n_threshold_leaves[i] = 0;
    bts->n_leaf_boxes[i] = 0;
    bts->mem_used[i] = 0;
    bts->mem_required[i] = 0;
  }
}


/**
 * \brief Update box-tree statistics.
 *
 * For most fields, we replace previous values with the current ones.
 *
 * For memory required, we are interested in the maximum values over time
 * (i.e. algorthm steps); this is the case even for the minimal memory
 * required, we is thus the time maximum of the rank minimum.
 *
 * \param [inout]   bts   Pointer to box tree statistics structure
 * \param [inout]   bt    Pointer to box tree structure
 *
 */

static void
_update_bt_statistics
(
_box_tree_stats_t     *bts,
const PDM_box_tree_t  *bt
)
{
  int dim;
  size_t i;
  size_t mem_required[3];

  assert(bts != NULL);

  dim = PDM_box_tree_get_stats (bt,
                                bts->depth,
                                bts->n_leaves,
                                bts->n_boxes,
                                bts->n_threshold_leaves,
                                bts->n_leaf_boxes,
                                bts->mem_used,
                                mem_required);

  bts->dim = dim;

  for (i = 0; i < 3; i++)
    bts->mem_required[i] = PDM_MAX(bts->mem_required[i], mem_required[i]);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier
 */

int
PDM_octree_create
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const double tolerance,
 const PDM_MPI_Comm comm
)
{
  if (_octrees == NULL) {
    _octrees = PDM_Handles_create (4);
  }

  _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

  int id = PDM_Handles_store (_octrees, octree);

  octree->octree_seq_id = PDM_octree_seq_create (n_point_cloud, depth_max,
                                                 points_in_leaf_max, tolerance);
  octree->comm = comm;
  octree->rankComm = PDM_MPI_COMM_NULL;

  //octree->extents_proc = NULL;
  octree->n_point_cloud = n_point_cloud; /*!< Number of point cloud */

  octree->n_points = (int *) malloc (sizeof(int) * n_point_cloud);

  octree->g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_point_cloud);

  for (int i = 0; i < n_point_cloud; i++) {
    octree->n_points[i] = 0;
    octree->g_num[i] = NULL;
  }

  octree->rankBoxes = NULL;  /*!< Rank Boxes */
  octree->usedRank = NULL;  /*!< Rank Boxes */
  octree->nUsedRank = 0;  /*!< Rank Boxes */
  octree->btShared = NULL;   /*!< Shared Boundary box tree */

  octree->maxTreeDepthShared = 10;
  octree->maxBoxesLeafShared = 6;
  octree->maxBoxRatioShared = 5;

  _init_bt_statistics (&(octree->btsShared));

  return id;
}



/**
 *
 * \brief Create an octree structure from a sequential octree
 *
 * \param [in]   octree_seq_id      Sequential octree identifier
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier
 */

int
PDM_octree_from_octree_seq_create
(
const int octree_seq_id,
const PDM_MPI_Comm comm
)
{
  if (_octrees == NULL) {
    _octrees = PDM_Handles_create (4);
  }

  _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

  int id = PDM_Handles_store (_octrees, octree);

  octree->octree_seq_id = octree_seq_id;

  octree->comm = comm;

  //octree->extents_proc = NULL;

  return id;
}


//void
//PROCF (pdm_octree_create, PDM_OCTREE_CREATE)
//(
// const int *n_point_cloud,
// const int *depth_max,
// const int *points_in_leaf_max,
// const double *tolerance,
// const PDM_MPI_Fint *fcomm,
// const int *id
//);

/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_octree_free
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);

  //free (octree->extents_proc);

  free (octree->n_points);
  free (octree->g_num);
  free (octree->usedRank);
  free (octree->usedRankExtents);

  PDM_box_set_destroy(&(octree->rankBoxes));

  PDM_box_tree_destroy(&(octree->btShared));

  PDM_octree_seq_free (octree->octree_seq_id);

  if (octree->rankComm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Comm_free (&(octree->rankComm));
  }

  free (octree);

  PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

  const int n_octrees = PDM_Handles_n_get (_octrees);

  if (n_octrees == 0) {
    _octrees = PDM_Handles_free (_octrees);
  }

}

//void
//PROCF (pdm_octree_free, PDM_OCTREE_FREE)
//(
// const int          *id
//);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id                 Identifier
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 * \param [in]   g_num              Point global number or NULL
 *
 */


void
PDM_octree_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords,
 const PDM_g_num_t *g_num
)
{
  _octree_t *octree = _get_from_id (id);

  const int idebug = 0;

  if (idebug) {

    printf ("*** PDM_octree_point_cloud_set d step 1  \n");
    for (int i = 0; i < n_points; i++) {
      printf ("     %d (%12.5e %12.5e %12.5e) : \n", i,
              coords[3*i], coords[3*i+1], coords[3*i+2]);
    }
    printf ("*** PDM_octree_point_cloud_set f step 1 \n");
  }

  octree->n_points[i_point_cloud] = n_points;
  octree->g_num[i_point_cloud] = (PDM_g_num_t *) g_num;
  PDM_octree_seq_point_cloud_set (octree->octree_seq_id, i_point_cloud,
                                  n_points, coords);

}

//void
//PROCF (pdm_octree_point_cloud_set, PDM_OCTREE_POINT_CLOUD_SET)
//(
// const int          *id
// const int          *i_point_cloud,
// const int          *n_points,
// const double       *coords
//);


/**
 *
 * \brief Build octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_octree_build
(
 const int          id
)
{

  _octree_t *octree = _get_from_id (id);

  const int nInfoLocation = 3;
  const int sExtents = 3 * 2;

  int myRank;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (octree->comm, &lComm);

  PDM_octree_seq_build (octree->octree_seq_id);

  double * extents = PDM_octree_seq_extents_get (octree->octree_seq_id);

  int n_proc;
  PDM_MPI_Comm_size (octree->comm, &n_proc);

  double *extents_proc = malloc (sizeof(double) * n_proc * 6);

  PDM_MPI_Allgather (extents, 6, PDM_MPI_DOUBLE,
                     extents_proc, 6, PDM_MPI_DOUBLE,
                     octree->comm);

  int root_id = PDM_octree_seq_root_node_id_get (octree->octree_seq_id);

  int n_pts = PDM_octree_seq_n_points_get(octree->octree_seq_id, root_id);

  int *n_pts_proc = (int *) malloc (sizeof(int) * lComm);
  PDM_MPI_Allgather (&n_pts, 1, PDM_MPI_INT,
                     n_pts_proc, 1, PDM_MPI_INT,
                     octree->comm);

  int nUsedRank = 0;
  for (int i = 0; i < lComm; i++) {
    if (n_pts_proc[i] > 0) {
      nUsedRank += 1;
    }
  }

  int *numProc = (int *) malloc (sizeof(int *) * nUsedRank);

  octree->usedRank = numProc;
  octree->nUsedRank = nUsedRank;

  PDM_g_num_t *gNumProc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nUsedRank);

  int idx = 0;

  for (int i = 0; i < lComm; i++) {
    if (n_pts_proc[i] > 0) {
      gNumProc[idx] = idx+1;
      numProc[idx] = i;

      for (int j = 0; j < sExtents; j++) {
        extents_proc[idx*sExtents + j] = extents_proc[i*sExtents + j];
      }
      idx += 1;
    }
  }

  free (n_pts_proc);

  extents_proc = (double *) realloc (extents_proc,
                                   sizeof(double) * sExtents * nUsedRank);

  int *initLocationProc = (int *) malloc (sizeof(int) * nInfoLocation * nUsedRank);
  for (int i = 0; i < nInfoLocation * nUsedRank; i++) {
    initLocationProc[i] = 0;
  }

  //PDM_MPI_Comm rankComm;
  PDM_MPI_Comm_split(octree->comm, myRank, 0, &(octree->rankComm));

  octree->rankBoxes = PDM_box_set_create(3,
                                         1,
                                         0,
                                         nUsedRank,
                                         gNumProc,
                                         extents_proc,
                                         1,
                                         &nUsedRank,
                                         initLocationProc,
                                         octree->rankComm);

  octree->btShared = PDM_box_tree_create (octree->maxTreeDepthShared,
                                          octree->maxBoxesLeafShared,
                                          octree->maxBoxRatioShared);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (octree->btShared,
                          octree->rankBoxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(octree->btsShared), octree->btShared);

  free (gNumProc);
  free (initLocationProc);

  octree->usedRankExtents = extents_proc;

}

//void
//PROCF (pdm_octree_build, PDM_OCTREE_BUILD)
//(
// const int          *id
//);

/**
 *
 * \brief Get root node id
 *
 * \param [in]   id                 Identifier
 *
 * \return     Root node identifier (-1 if octree is not built)
 *
 */

int
PDM_octree_root_node_id_get
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_root_node_id_get (octree->octree_seq_id);

}

//void
//PROCF (pdm_octree_root_node_id_get, PDM_OCTREE_ROOT_NODE_ID_GET)
//(
// const int          *id,
// int                *root_node_id
//);


/**
 *
 * \brief Get ancestor node id
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 *
 * \return     Ancestor node identifier
 *
 */

int
PDM_octree_ancestor_node_id_get
(
 const int          id,
 const int          node_id
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_ancestor_node_id_get(octree->octree_seq_id, node_id);
}

//void
//PROCF (pdm_octree_ancestor_node_id_get, PDM_OCTREE_ANCESTOR_NODE_ID_GET)
//(
// const int          *id,
// const int          *node_id,
// int                *ancestor_node_id
//);


/**
 *
 * \brief Get node extents
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 *
 * \return     Extents
 *
 */

const double *
PDM_octree_node_extents_get
(
 const int          id,
 const int          node_id
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_node_extents_get (octree->octree_seq_id, node_id);
}


/**
 *
 * \brief Get children of a node
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 * \param [in]   child              Children
 *
 * \return     Children node id
 *
 */

int
PDM_octree_children_get
(
 const int                id,
 const int                node_id,
 const PDM_octree_child_t child
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_children_get (octree->octree_seq_id, node_id,
                                      (PDM_octree_seq_child_t) child);
}


/**
 *
 * \brief Get Neighbor of node
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 * \param [in]   direction          Neighbor direction
 *
 * \return     Neighbor node id (-1 if no neighbor)
 *
 */

int
PDM_octree_neighbor_get
(
 const int                    id,
 const int                    node_id,
 const PDM_octree_direction_t direction
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_neighbor_get (octree->octree_seq_id, node_id,
                                      (PDM_octree_seq_direction_t) direction);
}

/**
 *
 * \brief Get the number of point inside a node
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 *
 * \return   Number of points
 *
 */

int
PDM_octree_n_points_get
(
 const int                id,
 const int                node_id
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_n_points_get (octree->octree_seq_id, node_id);

}


/**
 *
 * \brief Get indexes of points inside a node
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 * \param [out]  point_clouds_id    Point clouds number
 *                                  (size = Number of points inside the node)
 * \param [out]  point_indexes      Point indexes
 *                                  (size = Number of points inside the node)
 *
 */

void
PDM_octree_points_get
(
 const int                id,
 const int                node_id,
 int                    **point_clouds_id,
 int                    **point_indexes
)
{
  _octree_t *octree = _get_from_id (id);

  PDM_octree_seq_points_get (octree->octree_seq_id, node_id,
                             point_clouds_id, point_indexes);
}


/**
 *
 * \brief Is it a leaf
 *
 * \param [in]   id                 Identifier
 * \param [in]   node_id            Node identifier
 *
 * \return   1 or 0
 *
 */

int
PDM_octree_leaf_is
(
 const int                id,
 const int                node_id
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_leaf_is (octree->octree_seq_id, node_id);
}


/**
 *
 * \brief Get extents
 *
 * \param [in]   id                 Identifier
 *
 * \return     Extents
 *
 */

double *
PDM_octree_extents_get
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);

  return PDM_octree_seq_extents_get (octree->octree_seq_id);

}


/**
 *
 * \brief Used processes extents
 *
 * \param [in]   id                 Identifier
 * \param [out]  used_ranks         Used ranks
 * \param [out]  extents            Used ranks extents
 *
 * \return Number of used ranks
 */

int
PDM_octree_processes_extents_get
(
 const int          id,
 int              *used_ranks[],
 double           *extents[]
)
{
  _octree_t *octree = _get_from_id (id);

  *extents = octree->usedRankExtents;
  *used_ranks = octree->usedRank;

  return octree->nUsedRank;

}


/**
 *
 * Look for closest points stored inside an octree
 *
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest point in octree global number
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *
 */

void
PDM_octree_closest_point
(
const int    id,
const int    n_pts,
double      *pts,
PDM_g_num_t *pts_g_num,
PDM_g_num_t *closest_octree_pt_g_num,
double      *closest_octree_pt_dist2
)
{
  const int idebug = 0;

  _octree_t *octree = _get_from_id (id);

  int myRank;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (octree->comm, &lComm);

  for (int i = 0; i < n_pts; i++) {
    closest_octree_pt_g_num[i] = -1;
    closest_octree_pt_dist2[i] = HUGE_VAL;
  }

  /*******************************
   *
   * Look for the closest process
   *
   *******************************/

  int *rank_id = (int *) malloc (sizeof(int) * n_pts);
  double *rank_min_max_dist = (double *) malloc (sizeof(double) * n_pts);

  PDM_box_tree_min_dist_max_box (octree->btShared,
                                 n_pts,
                                 pts,
                                 rank_id,
                                 rank_min_max_dist);

  for (int i = 0; i < n_pts; i++) {
    if (rank_id[i] >= 0) {
      rank_id[i] = octree->usedRank[rank_id[i]];
    }
  }

  if (idebug == 1) {
    printf ("*** PDM_octree_closest_point d step 1 :"
            " the closest process \n");
    for (int i = 0; i < n_pts; i++) {
      //      if (i == ptprint)
      printf ("     %d (%12.5e %12.5e %12.5e) : %d %12.5e\n", i,
              pts[3*i], pts[3*i+1], pts[3*i+2],
              rank_id[i] ,rank_min_max_dist[i]);
    }
    printf ("*** PDM_octree_closest_point f step 1 :"
            " the closest process \n");
  }

  /***********************************
   *
   * Send points to closest processes
   *
   ***********************************/

  int *n_send_pts = (int *) malloc (sizeof(int) * lComm);

  for (int i = 0; i < lComm; i++) {
    n_send_pts[i] = 0;
  }

  for (int i = 0; i < n_pts; i++) {
    n_send_pts[rank_id[i]]++;
  }

  int *n_recv_pts = (int *) malloc (sizeof(int) * lComm);

  PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                    n_recv_pts, 1, PDM_MPI_INT,
                    octree->comm);

  int *i_send_pts = (int *) malloc (sizeof(int) * (lComm + 1));
  i_send_pts[0] = 0;

  int *i_recv_pts = (int *) malloc (sizeof(int) * (lComm + 1));
  i_recv_pts[0] = 0;

  for (int i = 0; i < lComm; i++) {
    i_send_pts[i+1] =  i_send_pts[i] + n_send_pts[i];
    n_send_pts[i] = 0;

    i_recv_pts[i+1] =  i_recv_pts[i] + n_recv_pts[i];
  }

  double *send_pts = malloc(sizeof(double) * 3 * i_send_pts[lComm]);
  double *recv_pts = malloc(sizeof(double) * 3 * i_recv_pts[lComm]);

  for (int i = 0; i < n_pts; i++) {
    int irank = rank_id[i];
    int idx = 3*(i_send_pts[irank] + n_send_pts[irank]);
    n_send_pts[irank] += 1;

    for (int j = 0; j < 3; j++) {
      send_pts [idx + j] = pts[3*i+j];
    }
  }

  for (int i = 0; i < lComm; i++) {
    n_send_pts[i] *= 3;
    i_send_pts[i] *= 3;
    n_recv_pts[i] *= 3;
    i_recv_pts[i] *= 3;
  }

  PDM_MPI_Alltoallv (send_pts, n_send_pts, i_send_pts, PDM_MPI_DOUBLE,
                     recv_pts, n_recv_pts, i_recv_pts, PDM_MPI_DOUBLE,
                     octree->comm);

  free (rank_min_max_dist);

  for (int i = 0; i < lComm; i++) {
    n_send_pts[i] = n_send_pts[i]/3;
    i_send_pts[i] = i_send_pts[i]/3;
    n_recv_pts[i] = n_recv_pts[i]/3;
    i_recv_pts[i] = i_recv_pts[i]/3;
  }

  /***************************************************
   *
   *  Look for the closest point in closest processes
   *
   ***************************************************/

  int *closest_pt = (int *) malloc(sizeof(int) * 2 * i_recv_pts[lComm]);
  double *closest_dist = (double *) malloc(sizeof(double) * i_recv_pts[lComm]);

  PDM_octree_seq_closest_point (octree->octree_seq_id, i_recv_pts[lComm],
                                recv_pts, closest_pt, closest_dist);

  if (idebug == 1) {
    printf ("*** PDM_octree_closest_point d step 2 :"
            " closest point in the closest process \n");
    for (int i = 0; i <  i_recv_pts[lComm]; i++) {

      printf ("     %d (%12.5e %12.5e %12.5e) : %d %d %12.5e\n", i,
              recv_pts[3*i], recv_pts[3*i+1], recv_pts[3*i+2],
              closest_pt[2*i] , closest_pt[2*i+1] ,closest_dist[i]);
    }
    printf ("*** PDM_octree_closest_point f step 2 :"
            " closest point in the closest process \n");
  }

  free (closest_pt);
  free (recv_pts);

  /************************************************************
   *
   * Receive distance to closest points from closest processes
   *
   ************************************************************/

  double *recv_dist = send_pts;

  PDM_MPI_Alltoallv (closest_dist, n_recv_pts, i_recv_pts, PDM_MPI_DOUBLE,
                     recv_dist, n_send_pts, i_send_pts, PDM_MPI_DOUBLE,
                     octree->comm);

  free (closest_dist);

  double *upper_bound_dist = (double *) malloc (sizeof(double) * n_pts);

  for (int i = 0; i < lComm; i++) {
    n_send_pts[i] = 0;
  }

  for (int i = 0; i < n_pts; i++) {
    int irank = rank_id[i];
    int idx = i_send_pts[irank] + n_send_pts[irank];
    n_send_pts[irank] += 1;

    upper_bound_dist[i] = recv_dist[idx] + 1e-6 * recv_dist[idx];
  }
  if (idebug == 1) {
    printf ("*** PDM_octree_closest_point d step 3 :"
            " exchange closest distance \n");
    for (int i = 0; i < n_pts; i++) {
      //if (i == ptprint)
      printf ("     %d (%12.5e %12.5e %12.5e) : %12.5e\n", i,
              pts[3*i], pts[3*i+1], pts[3*i+2], upper_bound_dist[i]);
    }
    printf ("*** PDM_octree_closest_point f step 3 :"
            " exchange closest distance \n");
  }

  free (recv_dist);
  free (rank_id);

  /****************************************************************************
   *
   *  Send points to processes that distance are inferior to computed distance
   *
   *   Be careful with number of processes ! Make several send !
   *
   ****************************************************************************/

  int *i_boxes = NULL;
  int *boxes = NULL;

  PDM_box_tree_closest_upper_bound_dist_boxes_get (octree->btShared,
                                                   n_pts,
                                                   pts,
                                                   upper_bound_dist,
                                                   &i_boxes,
                                                   &boxes);


  for (int i = 0; i < lComm; i++) {
    n_send_pts[i] = 0;
    n_recv_pts[i] = 0;
  }

  for (int i = 0; i < lComm+1; i++) {
    i_send_pts[i] = 0;
    i_recv_pts[i] = 0;
  }

  for (int i = 0; i < i_boxes[n_pts]; i++) {
    boxes[i] = octree->usedRank[boxes[i]];
    n_send_pts[boxes[i]]++;
  }

  if (idebug == 1) {
    printf ("*** PDM_octree_closest_point d step 4 :"
            " proc to send \n");
    for (int i = 0; i < n_pts; i++) {
      //if (i == ptprint) {
      printf ("     %d (%12.5e %12.5e %12.5e) %12.5e %d : ", i,
              pts[3*i], pts[3*i+1], pts[3*i+2], upper_bound_dist[i],
              i_boxes[i+1] - i_boxes[i]);
      for (int j = i_boxes[i]; j < i_boxes[i+1]; j++) {
        printf(" %d", boxes[j]);
      }
      printf("\n");
      //      }
    }
    printf ("*** PDM_octree_closest_point f step 4 :"
            " proc to send \n");
  }


  free (upper_bound_dist);
  PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                    n_recv_pts, 1, PDM_MPI_INT,
                    octree->comm);

  PDM_g_num_t n_sendrecv[2] = {0, 0};
  for (int i = 0; i < lComm; i++) {
    n_sendrecv[0] += n_send_pts[i];
    n_sendrecv[1] += n_recv_pts[i];
  }

  PDM_g_num_t max_n_exch[2];
  PDM_MPI_Allreduce (&n_sendrecv, &max_n_exch, 2,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, octree->comm);

  PDM_g_num_t max_max_n_exch = PDM_MAX (max_n_exch[0], max_n_exch[1]);

  PDM_g_num_t sum_npts;
  PDM_g_num_t _n_pts = n_pts;

  PDM_MPI_Allreduce (&_n_pts, &sum_npts, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, octree->comm);

  const int factor = 10;
  int n_data_exch_max = (int) (sum_npts/ (PDM_g_num_t) lComm) * factor;

  //FIXME: reactiver multi echanges : correction des instabilites

  n_data_exch_max = max_max_n_exch;  /* Pour obliger a faire un echange */

  // Fin correction des instabilites

  int n_exch = (int)(max_max_n_exch / n_data_exch_max);
  if ((int)(max_max_n_exch % n_data_exch_max) > 0) {
    n_exch += 1;
  }

  double *data_send_pts1 = malloc (sizeof(double) * 3 * n_data_exch_max); // Optimiser la taille
  double *data_recv_pts1 = malloc (sizeof(double) * 3 * n_data_exch_max);

  int *n_send_pts1 = NULL;
  int *i_send_pts1 = NULL;

  PDM_g_num_t *__closest_octree_pt_g_num = NULL;
  double      *__closest_octree_pt_dist2 = NULL;

  if (n_exch > 1) {
    __closest_octree_pt_g_num = malloc (sizeof(PDM_g_num_t) * n_pts);
    __closest_octree_pt_dist2 = malloc (sizeof(double) * n_pts);
  }
  else {
    __closest_octree_pt_g_num = closest_octree_pt_g_num;
    __closest_octree_pt_dist2 = closest_octree_pt_dist2;
  }

  if (n_exch == 1) {
    n_send_pts1 = n_send_pts;
    i_send_pts1 = i_send_pts;
  }
  else {
    n_send_pts1 = (int *) malloc (sizeof(int) * lComm);
    for (int i = 0; i < lComm; i++) {
      n_send_pts1[i] = 0;
    }
    i_send_pts1 = (int *) malloc (sizeof(int) * (lComm + 1));
    i_send_pts1[0] = 0;
  }
  int *n_recv_pts1 = n_recv_pts;
  int *i_recv_pts1 = i_recv_pts;

  PDM_g_num_t *data_send_gnum1 = malloc (sizeof(PDM_g_num_t) * n_data_exch_max); // Optimiser la taille
  PDM_g_num_t *data_recv_gnum1 = malloc (sizeof(PDM_g_num_t) * n_data_exch_max);

  int *n_send_gnum1 = (int *) malloc (sizeof(int) * lComm);
  int *n_recv_gnum1 = (int *) malloc (sizeof(int) * lComm);

  int *i_send_gnum1 = (int *) malloc (sizeof(int) * (lComm+1));
  int *i_recv_gnum1 = (int *) malloc (sizeof(int) * (lComm+1));

  for (int i = 0; i < lComm; i++) {
    n_send_gnum1[i] = 0;
    n_recv_gnum1[i] = 0;
  }
  for (int i = 0; i < lComm+1; i++) {
    i_send_gnum1[i] = 0;
    i_recv_gnum1[i] = 0;
  }
  double *data_send_pts2 = NULL;
  double *data_recv_pts2 = NULL;
  int *n_send_pts2 = NULL;
  int *n_recv_pts2 = NULL;
  int *i_send_pts2 = NULL;
  int *i_recv_pts2 = NULL;
  PDM_g_num_t *data_send_gnum2 = NULL;
  PDM_g_num_t *data_recv_gnum2 = NULL;
  int *n_send_gnum2 = NULL;
  int *n_recv_gnum2 = NULL;
  int *i_send_gnum2 = NULL;
  int *i_recv_gnum2 = NULL;
  int *send_bounds1 = NULL;
  int *send_bounds2 = NULL;
  int *send_counts = NULL;

  int *stride_ptb = malloc (sizeof(int) * n_data_exch_max);

  // Remplissage premier buffer et envoi

  if (n_exch > 1) {
    data_send_pts2 = malloc (sizeof(double) * 3 * n_data_exch_max);
    data_recv_pts2 = malloc (sizeof(double) * 3 * n_data_exch_max);
    n_send_pts2 = (int *) malloc (sizeof(int) * lComm);
    n_recv_pts2 = (int *) malloc (sizeof(int) * lComm);
    i_send_pts2 = (int *) malloc (sizeof(int) * (lComm+1));
    i_recv_pts2 = (int *) malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_send_pts2[i] = 0;
      n_recv_pts2[i] = 0;
    }
    for (int i = 0; i < lComm+1; i++) {
      i_send_pts2[i] = 0;
      i_recv_pts2[i] = 0;
    }

    data_send_gnum2 = malloc (sizeof(PDM_g_num_t) * n_data_exch_max); // Optimiser la taille
    data_recv_gnum2 = malloc (sizeof(PDM_g_num_t) * n_data_exch_max);
    n_send_gnum2 = (int *) malloc (sizeof(int) * lComm);
    n_recv_gnum2 = (int *) malloc (sizeof(int) * lComm);
    i_send_gnum2 = (int *) malloc (sizeof(int) * (lComm+1));
    i_recv_gnum2 = (int *) malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_send_gnum2[i] = 0;
      n_recv_gnum2[i] = 0;
    }
    for (int i = 0; i < lComm+1; i++) {
      i_send_gnum2[i] = 0;
      i_recv_gnum2[i] = 0;
    }

    send_bounds1 = (int *) malloc (sizeof(int) * 2 * lComm);
    send_bounds2 = (int *) malloc (sizeof(int) * 2 * lComm);

    send_counts = (int *) malloc (sizeof(int) * 2 * lComm);
    for (int i = 0; i < 2 * lComm; i++) {
      send_bounds1[i] = 0;
      send_bounds2[i] = 0;
      send_counts[i] = 0;
    }

  }

  if (n_exch == 1) {

    i_send_pts1[0]  = 0;
    i_send_gnum1[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_send_pts1[i+1]   = i_send_pts1[i] + 3 * n_send_pts1[i];
      i_send_gnum1[i+1] = i_send_gnum1[i] + n_send_pts1[i];
      n_send_pts1[i]    = 0;
      n_send_gnum1[i]   = 0;

      i_recv_pts1[i+1]   = i_recv_pts1[i] + 3 * n_recv_pts1[i];
      i_recv_gnum1[i+1] = i_recv_gnum1[i] + n_recv_pts1[i];
      n_recv_gnum1[i]   = n_recv_pts1[i];
      n_recv_pts1[i]    = 3 * n_recv_pts1[i];
    }
    for (int i = 0; i < n_pts; i++) {
      for (int j = i_boxes[i]; j < i_boxes[i+1]; j++) {
        int iproc     = boxes[j];
        int idx_pts1  = i_send_pts1[iproc] + n_send_pts1[iproc];
        int idx_gnum1 = i_send_gnum1[iproc] + n_send_gnum1[iproc];
        for (int k = 0; k < 3; k++) {
          data_send_pts1[idx_pts1+k] = pts[3*i+k];
        }
        data_send_gnum1[idx_gnum1] = pts_g_num[i];
        n_send_pts1[iproc] += 3;
        n_send_gnum1[iproc] += 1;
      }
    }
  }

  else {
    i_send_pts1[0]  = 0;
    i_send_gnum1[0] = 0;
    for (int i = 0; i < lComm; i++) {
      n_send_pts1[i]      = 0;
      n_send_gnum1[i]     = 0;
      send_counts[i]     = 0;
      send_bounds1[2*i]    = 0;
      send_bounds1[2*i+1]  = PDM_MIN (n_send_pts[i] / n_exch, n_data_exch_max);
    }

    for (int i = 0; i < n_pts; i++) {
      for (int j = i_boxes[i]; j < i_boxes[i+1]; j++) {
        int iproc     = boxes[j];
        if ((send_counts[iproc] >= send_bounds1[2*iproc]) &&
            (send_counts[iproc] < send_bounds1[2*iproc+1])) {
          n_send_pts1[iproc] += 3;
          n_send_gnum1[iproc] += 1;
        }
        send_counts[iproc]++;
      }
    }

    for (int i = 0; i < lComm; i++) {
      i_send_pts1[i+1]  = i_send_pts1[i] + n_send_pts1[i];
      i_send_gnum1[i+1] = i_send_gnum1[i] + n_send_gnum1[i];
      send_counts[i] = 0;
      n_send_pts1[i] = 0;
      n_send_gnum1[i] = 0;
    }

    for (int i = 0; i < n_pts; i++) {
      for (int j = i_boxes[i]; j < i_boxes[i+1]; j++) {
        int iproc     = boxes[j];
        if ((send_counts[iproc] >= send_bounds1[2*iproc]) &&
            (send_counts[iproc] < send_bounds1[2*iproc+1])) {
          int idx_pts1  = i_send_pts1[iproc] + n_send_pts1[iproc];
          int idx_gnum1 = i_send_gnum1[iproc] + n_send_gnum1[iproc];
          for (int k = 0; k < 3; k++) {
            data_send_pts1[idx_pts1+k] = pts[3*i+k];
          }
          data_send_gnum1[idx_gnum1] = pts_g_num[i];
          n_send_pts1[iproc] += 3;
          n_send_gnum1[iproc] += 1;
        }
        send_counts[iproc]++;
      }
    }

    PDM_MPI_Alltoall (n_send_gnum1, 1, PDM_MPI_INT,
                      n_recv_gnum1, 1, PDM_MPI_INT,
                      octree->comm);

    i_recv_pts1[0] = 0;
    i_recv_gnum1[0] = 0;

    for (int i = 0; i < lComm; i++) {
      i_recv_pts1[i+1]  = i_recv_pts1[i] + 3 * n_recv_gnum1[i];
      i_recv_gnum1[i+1] = i_recv_gnum1[i] + n_recv_gnum1[i];
      n_recv_pts1[i] = 3 * n_recv_gnum1[i];
    }

  }

  double      *data_send_pts  = data_send_pts1;
  double      *data_recv_pts  = data_recv_pts1;
               i_send_pts     = i_send_pts1;
               n_send_pts     = n_send_pts1;
               i_recv_pts     = i_recv_pts1;
               n_recv_pts     = n_recv_pts1;

  PDM_g_num_t *data_send_gnum = data_send_gnum1;
  PDM_g_num_t *data_recv_gnum = data_recv_gnum1;
  int         *i_send_gnum     = i_send_gnum1;
  int         *n_send_gnum     = n_send_gnum1;
  int         *i_recv_gnum     = i_recv_gnum1;
  int         *n_recv_gnum     = n_recv_gnum1;

  int         *send_bounds     = send_bounds1;

  double      *data_send_pts_next  = NULL;
  double      *data_recv_pts_next  = NULL;
  int         *i_send_pts_next     = NULL;
  int         *n_send_pts_next     = NULL;
  int         *i_recv_pts_next     = NULL;
  int         *n_recv_pts_next     = NULL;

  PDM_g_num_t *data_send_gnum_next  = NULL;
  PDM_g_num_t *data_recv_gnum_next  = NULL;
  int         *i_send_gnum_next     = NULL;
  int         *n_send_gnum_next     = NULL;
  int         *i_recv_gnum_next     = NULL;
  int         *n_recv_gnum_next     = NULL;

  int         *send_bounds_next     = NULL;

  PDM_MPI_Request Request_coord[2];
  PDM_MPI_Request Request_gnum[2];

  /* printf ("n_send_pts : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_send_pts[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_send_pts : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_send_pts[i]); */
  /* } */
  /* printf ("\n"); */


  /* printf ("n_recv_pts : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_recv_pts[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_recv_pts : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_recv_pts[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("data_send_pts : "); */
  /* for (int i = 0; i < i_send_pts[lComm]; i++) { */
  /*   printf(" %12.5e", data_send_pts[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("n_send_gnum : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_send_gnum[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_send_gnum : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_send_gnum[i]); */
  /* } */
  /* printf ("\n"); */


  /* printf ("n_recv_gnum : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_recv_gnum[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_recv_gnum : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_recv_gnum[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("data_send_gnum : "); */
  /* for (int i = 0; i < i_send_gnum[lComm]; i++) { */
  /*   printf(" %ld", data_send_gnum[i]); */
  /* } */
  /* printf ("\n"); */



  PDM_MPI_Ialltoallv (data_send_pts, n_send_pts, i_send_pts, PDM_MPI_DOUBLE,
                      data_recv_pts, n_recv_pts, i_recv_pts, PDM_MPI_DOUBLE,
                      octree->comm, &(Request_coord[0]));

  PDM_MPI_Ialltoallv (data_send_gnum, n_send_gnum, i_send_gnum, PDM__PDM_MPI_G_NUM,
                      data_recv_gnum, n_recv_gnum, i_recv_gnum, PDM__PDM_MPI_G_NUM,
                      octree->comm, &(Request_gnum[0]));

  int *_closest_octree_pt_id         = NULL;
  double *_closest_octree_pt_dist2   = NULL;
  PDM_g_num_t *_closest_octree_pt_g_num = NULL;

  int s_closest_octree_pt_dist2 = 0;

  for (int i = 0; i < n_exch; i++) {

    const int icurr = i%2;
    const int inext = 1 - icurr;

    if (icurr == 0) {
      n_send_pts    = n_send_pts1;
      n_recv_pts    = n_recv_pts1;
      i_send_pts    = i_send_pts1;
      i_recv_pts    = i_recv_pts1;
      data_send_pts = data_send_pts1;
      data_recv_pts = data_recv_pts1;

      n_send_gnum    = n_send_gnum1;
      n_recv_gnum    = n_recv_gnum1;
      i_send_gnum    = i_send_gnum1;
      i_recv_gnum    = i_recv_gnum1;
      data_send_gnum = data_send_gnum1;
      data_recv_gnum = data_recv_gnum1;

      send_bounds = send_bounds1;

      n_send_pts_next    = n_send_pts2;
      n_recv_pts_next    = n_recv_pts2;
      i_send_pts_next    = i_send_pts2;
      i_recv_pts_next    = i_recv_pts2;
      data_send_pts_next = data_send_pts2;
      data_recv_pts_next = data_recv_pts2;

      n_send_gnum_next    = n_send_gnum2;
      n_recv_gnum_next    = n_recv_gnum2;
      i_send_gnum_next    = i_send_gnum2;
      i_recv_gnum_next    = i_recv_gnum2;
      data_send_gnum_next = data_send_gnum2;
      data_recv_gnum_next = data_recv_gnum2;

      send_bounds_next = send_bounds2;
    }
    else {
      n_send_pts    = n_send_pts2;
      n_recv_pts    = n_recv_pts2;
      i_send_pts    = i_send_pts2;
      i_recv_pts    = i_recv_pts2;
      data_send_pts = data_send_pts2;
      data_recv_pts = data_recv_pts2;

      n_send_gnum    = n_send_gnum2;
      n_recv_gnum    = n_recv_gnum2;
      i_send_gnum    = i_send_gnum2;
      i_recv_gnum    = i_recv_gnum2;
      data_send_gnum = data_send_gnum2;
      data_recv_gnum = data_recv_gnum2;

      send_bounds = send_bounds2;

      n_send_pts_next    = n_send_pts1;
      n_recv_pts_next    = n_recv_pts1;
      i_send_pts_next    = i_send_pts1;
      i_recv_pts_next    = i_recv_pts1;
      data_send_pts_next = data_send_pts1;
      data_recv_pts_next = data_recv_pts1;

      n_send_gnum_next    = n_send_gnum1;
      n_recv_gnum_next    = n_recv_gnum1;
      i_send_gnum_next    = i_send_gnum1;
      i_recv_gnum_next    = i_recv_gnum1;
      data_send_gnum_next = data_send_gnum1;
      data_recv_gnum_next = data_recv_gnum1;

      send_bounds_next = send_bounds1;
    }

    // Remplissage prochain buffer et envoi

    if (i < (n_exch - 1)) {

      i_send_pts_next[0]  = 0;
      i_send_gnum_next[0] = 0;
      for (int j = 0; j < lComm; j++) {
        n_send_pts_next[j]      = 0;
        n_send_gnum_next[j]     = 0;
        send_counts[j]          = 0;
        send_bounds_next[2*j]    = send_bounds[2*j+1];
        send_bounds_next[2*j+1]  = send_bounds_next[2*j] +
                                   PDM_MIN (n_send_pts[j] / n_exch, n_data_exch_max);
      }

      for (int k = 0; k < n_pts; k++) {
        for (int j = i_boxes[k]; j < i_boxes[k+1]; j++) {
          int iproc     = boxes[j];
          if ((send_counts[iproc] >= send_bounds_next[2*iproc]) &&
              (send_counts[iproc] < send_bounds_next[2*iproc+1])) {
            n_send_pts_next[iproc] += 3;
            n_send_gnum_next[iproc] += 1;
          }
          send_counts[iproc]++;
        }
      }

      for (int j = 0; j < lComm; j++) {
        i_send_pts_next[j+1]  = i_send_pts_next[j] + n_send_pts_next[j];
        i_send_gnum_next[j+1] = i_send_gnum_next[j] + n_send_gnum_next[j];
        send_counts[j] = 0;
        n_send_pts_next[j] = 0;
        n_send_gnum_next[j] = 0;
      }

      for (int k = 0; k < n_pts; k++) {
        for (int j = i_boxes[k]; j < i_boxes[k+1]; j++) {
          int iproc     = boxes[j];
          if ((send_counts[iproc] >= send_bounds_next[2*iproc]) &&
              (send_counts[iproc] < send_bounds_next[2*iproc+1])) {
            int idx_pts_next  = i_send_pts_next[iproc] + n_send_pts_next[iproc];
            int idx_gnum_next = i_send_gnum_next[iproc] + n_send_gnum_next[iproc];
            for (int k1 = 0; k1 < 3; k1++) {
              data_send_pts_next[idx_pts_next+k1] = pts[3*k+k1];
            }
            data_send_gnum_next[idx_gnum_next] = pts_g_num[k];
            n_send_pts_next[iproc] += 3;
            n_send_gnum_next[iproc] += 1;
          }
          send_counts[iproc]++;
        }
      }

      PDM_MPI_Alltoall (n_send_gnum_next, 1, PDM_MPI_INT,
                        n_recv_gnum_next, 1, PDM_MPI_INT,
                        octree->comm);

      i_recv_pts_next[0] = 0;
      i_recv_gnum_next[0] = 0;

      for (int j = 0; j < lComm; j++) {
        i_recv_pts_next[j+1]  = i_recv_pts_next[j] + 3 * n_recv_gnum_next[j];
        i_recv_gnum_next[j+1] = i_recv_gnum_next[j] + n_recv_gnum_next[j];
        n_recv_pts_next[j] = 3 * n_recv_gnum_next[j];
      }

  /* printf ("n_send_pts : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_send_pts_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_send_pts_next : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_send_pts_next[i]); */
  /* } */
  /* printf ("\n"); */


  /* printf ("n_recv_pts_next : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_recv_pts_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_recv_pts_next : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_recv_pts_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("data_send_pts_next : "); */
  /* for (int i = 0; i < i_send_pts_next[lComm]; i++) { */
  /*   printf(" %12.5e", data_send_pts_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("n_send_gnum_next : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_send_gnum_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_send_gnum_next : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_send_gnum_next[i]); */
  /* } */
  /* printf ("\n"); */


  /* printf ("n_recv_gnum_next : "); */
  /* for (int i = 0; i < lComm; i++) { */
  /*   printf(" %d", n_recv_gnum_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("i_recv_gnum_next : "); */
  /* for (int i = 0; i < lComm +1; i++) { */
  /*   printf(" %d", i_recv_gnum_next[i]); */
  /* } */
  /* printf ("\n"); */

  /* printf ("data_send_gnum_next : "); */
  /* for (int i = 0; i < i_send_gnum_next[lComm]; i++) { */
  /*   printf(" %ld", data_send_gnum_next[i]); */
  /* } */
  /* printf ("\n"); */
      PDM_MPI_Ialltoallv (data_send_pts_next, n_send_pts_next,
                          i_send_pts_next, PDM_MPI_DOUBLE,
                          data_recv_pts_next, n_recv_pts_next,
                          i_recv_pts_next, PDM_MPI_DOUBLE,
                          octree->comm, &(Request_coord[inext]));

      PDM_MPI_Ialltoallv (data_send_gnum_next, n_send_gnum_next,
                          i_send_gnum_next, PDM__PDM_MPI_G_NUM,
                          data_recv_gnum_next, n_recv_gnum_next,
                          i_recv_gnum_next, PDM__PDM_MPI_G_NUM,
                          octree->comm, &(Request_gnum[inext]));

    }

    PDM_MPI_Wait (&(Request_coord[icurr]));
    PDM_MPI_Wait (&(Request_gnum[icurr]));

    // Attente reception buffer courant

    /* Look for the closest points for received points */

    if (s_closest_octree_pt_dist2 < i_recv_gnum[lComm]) {
      s_closest_octree_pt_dist2 = i_recv_gnum[lComm];
      _closest_octree_pt_id =
        realloc (_closest_octree_pt_id,
                 sizeof(int) * 2 * s_closest_octree_pt_dist2);

      _closest_octree_pt_dist2 =
        realloc (_closest_octree_pt_dist2,
                 sizeof(double) * s_closest_octree_pt_dist2);

      _closest_octree_pt_g_num =
        realloc (_closest_octree_pt_g_num,
                 sizeof(PDM_g_num_t) * s_closest_octree_pt_dist2);
    }

    /* for (int j1 = 0; j1 < i_recv_gnum[lComm]; j1++) { */
    /* /\*   if (j1 == 13) { *\/ */
    /*     /\* for (int j = i_recv_gnum[j1]; j < i_recv_gnum[j1+1]; j++) { *\/ */
    /*       /\* if (((ptprintc[0] - data_recv_pts[3*j]) *(ptprintc[0] - data_recv_pts[3*j]) + *\/ */
    /*       /\*     (ptprintc[1] - data_recv_pts[3*j+1]) *(ptprintc[1] - data_recv_pts[3*j+1]) + *\/ */
    /*       /\*     (ptprintc[2] - data_recv_pts[3*j+2]) *(ptprintc[2] - data_recv_pts[3*j+2])) < 1e-6) { *\/ */
    /*         printf ("pt recu du 13 %d : %12.5e  %12.5e  %12.5e\n", j1, data_recv_pts[3*j1], data_recv_pts[3*j1+1], data_recv_pts[3*j1+2]); */
    /*     /\* } *\/ */
    /*     } */
    /* /\*   } *\/ */
    /* /\* } *\/ */

    PDM_octree_seq_closest_point (octree->octree_seq_id,
                                  i_recv_gnum[lComm],
                                  data_recv_pts,
                                 _closest_octree_pt_id,
                                 _closest_octree_pt_dist2);

    for (int j = 0; j < i_recv_gnum[lComm]; j++) {
      _closest_octree_pt_g_num[j] = -1;
      if ((_closest_octree_pt_id[2*j]!= -1) &&
          (_closest_octree_pt_id[2*j + 1]!= -1)) {
        _closest_octree_pt_g_num[j] =
        octree->g_num[_closest_octree_pt_id[2*j]][_closest_octree_pt_id[2*j+1]];
      }
    }

    double * data_recv_dist2 = data_send_pts;

    PDM_MPI_Alltoallv (_closest_octree_pt_g_num, n_recv_gnum,
                       i_recv_gnum, PDM__PDM_MPI_G_NUM,
                       data_recv_gnum, n_send_gnum,
                       i_send_gnum, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    PDM_MPI_Alltoallv (_closest_octree_pt_dist2, n_recv_gnum,
                       i_recv_gnum, PDM_MPI_DOUBLE,
                       data_recv_dist2, n_send_gnum,
                       i_send_gnum, PDM_MPI_DOUBLE,
                       octree->comm);


    if (n_exch == 1) {

      for (int j = 0; j < lComm; j++) {
        n_send_gnum[j]     = 0;
      }

      for (int j = 0; j < n_pts; j++) {
        /*  if (825 == j) { */
        /*   printf ("Resultat calcul distance :\n"); */
        /* } */
        for (int l = i_boxes[j]; l < i_boxes[j+1]; l++) {
          int iproc = boxes[l];
          int idx   = i_send_gnum[iproc] + n_send_gnum[iproc];
        /* if (825 == j) { */
        /*   printf("iproc dist gnum : %d %12.5e %ld\n",iproc, data_recv_dist2[idx],  data_recv_gnum[idx]); */
        /* } */
          if (data_recv_dist2[idx] < closest_octree_pt_dist2[j]) {
            closest_octree_pt_dist2[j] = data_recv_dist2[idx];
            closest_octree_pt_g_num[j] = data_recv_gnum[idx];
          }
          n_send_gnum[iproc] += 1;
        }
      }

    }

    else {
      for (int j = 0; j < lComm; j++) {
        n_send_gnum[j]     = 0;
        send_counts[j]     = 0;
      }

      for (int j = 0; j < n_pts; j++) {
        for (int l = i_boxes[j]; l < i_boxes[j+1]; l++) {
          int iproc     = boxes[l];
          if ((send_counts[iproc] >= send_bounds[2*iproc]) &&
              (send_counts[iproc] < send_bounds[2*iproc+1])) {
            int idx   = i_send_gnum[iproc] + n_send_gnum[iproc];

            if (data_recv_dist2[idx] < closest_octree_pt_dist2[j]) {
              closest_octree_pt_dist2[j] = data_recv_dist2[idx];
              closest_octree_pt_g_num[j] = data_recv_gnum[idx];
            }
            n_send_gnum[iproc] += 1;
          }
          send_counts[iproc]++;
        }
      }
    }

  }

  if (n_exch > 1) {
    free (__closest_octree_pt_g_num);
    free (__closest_octree_pt_dist2);
  }

  free (stride_ptb);

  free (i_boxes);
  free (boxes);

  free (data_send_pts1);
  free (data_recv_pts1);
  free (n_send_pts1);
  free (n_recv_pts1);
  free (i_send_pts1);
  free (i_recv_pts1);
  free (data_send_gnum1);
  free (data_recv_gnum1);
  free (n_send_gnum1);
  free (n_recv_gnum1);
  free (i_send_gnum1);
  free (i_recv_gnum1);
  if (n_exch > 1) {
    free (data_send_pts2);
    free (data_recv_pts2);
    free (n_send_pts2);
    free (n_recv_pts2);
    free (i_send_pts2);
    free (i_recv_pts2);
    free (data_send_gnum2);
    free (data_recv_gnum2);
    free (n_send_gnum2);
    free (n_recv_gnum2);
    free (i_send_gnum2);
    free (i_recv_gnum2);
    free (send_counts);
    free (send_bounds1);
    free (send_bounds2);
  }

  free (_closest_octree_pt_id);
  free (_closest_octree_pt_dist2);
  free (_closest_octree_pt_g_num);

}


#ifdef	__cplusplus
}
#endif
