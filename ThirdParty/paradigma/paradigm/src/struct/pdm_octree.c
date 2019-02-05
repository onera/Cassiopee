
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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_octree.h"
#include "pdm_octree_seq.h"

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
//
///**
// * \struct _octant_t
// * \brief  Define an octant
// * 
// */
//
//typedef struct  {
//
//  int  ancestor_id; /*!< Ids of ancestor in octree array */
//  PDM_octree_child_t  location_in_ancestor; /*!< Location in ancestor */
//  int  depth;       /*!< Depth in the tree */
//  int  children_id[8]; /*!< Ids of children in octree array */
//  int  idx[9];         /*!< Start index of point list for each octant */
//  int  n_points;       /*!< Number of points in octant*/
//  double extents[6];   /*!< Extents of the node */
//  
//} _octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 * 
 */

typedef struct  {
  int    octree_seq_id;             /*!< Identifier of the associated octree seq */
//  double  extents[6];            /*!< Extents of current process */ 
  double *extents_proc;          /*!< Extents of processes */
//  int    depth_max;              /*!< Maximum depth of the three */
  PDM_MPI_Comm comm;             /*!< MPI communicator */
//  int points_in_leaf_max;        /*!< Maximum number of points in a leaf */
//  double tolerance;              /*!< Relative geometric tolerance */
//  int   n_nodes;                 /*!< Current number of nodes in octree */
//  int   n_nodes_max;             /*!< Maximum number of nodes in octree */
//  int   *n_points;               /*!< Number of points in each cloud */
//  int   t_n_points;              /*!< total number of points */
//  int   n_point_clouds;          /*!< Number of point cloud */
//  const double **point_clouds;         /*!< points cloud */
//  int *point_ids;                /*!< Id's of points in it cloud sorted by octree
//                                      (size: n_points + 1) */
//  int *point_icloud;             /*!< Cloud's of points sorted by octree
//                                      (size: n_points + 1) */
//  _octant_t   *nodes;            /*!< Array of octree nodes
//                                      (size: n_nodes_max) */
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
  
  octree->extents_proc = NULL;
  
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
  
  octree->extents_proc = NULL;
  
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

  free (octree->extents_proc);
  
  PDM_octree_seq_free (octree->octree_seq_id);
  
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
 * 
 */


void
PDM_octree_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords 
)
{
  _octree_t *octree = _get_from_id (id);
  
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

  PDM_octree_seq_build (octree->octree_seq_id);
  
  double * extents = PDM_octree_seq_extents_get (octree->octree_seq_id); 
 
  
  int n_proc;
  PDM_MPI_Comm_size (octree->comm, &n_proc);
  
  octree->extents_proc = malloc (sizeof(double)* n_proc * 6);
  
  PDM_MPI_Allgather (extents, 6, PDM_MPI_DOUBLE,
                     octree->extents_proc, 6, PDM_MPI_DOUBLE,
                     octree->comm);

  
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
 * \brief Processes extents  
 *
 * \param [in]   id                 Identifier 
 * \param [in]   i_proc             Process
 *
 */

const double *
PDM_octree_processes_extents_get
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);

  return octree->extents_proc;
}
