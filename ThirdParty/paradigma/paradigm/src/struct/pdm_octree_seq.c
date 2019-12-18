
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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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

/**
 * \struct _octant_t
 * \brief  Define an octant
 *
 */

typedef struct  {

  int  ancestor_id; /*!< Ids of ancestor in octree array */
  int  is_leaf; /*!< IS a leaf >*/
  PDM_octree_seq_child_t  location_in_ancestor; /*!< Location in ancestor */
  int  depth;       /*!< Depth in the tree */
  int  children_id[8]; /*!< Ids of children in octree array */
  int  range[2];       /*!< Ids of children in octree array */
  int  idx[9];         /*!< Start index of point list for each octant */
  int  n_points;       /*!< Number of points in octant*/
  double extents[6];   /*!< Extents of the node */

} _octant_t;


/**
 * \struct _octree_seq_t
 * \brief  Define an octree
 *
 */

typedef struct  {

  double  extents[6];            /*!< Extents of current process */
  int    depth_max;              /*!< Maximum depth of the three */
  int points_in_leaf_max;        /*!< Maximum number of points in a leaf */
  double tolerance;              /*!< Relative geometric tolerance */
  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */
  int   *n_points;               /*!< Number of points in each cloud */
  int   t_n_points;              /*!< total number of points */
  int   n_point_clouds;          /*!< Number of point cloud */
  const double **point_clouds;   /*!< points cloud */
  int *point_ids;                /*!< Id's of points in it cloud sorted by octree
                                      (size: n_points + 1) */
  int *point_icloud;             /*!< Cloud's of points sorted by octree
                                      (size: n_points + 1) */
  _octant_t   *nodes;            /*!< Array of octree nodes
                                      (size: n_nodes_max) */
} _octree_seq_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees   = NULL;

static const double _eps_default = 1.e-12;

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Compute distance to a box
 *
 * \param [in]   dim        Dimension
 * \param [in]   extents    Box extents
 * \param [in]   coords     Point coords
 * \param [out]  min_dist2  Square of minimum distance
 * \param [out]  max_dist2  Sqaure of maximum distance
 *
 * \return 1 if point is in the box, 0 otherwise
 *
 */

inline static int
_box_dist2
(
const int              dim,
const double          *restrict extents,
const double          *restrict coords,
double                *restrict min_dist2
//double                *restrict max_dist2
)
{

  int inbox = 0;

  double _min_dist2 = 0.;
  //  double _max_dist2 = 0.;
  for (int i = 0; i < dim; i++) {
    double x = coords[i];
    double xmin = extents[i];
    double xmax = extents[i+dim];

    if (x > xmax) {
      double diff_max = x - xmax;
      //double diff_min = x - xmin;
      _min_dist2 += diff_max * diff_max;
      //      _max_dist2 += diff_min * diff_min;
    }

    else if (x < xmin) {
      //double diff_max = x - xmax;
      double diff_min = x - xmin;
      _min_dist2 += diff_min * diff_min;
      //      _max_dist2 += diff_max * diff_max;
    }

    else {
      //double diff_max = x - xmax;
      //double diff_min = x - xmin;
      inbox += 1;
      // _max_dist2 += PDM_MAX (diff_min * diff_min,
      //                       diff_max * diff_max);
    }

  }

  *min_dist2 = _min_dist2;
//  *max_dist2 = _max_dist2;
  return inbox == dim;

}

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _octree_seq_t *
_get_from_id
(
 int  id
)
{
  _octree_seq_t *octree = (_octree_seq_t *) PDM_Handles_get (_octrees, id);

  if (octree == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return octree;
}

/**
 *
 * \brief Build a local octree's leaves.
 *
 * \param[in]  ancestor_id          Ancestor identifier
 * \param[in]  location_in_ancestor Location in ancestor
 * \param[in]  depth                Depth in the tree
 * \param[in]  extents              Extents associated with node:
 *                                  x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 * \param[in]  point_coords         Point coordinates
 * \param[in]  point_ids_tmp        Temporary point indexes
 * \param[in]  pos_tmp              Temporary point position in octree
 * \param[inout]  octree            Current octree structure
 * \param[inout]  point_range       Start and past-the end index in point_idx
 *                                  for current node (size: 2)
 */

static void
_build_octree_seq_leaves(const int       ancestor_id,
                     const PDM_octree_seq_child_t location_in_ancestor,
                     const int       depth,
                     const double    extents[],
                     const double   **point_coords,
                     int             *point_icloud_tmp,
                     int             *point_ids_tmp,
                     _octree_seq_t       *octree,
                     int             point_range[2])
{
  int i, j, k, _n_nodes, _n_points, tmp_size;

  int count[8], idx[9], octant_id[8];
  double mid[3], sub_extents[6];
  _octant_t  *_node;

  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  _n_nodes = octree->n_nodes;
  tmp_size = octree->n_nodes;

  /* Resize octree if necessary */

  if (octree->n_nodes >= octree->n_nodes_max) {
    if (octree->n_nodes == 0) {
      octree->n_nodes = 1;
      octree->n_nodes_max = 8;
    }
    octree->n_nodes_max *= 2;
    octree->nodes = realloc (octree->nodes,
                             octree->n_nodes_max * sizeof(_octant_t));
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  for (j = 0; j < 8; j++) {
    count[j] = 0;
    octant_id[j] = -1;
  }

  if (depth < octree->depth_max && _n_points > octree->points_in_leaf_max) {
    /* Extents center */

    for (j = 0; j < 3; j++) {
      mid[j]= (extents[j] + extents[j + 3]) * 0.5;
    }


    /* Count points in each octant */

    for (i = point_range[0]; i < point_range[1]; i++) {

      for (j = 0, k = 0; j < 3; j++) {
        if (point_coords[octree->point_icloud[i]][octree->point_ids[i]*3 + j] > mid[j])
          k += octant_mask[j];
      }

      count[k] += 1;
    }

    /* Build index */

    idx[0] = 0;
    for (j = 0; j < 8; j++) {
      idx[j+1] = idx[j] + count[j];
    }

    for (j = 0; j < 8; j++) {
      count[j] = 0;
    }

    for (i = point_range[0], j = 0; i < point_range[1]; i++) {

      for (j = 0, k = 0; j < 3; j++) {
        if (point_coords[octree->point_icloud[i]][octree->point_ids[i]*3 + j] > mid[j])
          k += octant_mask[j];
      }

      point_icloud_tmp[idx[k] + count[k]] = octree->point_icloud[i];
      point_ids_tmp[idx[k] + count[k]] = octree->point_ids[i];
      count[k] += 1;
    }

    /* Check if this subdivision is static
       and check coordinates to find multi point */

    for (i = point_range[0], j = 0; i < point_range[1]; i++, j++) {
      octree->point_icloud[i] = point_icloud_tmp[j];
      octree->point_ids[i] = point_ids_tmp[j];
    }

    for (i = 0; i < 9; i++)
      idx[i] = point_range[0] + idx[i];

    /* Build leaves recursively */

    for (i = 0; i < 8; i++) {

      if ((idx[i+1] - idx[i]) > 0) {

        tmp_size++;

        octant_id[i] = tmp_size;

        if (i < 4) {
          sub_extents[0] = extents[0];
          sub_extents[3] = mid[0];
        }
        else {
          sub_extents[0] = mid[0];
          sub_extents[3] = extents[3];
        }
        /* 1.0e-12 term in assert() used to allow for
           truncation error in for xmin = xmax case */
        assert(sub_extents[0] < sub_extents[3] + 1.0e-14);

        if (i%4 < 2) {
          sub_extents[1] = extents[1];
          sub_extents[4] = mid[1];
        }
        else {
          sub_extents[1] = mid[1];
          sub_extents[4] = extents[4];
        }
        assert(sub_extents[1] < sub_extents[4] + 1.0e-14);

        if (i%2 < 1) {
          sub_extents[2] = extents[2];
          sub_extents[5] = mid[2];
        }
        else {
          sub_extents[2] = mid[2];
          sub_extents[5] = extents[5];
        }
        assert(sub_extents[2] < sub_extents[5] + 1.0e-14);

        octree->n_nodes = tmp_size;

        _build_octree_seq_leaves(_n_nodes,
                             (PDM_octree_seq_child_t) i,
                             depth+1,
                             sub_extents,
                             point_coords,
                             point_icloud_tmp,
                             point_ids_tmp,
                             octree,
                             idx + i);

        tmp_size = octree->n_nodes;
      }

    }

  }
  /* Finalize node */


  _node = octree->nodes + _n_nodes;

  _node->range[0] = point_range[0];
  _node->range[1] = point_range[1];

  for (i = 0; i < 9; i++) {
    _node->idx[i] = idx[i];
  }

  for (i = 0; i < 6; i++) {
    _node->extents[i] = extents[i];
  }

  for (i = 0; i < 8; i++) {
    _node->children_id[i] = octant_id[i];
  }

  _node->is_leaf = (_node->children_id[0] == -1) &&
                (_node->children_id[1] == -1) &&
                (_node->children_id[2] == -1) &&
                (_node->children_id[3] == -1) &&
                (_node->children_id[4] == -1) &&
                (_node->children_id[5] == -1) &&
                (_node->children_id[6] == -1) &&
                (_node->children_id[7] == -1);


  _node->ancestor_id = ancestor_id;
  _node->depth = depth;

  _node->n_points = _n_points;
  _node->location_in_ancestor = location_in_ancestor;
}

/**
 *
 * \brief   Compute extents of a point set
 *
 *  \param [in] dim         Space dimension of points to locate_3d
 *  \param [in] n_points    Number of points to locate
 *  \param [in] point_index optional indirection array to point_coords
 *                          (1 to n_points numbering)
 *  \param [in] point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *
 */

static void
_point_extents(const int     dim,
               const int     n_points,
               const int     point_index[],
               const double  point_coords[],
               double        extents[])
{
  int i;
  int j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_index != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_index[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }

    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }
}


/**
 *
 * \brief Build an octree
 *
 * \param[in]  octree    Current octree
 * .
 */

static void
_build_octree
(
_octree_seq_t *octree
)
{

  int point_range[2];

  /* Initialization */

  octree->n_nodes = 0;
  octree->n_nodes_max = 0;
  octree->nodes = NULL;

  for (int i = 0; i < octree->n_point_clouds; i++) {
    octree->t_n_points += octree->n_points[i];
  };

  octree->point_ids = malloc (sizeof(int) * octree->t_n_points);
  octree->point_icloud = malloc (sizeof(int) * octree->t_n_points);

  int cpt = 0;
  for (int i = 0; i < octree->n_point_clouds; i++) {

    int n_points = octree->n_points[i];
    double extents[6];

    if (n_points > 0) {

      _point_extents(3,
                     n_points,
                     NULL,
                     octree->point_clouds[i],
                     extents);

      for (int i1 = 0; i1 < 3; i1++) {
        octree->extents[i1] = PDM_MIN (extents[i1], octree->extents[i1]);
        octree->extents[i1 + 3] = PDM_MAX (extents[i1 + 3], octree->extents[i1 + 3]);
      }

      for (int j = 0; j < n_points; j++) {
        octree->point_ids[cpt] = j;
        octree->point_icloud[cpt] = i;
        cpt +=1;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    double delta = PDM_MAX (octree->tolerance * (octree->extents[i + 3] - octree->extents[i]),
                            _eps_default);
    octree->extents[i] += -delta;
    octree->extents[i + 3] += delta;;
  }

  point_range[0] = 0;
  point_range[1] = octree->t_n_points;

  int *point_ids_tmp = malloc (sizeof(int) * octree->t_n_points);
  int *point_icloud_tmp = malloc (sizeof(int) * octree->t_n_points);

  _build_octree_seq_leaves(-1,
                       (PDM_octree_seq_child_t) 0,
                       -1,
                       octree->extents,
                       (const double **) octree->point_clouds,
                       point_icloud_tmp,
                       point_ids_tmp,
                       octree,
                       point_range);

  octree->n_nodes +=1;

  free (point_ids_tmp);
  free (point_icloud_tmp);

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
 *
 * \return     Identifier
 */

int
PDM_octree_seq_create
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const double tolerance
)
{

  if (_octrees == NULL) {
    _octrees = PDM_Handles_create (4);
  }

  _octree_seq_t *octree = (_octree_seq_t *) malloc(sizeof(_octree_seq_t));

  int id = PDM_Handles_store (_octrees, octree);

  octree->n_point_clouds = n_point_cloud;
  octree->depth_max = depth_max;
  octree->points_in_leaf_max = points_in_leaf_max;
  octree->tolerance = tolerance;

  octree->n_nodes = 0;
  octree->n_nodes_max = 0;

  octree->n_points = malloc (sizeof(double) * n_point_cloud);
  octree->point_clouds = malloc (sizeof(double *) * n_point_cloud);
  for (int i = 0; i < n_point_cloud; i++) {
    octree->n_points[i] = 0;
    octree->point_clouds[i] = NULL;
  }

  octree->point_icloud = NULL;
  octree->point_ids = NULL;
  octree->nodes = NULL;
  octree->t_n_points = 0;
  for (int i = 0; i < 3; i++) {
    octree->extents[i]     =  HUGE_VAL;
    octree->extents[i + 3] = -HUGE_VAL;
  }

  return id;
}


//void
//PROCF (pdm_octree_seq_create, PDM_OCTREE_SEQ_CREATE)
//(
// const int *n_point_cloud,
// const int *depth_max,
// const int *points_in_leaf_max,
// const double *tolerance,
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
PDM_octree_seq_free
(
 const int          id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  free (octree->n_points);
  free (octree->point_clouds);
  free (octree->point_ids);
  free (octree->nodes);
  free (octree->point_icloud);

  free (octree);

  PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

  const int n_octrees = PDM_Handles_n_get (_octrees);

  if (n_octrees == 0) {
    _octrees = PDM_Handles_free (_octrees);
  }

}

//void
//PROCF (pdm_octree_seq_free, PDM_OCTREE_SEQ_FREE)
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
PDM_octree_seq_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords
)
{
  _octree_seq_t *octree = _get_from_id (id);

  octree->n_points[i_point_cloud] = n_points;
  octree->point_clouds[i_point_cloud] = coords;
}

//void
//PROCF (pdm_octree_seq_point_cloud_set, PDM_OCTREE_SEQ_POINT_CLOUD_SET)
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
PDM_octree_seq_build
(
 const int          id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  if (octree->nodes == NULL) {
    _build_octree (octree);
  }

}

//void
//PROCF (pdm_octree_seq_build, PDM_OCTREE_SEQ_BUILD)
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
PDM_octree_seq_root_node_id_get
(
 const int          id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  if (octree->nodes == NULL) {
    return -1;
  }
  else {
    return 0;
  }
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
PDM_octree_seq_extents_get
(
 const int          id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  return octree->extents;
}

//void
//PROCF (pdm_octree_seq_root_node_id_get, PDM_OCTREE_SEQ_ROOT_NODE_ID_GET)
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
PDM_octree_seq_ancestor_node_id_get
(
 const int          id,
 const int          node_id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].ancestor_id;
}

//void
//PROCF (pdm_octree_seq_ancestor_node_id_get, PDM_OCTREE_SEQ_ANCESTOR_NODE_ID_GET)
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
PDM_octree_seq_node_extents_get
(
 const int          id,
 const int          node_id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].extents;
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
PDM_octree_seq_children_get
(
 const int                id,
 const int                node_id,
 const PDM_octree_seq_child_t child
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].children_id[child];
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
PDM_octree_seq_neighbor_get
(
 const int                    id,
 const int                    node_id,
 const PDM_octree_seq_direction_t direction
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  int neighbor_id = -1;

  _octant_t *octant = &(octree->nodes[node_id]);

  int isSameDiretion = 1;
  while (isSameDiretion && (octant->ancestor_id != 0)) {

    switch (direction) {
    case PDM_OCTREE_SEQ_NADIR:
    case PDM_OCTREE_SEQ_ZENITH:
      isSameDiretion = octant->location_in_ancestor%2  == direction;
      break;
    case PDM_OCTREE_SEQ_WEST:
    case PDM_OCTREE_SEQ_EAST:
      isSameDiretion = ((octant->location_in_ancestor%4 < 2) + 2) == direction;
      break;
    case PDM_OCTREE_SEQ_NORTH:
    case PDM_OCTREE_SEQ_SOUTH:
      isSameDiretion = ((octant->location_in_ancestor < 4) + 4) == direction;
      break;
    }

   octant = &(octree->nodes[octant->ancestor_id]);
  }

  if (octant->ancestor_id != 0) {
    if (direction < 2) {
      neighbor_id = (octant->location_in_ancestor/2) * 2 +
                    (octant->location_in_ancestor + 1)%2;
    }
    if (direction < 4) {
      neighbor_id = (octant->location_in_ancestor/4) * 4 +
                    (octant->location_in_ancestor + 2)%4;
    }
    else {
      neighbor_id = (octant->location_in_ancestor + 4)%8;
    }
  }

  return neighbor_id;
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
PDM_octree_seq_n_points_get
(
 const int                id,
 const int                node_id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].range[1] - octree->nodes[node_id].range[0];

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
PDM_octree_seq_points_get
(
 const int                id,
 const int                node_id,
 int                    **point_clouds_id,
 int                    **point_indexes
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  *point_clouds_id = octree->point_icloud + octree->nodes[node_id].range[0];

  *point_indexes = octree->point_ids + octree->nodes[node_id].range[0];
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
PDM_octree_seq_leaf_is
(
 const int                id,
 const int                node_id
)
{
  _octree_seq_t *octree = _get_from_id (id);

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].is_leaf;
}


/**
 *
 * Look for closest points stored inside an octree
 *
 * parameters:
 * \param [in]   id                     Identifier
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_octree_pt_id   Closest point in octree index
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *
 */

void
PDM_octree_seq_closest_point
(
const int          id,
const int        n_pts,
double          *pts,
int             *closest_octree_pt_id,
double          *closest_octree_pt_dist2
)
{

  const int n_children = 8;

  _octree_seq_t *octree = _get_from_id (id);

  int s_pt_stack = ((n_children - 1) * (octree->depth_max - 1) + n_children);
  int sort_child[n_children];
  double dist_child[n_children];
  int inbox_child[n_children];

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);

  int dim = 3;
  /* double ptprintc[3] = {5.83333e-01,  6.25000e-01,  5.00000e-01}; */
  /* double ptprintc2[3] = {5.83333e-01,  1.,  5.00000e-01}; */

  for (int i = 0; i < n_pts; i++) {

    int pos_stack = 0;
    const double *_pt = pts + dim * i;

      /* double t1 = (_pt[0] - ptprintc[0]) * (_pt[0] - ptprintc[0]); */
      /* double t2 = (_pt[1] - ptprintc[1]) * (_pt[1] - ptprintc[1]); */
      /* double t3 = (_pt[2] - ptprintc[2]) * (_pt[2] - ptprintc[2]); */
      /* t1=1; */
      /* if ((t1+t2+t3) < 1e-6) */
      /* printf ("\n ******** pt %d : %12.5e, %12.5e, %12.5e\n", i, _pt[0], _pt[1], _pt[2]); */

    /* Init stack */

    closest_octree_pt_id[2*i] = -1;
    closest_octree_pt_id[2*i+1] = -1;
    closest_octree_pt_dist2[i] = HUGE_VAL;

    stack[pos_stack] = 0; /* push root in th stack */

    double _min_dist2;
    _octant_t *root_node = &(octree->nodes[0]);
    int inbox1 =  _box_dist2 (dim,
                              root_node->extents,
                              _pt,
                              &_min_dist2);

    inbox_stack[pos_stack] = inbox1;
    min_dist2_stack[pos_stack] = _min_dist2;
    pos_stack++;

    while (pos_stack > 0) {

      int id_curr_node = stack[--pos_stack];
      _octant_t *curr_node = &(octree->nodes[id_curr_node]);

      double min_dist2 = min_dist2_stack[pos_stack];
      int inbox =  inbox_stack[pos_stack];

      /* int inbox =  _box_dist2 (dim, */
      /*                          curr_node->extents, */
      /*                          _pt, */
      /*                          &min_dist2); */

      if ((min_dist2 <= closest_octree_pt_dist2[i]) || (inbox == 1)) {

        if (!curr_node->is_leaf) {

          /* Sort children and store them into the stack */

          const int *_child_ids = curr_node->children_id;

          for (int j = 0; j < n_children; j++) {
            dist_child[j] = HUGE_VAL;
          }

          int n_selec = 0;
          for (int j = 0; j < n_children; j++) {


            int child_id = _child_ids[j];


            int child_inbox = 0;

            if (child_id != -1) {

              _octant_t *child_node = &(octree->nodes[child_id]);
              double child_min_dist2;

              child_inbox = _box_dist2 (dim,
                                        child_node->extents,
                                        _pt,
                                        &child_min_dist2);

              int i1 = 0;
              for (i1 = n_selec;
                   (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
                dist_child[i1] = dist_child[i1-1];
                sort_child[i1] = sort_child[i1-1];
                inbox_child[i1] = inbox_child[i1-1];
              }

              sort_child[i1] = child_id;
              dist_child[i1] = child_min_dist2;
              inbox_child[i1] = child_inbox;

              n_selec += 1;

            }
          }

          for (int j = 0; j < n_selec; j++) {
            int j1 = n_selec- 1 - j;
            int child_id = sort_child[j1];
            if (child_id != -1) {
              _octant_t *child_node = &(octree->nodes[child_id]);
              if ((dist_child[j1] < closest_octree_pt_dist2[i]) &&
                  (child_node->n_points > 0)) {

                min_dist2_stack[pos_stack] = dist_child[j1];
                inbox_stack[pos_stack] = inbox_child[j1];

                stack[pos_stack++] = child_id; /* push root in th stack */
              }
            }
          }
        }

        else {

          int *point_clouds_id = octree->point_icloud + curr_node->range[0];
          int *point_indexes = octree->point_ids + curr_node->range[0];

          for (int j = 0; j < curr_node->n_points; j++) {

            double point_dist2 = 0;
            const double *_coords = octree->point_clouds[point_clouds_id[j]]
                                    + dim * point_indexes[j];

            for (int k = 0; k < dim; k++) {
              point_dist2 += (_coords[k] - _pt[k]) *
                             (_coords[k] - _pt[k]);
            }

            if (point_dist2 < closest_octree_pt_dist2[i]) {
              closest_octree_pt_id[2*i] = point_clouds_id[j];
              closest_octree_pt_id[2*i+1] = point_indexes[j];
              closest_octree_pt_dist2[i] = point_dist2;
            }
          }
        }
      }
    }
      /*     if ((t1+t2+t3) < 1e-6) */
      /* printf ("\n ******** fin point ******************\n"); */

  }

  free (inbox_stack);
  free (min_dist2_stack);
  free (stack);

}

#ifdef	__cplusplus
}
#endif
