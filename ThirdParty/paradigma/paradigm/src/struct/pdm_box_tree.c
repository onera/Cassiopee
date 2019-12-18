/*============================================================================
 * Search octrees and quadtrees of boxes.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_box_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_box_tree.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro and Type definitions
 *============================================================================*/

#define PDM_BOX_TREE_MAX_BUILD_LOOPS 50

/* Structures for each octant or quadrant */
/*----------------------------------------*/

/* If the type is BOX_TREE_NODE, the ordering of children is defined as follows,
   using notation B: bottom, U: up, E: east, W: west, S: south,  N: north.

   octant:   0: BSW, 1: BSE, 2: BNW, 3: BNE, 4: USW, 5: USE, 6: UNW, 7: UNE
   quadrant: 0:  SW, 1:  SE, 2:  NW, 3:  NE
   segment:  0:   W, 1:   E
 */

typedef struct {

  _Bool              is_leaf;      /* True for leaf nodes */

  PDM_morton_code_t  morton_code;  /* Level and coordinates in the grid
                                      according to Morton encoding */

  int   n_boxes;             /* Number of associated bounding boxes */
  int   start_id;            /* Position of the first box_id */

} _node_t;

/* Structure used to manage statistics */

typedef struct {

  unsigned    max_level_reached;  /* Max level number reached */

  int   n_leaves;           /* Number of leaves in the tree */
  int   n_boxes;            /* Number of boxes to locate in the tree */
  int   n_linked_boxes;     /* Number of linked boxes in the tree */
  int   n_spill_leaves;     /* Number of leaves where n_boxes > threshold */

  int   min_linked_boxes;   /* Minimum number of boxes for a leaf */
  int   max_linked_boxes;   /* Maximum number of boxes for a leaf */

} PDM_box_tree_stats_t;

/* Main box tree structure */
/*-------------------------*/

struct _PDM_box_tree_t {

  int               n_children;      /* 8, 4, or 2 (2^dim) */

  int               max_level;       /* Max. possible level */
  int         threshold;       /* Max number of boxes linked to a
                                        node if max_level is not reached */
  float             max_box_ratio;   /* Max n_linked_boxes / n_boxes value */

  PDM_box_tree_stats_t stats;        /* Statistics related to the structure */

  int         n_max_nodes;     /* Current max. allocated nodes */
  int         n_nodes;         /* Number of nodes (including leaves) */

  _node_t          *nodes;           /* Array of nodes (root at index 0) */

  int        *child_ids;       /* Ids of associated children
                                        (size: 2^dim * n_max_nodes) */

  int        *box_ids;         /* List of associated box ids.
                                        size = stat.n_linked_boxes */

  int        *stack;           /* Stack for look for closest leaves */

  int        *pos_stack;       /* Current position in the stack */


  int     n_build_loops;             /* Number of loops required to build */

  PDM_MPI_Comm          comm;            /* Associated MPI communicator */

  const PDM_box_set_t  *boxes;       /* Associated boxes */
};

/*=============================================================================
 * Static global variables
 *============================================================================*/

static int iappel = 0;
/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Build minimal octree between two octants
 *
 * \param [in]     octree    Current octree
 * \param [in]     code      Morton code
 * \param [inout]  extents   Extents associated to the Morton code
 *
 */

static void
_extents
(
 const int dim,
 PDM_morton_code_t code,
 double    extents[]
)
{
  for (int i = 0; i < dim; i++) {
    extents[i] = (double) code.X[i]/(double) pow(2,code.L);
    extents[dim + i] = ((double) code.X[i] + 1)/(double) pow(2,code.L);
  }
}

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
_box_dist2_max
(
const int              dim,
const int              normalized,
const double          *restrict d,
const double          *restrict extents,
const double           *restrict coords,
double                *restrict max_dist2
)
{

  int inbox = 0;
  *max_dist2 = 0.;


  if (normalized) {

    for (int i = 0; i < dim; i++) {
      if (coords[i] > extents[i+dim]) {
        double _max_dist2 = d[i] * (coords[i] - extents[i]);
        *max_dist2 += _max_dist2 * _max_dist2;
      }

      else if (coords[i] < extents[i]) {
        double _max_dist2 = d[i] * (coords[i] - extents[dim+i]);
        *max_dist2 += _max_dist2 * _max_dist2;
      }

      else {
        inbox += 1;
        double val1 = d[i] * (coords[i] - extents[i]);
        double val2 = d[i] * (coords[i] - extents[dim+i]);
        *max_dist2 += PDM_MAX (val1 * val1, val2 * val2);
      }
    }

  }

  else {

    for (int i = 0; i < dim; i++) {
      if (coords[i] > extents[i+dim]) {
        double _max_dist2 = coords[i] - extents[i];
        *max_dist2 += _max_dist2 * _max_dist2;
      }

      else if (coords[i] < extents[i]) {
        double _max_dist2 = coords[i] - extents[dim+i];
        *max_dist2 += _max_dist2 * _max_dist2;
      }

      else {
        inbox += 1;
        double val1 = coords[i] - extents[i];
        double val2 = coords[i] - extents[dim+i];
        *max_dist2 += PDM_MAX (val1 * val1, val2 * val2);
      }
    }
  }

  return inbox == dim;

}


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
_box_dist2_min
(
const int              dim,
const int              normalized,
const double          *restrict d,
const double          *restrict extents,
const double           *restrict coords,
double                *restrict min_dist2
)
{

  int inbox = 0;
  *min_dist2 = 0.;


  if (normalized) {

    for (int i = 0; i < dim; i++) {
      if (coords[i] > extents[i+dim]) {
        double _min_dist2 = d[i] * (coords[i] - extents[dim+i]);
        *min_dist2 += _min_dist2 * _min_dist2;
      }

      else if (coords[i] < extents[i]) {
        double _min_dist2 = d[i] * (coords[i] - extents[i]);
        *min_dist2 += _min_dist2 * _min_dist2;
      }

      else {
        inbox += 1;
      }
    }

  }

  else {

    for (int i = 0; i < dim; i++) {
      if (coords[i] > extents[i+dim]) {
        double _min_dist2 = coords[i] - extents[dim+i];
        *min_dist2 += _min_dist2 * _min_dist2;
      }

      else if (coords[i] < extents[i]) {
        double _min_dist2 = coords[i] - extents[i];
        *min_dist2 += _min_dist2 * _min_dist2;
      }

      else {
        inbox += 1;
      }
    }
  }

  return inbox == dim;

}


/**
 *
 * \brief Add children nodes into stack
 *
 * \param [in]    bt            Box tree
 * \param [in]    dim           Dimension
 * \param [in]    id_curr_node  Identifier of current node
 * \param [in]    upper_bound   Upper_bound criteria to store a child in stack
 * \param [in]    pt            Distance to this point must be lesser van upper_bound
 * \param [inout] pos_stack     Position in the stack
 * \param [inout] stack         Stack
 *
 */

inline static void
_push_child_in_stack_v0
(
PDM_box_tree_t *bt,
const int       dim,
const int              normalized,
const double          *restrict d,
const int       id_curr_node,
const double    upper_bound,
const double    *restrict pt,
int             *restrict pos_stack,
int             *restrict stack,
int             *restrict inbox_stack,
double          *restrict min_dist2_stack,
int             flag,
int             sorted
)
{
  int sort_child[bt->n_children];
  double dist_child[bt->n_children];
  int inbox_child[bt->n_children];

  for (int i = 0; i < bt->n_children; i++) {
    dist_child[i] = HUGE_VAL;
  }

  /* Sort children and store them into the stack */

  const int *_child_ids = bt->child_ids + id_curr_node*bt->n_children;

  int _n_push = 0;

  double child_extents2[2*dim];

  for (int j = 0; j < bt->n_children; j++) {

    double child_min_dist2;

    int child_id = _child_ids[j];

    //const double *child_extents = bt->extents + dim * 2 * child_id;

    _node_t *curr_node = &(bt->nodes[child_id]);

    if (curr_node->n_boxes == 0) {
      continue;
    }

    _extents (dim, curr_node->morton_code, child_extents2);

    int inbox = _box_dist2_min (dim,
                                normalized,
                                d,
                                child_extents2,
                                pt,
                                &child_min_dist2);

    if (sorted) {
      int i1 = 0;
      for (i1 = _n_push; (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
        dist_child[i1] = dist_child[i1-1];
        sort_child[i1] = sort_child[i1-1];
        inbox_child[i1] = inbox_child[i1-1];
      }

      sort_child[i1] = _child_ids[j];
      dist_child[i1] =  child_min_dist2;
      inbox_child[i1] =  inbox;

      _n_push += 1;
    }
    else {
      if (((child_min_dist2 < upper_bound) || (inbox == 1))) {

        stack[*pos_stack]           = child_id; /* push child in th stack */
        inbox_stack[*pos_stack]     = inbox;
        min_dist2_stack[*pos_stack] = child_min_dist2;

      (*pos_stack)++;
      }
    }

  }

  if (sorted) {
    for (int j =  _n_push - 1; j >= 0; j--) {
      int child_id = sort_child[j];

      if (((dist_child[j] < upper_bound) || (inbox_child[j] == 1))) {

        stack[*pos_stack]           = child_id; /* push child in th stack */
        inbox_stack[*pos_stack]     = inbox_child[j];
        min_dist2_stack[*pos_stack] = dist_child[j];

        (*pos_stack)++;
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Get minimum coordinates for a given box.
 *
 * parameters:
 *   box_set <-- pointer to box set structure
 *   box_id  <-- id of box
 *
 * returns:
 *   pointer to minimum box coordinates
 *---------------------------------------------------------------------------*/

inline static const double *
_box_min(const PDM_box_set_t   *boxes,
         int              box_id)
{
  return boxes->extents + box_id*boxes->dim*2;
}

/*----------------------------------------------------------------------------
 * Get maxmum coordinates for a given box.
 *
 * parameters:
 *   box_set <-- pointer to box set structure
 *   box_id  <-- id of box
 *
 * returns:
 *   pointer to maximum box coordinates
 *---------------------------------------------------------------------------*/

inline static const double *
_box_max(const PDM_box_set_t   *boxes,
         int              box_id)
{
  return boxes->extents + box_id*boxes->dim*2 + boxes->dim;
}


/*----------------------------------------------------------------------------
 * Test for intersection between two bounding boxes.
 *
 * parameters:
 *   extents <-- array of box extents
 *   id_0    <-- id of first box in extents
 *   extentsB<-- box extents to intersect
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

inline static _Bool
_boxes_intersect_3d(const double  *extents,
                    int            id_0,
                    const double  *extentsB)
{
  const double *e0 = extents + id_0*6;

  if (   e0[0] > extentsB[3] || extentsB[0] > e0[3]
      || e0[1] > extentsB[4] || extentsB[1] > e0[4]
      || e0[2] > extentsB[5] || extentsB[2] > e0[5])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intersect_2d(const double  *extents,
                    int            id_0,
                    const double  *extentsB)
{
  const double *e0 = extents + id_0*4;

  if (   e0[0] > extentsB[2] || extentsB[0] > e0[2]
      || e0[1] > extentsB[3] || extentsB[1] > e0[3])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intersect_1d(const double  *extents,
                    int            id_0,
                    const double  *extentsB)
{
  const double *e0 = extents + id_0*2;

  if (   e0[0] > extentsB[1] || extentsB[0] > e0[1])
    return false;
  else
    return true;
}


/*----------------------------------------------------------------------------
 * Test for intersection between two bounding boxes.
 *
 * parameters:
 *   extents <-- array of box extents
 *   id_0    <-- id of first box
 *   id_1    <-- id of second box
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

inline static _Bool
_boxes_intern_intersect_3d(const double  *extents,
                           int          id_0,
                           int          id_1)
{
  const double *e0 = extents + id_0*6;
  const double *e1 = extents + id_1*6;

  if (   e0[0] > e1[3] || e1[0] > e0[3]
      || e0[1] > e1[4] || e1[1] > e0[4]
      || e0[2] > e1[5] || e1[2] > e0[5])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intern_intersect_2d(const double  *extents,
                           int          id_0,
                           int          id_1)
{
  const double *e0 = extents + id_0*4;
  const double *e1 = extents + id_1*4;

  if (   e0[0] > e1[2] || e1[0] > e0[2]
      || e0[1] > e1[3] || e1[1] > e0[3])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intern_intersect_1d(const double  *extents,
                           int          id_0,
                           int          id_1)
{
  const double *e0 = extents + id_0*2;
  const double *e1 = extents + id_1*2;

  if (   e0[0] > e1[1] || e1[0] > e0[1])
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------
 * Update octree stat structure (min, max, mean, box ratio, ...)
 *
 * parameters:
 *   bt      <-> pointer on the PDM_box_tree_t structure to deal with
 *   node_id <-- node on which we collect data
 *----------------------------------------------------------------------------*/

static void
_update_tree_stats(PDM_box_tree_t  *bt,
                   int        node_id)
{
  int  i;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {
    int n_children = bt->n_children;
    const int *_child_ids = bt->child_ids + node_id*bt->n_children;
    for (i = 0; i < n_children; i++)
      _update_tree_stats(bt, _child_ids[i]);
  }

  else { /* leaf node */

    PDM_box_tree_stats_t  s = bt->stats;

    s.n_leaves += 1;
    s.n_linked_boxes += node->n_boxes;

    if (node->n_boxes > bt->threshold)
      s.n_spill_leaves += 1;

    s.min_linked_boxes = PDM_MIN(s.min_linked_boxes, node->n_boxes);
    s.max_linked_boxes = PDM_MAX(s.max_linked_boxes, node->n_boxes);
    s.max_level_reached = PDM_MAX(s.max_level_reached, node->morton_code.L);

    bt->stats = s;
  }
}

/*----------------------------------------------------------------------------
 * Define box_tree->stat structure (min, max, mean, box ratio, ...)
 *
 * parameters:
 *   bt <-> pointer to the box-tree structure
 *----------------------------------------------------------------------------*/

static void
_get_box_tree_stats(PDM_box_tree_t  *bt)
{
  if (bt == NULL)
    return;

  /* Initialize statistics */

  bt->stats.max_level_reached = 0;

  bt->stats.n_leaves = 0;
  bt->stats.n_spill_leaves = 0;
  bt->stats.n_linked_boxes = 0;

  bt->stats.min_linked_boxes = INT_MAX;
  bt->stats.max_linked_boxes = 0;

  /* Recursively update stats, starting from root */

  if (bt->nodes != NULL)
    _update_tree_stats(bt, 0);
}

/*----------------------------------------------------------------------------
 * Get the coordinates in the grid for the current point at this level.
 *
 * parameters:
 *   level      <--   level on which we want the coordinates
 *   coords     <--   coords of the point to translate in the octree grid.
 *   XYZ        <->   pointer to the X, Y, Z coordinates in the grid
 *----------------------------------------------------------------------------*/

inline static void
_get_grid_coords_3d(PDM_morton_int_t  level,
                    const double      coords[3],
                    double            XYZ[])
{
  PDM_morton_int_t  refinement = 1 << level;

  XYZ[0] = coords[0] * refinement;
  XYZ[1] = coords[1] * refinement;
  XYZ[2] = coords[2] * refinement;
}

inline static void
_get_grid_coords_2d(PDM_morton_int_t  level,
                    const double      coords[2],
                    double            XY[])
{
  PDM_morton_int_t  refinement = 1 << level;

  XY[0] = coords[0] * refinement;
  XY[1] = coords[1] * refinement;
}

inline static void
_get_grid_coords_1d(PDM_morton_int_t  level,
                    const double      coords[1],
                    double            X[])
{
  PDM_morton_int_t  refinement = 1 << level;

  X[0] = coords[0] * refinement;
}

/*----------------------------------------------------------------------------
 * Return true if a leaf intersects a box, false otherwise.
 *
 * parameters:
 *   morton_code <-- Morton code of the leaf
 *   min_box     <-- coordinates of min. point of the bounding box
 *   max_box     <-- coordinates of max. point of the bounding box
 *
 * returns:
 *   true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_node_intersect_box_3d(PDM_morton_code_t  morton_code,
                       const double   min_box[3],
                       const double   max_box[3])
{
  int  i;
  double  min_oct[3], max_oct[3];

  for (i = 0; i < 3; i++) {
    min_oct[i] = (double)morton_code.X[i];
    max_oct[i] = (double)(morton_code.X[i] + 1);
  }

  /* printf ("min_oct : %12.5e %12.5e %12.5e\n", min_oct[0], min_oct[1], min_oct[2]); */
  /* printf ("max_oct : %12.5e %12.5e %12.5e\n", max_oct[0], max_oct[1], max_oct[2]); */

  /* printf ("min_box : %12.5e %12.5e %12.5e\n", min_box[0], min_box[1], min_box[2]); */
  /* printf ("max_box : %12.5e %12.5e %12.5e\n\n", max_box[0], max_box[1], max_box[2]); */


  if (   min_box[0] > max_oct[0] || min_oct[0] > max_box[0]
      || min_box[1] > max_oct[1] || min_oct[1] > max_box[1]
         || min_box[2] > max_oct[2] || min_oct[2] > max_box[2]){
    return false;
  }
  else{
    return true;
  }
}

inline static _Bool
_node_intersect_box_2d(PDM_morton_code_t  morton_code,
                       const double   min_box[2],
                       const double   max_box[2])
{
  int  i;
  double  min_oct[2], max_oct[2];

  for (i = 0; i < 2; i++) {
    min_oct[i] = (double)morton_code.X[i];
    max_oct[i] = (double)(morton_code.X[i] + 1);
  }

  if (   min_box[0] > max_oct[0] || min_oct[0] > max_box[0]
      || min_box[1] > max_oct[1] || min_oct[1] > max_box[1])
    return false;
  else
    return true;
}

inline static _Bool
_node_intersect_box_1d(PDM_morton_code_t  morton_code,
                       const double   min_box[1],
                       const double   max_box[1])
{
  double  min_oct, max_oct;

  min_oct = (double)morton_code.X[0];
  max_oct = (double)(morton_code.X[0] + 1);

  if (min_box[0] > max_oct || min_oct > max_box[0])
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------
 * Split a node into its children and evaluate the box distribution.
 *
 * parameters:
 *   bt       <-> pointer to the box tree being built
 *   boxes    <-- pointer to the associated box set structure
 *   node_id  <-- id of the node to split
 *
 * returns:
 *   the number of associations between nodes and their children
 *----------------------------------------------------------------------------*/

static int
_evaluate_splitting_3d(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       int             node_id)
{
  int   i, j;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[8];

  int  n_linked_boxes = 0;

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(boxes->dim == 3);

  /* Define a Morton code for each child */

  PDM_morton_get_children(3, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[3], max_grid_coord[3];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(3, next_level, box_min);
    max_code = PDM_morton_encode(3, next_level, box_max);

    if (   PDM_morton_compare(3, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {
        if (_node_intersect_box_3d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same octant */

      assert(   PDM_morton_compare(3, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {
        if (   PDM_morton_compare(3, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}

static int
_evaluate_splitting_2d(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       int             node_id)
{
  int   i, j;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[4];

  int  n_linked_boxes = 0;

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(boxes->dim == 2);

  /* Define a Morton code for each child */

  PDM_morton_get_children(2, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(2, next_level, box_min);
    max_code = PDM_morton_encode(2, next_level, box_max);

    if (   PDM_morton_compare(2, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {
        if (_node_intersect_box_2d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(2, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {
        if (   PDM_morton_compare(2, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}

static int
_evaluate_splitting_1d(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       int             node_id)
{
  int   i, j;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[2];

  int  n_linked_boxes = 0;

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  /* Define a Morton code for each child */

  PDM_morton_get_children(1, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[1], max_grid_coord[1];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(1, next_level, box_min);
    max_code = PDM_morton_encode(1, next_level, box_max);

    if (   PDM_morton_compare(1, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {
        if (_node_intersect_box_1d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(1, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {
        if (   PDM_morton_compare(1, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}


/*----------------------------------------------------------------------------
 * Evaluate the intersection between a box and a node tree
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxesB          <--  pointer boxes that intersect associated tree boxes
 *   boxB_id         <--  id of the current boxB in boxesB
 *   node_id         <--  id of the starting node
 *----------------------------------------------------------------------------*/

static _Bool
_boxes_intersect_node_3d (const PDM_box_tree_t     *bt,
                          const PDM_box_set_t      *boxesB,
                          int                       boxB_id,
                          int                       node_id)
{
  PDM_morton_code_t  min_code, max_code;

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  level = node.morton_code.L;

  assert(boxesB->dim == 3);

  double  min_grid_coord[3], max_grid_coord[3];

  const double  *box_min = _box_min(boxesB, boxB_id);
  const double  *box_max = _box_max(boxesB, boxB_id);

  min_code = PDM_morton_encode(3, level, box_min);
  max_code = PDM_morton_encode(3, level, box_max);

  if (   PDM_morton_compare(3, min_code, max_code)
         == PDM_MORTON_DIFFERENT_ID) {

    _get_grid_coords_3d(level, box_min, min_grid_coord);
    _get_grid_coords_3d(level, box_max, max_grid_coord);

    return _node_intersect_box_3d (node.morton_code, min_grid_coord, max_grid_coord);

  }

  else { /* Box is included in the same octant */

    assert(   PDM_morton_compare(3, max_code, min_code)
              == PDM_MORTON_EQUAL_ID);

    return (PDM_morton_compare(3, min_code, node.morton_code) == PDM_MORTON_EQUAL_ID);

  }
}


static _Bool
_boxes_intersect_node_2d (const PDM_box_tree_t     *bt,
                          const PDM_box_set_t      *boxesB,
                          int                       boxB_id,
                          int                       node_id)
{
  PDM_morton_code_t  min_code, max_code;

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  level = node.morton_code.L;

  assert(boxesB->dim == 2);

  double  min_grid_coord[2], max_grid_coord[2];

  const double  *box_min = _box_min(boxesB, boxB_id);
  const double  *box_max = _box_max(boxesB, boxB_id);

  min_code = PDM_morton_encode(2, level, box_min);
  max_code = PDM_morton_encode(2, level, box_max);

  if (   PDM_morton_compare(2, min_code, max_code)
         == PDM_MORTON_DIFFERENT_ID) {

    _get_grid_coords_2d(level, box_min, min_grid_coord);
    _get_grid_coords_2d(level, box_max, max_grid_coord);

    return _node_intersect_box_2d (node.morton_code, min_grid_coord, max_grid_coord);

  }

  else { /* Box is included in the same octant */

    assert(   PDM_morton_compare(2, max_code, min_code)
              == PDM_MORTON_EQUAL_ID);

    return (PDM_morton_compare(2, min_code, node.morton_code) == PDM_MORTON_EQUAL_ID);

  }
}


static _Bool
_boxes_intersect_node_1d (const PDM_box_tree_t     *bt,
                          const PDM_box_set_t      *boxesB,
                          int                       boxB_id,
                          int                       node_id)
{
  PDM_morton_code_t  min_code, max_code;

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  level = node.morton_code.L;

  assert(boxesB->dim == 1);

  double  min_grid_coord[1], max_grid_coord[1];

  const double  *box_min = _box_min(boxesB, boxB_id);
  const double  *box_max = _box_max(boxesB, boxB_id);

  min_code = PDM_morton_encode(1, level, box_min);
  max_code = PDM_morton_encode(1, level, box_max);

  if (   PDM_morton_compare(1, min_code, max_code)
         == PDM_MORTON_DIFFERENT_ID) {

    _get_grid_coords_1d(level, box_min, min_grid_coord);
    _get_grid_coords_1d(level, box_max, max_grid_coord);

    return _node_intersect_box_1d (node.morton_code, min_grid_coord, max_grid_coord);

  }

  else { /* Box is included in the same octant */

    assert(   PDM_morton_compare(1, max_code, min_code)
              == PDM_MORTON_EQUAL_ID);

    return (PDM_morton_compare(1, min_code, node.morton_code) == PDM_MORTON_EQUAL_ID);

  }
}


/*----------------------------------------------------------------------------
 * Evaluate the intersection between a box and a node tree
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxesB          <--  pointer boxes that intersect associated tree boxes
 *   boxB_id         <--  id of the current boxB in boxesB
 *   node_id         <--  id of the starting node
 *----------------------------------------------------------------------------*/

static _Bool
_boxes_intersect_node (const PDM_box_tree_t     *bt,
                       const PDM_box_set_t      *boxesB,
                       int                       boxB_id,
                       int                       node_id)
{
  _Bool is_intersect = false;
  if (boxesB->dim == 3) {
    is_intersect = _boxes_intersect_node_3d (bt, boxesB, boxB_id, node_id);
  }
  else if (boxesB->dim == 2) {
    is_intersect = _boxes_intersect_node_2d (bt, boxesB, boxB_id, node_id);
  }
  else if (boxesB->dim == 1) {
    is_intersect = _boxes_intersect_node_1d (bt, boxesB, boxB_id, node_id);
  }
  return is_intersect;
}



/*----------------------------------------------------------------------------
 * Evaluate the box distribution over the leaves of the box tree to help
 * determine if we should add a level to the tree structure.
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxes           <--  pointer to the associated box set structure
 *   node_id         <--  id of the starting node
 *   build_type      <--  layout variant for building the tree structure
 *   next_level_size -->  size of box_ids for the next level
 *----------------------------------------------------------------------------*/

static void
_count_next_level(PDM_box_tree_t           *bt,
                  const PDM_box_set_t      *boxes,
                  int                 node_id,
                  PDM_box_tree_sync_t       build_type,
                  int                *next_level_size)
{
  int   i;

  _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    assert(bt->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _count_next_level(bt,
                        boxes,
                        bt->child_ids[bt->n_children*node_id + i],
                        build_type,
                        next_level_size);

  }

  else { /* if (node->is_leaf == true) */
    if (   node->n_boxes < bt->threshold
        && node_id != 0                    /* Root node is always divided */
        && build_type == PDM_BOX_TREE_ASYNC_LEVEL)
      *next_level_size += node->n_boxes;

    else { /* Split node and evaluate box distribution between its children */
      if (boxes->dim == 3)
        *next_level_size += _evaluate_splitting_3d(bt, boxes, node_id);
      else if (boxes->dim == 2)
        *next_level_size += _evaluate_splitting_2d(bt, boxes, node_id);
      else if (boxes->dim == 1)
        *next_level_size += _evaluate_splitting_1d(bt, boxes, node_id);
    }
  }
}

/*----------------------------------------------------------------------------
 * Test if we have to continue the building of the box tree.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   build_type <--  layout variant for building the tree structure
 *   next_size  -->  size of box_ids for the next tree if required
 *
 * returns:
 *   true if we should continue, false otherwise.
 *----------------------------------------------------------------------------*/

static _Bool
_recurse_tree_build(PDM_box_tree_t       *bt,
                    const PDM_box_set_t  *boxes,
                    PDM_box_tree_sync_t   build_type,
                    int            *next_size)
{
  int  state = 0;
  int   _next_size = 0;

  _Bool retval = false;

  int  n_ranks = 1;
  PDM_MPI_Comm comm = boxes->comm;

  if (comm != PDM_MPI_COMM_NULL)
    PDM_MPI_Comm_size(comm, &n_ranks);

  bt->n_build_loops += 1;

  if (bt == NULL)
    state = 1;

  /* To avoid infinite loop on tree building */

  if (bt->n_build_loops > PDM_BOX_TREE_MAX_BUILD_LOOPS)
    state = 1;

  /* A sufficient accuracy has been reached */

  if ((int)(bt->stats.max_level_reached) == bt->max_level)
    state = 1;

  /* Algorithm is converged. No need to go further */

  if (   bt->stats.n_spill_leaves == 0
      && bt->stats.max_level_reached > 0)
    state = 1;

  if (n_ranks > 1 && build_type == PDM_BOX_TREE_SYNC_LEVEL) {
    int global_state;
    PDM_MPI_Allreduce(&state, &global_state, 1, PDM_MPI_INT, PDM_MPI_MIN, comm);
    state = global_state; /* Stop if all ranks require it */
  }

  if (state == 0) {

    float box_ratio;

    /* Limit, to avoid excessive memory usage */

    _count_next_level(bt,
                      boxes,
                      0,  /* Starts from root */
                      build_type,
                      &_next_size);

    if (bt->stats.n_boxes > 0)
      box_ratio = (float) ((_next_size*1.0)/bt->stats.n_boxes);
    else
      box_ratio = 0;

    if (bt->stats.max_level_reached > 0 && box_ratio > bt->max_box_ratio)
      state = 1;

  }

  if (n_ranks > 1 && build_type == PDM_BOX_TREE_SYNC_LEVEL) {
    int global_state;
    PDM_MPI_Allreduce(&state, &global_state, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
    state = global_state; /* Stop as as soon as any rank requires it */
  }

  /* If no condition is encoutered, we have to continue */

  *next_size = _next_size;

  if (state == 0)
    retval = true;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a box tree by copying from another.
 *
 * parameters:
 *   dest <-> pointer to destination box tree
 *   src  <-- pointer to source box tree
 *----------------------------------------------------------------------------*/

static void
_copy_tree(PDM_box_tree_t        *dest,
           const PDM_box_tree_t  *src)
{
  assert(dest != NULL && src != NULL);

  memcpy(dest, src, sizeof(PDM_box_tree_t));

  dest->nodes = (_node_t *) malloc(dest->n_max_nodes * sizeof(_node_t));
  dest->child_ids = (int *) malloc(dest->n_max_nodes*dest->n_children * sizeof(int));

  dest->box_ids = (int *) malloc((dest->stats).n_linked_boxes * sizeof(int));

  memcpy(dest->nodes, src->nodes, dest->n_nodes * sizeof(_node_t));
  memcpy(dest->child_ids,
         src->child_ids,
         dest->n_nodes * src->n_children * sizeof(int));

  memcpy(dest->box_ids,
         src->box_ids,
         (dest->stats).n_linked_boxes * sizeof(int));
}

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to PDM_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

static void
_free_tree_arrays(PDM_box_tree_t  *bt)
{
  assert(bt != NULL);

  free(bt->nodes);
  free(bt->child_ids);
  free(bt->box_ids);
}

/*----------------------------------------------------------------------------
 * Create a new node from a Morton code and associate it to a tree.
 *
 * parameters:
 *   bt          <->  pointer to the box tree being built
 *   morton_code <--  Morton identification number related to the node
 *   node_id     <--  id of the starting node
 *----------------------------------------------------------------------------*/

static inline void
_new_node(PDM_box_tree_t     *bt,
          PDM_morton_code_t   morton_code,
          int                 node_id)
{
  int  i;
  _node_t *node;

  assert(bt != NULL);

  node = bt->nodes + node_id;

  if ((int)(morton_code.L) > bt->max_level) {
    PDM_error(__FILE__, __LINE__, 0,
            "Error adding a new node in box tree (%p).\n"
            "Max level reached. Current level: %u and Max level: %d\n",
            (void *)bt, morton_code.L, bt->max_level);
    abort();
  }

  node->is_leaf = true;
  node->morton_code = morton_code;

  node->n_boxes = 0;
  node->start_id = -1; /* invalid value by default */

  for (i = 0; i < bt->n_children; i++) {
    bt->child_ids[node_id*bt->n_children + i] = -1;
  }

}

/*----------------------------------------------------------------------------
 * Split a node into its children and define the new box distribution.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   next_bt    <->  pointer to the next box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   node_id    <--  id of the starting node
 *   shift_ids  <->  first free position free in new box_ids
 *----------------------------------------------------------------------------*/

static void
_split_node_3d(PDM_box_tree_t       *bt,
               PDM_box_tree_t       *next_bt,
               const PDM_box_set_t  *boxes,
               int             node_id,
               int            *shift_ids)
{
  int j, i;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[8];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->n_nodes;
  _node_t  split_node = next_bt->nodes[node_id];

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 8);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 8 > next_bt->n_max_nodes) {
    assert(next_bt->n_max_nodes > 0);
    next_bt->n_max_nodes *= 2;
    next_bt->nodes = (_node_t *) realloc((void *) next_bt->nodes, next_bt->n_max_nodes * sizeof(_node_t));
    next_bt->child_ids = (int *) realloc((void *) next_bt->child_ids, next_bt->n_max_nodes*8 * sizeof(int));

  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(3, node.morton_code, children);

  for (i = 0; i < 8; i++) {

    const int   new_id = n_init_nodes + i;
    next_bt->child_ids[node_id*8 + i] = new_id;

    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = node.n_boxes;
  split_node.is_leaf = false;

  next_bt->nodes[node_id] = split_node;
  next_bt->n_nodes = n_init_nodes + 8;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {
    double  min_grid_coord[3], max_grid_coord[3];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(3, next_level, box_min);
    max_code = PDM_morton_encode(3, next_level, box_max);

    if (   PDM_morton_compare(3, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {
        if (_node_intersect_box_3d(children[i], min_grid_coord, max_grid_coord)) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;

        }
      }

    }
    else { /* Box is included in the same octant */

      assert(   PDM_morton_compare(3, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {
        if (   PDM_morton_compare(3, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 8; i++) {

    (next_bt->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 8; i++)
    (next_bt->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[3], max_grid_coord[3];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(3, next_level, box_min);
    max_code = PDM_morton_encode(3, next_level, box_max);

    if (   PDM_morton_compare(3, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {

        if (_node_intersect_box_3d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;
          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same octant */

      assert(   PDM_morton_compare(3, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {

        if (   PDM_morton_compare(3, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;

          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 8; i++) {
    _node_t  n1 = next_bt->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
         == (  next_bt->nodes[n_init_nodes + 8 - 1].start_id
             + next_bt->nodes[n_init_nodes + 8 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

static void
_split_node_2d(PDM_box_tree_t       *bt,
               PDM_box_tree_t       *next_bt,
               const PDM_box_set_t  *boxes,
               int             node_id,
               int            *shift_ids)
{
  int j, i;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[4];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->n_nodes;
  _node_t  split_node = next_bt->nodes[node_id];

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 4);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 4 > next_bt->n_max_nodes) {
    assert(next_bt->n_max_nodes > 0);
    next_bt->n_max_nodes *= 2;
    next_bt->nodes = (_node_t *) realloc((void *) next_bt->nodes, next_bt->n_max_nodes * sizeof(_node_t));
    next_bt->child_ids = (int *) realloc((void *) next_bt->child_ids, next_bt->n_max_nodes*4 * sizeof(int));
  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(2, node.morton_code, children);

  for (i = 0; i < 4; i++) {
    const int   new_id = n_init_nodes + i;
    next_bt->child_ids[node_id*4 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;

  next_bt->nodes[node_id] = split_node;
  next_bt->n_nodes = n_init_nodes + 4;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(2, next_level, box_min);
    max_code = PDM_morton_encode(2, next_level, box_max);

    if (   PDM_morton_compare(2, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {
        if (_node_intersect_box_2d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(2, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {
        if (   PDM_morton_compare(2, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 4; i++) {

    (next_bt->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 4; i++)
    (next_bt->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(2, next_level, box_min);
    max_code = PDM_morton_encode(2, next_level, box_max);

    if (   PDM_morton_compare(2, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {

        if (_node_intersect_box_2d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;
          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(2, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {

        if (   PDM_morton_compare(2, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;

          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 4; i++) {
    _node_t  n1 = next_bt->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
         == (  next_bt->nodes[n_init_nodes + 4 - 1].start_id
             + next_bt->nodes[n_init_nodes + 4 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

static void
_split_node_1d(PDM_box_tree_t       *bt,
               PDM_box_tree_t       *next_bt,
               const PDM_box_set_t  *boxes,
               int             node_id,
               int            *shift_ids)
{
  int j, i;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[2];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->n_nodes;
  _node_t  split_node = next_bt->nodes[node_id];

  const _node_t  node = bt->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 2);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 2 > next_bt->n_max_nodes) {
    assert(next_bt->n_max_nodes > 0);
    next_bt->n_max_nodes *= 2;
    next_bt->nodes = (_node_t *) realloc((void *) next_bt->nodes, next_bt->n_max_nodes * sizeof(_node_t));
    next_bt->child_ids = (int *) realloc((void *) next_bt->child_ids, next_bt->n_max_nodes*2 * sizeof(int));
  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(1, node.morton_code, children);

  for (i = 0; i < 2; i++) {

    const int   new_id = n_init_nodes + i;
    next_bt->child_ids[node_id*2 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;

  next_bt->nodes[node_id] = split_node;
  next_bt->n_nodes = n_init_nodes + 2;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(1, next_level, box_min);
    max_code = PDM_morton_encode(1, next_level, box_max);

    if (   PDM_morton_compare(1, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {
        if (_node_intersect_box_1d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same segment */

      assert(   PDM_morton_compare(1, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {
        if (   PDM_morton_compare(1, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 2; i++) {

    (next_bt->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 2; i++)
    (next_bt->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(1, next_level, box_min);
    max_code = PDM_morton_encode(1, next_level, box_max);

    if (   PDM_morton_compare(1, min_code, max_code)
        == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {

        if (_node_intersect_box_1d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;
          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same segment */

      assert(   PDM_morton_compare(1, max_code, min_code)
             == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {

        if (   PDM_morton_compare(1, min_code, children[i])
            == PDM_MORTON_EQUAL_ID) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;

          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 2; i++) {
    _node_t  n1 = next_bt->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
         == (  next_bt->nodes[n_init_nodes + 2 - 1].start_id
             + next_bt->nodes[n_init_nodes + 2 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

/*----------------------------------------------------------------------------
 * Evaluate the box distribution over the leaves of the box tree when adding
 * a level to the tree structure.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   next_bt    <->  pointer to the next box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   node_id    <--  id of the starting node
 *   build_type <--  layout variant for building the tree structure
 *   shift_ids  <->  first free position free in new box_ids
 *----------------------------------------------------------------------------*/

static void
_build_next_level(PDM_box_tree_t       *bt,
                  PDM_box_tree_t       *next_bt,
                  const PDM_box_set_t  *boxes,
                  int             node_id,
                  PDM_box_tree_sync_t   build_type,
                  int            *shift_ids)
{
  int   i;

  int   _shift_ids = *shift_ids;
  const _node_t  *cur_node = bt->nodes + node_id;

  if (cur_node->is_leaf == false) {

    assert(bt->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _build_next_level(bt,
                        next_bt,
                        boxes,
                        bt->child_ids[bt->n_children*node_id + i],
                        build_type,
                        &_shift_ids);
  }

  else { /* if (node->is_leaf == true) */

    if (   cur_node->n_boxes < bt->threshold
        && node_id != 0                    /* Root node is always divided */
        && build_type == PDM_BOX_TREE_ASYNC_LEVEL) {

      /* Copy related box_ids in the new next_ids */

      _node_t *next_node = next_bt->nodes + node_id;

      next_node->n_boxes = cur_node->n_boxes;
      next_node->start_id = _shift_ids;

      for (i = 0; i < cur_node->n_boxes; i++)
        next_bt->box_ids[_shift_ids++]
          = bt->box_ids[cur_node->start_id + i];
    }
    else {  /* Split node and evaluate box distribution between its children */

      if (boxes->dim == 3)
        _split_node_3d(bt,
                       next_bt,
                       boxes,
                       node_id,
                       &_shift_ids);
      else if (boxes->dim == 2)
        _split_node_2d(bt,
                       next_bt,
                       boxes,
                       node_id,
                       &_shift_ids);
      else if (boxes->dim == 1)
        _split_node_1d(bt,
                       next_bt,
                       boxes,
                       node_id,
                       &_shift_ids);
    }

  }

  /* Prepare return values */

  *shift_ids = _shift_ids;
}

/*----------------------------------------------------------------------------
 * Loop on all nodes of the box tree to define an array with Morton codes
 * and weights (= number of linked boxes) associated to each leaf
 *
 * parameters:
 *   bt         <-> pointer to PDM_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   node_id    <-- id of the current node (to traverse)
 *   n_leaves   <-> current number of leaves in the tree with n_boxes > 0
 *   leaf_codes <-> Morton code associated to each leaf
 *   weight     <-> number of boxes attached to each leaf
 *----------------------------------------------------------------------------*/

static void
_build_leaf_weight(const PDM_box_tree_t  *bt,
                   int              node_id,
                   int             *n_leaves,
                   PDM_morton_code_t     *leaf_codes,
                   int             *weight)
{
  int  i;

  int _n_leaves = *n_leaves;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false)
    for (i = 0; i < bt->n_children; i++)
      _build_leaf_weight(bt,
                         bt->child_ids[bt->n_children*node_id + i],
                         &_n_leaves,
                         leaf_codes,
                         weight);

  else { /* node is a leaf */

    if (node->n_boxes > 0) {
      leaf_codes[_n_leaves] = node->morton_code;
      weight[_n_leaves] = node->n_boxes;
      _n_leaves += 1;
    }
  }

  *n_leaves = _n_leaves;
}

/*----------------------------------------------------------------------------
 * Loop on all nodes of the tree to define an index of ranks related
 * to each box.
 *
 * parameters:
 *   bt           <-> pointer to PDM_box_tree_t structure.
 *   distrib      <-- structure holding box distribution data
 *   dim          <-- box tree layout dimension (1, 2, or 3)
 *   node_id      <-- id of the current node (to traverse)
 *   size         <-- size of index in which we search
 *   search_index <-- index on which box distribution is made
 *   id_rank      <-- relation between id and rank
 *----------------------------------------------------------------------------*/

static void
_build_rank_to_box_index(const PDM_box_tree_t  *bt,
                         PDM_box_distrib_t     *distrib,
                         int                    dim,
                         int              node_id,
                         size_t                 size,
                         PDM_morton_code_t      search_index[],
                         int                    id_rank[])
{
  int  i;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _build_rank_to_box_index(bt,
                               distrib,
                               dim,
                               bt->child_ids[bt->n_children*node_id + i],
                               size,
                               search_index,
                               id_rank);
  }
  else {

    if (node->n_boxes > 0) {

      int  id = PDM_morton_binary_search((int) size,
                                         node->morton_code,
                                         search_index);
      int  rank = id_rank[id];

      distrib->index[rank + 1] += node->n_boxes;

    }
  }

}

/*----------------------------------------------------------------------------
 * Loop on all nodes of the tree to define a list of ranks related
 * to each box.
 *
 * parameters:
 *   bt           <-> pointer to PDM_box_tree_t structure.
 *   distrib      <-- structure holding box distribution data
 *   dim          <-- box tree layout dimension (1, 2, or 3)
 *   node_id      <-- id of the current node (to traverse)
 *   counter      <-> counter array used to build the list
 *   size         <-- size of index in which we search
 *   search_index <-- index on which box distribution is made
 *   id_rank      <-- relation between id and rank
 *----------------------------------------------------------------------------*/

static void
_build_rank_to_box_list(const PDM_box_tree_t  *bt,
                        PDM_box_distrib_t     *distrib,
                        int                    dim,
                        int              node_id,
                        int              counter[],
                        size_t                 size,
                        PDM_morton_code_t      search_index[],
                        int                    id_rank[])
{
  int  i;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _build_rank_to_box_list(bt,
                              distrib,
                              dim,
                              bt->child_ids[bt->n_children*node_id + i],
                              counter,
                              size,
                              search_index,
                              id_rank);
  }
  else {

    if (node->n_boxes > 0) {

      int  id = PDM_morton_binary_search((int) size,
                                         node->morton_code,
                                         search_index);
      int  rank = id_rank[id];

      for (i = 0; i < node->n_boxes; i++) {

        int   box_id = bt->box_ids[node->start_id + i];
        int   shift = distrib->index[rank] + counter[rank];

        distrib->list[shift] = box_id;
        counter[rank] += 1;

      }
    }
  }

}


/*----------------------------------------------------------------------------
 * Recursively build an index on boxes which intersect.
 *
 * parameters:
 *   bt      <-> pointer to PDM_box_tree_t structure.
 *   boxesB  <-- pointer boxes that intersect associated tree boxes
 *   node_id <-- id of the current node (to traverse)
 *   count   <-> intersection count
 *----------------------------------------------------------------------------*/

static void
_count_boxes_intersections(const PDM_box_tree_t  *bt,
                           const PDM_box_set_t  *boxesB,
                           int                   boxB_id,
                           int                   node_id,
                           int                   count[])
{
  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->extents;
  const _node_t  *node = bt->nodes + node_id;

  const double *boxB_extents = boxesB->extents + 2 * boxes->dim * boxB_id;

  assert (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 node_id));

  if (node->is_leaf == false) {

    for (int i = 0; i < bt->n_children; i++) { /* traverse downwards */

      if (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 bt->child_ids[bt->n_children*node_id + i])) {

        _count_boxes_intersections(bt,
                                   boxesB,
                                   boxB_id,
                                   bt->child_ids[bt->n_children*node_id + i],
                                   count);
      }
    }
  }

  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = bt->box_ids[node->start_id + i];
        if (_boxes_intersect_3d (box_extents, id0, boxB_extents)) {
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 2) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = bt->box_ids[node->start_id + i];
        if (_boxes_intersect_2d (box_extents, id0, boxB_extents)) {
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 1) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = bt->box_ids[node->start_id + i];
        if (_boxes_intersect_1d(box_extents, id0, boxB_extents)) {
          count[id0] += 1;
        }
      }

    }
  }
}


/*----------------------------------------------------------------------------
 * Recursively build an index on boxes which intersect.
 *
 * parameters:
 *   bt      <-> pointer to PDM_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *   count   <-> intersection count
 *----------------------------------------------------------------------------*/

static void
_count_intern_intersections(const PDM_box_tree_t  *bt,
                            int                   node_id,
                            int                   count[])
{
  int   i, j;

  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->extents;
  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _count_intern_intersections(bt,
                                  bt->child_ids[bt->n_children*node_id + i],
                                  count);
  }
  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = bt->box_ids[node->start_id + i];
          int   id1 = bt->box_ids[node->start_id + j];
          if (_boxes_intern_intersect_3d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

    else if (boxes->dim == 2) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = bt->box_ids[node->start_id + i];
          int   id1 = bt->box_ids[node->start_id + j];
          if (_boxes_intern_intersect_2d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

    else if (boxes->dim == 1) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = bt->box_ids[node->start_id + i];
          int   id1 = bt->box_ids[node->start_id + j];
          if (_boxes_intern_intersect_1d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

  }

}


/*----------------------------------------------------------------------------
 * Recursively build a list on bounding boxes which intersect together.
 *
 * parameters:
 *   bt        <-> pointer to PDM_box_tree_t structure.
 *   boxesB    <-- pointer boxes that intersect associated tree boxes
 *   node_id   <-- id of the current node (to traverse)
 *   count     <-> intersection count (workin array)
 *   index     <-- index on intersections
 *   box_g_num <-> global number of intersection boxes
 *----------------------------------------------------------------------------*/

static void
_get_boxes_intersections(const PDM_box_tree_t  *bt,
                         const PDM_box_set_t  *boxesB,
                         int                   boxB_id,
                         int                   node_id,
                         int                   count[],
                         int                   box_index[],
                         int                   box_l_num[])
{
  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->extents;
  const _node_t  *node = bt->nodes + node_id;

  const double *boxB_extents = boxesB->extents + 2 * boxes->dim * boxB_id;

  assert (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 node_id));

  if (node->is_leaf == false) {

    for (int i = 0; i < bt->n_children; i++) { /* traverse downwards */

      if (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 bt->child_ids[bt->n_children*node_id + i])) {

        _get_boxes_intersections (bt,
                                  boxesB,
                                  boxB_id,
                                  bt->child_ids[bt->n_children*node_id + i],
                                  count,
                                  box_index,
                                  box_l_num);
      }
    }
  }

  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = bt->box_ids[node->start_id + i];
        if (_boxes_intersect_3d (box_extents, id0, boxB_extents)) {
          int   shift0 = box_index[id0] + count[id0];
          box_l_num[shift0] = boxB_id;
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 2) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = bt->box_ids[node->start_id + i];
        if (_boxes_intersect_2d (box_extents, id0, boxB_extents)) {
          int   shift0 = box_index[id0] + count[id0];
          box_l_num[shift0] = boxB_id;
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 1) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = bt->box_ids[node->start_id + i];
        if (_boxes_intersect_1d (box_extents, id0, boxB_extents)) {
          int   shift0 = box_index[id0] + count[id0];
          box_l_num[shift0] = boxB_id;
          count[id0] += 1;
        }
      }

    }
  } /* End if node is a leaf */
}


/*----------------------------------------------------------------------------
 * Recursively build a list on bounding boxes which intersect together.
 *
 * parameters:
 *   bt        <-> pointer to PDM_box_tree_t structure.
 *   node_id   <-- id of the current node (to traverse)
 *   count     <-> intersection count (workin array)
 *   index     <-- index on intersections
 *   box_g_num <-> global number of intersection boxes
 *----------------------------------------------------------------------------*/

static void
_get_intern_intersections(const PDM_box_tree_t  *bt,
                          int                    node_id,
                          int                    count[],
                          int                    box_index[],
                          PDM_g_num_t             box_g_num[])
{
  int  i, j;

  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->extents;
  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _get_intern_intersections(bt,
                                bt->child_ids[bt->n_children*node_id + i],
                                count,
                                box_index,
                                box_g_num);
  }
  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = bt->box_ids[node->start_id + i];
          int   id1 = bt->box_ids[node->start_id + j];

          if (_boxes_intern_intersect_3d(box_extents, id0, id1)) {
            int   shift0 = box_index[id0] + count[id0];
            int   shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->g_num[id1];
            box_g_num[shift1] = boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }
    else if (boxes->dim == 2) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = bt->box_ids[node->start_id + i];
          int   id1 = bt->box_ids[node->start_id + j];

          if (_boxes_intern_intersect_2d(box_extents, id0, id1)) {
            int   shift0 = box_index[id0] + count[id0];
            int   shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->g_num[id1];
            box_g_num[shift1] = boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }

    else if (boxes->dim == 1) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = bt->box_ids[node->start_id + i];
          int   id1 = bt->box_ids[node->start_id + j];

          if (_boxes_intern_intersect_1d(box_extents, id0, id1)) {
            int   shift0 = box_index[id0] + count[id0];
            int   shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->g_num[id1];
            box_g_num[shift1] = boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }

  } /* End if node is a leaf */
}

/*----------------------------------------------------------------------------
 * Recursively define a counter array on the number of bounding boxes
 * associated to a leaf.
 *
 * This will be used for displaying a histogram.
 *
 * parameters:
 *   bt      <-- pointer to PDM_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *   n_steps <-- number of steps in histogram
 *   step    <-- steps of the histogram
 *   h_min   <-- min. value of the histogram
 *   counter <-> counter (working array)
 *----------------------------------------------------------------------------*/

static void
_build_histogram(const PDM_box_tree_t  *bt,
                 int              node_id,
                 int              n_steps,
                 int              step,
                 int              h_min,
                 PDM_g_num_t              count[])
{
  int  i, j;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _build_histogram(bt,
                       bt->child_ids[bt->n_children*node_id + i],
                       n_steps,
                       step,
                       h_min,
                       count);
  }
  else {
    for (i = 0, j = 1; j < n_steps; i++, j++)
      if (node->n_boxes < h_min + j*step)
        break;
    count[i] += 1;
  }
}

/*----------------------------------------------------------------------------
 * Dump a box tree node.
 *
 * parameters:
 *   bt      <-- pointer to PDM_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *----------------------------------------------------------------------------*/

static void
_dump_node(const PDM_box_tree_t  *bt,
           int              node_id)
{
  int  i;

  const char *node_type[] = {"node", "leaf"};

  const _node_t  *node = bt->nodes + node_id;
  const PDM_morton_code_t  m_code = node->morton_code;

  PDM_printf("\n"
             "  node %10d (%s)\n"
             "    level:   %3u - anchor: [ %10u %10u %10u ]\n"
             "    n_boxes: %3d - start_id: %u\n"
             "    boxes:\n",
             node_id, node_type[(int)(node->is_leaf)],
	      m_code.L, m_code.X[0], m_code.X[1], m_code.X[2],
             node->n_boxes, node->start_id);

  for (i = 0; i < node->n_boxes; i++)
    PDM_printf("        %d\n", (int)(bt->box_ids[node->start_id + i]));

  if (node->is_leaf == false) {

    const int *c_id = bt->child_ids + bt->n_children*node_id;

    if (bt->n_children == 8) {
      PDM_printf("  children_id:  %d %d %d %d %d %d %d %d\n",
                 (int)c_id[0], (int)c_id[1], (int)c_id[2], (int)c_id[3],
                 (int)c_id[4], (int)c_id[5], (int)c_id[6], (int)c_id[7]);
    }
    else if (bt->n_children == 4) {
      PDM_printf("  children_id:  %d %d %d %d\n",
                 (int)c_id[0], (int)c_id[1], (int)c_id[2], (int)c_id[3]);
    }
    else if (bt->n_children == 2) {
      PDM_printf("  children_id:  %d %d\n",
                 (int)c_id[0], (int)c_id[1]);
    }

    for (i = 0; i < bt->n_children; i++)
      _dump_node(bt, c_id[i]);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a PDM_box_tree_t structure and initialize it.
 *
 * parameters:
 *  max_level     <-- max possible level
 *  threshold     <-- max number of  boxes linked to an octant if
 *                    max_level is not reached
 *  max_box_ratio <-- max n_linked_boxes / n_boxes ratio
 *
 * returns:
 *   pointer to an empty PDM_box_tree_t structure.
 *----------------------------------------------------------------------------*/

PDM_box_tree_t *
PDM_box_tree_create(int    max_level,
                    int    threshold,
                    float  max_box_ratio)
{
  PDM_box_tree_t  *bt = NULL;

  bt = (PDM_box_tree_t *) malloc(sizeof(PDM_box_tree_t));

  /* Sanity checks */

  if (max_level < 0) {
    PDM_error(__FILE__, __LINE__, 0,
            "  Forbidden max_level value (%d) in the tree structure\n",
            max_level);
    abort();
  }

  if (threshold < 1) {
    PDM_error(__FILE__, __LINE__, 0,
            "  Forbidden threshold value (%d) in the tree structure\n",
           threshold);
    abort();
  }

  if (max_box_ratio < 1.0) {
    PDM_error(__FILE__, __LINE__, 0,
            "  Forbidden max_box_ratio value (%f) in the tree structure\n",
            (double)max_box_ratio);
    abort();
  }

  /* Create and initialize tree structure according to its type */

  bt->max_level = max_level;
  bt->threshold = threshold;
  bt->max_box_ratio = max_box_ratio;

  bt->comm = PDM_MPI_COMM_NULL;

  /* Set stats */

  bt->stats.max_level_reached = 0;

  bt->stats.n_leaves = 0;
  bt->stats.n_spill_leaves = 0;
  bt->stats.n_linked_boxes = 0;

  bt->stats.min_linked_boxes = INT_MAX;
  bt->stats.max_linked_boxes = 0;

  /* Initialize nodes */

  bt->n_max_nodes = 0;
  bt->n_nodes = 0;

  bt->nodes = NULL;

  bt->box_ids = NULL;

  bt->n_build_loops = 0;

  bt->stack = NULL;
  bt->pos_stack = NULL;

  return bt;
}

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to PDM_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_destroy(PDM_box_tree_t  **bt)
{
  PDM_box_tree_t  *_bt = *bt;

  if (_bt != NULL) {

    free(_bt->nodes);
    free(_bt->child_ids);
    free(_bt->box_ids);

    if (_bt->stack != NULL) {
      free (_bt->stack);
    }
    if (_bt->pos_stack != NULL) {
      free (_bt->pos_stack);
    }

    free(_bt);
    *bt = _bt;
  }
}

/*----------------------------------------------------------------------------
 * Get the deepest level allowed by the tree structure.
 *
 * parameters:
 *   bt <-- pointer to PDM_box_tree_t structure.
 *
 * returns:
 *   deepest allowed level of the tree
 *----------------------------------------------------------------------------*/

int
PDM_box_tree_get_max_level(const PDM_box_tree_t  *bt)
{
  return bt->max_level;
}

/*----------------------------------------------------------------------------
 * Assign a set of boxes to an empty PDM_box_tree_t structure.
 *
 * The box tree structure must have been created using to PDM_tree_create().
 *
 * The depth of the tree is adjusted so that a maximum of max_n_elts boxes
 * will be assigned to each leaf, unless this would require going beyond
 * the tree's maximum level.
 *
 * If max_level = -1, the highest level reachable is PDM_TREE_MAX_LEVEL but
 * there is no defined target level.
 *
 * parameters:
 *   bt         <-> pointer to PDM_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   build_type <-- layout variant for building the tree structure
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_set_boxes(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       PDM_box_tree_sync_t   build_type)
{
  int   box_id;

  PDM_box_tree_t  tmp_bt;

  int   next_box_ids_size = 0, shift = 0;
  double anchor[3] = {0., 0., 0.};

  /* Initialization */

  assert(bt != NULL);

  bt->n_build_loops = 0;

  bt->comm = boxes->comm;
  bt->boxes = boxes;

  /* Preallocate for the two first levels of a tree */

  if (boxes->dim == 3) {
    bt->n_children = 8;
    bt->n_max_nodes = 73;
  }
  else if (boxes->dim == 2) {
    bt->n_children = 4;
    bt->n_max_nodes = 21;
  }
  else if (boxes->dim == 1) {
    bt->n_children = 2;
    bt->n_max_nodes = 7;
  }

  bt->n_nodes = 1;

  bt->nodes = (_node_t *) malloc(bt->n_max_nodes * sizeof(_node_t));
  bt->child_ids = (int *) malloc(bt->n_max_nodes*bt->n_children * sizeof(int));

  /* Define root node */

  _new_node(bt, PDM_morton_encode(boxes->dim, 0, anchor), 0);

  /* Initialize bt by assigning all boxes to the root leaf */

  bt->box_ids = (int *) malloc (boxes->n_boxes * sizeof(int));

  for (box_id = 0; box_id < boxes->n_boxes; box_id++)
    bt->box_ids[box_id] = box_id;

  (bt->nodes[0]).is_leaf = true;
  (bt->nodes[0]).n_boxes = boxes->n_boxes;
  (bt->nodes[0]).start_id = 0;

  bt->stats.n_boxes = boxes->n_boxes;

  _get_box_tree_stats(bt);

  /* Build local tree structure by adding boxes from the root */

  while (_recurse_tree_build(bt,
                             boxes,
                             build_type,
                             &next_box_ids_size)) {

    /* Initialize next_bt: copy of bt */

    _copy_tree(&tmp_bt, bt);

    /* Optimize memory usage */

    bt->n_max_nodes = bt->n_nodes;
    bt->nodes = (_node_t *) realloc((void *) bt->nodes, bt->n_nodes * sizeof(_node_t));
    bt->child_ids = (int *) realloc((void *) bt->child_ids,
                                    bt->n_max_nodes*bt->n_children * sizeof(int));

    /* Define a box ids list for the next level of the boxtree */

    tmp_bt.box_ids = (int*) realloc((void *) tmp_bt.box_ids, next_box_ids_size * sizeof(int));
    shift = 0;

    _build_next_level(bt,
                      &tmp_bt,
                      boxes,
                      0, /* Starts from root */
                      build_type,
                      &shift);

    assert(shift == next_box_ids_size);

    /* replace current tree by the tree computed at a higher level */

    _free_tree_arrays(bt);
    *bt = tmp_bt; /* Overwrite bt members with those of next_bt */

    _get_box_tree_stats(bt);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    PDM_printf("  - New box tree level -\n");
    PDM_box_tree_dump_statistics(bt);
#endif

  } /* While building should continue */
}

/*----------------------------------------------------------------------------
 * Compute an index based on Morton encoding to ensure a good distribution
 * of boxes among the participating ranks.
 *
 * parameters:
 *   bt         <-> pointer to PDM_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *
 * returns:
 *   pointer to newly created PDM_box_distrib_t structure.
 *----------------------------------------------------------------------------*/

PDM_box_distrib_t *
PDM_box_tree_get_distrib(PDM_box_tree_t        *bt,
                         const PDM_box_set_t   *boxes)
{
  int  i;

  int  reduce_size = 0;
  int   n_leaves = 0;
  int  *reduce_ids = NULL;
  PDM_morton_code_t  *leaf_codes = NULL, *reduce_index = NULL;
  int   *weight = NULL, *counter = NULL;

  PDM_box_distrib_t  *distrib = NULL;

  assert(bt != NULL);
  assert(boxes != NULL);

  /* Compute basic box distribution */

  distrib = PDM_box_distrib_create(boxes->n_boxes,
                                   boxes->n_g_boxes,
                                   (bt->stats).max_level_reached,
                                   boxes->comm);

  if (distrib == NULL)
    return NULL;

  leaf_codes = (PDM_morton_code_t *) malloc(bt->stats.n_leaves * sizeof(PDM_morton_code_t));
  weight = (int *) malloc(bt->stats.n_leaves * sizeof(int));

  /* Build index for boxes */

  _build_leaf_weight(bt,
                     0,
                     &n_leaves,
                     leaf_codes,
                     weight);

  assert(n_leaves <= bt->stats.n_leaves);

  leaf_codes = (PDM_morton_code_t *) realloc((void *) leaf_codes, n_leaves * sizeof(PDM_morton_code_t));
  weight = (int *) realloc((void *) weight, n_leaves * sizeof(int));

  /* Compute the resulting Morton index */

  PDM_box_set_build_morton_index(boxes,
                                 distrib,
                                 n_leaves,
                                 leaf_codes,
                                 weight);

  free(leaf_codes);
  free(weight);

  /* Compact Morton_index to get an array without "0 element" */

  for (i = 0; i < distrib->n_ranks; i++)
    if (PDM_morton_a_gt_b(distrib->morton_index[i+1],
                          distrib->morton_index[i]))
      reduce_size++;

  reduce_index = (PDM_morton_code_t *) malloc((reduce_size + 1) * sizeof(PDM_morton_code_t));
  reduce_ids  = (int *) malloc(reduce_size * sizeof(int));

  reduce_size = 0;
  reduce_index[0] = distrib->morton_index[0];

  for (i = 0; i < distrib->n_ranks; i++) {

    if (PDM_morton_a_gt_b(distrib->morton_index[i+1],
                          distrib->morton_index[i])) {

      reduce_index[reduce_size + 1] = distrib->morton_index[i+1];
      reduce_ids[reduce_size++] = i;

    }

  }

  /* Define a rank -> box indexed list */

  _build_rank_to_box_index(bt,
                           distrib,
                           boxes->dim,
                           0,  /* starts from root */
                           reduce_size,
                           reduce_index,
                           reduce_ids);

  for (i = 0; i < distrib->n_ranks; i++)
    distrib->index[i+1] += distrib->index[i];

  distrib->list = (int *) malloc(distrib->index[distrib->n_ranks] * sizeof(int));

  counter = (int *) malloc(distrib->n_ranks * sizeof(int));

  for (i = 0; i < distrib->n_ranks; i++)
    counter[i] = 0;

  _build_rank_to_box_list(bt,
                          distrib,
                          boxes->dim,
                          0,  /* starts from root */
                          counter,
                          reduce_size,
                          reduce_index,
                          reduce_ids);

  /* Free memory */

  free(counter);
  free(reduce_ids);
  free(reduce_index);

  /* Define the final index (without redundancies) and realloc list */

  PDM_box_distrib_clean(distrib);

  return distrib;
}


/*----------------------------------------------------------------------------
 * Build an indexed list on associated bt boxes to list boxes B intersections.
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * parameters:
 *   bt        <-- pointer to box tree structure to query
 *   boxesB    <-- pointer boxes that intersect associated tree boxes
 *   box_index --> pointer to the index array on bounding boxes
 *   box_g_num --> pointer to the list of intersecting bounding boxesB
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_get_boxes_intersects(PDM_box_tree_t       *bt,
                                  const PDM_box_set_t  *boxesB,
                                  int                  *box_index[],
                                  int                  *box_l_num[])
{
  int  i, list_size;

  int  *counter = NULL;
  int  *_index = NULL;
  int         *_l_num = NULL;
  const PDM_box_set_t  *boxes = bt->boxes;


  /* Build index */

  _index = (int *) malloc((boxes->n_boxes + 1) * sizeof(int));

  for (i = 0; i < boxes->n_boxes + 1; i++)
    _index[i] = 0;

  for (int k = 0; k < boxesB->n_boxes; k++) {

    _count_boxes_intersections(bt,
                               boxesB,
                               k,
                               0, /* start from root */
                               _index + 1);
  }

  /* Build index from counts */

  for (i = 0; i < boxes->n_boxes; i++)
    _index[i+1] += _index[i];

  list_size = _index[boxes->n_boxes];

  _l_num = (int *) malloc(list_size * sizeof(int));

  counter = (int *) malloc(boxes->n_boxes * sizeof(int));

  for (i = 0; i < boxes->n_boxes; i++)
    counter[i] = 0;

  /* Build list */

  for (int k = 0; k < boxesB->n_boxes; k++) {
    _get_boxes_intersections(bt,
                             boxesB,
                             k,
                             0, /* start from root */
                             counter,
                             _index,
                             _l_num);
  }

  free(counter);

  /* Return pointers */

  *box_index = _index;
  *box_l_num = _l_num;
}


/*----------------------------------------------------------------------------
 * Build an indexed list on boxes to list intersections.
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * parameters:
 *   bt        <-- pointer to box tree structure to query
 *   box_index --> pointer to the index array on bounding boxes
 *   box_g_num --> pointer to the list of intersecting bounding boxes
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_get_intern_intersects(PDM_box_tree_t       *bt,
                                   int                  *box_index[],
                                   PDM_g_num_t           *box_g_num[])
{
  int  i, list_size;

  int  *counter = NULL;
  int  *_index = NULL;
  PDM_g_num_t  *_g_num = NULL;
  const PDM_box_set_t  *boxes = bt->boxes;

  /* Build index */

  _index = (int *) malloc((boxes->n_boxes + 1) * sizeof(int));

  for (i = 0; i < boxes->n_boxes + 1; i++)
    _index[i] = 0;

  _count_intern_intersections(bt,
                              0, /* start from root */
                              _index + 1);

  /* Build index from counts */

  for (i = 0; i < boxes->n_boxes; i++)
    _index[i+1] += _index[i];

  list_size = _index[boxes->n_boxes];

  _g_num = (PDM_g_num_t *) malloc(list_size * sizeof(PDM_g_num_t));

  counter = (int *) malloc(boxes->n_boxes * sizeof(int));

  for (i = 0; i < boxes->n_boxes; i++)
    counter[i] = 0;

  /* Build list */

  _get_intern_intersections(bt,
                            0, /* start from root */
                            counter,
                            _index,
                            _g_num);

  free(counter);

  /* Return pointers */

  *box_index = _index;
  *box_g_num = _g_num;
}

/*----------------------------------------------------------------------------
 * Get global box tree statistics.
 *
 * All fields returned are optional: if their argument is set to NULL,
 * the corresponding information will not be returned.
 *
 * For each field not set to NULL, 3 values are always returned:
 * the mean on all ranks (rounded to the closest integer), the minimum,
 * and the maximum value respectively.
 *
 * In serial mode, the mean, minimum, and maximum will be identical for most
 * fields, but all 3 values are returned nonetheless.
 *
 * Note that the theoretical memory use includes that of the associated
 * box set.
 *
 * parameters:
 *   bt                 <-- pointer to box tree structure
 *   depth              --> tree depth (max level used)
 *   n_leaves           --> number of leaves in the tree
 *   n_boxes            --> number of boxes in the tree
 *   n_threshold_leaves --> number of leaves where n_boxes > threshold
 *   n_leaf_boxes       --> number of boxes for a leaf
 *   mem_used           --> theoretical used memory
 *   mem_allocated      --> theoretical allocated memory
 *
 * returns:
 *   the spatial dimension associated with the box tree layout (3, 2, or 1)
 *----------------------------------------------------------------------------*/

int
PDM_box_tree_get_stats(const PDM_box_tree_t  *bt,
                       int                    depth[3],
                       int              n_leaves[3],
                       int              n_boxes[3],
                       int              n_threshold_leaves[3],
                       int              n_leaf_boxes[3],
                       size_t                 mem_used[3],
                       size_t                 mem_allocated[3])
{
  int i;
  uint64_t mem_per_node;
  uint64_t s_mean[7], s_min[7], s_max[7];
  PDM_box_tree_stats_t s;

  int dim = 3;

  if (bt == NULL)
    return 0;

  s = bt->stats;

  if (bt->n_children == 4)
    dim = 2;
  else if (bt->n_children == 2)
    dim = 1;

  /* Prepare array of local values; prior to or in the absence of
     MPI communication, mean values are set to local values. */

  s_mean[0] = s.n_linked_boxes / s.n_leaves;
  /* Round to nearest integer, and not floor */
  if (s.n_linked_boxes % s.n_leaves >= s.n_leaves/2)
    s_mean[0] += 1;

  s_min[0] = s.min_linked_boxes;
  s_max[0] = s.max_linked_boxes;

  s_mean[1] = s.max_level_reached;
  s_mean[2] = s.n_leaves;
  s_mean[3] = s.n_boxes;
  s_mean[4] = s.n_spill_leaves;

  /* Estimate theoretical memory usage */

  mem_per_node = sizeof(_node_t) + bt->n_children*sizeof(int);

  s_mean[5] = sizeof(PDM_box_tree_t);
  s_mean[5] += bt->n_nodes * mem_per_node;
  s_mean[5] += s.n_linked_boxes * sizeof(int);

  s_mean[5] += sizeof(PDM_box_set_t);
  s_mean[5] += s.n_boxes * (  sizeof(PDM_g_num_t)
                           + (dim * 2 * sizeof(double)));

  s_mean[6] = s_mean[5] + (bt->n_max_nodes - bt->n_nodes)*mem_per_node;

 /* Pre-synchronize for serial cases (i = 0 already handled) */

  for (i = 1; i < 7; i++) {
    s_min[i] = s_mean[i];
    s_max[i] = s_mean[i];
  }

  /* In parallel mode, synchronize values */

  if (bt->comm != PDM_MPI_COMM_NULL) {

    int n_ranks;
    uint64_t s_l_sum[14], s_g_sum[14];

    PDM_MPI_Comm_size(bt->comm, &n_ranks);

    if (n_ranks > 1) { /* Should always be the case, bat play it safe) */

      /* Split value to avoid exceeding PDM_g_num_t limits
         (especially when it is a 32-bit value) */

      s_l_sum[0] = s.n_linked_boxes/n_ranks;
      s_l_sum[7] = s.n_linked_boxes%n_ranks;
      for (i = 1; i < 7; i++) {
        s_l_sum[i] = ((uint64_t) s_mean[i])/n_ranks;
        s_l_sum[i+7] = ((uint64_t) s_mean[i])%n_ranks;
      }

      PDM_MPI_Allreduce(s_l_sum, s_g_sum, 14, PDM_MPI_UINT64_T, PDM_MPI_SUM, bt->comm);

      s_mean[0] = s.min_linked_boxes;
      PDM_MPI_Allreduce(s_mean, s_min, 7, PDM_MPI_UINT64_T, PDM_MPI_MIN, bt->comm);
      s_mean[0] = s.max_linked_boxes;
      PDM_MPI_Allreduce(s_mean, s_max, 7, PDM_MPI_UINT64_T, PDM_MPI_MAX, bt->comm);

      /* Specific handling for linked boxes, so as to ensure correct
         total using large integers even if we do not know the
         corresponding MPI type */
      {
        uint64_t s_n = s_g_sum[0]*n_ranks + s_g_sum[7]; /* linked boxes */
        uint64_t s_d = s_g_sum[2]*n_ranks + s_g_sum[9]; /* leaves */
        uint64_t s_m = s_n / s_d;
        /* Round to nearest integer, and not floor */
        if (s_n % s_d >= s_d/2)
          s_m += 1;
        s_mean[0] = s_m;
      }
      for (i = 1; i < 7; i++) {
        s_mean[i] = s_g_sum[i] + s_g_sum[i+7]/n_ranks;
        /* Round to nearest integer, and not floor */
        if (s_g_sum[i+7]%n_ranks >= (PDM_g_num_t)n_ranks/2)
          s_mean[i] += 1;
      }
    }

  }

  /* Set values already in stats */

  if (depth != NULL) {
    depth[0] = (int) s_mean[1];
    depth[1] = (int) s_min[1];
    depth[2] = (int) s_max[1];
  }

  if (n_leaves != NULL) {
    n_leaves[0] = (int) s_mean[2];
    n_leaves[1] = (int) s_min[2];
    n_leaves[2] = (int) s_max[2];
  }

  if (n_boxes != NULL) {
    n_boxes[0] = (int) s_mean[3];
    n_boxes[1] = (int) s_min[3];
    n_boxes[2] = (int) s_max[3];
  }

  if (n_threshold_leaves != NULL) {
    n_threshold_leaves[0] = (int) s_mean[4];
    n_threshold_leaves[1] = (int) s_min[4];
    n_threshold_leaves[2] = (int) s_max[4];
  }

  if (n_leaf_boxes != NULL) {
    n_leaf_boxes[0] = (int) s_mean[0];
    n_leaf_boxes[1] = (int) s_min[0];
    n_leaf_boxes[2] = (int) s_max[0];
  }

  if (mem_used != NULL) {
    mem_used[0] = (int) s_mean[5];
    mem_used[1] = (int) s_min[5];
    mem_used[2] = (int) s_max[5];
  }

  if (mem_allocated != NULL) {
    mem_allocated[0] = (int) s_mean[6];
    mem_allocated[1] = (int) s_min[6];
    mem_allocated[2] = (int) s_max[6];
  }

  return dim;
}

/*----------------------------------------------------------------------------
 * Display local statistics about a PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_dump_statistics(const PDM_box_tree_t  *bt)
{
  int i, j;
  PDM_box_tree_stats_t s;
  unsigned g_max_level_reached;
  PDM_g_num_t n_g_leaves, n_g_boxes, n_g_linked_boxes, n_g_spill_leaves;
  int g_min_linked_boxes, g_max_linked_boxes;
  double mean_linked_boxes, box_ratio;

  PDM_g_num_t count[5];

  int step = 0, delta = 0;
  const int n_steps = 5;

  if (bt == NULL)
    return;

  s = bt->stats;

  g_max_level_reached = s.max_level_reached;
  n_g_leaves = s.n_leaves;
  n_g_boxes = s.n_boxes;
  n_g_linked_boxes = s.n_linked_boxes;
  n_g_spill_leaves = s.n_spill_leaves;
  g_min_linked_boxes = s.min_linked_boxes;
  g_max_linked_boxes = s.max_linked_boxes;

  if (bt->comm != PDM_MPI_COMM_NULL) {

    PDM_g_num_t l_min[1], g_min[1];
    PDM_g_num_t l_max[2], g_max[2];
    PDM_g_num_t l_sum[3], g_sum[3];

    l_sum[0] = n_g_leaves;
    l_sum[1] = n_g_spill_leaves;
    l_sum[2] = n_g_linked_boxes;

    l_min[0] = g_min_linked_boxes;
    l_max[0] = s.max_level_reached;
    l_max[1] = g_max_linked_boxes;

    PDM_MPI_Allreduce(l_sum, g_sum, 3, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, bt->comm);
    PDM_MPI_Allreduce(l_min, g_min, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MIN, bt->comm);
    PDM_MPI_Allreduce(l_max, g_max, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, bt->comm);

    n_g_leaves = l_sum[0];
    n_g_spill_leaves = l_sum[1];
    n_g_linked_boxes = l_sum[2];

    g_min_linked_boxes = (int) g_min[0];
    g_max_level_reached = (int) g_max[0];
    g_max_linked_boxes = (int) g_max[1];
  }

  /* Redefine final statistics */

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double _n_g_linked_boxes = (double) n_g_linked_boxes;
  double _n_g_leaves = (double) n_g_leaves;
  double _n_g_boxes = (double) n_g_boxes;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

  mean_linked_boxes = _n_g_linked_boxes / _n_g_leaves;
  box_ratio = _n_g_linked_boxes / _n_g_boxes;

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  delta = g_max_linked_boxes - g_min_linked_boxes;

  if (delta > 0) {

    step = delta/n_steps;

    _build_histogram(bt,
                     0, /* start from root */
                     n_steps,
                     step,
                     g_min_linked_boxes,
                     count);

  } /* max - min > 0 */

  /* Print statistics and bounding boxes histogram */

  PDM_printf("\n"
             "Box tree statistics:\n\n");
  PDM_printf("  Number of children per leaf:              %d\n"
             "  Max number of bounding boxes for a leaf:  %d\n"
             "  Max value for box ratio (final/init):     %f\n"
             "  Max level allowed:                        %d\n\n",
             bt->n_children, (int)(bt->threshold),
             (double)(bt->max_box_ratio), (int)(bt->max_level));

  PDM_printf("  Max level reached:                  %5u\n"
             "  Number of leaves:                   %10llu\n"
             "  Leaves with n_boxes > max_n_boxes:  %10llu\n"
             "  Initial number of boxes:            %10llu\n"
             "  Number of linked boxes:             %10llu\n"
             "  Mean number of leaves per box:      %10.4g\n\n",
             g_max_level_reached, (unsigned long long)n_g_leaves,
             (unsigned long long)n_g_spill_leaves, (unsigned long long)n_g_boxes,
             (unsigned long long)n_g_linked_boxes,  box_ratio);

  PDM_printf("Number of linked boxes per box tree leaf:\n"
             "  Mean value:         %10.4g\n"
             "  min. value:         %10llu\n"
             "  max. value:         %10llu\n\n",
             mean_linked_boxes,
             (unsigned long long)(s.min_linked_boxes),
             (unsigned long long)(s.max_linked_boxes));

  if (delta > 0) { /* Number of elements in each subdivision */

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      PDM_printf("    %3d : [ %10llu; %10llu [ = %10llu\n",
                 i+1,
                 (unsigned long long)(g_min_linked_boxes + i*step),
                 (unsigned long long)(g_min_linked_boxes + j*step),
                 (unsigned long long)(count[i]));

    PDM_printf("    %3d : [ %10llu; %10llu ] = %10llu\n",
               n_steps,
               (unsigned long long)(g_min_linked_boxes + (n_steps - 1)*step),
               (unsigned long long)(g_max_linked_boxes),
               (unsigned long long)(count[n_steps - 1]));

  }
}

/*----------------------------------------------------------------------------
 * Dump an PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_dump(PDM_box_tree_t  *bt)
{
  PDM_box_tree_stats_t s;

  if (bt == NULL) {
    PDM_printf("\nBox tree: nil\n");
    return;
  }

  PDM_printf("\nBox tree: %p\n\n", (void *)bt);

  PDM_printf("  n_max_nodes:  %d\n\n"
             "  n_nodes:      %d\n",
             (int)(bt->n_max_nodes), (int)(bt->n_nodes));

  s = bt->stats;

  /* Print statistics and bounding boxes histogram */

  PDM_printf("  Number of children per leaf:              %d\n"
             "  Max number of bounding boxes for a leaf:  %d\n"
             "  Max value for box ratio (linked/init):    %f\n"
             "  Max level allowed:                        %d\n\n",
             bt->n_children, (int)(bt->threshold),
             (double)(bt->max_box_ratio), (int)(bt->max_level));

  PDM_printf("  Max level reached:                  %5u\n"
             "  Number of leaves:                   %10llu\n"
             "  Leaves with n_boxes > max_n_boxes:  %10llu\n"
             "  Initial number of boxes:            %10llu\n"
             "  Number of linked boxes:             %10llu\n",
             s.max_level_reached,
             (unsigned long long)(s.n_leaves),
             (unsigned long long)(s.n_spill_leaves),
             (unsigned long long)(s.n_boxes),
             (unsigned long long)(s.n_linked_boxes));

  PDM_printf("Bounding boxes related to each leaf of the box tree.\n"
             "  min. value:         %10llu\n"
             "  max. value:         %10llu\n\n",
             (unsigned long long)(s.min_linked_boxes),
             (unsigned long long)(s.max_linked_boxes));

  _dump_node(bt, 0);
}

/*----------------------------------------------------------------------------
 * Get minimum of maximum distance of boxes
 *
 * parameters:
 *   bt           <-- pointer to box tree structure
 *   n_boxes      --> Number of boxes in the closest leaf
 *   box_g_num[]  --> Global number of boxes in the closest leaf
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_min_dist_max_box
(
PDM_box_tree_t  *bt,
const int        n_pts,
double          *pts,
int             *box_id,
double          *box_max_dist
)
{

  int s_pt_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);
  int pos_stack = 0;

  int dim = bt->boxes->dim;

  int normalized = bt->boxes->normalized;
  const double *d = bt->boxes->d;

  double *_pts = pts;
  if (normalized) {
    _pts = malloc (sizeof(double) * 3 * n_pts);
    for (int i = 0; i < n_pts; i++) {
      const double *_pt_origin =  pts + 3 * i;
      double *_pt        = _pts + 3 * i;
      PDM_box_set_normalize ((PDM_box_set_t *) bt->boxes, _pt_origin, _pt);
    }
  }

  double extents2[2*dim];

  for (int i = 0; i < n_pts; i++) {

    const double *_pt = _pts + 3 * i;

    /* Init stack */

    box_id[i] = -1;
    box_max_dist[i] = HUGE_VAL;

    pos_stack = 0;
    stack[pos_stack] = 0; /* push root in th stack */

    _extents (dim, bt->nodes[0].morton_code, extents2);

    inbox_stack[pos_stack] = _box_dist2_min (dim,
                                             normalized,
                                             d,
                                             extents2,
                                             _pt,
                                             min_dist2_stack);

    pos_stack++;

    while (pos_stack > 0) {

      int id_curr_node = stack[--pos_stack];

      _node_t *curr_node = &(bt->nodes[id_curr_node]);

      if (curr_node->n_boxes == 0)
        continue;

      _extents (dim, curr_node->morton_code, extents2);

      double max_dist2;

      _box_dist2_max (dim,
                      normalized,
                      d,
                      extents2,
                      _pt,
                      &max_dist2);

      if (max_dist2 <= box_max_dist[i]) {

        if (!curr_node->is_leaf) {

          _push_child_in_stack_v0 (bt,
                                dim,
                                normalized,
                                d,
                                id_curr_node,
                                box_max_dist[i],
                                _pt,
                                &pos_stack,
                                stack,
                                inbox_stack,
                                min_dist2_stack,
                                0,
                                1);

        }

        else {

          for (int j = 0; j < curr_node->n_boxes; j++) {

            double box_max_dist2;

            int   _box_id = bt->box_ids[curr_node->start_id + j];
            const double *_box_extents =  bt->boxes->extents + _box_id*dim*2;

            _box_dist2_max (dim,
                            normalized,
                            d,
                            _box_extents,
                            _pt,
                            &box_max_dist2);

            if (box_max_dist2 < box_max_dist[i]) {
              box_id[i] = _box_id;
              box_max_dist[i] = box_max_dist2;
            }
          }

        }
      }
    }
  }

  free (stack);
  free (inbox_stack);
  free (min_dist2_stack);

  if (_pts != pts) {
    free (_pts);
  }

}


/*----------------------------------------------------------------------------
 * Get minimum of maximum distance of boxes
 *
 * parameters:
 *   bt                <-- pointer to box tree structure
 *   n_pts             <-- Number of points
 *   pts               <-- Point coordinates (size = 3 * n_pts)
 *   upper_bound_dist2 <-- Upper bound of the square of the distance (size = n_pts)
 *   i_boxes           --> Index of boxes (size = n_pts + 1)
 *   boxes             --> Boxes (size = i_boxes[n_pts])
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_closest_upper_bound_dist_boxes_get
(
PDM_box_tree_t  *bt,
const int        n_pts,
double           pts[],
double           upper_bound_dist2[],
int             *i_boxes[],
int             *boxes[]
)
{
  //printf ("PDM_box_tree_closest_upper_bound_dist_boxes_get\n");

  int normalized = bt->boxes->normalized;
  const double *d = bt->boxes->d;

  double *_pts = pts;
  if (normalized) {
    _pts = malloc (sizeof(double) * 3 * n_pts);
    for (int i = 0; i < n_pts; i++) {
      const double *_pt_origin = pts + 3 * i;
      double *_pt        = _pts + 3 * i;
      PDM_box_set_normalize ((PDM_box_set_t *)bt->boxes, _pt_origin, _pt);
    }
  }

  int s_pt_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);

  *i_boxes = malloc (sizeof(int) * (n_pts + 1));
  int *_i_boxes = *i_boxes;

  for (int i = 0; i < n_pts + 1; i++) {
    _i_boxes[i] = 0;
  }

  int tmp_s_boxes = 4 * n_pts;
  *boxes = malloc (sizeof(int) * tmp_s_boxes);
  int *_boxes = *boxes;

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);

  int pos_stack = 0;

  int dim = bt->boxes->dim;

  int idx_box = 0;

  int n_boxes = bt->boxes->n_boxes;

  int *tag = malloc(sizeof(int) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    tag[i] = 0;
  }
  int n_visited_boxes = 0;

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  size_t n_node = 0;
  size_t n_node_vid = 0;

  double extents2[2*dim];

  for (int i = 0; i < n_pts; i++) {

    const double *_pt = _pts + 3 * i;
    int flag = 0;

    /* Init stack :  push root */

    pos_stack = 0;
    stack[pos_stack] = 0; /* push root in th stack */
    _extents (dim, bt->nodes[0].morton_code, extents2);
    inbox_stack[pos_stack] = _box_dist2_min (dim,
                                             normalized,
                                             d,
                                             extents2,
                                             _pt,
                                             min_dist2_stack);


    pos_stack++;
    n_visited_boxes = 0;

    while (pos_stack > 0) {

      int id_curr_node = stack[--pos_stack];

      _node_t *curr_node = &(bt->nodes[id_curr_node]);

      if (curr_node->n_boxes == 0)
        continue;

      n_node++;
      if (curr_node->n_boxes == 0)
        n_node_vid++;

      double min_dist2 = min_dist2_stack[pos_stack];

      int inbox = inbox_stack[pos_stack];

      if ((min_dist2 <= upper_bound_dist2[i]) || (inbox == 1)) {
        if (!curr_node->is_leaf) {

          _push_child_in_stack_v0 (bt,
                                dim,
                                normalized,
                                d,
                                id_curr_node,
                                upper_bound_dist2[i],
                                _pt,
                                &pos_stack,
                                stack,
                                inbox_stack,
                                min_dist2_stack,
                                flag,
                                0);

        }

        else {

          for (int j = 0; j < curr_node->n_boxes; j++) {

            double box_min_dist2;

            int   _box_id = bt->box_ids[curr_node->start_id + j];

            if (tag[_box_id] == 0) {

              const double *_box_extents =  bt->boxes->extents + _box_id*dim*2;

              inbox = _box_dist2_min (dim,
                                      normalized,
                                      d,
                                      _box_extents,
                                      _pt,
                                      &box_min_dist2);

              if ((box_min_dist2 <= upper_bound_dist2[i]) || (inbox == 1)) {
                if (idx_box >= tmp_s_boxes) {
                  tmp_s_boxes *= 2;
                  *boxes = realloc (*boxes, sizeof(int) * tmp_s_boxes);
                  _boxes = *boxes;
                }
                _boxes[idx_box++] = _box_id;
                _i_boxes[i+1]++;
              }
              visited_boxes[n_visited_boxes++] = _box_id;
              tag[_box_id] = 1;
            }
          }
        }
      }
    }

    for (int j = 0; j < n_visited_boxes; j++) {
      tag[visited_boxes[j]] = 0;
    }

  }

  //printf ("[%d] Parours arbre : %ld \n", iappel, n_node);
  iappel+=1;

  for (int i = 0; i < n_pts; i++) {
    _i_boxes[i+1] += _i_boxes[i];
  }

  *boxes = realloc (*boxes, sizeof(int) * _i_boxes[n_pts]);

  if (pts != _pts) {
    free (_pts);
  }

  free (tag);
  free (stack);
  free (inbox_stack);
  free (min_dist2_stack);
  free (visited_boxes);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
