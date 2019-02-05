#ifndef __PDM_BOX_H__
#define __PDM_BOX_H__

/*============================================================================
 * Handle boxes aligned with Cartesian axes.
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#include <stdbool.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_morton.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/* Collection of boxes */

typedef struct _PDM_box_set_t PDM_box_set_t;

/* Distribution on octree or quadtree */

typedef struct _PDM_box_distrib_t PDM_box_distrib_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a set of boxes and initialize it.
 *
 * parameters:
 *   dim              <-- spatial dimension
 *   normalize        <-- 1 if boxes are to be normalized, 0 otherwize
 *   allow_projection <-- if 1, project to lower dimension if all boxes
 *                        are cut by the median plane of the set.
 *   n_boxes          <-- number of elements to create
 *   box_gnum         <-- global numbering of boxes
 *   extents          <-- coordinate extents (size: n_boxes*dim*2, as
 *                        xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 *   origin   <--  initial location (size: n_boxes*3, as
 *                        iproc, ipart, local num, ...)
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   a new allocated pointer to a PDM_box_set_t structure.
 *---------------------------------------------------------------------------*/

PDM_box_set_t *
PDM_box_set_create(int                dim,
                   int                normalize,
                   int                allow_projection,
                   int                n_boxes,
                   const PDM_g_num_t  *box_gnum,
                   const double      *box_extents,
                   const int          n_part_orig,
                   const int         *n_boxes_orig,
                   const int         *origin,
                   PDM_MPI_Comm           comm);

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_set_t structure.
 *
 * parameters:
 *   boxes <-> pointer to pointer to the PDM_box_set_t structure to delete
 *----------------------------------------------------------------------------*/

void
PDM_box_set_destroy(PDM_box_set_t  **boxes);

/*----------------------------------------------------------------------------
 * Return the dimension associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   associated spatial dimension
 *---------------------------------------------------------------------------*/

int
PDM_box_set_get_dim(const PDM_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return the local number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   local number of boxes
 *---------------------------------------------------------------------------*/

int
PDM_box_set_get_size(const PDM_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return the global number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   local number of boxes
 *---------------------------------------------------------------------------*/

PDM_g_num_t
PDM_box_set_get_global_size(const PDM_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return extents associated with a set of boxes.
 *
 * The extents array is organized in the following fashion:
 * {x_min_0, y_min_0, ..., x_max_0, y_max_0, ...
 *  x_min_n, y_min_n, ..., x_max_n, y_max_n, ...}
 *
 * Its size is thus: n_boxes * dim * 2.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to extents array
 *---------------------------------------------------------------------------*/

const double *
PDM_box_set_get_extents(PDM_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return global numbers associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to global box numbers array
 *---------------------------------------------------------------------------*/

const PDM_g_num_t *
PDM_box_set_get_g_num(PDM_box_set_t  *boxes);


/*----------------------------------------------------------------------------
 * Return initial location associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to initial location array
 *---------------------------------------------------------------------------*/

const int *
PDM_box_set_origin_get(PDM_box_set_t  *boxes);


/*----------------------------------------------------------------------------
 * Build a Morton_index to get a well-balanced distribution of the boxes.
 *
 * parameters:
 *  boxes      <-- pointer to associated PDM_box_set_t structure
 *  distrib    <-> pointer to a PDM_box_distrib_t structure
 *  n_leaves   <-- number of leaves with weight > 0
 *  leaf_codes <-- Morton code for each leaf
 *  weight     <-- number of boxes related to each leaf
 *---------------------------------------------------------------------------*/

void
PDM_box_set_build_morton_index(const PDM_box_set_t  *boxes,
                               PDM_box_distrib_t    *distrib,
                               int             n_leaves,
                               PDM_morton_code_t    *leaf_codes,
                               int            *weight);

/*----------------------------------------------------------------------------
 * Redistribute boxes over the ranks according to the Morton index to
 * assume a better balanced distribution of the boxes.
 *
 * parameters:
 *  box_distrib <--  data structure on box distribution
 *  box_set     <->  pointer to the structure to redistribute
 *---------------------------------------------------------------------------*/

void
PDM_box_set_redistribute(const PDM_box_distrib_t  *box_distrib,
                         PDM_box_set_t            *boxes);

/*----------------------------------------------------------------------------
 * Dump a PDM_box_set_t structure.
 *
 * parameters:
 *   box_set   <-- pointer to the PDM_box_t structure
 *   verbosity <-- verbosity level (0 or 1)
 *----------------------------------------------------------------------------*/

void
PDM_box_set_dump(const PDM_box_set_t  *boxes,
                 int                   verbosity);


/*----------------------------------------------------------------------------
 * Receive data from origin for any box
 *
 * parameters:
 *   box_set                <-- pointer to the PDM_box_t structure
 *   t_stride               <-- Type of stride
 *   stride_cst             <-- Constant stride
 *   data_size              <-- Size of data
 *   origin_distrib_stride  <-- Origin stride distribution
 *   origin_distrib_data    <-- Origin data distribution
 *   current_distrib_stride <-> Current stride distribution (Allocate and compute if input is NULL, 
 *                              otherwise nothing)
 *   current_distrib_data   --> Current data distribution
 *
 * return : size of current_distrib_data
 *----------------------------------------------------------------------------*/

void
PDM_box_set_recv_data_from_origin_distrib 
(
 PDM_box_set_t  *boxes,
 PDM_stride_t   t_stride,
 int            stride_cst,
 size_t         data_size,
 int          **origin_distrib_stride,
 void         **origin_distrib_data, 
 int          **current_distrib_stride,
 void         **current_distrib_data 
 );


/*----------------------------------------------------------------------------
 * Send data to origin for any box
 *
 * parameters:
 *   box_set                <-- pointer to the PDM_box_t structure
 *   t_stride               <-- Type of stride
 *   stride_cst             <-- Constant stride
 *   data_size              <-- Size of data
 *   current_distrib_stride <-- Current stride distribution
 *   current_distrib_data   <-- Current data distribution
 *   origin_distrib_stride  --> Origin stride distribution
 *   origin_distrib_data    --> Origin data distribution
 *
 * return : size of origin_distrib_data
 *----------------------------------------------------------------------------*/

void
PDM_box_set_send_data_to_origin_distrib 
(
 PDM_box_set_t  *boxes,
 PDM_stride_t   t_stride,
 int            stride_cst,
 size_t         data_size,
 int           *current_distrib_stride,
 void          *current_distrib_data, 
 int          **origin_distrib_stride,
 void         **origin_distrib_data 
);

/*----------------------------------------------------------------------------
 * Create a PDM_box_distrib_t structure.
 *
 * parameters:
 *   n_boxes   <-- number of boxes
 *   n_g_boxes <-- global number of boxes
 *   max_level <-- max level reached locally in the related tree
 *   comm      <-- MPI comm. on which distribution takes place
 *
 * returns:
 *   a pointer to a new allocated PDM_box_distrib_t structure.
 *---------------------------------------------------------------------------*/

PDM_box_distrib_t *
PDM_box_distrib_create(int  n_boxes,
                       PDM_g_num_t  n_g_boxes,
                       int        max_level,
                       PDM_MPI_Comm   comm);

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_distrib_t structure.
 *
 * parameters:
 *   distrib <-> pointer to pointer to the structure to destroy
 *---------------------------------------------------------------------------*/

void
PDM_box_distrib_destroy(PDM_box_distrib_t  **distrib);

/*----------------------------------------------------------------------------
 * Delete redundancies in box distribution
 *
 * parameters:
 *   distrib <-> pointer to the PDM_box_distrib_t structure
 *---------------------------------------------------------------------------*/

void
PDM_box_distrib_clean(PDM_box_distrib_t  *distrib);

/*----------------------------------------------------------------------------
 * Display a histogramm on leaves associated to the boxes and
 * several other pieces of information (min, max, ...)
 *
 * parameters:
 *   distrib <-- pointer to the PDM_box_distrib_t structure
 *   comm    <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
PDM_box_distrib_dump_statistics(const PDM_box_distrib_t  *distrib,
                                PDM_MPI_Comm                  comm);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_BOX_H__ */
