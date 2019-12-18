/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dist_cloud_surf.h"
//#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
//#include "pdm_octree.h"
//#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"
//#include "pdm_hash_tab.h"
#include "pdm_gnum.h"

#include "pdm_dist_cellcenter_surf.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 4

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_Queue_t
 * \brief  Define a queue
 *
 */

typedef struct {

  int beg;       /*!< Beginning of the queue */
  int end;       /*!< End of the queue */

  int size;      /*!< Size of the queue */
  int n_elt;    /*!< Number of elements */

  int *queue;    /*!< queue array */

} _PDM_queue_t;


/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  FIRST_THICKNESS               = 1,
  PROPAGATION                   = 2,
  END                           = 3,

} _ol_timer_step_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {
  PDM_g_num_t          g_n_cell;      /*!< Global number of cells */
  PDM_g_num_t          g_n_face;      /*!< Global number of faces */
  PDM_g_num_t          g_n_vtx;       /*!< Global number of vertices */
  int                  n_part;        /*!< Number of partitions */
  int                 *n_cell;        /*!< Number of cells (size = \n_part) */
  int                 *n_face;        /*!< Number of faces (size = \n_part) */
  int                 *n_vtx;         /*!< Number of vertices (size = \n_part) */
  const int          **cell_face_idx; /*!< Cell -> face connectivity index */
        int          **cell_cell;     /*!< Cell -> cell connectivity */
  const int          **cell_face;     /*!< Cell -> face connectivity */
  const double       **cell_center;   /*!< Cell center or NULL */
  const PDM_g_num_t  **cell_ln_to_gn; /*!< Cell local numbering to global numbering */
  const int          **face_vtx_idx;  /*!< Face -> vtx connectivity index */
  const int          **face_vtx;      /*!< Face -> vtx connectivity */
        int          **face_cell;     /*!< Face -> cell connectivity */
  const PDM_g_num_t  **face_ln_to_gn; /*!< Face local numbering to global numbering */
  const double       **coords;        /*!< Vertex coordinates */
  const PDM_g_num_t  **vtx_ln_to_gn;  /*!< Vertex local numbering to global numbering */

  double      **closest_elt_distance;  /*!< Distance to the closest element */
  double      **closest_elt_projected; /*!< Projected point on the closest element */
  PDM_g_num_t **closest_elt_gnum;      /*!< Global numbering of the closest element */

} _PDM_vol_mesh_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  PDM_MPI_Comm comm;  /*!< MPI communicator */

  PDM_surf_mesh_t *surf_mesh;  /*!< Surface mesh pointer */

  _PDM_vol_mesh_t *vol_mesh; /*!< Volume mesh pointer */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */


} _PDM_dist_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dists   = NULL;

static int idebug = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return a new queue
 *
 * \param [in]   ini_size      Initial size of the queue
 *
 */

static _PDM_queue_t *
_PDM_queue_init
(
 const int ini_size
)
{
  _PDM_queue_t *q = malloc (sizeof(_PDM_queue_t));

  q->size = ini_size;

  q->queue = malloc(sizeof(int) * ini_size);
  q->beg = 0;
  q->end = -1;
  q->n_elt = 0;

  return q;

}


/**
 *
 * \brief Free a new queue
 *
 * \param [in]   q      Queue
 *
 */

static _PDM_queue_t *
_PDM_queue_free
(
 _PDM_queue_t *q
)
{
  free (q->queue);
  free (q);
  return NULL;
}


/**
 *
 * \brief Update size of the queue
 *
 * \param [in]   q      Queue
 *
 */

static void
_PDM_queue_update_size
(
 _PDM_queue_t *q,
 const int new_size
)
{
  assert (new_size > q->size);

  /* printf ("_PDM_queue_update_size 1 : %d\n",  q->size); */
  /* fflush(stdout); */

  q->queue = realloc (q->queue, sizeof(int) * new_size);

  const int _beg = q->beg;
  const int _end = q->end;

  assert (_beg != _end);

  const int step = new_size - q->size;

  if (_beg > _end) {

    for (int i = q->size - 1; i >= _beg; i--) {
      q->queue[i+step] = q->queue[i];
    }
    q->beg += step;

  }
  q->size = new_size;
  /* printf ("_PDM_queue_update_size 2 : %d\n",  q->size); */
  /* fflush(stdout); */
}


/**
 *
 * \brief Push value into the queue
 *
 * \param [in]   q      Queue
 * \param [in]   val    Value
 *
 */

static void
_PDM_queue_push
(
 _PDM_queue_t *q,
 const int val
)
{

  /* printf ("_PDM_queue_push val n_elt: %d %d\n", val, q->n_elt); */
  /* fflush(stdout); */

  if (q->n_elt >= q->size) {
    _PDM_queue_update_size (q, 2*q->size);
  }

  q->end = (q->end + 1) % q->size;

  q->queue[q->end] = val;

  q->n_elt += 1;

  /* printf ("_PDM_queue_push n_elt : %d\n", q->n_elt); */
  /* fflush(stdout); */
}


/**
 *
 * \brief Pull the firts value in the queue
 *
 * \param [in]   q      Queue
 * \param [out]  val    Value
 *
 * \return 1 if the queue is empty, 0 otherwise
 */

static int
_PDM_queue_pull
(
 _PDM_queue_t *q,
 int  *val
)
{

  int is_empty = 0;

  /* printf ("_PDM_queue_pull n_elt: %d\n", q->n_elt); */
  /* fflush(stdout); */

  if (q->n_elt <= 0) {
    is_empty = 1;
  /* printf ("_PDM_queue_pull val n_elt: empty \n"); */
  /* fflush(stdout); */
  }

  else {
    *val = q->queue[q->beg];
    q->beg = (q->beg + 1) % q->size;
    q->n_elt += -1;
  /* printf ("_PDM_queue_pull val n_elt: %d %d\n", *val, q->n_elt); */
  /* fflush(stdout); */
  }

  return is_empty;
}


/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _PDM_dist_t *
_get_from_id
(
 int  id
)
{
  _PDM_dist_t *dist = (_PDM_dist_t *) PDM_Handles_get (_dists, id);

  if (dist == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_dist error : Bad identifier\n");
  }

  return dist;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute distance of cell centers of a volume mesh
 * to a surface mesh
 *
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 *
 */

int
PDM_dist_cellcenter_surf_create
(
 const PDM_MPI_Comm comm
)
{
  if (_dists == NULL) {
    _dists = PDM_Handles_create (4);
  }

  _PDM_dist_t *dist = (_PDM_dist_t *) malloc(sizeof(_PDM_dist_t));

  int id = PDM_Handles_store (_dists, dist);

  dist->surf_mesh = NULL;
  dist->vol_mesh = NULL;
  dist->comm = comm;

  dist->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    dist->times_elapsed[i] = 0.;
    dist->times_cpu[i] = 0.;
    dist->times_cpu_u[i] = 0.;
    dist->times_cpu_s[i] = 0.;
  }

  return id;
}

void
PDM_dist_cellcenter_surf_create_cf
(
 const PDM_MPI_Fint comm,
 int *id
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_dist_cellcenter_surf_create (_comm);
}


/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_dist_cellcenter_surf_surf_mesh_global_data_set
(
 const int         id,
 const PDM_g_num_t n_g_face,
 const PDM_g_num_t n_g_vtx,
 const int         n_part
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  if (dist->surf_mesh != NULL) {
    dist->surf_mesh = PDM_surf_mesh_free (dist->surf_mesh);
  }

  dist->surf_mesh =
    PDM_surf_mesh_create (n_g_face, n_g_vtx, n_part, dist->comm);
}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_dist_cellcenter_surf_surf_mesh_part_set
(
 const int          id,
 const int          i_part,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
)
{

  _PDM_dist_t *dist = _get_from_id (id);

  PDM_surf_mesh_part_input (dist->surf_mesh,
                            i_part,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);
}

/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_g_face       Global number of cells
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_dist_cellcenter_surf_vol_mesh_global_data_set
(
 const int         id,
 const PDM_g_num_t n_g_cell,
 const PDM_g_num_t n_g_face,
 const PDM_g_num_t n_g_vtx,
 const int         n_part
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  assert (dist->vol_mesh == NULL);

  dist->vol_mesh = malloc (sizeof(_PDM_vol_mesh_t));

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  _vol_mesh->g_n_cell = n_g_cell;
  _vol_mesh->g_n_face = n_g_face;
  _vol_mesh->g_n_vtx = n_g_vtx;
  _vol_mesh->n_part = n_part;

  _vol_mesh->n_cell = malloc (sizeof(int) * n_part);
  _vol_mesh->n_face = malloc (sizeof(int) * n_part);
  _vol_mesh->n_vtx = malloc (sizeof(int) * n_part);

  _vol_mesh->cell_face_idx = malloc (sizeof(int *) * n_part);
  _vol_mesh->cell_face = malloc (sizeof(int *) * n_part);
  _vol_mesh->cell_cell = malloc (sizeof(int *) * n_part);
  _vol_mesh->cell_center = malloc (sizeof(double *) * n_part);
  _vol_mesh->cell_ln_to_gn = malloc (sizeof(PDM_g_num_t *) * n_part);

  _vol_mesh->face_vtx_idx = malloc (sizeof(int *) * n_part);
  _vol_mesh->face_vtx = malloc (sizeof(int *) * n_part);
  _vol_mesh->face_cell = malloc (sizeof(int *) * n_part);
  _vol_mesh->face_ln_to_gn = malloc (sizeof (PDM_g_num_t *) * n_part);

  _vol_mesh->coords = malloc (sizeof(double *) * n_part);
  _vol_mesh->vtx_ln_to_gn = malloc (sizeof(PDM_g_num_t *) * n_part);

  _vol_mesh->closest_elt_distance = malloc (sizeof(double *) * n_part);
  _vol_mesh->closest_elt_projected = malloc (sizeof(double *) * n_part);
  _vol_mesh->closest_elt_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of faces
 * \param [in]   cell_face_idx Index in the face -> vertex connectivity
 * \param [in]   cell_face     face -> vertex connectivity
 * \param [in]   cell_center   Cell center or NULL
 * \param [in]   cell_ln_to_gn Local cell numbering to global face numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_dist_cellcenter_surf_vol_mesh_part_set
(
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_face_idx,
 const int         *cell_face,
 const double      *cell_center,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  assert (dist->vol_mesh != NULL);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  _vol_mesh->n_cell[i_part] = n_cell;

  _vol_mesh->cell_face_idx[i_part] = cell_face_idx;
  _vol_mesh->cell_face[i_part] = cell_face;
  _vol_mesh->cell_center[i_part] = cell_center;
  _vol_mesh->cell_ln_to_gn[i_part] = cell_ln_to_gn;

  _vol_mesh->n_face[i_part] = n_face;
  _vol_mesh->face_vtx_idx[i_part] = face_vtx_idx;
  _vol_mesh->face_vtx[i_part] = face_vtx;
  _vol_mesh->face_cell[i_part] = malloc (sizeof(int) * 2 * n_face);

  for (int i = 0; i < 2*n_face; i++) {
    _vol_mesh->face_cell[i_part][i] = 0;
  }

  /* printf ("cell_face : \n"); */
  for (int i = 0; i < n_cell; i++) {
    for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
      int ifac = PDM_ABS(cell_face[j])-1;
      /* printf (" %d", ifac+1); */
      if (_vol_mesh->face_cell[i_part][2*ifac] == 0) {
        _vol_mesh->face_cell[i_part][2*ifac] = i + 1;
      }
      else{
        _vol_mesh->face_cell[i_part][2*ifac+1] = i + 1;
      }
    }
    /* printf ("\n"); */
  }

  _vol_mesh->cell_cell[i_part] = malloc(sizeof(int) * cell_face_idx[n_cell]);
  for (int i = 0; i < n_cell; i++) {
    for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
      int ifac = PDM_ABS (cell_face[j])-1;
      if (_vol_mesh->face_cell[i_part][2*ifac] == i+1) {
        int icell1 = _vol_mesh->face_cell[i_part][2*ifac+1];
        _vol_mesh->cell_cell[i_part][j] = icell1;
      }
      else if (_vol_mesh->face_cell[i_part][2*ifac+1] == i+1) {
        int icell1 = _vol_mesh->face_cell[i_part][2*ifac];
        _vol_mesh->cell_cell[i_part][j] = icell1;
      }
    }
  }

  /* printf ("cell_cell : \n"); */
  /* for (int i = 0; i < n_cell; i++) { */
  /*   if (i == 1069) { */
  /*   printf ("%d :", i+1); */
  /*   for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) { */
  /*     printf (" %d", _vol_mesh->cell_cell[i_part][j]); */
  /*   } */
  /*   printf("\n"); */
  /*   } */
  /* } */



  _vol_mesh->face_ln_to_gn[i_part] = face_ln_to_gn;

  _vol_mesh->n_vtx[i_part] = n_vtx;
  _vol_mesh->coords[i_part] = coords;
  _vol_mesh->vtx_ln_to_gn[i_part] = vtx_ln_to_gn;

}


/**
 *
 * \brief Compute distance
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_dist_cellcenter_surf_compute
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;
  PDM_surf_mesh_t *_surf_mesh = dist->surf_mesh;

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  dist->times_elapsed[BEGIN] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[BEGIN]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  /* First step : Look for boundary cells in the volume mesh */

  /* fflush(stdout); */
  /* printf(" --- First step\n"); */

  PDM_timer_hang_on(dist->timer);
  b_t_elapsed = PDM_timer_elapsed(dist->timer);
  b_t_cpu     = PDM_timer_cpu(dist->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  int _t_n_face = 0;
  for (int i = 0; i < _vol_mesh->n_part; i++) {
    _t_n_face += _vol_mesh->n_face[i];
  }

  int *bound_cell_idx = malloc (sizeof(int) * (_vol_mesh->n_part + 1));
  bound_cell_idx[0] = 0;

  int n_bound_cell = 0;
  int *bound_cell = malloc (sizeof(int) * _t_n_face);
  PDM_g_num_t *bound_parent_gnum = malloc (sizeof(PDM_g_num_t) * _t_n_face);
  double * bound_cell_center = malloc (sizeof(double) * 3 * _t_n_face);

  for (int i = 0; i < _vol_mesh->n_part; i++) {

    bound_cell_idx[i+1] = bound_cell_idx[i];

    const int *_face_cell = _vol_mesh->face_cell[i];
    const int _n_face = _vol_mesh->n_face[i];
    const int _n_cell = _vol_mesh->n_cell[i];
    const double *_cell_center = _vol_mesh->cell_center[i];
    const PDM_g_num_t *_gnum = _vol_mesh->cell_ln_to_gn[i];

    int *tag_cell = malloc(sizeof(int) * _n_cell);

    for (int j = 0; j < _n_cell; j++) {
      tag_cell[j] = 0;
    }

    /* printf ("face_cell :\n"); */
    for (int j = 0; j < _n_face; j++) {
      /* printf ("%d %d\n", _face_cell[2*j], _face_cell[2*j+1]); */
      if (_face_cell[2*j+1] == 0) {
        int icell = _face_cell[2*j] - 1;
        if (tag_cell[icell] == 0) {
          bound_cell[n_bound_cell] = icell;
          bound_parent_gnum[n_bound_cell] = _gnum[icell];
          bound_cell_idx[i+1] += 1;
          for (int k = 0; k < 3; k++) {
            bound_cell_center[3*n_bound_cell+k] = _cell_center[3*icell+k];
          }
          n_bound_cell += 1;
          tag_cell[icell] = 1;
        }
      }
    }

    free (tag_cell);

  }

  bound_cell = realloc (bound_cell, sizeof(int) * n_bound_cell);
  bound_cell_center = realloc (bound_cell_center, sizeof(double) * 3 * n_bound_cell);
  bound_parent_gnum = realloc (bound_parent_gnum, sizeof(PDM_g_num_t) * n_bound_cell);

  /* printf ("n_bound_cell : %d\n", n_bound_cell); */

  /* Second step : Compute distance to the surface mesh for the cell centers of boundary cells
                   Call PDM_mesh_dist */

  /* fflush(stdout); */
  /* printf(" --- Second step\n"); */

  int id_bound_dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_SURFACE_MESH,
                                            1,
                                            dist->comm);

  PDM_dist_cloud_surf_n_part_cloud_set (id_bound_dist, 0, 1);

  int id_gnum = PDM_gnum_create (3, 1, 0, 0.001, dist->comm);

  PDM_gnum_set_from_parents (id_gnum, 0,  n_bound_cell, bound_parent_gnum);

  PDM_gnum_compute (id_gnum);

  PDM_g_num_t *bound_gnum = PDM_gnum_get (id_gnum, 0);

  PDM_gnum_free (id_gnum, 1);

  PDM_dist_cloud_surf_cloud_set (id_bound_dist,
                           0,
                           0,
                           n_bound_cell,
                           bound_cell_center,
                           bound_gnum);

  PDM_dist_cloud_surf_surf_mesh_map (id_bound_dist, dist->surf_mesh);

  PDM_dist_cloud_surf_compute (id_bound_dist);

  double      *closest_elt_distance_bound = NULL;
  double      *closest_elt_projected_bound = NULL;
  PDM_g_num_t *closest_elt_gnum_bound = NULL;

  PDM_dist_cloud_surf_get (id_bound_dist,
                     0,
                     0,
                     &closest_elt_distance_bound,
                     &closest_elt_projected_bound,
                     &closest_elt_gnum_bound);

  PDM_dist_cloud_surf_free (id_bound_dist, 1);

  free (bound_gnum);
  free (bound_cell_center);
  free (bound_parent_gnum);

  PDM_timer_hang_on(dist->timer);
  e_t_elapsed = PDM_timer_elapsed(dist->timer);
  e_t_cpu     = PDM_timer_cpu(dist->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

  dist->times_elapsed[FIRST_THICKNESS] += e_t_elapsed - b_t_elapsed;
  dist->times_cpu[FIRST_THICKNESS]     += e_t_cpu - b_t_cpu;
  dist->times_cpu_u[FIRST_THICKNESS]   += e_t_cpu_u - b_t_cpu_u;
  dist->times_cpu_s[FIRST_THICKNESS]   += e_t_cpu_s - b_t_cpu_s;

  /* Third step : Get vertices about closest faces :
     - Build connectivity with coordinates
     - Part to block faces */

  /* fflush(stdout); */
  /* printf(" --- Third step\n"); */

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(dist->timer);

  int n_part_sm = PDM_surf_mesh_n_part_get (_surf_mesh);

  int *n_face_sm = malloc (sizeof(int) * n_part_sm);

  const PDM_g_num_t **gnum_sm = malloc (sizeof(PDM_g_num_t *) * n_part_sm);

  for (int i = 0; i < n_part_sm; i++) {
    n_face_sm[i] = PDM_surf_mesh_part_n_face_get (_surf_mesh, i);
    gnum_sm[i]   = PDM_surf_mesh_part_face_g_num_get (_surf_mesh, i);
  }

  PDM_part_to_block_t *ptb =
    PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                              1.,
                              (PDM_g_num_t **) gnum_sm,
                              NULL,
                              n_face_sm,
                              n_part_sm,
                              dist->comm);

  double **sm_face_coords = malloc (sizeof(double *) * n_part_sm);
  int **sm_face_coords_stride = malloc (sizeof(int *) * n_part_sm);

  for (int i = 0; i < n_part_sm; i++) {
    const int *_sm_face_vtx_idx =
      PDM_surf_mesh_part_face_vtx_idx_get (_surf_mesh, i);
    const int *_sm_face_vtx =
      PDM_surf_mesh_part_face_vtx_get (_surf_mesh, i);
    const double *_sm_vtx = PDM_surf_mesh_part_vtx_get (_surf_mesh, i);

    sm_face_coords[i] = malloc (sizeof(double) * 3 * _sm_face_vtx_idx[n_face_sm[i]]);
    sm_face_coords_stride[i] = malloc (sizeof(int) * n_face_sm[i]);

    for (int j = 0; j < n_face_sm[i]; j++) {
      sm_face_coords_stride[i][j] =
        3 * (_sm_face_vtx_idx[j+1] - _sm_face_vtx_idx[j]);
    }

    for (int j = 0; j < _sm_face_vtx_idx[n_face_sm[i]]; j++) {

      int iface = _sm_face_vtx[j]-1;
      const double *_coords = _sm_vtx + 3 * iface;

      for (int k = 0; k < 3; k++) {
        sm_face_coords[i][3*j+k] = _coords[k];
      }
    }

  }

  int     *sm_block_face_coords_stride;
  double  *sm_block_face_coords;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          0,
                          sm_face_coords_stride,
                          (void **) sm_face_coords,
                          &sm_block_face_coords_stride,
                          (void **) &sm_block_face_coords);

  const PDM_g_num_t  *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb);

  PDM_block_to_part_t *btp =
    PDM_block_to_part_create (blockDistribIdx,
                              (const PDM_g_num_t **) &closest_elt_gnum_bound,
                              &n_bound_cell,
                              1,
                              dist->comm);

  int  **sm_info_bound_cell_stride;
  double **sm_info_bound_cell;

  PDM_block_to_part_exch2 (btp,
                           sizeof(double),
                           PDM_STRIDE_VAR,
                           sm_block_face_coords_stride,
                           (void *) sm_block_face_coords,
                           &sm_info_bound_cell_stride,
                           (void ***) &sm_info_bound_cell);

  int *sm_info_bound_cell_idx = malloc(sizeof(int) * (n_bound_cell + 1));

  sm_info_bound_cell_idx[0] = 0;
  for (int i = 0; i < n_bound_cell; i++) {
    sm_info_bound_cell_idx[i+1] = sm_info_bound_cell_stride[0][i] + sm_info_bound_cell_idx[i];
  }

  PDM_part_to_block_free (ptb);
  PDM_block_to_part_free (btp);

  free (sm_block_face_coords_stride);
  free (sm_block_face_coords);

  /* Fourth step : Get connectivity with coordinates +
        Compute distance to the surface mesh for the other centers from the distance
        of the cell centers of boundary cells  */

  /* fflush(stdout); */
  /* printf(" --- Fourth step\n"); */

  for (int i = 0; i < _vol_mesh->n_part; i++) {

    // - Traitement element haut de pile :
    //       * Calcul distance à partir de toutes les cellules voisines ayant une distance
    //       * Si Changement de valeur on remet toutes les cellules voisines dans la pile sauf la source
    // - Arret losrque la pile est vide

    /* fflush(stdout); */
    /* printf(" --- Fourth step1\n"); */

    int          n_bound_cell_part = bound_cell_idx[i+1] - bound_cell_idx[i];

    const int    _n_cell = _vol_mesh->n_cell[i];

    double      *_closest_elt_distance = malloc (sizeof(double) * _n_cell);
    double      *_closest_elt_projected = malloc (sizeof(double) * 3 * _n_cell);
    PDM_g_num_t *_closest_elt_gnum = malloc (sizeof(PDM_g_num_t) * _n_cell);

    _vol_mesh->closest_elt_distance[i] = _closest_elt_distance;
    _vol_mesh->closest_elt_projected[i] = _closest_elt_projected;
    _vol_mesh->closest_elt_gnum[i] = _closest_elt_gnum;

    int         *_cell_cell = _vol_mesh->cell_cell[i];
    const int   *_cell_cell_idx = _vol_mesh->cell_face_idx[i];
    const double *_cell_center = _vol_mesh->cell_center[i];

    int *in_queue = malloc(sizeof(int) * _n_cell); /* -1 : boundary cell,
                                                       0 : outside the queue,
                                                       1 : inside the queue */


    int *cell_to_closest_bound_cell = malloc(sizeof(int) * _n_cell);
    for (int j = 0; j < _n_cell; j++) {
      cell_to_closest_bound_cell[j] = -1;
    }

    for (int j = bound_cell_idx[i]; j < bound_cell_idx[i+1]; j++) {
      cell_to_closest_bound_cell[bound_cell[j]] = j;
    }

    /* Initialize distances to a huge value */

    for (int j = 0; j < _n_cell; j++) {
      _closest_elt_distance[j] = HUGE_VAL;
      _closest_elt_gnum[j] = -1;
      in_queue[j] = 0;
      for (int k = 0; k < 3; k++) {
        _closest_elt_projected[3*j+k] = HUGE_VAL;
      }
    }

    /* Initialize the queue (add neighbours of the boundary cells) */


    /* fflush(stdout); */
    /* printf(" --- Fourth step2\n"); */

    _PDM_queue_t *q = _PDM_queue_init(2 * n_bound_cell_part);
    //_PDM_queue_t *q = _PDM_queue_init(2);

    for (int j = 0; j < n_bound_cell_part; j++) {
      int idx = bound_cell_idx[i] + j;
      int icell = bound_cell[idx];
      in_queue[icell] = -1;
      _closest_elt_distance[icell] = closest_elt_distance_bound[idx];
      _closest_elt_gnum[icell] = closest_elt_gnum_bound[idx];
      for (int k = 0; k < 3; k++) {
        _closest_elt_projected[3*icell+k] = closest_elt_projected_bound[3*idx+k];
      }
    }

    for (int j = 0; j < n_bound_cell_part; j++) {
      int idx = bound_cell_idx[i] + j;
      int icell = bound_cell[idx];
      for (int k = _cell_cell_idx[icell]; k < _cell_cell_idx[icell+1]; k++) {
        int icell1 = _cell_cell[k] - 1;
        if (icell1 > 0) {
          if (in_queue[icell1] == 0) {
            _PDM_queue_push (q, icell1);
            in_queue[icell1] = 1;
          }
        }
      }
    }

    /* fflush(stdout); */
    /* printf(" --- Fourth step3\n"); */

    int curr_cell;
    /* printf("\n\n**********************\n"); */
    while (!_PDM_queue_pull(q, &curr_cell)) {
      in_queue[curr_cell] = 0;
      double  _curr_closest_elt_distance = _closest_elt_distance[curr_cell];
      int _curr_closest_bound_cell = cell_to_closest_bound_cell[curr_cell];
      /* printf("*** queue : %d %d % 12.5e\n", curr_cell, _curr_closest_bound_cell, _curr_closest_elt_distance); */
      const double *_pt_coords = _cell_center + 3 * curr_cell;
      double new_closestPoint[3];
      int icell_src = -1;

      for (int k = _cell_cell_idx[curr_cell]; k < _cell_cell_idx[curr_cell+1]; k++) {
        int icell1 = _cell_cell[k] - 1;
        //        printf("icell 1 : %d\n", icell1+1);
        if (icell1 >= 0) {
          //if (in_queue[icell1] != 1) {
            int _icell1_closest_bound_cell = cell_to_closest_bound_cell[icell1];
            if (_icell1_closest_bound_cell != -1) {
              if (_curr_closest_bound_cell != _icell1_closest_bound_cell) {
                assert(_icell1_closest_bound_cell != -1);
                double * _coords_face_elt = sm_info_bound_cell[0] +
                  sm_info_bound_cell_idx[_icell1_closest_bound_cell];
                int    n_vtx_elt= sm_info_bound_cell_stride[0][_icell1_closest_bound_cell] / 3;

                double closestPoint[3];
                double minDist2;

                if (n_vtx_elt == 3) {

                  PDM_triangle_status_t status =
                    PDM_triangle_evaluate_position (_pt_coords,
                                                    _coords_face_elt,
                                                    closestPoint,
                                                    &minDist2,
                                                    NULL);

                  if (status == PDM_TRIANGLE_DEGENERATED) {
                    /* printf("degenerated tria\n"); */
                    continue;
                  }
                }

                else {

                  if (idebug) {
                    printf ("_pt_coords : %12.5e %12.5e %12.5e\n",
                            _pt_coords[0], _pt_coords[1], _pt_coords[2]);
                    printf ("_coords_face_elt %d %d %d :\n", n_vtx_elt,
                            sm_info_bound_cell_idx[_icell1_closest_bound_cell],
                            sm_info_bound_cell_idx[_icell1_closest_bound_cell]);

                    for (int k1 = 0; k1 < n_vtx_elt; k1++) {
                      printf (" / %12.5e %12.5e %12.5e /\n",
                              _coords_face_elt[3*k1],
                              _coords_face_elt[3*k1+1],
                              _coords_face_elt[3*k1+2]);
                    }
                    printf ("\n          *********\n");
                  }

                  PDM_polygon_status_t status =
                    PDM_polygon_evaluate_position (_pt_coords,
                                                   n_vtx_elt,
                                                   _coords_face_elt,
                                                   closestPoint,
                                                   &minDist2);

                  if (status == PDM_POLYGON_DEGENERATED) {
                    /* printf("degenerated poly\n"); */
                    continue;
                  }
                }
                if (minDist2 < _curr_closest_elt_distance) {
                  icell_src = icell1;
                  _curr_closest_elt_distance = minDist2;
                  _curr_closest_bound_cell = _icell1_closest_bound_cell;
                  for (int k1 = 0; k1 < 3; k1++) {
                    new_closestPoint[k1] = closestPoint[k1];
                  }
                }
              }
            }
            // }
        }
      }

      /* Update distance and add neighbours cell into the queue*/

      if (_curr_closest_bound_cell != cell_to_closest_bound_cell[curr_cell]) {
        /* printf ("Update distance\n"); */
        cell_to_closest_bound_cell[curr_cell] = _curr_closest_bound_cell;
        _closest_elt_distance[curr_cell] = _curr_closest_elt_distance;
        _closest_elt_gnum[curr_cell] = closest_elt_gnum_bound[_curr_closest_bound_cell];
        for (int k1 = 0; k1 < 3; k1++) {
          _closest_elt_projected[3*curr_cell+k1] = new_closestPoint[k1];
        }

        /* if (curr_cell == 1083) { */
        /*   printf("push pour %d :",curr_cell+1); */
        /* } */
        for (int k = _cell_cell_idx[curr_cell]; k < _cell_cell_idx[curr_cell+1]; k++) {

          int icell1 = _cell_cell[k] - 1;
          if (icell1 >= 0) {

            /* if (curr_cell == 1083) { */
            /*   printf(" %d", icell1+1); */
            /* } */
            if (icell1 != icell_src) {
              if (in_queue[icell1] == 0) {
                /* if (curr_cell == 1083) { */
                /*   printf("+"); */
                /* } */
                _PDM_queue_push (q, icell1);
                in_queue[icell1] = 1;
              }
            }
          }
        }
        /* if (curr_cell == 1083) { */
        /*   printf("\n"); */
        /* } */
      }
      /* printf("\n\n**********************\n"); */
    }

    /* fflush(stdout); */
    /* printf(" --- Fourth step4\n"); */

    free (in_queue);
    free (cell_to_closest_bound_cell);
    _PDM_queue_free (q);

  }

  free (n_face_sm);
  free (gnum_sm);
  free (sm_info_bound_cell[0]);
  free (sm_info_bound_cell_stride[0]);
  free (sm_info_bound_cell);
  free (sm_info_bound_cell_stride);
  free (sm_info_bound_cell_idx);

  for (int i = 0; i < n_part_sm; i++) {
    free (sm_face_coords[i]);
    free (sm_face_coords_stride[i]);
  }
  free (sm_face_coords);
  free (sm_face_coords_stride);

  free (closest_elt_distance_bound);
  free (closest_elt_projected_bound);
  free (closest_elt_gnum_bound);
  free (bound_cell);
  free (bound_cell_idx);

  PDM_timer_hang_on(dist->timer);
  e_t_elapsed = PDM_timer_elapsed(dist->timer);
  e_t_cpu     = PDM_timer_cpu(dist->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

  dist->times_elapsed[PROPAGATION] += e_t_elapsed - b_t_elapsed;
  dist->times_cpu[PROPAGATION]     += e_t_cpu - b_t_cpu;
  dist->times_cpu_u[PROPAGATION]   += e_t_cpu_u - b_t_cpu_u;
  dist->times_cpu_s[PROPAGATION]   += e_t_cpu_s - b_t_cpu_s;
  PDM_timer_resume(dist->timer);

  PDM_timer_hang_on(dist->timer);
  dist->times_elapsed[END] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[END]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[END]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[END]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  closest_elt_distance  Distance
 * \param [out]  closest_elt_projected Projected point coordinates
 * \param [out]  closest_elt_g_num     Global number of the closest element
 *
 */

void
PDM_dist_cellcenter_surf_get
(
 const int          id,
 const int          i_part,
       double      **closest_elt_distance,
       double      **closest_elt_projected,
       PDM_g_num_t **closest_elt_gnum
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  assert(_vol_mesh->closest_elt_distance != NULL);

  *closest_elt_distance = _vol_mesh->closest_elt_distance[i_part];
  *closest_elt_projected = _vol_mesh->closest_elt_projected[i_part];
  *closest_elt_gnum = _vol_mesh->closest_elt_gnum[i_part];
}


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_dist_cellcenter_surf_free
(
 const int id,
 const int partial
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  if (_vol_mesh->closest_elt_distance != NULL) {
    if (!partial) {
      for (int i = 0; i < _vol_mesh->n_part; i++) {
        if (_vol_mesh->closest_elt_distance[i] != NULL) {
          free (_vol_mesh->closest_elt_distance[i]);
        }
      }
    }
    free (_vol_mesh->closest_elt_distance);
  }

  if (_vol_mesh->closest_elt_projected != NULL) {
    if (!partial) {
      for (int i = 0; i < _vol_mesh->n_part; i++) {
        if (_vol_mesh->closest_elt_projected[i] != NULL) {
          free (_vol_mesh->closest_elt_projected[i]);
        }
      }
    }
    free (_vol_mesh->closest_elt_projected);
  }

  if (_vol_mesh->closest_elt_gnum != NULL) {
    if (!partial) {
      for (int i = 0; i < _vol_mesh->n_part; i++) {
        if (_vol_mesh->closest_elt_gnum[i] != NULL) {
          free (_vol_mesh->closest_elt_gnum[i]);
        }
      }
    }
    free (_vol_mesh->closest_elt_gnum);
  }

  if (_vol_mesh->n_cell != NULL) {
    free (_vol_mesh->n_cell);
  }

  if (_vol_mesh->n_face != NULL) {
    free (_vol_mesh->n_face);
  }

  if (_vol_mesh->n_vtx != NULL) {
    free (_vol_mesh->n_vtx);
  }

  if (_vol_mesh->cell_face_idx != NULL) {
    free (_vol_mesh->cell_face_idx);
  }

  if (_vol_mesh->cell_face != NULL) {
    free (_vol_mesh->cell_face);
  }

  if (_vol_mesh->cell_center != NULL) {
    free (_vol_mesh->cell_center);
  }

  if (_vol_mesh->cell_ln_to_gn != NULL) {
    free (_vol_mesh->cell_ln_to_gn);
  }

  if (_vol_mesh->face_vtx_idx != NULL) {
    free (_vol_mesh->face_vtx_idx);
  }

  if (_vol_mesh->face_vtx != NULL) {
    free (_vol_mesh->face_vtx);
  }

  if (_vol_mesh->face_cell != NULL) {
    for (int i = 0; i < _vol_mesh->n_part; i++) {
      if (_vol_mesh->face_cell[i] != NULL) {
        free (_vol_mesh->face_cell[i]);
      }
    }
    free (_vol_mesh->face_cell);
  }

  if (_vol_mesh->cell_cell != NULL) {
    for (int i = 0; i < _vol_mesh->n_part; i++) {
      if (_vol_mesh->cell_cell[i] != NULL) {
        free (_vol_mesh->cell_cell[i]);
      }
    }
    free (_vol_mesh->cell_cell);
  }

  if (_vol_mesh->face_ln_to_gn != NULL) {
    free (_vol_mesh->face_ln_to_gn);
  }

  if (_vol_mesh->vtx_ln_to_gn != NULL) {
    free (_vol_mesh->vtx_ln_to_gn);
  }

  if (_vol_mesh->coords != NULL) {
    free (_vol_mesh->coords);
  }

  free (_vol_mesh);

  if (dist->surf_mesh != NULL) {
    PDM_surf_mesh_free (dist->surf_mesh);
  }

  free (dist->timer);

  free (dist);

  PDM_Handles_handle_free (_dists, id, PDM_FALSE);

  const int n_dists = PDM_Handles_n_get (_dists);

  if (n_dists == 0) {
    _dists = PDM_Handles_free (_dists);
  }

}


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_dist_cellcenter_surf_dump_times
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  double t1 = dist->times_elapsed[END] - dist->times_elapsed[BEGIN];
  double t2 = dist->times_cpu[END] - dist->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (dist->times_elapsed, t_elaps_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (dist->times_cpu, t_cpu_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  int rank;
  PDM_MPI_Comm_rank (dist->comm, &rank);

  if (rank == 0) {

    PDM_printf( "Cell center distance timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "Cell center distance timer : Distance of the cell centers of the first thickness (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[FIRST_THICKNESS],
                t_cpu_max[FIRST_THICKNESS]);
    PDM_printf( "Cell center distance timer : Distance of the other cell centers (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[PROPAGATION],
                t_cpu_max[PROPAGATION]);
    PDM_printf_flush();
  }

}

#ifdef	__cplusplus
}
#endif
