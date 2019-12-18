
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_morton.h"
#include "pdm_handles.h"
#include "pdm_para_octree.h"
#include "pdm_timer.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define NTIMER 11

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  BUILD_ORDER_POINTS            = 1,
  BUILD_BLOCK_PARTITION         = 2,
  BUILD_LOCAL_NODES             = 3,
  BUILD_LOCAL_NEIGHBOURS_STEP1  = 4,
  BUILD_LOCAL_NEIGHBOURS_STEP2  = 5,
  BUILD_LOCAL_NEIGHBOURS_STEP3  = 6,
  BUILD_LOCAL_NEIGHBOURS        = 7,
  BUILD_DISTANT_NEIGHBOURS      = 8,
  BUILD_TOTAL                   = 9,
  END                           = 10,

} _ol_timer_step_t;


/**
 * \struct _heap_t
 * \brief  Heap used to recursively subdivide nodes
 *
 */

typedef struct  {

  int   top;                  /*!< Top of head  */
  int   size;                 /*!< Size of heap */
  PDM_morton_code_t *codes;   /*!< Morton codes */
  int *range;                 /*!< Points range */
  int *n_points;              /*!< Points number */
  int   max_top;

} _heap_t;


/**
 * \struct _l_octant_t
 * \brief  Define a list of octants
 *
 */

typedef struct  {

  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */

  PDM_morton_code_t *codes;        /*!< Morton codes */

  int  *n_points;          /*!< Number of points in octant*/
  int  *range;             /*!< Start index of point list for each octant */

  int   *neighbour_idx;
  int   *neighbours;               /*!< rank + id_node size = 2 * n_nodes */
  int   dim;

} _l_octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 *
 */

typedef struct  {

  double  global_extents[6];     /*!< Extents of current process */
  int     depth_max;             /*!< Maximum depth of the three */
  int     points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double      s[3];           /*!< Translation for the normalization */
  double      d[3];           /*!< Dilatation for the normalization */

  int     n_point_clouds;        /*!< Number of point cloud */

  PDM_g_num_t    t_n_points;         /*!< total number of points */
  int            n_points;           /*!< Number of points in each cloud */
  double *points;                    /*!< Point coordinates */
  int *points_icloud;                /*!< Point cloud */
  PDM_g_num_t *points_gnum;          /*!< Point global number */
  PDM_morton_code_t  *points_code;   /*!< Morton codes */

  _l_octant_t *octants;       /*!< list of octants */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */

  int n_part_boundary_elt;    /*!< Number of partitioning boundary element */
  int *part_boundary_elt_idx; /*!< Index for part_boundary_elt (size=\ref n_part_boundary_elt + 1 */
  int *part_boundary_elt;     /*!< Partitioning boundary elements description (proc number + element number) */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */

} _octree_t;



/**
 * \struct _neighbours_tmp_t
 * \brief  Define a temporary neighbour structure
 *
 */


typedef struct  {

  int n_neighbour[6];     /*!< Number of neighbours in the arrays  */
  int s_neighbour[6];     /*!< Size of arrays */
  int *neighbours[6];     /*!< Arrays */

} _neighbours_tmp_t;


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees    = NULL;

//static const double _eps_default  = 1.e-12;

static const PDM_morton_int_t max_morton_level = 15;
//static const int max_morton_level = 2;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Create a heap
 *
 * \param [in]  Size    Size of heap
 *
 * \return   a new heap
 *
 */

static _heap_t *
_heap_create
(
const int size
)
{
  _heap_t *heap = malloc(sizeof(_heap_t));
  heap->top = 0;
  heap->max_top = 0;
  heap->size = size;
  heap->codes = malloc(sizeof(PDM_morton_code_t) * size);
  heap->range =  malloc(sizeof(int) * size);
  heap->n_points =  malloc(sizeof(int) * size);
  return heap;
}


/**
 *
 * \brief Free a heap
 *
 * \param [in]  heap   Heap to free
 *
 * \return NULL
 *
 */

static _heap_t *
_heap_free
(
_heap_t *heap
)
{
  free (heap->codes);
  free (heap->range);
  free (heap->n_points);
  free (heap);
  return NULL;
}


/**
 *
 * \brief Push a new element in the heap
 *
 * \param [inout]  heap      Heap
 * \param [in]     code      Morton code
 * \param [in]     range     Range
 * \param [in]     n_points  Number of points
 *
 * \return  1 if pushed 0 otherwise
 *
 */

static int
_heap_push
(
 _heap_t *heap,
 const PDM_morton_code_t code,
 const int range,
 const int n_points
)
{
  if (heap->top >= heap->size) {
    return 0;
  }
  int idx = heap->top;
  PDM_morton_copy (code, &(heap->codes[idx]));
  heap->range[idx] = range;
  heap->n_points[idx] = n_points;
  heap->top++;
  heap->max_top = PDM_MAX(heap->top, heap->max_top);
  return 1;
}


/**
 *
 * \brief Pull top element of the heap
 *
 * \param [inout]  heap      Heap
 * \param [out]    code      Morton code
 * \param [out]    range     Range
 * \param [out]    n_points  Number of points
 *
 * \return  1 if pulled 0 otherwise
 *
 */

static int
_heap_pull
(
 _heap_t *heap,
 PDM_morton_code_t *code,
 int *range,
 int *n_points
)
{
  heap->top--;
  if (heap->top < 0) {
    return 0;
  }
  int idx = heap->top;
  PDM_morton_copy (heap->codes[idx], code);
  *range = heap->range[idx];
  *n_points = heap->n_points[idx];
  return 1;
}

/**
 *
 * \brief Neighbour
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static PDM_para_octree_direction_t
_inv_direction
(
 PDM_para_octree_direction_t direc
)
{
  if (direc == PDM_BOTTOM) {
    return PDM_UP;
  }
  else if (direc == PDM_UP) {
    return PDM_BOTTOM;
  }
  else if (direc == PDM_SOUTH) {
    return PDM_NORTH;
  }
  else if (direc == PDM_NORTH) {
    return PDM_SOUTH;
  }
  else if (direc == PDM_WEST) {
    return PDM_EAST;
  }
  else if (direc == PDM_EAST) {
    return PDM_WEST;
  }
  else {
    abort();
  }

}

/**
 *
 * \brief Neighbour
 *
 * \param [inout]   octants     Octants
 *
 * \return neighbour or NULL ()
 *
 */

static PDM_morton_code_t *
_neighbour
(
 PDM_morton_code_t code,
 PDM_para_octree_direction_t direction
)
{
  const int dim = direction / 2;
  const int _direction = 2 * (direction % 2) - 1;

  PDM_morton_code_t *neighbour = NULL;

  if (((_direction > 0) && (code.X[dim] < (pow(2,code.L) - 1))) ||
      ((_direction < 0) && (code.X[dim] > 0))) {

    neighbour = malloc(sizeof(PDM_morton_code_t));

    neighbour->L = code.L;
    neighbour->X[0] = code.X[0];
    neighbour->X[1] = code.X[1];
    neighbour->X[2] = code.X[2];

    neighbour->X[dim] = code.X[dim] + _direction;
  }

  return neighbour;
}


/**
 *
 * \brief Free octants
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static void
_octants_purge
(
 _l_octant_t *octants
)
{
  octants->n_nodes_max = 0;
  octants->n_nodes     = 0;

  if (octants->codes != NULL) {
    free (octants->codes);
  }

  if (octants->n_points != NULL) {
    free (octants->n_points);
  }

  if (octants->range != NULL) {
    free (octants->range);
  }

  if (octants->neighbour_idx != NULL) {
    free (octants->neighbour_idx);
  }

  if (octants->neighbours != NULL) {
    free (octants->neighbours);
  }
}

/**
 *
 * \brief Free octants
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static _l_octant_t *
_octants_free
(
 _l_octant_t *octants
)
{

  _octants_purge (octants);

  free(octants);
  return NULL;
}


/**
 *
 * \brief Initialize list of octants
 *
 * \param [inout]   octants     Octants
 * \param [in]      octant_dim  Dimension of an octant
 * \param [in]      init_size   Initial size of octants
 *
 */

static void
_octants_init
(
 _l_octant_t *octants,
 const int   octant_dim,
 const int   init_size
)
{
  octants->n_nodes_max = init_size;
  octants->n_nodes     = 0;

  octants->codes    = malloc (sizeof(PDM_morton_code_t) * octants->n_nodes_max);
  octants->n_points = malloc (sizeof(int) * octants->n_nodes_max);
  octants->range    = malloc (sizeof(int) * (octants->n_nodes_max+1));

  octants->neighbour_idx = NULL;
  octants->neighbours    = NULL;
  octants->dim = octant_dim;
}


/**
 *
 * \brief Check size of the size of a list of octants
 *
 * \param [in]   octants       Octants
 * \param [in]   n_free_node   Number of required fre nodes
 *
 */

static int
_octants_check_alloc
(
 _l_octant_t *octants,
 const int n_free_node
)
{
  int is_realloc = 0;
  if (octants->n_nodes + n_free_node > octants->n_nodes_max) {

    octants->n_nodes_max *= 2;

    octants->codes    = realloc (octants->codes,
                                 sizeof(PDM_morton_code_t) * octants->n_nodes_max);
    octants->n_points = realloc (octants->n_points,
                                 sizeof(int) * octants->n_nodes_max);
    octants->range = realloc (octants->range,
                              sizeof(int) * (octants->n_nodes_max+1));
    octants->neighbour_idx = NULL;
    octants->neighbours    = NULL;

    is_realloc = 1;
  }
  return is_realloc;
}


/**
 *
 * \brief Push back a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_push_back
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range
)
{

  _octants_check_alloc (octants, 1);

  const int idx = octants->n_nodes;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->n_nodes += 1;

}


/**
 *
 * \brief Push front a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_push_front
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range
)
{

  _octants_check_alloc (octants, 1);

  for (int i = octants->n_nodes; i > 0; i--) {

    PDM_morton_copy (octants->codes[i - 1], octants->codes + i);

    octants->n_points[i] =  octants->n_points[i-1];

    octants->range[i] = octants->range[i-1];
  }

  const int idx = 0;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->n_nodes += 1;

}


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
 *
 * \brief Build minimal octree between two octants
 *
 * \param [in]     octree    Current octree
 * \param [in]     code      Morton code
 * \param [inout]  extents   Extents associated to the Morton code
 *
 */

/* static void */
/* g_extents */
/* ( */
/*  _octree_t *octree, */
/*  PDM_morton_code_t code, */
/*  double    extents[] */
/* ) */
/* { */
/*   for (int i = 0; i < octree->dim; i++) { */
/*     extents[i] = */
/*       ((double) code.X[i]/((double) pow(2,code.L)))* octree->d[i] + octree->s[i]; */
/*     extents[octree->dim + i] = */
/*       (((double) code.X[i] + 1)/((double) pow(2,code.L))) * octree->d[i] + octree->s[i]; */
/*   } */
/* } */

/**
 *
 * \brief Remove duplicates
 *
 * \param [in]  octants List of octants
 *
 * \return octants without duplicates
 *
 */

static _l_octant_t *
_remove_duplicates
(
 _l_octant_t *octants
)
{
  PDM_morton_code_t *_codes = octants->codes;
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (r_octants, dim, octants->n_nodes);

  PDM_morton_code_t prev_code;

  prev_code.L = -1;
  prev_code.X[0] = 0;
  prev_code.X[1] = 0;
  prev_code.X[2] = 0;

  for (int i = 0; i < octants->n_nodes; i++) {

    if (_codes[i].L == prev_code.L) {
      if ((prev_code.X[0] == _codes[i].X[0]) &&
          (prev_code.X[1] == _codes[i].X[1]) &&
          (prev_code.X[2] == _codes[i].X[2])) {

        break;
      }
    }

    prev_code.L    = _codes[i].L;
    prev_code.X[0] = _codes[i].X[0];
    prev_code.X[1] = _codes[i].X[1];
    prev_code.X[2] = _codes[i].X[2];

    _octants_push_back (r_octants,
                        _codes[i],
                        0,
                        0);

  }

  return r_octants;
}

/**
 *
 * \brief Removing overlaps from a sorted lis of octants
 *
 * \param [inout]  octants A lis of octants
 *
 */

static _l_octant_t *
_linearize
(
 _l_octant_t *octants
)
{
  PDM_morton_code_t *_codes = octants->codes;
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (r_octants, dim, octants->n_nodes);

  if  (octants->n_nodes > 0) {
    for (int i = 0; i < octants->n_nodes - 1; i++) {

      if (!PDM_morton_ancestor_is (_codes[i], _codes[i+1])) {
        _octants_push_back (r_octants,
                            _codes[i],
                            0,
                            0);
      }

    }

    _octants_push_back (r_octants,
                        _codes[octants->n_nodes-1],
                        0,
                        0);
  }

  return r_octants;
}


/**
 *
 * \brief Constructing a minimal linear octree between two octants
 *
 * \param [in]  a     Morton code a
 * \param [in]  b     Morton code b
 *
 * \return octants The minimal linear octree between a and b or NULL if a >= b
 *
 */

static _l_octant_t *
_complete_region
(
 PDM_morton_code_t a,
 PDM_morton_code_t b
)
{
  const int dim = 3;

  _l_octant_t *r_octants = NULL;

  if (PDM_morton_a_gt_b (b, a)) {

    _l_octant_t *w_octants = malloc(sizeof(_l_octant_t));
    r_octants = malloc(sizeof(_l_octant_t));

    _octants_init (w_octants, dim, 4);
    _octants_init (r_octants, dim, 4);

    /* printf("_complete_region\n"); */

    /* printf("a_d\n"); */
    /* PDM_morton_dump(3, a); */
    /* printf("a_f\n"); */

    /* printf("b_d\n"); */
    /* PDM_morton_dump(3, b); */
    /* printf("b_f\n"); */

    PDM_morton_code_t nca;
    PDM_morton_nearest_common_ancestor (a, b, &nca);

    const int n_child = 8;

    int  size = PDM_morton_max_level * 8;
    _heap_t *heap = _heap_create (size);

    PDM_morton_code_t children[8];
    PDM_morton_get_children(dim,
                            nca,
                            children);

    for (int i = n_child - 1; i >= 0; i--) {
      int is_pushed = _heap_push (heap,
                                  children[i],
                                  0,
                                  0);
      if (!is_pushed) {
        printf ("Internal error PDM_para_octree 1 : heap is full\n");
        exit(1);
      }
    }

    PDM_morton_code_t code;
    int range;
    int n_points;

    while (_heap_pull (heap, &code, &range, &n_points)) {
      if (PDM_morton_a_gt_b (code, a) &&
          PDM_morton_a_gt_b (b, code) &&
          !PDM_morton_ancestor_is (code, b)) {

        _octants_push_back (r_octants,
                            code,
                            0,
                            0);
      }

      else if ((PDM_morton_ancestor_is (code, b) ||
                PDM_morton_ancestor_is (code, a)) &&
               !((code.X[0] == a.X[0]) &&
                 (code.X[1] == a.X[1]) &&
                 (code.X[2] == a.X[2]) &&
                 (code.L == a.L)) &&
               !((code.X[0] == b.X[0]) &&
                 (code.X[1] == b.X[1]) &&
                 (code.X[2] == b.X[2]) &&
                 (code.L == b.L))) {

        PDM_morton_get_children(dim,
                                code,
                                children);

        for (int i = n_child - 1; i >= 0; i--) {
          int is_pushed = _heap_push (heap,
                                      children[i],
                                      0,
                                      0);

          if (!is_pushed) {
            printf ("Internal error PDM_para_octree 2 : heap is full\n");
            exit(1);
          }
        }
      }
    }

    _octants_free (w_octants);

    _heap_free (heap);
  }

  return r_octants;
}



/**
 *
 * \brief Redistribute octants
 *
 * \param [inout]  L             Distributed list of octants
 * \param [in]     morton_index  Morton index
 * \param [in]     comm          MPI communicator
 *
 */

static void
_distribute_octants
(
 _l_octant_t       *L,
 PDM_morton_code_t *morton_index,
 PDM_MPI_Comm       comm
)
{
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int *send_count = malloc(sizeof(int) * n_ranks);
  size_t *send_shift = malloc(sizeof(size_t) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  size_t *recv_shift = malloc(sizeof(size_t) * (n_ranks+1));

  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  int irank = 0;
  for (int i = 0; i < L->n_nodes; i++) {
    if (PDM_morton_a_ge_b (L->codes[i], morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search (n_ranks - (irank + 1),
                                             L->codes[i],
                                             morton_index + irank + 1);
    }
    send_count[irank] += L->dim + 1;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  PDM_morton_int_t *send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < L->n_nodes; i++) {

    if (PDM_morton_a_ge_b (L->codes[i], morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search(n_ranks - (irank + 1),
                                            L->codes[i],
                                            morton_index + irank + 1);
    }

    int shift = send_shift[irank] + send_count[irank];
    send_codes[shift++] = L->codes[i].L;

    for (int j = 0; j < L->dim; j++) {
      send_codes[shift++] = L->codes[i].X[j];
    }

    send_count[irank] += L->dim + 1;
  }

  PDM_morton_int_t * recv_codes = malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */

  PDM_MPI_Alltoallv_l(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                      recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                      comm);

  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  /* - tri des codes recus */

  const int _dim = L->dim;
  _octants_purge (L);
  _octants_init (L, _dim, recv_shift[n_ranks]/(_dim + 1));

  int idx = 0;
  for (int i = 0; i < recv_shift[n_ranks]/4; i++) {
    PDM_morton_code_t _code;
    _code.L = recv_codes[idx++];
    for (int j = 0; j < L->dim; j++) {
     _code.X[j] = recv_codes[idx++];
    }
    _octants_push_back (L,
                        _code,
                        0,
                        0);
  }

  free (recv_shift);
  free (recv_codes);

  PDM_morton_local_sort (L->n_nodes, L->codes);

}


/**
 *
 * \brief Constructing a complete linear octree from partial set of octants
 *
 * \param [in]  L     Distributed list of octants
 * \param [in]  comm  MPI Communicator
 *
 * \return octants The complete linear octree
 *
 */

static _l_octant_t *
_complete_octree
(
 _l_octant_t *L,
 PDM_MPI_Comm comm
)
{
  const int dim = 3;

  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (comm, &rank);

  /* Remove duplicates */

  _l_octant_t *L1 = _remove_duplicates (L);

  /* Linearize */

  _l_octant_t *L2 = _linearize (L1);

  _octants_free (L1);

  PDM_g_num_t _n_nodes_global = 0;
  PDM_g_num_t _n_nodes_local =  L2->n_nodes;

  PDM_MPI_Allreduce (&_n_nodes_local, &_n_nodes_global, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  _l_octant_t *R = malloc(sizeof(_l_octant_t));
  _octants_init (R, dim, PDM_MAX(L2->n_nodes, 1));

  if (_n_nodes_global > 0) {

    PDM_morton_code_t *L2_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));

    int *order = malloc (sizeof(int) * L2->n_nodes);
    int *weight = malloc (sizeof(int) * L2->n_nodes);

    for (int i = 0; i < L2->n_nodes; i++) {
      weight[i] = 1;
      order[i] = i;
    }

    PDM_morton_int_t max_level = 0;
    for (int i = 0; i < L2->n_nodes; i++) {
    max_level = PDM_MAX (L2->codes[i].L, max_level);
    }

    PDM_morton_int_t max_max_level;
    PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                      PDM_MPI_UNSIGNED, PDM_MPI_MAX, comm);

    /* printf("n L2 : %d %d\n" , rank, L2->n_nodes); */
    /* printf("\n-------------\nL2 avant d\n"); */
    /* for (int i = 0; i < L2->n_nodes; i++) { */
    /*   PDM_morton_dump (3, L2->codes[i]); */
    /* } */
    PDM_morton_ordered_build_rank_index (dim,
                                         max_max_level,
                                         L2->n_nodes,
                                         L2->codes,
                                         weight,
                                         L2_morton_index,
                                         comm);

    /* printf("\nL2 morton index d\n"); */
    /* for (int i = 0; i < n_ranks + 1; i++) { */
    /*   PDM_morton_dump (3,  L2_morton_index[i]); */
    /* } */
    /* printf("L2 morton, index f\n"); */

    free(weight);
    free(order);

    _distribute_octants (L2, L2_morton_index, comm);

    free (L2_morton_index);

    /* printf("n L2 : %d %d\n" , rank, L2->n_nodes); */
    /* printf("\nL2 d\n"); */
    /* for (int i = 0; i < L2->n_nodes; i++) { */
    /*   PDM_morton_dump (3, L2->codes[i]); */
    /* } */
    /* printf("L2 f\n--------------\n"); */

  /* PDM_MPI_Barrier(comm); */
  /* exit(1); */

    int *rank_n_nodes = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Allgather(&L2->n_nodes, 1, PDM_MPI_INT,
                      rank_n_nodes, 1, PDM_MPI_INT,
                      comm);

    int first_rank = 0;
    while (first_rank < n_ranks-1
           && rank_n_nodes[first_rank] == 0) {
      first_rank++;
    }

    int last_rank = n_ranks-1;
    while (last_rank > 0
           && rank_n_nodes[last_rank] == 0) {
      last_rank--;
    }

    int next_rank = rank + 1;
    while (next_rank < n_ranks-1
           && rank_n_nodes[next_rank] == 0) {
      next_rank++;
    }

    int prev_rank = rank - 1;
    while (prev_rank > 0
           && rank_n_nodes[prev_rank] == 0) {
      prev_rank--;
    }

    /* printf ("[%d] first last next prev : %d %d %d %d\n", rank, */
    /*         first_rank, last_rank, next_rank, prev_rank); */

    if (rank == first_rank && rank_n_nodes[rank] > 0) {
      PDM_morton_code_t root_DFD;

      root_DFD.L = max_morton_level;
      root_DFD.X[0] = 0;
      root_DFD.X[1] = 0;
      root_DFD.X[2] = 0;

      PDM_morton_code_t FINA;

      PDM_morton_nearest_common_ancestor (root_DFD,
                                          L2->codes[0],
                                          &FINA);

      PDM_morton_code_t child[8];
      PDM_morton_get_children(dim,
                              FINA,
                              child);

      _octants_push_front (L2,
                           child[0],
                           0,
                           0);
    }


    if (rank == last_rank && rank_n_nodes[rank] > 0) {
      PDM_morton_code_t root_DLD;

      root_DLD.L = max_morton_level;
      root_DLD.X[0] = (1u << max_morton_level) - 1u;
      root_DLD.X[1] = (1u << max_morton_level) - 1u;
      root_DLD.X[2] = (1u << max_morton_level) - 1u;

      PDM_morton_code_t FINA;
      PDM_morton_nearest_common_ancestor (root_DLD,
                                          L2->codes[L2->n_nodes -1],
                                          &FINA);

      PDM_morton_code_t child[8];
      PDM_morton_get_children(dim,
                              FINA,
                              child);

      _octants_push_back (L2,
                          child[7],
                          0,
                          0);

    }

    unsigned int sbuff[4];
    unsigned int rbuff[4];
    PDM_MPI_Request srequest;
    PDM_MPI_Request rrequest;

    if (rank < last_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Irecv ((void *) rbuff,
                     4,
                     PDM_MPI_UNSIGNED,
                     next_rank,
                     0,
                     comm,
                     &rrequest);

    }

    if (rank > first_rank && rank_n_nodes[rank] > 0) {

      assert (L2->n_nodes > 0);
      sbuff[0] = L2->codes[0].L;
      sbuff[1] = L2->codes[0].X[0];
      sbuff[2] = L2->codes[0].X[1];
      sbuff[3] = L2->codes[0].X[2];

      PDM_MPI_Issend ((void *) sbuff,
                      4,
                      PDM_MPI_UNSIGNED,
                      prev_rank,
                      0,
                      comm,
                      &srequest);


    }

    if (rank < last_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Wait (&rrequest);
      PDM_morton_code_t code;

      code.L = rbuff[0];
      code.X[0] = rbuff[1];
      code.X[1] = rbuff[2];
      code.X[2] = rbuff[3];

      _octants_push_back (L2,
                          code,
                          0,
                          0);

    }

    if (rank > first_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Wait (&srequest);

    }

    for (int i = 0; i < L2->n_nodes - 1; i++) {
      _l_octant_t *A = _complete_region (L2->codes[i], L2->codes[i+1]);

      _octants_push_back (R,
                          L2->codes[i],
                          0,
                          0);

      if (A != NULL) {
        for (int j = 0; j < A->n_nodes; j++) {
          _octants_push_back (R,
                              A->codes[j],
                              0,
                              0);
        }
        _octants_free (A);
      }
    }

    if (rank == last_rank  && rank_n_nodes[rank] > 0) {
      _octants_push_back (R,
                          L2->codes[L2->n_nodes-1],
                          0,
                          0);
    }

    _octants_free (L2);

    free (rank_n_nodes);
  }

  else {
    if (rank == n_ranks - 1) {

      PDM_morton_code_t _code;
      _code.L = 0;
      _code.X[0] = 0;
      _code.X[1] = 0;
      _code.X[2] = 0;

      _octants_push_back (R,
                          _code,
                          0,
                          0);

    }

  }

  /* printf("fin complete_octree\n"); */
  /* fflush(stdout); */

  return R;
}



/**
 *
 * \brief Distribute points
 *
 * \param [in]   id                 Identifier
 *
 */

static void
_distribute_points
(
 int *n_points,
 double **points,
 int **points_icloud,
 PDM_g_num_t **points_gnum,
 PDM_morton_code_t **points_code,
 PDM_morton_code_t *morton_index,
 const PDM_MPI_Comm comm,
 const int dim,
 const PDM_morton_int_t max_level,
 const double *global_extents
)
{
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int _n_points = *n_points;

  double *__points = *points;
  int *__points_icloud = *points_icloud;
  PDM_g_num_t *__points_gnum = *points_gnum;
  PDM_morton_code_t *__points_code = *points_code;

  int *c_rank = malloc (_n_points * sizeof(int));

  for (int i = 0; i < _n_points; i++) {
    size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                __points_code[i],
                                                morton_index);
    c_rank[i] = (int) _c_rank;
  }

  int *send_count = malloc (n_ranks * sizeof (int));
  int *recv_count = malloc (n_ranks * sizeof (int));
  int *send_shift = malloc ((n_ranks + 1) * sizeof (int));
  int *recv_shift = malloc ((n_ranks + 1) * sizeof (int));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    send_count[c_rank[i]] += dim;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  double *send_coords = malloc (send_shift[n_ranks] * sizeof(double));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++)
    send_count[rank_id] = 0;

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    for (int j = 0; j < dim; j++)
      send_coords[shift + j] = __points[i*dim + j];
    send_count[rank_id] += dim;
  }

  double *recv_coords = malloc (recv_shift[n_ranks] * sizeof(double));

  /* Exchange coords between processes */

  PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                    recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                    comm);

  free(send_coords);

  /* Build send and receive buffers */

  for (int rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
    send_shift[rank_id] = send_shift[rank_id]/dim;
    recv_shift[rank_id] = recv_shift[rank_id]/dim;
  }

  int *send_points_icloud = malloc (send_shift[n_ranks] * sizeof(int));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    recv_count[rank_id] = recv_count[rank_id]/dim;
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    send_points_icloud[shift] = __points_icloud[i];
    send_count[rank_id] += 1;
  }

  int *recv_points_icloud = malloc (recv_shift[n_ranks] * sizeof(int));

  /* Exchange points_icloud between processes */

  PDM_MPI_Alltoallv(send_points_icloud, send_count, send_shift, PDM_MPI_INT,
                    recv_points_icloud, recv_count, recv_shift, PDM_MPI_INT,
                    comm);

  free(send_points_icloud);


  /* Build send and receive buffers : points_gnum*/

  PDM_g_num_t *send_points_gnum =
    malloc (send_shift[n_ranks] * sizeof(PDM_g_num_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    send_points_gnum[shift] = __points_gnum[i];
    send_count[rank_id] += 1;
  }

  free (c_rank);

  PDM_g_num_t *recv_points_gnum =
    malloc (recv_shift[n_ranks] * sizeof(PDM_g_num_t));

  /* Exchange points_gnum between processes */

  PDM_MPI_Alltoallv(send_points_gnum, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                    recv_points_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                    comm);

  free(send_points_gnum);

  _n_points = recv_shift[n_ranks];

  free (send_count);
  free (recv_count);
  free (send_shift);
  free (recv_shift);

  __points = realloc (__points, sizeof(double)  * 3 * _n_points);

  __points_icloud =
    realloc (__points_icloud, sizeof(int) * _n_points);

  __points_gnum =
    realloc (__points_gnum, sizeof(PDM_g_num_t) * _n_points);

  /* Re-encode points */

  __points_code = realloc (__points_code,
                           sizeof(PDM_morton_code_t) * _n_points);

  double d[3];
  double s[3];

  PDM_morton_encode_coords(dim,
                           max_level,
                           global_extents,
                           _n_points,
                           recv_coords,
                           __points_code,
                           d,
                           s);

  int *order = malloc (sizeof(int) * _n_points);

  for (int i = 0; i < _n_points; i++) {
    order[i] = i;
  }

  PDM_morton_local_order(_n_points, __points_code, order);

  for (int i = 0; i < _n_points; i++) {
    __points_icloud[i] = recv_points_icloud[order[i]];
    __points_gnum[i] = recv_points_gnum[order[i]];
    for (int j = 0; j < dim; j++) {
      __points[dim*i+j] = recv_coords[dim*order[i]+j];
    }
  }

  free (recv_points_icloud);
  free (recv_points_gnum);
  free (recv_coords);

  PDM_morton_code_t *_points_code =
    malloc (sizeof(PDM_morton_code_t) * _n_points);

  for (int i = 0; i < _n_points; i++) {
    _points_code[i].L = __points_code[order[i]].L;
    _points_code[i].X[0] = __points_code[order[i]].X[0];
    _points_code[i].X[1] = __points_code[order[i]].X[1];
    _points_code[i].X[2] = __points_code[order[i]].X[2];
  }

  free (__points_code);
  free (order);

  *points_code = _points_code;

  *points = __points;
  *points_icloud = __points_icloud;
  *points_gnum = __points_gnum;

  *n_points = _n_points;
}


/**
 *
 * \brief Partitioning octants into large contiguous blocks. The list of octants
 *        is redistributed
 *
 * \param [in]  octant_list  a list of distributed octants,
 *                           octant_list is not redistributed at the end
 *
 * \return block_octants  A list of distributed blocks
 *
 */

static _l_octant_t *
_block_partition
(
 _l_octant_t *octant_list,
 const PDM_MPI_Comm comm,
 PDM_morton_code_t **G_morton_index
)
{

  /* Complete region */

  _l_octant_t *T = NULL;

  int max_level = -1;
  int min_level = 32;

  int comm_rank;
  PDM_MPI_Comm_rank(comm, &comm_rank);

  /* printf("\n_block_partition : octant_list %d : d\n", comm_rank); */
  /* if (octant_list!=NULL) { */
  /*   printf("\n_nodes : %d\n", octant_list->n_nodes); */
  /*   for (int i = 0; i < octant_list->n_nodes; i++) { */
  /*     PDM_morton_dump (3, octant_list->codes[i]); */
  /*   } */
  /* } */
  /* else { */
  /*   printf ("octant_list NULL\n"); */
  /* } */

  /* printf("_block_partition :octant_list f\n\n"); */

  if (octant_list->n_nodes > 1 ) {

    T = _complete_region (octant_list->codes[0],
                          octant_list->codes[octant_list->n_nodes-1]);

    if (T !=NULL) {
      for (int i = 0; i < T->n_nodes; i++) {
        max_level = PDM_MAX (T->codes[i].L, max_level);
        min_level = PDM_MIN (T->codes[i].L, min_level);
      }
    }
  }

  /* printf("\n_block_partition : complete_region %d :d\n", comm_rank); */
  /* if (T!=NULL) { */
  /*   printf("\n_nodes : %d\n", T->n_nodes); */
  /*   for (int i = 0; i < T->n_nodes; i++) { */
  /*     PDM_morton_dump (3, T->codes[i]); */
  /*   } */
  /* } */
  /* else { */
  /*   printf ("T NULL\n"); */
  /* } */
  /* printf("_block_partition : complete_region f\n\n"); */

  int max_max_level;
  PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, comm);

  /* Intialize C with coarse octants */

  _l_octant_t *C = malloc(sizeof(_l_octant_t));

  if (T != NULL) {

    _octants_init (C, T->dim, T->n_nodes);

    for (int i = 0; i < T->n_nodes; i++) {

      if (T->codes[i].L <= min_level) {
        _octants_push_back (C,
                            T->codes[i],
                            T->n_points[i],
                            T->range[i]);
      }
    }
  }

  else {
    _octants_init (C, octant_list->dim, 1);
  }

  /* Complete octree */

  /* printf("\n_block_partition : before complete_octree %d %d : d\n", comm_rank, C->n_nodes); */
  /* for (int i = 0; i < C->n_nodes; i++) { */
  /*   PDM_morton_dump (3, C->codes[i]); */
  /* } */
  /* printf("_block_partition : before complete_octree f\n\n"); */

  _l_octant_t *G = _complete_octree (C, comm);


  double vol = 0;
  for (int i = 0; i < G->n_nodes; i++) {
    vol += ((1./pow(2, G->codes[i].L)) *
            (1./pow(2, G->codes[i].L)) *
            (1./pow(2, G->codes[i].L)));
    G->range[i+1] =
      G->range[i] +
      G->n_points[i];
  }
  double total_vol;
  PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  if ( (PDM_ABS(total_vol - 1.)>= 1e-15)) {
    printf("Erreur volume different de 1 apres complete_octree : %12.5e\n", total_vol);
    for (int i = 0; i < G->n_nodes; i++) {
      PDM_morton_dump (3, G->codes[i]);
    }
  }

  assert (PDM_ABS(total_vol - 1.) < 1e-15);


  /* printf("\n_block_partition : after complete_octree %d %d : d\n", comm_rank, C->n_nodes); */
  /* for (int i = 0; i < G->n_nodes; i++) { */
  /*   PDM_morton_dump (3, G->codes[i]); */
  /* } */
  /* printf("_block_partition : after complete_octree f\n\n"); */

  /* PDM_MPI_Barrier (comm); */
  /* exit(1); */
  _octants_free (C);

  if (T != NULL) {
    _octants_free (T);
  }

  /*
   * Compute weight
   */

  /* - exchange codes to ranks (weight per rank)*/

  int n_ranks;
  PDM_MPI_Comm_size(comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank(comm, &rank);

  PDM_morton_int_t *code_buff = malloc (sizeof(PDM_morton_int_t) * (G->dim + 1));
  PDM_morton_int_t *rank_buff = malloc (sizeof(PDM_morton_int_t) * n_ranks * (G->dim + 1));
  int *n_nodes_rank = malloc (sizeof(int) * n_ranks);

  for (int i = 0; i < G->dim + 1; i++) {
    code_buff[i] = 0;
  }

  if ( G->n_nodes > 0) {
    code_buff[0] = G->codes[0].L;
    for (int i = 0; i < G->dim; i++) {
      code_buff[i+1] =  G->codes[0].X[i];
    }
  }

  PDM_MPI_Allgather (&(G->n_nodes), 1, PDM_MPI_INT,
                     n_nodes_rank,  1, PDM_MPI_INT,
                     comm);

  PDM_MPI_Allgather (code_buff, G->dim + 1, PDM_MPI_UNSIGNED,
                     rank_buff, G->dim + 1, PDM_MPI_UNSIGNED,
                     comm);

  int n_active_ranks = 0;
  for (int i = 0; i < n_ranks; i++) {
    if (n_nodes_rank[i] > 0) {
      n_active_ranks++;
    }
  }

  //assert (n_active_ranks > 0);

  PDM_morton_code_t *rank_codes = malloc (sizeof(PDM_morton_code_t) * n_active_ranks);
  int *active_ranks = malloc (sizeof(int) * n_active_ranks);

  n_active_ranks = 0;
  for (int i = 0; i < n_ranks; i++) {
    if (n_nodes_rank[i] > 0) {
      active_ranks[n_active_ranks] = i;
      rank_codes[n_active_ranks].L = rank_buff[(G->dim + 1) * i];
      for (int j = 0; j < G->dim; j++) {
        rank_codes[n_active_ranks].X[j] = rank_buff[(G->dim + 1) * i + j + 1];
      }
      n_active_ranks++;
    }
  }

  free (code_buff);
  free (rank_buff);

  int *send_count = malloc(sizeof(int) * n_ranks);
  int *send_shift = malloc(sizeof(int) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  int *recv_shift = malloc(sizeof(int) * (n_ranks+1));

  int irank = 0;
  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  /* printf("rank codes deb\n"); */
  /* for (int i = 0; i < n_active_ranks; i++) { */
  /*   PDM_morton_dump(3, rank_codes[i]); */
  /* } */
  /* printf("rank codes fin\n"); */

  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_active_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_active_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }
    send_count[active_ranks[irank]] += octant_list->dim + 2;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  PDM_morton_int_t *send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_active_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_active_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }

    int shift = send_shift[active_ranks[irank]] + send_count[active_ranks[irank]];

    assert(octant_list->n_points[i] >= 0);

    send_codes[shift++] = octant_list->codes[i].L;

    for (int j = 0; j < octant_list->dim; j++) {
      send_codes[shift++] = octant_list->codes[i].X[j];
    }

    send_codes[shift++] = (PDM_morton_int_t) octant_list->n_points[i];

    send_count[active_ranks[irank]] += octant_list->dim + 2;
  }

  free (rank_codes);

  PDM_morton_int_t *recv_codes =
    malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */

  PDM_MPI_Alltoallv(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                    recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                    comm);


  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  const int _stride = octant_list->dim + 2;
  const int n_recv_codes = recv_shift[n_ranks] / _stride;

  free (recv_shift);

  int *weight = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i < G->n_nodes; i++) {
    weight[i] = 0;
  }

  /* - compute weight of each cell */

  for (int i = 0; i < n_recv_codes; i++) {

    PDM_morton_code_t code;

    code.L = recv_codes[i*_stride];

    for (int j = 0; j < octant_list->dim; j++) {
      code.X[j] = recv_codes[i*_stride+j+1];
    }

    int G_node =  PDM_morton_binary_search(G->n_nodes,
                                           code,
                                           G->codes);

    weight[G_node] +=  recv_codes[i*_stride + 1 + octant_list->dim];
  }

  free (recv_codes);

  /*
   * Load balancing G from weight
   */

  int *order = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i <  G->n_nodes; i++) {
    order[i] = i;
  }

  *G_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));
  PDM_morton_code_t *_G_morton_index = *G_morton_index;

  PDM_morton_ordered_build_rank_index (octant_list->dim,
                                       max_max_level,
                                       G->n_nodes,
                                       G->codes,
                                       weight,
                                       _G_morton_index,
                                       comm);

  free (order);
  free (weight);

  /* printf("\nblock_partition avant d\n"); */
  /* for (int i = 0; i <  G->n_nodes; i++) { */
  /*   PDM_morton_dump (3, G->codes[i]); */
  /* } */
  /* printf("block_partition avant f\n\n"); */

  _distribute_octants (G, _G_morton_index, comm);

  /*
   * Redistribute octant list from coarse load balancing
   */

  _distribute_octants (octant_list, _G_morton_index, comm);

  free (active_ranks);
  free (n_nodes_rank);
  return G;

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
PDM_para_octree_create
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const PDM_MPI_Comm comm
)
{

  if (_octrees == NULL) {
    _octrees = PDM_Handles_create (4);
  }

  _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

  int id = PDM_Handles_store (_octrees, octree);

  octree->dim = 3;

  for (int i = 0; i < octree->dim; i++) {
    octree->global_extents[i]   = -HUGE_VAL;
    octree->global_extents[octree->dim+i] =  HUGE_VAL;
    octree->s[i]         = 0.;
    octree->d[i]         = 0.;
  }

  octree->depth_max = depth_max;
  octree->points_in_leaf_max = points_in_leaf_max;

  octree->n_point_clouds = n_point_cloud;
  octree->t_n_points = 0;
  octree->n_points = 0;
  octree->points = NULL;
  octree->points_icloud = NULL;
  octree->points_gnum = NULL;
  octree->points_code = NULL;

  octree->octants = NULL;

  octree->n_part_boundary_elt = 0;
  octree->part_boundary_elt_idx = NULL;
  octree->part_boundary_elt = NULL;

  octree->comm = comm;

  octree->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    octree->times_elapsed[i] = 0.;
    octree->times_cpu[i] = 0.;
    octree->times_cpu_u[i] = 0.;
    octree->times_cpu_s[i] = 0.;
  }

  return id;

}


/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);

  if (octree->points != NULL) {
    free (octree->points);
  }

  if (octree->points_icloud != NULL) {
    free (octree->points_icloud);
  }

  if (octree->points_gnum != NULL) {
    free (octree->points_gnum);
  }

  if (octree->points_code != NULL) {
    free (octree->points_code);
  }

  if (octree->part_boundary_elt_idx != NULL) {
    free (octree->part_boundary_elt_idx);
  }

  if (octree->part_boundary_elt != NULL) {
    free (octree->part_boundary_elt);
  }

  if (octree->octants != NULL) {

    if (octree->octants->codes != NULL) {
      free (octree->octants->codes);
    }

    if (octree->octants->range != NULL) {
      free (octree->octants->range);
    }

    if (octree->octants->n_points != NULL) {
      free (octree->octants->n_points);
    }

    if (octree->octants->neighbour_idx != NULL) {
      free (octree->octants->neighbour_idx);
    }

    if (octree->octants->neighbours != NULL) {
      free (octree->octants->neighbours);
    }

    free (octree->octants);
  }

  PDM_timer_free (octree->timer);

  free (octree);

  PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

  const int n_octrees = PDM_Handles_n_get (_octrees);

  if (n_octrees == 0) {
    _octrees = PDM_Handles_free (_octrees);
  }
}


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
PDM_para_octree_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords,
 const PDM_g_num_t *g_num
)
{
  _octree_t *octree = _get_from_id (id);

  const int idx = octree->n_points;

  octree->n_points += n_points;
  octree->points =
    realloc (octree->points, octree->n_points * sizeof(double) * octree->dim);
  octree->points_icloud =
    realloc (octree->points_icloud, octree->n_points * sizeof(int));
  octree->points_gnum =
    realloc (octree->points_gnum, octree->n_points * sizeof(PDM_g_num_t));
  octree->points_code =
    realloc (octree->points_code, octree->n_points * sizeof(PDM_morton_code_t));

  for (int i = 0; i < octree->dim * n_points; i++) {
    octree->points[octree->dim*idx + i] = coords[i];
  }

  for (int i = 0; i < n_points; i++) {
    octree->points_gnum[idx + i] = g_num[i];
  }

  for (int i = 0; i < n_points; i++) {
    octree->points_icloud[idx + i] = i_point_cloud;
  }

}


/**
 *
 * \brief Build octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_build
(
 const int  id
)
{
  _octree_t *octree = _get_from_id (id);

  const int dim = octree->dim;
  //const PDM_morton_int_t max_level = 31u;
  const PDM_morton_int_t max_level = PDM_morton_max_level;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  octree->times_elapsed[BEGIN] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BEGIN]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(octree->timer);

  b_t_elapsed = octree->times_elapsed[BEGIN];
  b_t_cpu     = octree->times_cpu[BEGIN];
  b_t_cpu_u   = octree->times_cpu_u[BEGIN];
  b_t_cpu_s   = octree->times_cpu_s[BEGIN];
  PDM_timer_resume(octree->timer);

  /*
   * Get coord extents
   */

  PDM_morton_get_coord_extents(dim,
                               octree->n_points,
                               octree->points,
                               octree->global_extents,
                               octree->comm);

  /*
   * Encode coords
   */

  PDM_morton_encode_coords(dim,
                           max_level,
                           octree->global_extents,
                           octree->n_points,
                           octree->points,
                           octree->points_code,
                           octree->d,
                           octree->s);

  int *order = malloc (sizeof(int) * octree->n_points);

  for (int i = 0; i < octree->n_points; i++) {
    order[i] = i;
  }

  /**************************************
   *
   * Global order of codes and balancing
   *
   **************************************/

  PDM_morton_local_order (octree->n_points,
                          octree->points_code,
                          order);

  if (n_ranks > 1) {

    int *weight = malloc (sizeof(int) * octree->n_points);
    for (int i = 0; i < octree->n_points; i++) {
      weight[i] = 1;
    }

    PDM_morton_code_t *morton_index =
      malloc (sizeof(PDM_morton_code_t) * (n_ranks + 1));

    PDM_morton_build_rank_index(dim,
                                max_level,
                                octree->n_points,
                                octree->points_code,
                                weight,
                                order,
                                morton_index,
                                octree->comm);

    free (weight);
    free (order);

    /* distribute point from morton_index */

    _distribute_points (&octree->n_points,
                        &octree->points,
                        &octree->points_icloud,
                        &octree->points_gnum,
                        &octree->points_code,
                        morton_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);

    free(morton_index);

  }

  else {

    int *_points_icloud = malloc (sizeof(int) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_icloud[i] =  octree->points_icloud[order[i]];
    }

    free (octree->points_icloud);
    octree->points_icloud = _points_icloud;

    PDM_g_num_t *_points_gnum = malloc (sizeof(PDM_g_num_t) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_gnum[i] =  octree->points_gnum[order[i]];
    }

    free (octree->points_gnum);
    octree->points_gnum = _points_gnum;

    PDM_morton_code_t *_points_code =
      malloc (sizeof(PDM_morton_code_t) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_code[i].L = octree->points_code[order[i]].L;
      _points_code[i].X[0] = octree->points_code[order[i]].X[0];
      _points_code[i].X[1] = octree->points_code[order[i]].X[1];
      _points_code[i].X[2] = octree->points_code[order[i]].X[2];
    }

    free (octree->points_code);
    octree->points_code = _points_code;

    double *_points = malloc (sizeof(double) * dim * octree->n_points);
    for (int i = 0; i < octree->n_points; i++) {
      for (int j = 0; j < dim; j++) {
        _points[dim*i+j] = octree->points[dim*order[i]+j];
      }
    }

    free (octree->points);
    octree->points = _points;

    free (order);
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_ORDER_POINTS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_ORDER_POINTS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_ORDER_POINTS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_ORDER_POINTS]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  PDM_morton_code_t *block_octants_index = NULL;
  if (n_ranks > 1) {

    /*************************************************************************
     *
     * Store points in the octants (leaves) at the maximum depth of the octree
     * to build
     *
     *************************************************************************/

    int chg_code = 1;
    _l_octant_t *point_octants = malloc(sizeof(_l_octant_t));

    int curr_node = -1;

    _octants_init (point_octants, octree->dim, octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {

      PDM_morton_code_t _point_code;
      PDM_morton_copy (octree->points_code[i], &_point_code);

      PDM_morton_assign_level (&_point_code, octree->depth_max);

      if (curr_node != -1) {
        chg_code = !(PDM_morton_a_eq_b (point_octants->codes[curr_node],
                                        _point_code));
      }

      if (chg_code) {

        _octants_check_alloc (point_octants, 1);

        int idx = point_octants->n_nodes;

        curr_node = idx;

        PDM_morton_copy (octree->points_code[i], &(point_octants->codes[idx]));

        point_octants->n_points[idx] = 1;
        point_octants->range[idx] = i;

        point_octants->n_nodes += 1;
      }

      else {
        point_octants->n_points[curr_node] += 1;
      }

    }

    /* printf ("  - n_nodes before block partition : %d\n", point_octants->n_nodes); */
    /* for (int i = 0; i < point_octants->n_nodes; i++) { */
    /*   printf ("  %d : level %u code [%u, %u, %u], range %d , is_leaf %d , n_points %d\n", */
    /*           i, */
    /*           point_octants->codes[i].L, */
    /*           point_octants->codes[i].X[0], */
    /*           point_octants->codes[i].X[1], */
    /*           point_octants->codes[i].X[2], */
    /*           point_octants->range[i], */
    /*           point_octants->is_leaf[i], */
    /*           point_octants->n_points[i] */
    /*           ); */
    /* } */

    /* PDM_morton_code_t *point_octants_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1)); */

    /* PDM_morton_ordered_build_rank_index (dim, */
    /*                                      octree->depth_max, */
    /*                                      point_octants->n_nodes, */
    /*                                      point_octants->codes, */
    /*                                      point_octants->n_points, */
    /*                                      point_octants_morton_index, */
    /*                                      octree->comm); */

    /* _distribute_octants (point_octants, point_octants_morton_index, octree->comm); */

    /* _distribute_points (&octree->n_points, */
    /*                     &octree->points, */
    /*                     &octree->points_icloud, */
    /*                     &octree->points_gnum, */
    /*                     &octree->points_code, */
    /*                     point_octants_morton_index, */
    /*                     octree->comm, */
    /*                     octree->dim, */
    /*                     octree->depth_max, */
    /*                     octree->global_extents); */

    /* free (point_octants_morton_index); */

    /*************************************************************************
     *
     * Block partition (algo 2 sundar)
     *
     *************************************************************************/

    octree->octants = _block_partition (point_octants,
                                        octree->comm,
                                        &block_octants_index);

    _octants_free (point_octants);

    /* printf("\nblock_partition d\n"); */
    /* for (int i = 0; i <  octree->octants->n_nodes; i++) { */
    /*   PDM_morton_dump (3, octree->octants->codes[i]); */
    /* } */
    /* printf("block_partition f\n\n"); */

    /* printf("block_index d\n"); */
    /* for (int i = 0; i <  n_ranks + 1; i++) { */
    /*   PDM_morton_dump (3, block_octants_index[i]); */
    /* } */
    /* printf("block_index f\n\n"); */

    /*************************************************************************
     *
     * Redistribute points
     *
     *************************************************************************/

    _distribute_points (&octree->n_points,
                        &octree->points,
                        &octree->points_icloud,
                        &octree->points_gnum,
                        &octree->points_code,
                        block_octants_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);


    /* printf("Liste des points a ranger d\n"); */
    /* for (int i = 0; i < octree->n_points; i++) { */
    /*   PDM_morton_dump (3, octree->points_code[i]); */
    /* } */
    /* printf("Liste des points a ranger dpoints f\n\n"); */

    int iblock = 0;
    for (int i = 0; i < octree->n_points; i++) {
      /* printf("\n\n---- point d : %d %d %d\n",i, iblock, octree->octants->n_nodes); */
      /* PDM_morton_dump (3, octree->points_code[i]); */
      while (!PDM_morton_ancestor_is (octree->octants->codes[iblock],
                                      octree->points_code[i])) {
        iblock++;
        /* PDM_morton_dump (3, octree->octants->codes[iblock]); */
      }
      /* printf("---- point f : %d %d %d\n",i, iblock, octree->octants->n_nodes); */
      assert (iblock < octree->octants->n_nodes);
      octree->octants->n_points[iblock] += 1;
    }

    octree->octants->range[0] = 0;

    double vol = 0;
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      vol += ((1./pow(2, octree->octants->codes[i].L)) *
              (1./pow(2, octree->octants->codes[i].L)) *
              (1./pow(2, octree->octants->codes[i].L)));
      octree->octants->range[i+1] =
        octree->octants->range[i] +
        octree->octants->n_points[i];
    }
    double total_vol;
    PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, octree->comm);

    if ( (PDM_ABS(total_vol - 1.)>= 1e-15)) {
      printf("Erreur volume different de 1 : %12.5e\n", total_vol);
    }

    assert (PDM_ABS(total_vol - 1.) < 1e-15);

  }

  else {

    octree->octants = malloc(sizeof(_l_octant_t));

    _octants_init (octree->octants, octree->dim, octree->n_points);

    PDM_morton_code_t code;

    code.L = 0;
    code.X[0] = 0;
    code.X[1] = 0;
    code.X[2] = 0;

    _octants_push_back (octree->octants,
                        code,
                        octree->n_points,
                        0);

    octree->octants->range[0] = 0;
    octree->octants->range[1] = octree->n_points;
    octree->octants->n_points[0] = octree->n_points;

  }

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_BLOCK_PARTITION] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_BLOCK_PARTITION]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_BLOCK_PARTITION]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_BLOCK_PARTITION]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Build local octree
   *
   *************************************************************************/

  const int n_child = 8;
  const int n_direction = 6;

  int  size = octree->depth_max * 8;
  _heap_t *heap = _heap_create (size);
  for (int i = octree->octants->n_nodes - 1; i >= 0; i--) {
    int is_pushed = _heap_push (heap,
                                octree->octants->codes[i],
                                octree->octants->range[i],
                                octree->octants->n_points[i]);
    if (!is_pushed) {
      printf ("Internal error PDM_para_octree 3 : heap is full\n");
      exit(1);
    }
  }

  PDM_morton_code_t code;
  int range;
  int n_points;

  octree->octants->n_nodes = 0;

  while (_heap_pull (heap, &code, &range, &n_points)) {

    /* Add children into the heap*/

    if ((code.L < max_morton_level) && (code.L < max_level) &&
        (n_points > octree->points_in_leaf_max)) {

      PDM_morton_code_t children[n_child];
      PDM_morton_get_children(dim,
                              code,
                              children);

      int range_children[n_child];
      int n_points_children[n_child];

      for (int i = 0; i < n_child; i++) {
        n_points_children[i] = 0;
      }

      int ichild = 0;
      for (int i = 0; i < n_points; i++) {
        assert ((range + i) < octree->n_points);
        if (!PDM_morton_ancestor_is(code, octree->points_code[range + i])) {
          printf("Erreur : n'est pas un ancetre !!!!!\n");

        }
        assert (PDM_morton_ancestor_is(code, octree->points_code[range + i]));
        while (!PDM_morton_ancestor_is (children[ichild], octree->points_code[range + i])) {
          ichild += 1;
        }
        assert (ichild < n_child);
        n_points_children[ichild] += 1;
      }

      range_children[0] = 0;
      for (int i = 0; i < n_child - 1; i++) {
        range_children[i+1] = range_children[i] + n_points_children[i];
      }

      for (int i = n_child - 1; i >= 0; i--) {
        int is_pushed = _heap_push (heap,
                                    children[i],
                                    range + range_children[i],
                                    n_points_children[i]);
        if (!is_pushed) {
          printf ("Internal error PDM_para_octree 4 : heap is full\n");
          exit(1);
        }
      }

    }

    /* Store the leaf */

    else {
      _octants_push_back (octree->octants, code, n_points, range);
    }

  }


  double vol = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    vol += ((1./pow(2, octree->octants->codes[i].L)) *
            (1./pow(2, octree->octants->codes[i].L)) *
            (1./pow(2, octree->octants->codes[i].L)));
  }
  double total_vol;
  PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, octree->comm);

  assert (PDM_ABS(total_vol - 1.) < 1e-15);

  heap = _heap_free (heap);

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NODES] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NODES]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NODES]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NODES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /* printf("\noctant nodes d\n"); */
  /* for (int i = 0; i <  octree->octants->n_nodes; i++) { */
  /*   PDM_morton_dump (3, octree->octants->codes[i]); */
  /* } */
  /* printf("\noctant nodes f\n\n"); */

  /*************************************************************************
   *
   * Neighbours
   *
   *************************************************************************/

  //_compute_local_neighbours (octree->octants);

  _neighbours_tmp_t *neighbours_tmp = malloc (sizeof(_neighbours_tmp_t) * octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j =  0; j < n_direction; j++) {
      neighbours_tmp[i].n_neighbour[j] = 0;
      neighbours_tmp[i].s_neighbour[j] = 1;
      neighbours_tmp[i].neighbours[j] = malloc (sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
      neighbours_tmp[i].neighbours[j][0] = 0;
    }
  }

  /* Boucle sur les noeuds : */

  int *intersected_quantile = malloc(sizeof(int) * n_ranks);
  size_t  n_intersected_quantile = 0;

  int *intersect_nodes = malloc (sizeof(int) * octree->octants->n_nodes);
  size_t n_intersect_nodes = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 1; j < n_direction; j+=2) {

      PDM_morton_code_t *neighbour_code =
        _neighbour (octree->octants->codes[i], (PDM_para_octree_direction_t)j);
      PDM_para_octree_direction_t inv_j =
        _inv_direction((PDM_para_octree_direction_t) j);

      if (neighbour_code != NULL) {

        if (block_octants_index != NULL) {
          PDM_morton_quantile_intersect (n_ranks,
                                         *neighbour_code,
                                         block_octants_index,
                                         &n_intersected_quantile,
                                         intersected_quantile);
        }
        else {
          n_intersected_quantile = 1;
          intersected_quantile[0] = 0;
        }

        for (int i_inter = 0; i_inter <  n_intersected_quantile; i_inter++) {
          int neighbour_rank = intersected_quantile[i_inter];

          if (neighbour_rank == rank) {

            PDM_morton_quantile_intersect(octree->octants->n_nodes - (i+1),
                                          *neighbour_code,
                                          octree->octants->codes + i + 1,
                                          &n_intersect_nodes,
                                          intersect_nodes);

            for (int k = 0; k < n_intersect_nodes; k++) {
              int idx = intersect_nodes[k];
              if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
                neighbours_tmp[i].s_neighbour[j] *= 2;
                neighbours_tmp[i].neighbours[j] =
                  realloc (neighbours_tmp[i].neighbours[j],
                           sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
              }
              neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = idx;

              if (neighbours_tmp[idx].n_neighbour[inv_j] >= neighbours_tmp[idx].s_neighbour[inv_j]) {
                neighbours_tmp[idx].s_neighbour[inv_j] *= 2;
                neighbours_tmp[idx].neighbours[inv_j] =
                  realloc (neighbours_tmp[idx].neighbours[inv_j],
                           sizeof(int) * neighbours_tmp[idx].s_neighbour[inv_j]);
              }
              neighbours_tmp[idx].neighbours[inv_j][neighbours_tmp[idx].n_neighbour[inv_j]++] = i;
            }

          }

          else {
            if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
              neighbours_tmp[i].s_neighbour[j] *= 2;
              neighbours_tmp[i].neighbours[j] = realloc (neighbours_tmp[i].neighbours[j],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
            }
            neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = - (neighbour_rank + 1);
          }
        }
        free (neighbour_code);
      }
    }
  }

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j+=2) {
      PDM_morton_code_t *neighbour_code =
        _neighbour (octree->octants->codes[i], (PDM_para_octree_direction_t) j);

      if (neighbour_code != NULL) {

        if (block_octants_index != NULL) {

          PDM_morton_quantile_intersect (n_ranks,
                                         *neighbour_code,
                                         block_octants_index,
                                         &n_intersected_quantile,
                                         intersected_quantile);

        }
        else {
          n_intersected_quantile = 1;
          intersected_quantile[0] = 0;
        }

        for (int i_inter = 0; i_inter <  n_intersected_quantile; i_inter++) {
          int neighbour_rank = intersected_quantile[i_inter];

          if (neighbour_rank != rank) {

            if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
              neighbours_tmp[i].s_neighbour[j] *= 2;
              neighbours_tmp[i].neighbours[j] = realloc (neighbours_tmp[i].neighbours[j],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
            }
            neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = - (neighbour_rank + 1);
          }
        }
        free (neighbour_code);
      }
    }
  }

  free (intersected_quantile);

  if (block_octants_index != NULL) {
    free(block_octants_index);
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Copy temporary neighbours in the neighbour structure
   *
   *************************************************************************/

  octree->octants->neighbour_idx =
    malloc(sizeof(int) * (n_direction * octree->octants->n_nodes + 1));

  int idx = 0;
  octree->octants->neighbour_idx[0] = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      octree->octants->neighbour_idx[idx+1] =
        octree->octants->neighbour_idx[idx] + neighbours_tmp[i].n_neighbour[j];
      idx += 1;
    }
  }

  octree->octants->neighbours =
    malloc(sizeof(int) *
           octree->octants->neighbour_idx[n_direction * octree->octants->n_nodes]);

  idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
        octree->octants->neighbours[idx++] = neighbours_tmp[i].neighbours[j][k];
      }
    }
  }

  /* Free temporary arrays */
  /* printf("sortie 2 neighbours_tmp debut\n"); */
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      if (neighbours_tmp[i].neighbours[j] != NULL) {
        /* if (neighbours_tmp[i].neighbours[j][0] < 0) { */
          /* printf ("i neighbours[j] %d [%d / %d %d %d] %d : ", */
          /*         i, */
          /*         octree->octants->codes[i].L, */
          /*         octree->octants->codes[i].X[0], */
          /*         octree->octants->codes[i].X[1], */
          /*         octree->octants->codes[i].X[2], */
          /*         j); */
          /* for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) { */
          /*   printf (" %d",neighbours_tmp[i].neighbours[j][k]); */
          /* } */
          /* printf ("\n"); */
        /* } */
        free (neighbours_tmp[i].neighbours[j]);
      }
    }
  }
  /* printf("sortie 2 neighbours_tmp fin\n"); */

  free (neighbours_tmp);

  /* for (int i = 0; i < octree->octants->n_nodes; i++) { */
  /*   for (int j = 0; j < n_direction; j++) { */
  /*     idx = octree->octants->neighbour_idx[i*n_direction+j]; */
  /*     int n_b = octree->octants->neighbour_idx[i*n_direction+j+1] */
  /*             - octree->octants->neighbour_idx[i*n_direction+j]; */
  /*     if (n_b > 0) { */
  /*       printf ("i neighbours[j] %d [%d / %d %d %d] %d : ", */
  /*               i, */
  /*               octree->octants->codes[i].L, */
  /*               octree->octants->codes[i].X[0], */
  /*               octree->octants->codes[i].X[1], */
  /*               octree->octants->codes[i].X[2], */
  /*               j); */
  /*       for (int k = octree->octants->neighbour_idx[i*n_direction+j]; */
  /*            k < octree->octants->neighbour_idx[i*n_direction+j+1]; */
  /*            k++) { */

  /*         if (octree->octants->neighbours[k] < 0) { */
  /*           printf (" %d", octree->octants->neighbours[k]); */
  /*         } */
  /*       } */
  /*       printf ("\n"); */
  /*     } */
  /*   } */
  /* } */


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_s - b_t_cpu_s;


  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS] +=  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS] +=  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS] +=  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS] +=  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3];

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /* PDM_timer_hang_on(octree->timer); */
  /* e_t_elapsed = PDM_timer_elapsed(octree->timer); */
  /* e_t_cpu     = PDM_timer_cpu(octree->timer); */
  /* e_t_cpu_u   = PDM_timer_cpu_user(octree->timer); */
  /* e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer); */

  /* b_t_elapsed = e_t_elapsed; */
  /* b_t_cpu     = e_t_cpu; */
  /* b_t_cpu_u   = e_t_cpu_u; */
  /* b_t_cpu_s   = e_t_cpu_s; */

  /* PDM_timer_resume(octree->timer); */

  /*************************************************************************
   *
   * Build parallel partition boundary
   *
   *************************************************************************/

  if (n_ranks > 1) {

    const int n_quantile =  n_ranks * n_direction;

    int *neighbour_rank_n = malloc (sizeof(int) * n_quantile);
    int *neighbour_rank_idx = malloc (sizeof(int) * (n_quantile + 1));

    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_n[i] = 0;
    }

    /* Premiere boucle pour compter */

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = octree->octants->neighbour_idx[n_direction * i + j];
                 k < octree->octants->neighbour_idx[n_direction * i + j + 1];
                 k++) {
          if (octree->octants->neighbours[k] < 0) {
            neighbour_rank_n[(PDM_ABS(octree->octants->neighbours[k]) - 1)*n_direction +j]++;
          }
        }
      }
    }

    int max_node_dir = -1;
    neighbour_rank_idx[0] = 0;
    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_idx[i+1] = neighbour_rank_idx[i] + neighbour_rank_n[i];
      max_node_dir = PDM_MAX (max_node_dir, neighbour_rank_n[i]);
      neighbour_rank_n[i] = 0;
    }

    /* Allocation */

    /* printf("max_node_dir : %d\n", max_node_dir); */

    int *neighbour_rank_node_id = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *neighbour_rank_code = malloc (sizeof(PDM_morton_code_t) *
                                                     neighbour_rank_idx[n_quantile]);

    /* Deuxieme boucle pour stocker avec tri suivant la direction */

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = octree->octants->neighbour_idx[n_direction * i + j];
             k < octree->octants->neighbour_idx[n_direction * i + j + 1];
             k++) {
          if (octree->octants->neighbours[k] < 0) {
            int index = (PDM_ABS(octree->octants->neighbours[k]) - 1) * n_direction + j;
            int index2 = neighbour_rank_idx[index] + neighbour_rank_n[index];

            neighbour_rank_node_id[index2] = i;
            PDM_morton_copy (octree->octants->codes[i],
                             neighbour_rank_code + index2);

            neighbour_rank_n[index]++;
          }
        }
      }
    }

    /* printf("Avant tri d\n"); */
    /* for (int i = 0; i < n_ranks; i++) { */
    /*   for (int j = 0; j < n_direction; j++) { */
    /*     for (int k = neighbour_rank_idx[n_direction * i + j]; */
    /*          k < neighbour_rank_idx[n_direction * i + j + 1]; */
    /*          k++) { */
    /*       printf("rank dir node_id : %d %d %d\n", i, j, neighbour_rank_node_id[k]); */
    /*     } */
    /*   } */
    /* } */
    /* printf("Avant tri f\n"); */

    /* Tri des codes pour chaque direction de chaque rang */

    order = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_id = malloc (sizeof(int) * max_node_dir);
    PDM_morton_code_t *tmp_code = malloc (sizeof(PDM_morton_code_t) * max_node_dir);

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        PDM_morton_local_order (neighbour_rank_n[n_direction * i + j],
                                neighbour_rank_code + neighbour_rank_idx[n_direction * i + j],
                                order);
        int idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (neighbour_rank_code[k], tmp_code + idx1);
          tmp_node_id[idx1] = neighbour_rank_node_id[k];
          idx1 += 1;
        }

        idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (tmp_code[order[idx1]], neighbour_rank_code + k );
          neighbour_rank_node_id[k] = tmp_node_id[order[idx1]];
          idx1 += 1;
        }
      }
    }


    /* printf("Apres tri d\n"); */
    /* for (int i = 0; i < n_ranks; i++) { */
    /*   for (int j = 0; j < n_direction; j++) { */
    /*     for (int k = neighbour_rank_idx[n_direction * i + j]; */
    /*          k < neighbour_rank_idx[n_direction * i + j + 1]; */
    /*          k++) { */
    /*       printf("rank dir node_id : %d %d %d\n", i, j, neighbour_rank_node_id[k]); */
    /*     } */
    /*   } */
    /* } */
    /* printf("Apres tri f\n"); */

    free (tmp_code);
    free (order);
    free (tmp_node_id);

    /* Envoi / reception (Les donnees recues sont triees) */

    int *recv_neighbour_rank_n = malloc (sizeof(int) * n_quantile);

    for (int i = 0; i < n_quantile; i++) {
      recv_neighbour_rank_n[i] = 0;
    }

    PDM_MPI_Request *recv_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);
    PDM_MPI_Request *send_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);

    int *used_ranks = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Alltoall (neighbour_rank_n, n_direction, PDM_MPI_INT,
                      recv_neighbour_rank_n, n_direction, PDM_MPI_INT,
                      octree->comm);

    int *recv_neighbour_rank_idx = malloc (sizeof(int) * (n_direction * n_ranks + 1));
    recv_neighbour_rank_idx[0] = 0;

    /* printf("recv_neighbour_rank_idx : %d",recv_neighbour_rank_idx[0]); */
    for (int i = 0; i <  n_direction * n_ranks; i++) {
      recv_neighbour_rank_idx[i+1] = recv_neighbour_rank_idx[i] + recv_neighbour_rank_n[i];
      /* printf(" %d",recv_neighbour_rank_idx[i]); */
    }
    /* printf("\n"); */

    int *recv_neighbour_rank_node_id = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *recv_neighbour_rank_code =
      malloc (sizeof(PDM_morton_code_t) * recv_neighbour_rank_idx[n_quantile]);

    unsigned int *_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * neighbour_rank_idx[n_quantile]);
    unsigned int *_recv_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * recv_neighbour_rank_idx[n_quantile]);

    idx = 0;
    for (int i = 0; i < neighbour_rank_idx[n_quantile]; i++) {
      _neighbour_rank_code[idx++] = neighbour_rank_code[i].L;
      for (int j = 0; j < 3; j++) {
        _neighbour_rank_code[idx++] = neighbour_rank_code[i].X[j];
      }
    }

    int *rank_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));
    int *rank_recv_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_recv_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));

    rank_neighbour_rank_idx[0] = 0;
    rank_recv_neighbour_rank_idx[0] = 0;

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] = 0;
      rank_recv_neighbour_rank_n[i] = 0;
    }

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        rank_neighbour_rank_n[i] += neighbour_rank_n[i*n_direction+j];
        rank_recv_neighbour_rank_n[i] += recv_neighbour_rank_n[i*n_direction+j];
      }
      rank_neighbour_rank_idx[i+1] = rank_neighbour_rank_n[i] + rank_neighbour_rank_idx[i];
      rank_recv_neighbour_rank_idx[i+1] = rank_recv_neighbour_rank_n[i] + rank_recv_neighbour_rank_idx[i];
    }

    PDM_MPI_Alltoallv (neighbour_rank_node_id,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_id,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] *= 4;
      rank_recv_neighbour_rank_n[i] *= 4;
      rank_neighbour_rank_idx[i+1] *= 4;
      rank_recv_neighbour_rank_idx[i+1] *= 4;
    }


    PDM_MPI_Alltoallv (_neighbour_rank_code,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       _recv_neighbour_rank_code,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       octree->comm);


    free (_neighbour_rank_code);

    free (rank_neighbour_rank_n);
    free (rank_neighbour_rank_idx);
    free (rank_recv_neighbour_rank_n);
    free (rank_recv_neighbour_rank_idx);

    idx = 0;
    for (int i = 0; i < recv_neighbour_rank_idx[n_quantile]; i++) {
      recv_neighbour_rank_code[i].L = _recv_neighbour_rank_code[idx++];
      for (int j = 0; j < 3; j++) {
        recv_neighbour_rank_code[i].X[j] = _recv_neighbour_rank_code[idx++];
      }
    }
    free (_recv_neighbour_rank_code);

    free (recv_request);
    free (send_request);

    free (used_ranks);

    octree->n_part_boundary_elt = neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt_idx = malloc (sizeof(int) * (octree->n_part_boundary_elt + 1));

    int s_part_boundary_elt = 2 * 2 * neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt = malloc (sizeof(int) * s_part_boundary_elt);

    int n_part_boundary_elt = 0;

    idx = 0;
    int idx_part_boundary_elt = 0;
    octree->part_boundary_elt_idx[0] = 0;

    /* printf("\nnodes d\n"); */
    /* for (int k1 = 0; k1 < octree->octants->n_nodes; k1++) { */
    /*   PDM_morton_dump(3,  octree->octants->codes[k1]); */
    /* } */
    /* printf("nodes f\n"); */

    int *intersect = malloc(sizeof(int) *  recv_neighbour_rank_idx[n_quantile]);
    size_t n_intersect;

    for (int i = 0; i < n_ranks; i++ ) {

      for (PDM_para_octree_direction_t j = PDM_BOTTOM; j < n_direction; j++) {

        /* printf ("irank, idrection, n : %d %d %d \n", */
        /*         i, j, */
        /*         neighbour_rank_idx[i * n_direction+j+1] - neighbour_rank_idx[i * n_direction +j]); */

        int idx_recv = n_direction * i + _inv_direction(j);

        int idx_candidate = recv_neighbour_rank_idx[idx_recv];
        int n_candidate = recv_neighbour_rank_idx[idx_recv+1] - idx_candidate;

        if (n_candidate > 0) {

          for (int k = neighbour_rank_idx[i * n_direction + j];
                   k < neighbour_rank_idx[i * n_direction + j +1]; k++) {

            PDM_morton_code_t *neighbour_code = _neighbour (neighbour_rank_code[k], j);

            /* printf("\n*********************\nlocal code d\n"); */
            /* PDM_morton_dump(3, neighbour_rank_code[k]); */
            /* printf("local code f\n"); */

            /* printf("neighbour code d\n"); */
            /* PDM_morton_dump(3, *neighbour_code); */
            /* printf("neighbour code f\n"); */
            /* assert (neighbour_code != NULL); */

            /* printf("\ncandidates d\n"); */
            /* for (int k1 = 0; k1 < n_candidate; k1++) { */
            /*   PDM_morton_dump(3, (recv_neighbour_rank_code + idx_candidate)[k1]); */
            /* } */
            /* printf("candidates f\n"); */

            /* printf("\nintersect d\n"); */
            /* for (int k1 = 0; k1 < n_intersect; k1++) { */
            /*   PDM_morton_dump(3, (recv_neighbour_rank_code + idx_candidate)[ intersect[k1]]); */
            /* } */
            /* printf("intersect f\n"); */

            PDM_morton_quantile_intersect(n_candidate,
                                          *neighbour_code,
                                          recv_neighbour_rank_code + idx_candidate,
                                          &n_intersect,
                                          intersect);

            if (n_intersect > 0) {
              for (int k1 = 0; k1 < n_intersect; k1++) {
                if ((s_part_boundary_elt - n_part_boundary_elt) <= 2) {
                  s_part_boundary_elt *= 2;
                  octree->part_boundary_elt = realloc (octree->part_boundary_elt,
                                                       sizeof(int) * s_part_boundary_elt);
                }
                octree->part_boundary_elt_idx[n_part_boundary_elt+1]++;
                octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_id[idx_candidate + intersect[k1]];
              }
              octree->octants->neighbours[neighbour_rank_node_id[k]] = -(n_part_boundary_elt+1);
              n_part_boundary_elt++;
            }

            free (neighbour_code);

          }
        }
      }
    }

    free (intersect);

    if (octree->n_part_boundary_elt < n_part_boundary_elt) {

      printf("octree->n_part_boundary_elt == n_part_boundary_elt : %d %d\n", octree->n_part_boundary_elt, n_part_boundary_elt);
    }
    assert (octree->n_part_boundary_elt >= n_part_boundary_elt);
    octree->n_part_boundary_elt = n_part_boundary_elt;
    octree->part_boundary_elt_idx = realloc (octree->part_boundary_elt_idx, sizeof(int) * (n_part_boundary_elt+1));


    free (neighbour_rank_n);
    free (neighbour_rank_idx);
    free (neighbour_rank_node_id);
    free (neighbour_rank_code);

    free (recv_neighbour_rank_n);
    free (recv_neighbour_rank_idx);

    free (recv_neighbour_rank_node_id);
    free (recv_neighbour_rank_code);

  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_DISTANT_NEIGHBOURS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_DISTANT_NEIGHBOURS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  PDM_timer_hang_on(octree->timer);
  octree->times_elapsed[BUILD_TOTAL] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BUILD_TOTAL]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BUILD_TOTAL]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BUILD_TOTAL]   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  octree->times_elapsed[END] = octree->times_elapsed[BUILD_TOTAL];
  octree->times_cpu[END]     = octree->times_cpu[BUILD_TOTAL];
  octree->times_cpu_u[END]   = octree->times_cpu_u[BUILD_TOTAL];
  octree->times_cpu_s[END]   = octree->times_cpu_s[BUILD_TOTAL];

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
PDM_para_octree_extents_get
(
 const int  id
)
{
 _octree_t *octree = _get_from_id (id);

 return octree->global_extents;
}


/**
 *
 * \brief Dump octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_dump
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);
  PDM_printf("PDM_dump_para_octree : %d\n",id);

  PDM_printf("  - n_nodes : %d\n", octree->octants->n_nodes);
  printf("  - global_extents :");
  for (int i = 0; i < 2*octree->dim; i++) {
    printf(" %12.5e", octree->global_extents[i]);
  }
  PDM_printf("\n");
  PDM_printf("  - depth_max : %d\n", octree->depth_max);
  PDM_printf("  - points_in_leaf_max : %d\n", octree->points_in_leaf_max);

  PDM_printf("  - s : %12.5e %12.5e %12.5e\n", octree->s[0], octree->s[1], octree->s[2]);
  PDM_printf("  - d : %12.5e %12.5e %12.5e\n", octree->d[0], octree->d[1], octree->d[2]);

  PDM_printf("  - n_point_clouds : %d\n", octree->n_point_clouds);
  PDM_printf("  - t_n_points : "PDM_FMT_G_NUM"\n", octree->t_n_points);
  PDM_printf("  - n_points : %d\n", octree->n_points);
  for (int i = 0; i < octree->n_points; i++) {
    PDM_printf("  %d "PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e / level %u code [%u, %u, %u]\n",
               i, octree->points_gnum[i],
               octree->points[3*i], octree->points[3*i+1], octree->points[3*i+2],
               octree->points_code[i].L,
               octree->points_code[i].X[0],
               octree->points_code[i].X[1],
               octree->points_code[i].X[2]);
  }
  PDM_printf("  - n_nodes : %d\n", octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    PDM_printf("  %d : level %u code [%u, %u, %u], range %d, n_points %d\n",
               i,
               octree->octants->codes[i].L,
               octree->octants->codes[i].X[0],
               octree->octants->codes[i].X[1],
               octree->octants->codes[i].X[2],
               octree->octants->range[i],
               octree->octants->n_points[i]
               );
  }
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    PDM_printf("  %d : neighbors\n", i);
    for (int j = 0; j < 6; j++) {
      PDM_printf("    - direction %d : ", j);
      for (int k = octree->octants->neighbour_idx[6*i+j];
         k < octree->octants->neighbour_idx[6*i+j+1]; k++) {
        PDM_printf(" %d",octree->octants->neighbours[k]);
      }
      PDM_printf("\n");
    }
  }
}


/**
 *
 * Look for closest points stored inside an octree
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_closest_points       Number of closest points to find
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest points in octree global number
 * \param [out]  closest_octree_pt_dist Closest points in octree distance
 *
 */

void
PDM_para_octree_closest_point
(
const int    id,
const int    n_closest_points,
const int    n_pts,
double      *pts,
PDM_g_num_t *pts_g_num,
PDM_g_num_t *closest_octree_pt_g_num,
double      *closest_octree_pt_dist2
)
{
  // _octree_t *octree = _get_from_id (id);
}

/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_para_octree_dump_times
(
 const int id
)
{
  _octree_t *octree = _get_from_id (id);

  double t1 = octree->times_elapsed[END] - octree->times_elapsed[BEGIN];
  double t2 = octree->times_cpu[END] - octree->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (octree->times_elapsed, t_elaps_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (octree->times_cpu, t_cpu_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  if (rank == 0) {

    PDM_printf( "PDM_para_octree timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "PDM_para_octree timer : build octree : total (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_TOTAL],
                t_cpu_max[BUILD_TOTAL]);
    PDM_printf( "PDM_para_octree timer : build octree : step order points (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_ORDER_POINTS],
                t_cpu_max[BUILD_ORDER_POINTS]);
    PDM_printf( "PDM_para_octree timer : build octree : step block partition (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BLOCK_PARTITION],
                t_cpu_max[BUILD_BLOCK_PARTITION]);
    PDM_printf( "PDM_para_octree timer : build octree : step local nodes (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NODES],
                t_cpu_max[BUILD_LOCAL_NODES]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 1 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP1],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP1]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 2 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP2],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP2]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 3 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP3],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP3]);
    PDM_printf( "PDM_para_octree timer : build octree : step distant neighbours (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_DISTANT_NEIGHBOURS],
                t_cpu_max[BUILD_DISTANT_NEIGHBOURS]);

  }

}




#ifdef __cplusplus
}
#endif /* __cplusplus */
