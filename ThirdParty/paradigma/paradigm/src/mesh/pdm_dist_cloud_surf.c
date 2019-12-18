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
#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
#include "pdm_octree.h"
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"
#include "pdm_hash_tab.h"


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

#define NTIMER 8

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  UPPER_BOUND_DIST              = 1,
  CANDIDATE_SELECTION           = 2,
  LOAD_BALANCING_ELEM_DIST      = 3,
  COMPUTE_ELEM_DIST             = 4,
  RESULT_TRANSMISSION           = 5,
  END                           = 6,
  BBTREE_CREATE                 = 7,

} _ol_timer_step_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;
  double      **dist;
  double      **proj;
  PDM_g_num_t **closest_elt_gnum;

} _points_cloud_t;

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int  n_point_cloud; /*!< Number of point clouds */
  PDM_MPI_Comm comm;  /*!< MPI communicator */

  PDM_mesh_nature_t mesh_nature;  /*!< Nature of the mesh */

  PDM_surf_mesh_t *surf_mesh;  /*!< Surface mesh pointer */
  PDM_surf_mesh_t *_surf_mesh;  /*!< Surface mesh pointer */

  int  mesh_nodal_id;  /*!< Surface mesh identifier */

  _points_cloud_t *points_cloud; /*!< Point clouds */

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
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 */

int
PDM_dist_cloud_surf_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
)
{
  if (_dists == NULL) {
    _dists = PDM_Handles_create (4);
  }

  _PDM_dist_t *dist = (_PDM_dist_t *) malloc(sizeof(_PDM_dist_t));

  int id = PDM_Handles_store (_dists, dist);

  dist->mesh_nature = mesh_nature;
  dist->mesh_nodal_id = -1;
  dist->surf_mesh = NULL;
  dist->_surf_mesh = NULL;
  dist->n_point_cloud = n_point_cloud;
  dist->comm = comm;
  dist->points_cloud =
    (_points_cloud_t*) malloc (sizeof(_points_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    dist->points_cloud[i].n_part = -1;
    dist->points_cloud[i].n_points = NULL;
    dist->points_cloud[i].coords = NULL;
    dist->points_cloud[i].gnum = NULL;
    dist->points_cloud[i].dist = NULL;
    dist->points_cloud[i].proj = NULL;
    dist->points_cloud[i].closest_elt_gnum = NULL;
  }

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
PDM_dist_cloud_surf_create_cf
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Fint comm,
 int *id
)

{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_dist_cloud_surf_create (mesh_nature, n_point_cloud, _comm);

}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_dist_cloud_surf_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  dist->points_cloud[i_point_cloud].n_part = n_part;
  dist->points_cloud[i_point_cloud].n_points =
    realloc(dist->points_cloud[i_point_cloud].n_points, n_part * sizeof(int));
  dist->points_cloud[i_point_cloud].coords =
    realloc(dist->points_cloud[i_point_cloud].coords,
            n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].gnum =
    realloc(dist->points_cloud[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));
  dist->points_cloud[i_point_cloud].dist =
    realloc(dist->points_cloud[i_point_cloud].dist, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].proj =
    realloc(dist->points_cloud[i_point_cloud].proj, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].closest_elt_gnum =
    realloc(dist->points_cloud[i_point_cloud].closest_elt_gnum,
            n_part * sizeof(PDM_g_num_t * ));

  for (int i = 0; i < n_part; i++) {
    dist->points_cloud[i_point_cloud].n_points[i] = -1;
    dist->points_cloud[i_point_cloud].coords[i] = NULL;
    dist->points_cloud[i_point_cloud].gnum[i] = NULL;
    dist->points_cloud[i_point_cloud].dist[i] = NULL;
    dist->points_cloud[i_point_cloud].proj[i] = NULL;
    dist->points_cloud[i_point_cloud].closest_elt_gnum[i] = NULL;
  }
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_dist_cloud_surf_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  dist->points_cloud[i_point_cloud].n_points[i_part] = n_points;
  dist->points_cloud[i_point_cloud].coords[i_part] = coords;
  dist->points_cloud[i_point_cloud].gnum[i_part] = gnum;
}


/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Identifier
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 *
 */

void
PDM_dist_cloud_surf_nodal_mesh_set
(
 const int  id,
 const int  mesh_nodal_id
)
{
  _PDM_dist_t *dist = _get_from_id (id);
  dist->mesh_nodal_id = mesh_nodal_id;
}




/**
 *
 * \brief Map a surface mesh
 *
 * \param [in]   id         Identifier
 * \param [in]   surf_mesh  Surface mesh pointer
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_map
(
 const int  id,
 PDM_surf_mesh_t *surf_mesh
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  dist->_surf_mesh = surf_mesh;
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
PDM_dist_cloud_surf_surf_mesh_global_data_set
(
 const int         id,
 const PDM_g_num_t n_g_face,
 const PDM_g_num_t n_g_vtx,
 const int         n_part
)
{

  _PDM_dist_t *dist = _get_from_id (id);

  assert (dist->surf_mesh == NULL);

  dist->surf_mesh =
    PDM_surf_mesh_create (n_g_face, n_g_vtx, n_part, dist->comm);
  dist->_surf_mesh = dist->surf_mesh;
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
 * \param [in]   vtx_ln_to_gn  Local vertex numbering
 *                             to global vertex numbering
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_part_set
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
 * \brief Compute distance
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_dist_cloud_surf_compute
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  const int n_point_cloud = dist->n_point_cloud;
  const int mesh_id = dist->mesh_nodal_id;
  PDM_MPI_Comm comm = dist->comm;

  int rank;
  PDM_MPI_Comm_rank (comm, &rank);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  //PDM_timer_hang_on(dist->timer);
  dist->times_elapsed[BEGIN] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[BEGIN]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  /*
   * For each cloud
   */

  for (int i_point_cloud = 0; i_point_cloud < n_point_cloud; i_point_cloud++) {

    _points_cloud_t *pt_cloud = &(dist->points_cloud[i_point_cloud]);
    const int n_part = pt_cloud->n_part;

    /***************************************************************************
     *
     * Compute the upper bound distance. It is the distance from the closest
     * vertex
     *      - Store mesh vetices in a octree
     *      - Compute the closest distance of points to vertices stored in
     *        the octree
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    const double tolerance = 1e-4;
    // const int depth_max = 35;
    const int depth_max = 31;
    const int points_in_leaf_max = 4;

    int n_part_mesh = 0;
    if (dist->mesh_nodal_id != -1) {
      n_part_mesh = PDM_Mesh_nodal_n_part_get (mesh_id);
    }
    else if (dist->_surf_mesh != NULL) {
      n_part_mesh = PDM_surf_mesh_n_part_get (dist->_surf_mesh);
      /* PDM_surf_mesh_compute_faceExtentsMesh (dist->_surf_mesh, 1e-8); */

      /* double glob_extents[6]; */
      /* glob_extents[0] = HUGE_VAL; */
      /* glob_extents[1] = HUGE_VAL; */
      /* glob_extents[2] = HUGE_VAL; */
      /* glob_extents[3] = -HUGE_VAL; */
      /* glob_extents[4] = -HUGE_VAL; */
      /* glob_extents[5] = -HUGE_VAL; */

      /* for (int i = 0; i < n_part_mesh; i++) { */
      /*   const double *_extents = PDM_surf_mesh_part_extents_get (dist->_surf_mesh, i); */
      /*   glob_extents[0] = PDM_MIN (_extents[0], glob_extents[0]); */
      /*   glob_extents[1] = PDM_MIN (_extents[1], glob_extents[1]); */
      /*   glob_extents[2] = PDM_MIN (_extents[2], glob_extents[2]); */
      /*   glob_extents[3] = PDM_MAX (_extents[3], glob_extents[3]); */
      /*   glob_extents[4] = PDM_MAX (_extents[4], glob_extents[4]); */
      /*   glob_extents[5] = PDM_MAX (_extents[5], glob_extents[5]); */
      /* } */
      /* double min_dim = HUGE_VAL; */
      /* min_dim = glob_extents[3] - glob_extents[0]; */
      /* min_dim = PDM_MIN (min_dim, glob_extents[4] - glob_extents[1]); */
      /* min_dim = PDM_MIN (min_dim, glob_extents[5] - glob_extents[2]); */

      /* depth_max = 0; */
      /* while ((min_dim) > 1e-10) { */
      /*   min_dim = min_dim / 2.; */
      /*   depth_max += 1; */
      /* } */
    }
    else {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_dist_cloud_surf error : The surface mesh is not defined. "
                "To do that : \n"
                "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
                "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
                " PDM_dist_cloud_surf_surf_mesh_part_set\n");
    }

    int octree_id = PDM_octree_create (n_part_mesh,
                                       depth_max,
                                       points_in_leaf_max,
                                       tolerance,
                                       comm);

    for (int i_part = 0; i_part < n_part_mesh; i_part++) {

      int n_vertices = 0;
      const double *vertices_coords = NULL;
      const PDM_g_num_t *vertices_gnum = NULL;

      if (dist->mesh_nodal_id != -1) {
        n_vertices      = PDM_Mesh_nodal_n_vertices_get (mesh_id, i_part);
        vertices_coords = PDM_Mesh_nodal_vertices_get (mesh_id, i_part);
        vertices_gnum   = PDM_Mesh_nodal_vertices_g_num_get (mesh_id, i_part);
      }
      else if (dist->_surf_mesh != NULL) {
        n_vertices      = PDM_surf_mesh_part_n_vtx_get(dist->_surf_mesh, i_part);
        vertices_coords = PDM_surf_mesh_part_vtx_get (dist->_surf_mesh, i_part);
        vertices_gnum   = PDM_surf_mesh_part_vtx_g_num_get (dist->_surf_mesh,
                                                            i_part);
      }
      else {
        PDM_error(__FILE__, __LINE__, 0,
                  "PDM_dist_cloud_surf error : The surface mesh is not defined. "
                  "To do that : \n"
                  "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
                  "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
                  " PDM_dist_cloud_surf_surf_mesh_part_set\n");
      }

      PDM_octree_point_cloud_set (octree_id, i_part, n_vertices,
                                  vertices_coords, vertices_gnum);

    }

    /*
     * Build octree
     */

    PDM_octree_build (octree_id);

    /*
     * Concatenation of the partitions
     */

    int n_pts_rank = 0;

    for (int i_part = 0; i_part < n_part; i_part++) {
      n_pts_rank += pt_cloud->n_points[i_part];
    }

    double *pts_rank = malloc (sizeof(double) * n_pts_rank * 3);
    PDM_g_num_t *pts_g_num_rank = malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    n_pts_rank = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < pt_cloud->n_points[i_part]; i++) {
        for (int k = 0; k < 3; k++) {
          pts_rank[3*(n_pts_rank + i) + k] = pt_cloud->coords[i_part][3*i+k];
          pts_g_num_rank[n_pts_rank + i] = pt_cloud->gnum[i_part][i];
        }
      }
      n_pts_rank += pt_cloud->n_points[i_part];
    }

    /*
     * Look for closest surface mesh vertices
     */

    PDM_g_num_t *closest_vertices_gnum =
      malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    double *closest_vertices_dist2 =  malloc (sizeof(double) * n_pts_rank);

    PDM_octree_closest_point (octree_id,
                              n_pts_rank,
                              pts_rank,
                              pts_g_num_rank,
                              closest_vertices_gnum,
                              closest_vertices_dist2);

    //      debut test cube :

    /* int ierr = 0; */
    /* double xmin = 0.; */
    /* double ymin = 0.; */
    /* double zmin = 0.; */
    /* double xmax = 1.; */
    /* double ymax = 1.; */
    /* double zmax = 1.; */
    /* for (int i = 0; i < n_pts_rank; i++) { */
    /*   double d1 = PDM_MIN (PDM_ABS (pts_rank[3*i] - xmin), PDM_ABS (pts_rank[3*i] - xmax)); */
    /*   double d2 = PDM_MIN (PDM_ABS (pts_rank[3*i+1] - ymin), PDM_ABS (pts_rank[3*i+1] - ymax)); */
    /*   double d3 = PDM_MIN (PDM_ABS (pts_rank[3*i+2] - zmin), PDM_ABS (pts_rank[3*i+2] - zmax)); */
    /*   double d = PDM_MIN (PDM_MIN (d1,d2), d3); */
    /*   d = d * d; */
    /*   if (PDM_ABS(closest_vertices_dist2[i] - d) > 1e-6) { */
    /*     /\* printf ("Erreur distance plus proche somm %d %d %ld (%12.5e %12.5e %12.5e) : %12.5e %12.5e\n", rank, i,pts_g_num_rank[i], *\/ */
    /*     /\*         pts_rank[3*i], pts_rank[3*i+1], pts_rank[3*i+2], closest_vertices_dist2[i], d); *\/ */
    /*     ierr += 1; */
    /*   } */
    /*   /\* else { *\/ */
    /*   /\*   printf ("ok distance 1 %d (%12.5e %12.5e %12.5e) : %12.5e %12.5e\n", i, *\/ */
    /*   /\*           pts_rank[3*i], pts_rank[3*i+1], pts_rank[3*i+2], closest_vertices_dist2[i], d); *\/ */
    /*   /\* } *\/ */
    /* } */

    /* if (ierr > 0) { */
    /*   printf ("Erreur distance plus proche somm pour %d sommets\n", ierr); */
    /*   abort(); */
    /* } */

    //    fin test

    free (closest_vertices_gnum);

    PDM_octree_free (octree_id);

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[UPPER_BOUND_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[UPPER_BOUND_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[UPPER_BOUND_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[UPPER_BOUND_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);


    /***************************************************************************
     *
     * Compute bounding box structure to find candidates closest
     *     to the upper bound distance
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    PDM_dbbtree_t *dbbt = PDM_dbbtree_create (dist->comm, 3);

          int          *nElts   = malloc (sizeof(int) * n_part_mesh);
    const double      **extents = malloc (sizeof(double *) * n_part_mesh);
    const PDM_g_num_t **gNum    = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);

    if (dist->mesh_nodal_id != -1) {

/*     int PDM_Mesh_nodal_n_blocks_get */
/* ( */
/* const int   idx */
/* ); */

/* int * */
/* PDM_Mesh_nodal_blocks_id_get */
/* ( */
/* const int   idx */
/* ); */


/* PDM_Mesh_nodal_elt_t */
/* PDM_Mesh_nodal_block_type_get */
/* ( */
/* const int   idx, */
/* const int   id_block      */
/* ); */

/* void */
/* PDM_Mesh_nodal_block_std_get  */
/* (    */
/* const int            idx, */
/* const int            id_block,      */
/* const int            id_part,  */
/*       PDM_l_num_t  **connec    */
/* );  */

/* int */
/* PDM_Mesh_nodal_block_n_elt_get  */
/* (    */
/* const int            idx, */
/* const int            id_block,      */
/* const int            id_part  */
/* ); */

/* PDM_g_num_t * */
/* PDM_Mesh_nodal_block_g_num_get  */
/* (    */
/* const int            idx, */
/* const int            id_block,      */
/* const int            id_part  */
/* );  */


/* void */
/* PDM_Mesh_nodal_block_poly2d_get  */
/* ( */
/*  const int          idx, */
/*  const int          id_block,  */
/*  const int          id_part,  */
/*        PDM_l_num_t  **connec_idx,    */
/*        PDM_l_num_t  **connec */
/* );  */

    }
    else if (dist->_surf_mesh != NULL) {
      PDM_surf_mesh_compute_faceExtentsMesh (dist->_surf_mesh, 1e-8);
      for (int i_part = 0; i_part < n_part_mesh; i_part++) {
        nElts[i_part] = PDM_surf_mesh_part_n_face_get (dist->_surf_mesh,
                                                       i_part);

        gNum[i_part] = PDM_surf_mesh_part_face_g_num_get (dist->_surf_mesh,
                                                          i_part);

        extents[i_part] = PDM_surf_mesh_part_extents_get (dist->_surf_mesh,
                                                          i_part);

      }
    }

    else {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_dist_cloud_surf error : The surface mesh is not defined."
                " To do that : \n"
                "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
                "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
                " PDM_dist_cloud_surf_surf_mesh_part_set\n");
    }

    PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set (dbbt,
                                                             n_part_mesh,
                                                             nElts,
                                                             extents,
                                                             gNum);

    if (idebug) {
      printf ("surf_mesh_boxes->n_boxes : %d\n", PDM_box_set_get_size (surf_mesh_boxes));
      for (int i_part = 0; i_part < n_part_mesh; i_part++) {
        printf (" PDM_dbbtree_boxes_set nElts %d : %d\n", i_part, nElts[i_part]);
        for (int i = 0; i < nElts[i_part]; i++) {
          printf ("%d : extents %12.5e %12.5e %12.5e / %12.5e %12.5e %12.5e gnum "PDM_FMT_G_NUM"\n", i,
                  extents[i_part][6*i  ], extents[i_part][6*i+1], extents[i_part][6*i+2],
                  extents[i_part][6*i+3], extents[i_part][6*i+4], extents[i_part][6*i+5],
                  gNum[i_part][i]);
        }
      }
    }

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[BBTREE_CREATE] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[BBTREE_CREATE]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[BBTREE_CREATE]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[BBTREE_CREATE]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    /*
     * Find elements closer than closest_vertices_dist2 distance
     */

    int         *box_index;
    PDM_g_num_t *box_g_num;

    PDM_dbbtree_closest_upper_bound_dist_boxes_get (dbbt,
                                                    n_pts_rank,
                                                    pts_rank,
                                                    pts_g_num_rank,
                                                    closest_vertices_dist2,
                                                    &box_index,
                                                    &box_g_num);

    if (idebug) {
      printf (" PDM_dbbtree_closest_upper_bound_dist_boxes_get n_pts_rank : %d\n", n_pts_rank);
      for (int i = 0; i < n_pts_rank; i++) {
        printf (PDM_FMT_G_NUM" : (%12.5e %12.5e %12.5e) %12.5e\n", pts_g_num_rank[i],
                pts_rank[3*i], pts_rank[3*i+1], pts_rank[3*i+2],
                closest_vertices_dist2[i]);
        printf ("  boxes %d :" , box_index[i+1] - box_index[i]);
        for (int j = box_index[i]; j < box_index[i+1]; j++) {
          printf (" "PDM_FMT_G_NUM, box_g_num[j]);
        }
        printf ("\n");
      }
    }

    free (closest_vertices_dist2);

    PDM_dbbtree_free (dbbt);

    free (nElts);
    free (gNum);
    free (extents);

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[CANDIDATE_SELECTION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[CANDIDATE_SELECTION]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[CANDIDATE_SELECTION]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[CANDIDATE_SELECTION]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /***************************************************************************
     *
     * Load balancing of elementary computations
     * (distance from a point to an element)
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    double *box_n = malloc (sizeof(double) * n_pts_rank);
    int *i_box_n = malloc (sizeof(int) * n_pts_rank);

    for (int i = 0; i < n_pts_rank; i++) {
      box_n[i] = (double) (box_index[i+1] - box_index[i]);
    }

    for (int i = 0; i < n_pts_rank; i++) {
      i_box_n[i] = box_index[i+1] - box_index[i];
    }

    free (box_index);

    PDM_part_to_block_t *ptb_vtx =
      PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                PDM_PART_TO_BLOCK_POST_MERGE,
                                1.,
                                &pts_g_num_rank,
                                &box_n,
                                &n_pts_rank,
                                1,
                                comm);

    int *stride = malloc (sizeof(int) * n_pts_rank);
    for (int i = 0; i < n_pts_rank; i++) {
      stride[i] = 3;
    }

    int *block_stride;
    double *block_pts;
    PDM_part_to_block_exch (ptb_vtx,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &stride,
                            (void **) &pts_rank,
                            &block_stride,
                            (void **) &block_pts);

    free (pts_rank);
    free (stride);

    int *block_g_num_stride;
    PDM_g_num_t *block_g_num;
    PDM_part_to_block_exch (ptb_vtx,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &i_box_n,
                            (void **) &box_g_num,
                            &block_g_num_stride,
                            (void **) &block_g_num);

    free (i_box_n);
    free (box_n);
    free (box_g_num);

    const int n_block_vtx = PDM_part_to_block_n_elt_block_get (ptb_vtx);

    int *block_g_num_idx = malloc (sizeof(int) * (n_block_vtx+1));
    block_g_num_idx[0] = 0;

    /*
     * Remove double coordinates
     */

    int idx = 0;
    int idx2 = 0;
    for (int i = 0; i < n_block_vtx; i++) {
      for (int j = 0; j < 3; j++) {
        block_pts[idx2++] = block_pts[idx++];
      }
      for (int j = 3; j < block_stride[i]; j++) {
        idx += 1;
      }
    }
    free (block_stride);

    /*
     * Remove double boxes
     */


    for (int i = 0; i < n_block_vtx; i++) {
      block_g_num_idx[i+1] = block_g_num_stride[i] + block_g_num_idx[i];
      block_g_num_stride[i] = 0;
    }

    int *block_g_num_opt_idx = malloc (sizeof(int) * (n_block_vtx+1));
    int keyMax = 3 * n_block_vtx;
    PDM_hash_tab_t * ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                               &keyMax);

    block_g_num_opt_idx[0] = 0;
    idx2 = 0;
    for (int i = 0; i <  n_block_vtx; i++) {
      block_g_num_opt_idx[i+1] = block_g_num_opt_idx[i];
      for (int j = block_g_num_idx[i]; j < block_g_num_idx[i+1]; j++) {
        PDM_g_num_t curr_box = block_g_num[j];
        int key = (int) (curr_box % keyMax);
        int n_data = PDM_hash_tab_n_data_get (ht, &key);

        int found = 0;

        PDM_g_num_t **data = (PDM_g_num_t **) PDM_hash_tab_data_get (ht, &key);
        for (int k = 0; k < n_data; k++) {
          if (*(data[k]) == curr_box) {
            found = 1;
            break;
          }
        }

        if (!found) {
          PDM_hash_tab_data_add (ht, (void *) &key, block_g_num + idx2);
          block_g_num[idx2++] = curr_box;
          block_g_num_opt_idx[i+1] += 1;
          block_g_num_stride[i] += 1;
        }
      }
      //PDM_hash_tab_purge (ht, PDM_FALSE);
      PDM_hash_tab_purge (ht, PDM_FALSE);
    }

    PDM_hash_tab_free (ht);

    free (block_g_num_idx);

    block_g_num = realloc(block_g_num, sizeof(PDM_g_num_t) * block_g_num_opt_idx[n_block_vtx]);

    if (idebug) {
      PDM_g_num_t *block_vtx_gnum = PDM_part_to_block_block_gnum_get (ptb_vtx);

      printf ("\n\n **** vtx load balancing : %d\n", n_block_vtx);
      for (int i = 0; i < n_block_vtx; i++) {
        printf (PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e\n", block_vtx_gnum[i], block_pts[3*i], block_pts[3*i+1] , block_pts[3*i+2]);
        printf ("  boxes %d :" , block_g_num_opt_idx[i+1] - block_g_num_opt_idx[i]);
        for (int j = block_g_num_opt_idx[i]; j < block_g_num_opt_idx[i+1]; j++) {
          printf (" "PDM_FMT_G_NUM, block_g_num[j]);
        }
        printf ("\n");
      }
    }

    int block_g_num_n = block_g_num_opt_idx[n_block_vtx];

    /*
     * Receive of needed elements
     *    - PDM_part_to_part function have to be write
     *    - This step is realised by a couple of
     *         PDM_part_to_block and  PDM_block_to_part
     *
     */

    /* part to block */

    const PDM_g_num_t **gnum_face_mesh =
      malloc (sizeof(PDM_g_num_t *) * n_part_mesh);
    int *n_face_mesh = malloc (sizeof(int *) * n_part_mesh);

    for (int i = 0; i < n_part_mesh; i++) {
      n_face_mesh[i] = PDM_surf_mesh_part_n_face_get (dist->_surf_mesh, i);
      gnum_face_mesh[i] = PDM_surf_mesh_part_face_g_num_get(dist->_surf_mesh, i);
    }

    PDM_part_to_block_t *ptb_elt =
      PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                PDM_PART_TO_BLOCK_POST_CLEANUP,
                                1.,
                                (PDM_g_num_t **) gnum_face_mesh,
                                NULL,
                                n_face_mesh,
                                n_part_mesh,
                                comm);

    const int n_block_elt = PDM_part_to_block_n_elt_block_get (ptb_elt);

    double **coords_face_mesh = malloc (sizeof(double *) * n_part_mesh);
    int **coords_face_mesh_n = malloc (sizeof(int *) * n_part_mesh);

    for (int i = 0; i < n_part_mesh; i++) {
      coords_face_mesh[i] = NULL;
      const int *part_face_vtx     =
        PDM_surf_mesh_part_face_vtx_get (dist->_surf_mesh, i);
      const int *part_face_vtx_idx =
        PDM_surf_mesh_part_face_vtx_idx_get (dist->_surf_mesh, i);
      const double *part_vtx = PDM_surf_mesh_part_vtx_get (dist->_surf_mesh, i);

      coords_face_mesh_n[i] = malloc (sizeof(int) * n_face_mesh[i]);
      coords_face_mesh[i] =
        malloc(sizeof(double) * 3 * part_face_vtx_idx[n_face_mesh[i]]);

      for (int j = 0; j < n_face_mesh[i]; j++) {
        coords_face_mesh_n[i][j] =
          (part_face_vtx_idx[j+1] - part_face_vtx_idx[j]) * 3;
      }

      int idx3 = 0;
      for (int j = 0; j < part_face_vtx_idx[n_face_mesh[i]]; j++) {
        int _vtx = part_face_vtx[j] - 1;
        for (int k = 0; k < 3; k++) {
          coords_face_mesh[i][idx3++] = part_vtx[3*_vtx+k];
        }
      }

    }

    int *block_coords_face_mesh_n = NULL;
    double *block_coords_face_mesh = NULL;

    PDM_part_to_block_exch (ptb_elt,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            -1,
                            coords_face_mesh_n,
                            (void **) coords_face_mesh,
                            &block_coords_face_mesh_n,
                            (void **) &block_coords_face_mesh);

    if (idebug) {
      int idx3 = 0;
      PDM_g_num_t *block_elt_gnum = PDM_part_to_block_block_gnum_get (ptb_elt);
      printf ("\n\n **** part to block elt : %d\n", n_block_elt);
      for (int i = 0; i < n_block_elt; i++) {
        printf (PDM_FMT_G_NUM" :\n", block_elt_gnum[i]);
        for (int j = 0; j < block_coords_face_mesh_n[i]/3; j++) {
          printf ("/ %12.5e %12.5e %12.5e /\n",
                  block_coords_face_mesh[3*(idx3+j)],
                  block_coords_face_mesh[3*(idx3+j)+1],
                  block_coords_face_mesh[3*(idx3+j)+2]);
        }
        printf ("\n");
        idx3 += block_coords_face_mesh_n[i]/3;
      }
    }

    /* block to part */

    PDM_g_num_t *block_face_distrib_idx =
      PDM_part_to_block_distrib_index_get (ptb_elt);

    PDM_block_to_part_t *btp =
      PDM_block_to_part_create (block_face_distrib_idx,
                                (const PDM_g_num_t **) &block_g_num,
                                &block_g_num_n,
                                1,
                                comm);

    int un = 1;
    int *part_coords_vtx_face_n = malloc (sizeof(int) * block_g_num_n);
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &un,
                            block_coords_face_mesh_n,
                            NULL,
                            (void **) &part_coords_vtx_face_n);

    int *part_coords_vtx_face_idx = malloc (sizeof(int) * (block_g_num_n+1));
    part_coords_vtx_face_idx[0] = 0;

    for (int i = 0; i < block_g_num_n; i++) {
      part_coords_vtx_face_idx[i+1] = part_coords_vtx_face_idx[i] +
                                      part_coords_vtx_face_n[i]/3;
    }

    double *part_coords_vtx_face =
      malloc (sizeof(double) * 3 * part_coords_vtx_face_idx[block_g_num_n]);

    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_coords_face_mesh_n,
                            block_coords_face_mesh,
                            &part_coords_vtx_face_n,
                            (void **) &part_coords_vtx_face);

    /* free */

    for (int i = 0; i < n_part_mesh; i++) {
      free (coords_face_mesh[i]);
      free (coords_face_mesh_n[i]);
    }

    free (coords_face_mesh);
    free (coords_face_mesh_n);
    free (gnum_face_mesh);
    free (n_face_mesh);

    free (block_coords_face_mesh_n);
    free (block_coords_face_mesh);

    PDM_part_to_block_free (ptb_elt);

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[LOAD_BALANCING_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[LOAD_BALANCING_ELEM_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[LOAD_BALANCING_ELEM_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[LOAD_BALANCING_ELEM_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /***************************************************************************
     *
     * compute distance min per points
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    double *block_closest_dist = malloc (sizeof(double) * n_block_vtx);
    double *block_closest_proj = malloc (sizeof(double) * 3 * n_block_vtx);
    PDM_g_num_t *block_closest_gnum =
      malloc (sizeof(PDM_g_num_t) *  n_block_vtx);

    idx = 0;

    if (idebug) {
      printf ("\n\n****   compute distance min per points   ****\n");
    }

    PDM_g_num_t *block_vtx_gnum = PDM_part_to_block_block_gnum_get (ptb_vtx);
    for (int i = 0; i < n_block_vtx; i++) {
      double *_pt_coords = block_pts + 3*i;
      double *_block_closest_proj = block_closest_proj + 3*i;
      double *_block_closest_dist = block_closest_dist + i;
      _block_closest_dist[0] = HUGE_VAL;

      if (idebug) {
        printf ("    *** "PDM_FMT_G_NUM" : \n", block_vtx_gnum[i]);
      }

      PDM_g_num_t *_block_closest_gnum = block_closest_gnum + i;

      for (int j = 0; j < block_g_num_opt_idx[i+1] - block_g_num_opt_idx[i]; j++) {
        int n_vtx_elt = (part_coords_vtx_face_idx[idx+1] -
                         part_coords_vtx_face_idx[idx]);

        double *_coords_face_elt = part_coords_vtx_face + 3 *
                                   part_coords_vtx_face_idx[idx];

        PDM_g_num_t face_g_num = block_g_num[idx];

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
            idx += 1;
            continue;
          }
        }

        else {

          if (idebug) {
            printf ("_pt_coords : %12.5e %12.5e %12.5e\n",
                    _pt_coords[0], _pt_coords[1], _pt_coords[2]);

            printf ("_coords_face_elt %d %d %d :\n", n_vtx_elt, part_coords_vtx_face_idx[idx], part_coords_vtx_face_idx[idx+1]);
            for (int k = 0; k < n_vtx_elt; k++) {
              printf (" / %12.5e %12.5e %12.5e /\n",
                      _coords_face_elt[3*k],
                      _coords_face_elt[3*k+1],
                    _coords_face_elt[3*k+2]);
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
            idx += 1;
            continue;
          }

        }

        if (minDist2 < _block_closest_dist[0]) {
          _block_closest_dist[0] = minDist2;
          for (int k = 0; k < 3; k++) {
            _block_closest_proj[k] = closestPoint[k];
          }
          _block_closest_gnum[0] = face_g_num;
        }

        idx += 1;
      }
    }

    free (block_g_num_opt_idx);

    free (block_pts);
    free (block_g_num);
    free (block_g_num_stride);

    /* Free */

    free (part_coords_vtx_face_n);
    free (part_coords_vtx_face_idx);
    free (part_coords_vtx_face);

    PDM_block_to_part_free (btp);

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[COMPUTE_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[COMPUTE_ELEM_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[COMPUTE_ELEM_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[COMPUTE_ELEM_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /***************************************************************************
     *
     * Transfer results
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    PDM_g_num_t *block_vtx_distrib_idx =
      PDM_part_to_block_distrib_index_get (ptb_vtx);

    PDM_block_to_part_t *btp_vtx =
      PDM_block_to_part_create (block_vtx_distrib_idx,
                                (const PDM_g_num_t **) pt_cloud->gnum,
                                pt_cloud->n_points,
                                n_part,
                                comm);

    for (int i = 0; i < n_part; i++) {
      int npts = pt_cloud->n_points[i];
      pt_cloud->dist[i] = malloc (sizeof(double) * npts);
      pt_cloud->proj[i] = malloc (sizeof(double) * npts * 3);
      pt_cloud->closest_elt_gnum[i] = malloc (sizeof(PDM_g_num_t) * npts);
    }

    PDM_block_to_part_exch (btp_vtx,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &un,
                            block_closest_dist,
                            NULL,
                            (void **) pt_cloud->dist);

    int three = 3;
    PDM_block_to_part_exch (btp_vtx,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &three,
                            block_closest_proj,
                            NULL,
                            (void **) pt_cloud->proj);


    PDM_block_to_part_exch (btp_vtx,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            &un,
                            block_closest_gnum,
                            NULL,
                            (void **) pt_cloud->closest_elt_gnum);

    PDM_block_to_part_free (btp_vtx);
    PDM_part_to_block_free (ptb_vtx);

    free (block_closest_proj);
    free (block_closest_dist);
    free (block_closest_gnum);

    free (pts_g_num_rank);

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[RESULT_TRANSMISSION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[RESULT_TRANSMISSION]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[RESULT_TRANSMISSION]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[RESULT_TRANSMISSION]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

  }

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
 * \param [in]   id                Identifier
 * \param [in]   i_point_cloud     Current cloud
 * \param [in]   i_part            Index of partition of the cloud
 * \param [out]  distance          Distance
 * \param [out]  projected         Projected point coordinates
 * \param [out]  closest_elt_g_num Global number of the closest element
 *
 */

void
PDM_dist_cloud_surf_get
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
       double      **distance,
       double      **projected,
       PDM_g_num_t **closest_elt_gnum
)
{
 _PDM_dist_t *dist = _get_from_id (id);

 *distance = dist->points_cloud[i_point_cloud].dist[i_part];
 *projected = dist->points_cloud[i_point_cloud].proj[i_part];
 *closest_elt_gnum = dist->points_cloud[i_point_cloud].closest_elt_gnum[i_part];
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
PDM_dist_cloud_surf_free
(
 const int id,
 const int partial
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  if (!partial) {
    for (int i_point_cloud = 0;
         i_point_cloud < dist->n_point_cloud;
         i_point_cloud++) {

      for (int i = 0; i < (dist->points_cloud[i_point_cloud]).n_part; i++) {
        free (dist->points_cloud[i_point_cloud].dist[i]);
        free (dist->points_cloud[i_point_cloud].proj[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_gnum[i]);
      }
    }
  }

  for (int i_point_cloud = 0;
       i_point_cloud < dist->n_point_cloud;
       i_point_cloud++) {

    free (dist->points_cloud[i_point_cloud].n_points);
    free (dist->points_cloud[i_point_cloud].coords);
    free (dist->points_cloud[i_point_cloud].gnum);
    free (dist->points_cloud[i_point_cloud].dist);
    free (dist->points_cloud[i_point_cloud].proj);
    free (dist->points_cloud[i_point_cloud].closest_elt_gnum);
  }

  free (dist->points_cloud);

  PDM_timer_free(dist->timer);

  if (dist->_surf_mesh != NULL) {
    if (dist->surf_mesh != NULL) {
      PDM_surf_mesh_free (dist->surf_mesh);
    }
  }

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
PDM_dist_cloud_surf_dump_times
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


    PDM_printf( "distance timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "distance timer : Upper bound distance (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[UPPER_BOUND_DIST],
                t_cpu_max[UPPER_BOUND_DIST]);
    PDM_printf( "distance timer : Bbtree building (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BBTREE_CREATE],
                t_cpu_max[BBTREE_CREATE]);
    PDM_printf( "distance timer : Candidate selection (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[CANDIDATE_SELECTION],
                t_cpu_max[CANDIDATE_SELECTION]);
    PDM_printf( "distance timer : Load balacing of elementary computations of distance"
                " from the points to the candidates  (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[LOAD_BALANCING_ELEM_DIST],
                t_cpu_max[LOAD_BALANCING_ELEM_DIST]);
    PDM_printf( "distance timer : Computations of the distance"
                " from the points to the candidates  (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[COMPUTE_ELEM_DIST],
                t_cpu_max[COMPUTE_ELEM_DIST]);
    PDM_printf( "distance timer : Results exchange (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[RESULT_TRANSMISSION],
                t_cpu_max[RESULT_TRANSMISSION]);
    PDM_printf_flush();
  }
}


#ifdef	__cplusplus
}
#endif
