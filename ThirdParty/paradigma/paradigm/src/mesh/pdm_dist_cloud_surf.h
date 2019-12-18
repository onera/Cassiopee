#ifndef PDM_DIST_CLOUD_SURF_H
#define PDM_DIST_CLOUD_SURF_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_surf_mesh.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

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
 *
 */

int
PDM_dist_cloud_surf_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
);

void
PDM_dist_cloud_surf_create_cf
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Fint comm,
 int *id
);


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
 );


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
 );



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
 );

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
);


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
);


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
);


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
);


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  closest_elt_distance  Distance
 * \param [out]  closest_elt_projected Projected point coordinates
 * \param [out]  closest_elt_g_num     Global number of the closest element
 *
 */

void
PDM_dist_cloud_surf_get
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
       double      **closest_elt_distance,
       double      **closest_elt_projected,
       PDM_g_num_t **closest_elt_gnum
 );


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
 );


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
 );

#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_DIST_H
