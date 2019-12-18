#ifndef __PDM_MESH_NODAL_PRIV_H__
#define __PDM_MESH_NODAL_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_handles.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_Mesh_nodal_som_t
 * \brief  Vertices of a mesh partition
 *
 */

typedef struct PDM_Mesh_nodal_vtx_t PDM_Mesh_nodal_vtx_t;

struct PDM_Mesh_nodal_vtx_t {

  PDM_Mesh_nodal_vtx_t *parent;    /*!< Parent vertices if they are extraxted from an other mesh partition */
  PDM_l_num_t           n_vtx;     /*!< Number of vertices */
  double               *coords;    /*!< Coordinates (Locally allocated) */
  const double        *_coords;    /*!< Coordinates (Mapping) */
  const PDM_g_num_t   *_numabs;    /*!< Global numbering */
  const int           *_numparent; /*!< Numbering in the parent vertices (mapping) */
};

/**
 * \struct PDM_Mesh_nodal_block_std_t
 * \brief  Standard geometric block
 *
 */

typedef struct PDM_Mesh_nodal_block_std_t {

  PDM_Mesh_nodal_elt_t       t_elt;         /*!< Element type */
  PDM_bool_t                 st_free_data;  /*!< Release of memory during the destruction of the object or not */
  PDM_l_num_t                n_part;        /*!< Number of partitions */
  PDM_l_num_t               *n_elt;         /*!< Number elements */
  PDM_l_num_t              **_connec;       /*!< Connectivity (Memory mapping) */
  PDM_l_num_t              **_num_part;     /*!< Initial numbering int the partition (Memory mapping) */
  PDM_g_num_t              **_numabs;       /*!< Global numbering (Memory mapping) */
  PDM_g_num_t              **numabs_int;    /*!< Global numbering inside each block */
  PDM_l_num_t              **_parent_num;       /*!< Parent numbering or NULL */
  double                   **cell_centers;       /*!< Cell center coordinates */

} PDM_Mesh_nodal_block_std_t;


/**
 * \struct PDM_Mesh_nodal_block_poly2d_t
 * \brief  Polygon geometric block
 *
 */

typedef struct PDM_Mesh_nodal_block_poly2d_t {

  PDM_bool_t               st_free_data;  /*!< Release of memory during the destruction of the object or not */
  PDM_l_num_t              n_part;        /*!< Number of partitions */
  PDM_l_num_t             *n_elt;         /*!< Number of elements of each partition */
  PDM_l_num_t            **_connec_idx;   /*!< Index of elements connectivity of each partition (Memory mapping) */
  PDM_l_num_t            **_connec;       /*!< Elements connectivity of each partition (Memory mapping) */
  PDM_l_num_t            **_num_part;     /*!< Initial numbering in each partition (Memory mapping) */
  PDM_g_num_t            **_numabs;       /*!< Global numbering  (Memory mapping) */
  PDM_g_num_t            **numabs_int;    /*!< Global numbering inside each block */
  PDM_l_num_t            **_parent_num;       /*!< Parent numbering or NULL */
  double                 **cell_centers;       /*!< Cell center coordinates */

} PDM_Mesh_nodal_block_poly2d_t;


/**
 * \struct PDM_Mesh_nodal_block_poly3d_t
 * \brief  Polyhedron geometric block
 *
 */

typedef struct PDM_Mesh_nodal_block_poly3d_t{

  PDM_bool_t              st_free_data; /*!< Release of memory during the destruction of the object or not */
  PDM_l_num_t             n_part;       /*!< Number of partitions */
  PDM_l_num_t            *n_elt;        /*!< Number of elements of each partition */
  PDM_l_num_t            *n_face;       /*!< Number of face of each polyhedron of each partition */
  PDM_l_num_t           **_facvtx_idx;  /*!< Index of faces connectivity of each partition (Memory mapping) */

  PDM_l_num_t           **_facvtx;      /*!< Faces connectivity of each partition (Memory mapping) */
  PDM_l_num_t           **_cellfac_idx; /*!< Index of cell->face connectivity (Memory mapping) */

  PDM_l_num_t           **_cellfac;     /*!< cell->face connectivity (Memory mapping) */

  PDM_l_num_t           **_num_part;    /*!< Initial numbering in the partition (Memory mapping) */

  PDM_g_num_t           **_numabs;      /*!< Global numbering (Memory mapping) */

  PDM_g_num_t           **numabs_int;   /*!< Global numbering inside the block (Memory mapping) */
  PDM_l_num_t           **_parent_num;  /*!< Parent numbering or NULL */
  double                **cell_centers;       /*!< Cell center coordinates */

} PDM_Mesh_nodal_block_poly3d_t;


/**
 * \struct  PDM_Mesh_nodal_prepa_blocks_t
 *
 * \brief   Used to build blocks from cell to face face to edge connectivity
 *
 */

typedef struct PDM_Mesh_nodal_prepa_blocks_t {

  PDM_l_num_t  n_tria_proc;       /*!< Number of triangles per proc */
  PDM_l_num_t  n_quad_proc;       /*!< Number of quadrangles per proc */
  PDM_l_num_t  n_poly2d_proc;     /*!< Number of polygons per proc */
  PDM_l_num_t  n_tetra_proc;      /*!< Number of tetrahedra per proc */
  PDM_l_num_t  n_hexa_proc;       /*!< Number of hexahedra per proc */
  PDM_l_num_t  n_prism_proc;      /*!< Number of prisms per proc */
  PDM_l_num_t  n_pyramid_proc;    /*!< Number of pyramids per proc */
  PDM_l_num_t  n_poly3d_proc;     /*!< Number of poly3d per proc */
  PDM_l_num_t *add_etat;          /*!< Allows to check if all partitions are taking into account */
  PDM_l_num_t  t_add;             /*!< Type of input (1 : cell3d_cellface,
                                                      2 : cell2d_cellface,
                                                      3 : faces_facesvtx_add) */
  PDM_l_num_t  *n_tetra;          /*!< Number of tetrahedra per partition */
  PDM_l_num_t  *n_hexa;           /*!< Number of hexhedra per partition */
  PDM_l_num_t  *n_prism;          /*!< Number of prisms per partition */
  PDM_l_num_t  *n_pyramid;        /*!< Number of pyramids per partition */
  PDM_l_num_t  *n_poly3d;         /*!< Number of polyhedra per partition */
  PDM_l_num_t  *n_tria;           /*!< Number of triangles per partition */
  PDM_l_num_t  *n_quad;           /*!< Number of quadrangles per partition */
  PDM_l_num_t  *n_poly2d;         /*!< Number of polygons per partition */
  PDM_l_num_t  *l_connec_poly2d;  /*!< Size of polygon connectivity per partition */
  PDM_l_num_t  **face_vtx_idx;    /*!< Index of face vertex connectivity */
  PDM_l_num_t  **face_vtx_nb;     /*!< Number of vertex per face */
  PDM_l_num_t  **face_vtx;        /*!< Face vertex connectivity */
  PDM_l_num_t  *n_cell;           /*!< Number of cells */
  PDM_l_num_t  *n_face;           /*!< Number of faces */
  PDM_l_num_t  **cell_face_idx;   /*!< Index of cell face connectivity */
  PDM_l_num_t  **cell_face_nb;    /*!< Number of faces per cell */
  PDM_l_num_t  **cell_face;       /*!< Cell face connectivity */
  PDM_g_num_t  **numabs;          /*!< Global numbering per cell per partition */

} PDM_Mesh_nodal_prepa_blocks_t;


/**
 * \struct  PDM_Mesh_nodal_geom_prepa_blocks_t
 *
 * \brief   Used to build blocks from cell to face face to edge connectivity
 *
 */

struct _PDM_Mesh_nodal_t {

  /* PDM_g_num_t                         n_som_abs;                /\*!< Global number of vertices *\/ */
  /* PDM_g_num_t                         n_elt_abs;                /\*!< Global number of elements *\/ */
  int                                 n_part;                   /*!< Number of partitions */
  PDM_Mesh_nodal_vtx_t                **vtx;                    /*!< Description des sommmets de chaque partition */
  PDM_l_num_t                         *n_cell;                  /*!< Nombre de blocs d'elements standard */
  PDM_Handles_t                       *blocks_std;              /*!< Standard blocks */
  PDM_Handles_t                       *blocks_poly2d;           /*!< Polygon blocks */
  PDM_Handles_t                       *blocks_poly3d;           /*!< Polyhedron blocks */
  PDM_MPI_Comm                         pdm_mpi_comm;            /*!< MPI Communicator */
  PDM_Mesh_nodal_prepa_blocks_t  *prepa_blocks;            /*!< Blocks preparation */
  PDM_l_num_t                        **num_cell_parent_to_local;/*!< Initial local numbering to local numbering
                                                                 *   imposed by blocks */
  int                                 *blocks_id;               /*!< Blocks identifier */
  int                                  n_blocks;                /*!< Total number of blocks */
  int                      is_vtx_def_from_parent; /*< Are the points defined from parents */
} ;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_PRIV_H__ */
