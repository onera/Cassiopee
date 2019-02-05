#ifndef __PDM_MESH_NODAL_H__
#define __PDM_MESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Geometric element type
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_MESH_NODAL_POINT,     
  PDM_MESH_NODAL_BAR2,     
  PDM_MESH_NODAL_TRIA3,     
  PDM_MESH_NODAL_QUAD4,     
  PDM_MESH_NODAL_POLY_2D,     
  PDM_MESH_NODAL_TETRA4,     
  PDM_MESH_NODAL_PYRAMID5,     
  PDM_MESH_NODAL_PRISM6,     
  PDM_MESH_NODAL_HEXA8,     
  PDM_MESH_NODAL_POLY_3D     

} PDM_Mesh_nodal_elt_t;


typedef struct _PDM_Mesh_nodal_t PDM_Mesh_nodal_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

int 
PDM_Mesh_nodal_create
(
const int     n_part,
const PDM_MPI_Comm comm        
);

/**
 * \brief Free partially a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_partial_free
(
const int idx
);

/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_free
(
const int idx
);

/**
 * \brief Define partition vertices
 *
 * \param [in]  idx   Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  numabs    Global numbering
 *
 */

void
PDM_Mesh_nodal_coord_set
(
 const int          idx,
 const int          id_part, 
 const int          n_vtx,  
 const PDM_real_t  *coords,  
 const PDM_g_num_t *numabs
);


/**
 * \brief  Return number of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vertices_get
(
 const int          idx,
 const int          id_part 
);


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_Mesh_nodal_vertices_get
(
 const int          idx,
 const int          id_part 
);


/**
 * \brief  Return parent num of vertices
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Parent of vertices
 *
 */

const int *
PDM_Mesh_nodal_vertices_parent_get
(
 const int          idx,
 const int          id_part 
 );


/**
 * \brief  Return global numbering of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Global numbering of vertices
 *
 */

const PDM_g_num_t *
PDM_Mesh_nodal_vertices_g_num_get
(
 const int          idx,
 const int          id_part 
);


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return true if the vertices are defined from parents
 */

int
PDM_Mesh_nodal_is_set_coord_from_parent
(
 const int          idx
 );


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_vtx          Number of vertices
 * \param [in]  n_vtx_parent   Number of parent vertices
 * \param [in]  numabs         Global numbering (size = \ref n_vtx)
 * \param [in]  num_parent     Numbering in the parent numbering (size = \ref n_vtx)
 * \param [in]  coords_parent  Parent interlaced coordinates (size = 3 * \ref n_vtx_parent)
 * \param [in]  numabs_parent  Parent global numbering (size = \ref n_vtx_parent)
 *
 */

void
PDM_Mesh_nodal_coord_from_parent_set
(
 const int          idx,
 const int          id_part, 
 const int          n_vtx,  
 const int          n_vtx_parent,  
 const PDM_g_num_t *numabs,
 const int         *num_parent,
 const PDM_real_t  *coords_parent,  
 const PDM_g_num_t *numabs_parent
);


/**
 * \brief  Return number of blocks
 *
 * \param [in]  idx            Nodal mesh handle
 *
 * \return  Number of blocks
 *
 */

int
PDM_Mesh_nodal_n_blocks_get
(
const int   idx
);


/**
 * \brief  Return blocks identifier
 *
 * \param [in]  idx            Nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */

int *
PDM_Mesh_nodal_blocks_id_get
(
const int   idx
);


/**
 * \brief  Return number of partitions
 *
 * \param [in]  idx            Nodal mesh handle
 *
 * \return  Number of partitions
 *
 */

int
PDM_Mesh_nodal_n_part_get
(
const int   idx
);


/**
 * \brief  Return type of block
 *
 * \param [in]  idx        Nodal mesh handle
 * \param [in]  id_block   Block identifier
 *
 * \return  Type of block
 *
 */

PDM_Mesh_nodal_elt_t
PDM_Mesh_nodal_block_type_get
(
const int   idx,
const int   id_block     
);


/**
 * \brief  Add a new block to the current mesh
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory 
 *                             when the block is destroyed
 * \param [in]  id_block       Block identifier
 *
 * \return Block identifier     
 *
 */

int 
PDM_Mesh_nodal_block_add 
(
const int                    idx,
PDM_bool_t                   st_free_data,  
const PDM_Mesh_nodal_elt_t   t_elt
); 


/**
 * \brief Define a standard block
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_std_set 
(
const int          idx,
const int            id_block,     
const int            id_part, 
const int            n_elt,    
      PDM_l_num_t   *connec,   
      PDM_g_num_t   *numabs,
      PDM_l_num_t   *parent_num   
); 


/**
 * \brief Return standard block description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  connect        Connectivity
 *
 */

void
PDM_Mesh_nodal_block_std_get 
(   
const int            idx,
const int            id_block,     
const int            id_part, 
      PDM_l_num_t  **connec   
); 


/**
 * \brief Get number of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Number of elements
 *  
 */

int
PDM_Mesh_nodal_block_n_elt_get 
(   
const int            idx,
const int            id_block,     
const int            id_part 
); 


/**
 * \brief Get global numbering of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global numbering of block elements
 *  
 */

PDM_g_num_t *
PDM_Mesh_nodal_block_g_num_get 
(   
const int            idx,
const int            id_block,     
const int            id_part 
); 


/**
 * \brief Get global inside numbering of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global inside numbering of block elements
 *  
 */

PDM_g_num_t *
PDM_Mesh_nodal_block_inside_g_num_get 
(   
const int            idx,
const int            id_block,     
const int            id_part 
); 


/**
 * \brief Get parent numbering of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return parent numbering of block elements
 *  
 */

int *
PDM_Mesh_nodal_block_parent_num_get 
(   
const int            idx,
const int            id_block,     
const int            id_part 
); 


/**
 * \brief Define a polygon block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */
 
void
PDM_Mesh_nodal_block_poly2d_set 
(
const int            idx,
const int            id_block, 
const int            id_part, 
const PDM_l_num_t    n_elt,    
      PDM_l_num_t   *connec_idx,   
      PDM_l_num_t   *connec,
      PDM_g_num_t   *numabs,
      PDM_l_num_t   *parent_num
); 


/**
 * \brief Return a polygon block description
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] n_elt          Number of elements
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [out] numabs         Global numbering in the mesh
 * \param [out] numabs_block   Global numbering in the block or NULL (if not computed)
 * \param [out] parent_num     Parent numbering or NULL
 *
 */
 
void
PDM_Mesh_nodal_block_poly2d_get 
(
 const int          idx,
 const int          id_block, 
 const int          id_part, 
       PDM_l_num_t  **connec_idx,   
       PDM_l_num_t  **connec
); 


/**
 * \brief Define a polyhedra block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_poly3d_set 
(
const int            idx,
const int            id_block, 
const int            id_part, 
const PDM_l_num_t    n_elt,    
const PDM_l_num_t    n_face,   
      PDM_l_num_t   *facvtx_idx,   
      PDM_l_num_t   *facvtx,
      PDM_l_num_t   *cellfac_idx,   
      PDM_l_num_t   *cellfac,
      PDM_g_num_t   *numabs,
      PDM_l_num_t   *parent_num
); 


/**
 * \brief Define a polyhedra block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  n_elt          Number of polyhedra
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 * \param [out]  numabs         Global numbering
 * \param [out] numabs_block   Global numbering in the block or NULL (if not computed)
 * \param [out] parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_poly3d_get 
(
const int            idx,
const int            id_block, 
const int            id_part, 
      PDM_l_num_t   *n_face,   
      PDM_l_num_t  **facvtx_idx,   
      PDM_l_num_t  **facvtx,
      PDM_l_num_t  **cellfac_idx,   
      PDM_l_num_t  **cellfac
); 

/**
 * \brief  Add some 3D cells from cell face conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection 
 * between old and new numbering.
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  face_vtx_idx   Index of face vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each face
 * \param [in]  face_vtx       Face vertex connectivity
 * \param [in]  cell_face_idx  Index of cell face connectivity
 * \param [in]  cell_face_nb   Number of faces for each cell
 * \param [in]  cell_face      Cell face connectivity
 * \param [in]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_cell3d_cellface_add
(
const int         idx,
const int         id_part, 
const int         n_elt,
const int         n_face,
PDM_l_num_t      *face_vtx_idx,
PDM_l_num_t      *face_vtx_nb,
PDM_l_num_t      *face_vtx,
PDM_l_num_t      *cell_face_idx,
PDM_l_num_t      *cell_face_nb,
PDM_l_num_t      *cell_face,
PDM_g_num_t      *numabs
); 


/**
 * \brief  Add some 2D cells from cell edge conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection 
 * between old and new numbering.
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_edge         Number of edges used to describe polyhedra
 * \param [in]  edge_vtx_idx   Index of edge vertex connectivity
 * \param [in]  edge_vtx_nb    Number of vertices for each edge
 * \param [in]  edge_vtx       Edge vertex connectivity
 * \param [in]  cell_edge_idx  Index of cell edge connectivity
 * \param [in]  cell_edge_nb   Number of edges for each cell
 * \param [in]  cell_edge      Cell edge connectivity
 * \param [in]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_cell2d_celledge_add
(
const int          idx,
const int          id_part, 
const int          n_elt,
const int          n_edge,
PDM_l_num_t       *edge_vtx_idx,
PDM_l_num_t       *edge_vtx_nb,
PDM_l_num_t       *edge_vtx,
PDM_l_num_t       *cell_edge_idx,
PDM_l_num_t       *cell_edge_nb,
PDM_l_num_t       *cell_edge,
PDM_g_num_t       *numabs
); 


/**
 * \brief  Add some 2D cells from cell vertex connectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection 
 * between old and new numbering.
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygon
 * \param [in]  face_vtx_idx   Index of edge vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each edge
 * \param [in]  face_vtx       Edge vertex connectivity
 * \param [in]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_faces_facevtx_add
(
const int         idx,
const int         id_part, 
const int         n_face,
PDM_l_num_t      *face_vtx_idx,
PDM_l_num_t      *face_vtx_nb,
PDM_l_num_t      *face_vtx,
PDM_g_num_t      *numabs
); 


/**
 * \brief  Compute a global numbering in a block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 *
 */

void
PDM_Mesh_nodal_g_num_in_block_compute
(
const int         idx,
const int         id_block 
); 


/**
 * \brief  Return parent cell number to local number
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Parent cell number to local number
 * 
 */

int *
PDM_Mesh_nodal_num_cell_parent_to_local_get
(
const int  idx,
const int  id_part 
); 


/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Return number elements of a partition
 * 
 */

int
PDM_Mesh_nodal_n_cell_get
(
const int  idx,
const int  id_part 
); 


/**
 * \brief  Return parent  absolute number
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Parent of vertices
 *
 */

const PDM_g_num_t *
PDM_Mesh_nodal_vertices_g_num_parent_get
(
 const int          idx,
 const int          id_part 
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
