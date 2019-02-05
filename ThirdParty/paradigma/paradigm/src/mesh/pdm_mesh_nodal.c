
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
#include "pdm_mesh_nodal.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_BLOCK_ID_BLOCK_STD    = 0,    
  PDM_BLOCK_ID_BLOCK_POLY2D = 1000000,
  PDM_BLOCK_ID_BLOCK_POLY3D = 2000000

} PDM_block_id_block_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * \brief Storage of mesh handles
 *
 */

static PDM_Handles_t *mesh_handles = NULL; 


/*============================================================================
 * Private function definitions
 *============================================================================*/



/**
 * 
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 * 
 */

static 
PDM_Mesh_nodal_vtx_t *
_vtx_free
(
 PDM_Mesh_nodal_vtx_t *vtx
)
{
  
  if (vtx != NULL) {
    if (vtx->parent != NULL) {
      vtx->parent =_vtx_free (vtx->parent);
    }

    if (vtx->coords != NULL) {
      free (vtx->coords);
      vtx->coords = NULL;
    }

    free (vtx);
  }
  return NULL;
}


/**
 * 
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void   
_mesh_init
(
PDM_Mesh_nodal_t *mesh,
const int        n_part,
const PDM_MPI_Comm comm        
)
{

  mesh->n_som_abs                = 0;          
  mesh->n_elt_abs                = 0;          
  mesh->n_part                   = n_part;             

  mesh->vtx                      = malloc(n_part * sizeof(PDM_Mesh_nodal_vtx_t *));
  mesh->n_cell                   = malloc(n_part * sizeof(int));  
  for (int i = 0; i < n_part; i++) {
    mesh->vtx[i] = malloc(sizeof(PDM_Mesh_nodal_vtx_t));
    mesh->vtx[i]->_coords = NULL;
    mesh->vtx[i]->_numabs = NULL;
    mesh->vtx[i]->_numparent = NULL;
    mesh->vtx[i]->n_vtx   = 0;
    mesh->vtx[i]->parent = NULL;
    mesh->vtx[i]->coords = NULL;
    mesh->n_cell[i] = 0;
  }

  mesh->blocks_std               = NULL;  
  mesh->blocks_poly2d            = NULL;            
  mesh->blocks_poly3d            = NULL;           

  mesh->pdm_mpi_comm             = comm;
  mesh->prepa_blocks             = NULL;
  mesh->num_cell_parent_to_local = NULL;
  mesh->blocks_id                = NULL;
  mesh->n_blocks                 = 0;
  mesh->is_vtx_def_from_parent   = 0;
}

/**
 * 
 * \brief Update blocks identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static void   
_update_blocks_id
(
PDM_Mesh_nodal_t *mesh
)
{
  int n_blocks = 0;
  
  if (mesh->blocks_std != NULL) {
    n_blocks += PDM_Handles_n_get (mesh->blocks_std);
  }
    
  if (mesh->blocks_poly2d != NULL) {
    n_blocks += PDM_Handles_n_get (mesh->blocks_poly2d);
  }

  if (mesh->blocks_poly3d != NULL) {
    n_blocks += PDM_Handles_n_get (mesh->blocks_poly3d);
  }
  
  if (mesh->n_blocks < n_blocks) {
    mesh->blocks_id = (int *) realloc(mesh->blocks_id, sizeof(int) * n_blocks);
  }
  
  int k = 0;
  if (mesh->blocks_std != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->blocks_std);
    int n = PDM_Handles_n_get (mesh->blocks_std);
    for (int i = 0; i < n; i++) {
      mesh->blocks_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_STD;
    }
  }
    
  if (mesh->blocks_poly2d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->blocks_poly2d);
    int n = PDM_Handles_n_get (mesh->blocks_poly2d);
    for (int i = 0; i < n; i++) {
      mesh->blocks_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }
  
  if (mesh->blocks_poly3d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->blocks_poly3d);
    int n = PDM_Handles_n_get (mesh->blocks_poly3d);
    for (int i = 0; i < n; i++) {
      mesh->blocks_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }
  
  mesh->n_blocks = n_blocks;  

}


/**
 * 
 * \brief Free partially a standard block
 *
 * \param [inout]  _bloc_std    Standard block
 *  
 */

static
void
_block_std_free_partial
(
PDM_Mesh_nodal_block_std_t *_block_std
)
{

  if (_block_std == NULL) {
    return;
  }
  
  if (_block_std->_connec != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_connec[i] != NULL)
          free(_block_std->_connec[i]);
        _block_std->_connec[i] = NULL;
      }
    }
    free(_block_std->_connec);
    _block_std->_connec = NULL;
  }
  
  if (_block_std->_numabs != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_numabs[i] != NULL)
          free(_block_std->_numabs[i]);
        _block_std->_numabs[i] = NULL;
      }
    }
    free(_block_std->_numabs);
    _block_std->_numabs = NULL;
  }
  
  if (_block_std->_num_part != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_num_part[i] != NULL)
          free(_block_std->_num_part[i]);
        _block_std->_num_part[i] = NULL;
      }
    }
    free(_block_std->_num_part);
    _block_std->_num_part = NULL;
  }

}


/**
 * 
 * \brief Free a standard block
 *
 * \param [inout]  _bloc_std    Standard block
 *
 * \return         Null
 *   
 */

static
PDM_Mesh_nodal_block_std_t *
_block_std_free
(
PDM_Mesh_nodal_block_std_t *_block_std
)
{

  if (_block_std == NULL) {
    return NULL;
  }

  _block_std_free_partial(_block_std);

  if (_block_std->n_elt != NULL) {
    free(_block_std->n_elt);
    _block_std->n_elt = NULL;
  }

  if (_block_std->numabs_int != NULL) {
    for (int j = 0; j < _block_std->n_part; j++) {
      if (_block_std->numabs_int[j] != NULL) {
        free(_block_std->numabs_int[j]);
      }
    }
    free(_block_std->numabs_int);
    _block_std->numabs_int = NULL;
  }

  if (_block_std->_parent_num != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_parent_num[i] != NULL)
          free(_block_std->_parent_num[i]);
        _block_std->_parent_num[i] = NULL;
      }
    }
    free(_block_std->_parent_num);
    _block_std->_parent_num = NULL;
  }

  free(_block_std);
  return NULL;
}


/**
 * 
 * \brief Free partially a polygon block
 *
 * \param [inout]  _block_poly2d   polygon block
 *   
 */

static
void
_block_poly2d_free_partial
(
PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
)
{

  if (_block_poly2d->_connec_idx != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_connec_idx[i] != NULL)
          free(_block_poly2d->_connec_idx[i]);
        _block_poly2d->_connec_idx[i] = NULL;
      }
    }
    free(_block_poly2d->_connec_idx);
    _block_poly2d->_connec_idx = NULL;
  }

  if (_block_poly2d->_connec != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_connec[i] != NULL)
          free(_block_poly2d->_connec[i]);
        _block_poly2d->_connec[i] = NULL;
      }
    }
    free(_block_poly2d->_connec);
    _block_poly2d->_connec = NULL;
  }
  
  if (_block_poly2d->_num_part != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_num_part[i] != NULL)
          free(_block_poly2d->_num_part[i]);
        _block_poly2d->_num_part[i] = NULL;
      }
    }
    free(_block_poly2d->_num_part);
    _block_poly2d->_num_part = NULL;
  }
  
  if (_block_poly2d->_numabs != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_numabs[i] != NULL)
          free(_block_poly2d->_numabs[i]);
        _block_poly2d->_numabs[i] = NULL;
      }
    }
    free(_block_poly2d->_numabs);
    _block_poly2d->_numabs = NULL;
  }

}


/**
 * 
 * \brief Free a polygon block
 *
 * \param [inout]  _bloc_poly2d    Polygon block
 *
 * \return         Null
 *   
 */

static
PDM_Mesh_nodal_block_poly2d_t *
_block_poly2d_free
(
PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
)
{
  _block_poly2d_free_partial(_block_poly2d);
  
  if (_block_poly2d->n_elt != NULL) {
    free(_block_poly2d->n_elt);
    _block_poly2d->n_elt = NULL;
  }

  if (_block_poly2d->numabs_int != NULL) {
    for (int j = 0; j < _block_poly2d->n_part; j++) {
      if (_block_poly2d->numabs_int[j] != NULL) {
        free(_block_poly2d->numabs_int[j]);
      }
    }
    free(_block_poly2d->numabs_int);
    _block_poly2d->numabs_int = NULL;
  }

  if (_block_poly2d->_parent_num != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_parent_num[i] != NULL)
          free(_block_poly2d->_parent_num[i]);
        _block_poly2d->_parent_num[i] = NULL;
      }
    }
    free(_block_poly2d->_parent_num);
    _block_poly2d->_parent_num = NULL;
  }
  
  free(_block_poly2d);

  return NULL;
}


/**
 * 
 * \brief Free partially a polyhedron block
 *
 * \param [inout]  _block_poly3d   polyhedron block
 *   
 */

static
void
_block_poly3d_free_partial
(
PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
)
{
  
  if (_block_poly3d->_facvtx_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_facvtx_idx[i] != NULL)
          free(_block_poly3d->_facvtx_idx[i]);
        _block_poly3d->_facvtx_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_facvtx_idx);
    _block_poly3d->_facvtx_idx = NULL;
  }
  
  if (_block_poly3d->_facvtx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_facvtx[i] != NULL)
          free(_block_poly3d->_facvtx[i]);
        _block_poly3d->_facvtx[i] = NULL;
      }
    }
    free(_block_poly3d->_facvtx);
    _block_poly3d->_facvtx = NULL;
  }
  
  if (_block_poly3d->_cellfac_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellfac_idx[i] != NULL)
          free(_block_poly3d->_cellfac_idx[i]);
        _block_poly3d->_cellfac_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_cellfac_idx);
    _block_poly3d->_cellfac_idx = NULL;
  }
  
  if (_block_poly3d->_cellfac != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellfac[i] != NULL)
          free(_block_poly3d->_cellfac[i]);
        _block_poly3d->_cellfac[i] = NULL;
      }
    }
    free(_block_poly3d->_cellfac);
    _block_poly3d->_cellfac = NULL;
  }
  
  if (_block_poly3d->_numabs != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_numabs[i] != NULL)
          free(_block_poly3d->_numabs[i]);
        _block_poly3d->_numabs[i] = NULL;
      }
    }
    free(_block_poly3d->_numabs);
    _block_poly3d->_numabs = NULL;
  }
}


/**
 * 
 * \brief Free a polyhedron block
 *
 * \param [inout]  _block_poly3d    Polyhedron block
 *
 * \return         Null
 *   
 */

static
PDM_Mesh_nodal_block_poly3d_t *
_block_poly3d_free
(
PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
)
{
  _block_poly3d_free_partial(_block_poly3d);

  if (_block_poly3d->n_elt != NULL) {
    free(_block_poly3d->n_elt);
    _block_poly3d->n_elt = NULL;
  }

  if (_block_poly3d->n_face!= NULL) {
    free(_block_poly3d->n_face);
    _block_poly3d->n_face= NULL;
  }

  if (_block_poly3d->numabs_int != NULL) {
    for (int j = 0; j < _block_poly3d->n_part; j++) {
      if (_block_poly3d->numabs_int[j] != NULL) {
        free(_block_poly3d->numabs_int[j]);
      }
    }
    free(_block_poly3d->numabs_int);
    _block_poly3d->numabs_int = NULL;
  }

  if (_block_poly3d->_parent_num != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_parent_num[i] != NULL)
          free(_block_poly3d->_parent_num[i]);
        _block_poly3d->_parent_num[i] = NULL;
      }
    }
    free(_block_poly3d->_parent_num);
    _block_poly3d->_parent_num = NULL;
  }

  free(_block_poly3d);
  return NULL;
}


/**
 * 
 * \brief  Cross product
 *
 * \param[in]     a    First vector
 * \param[in]     b    Second vector
 * \param[inout]  c    \ref a X \ref b vector
 *   
 */

static inline void
_p_cross
(
 const double a[3],
 const double b[3],
       double c[3]
)
{
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = b[0] * a[2] - a[0] * b[2];
  c[2] = a[0] * b[1] - b[0] * a[1];
}


/**
 *
 * Dot product
 *
 * \param[in]     a    First vector
 * \param[in]     b    Second Vector
 *
 * \return    Dot product
 *   
 */

static inline double
_p_dot 
(
 const double a[3], 
 const double b[3] 
)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


/**
 * 
 * \brief Build tetrahedron nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[out] tetra_vtx   Tetrahedron connectivity
 *
 */

static void
_connec_tetra
(
PDM_Mesh_nodal_vtx_t *vtx,
PDM_l_num_t *tria_vtx,
PDM_l_num_t tetra_vtx[]
)
{

  /* Initialization */

  tetra_vtx[0] = tria_vtx[0];
  tetra_vtx[1] = tria_vtx[1];
  tetra_vtx[2] = tria_vtx[2];

  for (int i = 3; i < 11; i++) {
    if ((tria_vtx[i] != tetra_vtx[0]) &&
        (tria_vtx[i] != tetra_vtx[1]) &&
        (tria_vtx[i] != tetra_vtx[2]))
      tetra_vtx[3] = tria_vtx[i];
  } 

  /* Orientation */

  const PDM_real_t *_coords = vtx->_coords;
  double v1[3];
  double v2[3];
  double n[3];

  for (int i = 0; i < 3; i++) {
    v1[i] = _coords[3*(tria_vtx[1] - 1) + i] - _coords[3*(tria_vtx[0] - 1) + i];
    v2[i] = _coords[3*(tria_vtx[2] - 1) + i] - _coords[3*(tria_vtx[0] - 1) + i];
  }

  _p_cross(v1, v2, n);
  double orient = _p_dot(v1, n);

  if (orient < 0) {
    tetra_vtx[0] = tria_vtx[2];
    tetra_vtx[1] = tria_vtx[1];
    tetra_vtx[2] = tria_vtx[0];
  }
}


/**
 * 
 * \brief Build prism nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] prism_vtx   Prism connectivity
 *
 */

static void
_connec_prism
(
PDM_Mesh_nodal_vtx_t *vtx,
PDM_l_num_t *tria_vtx,
PDM_l_num_t *quad_vtx,
PDM_l_num_t prism_vtx[]
)
{

  /* Initialisation */

  for (int i = 0; i < 6; i++)
    prism_vtx[i] = tria_vtx[i];

  /* Orientation des faces */

  const PDM_real_t *_coords = vtx->_coords;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 3; j++) {
      int isom = prism_vtx[3*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += _coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 1.0/3.0;
    
    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;
    
    double v1[3];
    double v2[3];
    int isom3 = prism_vtx[3*i+2] - 1 ;
    int isom2 = prism_vtx[3*i+1] - 1;
    int isom1 = prism_vtx[3*i] - 1;

    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom2+k] - _coords[3*isom1+k]; 
      v2[k] = _coords[3*isom3+k] - _coords[3*isom1+k];
    } 
    _p_cross(v1, v2, n + 3*i);
  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  double orientation = _p_dot(cc, n);
  double orientation2 = _p_dot(cc, n+3);

  if (orientation < 0) {
    int tmp = prism_vtx[1];
    prism_vtx[1] = prism_vtx[2];
    prism_vtx[2] = tmp;
  } 

  if (orientation2 < 0) {
    int tmp = prism_vtx[4];
    prism_vtx[4] = prism_vtx[5];
    prism_vtx[5] = tmp;
  } 

  /* Permutation circulaire */

  int id1 = -1;
  for (int j = 0; j < 12; j++) {
    if (quad_vtx[j] == prism_vtx[0]) {
      id1 = j;
      break;
    }
  }

  int id2 = (id1 / 4) * 4 + (id1 + 1) % 4;
  if ((quad_vtx[id2] == prism_vtx[1]) ||
      (quad_vtx[id2] == prism_vtx[2]))
    id2 =  (id1 / 4) * 4 + (id1 + 3) % 4;

  int id_deb;
  for (int j = 0; j < 3; j++) {
    if (quad_vtx[id2] == prism_vtx[3+j]) {
      id_deb = j;
      break;
    }
  }

  int tmp[3];
  for (int j = 0; j < 3; j++)
    tmp[j] = prism_vtx[3+j];

  for (int j = 0; j < 3; j++) {
    int idx = (id_deb + j) % 3;
    prism_vtx[3+j] = tmp[idx];
  }

}


/**
 * 
 * \brief Build pyramid nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] pyramid_vtx Pyramid connectivity
 *
 */

static void
_connec_pyramid
(
PDM_Mesh_nodal_vtx_t *vtx,
PDM_l_num_t *tria_vtx,
PDM_l_num_t *quad_vtx,
PDM_l_num_t pyramid_vtx[]
)
{

  /* Initialisation */

  pyramid_vtx[0] = quad_vtx[0];
  pyramid_vtx[1] = quad_vtx[1];
  pyramid_vtx[2] = quad_vtx[2];
  pyramid_vtx[3] = quad_vtx[3];

  for (int i = 0; i < 9; i++) {
    if ((tria_vtx[i] != pyramid_vtx[0]) &&
        (tria_vtx[i] != pyramid_vtx[1]) &&
        (tria_vtx[i] != pyramid_vtx[2]) &&
        (tria_vtx[i] != pyramid_vtx[3])) {
      pyramid_vtx[4] = tria_vtx[i];
      break;
    }
  }

  /* Orientation */

  const PDM_real_t *_coords = vtx->_coords;

  double c[3];
  double n[3];

  for (int k = 0; k < 3; k++)
    c[k] = 0.;
  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    for (int k = 0; k < 3; k++)
      c[k] += _coords[3*isom+k];
  }
  for (int k = 0; k < 3; k++)
    c[k] *= 0.25;
    
  for (int k = 0; k < 3; k++)
    n[k] = 0.;
    
  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    int suiv = (j+1) % 4;
    int isom_suiv = pyramid_vtx[suiv] - 1;

    double v1[3];
    double v2[3];
    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom+k] -  c[k]; 
      v2[k] = _coords[3*isom_suiv+k] -  c[k]; 
    } 
      
    _p_cross(v1, v2, n);

  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = _coords[3*(pyramid_vtx[3] - 1) + k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_dot(cc, n);
  
  if (orientation < 0) {
    int tmp = pyramid_vtx[0];
    pyramid_vtx[0] = pyramid_vtx[3];
    pyramid_vtx[3] = tmp;
    tmp = pyramid_vtx[1];
    pyramid_vtx[1] = pyramid_vtx[2];
    pyramid_vtx[2] = tmp;
  } 

}


/**
 * 
 * \brief Build hexahedron nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] hexa_vtx    Hexahedron connectivity
 *
 */

static void
_connec_hexa
(
PDM_Mesh_nodal_vtx_t *vtx,
PDM_l_num_t *quad_vtx,
PDM_l_num_t hexa_vtx[]
)
{

  /* Initialization */

  hexa_vtx[0] = quad_vtx[0];
  hexa_vtx[1] = quad_vtx[1];
  hexa_vtx[2] = quad_vtx[2];
  hexa_vtx[3] = quad_vtx[3];

  PDM_l_num_t face_contact[4];

  for (int i = 1; i < 6; i++) {
    int cpt = 0;
    for (int j = 0; j < 4; j++) {
      PDM_l_num_t som_courant = quad_vtx[4*i+j];
      if ((som_courant != hexa_vtx[0]) &&
          (som_courant != hexa_vtx[1]) &&
          (som_courant != hexa_vtx[2]) &&
          (som_courant != hexa_vtx[3]))
        cpt += 1;
    }
    if (cpt == 4) {
      hexa_vtx[4] = quad_vtx[4*i];
      hexa_vtx[5] = quad_vtx[4*i+1];
      hexa_vtx[6] = quad_vtx[4*i+2];
      hexa_vtx[7] = quad_vtx[4*i+3];
    }
    if (cpt == 2) {
      face_contact[0] = quad_vtx[4*i];
      face_contact[1] = quad_vtx[4*i+1];
      face_contact[2] = quad_vtx[4*i+2];
      face_contact[3] = quad_vtx[4*i+3];
    }
  } 

  /* Calcul des centres et normales de la base et de la face opposee */

  const PDM_real_t *_coords = vtx->_coords;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += _coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 0.25;
    
    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;
    
    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      int suiv = (j+1) % 4;
      int isom_suiv = hexa_vtx[4*i+suiv] - 1;

      double v1[3];
      double v2[3];
      for (int k = 0; k < 3; k++) {
        v1[k] = _coords[3*isom+k] -  c[3*i+k]; 
        v2[k] = _coords[3*isom_suiv+k] -  c[3*i+k]; 
      } 
      
      _p_cross(v1, v2, n + 3*i);

    }

  }
  
  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_dot(cc, n);
  double orientation2 = _p_dot(cc, n+3);

  if (orientation < 0) {
    int tmp = hexa_vtx[0];
    hexa_vtx[0] = hexa_vtx[3];
    hexa_vtx[3] = tmp;
    tmp = hexa_vtx[1];
    hexa_vtx[1] = hexa_vtx[2];
    hexa_vtx[2] = tmp;
  } 

  if (orientation2 < 0) {
    int tmp = hexa_vtx[4];
    hexa_vtx[4] = hexa_vtx[7];
    hexa_vtx[7] = tmp;
    tmp = hexa_vtx[5];
    hexa_vtx[5] = hexa_vtx[6];
    hexa_vtx[6] = tmp;
  } 

  /* Permutation circulaire eventuelle de la face sup */

  int id1 = -1;
  int k1 = -1;
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      if (face_contact[j] == hexa_vtx[k]) {
        id1 = j;
        k1 = k;
        break;
      }
      if (id1 != -1)
        break;
    }
  }

  if (k1 == -1) {
    PDM_printf("Error connect_hexa : %d %d %d %d %d %d %d %d\n",
           hexa_vtx[0],
           hexa_vtx[1],
           hexa_vtx[2],
           hexa_vtx[3],
           hexa_vtx[4],
           hexa_vtx[5],
           hexa_vtx[6],
           hexa_vtx[7]);

    for (int i10 = 0; i10 < 4; i10++) {
      PDM_printf("   face %d : %d %d %d %d\n", i10+1, quad_vtx[4*i10],  
                                                  quad_vtx[4*i10+1],
                                                  quad_vtx[4*i10+2],
                                                  quad_vtx[4*i10+3]);
    }
    abort();
    
  }
    
  int id2 = (id1 + 1) % 4;
  int k2 = (k1 + 1) % 4;
  int k3 = (k1 + 3) % 4;
  
  if ((face_contact[id2] == hexa_vtx[k2]) ||
      (face_contact[id2] == hexa_vtx[k3]))
    id2 = (id1 + 3) % 4;

  int id_deb;
  for (int j = 0; j < 4; j++) {
    if (face_contact[id2] == hexa_vtx[4+j]) {
      id_deb = (j - k1);
      if (id_deb < 0) 
        id_deb += 4;
      id_deb = id_deb % 4;
      break;
    }
  }

  int tmp[4];
  for (int j = 0; j < 4; j++)
    tmp[j] = hexa_vtx[4+j];

  for (int j = 0; j < 4; j++) {
    int idx = (id_deb + j) % 4;
    hexa_vtx[4+j] = tmp[idx];
  }
}


/** 
 * \brief Get element type
 * 
 *   \param[in] n_face_cell     Number of faces in the current cell
 *   \param[in]  face_cell      Face to cell connectivity
 *   \param[in]  face_cell_idx  Face to vertex connectivity index 
 *   \param[in]  face_cell_n    Number of vertices for each face
 *   \param[in]  face_vtx       Face to vertex connectivity
 *   \param[out] tria_vtx       Quadrangles connectivity
 *   \param[out] quad_vtx       Quadrangles connectivity
 *
 *  \return   Cell type
 */

inline static
PDM_Mesh_nodal_elt_t
_type_cell_3D
(
 const int             n_face_cell,
 const PDM_l_num_t    *cell_face,
 const PDM_l_num_t    *face_vtx_idx,
 const PDM_l_num_t    *face_vtx_nb,
 const PDM_l_num_t    *face_vtx,
 PDM_l_num_t           tria_vtx[],
 PDM_l_num_t           quad_vtx[]
)
{
  int adjust = 0;
  if (face_vtx_idx[0] == 1) {
    adjust = 1;
  }

  
  PDM_l_num_t  n_trias = 0;
  PDM_l_num_t  n_quads = 0;

  if (n_face_cell > 6) {
    return PDM_MESH_NODAL_POLY_3D;
  }

  for (int i = 0; i < n_face_cell; i++) {

    const int face_id = PDM_ABS(cell_face[i]) - 1;
    const int n_som_face = face_vtx_nb[face_id];
    PDM_l_num_t idx = face_vtx_idx[face_id] - adjust;
 
    if (n_som_face == 3) {
      PDM_l_num_t *cell_som_tria_courant = tria_vtx + 3*n_trias;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_tria_courant[j-idx] = face_vtx[j];
      }
      n_trias += 1;
    }
    else if (n_som_face == 4) {
      PDM_l_num_t *cell_som_quad_courant = quad_vtx + 4*n_quads;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_quad_courant[j-idx] = face_vtx[j];
      }
      n_quads += 1;
    }
    else 
      return PDM_MESH_NODAL_POLY_3D;

  }
    
  PDM_Mesh_nodal_elt_t cell_type;

  if ((n_quads == 0) && (n_trias == 4))
    cell_type = PDM_MESH_NODAL_TETRA4;
  else if (n_quads == 6)
    cell_type = PDM_MESH_NODAL_HEXA8;
  else if ((n_quads == 1) && (n_trias == 4))
    cell_type = PDM_MESH_NODAL_PYRAMID5;
  else if ((n_quads == 3) && (n_trias == 2)) {
    int trias[6];
    n_trias = 0;
    for (int i = 0; i < n_face_cell; i++) {

      const int face_id = PDM_ABS(cell_face[i]) - 1;
      const int ideb = face_vtx_idx[face_id] - adjust;

      const int n_som_face = face_vtx_nb[face_id];
 
      if (n_som_face == 3) {
        for (int j = 0; j < 3; j++) {
          trias[3*n_trias+j] = face_vtx[ideb+j];
        }
        n_trias += 1;
      }
      if (n_trias >= 2)
        break;
    }

    cell_type = PDM_MESH_NODAL_PRISM6;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (trias[i] == trias[3+j]) {
          cell_type = PDM_MESH_NODAL_POLY_3D;
          break;
        }
      }
      if (cell_type == PDM_MESH_NODAL_POLY_3D)
        break;
    }
  }

  else {
    cell_type = PDM_MESH_NODAL_POLY_3D;
  }

  return cell_type;

}


/*=============================================================================
 * Public function definitions
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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) malloc (sizeof(PDM_Mesh_nodal_t));
  
  _mesh_init (mesh, n_part, comm);
  
  if (mesh_handles == NULL) {
    mesh_handles = PDM_Handles_create (4);
  }
  
  return PDM_Handles_store (mesh_handles, mesh);
}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->n_part;

}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
	
  if (mesh->blocks_std != NULL) {
    const int n_blocks_std = PDM_Handles_n_get (mesh->blocks_std);
    const int *list_ind = PDM_Handles_idx_get (mesh->blocks_std);

    for (int i = 0; i < n_blocks_std; i++) {
      PDM_Mesh_nodal_block_std_t *_block_std = 
              (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, list_ind[i]);
      _block_std_free_partial(_block_std);
    }
  }
	
  if (mesh->blocks_poly2d != NULL) {
    const int n_blocks_poly2d = PDM_Handles_n_get (mesh->blocks_poly2d);
    const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly2d);

    for (int i = 0; i < n_blocks_poly2d; i++) {
      PDM_Mesh_nodal_block_poly2d_t *_block_poly2d = 
              (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, list_ind[i]);
      _block_poly2d_free_partial(_block_poly2d);
    }
  }
	
  if (mesh->blocks_poly3d != NULL) {
    const int n_blocks_poly3d = PDM_Handles_n_get (mesh->blocks_poly3d);
    const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly3d);

    for (int i = 0; i < n_blocks_poly3d; i++) {
      PDM_Mesh_nodal_block_poly3d_t *_block_poly3d = 
              (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, list_ind[i]);
      _block_poly3d_free_partial(_block_poly3d);
    }
  }
}


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
)
{
  
  PDM_Mesh_nodal_partial_free (idx);
  
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh != NULL) {

    if (mesh->blocks_id != NULL) {
      free (mesh->blocks_id);
    }
    
    mesh->blocks_id = NULL;
    
    /* Free vertices */
  
    if (mesh->vtx != NULL) {
      for (int i = 0; i < mesh->n_part; i++) {
        mesh->vtx[i] = _vtx_free (mesh->vtx[i]);
      }

      free(mesh->vtx);
      mesh->vtx = NULL;
    }

    /* free standard blocks */

    if (mesh->blocks_std != NULL) {
      int n_blocks_std = PDM_Handles_n_get (mesh->blocks_std);
      const int *list_ind = PDM_Handles_idx_get (mesh->blocks_std);

      while (n_blocks_std > 0) {
        PDM_Mesh_nodal_block_std_t *_bloc_std = 
          (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, list_ind[0]);
        _block_std_free(_bloc_std);
        PDM_Handles_handle_free (mesh->blocks_std, list_ind[0], PDM_FALSE);
        n_blocks_std = PDM_Handles_n_get (mesh->blocks_std);
      }

      mesh->blocks_std = PDM_Handles_free (mesh->blocks_std); 
    }

    /* Free polygon blocks */ 

    if (mesh->blocks_poly2d != NULL) {
      int n_blocks_poly2d = PDM_Handles_n_get (mesh->blocks_poly2d);
      const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly2d);

      while (n_blocks_poly2d > 0) {
        PDM_Mesh_nodal_block_poly2d_t *_bloc_poly2d = 
          (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, list_ind[0]);
        _block_poly2d_free(_bloc_poly2d);
        PDM_Handles_handle_free (mesh->blocks_poly2d, list_ind[0], PDM_FALSE);
        n_blocks_poly2d = PDM_Handles_n_get (mesh->blocks_poly2d);
      }

      mesh->blocks_poly2d = PDM_Handles_free (mesh->blocks_poly2d); 
    }

    /* Free polyhedron blocks */ 

    if (mesh->blocks_poly3d != NULL) {
      int n_blocks_poly3d = PDM_Handles_n_get (mesh->blocks_poly3d);
      const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly3d);

      while (n_blocks_poly3d > 0) {
        PDM_Mesh_nodal_block_poly3d_t *_bloc_poly3d = 
          (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, list_ind[0]);
        _block_poly3d_free(_bloc_poly3d);
        PDM_Handles_handle_free (mesh->blocks_poly3d, list_ind[0], PDM_FALSE);
        n_blocks_poly3d = PDM_Handles_n_get (mesh->blocks_poly3d);
      }

      mesh->blocks_poly3d = PDM_Handles_free (mesh->blocks_poly3d); 
    }

    /* Free structure */ 

    if (mesh->num_cell_parent_to_local != NULL) {
      for (int ipart = 0; ipart < mesh->n_part; ipart++) {
        if (mesh->num_cell_parent_to_local[ipart] != NULL)
          free(mesh->num_cell_parent_to_local[ipart]);
      }
      free(mesh->num_cell_parent_to_local);
      mesh->num_cell_parent_to_local = NULL;
    }

    free(mesh->n_cell);
    mesh->n_cell = NULL;
    
    if (mesh->blocks_id != NULL) {
      free (mesh->blocks_id);
    }

    free(mesh);

    PDM_Handles_handle_free (mesh_handles, idx, PDM_FALSE);
  
    int n_mesh_array = PDM_Handles_n_get (mesh_handles); 

    if (n_mesh_array == 0) {
      mesh_handles = PDM_Handles_free (mesh_handles);
    }
  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  idx       Nodal mesh handle
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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "these partition vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;
  vtx->_numabs = numabs;

}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier %d %d\n", id_part, mesh->n_part);  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->n_vtx;
}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->_numparent;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_Mesh_nodal_vertices_get
(
 const int          idx,
 const int          id_part 
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->_coords;
}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->_numabs;
}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  return mesh->is_vtx_def_from_parent;

}

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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];
  
  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  vtx->parent = (PDM_Mesh_nodal_vtx_t *) malloc (sizeof (PDM_Mesh_nodal_vtx_t));
  PDM_Mesh_nodal_vtx_t *_parent = vtx->parent;
  _parent->parent = NULL;
  _parent->n_vtx = n_vtx_parent;
  _parent->coords = NULL;
  _parent->_coords = coords_parent ;
  _parent->_numabs = numabs_parent;
  _parent->_numparent = NULL;
  
  vtx->n_vtx      = n_vtx;
  vtx->coords     = malloc (sizeof(double) * 3 * n_vtx);
  vtx->_coords    = vtx->coords;
  vtx->_numabs    = numabs;
  vtx->_numparent = num_parent;

  for (int i = 0; i < n_vtx; i++) {
    int i_parent = num_parent[i] - 1;
    for (int j = 0; j < 3; j++) {
      vtx->coords[3*i+j] = _parent->_coords[3*i_parent+j];
    }
  }
  mesh->is_vtx_def_from_parent   = 1;

}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  
  return mesh->n_blocks;

}



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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  return mesh->blocks_id;

}


/**
 * \brief  Return type of block
 *
 * \param [in]  idx            Nodal mesh handle
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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  PDM_Mesh_nodal_elt_t t_elt;
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    t_elt = PDM_MESH_NODAL_POLY_3D;
    PDM_Mesh_nodal_block_std_t *block = 
            (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, id_block);

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad block identifier\n");  
    }

    t_elt = block->t_elt;
  }
  
  else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {
    
    t_elt = PDM_MESH_NODAL_POLY_2D;

  }
  
  else {
    
    t_elt = PDM_MESH_NODAL_POLY_3D;

  }
  
  return t_elt;
    
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int id_block;
  
  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :    
  case PDM_MESH_NODAL_BAR2     :   
  case PDM_MESH_NODAL_TRIA3    :    
  case PDM_MESH_NODAL_QUAD4    :    
  case PDM_MESH_NODAL_TETRA4   :     
  case PDM_MESH_NODAL_PYRAMID5 :     
  case PDM_MESH_NODAL_PRISM6   :     
  case PDM_MESH_NODAL_HEXA8    : 
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->blocks_std == NULL) {
        mesh->blocks_std = PDM_Handles_create (4);
      } 
      
      /* Allocation du bloc */
      
      PDM_Mesh_nodal_block_std_t *block_std = (PDM_Mesh_nodal_block_std_t *) malloc(sizeof(PDM_Mesh_nodal_block_std_t));

      id_block = PDM_Handles_store (mesh->blocks_std, block_std);

      /* Intialisation du bloc */

      block_std->t_elt = t_elt;
      block_std->st_free_data = st_free_data;
      block_std->n_part = mesh->n_part; 

      block_std->n_elt      = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t ) * block_std->n_part);
      block_std->_connec    = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * block_std->n_part);
      block_std->_num_part  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * block_std->n_part);
      block_std->_numabs    = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * block_std->n_part);
      block_std->numabs_int = NULL;
      block_std->_parent_num = NULL;

      for (int i = 0; i < block_std->n_part; i++) {
        block_std->n_elt[i]     = 0;
        block_std->_connec[i]   = NULL;
        block_std->_num_part[i] = NULL;
        block_std->_numabs[i]   = NULL;
      }
      
      id_block += PDM_BLOCK_ID_BLOCK_STD;
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard blocks must be less than %d\n", 
               PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }
    
    break;

  case PDM_MESH_NODAL_POLY_2D  :    
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->blocks_poly2d == NULL) {
        mesh->blocks_poly2d = PDM_Handles_create (4);
      }
      
      /* Allocation du bloc */
      
      PDM_Mesh_nodal_block_poly2d_t *block_poly2d =
              (PDM_Mesh_nodal_block_poly2d_t *) malloc(sizeof(PDM_Mesh_nodal_block_poly2d_t));

      id_block = PDM_Handles_store (mesh->blocks_poly2d, block_poly2d);
      
      /* Intialisation du bloc */

      block_poly2d->st_free_data = st_free_data;
      block_poly2d->n_part            = mesh->n_part; 

      block_poly2d->n_elt       = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t ) * block_poly2d->n_part);
      block_poly2d->_connec_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * block_poly2d->n_part);
      block_poly2d->_connec     = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * block_poly2d->n_part);
      block_poly2d->_num_part   = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * block_poly2d->n_part);
      block_poly2d->_numabs     = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * block_poly2d->n_part);
      block_poly2d->numabs_int = NULL;
      block_poly2d->_parent_num = NULL;

      for (int i = 0; i < block_poly2d->n_part; i++) {
        block_poly2d->n_elt[i]     = 0;
        block_poly2d->_connec_idx[i] = NULL;
        block_poly2d->_connec[i]     = NULL;
        block_poly2d->_num_part[i]= NULL;
        block_poly2d->_numabs[i]     = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_POLY2D;
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of polygon blocks must be less than %d\n",
               PDM_BLOCK_ID_BLOCK_POLY3D - PDM_BLOCK_ID_BLOCK_POLY2D);
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :     
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->blocks_poly3d == NULL) {
        mesh->blocks_poly3d = PDM_Handles_create (4);
      }
      
      /* Allocation du bloc */
      
      PDM_Mesh_nodal_block_poly3d_t *block_poly3d =
              (PDM_Mesh_nodal_block_poly3d_t *) malloc(sizeof(PDM_Mesh_nodal_block_poly3d_t));

      id_block = PDM_Handles_store (mesh->blocks_poly3d, block_poly3d);

      /* Intialisation du bloc */

      block_poly3d->n_part            = mesh->n_part; 
      block_poly3d->st_free_data = st_free_data;

      block_poly3d->n_elt        = (PDM_l_num_t *)   malloc(sizeof(PDM_l_num_t ) * block_poly3d->n_part);
      block_poly3d->n_face       = (PDM_l_num_t *)   malloc(sizeof(PDM_l_num_t ) * block_poly3d->n_part);
      block_poly3d->_facvtx_idx  = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * block_poly3d->n_part);
      block_poly3d->_facvtx      = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * block_poly3d->n_part);
      block_poly3d->_cellfac_idx = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * block_poly3d->n_part);
      block_poly3d->_cellfac     = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * block_poly3d->n_part);
      block_poly3d->_numabs      = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * block_poly3d->n_part);
      block_poly3d->numabs_int = NULL;
      block_poly3d->_parent_num = NULL;

      for (int i = 0; i < block_poly3d->n_part; i++) {
        block_poly3d->n_elt[i]        = 0;
        block_poly3d->n_face[i]       = 0;
        block_poly3d->_facvtx_idx[i]  = NULL;
        block_poly3d->_facvtx[i]      = NULL;
        block_poly3d->_cellfac_idx[i] = NULL;
        block_poly3d->_cellfac[i]     = NULL;
        block_poly3d->_numabs[i]      = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");

  }

  _update_blocks_id (mesh);
  return id_block ;
  
}


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
const int            idx,
const int            id_block,     
const int            id_part, 
const int            n_elt,    
      PDM_l_num_t   *connec,   
      PDM_g_num_t   *numabs,
      PDM_l_num_t   *parent_num   
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;
  
  PDM_Mesh_nodal_block_std_t *block = (PDM_Mesh_nodal_block_std_t *) 
     PDM_Handles_get (mesh->blocks_std, _id_block);
  
  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }
  
  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }
 
  /* Mapping */
  
  mesh->n_cell[id_part] += -block->n_elt[id_part]; 
  mesh->n_cell[id_part] += n_elt;
  block->n_elt[id_part] = n_elt;
  block->_connec[id_part] = connec;
  block->_numabs[id_part] = numabs;

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = parent_num;
  }
  
  for (int i = 0; i < n_elt; i++) {
    mesh->n_elt_abs = PDM_MAX (mesh->n_elt_abs, numabs[i]);
  }
 
}


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
 * \param [out]  n_elt          Number of elements
 * \param [out]  connect        Connectivity
 * \param [out]  numabs         Global numbering
 * \param [out] numabs_block   Global numbering in the block or NULL (if not computed)
 * \param [out] parent_num     Parent numbering or NULL
*
 */

void
PDM_Mesh_nodal_block_std_get 
(   
const int            idx,
const int            id_block,     
const int            id_part, 
      PDM_l_num_t  **connec   
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;
  
  const PDM_Mesh_nodal_block_std_t *block = (const PDM_Mesh_nodal_block_std_t *) 
     PDM_Handles_get (mesh->blocks_std, _id_block);
  
  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }
  
  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *connec = block->_connec[id_part];
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_block;
  
  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_Mesh_nodal_block_poly3d_t *block = (const PDM_Mesh_nodal_block_poly3d_t *) 
     PDM_Handles_get (mesh->blocks_poly3d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }
    
    return block->n_elt[id_part];
  }
  
  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_Mesh_nodal_block_poly2d_t *block = (const PDM_Mesh_nodal_block_poly2d_t *) 
     PDM_Handles_get (mesh->blocks_poly2d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->n_elt[id_part];
  }
  
  else {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_Mesh_nodal_block_std_t *block = (const PDM_Mesh_nodal_block_std_t *) 
     PDM_Handles_get (mesh->blocks_std, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->n_elt[id_part];
  }
  
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_block;
  
  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_Mesh_nodal_block_poly3d_t *block = (const PDM_Mesh_nodal_block_poly3d_t *) 
     PDM_Handles_get (mesh->blocks_poly3d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }
    
    return block->numabs_int[id_part];
  }
  
  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_Mesh_nodal_block_poly2d_t *block = (const PDM_Mesh_nodal_block_poly2d_t *) 
     PDM_Handles_get (mesh->blocks_poly2d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->numabs_int[id_part];
  }
  
  else {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_Mesh_nodal_block_std_t *block = (const PDM_Mesh_nodal_block_std_t *) 
     PDM_Handles_get (mesh->blocks_std, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->numabs_int[id_part];
  }
}


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
) 
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_block;
  
  PDM_g_num_t *_numabs_int = NULL;
  
  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_Mesh_nodal_block_poly3d_t *block = (const PDM_Mesh_nodal_block_poly3d_t *) 
     PDM_Handles_get (mesh->blocks_poly3d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }
    
    if (block->numabs_int != NULL) {      
      _numabs_int = block->numabs_int[id_part];
    }
  }
  
  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_Mesh_nodal_block_poly2d_t *block = (const PDM_Mesh_nodal_block_poly2d_t *) 
     PDM_Handles_get (mesh->blocks_poly2d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->numabs_int != NULL) {      
      _numabs_int = block->numabs_int[id_part];
    }
  }
  
  else {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_Mesh_nodal_block_std_t *block = (const PDM_Mesh_nodal_block_std_t *) 
     PDM_Handles_get (mesh->blocks_std, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->numabs_int != NULL) {      
      _numabs_int = block->numabs_int[id_part];
    }
  }
  
  return _numabs_int;
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_block;
  
  int *_parent_num = NULL;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_Mesh_nodal_block_poly3d_t *block = (const PDM_Mesh_nodal_block_poly3d_t *) 
     PDM_Handles_get (mesh->blocks_poly3d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }
    
    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }
  
  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_Mesh_nodal_block_poly2d_t *block = (const PDM_Mesh_nodal_block_poly2d_t *) 
     PDM_Handles_get (mesh->blocks_poly2d, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }
  
  else {
  
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_Mesh_nodal_block_std_t *block = (const PDM_Mesh_nodal_block_std_t *) 
     PDM_Handles_get (mesh->blocks_std, _id_block);
  
    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }
  
    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }
  
  return _parent_num;
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_Mesh_nodal_block_poly2d_t *block = 
          (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, _id_block);

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }
  
  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  /* Mapping */

  mesh->n_cell[id_part]      += -block->n_elt[id_part]; 
  mesh->n_cell[id_part]      += n_elt;
  block->n_elt[id_part]       = n_elt;
  block->_connec_idx[id_part] = connec_idx;
  block->_connec[id_part]     = connec;
  block->_numabs[id_part]     = numabs;

  for (int i = 0; i < n_elt; i++) {
    mesh->n_elt_abs = PDM_MAX(mesh->n_elt_abs, numabs[i]);
  }

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = parent_num;
  }
  
}



/**
 * \brief Return a polygon block description
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_Mesh_nodal_block_poly2d_t *block = 
          (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, _id_block);

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }
  
  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *connec_idx = block->_connec_idx[id_part];
  *connec     = block->_connec[id_part];
 
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
  PDM_Mesh_nodal_block_poly3d_t *block = 
          (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, _id_block);

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }
  
  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  mesh->n_cell[id_part]       += -block->n_elt[id_part]; 
  mesh->n_cell[id_part]       += n_elt;
  block->n_elt[id_part]        = n_elt;
  block->n_face[id_part]       = n_face;
  block->_facvtx_idx[id_part]  = facvtx_idx;
  block->_facvtx[id_part]      = facvtx;
  block->_cellfac_idx[id_part] = cellfac_idx;
  block->_cellfac[id_part]     = cellfac;
  block->_numabs[id_part]      = numabs;

  for (int i = 0; i < n_elt; i++) {
    mesh->n_elt_abs = PDM_MAX (mesh->n_elt_abs, numabs[i]);
  }

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = parent_num;
  }
  
}


/**
 * \brief Define a polyhedra block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
  PDM_Mesh_nodal_block_poly3d_t *block = 
          (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, _id_block);

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }
  
  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *n_face     = block->n_face[id_part];
  *facvtx_idx  = block->_facvtx_idx[id_part];
  *facvtx      = block->_facvtx[id_part];
  *cellfac_idx = block->_cellfac_idx[id_part];
  *cellfac     = block->_cellfac[id_part];

}


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
 * \param [out] ind_num        old numbering to new numbering
 *
 */

void
PDM_Mesh_nodal_cell3d_cellface_add
(
const int         idx,
const int         id_part, 
const int         n_cell,
const int         n_face,
PDM_l_num_t      *face_vtx_idx,
PDM_l_num_t      *face_vtx_nb,
PDM_l_num_t      *face_vtx,
PDM_l_num_t      *cell_face_idx,
PDM_l_num_t      *cell_face_nb,
PDM_l_num_t      *cell_face,
PDM_g_num_t      *numabs
)
{
  int adjust = 0;
  if (n_cell > 0) { 
    if (cell_face_idx[0] == 1) {
      adjust = 1;
    }
  }
  
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
    
  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  int n_part = 0;

  if (mesh->num_cell_parent_to_local == NULL) {
    mesh->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part); 
    for (int ipart = 0; ipart < mesh->n_part; ipart++) { 
      mesh->num_cell_parent_to_local[ipart] = NULL;
    }
  }

  mesh->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    mesh->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (mesh->prepa_blocks == NULL) {
    mesh->prepa_blocks = (PDM_Mesh_nodal_prepa_blocks_t *) malloc(sizeof(PDM_Mesh_nodal_prepa_blocks_t));
    mesh->prepa_blocks->t_add = 1;
    mesh->prepa_blocks->n_tria_proc = 0;    /* Nb de triangles par proc */
    mesh->prepa_blocks->n_quad_proc = 0;    /* Nb de quads par proc */
    mesh->prepa_blocks->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    mesh->prepa_blocks->n_tetra_proc = 0;   /* Nb de tetra par proc */
    mesh->prepa_blocks->n_hexa_proc = 0;    /* Nb d'hexa par proc */
    mesh->prepa_blocks->n_prism_proc = 0;   /* Nb de prisme par proc */
    mesh->prepa_blocks->n_pyramid_proc = 0; /* Nb de pyramide par proc */
    mesh->prepa_blocks->n_poly3d_proc = 0;  /* Nb de poly3d par proc */
    mesh->prepa_blocks->n_cell = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->face_vtx_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part); 
    mesh->prepa_blocks->face_vtx_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*mesh->n_part);
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (mesh->prepa_blocks->t_add != 1) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur Cs_geom_cell3d_cellface_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  /* Determination du type de chaque element */

  PDM_l_num_t cell_som_tria[18]; /* 6 triangles max in _type_cell_3D */
  PDM_l_num_t cell_som_quad[24]; /* 6 quadrangles max in _type_cell_3D */
  PDM_l_num_t n_tetra   = 0;
  PDM_l_num_t n_hexa    = 0;
  PDM_l_num_t n_prism   = 0;
  PDM_l_num_t n_pyramid = 0;
  PDM_l_num_t n_poly3d  = 0;
  
  if (1 == 0) {
    printf("cellface 2: \n");
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        printf(" %d", cell_face[j]);
      }
    printf("\n");
    }
      
    printf("facevtx 2: \n");
    for (int i = 0; i < n_face; i++) {
      for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
        printf(" %d", face_vtx[j]);
      }
    printf("\n");
    }
  }
    
  for (int i = 0; i < n_cell; i++) {
    
    PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_nb[i],
                                                    cell_face + cell_face_idx[i] - adjust,
                                                    face_vtx_idx,
                                                    face_vtx_nb,
                                                    face_vtx,
                                                    cell_som_tria,
                                                    cell_som_quad);
    switch(cell_type) {
    case PDM_MESH_NODAL_TETRA4 :
      n_tetra += 1;
      break;
    case PDM_MESH_NODAL_PYRAMID5 :
      n_pyramid += 1;
      break;
    case PDM_MESH_NODAL_PRISM6 :
      n_prism += 1;
      break;
    case PDM_MESH_NODAL_HEXA8 :
      n_hexa += 1;
      break;
    case PDM_MESH_NODAL_POLY_3D :
      n_poly3d += 1;
      break;
    default :
      break;
    }
  }

  mesh->prepa_blocks->n_tetra_proc          += n_tetra;
  mesh->prepa_blocks->n_hexa_proc           += n_hexa;
  mesh->prepa_blocks->n_prism_proc          += n_prism;
  mesh->prepa_blocks->n_pyramid_proc        += n_pyramid;
  mesh->prepa_blocks->n_poly3d_proc         += n_poly3d;
  mesh->prepa_blocks->n_tetra[id_part]       = n_tetra;
  mesh->prepa_blocks->n_hexa[id_part]        = n_hexa;
  mesh->prepa_blocks->n_prism[id_part]       = n_prism;
  mesh->prepa_blocks->n_pyramid[id_part]     = n_pyramid;
  mesh->prepa_blocks->n_poly3d[id_part]      = n_poly3d;
  mesh->prepa_blocks->face_vtx_idx[id_part]  = face_vtx_idx;
  mesh->prepa_blocks->face_vtx_nb[id_part]   = face_vtx_nb;
  mesh->prepa_blocks->face_vtx[id_part]      = face_vtx;
  mesh->prepa_blocks->cell_face_idx[id_part] = cell_face_idx;
  mesh->prepa_blocks->cell_face_nb[id_part]  = cell_face_nb;
  mesh->prepa_blocks->cell_face[id_part]     = cell_face;
  mesh->prepa_blocks->numabs[id_part]        = numabs;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_face[id_part]        = n_face;
  mesh->prepa_blocks->n_cell[id_part]          = n_cell;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < mesh->n_part; i++) {
    if (mesh->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (mesh->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[5];
    PDM_l_num_t som_elts[5];

    elts[0] = mesh->prepa_blocks->n_tetra_proc > 0;
    elts[1] = mesh->prepa_blocks->n_hexa_proc  > 0;
    elts[2] = mesh->prepa_blocks->n_prism_proc  > 0;
    elts[3] = mesh->prepa_blocks->n_pyramid_proc  > 0;
    elts[4] = mesh->prepa_blocks->n_poly3d_proc  > 0;
    
    PDM_MPI_Allreduce(elts, som_elts, 5, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);

    int id_bloc_tetra4;
    int id_bloc_hexa8;
    int id_bloc_prism6;
    int id_bloc_pyramid5;
    int id_bloc_poly_3d; 

    if (som_elts[0] > 0)
      id_bloc_tetra4 = PDM_Mesh_nodal_block_add(idx,
                                        PDM_TRUE,
                                        PDM_MESH_NODAL_TETRA4);

    if (som_elts[1] > 0)
      id_bloc_hexa8 = PDM_Mesh_nodal_block_add(idx,
                                       PDM_TRUE,
                                       PDM_MESH_NODAL_HEXA8);
    
    if (som_elts[2] > 0)
      id_bloc_prism6 = PDM_Mesh_nodal_block_add(idx,
                                        PDM_TRUE,
                                        PDM_MESH_NODAL_PRISM6);

    if (som_elts[3] > 0)
      id_bloc_pyramid5 = PDM_Mesh_nodal_block_add(idx,
                                          PDM_TRUE,
                                          PDM_MESH_NODAL_PYRAMID5);

    if (som_elts[4] > 0)
      id_bloc_poly_3d = PDM_Mesh_nodal_block_add(idx,
                                         PDM_TRUE,
                                         PDM_MESH_NODAL_POLY_3D);
                                                   
    /* Determination de la connectivite de chaque element */
    

    for (int ipart = 0; ipart < mesh->n_part; ipart++) {
      
      PDM_l_num_t n_cell_courant = mesh->prepa_blocks->n_cell[ipart];
      PDM_l_num_t *num_cell_parent_to_local_courant = mesh->num_cell_parent_to_local[ipart];
      PDM_l_num_t *face_som_idx_courant = mesh->prepa_blocks->face_vtx_idx[ipart];
      PDM_l_num_t *face_som_nb_courant = mesh->prepa_blocks->face_vtx_nb[ipart];
      PDM_l_num_t *face_som_courant = mesh->prepa_blocks->face_vtx[ipart];
      PDM_l_num_t *cell_face_idx_courant = mesh->prepa_blocks->cell_face_idx[ipart];
      PDM_l_num_t *cell_face_nb_courant = mesh->prepa_blocks->cell_face_nb[ipart];
      PDM_l_num_t *cell_face_courant = mesh->prepa_blocks->cell_face[ipart];
      PDM_g_num_t *numabs_courant = mesh->prepa_blocks->numabs[ipart];
  
      PDM_l_num_t n_face_part   = mesh->prepa_blocks->n_face[ipart];
 
      PDM_l_num_t n_tetra_part   = mesh->prepa_blocks->n_tetra[ipart];
      PDM_l_num_t n_hexa_part    = mesh->prepa_blocks->n_hexa[ipart];
      PDM_l_num_t n_prism_part   = mesh->prepa_blocks->n_prism[ipart];
      PDM_l_num_t n_pyramid_part = mesh->prepa_blocks->n_pyramid[ipart];
      PDM_l_num_t n_poly3d_part  = mesh->prepa_blocks->n_poly3d[ipart];
      
      PDM_l_num_t *connec_tetra = NULL;
      PDM_l_num_t *connec_hexa = NULL;
      PDM_l_num_t *connec_prism = NULL;
      PDM_l_num_t *connec_pyramid = NULL;
 
      PDM_g_num_t *numabs_tetra = NULL;
      PDM_g_num_t *numabs_hexa = NULL;
      PDM_g_num_t *numabs_prism = NULL;
      PDM_g_num_t *numabs_pyramid = NULL;
      PDM_g_num_t *numabs_poly3d = NULL;
    
      PDM_l_num_t *num_parent_tetra = NULL;
      PDM_l_num_t *num_parent_hexa = NULL;
      PDM_l_num_t *num_parent_prism = NULL;
      PDM_l_num_t *num_parent_pyramid = NULL;
      PDM_l_num_t *num_parent_poly3d = NULL;

      if (n_tetra_part > 0) {
        connec_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 *n_tetra_part);
        numabs_tetra = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tetra_part);
        num_parent_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_tetra_part);
      }

      if (n_hexa_part > 0) {
        connec_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 8 * n_hexa_part);
        numabs_hexa = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_hexa_part);
        num_parent_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_hexa_part);
      }

      if (n_prism_part > 0) {
        connec_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 6 * n_prism_part);
        numabs_prism = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_prism_part);
        num_parent_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_prism_part);
      }

      if (n_pyramid_part > 0) {
        connec_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 5 * n_pyramid_part);
        numabs_pyramid = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_pyramid_part);
        num_parent_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_pyramid_part);
      }

      if (n_poly3d_part > 0) {
        numabs_poly3d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly3d_part);
        num_parent_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_poly3d_part);
      }

      PDM_l_num_t *num_parent_tetra_courant = num_parent_tetra;
      PDM_l_num_t *num_parent_hexa_courant = num_parent_hexa;
      PDM_l_num_t *num_parent_prism_courant = num_parent_prism;
      PDM_l_num_t *num_parent_pyramid_courant = num_parent_pyramid;
      PDM_l_num_t *num_parent_poly3d_courant = num_parent_poly3d;

      PDM_l_num_t *connec_tetra_courant = connec_tetra;
      PDM_l_num_t *connec_hexa_courant = connec_hexa;
      PDM_l_num_t *connec_prism_courant = connec_prism;
      PDM_l_num_t *connec_pyramid_courant = connec_pyramid;

      PDM_g_num_t *numabs_tetra_courant = numabs_tetra;
      PDM_g_num_t *numabs_hexa_courant = numabs_hexa;
      PDM_g_num_t *numabs_prism_courant = numabs_prism;
      PDM_g_num_t *numabs_pyramid_courant = numabs_pyramid;
      PDM_g_num_t *numabs_poly3d_courant = numabs_poly3d;

      PDM_l_num_t *tag_face_poly3d = NULL;
      PDM_l_num_t  n_face_poly = 0;
      PDM_l_num_t *facsom_poly_idx = NULL;
      PDM_l_num_t *facsom_poly = NULL;
      PDM_l_num_t *cellfac_poly_idx = NULL;
      PDM_l_num_t *cellfac_poly = NULL;
      PDM_l_num_t l_cellfac_poly = 0;

      if (n_poly3d_part > 0) {
        tag_face_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face_part);
        for (int i = 0; i < n_face_part; i++) {
          tag_face_poly3d[i] = -1;
        }
        cellfac_poly_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly3d_part + 1));
        cellfac_poly_idx[0] = 0;
      }

      PDM_l_num_t idx_tetra = 0;
      PDM_l_num_t idx_hexa = n_tetra_part;
      PDM_l_num_t idx_prism = idx_hexa + n_hexa_part;
      PDM_l_num_t idx_pyramid = idx_prism + n_prism_part;
      PDM_l_num_t idx_poly3d = idx_pyramid + n_pyramid_part;

      n_poly3d_part = 0; 
      for (int i = 0; i < n_cell_courant; i++) {
        num_cell_parent_to_local_courant[i] = 0;
        PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_nb_courant[i],
                                                cell_face_courant + cell_face_idx_courant[i] - adjust,
                                                face_som_idx_courant,
                                                face_som_nb_courant,
                                                face_som_courant,
                                                cell_som_tria,
                                                cell_som_quad);
        
        switch(cell_type) {
        case PDM_MESH_NODAL_TETRA4 :
          _connec_tetra(mesh->vtx[ipart],
                        cell_som_tria,
                        connec_tetra_courant);
          *numabs_tetra_courant = numabs_courant[i];
          numabs_tetra_courant += 1;
          connec_tetra_courant += 4;
          *num_parent_tetra_courant = i;
          num_parent_tetra_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tetra++;
          break;
        case PDM_MESH_NODAL_HEXA8 :
          _connec_hexa(mesh->vtx[ipart],
                       cell_som_quad,
                       connec_hexa_courant);
          *numabs_hexa_courant = numabs_courant[i];
          numabs_hexa_courant += 1;
          connec_hexa_courant += 8;
          *num_parent_hexa_courant = i;
          num_parent_hexa_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_hexa++;
          break;
        case PDM_MESH_NODAL_PRISM6 :
          _connec_prism(mesh->vtx[ipart],
                        cell_som_tria,
                        cell_som_quad,
                        connec_prism_courant);
          *numabs_prism_courant = numabs_courant[i];
          numabs_prism_courant += 1;
          connec_prism_courant += 6;
          *num_parent_prism_courant = i;
          num_parent_prism_courant += 1;          
          num_cell_parent_to_local_courant[i] = idx_prism++;
          break;
        case PDM_MESH_NODAL_PYRAMID5 :
          _connec_pyramid(mesh->vtx[ipart],
                          cell_som_tria,
                          cell_som_quad,
                          connec_pyramid_courant);
          *numabs_pyramid_courant = numabs_courant[i];
          numabs_pyramid_courant += 1;
          connec_pyramid_courant += 5;
          *num_parent_pyramid_courant = i;
          num_parent_pyramid_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_pyramid++;
          break;
        case PDM_MESH_NODAL_POLY_3D : 
          {
            PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - adjust;
            for (int j = 0; j < cell_face_nb_courant[i]; j++) {
              tag_face_poly3d[cell_face_cell[j] - 1] = 0;
            }
            *numabs_poly3d_courant = numabs_courant[i];
            numabs_poly3d_courant += 1;
            l_cellfac_poly += cell_face_nb_courant[i];
            cellfac_poly_idx[n_poly3d_part+1] = l_cellfac_poly;
            n_poly3d_part += 1;
            *num_parent_poly3d_courant = i;
            num_parent_poly3d_courant += 1;
            num_cell_parent_to_local_courant[i] = idx_poly3d++;
            break;
          }
        default :
          break;
        }
      }
        
      if (n_poly3d_part > 0) {
        cellfac_poly = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_cellfac_poly);
        
        /* Stockage des faces du bloc */
        
        n_face_poly = 0;
        PDM_l_num_t l_facsom_poly = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] == 0) {
            tag_face_poly3d[i] = n_face_poly++;
            l_facsom_poly += face_som_nb_courant[i];
          }
        }
        
        facsom_poly_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_face_poly + 1));
        facsom_poly = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_facsom_poly);
        
        facsom_poly_idx[0] = 0;
        PDM_l_num_t idx_facsom_poly = 0;
        PDM_l_num_t idx_facsom = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] >= 0) {
            PDM_l_num_t ideb = face_som_idx_courant[i] - adjust;
            PDM_l_num_t ifin = ideb + face_som_nb_courant[i];
            facsom_poly_idx[idx_facsom+1] = facsom_poly_idx[idx_facsom] + face_som_nb_courant[i];
            idx_facsom += 1;
            for (int j = ideb; j < ifin; j++) {
              facsom_poly[idx_facsom_poly++] = face_som_courant[j];
            }
          }
        }

        /* Remplissage de la structure cellfac_poly */

        l_cellfac_poly = 0;
        for (int i = 0; i < n_cell_courant; i++) {
          PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_nb_courant[i],
                                                  cell_face_courant + cell_face_idx_courant[i] - adjust,
                                                  face_som_idx_courant,
                                                  face_som_nb_courant,
                                                  face_som_courant,
                                                  cell_som_tria,
                                                  cell_som_quad);
        
          switch(cell_type) {
            
          case PDM_MESH_NODAL_POLY_3D : 
            {
              PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - adjust;
              for (int j = 0; j < cell_face_nb_courant[i]; j++) {
                cellfac_poly[l_cellfac_poly++] = tag_face_poly3d[cell_face_cell[j] - 1] + 1;
              }
              break;
            }
         default:
            break;
          }
        }
        free(tag_face_poly3d);
      }

      if (som_elts[0] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                                     id_bloc_tetra4,
                                     ipart,
                                     n_tetra_part,
                                     connec_tetra,
                                     numabs_tetra,
                                     num_parent_tetra);

      if (som_elts[1] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                             id_bloc_hexa8,
                             ipart,
                             n_hexa_part,
                             connec_hexa,
                             numabs_hexa,
                             num_parent_hexa);
    
      if (som_elts[2] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                             id_bloc_prism6,
                             ipart,
                             n_prism_part,
                             connec_prism,
                             numabs_prism,
                             num_parent_prism);

      if (som_elts[3] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                             id_bloc_pyramid5,
                             ipart,
                             n_pyramid_part,
                             connec_pyramid,
                             numabs_pyramid,
                             num_parent_pyramid);

      if (som_elts[4] > 0)
        PDM_Mesh_nodal_block_poly3d_set(idx,
                                id_bloc_poly_3d,
                                ipart,
                                n_poly3d_part,
                                n_face_poly,
                                facsom_poly_idx,
                                facsom_poly,
                                cellfac_poly_idx,
                                cellfac_poly,
                                numabs_poly3d,
                                num_parent_poly3d);
    }

    if (mesh->prepa_blocks != NULL) {
      free(mesh->prepa_blocks->n_cell);
      free(mesh->prepa_blocks->n_face);
      free(mesh->prepa_blocks->n_tetra);
      free(mesh->prepa_blocks->n_hexa);
      free(mesh->prepa_blocks->n_prism);
      free(mesh->prepa_blocks->n_pyramid);
      free(mesh->prepa_blocks->n_poly3d);
      free(mesh->prepa_blocks->face_vtx_idx);
      free(mesh->prepa_blocks->face_vtx_nb);
      free(mesh->prepa_blocks->face_vtx);
      free(mesh->prepa_blocks->cell_face_idx);
      free(mesh->prepa_blocks->cell_face_nb);
      free(mesh->prepa_blocks->cell_face);
      free(mesh->prepa_blocks->add_etat);
      free(mesh->prepa_blocks->numabs);
      free(mesh->prepa_blocks);
      mesh->prepa_blocks = NULL;
    }
  }
  
}



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
const int          n_cell,
const int          n_edge,
PDM_l_num_t       *edge_vtx_idx,
PDM_l_num_t       *edge_vtx_nb,
PDM_l_num_t       *edge_vtx,
PDM_l_num_t       *cell_edge_idx,
PDM_l_num_t       *cell_edge_nb,
PDM_l_num_t       *cell_edge,
PDM_g_num_t       *numabs
) 
{
  int adjust = 0;
  if (n_cell > 0) {
    if (edge_vtx_idx[0] == 1) {
      adjust = 1;
    }
  }

  
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
    
  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  int n_part = 0;

  if (mesh->num_cell_parent_to_local == NULL) {
    mesh->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part); 
    for (int ipart = 0; ipart < mesh->n_part; ipart++) { 
      mesh->num_cell_parent_to_local[ipart] = NULL;
    }
  }

  mesh->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    mesh->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (mesh->prepa_blocks == NULL) {
    mesh->prepa_blocks = (PDM_Mesh_nodal_prepa_blocks_t *) malloc(sizeof(PDM_Mesh_nodal_prepa_blocks_t));
    mesh->prepa_blocks->t_add = 2;
    mesh->prepa_blocks->n_tria_proc = 0;    /* Nb de triangles par proc */
    mesh->prepa_blocks->n_quad_proc = 0;    /* Nb de quads par proc */
    mesh->prepa_blocks->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    mesh->prepa_blocks->n_cell = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->l_connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->face_vtx_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part); 
    mesh->prepa_blocks->face_vtx_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*mesh->n_part);
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (mesh->prepa_blocks->t_add != 2) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur Cs_geom_cell2d_cellface_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d = 0;

  for (int i = 0; i < n_cell; i++) {

    PDM_l_num_t n_face_cell = cell_edge_nb[i];
    if (n_face_cell == 3)
      n_tria += 1;
    else if (n_face_cell == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += cell_edge_nb[i]; 
    }
  }
  
  mesh->prepa_blocks->n_tria_proc           += n_tria;
  mesh->prepa_blocks->n_quad_proc           += n_quad;
  mesh->prepa_blocks->n_poly2d_proc         += n_poly2d;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_cell[id_part]        = n_cell;
  mesh->prepa_blocks->n_tria[id_part]        = n_tria; 
  mesh->prepa_blocks->n_quad[id_part]        = n_quad; 
  mesh->prepa_blocks->n_poly2d[id_part]      = n_poly2d;
  mesh->prepa_blocks->l_connec_poly2d[id_part] = l_connec_poly2d;
  mesh->prepa_blocks->face_vtx_idx[id_part]  = edge_vtx_idx;
  mesh->prepa_blocks->face_vtx_nb[id_part]   = edge_vtx_nb;
  mesh->prepa_blocks->face_vtx[id_part]      = edge_vtx;
  mesh->prepa_blocks->cell_face_idx[id_part] = cell_edge_idx;
  mesh->prepa_blocks->cell_face_nb[id_part]  = cell_edge_nb;
  mesh->prepa_blocks->cell_face[id_part]     = cell_edge;
  mesh->prepa_blocks->numabs[id_part]        = numabs;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_face[id_part]        = n_edge;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < mesh->n_part; i++) {
    if (mesh->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (mesh->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = mesh->prepa_blocks->n_tria_proc > 0;
    elts[1] = mesh->prepa_blocks->n_quad_proc > 0;
    elts[2] = mesh->prepa_blocks->n_poly2d_proc > 0;

    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);

    int id_bloc_tria3;
    int id_bloc_quad4;
    int id_bloc_poly_2d;

    if (som_elts[0] > 0)
      id_bloc_tria3 = PDM_Mesh_nodal_block_add (idx,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_TRIA3);

    if (som_elts[1] > 0)
      id_bloc_quad4 = PDM_Mesh_nodal_block_add (idx,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_QUAD4);
    
    if (som_elts[2] > 0)
      id_bloc_poly_2d = PDM_Mesh_nodal_block_add (idx,
                                                  PDM_TRUE,
                                                  PDM_MESH_NODAL_POLY_2D);

    /* Determination de la connectivite de chaque element */

    for (int ipart = 0; ipart < mesh->n_part; ipart++) {

      PDM_l_num_t n_cell_courant = mesh->prepa_blocks->n_cell[ipart];
      PDM_l_num_t *num_cell_parent_to_local_courant = mesh->num_cell_parent_to_local[ipart];
      PDM_l_num_t *face_som_courant = mesh->prepa_blocks->face_vtx[ipart];
      PDM_l_num_t *cell_face_idx_courant = mesh->prepa_blocks->cell_face_idx[ipart];
      PDM_l_num_t *cell_face_nb_courant = mesh->prepa_blocks->cell_face_nb[ipart];
      PDM_l_num_t *cell_face_courant = mesh->prepa_blocks->cell_face[ipart];
      PDM_g_num_t *numabs_courant = mesh->prepa_blocks->numabs[ipart];
   
      n_tria   = mesh->prepa_blocks->n_tria[ipart];
      n_quad    = mesh->prepa_blocks->n_quad[ipart];
      n_poly2d  = mesh->prepa_blocks->n_poly2d[ipart];
      l_connec_poly2d = mesh->prepa_blocks->l_connec_poly2d[ipart];

      PDM_l_num_t *connec_tria = NULL;
      PDM_l_num_t *connec_quad = NULL;
      PDM_l_num_t *connec_poly2d = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;
 
      PDM_g_num_t *numabs_tria = NULL;
      PDM_g_num_t *numabs_quad = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;
    
      PDM_l_num_t *num_parent_tria = NULL;
      PDM_l_num_t *num_parent_quad = NULL;
      PDM_l_num_t *num_parent_poly2d = NULL;

      if (n_tria > 0) {
        connec_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 3 *n_tria);
        numabs_tria = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tria);
        num_parent_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_tria);
      }

      if (n_quad > 0) {
        connec_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 * n_quad);
        numabs_quad = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_quad);
        num_parent_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_quad);
      }

      if (n_poly2d > 0) {
        connec_poly2d_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly2d + 1));
        connec_poly2d_idx[0] = 0;
        connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_connec_poly2d);
        numabs_poly2d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly2d);
        num_parent_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_poly2d);
      }

      PDM_l_num_t *connec_tria_courant = connec_tria;
      PDM_l_num_t *connec_quad_courant = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant = connec_poly2d;

      PDM_g_num_t *numabs_tria_courant = numabs_tria;
      PDM_g_num_t *numabs_quad_courant = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;
      
      PDM_l_num_t *num_parent_tria_courant = num_parent_tria; 
      PDM_l_num_t *num_parent_quad_courant = num_parent_quad; 
      PDM_l_num_t *num_parent_poly2d_courant = num_parent_poly2d; 

      /* Construction de la connectivit� sommet-> arrete */

      PDM_l_num_t *connec_som_are = 
              (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 2 * mesh->vtx[ipart]->n_vtx);

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      for (int j = 0; j < 2 * mesh->vtx[ipart]->n_vtx; j++) {
        connec_som_are[j] = -1;
      }

      for (int i = 0; i < n_cell_courant; i++) {

        PDM_l_num_t ideb = cell_face_idx_courant[i] - adjust;
        PDM_l_num_t n_face_cell = cell_face_nb_courant[i];
        PDM_l_num_t ifin = ideb + n_face_cell;
 
        for (int j = ideb; j < ifin; j++) {
          PDM_l_num_t ifac = PDM_ABS(cell_face_courant[j]) - 1;
          PDM_l_num_t isom1 = face_som_courant[2*ifac] - 1;
          PDM_l_num_t isom2 = face_som_courant[2*ifac+1] - 1;

          if (connec_som_are[2*isom1] == -1)
            connec_som_are[2*isom1] = ifac;
          else
            connec_som_are[2*isom1+1] = ifac;

          if (connec_som_are[2*isom2] == -1)
            connec_som_are[2*isom2] = ifac;
          else
            connec_som_are[2*isom2+1] = ifac;
        }      

        PDM_l_num_t *connec_courant;
        if (n_face_cell == 3) {
          *num_parent_tria_courant = i;
          num_parent_tria_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_face_cell;
        }
        else if (n_face_cell == 4) {
          *num_parent_quad_courant = i;
          num_parent_quad_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_face_cell;
        }
        else {
          *num_parent_poly2d_courant = i;
          num_parent_poly2d_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          connec_courant = connec_poly2d_courant;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) +  n_face_cell;
          connec_poly2d_idx_courant += 1;
          connec_poly2d_courant += n_face_cell;
        }

        /* Remplissage de la connectivite */
        
        PDM_l_num_t idx_som = 0;
        PDM_l_num_t face_courant = PDM_ABS(cell_face_courant[ideb]) - 1;
        PDM_l_num_t isom1 = face_som_courant[2*face_courant] - 1;
        PDM_l_num_t isom_suiv = face_som_courant[2*face_courant + 1] - 1;
        connec_courant[idx_som++] = isom1 + 1;

        while (isom1 != isom_suiv) {
          assert(idx_som <= n_face_cell);
          connec_courant[idx_som++] = isom_suiv + 1;
          
          /* Face suivante */
          
          PDM_l_num_t face_suiv = connec_som_are[2*isom_suiv];
          if (face_suiv == face_courant)
            face_suiv = connec_som_are[2*isom_suiv + 1];
          face_courant = face_suiv;

          /* Sommet suivant */

          PDM_l_num_t isom_tmp = face_som_courant[2*face_courant] - 1;
          if (isom_tmp == isom_suiv)
            isom_tmp = face_som_courant[2*face_courant + 1] - 1;
          isom_suiv = isom_tmp;
        }
        
        for (int j= 0; j < n_face_cell; j++) {
          connec_som_are[2*(connec_courant[j] -1)] = - 1;
          connec_som_are[2*(connec_courant[j] -1) + 1] = - 1;
        }
      }
      
      free(connec_som_are);

      if (som_elts[0] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                                     id_bloc_tria3,
                                     ipart,
                                     n_tria,
                                     connec_tria,
                                     numabs_tria,
                                     num_parent_tria);

      if (som_elts[1] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                                     id_bloc_quad4,
                                     ipart,
                                     n_quad,
                                     connec_quad,
                                     numabs_quad,
                                     num_parent_quad);
    
      if (som_elts[2] > 0)
        PDM_Mesh_nodal_block_poly2d_set(idx,
                                        id_bloc_poly_2d,
                                        ipart,
                                        n_poly2d,
                                        connec_poly2d_idx,
                                        connec_poly2d,
                                        numabs_poly2d,
                                        num_parent_poly2d);
    }
    if (mesh->prepa_blocks != NULL) {
      free(mesh->prepa_blocks->n_cell);
      free(mesh->prepa_blocks->n_face);
      free(mesh->prepa_blocks->n_tria);
      free(mesh->prepa_blocks->n_quad);
      free(mesh->prepa_blocks->n_poly2d);
      free(mesh->prepa_blocks->l_connec_poly2d);
      free(mesh->prepa_blocks->face_vtx_idx);
      free(mesh->prepa_blocks->face_vtx_nb);
      free(mesh->prepa_blocks->face_vtx);
      free(mesh->prepa_blocks->cell_face_idx);
      free(mesh->prepa_blocks->cell_face_nb);
      free(mesh->prepa_blocks->cell_face);
      free(mesh->prepa_blocks->add_etat);
      free(mesh->prepa_blocks->numabs);
      free(mesh->prepa_blocks);
      mesh->prepa_blocks = NULL;
    }
  }
}


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
)
{
  int adjust = 0;
  if (n_face > 0) {
    if (face_vtx_idx[0] == 1) {
      adjust = 1;
    }
  }

  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
    
  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  int n_part = 0;

  if (mesh->num_cell_parent_to_local == NULL) {
    mesh->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part); 
    for (int ipart = 0; ipart < mesh->n_part; ipart++) { 
      mesh->num_cell_parent_to_local[ipart] = NULL;
    }
  }

  mesh->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face);
  for (int i = 0; i < n_face; i++) {
    mesh->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (mesh->prepa_blocks == NULL) {
    mesh->prepa_blocks = (PDM_Mesh_nodal_prepa_blocks_t *) malloc(sizeof(PDM_Mesh_nodal_prepa_blocks_t));
    mesh->prepa_blocks->t_add = 3;
    mesh->prepa_blocks->n_tria_proc = 0;    /* Nb de triangles par proc */
    mesh->prepa_blocks->n_quad_proc = 0;    /* Nb de quads par proc */
    mesh->prepa_blocks->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    mesh->prepa_blocks->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->n_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->l_connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part); 
    mesh->prepa_blocks->face_vtx_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part); 
    mesh->prepa_blocks->face_vtx_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*mesh->n_part);
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (mesh->prepa_blocks->t_add != 3) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur Cs_geom_cell2d_cellface_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d  = 0;

  for (int i = 0; i < n_face; i++) {

    PDM_l_num_t n_som_face = face_vtx_nb[i];
    if (n_som_face == 3)
      n_tria += 1;
    else if (n_som_face == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += n_som_face;
    }
  }

  mesh->prepa_blocks->n_tria_proc           += n_tria;
  mesh->prepa_blocks->n_quad_proc           += n_quad;
  mesh->prepa_blocks->n_poly2d_proc         += n_poly2d;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_tria[id_part]        = n_tria; 
  mesh->prepa_blocks->n_quad[id_part]        = n_quad; 
  mesh->prepa_blocks->n_poly2d[id_part]      = n_poly2d;
  mesh->prepa_blocks->l_connec_poly2d[id_part] = l_connec_poly2d;
  mesh->prepa_blocks->face_vtx_idx[id_part]  = face_vtx_idx;
  mesh->prepa_blocks->face_vtx_nb[id_part]   = face_vtx_nb;
  mesh->prepa_blocks->face_vtx[id_part]      = face_vtx;
  mesh->prepa_blocks->numabs[id_part]        = numabs;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_face[id_part]        = n_face;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < mesh->n_part; i++) {
    if (mesh->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (mesh->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = mesh->prepa_blocks->n_tria_proc > 0;
    elts[1] = mesh->prepa_blocks->n_quad_proc > 0;
    elts[2] = mesh->prepa_blocks->n_poly2d_proc > 0;
    
    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);

    int id_bloc_tria3;
    int id_bloc_quad4;
    int id_bloc_poly_2d;

    if (som_elts[0] > 0) {
      id_bloc_tria3 = PDM_Mesh_nodal_block_add (idx,
                                               PDM_TRUE,
                                               PDM_MESH_NODAL_TRIA3);
    }
    
    if (som_elts[1] > 0) {
      id_bloc_quad4 = PDM_Mesh_nodal_block_add (idx,
                                               PDM_TRUE,
                                               PDM_MESH_NODAL_QUAD4);
    }
    
    if (som_elts[2] > 0) {
      id_bloc_poly_2d = PDM_Mesh_nodal_block_add (idx,
                                                 PDM_TRUE,
                                                 PDM_MESH_NODAL_POLY_2D);
    }
    
    /* Determination de la connectivite de chaque element */

    for (int ipart = 0; ipart < mesh->n_part; ipart++) {

      PDM_l_num_t n_face_courant = mesh->prepa_blocks->n_face[ipart];
      PDM_l_num_t *num_cell_parent_to_local_courant = mesh->num_cell_parent_to_local[ipart];
      PDM_l_num_t *face_som_idx_courant = mesh->prepa_blocks->face_vtx_idx[ipart];
      PDM_l_num_t *face_som_nb_courant = mesh->prepa_blocks->face_vtx_nb[ipart];
      PDM_l_num_t *face_som_courant = mesh->prepa_blocks->face_vtx[ipart];
      PDM_g_num_t *numabs_courant = mesh->prepa_blocks->numabs[ipart];
 
      n_tria   = mesh->prepa_blocks->n_tria[ipart];
      n_quad    = mesh->prepa_blocks->n_quad[ipart];
      n_poly2d  = mesh->prepa_blocks->n_poly2d[ipart];
      l_connec_poly2d  = mesh->prepa_blocks->l_connec_poly2d[ipart];

      PDM_l_num_t *connec_tria = NULL;
      PDM_l_num_t *connec_quad = NULL;
      PDM_l_num_t *connec_poly2d = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;
 
      PDM_g_num_t *numabs_tria = NULL;
      PDM_g_num_t *numabs_quad = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;
          
      PDM_l_num_t *num_parent_tria = NULL;
      PDM_l_num_t *num_parent_quad = NULL;
      PDM_l_num_t *num_parent_poly2d = NULL;


      if (n_tria > 0) {
        connec_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 3 *n_tria);
        numabs_tria = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tria);
        num_parent_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_tria);
      }

      if (n_quad > 0) {
        connec_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 * n_quad);
        numabs_quad = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_quad);
        num_parent_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_quad);
      }

      if (n_poly2d > 0) {
        connec_poly2d_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly2d + 1));
        connec_poly2d_idx[0] = 0;
        connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_connec_poly2d);
        numabs_poly2d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly2d);
        num_parent_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_poly2d);
      }

      PDM_l_num_t *connec_tria_courant = connec_tria;
      PDM_l_num_t *connec_quad_courant = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant = connec_poly2d;

      PDM_l_num_t *num_parent_tria_courant = num_parent_tria;
      PDM_l_num_t *num_parent_quad_courant = num_parent_quad;
      PDM_l_num_t *num_parent_poly2d_courant = num_parent_poly2d;

      PDM_g_num_t *numabs_tria_courant = numabs_tria;
      PDM_g_num_t *numabs_quad_courant = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      for (int i = 0; i < n_face_courant; i++) {
        PDM_l_num_t n_som_face = face_som_nb_courant[i];
        PDM_l_num_t idx_som_face = face_som_idx_courant[i] - adjust;
        PDM_l_num_t *connec_courant;

        if (n_som_face == 3) {
          *num_parent_tria_courant = i;
          num_parent_tria_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_som_face;
        }
        else if (n_som_face == 4) {
          *num_parent_quad_courant = i;
          num_parent_quad_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_som_face;
        }
        else {
          *num_parent_poly2d_courant = i;
          num_parent_poly2d_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) + n_som_face;
          connec_poly2d_idx_courant += 1;
          connec_courant = connec_poly2d_courant;
          connec_poly2d_courant += n_som_face;
        }

        /* Remplissage de la connectivite */

        for (int j = 0; j < n_som_face; j++)
          connec_courant[j] = face_som_courant[idx_som_face++];
      }

      if (som_elts[0] > 0) 
        PDM_Mesh_nodal_block_std_set(idx,
                             id_bloc_tria3,
                             ipart,
                             n_tria,
                             connec_tria,
                             numabs_tria,
                             num_parent_tria);

      if (som_elts[1] > 0)
        PDM_Mesh_nodal_block_std_set(idx,
                             id_bloc_quad4,
                             ipart,
                             n_quad,
                             connec_quad,
                             numabs_quad,
                             num_parent_quad);
    
      if (som_elts[2] > 0) 
        PDM_Mesh_nodal_block_poly2d_set(idx,
                                id_bloc_poly_2d,
                                ipart,
                                n_poly2d,
                                connec_poly2d_idx,
                                connec_poly2d,
                                numabs_poly2d,
                                num_parent_poly2d);
    }
    if (mesh->prepa_blocks != NULL) {
      free(mesh->prepa_blocks->n_face);
      free(mesh->prepa_blocks->n_tria);
      free(mesh->prepa_blocks->n_quad);
      free(mesh->prepa_blocks->n_poly2d);
      free(mesh->prepa_blocks->l_connec_poly2d);
      free(mesh->prepa_blocks->face_vtx_idx);
      free(mesh->prepa_blocks->face_vtx_nb);
      free(mesh->prepa_blocks->face_vtx);
      free(mesh->prepa_blocks->add_etat);
      free(mesh->prepa_blocks->numabs);
      free(mesh->prepa_blocks);
      mesh->prepa_blocks = NULL;
    }
  }
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
    
  int id_gnum = PDM_gnum_create (3, mesh->n_part, PDM_FALSE, 1e-3, mesh->pdm_mpi_comm);

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    PDM_Mesh_nodal_block_poly3d_t *block = 
            (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, _id_block);

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }    

    if (block->numabs_int == NULL) {
      block->numabs_int = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      return;
    }

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_gnum_set_from_parents (id_gnum, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    PDM_Mesh_nodal_block_poly2d_t *block = 
            (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, _id_block);

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }    

    if (block->numabs_int == NULL) {
      block->numabs_int = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      return;
    }

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_gnum_set_from_parents (id_gnum, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  else {
  
    int _id_block = id_block;
  
    PDM_Mesh_nodal_block_std_t *block = 
            (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, _id_block);

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }    

    if (block->numabs_int == NULL) {
      block->numabs_int = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      return;
    }

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_gnum_set_from_parents (id_gnum, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  PDM_gnum_compute (id_gnum);
  
  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    PDM_Mesh_nodal_block_poly3d_t *block = 
            (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, _id_block);
    
    for (int i = 0; i < mesh->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get (id_gnum, i);
    }
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    PDM_Mesh_nodal_block_poly2d_t *block = 
            (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, _id_block);

    for (int i = 0; i < mesh->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get (id_gnum, i);
    }
  }
  
  else {

    int _id_block = id_block;
  
    PDM_Mesh_nodal_block_std_t *block = 
            (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, _id_block);

    for (int i = 0; i < mesh->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get (id_gnum, i);
    }
  }

  PDM_gnum_free (id_gnum, PDM_TRUE);

}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  if (mesh->num_cell_parent_to_local != NULL)
    return mesh->num_cell_parent_to_local[id_part];
  else 
    return NULL;
  
}


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
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  return mesh->n_cell[id_part];

}


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
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  assert(vtx->parent != NULL);
  
  return vtx->parent->_numabs;
}



#ifdef __cplusplus
}
#endif /* __cplusplus */

