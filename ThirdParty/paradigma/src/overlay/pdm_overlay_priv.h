#ifndef __PDM_OVERLAY_PRIV_H__
#define __PDM_OVERLAY_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_dbbtree.h"
#include "pdm_surf_mesh.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"
#include "pdm_timer.h"
#include "pdm_overlay.h"
#include "pdm_handles.h"
#include "pdm_error.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 23

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum PDM_ol_mesh_t
 * \brief 3 Meshes
 *
 */

typedef enum {

    BEGIN                         = 0,
    INIT_DEF_DATA                 = 1,
    INIT_BUILD_EDGES              = 2,
    INIT_BUILD_EXCH_GRAPH         = 3,
    INIT_COMPUTE_CAR_LGTH         = 4,
    OL_FACE_EXTENTS               = 5,
    OL_IS_PLANE                   = 6,
    OL_BUILD_DBBTREE              = 7,
    OL_BOXES_INTERSECT            = 8,
    OL_DISTRIB_BOXESA_BLOCK       = 9,
    OL_EDGES_INTERSECTION         = 10,
    OL_EDGES_SYNCHRO              = 11,
    OL_CLIPPING                   = 12,
    OL_COMPUTE_ADD_SUB_FACESA     = 13,
    OL_DISTRIB_BOXESB_BLOCK       = 14,
    OL_COMPUTE_ADD_SUB_FACESB     = 15,
    OL_SEND_RESULTS_TO_INIT_PARTA = 16,
    OL_SEND_RESULTS_TO_INIT_PARTB = 17,
    OL_COMPUTE_LOCAL_CONNECTA     = 18,
    OL_COMPUTE_LOCAL_CONNECTB     = 19,
    OL_UPDATE_A_B_CONNECT_GRAPH   = 20,
    OL_SORT_A_B_CONNECT_GRAPH     = 21,
    END                           = 22,

} _ol_timer_step_t;


/**
 * \struct _ol_part_t
 * \brief  Ol mesh partition type
 *
 * \ref _part_t defines a overlay mesh partition structure
 *
 */

typedef struct _ol_part_t{

  PDM_surf_part_t  *initPart;/*!< Initial mesh partition */
  PDM_surf_part_t  *part;    /*!< Mesh partition */
  int       sInitToOlFace;    /*!< Size of \ref InitToOlFace */
  int      *initToOlFaceIdx;  /*!< Array adress of \ref initToOlFace index */
  int      *initToOlFace;     /*!< Array adress of initial to overlay faces */
  int       nLinkedFace;      /*!< Number of linked face */
  int      *linkedFacesProcIdx;/*!< Linked faces (size = nProc + 1) */
  int      *linkedFaces;      /*!< Linked faces (size = 4 * nLinkedFace) */
  int      *faceIniVtxIdx;    /*!< Intial faces with new vertices index */
  int      *faceIniVtx;       /*!< Intial faces with new vertices */

} _ol_part_t;


/**
 * \struct _ol_mesh_t
 * \brief  Mesh type
 *
 * \ref _ol_mesh_t defines a mesh structure
 *
 */

typedef struct _ol_mesh_t {

  PDM_g_num_t   nGFace;         /*!< Global number of faces     */
  PDM_g_num_t   nGVtx;          /*!< Global number of vertices  */
  int          nPart;          /*!< Number of local partitions */
  _ol_part_t **part;           /*!< Mesh partition             */
} _ol_mesh_t;

/**
 * \struct PDM_ol_t
 * \brief  Overlay type
 *
 * PDM_ol_t defines a overlaying structure
 *
 */


 typedef struct PDM_ol_t{

  double   projectCoeff;      /*!< Projection coefficient to define the overlay
                                   surface projection :
                                   If value == 0, the surface projection is MeshA
                                   If value == 1, the surface projection is MeshB
                                   If 0 < value < 1, the projection surface is an
                                   intermediate surface */
  double   vtxCarLengthTol;   /*!< Absolute tolerance used to define local geometric
                                   tolerance for vertex caracteristic lenght
                                   (tolerance > 0) */

  double   extentsTol;   /*!< Absolute tolerance used to define local geometric
                              tolerance for vertex caracteristic lenght
                              (tolerxance > 0) */

  double   samePlaneTol;   /*!< Absolute tolerance used to check if 2 surfaces
                             are the same plane surface */


  PDM_MPI_Comm comm;         /*!< MPI communicator */

  PDM_surf_mesh_t  *meshA;            /*!< Mesh A */
  PDM_surf_mesh_t  *meshB;            /*!< Mesh B */

  _ol_mesh_t  *olMeshA;       /*!< Overlay Mesh A */
  _ol_mesh_t  *olMeshB;       /*!< Overlay Mesh B */

  PDM_dbbtree_t * dbbtreeA;    /*!< Distributed boundary box tree on mesh A */

  PDM_timer_t *timer;


  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */

} PDM_ol_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

static PDM_Handles_t *olArray = NULL; /*!< Array to storage overlay identifiers */

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/**
 * \brief Return an intialized \ref _ol_part_t structure
 *
 * This function returns an initialized \ref _ol_part_t structure
 *
 * \param [in]  initPart        Initial part
 *
 * \return      A new initialized \ref _ol_part_t structure
 *
 */

static inline _ol_part_t *
_ol_part_create
(
PDM_surf_part_t          *initPart
)
{
  _ol_part_t *_ol_part = (_ol_part_t *) malloc(sizeof(_ol_part_t));

  _ol_part->initPart        = (PDM_surf_part_t *) initPart;
  _ol_part->part            = NULL;
  _ol_part->sInitToOlFace   = 0;
  _ol_part->initToOlFaceIdx = NULL;
  _ol_part->initToOlFace    = NULL;
  _ol_part->nLinkedFace     = 0;
  _ol_part->linkedFaces     = NULL;
  _ol_part->faceIniVtxIdx   = NULL;
  _ol_part->faceIniVtx      = NULL;

  return _ol_part;
}


/**
 * \brief Delete a \ref _ol_part_t structure
 *
 * This function deletes a \ref _ol_part_t structure
 *
 * \param [in]  ol_part      part to delete
 *
 * \return     Null pointer
 */

static inline _ol_part_t *
_ol_part_free
(
 _ol_part_t *ol_part
)
{
  if (ol_part != NULL) {
    ol_part->initPart = NULL;
    ol_part->part     = PDM_surf_part_free(ol_part->part);

    if (ol_part->initToOlFaceIdx != NULL) {
      free(ol_part->initToOlFaceIdx);
    }

    if (ol_part->initToOlFace != NULL) {
      free(ol_part->initToOlFace);
    }

    if (ol_part->linkedFacesProcIdx != NULL) {
      free (ol_part->linkedFacesProcIdx);
    }

    if (ol_part->linkedFaces != NULL) {
      free (ol_part->linkedFaces);
    }

    if (ol_part->faceIniVtxIdx != NULL) {
      free (ol_part->faceIniVtxIdx);
    }

    if (ol_part->faceIniVtx != NULL) {
      free (ol_part->faceIniVtx);
    }

    free(ol_part);

  }
  return NULL;
}


/**
 * \brief Delete a \ref _mesh_t structure
 *
 * This function returns an initialized \ref _mesh_t structure
 *
 * \param [in]  mesh         Mesh to delete
 *
 * \return     Null pointer
 *
 */

static inline _ol_mesh_t *
_ol_mesh_free
(
_ol_mesh_t *mesh
)
{

  if (mesh != NULL) {
    if (mesh->part != NULL) {
      for (int i = 0; i < mesh->nPart; i++)
        mesh->part[i] = _ol_part_free(mesh->part[i]);
      free(mesh->part);
      mesh->part = NULL;
    }
    free (mesh);
  }

  return NULL;
}


/**
 * \brief Return an intialized \ref _ol_mesh_t structure
 *
 * This function returns an initialized \ref _ol_mesh_t structure
 *
 * \param [in]  nGface       Number of global faces
 * \param [in]  nGVtx        Number of global vertices
 * \param [in]  nPart        Number of partition
 *
 * \return      A new initialized \ref _ol_mesh_t structure
 *
 */

static inline _ol_mesh_t *
_ol_mesh_create
(
const PDM_g_num_t  nGFace,
const PDM_g_num_t  nGVtx,
const int         nPart
)
{
  _ol_mesh_t *_ol_mesh = (_ol_mesh_t *) malloc(sizeof(_ol_mesh_t));

  _ol_mesh->nGFace  = nGFace;
  _ol_mesh->nGVtx   = nGVtx;
  _ol_mesh->nPart   = nPart;
  _ol_mesh->part = (_ol_part_t **) malloc(nPart * sizeof(_ol_part_t *));


  return _ol_mesh;
}


/**
 * \brief Return an overlay object from its identifier
 *
 * This function returns an overlay object from its identifier
 *
 * \param [in]  id          ol identifier
 *
 */

static inline PDM_ol_t *
_ol_get
(
const int  id
)
{
  PDM_ol_t *ol = (PDM_ol_t *) PDM_Handles_get (olArray, id);

  if (ol == NULL) {
    PDM_error (__FILE__, __LINE__, 0,  "Error _ol_get : '%i' Undefined identifier\n", id);
  }

  return ol;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_OVERLAY_PRIV_H__ */
