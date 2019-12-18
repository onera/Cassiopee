#ifndef __PDM_OVERLAY_H__
#define __PDM_OVERLAY_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Macro pour la manipulation des noms des fonctions entre C et Fortran pour
 * le link.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

/*----------------------------------------------------------------------------
 * Macro permettant de prendre en compte les arguments caches des longueurs
 * de chaine du Fortran (Evite les erreurs de compilation)
 *----------------------------------------------------------------------------*/

#if defined (__uxpv__)
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_ol_mesh_t
 * \brief 3 Meshes
 *
 */

typedef enum {

  PDM_OL_MESH_A       = 0,  /*!< First mesh to overlay */
  PDM_OL_MESH_B       = 1,  /*!< Second mesh to overlay */

} PDM_ol_mesh_t;


/**
 * \enum PDM_ol_parameter_t
 * \brief Parameters for ovelay meshes building
 *
 */

typedef enum {

  PDM_OL_CAR_LENGTH_TOL = 0,  /*!< Absolute tolerance for caracteristic length */
  PDM_OL_EXTENTS_TOL    = 1,  /*!< Absolute tolerance for extents */
  PDM_OL_SAME_PLANE_TOL = 2,  /*!< Absolute tolerance for check if 2 surfaces are
                                   the same plane surface*/

} PDM_ol_parameter_t;

/**
 * \enum PDM_ol_mv_t
 * \brief Type of moving mesh
 *
 */

typedef enum {

  PDM_OL_MV_TRANSFORMATION  = 0,  /*!< Moving with combination of geometric transformations */
  PDM_OL_MV_UNKNOWN         = 1,  /*!< Unknown moving type */

} PDM_ol_mv_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Build and initialize an overlaying object
 *
 * This function builds an initializes an overlaying surface meshes object
 *
 * \param [in]  nPartMeshA   Number of local partitions of the meshA input
 * \param [in]  nGFaceA      Number of global faces of the meshA input
 * \param [in]  nGVtxA       Number of global vertices of the meshA input
 * \param [in]  nPartMeshB   Number of local partitions of the meshB input
 * \param [in]  nGFaceB      Number of global faces of the meshB input
 * \param [in]  nGVtxB       Number of global vertices of the meshB input
 * \param [in]  projectCoeff Projection coefficient to define the overlay surface projection
 *                           If value == 0, the surface projection is MeshA
 *                           If value == 1, the surface projection is MeshB
 *                           If 0 < value < 1 , the projection surface is an intermediate surface
 * \param [in]  comm         MPI communicator.
 *
 * \return      id           Overlay object identifier.
 *
 */

int
PDM_ol_create
(
 const int         nPartMeshA,
 const PDM_g_num_t  nGFaceMeshA,
 const PDM_g_num_t  nGVtxMeshA,
 const int         nPartMeshB,
 const PDM_g_num_t  nGFaceMeshB,
 const PDM_g_num_t  nGVtxMeshB,
 const double      projectCoeff,
 const PDM_MPI_Comm    comm
);

void
PROCF (pdm_ol_create, PDM_OL_CREATE)
(
 int        *nPartMeshA,
 PDM_g_num_t *nGFaceMeshA,
 PDM_g_num_t *nGVtxMeshA,
 int        *nPartMeshB,
 PDM_g_num_t *nGFaceMeshB,
 PDM_g_num_t *nGVtxMeshB,
 double     *projectCoeff,
 PDM_MPI_Fint   *comm,
 int        *id
);


/**
 * \brief Set an overlay parameter
 *
 * This function sets en overlay parameter
 *
 * \param [in]  id          PDM_ol identifier
 * \param [in]  parameter   Parameter to define
 * \param [in]  value       Parameter value
 *
 */

void
PDM_ol_parameter_set
(
 const int                id,
 const PDM_ol_parameter_t parameter,
 const double             value
);

void
PROCF (pdm_ol_parameter_set, PDM_OL_PARAMETER_SET)
(
 int                *id,
 PDM_ol_parameter_t *parameter,
 double             *value
);

/**
 * \brief Define input meshes properties
 *
 * This function defines the input meshes properties
 *
 * \param [in]  id          PDM_ol identifier
 * \param [in]  mesh        Input mesh to define
 *                          (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  ipart       Partition to define
 * \param [in]  nFace       Number of faces
 * \param [in]  faceVtxIdx  Index in the face -> vertex connectivity
 * \param [in]  faceVtxIdx  face -> vertex connectivity
 * \param [in]  faceLnToGn  Local face numbering to global face numbering
 * \param [in]  nVtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtxLnToGn   Local vertex numbering to global vertex numbering
 *
 */

void
PDM_ol_input_mesh_set
(
 const int           id,
 const PDM_ol_mesh_t mesh,
 const int           ipart,
 const int           nFace,
 const int          *faceVtxIdx,
 const int          *faceVtx,
 const PDM_g_num_t   *faceLnToGn,
 const int           nVtx,
 const double       *coords,
 const PDM_g_num_t   *vtxLnToGn
);

void
PROCF (pdm_ol_input_mesh_set, PDM_OL_INPUT_MESH_SET)
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 int           *ipart,
 int           *nFace,
 int           *faceVtxIdx,
 int           *faceVtx,
 PDM_g_num_t    *faceLnToGn,
 int           *nVtx,
 double        *coords,
 PDM_g_num_t    *vtxLnToGn
);


/**
 * \brief Define the type of a mesh moving
 *
 * This function defines the type of a mesh moving.
 * Only a mesh can move
 *
 * \param [in]  id       PDM_ol identifier
 * \param [in]  mesh     Moving mesh
 * \param [in]  mv       Type of moving
 *
 */

void
PDM_ol_moving_type_set
(
 const int           id,
 const PDM_ol_mesh_t mesh,
 const PDM_ol_mv_t   mv
);


void
PROCF (pdm_ol_moving_type_set, PDM_OL_MOVING_TYPE_SET)
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 PDM_ol_mv_t   *mv
);


/**
 * \brief Define a translation
 *
 * This function defines a translation for the moving mesh
 *
 * \param [in]  id       PDM_overlay identifier
 * \param [in]  vect     Translation vector
 * \param [in]  center   Translation center
 *
 */

void
PDM_ol_translation_set
(
 const int           id,
 const double       *vect,
 const double       *center
 );


void
PROCF (pdm_ol_translation_set, PDM_OL_TRANSLATION_SET)
(
 int          *id,
 double       *vect,
 double       *center
);


/**
 * \brief Define a rotation
 *
 * This function defines a rotation for the moving mesh
 *
 * \param [in]  id        PDM_ol identifier
 * \param [in]  direction Rotation direction
 * \param [in]  center    Rotation center
 * \param [in]  angle     Rotation center (degrees)
 *
 */

void
PDM_ol_rotation_set
(
 const int      id,
 const double  *direction,
 const double  *center,
 const double   angle
);


void
PROCF (pdm_ol_rotation_set, PDM_OL_ROTATION_SET)
(
 int     *id,
 double  *direction,
 double  *center,
 double  *angle
);

/**
 * \brief Overlaying the input surface meshes
 *
 * This function overlays the input surface meshes
 *
 * \param [in]  id       PDM_ol identifier
 *
 */

void
PDM_ol_compute
(
 const int          id
);


void
PROCF (pdm_ol_compute, PDM_OL_COMPUTE)
(
 int         *id
);


/**
 * \brief Return the entitie sizes of the overlay mesh
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id        PDM_ol identifier
 * \param [in]  mesh      Input mesh
 *                        (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [out] nGOlFace  Global number of faces of the overlay mesh
 * \param [out] nGOlVtx   Global number of vertices of the overlay mesh
 *
 */


void
PDM_ol_mesh_dim_get
(
 const int            id,
 const PDM_ol_mesh_t  mesh,
       PDM_g_num_t    *nGOlFace,
       PDM_g_num_t    *nGOlVtx
);


void
PROCF (pdm_ol_mesh_dim_get, PDM_OL_MESH_DIM_GET)
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 PDM_g_num_t    *nGOlFace,
 PDM_g_num_t    *nGOlVtx
);


/**
 * \brief Return the entitie sizes of the overlay mesh
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id            PDM_ol identifier
 * \param [in]  mesh          Input mesh
 *                            (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  ipart         Partition to define
 * \param [out] nOlFace       Number of faces of the overlay mesh
 * \param [out] nOlLinkedFace Number of linked faces
 * \param [out] nOlVtx        Number of vertices of the overlay mesh
 * \param [out] sOlFaceVtx    Size of olFaceVtx for each partition
 * \param [out] sInitToOlFace Size of initToOlFace for each partition
 *
 */

void
PDM_ol_part_mesh_dim_get
(
 const int            id,
 const PDM_ol_mesh_t  mesh,
 const int            ipart,
       int           *nOlFace,
       int           *nOlLinkedFace,
       int           *nOlVtx,
       int           *sOlFaceVtx,
       int           *sInitToOlFace
);


void
PROCF (pdm_ol_part_mesh_dim_get, PDM_OL_PART_MESH_DIM_GET)
(
 int           *id,
 PDM_ol_mesh_t *mesh,
 int           *ipart,
 int           *nOlFace,
 int           *nOlLinkedFace,
 int           *nOlVtx,
 int           *sOlFaceVtx,
 int           *sInitToOlFace
);


/**
 * \brief Return the entitie of the overlay mesh
 *
 * This function returns the entities of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id              PDM_ol identifier
 * \param [in]  mesh            Input mesh
 *                              (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  ipart           Mesh partition identifier
 * \param [out] olFaceVtxIdx    Array adress of \ref olFaceVtx index
 *                              (size : \ref nOlFace + 1)
 * \param [out] olFaceVtx       Array adress of face vertex connectivity
 *                              (size : \ref sOlFaceVtx[\ref ipart])
 * \param [out] olLinkedFaceProcIdx olLinkedFace Index (size = nProc + 1)
 * \param [out] olLinkedFace    Array adress of linked face in other mesh
 *                              For each face, 4 link properties :
 *                                    - local face number
 *                                    - linked process,
 *                                    - linked part number,
 *                                    - linked local face number
 *                              (size : \ref 4 * nOlLinkedFace)
 * \param [out] olFaceLnToGn    Array adress of local to global face numbering
 *                              (size : \ref nOlFace)
 * \param [out] olCoords        Array adress of vertex coodinates
 *                              (size : 3 * \ref nOlVtx)
 * \param [out] olVtxLnToGn     Array adress of local to global vertex numbering array
 *                              (size : \ref nOlVtx)
 * \param [out] initToOlFaceIdx Array adress of \ref initToOlFace index
 *                              (size : \ref nOlVtx + 1)
 * \param [out] initToOlFace    Array adress of initial to overlay faces
 * \param [out] initToOlVtx     Array adress of initial to overlay vertices
 *
 */

void
PDM_ol_mesh_entities_get
(
 const int              id,
 const PDM_ol_mesh_t    mesh,
 const int              iPart,
 const int            **olFaceIniVtxIdx,
 const int            **olFaceIniVtx,
 const int            **olFaceVtxIdx,
 const int            **olFaceVtx,
 const int            **olLinkedFaceProcIdx,
 const int            **olLinkedFace,
 const PDM_g_num_t     **olFaceLnToGn,
 const double         **olCoords,
 const PDM_g_num_t     **olVtxLnToGn,
 const int            **initToOlFaceIdx,
 const int            **initToOlFace
);

void
PROCF (pdm_ol_mesh_entities_get, PDM_OL_MESH_ENTITIES_GET)
(
 int            *id,
 PDM_ol_mesh_t  *mesh,
 int            *ipart,
 int            *olFaceIniVtxIdx,
 int            *olFaceIniVtx,
 int            *olFaceVtxIdx,
 int            *olFaceVtx,
 int            *olLinkedFaceProcIdx,
 int            *olLinkedFace,
 PDM_g_num_t     *olFaceLnToGn,
 double         *olCoords,
 PDM_g_num_t     *olVtxLnToGn,
 int            *initToOlFaceIdx,
 int            *initToOlFace
);

/**
 * \brief Delete an overlay object
 *
 * This function deletes an overlay object
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

void
PDM_ol_del
(
 const int     id
);

void
PROCF (pdm_ol_del, PDM_OL_DEL)
(
 int     *id
);


/**
 * \brief Dump elapsed an CPU time
 *
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

void
PDM_ol_dump_times
(
 const int     id
);

void
PROCF (pdm_ol_dump_times, PDM_OL_DUMP_TIMES)
(
 int     *id
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_OVERLAY_H__ */
