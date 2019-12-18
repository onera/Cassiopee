#ifndef __PDM_EDGES_INTERSECT_H__
#define __PDM_EDGES_INTERSECT_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_hash_tab.h"
#include "pdm_line.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

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
 * \enum PDM_edges_get_t
 * \brief Type of return
 *
 */

typedef enum {

  PDM_EDGES_GET_FROM_A,   /*!< Mesh A  */
  PDM_EDGES_GET_FROM_B,   /*!< Mesh B  */
  PDM_EDGES_GET_FROM_AB,  /*!< Mesh AB  */

} PDM_edges_get_t;

/**
 * \enum PDM_line_status_t
 * \brief Polygon status type
 *
 */

typedef enum {

  PDM_EDGES_INTERSECT_MESHA,  /*!< MEsh A  */
  PDM_EDGES_INTERSECT_MESHB,  /*!< Mesh B  */

} PDM_edges_intersect_mesh_t;

/**
 * \enum PDM_edges_intersect_point_t
 * \brief Type of intersection point
 *
 */

typedef enum {

  PDM_EDGES_INTERSECT_POINT_NEW,     /*!< New point : True intersection */
  PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB,  /*!< Point from A  */
  PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA,  /*!< Point from B  */
  PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB,/*!< Point from A and B <=> A == B */

} PDM_edges_intersect_point_t;

/**
 * \struct PDM_edges_intersect_t
 * \brief Manage a set of edge intersections
 *
 *  PDM_edges_intersect_t manages a set of edge intersections
 *
 */

typedef struct _edges_intersect_t PDM_edges_intersect_t;

/**
 * \struct PDM_edges_intersect_res_t
 * \brief Result of the intersection between two edges
 *
 *
 *  PDM_edges_intersect_res_t describes the result of the intersection
 *  between two edges
 *
 */

typedef struct _edges_intersect_res_t PDM_edges_intersect_res_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Create a new \ref PDM_edges_intersect_t object
 *
 * \param [in]   maxGNEdgeA      max global number of edges in mesh A
 * \param [in]   maxGNEdgeN      max global number of edges in mesh B
 * \param [in]   vtxCarLengthTol Absolute tolerance for characteristic length
 * \param [in]   sMSGComm        size of mpicomm
 *
 * \return      A new \ref PDM_edges_intersect_t
 */

PDM_edges_intersect_t *
PDM_edges_intersect_create
(
const PDM_g_num_t maxGNEdgeA,
const PDM_g_num_t maxGNEdgeB,
const double     vtxCarLengthTol,
const PDM_MPI_Comm   comm
);


/**
 *
 * \brief Get result of the intersection
 *
 * \param [in]   ei               Current edges intersection pointer
 * \param [in]   get_t            Type of key to return data
 * \param [in]   nGEdgeA          Global number of meshA edge
 * \param [in]   nGEdgeB          Global number of meshB edge
 * \param [out]  n_intersect      Number of intersections
 *
 * \return    Result of the intersection
 */

PDM_edges_intersect_res_t **
PDM_edges_intersect_get
(
PDM_edges_intersect_t       *ei,
PDM_edges_get_t              get_t,
const PDM_g_num_t             nGEdgeA,
const PDM_g_num_t             nGEdgeB,
int                         *n_intersect
);


/**
 *
 * \brief Perform an intersection between a meshA edge and a meshB edge
 *
 * \param [in]   ei               Current edges intersection pointer
 * \param [in]   nGEdgeA          Global number of meshA edge
 * \param [in]   nGVtxA           Global number of edgeA vertices
 * \param [in]   charLgthVtxA     Characteristic length of edgeA vertices
 * \param [in]   coordsVtxA       Coordinates of edgeA vertices
 * \param [in]   nGEdgeB          Global number of meshB edge
 * \param [in]   nGVtxB           Global number of edgeB vertices
 * \param [in]   charLgthVtxB     Characteristic length of edgeB vertices
 * \param [in]   coordsVtxB       Coordinates of edgeB vertices
 *
 * \return    Result of the intersection
 */

PDM_edges_intersect_res_t *
PDM_edges_intersect_add
(
PDM_edges_intersect_t       *ei,
const PDM_g_num_t             nGEdgeA,
const PDM_g_num_t             nGVtxA[2],
const double                 charLgthVtxA[2],
const double                 coordsVtxA[6],
const PDM_g_num_t             nGEdgeB,
const PDM_g_num_t             nGVtxB[2],
const double                 charLgthVtxB[2],
const double                 coordsVtxB[6]
);


/**
 *
 * \brief Free \ref PDM_edges_intersect_t object
 *
 * \param [in]   ei   Current edges intersection pointer
 *
 * \return     NULL
 */

PDM_edges_intersect_t *
PDM_edges_intersect_free
(
PDM_edges_intersect_t *ei
);

/**
 *
 * \brief Get data intersection
 *
 * \param [in]  eir             Current edges intersection result pointer
 * \param [in]  mesh            Origin mesh \ref PDM_edges_intersect_MESHA or
 *                                          \ref PDM_edges_intersect_MESHB
 * \param [out] nGEdge          Global number of meshA edge
 * \param [out] originEdge      Global number of edge origin
 * \param [out] tIntersect      Intersection type
 * \param [out] nNewPoints      Number of intersection points
 * \param [out] oNewPoints      Origin of intersection points
 * \param [out] link            Linked vertex in linked mesh
 * \param [out] gNum            Global number in overlay mesh
 * \param [out] coords          Coordinates of intersection point
 * \param [out] u               Parameter of the intersections in edges
 *                              parametric coordinates
 *
 */

void
PDM_edges_intersect_res_data_get
(
PDM_edges_intersect_res_t   *eir,
PDM_edges_intersect_mesh_t   mesh,
PDM_g_num_t                  *nGEdge,
PDM_g_num_t                  *originEdge,
PDM_line_intersect_t        *tIntersect,
int                         *nNewPoints,
PDM_edges_intersect_point_t **oNewPoints,
PDM_g_num_t                  **link,
PDM_g_num_t                  **gNum,
double                      **coords,
double                      **u
);


/**
 *
 * \brief Perform edges intersection from two polygons
 *
 * \param [in]    intersect            Edges intersection management
 * \param [in]    nVtxA                Number of polygon A vertices
 * \param [in]    faceToEdgeA          Polygon A face to edge connectivity
 * \param [in]    faceToVtxA           Polygon A face to vertex connectivity
 * \param [in]    faceVtxCooA          Polygon A vertex coordinates
 * \param [in]    faceVtxEpsA          Polygon A vertex characteristic length
 * \param [in]    nVtxB                Number of polygon B vertices
 * \param [in]    faceToEdgeB          Polygon B face to edge connectivity
 * \param [in]    faceToVtxB           Polygon B face to vertex connectivity
 * \param [in]    faceVtxCooB          Polygon B vertex coordinates
 * \param [in]    faceVtxEpsB          Polygon B vertex characteristic length
 *
 */

void
PDM_edges_intersect_poly_add
(
PDM_edges_intersect_t  *ei,
const int               nVtxA,
PDM_g_num_t             *faceToEdgeA,
PDM_g_num_t             *faceToVtxA,
double                 *faceVtxCooA,
double                 *faceVtxEpsA,
const int               nVtxB,
PDM_g_num_t             *faceToEdgeB,
PDM_g_num_t             *faceToVtxB,
double                 *faceVtxCooB,
double                 *faceVtxEpsB
);


/**
 *
 * \brief Remove inconsistencies between processes
 *
 * \param [in]   ei           Current edges intersection pointer
 * \param [in]   nAbsVtxA     Absolute number of vertices in initial A mesh
 * \param [in]   nAbsVtxB     Absolute number of vertices in initial B mesh
 * \param [out]  nAbsNewVtxA  Absolute number of vertices in A mesh after intersections
 * \param [out]  nAbsNewVtxB  Absolute number of vertices in B mesh after intersections
 *
 */

void
PDM_edges_intersect_synchronize
(
PDM_edges_intersect_t       *ei,
PDM_g_num_t             nAbsVtxA,
PDM_g_num_t             nAbsVtxB,
PDM_g_num_t            *nAbsNewVtxA,
PDM_g_num_t            *nAbsNewVtxB
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_EDGES_INTERSECTION_H__ */
