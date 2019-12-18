/*
 * This file is an implementation of Greiner & Hormann algorithm with Foster
 * Overleft extension to remove degenerate cases
 *
 */

#ifndef __PDM_POLYGON_CLIPP_H__
#define __PDM_POLYGON_CLIPP_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_edges_intersect.h"

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
 * \enum PDM_edges_intersect_point_t
 * \brief Type of intersection point
 *
 */

typedef enum {

  PDM_POLY_CLIPP_CLIP,            /*!< Perform clipped polygons */
  PDM_POLY_CLIPP_REVERSE,         /*!< Perform opposit polygons */

} PDM_poly_clipp_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Perform polygon clipping
 *
 * \param [in]    ei                   Edges intersection management
 * \param [in]    gNumA                Polygon A global number
 * \param [in]    nVtxA                Number of polygon A vertices
 * \param [in]    faceToEdgeA          Polygon A face to edge connectivity
 * \param [in]    faceToVtxA           Polygon A face to vertex connectivity
 * \param [in]    faceVtxCooA          Polygon A vertex coordinates
 * \param [in]    faceVtxEpsA          Polygon A vertex characteristic length
 * \param [in]    gNumB                Polygon A global number
 * \param [in]    nVtxB                Number of polygon B vertices
 * \param [in]    faceToEdgeB          Polygon B face to edge connectivity
 * \param [in]    faceToVtxB           Polygon B face to vertex connectivity
 * \param [in]    faceVtxCooB          Polygon B vertex coordinates
 * \param [in]    faceVtxEpsB          Polygon B vertex characteristic length
 * \param [in]    performed_t          Type of performed polygon
 * \param [out]   nPolyClippA           Number of clipped polygon
 * \param [out]   polyClippIdxA         Connectivity index for each polygon
 *                                     size = nPolyClipp + 1
 * \param [out]   polyClippConnecA     Connectivity of each clipped polygon
 *                                     size = polyClippIdx[nPolyClipp]
 * \param [out]   polyClippCoordsA     Vertices coordinates of clipping polygon
 * \param [out]   nPolyClippB          Number of clipped polygon
 * \param [out]   polyClippIdxB         Connectivity index for each polygon
 *                                     size = nPolyClipp + 1
 * \param [out]   polyClippConnecB     Connectivity of each clipped polygon
 *                                     size = polyClippIdx[nPolyClipp]
 * \param [out]   polyClippCoordsB     Vertices coordinates of clipping polygon
 *
 */

void
PDM_poly_clipp
(
PDM_edges_intersect_t  *ei,
const int               nVtxA,
PDM_g_num_t             *faceToEdgeA,
PDM_g_num_t             *faceToVtxA,
double                 *faceVtxCooA,
const int               nVtxB,
PDM_g_num_t             *faceToEdgeB,
PDM_g_num_t             *faceToVtxB,
double                 *faceVtxCooB,
PDM_poly_clipp_t        performed_t,
int                    *nPolyClippA,
int                   **polyClippIdxA,
PDM_g_num_t            **polyClippConnecA,
double                **polyClippCoordsA,
int                    *nPolyClippB,
int                   **polyClippIdxB,
PDM_g_num_t            **polyClippConnecB,
double                **polyClippCoordsB
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POLY_CLIPP_H__ */
