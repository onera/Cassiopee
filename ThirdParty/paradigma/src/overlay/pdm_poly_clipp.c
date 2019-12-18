/*
 * This file is an implementation of Greiner & Hormann algorithm with Foster
 * Overleft extension to remove degenerate cases
 *
 */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_poly_clipp.h"
#include "pdm_poly_clipp_priv.h"
#include "pdm_binary_search.h"
#include "pdm_plane.h"
#include "pdm_polygon.h"
#include "pdm_priv.h"
#include "pdm_edges_intersect.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/


/**
 *
 * \brief Create a new vertex linked list
 *
 * \param [in]    coords   Vertex coordinates
 * \param [in]    charLgthVtxA      Characteristic length vertex
 * \param [in]    type         Type of point
 * \param [in]    gN       Global number of vertex
 * \param [in]    gNEdge   Global number of edge
 *
 * \return new \ref _vertex_poly_t object
 *
 */

static _vertex_poly_t *
_poly_clipp_new
(
const double *coords,
const PDM_g_num_t gN,
const PDM_g_num_t gNEdge
)
{
  _vertex_poly_t *vtxp = malloc (sizeof(_vertex_poly_t));

  vtxp->coords       = coords;
  vtxp->u            = 0.;
  vtxp->gN           = gN;
  vtxp->gNEdge       = gNEdge;
  vtxp->next         = vtxp;
  vtxp->previous     = vtxp;
  vtxp->used         = false;
  vtxp->isect        = false;
  vtxp->tag          = false;
  vtxp->neighbor     = NULL;
  vtxp->first        = vtxp;
  return vtxp;
}


/**
 *
 * \brief Add a new vertex in linked list
 *
 * \param [in]    coords       Vertex coordinates
 * \param [in]    charLgthVtxA Characteristic length vertex
 * \param [in]    type         Type of point
 * \param [in]    gN           Global number of vertex
 * \param [in]    gNEdge       Global number of edge
 * \param [in]    loc          Location of new vertex
 * \param [in]    linked_vtxp  linked element
 *
 *
 * \return new added \ref _vertex_poly_t object
 *
 */

static _vertex_poly_t *
_poly_clipp_add
(
const double      *coords,
const PDM_g_num_t   gN,
const PDM_g_num_t   gNEdge,
const _poly_clipp_loc_t  loc,
_vertex_poly_t    *linked_vtxp
)
{
  _vertex_poly_t *vtxp = malloc (sizeof(_vertex_poly_t));

  vtxp->coords       = coords;
  vtxp->u            = 0.;
  vtxp->gN           = gN;
  vtxp->gNEdge       = gNEdge;
  vtxp->used         = false;
  vtxp->first        = linked_vtxp->first;

  if (loc == POLY_CLIPP_LOC_NEXT) {
    vtxp->previous       = linked_vtxp;
    vtxp->next           = linked_vtxp->next;
    linked_vtxp->next        = vtxp;
    vtxp->next->previous = vtxp;
  }
  else {
    vtxp->next            = linked_vtxp;
    vtxp->previous        = linked_vtxp->previous;
    linked_vtxp->previous = vtxp;
    vtxp->previous->next  = vtxp;
  }

  vtxp->isect    = false;
  vtxp->tag      = false;
  vtxp->neighbor = NULL;
  return vtxp;
}


/**
 *
 * \brief Add a new intersection in linked list
 *
 * \param [in]    u            Parameterization in current edge
 * \param [in]    i            New point index
 * \param [in]    type         Type of point
 * \param [in]    loc          Location of new vertex
 * \param [in]    linked_vtxp  Linked element
 *
 *
 * \return new added \ref _vertex_poly_t object
 *
 */

static _vertex_poly_t *
_poly_clipp_intersect_add
(
const double             u,
const PDM_g_num_t         gN,
const _poly_clipp_loc_t  loc,
_vertex_poly_t          *linked_vtxp
)
{
  _vertex_poly_t *vtxp = malloc (sizeof(_vertex_poly_t));

  vtxp->coords       = NULL;
  vtxp->u            = u;
  vtxp->gN           = gN;
  vtxp->gNEdge       = -1;
  vtxp->used         = false;
  vtxp->first        = linked_vtxp->first;

  if (loc == POLY_CLIPP_LOC_NEXT) {
    vtxp->previous        = linked_vtxp;
    vtxp->next            = linked_vtxp->next;
    linked_vtxp->next     = vtxp;
    vtxp->next->previous  = vtxp;
  }
  else {
    vtxp->next            = linked_vtxp;
    vtxp->previous        = linked_vtxp->previous;
    linked_vtxp->previous = vtxp;
    vtxp->previous->next  = vtxp;
  }

  vtxp->isect    = true;
  vtxp->tag      = true;
  vtxp->neighbor = NULL;
  return vtxp;
}


/**
 *
 * \brief Link neighbors
 *
 * \param [in] vtxpA  A intersection
 * \param [in] vtxpB  B intersection
 *
 */

static void
_poly_clipp_link_neighbors
(
 _vertex_poly_t * vtxpA,
 _vertex_poly_t * vtxpB
)
{
  vtxpA->neighbor = vtxpB;
  vtxpB->neighbor = vtxpA;
  vtxpA->isect    = true;
  vtxpB->isect    = true;
  vtxpA->tag      = true;
  vtxpB->tag      = true;
}


/**
 *
 * \brief Remove a vertex in linked list
 *
 * \param [in]  vtxp   Vertex to remove
 *
 * \return NULL
 *
 */

//static _vertex_poly_t *
//_poly_clipp_remove
//(
// _vertex_poly_t * vtxp
//)
//{
//  if (vtxp != NULL) {
//    vtxp->previous->next = vtxp->next;
//    vtxp->next->previous = vtxp->previous;
//  }
//  free (vtxp);
//  return NULL;
//}



/**
 *
 * \brief Unset intersection vertex (Keep old neighbor in memory)
 *
 * \param [inout] current Current intersection
 *
 */

static void
_poly_clipp_unset_intersect
(
 _vertex_poly_t *current
)
{
  current->isect           = false;
  current->neighbor->isect = false;
}


/**
 *
 * \brief Free a linked list
 *
 * \return NULL
 *
 */


//FIXME: Voir l'appel ou faire l'appel a _poly_clipp_free

//static _vertex_poly_t *
//_poly_clipp_free
//(
// _vertex_poly_t * vtxp
//)
//{
//  _vertex_poly_t *init = vtxp;
//  _vertex_poly_t *next = vtxp->next;
//
//  while (next != init) {
//    _poly_clipp_remove (next);
//    next = next->next;
//
//  }
//  return _poly_clipp_remove (init);
//}


/**
 *
 * \brief Perform vertex location (in, out) for no 'on' vertex
 *                (ray tracing algorithm)
 *
 * \param [inout] vtxpA        Subject polygon
 * \param [inout] vtxpB        Constraint polygon
 * \param [in]    faceVtxCooA  A vertex coordinates
 * \param [in]    faceVtxCooB  B vertex coordinates
 * \param [in]    nVtxA        Number of A vertices
 * \param [in]    nVtxB        Number of B vertices
 * \param [in]    nA           A normal
 * \param [in]    nB           B normal
 *
 */

static void
_location
(
_vertex_poly_t *vtxpA,
_vertex_poly_t *vtxpB,
double *faceVtxCooA,
double *faceVtxCooB,
int nVtxA,
int nVtxB,
double nA[3],
double nB[3]
)
{
  double *boundsA = PDM_bounds_get (nVtxA, faceVtxCooA);
  double *boundsB = PDM_bounds_get (nVtxB, faceVtxCooB);

  _vertex_poly_t *vtx_currB = vtxpB;
  for (int i = 0; i < nVtxB; i++) {
    const double *_coo = faceVtxCooB + 3*i;
    PDM_polygon_status_t stat = PDM_polygon_point_in (_coo,
                                                      nVtxA,
                                                      faceVtxCooA,
                                                      boundsA,
                                                      nA);
    if (stat == PDM_POLYGON_INSIDE) {
      vtx_currB->tag = true;
    }
    else if (stat == PDM_POLYGON_OUTSIDE) {
      vtx_currB->tag = false;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : degenerated polygon\n");
      abort();
    }
    vtx_currB = vtx_currB->next;
  }

  _vertex_poly_t *vtx_currA = vtxpA;
  for (int i = 0; i < nVtxA; i++) {
    const double *_coo = faceVtxCooA + 3*i;
    PDM_polygon_status_t stat = PDM_polygon_point_in (_coo,
                                                      nVtxB,
                                                      faceVtxCooB,
                                                      boundsB,
                                                      nB);
    if (stat == PDM_POLYGON_INSIDE) {
      vtx_currA->tag = true;
    }
    else if (stat == PDM_POLYGON_OUTSIDE) {
      vtx_currA->tag = false;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : degenerated polygon\n");
      abort();
    }
    vtx_currA = vtx_currA->next;
  }

  free (boundsA);
  free (boundsB);

}


/**
 *
 * \brief Tag intersection vertices
 *
 * \param [inout] vtxpA        Subject polygon
 * \param [inout] vtxpB        Constraint polygon
 *
 */

static void
_tag_current
(
_vertex_poly_t *current
)
{

  _vertex_poly_t *neighbor = current->neighbor;

  /*
   * on/on
   */

  if (current->previous->isect && current->next->isect) {

    /*
     * Determine what to do based on the neighbor
     * en tag is the opposite of the neighbor's tag
     */

    if (current->previous->isect && current->next->isect) {
      _poly_clipp_unset_intersect (current);
      current->tag = true;
      neighbor->tag = true;
    }

    else if (neighbor->previous->isect && !neighbor->next->tag) {
      current->tag = false;
    }

    else if (neighbor->previous->isect && neighbor->next->tag) {
      current->tag = true;
    }

    else if (!neighbor->previous->tag && neighbor->next->isect) {
      current->tag = false;
    }

    else if (!(neighbor->previous->tag || neighbor->next->tag)) {
      _poly_clipp_unset_intersect (current);
      current->tag = true;
      neighbor->tag = false;
    }

    else if (!neighbor->previous->tag && neighbor->next->tag) {
      current->tag = false;
    }

    else if (neighbor->previous->tag && neighbor->next->isect) {
      current->tag = true;
    }

    else if (neighbor->previous->tag && !neighbor->next->tag) {
      current->tag = true;
    }

    else if (neighbor->previous->tag && neighbor->next->tag) {
      _poly_clipp_unset_intersect (current);
      current->tag = false;
      neighbor->tag = true;
    }
  }

  /*
   * on/out
   */

  else if (current->previous->isect && !current->next->tag) {
    current->tag = false;
  }

  /*
   * on/in
   */

  else if (current->previous->isect && current->next->tag) {
    current->tag = true;
  }

  /*
   * out/on
   */

  else if (!current->previous->tag && current->next->isect) {
    current->tag = true;
  }

  /*
   * out/out
   */

  else if (!(current->previous->tag || current->next->tag)) {
    if (neighbor->previous->isect && neighbor->next->isect) {
      _poly_clipp_unset_intersect (current);
      neighbor->tag = true;
    }

    else if (neighbor->previous->tag && neighbor->next->tag) {
      _poly_clipp_unset_intersect (current);
    }

    else {
      if (neighbor->previous->tag && !neighbor->next->tag) {
        current->tag = true;
      }
      else {
        current->tag = false;
      }
    }
  }

  /*
   * out/in
   */

  else if (!current->previous->tag && current->next->tag) {
    current->tag = true;
  }

  /*
   * in/on
   */

  else if (current->previous->tag && current->next->isect) {
    current->tag = false;
  }

  /*
   * in/out
   */

  else if (current->previous->tag && !current->next->tag) {
    current->tag = false;
  }

  /*
   * in/in
   */

  else if (current->previous->tag && current->next->tag) {
    if (neighbor->previous->isect && neighbor->next->isect) {
      _poly_clipp_unset_intersect (current);
      neighbor->tag = false;
    }
    else if (neighbor->previous->tag == neighbor->next->tag) {
      _poly_clipp_unset_intersect (current);
    }
    else {
      if (neighbor->previous->tag && !neighbor->next->tag) {
        current->tag = true;
      }
      else {
        current->tag = false;
      }
    }
  }

}


/**
 *
 * \brief Tag intersection vertices
 *
 * \param [inout] vtxpA   Subject polygon
 * \param [inout] vtxpB   Constraint polygon
 *
 */


//FIXME :Fonction _tag : Voir s'il faut faire le travail pour vtxpB

static void
_tag
(
_vertex_poly_t *vtxpA,
_vertex_poly_t *vtxpB
)
{
  vtxpB;

  _vertex_poly_t *vtx_currA = vtxpA;

  do {

    _vertex_poly_t *nextA = vtx_currA->next;
    if (vtx_currA->isect) {
      _tag_current (vtx_currA);
      if (vtx_currA->isect) {
        _tag_current (vtx_currA->neighbor);

        if (vtx_currA->neighbor->tag == vtx_currA->tag) {
          _poly_clipp_unset_intersect (vtx_currA);
        }
      }
    }

    vtx_currA = nextA;

  } while (vtx_currA != vtxpA);

}

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
)
{
  PDM_g_num_t *_faceToEdgeA = faceToEdgeA;
  PDM_g_num_t *_faceToVtxA  = faceToVtxA;
  double     *_faceVtxCooA = faceVtxCooA;

  PDM_g_num_t *_faceToEdgeB = faceToEdgeB;
  PDM_g_num_t *_faceToVtxB  = faceToVtxB;
  double     *_faceVtxCooB = faceVtxCooB;

  /*
   * Compute Normal
   *
   */

  double nA[3];
  PDM_plane_normal (nVtxA, faceVtxCooA, nA);

  double nB[3];
  PDM_plane_normal (nVtxB, faceVtxCooB, nB);

  double dot = PDM_DOT_PRODUCT (nA, nB);

  bool revert = false;

  if (dot < 0) {
    revert = true;
  }

  /*
   * Reorient if necessary
   *
   */

  if (revert) {

    _faceToEdgeB = malloc (sizeof(PDM_g_num_t) * nVtxB);
    _faceToVtxB  = malloc (sizeof(PDM_g_num_t) * nVtxB);
    _faceVtxCooB = malloc (sizeof(double) * 3 * nVtxB);

    int j = nVtxB - 1;
    for (int i = 0; i < nVtxB; i++) {
      _faceToEdgeB[i] = faceToEdgeB[j];
      _faceToVtxB[i] = faceToVtxB[j];
      for (int k = 0; k < 3; k++) {
        _faceVtxCooB[3*i+k] = faceVtxCooB[3*j+k];
      }
      j += -1;
    }

  }

  /*
   * Create double linked list vertex structures
   *
   */

  _vertex_poly_t *vtxA = _poly_clipp_new (_faceVtxCooA,
                                          *_faceToVtxA,
                                          *_faceToEdgeA);

  _vertex_poly_t *vtxB = _poly_clipp_new (_faceVtxCooB,
                                          *_faceToVtxB,
                                          *_faceToEdgeB);

  for (int i = 1; i < nVtxA; i++) {
    const double *_coo   = _faceVtxCooA + 3*i;
    PDM_g_num_t   *_gN     = _faceToVtxA +i;
    PDM_g_num_t   *_gNEdge = _faceToEdgeA +i;

    _poly_clipp_add (_coo,
                     *_gN,
                     *_gNEdge,
                     POLY_CLIPP_LOC_NEXT,
                     vtxA);
  }

  for (int i = 1; i < nVtxB; i++) {
    const double *_coo = _faceVtxCooB + 3*i;
    PDM_g_num_t   *_gN     = _faceToVtxB +i;
    PDM_g_num_t   *_gNEdge = _faceToEdgeB +i;

    _poly_clipp_add (_coo,
                     *_gN,
                     *_gNEdge,
                     POLY_CLIPP_LOC_NEXT,
                     vtxB);
  }

  /*
   *   - Add intersection into linked list
   *   - Move to Intersection tag for Vertices located on Polygon
   *
   */

  _vertex_poly_t *vtx_currB;

  _vertex_poly_t *vtx_currA = vtxA;

  do {

    _vertex_poly_t *nextA = vtx_currA->next;
    if (!vtx_currA->isect) {

      while (!nextA->isect) {
        nextA = nextA->next;
      }

      vtx_currB = vtxB;

      do {

        _vertex_poly_t *nextB = vtx_currB->next;
        if (!vtx_currB->isect) {

          while (!nextB->isect) {
            nextB = nextB->next;
          }

          /*
           * Look for intersection : (compute if not already stored)
           */

          int nData = 0;
          PDM_edges_intersect_res_t **_eir =
                  PDM_edges_intersect_get (ei,
                                           PDM_EDGES_GET_FROM_AB,
                                           vtx_currA->gNEdge,
                                           vtx_currB->gNEdge,
                                           &nData);
          PDM_edges_intersect_res_t *eir = *_eir;
          /*
           * Get new points and modified vertices coming from intersections
           */

          PDM_line_intersect_t         tIntersect;

          PDM_g_num_t                   nGEdgeA;
          PDM_g_num_t                   originEdgeA;
          int                          nNewPointsA;
          PDM_edges_intersect_point_t *oNewPointsA;
          PDM_g_num_t                  *linkA;
          PDM_g_num_t                  *gNumVtxA;
          double                      *coordsA;
          double                      *uA;

          /*
           * Get intersection properties
           */

          PDM_edges_intersect_res_data_get (eir,
                                            PDM_EDGES_INTERSECT_MESHA,
                                            &nGEdgeA,
																						&originEdgeA,
                                            &tIntersect,
                                            &nNewPointsA,
                                            &oNewPointsA,
                                            &linkA,
                                            &gNumVtxA,
                                            &coordsA,
                                            &uA);

          PDM_g_num_t                   nGEdgeB;
          PDM_g_num_t                   originEdgeB;
          int                          nNewPointsB;
          PDM_edges_intersect_point_t *oNewPointsB;
          PDM_g_num_t                  *linkB;
          PDM_g_num_t                  *gNumVtxB;
          double                      *coordsB;
          double                      *uB;

          PDM_edges_intersect_res_data_get (eir,
                                            PDM_EDGES_INTERSECT_MESHB,
                                            &nGEdgeB,
																						&originEdgeB,
                                            &tIntersect,
                                            &nNewPointsB,
                                            &oNewPointsB,
                                            &linkB,
                                            &gNumVtxB,
                                            &coordsB,
                                            &uB);

          free (eir);

          /*
           * Add new intersections vertex or switch vertex to intersection
           * if vertex is on subject polygon
           */

          if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {

            if (nNewPointsA == 2) {

              int iA = ((linkA[0] != vtx_currB->gN) ? 0 : 1);

              if (oNewPointsA[iA] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {

                _vertex_poly_t *interVtxA = _poly_clipp_intersect_add (uA[iA],
                                                                       gNumVtxA[iA],
                                                                       POLY_CLIPP_LOC_NEXT,
                                                                       vtx_currA);

                _poly_clipp_link_neighbors (interVtxA, vtx_currB);

              }

              else if (oNewPointsA[iA] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
                if (   ((originEdgeA == vtx_currA->gN) && (uA[iA] <=  0.1))
                    || ((originEdgeA != vtx_currA->gN) && ((1 - uA[iA]) <=  0.1))) {

                  _poly_clipp_link_neighbors (vtx_currA, vtx_currB);

                }
              }
            }

            else if (nNewPointsA == 1) {

              if (linkA[0] == vtx_currB->gN) {

                if (oNewPointsA[0] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
                  if (   ((originEdgeA == vtx_currA->gN) && (uA[0] <=  0.1))
                      || ((originEdgeA != vtx_currA->gN) && ((1 - uA[0]) <=  0.1))) {

                    _poly_clipp_link_neighbors (vtx_currA, vtx_currB);

                  }
                }

                else if (oNewPointsA[0] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
                  _vertex_poly_t *interVtxA = _poly_clipp_intersect_add (uA[0],
                                                                         gNumVtxA[0],
                                                                         POLY_CLIPP_LOC_NEXT,
                                                                         vtx_currA);

                  _poly_clipp_link_neighbors (interVtxA, vtx_currB);

                }
              }
            }

            if (nNewPointsB == 2) {

              int iB = ((linkB[0] != vtx_currA->gN) ? 0 : 1);

              if (oNewPointsB[iB] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {

                _vertex_poly_t *interVtxB = _poly_clipp_intersect_add (uB[iB],
                                                                       gNumVtxB[iB],
                                                                       POLY_CLIPP_LOC_NEXT,
                                                                       vtx_currB);

                _poly_clipp_link_neighbors (interVtxB, vtx_currA);

              }

              else if (oNewPointsB[iB] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
                if (   ((originEdgeB == vtx_currB->gN) && (uB[iB] <=  0.1))
                    || ((originEdgeB != vtx_currB->gN) && ((1 - uB[iB]) <=  0.1))) {

                  _poly_clipp_link_neighbors (vtx_currB, vtx_currA);

                }
              }
            }

            else if (nNewPointsB == 1) {

              if (linkB[0] == vtx_currA->gN) {

                if (oNewPointsB[0] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
                  if (   ((originEdgeB == vtx_currB->gN) && (uB[0] <=  0.1))
                      || ((originEdgeB != vtx_currB->gN) && ((1 - uB[0]) <=  0.1))) {

                    _poly_clipp_link_neighbors (vtx_currB, vtx_currA);

                  }
                }

                else if (oNewPointsB[0] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
                  _vertex_poly_t *interVtxB = _poly_clipp_intersect_add (uB[0],
                                                                         gNumVtxB[0],
                                                                         POLY_CLIPP_LOC_NEXT,
                                                                         vtx_currB);

                  _poly_clipp_link_neighbors (interVtxB, vtx_currA);

                }
              }
            }
          }

          else if (tIntersect == PDM_LINE_INTERSECT_YES) {

            /*
             * Look for if A vertex is on B edge
             */

            if (oNewPointsB[0] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {

              _vertex_poly_t *_currA;
              if (linkB[0] == vtx_currA->gN) {
                _currA = vtx_currA;

                if (_currA->neighbor == NULL) {

                  _vertex_poly_t *_currB = vtx_currB;
                  _vertex_poly_t *_next_currB = _currB->next;
                  while (_next_currB->isect) {
                    if (uB[0] < _next_currB->u) {
                      break;
                    }
                    _currB = _next_currB;
                    _next_currB = _currB->next;
                  };

                  _vertex_poly_t *interVtxB = _poly_clipp_intersect_add (uB[0],
                                                                        gNumVtxB[0],
                                                                        POLY_CLIPP_LOC_NEXT,
                                                                        _currB);


                  _poly_clipp_link_neighbors (interVtxB, _currA);
                }

                else {

                  _vertex_poly_t *firstInterVtxB = _currA->neighbor;

                  _vertex_poly_t *_currB = firstInterVtxB;

                  while (_currB->isect) {
                    _currB = _currB->previous;
                  }

                  if ((vtx_currB != _currB) && (nextB != _currB)) {

                    PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : Internal error 3\n"
                                    " B vertex is on 2 different A edges anisotropy ?\n"
                                    " ptB Origin1 origin2 : "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM" "PDM_FMT_G_NUM")\n",
                                    _currA->gN, _currB->gN, vtx_currB->gN);
                    abort();

                  }
                }
              }
            }

            /*
             * Look for if B vertex is on A edge
             */

            if (oNewPointsA[0] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
              _vertex_poly_t *_currB;
              if (linkA[0] == vtx_currB->gN) {
                _currB = vtx_currB;

                if (_currB->neighbor == NULL) {

                  _vertex_poly_t *_currA = vtx_currA;
                  _vertex_poly_t *_next_currA = _currA->next;
                  while (_next_currA->isect) {
                    if (uA[0] < _next_currA->u) {
                      break;
                    }
                    _currA = _next_currA;
                    _next_currA = _currA->next;
                  };

                  _vertex_poly_t *interVtxA = _poly_clipp_intersect_add (uA[0],
                                                                         gNumVtxA[0],
                                                                         POLY_CLIPP_LOC_NEXT,
                                                                         _currA);


                  _poly_clipp_link_neighbors (interVtxA, _currB);
                }

                else {

                  _vertex_poly_t *firstInterVtxA = _currB->neighbor;

                  _vertex_poly_t *_currA = firstInterVtxA;

                  while (_currA->isect) {
                    _currA = _currA->previous;
                  }

                  if ((vtx_currA != _currA) && (nextA != _currA)) {

                    PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : Internal error 3\n"
                                    " B vertex is on 2 different A edges anisotropy ?\n"
                                    " ptB Origin1 origin2 : "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM" "PDM_FMT_G_NUM")\n",
                                    _currB->gN, _currA->gN, vtx_currA->gN);
                    abort();

                  }
                }
              }
            }

            /*
             * Look for if A and B are the same vertex (link vertices if same)
             */

            if ((oNewPointsA[0] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) &&
                (oNewPointsB[0] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {

              _vertex_poly_t *_currA;
              _vertex_poly_t *_currB;
              if ((linkB[0] == vtx_currA->gN) && (linkA[0] == vtx_currB->gN)) {
                _currA = vtx_currA;
                _currB = vtx_currB;

                if ((_currA->neighbor != NULL) && (_currA->neighbor != _currB)) {

                  PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : Internal error 3\n"
                                  " B vertex is on 2 different A edges anisotropy ?\n"
                                  " ptA ptB1 ptB2 : "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM" "PDM_FMT_G_NUM")\n",
                                  _currA->gN, _currB->gN, _currA->neighbor->gN);
                  abort();

                }

                _poly_clipp_link_neighbors (_currA, _currB);

              }
            }

            /*
             * General case : insert a new intersection in linked lists
             */

            if ((oNewPointsA[0] == PDM_EDGES_INTERSECT_POINT_NEW) &&
                (oNewPointsB[0] == PDM_EDGES_INTERSECT_POINT_NEW)) {

              _vertex_poly_t *_currA = vtx_currA;
              _vertex_poly_t *_next_currA = _currA->next;
              while (_next_currA->isect) {
                if (uA[0] < _next_currA->u) {
                  break;
                }
                _currA = _next_currA;
                _next_currA = _currA->next;
              };

              _vertex_poly_t *interVtxA = _poly_clipp_intersect_add (uA[0],
                                                                     gNumVtxA[0],
                                                                     POLY_CLIPP_LOC_NEXT,
                                                                     _currA);

              _vertex_poly_t *_currB = vtx_currB;
              _vertex_poly_t *_next_currB = _currB->next;
              while (_next_currB->isect) {
                if (uB[0] < _next_currB->u) {
                  break;
                }
                _currB = _next_currB;
                _next_currB = _currB->next;
              };

              _vertex_poly_t *interVtxB = _poly_clipp_intersect_add (uB[0],
                                                                     gNumVtxB[0],
                                                                     POLY_CLIPP_LOC_NEXT,
                                                                     _currB);

              _poly_clipp_link_neighbors (interVtxA, interVtxB);

            }

          }

        }

        vtx_currB = vtx_currB->next;

      } while (vtx_currB != vtxB);

    }

    vtx_currA = nextA;

  } while (vtx_currA != vtxA);

  /*
   * Perform vertex location (in, out) for no 'on' vertex
   *    - First step  : In or out (ray tracing algorithm)
   *
   */

  _location (vtxA, vtxB, _faceVtxCooA, _faceVtxCooB, nVtxA, nVtxB, nA, nB);

  /*
   * Tag Intersection points
   *
   */

  _tag (vtxA, vtxB);

  /*
   * Clipping
   */

  int nPolyPredicA = 2;
  int sPolyConnecA = nVtxA + nVtxB;

  *nPolyClippA = 0;
  *polyClippConnecA = NULL;
  *polyClippIdxA = NULL;
  int idx1 = 0;
  int idx2 = 0;

  int nPolyPredicB = 2;
  int sPolyConnecB = nVtxA + nVtxB;

  *nPolyClippB = 0;
  *polyClippConnecB = NULL;
  *polyClippIdxB = NULL;

  int sPolyCoordA = 3 * sPolyConnecA;
  int sPolyCoordB = 3 * sPolyConnecB;

  if (performed_t == PDM_POLY_CLIPP_CLIP) {

    *polyClippIdxA = malloc (sizeof(int) * (nPolyPredicA + 1));
    (*polyClippIdxA)[0] = 0;

    *polyClippIdxB = *polyClippIdxA;

    *polyClippConnecA = malloc (sizeof(PDM_g_num_t) * sPolyConnecA);
    *polyClippConnecB = malloc (sizeof(PDM_g_num_t) * sPolyConnecB);

    *polyClippCoordsA = malloc (sizeof(double) * sPolyCoordA);
    *polyClippCoordsB = *polyClippCoordsA;

    sPolyCoordB = sPolyCoordA;

    /* Look for First intersection */

    _vertex_poly_t *first = vtxA;
    while (!first->isect && first->next != vtxA) {
      first = first->next;
    }

    /* If intersection found, clip  */

    if (first->isect) {

      int s_clipped_vtx = nVtxA + nVtxB;

      _vertex_poly_t **clipped_vtx = malloc (sizeof(_vertex_poly_t*) * s_clipped_vtx);

      _vertex_poly_t *curr = first;

      do {

        /* Clip */

        int nVtxClipp = 0;

        _vertex_poly_t *first_poly = curr;


        if (!curr->used) {

          clipped_vtx[nVtxClipp++] = first_poly;
          first_poly->used = true;

          bool direction = curr->tag;

          do {

            do {
              curr = direction ? curr->next : curr->previous;
              curr->used = true;

              if (nVtxClipp >= s_clipped_vtx) {
                while (nVtxClipp >= s_clipped_vtx) {
                  s_clipped_vtx *= 2;
                }
                clipped_vtx =
                realloc (clipped_vtx, sizeof(_vertex_poly_t*) * s_clipped_vtx);
              }

              clipped_vtx[nVtxClipp++] = curr;
              if (curr->neighbor != NULL) {
                bool isBreak = true;
                if (!curr->isect) {
                  if (!(curr->neighbor->tag && curr->tag)) {
                    isBreak = false;
                  }
                }
                if (isBreak) {
                  break;
                }
              }
            } while (true);

            /*
             * Change direction only for intersection vertex for not taken into
             * account contact area
             */

            if (curr->isect) {
              direction = curr->neighbor->tag;
            }

            curr = curr->neighbor;

          } while (curr != first_poly);


          /* Check if the polygon intersect at a point */

          if (nVtxClipp > 1) {
            nVtxClipp += -1;
          }

          /* Edges not taken into account */

          if (nVtxClipp > 2) {

            if (*nPolyClippA >= nPolyPredicA) {
              while (*nPolyClippA >= nPolyPredicA) {
                nPolyPredicA *= 2;
              }
              *polyClippIdxA = realloc (*polyClippIdxA, sizeof(int) * (nPolyPredicA + 1));
              *polyClippIdxB = *polyClippIdxA;
            }

            *nPolyClippA += 1;
            *nPolyClippB = *nPolyClippA;
            *polyClippIdxA[*nPolyClippA] =
                    *polyClippIdxA[*nPolyClippA - 1] + nVtxClipp;

            if (*polyClippIdxA[*nPolyClippA] >= sPolyConnecA) {
              while (*polyClippIdxA[*nPolyClippA] >= sPolyConnecA) {
                sPolyConnecA *= 2;
              }
              *polyClippConnecA = realloc (*polyClippConnecA, sizeof(PDM_g_num_t) *sPolyConnecA);
              sPolyConnecB = sPolyConnecA;
              *polyClippConnecB = realloc (*polyClippConnecB, sizeof(PDM_g_num_t) *sPolyConnecB);
            }

            if ((3 * (*polyClippIdxA[*nPolyClippA])) >= sPolyCoordA) {
              while ((3 * (*polyClippIdxA[*nPolyClippA])) >= sPolyCoordA) {
                sPolyCoordA *= 2;
              }
              *polyClippCoordsA = realloc (*polyClippCoordsA, sizeof(double) *sPolyCoordA);
              sPolyCoordB = sPolyCoordA;
              *polyClippCoordsB = *polyClippCoordsA;
            }

            for (int i = 0; i < nVtxClipp; i++) {
              if (clipped_vtx[i]->first == vtxA) {
                *polyClippConnecA[idx1  ] = clipped_vtx[i]->gN;
                *polyClippConnecB[idx1++] = clipped_vtx[i]->neighbor->gN;
              }
              else {
                *polyClippConnecB[idx1  ] = clipped_vtx[i]->gN;
                *polyClippConnecA[idx1++] = clipped_vtx[i]->neighbor->gN;
              }
              for (int j = 0; j < 3; j++) {
                *polyClippCoordsA[idx2++] = clipped_vtx[i]->coords[j];
              }
            }
          }
        }

        /* Look for next polygon */

        do {
          curr = curr->next;
        } while (!curr->isect);

      } while (curr != first);

      free (clipped_vtx);

    }

    else {

      /*
       * Check if the constraint polygon is inside the subject polygon
       */

      bool inside = true;

      _vertex_poly_t *curr = first;

      do {
        if (!curr->tag) {
          inside = false;
          break;
        }

        *polyClippConnecA[idx1  ] = curr->gN;
        *polyClippConnecB[idx1++] = curr->neighbor->gN;

        curr = curr->next;

      } while (curr != first);


      if (inside) {
        *nPolyClippA = 1;
        *nPolyClippB = 1;
        (*polyClippIdxA)[1] = nVtxA;
      }

      /*
       * Check if the subject polygon is inside the constraint polygon
       */

      else {

        curr = vtxB;

        idx1 = 0;

        do {
          if (!curr->tag) {
            inside = false;
            break;
          }

          *polyClippConnecB[idx1  ] = curr->gN;
          *polyClippConnecA[idx1++] = curr->neighbor->gN;

          curr = curr->next;
        } while (curr != vtxB);


        if (inside) {
          *nPolyClippA = 1;
          *nPolyClippB = 1;
          (*polyClippIdxA)[1] = nVtxB;
        }

      }

    }

    /*
     * Update size
     */

    *polyClippIdxA =
            realloc (*polyClippIdxA, (sizeof(int) * (*nPolyClippA + 1)));

    *polyClippConnecA =
            realloc (*polyClippConnecA, (sizeof(PDM_g_num_t) * (*polyClippIdxA)[*nPolyClippA]));

    *polyClippCoordsA =
            realloc (*polyClippCoordsA, (sizeof(double) * 3 * (*polyClippIdxA)[*nPolyClippA]));

    *polyClippIdxB = *polyClippIdxA;
    *polyClippConnecB =
            realloc (*polyClippConnecB, (sizeof(PDM_g_num_t) * (*polyClippIdxB)[*nPolyClippB]));

  }

  /*
   * Reverse mod
   */

  else if (performed_t == PDM_POLY_CLIPP_REVERSE) {

    *polyClippIdxA = malloc (sizeof(int) * (nPolyPredicA + 1));
    (*polyClippIdxA)[0] = 0;

    *polyClippIdxB = malloc (sizeof(int) * (nPolyPredicB + 1));
    (*polyClippIdxB)[0] = 0;

    *polyClippCoordsA = malloc (sizeof(double) * sPolyCoordA);
    *polyClippCoordsB = malloc (sizeof(double) * sPolyCoordB);

    /*
     * A Reverse clipping
     */

    /* Look for first out vertex*/

    _vertex_poly_t *first_out = vtxA;
    while (first_out->isect || first_out->tag) {
      first_out = first_out->next;
    }

    /* Fist polygon */

    idx1 = 0;

    if (!first_out->isect && !first_out->tag) {

      int s_clipped_vtx = nVtxA + nVtxB;

      _vertex_poly_t **clipped_vtx = malloc (sizeof(_vertex_poly_t*) * s_clipped_vtx);

      _vertex_poly_t *curr = first_out;

      do {

        /* Clip */

        int nVtxClipp = 0;

        _vertex_poly_t *first_poly_out = curr;

        if (!curr->used) {

          curr->used = true;

          bool direction = true;

          do {

            if (curr->neighbor != NULL) {

              if (curr->isect) {

                direction = !direction;
                curr = curr->neighbor;

              }

              else {

                if ((direction && curr->neighbor->tag) ||
                    (!direction && !curr->neighbor->tag)) {
                  direction = !direction;
                  curr = curr->neighbor;
                }

              }

            }

            do {

              if (nVtxClipp >= s_clipped_vtx) {
                while (nVtxClipp >= s_clipped_vtx) {
                  s_clipped_vtx *= 2;
                }
                clipped_vtx =
                realloc (clipped_vtx, sizeof(_vertex_poly_t*) * s_clipped_vtx);
              }

              clipped_vtx[nVtxClipp++] = curr;
              curr->used = true;

              curr = direction ? curr->next : curr->previous;

            } while ((curr != first_poly_out) && (curr->neighbor == NULL));

          } while (curr != first_poly_out);

          if (nVtxClipp > 2) {

            if (*nPolyClippA >= nPolyPredicA) {
              while (*nPolyClippA >= nPolyPredicA) {
                nPolyPredicA *= 2;
              }
              *polyClippIdxA = realloc (*polyClippIdxA, sizeof(int) * (nPolyPredicA + 1));
            }

            *nPolyClippA += 1;
            *polyClippIdxA[*nPolyClippA] =
                    *polyClippIdxA[*nPolyClippA - 1] + nVtxClipp;

            if (*polyClippIdxA[*nPolyClippA] >= sPolyConnecA) {
              while (*polyClippIdxA[*nPolyClippA] >= sPolyConnecA) {
                sPolyConnecA *= 2;
              }
              *polyClippConnecA = realloc (*polyClippConnecA, sizeof(PDM_g_num_t) *sPolyConnecA);
            }

            if ((3 * (*polyClippIdxA[*nPolyClippA])) >= sPolyCoordA) {
              while ((3 * (*polyClippIdxA[*nPolyClippA])) >= sPolyCoordA) {
                sPolyCoordA *= 2;
              }
              *polyClippCoordsA = realloc (*polyClippCoordsA, sizeof(double) *sPolyCoordA);
            }

            for (int i = 0; i < nVtxClipp; i++) {
              if (clipped_vtx[i]->first == vtxA) {
                *polyClippConnecA[idx1++] = clipped_vtx[i]->gN;
              }
              else {
                *polyClippConnecA[idx1++] = clipped_vtx[i]->neighbor->gN;
              }
              for (int j = 0; j < 3; j++) {
                *polyClippCoordsA[idx2++] = clipped_vtx[i]->coords[j];
              }
            }
          }
        }

        do {
          curr = curr->next;
        } while (curr->isect && (curr != first_out));

      } while (curr != first_out);

      free (clipped_vtx);

    }

    /*
     * B Reverse clipping
     */

    /* Look for first out vertex*/

    first_out = vtxB;
    while (first_out->isect || first_out->tag) {
      first_out = first_out->next;
    }

    idx1 = 0;
    idx2 = 0;

    /* Fist polygon */

    if (!first_out->isect && !first_out->tag) {

      int s_clipped_vtx = nVtxA + nVtxB;

      _vertex_poly_t **clipped_vtx = malloc (sizeof(_vertex_poly_t*) * s_clipped_vtx);

      _vertex_poly_t *curr = first_out;

      do {

        /* Clip */

        int nVtxClipp = 0;

        _vertex_poly_t *first_poly_out = curr;

        if (!curr->used) {

          curr->used = true;

          bool direction = true;

          do {

            if (curr->neighbor != NULL) {

              if (curr->isect) {
                direction = !direction;
                curr = curr->neighbor;
              }

              else {
                if ((direction && curr->neighbor->tag) ||
                    (!direction && !curr->neighbor->tag)) {
                  direction = !direction;
                  curr = curr->neighbor;
                }
              }
            }

            do {

              if (nVtxClipp >= s_clipped_vtx) {
                while (nVtxClipp >= s_clipped_vtx) {
                  s_clipped_vtx *= 2;
                }
                clipped_vtx =
                realloc (clipped_vtx, sizeof(_vertex_poly_t*) * s_clipped_vtx);
              }

              clipped_vtx[nVtxClipp++] = curr;
              curr->used = true;

              curr = direction ? curr->next : curr->previous;

            } while ((curr != first_poly_out) && (curr->neighbor == NULL));

          } while (curr != first_poly_out);

          if (nVtxClipp > 2) {

            if (*nPolyClippB >= nPolyPredicB) {
              while (*nPolyClippB >= nPolyPredicB) {
                nPolyPredicB *= 2;
              }
              *polyClippIdxB = realloc (*polyClippIdxB, sizeof(int) * (nPolyPredicB + 1));
            }

            *nPolyClippB += 1;
            *polyClippIdxB[*nPolyClippB] =
                    *polyClippIdxB[*nPolyClippB - 1] + nVtxClipp;

            if (*polyClippIdxB[*nPolyClippB] >= sPolyConnecB) {
              while (*polyClippIdxB[*nPolyClippB] >= sPolyConnecB) {
                sPolyConnecB *= 2;
              }
              *polyClippConnecB = realloc (*polyClippConnecB, sizeof(PDM_g_num_t) *sPolyConnecB);
            }

            if ((3 * (*polyClippIdxB[*nPolyClippB])) >= sPolyCoordB) {
              while ((3 * (*polyClippIdxB[*nPolyClippB])) >= sPolyCoordB) {
                sPolyCoordB *= 2;
              }
              *polyClippCoordsB = realloc (*polyClippCoordsB, sizeof(double) *sPolyCoordB);
            }

            for (int i = 0; i < nVtxClipp; i++) {
              if (clipped_vtx[i]->first == vtxB) {
                *polyClippConnecB[idx1++] = clipped_vtx[i]->gN;
              }
              else {
                *polyClippConnecB[idx1++] = clipped_vtx[i]->neighbor->gN;
              }
              for (int j = 0; j < 3; j++) {
                *polyClippCoordsB[idx2++] = clipped_vtx[i]->coords[j];
              }
            }
          }
        }

        do {
          curr = curr->next;
        } while (curr->isect && (curr != first_out));

      } while (curr != first_out);

      free (clipped_vtx);
    }

    /*
     * Update size
     */

    *polyClippIdxA =
            realloc (*polyClippIdxA, (sizeof(int) * (*nPolyClippA + 1)));
    *polyClippConnecA =
            realloc (*polyClippConnecA, (sizeof(PDM_g_num_t) * (*polyClippIdxA)[*nPolyClippA]));
    *polyClippCoordsA =
            realloc (*polyClippCoordsA, (sizeof(double)* 3 * (*polyClippIdxA)[*nPolyClippA]));

    *polyClippIdxB =
            realloc (*polyClippIdxB, (sizeof(int) * (*nPolyClippB + 1)));
    *polyClippConnecB =
            realloc (*polyClippConnecB, (sizeof(PDM_g_num_t) * (*polyClippIdxB)[*nPolyClippB]));

    *polyClippCoordsB =
            realloc (*polyClippCoordsB, (sizeof(double)* 3 * (*polyClippIdxB)[*nPolyClippB]));
  }


  /*
   * Free local memory
   */

  if (_faceToEdgeA != faceToEdgeA) {
    free (_faceToEdgeA);
  }

  if (_faceToVtxA  != faceToVtxA) {
    free (_faceToVtxA);
  }

  if (_faceVtxCooA != faceVtxCooA) {
    free (_faceVtxCooA);
  }

  if (_faceToEdgeB != faceToEdgeB) {
    free (_faceToEdgeB);
  }

  if (_faceToVtxB != faceToVtxB) {
    free (_faceToVtxB);
  }

  if (_faceVtxCooB != faceVtxCooB) {
    free (_faceVtxCooB);
  }

}

#ifdef	__cplusplus
}
#endif


