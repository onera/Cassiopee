#ifndef __PDM_GEOM_ELEM_H__
#define __PDM_GEOM_ELEM_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2012  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *  \brief Compute a dynamic geometric epsilon from a characteristic length
 *
 *    @param [in]  characteristicLength  Characteristic length
 *    @param [in]  consEpsilon           Constant part
 *    @return                            Geometric epsilon
 */

double
PDM_geom_elem_geometric_epsilon
(
 const double characteristicLength,
 const double constEpsilon
 );

/**
 *  \brief Triangle surface vector
 *
 *  @param [in]  nTriangle      Number of triangles
 *  @param [in]  connectivity   Connectivity
 *  @param [in]  coords         Vertice coordinates
 *  @param [out] surfaceVector  Surface Vector
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tria_surface_vector
(
 const int     nTriangle,
 const int    *connectivity,
 const double *coords,
 double       *surfaceVector,
 double       *characteristicLength,
 int          *isDegenerated
 );

/**
 *  \brief Triangle area
 *
 *  @param [in]  nTriangle      Number of triangles
 *  @param [in]  surfaceVector         SurfaceVector vectors
 *  @param [out] area           Area
 */

void
PDM_geom_elem_tria_area
(
 const int     nTriangle,
 const double *surfaceVector,
 double *area
);


/**
 *  \brief Triangle center
 *
 *  @param [in]  nTriangle      Number of triangles
 *  @param [in]  connectivity   Connectivity
 *  @param [in]  coords         Vertice coordinates
 *  @param [out] center         center
 */

void
PDM_geom_elem_tria_center
(
 const int     nTriangle,
 const int    *connectivity,
 const double *coords,
 double *center
 );


/**
 *  \brief Tetrahedra oriented volume
 *
 *  @param [in]  nTetrahedra           Number of tetrahedra
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  coords                Vertice coordinates
 *  @param [out] volume                Volume
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tetra_oriented_volume
(
 const int     nTetrahedra,
 const int    *connectivity,
 const double *coords,
 double       *volume,
 double       *characteristicLength,
 int         *isDegenerated
 );

/**
 *  \brief Tetrahedra center
 *
 *  @param [in]  nTetrahedra    Number of tetrahedra
 *  @param [in]  connectivity   Connectivity
 *  @param [in]  coords         Vertice coordinates
 *  @param [out] center         center
 */

void
PDM_geom_elem_tetra_center
(
 const int     nTetrahedra,
 const int    *connectivity,
 const double *coords,
 double *center
 );

/**
 *  \brief Tetrahedra Faces
 *
 *  @param [in]  nTetrahedra       Number of tetrahedra
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] faceConnectivity  Face connectivity
 */

void
PDM_geom_elem_tetra_faces
(
 const int     nTetrahedra,
 const int     orientation,
 const int    *connectivity,
 int          *faceConnectivityIndex,
 int          *faceConnectivity
);


/**
 *  \brief HexahedraFaces
 *
 *  @param [in]  nHexahedra        Number of hexahedra
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] faceConnectivity  Face connectivity
 */

void
PDM_geom_elem_hexa_faces
(
 const int     nHexahedra,
 const int     orientation,
 const int    *connectivity,
 int          *faceConnectivityIndex,
 int          *faceConnectivity
 );


/**
 *  \brief Prism Faces
 *
 *  @param [in]  nPrism            Number of Prism
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] faceConnectivity  Face connectivity
 */

void
PDM_geom_elem_prism_faces
(
 const int     nPrism,
 const int     orientation,
 const int    *connectivity,
 int          *faceConnectivityIndex,
 int          *faceConnectivity
 );

/**
 *  \brief Pyramid Faces
 *
 *  @param [in]  nPyramid          Number of pyramid
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] faceConnectivity  Face connectivity
 */

void
PDM_geom_elem_pyramid_faces
(
 const int   nPyramid,
 const int   orientation,
 const int  *connectivity,
 int        *faceConnectivityIndex,
 int        *faceConnectivity
 );


/**
 *  \brief Edges properties
 *
 *  @param [in]  nEdges                Number of edges
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] center                Center
 *  @param [out] length                Length
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_edges_properties
(
 const int     nEdges,
 const int    *connectivity,
 const double *coords,
 double       *center,
 double       *length,
 double       *characteristicLength,
 int         *isDegenerated
 );


/**
 *  \brief Triangle properties
 *
 *  @param [in]  nTriangle             Number of triangles
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surfaceVector         Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tria_properties
(
 const int     nTriangle,
 const int    *connectivity,
 const double *coords,
 double       *surfaceVector,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
 );


/**
 * \brief Quadrangle properties
 *
 *  @param [in]  nTriangle             Number of quadrangles
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surfaceVector         Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 *
 *  @return                     The status of properties computation convergence
 */

int
PDM_geom_elem_quad_properties
(
 const int     nQuadrangle,
 const int    *connectivity,
 const double *coords,
 double       *surfaceVector,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
);


/**
 * \brief Compute the barycentric coordinates of a set of points inside
          their belonging polygons.
 *
 *  @param [in]  nPoints               Number of points
 *  @param [in]  ptsLocations          Numbering of the belonging polygons inside the connectivityIndex
 *  @param [in]  connectivityIndex     Mesh connectivity Index
 *  @param [in]  connectivity          Mesh connectivity
 *  @param [in]  coords                Mesh coordinates
 *  @param [out] barCoordsIndex        Pointer to the barycentric coordinates index
 *  @param [out] barCoordsIndex        Pointer to the barycentric coordinates
 *
 *  @return                     The status of properties computation convergence
 */

int
PDM_geom_elem_compute_polygon_barycentric_coordinates
(
 const int           n_points,
 const int          *pts_locations,
 const double       *pts_coords,
 const int          *connectivityIndex,
 const int          *connectivity,
 const double       *coords,
 int               **barCoordsIndex,
 double            **barCoords
);

/**
 *  \brief Polygon properties
 *
 *  @param [in]  nPolygon              Number of polygon
 *  @param [in]  connectivityIndex     Connectivity Index
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surfaceVector         Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 *
 *  @return                        The status of properties computation convergence
 */

int
PDM_geom_elem_polygon_properties
(
 const int     nPolygon,
 const int    *connectivityIndex,
 const int    *connectivity,
 const double *coords,
 double       *surfaceVector,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
);


/**
 *  \brief Tetrahedra properties
 *
 *  @param [in]  nTetrahedra           Number of tetrahedra
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tetra_properties
(
 const int     nTetrahedra,
 const int    *connectivity,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
 );


/**
 *  \brief Hexahedra properties
 *
 *  @param [in]  nHexahedra            Number of hexahedra
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_hexa_properties
(
 const int     nHexahedra,
 const int    *connectivity,
 const int     nVertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
);


/**
 *  \brief Prism properties
 *
 *  @param [in]  nPrism                Number of prism
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_prism_properties
(
 const int     nPrism,
 const int    *connectivity,
 const int     nVertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
);


/**
 *  \brief Pyramid properties
 *
 *  @param [in]  nPyramid              Number of pyramid
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  nVertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristicLength  Characteristic length (active if != NULL)
 *  @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_pyramid_properties
(
 const int     nPyramid,
 const int    *connectivity,
 const int     nVertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
);


/**
 *  \brief Polyhedra properties
 *
 *  @param [in]  nPolyhedra                 Number of polyhedra
 *  @param [in]  nFace                      Number of faces
 *  @param [in]  faceConnectivityIdx        Face connectivity index
 *  @param [in]  faceConnectivity           Face connectivity
 *  @param [in]  cellToFaceConnectivityIdx  Cell to face connectivity index
 *  @param [in]  cellToFaceConnectivity     Cell to face connectivity
 *  @param [in]  nVertices                  Number of vertices
 *  @param [in]  coords                     Vertices coordinates
 *  @param [out] volume                     Volume
 *  @param [out] center                     Center
 *  @param [out] characteristicLength       Characteristic length (active if != NULL)
 *  @param [out] isDegenerated              Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_polyhedra_properties
(
 const  int    isOriented,
 const int     nPolyhedra,
 const int     nFace,
 const int    *faceConnectivityIdx,
 const int    *faceConnectivity,
 const int    *cellToFaceConnectivityIdx,
       int    *cellToFaceConnectivity,
 const int     nVertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristicLength,
 int         *isDegenerated
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ELEM_GEOM__ */
