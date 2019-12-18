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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_polygon.h"
#include "pdm_plane.h"
#include "pdm_geom_elem.h"
#include "pdm_hash_tab.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

enum {false, true};

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double GEOM_EPS_MIN  = 1e-30; /*!< Minimum value allowed for geometric computation */

static const double GEOM_EPS_VOL  = 1e-9; /*!< Constant value used to compute geomtric epsilon for volume */

static const double GEOM_EPS_SURF = 1e-9; /*!< Constant value used to compute geomtric epsilon for surface */

static const double GEOM_EPS_DIST = 1e-9; /*!< Minimum distance between two vertices */

enum {
  NOT_DEFINE,
  IN_STACK,
  UNCHANGED_CYCLE,
  CHANGED_CYCLE,
};

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
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
 )
{
  return PDM_MAX(constEpsilon * characteristicLength, GEOM_EPS_MIN);
}

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
 )
{
  for (int itri = 0; itri < nTriangle; itri++) {

    const int i = connectivity[3*itri    ] - 1;
    const int j = connectivity[3*itri + 1] - 1;
    const int k = connectivity[3*itri + 2] - 1;

    double v1[3];
    double v2[3];
    double *surfaceVectorTri = surfaceVector + 3*itri;

    v1[0] = coords[3*j    ] - coords[3*i    ];
    v1[1] = coords[3*j + 1] - coords[3*i + 1];
    v1[2] = coords[3*j + 2] - coords[3*i + 2];

    v2[0] = coords[3*k    ] - coords[3*i    ];
    v2[1] = coords[3*k + 1] - coords[3*i + 1];
    v2[2] = coords[3*k + 2] - coords[3*i + 2];

    if (characteristicLength != NULL) {

      double v3[3];

      v3[0] = coords[3*k    ] - coords[3*j    ];
      v3[1] = coords[3*k + 1] - coords[3*j + 1];
      v3[2] = coords[3*k + 2] - coords[3*j + 2];

      double normV1 = PDM_MODULE(v1);
      double normV2 = PDM_MODULE(v2);
      double normV3 = PDM_MODULE(v3);

      characteristicLength[itri] = PDM_MIN(PDM_MIN(normV1, normV2), normV3);

    }

    PDM_CROSS_PRODUCT(v1, v2, surfaceVectorTri);

    surfaceVectorTri[0] *= 0.5;
    surfaceVectorTri[1] *= 0.5;
    surfaceVectorTri[2] *= 0.5;

    if ((characteristicLength != NULL) && (isDegenerated != NULL)) {

      double normSurfaceVectorTri = PDM_MODULE(surfaceVectorTri);
      double eps_loc = PDM_geom_elem_geometric_epsilon(characteristicLength[itri], GEOM_EPS_SURF);
      isDegenerated[itri] = 0;
      if (normSurfaceVectorTri <= eps_loc)
        isDegenerated[itri] = 1;

    }
  }
}

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
 )

{
  for (int itri = 0; itri < nTriangle; itri++) {
    const double *surfaceVectorTri = surfaceVector + 3*itri;
    area[itri] = PDM_MODULE(surfaceVectorTri);
  }
}

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
 )
{
  for (int itri = 0; itri < nTriangle; itri++) {

    const int i = connectivity[3*itri    ] - 1;
    const int j = connectivity[3*itri + 1] - 1;
    const int k = connectivity[3*itri + 2] - 1;

    double *centerTri = center + 3*itri;

    centerTri[0] = (coords[3*i    ] + coords[3*j    ] + coords[3*k    ])/3;
    centerTri[1] = (coords[3*i + 1] + coords[3*j + 1] + coords[3*k + 1])/3;
    centerTri[2] = (coords[3*i + 2] + coords[3*j + 2] + coords[3*k + 2])/3;
  }
}


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
 )
{
  for (int itet = 0; itet < nTetrahedra; itet++) {

    const int *connectivityTet = connectivity + 4*itet;

    const int i1 = connectivityTet[0] - 1;
    const int i2 = connectivityTet[1] - 1;
    const int i3 = connectivityTet[2] - 1;
    const int i4 = connectivityTet[3] - 1;

    const int nTriangle = 1;
    double surfaceVector[3];
    double faceCenter[3];
    double fC_i4[3];

    PDM_geom_elem_tria_surface_vector(nTriangle,
                                      connectivityTet,
                                      coords,
                                      surfaceVector,
                                      NULL,
                                      NULL);

    if (characteristicLength != NULL) {

      double vectV1V2[3] =
        {coords[3*i2    ] - coords[3*i1    ],
         coords[3*i2 + 1] - coords[3*i1 + 1],
         coords[3*i2 + 2] - coords[3*i1 + 2]};

      double vectV1V3[3] =
        {coords[3*i3    ] - coords[3*i1    ],
         coords[3*i3 + 1] - coords[3*i1 + 1],
         coords[3*i3 + 2] - coords[3*i1 + 2]};

      double vectV2V3[3] =
        {coords[3*i3    ] - coords[3*i2    ],
         coords[3*i3 + 1] - coords[3*i2 + 1],
         coords[3*i3 + 2] - coords[3*i2 + 2]};

      double vectV1V4[3] =
        {coords[3*i4    ] - coords[3*i1    ],
         coords[3*i4 + 1] - coords[3*i1 + 1],
         coords[3*i4 + 2] - coords[3*i1 + 2]};

      double normV1V2 = PDM_MODULE(vectV1V2);
      double normV1V3 = PDM_MODULE(vectV1V3);
      double normV2V3 = PDM_MODULE(vectV2V3);
      double normV1V4 = PDM_MODULE(vectV1V4);

      characteristicLength[itet] = PDM_MIN(normV1V2,                    normV1V3);
      characteristicLength[itet] = PDM_MIN(characteristicLength[itet],  normV2V3);
      characteristicLength[itet] = PDM_MIN(characteristicLength[itet],  normV1V2);
      characteristicLength[itet] = PDM_MIN(characteristicLength[itet],  normV1V3);
      characteristicLength[itet] = PDM_MIN(characteristicLength[itet],  normV1V4);

    }

    PDM_geom_elem_tria_center (nTriangle,
                               connectivity,
                               coords,
                               faceCenter);

    fC_i4[0] = coords[3*i4    ] - faceCenter[0];
    fC_i4[1] = coords[3*i4 + 1] - faceCenter[1];
    fC_i4[2] = coords[3*i4 + 2] - faceCenter[2];

    volume[itet] = 1./3. * PDM_DOT_PRODUCT(fC_i4, surfaceVector);

    if ((characteristicLength != NULL) && (isDegenerated != NULL)) {

      double eps_loc = PDM_geom_elem_geometric_epsilon(characteristicLength[itet], GEOM_EPS_VOL);
      isDegenerated[itet] = 0;
      if (volume[itet]  <= eps_loc)
        isDegenerated[itet] = 1;

    }
  }
}

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
 )

{
  for (int itet = 0; itet < nTetrahedra; itet++) {

    const int i1 = connectivity[4*itet    ] - 1;
    const int i2 = connectivity[4*itet + 1] - 1;
    const int i3 = connectivity[4*itet + 2] - 1;
    const int i4 = connectivity[4*itet + 3] - 1;

    double *centerTet = center + 3*itet;

    centerTet[0] = (coords[3*i1    ] + coords[3*i2    ]
                    + coords[3*i3    ] + coords[3*i4    ])/4;
    centerTet[1] = (coords[3*i1 + 1] + coords[3*i2 + 1]
                    + coords[3*i3 + 1] + coords[3*i4 + 1])/4;
    centerTet[2] = (coords[3*i1 + 2] + coords[3*i2 + 2]
                    + coords[3*i3 + 2] + coords[3*i4 + 2])/4;
  }
}


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
 )
{
  const int nFace       = 4;         /* 4 triangles */
  const int nVertex     = nFace * 3; /* 4 * 3 vertices */
  const int nVertexElt  = 4;         /* 4 vertices */

  faceConnectivityIndex[0] = 0;

  for (int ielt = 0; ielt < nTetrahedra; ielt++) {

    for (int iface = 0; iface < nFace; iface++)
      faceConnectivityIndex[ielt * nFace + iface + 1] =
        faceConnectivityIndex[ielt * nFace + iface] + 3;

    if (orientation == 0) {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];

      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 1];

      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 3];

      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 2];

    }

    else {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 1];

    }
  }
}


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
 )
{
  const int nFace       = 6;         /* 6 quadrangles */
  const int nVertex     = nFace * 4; /* 6 * 4 vertices */
  const int nVertexElt  = 8;         /* 8 vertices */

  faceConnectivityIndex[0] = 0;

  for (int ielt = 0; ielt < nHexahedra; ielt++) {

    for (int iface = 0; iface < nFace; iface++)
      faceConnectivityIndex[ielt * nFace + iface + 1] =
        faceConnectivityIndex[ielt * nFace + iface] + 4;

    if (orientation == 0) {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 3];

      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 7];
      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 6];

      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 7];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 4];

      faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 6];
      faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 7];

      faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 18] = connectivity[nVertexElt * ielt + 6];
      faceConnectivity[nVertex * ielt + 19] = connectivity[nVertexElt * ielt + 2];

      faceConnectivity[nVertex * ielt + 20] = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 21] = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 22] = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 23] = connectivity[nVertexElt * ielt + 1];

    }

    else {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 6];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 7];
      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 5];

      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 7];
      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 7];
      faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 6];
      faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 3];

      faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt + 6];
      faceConnectivity[nVertex * ielt + 18] = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 19] = connectivity[nVertexElt * ielt + 1];

      faceConnectivity[nVertex * ielt + 20] = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 21] = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 22] = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 23] = connectivity[nVertexElt * ielt + 0];

    }
  }
}


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
 )
{
  const int nFace       = 5;         /* 3 quadrangles + 2 triangles */
  const int nVertex     = 3*4 + 2*3; /* */
  const int nVertexElt  = 6;         /* 6 vertices */

  faceConnectivityIndex[0] = 0;

  for (int ielt = 0; ielt < nPrism; ielt++) {

    faceConnectivityIndex[ielt * nFace + 1] = faceConnectivityIndex[ielt * nFace    ] + 3;
    faceConnectivityIndex[ielt * nFace + 2] = faceConnectivityIndex[ielt * nFace + 1] + 3;
    faceConnectivityIndex[ielt * nFace + 3] = faceConnectivityIndex[ielt * nFace + 2] + 4;
    faceConnectivityIndex[ielt * nFace + 4] = faceConnectivityIndex[ielt * nFace + 3] + 4;
    faceConnectivityIndex[ielt * nFace + 5] = faceConnectivityIndex[ielt * nFace + 4] + 4;

    if (orientation == 0) {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];

      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 4];

      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 5];

      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 4];

      faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt + 3];

    }

    else {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 3];

      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 2];

      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 1];

      faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 5];
      faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt    ];

    }
  }
}

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
 )
{
  const int nFace       = 5;         /* 1 quadrangle + 4 triangles */
  const int nVertex     = 1*4 + 4*3; /* */
  const int nVertexElt  = 5;         /* 5 vertices */

  faceConnectivityIndex[0] = 0;

  for (int ielt = 0; ielt < nPyramid; ielt++) {

    faceConnectivityIndex[ielt * nFace + 1] = faceConnectivityIndex[ielt * nFace    ] + 4;
    faceConnectivityIndex[ielt * nFace + 2] = faceConnectivityIndex[ielt * nFace + 1] + 3;
    faceConnectivityIndex[ielt * nFace + 3] = faceConnectivityIndex[ielt * nFace + 2] + 3;
    faceConnectivityIndex[ielt * nFace + 4] = faceConnectivityIndex[ielt * nFace + 3] + 3;
    faceConnectivityIndex[ielt * nFace + 5] = faceConnectivityIndex[ielt * nFace + 4] + 3;

    if (orientation == 0) {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 3];

      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 4];

      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 4];

      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 4];

      faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 4];

    }

    else {

      faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt    ];

      faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt    ];
      faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 1];

      faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 1];
      faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 2];

      faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 2];
      faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 3];

      faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 4];
      faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 3];
      faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt    ];

    }
  }
}


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
 )
{

  for (int iedge = 0; iedge < nEdges; iedge++) {

    const int *connectivityEdge = connectivity + 2*iedge;

    const int i1 = connectivityEdge[0] - 1;
    const int i2 = connectivityEdge[1] - 1;

    double *centerEdge = center + 3*iedge;

    for (int i = 0; i < 3; i++)
      centerEdge[i] = 0.5 * (coords[3*i1 + i] + coords[3*i2 + i]);

    length[iedge] = sqrt((coords[3*i2 + 0] - coords[3*i1 + 0])
                         * (coords[3*i2 + 0] - coords[3*i1 + 0])
                         + (coords[3*i2 + 1] - coords[3*i1 + 1])
                         * (coords[3*i2 + 1] - coords[3*i1 + 1])
                         + (coords[3*i2 + 2] - coords[3*i1 + 2])
                         * (coords[3*i2 + 2] - coords[3*i1 + 2]));

    if (characteristicLength != NULL)
      characteristicLength[iedge] = length[iedge];

    if (isDegenerated != NULL) {
      isDegenerated[iedge] = 0;
      if (length[iedge] < GEOM_EPS_DIST)
        isDegenerated[iedge] = 1;
    }
  }
}


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
 )
{
  PDM_geom_elem_tria_surface_vector (nTriangle,
                                     connectivity,
                                     coords,
                                     surfaceVector,
                                     characteristicLength,
                                     isDegenerated);

  PDM_geom_elem_tria_center (nTriangle,
                             connectivity,
                             coords,
                             center);

}


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
 )

{
  PDM_geom_elem_tetra_oriented_volume (nTetrahedra,
                                       connectivity,
                                       coords,
                                       volume,
                                       characteristicLength,
                                       isDegenerated);

  PDM_geom_elem_tetra_center (nTetrahedra,
                              connectivity,
                              coords,
                              center);
}


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
)
{
  int *connectivityIndex = (int *) malloc (sizeof(int) * (nQuadrangle + 1));
  int convergence;

  connectivityIndex[0] = 0;
  for (int i = 1; i < nQuadrangle + 1; i++) {
    connectivityIndex[i] = connectivityIndex[i-1] + 4;
  }

  convergence = PDM_geom_elem_polygon_properties(nQuadrangle,
                                                 connectivityIndex,
                                                 connectivity,
                                                 coords,
                                                 surfaceVector,
                                                 center,
                                                 characteristicLength,
                                                 isDegenerated);

  free (connectivityIndex);

  return convergence;
}

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

int PDM_geom_elem_compute_polygon_barycentric_coordinates(const int           nPoints,
                                                           const int          *ptsLocations,
                                                           const double       *pts_coords,
                                                           const int          *connectivityIndex,
                                                           const int          *connectivity,
                                                           const double       *coords,
                                                           int               **barCoordsIndex,
                                                           double            **barCoords
)
{

  int convergence = 1;
  double local_pts[3];

  /* Tableaux locaux */

  const double eps_base = 1e-10;
  double* coords_sommets = NULL;
  double* s = NULL;
  double* dist = NULL;
  double* aire = NULL;
  double* proScal = NULL;

  *barCoordsIndex = (int*) malloc(sizeof(int) * (nPoints+1) );
  int* _barCoordsIndex = *barCoordsIndex;

  int prev_n_sommets = 0;
  int n_sommets = 0;

  _barCoordsIndex[0] = 0;
  /* Boucle sur les points distants */
  for (int ipoint =  0; ipoint < nPoints; ipoint++ ) {

    /* Initialisation - Copie locale */

    int isOnEdge = 0;
    int isVertex = 0;
    int ielt = ptsLocations[ipoint] - 1;
    prev_n_sommets = n_sommets;
    n_sommets =  connectivityIndex[ielt+1] -
                       connectivityIndex[ielt];

    local_pts[0] = pts_coords[3*ipoint];
    local_pts[1] = pts_coords[3*ipoint + 1];
    local_pts[2] = pts_coords[3*ipoint + 2];

    if (ipoint == 0) {
      coords_sommets = (double*)malloc(sizeof(double)* 3 * n_sommets);
      s              = (double*)malloc(sizeof(double)* 3 * n_sommets);
      dist           = (double*)malloc(sizeof(double)* n_sommets);
      aire           = (double*)malloc(sizeof(double)* n_sommets);
      proScal        = (double*)malloc(sizeof(double)* n_sommets);
    }
    else {
      if (prev_n_sommets < n_sommets) {
        coords_sommets = (double*)realloc(coords_sommets,sizeof(double)* 3 * n_sommets);
        s              = (double*)realloc(s,      sizeof(double)* 3 * n_sommets);
        dist           = (double*)realloc(dist,   sizeof(double)* n_sommets);
        aire           = (double*)realloc(aire,   sizeof(double)* n_sommets);
        proScal        = (double*)realloc(proScal,sizeof(double)* n_sommets);
      }
    }

    for (int isom = 0; isom < n_sommets; isom++) {
      coords_sommets[3*isom]   =
        coords[3*(connectivity[connectivityIndex[ielt]+isom]-1)];

      coords_sommets[3*isom+1] =
        coords[3*(connectivity[connectivityIndex[ielt]+isom]-1)+1];

      coords_sommets[3*isom+2] =
        coords[3*(connectivity[connectivityIndex[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    double bary[3];

    PDM_polygon_compute_barycenter (n_sommets, &(coords_sommets[0]), bary);

    double n[3]   = {0, 0, 1};
    double p0 [3] = {0 ,0, 0};
    double p10[3] = {0 ,0, 0};
    double l10[3] = {0 ,0, 0};
    double p20[3] = {0 ,0, 0};
    double l20[3] = {0 ,0, 0};

    /*Compute polygon normal*/
    PDM_polygon_parameterize (n_sommets, &(coords_sommets[0]),p0,p10,l10,p20,l20, n);

    PDM_plane_projection2 (local_pts, bary, n, local_pts);

    for (int isom = 0; isom < n_sommets; isom++) {

      double *pt1 = &(coords_sommets[0]) + 3 *isom;
      PDM_plane_projection2 (pt1, bary, n, pt1);

    }

    double bounds[6] = {DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX};

    for (int isom = 0; isom < n_sommets; isom++) {
      bounds[0] = PDM_MIN(bounds[0], coords_sommets[3*isom]);
      bounds[1] = PDM_MAX(bounds[1], coords_sommets[3*isom]);

      bounds[2] = PDM_MIN(bounds[2], coords_sommets[3*isom + 1]);
      bounds[3] = PDM_MAX(bounds[3], coords_sommets[3*isom + 1]);

      bounds[4] = PDM_MIN(bounds[4], coords_sommets[3*isom + 2]);
      bounds[5] = PDM_MAX(bounds[5], coords_sommets[3*isom + 2]);
    }


    /* Verification que le point est dans l'element */
    double closest[3];
    double dist_min = DBL_MAX;

    PDM_polygon_status_t position_inout
      = PDM_polygon_evaluate_position(local_pts, n_sommets, &(coords_sommets[0]), closest, &dist_min);

    if (position_inout == PDM_POLYGON_OUTSIDE) {
      local_pts[0] = closest[0];
      local_pts[1] = closest[1];
      local_pts[2] = closest[2];
    }

    /* Calcul des coordonnnees barycentriques */

    double min_dist = DBL_MAX;
    for (int isom = 0; isom < n_sommets; isom++) {

      int inext = (isom + 1) % n_sommets;
      double *vect = &s[0] + 3*isom;
      double l_edge;
      vect[0] = coords_sommets[3*inext]   - coords_sommets[3*isom];
      vect[1] = coords_sommets[3*inext+1] - coords_sommets[3*isom+1];
      vect[2] = coords_sommets[3*inext+2] - coords_sommets[3*isom+2];
      l_edge  = PDM_MODULE (vect);
      min_dist = PDM_MIN(l_edge, min_dist);
    }
    double eps = PDM_MAX(min_dist * eps_base, 1.e-30);

    for (int isom = 0; isom < n_sommets; isom++) {

      double *vect = &s[0] + 3*isom;
      vect[0] = coords_sommets[3*isom]   - local_pts[0];
      vect[1] = coords_sommets[3*isom+1] - local_pts[1];
      vect[2] = coords_sommets[3*isom+2] - local_pts[2];
      dist[isom] = PDM_MODULE (vect);
    }

    int currentVertex;
    for (int isom = 0; isom < n_sommets; isom++) {
      int inext = (isom + 1) % n_sommets;
      double *vect1 = &s[0] + 3 * isom;
      double *vect2 = &s[0] + 3 * inext;
      double pvect[3];

      proScal[isom] = PDM_DOT_PRODUCT (vect1, vect2);
      PDM_CROSS_PRODUCT(pvect, vect1, vect2);

      double sign = PDM_DOT_PRODUCT (pvect, n);
      aire[isom] = PDM_MODULE(pvect);

      if (sign < 0) {
        aire[isom] = -aire[isom];
      }

      if (dist[isom] <= eps) {

        isVertex = 1;
        currentVertex = isom;
        break;
      }

      else if ((fabs(aire[isom]) <= eps)  && (proScal[isom] < 0)) {

        isOnEdge = 1;
        currentVertex = isom;
        break;

      }

    }

    _barCoordsIndex[ipoint+1] = _barCoordsIndex[ipoint] + n_sommets;
    //Vector/Pointer containing Barycentric coordinates

    *barCoords = (double*)realloc( *barCoords, sizeof(double)*(_barCoordsIndex[ipoint+1]) );

    double *_barCoords = *barCoords;

    double* _localBaryCoords  = &(_barCoords[ _barCoordsIndex[ipoint] ]);

    /* Le point distant est un sommet */

    if (isVertex) {
      for (int isom = 0; isom < n_sommets; isom++)
        _localBaryCoords[isom] = 0.;
      _localBaryCoords[currentVertex] = 1.;
    }
    else if (isOnEdge) {
      /* Le point distant est sur arete */
      for (int isom = 0; isom < n_sommets; isom++)
        _localBaryCoords[isom] = 0.;

      int nextPoint = (currentVertex + 1) % n_sommets;

      _localBaryCoords[currentVertex] =
        dist[nextPoint]     / (dist[nextPoint]+dist[currentVertex]);
      _localBaryCoords[nextPoint]     =
        dist[currentVertex] / (dist[nextPoint]+dist[currentVertex]);
    }
    else {
      /* Cas general */
      double sigma = 0;
      for (int isom = 0; isom < n_sommets; isom++) {
        double coef = 0.;
        int previousVertex = (isom - 1 + n_sommets) % n_sommets;
        int nextVertex = (isom + 1) % n_sommets;

        if (fabs(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - proScal[previousVertex]/dist[isom]) / aire[previousVertex];
        if (fabs(aire[isom]) > eps)
          coef += (dist[nextVertex]     - proScal[isom]/dist[isom])           / aire[isom];
        sigma += coef;
        _localBaryCoords[isom] = coef;
      }

      if (PDM_ABS(sigma) >= eps ) {
        for (int isom = 0; isom < n_sommets; isom++) {
          _localBaryCoords[isom] /= sigma;
        }
      }
      else {
        double abs_sigma = fabs(sigma);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int isom = 0; isom < n_sommets; isom++) {
          _localBaryCoords[isom] = NAN;
        }
      }

      /* Check Result */

      for (int isom = 0; isom <  n_sommets; isom++) {
        if ( _localBaryCoords[isom] != _localBaryCoords[isom] ||
             _localBaryCoords[isom] < 0. ||
             _localBaryCoords[isom] > 1. ) {

          convergence = 0;
  /*        double dist_min = DBL_MAX;
          int k_min = 0;
          double t_min;

          for (int k = 0; k < n_sommets; k++) {
            _localBaryCoords[k] = 0.0;
          }

          for (int k = 0; k < n_sommets; k++) {
            double *p1 = &(coords_sommets[3 * k]);
            double *p2 = &(coords_sommets[3 * ((k+1) % n_sommets)]);
            double closest[3];
            double t;

            double dist2 = fvmc_distance_to_line (local_pts,
                                                 p1,
                                                 p2,
                                                 &t,
                                                 closest);
            if (dist2 < dist_min) {
              t_min = t;
              k_min = k;
            }
          }

          _localBaryCoords[k_min] = 1 - t_min;
          _localBaryCoords[(k_min + 1) % n_sommets] = t_min;
*/
          break;

        }

      }

    }

    if (0 == 1) {
      if ((nPoints == 1) && (ptsLocations[0] == 1)) {

        PDM_printf("coord %i %i :", ipoint+1, ielt+1);
        PDM_printf(" %12.5e %12.5e %12.5e", pts_coords[3*ipoint],
                    pts_coords[3*ipoint+1],
                    pts_coords[3*ipoint+2] );
        PDM_printf("\n");

        PDM_printf("coo b %i :", ipoint+1);
        for (int isom = 0; isom < n_sommets; isom++) {
          PDM_printf(" %f", _localBaryCoords[isom]);
        }
        PDM_printf("\n");
      }
    }
  }

  free(coords_sommets);
  free(s);
  free(aire);
  free(dist);
  free(proScal);

  return convergence;
}


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
)
{

  int convergence = 1;

  const double dispMin = 1e-9; /* Minimum displacement */
  const double big = 1e30;     /* Big value */

  const int nIterMax = 100;    /* Maximum iteration number */

  for (int ifac = 0; ifac < nPolygon; ifac++) {

    /*
     * Define local pointer
     * --------------------
     */

    const int nVerticesFace = connectivityIndex[ifac + 1]
      - connectivityIndex[ifac    ];

    const int *connectivityFace = connectivity + connectivityIndex[ifac];

    double *surfaceVectorFace = surfaceVector + 3*ifac; /* Face surface vector */

    double *centerFace = center + 3*ifac; /* Face Center */

    if (characteristicLength != NULL) {
      characteristicLength[ifac] = big;
    }

    /*
     * Initialization
     * --------------
     */

    int nIter = 0;       /* Number of iteration */

    for (int i = 0; i < 3; i++) {
      surfaceVectorFace[i] = 0;
      centerFace[i]        = 0;
    }

    /* Initialize face center to the barycenter */

    for (int ivert = 0; ivert < nVerticesFace; ivert++) {
      const int vert = connectivityFace[ivert] - 1;
      for (int i = 0; i < 3; i++)
        centerFace[i] += coords[3*vert + i];
    }

    for (int i = 0; i < 3; i++)
      centerFace[i] /= nVerticesFace;

    /*
     * Compute cell center and surface vector
     * --------------------------------------
     */

    while (1) {

      double displacement[3] = {0, 0, 0}; /* Cell center displacement */

      nIter += 1;

      for (int i = 0; i < 3; i++)
        surfaceVectorFace[i] = 0;

      double areaFace = 0; /* Face area */

      for (int ivert = 0; ivert < nVerticesFace; ivert++) {

        const int vert1 = connectivityFace[ivert] - 1;
        const int vert2 = connectivityFace[(ivert + 1) % nVerticesFace] - 1;

        /* Edge center */

        double edgeCenter[3] = {0.5 * (coords[3*vert1    ] + coords[3*vert2    ]),
                                0.5 * (coords[3*vert1 + 1] + coords[3*vert2 + 1]),
                                0.5 * (coords[3*vert1 + 2] + coords[3*vert2 + 2])};

        /* Edge vector */

        double vectV1V2[3] = {coords[3*vert2    ] - coords[3*vert1    ],
                              coords[3*vert2 + 1] - coords[3*vert1 + 1],
                              coords[3*vert2 + 2] - coords[3*vert1 + 2]};

        if (characteristicLength != NULL) {
          double norm_V1V2 = PDM_MODULE(vectV1V2);
          characteristicLength[ifac] = PDM_MIN(characteristicLength[ifac], norm_V1V2);
        }

        /* Vector face center -> edge center */

        double vectFECenter[3] = {edgeCenter[0] - centerFace[0],
                                  edgeCenter[1] - centerFace[1],
                                  edgeCenter[2] - centerFace[2]};

        /* Compute area of the triangle (face center, vert1, vert2) */

        double surfaceVectorTria[3];
        PDM_CROSS_PRODUCT(surfaceVectorTria, vectFECenter, vectV1V2);

        for (int i = 0; i < 3; i++)
          surfaceVectorTria[i] *= 0.5;

        const double areaTri = PDM_MODULE(surfaceVectorTria);

        areaFace += areaTri;
        for (int i = 0; i < 3; i++) {
          surfaceVectorFace[i] += surfaceVectorTria[i];
          displacement[i] += areaTri * vectFECenter[i];
        }

        areaFace += areaTri;

      }

      double denomAreaFace = 1. / PDM_MAX(fabs(areaFace), GEOM_EPS_MIN);

      for (int i = 0; i < 3; i++) {
        displacement[i] = 2./3. * denomAreaFace * displacement[i];
        centerFace[i] += displacement[i];
      }

      /*
       * Check convergence
       */

      const double normDisp = PDM_MODULE(displacement);

      if (normDisp < dispMin) {
        break;
      }

      /*
       * Check Number of iteration
       */

      else if (nIterMax < nIter) {
        convergence = false;
        break;
      }
    } /* while (1) */

    if ((characteristicLength != NULL) && (isDegenerated != NULL)) {

      double normSurfaceVector = PDM_MODULE(surfaceVectorFace);
      double eps_loc = PDM_geom_elem_geometric_epsilon(characteristicLength[ifac], GEOM_EPS_SURF);
      isDegenerated[ifac] = 0;
      if (normSurfaceVector <= eps_loc)
        isDegenerated[ifac] = 1;
    }
  } /* for (int ifac = 0; ifac < nPolygon; ifac++) */

  if (0 == 1) {
    PDM_printf( "surfacevector : ");
    for (int ifac = 0; ifac < 3*nPolygon; ifac++) {
      PDM_printf( "%12.5e ",surfaceVector[ifac]);
    }
    PDM_printf( "\n");

    PDM_printf( "center : ");
    for (int ifac = 0; ifac < 3*nPolygon; ifac++) {
      PDM_printf( "%12.5e ",center[ifac]);
    }
    PDM_printf( "\n");
  }
  return convergence;
}


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
)
{

  const int orientation = 1; /*  Surface vector oriented towards inside cell outside */

  const int nQuadrangle = 6;
  const int nTriangle = 0;
  const int nHexahedraFaces = nQuadrangle + nTriangle;
  const int nFaces = nHexahedraFaces * nHexahedra;

  int *faceConnectivity          = (int *) malloc (sizeof(int) * ((nQuadrangle*4 + nTriangle*3) * nHexahedra));
  int *faceConnectivityIdx       = (int *) malloc (sizeof(int) * (nFaces + 1));
  int *cellToFaceConnectivityIdx = (int *) malloc (sizeof(int) * (nHexahedra + 1));
  int *cellToFaceConnectivity    = (int *) malloc (sizeof(int) * (nFaces));

  /*
   * Get hexahedra faces
   */

  PDM_geom_elem_hexa_faces (nHexahedra,
                            orientation,
                            connectivity,
                            faceConnectivityIdx,
                            faceConnectivity);

  /*
   * Define cell to face connectivity
   */

  for (int i = 0; i < nFaces; i++)
    cellToFaceConnectivity[i] = i + 1;

  cellToFaceConnectivityIdx[0] = 0;
  for (int i = 1; i < nHexahedra + 1; i++)
    cellToFaceConnectivityIdx[i] = cellToFaceConnectivityIdx[i-1] + nHexahedraFaces;

  /*
   * Compute Volume and center
   */

  PDM_geom_elem_polyhedra_properties (1,
                                      nHexahedra,
                                      nFaces,
                                      faceConnectivityIdx,
                                      faceConnectivity,
                                      cellToFaceConnectivityIdx,
                                      cellToFaceConnectivity,
                                      nVertices,
                                      coords,
                                      volume,
                                      center,
                                      characteristicLength,
                                      isDegenerated);

  /*
   * Free
   */

  free (faceConnectivity);
  free (faceConnectivityIdx);
  free (cellToFaceConnectivity);
  free (cellToFaceConnectivityIdx);

}


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
)
{

  const int orientation = 1; //  Surface vector oriented towards inside cell outside

  const int nQuadrangle = 3;
  const int nTriangle = 2;
  const int nPrismFaces = nQuadrangle + nTriangle;
  const int nFaces = nPrismFaces * nPrism;

  int *faceConnectivity          = (int *) malloc (sizeof(int) * ((nQuadrangle*4 + nTriangle*3) * nPrism));
  int *faceConnectivityIdx       = (int *) malloc (sizeof(int) * (nFaces + 1));
  int *cellToFaceConnectivityIdx = (int *) malloc (sizeof(int) * (nPrism + 1));
  int *cellToFaceConnectivity    = (int *) malloc (sizeof(int) * (nFaces));

  /*
   * Get prism faces
   */

  PDM_geom_elem_prism_faces (nPrism,
                             orientation,
                             connectivity,
                             faceConnectivityIdx,
                             faceConnectivity);

  /*
   * Define cell to face connectivity
   */

  for (int i = 0; i < nFaces; i++)
    cellToFaceConnectivity[i] = i + 1;

  cellToFaceConnectivityIdx[0] = 0;
  for (int i = 1; i < nPrism + 1; i++)
    cellToFaceConnectivityIdx[i] = cellToFaceConnectivityIdx[i-1] + nPrismFaces;

  /*
   * Compute Volume and center
   */

  PDM_geom_elem_polyhedra_properties (1,
                                      nPrism,
                                      nFaces,
                                      faceConnectivityIdx,
                                      faceConnectivity,
                                      cellToFaceConnectivityIdx,
                                      cellToFaceConnectivity,
                                      nVertices,
                                      coords,
                                      volume,
                                      center,
                                      characteristicLength,
                                      isDegenerated);


  /*
   * Free
   */

  free (faceConnectivity);
  free (faceConnectivityIdx);
  free (cellToFaceConnectivity);
  free (cellToFaceConnectivityIdx);

}


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
)
{
  const int orientation = 1; /*  Surface vector oriented towards inside cell outside */

  const int nQuadrangle = 1;
  const int nTriangle = 4;
  const int nPyramidFaces = nQuadrangle + nTriangle;
  const int nFaces = nPyramidFaces * nPyramid;

  int *faceConnectivity          = (int *) malloc (sizeof(int) * ((nQuadrangle*4 + nTriangle*3) * nPyramid));
  int *faceConnectivityIdx       = (int *) malloc (sizeof(int) * (nFaces + 1));
  int *cellToFaceConnectivityIdx = (int *) malloc (sizeof(int) * (nPyramid + 1));
  int *cellToFaceConnectivity    = (int *) malloc (sizeof(int) * (nFaces));

  /*
   * Get pyramid faces
   */

  PDM_geom_elem_pyramid_faces (nPyramid,
                               orientation,
                               connectivity,
                               faceConnectivityIdx,
                               faceConnectivity);

  /*
   * Define cell to face connectivity
   */

  for (int i = 0; i < nFaces; i++)
    cellToFaceConnectivity[i] = i + 1;

  cellToFaceConnectivityIdx[0] = 0;
  for (int i = 1; i < nPyramid + 1; i++)
    cellToFaceConnectivityIdx[i] = cellToFaceConnectivityIdx[i-1] + nPyramidFaces;

  /*
   * Compute Volume and center
   */

  PDM_geom_elem_polyhedra_properties (1,
                                      nPyramid,
                                      nFaces,
                                      faceConnectivityIdx,
                                      faceConnectivity,
                                      cellToFaceConnectivityIdx,
                                      cellToFaceConnectivity,
                                      nVertices,
                                      coords,
                                      volume,
                                      center,
                                      characteristicLength,
                                      isDegenerated);

  /*
   * Free
   */

  free (faceConnectivity);
  free (faceConnectivityIdx);
  free (cellToFaceConnectivity);
  free (cellToFaceConnectivityIdx);

}


/**
 *  \brief Polyhedra properties
 *
 *  Compute cellCenter and volume. Reorient cellToFaceConnectivity
 *
 *  @param [in]  isOriented                 1 if cellToFaceConnectivity is already oriented, 0 otherwise
 *  @param [in]  nPolyhedra                 Number of polyhedra
 *  @param [in]  nFace                      Number of faces
 *  @param [in]  faceConnectivityIdx        Face connectivity index
 *  @param [in]  faceConnectivity           Face connectivity
 *  @param [in,out]  cellToFaceConnectivityIdx  Cell to face connectivity index
 *  @param [in,out]  cellToFaceConnectivity     Cell to face connectivity
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
 const int     isOriented,
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
 int          *isDegenerated
)
{
  const double big = 1e30;
  int convergence = 1;
  int *colorVertice = (int *) malloc (sizeof(int) * nVertices);

  int  lPolyhedraVertices = 24;
  int *polyhedraVertices = (int *) malloc (sizeof(int) *lPolyhedraVertices); //First allocation

  for (int i = 0; i < nVertices; i++)
    colorVertice[i] = false;

  /*
   * Compute face properties
   */

  double *surfaceVector = (double *) malloc (sizeof(double) * 3 *nFace);
  double *faceCenter    = (double *) malloc (sizeof(double) * 3 *nFace);

  // PDM_printf( "faceConnectivity : \n");
  // for (int ipoly = 0; ipoly < nFace; ipoly++) {
  //   PDM_printf( "  - face %i : ", ipoly+1);
  //   for (int j = faceConnectivityIdx[ipoly]; j < faceConnectivityIdx[ipoly+1]; j++) {
  //     PDM_printf( "%i ",faceConnectivity[j]);
  //   }
  //   PDM_printf( "\n");
  // }



  int convergenceFace = PDM_geom_elem_polygon_properties (nFace,
                                                          faceConnectivityIdx,
                                                          faceConnectivity,
                                                          coords,
                                                          surfaceVector,
                                                          faceCenter,
                                                          NULL,
                                                          NULL);

  if (0 == 1) {

    PDM_printf( "faceConnectivity : \n");
    for (int ipoly = 0; ipoly < nFace; ipoly++) {
      PDM_printf( "  - face %i : ", ipoly+1);
      for (int j = faceConnectivityIdx[ipoly]; j < faceConnectivityIdx[ipoly+1]; j++) {
        PDM_printf( "%i ",faceConnectivity[j]);
      }
      PDM_printf( "\n");
    }

    PDM_printf( "surfacevector : ");
    for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
      PDM_printf( "%12.5e ",surfaceVector[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "facecenter : ");
    for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
      PDM_printf( "%12.5e ",faceCenter[ipoly]);
    }
    PDM_printf( "\n");
  }

  /*
   * Loop on polyhedra
   */

  int keyMax = 2 * nVertices;

  PDM_hash_tab_t *hashOrient = NULL;
  int nKeyPoly = 0;
  int sKeyPoly = 10;

  int *keyPoly = NULL;

  double volume_t = 0;
  int maxNPolyFace = -1;

  int *stack = NULL;
  int *tagFace = NULL;

  if (!isOriented) {

    hashOrient = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT, &keyMax);

    keyPoly = (int *) malloc (sizeof(int) * sKeyPoly);

    for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
      const int polyIdx   = cellToFaceConnectivityIdx[ipoly];
      const int nPolyFace = cellToFaceConnectivityIdx[ipoly + 1] - polyIdx;
      maxNPolyFace = PDM_MAX (maxNPolyFace, nPolyFace);
    }

    stack = (int *) malloc (sizeof(int) * maxNPolyFace);
    tagFace = (int *) malloc (sizeof(int) * maxNPolyFace);

  }

  for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {

    double *polyCenter = center + 3*ipoly;
    double  disp[3];

    const int polyIdx   = cellToFaceConnectivityIdx[ipoly];
    const int nPolyFace = cellToFaceConnectivityIdx[ipoly + 1] - polyIdx;

    if (characteristicLength != NULL)
      characteristicLength[ipoly] = big;

    /*
     * Intialize cell center to the barycenter
     */

    nKeyPoly = 0;

    /* Search polyhedra vertices */

    int nPolyhedraVertices = 0;

    for (int iface = 0; iface < nPolyFace; iface++) {

      if (!isOriented) {
        tagFace[iface] = NOT_DEFINE;
      }

      const int face          = abs(cellToFaceConnectivity[polyIdx + iface]) - 1;
      const int faceIdx       = faceConnectivityIdx[face];
      const int nFaceVertices = faceConnectivityIdx[face+1] - faceIdx;

      for (int ivert = 0; ivert < nFaceVertices; ivert++) {
        const int vertex = faceConnectivity[faceIdx + ivert] - 1;

        if (!isOriented) {
          const int inext = (ivert + 1) % nFaceVertices;
          const int vertexNext = faceConnectivity[faceIdx + inext] - 1;
          const int key = vertex + vertexNext;

          if (nKeyPoly >= sKeyPoly) {
            sKeyPoly *= 2;
            keyPoly = (int *) realloc (keyPoly, sizeof(int) * sKeyPoly);
          }

          keyPoly[nKeyPoly++] = key;

          int *edge = malloc (sizeof(int)*3);
          edge[0] = vertex;
          edge[1] = vertexNext;
          edge[2] = iface;

          PDM_hash_tab_data_add (hashOrient, (void *) &key, edge);
        }

        if (!colorVertice[vertex]) {
          colorVertice[vertex] = 1;

          if (nPolyhedraVertices >= lPolyhedraVertices) {
            lPolyhedraVertices *= 2;
            polyhedraVertices =
                    (int *) realloc(polyhedraVertices, lPolyhedraVertices*sizeof(int));
          }
          polyhedraVertices[nPolyhedraVertices++] = vertex;
        }
      }
    }

    if (!isOriented) {
      int nStack = -1;
      stack[++nStack] = 0;
      tagFace[0] = IN_STACK;

      while (nStack >= 0) {

        int iFace = stack[nStack--];

        const int face          = abs(cellToFaceConnectivity[polyIdx + iFace]) - 1;
        const int faceIdx       = faceConnectivityIdx[face];
        const int nFaceVertices = faceConnectivityIdx[face+1] - faceIdx;

        for (int ivert = 0; ivert < nFaceVertices; ivert++) {
          const int inext = (ivert + 1) % nFaceVertices;

          const int vertex = faceConnectivity[faceIdx + ivert] - 1;
          const int vertexNext = faceConnectivity[faceIdx + inext] - 1;
          int key = vertex + vertexNext;

          int nData = PDM_hash_tab_n_data_get (hashOrient, &key);
          int **data = (int **) PDM_hash_tab_data_get (hashOrient, &key);

          int jCurrentEdge = -1;
          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
              int isSameFace    = iFace == _edge[2];
              if (isSameEdge && isSameFace) {
                jCurrentEdge = j;
                break;
              }
            }
          }

          assert (jCurrentEdge > -1);

          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isInverseEdge = (vertex == _edge[1]) && (vertexNext == _edge[0]);
              int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
              int isSameFace    = iFace == _edge[2];

              int neighbour = _edge[2];

              if (!isSameFace) {
                if (isSameEdge || isInverseEdge) {

                  if (tagFace[iFace] < UNCHANGED_CYCLE) {

                    if (tagFace[neighbour] >= UNCHANGED_CYCLE) {
                      if (tagFace[neighbour] == UNCHANGED_CYCLE) {
                        if (isSameEdge) {
                          tagFace[iFace] = CHANGED_CYCLE;
                        }
                        else  {
                          tagFace[iFace] = UNCHANGED_CYCLE;
                        }
                      }
                      else {
                        if (isSameEdge) {
                          tagFace[iFace] = UNCHANGED_CYCLE;
                        }
                        else  {
                          tagFace[iFace] = CHANGED_CYCLE;
                        }
                      }

                      free (data[j]);
                      free (data[jCurrentEdge]);
                      data[j] = NULL;
                      data[jCurrentEdge] = NULL;

                    }
                  }

                  if (tagFace[neighbour] == NOT_DEFINE) {
                    stack[++nStack] = neighbour;
                    tagFace[neighbour] = IN_STACK;
                  }

                  break;

                }
              }
            }
          }
        }

        if (tagFace[iFace] == IN_STACK) {
          tagFace[iFace] = UNCHANGED_CYCLE;
        }
      }


      for (int iface = 0; iface < nPolyFace; iface++) {

        if (tagFace[iface] == CHANGED_CYCLE) {
          cellToFaceConnectivity[polyIdx + iface] = -cellToFaceConnectivity[polyIdx + iface];
        }
      }

    }

    /* Compute barycenter */

    for (int i = 0; i < 3; i++)
      polyCenter[i] = 0.;

    for (int j = 0; j < nPolyhedraVertices; j++) {
      const int vertex = polyhedraVertices[j];
      colorVertice[vertex] = false;
      for (int i = 0; i < 3; i++)
        polyCenter[i] += coords[3*vertex + i];
    }

    for (int i = 0; i < 3; i++)
      polyCenter[i] /= nPolyhedraVertices;

    nPolyhedraVertices = 0;

    /*
     * Intialize volume
     */

    volume[ipoly] = 0.;

    /*
     * Intialize cell center displacement
     */

    for (int i = 0; i < 3; i++)
      disp[i] = 0.;

    /*
     * Loop on faces
     */

    for (int iface = 0; iface < nPolyFace; iface++) {

      const int face          = abs(cellToFaceConnectivity[polyIdx + iface]) - 1;
      const int direction     = (cellToFaceConnectivity[polyIdx + iface] < 0) ? -1 : 1;

      const int faceIdx       = faceConnectivityIdx[face];
      const int nFaceVertices = faceConnectivityIdx[face+1] - faceIdx;

      /*
       * Loop on vertices
       */

      for (int ivert = 0; ivert < nFaceVertices; ivert++) {

        int vert1;
        int vert2;
        int ivert1 = ivert;

        if (direction > 0) {
          vert1 = faceConnectivity[faceIdx + ivert1] - 1;
          vert2 = faceConnectivity[faceIdx + (ivert1 + 1) % nFaceVertices] - 1;
        }
        else {
          ivert1 = nFaceVertices - 1 - ivert;
          vert1 = faceConnectivity[faceIdx + (ivert1 + 1) % nFaceVertices] - 1;
          vert2 = faceConnectivity[faceIdx + ivert1] - 1;
        }

        if (characteristicLength != NULL)  {

          /* Vector vert1 -> vert2 */

          const double vectV1V2[3] =
            {coords[3*vert2    ] - coords[3*vert1    ],
             coords[3*vert2 + 1] - coords[3*vert1 + 1],
             coords[3*vert2 + 2] - coords[3*vert1 + 2]};

          double normV1V2 = PDM_MODULE (vectV1V2);

          characteristicLength[ipoly] = PDM_MIN(characteristicLength[ipoly], normV1V2);

        }

        /* Vector face center -> vert1 */

        double vectFCV1[3] =
          {coords[3*vert1    ] - faceCenter[3 * face    ],
           coords[3*vert1 + 1] - faceCenter[3 * face + 1],
           coords[3*vert1 + 2] - faceCenter[3 * face + 2]};

        /* Vector face center -> vert2 */

        double vectFCV2[3] =
          {coords[3*vert2    ] - faceCenter[3 * face    ],
           coords[3*vert2 + 1] - faceCenter[3 * face + 1],
           coords[3*vert2 + 2] - faceCenter[3 * face + 2]};

        /* Vector cell center -> face center */

        double vectCCFC[3] =
          {faceCenter[3 * face    ] - polyCenter[0],
           faceCenter[3 * face + 1] - polyCenter[1],
           faceCenter[3 * face + 2] - polyCenter[2]};

        double surfaceVectorTri[3];
        PDM_CROSS_PRODUCT (surfaceVectorTri, vectFCV1, vectFCV2);

        for (int i = 0; i < 3; i++)
          surfaceVectorTri[i] *= 0.5;

        /* Oriented volume */

        double volumeTet = 1./3 * PDM_DOT_PRODUCT (surfaceVectorTri, vectCCFC);

        volume[ipoly] += volumeTet;

        for (int i = 0; i < 3; i++)
          disp[i] = disp[i] + volumeTet * vectCCFC[i];

      }

    }

    int signeVol = (volume[ipoly] < 0.) ? -1 : 1 ;

    if (signeVol == -1) {

      if (!isOriented) {
        volume[ipoly] = -volume[ipoly];
        for (int iface = 0; iface < nPolyFace; iface++) {
          cellToFaceConnectivity[polyIdx + iface] = - cellToFaceConnectivity[polyIdx + iface];
        }
      }

      else {
        PDM_printf( "Warning polyhedraProperties : volume < 0 for polyhedron '%i'\n",
                    ipoly + 1);
      }
    }

    double denomVol = 1 / PDM_MAX(fabs(volume[ipoly]), GEOM_EPS_MIN);

    for (int i = 0; i < 3; i++)
      polyCenter[i] =  polyCenter[i] + signeVol * denomVol * disp[i];


    volume_t += volume[ipoly];

    /*
     * Check convergence
     */

    if (!convergenceFace)
      convergence = false;

    if ((characteristicLength != NULL) && (isDegenerated != NULL)) {

      double eps_loc = PDM_geom_elem_geometric_epsilon (characteristicLength[ipoly], GEOM_EPS_VOL);
      isDegenerated[ipoly] = 0;
      if (fabs(volume[ipoly]) <= eps_loc)
        isDegenerated[ipoly] = 1;
    }

    if (!isOriented) {
      PDM_hash_tab_purge (hashOrient, PDM_TRUE);
    }

  }

  free (polyhedraVertices);

  if (!isOriented) {
    free (keyPoly);
    free (tagFace);
    free (stack);
    PDM_hash_tab_free (hashOrient);
  }

  if (0 == 1) {

    PDM_printf( "surfacevector : ");
    for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
      PDM_printf( "%12.5e ",surfaceVector[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "facecenter : ");
    for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
      PDM_printf( "%12.5e ",faceCenter[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "isDegenrated : ");
    for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
      PDM_printf( "%i ",isDegenerated[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "characteristicLength  : ");
    for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
      PDM_printf( "%12.5e ",characteristicLength[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "volume  : ");
    for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
      PDM_printf( "%12.5e ",volume[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "center  : ");
    for (int ipoly = 0; ipoly < 3 * nPolyhedra; ipoly++) {
      PDM_printf( "%12.5e ",center[ipoly]);
    }
    PDM_printf( "\n");

  }

  free (surfaceVector);
  free (faceCenter);
  free (colorVertice);

  if (!convergence)
    PDM_printf( "Warning polyhedraProperties : some polyhedra faces are not planar\n");

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
