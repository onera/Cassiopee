#ifndef __PDM_POLY_SURF_GEN_H__
#define __PDM_POLY_SURF_GEN_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

//#define PROCF(x, y) x
#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

void PDM_poly_surf_gen
(
PDM_MPI_Comm     localComm,
double       xmin,
double       xmax,
double       ymin,
double       ymax,
int          haveRandom,
int          initRandom,
PDM_g_num_t   nx,
PDM_g_num_t   ny,
PDM_g_num_t  *nGFace,
PDM_g_num_t  *nGVtx,
PDM_g_num_t  *nGEdge,
int         *dNVtx,
double     **dVtxCoord,
int         *dNFace,
int        **dFaceVtxIdx,
PDM_g_num_t **dFaceVtx,
PDM_g_num_t **dFaceEdge,
int         *dNEdge,
PDM_g_num_t **dEdgeVtx,
PDM_g_num_t **dEdgeFace,
int         *nEdgeGroup,
int        **dEdgeGroupIdx,
PDM_g_num_t **dEdgeGroup
);

/* void PROCF(creemaillagepolygone2d_f, CREEMAILLAGEPOLYGONE2D_F) */
/* ( */
/* PDM_MPI_Fint* localFComm, */
/* double     *xmin, */
/* double     *xmax, */
/* double     *ymin, */
/* double     *ymax, */
/* int        *initRandom, */
/* PDM_g_num_t *nx, */
/* PDM_g_num_t *ny, */
/* int        *nVertex_f, */
/* double     *coords_f, */
/* int        *nFace_f, */
/* int        *faceVertexIdx_f, */
/* PDM_g_num_t *faceVertex_f, */
/* PDM_g_num_t *faceEdge_f,    */
/* int        *nEdge_f, */
/* PDM_g_num_t *edgeVertex_f, */
/* PDM_g_num_t *edgeFace_f */
/* ); */

#ifdef __cplusplus
}
#endif
#endif
