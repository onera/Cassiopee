#ifndef __PDM_MESH_CHECK_H__
#define __PDM_MESH_CHECK_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Remove unconnected vertices in a mesh connectivity
 *
 * If a vertex is not cited in the face->vtx connectivity, the function
 * removes it from the mesh to ensure contiguity
 *
 * /TODO : PDM_Mesh_check_unconnected_vertex complexity is n^2. This function must be optimized
 *
 * \param [in, out] n_vtx       Number of vertices
 * \param [in, out] l_face_vtx  Size of face->vtx connectivity
 * \param [in, out] face_vtx    Face->vtx connectivity
 * \param [in, out] coords      Vertices coordinates
 * \param [in, out] n_holes     Number of holes
 *
 */

void PROCF (pdm_mesh_check_unconnected_vertex, PDM_MESH_CHECK_UNCONNECTED_VERTEX)
(
PDM_g_num_t *nb_vtx,
PDM_g_num_t *l_face_vtx,
PDM_g_num_t *face_vtx,
double      *coords,
int         *nb_holes
);

void PDM_Mesh_check_unconnected_vertex
(
PDM_g_num_t* nb_vtx,
PDM_g_num_t* face_vtx,
PDM_g_num_t* l_face_vtx,
double* coords,
int* nb_holes
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GNUM_H__ */
