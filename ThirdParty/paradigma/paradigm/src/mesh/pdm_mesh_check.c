/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2018       ONERA

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_check.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \brief Find unconnected vertex in a mesh connectivity
 *
 * If a vertex is not cited in the face->vtx connectivity, the function
 * removes it from the mesh to ensure contiguity
 *
 * \param [in, out] n_vtx                 Number of vertices
 * \param [in, out] l_face_vtx            Size of face->vtx connectivity
 * \param [in, out] face_vtx              Face->vtx connectivity
 * \param [in     ] nb_vtx_last_problem   last unconnected vertex
 *
 * \return the first unconnected vertex behind \ref nb_vtx_last_problem
 */


static PDM_g_num_t
_unconnected_vertex_find
(
PDM_g_num_t* n_vtx,
PDM_g_num_t* l_face_vtx,
PDM_g_num_t* face_vtx,
PDM_g_num_t  nb_vtx_last_problem
)
{

  // Déclarations
  PDM_g_num_t *check_som = NULL;
  PDM_g_num_t index_courant ;
  PDM_g_num_t nb_som_problem = 0 ;

  // Allocation
  check_som = (PDM_g_num_t*) calloc((*n_vtx +1), sizeof(PDM_g_num_t));

  // check_som à 1 pour tous les sommets cités dans la connectivité face->som
  for (PDM_g_num_t i=1; i < *l_face_vtx; i++) {
    index_courant = face_vtx[i] ;
    check_som[index_courant] = 1 ;
  }

  // si check_som est 0 quelquepart, c'est qu'il y a un trou dans la connectivité
  for (PDM_g_num_t i=nb_vtx_last_problem; i > 0; i--) {
    if (check_som[i] == 0) {
      nb_som_problem = i ;
      break ;
    }
  }

  free(check_som) ;

  return nb_som_problem ;
}


/**
 * \brief Remove unconnected vertex in a mesh connectivity
 *
 * If a vertex is not cited in the face->vtx connectivity, the function
 * removes it from the mesh to ensure contiguity
 *
 * \param [in, out] n_vtx                 Number of vertices
 * \param [in     ] nb_vtx_problem        unconnected vertex
 * \param [in     ] coords                Coordinates
 * \param [in, out] l_face_vtx            Size of face->vtx connectivity
 * \param [in, out] face_vtx              Face->vtx connectivity
 *
 */

static void
_unconnected_vertex_remove
(
PDM_g_num_t* n_vtx,
PDM_g_num_t nb_vtx_problem,
double* coords,
PDM_g_num_t* l_face_vtx,
PDM_g_num_t* face_vtx
)
{

  PDM_g_num_t index_courant ;

  // Coordonnees
  for (PDM_g_num_t i=nb_vtx_problem; i < *n_vtx; i++) {
    coords[3*(i-1)] = coords[3*(i-1)+3] ; //x(n-1) <- x(n)
    coords[3*(i-1)+1] = coords[3*(i-1)+1+3] ; //y(n-1) <- y(n)
    coords[3*(i-1)+2] = coords[3*(i-1)+2+3] ; //z(n-1) <- z(n)
  }

  // Connectivité faces->sommets
  for (PDM_g_num_t i=1; i < *l_face_vtx; i++) {
    index_courant = face_vtx[i] ;
    if (index_courant > nb_vtx_problem) {
      face_vtx[i] = index_courant - 1 ;
    }
  }
}

/*=============================================================================
 * Public function definitions
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
)
{

  PDM_Mesh_check_unconnected_vertex(nb_vtx,
                     l_face_vtx,
                     face_vtx,
                     coords,
                     nb_holes);
}


void PDM_Mesh_check_unconnected_vertex
(
PDM_g_num_t* nb_vtx,
PDM_g_num_t* face_vtx,
PDM_g_num_t* l_face_vtx,
double* coords,
int* nb_holes
)
{

  PDM_g_num_t nb_som_problem = 0 ; // (absolute) number of problematic vertex
  PDM_g_num_t nb_som_last_problem = *nb_vtx ; // problematic vertex from the previous iteration
  *nb_holes = 0 ;               // total number of holes in the mesh

  do {

    nb_som_problem = _unconnected_vertex_find // Find the number of a problematic vertex
                     (nb_vtx,
                      l_face_vtx,
                      face_vtx,
                      nb_som_last_problem) ;

    if (nb_som_problem != 0) {
      _unconnected_vertex_remove                // if found, fix it
        (nb_vtx,
         nb_som_problem,
         coords,
         l_face_vtx,
         face_vtx) ;

      (*nb_holes)++ ;
      nb_som_last_problem = nb_som_problem ;
    }

  } while (nb_som_problem != 0) ;

}



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
