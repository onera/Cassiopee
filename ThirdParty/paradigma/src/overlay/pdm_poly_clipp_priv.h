/*
 * This file is an implementation of Greiner & Hormann algorithm with Foster
 * Overleft extension to remove degenerate cases
 *
 */ 

#ifndef __PDM_POLYGON_CLIPP_PRIV_H__
#define __PDM_POLYGON_CLIPP_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_poly_clipp.h"

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
 * \enum _poly_clipp_tag_t
 * \brief Vertex location                
 *
 */
 
typedef enum {

    POLY_CLIPP_LOC_NEXT,           /*!< Next */
    POLY_CLIPP_LOC_PREVIOUS,       /*!< Previous */
    
} _poly_clipp_loc_t;

/**
 * \struct _vertex_poly_t
 * \brief  Vertex data structure
 * 
 */

typedef struct _vertex_poly_t _vertex_poly_t;

struct _vertex_poly_t {

  const double           *coords;    /*!< Coordinates */
  double                  u;         /*!< Paretrization of intersection in Edge */
  struct _vertex_poly_t  *first;    /*!< First vertex */   
  
  PDM_g_num_t              gN;        /*!< Global number */
  PDM_g_num_t              gNEdge;    /*!< Edge global number only for vertex*/
  struct _vertex_poly_t  *next;      /*!< Link to next vertex or intersection */
  struct _vertex_poly_t  *previous;  /*!< Link to previous vertex or intersection */
  bool                    isect;     /*!< Intersection (True  : intersection,
                                      *                 False : vertex) */
  bool                    tag;       /*!< Entry or exit if intersection
                                      *  (True  : Entry,
                                      *   False : Exit)
                                      *   In or out if vertex
                                      *  (True  : In, 
                                      *   False : Out) */
  struct _vertex_poly_t  *neighbor;  /*!< Pointer to intersection in adjacent
                                      *   polygon */
 
  bool                    used;      /*!< Is used in su polygon */
    
};

/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POLY_CLIPP_H__ */
