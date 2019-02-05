/* 
 * File:   pdm_edges_intersect_priv.h
 * Author: equemera
 *
 * Created on April 12, 2016, 11:23 AM
 */

#ifndef PDM_EDGES_INTERSECT_PRIV_H
#define	PDM_EDGES_INTERSECT_PRIV_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"

#ifdef	__cplusplus
extern "C" {
#endif


/**
 * \struct _edges_intersect_t
 * \brief  Manage a set of edges
 * 
 *  PDM_edges_intersect_t manages a set of edge intersections
 *
 */

typedef struct {

  PDM_hash_tab_t *ht;          /*!< Hash table to store intersections key = nGEdgeA + nGEdgeB*/
  PDM_hash_tab_t *htA;         /*!< Hash table to store intersections key = nGEdgeA */
  PDM_hash_tab_t *htB;         /*!< Hash table to store intersections key = nGEdgeB */
  PDM_g_num_t      maxGNEdgeA;  /*!< Max global number of edges in mesh A */
  PDM_g_num_t      maxGNEdgeB;  /*!< Max global number of edges in mesh B */      
  double         vtxCarLengthTol; /*< Absolute tolerance for characteristic length */
  PDM_MPI_Comm        comm;        /*!< MSG Comm */
  int             sMSGComm;    /*!< Size of MSG Comm */

} _edges_intersect_t;


/**
 * \struct PDM_edges_intersect_res_t
 * \brief Result of the intersection between two edges
 *
 * 
 *  PDM_edges_intersect_res_t describes the result of the intersection 
 *  between two edges
 * 
 */

typedef struct {

  PDM_g_num_t              nGEdgeA;      /*!< Global number of meshA edge */
  PDM_g_num_t              nGEdgeB;      /*!< Global number of meshB edge */
  PDM_line_intersect_t    tIntersect;   /*!< Intersection type */
  PDM_g_num_t              originEdgeA;  /*!< Global number of meshA edge */
  PDM_g_num_t              originEdgeB;  /*!< Global number of meshB edge */

  int                     nNewPointsA;  /*!< Number of A intersection points */
  PDM_edges_intersect_point_t *oNewPointsA; /*!< Type of intersection */                               
  PDM_g_num_t              *gNumA;        /*!< Global number in A overlay mesh */
  PDM_g_num_t              *linkA;        /*!< Linked vertex in B Mesh */
  double                  *uA;           /*!< Location of intersection on edge A */
  double                  *coordsA;      /*!< Coordinates for A new points */
  
  int                     nNewPointsB;  /*!< Number of B intersection points */
  PDM_edges_intersect_point_t *oNewPointsB; /*!< Origin of B intersection points */                               
  PDM_g_num_t              *gNumB;        /*!< Global number in B overlay mesh */
  PDM_g_num_t              *linkB;      /*!< Linked vertex in A Mesh */
  double                  *uB;           /*!< Location of intersection on edge B */
  double                  *coordsB;      /*!< Coordinates for B new points */

} _edges_intersect_res_t;

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_EDGES_INTERSECT_PRIV_H */

