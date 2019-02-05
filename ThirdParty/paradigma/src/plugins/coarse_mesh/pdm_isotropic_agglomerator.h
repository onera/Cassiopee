/* 
 * File:   pdm_isotropic_agglomerator.h
 * Author: Nicolas Lantos
 *
 * Created on November 18, 2017, 1:24 PM
 */

#ifndef __TEST_AGGLOMERATOR_ISOTROPIC_H__
#define __TEST_AGGLOMERATOR_ISOTROPIC_H__

#include "pdm.h"

#if !defined (__hpux) && !defined (_AIX) 
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif


#ifdef  __cplusplus
extern "C" {
#endif
  
/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Orient cell->face connectivity
 * 
 * At the output of th function, a face number in \ref cellFace is positive 
 * if surface normal is inside the cell, negative otherwise. \ref faceCell is 
 * oriented in the same way  
 *
 * \param [in]     nCell       Number of cells
 * \param [in]     nFace       Number of faces
 * \param [in]     nVtx        Number of vertices
 * \param [in]     coords      Vertices coordinates
 * \param [in]     cellFaceIdx Cell to face connectivity index (size = \ref nCell + 1)
 * \param [in, out]cellFace    Cell to face connectivity (size = cellFaceIdx[nCell])
 * \param [in, out]faceCell    face to cell connectivity (size = 2 * \ref nFace) or NULL
 * \param [in]     faceVtxIdx  face to vertex connectivity index (size = \ref nFace + 1)
 * \param [in]     faceVtx     face to vertex connectivity (size = faceVtxIdx[nFace])
 *
 */

void agglomerateOneLevel(int *sizes,

                         int *adjMatrix_row_ptr,
                         int *adjMatrix_col_ind,
                         double *adjMatrix_areaValues,
                         double *volumes,

                         int *arrayOfFineAnisotropicCompliantCells,

                         int *isOnFineBnd,
                         int *array_isOnValley,
                         int *array_isOnRidge,
                         int *array_isOnCorner,

                         int isFirstAgglomeration,
                         int isAnisotropic,

                         int *fineCellToCoarseCell,

                         int *agglomerationLines_Idx,
                         int *agglomerationLines,

                         int dimension,
                         int goalCard1,
                         int minCard1,
                         int maxCard3,
                         int checks,
                         int verbose);

void agglomerateOneLevel_v_Paradigma(int *sizes,

                                     int *adjMatrix_row_ptr,
                                     int *adjMatrix_col_ind,
                                     double *volumes,

                                     int *arrayOfFineAnisotropicCompliantCells,

                                     int *isOnFineBnd_l,
                                     int * faceCell,
                                     double *Face_area,

                                     int isFirstAgglomeration_int,
                                     int isAnisotropic_int,

                                     int *fineCellToCoarseCell,

                                     int *agglomerationLines_Idx,
                                     int *agglomerationLines,

                                     int dimension,
                                     int goalCard1,
                                     int minCard1,
                                     int maxCard3,
                                     int checks,
                                     int verbose);

#ifdef  __cplusplus
}
#endif

#endif /*__TEST_AGGLOMERATOR_ISOTROPIC_H__ */
