#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _dcube_t
 * \brief  Distributed cube
 * 
 * _dcube_t define a distributed mesh of a cube 
 *
 */

typedef struct  {
  PDM_MPI_Comm       comm;          /*!< MPI communicator                          */
  PDM_g_num_t   nVtxSeg;       /*!< Number of vertices in segments            */
  double         length;        /*!< Segment length                            */
  double         zero_x;          /*!< Coordinates of the origin                 */
  double         zero_y;          /*!< Coordinates of the origin                 */
  double         zero_z;          /*!< Coordinates of the origin                 */
  int            nFaceGroup;    /*!< Number of faces groups                    */
  int            dNCell;        /*!< Number of cells stored in this process    */
  int            dNFace;        /*!< Number of faces stored in this process    */
  int            dNVtx;         /*!< Number of vertices stored in this process */
  PDM_g_num_t  *dFaceCell;     /*!< Faces from cells connectivity             */
  int           *dFaceVtxIdx;   /*!< Faces from vertices connectivity index    */
  PDM_g_num_t  *dFaceVtx;      /*!< Faces from vertices connectivity          */
  double        *dVtxCoord;     /*!< Vertices coordinates                      */
  int           *dFaceGroupIdx; /*!< Faces groups index                        */
  PDM_g_num_t  *dFaceGroup;    /*!< Faces groups                              */
} _dcube_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dcubes  = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _dcube_t *
_get_from_id
(
 int  id
)
{
  _dcube_t *dcube = (_dcube_t *) PDM_Handles_get (_dcubes, id); 
    
  if (dcube == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "PDM_part_dcube error : Bad dcube identifier\n");
  }

  return dcube;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a distributed cube
 *
 * \param [out]  id             dcube identifier
 * \param [in]   comm           Communicator
 * \param [in]   nVtxSeg        Number of vertices in segments
 * \param [in]   length         Segment length
 * \param [in]   zero_x         Coordinates of the origin 
 * \param [in]   zero_y         Coordinates of the origin 
 * \param [in]   zero_z         Coordinates of the origin 
 *
 */

void
PDM_dcube_gen_init 
(
 int                *id,
 PDM_MPI_Comm        comm, 
 const PDM_g_num_t  nVtxSeg, 
 const double        length,
 const double        zero_x,
 const double        zero_y,
 const double        zero_z
)
{

  int nRank;
  int myRank;

  PDM_MPI_Comm_size(comm, &nRank);
  PDM_MPI_Comm_rank(comm, &myRank);
  
  /*
   * Search a dcube free id
   */

  if (_dcubes == NULL) {
    _dcubes = PDM_Handles_create (4);
  }

  _dcube_t *dcube = (_dcube_t *) malloc(sizeof(_dcube_t));
  
  *id = PDM_Handles_store (_dcubes, dcube);

  /*
   * Build dcube structure
   */

  dcube->comm    = comm;
  dcube->nVtxSeg = nVtxSeg;
  dcube->length  = length;
  dcube->zero_x  = zero_x;
  dcube->zero_y  = zero_y;
  dcube->zero_z  = zero_z;

  PDM_g_num_t nVtx      = nVtxSeg * nVtxSeg * nVtxSeg;
  PDM_g_num_t nFaceSeg  = nVtxSeg - 1;
  PDM_g_num_t nFace     = 3 * nFaceSeg * nFaceSeg * nVtxSeg;
  PDM_g_num_t nCell     = nFaceSeg * nFaceSeg * nFaceSeg; 
  PDM_g_num_t nFaceFace = nFaceSeg * nFaceSeg;
  PDM_g_num_t nVtxFace  = nVtxSeg * nVtxSeg;
  PDM_g_num_t nFaceLim  = 6 * nFaceFace;
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double step = length / (double) nFaceSeg;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif  
  PDM_g_num_t *distribVtx     = (PDM_g_num_t *) malloc((nRank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distribFace    = (PDM_g_num_t *) malloc((nRank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distribCell    = (PDM_g_num_t *) malloc((nRank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distribFaceLim = (PDM_g_num_t *) malloc((nRank + 1) * sizeof(PDM_g_num_t));

  //
  // Define distribution

  distribVtx[0]     = 0;
  distribFace[0]    = 0;
  distribCell[0]    = 0;
  distribFaceLim[0] = 0;

  PDM_g_num_t stepVtx = nVtx / nRank;
  PDM_g_num_t remainderVtx = nVtx % nRank;

  PDM_g_num_t stepFace = nFace / nRank;
  PDM_g_num_t remainderFace = nFace % nRank;

  PDM_g_num_t stepCell = nCell / nRank;
  PDM_g_num_t remainderCell = nCell % nRank;

  PDM_g_num_t stepFaceLim = nFaceLim / nRank;
  PDM_g_num_t remainderFaceLim = nFaceLim % nRank;

  for (int i = 1; i < nRank + 1; i++) {
    distribVtx[i]     = stepVtx;
    distribFace[i]    = stepFace;
    distribCell[i]    = stepCell;
    distribFaceLim[i] = stepFaceLim;
    const int i1 = i - 1;
    if (i1 < remainderVtx)
      distribVtx[i]  += 1;
    if (i1 < remainderFace)
      distribFace[i]  += 1;
    if (i1 < remainderCell)
      distribCell[i]  += 1;
    if (i1 < remainderFaceLim)
      distribFaceLim[i]  += 1;
  }

  for (int i = 1; i < nRank + 1; i++) {
    distribVtx[i]  += distribVtx[i-1];
    distribFace[i] += distribFace[i-1];
    distribCell[i] += distribCell[i-1];
    distribFaceLim[i] += distribFaceLim[i-1];
  }

  dcube->nFaceGroup = 6;
  PDM_g_num_t _dNCell = distribCell[myRank+1] - distribCell[myRank];
  dcube->dNCell       = (int) _dNCell;
  PDM_g_num_t _dNFace = distribFace[myRank+1]    - distribFace[myRank];
  dcube->dNFace       = (int) _dNFace;
  PDM_g_num_t _dNVtx  = distribVtx[myRank+1]     - distribVtx[myRank];
  dcube->dNVtx        = (int) _dNVtx;
  PDM_g_num_t _dNFaceLim = distribFaceLim[myRank+1] - distribFaceLim[myRank];
  int dNFaceLim = (int) _dNFaceLim;

  dcube->dFaceCell     = (PDM_g_num_t *) malloc(2*(dcube->dNFace) * sizeof(PDM_g_num_t *));
  dcube->dFaceVtxIdx   = (int *)          malloc((dcube->dNFace + 1) * sizeof(int *));
  dcube->dFaceVtx      = (PDM_g_num_t *) malloc(4*(dcube->dNFace)  * sizeof(PDM_g_num_t *));
  dcube->dVtxCoord     = (double *)       malloc(3*(dcube->dNVtx)  * sizeof(double *));
  dcube->dFaceGroupIdx = (int *)          malloc((dcube->nFaceGroup + 1)  * sizeof(int *));
  dcube->dFaceGroup    = (PDM_g_num_t *) malloc(dNFaceLim * sizeof(PDM_g_num_t *));

  PDM_g_num_t  *_dFaceCell     = dcube->dFaceCell;
  int           *_dFaceVtxIdx   = dcube->dFaceVtxIdx; 
  PDM_g_num_t  *_dFaceVtx      = dcube->dFaceVtx;
  double        *_dVtxCoord     = dcube->dVtxCoord;
  int           *_dFaceGroupIdx = dcube->dFaceGroupIdx;
  PDM_g_num_t  *_dFaceGroup    = dcube->dFaceGroup;

  _dFaceVtxIdx[0] = 0;
  for (int i = 1; i < dcube->dNFace + 1; i++) {
    _dFaceVtxIdx[i] = 4 + _dFaceVtxIdx[i-1];
  } 

  //
  // Coordinates
  
  const PDM_g_num_t bVtxZ = distribVtx[myRank] / nVtxFace;
  const PDM_g_num_t rVtxZ = distribVtx[myRank] % nVtxFace;

  const PDM_g_num_t bVtxY = rVtxZ / nVtxSeg;
  const PDM_g_num_t bVtxX = rVtxZ % nVtxSeg;
  
  int iVtx = 0;
  int cpt  = 0;

  for(PDM_g_num_t k = bVtxZ; k < nVtxSeg; k++) {
    PDM_g_num_t _bVtxY = 0;
    if (k == bVtxZ)
      _bVtxY = bVtxY; 
    for(PDM_g_num_t j = _bVtxY; j < nVtxSeg; j++) {
      PDM_g_num_t _bVtxX = 0;
      if ((k == bVtxZ) && (j == bVtxY))
        _bVtxX = bVtxX; 
      for(PDM_g_num_t i = _bVtxX; i < nVtxSeg; i++) {
        _dVtxCoord[3 * iVtx    ] = i * step + zero_x; 
        _dVtxCoord[3 * iVtx + 1] = j * step + zero_y;
        _dVtxCoord[3 * iVtx + 2] = k * step + zero_z;
        cpt += 1;
        iVtx += 1;
        if (cpt == dcube->dNVtx)
          break;
      }
      if (cpt == dcube->dNVtx)
        break;
    }
    if (cpt == dcube->dNVtx)
      break;
  }
  
  //
  // faceVtx et faceCell
  
  cpt = 0;
  
  PDM_g_num_t serie  = nFace / 3;
  PDM_g_num_t iSerie = distribFace[myRank] / serie;
  PDM_g_num_t rSerie = distribFace[myRank] % serie;

  PDM_g_num_t b1 = 0;
  PDM_g_num_t r1 = 0;
    
  PDM_g_num_t b2 = 0;
  PDM_g_num_t b3 = 0;

  b1 = rSerie / nFaceFace;
  r1 = rSerie % nFaceFace;
    
  b2 = r1 / nFaceSeg;
  b3 = r1 % nFaceSeg;

  switch (iSerie) {

  case 0 :

    //
    // Faces zmin -> zmax
  
    for(PDM_g_num_t k = b1; k < nVtxSeg; k++) {
      PDM_g_num_t _b2 = 0;
      if (k == b1)
        _b2 = b2; 
      for(PDM_g_num_t j = _b2; j < nFaceSeg; j++) {
        PDM_g_num_t _b3 = 0;
        if ((k == b1) && (j == b2))
          _b3 = b3; 
        for(PDM_g_num_t i = _b3; i < nFaceSeg; i++) {
          _dFaceVtx[cpt * 4    ] = k * nVtxSeg * nVtxSeg + (    j * nVtxSeg + i + 1);
          _dFaceVtx[cpt * 4 + 1] = k * nVtxSeg * nVtxSeg + ((j+1) * nVtxSeg + i + 1);
          _dFaceVtx[cpt * 4 + 2] = k * nVtxSeg * nVtxSeg + ((j+1) * nVtxSeg + i + 2); 
          _dFaceVtx[cpt * 4 + 3] = k * nVtxSeg * nVtxSeg + (    j * nVtxSeg + i + 2);
          if (k == 0) {
            _dFaceCell[2*cpt + 0] = j * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] = 0;
          }
          else if (k == nFaceSeg) {
            _dFaceCell[2*cpt + 0] = (k-1) * nFaceSeg * nFaceSeg + j * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] = 0;
          }
          else {
            _dFaceCell[2*cpt + 0] = (k-1) * nFaceSeg * nFaceSeg + j * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] =     k * nFaceSeg * nFaceSeg + j * nFaceSeg + i + 1;
          }
          cpt += 1;
          if (cpt == dcube->dNFace)
            break;
        }
        if (cpt == dcube->dNFace)
          break;
      }
      if (cpt == dcube->dNFace)
        break;
    }

    b1 = 0;
    b2 = 0;
    b3 = 0;

    if (cpt == dcube->dNFace)
      break;

  case 1 :
  
    //
    // Faces xmin -> xmax
  
    for(PDM_g_num_t i = b1; i < nVtxSeg; i++) {
      PDM_g_num_t _b2 = 0;
      if (i == b1)
        _b2 = b2; 
      for(PDM_g_num_t k = _b2; k < nFaceSeg; k++) {
        PDM_g_num_t _b3 = 0;
        if ((i == b1) && (k == b2))
          _b3 = b3; 
        for(PDM_g_num_t j = _b3; j < nFaceSeg; j++) {
          _dFaceVtx[cpt * 4    ] =     k * nVtxSeg * nVtxSeg +     j * nVtxSeg + i + 1;
          _dFaceVtx[cpt * 4 + 1] =     k * nVtxSeg * nVtxSeg + (j+1) * nVtxSeg + i + 1;
          _dFaceVtx[cpt * 4 + 2] = (k+1) * nVtxSeg * nVtxSeg + (j+1) * nVtxSeg + i + 1;
          _dFaceVtx[cpt * 4 + 3] = (k+1) * nVtxSeg * nVtxSeg +     j * nVtxSeg + i + 1;
          if (i == 0) {
            _dFaceCell[2*cpt + 0] = k * nFaceSeg * nFaceSeg + j * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] = 0;
          }
          
          else if (i == nFaceSeg) {
            _dFaceCell[2*cpt + 0] = k * nFaceSeg * nFaceSeg + j * nFaceSeg + i;
            _dFaceCell[2*cpt + 1] = 0;
          }
          
          else {
            _dFaceCell[2*cpt + 0] = k * nFaceSeg * nFaceSeg + j * nFaceSeg + i ;
            _dFaceCell[2*cpt + 1] = k * nFaceSeg * nFaceSeg + j * nFaceSeg + i + 1;
          }
          cpt += 1;
          if (cpt == dcube->dNFace)
            break;
        }
        if (cpt == dcube->dNFace)
          break;
      }
      if (cpt == dcube->dNFace)
        break;
    }

    b1 = 0;
    b2 = 0;
    b3 = 0;

    if (cpt == dcube->dNFace)
      break;

  case 2 :

    //
    // Faces ymin -> ymax

    for(PDM_g_num_t j = b1; j < nVtxSeg; j++) {
      PDM_g_num_t _b2 = 0;
      if (j == b1)
        _b2 = b2; 
      for(PDM_g_num_t i = _b2; i < nFaceSeg; i++) {
        PDM_g_num_t _b3 = 0;
        if ((j == b1) && (i == b2))
          _b3 = b3; 
        for(PDM_g_num_t k = _b3; k < nFaceSeg; k++) {
          _dFaceVtx[cpt * 4    ] =     k * nVtxSeg * nVtxSeg + j * nVtxSeg + i + 1    ;
          _dFaceVtx[cpt * 4 + 1] =     k * nVtxSeg * nVtxSeg + j * nVtxSeg + i + 1 + 1;
          _dFaceVtx[cpt * 4 + 2] = (k+1) * nVtxSeg * nVtxSeg + j * nVtxSeg + i + 1 + 1;
          _dFaceVtx[cpt * 4 + 3] = (k+1) * nVtxSeg * nVtxSeg + j * nVtxSeg + i + 1    ;
          if (j == 0) {
            _dFaceCell[2*cpt + 0] = k * nFaceSeg * nFaceSeg + j * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] = 0;
          }
          
          else if (j == nFaceSeg) {
            _dFaceCell[2*cpt + 0] =  k * nFaceSeg * nFaceSeg + (j-1) * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] = 0;
          }
          
          else {
            _dFaceCell[2*cpt + 0] = k * nFaceSeg * nFaceSeg + (j-1) * nFaceSeg + i + 1;
            _dFaceCell[2*cpt + 1] = k * nFaceSeg * nFaceSeg +     j * nFaceSeg + i + 1;
          }
          cpt += 1;
          if (cpt == dcube->dNFace)
            break;
        }
        if (cpt == dcube->dNFace)
          break;
      }
      if (cpt == dcube->dNFace)
        break;
    }
  }

  //
  // Faces limite

  cpt = 0;
  PDM_g_num_t bFace;
  int cpt1 = 0;
  int cpt3 = 0;
  int firstGroup = 0;
  
  serie  = nFaceLim / dcube->nFaceGroup;
  iSerie = distribFaceLim[myRank] / serie;
  rSerie = distribFaceLim[myRank] % serie;

  for (int i = 0; i < dcube->nFaceGroup + 1; i++)
    _dFaceGroupIdx[i] = 0;

  switch (iSerie) {

  case 0 :

    //
    // Faces zmin 

    if (cpt == 0)
      firstGroup = 1;
    
    cpt1 = cpt;

    bFace = 0;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < nFaceSeg; j++) {
      for(PDM_g_num_t i = 0; i < nFaceSeg; i++) {
        cpt3 += 1;
        if (!firstGroup || (firstGroup && ((cpt3 - 1)  >= rSerie))) {
          _dFaceGroup[cpt] = bFace + j * nFaceSeg + i + 1;
          cpt += 1;
	  if (cpt == dNFaceLim)
	    break;
	}
      }
      if (cpt == dNFaceLim)
        break;
    }
    
    _dFaceGroupIdx[1] = cpt - cpt1;

    if (cpt == dNFaceLim)
      break;

    firstGroup = 0;

  case 1 :

    //
    // Faces zmax 

    if (cpt == 0)
      firstGroup = 1;
    
    cpt1 = cpt;

    bFace = nFaceSeg * nFaceSeg * nFaceSeg;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < nFaceSeg; j++) {
      for(PDM_g_num_t i = 0; i < nFaceSeg; i++) {
        cpt3 += 1;
        if (!firstGroup || (firstGroup && ((cpt3 - 1)  >= rSerie))) {
          _dFaceGroup[cpt] = bFace + j * nFaceSeg + i + 1;
          cpt += 1;
          if (cpt == dNFaceLim)
            break;
        }
      }
      if (cpt == dNFaceLim)
        break;
    }

    _dFaceGroupIdx[2] = cpt - cpt1;

    if (cpt == dNFaceLim)
      break;

    firstGroup = 0;
          
  case 2 :

    //
    // Faces xmin 
  
    if (cpt == 0)
      firstGroup = 1;

    cpt1 = cpt;

    bFace = nFaceSeg * nFaceSeg * nVtxSeg;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < nFaceSeg; j++) {
      for(PDM_g_num_t i = 0; i < nFaceSeg; i++) {
        cpt3 += 1;
        if (!firstGroup || (firstGroup && ((cpt3 - 1)  >= rSerie))) {
          _dFaceGroup[cpt] = bFace + j * nFaceSeg + i + 1;
          cpt += 1;
          if (cpt == dNFaceLim)
            break;
        }
      }
      if (cpt == dNFaceLim)
        break;
    }

    _dFaceGroupIdx[3] = cpt - cpt1;

    if (cpt == dNFaceLim)
      break;

    firstGroup = 0;

  case 3 :

    //
    // Faces xmax 
  
    if (cpt == 0)
      firstGroup = 1;

    cpt1 = cpt;

    bFace = nFaceSeg * nFaceSeg * (nVtxSeg + nFaceSeg);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < nFaceSeg; j++) {
      for(PDM_g_num_t i = 0; i < nFaceSeg; i++) {
        cpt3 += 1;
        if (!firstGroup || (firstGroup && ((cpt3 - 1)  >= rSerie))) {
          _dFaceGroup[cpt] = bFace + j * nFaceSeg + i + 1;
          cpt += 1;
          if (cpt == dNFaceLim)
            break;
        }
      }
      if (cpt == dNFaceLim)
        break;
    }

    _dFaceGroupIdx[4] = cpt - cpt1;

    if (cpt == dNFaceLim)
      break;

    firstGroup = 0;

  case 4 :

    //
    // Faces ymin 
  
    if (cpt == 0)
      firstGroup = 1;

    cpt1 = cpt;

    bFace = nFaceSeg * nFaceSeg * (nVtxSeg + nVtxSeg);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < nFaceSeg; j++) {
      for(PDM_g_num_t i = 0; i < nFaceSeg; i++) {
        cpt3 += 1;
        if (!firstGroup || (firstGroup && ((cpt3 - 1)  >= rSerie))) {
          _dFaceGroup[cpt] = bFace + j * nFaceSeg + i + 1;
          cpt += 1;
          if (cpt == dNFaceLim)
            break;
        }
      }
      if (cpt == dNFaceLim)
        break;
    }

    _dFaceGroupIdx[5] = cpt - cpt1;

    if (cpt == dNFaceLim)
      break;

    firstGroup = 0;

  case 5 :

    //
    // Faces ymax 
  
    if (cpt == 0)
      firstGroup = 1;

    cpt1 = cpt;

    bFace = nFaceSeg * nFaceSeg * (nVtxSeg + nVtxSeg + nFaceSeg);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < nFaceSeg; j++) {
      for(PDM_g_num_t i = 0; i < nFaceSeg; i++) {
        cpt3 += 1;
        if (!firstGroup || (firstGroup && ((cpt3 - 1)  >= rSerie))) {
          _dFaceGroup[cpt] = bFace + j * nFaceSeg + i + 1;
          cpt += 1;
          if (cpt == dNFaceLim)
            break;
        }
      }
      if (cpt == dNFaceLim)
        break;
    }

    _dFaceGroupIdx[6] = cpt - cpt1;

    if (cpt == dNFaceLim)
      break;

    firstGroup = 0;

  }

  for (int i = 1; i < dcube->nFaceGroup + 1; i++)
    _dFaceGroupIdx[i] += _dFaceGroupIdx[i-1];

  free(distribVtx);
  free(distribFace);
  free(distribCell);
  free(distribFaceLim);

}

void
PROCF (pdm_dcube_gen_init, PDM_DCUBE_GEN_INIT)  
(
 int                *id,
 const PDM_MPI_Fint *comm,
 const PDM_g_num_t  *nVtxSeg, 
 const double       *length,
 const double       *zero_x,
 const double       *zero_y,
 const double       *zero_z
)
{

  PDM_MPI_Fint comm1 = *((PDM_MPI_Fint *) comm);

  PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c(comm1);

  PDM_dcube_gen_init (id,
                      c_comm,
                      *nVtxSeg, 
                      *length,
		      *zero_x,
		      *zero_y,
		      *zero_z);
}


/**
 *
 * \brief Return distributed cube size
 *
 * \param [in]   id          dcube identifier
 * \param [out]  NFaceGroup  Number of faces groups
 * \param [out]  dNCell      Number of cells stored in this process 
 * \param [out]  dNFace      Number of faces stored in this process
 * \param [out]  dNVtx       Number of vertices stored in this process
 * \param [out]  dFaceVtxL   Length of dFaceVtx array
 * \param [out]  dFacegroupL Length of dFacegroup array
 *
 */

void
PDM_dcube_gen_dim_get 
(
 int                id,
 int                *nFaceGroup,
 int                *dNCell,
 int                *dNFace,
 int                *dNVtx,
 int                *dFaceVtxL,
 int                *dFacegroupL
)
{
  _dcube_t *dcube = _get_from_id(id);

  *nFaceGroup = dcube->nFaceGroup;
  *dNCell     = dcube->dNCell;
  *dNFace     = dcube->dNFace;
  *dNVtx      = dcube->dNVtx;
  *dFaceVtxL  = dcube->dFaceVtxIdx[dcube->dNFace]; 
  *dFacegroupL= dcube->dFaceGroupIdx[dcube->nFaceGroup];
}


void
PROCF(pdm_dcube_gen_dim_get, PDM_DCUBE_GEN_DIM_GET)
(
 int                *id,
 int                *nFaceGroup,
 int                *dNCell,
 int                *dNFace,
 int                *dNVtx,
 int                *dFaceVtxL,
 int                *dFacegroupL
)
{
  PDM_dcube_gen_dim_get (*id,
                         nFaceGroup,
                         dNCell,
                         dNFace,
                         dNVtx,
                         dFaceVtxL,
                         dFacegroupL);

}

/**
 *
 * \brief Return distributed cube data
 *
 * \param [in]  id            dcube identifier
 * \param [out] dFaceCell     Faces from cells connectivity (size = 2 * dNFace)
 * \param [out] dFaceVtxIdx   Faces from vertices connectivity index (size = dNface + 1)
 * \param [out] dFaceVtx      Faces from vertices connectivity (size = dFaceVtxL)
 * \param [out] dVtxCoord     Vertices coordinates (size = 3 * dNVtx)
 * \param [out] dFaceGroupIdx Faces groups index (size = NFaceGroup + 1)
 * \param [out] dFaceGroup    Faces groups (size = dFacegroupL)
 *
 */

void
PDM_dcube_gen_data_get 
(
 int                 id,
 PDM_g_num_t      **dFaceCell,
 int               **dFaceVtxIdx, 
 PDM_g_num_t      **dFaceVtx,
 double            **dVtxCoord,
 int               **dFaceGroupIdx,
 PDM_g_num_t      **dFaceGroup 
)
{
  _dcube_t *dcube = _get_from_id(id);

  *dFaceCell     = dcube->dFaceCell;
  *dFaceVtxIdx   = dcube->dFaceVtxIdx;
  *dFaceVtx      = dcube->dFaceVtx;
  *dVtxCoord     = dcube->dVtxCoord;
  *dFaceGroupIdx = dcube->dFaceGroupIdx;
  *dFaceGroup    = dcube->dFaceGroup;
} 


void 
PROCF (pdm_dcube_gen_data_get, PDM_DCUBE_GEN_DATA_GET)
(
 int               *id,
 PDM_g_num_t      *dFaceCell,
 int               *dFaceVtxIdx, 
 PDM_g_num_t      *dFaceVtx,
 double            *dVtxCoord,
 int               *dFaceGroupIdx,
 PDM_g_num_t      *dFaceGroup 
)
{
  _dcube_t *dcube = _get_from_id(*id);

  for (int i = 0; i < 2 * dcube->dNFace; i++)
    dFaceCell[i] = dcube->dFaceCell[i];

  for (int i = 0; i < dcube->dNFace + 1 ; i++)
    dFaceVtxIdx[i] = dcube->dFaceVtxIdx[i]; 

  for (int i = 0; i < dcube->dFaceVtxIdx[dcube->dNFace] ; i++)
    dFaceVtx[i] = dcube->dFaceVtx[i]; 
    
  for (int i = 0; i < 3 * dcube->dNVtx; i++)
    dVtxCoord[i] = dcube->dVtxCoord[i];

  for (int i = 0; i < dcube->nFaceGroup + 1; i++)
    dFaceGroupIdx[i] = dcube->dFaceGroupIdx[i]; 

  for (int i = 0; i < dcube->dFaceGroupIdx[dcube->nFaceGroup]; i++)
    dFaceGroup[i] = dcube->dFaceGroup[i]; 
} 


/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  id            dcube identifier
 *
 */

void
PDM_dcube_gen_free
(
 int id
 )
{
  _dcube_t *dcube = _get_from_id(id);

  if (dcube->dFaceCell  != NULL)
    free(dcube->dFaceCell);

  if (dcube->dFaceVtxIdx  != NULL)
    free(dcube->dFaceVtxIdx);

  if (dcube->dFaceVtx  != NULL)
    free(dcube->dFaceVtx);

  if (dcube->dVtxCoord  != NULL)
    free(dcube->dVtxCoord);

  if (dcube->dFaceGroupIdx  != NULL)
    free(dcube->dFaceGroupIdx);

  if (dcube->dFaceGroup  != NULL)
    free(dcube->dFaceGroup);

  free(dcube);
  
  PDM_Handles_handle_free (_dcubes, id, PDM_FALSE);
  
  const int n_dcube = PDM_Handles_n_get (_dcubes);
  
  if (n_dcube == 0) {
    PDM_Handles_free (_dcubes);
    
  }
}
 

void 
PROCF (pdm_dcube_gen_free, PDM_DCUBE_GEN_FREE)
(
 int *id
 )
{
  PDM_dcube_gen_free (*id);
} 
