
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part_geom.h"
#include "pdm_hilbert.h"
#include "pdm_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


static void
_prepare_connectivity
(
const PDM_MPI_Comm    comm,
const PDM_g_num_t     dNCell,
const int            *dCellFaceIdx,
const PDM_g_num_t    *dCellFace,
const int            *dFaceVtxIdx,
const PDM_g_num_t    *dFaceVtx,
const PDM_g_num_t    *dFaceProc,
int                **faceVtx,
int                **faceVtxIdx,
int                  *sizeFaceVtxIdx,
const PDM_g_num_t    *dVtxProc,
const double         *dVtxCoord,
      double       **lVtxCoord
)
{
  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(comm, &myRank);
  PDM_MPI_Comm_size(comm, &nRank);

  /* Offset Distribution to work with block_to_part */
  PDM_g_num_t* dFaceProcLoc = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * nRank + 1);
  PDM_g_num_t* dVtxProcLoc  = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * nRank + 1);
  for (int i = 0; i < nRank+1; i++) {
    dFaceProcLoc[i] = dFaceProc[i]-1;
    dVtxProcLoc[i]  = dVtxProc[i]-1;
  }


  /* Verbose */
  if(0 == 1){
    PDM_printf("dFaceProc : \n");
    for (int i = 0; i < nRank+1; i++) {
      PDM_printf("[%i] -> %i \n ",i,  dFaceProc[i]);
      PDM_printf("[%i] -> %i \n ",i,  dFaceProcLoc[i]);
    }
    for (int i = 0; i < dNCell; i++) {
      int nFac = dCellFaceIdx[i+1] - dCellFaceIdx[i];
      printf("nFac[%i] :  %i \n ", i, nFac);
      for (int iVtx = dCellFaceIdx[i]; iVtx < dCellFaceIdx[i+1]; iVtx++) {
        printf("dCellFace[%i] :  "PDM_FMT_G_NUM" \n ", iVtx, dCellFace[iVtx]);
      }
    }
  }


  /* -------------------------------------------------------------- */
  PDM_block_to_part_t *ptb = PDM_block_to_part_create (dFaceProcLoc,
                                                       &dCellFace,
                                                       &dCellFaceIdx[dNCell],
                                                        1,
                                                        comm);

  /* Echange sur les nombre de sommet de chaque face */
  int dNFace = dFaceProc[myRank+1] - dFaceProc[myRank];

  int* nVtxFaceLoc = (int *) malloc (sizeof(int) * dNFace              );
  int* nVtxFace    = (int *) malloc (sizeof(int) * dCellFaceIdx[dNCell]);

  for (int i = 0; i < dNFace; i++) {
    nVtxFaceLoc[i] = dFaceVtxIdx[i+1] - dFaceVtxIdx[i];
    for (int iVtx = dFaceVtxIdx[i]; iVtx < dFaceVtxIdx[i+1]; iVtx++) {
    }
  }

  int strideOne = 1;

  PDM_block_to_part_exch (ptb,
                          sizeof(PDM_l_num_t),
                          PDM_STRIDE_CST,
                          &strideOne,
                          (void *) nVtxFaceLoc,
                          NULL,
                          (void **) &nVtxFace);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < dCellFaceIdx[dNCell]; i++) {
      printf("nVtxFace :  %i \n ", nVtxFace[i]);
    }
  }

  int ndCellFaceTot = 0;
  for (int i = 0; i < dCellFaceIdx[dNCell]; i++) {
    ndCellFaceTot += nVtxFace[i];
  }
  // printf("ndCellFaceTot : %i \n", ndCellFaceTot);

  /* Alloc field */
  PDM_g_num_t* lFaceVtx    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * ndCellFaceTot);

  PDM_block_to_part_exch (           ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     nVtxFaceLoc,
                          (void *)   dFaceVtx,
                                    &nVtxFace,
                          (void **) &lFaceVtx);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < dCellFaceIdx[dNCell]; i++) {
      printf("nVtxFace :  %i \n ", nVtxFace[i]);
    }
    for (int i = 0; i < ndCellFaceTot; i++) {
      printf("lFaceVtx[%i] :  "PDM_FMT_G_NUM" \n ",i, lFaceVtx[i]);
    }
  }

  PDM_block_to_part_free(ptb);
  /* -------------------------------------------------------------- */


  /* -------------------------------------------------------------- */
  ptb = PDM_block_to_part_create (dVtxProcLoc,
                                  (const PDM_g_num_t **) &lFaceVtx,
                                  &ndCellFaceTot,
                                  1,
                                  comm);

  int stride = 3;


  // int dNVtx = dVtxProc[myRank+1] - dVtxProc[myRank];
  // for (int i = 0; i < dNVtx; i++) {
  //   printf("dVtxCoord[%i] :  %12.5e \n ",i, dVtxCoord[i]);
  // }

  *lVtxCoord = (double *) malloc (sizeof(double) * 3 * ndCellFaceTot);
  PDM_block_to_part_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &stride,
                          (void *) dVtxCoord,
                          NULL,
                          (void **) lVtxCoord);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < ndCellFaceTot; i++) {
      printf("lVtxCoord[%i] :  %12.5e \n ",i, (*lVtxCoord)[i]);
    }
  }

  PDM_block_to_part_free(ptb);
  /* -------------------------------------------------------------- */

  *sizeFaceVtxIdx = dCellFaceIdx[dNCell] + 1 ;

  *faceVtx    = (int *) malloc( sizeof(int *) *   ndCellFaceTot              );
  *faceVtxIdx = (int *) malloc( sizeof(int *) * ( dCellFaceIdx[dNCell] + 1 ) );

  for (int i = 0; i < ndCellFaceTot; i++) {
    (*faceVtx)[i] = i+1;
  }

  (*faceVtxIdx)[0] = 0;
  for (int i = 0; i < dCellFaceIdx[dNCell]; i++) {
    (*faceVtxIdx)[i+1] = (*faceVtxIdx)[i] + nVtxFace[i];
  }

  /* Free memory */
  free(nVtxFace);
  free(nVtxFaceLoc);
  free(lFaceVtx);

}
/**
 *
 * \brief Compute cell center of elements
 *
 * \param [in]  comm         MPI Communicator
 * \param [in]  dNCell       Number of cells in the current process
 * \param [in]  dNFace       Number of faces in the current process
 * \param [in]  dNVtx        Number of vertices in the current process
 * \param [in]  dCellFaceIdx Index of cellFace
 * \param [in]  dCellFace    cell face connectivity in the current process
 * \param [in]  dFaceVtxIdx  Index of faceVtx
 * \param [in]  dFaceVtx     face vertex connectivity in the current process
 * \param [in]  dFaceProc    face distribution
 * \param [in]  dVtxCoord    coordinates of vertices
 * \param [in]  dVtxProc     Vertex distribution
 *
 * \param [out] cellCenter   Cell centers
 *
 */

static void
_compute_cellCenter
(
  const PDM_MPI_Comm  comm,
  const int           dNCell,
  const int          *dCellFaceIdx,
  const PDM_g_num_t  *dCellFace,
  const int          *dFaceVtxIdx,
  const PDM_g_num_t  *dFaceVtx,
  const PDM_g_num_t  *dFaceProc,
  const double       *dVtxCoord,
  const PDM_g_num_t  *dVtxProc,
  double             *cellCenter
)
{

  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(comm, &myRank);
  PDM_MPI_Comm_size(comm, &nRank);

  /* Local connectivity (Allocate in _prepareconnectivity */
  int sizeFaceVtxIdx = 0;

  int *faceVtx;
  int *faceVtxIdx;

  /* Local coordinates */
  double *lVtxCoord;

  /*
   * Get face connectivities from other process
   */
  _prepare_connectivity( comm,
                         dNCell,
                         dCellFaceIdx,
                         dCellFace,
                         dFaceVtxIdx,
                         dFaceVtx,
                         dFaceProc,
                        &faceVtx,
                        &faceVtxIdx,
                        &sizeFaceVtxIdx,
                         dVtxProc,
                         dVtxCoord,
                        &lVtxCoord);


  /* Verbose */
  if(0 == 1 ){
    PDM_printf("faceVtx : \n");
    for (int i1 = 0; i1 < sizeFaceVtxIdx; i1++) {
      for (int j = faceVtxIdx[i1]; j < faceVtxIdx[i1+1]; j++)
        PDM_printf(" %i", faceVtx[j]);
      PDM_printf("\n");
    }
  }

  /*
   * Compute cell centers
   */
  int iFace = 0;
  for (int iCell = 0; iCell < dNCell; iCell++) {

    int iBeg = dCellFaceIdx[iCell  ];
    int iEnd = dCellFaceIdx[iCell+1];

    cellCenter[3*iCell    ] = 0.;
    cellCenter[3*iCell + 1] = 0.;
    cellCenter[3*iCell + 2] = 0.;

    int nVtxOnCell = 0;

    /* Loop on all faces */
    for (int iFaceG = iBeg; iFaceG < iEnd; iFaceG++) {

      int iBegVtx = faceVtxIdx[iFace  ];
      int iEndVtx = faceVtxIdx[iFace+1];
      nVtxOnCell += iEndVtx - iBegVtx;

      /* Loop on all Vtx of current faces */
      for (int iVtx = iBegVtx; iVtx < iEndVtx; iVtx++){

        cellCenter[3*iCell    ] += lVtxCoord[3*iVtx    ];
        cellCenter[3*iCell + 1] += lVtxCoord[3*iVtx + 1];
        cellCenter[3*iCell + 2] += lVtxCoord[3*iVtx + 2];

      }

      /* Go to next face in local numerotation */
      iFace++;

    }

    /* Ponderate */
    cellCenter[3*iCell    ] = cellCenter[3*iCell    ]/nVtxOnCell;
    cellCenter[3*iCell + 1] = cellCenter[3*iCell + 1]/nVtxOnCell;
    cellCenter[3*iCell + 2] = cellCenter[3*iCell + 2]/nVtxOnCell;
  }

  /* Verbose */
  if(0 == 1){
    for (int iCell = 0; iCell < dNCell; iCell++) {
      PDM_printf(" --------- %i / %i \n ", iCell);
      PDM_printf(" %12.5e %12.5e %12.5e \n", cellCenter[3*iCell    ], cellCenter[3*iCell + 1], cellCenter[3*iCell + 2]);
    }
  }

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Perform geomtric patitionning
 *
 * \param [in]   method         Geometric method
 * \param [in]   nPart          Number of partition to build on this process
 * \param [in]   comm           Communicator
 * \param [in]   dNCell         Number of distributed cells
 * \param [in]   dNFace         Number of distributed faces
 * \param [in]   dNVtx          Number of distributed vertices
 * \param [in]   dCellFaceIdx   Distributed cell face connectivity index or NULL
 *                              (size : dNCell + 1, numbering : 0 to n-1)
 * \param [in]   dCellFace      Distributed cell face connectivity or NULL
 *                              (size : dFaceVtxIdx[dNCell], numbering : 1 to n)
 * \param [in]   dCellWeight    Cell weight (size : nCell) or NULL
 * \param [in]   dFaceVtxIdx    Distributed face to vertex connectivity index
 *                              (size : dNFace + 1, numbering : 0 to n-1)
 * \param [in]   dFaceVtx       Distributed face to vertex connectivity
 *                              (size : dFaceVtxIdx[dNFace], numbering : 1 to n)
 * \param [in]   dVtxCoord      Distributed vertex coordinates
 *                              (size : 3*dNVtx)
 * \param [inout]   dCellPart      Distributed cell partitioning
 *                              (size = dNCell)
 *
 */

void
PDM_part_geom
(
 PDM_part_geom_t     method,
 const int           nPart,
 const PDM_MPI_Comm  comm,
 const int           dNCell,
 const int          *dCellFaceIdx,
 const PDM_g_num_t  *dCellFace,
 const int          *dCellWeight,
 const int          *dFaceVtxIdx,
 const PDM_g_num_t  *dFaceVtx,
 const PDM_g_num_t  *dFaceProc,
 const double        *dVtxCoord,
 const PDM_g_num_t  *dVtxProc,
 int                *dCellPart
)
{

  assert (method == PDM_PART_GEOM_HILBERT);

  const int dim = 3;

  double *barycenterCoords = (double *) malloc (dNCell * 3 * sizeof(double ));

  PDM_hilbert_code_t *hilbertCodes     = (PDM_hilbert_code_t *) malloc (dNCell * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_code_t *tmp_hilbertCodes = (PDM_hilbert_code_t *) malloc (dNCell * sizeof(PDM_hilbert_code_t));

  /*
   * cell center computation
   */
  _compute_cellCenter (comm,
                       dNCell,
                       dCellFaceIdx,
                       dCellFace,
                       dFaceVtxIdx,
                       dFaceVtx,
                       dFaceProc,
                       dVtxCoord,
                       dVtxProc,
                       barycenterCoords);


  /** TRAITEMENT HILBERT FVM **/

	/** Initialisation **/

  double extents[2*dim]; /** DIM x 2**/

	/** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dNCell, barycenterCoords, extents, comm);

	/** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dNCell, barycenterCoords, hilbertCodes);

  for (int i = 0; i < dNCell; ++i) {
    tmp_hilbertCodes [i] = hilbertCodes [i];
	}

  ///** Calcul des index des codes Hilbert **/

  int * hilbert_order = (int * ) malloc (dNCell * sizeof(int));

  for (int i = 0; i < dNCell; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbertCodes, hilbert_order, dNCell);

  free(tmp_hilbertCodes);

  int nRank;
  PDM_MPI_Comm_size (comm, &nRank);

  PDM_hilbert_code_t *hilbertCodesIdx = (PDM_hilbert_code_t *) malloc ((nRank+1)*nPart * sizeof(PDM_hilbert_code_t));

  int * weight = (int *) malloc (dNCell * sizeof(int));
  if (dCellWeight != NULL) {
    for(int i = 0; i < dNCell; ++i) {
		  weight [i] = dCellWeight [i];
    }
  }
  else {
    for(int i = 0; i < dNCell; ++i) {
		  weight [i] = 1;
    }
  }

  int nTPart;
  PDM_MPI_Allreduce ((void *) &nPart, &nTPart, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

  PDM_hilbert_build_rank_index (dim,
                                nTPart,
                                dNCell,
                                hilbertCodes,
                                weight,
                                hilbert_order,
                                hilbertCodesIdx,
                                comm);

  free(weight);

  /** Remplissage de cellParts -> en fct des codes Hilbert **/

  for(int i = 0; i < dNCell; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(nTPart,
                                                hilbertCodes[i],
                                                hilbertCodesIdx);
    dCellPart [i] = (int) quantile;

  }

  free(barycenterCoords);
  free(hilbertCodesIdx);
  free(hilbert_order);
  free(hilbertCodes);

}



// /**
//  * \brief Remove duplicates in a array of long
//  *
//  * \param [in,out] array  Array to clean
//  * \param [in,out] size   Size of array to clean
//  *
//  */

// static void
// _remove_duplicate_pdm_part_long_t
// (
//  PDM_g_num_t *array,
//  int          *size
// )
// {

//   PDM_g_num_t * oldArray = (PDM_g_num_t * ) malloc( (*size) * sizeof(PDM_g_num_t));

//   for(int i = 0; i < *size; ++i ){
//     oldArray[i] = array[i];
//     array[i] = -1;
//   }
//   int * initialIdx = (int *) malloc( (*size) * sizeof(int));

//   for (int i = 0; i < *size; ++i ){
//     initialIdx[i] = i;
//   }

//   int k = 1;
//   int finalSize = 0;

//   array[0] = oldArray[0];
//   PDM_g_num_t last = array[0];
//   for (int i = 1; i < *size; ++i) {
//     PDM_g_num_t elem = oldArray[k];
//     if(oldArray[i] != elem  && oldArray[i] > last) {
//       k ++;
//       array[initialIdx[i]] = oldArray[i];
//       last = oldArray[i];
//     }
//   }

//   finalSize = k;

//   k = 1;
//   for (int i = 1; i < *size; ++i ){
//     if(array[i] != -1) {
//       array[k++] = array [i];
//     }
//   }

//   *size = finalSize;
//   free(oldArray);
//   free(initialIdx);
// }


// *
//  * \brief Sort and remove duplicate elements in an array of coordinates
//  *        coordinates are sorted in the local numbering order
//  *
//  * \param [in,out] array      Array to clean
//  * \param [in,out] size       Size of array to clean
//  * \param [in,out] orderArray Array that indicates the original numbering
//  *                            orderArray must be sorted and with duplicate elements
//  * \param [in,out] initialIdx Array that indicates the original numbering
//  *                            of orderArray after a quickSort


// static void
// _remove_duplicate_coords
// (
// double        *array,
// int            size,
// PDM_g_num_t *orderArray,
// int          *initialIdx
// )
// {

//   PDM_g_num_t * oldOrderArray = (PDM_g_num_t * ) malloc( (size) * sizeof(PDM_g_num_t));
//   double * oldArray = (double * ) malloc( (size) * 3 * sizeof(double));
//   for(int i = 0; i < size; ++i ){
//     oldOrderArray[i] = orderArray[i];

//     for(int j = 0; j < 3; ++j) {
//       oldArray[i * 3 +j ] = array [i * 3 +j ];
//     }
//   }

//   /*
//    * Sort
//    */

//   for (int i = 0; i < size; ++i ){
//     for(int j = 0; j < 3; ++j) {
//       array[i * 3 + j] = oldArray[initialIdx[i] * 3 + j];
//     }
//   }

//   int *tab_doublons = (int *) malloc ( size * sizeof(int) );
//   for (int i = 0; i < size; ++i) {
//     tab_doublons[i] = 0;
//   }

//   int k = 0;
//   PDM_g_num_t elem = oldOrderArray[0];
//   for (int i = 1; i < size; ++i) {
//     if (elem != oldOrderArray[i]) {
//       k++;
//       elem = oldOrderArray[i];
//     }
//     else {
//       tab_doublons[k]++;
//     }
//   }

//   /*
//    * Remove duplicate elements
//    */

//   int sizeTabDblIdx = k+2;
//   int * tabDblIdx = (int *) malloc (sizeTabDblIdx * sizeof(int));
//   tabDblIdx[0] = 0;
//   for (int i = 1; i < sizeTabDblIdx; ++i) {
//     tabDblIdx[i] = tabDblIdx[i-1] + tab_doublons[i-1]+1;
//   }
//   free(tab_doublons);


//   for(int i = 0; i < sizeTabDblIdx-1; ++i) {
//     for(int j = 0; j < 3; ++j) {
//       array [i * 3 +j ] = array [(tabDblIdx[i]) * 3 + j ];
//     }
//   }


//   free(oldArray);
//   free(oldOrderArray);
//   free(tabDblIdx);

// }
// /**
//  *
//  * \brief Distributes face arrays
//  *
//  * \param [in]  comm           MPI Communicator
//  * \param [in]  dNCell         Number of cells in the current process
//  * \param [in]  dCellFaceIdx   Index of cellFace
//  * \param [in]  dCellFace      cell face connectivity in the current process
//  * \param [in]  dFaceVtxIdx    Index of faceVtx
//  * \param [in]  dFaceVtx       face vertex connectivity in the current process
//  * \param [in]  dFaceProc      face distribution
//  * \param [out] faceVtx        faceVtx connectivity for all necessary faces
//  * \param [out] faceVtxIdx     Index of faceVtx connectivity
//  * \param [out] sizeFaceVtx    size of faceVtx
//  * \param [out] sizeFaceVtxIdx size of faceVtxIdx
//  *
//  */

// static void
// _distrib_face
// (
// const PDM_MPI_Comm  comm,
// const PDM_g_num_t   dNCell,
// const int          *dCellFaceIdx,
// const PDM_g_num_t  *dCellFace,
// const int          *dFaceVtxIdx,
// const PDM_g_num_t  *dFaceVtx,
// const PDM_g_num_t  *dFaceProc,
// PDM_g_num_t       **faceVtx,
// int                *faceVtxIdx,
// int                *sizeFaceVtx,
// int                *sizeFaceVtxIdx
// )
// {
//   int myRank;
//   int nRank;

//   PDM_MPI_Comm_rank(comm, &myRank);
//   PDM_MPI_Comm_size(comm, &nRank);

//   const int nData     = 1;
//   // int       nDataFace = 2;
//   int       nDataFace = 1;

//   int *faceToSendIdx = (int *) calloc((nRank + 1), sizeof(int));
//   int *faceToSendN   = (int *) calloc( nRank     , sizeof(int));

//   int *requestedFaceN   = (int *) malloc( nRank      * sizeof(int));
//   int *requestedFaceIdx = (int *) malloc((nRank + 1) * sizeof(int));

//   int size_dCellFace = dCellFaceIdx[dNCell];

//   if (1 == 1) {
//     // PDM_printf("dCellFace : \n");
//     // for (int i = 0; i < dNCell; i++) {
//     //   for (int j = dCellFaceIdx[i]; j < dCellFaceIdx[i+1]; j++)
//     //     PDM_printf(" "PDM_FMT_G_NUM, dCellFace[j]);
//     //   PDM_printf("\n");
//     // }
//     PDM_printf("dFaceProc : \n");
//     for (int i = 0; i < nRank+1; i++) {
//       PDM_printf("[%i] -> %i \n ",i,  dFaceProc[i]);
//     }

//   }

//   int * allToallNToLN = (int *) malloc(size_dCellFace * sizeof(int));

//   const PDM_g_num_t *cellFace = dCellFace;

//   for(int iFace = 0; iFace < size_dCellFace; iFace ++) {
//     int irank = PDM_binary_search_gap_long(PDM_ABS (cellFace[iFace]), dFaceProc, nRank+1);
//     // printf(" DENUG faceToSendIdx [%i] = %i - %i \n", irank, iFace, cellFace[iFace]);
//     faceToSendIdx[irank+1] += nData;
//   }

//   for (int i = 0; i < nRank; i++) {
//       faceToSendIdx[i+1] += faceToSendIdx[i];
//   }

//   int sizeFaceToSend = faceToSendIdx[nRank];
//   printf("RANK%d -- sizeFaceToSend = %d\n", myRank, sizeFaceToSend);
//   printf("RANK%d -- size_dCellFace = %d\n", myRank, size_dCellFace);
//   PDM_g_num_t * faceToSend = (PDM_g_num_t *) malloc(  sizeFaceToSend * sizeof(PDM_g_num_t));

//   for(int iFace = 0; iFace < size_dCellFace; iFace ++) {
//     int irank = PDM_binary_search_gap_long (PDM_ABS (cellFace[iFace]), dFaceProc, nRank+1);
//     int idx = faceToSendIdx[irank] + faceToSendN[irank];

//     // printf(" [%i] = %d\n", idx, iFace);
//     allToallNToLN[idx/nData] = iFace;
//     faceToSend[idx++]   = PDM_ABS (cellFace[iFace]);        // Face global numbering
//     faceToSendN[irank] += nData;
//   }

//   PDM_g_num_t *requestedFace    = NULL;

//   _alltoall(            faceToSend,
//                         faceToSendN,
//                         faceToSendIdx,
//              (void **) &requestedFace,
//                         requestedFaceN,
//                         requestedFaceIdx,
//                         PDM__PDM_MPI_G_NUM,
//                         sizeof(PDM_g_num_t),
//                         comm);

//   free(faceToSend);

//   /*
//    *  Processes exchange information about requested faces
//    *  For each face, information contains :
//    *     - tag (if ppart->_dFaceTag != NULL)
//    *     - Number of vertices
//    *     - Vertices
//    *
//    */

//   int *sFaceInfoIdx = faceToSendIdx; for(int i = 0; i < nRank+1; ++i) sFaceInfoIdx[i] = 0;
//   int *sFaceInfoN   = faceToSendN; for(int i = 0; i < nRank; ++i) sFaceInfoN[i] = 0;

//   for (int i = 0; i < nRank; i++) {
//     for (int k = requestedFaceIdx[i]; k < requestedFaceIdx[i+1]; k+=nData) {

//       PDM_g_num_t gFace     = requestedFace[k];

//       PDM_g_num_t _lFace = gFace - dFaceProc[myRank];
//       int          lFace = (int) _lFace;

//       int          nbVtxFace = (int) (dFaceVtxIdx[lFace+1]
//                                     - dFaceVtxIdx[lFace]);
//       sFaceInfoIdx[i+1] += nDataFace + nbVtxFace;
//     }
//   }

//   for (int i = 0; i < nRank; i++) {
//     sFaceInfoIdx[i+1] += sFaceInfoIdx[i];
//   }

//   PDM_g_num_t *sFaceInfo = (PDM_g_num_t *) calloc(sFaceInfoIdx[nRank], sizeof(PDM_g_num_t));

//   for (int i = 0; i < nRank; i++) {
//     for (int k = requestedFaceIdx[i]; k < requestedFaceIdx[i+1]; k+=nData) {
//       PDM_g_num_t gFace      = requestedFace[k];
//       PDM_g_num_t _lFace     = gFace - dFaceProc[myRank];
//       int          lFace     = (int) _lFace;
//       int          nbVtxFace = (int) (dFaceVtxIdx[lFace+1]
//                              - dFaceVtxIdx[lFace]);
//       int idx = sFaceInfoIdx[i] + sFaceInfoN[i];

//       sFaceInfo[idx++] = nbVtxFace;                   // Number of vertices
//       sFaceInfoN[i] += 1;

//       for(int j = dFaceVtxIdx[lFace]; j < dFaceVtxIdx[lFace+1]; j++) {
//         sFaceInfo[idx++] =  dFaceVtx[j];  //numero global du sommet qui compose la face
//         sFaceInfoN[i] += 1;
//       }
//     }
//   }

//   free(requestedFace);

//   PDM_g_num_t  *rFaceInfo    = NULL;
//   int          *rFaceInfoN   = requestedFaceN;
//   int          *rFaceInfoIdx = requestedFaceIdx;

//   _alltoall(           sFaceInfo,
//                        sFaceInfoN,
//                        sFaceInfoIdx,
//             (void **) &rFaceInfo,
//                        rFaceInfoN,
//                        rFaceInfoIdx,
//                        PDM__PDM_MPI_G_NUM,
//                        sizeof(PDM_g_num_t),
//                        comm);

//   free(sFaceInfo);
//   sFaceInfo = NULL;

//   int k = 0;
//   for (int i = 0; i < *sizeFaceVtxIdx - 1 ; i++) {
//     int nVtx = (int) rFaceInfo[k++];

//     faceVtxIdx[allToallNToLN[i]+1] = nVtx;
//     k += nVtx;
//   }

//   faceVtxIdx[0] = 0;
//   for (int i = 0; i < size_dCellFace; i++) {
//     faceVtxIdx[i+1] += faceVtxIdx[i];
//   }

//   *sizeFaceVtx = faceVtxIdx[*sizeFaceVtxIdx - 1];

//   (*faceVtx) = (PDM_g_num_t *) malloc ( *sizeFaceVtx * sizeof(PDM_g_num_t));

//   k = 0;
//   for (int i = 0; i < (*sizeFaceVtxIdx) -1 ; i++) {

//     int nVtx = (int) rFaceInfo[k++];
//     int idx = faceVtxIdx[allToallNToLN[i]];

//     for (int j = 0; j < nVtx; j++) {
//       (*faceVtx)[idx + j] = rFaceInfo[k++];
//     }
//   }

//   if (rFaceInfo != NULL) {
//     free(rFaceInfo);
//     rFaceInfo = NULL;
//   }

//   if (allToallNToLN != NULL)
//     free(allToallNToLN);

//   free(faceToSendN);
//   free(faceToSendIdx);
//   free(requestedFaceN);
//   free(requestedFaceIdx);

// }


// /**
//  *
//  * \brief Distributes vertex arrays
//  *
//  * \param [in]  comm         MPI Communicator
//  * \param [in]  dVtxCoord    coordinates of vertices
//  * \param [in]  dVtxProc     Vertex distribution
//  * \param [in]  sizeFaceVtx  Size of face vertex connectivity
//  * \param [in]  faceVtx      Face vertex connectivity
//  * \param [out] necessaryVtx Necessary vertices
//  *
//  */

// static void
// _distrib_vtx
// (
// const PDM_MPI_Comm      comm,
// const double       *dVtxCoord,
// const PDM_g_num_t *dVtxProc,
// const int           sizeFaceVtx,
// PDM_g_num_t       *faceVtx,
// double             *newVtx
// )
// {

//   const int nData    = 1;
//   int nDataVtx = 0;

//   int myRank;
//   int nRank;

//   PDM_MPI_Comm_rank(comm, &myRank);
//   PDM_MPI_Comm_size(comm, &nRank);

//   int  *vtxToSendIdx = (int *) calloc((nRank + 1) , sizeof(int));
//   int  *vtxToSendN   = (int *) calloc(nRank , sizeof(int));

//   int  *requestedVtxN   = (int *) calloc(nRank, sizeof(int));
//   int  *requestedVtxIdx = (int *) calloc((nRank + 1), sizeof(int));

//   int *allToallNToLN = (int *) calloc(sizeFaceVtx, sizeof(int));

//   for (int iVtx = 0; iVtx < sizeFaceVtx; ++iVtx) {
//     int irank =  PDM_binary_search_gap_long(faceVtx[iVtx], dVtxProc, nRank+1);
//     vtxToSendIdx[irank+1] += nData;
//   }

//   for (int i = 0; i < nRank; i++) {
//     vtxToSendIdx[i+1] += vtxToSendIdx[i] ;
//   }

//   PDM_g_num_t * vtxToSend = (PDM_g_num_t *) malloc(vtxToSendIdx[nRank] * sizeof(PDM_g_num_t));

//   for (int iVtx = 0; iVtx < sizeFaceVtx; iVtx++) {
//     int irank =  PDM_binary_search_gap_long (faceVtx[iVtx], dVtxProc, nRank+1);
//     int idx = vtxToSendIdx[irank] + vtxToSendN[irank];

//     allToallNToLN[idx/nData] = iVtx;
//     vtxToSend[idx++]   = faceVtx[iVtx];        /* Vtx global numbering */
//     vtxToSendN[irank] += nData;
//   }


//   PDM_g_num_t *requestedVtx    = NULL;

//   _alltoall(vtxToSend,
//             vtxToSendN,
//             vtxToSendIdx,
//             (void **) &requestedVtx,
//             requestedVtxN,
//             requestedVtxIdx,
//             PDM__PDM_MPI_G_NUM,
//             sizeof(PDM_g_num_t),
//             comm);
//   free(vtxToSend);

//   /*
//    *  Processes exchange information about requested vtxs
//    *  For each vtx, information contains :
//    *     - Tag (if ppart->_dVtxTag != NULL)
//    *     - Coordinates
//    *
//    */

//   int *sVtxInfoIdx = vtxToSendIdx; for( int i = 0; i < nRank+1; ++i) sVtxInfoIdx[i] = 0;
//   int *sVtxInfoN   = vtxToSendN; for (int i = 0; i < nRank; ++i)  sVtxInfoN[i] = 0;

//   for (int i = 0; i < nRank; i++) {
//     for (int k = requestedVtxIdx[i]; k < requestedVtxIdx[i+1]; k += nData) {
//       sVtxInfoIdx[i+1] += nDataVtx * (int) sizeof(int) + 3 * (int) sizeof(double);
//     }
//   }

//   for (int i = 0; i < nRank; i++) {
//     sVtxInfoIdx[i+1] += sVtxInfoIdx[i];
//   }

//   unsigned char *sVtxInfo = (unsigned char *)
//     malloc(sVtxInfoIdx[nRank] * sizeof(unsigned char));

//   for (int i = 0; i < nRank; i++) {
//     for (int k = requestedVtxIdx[i]; k < requestedVtxIdx[i+1]; k+=nData) {
//       PDM_g_num_t gVtx     = requestedVtx[k];
//       PDM_g_num_t _lVtx    = gVtx - dVtxProc[myRank];
//       int          lVtx    = (int) _lVtx;

//       int idx = sVtxInfoIdx[i] + sVtxInfoN[i];

//       double *_d_sVtxInfo = (double *) (sVtxInfo + idx);

//       for(int j = 0; j < 3; ++j) {
//         _d_sVtxInfo[j] = dVtxCoord[3*lVtx + j];
//       }

//       sVtxInfoN[i] += 3 * sizeof(double);
//     }
//   }

//   free(requestedVtx);

//   unsigned char *rVtxInfo    = NULL;
//   int           *rVtxInfoN   = requestedVtxN;
//   int           *rVtxInfoIdx = requestedVtxIdx;

//   _alltoall(sVtxInfo,
//       sVtxInfoN,
//       sVtxInfoIdx,
//       (void **) &rVtxInfo,
//       rVtxInfoN,
//       rVtxInfoIdx,
//       PDM_MPI_UNSIGNED_CHAR,
//       sizeof(PDM_MPI_UNSIGNED_CHAR),
//       comm);

//   if (sVtxInfo != NULL)
//     free(sVtxInfo);

//   /* Complete vtx */

//   assert(newVtx != NULL);

//   int k = 0;
//   for (int i = 0; i < sizeFaceVtx; i++) {
//     double *_d_rVtxInfo = (double *) (rVtxInfo + k);
//     for(int j = 0; j < 3; ++j)
//       newVtx[3*allToallNToLN[i] + j ] = _d_rVtxInfo[j];
//     k += 3*sizeof(double);
//   }

// /***********************************************/

//   if (allToallNToLN != NULL)
//     free(allToallNToLN);

//   if (rVtxInfo != NULL)
//     free(rVtxInfo);
//   rVtxInfo = NULL;

//   free(vtxToSendN);
//   free(vtxToSendIdx);
//   free(requestedVtxN);
//   free(requestedVtxIdx);

// }


/**
 *
 * \brief Call a couple MPI_Alltoall MPI_Alltoallv
 *
 * \param [in]   sendBuff,            Sending buffer
 * \param [in]   sendBuffN,           Number of data to send to each process
 *                                    (size : communicator size)
 * \param [in]   sendBuffIdx,         Index in sendBuff for each process
 *                                    (size : communicator size)
 * \param [out]  recvBuff,            Receiving buffer
 * \param [out]  recvBuffSize         Receiving buffer size
 * \param [in]   exch_mpi_data_type   Data type to exchange
 * \param [in]   type_exch_size       Size of data type
 * \param [in]   comm                 Communicator
 *
 */

// static void
// _alltoall
// (
//  void              *sendBuff,
//  int               *sendBuffN,
//  int               *sendBuffIdx,
//  void             **recvBuff,
//  int               *recvBuffN,
//  int               *recvBuffIdx,
//  PDM_MPI_Datatype   MPIDataType,
//  size_t             MPIDataTypeSize,
//  PDM_MPI_Comm       comm
// )
// {
//   int nRank = 0;
//   PDM_MPI_Comm_size(comm, &nRank);

//   /* Get number data to receive from each process */
//   PDM_MPI_Alltoall(sendBuffN,
//                    1,
//                    PDM_MPI_INT,
//                    recvBuffN,
//                    1,
//                    PDM_MPI_INT,
//                    comm);

//   recvBuffIdx[0] = 0;
//   for(int i = 0; i < nRank; i++) {
//       recvBuffIdx[i+1] = recvBuffIdx[i] + recvBuffN[i];
//   }

//   *recvBuff = malloc(recvBuffIdx[nRank] * MPIDataTypeSize);

//   /* Receive data from each process */
//   PDM_MPI_Alltoallv(sendBuff,
//                     sendBuffN,
//                     sendBuffIdx,
//                     MPIDataType,
//                     *recvBuff,
//                     recvBuffN,
//                     recvBuffIdx,
//                     MPIDataType,
//                     comm);

// }



#ifdef __cplusplus
}
#endif /* __cplusplus */
