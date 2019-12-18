/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_part_graph.h"
#include "pdm_hilbert.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_ext_wrapper.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */

static void
_quickSort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Splits the graph
 *
 * \param [in]  method       Method to be used: choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]  nPart        Number of partitions
 * \param [in]  part_ini     Part object fine mesh
 *
 * \param [in]  cellCell                  Dual graph (size : cellCellIdx[nCell])
 * \param [in]  cellCellIdx               Array of indexes of the dual graph (size : nCell + 1)
 * \param [in]  cellWeight         Cell weight (size = nCell)
 * \param [in]  faceWeight         Face weight (size = nFace)
 *
 * \param [inout] cellPart  Cell partitioning (size : nCell)
 *
 */

void
PDM_part_graph_split
(
 int         method,
 int         nPart,
 _part_t    *part_ini,
 int        *cellCellIdx,
 int        *cellCell,
 int        *cellWeight,
 int        *faceWeight,
 int       **cellPart
)
{
  *cellPart = (int *) malloc(part_ini->nCell * sizeof(int));

  for (int i = 0; i < part_ini->nCell; i++){
    (*cellPart)[i] = 0;
  }

  switch(method) {
  case 1:
    {
#ifdef PDM_HAVE_PARMETIS
      //            Define Metis properties

      //          int flag_weights = 0; //0 = False -> weights are unused

      int flag_weights = 0; //0 = False -> weights are unused

      int ncon = 1; //The number of balancing constraints

      int *vwgt = cellWeight; //Weights of the vertices of the graph (NULL if unused)

      int *adjwgt = faceWeight; //Weights of the edges of the graph (NULL if unused)

      double *tpwgts = NULL;
      if (flag_weights != 0) {
        tpwgts = (double *) malloc(ncon * nPart * sizeof(double));
        for (int i = 0; i < ncon * nPart; i++){
          tpwgts[i] = (double) (1./nPart);
        }
      }

      double *ubvec = NULL;
      if (flag_weights != 0) {
        ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }
      }

      //TO ADD: USE OF ADJWGT IN AN IF STATEMENT

      //This value is solely a memory space to be filled by METIS

      int edgecut;
      printf("PDM_part_graph_split \n");
      if (nPart < 8) {

        PDM_METIS_PartGraphRecursive (&(part_ini->nCell),
                                      &ncon,
                                      cellCellIdx,
                                      cellCell,
                                      vwgt,
                                      adjwgt,
                                      &nPart,
                                      tpwgts,
                                      ubvec,
                                      &edgecut,
                                      *cellPart);
      }

      else {

        PDM_METIS_PartGraphKway (&(part_ini->nCell),
                                 &ncon,
                                 cellCellIdx,
                                 cellCell,
                                 vwgt,
                                 adjwgt,
                                 &nPart,
                                 tpwgts,
                                 ubvec,
                                 &edgecut,
                                 *cellPart);
      }
      // double inbalance = 0.03;
      // double balance   = 0;
      // edgecut   = 0;
      // // bool   suppress_output = False;
      // // bool   graph_partitioned = False;
      // int time_limit = 0;
      // int seed  = 0;
      // int mode = 2;
      // PDM_kaffpa(&(part_ini->nCell),
      //            NULL,
      //            cellCellIdx,
      //            NULL,
      //            cellCell,
      //            &nPart,
      //             &inbalance,
      //             seed,
      //             mode,
      //             &edgecut,
      //             *cellPart );

      if (0 == 1) {
        PDM_printf("\n Contenu de cellPart : \n");
        for (int i = 0; i < part_ini->nCell; i++) {
          PDM_printf(" %d ", (*cellPart)[i]);
        }
        PDM_printf("\n");
      }

      if (flag_weights != 0) {
        if(ubvec!= NULL)
          free(ubvec);
        if(tpwgts!= NULL)
          free(tpwgts);
        // if(adjwgt!= NULL)
        //   free(adjwgt);
      }

#else
      PDM_printf("PDM_part error : METIS unavailable\n");
      exit(1);

#endif
      break;
    }
  case 2:
    {
#ifdef PDM_HAVE_PTSCOTCH

      int check = 0;

      PDM_SCOTCH_part (part_ini->nCell,
                       cellCellIdx,
                       cellCell,
                       cellWeight,
                       faceWeight,
                       check,
                       nPart,
                       *cellPart);

#else
      PDM_printf("PDM_part error : Scotch unavailable\n");
      exit(1);
#endif

      break;
    }
    case 3:
    {

      // To see with eric ...
      // abort();
      /* Allocation */
      double *cellCenter = (double *) malloc (part_ini->nCell * 3 * sizeof(double ));

      PDM_hilbert_code_t *hilbertCodes = (PDM_hilbert_code_t *) malloc (part_ini->nCell * sizeof(PDM_hilbert_code_t));

      /** Barycentre computation **/

      /* Allocate */
      double *cellPond = (double *) malloc (part_ini->nCell * sizeof(double));

      /* Nulliffy cellCenterArray */
      for(int iCell = 0; iCell < part_ini->nCell; iCell++) {
        cellCenter[3*iCell  ] = 0.;
        cellCenter[3*iCell+1] = 0.;
        cellCenter[3*iCell+2] = 0.;
        cellPond[iCell]     = 0.;
      }

      /* Compute */
      for(int iCell = 0; iCell < part_ini->nCell; iCell++) {

        /* Cellule composé de nFace */
        int aFac = part_ini->cellFaceIdx[iCell];
        int nFac = part_ini->cellFaceIdx[iCell+1] - aFac;

        for(int iFac = 0; iFac < nFac; iFac++) {

          /* Face composé de nVtx */
          int lFac = PDM_ABS(part_ini->cellFace[aFac + iFac]) - 1;

          int aVtx = part_ini->faceVtxIdx[lFac];
          int nVtx = part_ini->faceVtxIdx[lFac+1] - aVtx;

          for(int iVtx = 0; iVtx < nVtx; iVtx++) {

            /* Face composé de nVtx */
            int lVtx = part_ini->faceVtx[aVtx + iVtx] - 1;

            /* Add to current cell and stack weight */
            cellCenter[3*iCell  ] += part_ini->vtx[3*lVtx  ];
            cellCenter[3*iCell+1] += part_ini->vtx[3*lVtx+1];
            cellCenter[3*iCell+2] += part_ini->vtx[3*lVtx+2];

            cellPond[iCell] += 1.;
          }
        }
      }

      /* Nulliffy cellCenterArray */
      for(int iCell = 0; iCell < part_ini->nCell; iCell++) {
        cellCenter[3*iCell  ] = cellCenter[3*iCell  ]/cellPond[iCell];
        cellCenter[3*iCell+1] = cellCenter[3*iCell+1]/cellPond[iCell];
        cellCenter[3*iCell+2] = cellCenter[3*iCell+2]/cellPond[iCell];
      }


      double extents[3 * 2];

      /** Get EXTENTS LOCAL **/

      PDM_hilbert_get_coord_extents_seq(3, part_ini->nCell, cellCenter, extents);

      /** Hilbert Coordinates Computation **/

      PDM_hilbert_encode_coords(3, PDM_HILBERT_CS, extents, part_ini->nCell, cellCenter, hilbertCodes);

      /** CHECK H_CODES **/

      free(cellCenter);
      free(cellPond);

      int *newToOldOrder = (int *) malloc (part_ini->nCell * sizeof(int));
      for(int i = 0; i < part_ini->nCell; ++i) {
        newToOldOrder [i] = i;
      }

      PDM_sort_double (hilbertCodes, *cellPart, part_ini->nCell);

      /* Free */
      free (hilbertCodes);
      free (newToOldOrder);

      break;
    }
  default:
    PDM_printf("PART error : '%i' unknown partitioning method\n", method);
    exit(1);
  }

}


/**
 *
 * \brief Builds dual graph face cell connectivity
 *
 * \param [inout] part_ini                 Part object - fine mesh partition
 *
 * \param [inout] cellCellIdxCompressed    Array of indexes of the dual graph
 * \param [inout] cellCellCompressed       Dual graph
 *
 */
void
PDM_part_graph_compute_from_face_cell
(
  _part_t        *part_ini,
  int           **cellCellIdxCompressed,
  int           **cellCellCompressed
)
{
  int *cellCellN = (int *) malloc(part_ini->nCell * sizeof(int));
  for (int i = 0; i < part_ini->nCell; i++) {
    cellCellN[i] = 0;
  }

  int *cellCell = (int *) malloc(part_ini->cellFaceIdx[part_ini->nCell] * sizeof(int));
  for (int i = 0; i < part_ini->cellFaceIdx[part_ini->nCell]; i++) {
    cellCell[i] = -1;
  }

  int *cellCellIdx = (int *) malloc((part_ini->nCell + 1) * sizeof(int));
  for(int i = 0; i < part_ini->nCell + 1; i++) {
    cellCellIdx[i] = part_ini->cellFaceIdx[i];
  }

  for (int i = 0; i < part_ini->nFace; i++) {
    int iCell1 = PDM_ABS (part_ini->faceCell[2*i    ]);
    int iCell2 = PDM_ABS (part_ini->faceCell[2*i + 1]);
    //Only the non-boundary faces are stored
    if (iCell2 > 0) {
      int idx1 = cellCellIdx[iCell1-1] + cellCellN[iCell1-1];
      cellCell[idx1] = iCell2;
      cellCellN[iCell1-1] += 1;

      int idx2 = cellCellIdx[iCell2-1] + cellCellN[iCell2-1];
      cellCell[idx2] = iCell1;
      cellCellN[iCell2-1] += 1;
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cellCellN after looping over cellFace: ");
    for(int i = 0; i < part_ini->nCell; i++) {
      PDM_printf(" %d ", cellCellN[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cellCell after looping over cellFace: ");
    for(int i = 0; i < part_ini->cellFaceIdx[part_ini->nCell]; i++) {
      PDM_printf(" %d ", cellCell[i]);
    }
    PDM_printf("\n");
  }

  //cellCellIdx is rebuilt
  assert((*cellCellIdxCompressed) == NULL);
  (*cellCellIdxCompressed) = (int *) malloc((part_ini->nCell + 1) * sizeof(int));

  (*cellCellIdxCompressed)[0] = 0;
  for(int i = 0; i < part_ini->nCell; i++) {
    (*cellCellIdxCompressed)[i + 1] = (*cellCellIdxCompressed)[i] + cellCellN[i];
  }

  //We compress the dual graph since cellCellIdx was built from cellFaceIdx
  //We have then nFace elements in cellCell whereas it needs to be composed of nCell elements

  //    PDM_printf("(*cellCellIdxCompressed)[part_ini->nCell] : %d \n", (*cellCellIdxCompressed)[part_ini->nCell]);
  //
  assert( (*cellCellCompressed) == NULL);
  (*cellCellCompressed) = (int *) malloc((*cellCellIdxCompressed)[part_ini->nCell] * sizeof(int));

  int cpt_cellCellCompressed = 0;
  for(int i = 0; i < part_ini->cellFaceIdx[part_ini->nCell]; i++) {
    //        PDM_printf("I am testing a value for the %d time! \n", i);

    //We have an information to store when a neighboring cell exists
    if(cellCell[i] > -1){
      //            PDM_printf("I am storing a value for the %d time! \n", i);
      //We add a -1 to have the graph vertices numbered from 0 to n (C numbering)
      (*cellCellCompressed)[cpt_cellCellCompressed++] = cellCell[i] - 1;
      //            PDM_printf("Valeur stockee : %d \n ", (*cellCellCompressed)[cpt_cellCellCompressed - 1]);
    }
  }

  if( 0 == 1) {
    PDM_printf("Content of cellCellCompressed after compression and renumbering: ");
    for(int i = 0; i < (*cellCellIdxCompressed)[part_ini->nCell]; i++) {
      PDM_printf(" %d ", (*cellCellCompressed)[i]);
    }
    PDM_printf("\n");
  }

  /* Free temporary arrays*/

  free(cellCellN);
  free(cellCell);
  free(cellCellIdx);

  //Remove duplicate cells of the dual graph
  //We use the following scheme:
  //We loop over the indexes for the whole array to subdivide it into subarrays
  //We sort locally each subarray (determined thanks to cellCellIdxCompressed)
  //We loop over each subarray
  //We store the first value of each subarray anyway
  //We store each non-duplicated value and increment the writing index
  //We update the index array at each iteration

  int idx_write = 0;
  int tabIdxTemp = 0;

  for (int i = 0; i < part_ini->nCell; i++) {
    _quickSort_int((*cellCellCompressed), tabIdxTemp, (*cellCellIdxCompressed)[i + 1] - 1);

    int last_value = -1;

    for (int j = tabIdxTemp; j < (*cellCellIdxCompressed)[i + 1]; j++) {
      //We need to have a local index (between 0 and nFace)
      //If the value is different from the previous one (higher than is the same as different since the array is sorted)

      if(last_value != (*cellCellCompressed)[j]) {
        (*cellCellCompressed)[idx_write++] = (*cellCellCompressed)[j];
        last_value = (*cellCellCompressed)[j];
      }
    }

    if (0 == 1) {
      PDM_printf("\n Contenu de cellCellCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < (*cellCellIdxCompressed)[part_ini->nCell]; i1++) {
        PDM_printf(" %d ", (*cellCellCompressed)[i1]);
      }
      PDM_printf("\n");
    }

    tabIdxTemp = (*cellCellIdxCompressed)[i + 1];
    (*cellCellIdxCompressed)[i + 1] = idx_write;

    if (0 == 1) {
      PDM_printf("\n Contenu de cellCellIdxCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < part_ini->nCell + 1; i1++) {
        PDM_printf(" %d ", (*cellCellIdxCompressed)[i1]);
      }
      PDM_printf("\n");
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cellCellIdxCompressed after compression: ");
    for(int i1 = 0; i1 < part_ini->nCell + 1; i1++) {
      PDM_printf(" %d ", (*cellCellIdxCompressed)[i1]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cellCellCompressed after compression: ");
    for(int i1 = 0; i1 < (*cellCellIdxCompressed)[part_ini->nCell]; i1++) {
      PDM_printf(" %d ", (*cellCellCompressed)[i1]);
    }
    PDM_printf("\n");
  }

  //We reallocate the memory in case of duplicated values removed
  //The new array size is idx_write (stored in (*cellCellIdxCompressed)[part_ini->nCell])
  *cellCellCompressed = realloc(*cellCellCompressed,
                                (*cellCellIdxCompressed)[part_ini->nCell] * sizeof(int));

}




void
PDM_part_graph_split_bis
(
 int         method,
 int         nPart,
 int         graphSize,
 int        *cellCellIdx,
 int        *cellCell,
 int        *cellWeight,
 int        *faceWeight,
 int       **cellPart
)
{
  *cellPart = (int *) malloc(graphSize * sizeof(int));

  for (int i = 0; i < graphSize; i++){
    (*cellPart)[i] = 0;
  }

  switch(method) {
  case 1:
    {
#ifdef PDM_HAVE_PARMETIS
      //            Define Metis properties

      //          int flag_weights = 0; //0 = False -> weights are unused

      int flag_weights = 1; //0 = False -> weights are unused

      int ncon = 1; //The number of balancing constraints

      int *vwgt = cellWeight; //Weights of the vertices of the graph (NULL if unused)

      int *adjwgt = faceWeight; //Weights of the edges of the graph (NULL if unused)

      double *tpwgts = NULL;
      if (flag_weights != 0) {
        tpwgts = (double *) malloc(ncon * nPart * sizeof(double));
        for (int i = 0; i < ncon * nPart; i++){
          tpwgts[i] = (double) (1./nPart);
        }
      }

      double *ubvec = NULL;
      if (flag_weights != 0) {
        ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }
      }

      //TO ADD: USE OF ADJWGT IN AN IF STATEMENT

      //This value is solely a memory space to be filled by METIS

      int edgecut;

      if (nPart < 8) {

        PDM_METIS_PartGraphRecursive (&(graphSize),
                                      &ncon,
                                      cellCellIdx,
                                      cellCell,
                                      vwgt,
                                      adjwgt,
                                      &nPart,
                                      tpwgts,
                                      ubvec,
                                      &edgecut,
                                      *cellPart);
      }

      else {

        PDM_METIS_PartGraphKway (&(graphSize),
                                 &ncon,
                                 cellCellIdx,
                                 cellCell,
                                 vwgt,
                                 adjwgt,
                                 &nPart,
                                 tpwgts,
                                 ubvec,
                                 &edgecut,
                                 *cellPart);
      }

      // double inbalance = 0.03;
      // double balance   = 0;
      // edgecut   = 0;
      // // bool   suppress_output = False;
      // // bool   graph_partitioned = False;
      // int time_limit = 0;
      // int seed  = 0;
      // int mode = 2;
      // PDM_kaffpa(&(graphSize),
      //            NULL,
      //            cellCellIdx,
      //            NULL,
      //            cellCell,
      //            &nPart,
      //             &inbalance,
      //             seed,
      //             mode,
      //             &edgecut,
      //             *cellPart );

      if (0 == 1) {
        PDM_printf("\n Contenu de cellPart : \n");
        for (int i = 0; i < graphSize; i++) {
          PDM_printf(" %d ", (*cellPart)[i]);
        }
        PDM_printf("\n");
      }

      if (flag_weights != 0) {
        if(ubvec!= NULL)
          free(ubvec);
        if(tpwgts!= NULL)
          free(tpwgts);
        // if(adjwgt!= NULL)
        //   free(adjwgt);
      }

#else
      PDM_printf("PDM_part error : METIS unavailable\n");
      exit(1);

#endif
      break;
    }
  case 2:
    {
#ifdef PDM_HAVE_PTSCOTCH

      int check = 0;

      PDM_SCOTCH_part (graphSize,
                       cellCellIdx,
                       cellCell,
                       cellWeight,
                       faceWeight,
                       check,
                       nPart,
                       *cellPart);

#else
      PDM_printf("PDM_part error : Scotch unavailable\n");
      exit(1);
#endif

      break;
    }
    case 3:
    {

      // To see with eric ...
      abort();
      break;
    }
  default:
    PDM_printf("PART error : '%i' unknown partitioning method\n", method);
    exit(1);
  }

}


#ifdef  __cplusplus
}
#endif
