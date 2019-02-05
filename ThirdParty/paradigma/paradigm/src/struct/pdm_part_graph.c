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

      int flag_weights = 1; //0 = False -> weights are unused
          
      int ncon = 1; //The number of balancing constraints
            
      int *vwgt = cellWeight; //Weights of the vertices of the graph (NULL if unused)
                    
      int *adjwgt = faceWeight; //Weights of the edges of the graph (NULL if unused)

      if (flag_weights != 0) {
        double *tpwgts = (double *) malloc(ncon * nPart * sizeof(double));
        for (int i = 0; i < ncon * nPart; i++){
          tpwgts[i] = (double) (1./nPart);
        }
      }
          
      // double *tpwgts = NULL;
      double *tpwgts = (double *) malloc(ncon * nPart * sizeof(double));
          
      if (flag_weights != 0) {
        double *ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }
      }
          
      // double *ubvec = NULL;
      double *ubvec = (double *) malloc(ncon * sizeof(double));
          
      //TO ADD: USE OF ADJWGT IN AN IF STATEMENT                
        
      //This value is solely a memory space to be filled by METIS

      int edgecut;
          
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
          
      if (0 == 1) {
        PDM_printf("\n Contenu de cellPart : \n");            
        for (int i = 0; i < part_ini->nCell; i++) {
          PDM_printf(" %d ", (*cellPart)[i]);
        }
        PDM_printf("\n");
      }
        
      if (flag_weights != 0) {
          free(ubvec);
          free(tpwgts);
          free(adjwgt);
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
  *cellCellIdxCompressed = malloc((part_ini->nCell + 1) * sizeof(int)); 
    
  (*cellCellIdxCompressed)[0] = 0;
  for(int i = 0; i < part_ini->nCell; i++) {
    (*cellCellIdxCompressed)[i + 1] = (*cellCellIdxCompressed)[i] + cellCellN[i];
  }    
    
  //We compress the dual graph since cellCellIdx was built from cellFaceIdx
  //We have then nFace elements in cellCell whereas it needs to be composed of nCell elements
 
  //    PDM_printf("(*cellCellIdxCompressed)[part_ini->nCell] : %d \n", (*cellCellIdxCompressed)[part_ini->nCell]);
  *cellCellCompressed = malloc((*cellCellIdxCompressed)[part_ini->nCell] * sizeof(int));
    
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

#ifdef  __cplusplus
}
#endif
