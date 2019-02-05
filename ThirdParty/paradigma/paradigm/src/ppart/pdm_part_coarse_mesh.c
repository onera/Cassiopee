#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_geom_elem.h"
#include "pdm_part_coarse_mesh.h"
#include "pdm_part_coarse_mesh_priv.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"

#include "pdm_part.h"
#include "pdm_mpi.h"

#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_ext_wrapper.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"

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
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/


/**
 * \def _PDM_part_MIN(a,b)
 * Computes the minimum of \a x and \a y. 
 *
 */

#define _PDM_part_MIN(a,b) ((a) > (b) ? (b) : (a))

/**
 * \def _PDM_part_MAX(a,b)
 * Computes the maximum of \a x and \a y. 
 *
 */

#define _PDM_part_MAX(a,b) ((a) < (b) ? (b) : (a))

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

static PDM_Handles_t *_cm  = NULL;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * Storage of face renumbering methods
 */

static PDM_Handles_t *_coarse_mesh_methods = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/
  

    
/**
 *
 * \brief Return an initialized coarse part object
 * 
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   nPart             Number of partitions
 * \param [in]   nTPart            Total number of partitions
 * \param [in]   nFaceGroup        Number of boundaries
 * \param [in]   have_cellTag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_faceTag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtxTag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cellWeight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_faceWeight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_faceGroup    Presence des tableaux de groupes de faces
 */

static _coarse_mesh_t * 
_coarse_mesh_create
(
 const PDM_MPI_Comm  comm,        
 const char         *method,
 const int           nPart,
 const int           nTPart,
 const int           nFaceGroup,
 const int           have_cellTag,
 const int           have_faceTag,
 const int           have_vtxTag,
 const int           have_cellWeight,
 const int           have_faceWeight,
 const int           have_faceGroup

 )
{     
   _coarse_mesh_t *cm = (_coarse_mesh_t *) malloc(sizeof(_coarse_mesh_t));

   cm->nPart = nPart;
   cm->comm  = comm; 
   
   int _method = PDM_coarse_mesh_method_idx_get(method);
  
   if (_method == -1) {
     PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown coarse mesh method\n", method);
   }
   
   cm->method = _method;

   cm->nTPart = nTPart;
   
   cm->nFaceGroup = nFaceGroup;

   cm->have_cellTag    = have_cellTag;
   cm->have_faceTag    = have_faceTag;
   cm->have_vtxTag     = have_vtxTag;
   cm->have_cellWeight = have_cellWeight;
   cm->have_faceWeight = have_faceWeight;
   cm->have_faceGroup  = have_faceGroup;
   
   cm->part_ini = malloc(sizeof(_part_t *) * nPart); //On déclare un tableau de partitions
   
   cm->part_res = malloc(sizeof(_coarse_part_t *) * nPart);

   cm->specific_data = NULL;
   
   for (int i = 0; i < nPart; i++) {
     cm->part_ini[i] = _part_create(); 
     
     cm->part_res[i] = _coarse_part_create();     
     
   }   
    
   return cm;
}

/**
 *
 * \brief Perform the coarse mesh from the SCOTCH graph method
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_coarse_from_scotch
(
_coarse_mesh_t* cm,
const int       iPart,
int            *nCoarseCellComputed,
int            *cellCellIdx,
int            *cellCell,
int            *cellPart
)
{
#ifdef PDM_HAVE_PTSCOTCH

      int check = 0;
      int nPart  = cm->part_res[iPart]->nCoarseCellWanted;

      PDM_SCOTCH_part (cm->part_ini[iPart]->nCell,
                       cellCellIdx,
                       cellCell,
                       (int *) cm->part_ini[iPart]->cellWeight,
                       (int *) cm->part_ini[iPart]->faceWeight,
                       check,
                       nPart,
                       cellPart);

      (*nCoarseCellComputed) = nPart;

      if (0 == 1) {
        PDM_printf("\nContent of cellPart\n");
        for(int i = 0; i < cm->part_ini[iPart]->nCell ; i++) {
          PDM_printf(" %d ", cellPart[i]);
        }
        PDM_printf("\n");
      }

#else
      PDM_printf("PDM_part error : Scotch unavailable\n");
      exit(1);
#endif
}


/**
 *
 * \brief Perform the coarse mesh from the SCOTCH graph method
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_coarse_from_metis
(
_coarse_mesh_t* cm,
const int       iPart,
int            *nCoarseCellComputed,
int            *cellCellIdx,
int            *cellCell,
int            *cellPart
)
{
#ifdef PDM_HAVE_PARMETIS
      _part_t * part_ini       = cm->part_ini[iPart];
      _coarse_part_t *part_res = cm->part_res[iPart];

      int *cellWeight = (int *) part_ini->cellWeight;
      int *faceWeight = (int *) part_ini->faceWeight;

      int nPart  = part_res->nCoarseCellWanted;
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

      double *tpwgts = NULL;

      if (flag_weights != 0) {
        double *ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }
      }

      double *ubvec = NULL;

      //TO ADD: USE OF ADJWGT IN AN IF STATEMENT

      //This value is solely a memory space to be filled by METIS

      int edgecut;

      if (nPart < 8) {

        PDM_printf("\n \t\t\t\t PDM_METIS_PartGraphRecursive\n");
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
                                      cellPart);
      }

      else {

        PDM_printf("\n \t\t\t\tPDM_METIS_PartGraphKway \n");
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
                                 cellPart);
      }

      if (1 == 0) {
        PDM_printf("\n Contenu de cellPart : \n");
        for (int i = 0; i < part_ini->nCell; i++) {
          PDM_printf(" %d ", cellPart[i]);
        }
        PDM_printf("\n");
      }

      if (flag_weights != 0) {
          free(ubvec);
          free(tpwgts);
          free(adjwgt);
      }

      (*nCoarseCellComputed) = nPart;

#else
      PDM_printf("PDM_part error : METIS unavailable\n");
      exit(1);

#endif
}

/**
 *
 * \brief Return coarse mesh object from its identifier
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

static _coarse_mesh_t *
_get_from_id
(
 int  cmId
)
{
  _coarse_mesh_t *cm = (_coarse_mesh_t *) PDM_Handles_get (_cm, cmId);
      
  if (cm == NULL) {
    PDM_printf("PPART error : Bad cm identifier\n");
    exit(1);
  }

  return cm;
}


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

static void 
_dual_graph_from_face_cell
(
  _part_t        *part_ini,
  int           **cellCellIdxCompressed,
  int           **cellCellCompressed
)
{
  //cellCellN: array of counters of the numbers of connectivities
  //cellCell: dual graph to be built
  //cellCellIdx: array of indexes of the dual graph (same as cellFaceIdx)
    
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

static void 
_split
(
_coarse_mesh_t *cm,
const int       iPart,
int            *nCoarseCellComputed,
int            *cellCellIdx,
int            *cellCell,
int            **cellPart)
{    

  _part_t * part_ini       = cm->part_ini[iPart];
  
  int method = cm->method;

  *cellPart = (int *) malloc(part_ini->nCell * sizeof(int));
    
  for (int i = 0; i < part_ini->nCell; i++){
    (*cellPart)[i] = 0;    
  }
  
  const _coarse_mesh_method_t *method_ptr = (const _coarse_mesh_method_t *) 
                                      PDM_Handles_get (_coarse_mesh_methods, method);
  
  PDM_coarse_mesh_fct_t fct = method_ptr->fct;
  
  if (fct != NULL) {
    (fct) (cm,
           iPart,
           nCoarseCellComputed,
           cellCellIdx,
           cellCell,
           *cellPart);
  }

}

/**
 *
 * \brief Obtains the cells per processor from the partition of each cell
 * 
 * \param [in]  nCoarseCell  Number of partitions ( = number of coarse cells) 
 * \param [in]  nCell        Number of cells before refining
 * \param [in]  cellPart     Cell partitioning (size : nCell)
 * 
 * \param [inout] partCellIdx  Array of indexes of the partitioning array (size : nCoarseCell + 1)
 * \param [inout] partCell     Partitioning array (size : partCellIdx[nCoarseCell] = nCell)
 *
 */

static void 
_partCell_from_cellPart
(
 int            nCoarseCell,
 int            nCell,       
 int           *cellPart,
 int          **partCellIdx,
 int          **partCell
)
{
  //Allocation of an array to count the number of cells per partition
  int * cptCellsPerPartitions = (int *) malloc(nCoarseCell * sizeof(int));
  for (int i = 0; i < nCoarseCell; i++){
    cptCellsPerPartitions[i] = 0;
  }
    
  for (int i = 0; i < nCell; i++){
    int color = cellPart[i]; //A color is a number of partition (output of Metis or Scotch)
    cptCellsPerPartitions[color]++;
  }
    
  if(0 == 1) {
    PDM_printf("\n Contenu de cptCellsPerPartitions : \n");
    for(int i = 0; i < nCoarseCell; i++) {
      PDM_printf(" %d ", cptCellsPerPartitions[i]);
    }  
    PDM_printf("\n");
  }
    
  //Allocation of an array for counter indexes    
  *partCellIdx = (int *) malloc((nCoarseCell + 1) * sizeof(int));
  (*partCellIdx)[0] = 0;
  for (int i = 0; i < nCoarseCell; i++){
    (*partCellIdx)[i + 1] = (*partCellIdx)[i] + cptCellsPerPartitions[i];
  }
    
  if (0 == 1) {
    PDM_printf("\n Contenu de partCellIdx : \n");
    for(int i = 0; i < nCoarseCell + 1; i++) {
      PDM_printf(" %d ", (*partCellIdx)[i]);
    }
    PDM_printf("\n");
  }
    
  *partCell = (int *) malloc((*partCellIdx)[nCoarseCell] * sizeof(int));
    
  //cptCellsPerPartitions is reused for building partCell
  for (int i = 0; i < nCoarseCell; i++){
    cptCellsPerPartitions[i] = 0;
  }
    
  //We store each cell in partCell by means of (*partCellIdx)
  for (int i = 0; i < nCell; i++){
    int color = cellPart[i]; //A color is a number of partition (output of Metis or Scotch)
    int idx = (*partCellIdx)[color] + cptCellsPerPartitions[color];
    (*partCell)[idx] = i;
    cptCellsPerPartitions[color]++;
  }
    
  if (0 == 1) {
    PDM_printf("\nContenu de partCell \n");
    for (int i = 0; i < nCoarseCell; i++){
      PDM_printf("Valeur de i + 1 : %d \n", i + 1);
      for (int j = (*partCellIdx)[i]; j < (*partCellIdx)[i + 1]; j++){
        PDM_printf("%d " ,(*partCell)[j]);              
      }
      PDM_printf("\n");
    }
    PDM_printf("\n");
  }
        
  //Free
  free(cptCellsPerPartitions);
    
}

/**
 *
 * \brief Checks the neighboring cell of each studied cell by filling an array of tags
 *        A tag is set to -1 by default and set to numberGlobalPartition if the studied cell has valid neighbors
 * 
 * \param [in]  cellNumber                Number of the cell studied
 * \param [in]  cellPart                  Cell partitioning (size : nCell)
 * \param [in]  cellCell                  Dual graph (size : cellCellIdx[nCell])
 * \param [in]  cellCellIdx               Array of indexes of the dual graph (size : nCell + 1)
 * \param [in]  numberGlobalPartition     Number of the coarse cell to be applied
 * 
 * \param [inout] cellCoarseCell          Cell partitioning (size : nCell) (partitions are equal to coarse cells for a good partitioning)
 * \param [inout] cptCellConnectedLocal   Number of cells that have been tagged
 */

static void 
_fill_Neighboring
(
 int            cellNumber,      
 int           *cellPart,
 int           *cellCell,
 int           *cellCellIdx,
 int          **cellCoarseCell,
 int            numberGlobalPartition,
 int           *cptCellConnectedLocal
)
{
  /*
   * If the studied cell is not tagged, it is tagged
   */
  
  if ((*cellCoarseCell)[cellNumber] == -1) {
    (*cellCoarseCell)[cellNumber] = numberGlobalPartition;
    (*cptCellConnectedLocal)++;
  }
  
  /*
   * If the cell is tagged, I break
   */
  
  else {
    return;
  }

  /*
   * Loop over the partition cells (k is a neighbor cell number)
   */
  
  for (int k = cellCellIdx[cellNumber]; k < cellCellIdx[cellNumber + 1]; k++) {

    /*
     * if the neighboring cell is part of the same partition
     */
    
    if (cellPart[cellNumber] == cellPart[cellCell[k]]) {
      _fill_Neighboring (cellCell[k],
                         cellPart,
                         cellCell,
                         cellCellIdx,
                         cellCoarseCell,
                         numberGlobalPartition,
                         cptCellConnectedLocal);  
    }
  }
}


/**
 *
 * \brief Checks that the partitioning is fully connected
 * Otherwise, the cells that are not connected are removed for their initial partition to be part of a new one
 * partCell and partCellIdx are updated
 * 
 * \param [in]  nCoarseCellChecked  Number of partitions checked ( >= number of coarse cells wanted by the user) 
 * \param [in]  nCell               Number of cells before refining
 * \param [in]  cellPart            Cell partitioning (size : nCell) * 
 * \param [in]  cellCell            Dual graph (size : cellCellIdx[nCell])
 * \param [in]  cellCellIdx         Array of indexes of the dual graph (size : nCell + 1)
 * \param [in]  partCell            Partitioning array (size : partCellIdx[nCoarseCell] = nCell)
 * \param [in]  partCellIdx         Array of indexes of the partitions (size : nCoarseCellWanted + 1)
 *
 * \param [inout] coarseCellCellIdx  Array of indexes of the connected partitions (size : nCoarseCellChecked + 1)
 * \param [inout] coarseCellCell     Partitioning array (size : CoarseCellCellIdx[nCoarseCellChecked])
 * \param [inout] cellCoarseCell     Cell partitioning with coarse cells (size : nCoarseCellChecked)
 */

static void 
_adapt_Connectedness
(
 int           *nCoarseCellChecked,
 int            nCell,       
 int           *cellPart,
 int          **cellCoarseCell,
 int           *cellCell,
 int           *cellCellIdx,
 int           *partCell,
 int           *partCellIdx,
 int          **coarseCellCell,
 int          **coarseCellCellIdx 
)
{
  int numberGlobalPartition = 1;
    
  *cellCoarseCell = (int *) malloc(nCell * sizeof(int));
    
  for (int i = 0; i < nCell; i++) {
    (*cellCoarseCell)[i] = -1;
  }
    
  /*
   *  We store the initial number of coarse cells wanted by the user
   */

  int nCoarseCellWanted = (*nCoarseCellChecked);
 
  if (0 == 1) {
    PDM_printf("Valeur finale de (*nCoarseCellChecked) : %d %d\n", (*nCoarseCellChecked), nCell);

    PDM_printf("partCell : \n");
    for (int i = 0; i < nCoarseCellWanted; i++) {
      for (int j = partCellIdx[i]; j < partCellIdx[i+1]; j++) {
        PDM_printf(" %d ", partCell[j]);
      }
      PDM_printf("\n");
    }

    PDM_printf("cellCell : \n");
    for (int i = 0; i < nCell; i++) {
      for (int j = cellCellIdx[i]; j < cellCellIdx[i+1]; j++) {
        PDM_printf(" %d ", cellCell[j]);
      }
      PDM_printf("\n");
    }

  }

    
  /*
   * cellNumber will be replaced by the first cell of the first partition at the beginning of the loop
   */
  
  int cellNumber = -1;
    
  /* 
   * Loop over the partitions (i is a partition number)
   */
  
  for (int i = 0; i < nCoarseCellWanted; i++) {
    /*
     * We study the first cell of the partition
     */
    
    cellNumber = partCell[partCellIdx[i]];
    int cptCellConnectedLocal = 0;
        
    /*
     * We tag all the neighboring cells of the cell cellNumber of the partition
     */

    _fill_Neighboring(cellNumber, cellPart, cellCell, cellCellIdx, &(*cellCoarseCell), numberGlobalPartition, &cptCellConnectedLocal);
       
    numberGlobalPartition++;        
        
    /*
     * If the size of array indexes is too low
     */
        
    int nCellLocal = partCellIdx[i + 1] - partCellIdx[i];
    int nCellLocalRemaining = nCellLocal - cptCellConnectedLocal;
        
    /* 
     * If the partition has not been fully looped over, we will have to create an extra coarse cell
     */

    if (cptCellConnectedLocal < nCellLocal) {          
      /* 
       * We reinitialize cptCellConnectedLocal since we have a new coarse cell
       */
      cptCellConnectedLocal = 0;
    }
        
    /* 
     *  As long as the partition has not been fully looped over, we call the recursive function
     */

    while (cptCellConnectedLocal < nCellLocalRemaining) {
      for (int j = 1; j < nCellLocal; j++) {
        cellNumber = partCell[partCellIdx[i] + j];
        if ((*cellCoarseCell)[cellNumber] == -1) {
          break;
        }
      }
            
      _fill_Neighboring(cellNumber, cellPart, cellCell, cellCellIdx, &(*cellCoarseCell), numberGlobalPartition, &cptCellConnectedLocal);

      numberGlobalPartition++;

      /* 
       * If the size of array indexes is too low
       */

    }
         
  }
    
  (*nCoarseCellChecked) = numberGlobalPartition - 1;    
    
  /* 
   * Size of *coarseCellCellIdx may be dynamic
   */
      
  *coarseCellCellIdx = malloc(((*nCoarseCellChecked) + 1) * sizeof(int));
  for (int i = 0; i < (*nCoarseCellChecked) + 1; i++) {
    (*coarseCellCellIdx)[i] = 0;
  }

  for (int i = 0; i < nCell; i++) {
    (*coarseCellCellIdx)[(*cellCoarseCell)[i]]++;
  }

  for (int i = 0; i < (*nCoarseCellChecked); i++) {
    (*coarseCellCellIdx)[i+1] += (*coarseCellCellIdx)[i];
  }
  
  if (0 == 1) {
    PDM_printf("Valeur finale de (*nCoarseCellChecked) : %d %d\n", (*nCoarseCellChecked), nCell);
        
    PDM_printf("Affichage de *coarseCellCellIdx");
    for (int i = 0; i < (*nCoarseCellChecked) + 1; i++) {
      PDM_printf(" %d ", (*coarseCellCellIdx)[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cellCoarseCell: ");
    for (int i = 0; i < nCell; i++) {
      PDM_printf(" %d ", (*cellCoarseCell)[i]);
    }
    PDM_printf("\n");
  }
    
  /*
   * Creation of coarseCellCell from cellCoarseCell and cellCoarseCellIdx    
   */

  *coarseCellCell = (int *) malloc(nCell * sizeof(int));
    
  int * cptCellsPerPartitions = (int *) malloc((*nCoarseCellChecked) * sizeof(int));
  for (int i = 0; i < (*nCoarseCellChecked); i++){
    cptCellsPerPartitions[i] = 0;
  }      
    
  /*
   * We store each cell in partCell by means of (*partCellIdx)
   */

  for (int i = 0; i < nCell; i++){
    int color = (*cellCoarseCell)[i] - 1; //A color is a number of partition (output of Metis or Scotch)
    int idx = (*coarseCellCellIdx)[color] + cptCellsPerPartitions[color];
    (*coarseCellCell)[idx] = i + 1;
    cptCellsPerPartitions[color]++;
  }
    
  free(cptCellsPerPartitions);    
}

/**
 *
 * \brief Builds the array faceCoarseCell with all the inner faces removed
 * 
 * \param [in]  nFaceChecked                Number of faces after refining ( <= nFace)
 * \param [in]  faceCell                    Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [in]  cellCoarseCell              Cell partitioning with coarse cells (size : nCoarseCellChecked)
 *
 * \param [inout] faceCoarseCell            Face to coarse cell connectivity  (size = 2 * nFaceChecked, numbering : 1 to n)
 * \param [inout] fineFaceToCoarseFace      Fine face - coarse face connectivity (size = nFace)
 * \param [inout] coarseFaceToFineFace      Coarse face - fine face connectivity (size = nFaceChecked)
 *
 */

static void 
_build_faceCoarseCell
(
 int           *nFaceChecked,       
 int           *faceCell,
 int           *cellCoarseCell,
 int          **faceCoarseCell,
 int          **fineFaceToCoarseFace,
 int          **coarseFaceToFineFace
)
{
  /*
   * nFaceChecked = nFace at the beginning of the method 
   */
  
  int nFace = (*nFaceChecked);
    
  /*
   * Fine face - coarse face connectivity (size = nFace)
   */
  
  *fineFaceToCoarseFace = malloc(nFace * sizeof(int));
    
  int *faceCellTemp = (int *) malloc(2 * nFace * sizeof(int));
    
  for (int i = 0; i < nFace; i++) {
    (*fineFaceToCoarseFace)[i] = -1;
  }
    
  for (int i = 0; i < 2 * nFace; i++) {
    faceCellTemp[i] = 0;
  }
    
  /*
   * Loop over faceCell. i = number of face, faceCell[i] = cell number
   * We fill faceCellTemp with the coarse cells associated
   */
  
  for (int i = 0; i < 2 * nFace; i++) {

    /*
     * If we have a boarding face, we have a -1 -> nothing to do
     * If we have a "real" neighboring cell, we store its coarse cell
     */
    
    if (faceCell[i] != 0) {
      faceCellTemp[i] = cellCoarseCell[PDM_ABS (faceCell[i]) - 1];
    }
  }
    
  if (0 == 1) {
    PDM_printf("Content of faceCellTemp: |");
    for (int i = 0; i < 2 * nFace; i++) {
      PDM_printf(" %d ", faceCellTemp[i]);
      if (i % 2 == 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }
    
  *faceCoarseCell = (int *) malloc(2 * nFace * sizeof(int));
    
  /*
   * Loop over faceCellTemp which is to be compressed. i = face number    
   */

  int idx = 0;
  for (int i = 0; i < nFace; i++) {
    int iCell1 = PDM_ABS (faceCellTemp[2 * i    ]);
    int iCell2 = PDM_ABS (faceCellTemp[2 * i + 1]);        
        
    /* 
     * If a face is surrounded by the same coarse cell, it is not stored
     */
    if (iCell1 != iCell2) {
      (*faceCoarseCell)[2 * idx]     = iCell1;
      (*faceCoarseCell)[2 * idx + 1] = iCell2; 
            
      (*fineFaceToCoarseFace)[i] = idx + 1;            
      idx++;
    }
  }
        
  (*nFaceChecked) = idx;
  /* 
   * realloc of the correct size
   */

  *faceCoarseCell = realloc((*faceCoarseCell), 2 * (*nFaceChecked) * sizeof(int));
    
  /*
   * Fine face - coarse face connectivity (size = nFace)
   */

  *coarseFaceToFineFace = malloc(nFace * sizeof(int));
    
  int idx_coarseFaceToFineFace = 0;
    
  /*
   *  Loop over fineFaceToCoarseFace
   */

  for (int i = 0; i < nFace; i++) {

    /*
     * If the fine face has not been removed, I store it
     */
      
    if((*fineFaceToCoarseFace)[i] != -1) {
      (*coarseFaceToFineFace)[idx_coarseFaceToFineFace++] = i + 1;
    }
  }
    
  /*
   * At the end of the loop, idx_coarseFaceToFineFace must be equal to nFaceChecked
   */
    
  assert(idx_coarseFaceToFineFace == (*nFaceChecked));
    
  *coarseFaceToFineFace = realloc((*coarseFaceToFineFace), idx_coarseFaceToFineFace * sizeof(int));
    
  if(0 == 1) {
    PDM_printf("Valeur finale de (*nFaceChecked) : %d \n", (*nFaceChecked));

    PDM_printf("Final content of faceCoarseCell: |");
    for (int i = 0; i < 2 * (*nFaceChecked); i++) {
      PDM_printf(" %d ", (*faceCoarseCell)[i]);
      if (i % 2 == 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");

    PDM_printf("Final content of fineFaceToCoarseFace: \n");
    for (int i = 0; i < nFace; i++) {
      PDM_printf(" %d ", (*fineFaceToCoarseFace)[i]);
    }
    PDM_printf("\n");


    PDM_printf("Affichage final de (*coarseFaceToFineFace) \n");
    for (int i = 0; i < (*nFaceChecked); i++) {
      PDM_printf(" %d ", (*coarseFaceToFineFace)[i]);
    }
    PDM_printf("\n");
  }
        
  free(faceCellTemp);
}

/**
 *
 * \brief Obtains the faces per coarse cell from the coarse cell of each face
 * 
 * \param [in]  nCoarseCellChecked    Number of coarse cells after the connectedness check
 * \param [in]  nFaceChecked          Number of faces obtained after the creation of faceCoarseCell
 * \param [in]  faceCoarseCell        Face to coarse cell connectivity  (size = 2 * nFaceChecked, numbering : 1 to n)
 
 * \param [inout] coarseCellFaceIdx   Array of indexes of the coarse cell to face connectivity (size = nCoarseCellChecked + 1, numbering : 1 to n)
 * \param [inout] coarseCellFace      Coarse cell to face connectivity  (size = coarseCellFaceIdx[nCoarseCellChecked], numbering : 1 to n)
 *
 */

static void 
_coarseCellFace_from_faceCoarseCell
(
 int            nCoarseCellChecked,       
 int            nFaceChecked,
 int           *faceCoarseCell,
 int          **coarseCellFaceIdx,
 int          **coarseCellFace
)
{
  /* 
   *  Allocation of an array to count the number of faces per coarse cell
   */

  int * cptFacesPerCoarseCell = (int *) malloc(nCoarseCellChecked * sizeof(int));
  for (int i = 0; i < nCoarseCellChecked; i++) {
    cptFacesPerCoarseCell[i] = 0;
  }
    
  /*
   * Loop over faceCoarseCell. i = number of face
   */

  for (int i = 0; i < nFaceChecked; i++) {
    int coarseCell1 = faceCoarseCell[2 * i    ];
    int coarseCell2 = faceCoarseCell[2 * i + 1];
    cptFacesPerCoarseCell[coarseCell1 - 1]++;
    
    /* 
     * If coarseCell2 != -1, it is not a boarder cell
     * A non-boarder cell touches two coarse cells
     */

    if(coarseCell2 != 0) {
      cptFacesPerCoarseCell[coarseCell2 - 1]++;
    }
  }
    
  if(0 == 1) {
    PDM_printf("\n Contenu de cptFacesPerCoarseCell : \n");
    for(int i = 0; i < nCoarseCellChecked; i++) {
      PDM_printf(" %d ", cptFacesPerCoarseCell[i]);
    }  
    PDM_printf("\n");
  }
    
  /*
   * Allocation of an array for counter indexes
   */

  *coarseCellFaceIdx = (int *)malloc((nCoarseCellChecked + 1) * sizeof(int));
  (*coarseCellFaceIdx)[0] = 0;
  for (int i = 0; i < nCoarseCellChecked; i++) {
    (*coarseCellFaceIdx)[i + 1] = (*coarseCellFaceIdx)[i] + cptFacesPerCoarseCell[i];
  }
    
  *coarseCellFace = (int *) malloc((*coarseCellFaceIdx)[nCoarseCellChecked] * sizeof(int));

  /* 
   *  cptFacesPerCoarseCell is reused for building coarseCellFace
   */

  for (int i = 0; i < nCoarseCellChecked; i++){
    cptFacesPerCoarseCell[i] = 0;
  }
    
  /*
   * We store each face in coarseCellFace by means of (*coarseCellFaceIdx)
   * Loop over faceCoarseCell. i = number of face
   */
  
  for (int i = 0; i < nFaceChecked; i++) {
    int coarseCell1 = faceCoarseCell[2 * i]; 
    int coarseCell2 = faceCoarseCell[2 * i + 1]; 
        
    int idx1 = (*coarseCellFaceIdx)[coarseCell1 - 1] + cptFacesPerCoarseCell[coarseCell1 - 1];
    int idx2 = -1;  
        
    /*
     * If the face is not on the boarder, we store it
     */

    if (coarseCell2 != 0) {
      idx2 = (*coarseCellFaceIdx)[coarseCell2 - 1] + cptFacesPerCoarseCell[coarseCell2 - 1];          
    }
        
    (*coarseCellFace)[idx1] = i + 1;
    cptFacesPerCoarseCell[coarseCell1 - 1]++;
        
    /*
     * If idx2 is higher than -1, it means that the face is not on the boarder
     */

    if (idx2 > -1) {
      (*coarseCellFace)[idx2] = i + 1;
      cptFacesPerCoarseCell[coarseCell2 - 1]++;
    }
  }
    
  if(0 == 1) {
    PDM_printf("Contenu de (*coarseCellFace) \n");
    for (int i = 0; i < (*coarseCellFaceIdx)[nCoarseCellChecked]; i++) {
      PDM_printf(" %d ", (*coarseCellFace)[i]);
      if (i % (*coarseCellFaceIdx)[1] == (*coarseCellFaceIdx)[1] - 1) {
        PDM_printf("|");
      }
    }       

    PDM_printf("\n Contenu de (*coarseCellFaceIdx) : \n");
    for (int i = 0; i < nCoarseCellChecked + 1; i++) {
      PDM_printf(" %d ", (*coarseCellFaceIdx)[i]);
    }
    PDM_printf("\n");
  }
    
  free(cptFacesPerCoarseCell);
}

/**
 *
 * \brief Builds the array faceVtx with all the inner vertices removed
 * 
 * \param [in] nFace                 Number of faces before refining
 * \param [in] nFaceChecked          Number of faces after refining ( <= nFace)
 * \param [in] nVtx                  Number of vertices before refining
 * \param [in] fineFaceToCoarseFace  Fine face - coarse face connectivity (size = nFace)
 * 
 * \param [inout] faceVtxIdx         Face vertex connectivity index (final size = nFaceChecked + 1) 
 * \param [inout] faceVtx            Face vertex connectivity (final size = faceVtxIdx[nFaceChecked])
 * \param [inout] nVtxChecked        Number of vertices before refining becoming the number of vertices after refining
 * \param [inout] fineVtxToCoarseVtx Fine vertex - coarse vertex connectivity (size = nVtx)
 * \param [inout] coarseVtxToFineVtx Coarse vertex - fine vertex connectivity (size = nVtxChecked)
 *
 */

static void 
_build_faceVtx
(
 int            nFace,       
 int            nFaceChecked,
 int            nVtx, 
 int           *fineFaceToCoarseFace,
 int          **faceVtxIdx,
 int          **faceVtx,
 int           *nVtxChecked,
 int          **fineVtxToCoarseVtx,
 int          **coarseVtxToFineVtx
)
{    
  int idx_write_faceVtx = 0;
    
  (*faceVtxIdx)[0] = 0;
  int idx_write_faceVtxIdx = 1;    
    
  /*
   * Loop over the old faceVtxIdx, i = face number
   */
  
  for (int i = 0; i < nFace; i++) {
    //Loop over the old faceVtx, j = vertex number
    if (fineFaceToCoarseFace[i] != - 1) {

      for (int j = (*faceVtxIdx)[i]; j < (*faceVtxIdx)[i + 1]; j++) {            
        //If the face studied has been removed, I skip it
        int vtx = (*faceVtx)[j];                
        (*faceVtx)[idx_write_faceVtx++] = vtx; 
      }
        
      (*faceVtxIdx)[idx_write_faceVtxIdx] = (*faceVtxIdx)[i + 1] - (*faceVtxIdx)[i] + (*faceVtxIdx)[idx_write_faceVtxIdx - 1];
      idx_write_faceVtxIdx++;    
    }
  }
  
  *faceVtxIdx = realloc((*faceVtxIdx), (nFaceChecked + 1) * sizeof(int));
  *faceVtx = realloc((*faceVtx), (*faceVtxIdx)[nFaceChecked] * sizeof(int));     
  
  if (0 == 1) {
    PDM_printf("Valeur de (*faceVtxIdx)[nFaceChecked] : %d \n", (*faceVtxIdx)[nFaceChecked]);

    for (int i = 0; i < nFaceChecked; i++) {
      for (int j = (*faceVtxIdx)[i]; j < (*faceVtxIdx)[i + 1]; j++) {            
        //If the face studied has been removed, I skip it
        int vtx = (*faceVtx)[j];                
        // A supprimer
        for (int j1 = (*faceVtxIdx)[i]; j1 < (*faceVtxIdx)[i + 1]; j1++) {            
          //If the face studied has been removed, I skip it
          int vtx1 = (*faceVtx)[j1];                
          if (j != j1 && vtx == vtx1) {
            PDM_printf("Error multiple vertex in a face\n");
            abort();
          }
        }
      }
    }
    PDM_printf("\n");
  }
    
  /* 
   * Creation of a correspondence table coarse vertex to fine vertex
   */
    
  *coarseVtxToFineVtx = malloc((*faceVtxIdx)[nFaceChecked] * sizeof(int));
    
  int idx_write_coarseVtxToFineVtx = 0;
    
  /*
   * It is a copy of faceVtx at first
   * Then, it is sorted
   * All the double vertices from the sorted array are removed
   * We have our correspondence table
   */

  for (int i = 0; i < (*faceVtxIdx)[nFaceChecked]; i++) {
    (*coarseVtxToFineVtx)[i] = (*faceVtx)[i];
  }
    
  _quickSort_int((*coarseVtxToFineVtx), 0, (*faceVtxIdx)[nFaceChecked] - 1);  
    
  int last_value = -1;
    
  /*
   * Loop over (*coarseVtxToFineVtx)
   * Each vertex is stored only once
   */

  for (int i = 0; i < (*faceVtxIdx)[nFaceChecked]; i++) {
    if (last_value != (*coarseVtxToFineVtx)[i]) {
      (*coarseVtxToFineVtx)[idx_write_coarseVtxToFineVtx++] = (*coarseVtxToFineVtx)[i];
      last_value = (*coarseVtxToFineVtx)[i];
    }
  }
    
  (*nVtxChecked) = idx_write_coarseVtxToFineVtx;
    
  (*coarseVtxToFineVtx) = realloc((*coarseVtxToFineVtx), (*nVtxChecked) * sizeof(int));
    
  if (0 == 1) {
    PDM_printf("\nFinal content of coarseVtxToFineVtx: ");
    for (int i = 0; i < (*nVtxChecked); i++) {
      PDM_printf(" %d ", (*coarseVtxToFineVtx)[i]);     
    }
    PDM_printf("\n");
  }
    
  /*
   * Creation of a correspondence table fine vertex to coarse vertex
   */

  *fineVtxToCoarseVtx = malloc(nVtx * sizeof(int));
    
  for (int i = 0; i < nVtx; i++) {
    (*fineVtxToCoarseVtx)[i] = -1;
  }
        
  /*
   * Loop over (*coarseVtxToFineVtx)
   */

  for (int i = 0; i < (*nVtxChecked); i++) {
    int fineVtx = (*coarseVtxToFineVtx)[i];
    //        PDM_printf("Valeur de fineVtx : %d \n", fineVtx);
    (*fineVtxToCoarseVtx)[fineVtx - 1] = i + 1;
  }
    
  if(0 == 1) {
    PDM_printf("Content of fineVtxToCoarseVtx: ");
    for (int i = 0; i < nVtx; i++) {
      PDM_printf(" %d ", (*fineVtxToCoarseVtx)[i]);        
    }
    PDM_printf("\n");      
  }
        
  //Loop over faceVtx to re-number faceVtx thanks to (*coarseVtxToFineVtx)
  for (int i = 0; i < (*faceVtxIdx)[nFaceChecked]; i++) {
    (*faceVtx)[i] = (*fineVtxToCoarseVtx)[(*faceVtx)[i] - 1];        
  }   
    
  if (0 == 1) {
    PDM_printf("Valeur de (*nVtxChecked) : %d \n", (*nVtxChecked));    
    PDM_printf("Valeur de idx_write_faceVtx : %d \n", idx_write_faceVtx);

    PDM_printf("Final content of faceVtxIdx: ");
    for(int i = 0; i < nFaceChecked + 1; i++) {
      PDM_printf(" %d ", (*faceVtxIdx)[i]);
    }  
    PDM_printf("\n");

    PDM_printf("Final content of faceVtx: |");
    for (int i = 0; i < (*faceVtxIdx)[nFaceChecked]; i++) {
      PDM_printf(" %d ", (*faceVtx)[i]);
      if (i % (*faceVtxIdx)[1] == (*faceVtxIdx)[1] - 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }
        
}

/**
 *
 * \brief Builds the array vtx with all the coordinates of the inner vertices removed
 * 
 * \param [in]  nVtx                 Number of vertices before refining
 * \param [in]  nVtxChecked          Number of vertices after refining
 * \param [in]  fineVtxToCoarseVtx   Fine vertex - coarse vertex connectivity (size = nVtx)
 *
 * \param [inout] vtx                Vertex coordinates (size = nVtxChecked)
 * 
 */

static void 
_build_vtx
(
 int            nVtx,       
 int            nVtxChecked,
 int           *fineVtxToCoarseVtx,
 double       **vtx 
)
{
  //If no vertex has been removed, nothing to do!
  if (nVtx == nVtxChecked) {
      return;
  }

  int idx_write = 0;    

  //Loop over fineVtxToCoarseVtx, i = index of a vertex number (vertex number - 1)
  for (int i = 0; i < nVtx; i++) {
    //We store each vertex that has not been removed
    if (fineVtxToCoarseVtx[i] != -1) {
      double coord1 = (*vtx)[3 * i    ];
      double coord2 = (*vtx)[3 * i + 1];
      double coord3 = (*vtx)[3 * i + 2];

      (*vtx)[idx_write++] = coord1;
      (*vtx)[idx_write++] = coord2;
      (*vtx)[idx_write++] = coord3;
    }
  }

  //Reallocation of vtx at the suitable size
  *vtx = realloc((*vtx), 3 * nVtxChecked * sizeof(double));

  assert(3 * nVtxChecked == idx_write);

  if(0 == 1) {
    PDM_printf("Contenu final de vtx\n");
    for (int i = 0; i < 3 * nVtxChecked; i++) {        
      PDM_printf(" %.1f ", (*vtx)[i]);
      if (i % 3 == 2) {
          PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Updates the array cellTag in an array called coarseCellTag
 * 
 * \param [in] nCoarseCellChecked Number of partitions checked ( >= number of coarse cells wanted by the user) 
 * \param [in] coarseCellCellIdx  Array of indexes of the connected partitions (size : nCoarseCellChecked + 1)
 * \param [in] coarseCellCell     Partitioning array (size : coarseCellCellIdx[nCoarseCellChecked])
 * \param [in] cellTag            Cell tag (size = nCell)
 *                                  
 * \param [inout] coarseCellTag   Tag coarse cell connectivity index (size = nCoarseCellChecked)
 * 
 */

static void 
_build_coarseCellTag
(
 int            nCoarseCellChecked,       
 int           *coarseCellCellIdx,
 int           *coarseCellCell,
 int           *cellTag,
 int          **coarseCellTag
)
{
  if (cellTag == NULL) {
      return;
  }

  *coarseCellTag = (int *) malloc(nCoarseCellChecked * sizeof(int));

  //Loop over coarseCellCellIdx, i = index of coarse cell
  for (int i = 0; i < nCoarseCellChecked; i++) {
    //This should be the tag of all the fine cells of the coarse cell j
    int tag = cellTag[coarseCellCell[coarseCellCellIdx[i]] - 1];

    //Loop over coarseCellCell, j = coarse cell number
    for (int j = coarseCellCellIdx[i]; j < coarseCellCellIdx[i + 1]; j++) {
      //If any fine cell does not have the same tag as the previous one, the cellTag array is incorrect
      if (cellTag[coarseCellCell[j] - 1] != tag) {
        PDM_printf("Incorrect cellTag array provided!\n");
        PDM_printf("Please check the fine cell %d\n", j + 1);
        PDM_printf("A default tag of 0 will be written in the coarse cell %d\n", i + 1);
        tag = 0;
        break;
      }
      tag = cellTag[coarseCellCell[j] - 1];
    }
    (*coarseCellTag)[i] = tag;        
  }

  if(0 == 1) {
    PDM_printf("Affichage de (*coarseCellTag)\n");
    for (int i = 0; i < nCoarseCellChecked; i++) {
      PDM_printf(" %d ", (*coarseCellTag)[i]);
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Updates the array faceTag
 * 
 * \param [in]    nFaceChecked          Number of faces after refining ( <= nFace)
 * \param [in]    coarseFaceToFineFace  Coarse face - fine face connectivity (size = nFaceChecked)
 *
 * \param [inout] faceTag               Tag face connectivity index (size = nFace at first and nFaceChecked at the end)
 * 
 */

static void 
_build_faceTag
(
 int            nFaceChecked,       
 int           *coarseFaceToFineFace,
 int          **faceTag
)
{
  if(*faceTag == NULL) {
    return;
  }

  //Loop over coarseFaceToFineFace, i = number of a face after refinement
  for (int i = 0; i < nFaceChecked; i++) {
    (*faceTag)[i] = (*faceTag)[coarseFaceToFineFace[i] - 1];
  }

  (*faceTag) = realloc((*faceTag), nFaceChecked * sizeof(int));

  if(0 == 1) {
    PDM_printf("Contenu de (*faceTag)\n");
    for (int i = 0; i < nFaceChecked; i++) {
        PDM_printf(" %d ", (*faceTag)[i]);
    }
    PDM_printf("\n");
  }    
}

/**
 *
 * \brief Updates the array vtxTag
 * 
 * \param [in]    nVtxChecked          Number of vertices after refining
 * \param [in]    coarseVtxToFineVtx   Coarse vertex - fine vertex connectivity (size = nVtxChecked)
 *
 * \param [inout] vtxTag               Tag vertex connectivity index (size = nVtx at first and nVtxChecked at the end)
 * 
 */

static void 
_build_vtxTag
(
 int            nVtxChecked,       
 int           *coarseVtxToFineVtx,
 int          **vtxTag
)
{    
  if(*vtxTag == NULL) {
    return;
  }

  //Loop over coarseFaceToFineFace, i = number of a face after refinement
  for (int i = 0; i < nVtxChecked; i++) {
    (*vtxTag)[i] = (*vtxTag)[coarseVtxToFineVtx[i] - 1];
  }

  (*vtxTag) = realloc((*vtxTag), nVtxChecked * sizeof(int));

  if(0 == 1) {
    PDM_printf("Contenu de (*vtxTag)\n");
    for (int i = 0; i < nVtxChecked; i++) {
      PDM_printf(" %d ", (*vtxTag)[i]);
    }
    PDM_printf("\n");
  }    
}

/**
 *
 * \brief Updates the array faceGroup by renumbering the faces and removing the removed faces
 * 
 * \param [in] nFaceGroup            Number of groups of faces
 * \param [in] fineFaceToCoarseFace  Fine face - coarse face connectivity (size = nFace)
 * 
 * \param [inout] faceGroup          Face group index (size = faceGroupIdx[nFaceGroup]) 
 * \param [inout] faceGroupIdx       Face group index (size = nFaceGroup + 1)
 *
 */

static void 
_build_faceGroup
(
 int            nFaceGroup,       
 int          **faceGroup, 
 int          **faceGroupIdx,
 int          **coarseFaceGroupToFineFaceGroup
)
{ 
  
  if(*faceGroup == NULL || *faceGroupIdx == NULL || nFaceGroup == 0) {
    return;
  }
  *coarseFaceGroupToFineFaceGroup = malloc((*faceGroupIdx)[nFaceGroup] * sizeof(int));
    
  //Renumbering of partGroup from the fine numbering to the coarse one
  //Loop over faceGroup, i = face number
//  for (int i = 0; i < (*faceGroupIdx)[nFaceGroup]; i++) {
//      (*faceGroup)[i] = fineFaceToCoarseFace[(*faceGroup)[i] - 1];
//  }
    
  if (0 == 1) {
    PDM_printf("Content of faceGroup after renumbering: |");
    for (int i = 0; i < (*faceGroupIdx)[nFaceGroup]; i++) {
      PDM_printf(" %d ", (*faceGroup)[i]);
    }
    PDM_printf("\n");
  }
    
  int idx = 0;
    
  //Counter of faces per group
  int *cptFacesPerGroup = malloc(nFaceGroup * sizeof(int));
    
  for (int i = 0; i < nFaceGroup; i++) {
    cptFacesPerGroup[i] = 0;
  }
    
  //faceGroupIdx is rebuilt
  //Loop over faceGroupIdx, i = group number
  for (int i = 0; i < nFaceGroup; i++) {
    int startNumberingFace = 1;
    //Loop over faceGroup, j = face number
    for (int j = (*faceGroupIdx)[i]; j < (*faceGroupIdx)[i + 1]; j++) {
      //If we do not have a -1, the face has not been removed and is saved
      if ((*faceGroup)[j] != -1) {
        cptFacesPerGroup[i]++;
        (*faceGroup)[idx] = (*faceGroup)[j];              
        (*coarseFaceGroupToFineFaceGroup)[idx++] = startNumberingFace;               
      }
      startNumberingFace++;
    }
  }

  if (0 == 1) {
    PDM_printf("Contenu de cptFacesPerGroup apres remplissage\n");
    for (int i = 0; i < nFaceGroup; i++) {
      PDM_printf(" %d ", cptFacesPerGroup[i]);
    }
    PDM_printf("\n");
  }
    
  //Update of faceGroupIdx
  (*faceGroupIdx)[0] = 0;
  
  //Loop over cptFacesPerGroup, i = group number
  for (int i = 0; i < nFaceGroup; i++) {
    (*faceGroupIdx)[i + 1] = (*faceGroupIdx)[i] + cptFacesPerGroup[i];
  }
    
  (*faceGroup) = realloc((*faceGroup), (*faceGroupIdx)[nFaceGroup] * sizeof(int));
  (*coarseFaceGroupToFineFaceGroup) = realloc((*coarseFaceGroupToFineFaceGroup), (*faceGroupIdx)[nFaceGroup] * sizeof(int));
    
  if (0 == 1) {
    PDM_printf("Final content of faceGroupIdx: ");
    for (int i = 0; i < nFaceGroup + 1; i++) {
      PDM_printf(" %d ", (*faceGroupIdx)[i]);
    }
    PDM_printf("\n");
    
    PDM_printf("Final content of faceGroup: |");
    for (int i = 0; i < (*faceGroupIdx)[nFaceGroup]; i++){
      PDM_printf(" %d ", (*faceGroup)[i]);        
    }
    PDM_printf("\n");

    PDM_printf("Final content of coarseFaceGroupToFineFaceGroup: ");
    for (int i = 0; i < (*faceGroupIdx)[nFaceGroup]; i++) {
      PDM_printf(" %d ", (*coarseFaceGroupToFineFaceGroup)[i]);        
    }
    PDM_printf("\n\n");
  }
  free(cptFacesPerGroup);
}

/**
 *
 * \brief SetUp data array fo r coarse mesh
 *
 * \param [out] cgId              Coarse grid identifier
 * 
 * \param [in]  iPart              Partition identifier
 * \param [in]  nCoarseCellWanted  Number of cells in the coarse grid wanted by the user
 * \param [in]  nCell              Number of cells
 * \param [in]  nFace              Number of faces
 * \param [in]  nVtx               Number of vertices
 * \param [in]  nFaceGroup         Number of face groups
 * \param [in]  nFacePartBound     Number of partitioning boundary faces            
 * \param [in]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1, numbering : 0 to n-1)
 * \param [in]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace
 *                                                             numbering : 1 to n)
 * \param [in]  cellTag            Cell tag (size = nCell)
 * \param [in]  cellWeight         Cell weight (size = nCell)
 * \param [in]  faceWeight         Face weight (size = nFace)
 * \param [in]  cellLNToGN         Cell local numbering to global numbering (size = nCell, numbering : 1 to n)
 * \param [in]  faceCell           Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [in]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1, numbering : 0 to n-1)
 * \param [in]  faceVtx            Face to Vertex connectivity (size = faceVertexIdx[nFace], numbering : 1 to n)
 * \param [in]  faceTag            Face tag (size = nFace)
 * \param [in]  faceLNToGN         Face local numbering to global numbering (size = nFace, numbering : 1 to n)
 * \param [in]  vtxCoord           Vertex coordinates (size = 3 * nVertex)
 * \param [in]  vtxTag             Vertex tag (size = nVertex)
 * \param [in]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx, numbering : 1 to n)
 * \param [in]  faceGroupIdx       Face group index (size = nFaceGroup + 1, numbering : 1 to n-1)
 * \param [in]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [in]  faceGroupLNToGN    Faces global numbering for each group 
 *                                  (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [in]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
 * \param [in]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
 * \param [in]  facePartBound      Partitioning boundary faces (size = 4 * nFacePartBound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition 
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number 
 *                                          in the connected partition (numbering :1 to n)
 */

static void 
_coarse_grid_mesh_input
( 
 _coarse_mesh_t     *cm,
 const int           iPart,       
 const int           nCoarseCellWanted,
 const int           nCell,
 const int           nFace,
 const int           nVtx,
 const int           nFaceGroup,
 const int           nFacePartBound,
 const int          *cellFaceIdx,
 const int          *cellFace,
 const int          *cellTag,
 const int          *cellWeight,
 const int          *faceWeight,
 const PDM_g_num_t  *cellLNToGN,
 const int          *faceCell,
 const int          *faceVtxIdx,
 const int          *faceVtx,
 const int          *faceTag,
 const PDM_g_num_t  *faceLNToGN,
 const double       *vtxCoord,
 const int          *vtxTag,
 const PDM_g_num_t  *vtxLNToGN,
 const int          *faceGroupIdx,
 const int          *faceGroup,
 const PDM_g_num_t  *faceGroupLNToGN,
 const int          *facePartBoundProcIdx,       
 const int          *facePartBoundPartIdx,
 const int          *facePartBound      
)
{
  _part_t * part_ini = cm->part_ini[iPart];
  _coarse_part_t *part_res = cm->part_res[iPart];

  // const int *_cellWeight = cellWeight;
  // const int *_faceWeight = faceWeight;
  
  part_ini->cellWeight = cellWeight;
  part_ini->faceWeight = faceWeight;
  
  part_ini->nVtx = nVtx;
  part_ini->nCell = nCell;
  part_ini->nFace = nFace;
  part_ini->nFaceGroup = nFaceGroup;
  part_ini->nFacePartBound = nFacePartBound;
  part_ini->cellFaceIdx = (int *) cellFaceIdx;
  part_ini->cellFace = (int *) cellFace;
  part_ini->cellTag = (int *) cellTag;
  part_ini->faceCell = (int *) faceCell;
  part_ini->faceVtxIdx = (int *) faceVtxIdx; 
  part_ini->faceVtx = (int *) faceVtx;
  part_ini->faceTag = (int *) faceTag;
  part_ini->facePartBoundProcIdx = (int *) facePartBoundProcIdx;
  part_ini->facePartBoundPartIdx = (int *) facePartBoundPartIdx;
  part_ini->facePartBound = (int *) facePartBound;
  part_ini->faceGroupIdx = (int *) faceGroupIdx;
  part_ini->faceGroup = (int *) faceGroup;
  part_ini->vtx = (double *) vtxCoord;
  part_ini->vtxTag = (int *) vtxTag;    
  
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
  part_ini->cellLNToGN = (PDM_g_num_t *) cellLNToGN;
  part_ini->faceLNToGN = (PDM_g_num_t *) faceLNToGN;
  part_ini->faceGroupLNToGN = (PDM_g_num_t *) faceGroupLNToGN;
  part_ini->vtxLNToGN = (PDM_g_num_t *) vtxLNToGN;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
  
  part_res->nCoarseCellWanted = nCoarseCellWanted;
  
}



/**
 *
 * \brief Build a coarse grid prealably setUp with _coarse_grid_mesh_input
 *
 * \param [out] cgId              Coarse grid identifier
 * 
 * \param [in]  iPart              Partition identifier
 * \param [in]  nCoarseCellWanted  Number of cells in the coarse grid wanted by the user
 */

static void 
_coarse_grid_compute
( 
 _coarse_mesh_t     *cm,
 const int           iPart
)
{

  _part_t * part_ini = cm->part_ini[iPart];
  _coarse_part_t *part_res = cm->part_res[iPart];

  

  cm->timer = PDM_timer_create();
  for (int i = 0; i < 18; i++) {
    cm->times_elapsed[i] = 0.;
    cm->times_cpu[i] = 0.;
    cm->times_cpu_u[i] = 0.;
    cm->times_cpu_s[i] = 0.;
  }
  
  PDM_timer_resume(cm->timer);  
   
  int *dualGraphIdx = NULL;
  int *dualGraph    = NULL;  
  
  _dual_graph_from_face_cell(part_ini,
                             (int **) &dualGraphIdx,
                             (int **) &dualGraph);  
  
  int itime = 1;
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  //Call Metis or Scotch to get the cellPart array
  //cellPart must be allocated before proceeding (the initialization is part of the split method)
  
  PDM_timer_resume(cm->timer);
  
  int *cellPart = NULL;   
  
  int nCoarseCellComputed;

  _split( cm, 
          iPart,
         &nCoarseCellComputed,
         dualGraphIdx, 
         dualGraph, 
         (int **) &cellPart);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  
  /* Assign size of multigrid */
  part_res->part->nCell = nCoarseCellComputed;
  
  //  From the cellPart array, get the partCell
  int *partCellIdx = NULL;
  int *partCell = NULL;  
  
  // _partCell_from_cellPart(nCoarseCellWanted,
  _partCell_from_cellPart(part_res->part->nCell,
                          part_ini->nCell, 
                          cellPart, 
                          (int **) &partCellIdx,
                          (int **) &partCell);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  //Check that all partitions are correctly connected

  PDM_timer_resume(cm->timer);
  
  int *cellCoarseCell = NULL; 
  
  // part_res->part->nCell = nCoarseCellWanted;
    
  _adapt_Connectedness(&(part_res->part->nCell),
                       part_ini->nCell,
                       cellPart,
                       (int **) &cellCoarseCell,
                       dualGraph,
                       dualGraphIdx, 
                       partCell, 
                       partCellIdx,
                       (int **) &(part_res->coarseCellCell),
                       (int **) &(part_res->coarseCellCellIdx));
  
  free(partCellIdx);
  free(partCell);
  
  free(cellPart);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  //Compress the faceCell array to create the faceCoarseCell array
  
  PDM_timer_resume(cm->timer);  
  int *fineFaceToCoarseFace = NULL;
  
  //Temporary storage of the data of part_ini
  part_res->part->faceCell = part_ini->faceCell;
  part_res->part->nFace = part_ini->nFace;
  
  _build_faceCoarseCell(&(part_res->part->nFace),
                        part_ini->faceCell,
                        cellCoarseCell,
                        &(part_res->part->faceCell),
                        (int **) &fineFaceToCoarseFace,
                        &(part_res->coarseFaceToFineFace));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
    
  //Updates the faceGroupIdx and faceGroup arrays  
  PDM_timer_resume(cm->timer);
  
  part_res->part->faceGroupIdx = NULL;
  part_res->part->faceGroup = NULL;
  part_res->part->faceGroupLNToGN = NULL;

  if (part_ini->nFaceGroup > 0) {
    part_res->part->faceGroupIdx = malloc((part_ini->nFaceGroup + 1) * sizeof(int));
    for (int i = 0; i < (part_ini->nFaceGroup + 1); i++) {
      part_res->part->faceGroupIdx[i] = part_ini->faceGroupIdx[i];
    }
    part_res->part->faceGroup = malloc(part_res->part->faceGroupIdx[part_ini->nFaceGroup] * sizeof(int));
    for (int i = 0; i < part_ini->faceGroupIdx[part_ini->nFaceGroup]; i++) {
      part_res->part->faceGroup[i] = fineFaceToCoarseFace[part_ini->faceGroup[i] - 1]; 
    }
  }
  
  _build_faceGroup(part_ini->nFaceGroup,
                   &(part_res->part->faceGroup),
                   &(part_res->part->faceGroupIdx),
                   &(part_res->coarseFaceGroupToFineFaceGroup));
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Compress the cellFaceIdx and cellFace arrays
  
  PDM_timer_resume(cm->timer);
  
  _coarseCellFace_from_faceCoarseCell(part_res->part->nCell,
                                      part_res->part->nFace, 
                                      part_res->part->faceCell, 
                                      &(part_res->part->cellFaceIdx),
                                      &(part_res->part->cellFace));
   
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Compress the faceVtxIdx and faceVtx arrays
  
  PDM_timer_resume(cm->timer);
  
  part_res->part->nVtx = part_ini->nVtx;
  
  part_res->part->faceVtxIdx = malloc((part_ini->nFace + 1) * sizeof(int));
  for (int i = 0; i < (part_ini->nFace + 1); i++) {
    part_res->part->faceVtxIdx[i] = part_ini->faceVtxIdx[i];
  }
  
  part_res->part->faceVtx = malloc(part_res->part->faceVtxIdx[part_ini->nFace] * sizeof(int));
  for (int i = 0; i < part_res->part->faceVtxIdx[part_ini->nFace]; i++) {
    part_res->part->faceVtx[i] = part_ini->faceVtx[i];
  }
    
  int *fineVtxToCoarseVtx = NULL;
  
  _build_faceVtx(part_ini->nFace,
                 part_res->part->nFace, 
                 part_ini->nVtx, 
                 fineFaceToCoarseFace, 
                 &(part_res->part->faceVtxIdx),
                 &(part_res->part->faceVtx), 
                 &(part_res->part->nVtx), 
                 (int **) &fineVtxToCoarseVtx,
                 (int **) &(part_res->coarseVtxToFineVtx));

  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  //  Compress the vtxCoord array

  PDM_timer_resume(cm->timer);
  
  part_res->part->vtx = malloc(3 * part_ini->nVtx * sizeof(double));
  for (int i = 0; i < 3 * part_ini->nVtx; i++) {
    part_res->part->vtx[i] = part_ini->vtx[i];
  }
  
  _build_vtx(part_ini->nVtx, 
             part_res->part->nVtx, 
             fineVtxToCoarseVtx, 
             &(part_res->part->vtx));
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  //Update the tag arrays
  
  PDM_timer_resume(cm->timer);
  
  _build_coarseCellTag(part_res->part->nCell, 
                       part_res->coarseCellCellIdx,
                       part_res->coarseCellCell, 
                       part_ini->cellTag, 
                       &(part_res->part->cellTag));
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  
  if (part_ini->faceTag != NULL) {
    part_res->part->faceTag = malloc(part_ini->nFace * sizeof(int));
    for (int i = 0; i < part_ini->nFace; i++) {
      part_res->part->faceTag[i] = part_ini->faceTag[i];
    }
  }
  else {
    part_res->part->faceTag = NULL;
  }
  
  _build_faceTag(part_res->part->nFace, 
                 part_res->coarseFaceToFineFace, 
                 &(part_res->part->faceTag));
 
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  if (part_ini->vtxTag != NULL) {
    part_res->part->vtxTag = malloc(part_ini->nVtx * sizeof(int));
    for (int i = 0; i < part_ini->nVtx; i++) {
      part_res->part->vtxTag[i] = part_ini->vtxTag[i];
    }
  }
  else {
    part_res->part->vtxTag = NULL;
  }  
  
  _build_vtxTag(part_res->part->nVtx,
                part_res->coarseVtxToFineVtx, 
                &(part_res->part->vtxTag));
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  free(cellCoarseCell);
  
  free(dualGraphIdx);
  free(dualGraph);
  
  free(fineFaceToCoarseFace);

  free(fineVtxToCoarseVtx);
}

/**
 *
 * \brief Updates the cellLNToGN array for the coarse mesh are save it in the coarse mesh partition
 * 
 * \param [in]   cm                 Coarse mesh
 */

static void
_build_coarseCellLNToGN
(
 _coarse_mesh_t * cm
 )
{    
  //Calculation of the number of cells on the processor
  PDM_g_num_t nCellProc = 0;
  
  //Loop over the partition numbers, i = partition number
  for (int i = 0; i < cm->nPart; i++) {        
    nCellProc += cm->part_res[i]->part->nCell;
  }    
  
  //    PDM_printf("\nValeur de nCellProc : %d \n", nCellProc);
  
  //Global numbering of the cells
  PDM_g_num_t beg_NumAbs;
  
  PDM_MPI_Scan(&nCellProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);
  
  //Index to position the local cells
  beg_NumAbs -= nCellProc;
  
  int idx_write = 0;    
  
  //Loop over the partition numbers, i = partition number
  for (int i = 0; i < cm->nPart; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    int nCell = cm->part_res[i]->part->nCell;
    cmp->cellLNToGN = (PDM_g_num_t *) malloc(nCell * sizeof(PDM_g_num_t));        
    //Loop over the partition cells, j = cell number
    for (int j = 0; j < nCell; j++) {
      cmp->cellLNToGN[j] = beg_NumAbs + idx_write + 1;
      idx_write++;
    }
    
  }  
  
  if(0 == 1) {
    for (int iPart = 0; iPart < cm->nPart; iPart++) {
      PDM_printf("\nContenu de cm->part_res[%d]->part->cellLNToGN\n", iPart);
      for (int j = 0; j < cm->part_res[iPart]->part->nCell; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cm->part_res[iPart]->part->cellLNToGN[j]);
      }
      PDM_printf("\n\n");   
    }
  }
}

/**
 *
 * \brief Updates the faceLNToGN array for the coarse mesh are save it in the coarse mesh partition
 * 
 * \param [in]   cm                 Coarse mesh 
 *
 */

static void
_build_faceLNToGN
(
_coarse_mesh_t * cm
)
{
  PDM_g_num_t **faceLNToGNPart = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));
  int *nFacePart = (int *) malloc(cm->nPart * sizeof(int));

  for (int i = 0; i < cm->nPart; i++) {
    faceLNToGNPart[i] = cm->part_ini[i]->faceLNToGN;
    nFacePart[i] = cm->part_ini[i]->nFace;
  }

  if(0 == 1) {
    PDM_printf("Contenu de faceLNToGNPart\n");
    for (int i = 0; i < cm->nPart; i++) {
      for (int j = 0; j < nFacePart[i]; j++) {
         PDM_printf(" "PDM_FMT_G_NUM" ", faceLNToGNPart[i][j]);
      }
    PDM_printf("\n");
    }

  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCk_POST_CLEANUP,
                                                     1.,
                                                     (PDM_g_num_t **) faceLNToGNPart,
                                                     nFacePart,
                                                     cm->nPart,
                                                     cm->comm);    

  PDM_g_num_t **faceLNToGNTag = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));

  int idx_write = 0;

  for (int i = 0; i < cm->nPart; i++) {
    idx_write = 0;
    faceLNToGNTag[i] = (PDM_g_num_t *) malloc(cm->part_ini[i]->nFace * sizeof(PDM_g_num_t));
      //Loop over coarseFaceToFineFace, i = index of coarseFaceToFineFace (from 0 to cm->part_res[iPart]->part->nFace)

    for (int j = 0; j < cm->part_ini[i]->nFace; j++) {
      faceLNToGNTag[i][j] = -1;
    }

    for (int j = 0; j < cm->part_res[i]->part->nFace; j++) {
          //If the vertex studied is the same as in coarseFaceToFineFace, it is to be stored
     int k =  cm->part_res[i]->coarseFaceToFineFace[j] - 1;
     faceLNToGNTag[i][k] = 0;
    }
  }

  if(0 == 1) {
    PDM_printf("Contenu de faceLNToGNTag\n");    
    for (int i = 0; i < cm->nPart; i++) {
      for (int j = 0; j < cm->part_res[i]->part->nFace; j++) {
          PDM_printf(" "PDM_FMT_G_NUM" ", faceLNToGNTag[i][j]);
      }
    }
    PDM_printf("\n");
  }

  PDM_g_num_t *b_tIntersects = NULL;
  int *b_strideOne = NULL;
  int *part_stride = NULL;

  PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            1,
                            &part_stride,
                            (void **) faceLNToGNTag,                                   
                            &b_strideOne,
                            (void **) &b_tIntersects);

  //Calculation of the number of faces on the processor
  PDM_g_num_t nFaceProc = 0;

  int size_block = PDM_part_to_block_n_elt_block_get(ptb);
  for (int i = 0; i < size_block; i++) {
    //If the face has not been removed
     if(b_tIntersects[i] == 0) {
         nFaceProc++;
     }        
  }    

  //Global numbering of the faces
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&nFaceProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

  //Index to position the local vertices
  beg_NumAbs -= nFaceProc;

  idx_write = 0;

  //Loop over the partition numbers, i = partition number

  for (int i = 0; i < size_block; i++) {
    //If the vertex has not been removed
    if(b_tIntersects[i] == 0) {
        b_tIntersects[i] = beg_NumAbs + (idx_write++) + 1;

    }
    else {
        b_tIntersects[i] = -1;
    }        
  }    

  PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb); 

  PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                       (PDM_g_num_t **) faceLNToGNPart,
                                                       nFacePart,
                                                       cm->nPart,
                                                       cm->comm);

  PDM_g_num_t  **faceLNToGNFine = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));

  for (int i = 0; i < cm->nPart; i++) {        
    faceLNToGNFine[i] = (PDM_g_num_t *) malloc(cm->part_ini[i]->nFace * sizeof(PDM_g_num_t));        
  }

  int strideOne = 1;

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &strideOne, 
                          (void *) b_tIntersects,
                          &part_stride,
                          (void **) faceLNToGNFine);

  if(0 == 1) {
    PDM_printf("\nContenu de faceLNToGNFine\n");
    for (int i = 0; i < cm->nPart; i++) {        
      //Loop over the partition faces, j = face number
      for (int j = 0; j < cm->part_ini[i]->nFace; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", faceLNToGNFine[i][j]);
      }
    }
    PDM_printf("\n");
  }

  for (int i = 0; i < cm->nPart; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    int nFace = cm->part_res[i]->part->nFace;
    cmp->faceLNToGN = (PDM_g_num_t *) malloc(nFace * sizeof(PDM_g_num_t));
    for (int j = 0; j < nFace; j++) {            
      cmp->faceLNToGN[j] = (PDM_g_num_t) faceLNToGNFine[i][cm->part_res[i]->coarseFaceToFineFace[j] - 1];
    }
  }

  if(0 == 1) {
    PDM_printf("\nContenu de faceLNToGN de la structure\n");
    for (int i = 0; i < cm->nPart; i++) {        
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cm->part_res[i]->part->nFace; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cm->part_res[i]->part->faceLNToGN[j]);
      }
    }
    PDM_printf("\n");
  }

  free(faceLNToGNPart);
  free(nFacePart);

  for (int i = 0; i < cm->nPart; i++) {
    free(faceLNToGNTag[i]);
  }
  free(faceLNToGNTag);

  for (int i = 0; i < cm->nPart; i++) {
    free(faceLNToGNFine[i]);
  }
  free(faceLNToGNFine);

  free (b_strideOne);
  free (part_stride);
  free (b_tIntersects);

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);    
}

/**
 *
 * \brief Updates the vtxLNToGN array for the coarse mesh are save it in the coarse mesh partition
 * 
 * \param [in]   cm                 Coarse mesh 
 *
 */

static void
_build_vtxLNToGN
(
_coarse_mesh_t * cm
)
{    
  PDM_g_num_t **vtxLNToGNPart = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));
  int *nVtxPart = (int *) malloc(cm->nPart * sizeof(int));

  for (int i = 0; i < cm->nPart; i++) {
    vtxLNToGNPart[i] = cm->part_ini[i]->vtxLNToGN;
    nVtxPart[i] = cm->part_ini[i]->nVtx;
  }

  if(0 == 1) {
    PDM_printf("Contenu de vtxLNToGNPart\n");
    for (int i = 0; i < cm->nPart; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", *(vtxLNToGNPart[i]));
    }
    PDM_printf("\n");

    PDM_printf("Contenu de nVtxPart\n");
    for (int i = 0; i < cm->nPart; i++) {
      PDM_printf(" %d ", nVtxPart[i]);
    }
    PDM_printf("\n");
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCk_POST_CLEANUP,
                                                     1.,
                                                     (PDM_g_num_t **) vtxLNToGNPart,
                                                     nVtxPart,
                                                     cm->nPart,
                                                     cm->comm);    

  PDM_g_num_t **vtxLNToGNTag = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));

  int idx_write = 0;

  for (int i = 0; i < cm->nPart; i++) {
    int nFineVtx = cm->part_ini[i]->nVtx;
    int nCoarseVtx = cm->part_res[i]->part->nVtx;
    idx_write = 0;
    vtxLNToGNTag[i] = (PDM_g_num_t *) malloc(nFineVtx * sizeof(PDM_g_num_t));
    //Loop over coarseFaceToFineFace, i = index of coarseVtxToFineVtx (from 0 to cm->part_res[iPart]->part->nVtx)
    for (int j = 0; j < nFineVtx; j++) {
      vtxLNToGNTag[i][j] = -1;
    }

    for (int j = 0; j < nCoarseVtx; j++) {
        //If the vertex studied is the same as in coarseVtxToFineVtx, it is to be stored
      int k = cm->part_res[i]->coarseVtxToFineVtx[j] - 1;
      vtxLNToGNTag[i][k] = 0;
    }
  }

  if (0 == 1) {
    PDM_printf("Contenu de vtxLNToGNTag\n");    
    for (int i = 0; i < cm->nPart; i++) {
      for (int j = 0; j < cm->part_res[i]->part->nVtx; j++) {
        PDM_printf(" "PDM_FMT_G_NUM, vtxLNToGNTag[i][j]);
      }
    }
    PDM_printf("\n");
  }

  PDM_g_num_t *b_tIntersects = NULL;
  int *b_strideOne = NULL;
  int *part_stride = NULL;

  PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            1,
                            &part_stride,
                            (void **) vtxLNToGNTag,                                               
                            &b_strideOne,
                            (void **) &b_tIntersects);

  //Calculation of the number of vertices on the processor
  PDM_g_num_t nVtxProc = 0;

  int size_block = PDM_part_to_block_n_elt_block_get(ptb);
  for (int i = 0; i < size_block; i++) {
    //If the vertex has not been removed
    if(b_tIntersects[i] == 0) {
      nVtxProc++;
    }
  }

  //Global numbering of the vertices
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&nVtxProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

  //Index to position the local vertices
  beg_NumAbs -= nVtxProc;

  idx_write = 0;

  //Loop over the partition numbers, i = partition number    
  for (int i = 0; i < size_block; i++) {
    //If the vertex has not been removed
    if (b_tIntersects[i] == 0) {
      b_tIntersects[i] = beg_NumAbs + (idx_write++) + 1;
    }
    else {
      b_tIntersects[i] = -1;
    }           
  }

  PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb); 

  PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                       (PDM_g_num_t **) vtxLNToGNPart,
                                                       nVtxPart,
                                                       cm->nPart,
                                                       cm->comm);

  PDM_g_num_t  **vtxLNToGNFine = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));

  for (int i = 0; i < cm->nPart; i++) {        
    vtxLNToGNFine[i] = (PDM_g_num_t *) malloc(cm->part_ini[i]->nVtx * sizeof(PDM_g_num_t));        
  }

  int strideOne = 1;

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &strideOne, 
                          (void *) b_tIntersects,
                          &part_stride,
                          (void **) vtxLNToGNFine);

  if(0 == 1) {
    PDM_printf("\nContenu de vtxLNToGNFine\n");
    for (int i = 0; i < cm->nPart; i++) {        
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cm->part_ini[i]->nVtx; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", vtxLNToGNFine[i][j]);
      }
    }
    PDM_printf("\n");
  }

  for (int i = 0; i < cm->nPart; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    int nVtx = cm->part_res[i]->part->nVtx;
    cmp->vtxLNToGN = (PDM_g_num_t *) malloc(nVtx * sizeof(PDM_g_num_t));
    for (int j = 0; j < nVtx; j++) {            
      cmp->vtxLNToGN[j] = (PDM_g_num_t) vtxLNToGNFine[i][cm->part_res[i]->coarseVtxToFineVtx[j] - 1];
    }
  }

  if(0 == 1) {
    PDM_printf("\nContenu de vtxLNToGN de la structure\n");
    for (int i = 0; i < cm->nPart; i++) {        
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cm->part_res[i]->part->nVtx; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cm->part_res[i]->part->vtxLNToGN[j]);
      }
    }
    PDM_printf("\n");
  }

  free(vtxLNToGNPart);
  free(nVtxPart);

  for (int i = 0; i < cm->nPart; i++) {
    free(vtxLNToGNTag[i]);
  }
  free(vtxLNToGNTag);

  for (int i = 0; i < cm->nPart; i++) {
    free(vtxLNToGNFine[i]);
  }
  free(vtxLNToGNFine);

  free (b_strideOne);
  free (part_stride);
  free (b_tIntersects);

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
}

/**
 *
 * \brief Updates the faceGroupLNToGN array for the coarse mesh are save it in the coarse mesh partition
 * 
 * \param [in]   cm                 Coarse mesh
 *
 */

static void 
_build_faceGroupLNToGN
(
_coarse_mesh_t * cm
)
{        
  for (int i = 0; i < cm->nPart; i++) {
    //Si un des faceGroupIdx est NULL, on n'utilise pas les groupes
    //On quitte donc la boucle et la fonction !
    if(cm->part_ini[i]->faceGroupIdx == NULL)  {
      return;
    }
  }

  int **faceLNToGNTag = (int **) malloc(cm->nPart * sizeof(int *));
    
  int idx_write = 0;
    
  for (int i = 0; i < cm->nPart; i++) {
    idx_write = 0;
    faceLNToGNTag[i] = (int *) malloc(cm->part_ini[i]->nFace * sizeof(int));
    int nFace = cm->part_res[i]->part->nFace;
    for (int j = 0; j < cm->part_ini[i]->nFace; j++) {
      faceLNToGNTag[i][j] = -1;
    }
    //Loop over coarseFaceToFineFace, i = index of coarseFaceToFineFace (from 0 to cm->part_res[iPart]->part->nFace)
    for (int j = 0; j < nFace; j++) {
      //If the face studied is the same as in coarseFaceToFineFace, it is to be stored
      int k =  cm->part_res[i]->coarseFaceToFineFace[j] - 1;
      faceLNToGNTag[i][k] = 0;
    }
  }

  if (0 == 1) {
    PDM_printf("Contenu de faceLNToGNTag\n");    
    for (int i = 0; i < cm->nPart; i++) {
      for (int j = 0; j < cm->part_ini[i]->nFace; j++) {
        PDM_printf(" %d ", faceLNToGNTag[i][j]);
      }
    }
    PDM_printf("\n");
  }
    
  int **faceGroupLNToGNTag = (int **) malloc(cm->nPart * sizeof(int *));
    
  idx_write = 0;
  
  int _rank;
  
  PDM_MPI_Comm_rank (cm->comm, &_rank);
  
  for (int i = 0; i < cm->nPart; i++) {
    _part_t *cmp_coarse = cm->part_res[i]->part;
    _part_t *cmp_fine = cm->part_ini[i];
    fflush(stdout);
    faceGroupLNToGNTag[i] = (int *) malloc(cmp_coarse->faceGroupIdx[cm->nFaceGroup] * sizeof(int));
        
    //Loop over faceGroupIdx, i = index of group of faces (from 0 to cm->part_res[iPart]->part->nFaceGroup)
    for (int j = 0; j < cmp_fine->faceGroupIdx[cm->nFaceGroup]; j++) {
      faceGroupLNToGNTag[i][j] = faceLNToGNTag[i][cmp_fine->faceGroup[j] - 1];
    }
  }

  if(0 == 1) {
    PDM_printf("Contenu de faceGroupLNToGNTag\n");    
    for (int i = 0; i < cm->nPart; i++) {
      for (int j = 0; j < cm->part_res[i]->part->faceGroupIdx[cm->nFaceGroup]; j++) {
        PDM_printf(" %d ", faceGroupLNToGNTag[i][j]);
      }
    }
    PDM_printf("\n");
  }
    
  for (int i = 0; i < cm->nPart; i++) {
    free(faceLNToGNTag[i]);
  }
  free(faceLNToGNTag);
    
  for (int i = 0; i < cm->nPart; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    cmp->faceGroupLNToGN = (PDM_g_num_t *) malloc(cmp->faceGroupIdx[cm->nFaceGroup] * sizeof(PDM_g_num_t));
  }
    
  for(int iGroup = 0; iGroup < cm->nFaceGroup; iGroup++) {
    
    PDM_g_num_t **faceGroupLNToGNPart = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));
    int *nFaceGroupPart = (int *) malloc(cm->nPart * sizeof(int));
    
    for (int i = 0; i < cm->nPart; i++) {
      _part_t *cmp = cm->part_ini[i];
      faceGroupLNToGNPart[i] = &(cmp->faceGroupLNToGN[cmp->faceGroupIdx[iGroup]]);
      nFaceGroupPart[i] = cmp->faceGroupIdx[iGroup + 1] - cmp->faceGroupIdx[iGroup];
    }

    if(0 == 1) {
      PDM_printf("Contenu de faceGroupLNToGNPart\n");
      for (int i = 0; i < cm->nPart; i++) {
        int nFaceCurrentGroup = cm->part_ini[i]->faceGroupIdx[iGroup+1] - cm->part_ini[i]->faceGroupIdx[iGroup];
        for (int j = 0; j < nFaceCurrentGroup; j++)
          PDM_printf(" "PDM_FMT_G_NUM" ", faceGroupLNToGNPart[i][j]);
        PDM_printf("\n");
      }

      PDM_printf("Contenu de nFaceGroupPart\n");
      for (int i = 0; i < cm->nPart; i++) {
        PDM_printf(" %d ", nFaceGroupPart[i]);
      }
      PDM_printf("\n");
    }
        
    int rank;
    PDM_MPI_Comm_rank(cm->comm, &rank);

    PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCk_POST_CLEANUP,
                                                         1.,
                                                         (PDM_g_num_t **) faceGroupLNToGNPart,
                                                         nFaceGroupPart,
                                                         cm->nPart,
                                                         cm->comm);    

    PDM_g_num_t *b_tIntersects = NULL;
    int *b_strideOne = NULL;
    int *part_stride = NULL;
    
    PDM_g_num_t **faceGroupLNToGNTagGroup = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));

    for (int i = 0; i < cm->nPart; i++) {
      _part_t *cmp = cm->part_res[i]->part;
      int nFacePerGroup = cmp->faceGroupIdx[iGroup + 1] - cmp->faceGroupIdx[iGroup];
      faceGroupLNToGNTagGroup[i] = (PDM_g_num_t *) malloc(nFacePerGroup * sizeof(PDM_g_num_t));
      
      idx_write = 0;
      //Copy of the sub-array faceGroupLNToGN for each group
      for (int j = cmp->faceGroupIdx[iGroup]; j < cmp->faceGroupIdx[iGroup + 1]; j++) {
        faceGroupLNToGNTagGroup[i][idx_write++] = faceGroupLNToGNTag[i][j];
      }
    }

    if (0 == 1) {
      PDM_printf("Contenu de faceGroupLNToGNTagGroup\n");    
      for (int i = 0; i < cm->nPart; i++) {
        int nFacePerGroup = cm->part_res[i]->part->faceGroupIdx[iGroup + 1] - cm->part_res[i]->part->faceGroupIdx[iGroup];
        for (int j = 0; j < nFacePerGroup; j++) {
          PDM_printf(" %d ", faceGroupLNToGNTag[i][j]);
        }
      }
      PDM_printf("\n");
    }
        
    PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            1,
                            &part_stride,
                            (void **) faceGroupLNToGNTagGroup,
                            &b_strideOne,
                            (void **) &b_tIntersects);

    //Calculation of the number of faces for all the groups on the processor
    PDM_g_num_t nFaceGroupProc = 0;
    
    int size_block = PDM_part_to_block_n_elt_block_get(ptb);

    for (int i = 0; i < size_block; i++) {
      //If the face of the group has not been removed
      if (b_tIntersects[i] == 0) {
        nFaceGroupProc++;
      }

    }
        
    //Global numbering of the vertices
    PDM_g_num_t beg_NumAbs;

    PDM_MPI_Scan(&nFaceGroupProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

    //Index to position the local vertices
    beg_NumAbs -= nFaceGroupProc;
    
    idx_write = 0;
    
    //Loop over the partition numbers, i = partition number

    for (int i = 0; i < size_block; i++) {
      //If the vertex has not been removed
      if(b_tIntersects[i] == 0) {
        b_tIntersects[i] = beg_NumAbs + (idx_write++) + 1;
        
      }
      else {
        b_tIntersects[i] = -1;
      }           

    }

    PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb);
        
    //        PDM_printf("assert : [%d] %d %d\n", rank, size_block, blockDistribIdx[rank+1] - blockDistribIdx[rank] );
    //        assert(blockDistribIdx[rank+1] - blockDistribIdx[rank] ==  size_block);
    
    PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                         (PDM_g_num_t **) faceGroupLNToGNPart,
                                                         nFaceGroupPart,
                                                         cm->nPart,
                                                         cm->comm);

    PDM_g_num_t  **faceGroupLNToGNFine = (PDM_g_num_t **) malloc(cm->nPart * sizeof(PDM_g_num_t *));

    for (int i = 0; i < cm->nPart; i++) {
      _part_t *cmp = cm->part_ini[i];
      int nFacePerGroup = cmp->faceGroupIdx[iGroup + 1] - cmp->faceGroupIdx[iGroup];
      faceGroupLNToGNFine[i] = (PDM_g_num_t *) malloc(nFacePerGroup * sizeof(PDM_g_num_t));        
    }

    int strideOne = 1;
    
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &strideOne, 
                            (void *) b_tIntersects,
                            &part_stride,
                            (void **) faceGroupLNToGNFine);

    if(0 == 1) {
      PDM_printf("\nContenu de faceGroupLNToGNFine\n");
      for (int i = 0; i < cm->nPart; i++) {
        _part_t *cmp = cm->part_ini[i];
        //Loop over the partition vertices, j = vertex number
        int nFacePerGroup = cmp->faceGroupIdx[iGroup + 1] - cmp->faceGroupIdx[iGroup];
        for (int j = 0; j < nFacePerGroup; j++) {
          PDM_printf(" "PDM_FMT_G_NUM" ", faceGroupLNToGNFine[i][j]);
        }
      }
      PDM_printf("\n");  
    }

    for (int i = 0; i < cm->nPart; i++) {
      _part_t *cmp = cm->part_res[i]->part;int nFacePerGroupCoarse = cmp->faceGroupIdx[iGroup + 1] - cmp->faceGroupIdx[iGroup];
      for (int j = 0; j < nFacePerGroupCoarse; j++) {
        int idxConcatenation = cmp->faceGroupIdx[iGroup];
        cmp->faceGroupLNToGN[idxConcatenation + j] = 
                faceGroupLNToGNFine[i][cm->part_res[i]->coarseFaceGroupToFineFaceGroup[idxConcatenation + j] - 1];
      }
    }
        
    free(faceGroupLNToGNPart);
    free(nFaceGroupPart);

    for (int i = 0; i < cm->nPart; i++) {
      free(faceGroupLNToGNTagGroup[i]);
    }
    free(faceGroupLNToGNTagGroup);        
    
    for (int i = 0; i < cm->nPart; i++) {
      free(faceGroupLNToGNFine[i]);
    }
    free(faceGroupLNToGNFine);
    
    free (b_strideOne);
    free (part_stride);
    free (b_tIntersects);
    
    PDM_part_to_block_free(ptb);
    PDM_block_to_part_free(btp);
    
  }
    
  for (int i = 0; i < cm->nPart; i++) {
    free(faceGroupLNToGNTag[i]);
  }
  free(faceGroupLNToGNTag);

  if (0 == 1) {
    for (int i = 0; i < cm->nPart; i++) {
      PDM_printf("\nContenu de faceGroupLNToGN de la structure %d\n", i);
      _part_t *cmp = cm->part_res[i]->part;
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cmp->faceGroupIdx[cm->nFaceGroup]; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cmp->faceGroupLNToGN[j]);
      }
    }
    PDM_printf("\n"); 
  }

}

/**
 *
 * \brief Updates the facePartBound array for the coarse mesh are save it in the coarse mesh partition
 *  
 * \param [in]   cm                 Coarse mesh 
 *
 */

static void 
_build_facePartBound
(
_coarse_mesh_t * cm
)
{
  PDM_printf("_build_facePartBound \n");
  //Number of processors
  int nProc;
  PDM_MPI_Comm_size(cm->comm, &nProc);
  
  //Loop over the partitions of part_ini and part_res
  for (int iPart = 0; iPart < cm->nPart; iPart++) {
    _part_t *cmp_fine = cm->part_ini[iPart];
    _part_t *cmp_coarse = cm->part_res[iPart]->part;
    
    //Copy of facePartBoundPartIdx for the coarse mesh
    int *coarseFacePartBoundPartIdx = (int *) malloc((cm->nTPart + 1) * sizeof(int));
    for (int i = 0; i < cm->nTPart + 1; i++) {
      coarseFacePartBoundPartIdx[i] = cmp_fine->facePartBoundPartIdx[i];
    }
    cmp_coarse->facePartBoundPartIdx = coarseFacePartBoundPartIdx;
    
    //Copy of facePartBoundProcIdx for the coarse mesh
    int *coarseFacePartBoundProcIdx = (int *) malloc((nProc + 1) * sizeof(int));
    for (int i = 0; i < nProc + 1; i++) {
      coarseFacePartBoundProcIdx[i] = cmp_fine->facePartBoundProcIdx[i];
    }
    cmp_coarse->facePartBoundProcIdx = coarseFacePartBoundProcIdx;
  }
  
  int **fineFaceToCoarseFace = (int **) malloc(cm->nPart * sizeof(int *));
  
  for (int iPart = 0; iPart < cm->nPart; iPart++) {
    _part_t *cmp_fine = cm->part_ini[iPart];
    _part_t *cmp_coarse = cm->part_res[iPart]->part;
    
    //Creation of fineFaceToCoarseFace
    fineFaceToCoarseFace[iPart] = (int *) malloc(cmp_fine->nFace * sizeof(int));
    
    //Initialization to -1
    for (int i = 0; i < cmp_fine->nFace; i++) {
      fineFaceToCoarseFace[iPart][i] = -1;
    }
    
    //Loop over coarseFaceToFineFace
    for (int i = 0; i < cmp_coarse->nFace; i++) {
      int fineFace = cm->part_res[iPart]->coarseFaceToFineFace[i] - 1;
      fineFaceToCoarseFace[iPart][fineFace] = i + 1;
    }
    
    if(0 == 1) {
      PDM_printf("Final content of fineFaceToCoarseFace[%d]: \n",iPart);
      PDM_printf("Valeur de cm->part_ini[iPart]->nFace : %d \n", cm->part_ini[iPart]->nFace);
      for (int i = 0; i < cm->part_ini[iPart]->nFace; i++) {
        PDM_printf(" %d ", fineFaceToCoarseFace[iPart][i]);
      }
      PDM_printf("\n");
      PDM_printf("------------------------------------------\n\n");
    }

  }
        
  int *sendIdx = malloc(sizeof(int) * nProc);
  int *sendN = malloc(sizeof(int) * nProc); 
  PDM_g_num_t *sendBuff = NULL;
  
  int *recvIdx = malloc(sizeof(int) * nProc); 
  int *recvN = malloc(sizeof(int) * nProc); 
  PDM_g_num_t *recvBuff = NULL;
  
  for (int i = 0; i < nProc; i++) {
    sendN[i] = 0;
  }
  
  int n_t_send = 0;
  for (int i = 0; i < cm->nPart; i++)  {
    _part_t *cmp = cm->part_ini[i];         
    n_t_send += cmp->nFacePartBound; 
    
    for (int j = 0; j < cmp->nFacePartBound; j++) {
      int iProc    = cmp->facePartBound[4*j + 1]; 
      /* int iPart    = cmp->facePartBound[4*j + 2]; */ 
      /* int iFacDist = cmp->facePartBound[4*j + 3]; */
      sendN[iProc] += 1;
    }  
  }

  int n_t_recv = 0;
  for (int i = 0; i < cm->nPart; i++)  {
    _part_t *cmp_coarse = cm->part_res[i]->part; 
    int nFacePartBound = cm->part_ini[i]->nFacePartBound;
    cmp_coarse->nFacePartBound = nFacePartBound;
    
    n_t_recv += cmp_coarse->nFacePartBound; 
  }    
    
  sendIdx[0] = 0;
  for (int i = 1; i < nProc; i++) {
    sendIdx[i] = sendIdx[i - 1] + sendN[i - 1];
    sendN[i - 1] = 0;
  }    

  sendN[nProc - 1] = 0;
  
  sendBuff = malloc(sizeof(PDM_g_num_t) * n_t_send * 3);
  recvBuff = malloc(sizeof(PDM_g_num_t) * n_t_recv * 3);
  
  int **iFaceLocToIPartBound = (int **) malloc(cm->nPart * sizeof(int *));
  
  for (int i = 0; i < cm->nPart; i++) {
    _part_t *cmp = cm->part_ini[i]; 
    //Creation of iFaceLocToIPartBound
    iFaceLocToIPartBound[i] = (int *) malloc(cmp->nFace * sizeof(int));
    
    //Initialization to -1
    for (int cpt = 0; cpt < cmp->nFace; cpt++) {
      iFaceLocToIPartBound[i][cpt] = -1;
    }
    
    if(0 == 1) {
      PDM_printf("Valeur de nFacePartBound : %d \n", cmp->nFacePartBound);
      
      PDM_printf("Contenu de facePartBound initial de la partition %d\n", i);
      for (int cpt = 0; cpt < 4* cmp->nFacePartBound; cpt++) {
        if (cpt % 4 == 0)
          PDM_printf("|");
        PDM_printf(" %d ", cmp->facePartBound[cpt]);    
      }
      PDM_printf("\n");
    }
         
    for (int j = 0; j < cmp->nFacePartBound; j++) {             
      int iFacLoc  = cmp->facePartBound[4*j    ];
      int iProc    = cmp->facePartBound[4*j + 1];
      int iPart    = cmp->facePartBound[4*j + 2];
      int iFacDist = cmp->facePartBound[4*j + 3];
      
      iFaceLocToIPartBound[i][iFacLoc - 1] = j;             
      
      int id = sendIdx[iProc] + sendN[iProc];
             
      ++sendN[iProc];
//             
      sendBuff[3*id    ] = iPart;
      sendBuff[3*id + 1] = iFacDist;
      sendBuff[3*id + 2] = fineFaceToCoarseFace[i][iFacLoc - 1]; 
             
    }
        
    if (0 == 1) {
      PDM_printf("Contenu de iPartBoundToIFacLoc \n");
      for (int j = 0; j < cmp->nFacePartBound; j++) { //Modif : olp->nLinkedFace = cmp->nFacePartBound
        PDM_printf(" %d ", cmp->facePartBound[4*j    ]);
      }
      PDM_printf("\n");
      
      PDM_printf("Final content of iFaceLocToIPartBound[i]: \n");
      for (int cpt = 0; cpt < cm->part_ini[i]->nFace; cpt++) {
        PDM_printf(" %d ", iFaceLocToIPartBound[i][cpt]);
      }
      PDM_printf("\n");
    }
    
  }
    
  for (int i = 0; i < nProc; i++) {
    sendN[i] *= 3;
    sendIdx[i] *= 3;
  }
  
  PDM_MPI_Alltoall (sendN, 1, PDM_MPI_INT, 
                recvN, 1, PDM_MPI_INT, 
                cm->comm);
  recvIdx[0] = 0;
    
  for (int i = 1; i < nProc; i++)  {
    recvIdx[i] = recvIdx[i - 1] + recvN[i - 1];
  }

  PDM_MPI_Alltoallv(sendBuff, sendN, sendIdx, PDM__PDM_MPI_G_NUM,
                recvBuff, recvN, recvIdx, PDM__PDM_MPI_G_NUM, cm->comm);
  
  //Loop over the partitions
  for (int iPart = 0; iPart < cm->nPart; iPart++) {
    _part_t *cmp = cm->part_ini[iPart];
    _part_t *cmp_coarse = cm->part_res[iPart]->part;
    //Memory allocation of coarseFacePartBound (linked to part_res after)
    int *coarseFacePartBound = (int  *) malloc(4 * cmp->nFacePartBound * sizeof(int));
        
    //Copy of the facePartBound of part_ini
    for (int i = 0; i < 4 * cmp->nFacePartBound; i++) {
      coarseFacePartBound[i] = cmp->facePartBound[i];
    }
    cmp_coarse->facePartBound = coarseFacePartBound;
  }    
        
  //Loop over recvBuff
  for (int i = 0; i < (recvIdx[nProc - 1] + recvN[nProc - 1]) / 3; i++) {
    int iPartLoc       = (int) recvBuff[3 * i    ];
    int iFacLocFine    = (int) recvBuff[3 * i + 1];
    int iFacDistCoarse = (int) recvBuff[3 * i + 2];
    
    int posFacePartBoundCoarse = iFaceLocToIPartBound[iPartLoc - 1][iFacLocFine - 1];
    
    //Update of the data about faces in the facePartBound of part_res
    _part_t *cmp = cm->part_res[iPartLoc - 1]->part;
    cmp->facePartBound[4 * posFacePartBoundCoarse    ] = fineFaceToCoarseFace[iPartLoc - 1][iFacLocFine - 1]; 
    cmp->facePartBound[4 * posFacePartBoundCoarse + 3] = iFacDistCoarse;
  }
    
  if(0 == 1) {
    for (int iPart = 0; iPart < cm->nPart; iPart++) {
      PDM_printf("\nContent of facePartBound of part_res[%d]\n", iPart);
      PDM_printf("Valeur de cm->part_res[%d]->part->nFace : %d\n",iPart, cm->part_res[iPart]->part->nFace);
      for (int i = 0; i < 4 * cm->part_res[iPart]->part->nFacePartBound; i++) {
        if (i % 4 == 0)
          PDM_printf("|");
        PDM_printf(" %d ", cm->part_res[iPart]->part->facePartBound[i]);
      }
      PDM_printf("\n");
    }
    PDM_printf("\n");
  }

  for (int iPart = 0; iPart < cm->nPart; iPart++) { //Modif : nPartB => cm->nPart
  
    free(fineFaceToCoarseFace[iPart]);
    free(iFaceLocToIPartBound[iPart]);
  }
    
  free(fineFaceToCoarseFace);
  free(iFaceLocToIPartBound);
  
  free(sendN);
  free(sendIdx);
  free(sendBuff);
  
  free(recvN);
  free(recvIdx);
  free(recvBuff);
}

/**
 *
 * \brief Displays all the arrays of a partition of type _part_t
 * 
 * \param [in]  nPart        Number of partitions to define on this process
 * \param [in]  nTPart       Total number of partitions
 * \param [in]  nFaceGroup   Number of boundaries
 * 
 */

static void 
_part_display
(
  _part_t *part,
  int nPart,
  int nTPart,
  int nFaceGroup
)
{
  if (part == NULL) {
    PDM_printf("Incorrect part to display\n");
    return;        
  }
  
  PDM_printf("Value of nVtx : %d \n", part->nVtx);
  PDM_printf("Value of nCell : %d \n", part->nCell);
  PDM_printf("Value of nFace : %d \n", part->nFace);
  PDM_printf("Value of nFacePartBound : %d \n", part->nFacePartBound);
    
  if (part->cellFaceIdx != NULL) {
    PDM_printf("\nContent of cellFaceIdx\n");    
    for(int i = 0; i < part->nCell + 1; i++) {
      PDM_printf(" %d ", part->cellFaceIdx[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->gCellFace != NULL) {
    PDM_printf("\nContent of gCellFace\n");    
    for(int i = 0; i < part->cellFaceIdx[part->nCell]; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->gCellFace[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->cellFace != NULL) {
    PDM_printf("\nContent of cellFace\n");    
    for(int i = 0; i < part->cellFaceIdx[part->nCell]; i++) {
      PDM_printf(" %d ", part->cellFace[i]);
      if (i % (part->cellFaceIdx)[1] == (part->cellFaceIdx)[1] - 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }
    
  if (part->cellLNToGN != NULL) {
    PDM_printf("\nContent of cellLNToGN\n");    
    for(int i = 0; i < part->nCell; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->cellLNToGN[i]);
    }
    PDM_printf("\n");
  }
     
  if (part->cellTag != NULL) {
    PDM_printf("\nContent of cellTag\n");    
    for(int i = 0; i < part->nCell; i++) {
      PDM_printf(" %d ", part->cellTag[i]);
    }
    PDM_printf("\n");
  }    
    
  if (part->faceCell != NULL) {
    PDM_printf("\nContent of faceCell\n");    
    for(int i = 0; i < 2 * part->nFace; i++) {
      PDM_printf(" %d ", part->faceCell[i]);
      if (i % 2 == 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }    
    
  if (part->faceVtxIdx != NULL) {        
    PDM_printf("\nContent of faceVtxIdx\n");    
    for(int i = 0; i < part->nFace + 1; i++) {
      PDM_printf(" %d ", part->faceVtxIdx[i]);
    }
    PDM_printf("\n");
  }    
    
  if (part->gFaceVtx != NULL) {
    PDM_printf("\nContent of gFaceVtx\n");    
    for(int i = 0; i < part->faceVtxIdx[part->nFace]; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->gFaceVtx[i]);
    }
    PDM_printf("\n");
  }    
    
  if (part->faceVtx != NULL) {        
    PDM_printf("\nContent of faceVtx\n");            
    for(int i = 0; i < part->faceVtxIdx[part->nFace]; i++) {
      PDM_printf(" %d ", part->faceVtx[i]);
      if (i % (part->faceVtxIdx)[1] == (part->faceVtxIdx)[1] - 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }    
  
  if (part->faceLNToGN != NULL) {
    PDM_printf("\nContent of faceLNToGN\n");    
    for(int i = 0; i < part->nFace; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->faceLNToGN[i]);
    }
    PDM_printf("\n");
  }    
    
  if (part->faceTag != NULL) {
    PDM_printf("\nContent of faceTag\n");    
    for(int i = 0; i < part->nFace; i++) {
      PDM_printf(" %d ", (part->faceTag)[i]);
    }
    PDM_printf("\n");        
  }
    
  if (part->facePartBoundPartIdx != NULL) {
    PDM_printf("\nContent of facePartBoundPartIdx\n");    
    for(int i = 0; i < nPart + 1; i++) {
      PDM_printf(" %d ", part->facePartBoundPartIdx[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->facePartBoundProcIdx != NULL) {
    PDM_printf("\nContent of facePartBoundProcIdx\n");    
    for(int i = 0; i < nTPart + 1; i++) {
      PDM_printf(" %d ", part->facePartBoundProcIdx[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->facePartBound != NULL) {
    PDM_printf("\nContent of facePartBound\n");    
    for(int i = 0; i < 4 * part->nFacePartBound; i++) {
      PDM_printf(" %d ", part->facePartBound[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->faceGroupIdx != NULL) {
    PDM_printf("\nContent of faceGroupIdx\n");    
    for(int i = 0; i < nFaceGroup + 1; i++) {
      PDM_printf(" %d ", part->faceGroupIdx[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->faceGroup != NULL) {       
    PDM_printf("\nContent of faceGroup\n");    
    for(int i = 0; i < part->faceGroupIdx[nFaceGroup]; i++) {
      PDM_printf(" %d ", part->faceGroup[i]);
    }
    PDM_printf("\n");
  }
    
  if (part->faceGroupLNToGN != NULL) {
    PDM_printf("\nContent of faceGroupLNToGN\n");    
    for(int i = 0; i < part->faceGroupIdx[nFaceGroup]; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->faceGroupLNToGN[i]);
    }
    PDM_printf("\n");
  }
  
  if (part->vtx != NULL) {
    PDM_printf("\nContent of vtx\n");    
    for(int i = 0; i < 3 * part->nVtx; i++) {
      PDM_printf(" %.1f ", part->vtx[i]);
      if (i % 3 == 2) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }
  
  if (part->vtxLNToGN != NULL) {
    PDM_printf("\nContent of vtxLNToGN\n");    
    for(int i = 0; i < part->nVtx; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->vtxLNToGN[i]);
    }
    PDM_printf("\n");
  }
    
    
  if (part->vtxTag != NULL) {
    PDM_printf("\nContent of vtxTag\n");    
    for(int i = 0; i < part->nVtx; i++) {
      PDM_printf(" %d ", part->vtxTag[i]);
    }
    PDM_printf("\n");
  }
    
}

/**
 *
 * \brief Displays all the arrays of a coarse partition of type _coarse_part_t
 * 
 * \param [in]  nPart        Number of partitions to define on this process
 * \param [in]  nTPart       Total number of partitions
 * \param [in]  nFaceGroup   Number of boundaries
 * 
 */

static void 
_coarse_part_display
(
  _coarse_part_t *coarse_part,
  int nPart,
  int nTPart,
  int nFaceGroup
)
{
  if (coarse_part == NULL) {
    PDM_printf("Incorrect coarse part to display\n");
    return;        
  }
    
  _part_display(coarse_part->part, nPart, nTPart, nFaceGroup);
        
  if (coarse_part->coarseCellCellIdx != NULL) {
    PDM_printf("\nContent of coarseCellCellIdx\n");    
    for(int i = 0; i < coarse_part->part->nCell + 1; i++) {
      PDM_printf(" %d ", coarse_part->coarseCellCellIdx[i]);
    }
    PDM_printf("\n");
  }
    
  if (coarse_part->coarseCellCell != NULL) {
    PDM_printf("\nContent of coarseCellCell\n");    
    for(int i = 0; i < coarse_part->coarseCellCellIdx[coarse_part->part->nCell]; i++) {
      PDM_printf(" %d ", coarse_part->coarseCellCell[i]);
    }
    PDM_printf("\n");
  }
    
  if (coarse_part->coarseFaceGroupToFineFaceGroup != NULL) {
    PDM_printf("\nContent of coarseFaceGroupToFineFaceGroup\n");    
    for(int i = 0; i < coarse_part->part->faceGroupIdx[nFaceGroup]; i++) {
      PDM_printf(" %d ", coarse_part->coarseFaceGroupToFineFaceGroup[i]);
    }
    PDM_printf("\n");
  }
    
  if (coarse_part->coarseFaceToFineFace != NULL) {
    PDM_printf("\nContent of coarseFaceToFineFace\n");    
    for(int i = 0; i < coarse_part->part->nFace; i++) {
      PDM_printf(" %d ", coarse_part->coarseFaceToFineFace[i]);
    }
    PDM_printf("\n");
  }
    
  if (coarse_part->coarseVtxToFineVtx != NULL) {
    PDM_printf("\nContent of coarseVtxToFineVtx\n");    
    for(int i = 0; i < coarse_part->part->nVtx; i++) {
      PDM_printf(" %d ", coarse_part->coarseVtxToFineVtx[i]);
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Free partition
 *
 * \param [in]   part      partition
 *
 */

static void 
_part_free
(
 _part_t *part
)
{
  if (part == NULL) {
    PDM_printf("Incorrect part to free");
    return;        
  }
    
  if (part->cellFaceIdx != NULL)
    free(part->cellFaceIdx);
  part->cellFaceIdx = NULL;    
  
  if (part->gCellFace != NULL)
    free(part->gCellFace);
  part->gCellFace = NULL;    
  
  if (part->cellFace != NULL)
    free(part->cellFace);
  part->cellFace = NULL;    
  
  if (part->cellLNToGN != NULL)
    free(part->cellLNToGN);
  part->cellLNToGN = NULL;    
  
  if (part->cellTag != NULL)
    free(part->cellTag);
  part->cellTag = NULL;    
  
  if (part->faceCell != NULL)
    free(part->faceCell);
  part->faceCell = NULL;    

  if (part->faceVtxIdx != NULL)    
    free(part->faceVtxIdx);
  part->faceVtxIdx = NULL;    

  if (part->gFaceVtx != NULL)
    free(part->gFaceVtx);
  part->gFaceVtx = NULL;    

  if (part->faceVtx != NULL)
    free(part->faceVtx);
  part->faceVtx = NULL;    

  if (part->faceLNToGN != NULL)
    free(part->faceLNToGN);
  part->faceLNToGN = NULL;    

  if (part->faceTag != NULL)      
    free(part->faceTag);
  part->faceTag = NULL;    

  if (part->facePartBoundProcIdx != NULL)
    free(part->facePartBoundProcIdx);
  part->facePartBoundProcIdx = NULL;    

  if (part->facePartBoundPartIdx != NULL)
    free(part->facePartBoundPartIdx);
  part->facePartBoundPartIdx = NULL;    

  if (part->facePartBound != NULL)
    free(part->facePartBound);
  part->facePartBound = NULL;    

  if (part->faceGroupIdx != NULL)
    free(part->faceGroupIdx);
  part->faceGroupIdx = NULL;    

  if (part->faceGroup != NULL)
    free(part->faceGroup);
  part->faceGroup = NULL;    

  if (part->faceGroupLNToGN != NULL)
    free(part->faceGroupLNToGN);
  part->faceGroupLNToGN = NULL;    
          
  if (part->vtx != NULL)
    free(part->vtx);
  part->vtx = NULL;    

  if (part->vtxLNToGN != NULL)
    free(part->vtxLNToGN);
  part->vtxLNToGN = NULL;    

  if (part->vtxTag != NULL)
    free(part->vtxTag);
  part->vtxTag = NULL;    

  free(part);
}

/**
 *
 * \brief Free coarse partition
 *
 * \param [in]   coarse_part     coarse partition
 *
 */

static void 
_coarse_part_free
(
 _coarse_part_t *coarse_part
)
{
  _part_free(coarse_part->part);

  if (coarse_part->specific_data != NULL) {
    free (coarse_part->specific_data);
  }
  
  if (coarse_part->coarseCellCell != NULL)
    free(coarse_part->coarseCellCell);
  coarse_part->coarseCellCell = NULL;    
  
  if (coarse_part->coarseCellCellIdx != NULL)
    free(coarse_part->coarseCellCellIdx);
  coarse_part->coarseCellCellIdx = NULL;  
  
  if (coarse_part->coarseFaceGroupToFineFaceGroup != NULL)
    free(coarse_part->coarseFaceGroupToFineFaceGroup);
  coarse_part->coarseFaceGroupToFineFaceGroup = NULL;   
  
  if (coarse_part->coarseFaceToFineFace != NULL)
    free(coarse_part->coarseFaceToFineFace);
  coarse_part->coarseFaceToFineFace = NULL;   
  
  if (coarse_part->coarseVtxToFineVtx != NULL)
    free(coarse_part->coarseVtxToFineVtx);
  coarse_part->coarseVtxToFineVtx = NULL;  
    
  free(coarse_part);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Return an initialized coarse mesh object
 *
 * \param [out]  cmId              Coarse mesh identifier
 * 
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   nPart             Number of partitions
 * \param [in]   nTPart            Total number of partitions
 * \param [in]   nFaceGroup        Total number of groups
 * \param [in]   have_cellTag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_faceTag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtxTag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cellWeight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_faceWeight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_faceGroup    Presence des tableaux de groupes de faces
 */

void 
PDM_part_coarse_mesh_create
(
 int                *cmId,
 PDM_MPI_Comm        comm, 
 const char*         method,
 const int           nPart,
 const int           nTPart,
 const int           nFaceGroup,
 const int           have_cellTag,
 const int           have_faceTag,
 const int           have_vtxTag,
 const int           have_cellWeight,
 const int           have_faceWeight,
 const int           have_faceGroup
)
{
  if (_cm == NULL) {
    _cm = PDM_Handles_create (4);
  }

  _coarse_mesh_t *cm  = _coarse_mesh_create (comm,
                                             method,
                                             nPart,
                                             nTPart,
                                             nFaceGroup,
                                             have_cellTag,
                                             have_faceTag,
                                             have_vtxTag,
                                             have_cellWeight,
                                             have_faceWeight,
                                             have_faceGroup);

  *cmId = PDM_Handles_store (_cm, cm);
}

void
PROCF (pdm_part_coarse_mesh_create_cf, PDM_PART_COARSE_MESH_CREATE_CF)
(
 int                *cmId,
 PDM_MPI_Fint       *fcomm,        
 const char         *method,
 const int          *l_method,
 const int          *nPart, 
 const int          *nTPart, 
 const int          *nFaceGroup,
 const int          *have_cellTag,
 const int          *have_faceTag,
 const int          *have_vtxTag,
 const int          *have_cellWeight,
 const int          *have_faceWeight,
 const int          *have_faceGroup
)
{
  
  PDM_MPI_Comm comm = PDM_MPI_Comm_f2c (*fcomm);

  char *_method = PDM_fortran_to_c_string (method, *l_method); 
  
  PDM_part_coarse_mesh_create (cmId,
                               comm,
                               _method,
                               *nPart,
                               *nTPart,
                               *nFaceGroup,
                               *have_cellTag,
                               *have_faceTag,
                               *have_vtxTag,
                               *have_cellWeight,
                               *have_faceWeight,
                               *have_faceGroup);

  free (_method);
  
}

/**
 *
 * \brief Build a coarse mesh
 *
 * \param [in]  cmId               Coarse mesh identifier 
 * \param [in]  iPart              Partition identifier
 * \param [in]  nCoarseCell        Number of cells in the coarse grid
 * \param [in]  nCell              Number of cells
 * \param [in]  nFace              Number of faces
 * \param [in]  nFacePartBound     Number of partitioning boundary faces
 * \param [in]  nVtx               Number of vertices 
 * \param [in]  nFaceGroup         Number of face groups             
 * \param [in]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1, numbering : 0 to n-1)
 * \param [in]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace
 *                                                             numbering : 1 to n)
 * \param [in]  cellTag            Cell tag (size = nCell)
 * \param [in]  cellLNToGN         Cell local numbering to global numbering (size = nCell, numbering : 1 to n)
 * \param [in]  cellWeight         Cell weight (size = nCell)
 * \param [in]  faceWeight         Face weight (size = nFace)
 * \param [in]  faceCell           Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [in]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1, numbering : 0 to n-1)
 * \param [in]  faceVtx            Face to Vertex connectivity (size = faceVertexIdx[nFace], numbering : 1 to n)
 * \param [in]  faceTag            Face tag (size = nFace)
 * \param [in]  faceLNToGN         Face local numbering to global numbering (size = nFace, numbering : 1 to n)
 * \param [in]  vtxCoord           Vertex coordinates (size = 3 * nVertex)
 * \param [in]  vtxTag             Vertex tag (size = nVertex)
 * \param [in]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx, numbering : 1 to n)
 * \param [in]  faceGroupIdx       Face group index (size = nFaceGroup + 1, numbering : 1 to n-1)
 * \param [in]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 */

void 
PDM_part_coarse_mesh_input
(
 int                 cmId,
 int                 iPart,
 const int           nCoarseCellWanted,
 const int           nCell,
 const int           nFace,
 const int           nVtx,
 const int           nFaceGroup, //FIXME: Argument a eliminer : Information deja donnee
 const int           nFacePartBound,
 const int          *cellFaceIdx,
 const int          *cellFace,
 const int          *cellTag,
 const int          *cellWeight,
 const int          *faceWeight,
 const PDM_g_num_t  *cellLNToGN,       
 const int          *faceCell,
 const int          *faceVtxIdx,
 const int          *faceVtx,
 const int          *faceTag,       
 const PDM_g_num_t  *faceLNToGN,       
 const double       *vtxCoord,
 const int          *vtxTag,
 const PDM_g_num_t  *vtxLNToGN,       
 const int          *faceGroupIdx,
 const int          *faceGroup,
 const PDM_g_num_t  *faceGroupLNToGN,
 const int          *facePartBoundProcIdx,       
 const int          *facePartBoundPartIdx,
 const int          *facePartBound               
)
{   
  _coarse_mesh_t * cm = _get_from_id (cmId);  
    
  _coarse_grid_mesh_input (cm,
                           iPart,
                           nCoarseCellWanted,
                           nCell,
                           nFace,
                           nVtx,
                           nFaceGroup,
                           nFacePartBound,
                           cellFaceIdx,
                           cellFace,
                           cellTag,
                           cellWeight,
                           faceWeight,
                           cellLNToGN,
                           faceCell,
                           faceVtxIdx,
                           faceVtx,
                           faceTag,
                           faceLNToGN,
                           vtxCoord,
                           vtxTag,
                           vtxLNToGN,
                           faceGroupIdx,
                           faceGroup,
                           faceGroupLNToGN,
                           facePartBoundProcIdx,
                           facePartBoundPartIdx,
                           facePartBound);   
}

void 
PROCF (pdm_part_coarse_mesh_input, PDM_PART_COARSE_MESH_INPUT)
(
 int                *cmId,
 int                *iPart,
 const int          *nCoarseCellWanted,
 const int          *nCell,
 const int          *nFace,
 const int          *nVtx,
 const int          *nFaceGroup,
 const int          *nFacePartBound,
 const int          *cellFaceIdx,
 const int          *cellFace,
 const int          *have_cellTag,
 const int          *cellTag,
 const int          *have_cellWeight,
 const int          *cellWeight,
 const int          *have_faceWeight,
 const int          *faceWeight,
 const PDM_g_num_t *cellLNToGN,        
 const int          *faceCell,
 const int          *faceVtxIdx,
 const int          *faceVtx,
 const int          *have_faceTag,
 const int          *faceTag,
 const PDM_g_num_t *faceLNToGN,        
 const double       *vtxCoord,
 const int          *have_vtxTag,
 const int          *vtxTag,
 const PDM_g_num_t *vtxLNToGN,
 const int          *have_faceGroup,
 const int          *faceGroupIdx,
 const int          *faceGroup,
 const PDM_g_num_t *faceGroupLNToGN,
 const int          *facePartBoundProcIdx,
 const int          *facePartBoundPartIdx,
 const int          *facePartBound        
)
{

  int *_cellTag = (int *) cellTag;
  int *_faceTag = (int *) faceTag;
  int *_vtxTag = (int *) vtxTag;
  int *_cellWeight = (int *) cellWeight;
  int *_faceWeight = (int *) faceWeight;
  int *_faceGroupIdx = (int *) faceGroupIdx;
  int *_faceGroup = (int *) faceGroup;
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
  PDM_g_num_t *_faceGroupLNToGN = (PDM_g_num_t *) faceGroupLNToGN;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
  
  if (*have_cellTag == 0) {
    _cellTag = NULL;
  }

  if (*have_faceTag == 0) {
    _faceTag = NULL;
  }

  if (*have_vtxTag == 0) {
    _vtxTag = NULL;
  }

  if (*have_cellWeight == 0) {
    _cellWeight = NULL;
  }

  if (*have_faceWeight == 0) {
    _faceWeight = NULL;
  }

  if (*have_faceGroup == 0) {
    _faceGroupIdx = NULL;
    _faceGroup = NULL;
    _faceGroupLNToGN = NULL;
  }

  PDM_part_coarse_mesh_input(*cmId,
                             *iPart,
                             *nCoarseCellWanted,
                             *nCell,
                             *nFace,
                             *nVtx,
                             *nFaceGroup,
                             *nFacePartBound,
                             cellFaceIdx,
                             cellFace,
                             _cellTag,
                             _cellWeight,
                             _faceWeight,
                             cellLNToGN,
                             faceCell,
                             faceVtxIdx,
                             faceVtx,
                             _faceTag,  
                             faceLNToGN,
                             vtxCoord,
                             _vtxTag,
                             vtxLNToGN,
                             _faceGroupIdx,
                             _faceGroup,
                             _faceGroupLNToGN,
                             facePartBoundProcIdx,
                             facePartBoundPartIdx,
                             facePartBound);      
}

/**
 *
 * \brief Updates all the arrays dealing with MPI exchanges
 *
 * \param [in] cmId               Coarse mesh identifier
 */

void 
PDM_part_coarse_mesh_compute
(
                                int cmId
)
{    
  _coarse_mesh_t * cm = _get_from_id (cmId);
  
  /* First step : Manage independently coarse grid generation */
  
  for (int iPart = 0; iPart < cm->nPart; iPart++) {
    _coarse_grid_compute(cm, iPart);
  }
  
  /* Second step : Manage MPI */
  int itime = 13;
  
  //    PDM_part_coarse_mesh_display(cmId); 
  PDM_timer_resume(cm->timer);
  
  _build_coarseCellLNToGN(cm);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  
  _build_faceLNToGN(cm);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  
  _build_vtxLNToGN(cm);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  
  _build_faceGroupLNToGN(cm);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;
  
  PDM_timer_resume(cm->timer);
  
  _build_facePartBound(cm);
  
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);  
  
  cm->times_elapsed[0]     = cm->times_elapsed[itime];
  cm->times_cpu[0]         = cm->times_cpu[itime];
  cm->times_cpu_u[0]       = cm->times_cpu_u[itime];
  cm->times_cpu_s[0]       = cm->times_cpu_s[itime];
  
  for (int i = itime; i > 1; i--) {
    cm->times_elapsed[i] -= cm->times_elapsed[i-1];
    cm->times_cpu[i]     -= cm->times_cpu[i-1];
    cm->times_cpu_u[i]   -= cm->times_cpu_u[i-1];
    cm->times_cpu_s[i]   -= cm->times_cpu_s[i-1];
  }  
}

void
PROCF (pdm_part_coarse_mesh_compute, PDM_PART_COARSE_MESH_COMPUTE)
(
 int *cmId
)
{
  PDM_part_coarse_mesh_compute(*cmId);
}

/**
 *
 * \brief Return a coarse mesh partition dimensions
 * 
 * \param [in]   ppartId            Coarse mesh identifier
 * \param [in]   iPart              Current partition
 * 
 * \param [out]  nCell              Number of cells
 * \param [out]  nFace              Number of faces
 * \param [out]  nFacePartBound     Number of partitioning boundary faces
 * \param [out]  nVtx               Number of vertices 
 * \param [out]  nProc              Number of processus
 * \param [out]  nTPart             Number of partitions
 * \param [out]  sCellFace          Size of cell-face connectivity 
 * \param [out]  sFaceVtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound     Size of facePartBound array
 * \param [out]  sFaceGroup         Size of faceGroup array 
 * \param [out]  sCoarseCellToFineCell  Size of coarseCellToFineCell array
 *
 */

void 
PDM_part_coarse_mesh_part_dim_get
(
 const int      cmId,
 int            iPart,
 int           *nCell,
 int           *nFace,
 int           *nFacePartBound,
 int           *nVtx,
 int           *nProc,
 int           *nTPart,
 int           *nFaceGroup,
 int           *sCellFace,
 int           *sFaceVtx,
 int           *sFaceGroup,
 int           *sCoarseCellToFineCell
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId); 
  
  _coarse_part_t *part_res = NULL;   
  
  int numProcs;
  PDM_MPI_Comm_size(cm->comm, &numProcs);

  if (iPart < cm->nPart) {        
    part_res = cm->part_res[iPart]; 
  }
  
  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_dim_get error : unknown partition\n");
    exit(1);
  }
  
  *nFaceGroup = cm->nFaceGroup;
  
  *nCell = part_res->part->nCell;
  *nFace = part_res->part->nFace;
  *nVtx  = part_res->part->nVtx;    

  *nFacePartBound  = part_res->part->nFacePartBound;
  *nProc           = numProcs;  
  *nTPart          = cm->nTPart;
  *sCellFace       = part_res->part->cellFaceIdx[*nCell];
  *sFaceVtx        = part_res->part->faceVtxIdx[*nFace];
  *sCoarseCellToFineCell = part_res->coarseCellCellIdx[*nCell];
  *sFaceGroup      = 0;
  if (cm->nFaceGroup > 0 && part_res->part->faceGroupIdx != NULL) {
    *sFaceGroup    = part_res->part->faceGroupIdx[cm->nFaceGroup];
  }
}

void
PROCF (pdm_part_coarse_mesh_part_dim_get, PDM_PART_COARSE_MESH_PART_DIM_GET)
(
 const int     *cmId,
 int           *iPart,
 int           *nCell,
 int           *nFace,
 int           *nFacePartBound,
 int           *nVtx,
 int           *nProc,
 int           *nTPart,
 int           *nFaceGroup,
 int           *sCellFace,
 int           *sFaceVtx,
 int           *sFaceGroup,
 int           *sCoarseCellToFineCell
)
{
  PDM_part_coarse_mesh_part_dim_get(*cmId,
                                    *iPart,
                                    nCell,
                                    nFace,
                                    nFacePartBound,
                                    nVtx,
                                    nProc,
                                    nTPart,
                                    nFaceGroup,
                                    sCellFace,
                                    sFaceVtx, 
                                    sFaceGroup,
                                    sCoarseCellToFineCell);
}

/**
 *
 * \brief Return a mesh partition
 * 
 * \param [in]   cmId               Coarse mesh identifier
 * \param [in]   iPart              Current partition
 * 
 * \param [out]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1, numbering : 0 to n-1)
 * \param [out]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace
 *                                                             numbering : 1 to n)
 * \param [out]  cellTag            Cell tag (size = nCell)
 * \param [out]  cellLNToGN         Cell local numbering to global numbering (size = nCell, numbering : 1 to n)
 * \param [out]  cellInitCellIdx    Array of indexes of the connected partitions (size : nCoarseCell + 1)
 * \param [out]  cellInitCell       Partitioning array (size : cellInitCellIdx[nCoarseCell]) 
 * 
 * \param [out]  faceCell           Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [out]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1, numbering : 0 to n-1)
 * \param [out]  faceVtx            Face to Vertex connectivity (size = faceVertexIdx[nFace], numbering : 1 to n)
 * \param [out]  faceTag            Face tag (size = nFace)
 * \param [out]  faceLNToGN         Face local numbering to global numbering (size = nFace, numbering : 1 to n)
 * \param [out]  faceInitFace       Coarse face - fine face connectivity (size = nCoarseFace)
 * 
 * \param [out]  vtxCoord           Vertex coordinates (size = 3 * nVtx)
 * \param [out]  vtxTag             Vertex tag (size = nVtx)
 * \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx, numbering : 1 to n)
 * \param [out]  vtxInitVtx         Coarse vertex - fine vertex connectivity (size = nCoarseVtx)
 * 
 * \param [out]  faceGroupIdx       Face group index (size = nFaceGroup + 1, numbering : 1 to n-1)
 * \param [out]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [out]  faceGroupLNToGN    Faces global numbering for each group 
 *                                  (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * 
 * \param [out]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
 * \param [out]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
 * \param [out]  facePartBound      Partitioning boundary faces (size = 4 * nFacePartBound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition 
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number 
 *                                          in the connected partition (numbering :1 to n)
 * 
 */

void 
PDM_part_coarse_mesh_part_get
(
 const int    cmId,
 const int    iPart,       
 int          **cellFaceIdx,
 int          **cellFace,
 int          **cellTag,
 PDM_g_num_t **cellLNToGN,
 int          **cellInitCellIdx,                  
 int          **cellInitCell,          
 int          **faceCell,
 int          **faceVtxIdx,
 int          **faceVtx,
 int          **faceTag,
 PDM_g_num_t **faceLNToGN,        
 int          **faceGroupInitFaceGroup,          
 int          **faceInitFace,          
 double       **vtxCoord,
 int          **vtxTag,
 PDM_g_num_t **vtxLNToGN,        
 int          **vtxInitVtx,          
 int          **faceGroupIdx,
 int          **faceGroup,
 PDM_g_num_t **faceGroupLNToGN,
 int          **facePartBoundProcIdx,
 int          **facePartBoundPartIdx,
 int          **facePartBound        
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId); 
  
  _coarse_part_t *part_res = NULL;   
  
  if (iPart < cm->nPart) {        
    part_res = cm->part_res[iPart]; 
  }
  
  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_get error : unknown partition\n");
    exit(1);
  }
  
  *cellInitCellIdx        = part_res->coarseCellCellIdx;
  *cellInitCell           = part_res->coarseCellCell;
  *faceGroupInitFaceGroup = part_res->coarseFaceGroupToFineFaceGroup;
  *faceInitFace           = part_res->coarseFaceToFineFace;
  *vtxInitVtx             = part_res->coarseVtxToFineVtx;
    
  *cellFaceIdx          = part_res->part->cellFaceIdx;
  *cellFace             = part_res->part->cellFace;
  *cellTag              = part_res->part->cellTag;
  *cellLNToGN           = part_res->part->cellLNToGN;
  *faceCell             = part_res->part->faceCell;
  *faceVtxIdx           = part_res->part->faceVtxIdx;
  *faceVtx              = part_res->part->faceVtx;
  *faceTag              = part_res->part->faceTag;
  *faceLNToGN           = part_res->part->faceLNToGN;
  *vtxCoord             = part_res->part->vtx;
  *vtxTag               = part_res->part->vtxTag;
  *vtxLNToGN            = part_res->part->vtxLNToGN;
  *faceGroupIdx         = part_res->part->faceGroupIdx;
  *faceGroup            = part_res->part->faceGroup;
  *faceGroupLNToGN      = part_res->part->faceGroupLNToGN;
  *facePartBoundProcIdx = part_res->part->facePartBoundProcIdx;
  *facePartBoundPartIdx = part_res->part->facePartBoundPartIdx;
  *facePartBound        = part_res->part->facePartBound;
  
}


void
PROCF (pdm_part_coarse_mesh_part_get, PDM_PART_COARSE_MESH_PART_GET)
(
 int          *cmId,
 int          *iPart,       
 int          *cellFaceIdx,
 int          *cellFace,
 int          *cellTag,
 PDM_g_num_t *cellLNToGN,
 int          *cellInitCellIdx,                  
 int          *cellInitCell,          
 int          *faceCell,
 int          *faceVtxIdx,
 int          *faceVtx,
 int          *faceTag,
 PDM_g_num_t *faceLNToGN, 
 int          *faceGroupInitFaceGroup,
 int          *faceInitFace,          
 double       *vtxCoord,
 int          *vtxTag,
 PDM_g_num_t *vtxLNToGN,        
 int          *vtxInitVtx,          
 int          *faceGroupIdx,
 int          *faceGroup,
 PDM_g_num_t *faceGroupLNToGN,
 int          *facePartBoundProcIdx,
 int          *facePartBoundPartIdx,
 int          *facePartBound
)
{
  _coarse_mesh_t * cm = _get_from_id (*cmId); 
  
  int numProcs;
  PDM_MPI_Comm_size(cm->comm, &numProcs);
  
  _coarse_part_t *part_res = NULL;   
    
  if (*iPart < cm->nPart) {        
    part_res = cm->part_res[*iPart]; 
  }
  
  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_get error : unknown partition\n");
    exit(1);
  }
    
  for (int i = 0; i < part_res->part->nCell + 1; i++)
    cellFaceIdx[i] = part_res->part->cellFaceIdx[i];

  for (int i = 0; i < part_res->part->cellFaceIdx[part_res->part->nCell]; i++)
    cellFace[i] = part_res->part->cellFace[i];
  
  for (int i = 0; i < part_res->part->nCell; i++){        
    cellLNToGN[i] = part_res->part->cellLNToGN[i];
    
    if (part_res->part->cellTag != NULL)          
        cellTag[i] = part_res->part->cellTag[i];    
  }
    
  for (int i = 0; i < part_res->part->nCell + 1; i++)
    cellInitCellIdx[i] = part_res->coarseCellCellIdx[i];
  
  for (int i = 0; i < part_res->coarseCellCellIdx[part_res->part->nCell]; i++)
    cellInitCell[i] = part_res->coarseCellCell[i];
  
  for (int i = 0; i < 2 * part_res->part->nFace; i++){
    faceCell[i] = part_res->part->faceCell[i];
  }
  
  for (int i = 0; i < part_res->part->nFace + 1; i++)
    faceVtxIdx[i] = part_res->part->faceVtxIdx[i];
  
  for (int i = 0; i < part_res->part->faceVtxIdx[part_res->part->nFace]; i++)
    faceVtx[i] = part_res->part->faceVtx[i];
  
  for (int i = 0; i < part_res->part->nFace; i++){        
    faceLNToGN[i] = part_res->part->faceLNToGN[i];
    
    if (part_res->part->faceTag != NULL) 
      faceTag[i] = part_res->part->faceTag[i];      
  }
  
  if (part_res->part->faceGroupIdx != NULL) {
    for (int i = 0; i < part_res->part->faceGroupIdx[cm->nFaceGroup]; i++) {
      faceGroupInitFaceGroup[i] = part_res->coarseFaceGroupToFineFaceGroup[i];
    }
  }
  
  for (int i = 0; i < part_res->part->nFace; i++)
    faceInitFace[i] = part_res->coarseFaceToFineFace[i];
  
  for (int i = 0; i < 3 * part_res->part->nVtx; i++){
    vtxCoord[i] = part_res->part->vtx[i];
  }
  
  for (int i = 0; i < part_res->part->nVtx; i++){        
    vtxLNToGN[i] = part_res->part->vtxLNToGN[i];
    
    if (part_res->part->vtxTag != NULL) 
      vtxTag[i] = part_res->part->vtxTag[i];
    
  }
  
  for (int i = 0; i < part_res->part->nVtx; i++)
    vtxInitVtx[i] = part_res->coarseVtxToFineVtx[i];
  
  if (part_res->part->faceGroupIdx != NULL) {
    for (int i = 0; i < cm->nFaceGroup + 1; i++) {
      faceGroupIdx[i] = part_res->part->faceGroupIdx[i];
    }
  
    for (int i = 0; i < part_res->part->faceGroupIdx[cm->nFaceGroup]; i++) {
      faceGroup[i]       = part_res->part->faceGroup[i];
      faceGroupLNToGN[i] = part_res->part->faceGroupLNToGN[i];
    }
  }
  
  for (int i = 0; i < 4 * part_res->part->nFacePartBound; i++)
    facePartBound[i] = part_res->part->facePartBound[i];
  
  for (int i = 0; i < numProcs + 1; i++)
    facePartBoundProcIdx[i] = part_res->part->facePartBoundProcIdx[i];
  
  for (int i = 0; i < cm->nTPart + 1; i++)
    facePartBoundPartIdx[i] = part_res->part->facePartBoundPartIdx[i];
  
}


/**
 *
 * \brief Free coarse mesh
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

void 
PDM_part_coarse_mesh_free
(
 const int    cmId
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);   
  
  for (int i = 0; i < cm->nPart; i++) { 
    free(cm->part_ini[i]);   
    _coarse_part_free(cm->part_res[i]);
    cm->part_ini[i] = NULL;
    cm->part_res[i] = NULL;
  }

  if (cm->specific_data != NULL) {
    free (cm->specific_data);
  }
  
  free(cm->part_ini);
  free(cm->part_res);

  cm->part_ini = NULL;
  cm->part_res = NULL;
  
  PDM_timer_free(cm->timer);
  cm->timer = NULL;

  free(cm);
  
  PDM_Handles_handle_free (_cm, cmId, PDM_FALSE);

  const int n_cm = PDM_Handles_n_get (_cm);
  
  if (n_cm == 0) {
    _cm = PDM_Handles_free (_cm);
  }

}

void
PROCF (pdm_part_coarse_mesh_free, PDM_PART_COARSE_MESH_FREE)
(
  int *cmId
)
{
  PDM_part_coarse_mesh_free(*cmId);
}

/**
 *
 * \brief Return times
 * 
 * \param [in]   cmId        coarse mesh identifier
 * \param [out]  elapsed     elapsed times (size = 18)
 * \param [out]  cpu         cpu times (size = 18)
 * \param [out]  cpu_user    user cpu times (size = 18)
 * \param [out]  cpu_sys     system cpu times (size = 18)
 *
 */

void PDM_part_coarse_mesh_time_get
(
 int       cmId,
 double  **elapsed,
 double  **cpu,
 double  **cpu_user,
 double  **cpu_sys
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);   

  *elapsed  = cm->times_elapsed;
  *cpu      = cm->times_cpu;
  *cpu_user = cm->times_cpu_u;
  *cpu_sys  = cm->times_cpu_s;
}

void 
PROCF (pdm_part_coarse_mesh_time_get, PDM_PART_COARSE_MESH_TIME_GET)
(
 int      *cmId,
 double   *elapsed,
 double   *cpu,
 double   *cpu_user,
 double   *cpu_sys
 )
{
  _coarse_mesh_t * cm = _get_from_id (*cmId); 

  for (int i = 0; i < 18; i++) {
    elapsed[i]  = cm->times_elapsed[i];
    cpu[i]      = cm->times_cpu[i];
    cpu_user[i] = cm->times_cpu_u[i];
    cpu_sys[i]  = cm->times_cpu_s[i];
  }
}


/**
 *
 * \brief Displays all the arrays of a coarse mesh
 * 
 * \param [in]   cmId        Coarse mesh identifier
 * 
 */

void 
PDM_part_coarse_mesh_display
(
 const int    cmId
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);
  
  //Display all the elements of the structure (not part of the other structures)
  
  PDM_printf("\n");
  PDM_printf("Value of nPart : %d \n", cm->nPart);
  PDM_printf("Value of comm : %d \n", cm->comm);
  PDM_printf("Value of method : %d \n", cm->method);
  PDM_printf("Value of nTPart : %d \n", cm->nTPart);
  PDM_printf("Value of nFaceGroup : %d \n", cm->nFaceGroup);    
  
  for (int i = 0; i < cm->nPart; i++) { 
    PDM_printf("\n=============================================\n");
    PDM_printf("Valeur de i : %d \n", i);
    PDM_printf("\n=============================================\n");
    
    PDM_printf("\n----------Affichage de part_ini-----------------\n");
    _part_display(cm->part_ini[i],cm->nPart,cm->nTPart,cm->nFaceGroup);
    
    PDM_printf("\n----------Affichage de part_res-----------------\n");
    _coarse_part_display(cm->part_res[i], cm->nPart,cm->nTPart,cm->nFaceGroup);
    PDM_printf("\n-----------------------------------------\n");
  }
}

void
PROCF (pdm_part_coarse_mesh_display, PDM_PART_COARSE_MESH_DISPLAY)
(
  int *cmId
)
{
  PDM_part_coarse_mesh_display (*cmId);
}


/**
 *
 * \brief Add a new coarse mesh method
 *
 * \param [in]      name          Mesh entity to renumber
 * \param [in]      fct           Function
 *
 */

int
PDM_coarse_mesh_method_add
(
 const char                 *name,     /*!< Name          */
 PDM_coarse_mesh_fct_t       fct       /*!< Function      */
)
{
  if (_coarse_mesh_methods == NULL) {
      PDM_coarse_mesh_method_load_local();
  }

  _coarse_mesh_method_t *method_ptr = malloc (sizeof(_coarse_mesh_method_t));

  int idx = PDM_Handles_store  (_coarse_mesh_methods, method_ptr);

  method_ptr->name = malloc (sizeof(char) * (strlen(name) + 1));
  strcpy (method_ptr->name, name);

  method_ptr->fct = fct;

  return idx;
}


/**
 *
 * \brief Get index of a coarse mesh method from it's name
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PROCF (pdm_coarse_mesh_method_idx_get_cf, PDM_COARSE_MESH_METHOD_IDX_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  char *_name = PDM_fortran_to_c_string (name, *l_name); 

  *idx = PDM_coarse_mesh_method_idx_get (_name);
  
  free (_name);

}

int
PDM_coarse_mesh_method_idx_get
(
const char *name
)
{
  if (_coarse_mesh_methods == NULL) {
    PDM_coarse_mesh_method_load_local();
  }
  int idx = -1;
  
  if (_coarse_mesh_methods != NULL) {
    int n_methods = PDM_Handles_n_get (_coarse_mesh_methods);
    const int *index =  PDM_Handles_idx_get (_coarse_mesh_methods);

    for (int i = 0; i < n_methods; i++) {
      _coarse_mesh_method_t *method_ptr = 
              (_coarse_mesh_method_t *) PDM_Handles_get (_coarse_mesh_methods, index[i]);
      if (!strcmp(method_ptr->name, name)) {
        idx = index[i];
        break;
      }      
    }
  }
  return idx;
}


/**
 *
 * \brief Get name of a coarse mesh method from it's index
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PROCF (pdm_coarse_mesh_method_name_get_cf, PDM_COARSE_MESH_METHOD_NAME_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  const char *_name = PDM_coarse_mesh_method_name_get (*idx);

  const int _l_name = strlen(_name);

  *l_name = PDM_MAX (_l_name, PDM_MAX_CHAR_LENGTH);

  strncpy (name, _name, *l_name);
}

char *
PDM_coarse_mesh_method_name_get
(
const int id
)
{
  if (_coarse_mesh_methods == NULL) {
    PDM_coarse_mesh_method_load_local();
  }

  int n_methods = PDM_Handles_n_get (_coarse_mesh_methods);

  if (id >= n_methods) {
    return NULL;
  }

  const int *index =  PDM_Handles_idx_get (_coarse_mesh_methods);

  _coarse_mesh_method_t *method_ptr = 
            (_coarse_mesh_method_t *) PDM_Handles_get (_coarse_mesh_methods, index[id]);

  return method_ptr->name;
}


/**
 *
 * \brief Get the number of coarse mesh method
 *
 * \return Number of methods
 *
 */

void
PROCF (pdm_coarse_mesh_method_n_get, PDM_COARSE_MESH_METHOD_N_GET)
(
 int  *n_method
 )
{
  *n_method = PDM_coarse_mesh_method_n_get ();
}

int
PDM_coarse_mesh_method_n_get
(
void
)
{
  if (_coarse_mesh_methods == NULL) {
    PDM_coarse_mesh_method_load_local();
  }
  
  return PDM_Handles_n_get (_coarse_mesh_methods);
  
}

/**
 *
 * \brief Purge coarse mesh methods catalog
 *
 */

void
PDM_coarse_mesh_method_purge
(
void
)
{
  if (_coarse_mesh_methods != NULL) {

    const int *index =  PDM_Handles_idx_get (_coarse_mesh_methods);
    int n_methods = PDM_Handles_n_get (_coarse_mesh_methods);
    
    while (n_methods > 0) {
      int idx = index[0];
      _coarse_mesh_method_t *method_ptr = 
              (_coarse_mesh_method_t *) PDM_Handles_get (_coarse_mesh_methods, idx);
      free (method_ptr->name);
      PDM_Handles_handle_free (_coarse_mesh_methods, idx, PDM_TRUE);
      n_methods = PDM_Handles_n_get (_coarse_mesh_methods);
    }

    _coarse_mesh_methods = PDM_Handles_free (_coarse_mesh_methods);
    
  }
}

/**
 *
 * \brief Load local coarse mesh methods
 *
 */

void
PDM_coarse_mesh_method_load_local
(
void
)
{
  if (_coarse_mesh_methods == NULL)  {
    
    const int n_default_methods = 2;
    _coarse_mesh_methods = PDM_Handles_create (n_default_methods);
    
    PDM_coarse_mesh_method_add ("PDM_COARSE_MESH_SCOTCH",
                             _coarse_from_scotch);
    PDM_coarse_mesh_method_add ("PDM_COARSE_MESH_METIS",
                             _coarse_from_metis);

  }

}



/**
 *
 * \brief Return coarse mesh object from its identifier
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

_coarse_mesh_t *
PDM_part_coarse_mesh_get_from_id
(
 int  cmId
 )
{
  return _get_from_id (cmId);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

