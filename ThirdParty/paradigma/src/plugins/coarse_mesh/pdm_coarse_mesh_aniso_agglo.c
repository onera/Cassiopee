
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

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_geom_elem.h"
#include "pdm_part_coarse_mesh.h"
#include "pdm_part_coarse_mesh_priv.h"
#include "pdm_part_renum.h"
#include "pdm_order.h"
#include "pdm_isotropic_agglomerator.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_coarse_mesh_aniso_agglo.h"

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

/*=============================================================================
 * Static global variables
 *============================================================================*/


/**
 * \struct _part_aniso_agglo_data_t
 * \brief Specific data for aniso_agglo method
 *
 */

typedef struct  {

  int *anisotropicOption;

} _aniso_agglo_data_t;

/**
 * \struct _part_aniso_agglo_data_t
 * \brief Specific data for aniso_agglo method
 *
 */

typedef struct  {

  /* Array specific to anisotropic agglomeration */
  int *agglomerationLines;
  int *agglomerationLinesIdx;
  int  agglomerationLinesIdx_size;
  int *isOnFineBnd;

  /* Array specific to anisotropic agglomeration if Initialise from a finer grid */
  int *agglomerationLinesInit;
  int *agglomerationLinesInitIdx;
  int  agglomerationLinesInitIdx_size;
  int *isOnFineBndInit;

} _part_aniso_agglo_data_t;

/*============================================================================
 * Private function definitions
 *============================================================================*/

#ifdef PDM_HAVE_ANISO_AGGLO

static
void
_coarse_aniso_agglo_renum
(
  _coarse_part_t *cp
)
{
  // PDM_printf("_coarse_aniso_agglo_renum : \n");

  _part_aniso_agglo_data_t *part_aniso_agglo_data = (_part_aniso_agglo_data_t *) cp->specific_data;

  // Je pense qu'il faut faire le Init aussi

  if(part_aniso_agglo_data->agglomerationLines != NULL){

    int nAgglo = part_aniso_agglo_data->agglomerationLinesIdx[part_aniso_agglo_data->agglomerationLinesIdx_size-1];

    int *oldToNewOrder = (int *) malloc (cp->part->nCell * sizeof(int));
    for(int i = 0; i < cp->part->nCell; i++) {
     oldToNewOrder[cp->part->newToOldOrderCell[i]] = i;
    }

    for(int i = 0; i < nAgglo; ++i) {
      part_aniso_agglo_data->agglomerationLines[i] = oldToNewOrder[part_aniso_agglo_data->agglomerationLines[i]];
    }

    free(oldToNewOrder);
    // PDM_printf("_coarse_aniso_agglo_renum flag 2: \n");

  }

  if(part_aniso_agglo_data->isOnFineBnd != NULL){
    /* isOnFineBnd -> 0,1,2,3 */
    PDM_order_array(cp->part->nCell,
                    sizeof(int),
                    cp->part->newToOldOrderCell,
                    part_aniso_agglo_data->isOnFineBnd);
  }

}

/**
 *
 * \brief Perform the coarse mesh from the SCOTCH graph method
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_coarse_from_aniso_agglo
(
_coarse_mesh_t *cm,
const int       iPart,
int            *nCoarseCellComputed,
int            *cellCellIdx,
int            *cellCell,
int            *cellPart
)
{
  _part_t        *part_ini = cm->part_ini[iPart];
  _coarse_part_t *part_res = cm->part_res[iPart];

  /* Setup function ptr for renumbering */
  cm->specific_func = &_coarse_aniso_agglo_renum;

  /* TODO : See with Eric */
  if(part_res->specific_data == NULL ){
    part_res->specific_data = malloc (sizeof(_part_aniso_agglo_data_t));
  }
  _part_aniso_agglo_data_t *part_aniso_agglo_data = (_part_aniso_agglo_data_t *) part_res->specific_data;


  if ( cm->specific_data == NULL){
    cm->specific_data = malloc (sizeof(_aniso_agglo_data_t));
  }
  _aniso_agglo_data_t *aniso_agglo_data = (_aniso_agglo_data_t *) cm->specific_data;

   /** Allocate  **/
  int    isOriented  = 0;
  double *volume     = (double *) malloc(  part_ini->nCell * sizeof(double));
  double *centercell = (double *) malloc(3*part_ini->nCell * sizeof(double));

  double *surfaceArea   = (double *) malloc(  part_ini->nFace * sizeof(double));
  double *surfaceVector = (double *) malloc(3*part_ini->nFace * sizeof(double));
  double *surfaceCenter = (double *) malloc(3*part_ini->nFace * sizeof(double));

  /** Build up volume of all cells of Fine partition **/
  PDM_geom_elem_polyhedra_properties (isOriented,
                                      part_ini->nCell,
                                      part_ini->nFace,
                                      part_ini->faceVtxIdx,
                                      part_ini->faceVtx,
                                      part_ini->cellFaceIdx,
                                      part_ini->cellFace,
                                      part_ini->nVtx,
                                      part_ini->vtx,
                                      volume,
                                      centercell,
                                      NULL,
                                      NULL);

  /** Build up Area of all faces of Fine partition **/
  PDM_geom_elem_polygon_properties (part_ini->nFace,
                                    part_ini->faceVtxIdx,
                                    part_ini->faceVtx,
                                    part_ini->vtx,
                                    surfaceVector,
                                    surfaceCenter,
                                    NULL,
                                    NULL);

  /* Deduce surface area */
  for (int iface = 0; iface < part_ini->nFace; iface++) {
    surfaceArea[iface] = sqrt(  surfaceVector[3*iface  ]*surfaceVector[3*iface  ]
                              + surfaceVector[3*iface+1]*surfaceVector[3*iface+1]
                              + surfaceVector[3*iface+2]*surfaceVector[3*iface+2] );
  }


  int *agglomerationLinesInit         = part_aniso_agglo_data->agglomerationLinesInit;
  int *agglomerationLinesInitIdx      = part_aniso_agglo_data->agglomerationLinesInitIdx;
  int  agglomerationLinesInitIdx_size = part_aniso_agglo_data->agglomerationLinesInitIdx_size;
  int *isOnFineBndInit                = part_aniso_agglo_data->isOnFineBndInit;
  int *anisotropicOption              = aniso_agglo_data->anisotropicOption;

  /* Dramatic verbose */
  if(0 == 1)
  {
    PDM_printf("Address of agglomerationLinesInit    : %p \n", agglomerationLinesInit   );
    PDM_printf("Address of agglomerationLinesInitIdx : %p \n", agglomerationLinesInitIdx);
    PDM_printf("Address of isOnFineBndInit           : %p \n", isOnFineBndInit          );
    PDM_printf("Address of anisotropicOption         : %p \n", anisotropicOption        );
  }

  if (anisotropicOption[0] == 1){
    // first agglomeration
    // some array need to be allocated:
    // part_res->isOnFineBndInit, part_aniso_agglo_data->agglomerationLinesInit and part_aniso_agglo_data->agglomerationLinesInitIdx are supposed to be NULL!
    // 1) Create isOnFineBnd array to describe for every cells if it is on boundary

    isOnFineBndInit = malloc(( part_ini->nCell  )* sizeof(int));

    for (int icell = 0; icell<part_ini->nCell; icell++){
      isOnFineBndInit[icell] = 0;
    }

    /* Compute isOnFineBndInit */
    for (int iFace = 0; iFace<part_ini->nFace; iFace++){
      if(part_ini->faceCell[2*iFace+1] == 0){
        if(isOnFineBndInit[part_ini->faceCell[2*iFace]- 1] < 3){
          isOnFineBndInit[part_ini->faceCell[2*iFace]-1] += 1;
        }
      }
    }

    // C'est pas la meme condition que en haut ???
    // if (anisotropicOption[0] == 1){
    if (anisotropicOption[1] == 1){
      // Anisotropic agglomeration
      // we need to allocated
      PDM_printf("\n Anisotropic agglomeration\n");
      agglomerationLinesInit         = malloc( ( part_ini->nCell ) * sizeof(int));
      agglomerationLinesInitIdx      = malloc( ( part_ini->nCell ) * sizeof(int));
      agglomerationLinesInitIdx_size = part_ini->nCell ;

      for (int icell = 0; icell < part_ini->nCell; icell++){
        agglomerationLinesInit[icell]    = 0;
        agglomerationLinesInitIdx[icell] = 0;
      }
    }
  } /* End anisotropicOption */

  // Computation of countOnBndCells to know the number of cells (coarse or fine) that are on the boundary of the domain
  int countOnBndCells = 0;
  for (int icell = 0; icell < part_ini->nCell; icell++){
    // PDM_printf(" isOnFineBndInit[%i] : %i \n ",icell, isOnFineBndInit[icell]);
    if (isOnFineBndInit[icell]>0){
      countOnBndCells = countOnBndCells+1;
    }
  }

  if(0 == 1){
   PDM_printf(" countOnBndCells : %i \n ",countOnBndCells);
  }


  /* ---------------------------------------- */
  /* 2) Create adjacency with current cell iff on boundary */

  int *cellCellTmpIdx  = malloc(( part_ini->nCell + 1                           )* sizeof(int));
  int *cellCellTmp     = malloc(( cellCellIdx[part_ini->nCell] + countOnBndCells)* sizeof(int));

  cellCellTmpIdx[0] = 0;
  for (int icell = 0; icell < part_ini->nCell; icell++){

    if(isOnFineBndInit[icell] >0){
      int beg = cellCellIdx[icell  ];
      int end = cellCellIdx[icell+1];

      int nbc = end-beg;

      /* Add One slot for current cell */
      cellCellTmpIdx[icell+1] = cellCellTmpIdx[icell]+nbc+1;

      int beg2 = cellCellTmpIdx[icell  ];
      int end2 = cellCellTmpIdx[icell+1];
      // int nbc2 = end2-beg2;

      cellCellTmp[beg2] = icell;
      int cpt = 0;
      for (int icell2 = beg2+1; icell2 < end2; icell2++){
        cellCellTmp[icell2] = cellCell[beg+cpt++];
      }
    }
    else{
      int beg = cellCellIdx[icell  ];
      int end = cellCellIdx[icell+1];

      int nbc = end-beg;

      /* Add One slot for current cell */
      cellCellTmpIdx[icell+1] = cellCellTmpIdx[icell]+nbc;

      int beg2 = cellCellTmpIdx[icell  ];
      int end2 = cellCellTmpIdx[icell+1];

      int cpt = 0;
      for (int icell2 = beg2; icell2 < end2; icell2++){
        cellCellTmp[icell2] = cellCell[beg+cpt++];
      }
    }
  }

  //
  if(0 == 1){
    PDM_printf("_split \n ");
    for (int i=0; i<8; i++){
      PDM_printf("_split anisotropicOption[%i] : % i\n ", i, anisotropicOption[i]);
    }
  }
  /* ---------------------------------------- */
  /* Parameters */
  int *sizes = malloc(( 12  )* sizeof(int));

  sizes[0 ] = part_ini->nCell;
  sizes[1 ] = cellCellTmpIdx[part_ini->nCell];
  sizes[2 ] = 0;
  sizes[3 ] = 0;
  sizes[4 ] = 0;
  sizes[5 ] = 0;
  sizes[6 ] = 0;
  sizes[7 ] =   part_ini->nCell;
  sizes[10] = 2*part_ini->nFace;
  sizes[11] =   part_ini->nFace;

  if ( (anisotropicOption[0]==1) || (agglomerationLinesInitIdx_size == 0) )
  {
    sizes[8 ] = 0;
    sizes[9 ] = 0;
  }
  else
  {
    sizes[8 ] = agglomerationLinesInitIdx_size;
    sizes[9 ] = agglomerationLinesInitIdx[agglomerationLinesInitIdx_size] ;
  }

  int *arrayOfFineAnisotropicCompliantCells = malloc( (part_ini->nCell) * sizeof(int));
  for (int icell = 0; icell < part_ini->nCell; icell++){
    arrayOfFineAnisotropicCompliantCells[icell] = icell;
  }

  for (int i = 0; i < part_ini->nCell; i++){
    cellPart[i] = -1;
  }

  /* ---------------------------------------- */
  if(0 == 1 ){
    PDM_printf(" Size of cellCellIdx : %i \n ", cellCellIdx[part_ini->nCell]);
    PDM_printf("Sizes Before: \n ");
    for (int i=8; i<10; i++){
      PDM_printf(" sizes[%i] : % i\n ", i, sizes[i]);
    }
  }

  /* ---------------------------------------- */
  agglomerateOneLevel_v_Paradigma(sizes,
                                  cellCellTmpIdx,
                                  cellCellTmp,
                                  volume,
                                  arrayOfFineAnisotropicCompliantCells,
                                  isOnFineBndInit,
                                  part_ini->faceCell,
                                  surfaceArea,
                                  anisotropicOption[0],      // isFirstAgglomeration_int
                                  anisotropicOption[1],      // isAnisotropic_int
                                  cellPart,
                                  agglomerationLinesInitIdx,
                                  agglomerationLinesInit,
                                  anisotropicOption[2],      // dimension
                                  anisotropicOption[3],      // goalCard
                                  anisotropicOption[4],      // minCard
                                  anisotropicOption[5],      // maxCard
                                  anisotropicOption[6],      // checks_int
                                  anisotropicOption[7] );    // verbose_int

  if(0 == 1 ){
    PDM_printf("sizes After \n ");
    for (int i=8; i<10; i++){
      PDM_printf(" sizes[%i] : % i\n ", i, sizes[i]);
    }
    PDM_printf("After agglomerateOneLevel_v_Paradigma \n ");
    PDM_printf("Nb Coarse cell %i", sizes[2]);
    PDM_printf(" cellPart");
    for (int i = 0; i < part_ini->nCell; i++){
      PDM_printf(" %i ", cellPart[i]);
    }
  }

  (*nCoarseCellComputed) = sizes[2];

  if(0 == 1){
    PDM_printf("NbCoarse  : %i \n ",(*nCoarseCellComputed));
  }

  // Update of  Aniso_Agglo data:
  // first: isOnCoarseBnd!
  // int   *isOnCoarseBnd = part_res->isOnFineBnd;  //useless
  int *isOnCoarseBnd = malloc(( (*nCoarseCellComputed)  )* sizeof(int));
  for (int iCC = 0; iCC<(*nCoarseCellComputed); iCC++){
    isOnCoarseBnd[iCC]=0;
  }

  int numberOfFineCell = sizes[0];
  // We process every cell of level iLevel -1
  for(int iFineCell=0; iFineCell<numberOfFineCell; iFineCell++){
    short bndFC = (short) isOnFineBndInit[iFineCell];
    int iCoarseCell =cellPart[iFineCell];
    if (bndFC > isOnCoarseBnd[iCoarseCell]){
      isOnCoarseBnd[iCoarseCell] = bndFC;
    }
  }
  part_aniso_agglo_data->isOnFineBnd = isOnCoarseBnd;

  // then: AgglomerationLines
  if( sizes[8] > 0){

    /* Allocation agglomeration lines */
    int *agglomerationLinesIdx = malloc( ( sizes[8] ) * sizeof(int));
    int *agglomerationLines    = malloc( ( sizes[9] ) * sizeof(int));

    for (int iCC =0; iCC<sizes[8]; iCC++){
      // printf("agglomerationLinesIdx[%i] = %i \n", iCC, agglomerationLinesInitIdx[iCC] );
      agglomerationLinesIdx[iCC] = agglomerationLinesInitIdx[iCC];
    }

    for (int iCC =0; iCC<sizes[9]; iCC++){
      // printf("agglomerationLines[%i] = %i \n", iCC, agglomerationLinesInit[iCC] );
      agglomerationLines[iCC] = agglomerationLinesInit[iCC];
    }

    part_aniso_agglo_data->agglomerationLinesIdx      = agglomerationLinesIdx;
    part_aniso_agglo_data->agglomerationLines         = agglomerationLines;
    part_aniso_agglo_data->agglomerationLinesIdx_size = sizes[8];

    if( 0 == 1 ){
      // PDM_printf(" agglomerationLinesInitIdx_size: %i \n", part_aniso_agglo_data->agglomerationLinesInitIdx_size);
      PDM_printf(" agglomerationLinesIdx_size    : %i \n", part_aniso_agglo_data->agglomerationLinesIdx_size);
    }
  }
  else
  {
    part_aniso_agglo_data->agglomerationLinesIdx      = NULL;
    part_aniso_agglo_data->agglomerationLines         = NULL;
    part_aniso_agglo_data->agglomerationLinesIdx_size = 0;
  }

  /* Free array */
  free(cellCellTmpIdx);
  free(cellCellTmp);
  free(sizes);
  free(arrayOfFineAnisotropicCompliantCells);

  free(isOnFineBndInit);

  if(anisotropicOption[1] == 1){
    free(agglomerationLinesInit);
    free(agglomerationLinesInitIdx);
  }

  free(volume);
  free(centercell);
  free(surfaceArea);
  free(surfaceVector);
  free(surfaceCenter);

  free (cm->specific_data);
  cm->specific_data = NULL;

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Add aniso_agglo coarse mesh method
 *
 */

void
PDM_coarse_mesh_aniso_agglo_add
(
void
)
{
  PDM_coarse_mesh_method_add ("PDM_COARSE_MESH_ANISO_AGGLO", _coarse_from_aniso_agglo);
}


/**
 *
 * \brief Add isotropic array to current coarse mesh
 *
 * \param [in]   cmId                      Coarse mesh identifier
 * \param [in]   iPart                     Current partition
 *
 * \param [in]  agglomerationLines
 * \param [in]  agglomerationLinesIdx
 * \param [in]  isOnFineBnd
 *
 */

void
PDM_part_coarse_mesh_part_set_anisotropic_info
(
 const int     cmId,
 const int     iPart,
 const int    *agglomerationLinesInit,
 const int    *agglomerationLinesInitIdx,
 const int     agglomerationLinesInitIdx_size,
 const int    *isOnFineBndInit
)
{
  _coarse_mesh_t * cm = PDM_part_coarse_mesh_get_from_id (cmId);

  _coarse_part_t *part_res = cm->part_res[iPart];

  assert(part_res->specific_data == NULL);
  part_res->specific_data = malloc (sizeof(_part_aniso_agglo_data_t));
  _part_aniso_agglo_data_t *part_aniso_agglo_data  = (_part_aniso_agglo_data_t *) part_res->specific_data;

  part_aniso_agglo_data->agglomerationLinesInit         = (int *) agglomerationLinesInit;
  part_aniso_agglo_data->agglomerationLinesInitIdx      = (int *) agglomerationLinesInitIdx;
  part_aniso_agglo_data->agglomerationLinesInitIdx_size = agglomerationLinesInitIdx_size;
  part_aniso_agglo_data->isOnFineBndInit                = (int *) isOnFineBndInit;

  if(0 == 1){
    PDM_printf("Inside PDM_part_coarse_mesh_part_set_anisotropic_info(...),part_aniso_agglo_data->agglomerationLinesInitIdx_size: %i \n", part_aniso_agglo_data->agglomerationLinesInitIdx_size);
  }
}

void
PROCF (pdm_part_coarse_mesh_part_set_anisotropic_info, PDM_PART_COARSE_MESH_PART_SET_ANISOTROPIC_INFO)
(
 int          *cmId,
 int          *iPart,
 int          *agglomerationLinesInit,
 int          *agglomerationLinesInitIdx,
 int          *agglomerationLinesInitIdx_size,
 int          *isOnFineBndInit
)
{
  PDM_part_coarse_mesh_part_set_anisotropic_info(*cmId,
                                                 *iPart,
                                                 agglomerationLinesInit,
                                                 agglomerationLinesInitIdx,
                                                 *agglomerationLinesInitIdx_size,
                                                 isOnFineBndInit);

}



/**
 *
 * \brief Add option for anisotropic mesh agglomeration
 *
 * \param [out]  cmId              Coarse mesh identifier
 * \param [in]   anisotropicOption
 */


void
PDM_part_coarse_mesh_add_option_anisotropic
(
 int        cmId,
 const int* anisotropicOption
)
{
  _coarse_mesh_t * cm       = PDM_part_coarse_mesh_get_from_id (cmId);

  assert(cm->specific_data == NULL);
  cm->specific_data = malloc (sizeof(_aniso_agglo_data_t));
  _aniso_agglo_data_t *aniso_agglo_data = (_aniso_agglo_data_t *) cm->specific_data;

  aniso_agglo_data->anisotropicOption = (int *) anisotropicOption;

  if(0 == 1){
    PDM_printf("PDM_part_coarse_mesh_add_option_anisotropic \n ");
    for (int i=0; i<8; i++){
      PDM_printf("PDM_part_coarse_mesh_add_option_anisotropic[%i] : % i\n ", i, aniso_agglo_data->anisotropicOption[i]);
    }
  }
}

void
PROCF (pdm_part_coarse_mesh_add_option_anisotropic, PDM_PART_COARSE_MESH_ADD_OPTION_ANISOTROPIC)
(
 int        *cmId,
 const int  *anisotropicOption
)
{
  PDM_part_coarse_mesh_add_option_anisotropic(*cmId,
                                              anisotropicOption);
}

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   cmId                      Coarse mesh identifier
 * \param [in]   iPart                     Current partition
 *
 * \param [out]  agglomerationLines
 * \param [out]  agglomerationLinesIdx
 * \param [out]  isOnFineBnd
 *
 */

void
PDM_part_coarse_mesh_part_get_anisotropic_info
(
 const int   cmId,
 const int   iPart,
 int       **agglomerationLines,
 int       **agglomerationLinesIdx,
 int        *agglomerationLinesIdx_size,
 int       **isOnFineBnd
)
{
  PDM_printf("Call of PDM_part_coarse_mesh_part_get_anisotropic_info\n");
  _coarse_mesh_t * cm = PDM_part_coarse_mesh_get_from_id (cmId);

  _coarse_part_t *part_res = NULL;

  if (iPart < cm->nPart) {
    part_res = cm->part_res[iPart];
  }

  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_get error : unknown partition\n");
    exit(1);
  }

  _part_aniso_agglo_data_t *part_aniso_agglo_data                = (_part_aniso_agglo_data_t *) part_res->specific_data;

  if(0 == 1 ){
    PDM_printf("Call of PDM_part_coarse_mesh_part_get_anisotropic_info %i \n", cm->method);
  }

  // Bruno : A comprendre ici pourquoi on a pas pointeur null initiliser normalement dans priv.h
  if(cm->method == 2){

    PDM_printf("get of  agglomerationLines\n");
     *agglomerationLines           = part_aniso_agglo_data->agglomerationLines;
     *agglomerationLinesIdx        = part_aniso_agglo_data->agglomerationLinesIdx;
     (*agglomerationLinesIdx_size) = part_aniso_agglo_data->agglomerationLinesIdx_size;
     *isOnFineBnd                  = part_aniso_agglo_data->isOnFineBnd;
     // PDM_printf("\t\tpart_aniso_agglo_data->agglomerationLinesIdx: %i ,%i, %i, %i\n", part_aniso_agglo_data->agglomerationLinesIdx[0] , part_aniso_agglo_data->agglomerationLinesIdx[1], part_aniso_agglo_data->agglomerationLinesIdx[2], part_aniso_agglo_data->agglomerationLinesIdx[3]);
     // PDM_printf("\t\tpart_aniso_agglo_data->agglomerationLinesIdx_size: %i \n", part_aniso_agglo_data->agglomerationLinesIdx_size);
     // PDM_printf("\tEnd of PDM_part_coarse_mesh_part_get_anisotropic_info\n");
   }
   else
   {
     *agglomerationLines           = NULL;
     *agglomerationLinesIdx        = NULL;
     (*agglomerationLinesIdx_size) = 0;
     *isOnFineBnd                  = NULL;
   }
}

void
PROCF (pdm_part_coarse_mesh_part_get_anisotropic_info, PDM_PART_COARSE_MESH_PART_GET_ANISOTROPIC_INFO)
(
 int          *cmId,
 int          *iPart,
 int          *agglomerationLines,
 int          *agglomerationLinesIdx,
 int          *agglomerationLinesIdx_size,
 int          *isOnFineBnd
)
{
  _coarse_mesh_t * cm = PDM_part_coarse_mesh_get_from_id (*cmId);

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

  _part_aniso_agglo_data_t *part_aniso_agglo_data = (_part_aniso_agglo_data_t *) part_res->specific_data;

  for (int i = 0; i < part_aniso_agglo_data->agglomerationLinesIdx_size; i++)
    agglomerationLinesIdx[i] = part_aniso_agglo_data->agglomerationLinesIdx[i];

  for (int i = 0; i < part_aniso_agglo_data->agglomerationLinesIdx[part_aniso_agglo_data->agglomerationLinesIdx_size-1]; i++)
    agglomerationLines[i] = part_aniso_agglo_data->agglomerationLines[i];

  for (int i = 0; i < part_res->part->nFace; i++)
    isOnFineBnd[i] = part_aniso_agglo_data->isOnFineBnd[i];

  (*agglomerationLinesIdx_size) = part_aniso_agglo_data->agglomerationLinesIdx_size;

}

#endif

#ifdef  __cplusplus
}
#endif
