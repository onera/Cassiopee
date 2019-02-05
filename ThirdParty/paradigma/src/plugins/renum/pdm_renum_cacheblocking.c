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
#include "pdm_renum_cacheblocking.h"
#include "pdm_part_graph.h"
#include "pdm_part_renum.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
 *
 * \brief Perform a cells renumbering with cache blocking
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_renum_cells_cacheblocking
(
_PDM_part_t* ppart
)
{
  int methodface = ppart->renum_face_method;
  
  if(ppart->nPropertyCell != 3) {
    PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : You need to specifie [ nCellPerCacheWanted, isAsynchrone, isVectorisation ] in  renum_properties_cell \n");
  }
  
  int nCellPerCacheWanted = ppart->renum_properties_cell[0];
  int isAsynchrone        = ppart->renum_properties_cell[1];
  int isVectorisation     = ppart->renum_properties_cell[2];
  
  const char *_name = PDM_part_renum_method_face_name_get (methodface);
  
  if (strcmp (_name, "PDM_PART_RENUM_FACE_NONE")) {
   PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : face numbering for cacheblocking need to be set to PDM_PART_RENUM_FACE_NONE \n");
  }
  
  /* Loop over all part of the current process */
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    /* Get current part id */
    _part_t *part = ppart->meshParts[ipart];
    
    PDM_renum_cacheblocking(part, 
                            ppart->split_method, 
                            nCellPerCacheWanted, 
                            isAsynchrone, 
                            isVectorisation);
    
  }
}


/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \brief Add renumbering cacheblocking method in the ppart methods
 *
 */

static void 
_computeIdx
(
int* arrayIdx, 
int  size  
)
{
  int tmp1 = arrayIdx[0];
  int tmp2 = -1;
  arrayIdx[0] = 0;
  for (int i = 0; i < size; i++){
    tmp2 = arrayIdx[i+1];
    arrayIdx[i+1] = arrayIdx[i] + tmp1;
    tmp1 = tmp2;
  }
  
}


/**
 *
 * \brief Prepare subdomain
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout  
 * 
 */
  
static void 
_prepare_subdomain
(
 _part_t          *part, 
 _subpartlayout_t *subpartlayout 
)
{  
  int tmp1     = -1;
  int tmp2     = -1;
  
  /* Init to zero */
  subpartlayout->nFaceInt =  0;
  subpartlayout->nFaceExt =  0;
  
  /* Identify max number of subdomaie for faces */
  for (int iface = 0; iface < part->nFace; iface++){
    tmp1 = PDM_MAX(tmp1, part->faceColor[iface]);
    
    int iCell2 = part->faceCell[2*iface + 1];
    if(iCell2 < 1){
      subpartlayout->nFaceExt = subpartlayout->nFaceExt + 1;
    }
    else{
      subpartlayout->nFaceInt = subpartlayout->nFaceInt + 1;
    }
  }
  
  /* Identify max number of subdomaie for faces */
  for (int iface = 0; iface < part->nFace; iface++){
    tmp1 = PDM_MAX(tmp1, part->faceColor[iface]);
  }
  
  /* Identify max number of subdomain for cells  */
  for (int icell = 0; icell < part->nFace; icell++){
    tmp2 = PDM_MAX(tmp2, part->cellColor[icell]);
  }
  
  /* Assert that subdomain are the same size for face and cell */
  if(tmp1 != tmp2)
  {
    PDM_error(__FILE__, __LINE__, 0, "Inconsistency betwenn cell and face tag \n");
  }
  
  subpartlayout->nSdom = tmp1 + 1;
  
  /* Verbose */
  if(0 == 1)
  {
    PDM_printf("nSdom     : " PDM_FMT_L_NUM" \n ", subpartlayout->nSdom);
    PDM_printf("nFaceExt  : " PDM_FMT_L_NUM" \n ", subpartlayout->nFaceExt);
    PDM_printf("nFaceInt  : " PDM_FMT_L_NUM" \n ", subpartlayout->nFaceInt);
  }
  
}


/**
 *
 * \brief Compute subdomain
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout  
 * 
 */

static void 
_compute_subdomain
(
 _part_t          *part, 
 _subpartlayout_t *subpartlayout 
)
{
  int nSdom = subpartlayout->nSdom;
  
  /* Allocate */
  subpartlayout->cellTileIdx    = (int *) malloc (sizeof(int) * (nSdom + 1) );
  subpartlayout->faceTileIdx    = (int *) malloc (sizeof(int) * (nSdom + 1) );
  subpartlayout->faceBndTileIdx = (int *) malloc (sizeof(int) * (nSdom + 1) );
  
  /* Initialise the Idx array  */
  for (int iSub = 0; iSub < nSdom+1; iSub++){
    subpartlayout->cellTileIdx   [iSub] = 0;
    subpartlayout->faceTileIdx   [iSub] = 0;
    subpartlayout->faceBndTileIdx[iSub] = 0;
  }
  
  /* Loop on internal faces */
  for (int iface = 0; iface < subpartlayout->nFaceInt; iface++){
    int color = part->faceColor[iface];
    subpartlayout->faceTileIdx[color]++;
  }
  
  /* Loop on external faces */
  for (int iface = subpartlayout->nFaceInt; iface < part->nFace; iface++){
    int color = part->faceColor[iface];
    subpartlayout->faceBndTileIdx[color]++;
  }
  
  /* Loop on cell */
  for (int icell = 0; icell < part->nCell; icell++){
    int color = part->cellColor[icell];
    subpartlayout->cellTileIdx[color]++;
  }
  
  /* Switch to Idx array */
  _computeIdx(subpartlayout->faceTileIdx   , nSdom+1);
  _computeIdx(subpartlayout->faceBndTileIdx, nSdom+1);
  _computeIdx(subpartlayout->cellTileIdx   , nSdom+1);
  
  /* Panic verbose */
  if(0 == 1)
  {
    /* Debug faceTileIdx */
    PDM_printf(" --------------------- faceTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int iface = subpartlayout->faceTileIdx[iSub]; iface < subpartlayout->faceTileIdx[iSub+1]; iface++){
        printf("iface : %i \n", iface);
      }
    }
    
    /* Debug faceBndTileIdx */
    PDM_printf(" --------------------- faceBndTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int iface = subpartlayout->faceBndTileIdx[iSub]; iface < subpartlayout->faceBndTileIdx[iSub+1]; iface++){
        printf("iface : %i \n", iface);
      }
    }
    
    /* Debug faceBndTileIdx */
    PDM_printf(" --------------------- cellTileIdx on nSdom = %i \n", nSdom);
    for (int iSub = 0; iSub < nSdom; iSub++){
      PDM_printf("iSub         : %i  \n", iSub);
      for (int icell = subpartlayout->cellTileIdx[iSub]; icell < subpartlayout->cellTileIdx[iSub+1]; icell++){
        printf("icell : %i \n", icell);
      }
    }
    
  }
  
}

/**
 *
 * \brief Compute maskage
 *
 * \param [in]      ppart          Ppart object
 * \param [in]      subpartlayout  
 * 
 */

static void 
_compute_mask
(
 _part_t          *part, 
 _subpartlayout_t *subpartlayout 
)
{
  int nSdom = subpartlayout->nSdom;
  
  /* Allocate */
  subpartlayout->maskTileIdx     = (int *) malloc (sizeof(int) * (nSdom + 1   ) );
  subpartlayout->cellVectTileIdx = (int *) malloc (sizeof(int) * (nSdom + 1   ) );
  
  subpartlayout->maskTileN       = (int *) malloc (sizeof(int) * (nSdom       ) );
  subpartlayout->cellVectTileN   = (int *) malloc (sizeof(int) * (nSdom       ) );
  
  subpartlayout->maskTile        = (int *) malloc (sizeof(int) * (part->nCell ) );
  
  int* alreadyInit               = (int *) malloc (sizeof(int) * (part->nCell ) );
  
  /* Initialize */
  for (int iSub = 0; iSub < nSdom; iSub++){
    subpartlayout->maskTileN[iSub]     = 0;
    subpartlayout->cellVectTileN[iSub] = 0;
  }
  
  for (int iSub = 0; iSub < nSdom+1; iSub++){
    subpartlayout->maskTileIdx[iSub]     = 0;
  }
  
  for (int icell = 0; icell < part->nCell; icell++){
    alreadyInit[icell]     = -1;
  }
  
  free(alreadyInit);
  
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Add renumbering cacheblocking method in the ppart methods
 *
 */

void 
PDM_renum_cacheblocking_ppart_add
(
void
)
{
  PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_CACHEBLOCKING", _renum_cells_cacheblocking);
}

/**
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition 
 * \param [in]   split_method        Split method
 * \param [in]   nCellPerCacheWanted Approximate number of cells on each cache
 * \param [in]   isAsynchrone        [0 : Synchrone/1 : Asynchrone]
 * \param [in]   isVectorisation     [0 : No/1 : Yes]
 *
 */

void 
PDM_renum_cacheblocking
(
 _part_t     *part,
int           split_method,  
int           nCellPerCacheWanted, 
int           isAsynchrone,  
int           isVectorisation 
)
{      
  /* Get nFac and nCel */
  const int nCell = part->nCell;
  const int nFace = part->nFace;
  
  /** Allocate reoerdering/permutation array **/
  int *CellOrder = (int *) malloc (sizeof(int) * nCell);
  int *FaceOrder = (int *) malloc (sizeof(int) * nFace);
  
  /*
   * I/ Coloring the cells with classic graph library
   */
  int *CellCellIdx = NULL;
  int *CellCell    = NULL;   

  /* Compute graph associate to mesh */
  PDM_part_graph_compute_from_face_cell(part,
                                   (int **) &CellCellIdx,
                                   (int **) &CellCell); 
  
  /* 
   * Determine the optimal size for cache blocking 
   *   -> The user specify the desired number of cells on each subdomain
   *   -> An other idea is to specify directly the number of block he want
   */
  int nBlkCacheWanted;
  if(nCellPerCacheWanted == 0){nBlkCacheWanted = 1;}
  else                        {nBlkCacheWanted = PDM_MAX(nCell/nCellPerCacheWanted,1);}
  
  /* Split the graph */
  PDM_part_graph_split(split_method,
                  nBlkCacheWanted,
                  part,
                  CellCellIdx, 
                  CellCell, 
                  (int *) NULL,
                  (int *) NULL,
                  (int **) &part->cellColor);
  
  /*
   * II/ Create a proper cells order :
   *       -> [ [SubDom1] , [SubDom2], ..., [SubDomN] ]
   */
  
  /* Allocate */
  int *nCellPerCache    = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int *nCellPerCacheBeg = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int* partCellIdx      = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  
  /* II-Bis - Asynchrone : We change only number of each sub-domain 
   * in order to have interior sub-domain first and exterior after 
   */
  if(isAsynchrone == 1)
  {
    /* Allocate */
    int *flagHaveBnd     = (int *) malloc(nBlkCacheWanted * sizeof(int));
    int *syncToAsynColor = (int *) malloc(nBlkCacheWanted * sizeof(int));
    
    /* All sub-domain is initialize */
    for(int i = 0; i < nBlkCacheWanted; i++) {
      flagHaveBnd[i] = 0;
    }
    
    /* 
     * Loop on faces : if at least one is border, flagHaveBnd is set to 1 
     */
    for (int i = 0; i < part->nFace; i++) {
      int iCell1 = part->faceCell[2*i    ];
      int iCell2 = part->faceCell[2*i + 1];
      int color1 = part->cellColor[iCell1-1];

      /* Solution 1 */        
      // if(iCell2 > 0 ){
      //   int color2 = part->cellColor[iCell2-1];
      //   int color  = _PDM_part_MIN(color1, color2);
      // }
      // else
      // {
      //   flagHaveBnd[color1] = 1;
      // }
      
      /* Solution 2 */     
      if(iCell2 == 0){flagHaveBnd[color1] = 1;}
      
    }
    
    /*
     * Built the syncToAsynColor array
     */
    int nSdomB = 0;                     /* Number of   Blocking domain */
    int nSdomU = 0;                     /* Number of Unblocking domain */
    int nextColorB = nBlkCacheWanted-1; /* Next Blocking   color       */
    int nextColorU = 0;                 /* Next UnBlocking color       */ 
    
    for(int i = 0; i < nBlkCacheWanted; i++) {
      if(flagHaveBnd[i] == 0)
      {
        syncToAsynColor[i] = nextColorU;
        nextColorU++;
        nSdomU++;
      }
      else
      {
        syncToAsynColor[i] = nextColorB;
        nextColorB--;
        nSdomB++;
      }
    }
    
    /* Dramatic verbose */
    if(0 == 1)
    {
      printf(" ----------- : %i \n", nBlkCacheWanted);
      for (int i = 0; i < nBlkCacheWanted; i++){
        printf("SyncToAsynColor[%i] = %i\n", i, syncToAsynColor[i]);
      }
    }
    
    /*
     * Apply the computed array on cellPart 
     */
    for(int i = 0; i < part->nCell; i++) {
      part->cellColor[i] = syncToAsynColor[part->cellColor[i]];
    }
    
    
    /* Free memory */
    free(flagHaveBnd);
    free(syncToAsynColor);
    
  } /* End asynchrone */

  
  /* Init to Zero */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nCellPerCache[i] = 0;
  }
  
  /* First loop to identify the size of each color */
  for (int i = 0; i < nCell; i++){
    int color = part->cellColor[i]; // A color is a number of partition (output of Metis or Scotch)
    nCellPerCache[color]++;
  }

  /* Compute the index of each Sub-Domain **/
  partCellIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partCellIdx[i + 1] = partCellIdx[i] + nCellPerCache[i];
  }
  

  if(isVectorisation == 0){
    /** First method - Standard **/
    
    /* Reset array */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCache[i] = 0;
    }
    
    for (int i = 0; i < nCell; i++){
      int color  = part->cellColor[i]; 
      int idx    = partCellIdx[color] + nCellPerCache[color];
      
      CellOrder[idx] = i;
      nCellPerCache[color]++;
    }
  }
  else
  {
    /** Second method - Vectorisation on cell on current domain **/
    
    /* Store the current Index to add a vectorisable cell */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCacheBeg[i] = 0;
    }
    
    /* Loop on cells : 
     *     -> If current cell is adjacent to another color : add to end 
     *     -> Else add to begin
     * In fact this function sort the interior cell and exterior cell leads to the following numbering : 
     *         [ [ Interior Cells] , [ Exterior Cells]]
     */
    for (int i = 0; i < nCell; i++){
      
      /* Get color and prepare flag */
      int color = part->cellColor[i];
      int flag  = -1;
      for(int j = CellCellIdx[i]; j < CellCellIdx[i+1]; j++){
        int iCell = CellCell[j];
        if(part->cellColor[iCell] != color){flag = 1;} // Alors on est au bord d'un sous domaine !
      }
      
      if(flag == -1){  // Cell is interior : add to begin
        
        int idx = partCellIdx[color] + nCellPerCacheBeg[color];
        
        CellOrder[idx] = i;
        nCellPerCacheBeg[color]++;
        
      }
      else{ // Cell is exterior : add to end
        
        int idx = partCellIdx[color] + nCellPerCache[color]-1;
        
        CellOrder[idx] = i;   
        nCellPerCache[color]--;
        
      }
      
      /* Panic verbose */
      // if(0 == 1){
      //   PDM_printf("Begin : %i %i %i \n", color, idx, nCellPerCacheBeg[color]);
      // }
      
    }
  }
  
  /* 
   * Verbose 
   */
  if(0 == 1)
  {
    printf(" ----------- : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> CellOrder[%i]  = %i\n", i, CellOrder[i]);
    }
  }
  
  /*
   * We apply renumbering here because we need it for faces renumbering 
   */
  PDM_part_reorder_cell(part, CellOrder);
  
  /*
   * III/ Create a proper faces order :
   *       -> [ [SubDom1Int/SubDom1Ext] , [SubDom2Int/SubDom2Ext], ..., [SubDomNInt/SubDomNExt],  ]
   *      Renumbering the faces according to the sub block new ordering 
   */

  /* Allocate */    
  int *nFacePerCache    = (int *) malloc(nBlkCacheWanted * sizeof(int));
  int *nFaceBndPerCache = (int *) malloc(nBlkCacheWanted * sizeof(int));
  
  /* Init array */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }
  
  /* 
   * First pass : 
   *      -> the face is associate with the subdomain that have the lowest color number 
   */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];
    
    int color1 = part->cellColor[iCell1-1];
    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      nFacePerCache[color]++;
    }
    else
    {
      nFaceBndPerCache[color1]++;
    }
    
    /* Dramatic test */
    if(iCell1 < 0 ){ 
      printf("PPART internal error \n");
      exit(1); 
    }
  }
   
  /* Second pass : 
   *      -> the face is associate with the subdomain that have the lowest color number 
   */
  
  /* Allocate */
  int* partFaceIdx    = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  int* partFaceBndIdx = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  
  /* Initialise */
  partFaceIdx[0]    = 0;
  partFaceBndIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partFaceIdx[i + 1]    = partFaceIdx[i]    + nFacePerCache[i];
    partFaceBndIdx[i + 1] = partFaceBndIdx[i] + nFaceBndPerCache[i];
  }
  
  for (int i = 0; i < nBlkCacheWanted+1; i++){
    partFaceBndIdx[i] = partFaceBndIdx[i] + partFaceIdx[nBlkCacheWanted];
  }
  
  /* 
   * Verbose 
   */
  if(0 == 1)
  {
    printf(" ----------- : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nBlkCacheWanted; i++){
      printf("~> %i\n", i);
      printf(" =====> partFaceIdx    %i\n", partFaceIdx[i + 1]);
      printf(" =====> partFaceBndIdx %i\n", partFaceBndIdx[i + 1]);
      
      printf(" =====> nFacePerCache    %i\n", nFacePerCache[i]);
      printf(" =====> nFaceBndPerCache %i\n", nFaceBndPerCache[i]);
    }
  }
      
  /* Reset Idx */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }

  /* Determine for each subDomain the associate faces */    
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];
    
    int color1 = part->cellColor[iCell1-1];
    
    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      int idx    = partFaceIdx[color]    + nFacePerCache[color];
      
      FaceOrder[idx] = i;
      nFacePerCache[color]++;
    }
    else
    {
      int color  = color1;
      int idxBnd = partFaceBndIdx[color] + nFaceBndPerCache[color];
      
      FaceOrder[idxBnd] = i;
      nFaceBndPerCache[color]++;
    }
    
    if(iCell1 < 0 ){ 
      printf("PPART internal error \n");
      exit(1); 
    }
    
  } /* End file FaceOrder */
  
  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrder[%i] = %i \n", i, FaceOrder[i]);
    }
  }
  
  /*
   * Apply renumbering
   */
  PDM_part_reorder_face(part, FaceOrder);
  
  /*
   * Save in array the color of each faces 
   */
  part->faceColor = (int *) malloc (sizeof(int) * nFace);
  
  for (int i = 0; i < nBlkCacheWanted; i++){
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
    for (int iface = partFaceBndIdx[i]; iface < partFaceBndIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
  }
  
      
  /*
   * IV/ Create a proper faces order for vectorisation :
   *         -> The idea is to found a pool of faces that not pointed same cells : "independant faces"
   *         Rmk : Boundary faces are not vecotrised
   */
  for (int i = 0; i < part->nFace; i++) {
    FaceOrder[i] = i;
  }
  
  /* Allocate */
  int *flagCell = (int *) malloc (sizeof(int) * nCell);
  int *flagFace = (int *) malloc (sizeof(int) * nFace);
  
  /* Core loop of reordenencing */
  for (int i = 0; i < nBlkCacheWanted; i++){
    
    /* Get local size of SubDomain */
    int nFacSub = nFacePerCache[i];
    // int nCelSub = nCellPerCache[i];
    
    /* Reset for one pacquet all flag */
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      flagFace[iface] = -1;
    }
    
    int nFacTreated = 0;
    // printf("nFacSub : %i \n", nFacSub);
    while(nFacTreated != nFacSub){
      
      /* Reset for one pacquet all flag */
      // for (int iCel = 0; iCel < nCell; iCel++){
      //   flagCell[iCel] = -1;
      // }
      // TEST : A remplacer par : (Beacoup moins couteux)
      for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
        int iCell1 = part->faceCell[2*iface  ];
        int iCell2 = part->faceCell[2*iface+1];
        flagCell[iCell1-1] = -1;
        flagCell[iCell2-1] = -1;
      }
      
      /* Dramatic verbose */
      if(0 == 1){
        printf("nFacTreated/nFacSub : %i -> %i \n", nFacTreated, nFacSub);
      }
      
      /* Loop on face */
      for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
        
        if(flagFace[iface] == -1){
          int iCell1 = part->faceCell[2*iface  ];
          int iCell2 = part->faceCell[2*iface+1];
          
          int t1 = flagCell[iCell1-1];
          int t2 = flagCell[iCell2-1];
          
          // printf("t1/t2 : %i/%i \n", t1, t2);
          if( (t1 == -1) && (t2 == -1)){
            flagCell[iCell1-1] = 1;
            flagCell[iCell2-1] = 1;
            
            int bFac = partFaceIdx[i];
            
            FaceOrder[bFac+nFacTreated] = iface;
            
            flagFace[iface] = 1;
            
            nFacTreated++;
          }
        } /* End Face loop */
      } /* End While */
    } /* End SubDom loop */
     
  }
  
  /*
   * Apply renumbering -> Attention au tableau en plus à reordencer !
   */
  PDM_part_reorder_face(part, FaceOrder);
  
  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrderVect[%i] = %i \n", i, FaceOrder[i]);
    }
  }
  
  /* Free memory */
  free(CellOrder);
  free(FaceOrder);
  free(CellCellIdx);
  free(CellCell);
  free(nCellPerCacheBeg);
  free(nCellPerCache);
  free(nFacePerCache);
  free(nFaceBndPerCache);
  free(partFaceIdx);
  free(partFaceBndIdx);
  free(flagCell);
  free(flagFace);

}

/*
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]   part                Mesh Partition 
 *
 */
void 
PDM_renum_cacheblocking_compute_loop_array
(
 _part_t     *part
)
{  
  if(part->faceColor == NULL){PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : You need to faceColor \n");}
  if(part->cellColor == NULL){PDM_error(__FILE__, __LINE__, 0, "_renum_cells_cacheblocking Error : You need to cellColor \n");}
  
  /* Prepare sub-domain */
  _prepare_subdomain(part, part->subpartlayout);
  
  /* Compute sub-domain */
  _compute_subdomain(part, part->subpartlayout);
  
  /* Compute mask */
  _compute_mask(part, part->subpartlayout);
  
  
}

#ifdef  __cplusplus
}
#endif
