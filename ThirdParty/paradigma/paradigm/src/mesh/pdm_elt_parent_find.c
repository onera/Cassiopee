/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

//FIXME: pdm_elt_parent_find ne fonctionne pas en 64bit !!!!

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_elt_parent_find.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_block.h"
#include "pdm_sort.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

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


/**
 * \def _find_pairs
 * Search common faces in a distribution
 *
 */
static PDM_g_num_t
_find_pairs
(
const int         *IdxFace,
const int         *data,
const int          nFac,
      PDM_g_num_t *connect,
      PDM_g_num_t  iAbsFace
)
{
  /*
   * Make array of already treated face
   */
  int *AlreadyTreat = (int * ) malloc( sizeof(int *) * nFac);
  for(int i=0; i < nFac; i++){
    AlreadyTreat[i] = -1;
  }

  /*
   * Begin with the first face of the array
   */
  for(int iPos=0; iPos < nFac; iPos++){

    /*
     * Verbose
     */
    // printf("---------------------------- : %d - %d - %d\n", iPos, IdxFace[iPos ], data[IdxFace[iPos ]]);

    /*
     * Face deja traité ??
     */
    if(AlreadyTreat[iPos] != 1){

      int curFac = IdxFace[iPos ];
      int nVtx1  = data[curFac+2];

      /** Prepare ElemCon **/
      int *ElmCon1 = (int *) malloc( sizeof(int *) * nVtx1); /* To sort out of loop */
      int  tKey1   = 0;
      for(int iVtx=0; iVtx < nVtx1; iVtx++){
        ElmCon1[iVtx] = data[curFac+3+iVtx];
        tKey1        += data[curFac+3+iVtx];
      }
      // Normalement c'est 64 bit ici !!!
      // PDM_sort_long(ElmCon1, 0, nVtx1-1);
      _quickSort_int(ElmCon1, 0, nVtx1-1);

      if(0 == 1){
        for(int iVtx=0; iVtx<nVtx1; iVtx++){
          printf("ElmCon1[%d] : %d \n", iVtx, ElmCon1[iVtx]);
        }
      }

      /** Loop on candidates faces **/
      for(int iNex=iPos+1; iNex < nFac; iNex++){

        /*
         * Verbose
         */
        // printf("++++++++++++++++++++++++++++++++++ : %d - %d - %d \n", iNex, IdxFace[iNex ], data[IdxFace[iNex ]]);

        /*
         * On sait que la prochaine a traiter est forcement plus grande
         */
        int nexFac = IdxFace[iNex ];
        int nVtx2  = data[nexFac+2];

        // printf("+++++++++++++++++++ : %d - %d \n", nVtx1, nVtx2);
        /*
         * First sort - Number of Vertex is different -> Si Vtx different c'est pas la meme
         */
        if(nVtx1 == nVtx2){

          /*
           * Allocate memory and copy
           */
          int *ElmCon2 = (int *) malloc( sizeof(int *) * nVtx1);
          int  tKey2   = 0;
          for(int iVtx=0; iVtx < nVtx1; iVtx++){
            ElmCon2[iVtx] = data[nexFac+3+iVtx];
            tKey2        += data[nexFac+3+iVtx];
          }

          /*
           * Sort ElemCon2
           */
          _quickSort_int(ElmCon2, 0, nVtx2-1);

          if(0 == 1){
            for(int iVtx=0; iVtx < nVtx1; iVtx++){
              printf("ElmCon2[%d] : %d \n", iVtx, ElmCon2[iVtx]);
            }
          }

          /** Compare **/
          int isSameFace = 1;
          for(int iVtx=0; iVtx < nVtx1; iVtx++){
             if(ElmCon1[iVtx] != ElmCon2[iVtx]){
                isSameFace = -1;
                break;
             }
          }

          /*
           * Verbose
           */
          if(0 == 1){
            printf("isSameFace : %d \n", isSameFace);
            printf("tKey : %d / %d \n", tKey1, tKey2);
            assert(tKey1 == tKey2);
          }

          /*
           * Fill the Face data if Face is the same
           */
          if(isSameFace == 1){

            int flagCur = data[curFac+1];
            int flagNex = data[nexFac+1];

            assert(flagCur != flagNex);

            if(flagCur == -1){
              // connect[data[curFac]] = data[nexFac];
              // connect[data[curFac]-offsetconnect] = data[nexFac];
              connect[iAbsFace] = data[nexFac];
              // printf("Cas 1 : p[%d] = %d  ---- %d \n", data[curFac]-1, data[nexFac], offsetconnect);
            }
            else if(flagNex == -1){
              // connect[data[nexFac]] = data[curFac];
              // connect[data[nexFac]-offsetconnect] = data[curFac];
              connect[iAbsFace] = data[curFac];
              // printf("Cas 2 : p[%d] = %d ---- %d \n", data[nexFac]-1, data[curFac], offsetconnect);
            }
            else
            {
              printf("Something strange\n");
              exit(2);
            }

            iAbsFace++;

            // if(nFacApprox[0] <= iAbsFace ){
            //     printf("Realloc :  %d - %d \n", nFacApprox[0], iAbsFace);

            //     nFacApprox[0] *= 2;
            //     // printf("Realloc :  %d - %d - %d \n", nFacApprox[0], nFacApprox2, iAbsFace);
            //     connect   = (PDM_g_num_t *) realloc(connect,  sizeof(PDM_g_num_t) * nFacApprox[0] );
            // }

            /*
             * Flags the two faces as treated
             */
            AlreadyTreat[iPos] = 1;
            AlreadyTreat[iNex] = 1;

          } /** End if (isSameFace) **/

          /** Free Cendidate ElemCon **/
          free(ElmCon2);

        } /* Enf if Same Number of Vertex */

      }
      /** Free current ElemCon **/
      free(ElmCon1);

    } /** End If alreadyTreat **/

    /* Boundary management **/
    if(AlreadyTreat[iPos] != 1){
      // printf("Alone face nothing to do \n");

      int curFac  = IdxFace[iPos ];
      int flagCur = data[curFac+1];
      int nVtx1   = data[curFac+2];
      if(flagCur == -1){
        printf("flagCur        : %d \n", flagCur);
        printf("curFac         : %d \n", curFac);
        printf("data[curFac];  : %d \n", data[curFac]);

        for(int iVtx=0; iVtx < nVtx1; iVtx++){
          printf("iVtx[%d] : %d \n", iVtx, data[curFac+3+iVtx]);
        }

      }

      assert(flagCur == 1);
    }

  }

  /** Boundary management - Temporary **/


  /** Free **/
  free(AlreadyTreat);

  return iAbsFace;

}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Find parent in a set of elements
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     elt_def_idx          Element definition index
 * \param [in]     elt_def              Element definition
 * \param [in]     n_elt_to_find        Number of elements to find
 * \param [in]     elt_to_find_def_idx  Element to find definition index
 * \param [in]     elt_to_find_def      Element to find definition
 * \param [in]     comm                 MPI Communicator
 * \param [inout]  parent               Parent element of found element, 0 otherwise

 */

void
PDM_elt_parent_find
(
 const int           dnelt,
 const int          *elt_def_idx,
 const PDM_g_num_t  *elt_def,
 const int           dnelt_to_find,
 const int          *elt_to_find_def_idx,
 const PDM_g_num_t  *elt_to_find_def,
 const PDM_MPI_Comm  comm,
       PDM_g_num_t  *parent
)
{
  if (sizeof(int) != sizeof(PDM_g_num_t)) {


    printf("pdm_elt_parent_find : Erreur : Cette fontion ne fonctionne pas en 64bit\n");
    exit(1);

  }

  //FIXME: dNFace

  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(comm, &myRank);
  PDM_MPI_Comm_size(comm, &nRank);

  /* Compute distribution for element */

  PDM_g_num_t* elt_distrib = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t  _dnelt      = (PDM_g_num_t) dnelt;

  PDM_MPI_Allgather((void *) &_dnelt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&elt_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  elt_distrib[0] = 1;

  for (int i = 1; i < nRank+1; i++) {
    elt_distrib[i] +=  elt_distrib[i-1];
  }

  /* Verbose */
  if (1 == 0) {
    PDM_printf("elt_distrib : "PDM_FMT_G_NUM,  elt_distrib[0]);
    for (int i = 1; i < nRank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, elt_distrib[i]);
    }
    PDM_printf("\n");
  }

  /* Compute distribution for element to find */
  PDM_g_num_t* elt_to_find_distrib = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t  _dnelt_to_find      = (PDM_g_num_t) dnelt_to_find;

  PDM_MPI_Allgather((void *) &_dnelt_to_find,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&elt_to_find_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  elt_to_find_distrib[0] = 1;

  for (int i = 1; i < nRank+1; i++) {
    elt_to_find_distrib[i] +=  elt_to_find_distrib[i-1];
  }

  /* Verbose */
  if (1 == 0) {
    PDM_printf("elt_distrib : "PDM_FMT_G_NUM,  elt_to_find_distrib[0]);
    for (int i = 1; i < nRank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, elt_to_find_distrib[i]);
    }
    PDM_printf("\n");
  }

  /* Compute */
  PDM_elt_parent_find_from_distrib(elt_distrib,
                                   elt_def_idx,
                                   elt_def,
                                   elt_to_find_distrib,
                                   elt_to_find_def_idx,
                                   elt_to_find_def,
                                   comm,
                                   parent);

  /* Free */
  free(elt_distrib);
  free(elt_to_find_distrib)  ;

}

/**
 * \brief Find parent in a set of elements
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     elt_def_idx          Element definition index
 * \param [in]     elt_def              Element definition
 * \param [in]     n_elt_to_find        Number of elements to find
 * \param [in]     elt_to_find_def_idx  Element to find definition index
 * \param [in]     elt_to_find_def      Element to find definition
 * \param [in]     comm                 MPI Communicator
 * \param [inout]  parent               Parent element of found element, 0 otherwise

 */

void
PDM_elt_parent_find_from_distrib
(
 const PDM_g_num_t  *elt_distrib,
 const int          *elt_def_idx,
 const PDM_g_num_t  *elt_def,
 const PDM_g_num_t  *elt_to_find_distrib,
 const int          *elt_to_find_def_idx,
 const PDM_g_num_t  *elt_to_find_def,
 const PDM_MPI_Comm  comm,
       PDM_g_num_t  *parent
)
{
  if (sizeof(int) != sizeof(PDM_g_num_t)) {


    printf("pdm_elt_parent_find : Erreur : Cette fontion ne fonctionne pas en 64bit\n");
    exit(1);

  }
  /* Question Eric :    */
  /* 1)  elt_to_find_distrib -> pas en gnum ???? */
  /* 2)  Distribution commence a 1 ou 0 ????, */

  /* elt_distrib = Distribution de face */
  /* elt_def_idx ~= dFaceVtxIdx         */
  /* elt_def     ~= dFace               */

  printf("PDM_elt_parent_find_from_distrib \n");

  int myRank;
  int nRank;

  PDM_MPI_Comm_rank(comm, &myRank);
  PDM_MPI_Comm_size(comm, &nRank);

  /* Build distributed hash tab (key = sum of integer used for element definition */

  /*
   * Allocate Disctribute Hash Key
   */
  PDM_g_num_t nFacApprox = (elt_distrib        [myRank+1] - elt_distrib        [myRank]);
  nFacApprox            += (elt_to_find_distrib[myRank+1] - elt_to_find_distrib[myRank]);

  // nFacApprox *= 2;

  int nDataApprox = nFacApprox*8*2;

  /** Allocate (SurDim the part_data and part_stride ) **/
  int*         part_data = (int *) malloc( sizeof(int) * nDataApprox );
  int*         part_stri = (int *) malloc( sizeof(int) * nFacApprox  );
  PDM_g_num_t* LNToGN    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * nFacApprox );
  PDM_g_num_t* connect   = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * nFacApprox );
  // PDM_g_num_t* connect   = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 18 );

  /** Initialisation of some data **/
  int nFace  = 0;
  int nData  = 0;

  /* Prepare hash table with current element */
  // for (PDM_g_num_t i = elt_distrib[myRank]; i < elt_distrib[myRank+1]; i++) {

  int dNElmts = (int) (elt_distrib[myRank+1]-elt_distrib[myRank]);

  /* Loop over all elements on current process */
  for(PDM_g_num_t iElmt = 0 ; iElmt < dNElmts; iElmt++) {

    if(0 == 1){
      printf("-----> iElmt : "PDM_FMT_G_NUM" \n", iElmt);
      printf("%d/"PDM_FMT_G_NUM" : ", nFace, nFacApprox);
    }

    part_data[nData++] = iElmt+elt_distrib[myRank]; //+1;
    part_data[nData++] = 1;

    /* Setup nVtx */
    int nVtx = elt_def_idx[iElmt+1]-elt_def_idx[iElmt] ;
    part_data[nData++] = nVtx;

    PDM_g_num_t iKey = 0;
    for(PDM_g_num_t idxElmt = elt_def_idx[iElmt] ; idxElmt < elt_def_idx[iElmt+1]; idxElmt++) {
      iKey += elt_def[idxElmt];
      part_data[nData++] = elt_def[idxElmt];
    }

    LNToGN[nFace]    = iKey;

    /** Prepare part_data **/
    part_stri[nFace] = 1+1+1+nVtx;

    /** Go to next face **/
    nFace++;

  }

  /* Prepare hash table with current element to find */

  int dNElmtsToFind = (int) (elt_to_find_distrib[myRank+1]-elt_to_find_distrib[myRank]);

  /* Loop over all elements on current process */
  for(PDM_g_num_t iElmt = 0 ; iElmt < dNElmtsToFind; iElmt++) {

    if(0 == 1){
      printf("-----> iElmt : "PDM_FMT_G_NUM" \n", iElmt);
      printf("%d/"PDM_FMT_G_NUM" : ", nFace, nFacApprox);
    }


    // En parallele on va core dans find_pairs
    part_data[nData++] = iElmt+elt_to_find_distrib[myRank];
    // part_data[nData++] = iElmt+1;
    part_data[nData++] = -1;

    /* Setup nVtx */
    int nVtx = elt_to_find_def_idx[iElmt+1]-elt_to_find_def_idx[iElmt] ;
    part_data[nData++] = nVtx;

    PDM_g_num_t iKey = 0;
    for(PDM_g_num_t idxElmt = elt_to_find_def_idx[iElmt] ; idxElmt < elt_to_find_def_idx[iElmt+1]; idxElmt++) {
      iKey += elt_to_find_def[idxElmt];
      part_data[nData++] = elt_to_find_def[idxElmt];
    }

    LNToGN[nFace]    = iKey;

    /** Prepare part_data **/
    part_stri[nFace] = 1+1+1+nVtx;

    /** Go to next face **/
    nFace++;

  }


  /*
   * Create PartToBlock Structure
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &LNToGN,
                                                      NULL,
                                                      &nFace,
                                                      1,
                                                      comm);

  /*
   * Exchange data in one call ( better do multiple -> Eric ??)
   */
  int *BlkStri = NULL;
  int *BlkData = NULL;

  PDM_part_to_block_exch(          ptb,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_VAR,
                                  -1,
                                   &part_stri,
                         (void **) &part_data,
                                   &BlkStri,
                         (void **) &BlkData);

  /*
   *  Get the size of the current process bloc
   */
  int BlkSize = PDM_part_to_block_n_elt_block_get(ptb);


  /*
   * Verbose
   */
  if(0 == 1){
    printf("BlkSize : %d\n", BlkSize);
    // for(int i = 0; i < BlkSize; i++) {
    //   // printf("BlockData[%d]    : %d \n", i, BlkData[i]);
    //   printf("BlkStri[%d] : %d \n", i, BlkStri[i]);
    // }
  }


  /*
   *  Creation of array of diplacement
   */
  int* BlkStriIdx = (int *) malloc( sizeof(int *) * (BlkSize+1) ) ;
  BlkStriIdx[0] = 0;
  for (int i = 0; i < BlkSize; i++) {
    BlkStriIdx[i+1] = BlkStriIdx[i] + BlkStri[i];
  }
  free(BlkStri);

  /* Find parent in distributed hash table */

  PDM_g_num_t dNFace = 0; //FIXME: Attention : dNFace doit probabement etre declare en entier puis caster en long

  int nFacLocApprox = 10;
  int         *IdxFace      = (int         *) malloc( sizeof(int         *) * nFacLocApprox + 1);
  PDM_g_num_t *connectLocal = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t *) * nFacLocApprox    );

  /*
   * Loop over BlkData : Each block correspond to a key
   *              Each key can identify the same faces
   *              Need to solve all conflict -> See _find_pairs
   */
  for(int iBlk=0; iBlk < BlkSize; iBlk++){

    /*
     * Traitement des groupes de données
     * Reminder Blk = [ [iCell, iType, nVtx, i1, i2, ...], ... ]
     */
    int nEntryFace = BlkStriIdx[iBlk+1] - BlkStriIdx[iBlk];

    /*
     * Verbose
     */
    if(0 == 1){
      printf("BlkStri[%d] : %d -> %d \n", iBlk, BlkStriIdx[iBlk], nEntryFace);
      for(int j=BlkStriIdx[iBlk]; j < BlkStriIdx[iBlk+1]; j++){
        printf("BlkData[%d] : %d \n", j, BlkData[j]);
      }
    }

    /*
     * First Pass - Identify face conflict
     */

    int nFac = 0;
    int iFacData = BlkStriIdx[iBlk];
    while(iFacData < BlkStriIdx[iBlk+1]){
      // printf("Number of Vextex for Face : %d \n", BlkData[iFacData+2]);
      iFacData += 3+BlkData[iFacData+2];
      nFac += 1;
    }

    /*
     *  Allocate array of displacement for all faces in conflict
     *         -> Set the adress of the block data of current faces
     *         -> Make a surdim before and not realloc ...
     */
    if(nFac >= nFacLocApprox){
        printf("Realloc IdxFace -> nFac %d // nFacLocApprox : %d \n", nFac, nFacLocApprox);
        nFacLocApprox = nFac;
        IdxFace       = (int *)         realloc( IdxFace     , sizeof(int *        ) * nFac + 1);
        connectLocal  = (PDM_g_num_t *) realloc( connectLocal, sizeof(PDM_g_num_t *) * nFac    );
    }

    /** New **/
    iFacData   = BlkStriIdx[iBlk];
    IdxFace[0] = BlkStriIdx[iBlk];
    nFac       = 0;

    /** All Faces **/
    while(iFacData < BlkStriIdx[iBlk+1]){
      IdxFace[nFac+1] = IdxFace[nFac] + BlkData[iFacData+2] + 3;
      // printf("BlkData : %d \n ", BlkData[iFacData+2]);
      // printf("IdxFace[%d] : %d \n ", nFac+1, IdxFace[nFac+1]);
      iFacData     += BlkData[iFacData+2] + 3;
      nFac         += 1;
    }
    /* ************************************************** */


    /*
     * Verbose
     */
    if(0 == 1){
      for(int j=0; j < nFac+1; j++){
        printf("IdxFace[%d] : %d \n", j, IdxFace[j]);
        // printf("DeltaIdxFace[%d] : %d \n", j, IdxFace[j+1]-IdxFace[j]);
      }
    }


    /*
     * Solve conflict between all faces and build connectivity array
     *   WATCH OUT -> iAbsFace is changed
     */


    // dNFace = _find_pairs(IdxFace, BlkData, nFac, connect, elt_to_find_distrib[myRank], dNFace, &nFacApprox);

    PDM_g_num_t iAbsFace = 0;
    iAbsFace = _find_pairs(IdxFace, BlkData, nFac, connectLocal, iAbsFace);


    if(dNFace + iAbsFace > nFacApprox){
        printf("[%d/%d] - Realloc :  "PDM_FMT_G_NUM" - "PDM_FMT_G_NUM" \n", myRank, nRank, nFacApprox, dNFace + iAbsFace);
        nFacApprox *= 2;
        connect   = (PDM_g_num_t *) realloc(connect,  sizeof(PDM_g_num_t) * nFacApprox );
    }

    /* Copy on current array */
    for(int ilocface = 0; ilocface < iAbsFace; ilocface++) {
      connect[dNFace+ilocface] = connectLocal[ilocface];
    }

    dNFace = dNFace + iAbsFace;




  }

  /*
   * Free
   */
  free(IdxFace);
  free(connectLocal);


  /* Verbose */
  if(0 == 1){
    // printf("dElmtTot : %d\n", dNFace);
    printf("[%d/%d] - dNFace : "PDM_FMT_G_NUM" / nFaceApprox : "PDM_FMT_G_NUM" \n", myRank, nRank, dNFace, nFacApprox);

    // for(int i = 0; i < dNFace; i++) {
    //   printf("connect[%d]    : %d \n", i, connect[i]);
      // printf("BlkStri2[%d] : %d \n", i, BlkStri2[i]);
    // }
  }

  /* Realloc the connect array */
  if(dNFace > 0){
    connect   = (PDM_g_num_t *) realloc(connect,  sizeof(PDM_g_num_t) * dNFace );
  }

  /*
   * Generate absolute numerotation of faces
   */
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&dNFace, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  beg_NumAbs -= dNFace;

  /** Compute the distribution of elements amont proc **/
  PDM_g_num_t* FaceDistrib = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dNFace = (PDM_g_num_t) dNFace;
  PDM_MPI_Allgather((void *) &_dNFace,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&FaceDistrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  FaceDistrib[0] = 0;

  for (int i = 1; i < nRank+1; i++) {
    FaceDistrib[i] +=  FaceDistrib[i-1];
  }

  /* Verbose */
  if (0 == 1) {
    printf("beg_NumAbs::Face : "PDM_FMT_G_NUM" \n", beg_NumAbs);
    printf("mesh->face_distrib : "PDM_FMT_G_NUM,  FaceDistrib[0]);
    for (int i = 1; i < nRank+1; i++) {
      printf(" "PDM_FMT_G_NUM, FaceDistrib[i]);
    }
    printf("\n");
  }

  /* Return result in the same order of elt_to_find_def array */
  PDM_g_num_t *LNToGNElem = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * dNFace );
  for (PDM_g_num_t i = FaceDistrib[myRank]; i < FaceDistrib[myRank+1]; i++) {
    int idx = (int) (i-FaceDistrib[myRank]);
    LNToGNElem[idx] = i+1;
  }

  int __dNFace = (int) dNFace;

  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                       &LNToGNElem,
                                                       NULL,
                                                       &__dNFace,
                                                       1,
                                                       comm);


  int *BlkStri2 = NULL;
  int *BlkData2 = NULL;

  /** Exchange connect  **/
  PDM_part_to_block_exch(          ptb2,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_CST,
                                   1,
                                   NULL,
                         (void **) &connect,
                                   &BlkStri2,
                         (void **) &BlkData2);

  /*
   *  Get the size of the current process bloc
   */
  int dElmtTot = PDM_part_to_block_n_elt_block_get(ptb2);


  /*
   * Verbose
   */
  if(0 == 1){
    printf("[%d/%d] - ElmtTot : %d\n", myRank, nRank, dElmtTot);
    // for(int i = 0; i < dElmtTot; i++) {
    //   printf("BlockData[%d]    : %d \n", i, BlkData2[i]);
    //   // printf("BlkStri2[%d] : %d \n", i, BlkStri2[i]);
    // }
  }


  /* Push result into parent array */


  for(int i = 0; i < dElmtTot; i++) {
    parent[i] = BlkData2[i];
  }


  PDM_g_num_t* FaceDistribNew = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));


  FaceDistribNew[0] = 0;
  FaceDistribNew[1] = FaceDistrib[nRank];
  for (int i = 2; i < nRank+1; i++) {
    FaceDistribNew[i] = 0;
  }
  for (int i = 1; i < nRank+1; i++) {
    FaceDistribNew[i] +=  FaceDistribNew[i-1];
  }

  if (0 == 1) {
    printf("FaceDistribNew : "PDM_FMT_G_NUM,  FaceDistribNew[0]);
    for (int i = 1; i < nRank+1; i++) {
      printf(" "PDM_FMT_G_NUM, FaceDistribNew[i]);
    }
    printf("\n");
  }

  PDM_part_to_block_free(ptb );

  PDM_part_to_block_free(ptb2);

  free(LNToGN);
  free(LNToGNElem);
  free(part_data);
  free(part_stri);
  free(BlkData   );
  free(BlkStriIdx);
  free(BlkData2   );
  free(FaceDistrib);
  free(connect);

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
