/* 
 * File:   pdm_isotropic_agglomerator.h
 * Author: Nicolas Lantos
 *
 * Created on November 18, 2017, 1:24 PM
 */

#ifndef __TEST_AGGLOMERATOR_ISOTROPIC_PRIV_H__
#define __TEST_AGGLOMERATOR_ISOTROPIC_PRIV_H__

#include "pdm.h"
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <queue>
using namespace std;

//#if !defined (__hpux) && !defined (_AIX)
//#define PROCF(x, y) x##_
//#else
//#define PROCF(x, y) x
//#endif


//#ifdef  __cplusplus
//extern "C" {
//#endif
  
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

// void agglomerateOneLevel(int *sizes,

//                          int *adjMatrix_row_ptr,
//                          int *adjMatrix_col_ind,
//                          double *adjMatrix_areaValues,
//                          double *volumes,

//                          int *arrayOfFineAnisotropicCompliantCells,

//                          int *isOnFineBnd,
//                          int *array_isOnValley,
//                          int *array_isOnRidge,
//                          int *array_isOnCorner,

//                          int isFirstAgglomeration,
//                          int isAnisotropic,

//                          int *fineCellToCoarseCell,

//                          int *agglomerationLines_Idx,
//                          int *agglomerationLines,

//                          int dimension = 3,
//                          int goalCard = -1,
//                          int minCard = -1,
//                          int maxCard = -3,
//                          int checks = 0,
//                          int verbose = 0);

// void agglomerateOneLevel_v_Paradigma(int *sizes,

//                                      int *adjMatrix_row_ptr,
//                                      int *adjMatrix_col_ind,
//                                      double *volumes,

//                                      int *arrayOfFineAnisotropicCompliantCells,

//                                      int *isOnFineBnd_l,
//                                      int * faceCell,
//                                      double *Face_area,

//                                      int isFirstAgglomeration_int,
//                                      int isAnisotropic_int,

//                                      int *fineCellToCoarseCell,

//                                      int *agglomerationLines_Idx,
//                                      int *agglomerationLines,

//                                      int dimension = 3,
//                                      int goalCard = -1,
//                                      int minCard = -1,
//                                      int maxCard = -3,
//                                      int checks = 0,
//                                      int verbose = 0);

void agglomerate_Isotropic_One_Level_v_2(int *sizes,
                                        int *matrixAdj_CRS_row_ptr,
                                        int *matrixAdj_CRS_col_ind,
                                        double *matrixAdj_CRS_values,
                                        double *volumes,
                                        int *fineCellToCoarseCell,
                                        bool *isFineCellAgglomerated,
                                        unordered_set<int>& isOnValley,
                                        unordered_set<int>& isOnRidge,
                                        unordered_set<int>& isOnCorner,
                                        int *isOnFineBnd,  //size numberOfFineCells
                                        int minCard,
                                        int goalCard,
                                        int maxCard,
                                        int thresholdCard,
                                        bool checks,
                                        bool verbose);

void agglomerate_Isotropic_First_Step(int *sizes,
                                      int *matrixAdj_CRS_row_ptr,
                                      int *matrixAdj_CRS_col_ind,
                                      double *matrixAdj_CRS_values,
                                      double *volumes,
                                      unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseCells,
                                      unordered_map<int, unordered_set<int>>& dict_Coarse_Cells,
                                      unordered_map<int, unordered_set<int>>& dict_Card_Coarse_Cells,
        // numberOfFineCells,
                                      unordered_set<int>& isOnValley,
                                      unordered_set<int>& isOnRidge,
                                      unordered_set<int>& isOnCorner,
        // indCoarseCell,
        // numberOfFineAgglomeratedCells,
                                      bool *isFineCellAgglomerated,
                                      int *isOnFineBnd,
                                      int *fineCellToCoarseCell,
                                      list<unordered_set<int>>& delayedCoarseCells,
                                      int goalCard,
                                      int thresholdCard,
                                      int maxCard);

void agglomerate_Isotropic_CreateDelayedCoarseCells(unordered_map<int, unordered_set<int>> & dict_Coarse_Cells,
                                                    unordered_map<int, unordered_set<int>> & dict_Card_Coarse_Cells,
                                                    list<unordered_set<int>> delayedCoarseCells,
                                                    int &indCoarseCell,
                                                    int *fineCellToCoarseCell);


unordered_map<int, queue<int>*> findSeedViaFrontalMethod(int numberOfInts, int* sizes,
                                                         vector<int> listOfFineCells,
                                                         int* matrixAdj_CRS_row_ptr,
                                                         int* matrixAdj_CRS_col_ind);

int computeDistanceFromSeedTofirstVertexOfDegree2(int seed, unordered_map<int, queue<int>*> dict_Covering_Tree);

unordered_map<int, int> computation_Of_Neighbourhood(int seed, int numberOfOrderOfNeighbourhood,
                                                      int *matrixAdj_CRS_row_ptr, int *matrixAdj_CRS_col_ind,
                                                      int maxCard,
                                                      bool *isFineCellAgglomerated_tmp,
                                                      unordered_set<int>* setOfFineCells = nullptr);
int computeNumberOfCommonFaces(int iFine, int iCoarse,
                               int* matrixAdj_CRS_row_ptr,
                               int* matrixAdj_CRS_col_ind,
                               int* fine_Cell_indices_To_Coarse_Cell_Indices);

int removeSeparatingVertex(int seed, unordered_map<int, queue<int>*> dict_ConnectivityTree,
                           unordered_set<int>& setOfFineCells,
                           int* matrixAdj_CRS_row_ptr,
                           int* matrixAdj_CRS_col_ind, int verbose=0);
//bool checkConnectivity(vector<int> listFineCells, int* matrixAdj_CRS_row_ptr, int* matrixAdj_CRS_col_ind, int verbose = 0);
bool checkConnectivity_w_set(unordered_set<int> listFineCells, int *matrixAdj_CRS_row_ptr, int *matrixAdj_CRS_col_ind, int verbose = 0);



list<unordered_set<int>> partsList(vector<int> seq, int length=0);

unordered_set<int> swapFineCell(int iFineCell, int iOrigineCoarseCell, int iDestinationCoarseCell,
                                 unordered_map<int, unordered_set<int>>& dict_Coarse_Elem,
                                 unordered_map<int, unordered_set<int>>& dict_Card_Coarse_Cells,
                                 unordered_map<int, int>& dict_DistributionOfCardinalOfCoarseElements,
                                 int *fineCellToCoarseCell);


void splitNonConnectedCoarseCell(int& indCoarseElement,
                                 int& numberOfFineAgglomeratedCells,
                                 int& iCoarseCell,
                                 unordered_map<int, unordered_set<int>>& dict_Coarse_Cells,
                                 unordered_map<int, unordered_set<int>>& dict_Card_Coarse_Cells,
                                 unordered_map<int, int>& dict_DistributionOfCardinalOfCoarseElements,
                                 int * matrixAdj_CRS_row_ptr,
                                 int * matrixAdj_CRS_col_ind,
                                 bool * isFineCellAgglomerated,
                                 int* fine_Cell_indices_To_Coarse_Cell_Indices);
list<unordered_set<int>> computeConnectedComponent(
        unordered_set<int>listInitialCoarseCell, int* matrixAdj_CRS_row_ptr, int* matrixAdj_CRS_col_ind);

void createCoarseCell(unordered_set<int> l,
                      unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                      unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                      unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements,
                      int &indCoarseElement,
                      int &numberOfFineAgglomeratedCells_tmp,
                      bool *isFineCellAgglomerated_tmp,
                      int *Fine_Cell_indices_To_Coarse_Cell_Indices,
                      bool isMutable = true,
                      bool isCreationDelayed = false);

void createADelayedCoarseCell(unordered_set<int> l,
                              unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                              unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                              int & indCoarseElement,
                              int* Fine_Cell_indices_To_Coarse_Cell_Indices);

void removeDeletedCoarseCells_v3(unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                                 unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                                 int * fine_Cell_indices_To_Coarse_Cell_Indices,
                                 unordered_set<int> set_removedCoarseCells,
                                 int& numberOfCoarseCells);

unordered_set<int> choice_Of_Agglomerated_Cells(int seed,
                                                 vector<queue<int>>& listOfSeeds,
                                                 unordered_map<int, int>& dict_Neighbours_Of_Seed,
                                                 int *matrixAdj_CRS_row_ptr,
                                                 int *matrixAdj_CRS_col_ind,
                                                 double *matrixAdj_CRS_values,
                                                 double *volumes,
                                                 int goalCard,
                                                 int maxCard,
                                                 bool *isFineCellAgglomerated_tmp,
                                                 int *isOnFineBnd,
                                                 int &numberOfFineAgglomeratedCells_tmp,
                                                 bool isOrderPrimary = false);

void agglomerate_Isotropic_createCoarseCell(unordered_set<int> l,
                                            unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                                            unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                                            unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements,
                                            int &indCoarseCell,
                                            int &numberOfFineAgglomeratedCells,
                                            bool *isFineCellAgglomerated,
                                            int *fineCellToCoarseCell,
                                            list<unordered_set<int>> &delayedCoarseCells,
                                            bool isMutable=true,
                                            bool isCreationDelayed=false);

void agglomerate_Isotropic_Choice_Of_Agglomerated_Cells(int seed,
                                                        vector<queue<int>>& listOfSeeds,
                                                        unordered_map<int, int>&dict_Neighbours_Of_Seed,
                                                        int* matrixAdj_CRS_row_ptr,
                                                        int* matrixAdj_CRS_col_ind,
                                                        double* matrixAdj_CRS_values,
                                                        double* volumes,
                                                        unordered_map<int, unordered_set<int>>& dict_Coarse_Elem,
                                                        unordered_map<int, unordered_set<int>>& dict_Card_Coarse_Cells,
                                                        unordered_map<int, int>& dict_DistributionOfCardinalOfCoarseElements,
                                                        int& indCoarseCell,
                                                        int& numberOfFineAgglomeratedCells,
                                                        bool* isFineCellAgglomerated,
                                                        int* fineCellToCoarseCell,
                                                        list<unordered_set<int>> & delayedCoarseCells,
                                                        int* isOnFineBnd,
                                                        int goalCard,
                                                        int thresholdCard,
                                                        int maxCard);
int agglomerate_Isotropic_Choice_Of_Seed(vector<queue<int>> &listOfSeeds,
                                         int numberOfFineCells,
                                         const bool* isFineCellAgglomerated,
                                         unordered_set<int> isOnRidge,
                                         unordered_set<int> isOnValley);

void remove_Too_Small_Cells_v2(int thresholdCard,
                               int *fineCellIndicesToCoarseCellIndices,
                               int &indCoarseCell,
                               int *matrixAdj_CRS_row_ptr,
                               int *matrixAdj_CRS_col_ind,
                               double *matrixAdj_CRS_values,
                               unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                               unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                               unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements);

void makeSmallCellBigger(unordered_map<int, unordered_set<int>> & dict_Coarse_Elem,
                         unordered_map<int, unordered_set<int>> & dict_Card_Coarse_Cells,
                         int * matrixAdj_CRS_row_ptr,
                         int * matrixAdj_CRS_col_ind,
                         unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements,
                         int& indCoarseCell,
                         int& numberOfFineAgglomeratedCells,
                         bool* isFineCellAgglomerated,
                         int * fineCellToCoarseCell,
                         int minCard,
                         int goalCard,
                         int thresholdCard,
                         bool verbose);

void agglomerate_Isotropic_Correction_Swap(unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                                           unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                                           int *matrixAdj_CRS_row_ptr,
                                           int *matrixAdj_CRS_col_ind,
                                           unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements,
                                           int &indCoarseCell,
                                           int numberOfFineCells,
                                           int *fineCellToCoarseCell,
                                           bool verbose);

void agglomerate_Isotropic_Correction_Too_Big_Cells(unordered_map<int, unordered_set<int>> & dict_Coarse_Elem,
                                                    unordered_map<int, unordered_set<int>> & dict_Card_Coarse_Cells,
                                                    int * matrixAdj_CRS_row_ptr,
                                                    int * matrixAdj_CRS_col_ind,
                                                    unordered_map<int, int> & dict_DistributionOfCardinalOfCoarseElements,
                                                    int * fineCellToCoarseCell,
                                                    int & indCoarseCell,
                                                    int goalCard,
                                                    bool verbose);

void compute_Dicts_From_FineCellIndicesToCoarseCellIndices( int nbOfFineCells,
                                                            int * fineCellIndicesToCoarseCellIndices,
                                                            unordered_map<int, unordered_set<int>> &dict_Coarse_Cells,
                                                            unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                                                            unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements);

void agglomerate_Isotropic_Correction_SplitTooBigCoarseCellInTwo(int Nbsizes,
                                                                 int *sizes,
                                                                 vector<queue<int>> &listOfSeeds,
                                                                 unordered_map<int, unordered_set<int>> &dict_Coarse_Elem,
                                                                 unordered_map<int, unordered_set<int>> &dict_Card_Coarse_Cells,
                                                                 int *matrixAdj_CRS_row_ptr,
                                                                 int *matrixAdj_CRS_col_ind,
                                                                 double *matrixAdj_CRS_values,
                                                                 double *volumes,
                                                                 unordered_map<int, int> &dict_DistributionOfCardinalOfCoarseElements,
                                                                 int &indCoarseCell,
                                                                 int *fineCellToCoarseCell,
                                                                 int numberOfFineAgglomeratedCells,
                                                                 bool *isFineCellAgglomerated,
                                                                 int *isOnFineBnd,
                                                                 int minCard,
                                                                 int maxCard,
                                                                 bool checks,
                                                                 bool verbose);


void agglomerate_Isotropic_CheckConsistancyDictCoarseCells(unordered_map<int, unordered_set<int>>& dict_Coarse_Cells,
                                                           unordered_map<int, unordered_set<int>>& dict_Card_Coarse_Cells,
                                                           int fineCellIndicesToCoarseCellIndices_size,
                                                           int *fineCellIndicesToCoarseCellIndices);

void agglomerate_Isotropic_Second_Step_Correction(int numberOfInts, int *sizes,
                                                  int* matrixAdj_CRS_row_ptr,
                                                  int* matrixAdj_CRS_col_ind,
                                                  double* matrixAdj_CRS_values,
                                                  double* volumes,
                                                  unordered_map<int, unordered_set<int>>& dict_CoarseCells,
                                                  unordered_map<int, unordered_set<int>>& dict_CardCoarseCells,
                                                  unordered_map<int, int> & dict_DistribOfCardOfCoarseCells,

                                                  bool* isFineCellAgglomerated,
                                                  int* fineCellToCoarseCell,
                                                  int* isOnFineBnd,
                                                  int minCard,
                                                  int goalCard,
                                                  int maxCard,
                                                  int thresholdCard,
                                                  bool checks,
                                                  bool verbose);

//#ifdef  __cplusplus
//}
//#endif

#endif // __TEST_AGGLOMERATOR_ISOTROPIC_PRIV_H__
