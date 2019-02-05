/* 
 * File:   pdm_anisotropic_agglomerator.cxx
 * Author: Nicolas Lantos
 *
 * Created on November 18, 2017, 1:24 PM
 */

#include "pdm_anisotropic_agglomerator.h"
#include <unordered_map>
#include <map>
#include <limits>
#include <deque>
#include <forward_list>
#include <iostream>
#include <assert.h>
#include <set>
#include <algorithm>

using namespace std;

template<typename A, typename B> std::pair<B,A> flip_pair(const pair<A,B> &p)
{
    return pair<B,A>(p.second, p.first);
};

template<typename A, typename B, template<class, class, class...> class M, class... Args> multimap<B, A> flip_map(const M<A,B,Args...> &src)
{
    multimap<B,A> dst;
    transform(src.begin(), src.end(), inserter(dst, dst.begin()), flip_pair<A,B>);
    return dst;
};


bool computeAnisotropicLine(int* sizes,
                           int * adjMatrix_row_ptr,
                           int *adjMatrix_col_ind,
                           double * adjMatrix_areaValues,
                           int * arrayOfFineAnisotropicCompliantCells,
                           int * agglomerationLines_Idx,
                           int * agglomerationLines,
                           bool verbose){

    // Rmk: fonction couteuse Il y a un tri d'un dictionnaire!

    int numberOfFineCells = sizes[0];
//    int adjMatrix_row_ptr_size = numberOfFineCells+1;
//    int adjMatrix_col_ind_size = sizes[1];
//    int adjMatrix_areaValues_size = sizes[1];

    int numberOfFineAnisotropicCompliantCells = sizes[7];
    int numberOfAnisotropicLinesPOne_size = sizes[8];  // numberOfFineAnisotropicCompliantCells at the  beginning
    int agglomerationLines_size = sizes[9];            // numberOfFineAnisotropicCompliantCells at the  beginning

    bool isAnisotropicLines = false;
//    cout<<"numberOfFineCells "<<numberOfFineCells<<endl;
//    cout<<"adjMatrix_row_ptr_size "<<adjMatrix_row_ptr_size<<endl;
//    cout<<"adjMatrix_col_ind_size "<<adjMatrix_col_ind_size<<endl;
//    cout<<"numberOfFineAnisotropicCompliantCells "<<numberOfFineAnisotropicCompliantCells<<endl;
//    cout<<"numberOfAnisotropicLinesPOne_size "<<numberOfAnisotropicLinesPOne_size<<endl;
//    cout<<"agglomerationLines_size "<<agglomerationLines_size<<endl;
    // self._listOfSetAnisotropicCompliant[iLevel] contains the set of anisotropic compliant cells (i.e. Hexa and
    //                                                                                                prism)
    // dictAnisotropicCell
    // isAnisotropic

    // TODO Au lieu de faire un tableau avec toutes les cellules, les Hexa et les prismes sont suffisants!
    double *maxArray = new double[numberOfFineCells];
    unordered_map<int, double> dictAnisotropicCell;// keys are the ratio Max to average (ratioArray[iCell]) and value
    //                                                 the (global) index of the cell.
    //dictAnisotropicCell = dict()

    // Process of every fine cells:
    // TODO a priori on pourrait faire le calcul uniquement pour les cellules hexe/prism
    // for iCell in xrange(numberOfFineCells):
    int cellNumber;
    int ind, indPOne, indNeighborCell;
    double minWeight, maxWeight, averageWeight, weight;

    for (int iCell = 0; iCell<numberOfFineAnisotropicCompliantCells; ++iCell){
        cellNumber = arrayOfFineAnisotropicCompliantCells[iCell];
        ind = adjMatrix_row_ptr[cellNumber];
        indPOne = adjMatrix_row_ptr[cellNumber + 1];

        minWeight = numeric_limits<double>::max();
        maxWeight = 0.0;
        averageWeight = 0.0;

        // computation of minWeight, maxWeight and averageWeight for the current cell iCell
        // Process of every faces/Neighbours
        for (int iN=ind; iN<indPOne; ++iN) {
            indNeighborCell = adjMatrix_col_ind[iN];
            if( indNeighborCell != cellNumber) {  // to avoid special case where the boundary value are stored
                weight = adjMatrix_areaValues[iN];
                if (maxWeight < weight) {
                    maxWeight = weight;
                }
                if (minWeight > weight) {
                    minWeight = weight;
                }
            }

            averageWeight += adjMatrix_areaValues[iN] / (indPOne - ind);
        }
        maxArray[cellNumber] = maxWeight;

        // TODO comprendre la raison d'etre des 2 criteres: present chez Mavriplis et DLR
        // Anisotropy criteria for the initial sort
        // if ratioArray[iCell] >= 2.0:
        //     dictAnisotropicCell[ratioArray[iCell]] = iCell
        //
        //     count += 1

        // Anisotropy criteria for the line Admissibility
        if (maxWeight / minWeight >= 4.0) {
            // if iCell in setOfAnisotropicCompliantFineCells:
            // we check if the cell is for lvl 0 hexa or prism
            // or for lvl >0 if is build from hexa or prism.
            // isAnisotropic[iCell] = True
            dictAnisotropicCell[cellNumber] = maxWeight / averageWeight;
        }
    }

//    cout<<"dictAnisotropicCell.empty() "<<dictAnisotropicCell.empty()<<endl;
//    if (!dictAnisotropicCell.empty()) {
//    unordered_map<int, double>::iterator uMapIt;
//    for (uMapIt=dictAnisotropicCell.begin(); uMapIt!=dictAnisotropicCell.end(); uMapIt++)
//    {
//        cout<<uMapIt->first<<" "<<uMapIt->second<<endl;
//    }
//    }
    if (!dictAnisotropicCell.empty()) {
        // There are anisotropic cells

        // We sort the dict w.r.t. the ratio  ????
        // From the biggest "anisotropy" to the smallest
        // orderedDictAnisotropyCell = OrderedDict(sorted(dictAnisotropicCell.items(), key=lambda t: t[1],
        //                                                reverse=True))
        //sortedList_anisotropicCells = sorted(dictAnisotropicCell, key = dictAnisotropicCell.get, reverse = True)

        multimap<double, int> sorted_anisotropicCells = flip_map(dictAnisotropicCell);  // Sorting!
        multimap<double, int>::reverse_iterator rIt;

        int numberOfAnisotropicCells = dictAnisotropicCell.size();
//        cout<<"numberOfAnisotropicCells "<<numberOfAnisotropicCells<<endl;
        int* indexLocalToGlobal = new int[numberOfAnisotropicCells];
        unordered_map<int, int> dictGlobalToLocal;

        // Initialisation of indexLocalToGlobal
        // To work only on anisotropic cells.
        int i =0, index;
        for (rIt = sorted_anisotropicCells.rbegin();rIt!=sorted_anisotropicCells.rend(); rIt++)
        {
            index = rIt->second;
            indexLocalToGlobal[i] = index;
            dictGlobalToLocal[index] = i;
            ++i;
        }

        // isAnisotropic_local: to check if a cell has already been treated
        bool* isAnisotropic_local = new bool[numberOfAnisotropicCells];
        for (int iLoc =0; iLoc< numberOfAnisotropicCells; iLoc++)
        {
            isAnisotropic_local[iLoc]=true;
        }

        // TODO Think of a more pertinent cell???
        int iLocalSeed = 0;

        int iGlobalSeed = indexLocalToGlobal[iLocalSeed];

//        cout<<"iGlobalSeed "<<iGlobalSeed<<endl;

        int lines_size = 0;
        forward_list<deque<int>*> lines;

        int currentCell;

        bool searchOppositeDirection= false;
        int cellAtTheLimitIsoAniso;
        int previousCurrentCell;

        // We try every possible cell as seed:
        while (iLocalSeed < numberOfAnisotropicCells) {
//            cout<<"iLocalSeed "<< iLocalSeed<<endl;
//            cout<<"lines_size "<<lines_size<<endl;
            currentCell = iGlobalSeed;
            deque<int>* dQue;
            dQue = new deque<int>();
            (*dQue).push_back(iGlobalSeed);
            deque<int>::iterator iter;
//            for (iter = (*dQue).begin();iter!=(*dQue).end(); iter++) {
//                cout<<"iter1 "<<*iter<<endl;
//            }
            // Line contains the global numbering of cells
            searchOppositeDirection = false;

            isAnisotropic_local[iLocalSeed] = false;
            cellAtTheLimitIsoAniso = -1;

            while (currentCell != -1) {
                // Ca c'est pour verifier si on part ou non d'une extremite d'une ligne.
                // voir if currentCell == previousCurrentCell:
                previousCurrentCell = currentCell;
//                cout<<"CurrentCell "<<currentCell <<endl;
                ind = adjMatrix_row_ptr[currentCell];
                indPOne = adjMatrix_row_ptr[currentCell + 1];

                // Process neighbours of current cell
                for (int iN = ind; iN < indPOne; ++iN) {
                    indNeighborCell = adjMatrix_col_ind[iN];
//                    cout<<"currentCell "<< currentCell<< " indNeighborCell "<<indNeighborCell<<endl;
                    if (indNeighborCell != currentCell) {
                        // Reminder: if indNeighborCell == currentCell, we are at the boundary of the domain
                        if (adjMatrix_areaValues[iN] > 0.75 * maxArray[currentCell]) {
                            // Search faces with biggest areas
//                            cout<<"(adjMatrix_areaValues[iN] > 0.75 * maxArray[currentCell])"<<endl;
                            if (dictGlobalToLocal.count(indNeighborCell) == 1 and
                                isAnisotropic_local[dictGlobalToLocal[indNeighborCell]]) {

                                // We find an anisotropic neighbour

                                // On veut que la ligne soit definie de maniere contigue:
                                // c'est a  dire qu'on parte d'un bout vers l'autre et
                                // pas du milieu, jusqu'a un bout et puis on repart dans l'autre sens.
                                // TODO le cas test en ne triant pas initialement devra me permettre de tester ca!
                                if (!searchOppositeDirection) {
                                    // General case:
                                    // Correct direction from wall to farfield
                                    (*dQue).push_back(indNeighborCell);
                                    currentCell = indNeighborCell;
                                } else {
                                    (*dQue).push_front(indNeighborCell);
                                    currentCell = indNeighborCell;
                                }
                                // TODO Est-ce le bon moyen de marquer les cellules comme visitees????
                                // isAnisotropic[indNeighborCell] = False
                                isAnisotropic_local[dictGlobalToLocal[indNeighborCell]] = false;
                                break;
                            } else {
                                // We find an isotropic neighbour! (and we are not at the boundary

                                // assert cellAtTheLimitIsoAniso ==-1, "Problem cellAtTheLimitIsoAniso
                                // is overwritten"
                                bool isInDeque = false;
                                for (int index2: (*dQue)) {
                                    if (indNeighborCell == index2) {
                                        isInDeque = true;
                                        break;
                                    }
                                }
                                if (!isInDeque) {
                                    // As isAnisotropic[indNeighborCell] is used to mark/color cells already treated
                                    // we check is the neighbour cell has been treated in the current line
                                    // TODO what happend if the neighbour is in an another agglomerationLine???
                                    cellAtTheLimitIsoAniso = currentCell;
                                }
                            }
                        }
                    }
                }


                // TODO Could it be 3 anisotropic neighbour with area bigger than the criteria?
                if (currentCell == previousCurrentCell) {
                    // in this case, the seed was not at one end of the line.
                    // We have done a part of the line and we now try the other direction
                    if (!searchOppositeDirection) {
                        //assert(dQue.size() > 1); // "Trouble"
                        searchOppositeDirection = true;
                        currentCell = iGlobalSeed;
                    } else {
                        currentCell = -1;
                    }
                }
            }
            //cout<<"Line found: ";
            //for (iter = (*dQue).begin();iter!=(*dQue).end(); iter++) {
            //    cout<<*iter<<" ";
            //}
            //cout<<endl;

            // The current line contains more cells than the seed.
            // We keep it
            // otherwise we rewrite on it.
            if ((*dQue).size() > 1) {
//                cout<<"(*dQue).size() > 1"<<endl;
                //Si on recontre une cellulle isotrope, on inverse l'ordre!
                if ((*dQue).front() == cellAtTheLimitIsoAniso) {
                    //lines[lineNumber] = deque(reversed(lines[lineNumber]));
                    deque<int>* dQ2;
                    dQ2 = new deque<int>();
                    for (auto iFC :(*dQue)){
//                        cout<<"(*dQue) "<<i<<endl;
                        (*dQ2).push_front(iFC);
                    }
                    assert((*dQ2).front() != cellAtTheLimitIsoAniso) ;// "Problem both ends are at the limit"
                    // deletion of dQue which is useless now
                    delete dQue;

                    // Creation of a new line
                    lines.push_front(dQ2);
                    lines_size ++;
                }
                else{
                    // Creation of a new line
                    lines.push_front(dQue);
                    lines_size ++;
                }
            }else{
                delete dQue;
            }

            // On pourrait prendre directement le suivant i.e. iSeed+=1
            // Sauf que potentiellement ce suivant pourrait deja appartenir a une ligne.
            while ((iLocalSeed < numberOfAnisotropicCells) && (! isAnisotropic_local[iLocalSeed])) {
                // while iLocalSeed < numberOfAnisotropicCells and not isAnisotropic[indexLocalToGlobal[iLocalSeed]]:
                iLocalSeed += 1;
            }
            if (iLocalSeed == numberOfAnisotropicCells) {
                continue;
            }
            iGlobalSeed = indexLocalToGlobal[iLocalSeed];
        }
//        cout<<"End of While Loop"<<endl;
//        forward_list<deque<int>*>::iterator fLIter;
//        for (fLIter = lines.begin(); fLIter!=lines.end(); fLIter++)
//        {
//            cout<<(*(*fLIter)).size()<<endl;
//        }
//        cout<<"lines.front().size() "<< (*lines.front()).size()<<endl;

        // The last line may be of size 1 and will not be overwritten.
        if (!lines.empty()) {
            if ((*lines.front()).size() == 1) {

                if (verbose) {
                    cout << "lines.pop() " << (*lines.front()).front() << endl;
                }
                lines.pop_front();
            }
        }
//        cout<<"Conversion of deque to array"<<endl;
        int numberOfAgglomerationLines = lines_size;

        if(numberOfAgglomerationLines == 0) {
//            numberOfAgglomerationLines = 1;
            agglomerationLines_Idx[0] = 0;
            agglomerationLines_Idx[1] = 0;
            agglomerationLines_size =0;
            numberOfAnisotropicLinesPOne_size =2;
        }
        else{

            int numberOfFCellsInAgglomerationLines = 0;

            agglomerationLines_Idx[0] = 0;
            forward_list<deque<int>*>::iterator fLIt;
            int iLines =1;
            for (fLIt = lines.begin();fLIt!=lines.end(); fLIt++)
            {

                agglomerationLines_Idx[iLines] = (*(*fLIt)).size()+numberOfFCellsInAgglomerationLines;
                int jCount=0;

                for (auto i2 :(*(*fLIt))){
                    agglomerationLines[jCount + numberOfFCellsInAgglomerationLines]=i2;
                    jCount++;
                }
                numberOfFCellsInAgglomerationLines +=(*(*fLIt)).size();
                iLines++;
            }

            // Deletion of pointer to deque:
            for (fLIt = lines.begin();fLIt!=lines.end(); fLIt++)
            {
                delete (*fLIt);
            }

            numberOfAnisotropicLinesPOne_size = iLines;
            agglomerationLines_size = numberOfFCellsInAgglomerationLines;
            isAnisotropicLines = true;
        }


        delete[] indexLocalToGlobal;
        delete[] isAnisotropic_local;
    }else{
        agglomerationLines_Idx[0] = 0;
        agglomerationLines_Idx[1] = 0;
        agglomerationLines_size =0;
        numberOfAnisotropicLinesPOne_size =2;

    }

//    sizes[0] = numberOfFineCells;  // not modified!
//    sizes[1] = adjMatrix_col_ind_size;  // not modified!
//    sizes[2] = numberOfFineAnisotropicCompliantCells;  // not modified!
    sizes[8] = numberOfAnisotropicLinesPOne_size;
    sizes[9] = agglomerationLines_size;

    delete[] maxArray;
    return isAnisotropicLines;
}


void agglomerate_Anisotropic_One_Level_without_list_lines(int *sizes,
                                                          int *fineAgglomerationLines_array_Idx, int *fineAgglomerationLines_array,
                                                          int *fineCellToCoarseCell,
                                                          bool *isFineCellAgglomerated,
                                                          int *AnisotropicCompliantCoarseCells_array) {

    int fineAgglomerationLines_array_Idx_size = sizes[0];
    int numberOfFineCells = sizes[1];
    int numberOfFineAgglomeratedCells = sizes[2];
    int indCoarseCell = sizes[3];
    int AnisotropicCompliantCoarseCells_size = sizes[4];

    set<int> setOfAnisotropicCompliantCoarseCells;

    // Assert that lines exists:
    if ( (fineAgglomerationLines_array_Idx_size == 2) && (fineAgglomerationLines_array_Idx[1] == 0)){
        return ;
    }

    //new_Lines_index = np.zeros(len(fineAgglomerationLines_array_Idx), dtype=np.int32)
    int* new_Lines_index = new int[fineAgglomerationLines_array_Idx_size];
    for(int iNLI=0; iNLI<fineAgglomerationLines_array_Idx_size; iNLI++)
    {
        new_Lines_index[iNLI]=0;
    }
//    new_Lines_index[0] = 0;


    int fineAgglomerationLines_array_size = fineAgglomerationLines_array_Idx[fineAgglomerationLines_array_Idx_size-1]-1;
    //cout<<"fineAgglomerationLines_array_size "<<fineAgglomerationLines_array_size<<endl;

    //new_Lines_array = np.zeros(fineAgglomerationLines_array_Idx[-1]-1, dtype=np.int32)
    int* new_Lines_array = new int[fineAgglomerationLines_array_size];

    int iCountCoarseCell = 0;
    int i_New_numberOfLines = 0;

    int ind, indPOne, lineSize;
    int iCount;
    // Process of every agglomeration lines:
    for (int i=0; i<fineAgglomerationLines_array_Idx_size -1; i++){

        ind = fineAgglomerationLines_array_Idx[i];
        indPOne = fineAgglomerationLines_array_Idx[i + 1];
        lineSize = indPOne-ind;
        if (lineSize <= 1){
            continue;
        }
        iCount = 0;

        while (iCount + 2 <= lineSize){

            // Update of the _Fine_Cell_indices_To_Coarse_Cell_Indices table
            // TODO not check this in production
            assert(fineCellToCoarseCell[fineAgglomerationLines_array[ind+iCount]] == -1);  // "Already assigned fine cell"
            assert(fineCellToCoarseCell[fineAgglomerationLines_array[ind+iCount + 1]] == -1);  // "Already assigned fine cell"

            fineCellToCoarseCell[fineAgglomerationLines_array[ind+iCount]] = indCoarseCell;
            fineCellToCoarseCell[fineAgglomerationLines_array[ind+iCount + 1]] = indCoarseCell;
            numberOfFineAgglomeratedCells += 2;
            isFineCellAgglomerated[fineAgglomerationLines_array[ind+iCount]] = 1;
            isFineCellAgglomerated[fineAgglomerationLines_array[ind+iCount + 1]] = 1;

            new_Lines_array[iCountCoarseCell] = indCoarseCell;

            setOfAnisotropicCompliantCoarseCells.insert(indCoarseCell);

            iCount += 2;
            indCoarseCell += 1;
            iCountCoarseCell += 1;

            // if iCount < len(line): il reste une cellule toute seule!
            // i.e. la ligne etait de taile impaire.
            // 2 situations: si la cellule abandonnee est contigue a la zone Euler, c'est OK.
            //               sinon, c'est la merde!

        }  //End: while (iCount + 2 <= lineSize)

        if (iCount < lineSize){
            // Problematic Check! it may happen that the line touches the farfield
            // This is correct!
            // for example: RAE 2D case, Boxes iso_and_aniso
            // # check
            // ind = matrixAdj_CRS_row_ptr[cell]
            // ind_p_one = matrixAdj_CRS_row_ptr[cell + 1]
            // isOnBoundary = False
            // for i in xrange(ind, ind_p_one):
            //     indNeighbourCell = matrixAdj_CRS_col_ind[i]
            //     if indNeighbourCell == cell:
            //         isOnBoundary = True
            // assert not isOnBoundary, "The left alone anisotropic cell to agglomerate is on boundary"
            // # End Check

            numberOfFineAgglomeratedCells += 1;
            isFineCellAgglomerated[fineAgglomerationLines_array[ind+iCount]] = 1;  // the extra cell is agglomerated and
            fineCellToCoarseCell[fineAgglomerationLines_array[ind+iCount]] = indCoarseCell - 1;
        }

        i_New_numberOfLines += 1;
        new_Lines_index[i_New_numberOfLines] = iCountCoarseCell;

    }  // end: for (int i=0; i<fineAgglomerationLines_array_Idx_size -1; i++){
    cout<<"End of loop"<<endl;

    // update of arrays fineAgglomerationLines_array_Idx, fineAgglomerationLines_array for the new coarse level!
    for(int i=0; i<i_New_numberOfLines; i++){
        //cout<<"i= "<< i;
        fineAgglomerationLines_array_Idx[i] = new_Lines_index[i];
        ind = new_Lines_index[i];
        indPOne = new_Lines_index[i + 1];
        //cout<<" ind "<<ind<< " indPOne "<<indPOne<<endl;
        for(int j=ind; j<indPOne; j++){

            fineAgglomerationLines_array[j] = new_Lines_array[j];
        }
    }

    fineAgglomerationLines_array_Idx[i_New_numberOfLines] = new_Lines_index[i_New_numberOfLines];
    cout<<"End of update"<<endl;

    fineAgglomerationLines_array_Idx_size = i_New_numberOfLines+1;
    // fineAgglomerationLines_array_Idx.resize((i_New_numberOfLines+1,), refcheck=False)
    // fineAgglomerationLines_array.resize((fineAgglomerationLines_array_Idx[i_New_numberOfLines],),refcheck=False)

    // array_AnisotropicCompliantCoarseCells = None
    if (!setOfAnisotropicCompliantCoarseCells.empty()){
        // array_AnisotropicCompliantCoarseCells = np.zeros((len(setOfAnisotropicCompliantCoarseCells),), dtype=int);
        set<int>::const_iterator sit  = setOfAnisotropicCompliantCoarseCells.begin ();
        set<int>::const_iterator send = setOfAnisotropicCompliantCoarseCells.end ();
        int i=0;
        for(; sit != send; sit++){
            AnisotropicCompliantCoarseCells_array[i] = * sit ;
            ++i;
        }
        AnisotropicCompliantCoarseCells_size = setOfAnisotropicCompliantCoarseCells.size();

    }
    cout<<"numberOfFineAgglomeratedCells "<<numberOfFineAgglomeratedCells<<endl;
    cout<<"indCoarseCell "<<indCoarseCell<<endl;

    //return numberOfFineAgglomeratedCells, indCoarseCell, array_AnisotropicCompliantCoarseCells
    delete[] new_Lines_index;
    delete[] new_Lines_array;

    sizes[0] = fineAgglomerationLines_array_Idx_size;
    sizes[1] = numberOfFineCells;
    sizes[2] = numberOfFineAgglomeratedCells;
    sizes[3] = indCoarseCell;
    sizes[4] = AnisotropicCompliantCoarseCells_size;
}
