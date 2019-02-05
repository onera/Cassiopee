/* 
 * File:   pdm_anisotropic_agglomerator.h
 * Author: Nicolas Lantos
 *
 * Created on November 18, 2017, 1:24 PM
 */

#ifndef UNTITLED_AGGLOMERATOR_H
#define UNTITLED_AGGLOMERATOR_H


bool computeAnisotropicLine(int* sizes,
                           int * adjMatrix_row_ptr,
                           int *adjMatrix_col_ind,
                           double * adjMatrix_areaValues,
                           int * arrayOfFineAnisotropicCompliantCells,
                           int * agglomerationLines_Idx,
                           int * agglomerationLines,
                           bool verbose);

void agglomerate_Anisotropic_One_Level_without_list_lines(int *sizes,
                                                          int *fineAgglomerationLines_array_Idx, int *fineAgglomerationLines_array,
                                                          int *fineCellToCoarseCell,
                                                          bool *isFineCellAgglomerated,
                                                          int *AnisotropicCompliantCoarseCells_array);


class Agglomerator {

};


#endif //UNTITLED_AGGLOMERATOR_H
