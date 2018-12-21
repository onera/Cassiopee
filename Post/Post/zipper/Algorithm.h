/*    
    Copyright 2013-2019 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _POST_ZIPPER_ALGORITHM_H
#define _POST_ZIPPER_ALGORITHM_H

# include "kcore.h"
# include <vector>
# include "StructBlock.h"
# include "CString.h"
# include "SegmentPair.h"
# include "TriangleZ.h"
# include "SingleSegment.h"
# include "Pocket.h"
# include "ZipLib.h"

# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI
# define FldArrayB K_FLD::FldArrayB

void computeIBlank(std::vector<StructBlock*>& structBlocks);
void computeStringing(std::vector<StructBlock*>& structBlocks,
                      std::vector<CString*>& strings);
void computeMatchingSegments(std::vector<CString*>& strings,
                             std::vector<SegmentPair*>& segPairs);
void closePockets(std::vector<CString*>& strings, 
                  std::vector<SegmentPair*>& segPairs,
                  std::vector<FldArrayF*>& field,
                  std::vector<FldArrayI*>& triConnect);
void mergeAllZonesInOne(std::vector<StructBlock*>& structBlocks,
                        std::vector<FldArrayF*>& field,
                        std::vector<FldArrayI*>& idgT, 
                        FldArrayI& idgG, FldArrayF& fieldG,
                        FldArrayI& FirstPt);
void computeConnectivity(std::vector<StructBlock*>& structBlocks,
                         std::vector<FldArrayI*>& triConnect,
                         std::vector<FldArrayI*>& idgT,
                         FldArrayI& idgG, 
                         FldArrayI& unsConnectENG,
                         FldArrayI& FirstPt);
void deletingUnusedNodes(FldArrayI& unsConnectENG, 
                         FldArrayF& fieldG);
// void computeIntegrals(FldArrayI& unsConnectENG,
//                       FldArrayF& fieldG);x
void zipInbetweenPairs(std::vector<SegmentPair*>& segPairs,
                       std::vector<FldArrayF*>& field,
                       std::vector<FldArrayI*>& triConnect,
                       FldArrayB& isZipped);
// void writeUnstructMesh(char* outfile, char* varString, 
//                        FldArrayF& fieldG, 
//                        FldArrayI& unsConnectENG);
// void zipper(char* outfile, char* varString,
//             vector<StructBlock*>& structBlocks);
# undef FldArrayF
# undef FldArrayI
# undef FldArrayB

#endif
