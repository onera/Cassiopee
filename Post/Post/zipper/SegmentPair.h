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
#ifndef _POST_ZIPPER_SEGMENT_PAIR_H_
#define _POST_ZIPPER_SEGMENT_PAIR_H_

# include "kcore.h"
# include "StructBlock.h"
# include "CString.h"
# include "TriangleZ.h"
# include "SingleSegment.h"
# include "ZipLib.h"
# include <vector>
# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI

//=============================================================================
/* Class defining a segment pair  */
//=============================================================================
class SegmentPair
{
  public: 
    ///+ 1- Constructors / Destructor
    
    /** Constructor */
    SegmentPair(CString* st1, CString* st2);
    
    /** Destructor */
    ~SegmentPair();

    /** Return the array for first string in pair */
    FldArrayI& getIndArray1();
    /** Return the array for second string in pair */
    FldArrayI& getIndArray2();
    /** Return block of first string */
    StructBlock* getBlock1();
    /** Return block of second string */
    StructBlock* getBlock2();

    /** Computes zipper for the pair of segment.
        out : field : coordinates of triangles 
        out : connect : connectivity
        E_Boolean = true if the zipping is done */
    E_Boolean computeZipper(const E_Int nfieldTot,
                            FldArrayF& field,
                            FldArrayI& connect);

    /** Modify the flag array: 
       flag = 0: point not set in a segment pair
       flag = 1: point = extremity of a segment pair element
       flag = 2: point = internal pt of a segment pair element 
    */
    void updateFlagForStrings(std::vector<SegmentPair*>& segPairs);    
    
    /** Add extremities of segment pairs to the list of single segments to 
        be treated */
    void identifyLastEdges(std::vector<FldArrayF*>& lastEdges);

  private:

    /* Computes the  distance between ind1 of seg1 and ind2 of seg2 */
    E_Float compSegmentPtsDistance(E_Int ind1, E_Int ind2);

    /* Given point A of index iprev and point B of index inext of seg1
             point D of index jnext and point C of index jprev of seg2
             build the triangulation if it is possible. Else exit */
    void compDelaunay(const E_Int nfieldTot,
                      E_Int& iprev, E_Int& inext,
                      E_Int& jprev, E_Int& jnext,
                      std::vector<TriangleZ*>& triangles);

    /* Test the validity of the triangulation. Returns :
       0 if both diagonals of ABCD are valid
       1 if only AC is valid 
       2 if only Bd is valid*/ 
    E_Int checkTriangles(FldArrayF& fieldA, FldArrayF& fieldB,
                         FldArrayF& fieldC, FldArrayF& fieldD);

    /* Connect the remaining points of the bigger segment to the extremum of
       the other segment */
    void connectRemainingPtsToEnd(const E_Int nfieldTot,
                                  E_Int smallestseg, E_Int iend,
                                  E_Int iflag, E_Int imax, 
                                  std::vector<TriangleZ*>& triangles);
    
    /* Computes connectivity for list of triangles */
    void compConnectivity(const E_Int nfieldTot,
                          std::vector<TriangleZ*>& triangles,
                          FldArrayF& field,
                          FldArrayI& connect);
    
//     /* Reverse segment pairs ind arrays */
//     void reverseSegmentIndArray();

    /* Search for corresponding segment and return :
       True:  if extremities of the current segment matches with extremities of
       another segment pair
       False: if extremities of the current segment matches with extremities of
       segments that come from 2 different pairs or if not found*/
    E_Boolean searchForCorrespondingSegments(
      E_Int ind1, E_Int ind2, 
      std::vector<SegmentPair*>& segPairs);
  private:
    // Indices of points of the 1st matching segment 
    FldArrayI _indseg1;
    // blk of the 1st matching segment
    StructBlock* _blk1;
    // Original string from which derives the segment s1
    CString* _origSt1;

    // Indices of points of the 2nd matching segment 
    FldArrayI _indseg2;
    // blk of the 2nd matching segment
    StructBlock* _blk2;
    // Original string from which derives the segment s2
    CString* _origSt2;
 };
# undef FldArrayF
# undef FldArrayI

#endif
