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
# ifndef _POST_ZIPPER_POCKET_H_
# define _POST_ZIPPER_POCKET_H_

# include "kcore.h"
# include "CString.h"
# include "StructBlock.h"
# include "TriangleZ.h"
# include<vector>
# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI

//=============================================================================
/* Class defining the list of segments set in a pocket   */
//=============================================================================
class Pocket
{ 
  public: 
    ///+ 1- Constructors / Destructor
    
    /** Constructor */
    Pocket(FldArrayF& field);

    /** Destructor */
    ~Pocket();
    
    /** Write the pocket contours and store it in file fileName*/
    void writeLine( char* fileName, E_Boolean add);  
    
    /** Close pocket by building triangles. 
        Return true if triangles are created*/
    E_Boolean closePocket( FldArrayF& field, FldArrayI& connect);

  private :
    /* Close pocket when a set of points are aligned*/
    void closeSpecificPocket(FldArrayF& field1, 
                             FldArrayF& field2,
                             std::vector<TriangleZ*>& triangles);
    /* Check if a set of points are aligned. If true, then return field1 
       the list of aligned points, and field2 opposite points */
    E_Boolean checkIfPtsAreAligned(FldArrayF& field1, 
                                   FldArrayF& field2);

    /* Computes triangles  starting from istart */
    E_Boolean computeTriangulation(std::vector<TriangleZ*>& triangles);
    
    /* Compute the connectivity of triangles */
    void compConnectivity(std::vector<TriangleZ*>& triangles,
                          FldArrayF& field,
                          FldArrayI& connect);
    
    /* Compute triangulation for the quad PiPi+1PjPj+1 of the pocket */
    E_Boolean compDelaunay(E_Int& iA, E_Int& iB,
                           E_Int& iC, E_Int& iD,
                           FldArrayF& field1, FldArrayF& field2,
                           std::vector<TriangleZ*>& triangles);

   /* Compute the  distance between ind1 and ind2 */
    E_Float compDistanceBetweenPoints(FldArrayF& field1, 
                                      FldArrayF& field2);

    /* Test the validity of the triangulation done with respect to the given
       diagonal:
       diagType = 1 : (AC)
       diagType = 2 : (BD)
       Return true if the normals of triangles are in the same direction */
    E_Boolean checkTriangles( E_Int diagType,
                              FldArrayF& fieldA, 
                              FldArrayF& fieldB,
                              FldArrayF& fieldC, 
                              FldArrayF& fieldD);
 private:
    E_Int _size;
    E_Int _nfld;
    FldArrayF _field;
};
# undef FldArrayF
# undef FldArrayI

#endif
