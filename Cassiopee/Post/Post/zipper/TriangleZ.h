/*    
    Copyright 2013-2026 ONERA.

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
#ifndef _POST_ZIPPER_TRIANGLEZ_H_
#define _POST_ZIPPER_TRIANGLEZ_H_

#include "StructBlock.h"
# define FldArrayF K_FLD::FldArrayF
//=============================================================================
/* Class defining a triangle  */
//=============================================================================
class TriangleZ
{
  public:
    
    ///+ 1- Constructors / Destructor
    /** Constructor :
     IN: field1: coord + solution of point 1 
     IN: field2: for point 2
     IN: field3: for point 3 */
    TriangleZ(FldArrayF& field1, FldArrayF& field2, 
              FldArrayF& field3);
    
    /** Destructor */
    ~TriangleZ();

    /** Return the coordinates and solution of the triangle vertex of 
        index 1 */
    FldArrayF& getField1();

    /** Return the coordinates and solution of the triangle vertex of 
        index 2 */
    FldArrayF& getField2();

    /** Return the coordinates and solution of the triangle vertex of 
        index 3 */
    FldArrayF& getField3();
    
    
  private:
    FldArrayF _fieldA; // coord + solution of pt 1
    FldArrayF _fieldB;
    FldArrayF _fieldC;
# undef FldArrayF
};

#endif
