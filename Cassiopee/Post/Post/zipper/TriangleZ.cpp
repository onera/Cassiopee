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
#include "TriangleZ.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
TriangleZ::TriangleZ(FldArrayF& field1, FldArrayF& field2, FldArrayF& field3)
{  
  _fieldA = field1;
  _fieldB = field2;
  _fieldC = field3;
}

//=============================================================================
TriangleZ::~TriangleZ()
{
}
//=============================================================================
/* Return the coordinates of the triangle vertex of ind1 */
//=============================================================================
FldArrayF& TriangleZ::getField1()
{
  return _fieldA;
}
//=============================================================================
/* Return the coordinates of the triangle vertex of ind2 */
//=============================================================================
FldArrayF& TriangleZ::getField2()
{
  return _fieldB;
}
//=============================================================================
/* Return the coordinates of the triangle vertex of ind3 */
//=============================================================================
FldArrayF& TriangleZ::getField3()
{
  return _fieldC;
}

