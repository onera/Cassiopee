/*    
    Copyright 2013-2025 Onera.

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

#include "SingleSegment.h"
using namespace K_FLD;

//=============================================================================
SingleSegment::SingleSegment(FldArrayI& indir, StructBlock* blk)
{
  E_Int size = indir.getSize();
  FldArrayF& field = blk->getGlobalField();
  E_Int nfld = field.getNfld();
  E_Int ind;

  _field.malloc(size,nfld);
  
  for (E_Int i = 0; i < size; i++)
  {
    ind = indir[i];
    for (E_Int j = 1; j <= nfld; j++)
      _field(i,j) = field(ind,j);
  }
}
//=============================================================================
SingleSegment::SingleSegment(FldArrayF& field)
{
  _field = field;
}
//=============================================================================
SingleSegment::~SingleSegment()
{
}
//=============================================================================
/* Return the indirection array  */
//=============================================================================
FldArrayF& SingleSegment::getGlobalField()
{
  return _field;
}

