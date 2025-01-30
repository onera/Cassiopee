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
#include "Array.h"

//=============================================================================
void K_ARRAY::cleanStructFields(std::vector<K_FLD::FldArrayF*>& structF)
{
  E_Int structFSize = structF.size();
  for (E_Int v = 0; v < structFSize; v++) delete structF[v];
  structF.clear();
}
//============================================================================
void K_ARRAY::cleanUnstrFields(std::vector<K_FLD::FldArrayF*>& unstrF, 
                               std::vector<K_FLD::FldArrayI*> cnt,
                               std::vector<char*>& eltType)
{
  E_Int unstrFSize = unstrF.size();
  E_Int cntSize = cnt.size();
  for (E_Int v = 0; v < unstrFSize; v++) delete unstrF[v];
  for (E_Int v = 0; v < cntSize; v++) delete cnt[v];

  unstrF.clear(); cnt.clear(); eltType.clear();
}
