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
# include <string.h>
# include "Array/Array.h"

//=============================================================================
/* Retourne le nombre de variables dans varString (style "a,b,c")
   IN: varString: la chaine de variables. */
//=============================================================================
E_Int K_ARRAY::getNumberOfVariables(const char* varString)
{
  E_Int lenString = strlen(varString);

  // Compte les virgules
  E_Int cnt = 0;
  if (lenString == 0) return 0;
  for (E_Int p = 0; p < lenString; p++)
    if (varString[p] == ',') cnt++;
  return cnt+1;
}
