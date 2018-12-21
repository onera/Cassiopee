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
#include <string.h>
#include "Array/Array.h"
#include "String/kstring.h"

using namespace std;
using namespace K_FLD;

//============================================================================
// Compare deux varstrings (du style "a,b,c")
// Retourne -1 si les deux varStrings sont differentes
// sinon retourne 0
// Les variables doivent etre positionnees au meme endroit.
// Les noms x (yz) et CoordinateX (YZ) sont equivalents.
//============================================================================
E_Int K_ARRAY::compareVarStrings(char* varString1, char* varString2)
{
  vector<char*> vars1; vector<char*> vars2;
  extractVars(varString1, vars1);
  extractVars(varString2, vars2);
  
  E_Int size1 = vars1.size();
  E_Int size2 = vars2.size();
  
  if (size1 != size2)
  {
    for (E_Int v = 0; v < size1; v++) delete [] vars1[v];
    for (E_Int v = 0; v < size2; v++) delete [] vars2[v];
    return -1;
  }
  for (E_Int v1 = 0; v1 < size1; v1++)
  {
    if (K_STRING::cmp(vars1[v1], vars2[v1]) != 0)
    {
      // Il existe des variables identiques mais de nom different?
      if (isCoordinateXPresent(vars1[v1]) == 0 &&
          isCoordinateXPresent(vars2[v1]) == 0) goto skip;
      if (isCoordinateYPresent(vars1[v1]) == 0 &&
          isCoordinateYPresent(vars2[v1]) == 0) goto skip;
      if (isCoordinateZPresent(vars1[v1]) == 0 &&
          isCoordinateZPresent(vars2[v1]) == 0) goto skip;

      // Delete vars1, vars2
      for (E_Int v = 0; v < size1; v++) delete [] vars1[v];
      for (E_Int v = 0; v < size2; v++) delete [] vars2[v];
      return -1;
    }
    skip: ;
  }

  // Delete vars1, vars2
  for (E_Int v = 0; v < size1; v++) delete [] vars1[v];
  for (E_Int v = 0; v < size2; v++) delete [] vars2[v];
  return 0;
}
