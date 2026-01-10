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
#include <string.h>
#include "Array/Array.h"
#include "String/kstring.h"

using namespace std;

//=============================================================================
/* A partir de 2 varStrings, recherche si des variables communes existent et 
   retourne les vecteurs des positions des variables correspondantes dans les
   arrays (1er element: i = 1).
   Les noms des variables communes sont dans varString.
   IN: varString1, varString2: les deux varStrings
   OUT: pos1: position des variables communes dans varString1
   OUT: pos2: position des variables communes dans varString2
   OUT: varString: nom des variables communes style ("a,b,c")
   Retourne:
   0 si les chaines varString1 et varString2 sont differentes
   1 si les deux chaines sont identiques
 */
//=============================================================================
E_Int K_ARRAY::getPosition(char* varString1, char* varString2,
                           vector<E_Int>& pos1, vector<E_Int>& pos2,
                           char* varString)
{
  E_Int equal = 1; // vaut 0 s'il existe une variable differente au moins

  vector<char*> vars1;
  extractVars(varString1, vars1);
  E_Int c1, c2;
  E_Int isfound = 0;

  // Recherche de variables communes
  E_Int vars1Size = vars1.size();
  for (c1 = 0; c1 < vars1Size; c1++)
  {
    if (K_STRING::cmp(vars1[c1], ",") != 0)
    {
      c2 = isNamePresent(vars1[c1], varString2);
      if (c2 == -1) equal = 0; // variable non presente
      else
      {
        if (isfound == 0)
        {
          strcpy(varString, vars1[c1]);
          isfound = 1;
        }
        else 
        {
          strcat(varString, ",");
          strcat(varString, vars1[c1]);
        }
        pos1.push_back(c1+1); pos2.push_back(c2+1);
      }
    } 
  }
    
  // Check if all the variables in string1 are present in the final string
  E_Int pos1Size = pos1.size();
  E_Int pos2Size = pos2.size();
  if (pos1Size != vars1Size) equal = 0; 

  // Check if the common list of variables is not empty
  if (pos1Size == 0 || pos2Size == 0) equal = -1;
      
  for (c1 = 0; c1 < vars1Size; c1++) delete [] vars1[c1];

  return equal;
}
