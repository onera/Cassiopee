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
# include <string.h>
# include "Array/Array.h"

using namespace std;

//=============================================================================
// Construit une liste des variables a partir de la chaine des variables
// C'est la responsabilite de l'appelant de liberer la memoire des vars.
// Note: supprime les blancs.
// IN: varString: de type "x,y,z,ro"
// OUT: un vector <"x", "y", "z", "ro">
//=============================================================================
void K_ARRAY::extractVars(const char* varString, vector<char*>& vars)
{
  E_Int l = strlen(varString);
  E_Int c = 0;
  
  E_Int p = 0;
  char* temp = new char [l+1];
  
  while (c < l)
  {
    if (varString[c] != ' ')
    {
      if (varString[c] == ',')
      {
        char* m = new char[VARNAMELENGTH];
        temp[p] = '\0';
        for (E_Int i = 0; i <= p; i++) m[i] = temp[i];
        vars.push_back(m);
        p = 0;
      }
      else
      {
        temp[p] = varString[c]; p++;
      }
    }
    c++;
  }

  if (p > 0)
  {
    char* m = new char[VARNAMELENGTH];
    temp[p] = '\0';
    for (E_Int i = 0; i <= p; i++) m[i] = temp[i];  
    vars.push_back(m);
  }

  delete [] temp;
}
