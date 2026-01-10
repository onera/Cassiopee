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
# include "Array/Array.h"
# include <string.h>

using namespace std;

//=============================================================================
/* Ajoute un prefixe a tous les noms de variables contenus dans varString. 
   IN: varString: la chaine des variables style "a,b,c"
   IN: prefix: la chaine de prefix a ajouter a chaque nom de variable
   OUT: varStringOut: la chaine resultat (doit etre deja allouee). */
//=============================================================================
void K_ARRAY::addPrefixInVarString(char* varString, const char* prefix,
                                   char* varStringOut)
{
  vector<char*> vars;
  extractVars(varString, vars);
  E_Int size = vars.size();

  strcpy(varStringOut, "");
  for (E_Int v = 0; v < size; v++)
  {
    strcat(varStringOut, prefix);
    strcat(varStringOut, vars[v]);
    if (v < size-1) strcat(varStringOut, ",");
  } 

  for (E_Int v = 0; v < size; v++) delete [] vars[v];
}
