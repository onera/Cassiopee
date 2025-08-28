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
#include "Data.h"

// Construit les chaines de toutes les variables contenues dans structVarString
// et unstrVarString
// il faut liberer la memoire de referenceVarNames
void Data::getAllVars(std::vector<char*>& structVarString,
                      std::vector<char*>& unstrVarString,
                      E_Int& referenceNfield,
                      char**& referenceVarNames)
{
  std::vector<char*> allvars;
  for (size_t i = 0; i < structVarString.size(); i++)
  {
    char* v = structVarString[i];
    std::vector<char*> vars;
    K_ARRAY::extractVars(v, vars);
    for (size_t j = 0; j < vars.size(); j++)
    {
      E_Bool found = false;
      char* v1 = vars[j];
      for (size_t k = 0; k < allvars.size(); k++)
      {
        char* v2 = allvars[k];
        if (strcmp(v1, v2) == 0) { found = true; break; }
      }
      if (found == false) 
      {
        if (strcmp(v1, "x") != 0 && strcmp(v1, "y") != 0 && strcmp(v1, "z") != 0 &&
            strcmp(v1, "CoordinateX") != 0 && strcmp(v1, "CoordinateY") != 0 && strcmp(v1, "CoordinateZ") != 0)
        {
          char* s = new char [strlen(v1)+1];
          strcpy(s, v1);
          allvars.push_back(s);
        }
      }
    }
    for (size_t j = 0; j < vars.size(); j++) delete [] vars[j];
  }
  for (size_t i = 0; i < unstrVarString.size(); i++)
  {
    char* v = unstrVarString[i];
    std::vector<char*> vars;
    K_ARRAY::extractVars(v, vars);
    for (size_t j = 0; j < vars.size(); j++)
    {
      E_Bool found = false;
      char* v1 = vars[j];
      for (size_t k = 0; k < allvars.size(); k++)
      {
        char* v2 = allvars[k];
        if (K_STRING::cmp(v1, v2) == 0) { found = true; break; }
      }
      if (found == false) 
      {
        if (K_STRING::cmp(v1, "x") != 0 && K_STRING::cmp(v1, "y") != 0 && K_STRING::cmp(v1, "z") != 0 &&
            K_STRING::cmp(v1, "CoordinateX") != 0 && K_STRING::cmp(v1, "CoordinateY") != 0 && K_STRING::cmp(v1, "CoordinateZ") != 0)
        {
          char* s = new char [strlen(v1)+1];
          strcpy(s, v1);
          allvars.push_back(s);
        }
      }
    }
    for (size_t j = 0; j < vars.size(); j++) delete [] vars[j];
  }
  referenceNfield = allvars.size();
  referenceVarNames = new char* [referenceNfield];
  for (E_Int i = 0; i < referenceNfield; i++) referenceVarNames[i] = allvars[i];

  //for (size_t i = 0; i < allvars.size(); i++) printf("%s\n", allvars[i]);
}
