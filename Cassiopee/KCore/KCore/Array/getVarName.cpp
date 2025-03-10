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

# include "Array/Array.h"
# include <string.h>
# include <stdio.h>

//=============================================================================
/* Construit la varString a partir d'un object python liste de strings.
   Retourne le nombre de variables.
   IN: varNames: objet liste python ['a','b','c']
   OUT: varString: "a,b,c"
*/
//=============================================================================
E_Int K_ARRAY::getVarName(PyObject* varNames, char* varString)
{
  PyObject* tpl;
  char* s;
 
  E_Int nvars = -1;
  if (PyList_Check(varNames) == 0)
  {
    printf("Error: getVarName: unrecognised variables list.\n"); 
    printf(" Set to default: x,y,z,ro,rou,rov,row,roE.\n");
    strcpy(varString, "x,y,z,ro,rou,rov,row,roE");
    return 5;
  }
  
  nvars = PyList_Size(varNames);
  strcpy(varString, "x,y,z");

  for (int i = 0; i < nvars; i++)
  {
    tpl = PyList_GetItem(varNames, i);
    if (PyString_Check(tpl))
    {
      s = PyString_AsString(tpl);
      strcat(varString, ",");
      strcat(varString, s); 
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl))
    {
      s = (char*)PyUnicode_AsUTF8(tpl); 
      strcat(varString, ",");
      strcat(varString, s);
    }
#endif
    else
    {
      printf("Error: getVarName: variables must be strings.\n");
      return -1;
    }
    
  }
  return nvars;
}
