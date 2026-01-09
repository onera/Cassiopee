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

// Formated IGES file support (not usable, use OCC module instead)

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* igesread */
// Sont reconnus :
// - Point (116)
// - Line
// - splines
// - nurbs
// - surface splines
// - surface nurbs
// - trim surfaces 
//=============================================================================
E_Int K_IO::GenIO::igesread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  E_Float density = 1.;

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: igesread: cannot open file %s.\n", file);
    return 1;
  }

  char l[81];
  char type;

  // Boucle de lecture
  printf("\n");
  E_Int ret = 0;
  E_Int pos = 0;
  E_Int ival;
  while (ret != -1)
  {
    ret = readline(ptrFile, l, 82); pos = 0;
    if (ret == 1)
    {
      printf("Warning: igesread: strangely formatted file.\n");
    }
    //printf("%s\n", l);
    type = l[72];
    printf("line type %c\n", type);
    if (type == 'D') // Directory entry
    {
      // Type (1)
      ret = readInt(l, 73, pos, ival); 
      if (ret == 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" type=" SF_D_ "\n", ival);

      // Parameter datas (2)
      ret = readInt(l, 73, pos, ival); 
      if (ret == 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" params=" SF_D_ "\n", ival);

      // Structure (3)
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" struct=" SF_D_ "\n", ival);
      
      
      // Line font pattern (4)
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" font=" SF_D_ "\n", ival);

      // Level
      ret = readInt(l, 73, pos, ival);
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" level=" SF_D_ "\n", ival);

      // View
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" view=" SF_D_ "\n", ival);

      // Matrix
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0)
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" matrix=" SF_D_ "\n", ival);

      
      // Label display
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" label=" SF_D_ "\n", ival);
      
      // Line weight
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" line weight=" SF_D_ "\n", ival);
      
      // Color number
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" color=" SF_D_ "\n", ival);
      
      // Line count
      ret = readInt(l, 73, pos, ival);
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" count=" SF_D_ "\n", ival);
      
      // Form number
      ret = readInt(l, 73, pos, ival); 
      if (ret < 0) 
      { pos = 0; ret = readline(ptrFile, l, 82); 
        ret = readInt(l, 73, pos, ival); ret = readInt(l, 73, pos, ival); }
      printf(" form=" SF_D_ "\n", ival);
    }
    else if (type == 'P') // Parameters
    {
      ret = readInt(l, 73, pos, ival);
      switch (ival)
      {
        case 100: // circular arc
        {
          E_Float zt, x1, y1, x2, y2, x3, y3;
          ret = readDouble(l, 73, pos, zt);
          ret = readDouble(l, 73, pos, x1);
          ret = readDouble(l, 73, pos, y1);
          ret = readDouble(l, 73, pos, x2);
          ret = readDouble(l, 73, pos, y2);
          ret = readDouble(l, 73, pos, x3);
          ret = readDouble(l, 73, pos, y3);
        }

        case 110: // line
        {
          E_Float x1, y1, z1, x2, y2, z2;
          ret = readDouble(l, 73, pos, x1);
          ret = readDouble(l, 73, pos, y1);
          ret = readDouble(l, 73, pos, z1);
          ret = readDouble(l, 73, pos, x2);
          ret = readDouble(l, 73, pos, y2);
          ret = readDouble(l, 73, pos, z2);
          createLine(structField, ni, nj, nk, density, x1, y1, z1,
                     x2, y2, z2);          
          break;
        }

        case 112: // spline curve
        {
          E_Int ctype, ndim;
          E_Float h, x1, y1, z1, x2, y2, z2;
          ret = readInt(l, 73, pos, ctype);
          ret = readDouble(l, 73, pos, h);
          ret = readInt(l, 73, pos, ndim);

          ret = readDouble(l, 73, pos, x1);
          ret = readDouble(l, 73, pos, y1);
          ret = readDouble(l, 73, pos, z1);
          ret = readDouble(l, 73, pos, x2);
          ret = readDouble(l, 73, pos, y2);
          ret = readDouble(l, 73, pos, z2);
          createLine(structField, ni, nj, nk, density, x1, y1, z1,
                     x2, y2, z2);
          break;
        }


        case 116: // point
        {
          E_Float x = 0;
          E_Float y = 0;
          E_Float z = 0;
          ret = readDouble(l, 73, pos, x);
          ret = readDouble(l, 73, pos, y);
          ret = readDouble(l, 73, pos, z);
          createPoint(unstructField, connect, eltType, density, x, y, z);
          break;
        }

        default:
          printf("Warning: I got an unknown entity (" SF_D_ ")\n", ival);
      }
    }
  }

  // Cree les noms des zones structurees
  E_Int structSize = structField.size();
  for (E_Int i=0; i<structSize;i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, i);
    zoneNames.push_back(zoneName);
  }

  // Cree les noms des zones non structurees
  for (unsigned int i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, i+structSize);
    zoneNames.push_back(zoneName);
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  return 0;
}

