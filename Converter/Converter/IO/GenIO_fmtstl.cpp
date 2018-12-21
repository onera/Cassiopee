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

// Formated Stl file support

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
/* fstlread */
// Lit le fichier comme une seule zone TRI.
//=============================================================================
E_Int K_IO::GenIO::fstlread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t;
  E_Int i, nv;
  FldArrayF* f;

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: fstlread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de tous les vertex
  nv = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "VERTEX "); nv++;
  }
  nv = nv-1;
  //printf("I found %d vertices\n", nv);
  if (nv == 0) { fclose(ptrFile); return 1; } // FAILED

  f = new FldArrayF(nv, 3);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);
  fseek(ptrFile, 0, SEEK_SET);
  i = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "VERTEX ");
    if (res == 1)
    {
      res = readDouble(ptrFile, t, -1); fx[i] = t; //printf("%f ", t);
      res = readDouble(ptrFile, t, -1); fy[i] = t; //printf("%f ", t);
      res = readDouble(ptrFile, t, -1); fz[i] = t; i++; //printf("%f\n", t);
      if (res == 0) res = 1; else res = 0;
    }
  }

  // facettes
  E_Int ne = nv/3;
  FldArrayI* cn = new FldArrayI(ne, 3);
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);

  for (i = 0; i < ne; i++)
  {
    cn1[i] = 3*i+1; cn2[i] = 3*i+2; cn3[i] = 3*i+3;
  }

  unstructField.push_back(f);
  eltType.push_back(2);
  connect.push_back(cn);

  K_CONNECT::cleanConnectivity(1, 2, 3, 
                               1.e-14,  "TRI",
                               *unstructField[0], *connect[0]);

  // Cree le nom de zone
  char* zoneName = new char [128];
  sprintf(zoneName, "Zone0");
  zoneNames.push_back(zoneName);

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  
  return 0;
}

//=============================================================================
// Ecrit seulement une seule zone TRI.
//=============================================================================
E_Int K_IO::GenIO::fstlwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  E_Int no = -1;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltType[zone] == 2) // triangles
    { nvalidZones++; if (no == -1) no = zone; }
    else
      printf("Warning: fstlwrite: zone %d not written (not a triangle zone).", zone);
  } 

  if (nvalidZones == 0) return 1;
  if (nvalidZones > 1) printf("Warning: fstlwrite: only first zone will be written.");
  
  // Zone must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: fstlwrite: zone do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  // Open file
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: fstlwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Build writing data format
  char format1[134]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt);
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format of data
  sprintf(format1,"      vertex %s%s%s\n", dataFmt, dataFmt, dataFmtl);

  // Header
  fprintf(ptrFile, "solid Cassiopee STL\n");
  
  FldArrayF* a = unstructField[no];
  FldArrayI& c = *connect[no];
  int ne = c.getSize();
  E_Float* fx = a->begin(posx);
  E_Float* fy = a->begin(posy);
  E_Float* fz = a->begin(posz);
  
  E_Int* c1 = c.begin(1);
  E_Int* c2 = c.begin(2);
  E_Int* c3 = c.begin(3);
  E_Int p;
  for (E_Int i = 0; i < ne; i++)
  {
    fprintf(ptrFile, "  facet normal 0. 0. 0.\n");
    fprintf(ptrFile, "    outer loop\n");

    p = c1[i]-1;
    fprintf(ptrFile, format1, fx[p], fy[p], fz[p]);
    p = c2[i]-1;
    fprintf(ptrFile, format1, fx[p], fy[p], fz[p]);
    p = c3[i]-1;
    fprintf(ptrFile, format1, fx[p], fy[p], fz[p]);
    fprintf(ptrFile, "    endloop\n");
    fprintf(ptrFile, "  endfacet\n");
  }

  fprintf(ptrFile, "endsolid Cassiopee STL\n");
  fclose(ptrFile);
  return 0;
}
