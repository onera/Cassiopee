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

// Formated selig (1D airfoils) file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* seligread */
// Lit le fichier selig comme une seule zone STRUCT 1DE.
//=============================================================================
E_Int K_IO::GenIO::seligread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t;
  E_Int nv, l;
  FldArrayF* f;

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");
  
  if (ptrFile == NULL)
  {
    printf("Warning: seligread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture premiere ligne -> nom +
  char* zoneName = new char [128];
  res = readline(ptrFile, zoneName, 128);
  
  if (res == -1)
  {
    printf("Warning: seligread: cannot read first line.\n");
    return 1;
  }

  // Lecture de tous les vertex pour dimensionner
  nv = 0; res = 1;
  while (res == 1 || res == 0)
  {
    res = readDouble(ptrFile, t, -1);
    res = readDouble(ptrFile, t, -1);
    nv += 1;
  }
  nv = nv-1;
  
  //printf("I found " SF_D_ " vertices\n", nv); fflush(stdout);
  if (nv == 0) { fclose(ptrFile); return 1; } // FAILED

  f = new FldArrayF(nv, 3);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);
  fseek(ptrFile, 0, SEEK_SET);
  readline(ptrFile, zoneName, 1023);
  l = strlen(zoneName);
  // rip blank et cr
  E_Int i;
  for (i = l-1; i > 0; i--)
  {
    if (zoneName[i] != '\n' && zoneName[i] != '\r' && zoneName[i] != ' ') break;
  }
  zoneName[i+1] = '\0';
  //printf(SF_D_ " %s\n", i, zoneName);
  
  for (E_Int i = 0; i < nv; i++)
  {
    res = readDouble(ptrFile, t, -1); fx[i] = t; //printf(SF_F_ " ", t);
    res = readDouble(ptrFile, t, -1); fy[i] = t; //printf(SF_F_ " ", t);
    fz[i] = 0.; //printf(SF_F_ "\n", t);
  }

  structField.push_back(f);
  ni.push_back(nv); nj.push_back(1); nk.push_back(1);

  // Cree le nom de zone
  zoneNames.push_back(zoneName);

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  
  return 0;
}

//=============================================================================
// Ecrit seulement une seule zone STRUCT 1D.
//=============================================================================
E_Int K_IO::GenIO::seligwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = structField.size();
  E_Int nvalidZones = 0;
  E_Int no = -1;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (nj[zone] == 1 && nj[zone] == 1)
    { nvalidZones++; if (no == -1) no = zone; }
    else
      printf("Warning: seligwrite: zone " SF_D_ " not written (not a 1D STRUCT).", zone);
  }

  if (nvalidZones == 0) return 1;
  if (nvalidZones > 1) printf("Warning: seligwrite: only first zone will be written.");
  
  // Zone must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: seligwrite: zone do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  // Open file
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: seligwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Build writing data format
  char format1[134]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  E_Int l = strlen(dataFmt);
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format of data
  sprintf(format1," %s %s\n", dataFmt, dataFmtl);

  // Header
  fprintf(ptrFile, "%s\n", zoneNames[no]);
  
  FldArrayF* a = structField[no];
  E_Int niz = ni[no];
  E_Float* fx = a->begin(posx);
  E_Float* fy = a->begin(posy);
    
  for (E_Int i = 0; i < niz; i++)
  {
    fprintf(ptrFile, format1, fx[i], fy[i]);
  }

  fclose(ptrFile);
  return 0;
}
