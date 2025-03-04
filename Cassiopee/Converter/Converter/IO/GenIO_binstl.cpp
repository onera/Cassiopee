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

// Binary STL (StereoLithography) file support

# include "GenIO.h"
# include "Array/Array.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   stlread 
   Un seul domaine peut-etre stocke dans ce type de fichier.
*/
//=============================================================================
E_Int K_IO::GenIO::stlread( 
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: stlread: cannot open file %s.\n", file);
    return 1;
  }
  
  // Header
  char dummy[81]; E_Int ret;
  ret = fread(dummy, sizeof(char), 80, ptrFile);
  if (ret < 80) { fclose(ptrFile); return 1; }

  // Nombre de noeuds
  int nd = 0;
  ret = fread(&nd, sizeof(int), 1, ptrFile);
  if (ret != 1) { fclose(ptrFile); return 1; }
  if (nd < 0) 
  { fclose(ptrFile); return 1;} // cette valeur est arbitraire
  

  // Champ des vertex
  float buf[12];
  short count;
  FldArrayF* f = new FldArrayF(3*nd, 3);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);
  unstructField.push_back(f);

  E_Int sizeoffloat = sizeof(float);
  E_Int sizeofshort = sizeof(short);
  for (E_Int i = 0; i < nd; i++)
  {
    fread(buf, sizeoffloat, 12, ptrFile);
    fread(&count, sizeofshort, 1, ptrFile);
    fx[3*i] = buf[3];
    fy[3*i] = buf[4];
    fz[3*i] = buf[5];
    fx[3*i+1] = buf[6];
    fy[3*i+1] = buf[7];
    fz[3*i+1] = buf[8];
    fx[3*i+2] = buf[9];
    fy[3*i+2] = buf[10];
    fz[3*i+2] = buf[11];
  }

  // Connectivite
  eltType.push_back(2);
  FldArrayI* cp = new FldArrayI(nd, 3);
  connect.push_back(cp);
  E_Int* c1 = cp->begin(1);
  E_Int* c2 = cp->begin(2);
  E_Int* c3 = cp->begin(3);
  for (E_Int i = 0; i < nd; i++)
  {
    c1[i] = 3*i+1;
    c2[i] = 3*i+2;
    c3[i] = 3*i+3;
  }

  K_CONNECT::cleanConnectivity(1, 2, 3, 
                               1.e-14,  "TRI",
                               *unstructField[0], *connect[0]);

  // Cree les noms de zones
  char* zoneName = new char [128];
  sprintf(zoneName, "Zone0");
  zoneNames.push_back(zoneName);

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Only write ONE triangle array. Others are discarded.
//=============================================================================
E_Int K_IO::GenIO::stlwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  E_Int no = -1;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltTypes[zone][0] == 2) // triangles
    { nvalidZones++; if (no == -1) no = zone; }
    else
      printf("Warning: stlwrite: zone " SF_D_ " not written (not a triangle zone).", zone);
  } 

  if (nvalidZones == 0) return 1;
  if (nvalidZones > 1) printf("Warning: stlwrite: only first zone will be written.");
  
  // Zone must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: stlwrite: zone do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  // Open file
  FILE* ptrFile = fopen(file, "wb");
  if (ptrFile == NULL) 
  {
    printf("Warning: stlwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Header
  char dummy[80];
  for (E_Int i = 0; i < 80; i++) dummy[i] = ' ';
  strcpy(dummy, "Written by Cassiopee");
  fwrite(dummy, sizeof(char), 80, ptrFile);
  
  // Nombre de facettes
  FldArrayF* a = unstructField[no];
  FldArrayI& c = *connect[no];
  int nd = c.getSize();
  E_Float* fx = a->begin(posx);
  E_Float* fy = a->begin(posy);
  E_Float* fz = a->begin(posz);
  float buf[12];
  E_Int p;
  short count = 0;
  fwrite(&nd, sizeof(int), 1, ptrFile);
  
  E_Int sizeoffloat = sizeof(float);
  E_Int sizeofshort = sizeof(short);
  
  for (E_Int i = 0; i < nd; i++)
  {
    buf[0] = 0.; // normales a detecter
    buf[1] = 0.;
    buf[2] = 0.;
    p = c(i,1)-1;
    buf[3] = fx[p];
    buf[4] = fy[p];
    buf[5] = fz[p];
    p = c(i,2)-1;
    buf[6] = fx[p];
    buf[7] = fy[p];
    buf[8] = fz[p];
    p = c(i,3)-1;
    buf[9] = fx[p];
    buf[10] = fy[p];
    buf[11] = fz[p];

    fwrite(buf, sizeoffloat, 12, ptrFile);
    fwrite(&count, sizeofshort, 1, ptrFile);
  }
  fclose(ptrFile);
  return 0;
}
