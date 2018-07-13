/*    
    Copyright 2013-2018 Onera.

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

// Formated Obj file support

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
/* objread */
// Manque: commentaires
// Nettoyage de f a la fin (vertex reellement utilises)
//=============================================================================
E_Int K_IO::GenIO::objread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t;
  E_Int ti, i, j, nv, nf, nt, nq, k;
  char buf[256];
  FldArrayF* f;
  FldArrayI* cnq;
  FldArrayI* cnt;

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: objread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de tous les vertex
  nv = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "V "); nv++;
  }
  nv = nv-1;
  //printf("I found %d vertices\n", nv);

  f = new FldArrayF(nv, 3);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);
  KFSEEK(ptrFile, 0, SEEK_SET);
  i = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "V ");
    if (res == 1)
    {
      res = readDouble(ptrFile, t, -1); fx[i] = t; //printf("%f ", t);
      res = readDouble(ptrFile, t, -1); fy[i] = t; //printf("%f ", t);
      res = readDouble(ptrFile, t, -1); fz[i] = t; i++; //printf("%f\n", t);
      if (res == 0) res = 1; else res = 0;
    }
  }

  // Lecture des faces
  KFSEEK(ptrFile, 0, SEEK_SET);
  nf = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "F "); nf++;
  }
  nf = nf-1;
  //printf("I found %d faces\n", nf);

  // Type des faces
  KFSEEK(ptrFile, 0, SEEK_SET);
  FldArrayI elts(nf);
  res = 1;
  res = readGivenKeyword(ptrFile, "F ");
  i = 0; nt = 0; nq = 0;
  while (res >= 1)
  {
    for (j = 0; j < 4; j++)
    {
      res = readWord(ptrFile, buf); 
    }
    if (res >= 1 && strcmp(buf, "f") == 0)
    {elts[i] = 3; nt++; }
    else if ((res >= 1 || res == 0) && (buf[0] < 48 || buf[0] > 57))
    {elts[i] = 3; nt++; res = readGivenKeyword(ptrFile, "F "); }
    else if (res == 0 && buf[0] >= 48 && buf[0] <= 57)
    {elts[i] = 4; nq++; res = readGivenKeyword(ptrFile, "F ");}
    else if (res == -1)
    {elts[i] = 3; nt++; res = readGivenKeyword(ptrFile, "F "); }
    else { elts[i] = 4; nq++; res = readGivenKeyword(ptrFile, "F ");}
    i++;
  }
  
  //printf("I found %d tri %d quads\n", nt, nq);

  // Look for quads
  if (nq > 0)
  {
    cnq = new FldArrayI(nq, 4);
    KFSEEK(ptrFile, 0, SEEK_SET);
    res = 1; i = 0; k = 0;
    res = readGivenKeyword(ptrFile, "F ");
    while (res >= 1)
    {
      if (elts[k] == 4)
      {
        for (j = 1; j <= 4; j++)
        {
          res = readIntTuple(ptrFile, ti); (*cnq)(i,j) = ti;
        }
        i++;
      }
      k++;
      res = readGivenKeyword(ptrFile, "F ");
    }
  }

  // Look for tri
  if (nt > 0)
  {
    cnt = new FldArrayI(nt, 3);
    KFSEEK(ptrFile, 0, SEEK_SET);
    res = 1; i = 0; k = 0;
    res = readGivenKeyword(ptrFile, "F ");
    while (res >= 1)
    {
      if (elts[k] == 3)
      {
        for (j = 1; j <= 3; j++)
        {
          res = readIntTuple(ptrFile, ti); (*cnt)(i,j) = ti;
        }
        //printf("%d %d %d\n", (*cnt)(i,1), (*cnt)(i,2), (*cnt)(i,3));
        i++;
      }
      k++;
      res = readGivenKeyword(ptrFile, "F ");
    }
  }

  if (nt > 0)
  {
    unstructField.push_back(f);
    eltType.push_back(2);
    connect.push_back(cnt);
  }
  if (nq > 0)
  {
    if (nt > 0)
    {
      FldArrayF* f2 = new FldArrayF(*f);
      unstructField.push_back(f2);
    }
    else unstructField.push_back(f);
    eltType.push_back(3);
    connect.push_back(cnq);
  }
  // Cree les noms des zones
  for (unsigned int i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%d",i);
    zoneNames.push_back(zoneName);
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  
  return 0;
}

//=============================================================================
// Only write triangles, quads, tetra, hexa meshes. Others are discarded.
//=============================================================================
E_Int K_IO::GenIO::objwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  E_Int zone;
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (zone = 0; zone < nzone; zone++)
  {
    // triangles, quads supported
    if (eltType[zone] == 2 || eltType[zone] == 3)
      nvalidZones++;
    else
      printf("Warning: objwrite: zone %d not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: objwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Ouverture fichier
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: objwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Build writing data format
  char format1[40]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt);
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format of data
  sprintf(format1,"v %s%s%s\n", dataFmt, dataFmt, dataFmtl);

  // Write zones
  E_Int ng = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    // vertices
    FldArrayF& field = *unstructField[i];
    E_Float* fx = field.begin(posx);
    E_Float* fy = field.begin(posy);
    E_Float* fz = field.begin(posz);
    for (E_Int n = 0; n < field.getSize(); n++)
    {
      fprintf(ptrFile, format1, fx[n], fy[n], fz[n]);
    }
    if (eltType[i] == 2)
    {
      // faces
      FldArrayI& cn = *connect[i];
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        fprintf(ptrFile, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", 
                cn1[n]+ng, cn1[n]+ng, cn1[n]+ng,
                cn2[n]+ng, cn2[n]+ng, cn2[n]+ng,
                cn3[n]+ng, cn3[n]+ng, cn3[n]+ng);
      }
    }
    if (eltType[i] == 3)
    {
      // faces
      FldArrayI& cn = *connect[i];
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      E_Int* cn4 = cn.begin(4);
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        fprintf(ptrFile, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n", 
                cn1[n]+ng, cn1[n]+ng, cn1[n]+ng,
                cn2[n]+ng, cn2[n]+ng, cn2[n]+ng,
                cn3[n]+ng, cn3[n]+ng, cn3[n]+ng,
                cn4[n]+ng, cn4[n]+ng, cn4[n]+ng);
      }
    }
    ng += field.getSize();
  }

  fclose(ptrFile);
  return 0;
}
