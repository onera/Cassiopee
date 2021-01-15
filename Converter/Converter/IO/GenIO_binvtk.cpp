/*    
    Copyright 2013-2021 Onera.

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

// Binary vtk (legacy) file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

// retourne 0 si ok, 1 sinon
int readLine(FILE* ptrFile, char* buf)
{
  char c = '1'; E_Int i = 0;
  while (c != '\n')
  { 
    if (feof(ptrFile)) return 1;
    c = fgetc(ptrFile);
    buf[i] = c; i++; 
  }
  if (i > 0) buf[i-1] = '\0';
  else buf[0] = '\0';
  printf("LINE: %s\n", buf);
  return 0;
}

void getFirstWord(char* buf, char* word, E_Int& length)
{
  E_Int i = 0;
  while (buf[i] != '\n')
  {
    word[i] = buf[i]; i++;
  }
  if (i > 0) { word[i-1] = '\0'; length = i-1; }
  else { word[0] = '\0'; length = 0; }
  printf("WORD:%s\n", word);
}

void readPoints(FILE* ptrFile, E_Boolean changeEndian, E_Int& npts, float*& b)
{
  char buf[256];
  readLine(ptrFile, buf);
  E_Int l = strlen(buf);
  E_Int count = 0;
  while (buf[count] != ' ' && buf[count] != '\n') count++;
  for (E_Int i = 0; i <= l-count; i++) buf[i] = buf[i+count];
  npts = 0;
  sscanf(buf, "%d", &npts);
  printf("npts: %d\n", npts);
  b = new float [npts*3];
  fread(b, sizeof(float), 3*npts, ptrFile);
  if (changeEndian)
  {
    for (E_Int i = 0; i < 3*npts; i++) b[i] = FBE(b[i]);
  }
  printf("point %f\n", b[0]);
  fgetc(ptrFile); // avoid \n
}

void readCells(FILE* ptrFile, E_Boolean changeEndian, E_Int& ncells, E_Int& size, int*& cells)
{
  char buf[256]; char buf2[256];
  readLine(ptrFile, buf);
  // Recupere le nombre de cellules
  E_Int i = 6;
  while (buf[i] != '\0' && buf[i] != ' ') { buf2[i-6] = buf[i]; i++; }
  buf2[i-6] = '\0';
  ncells = 0;
  sscanf(buf2, "%d", &ncells);
  printf("ncells: %d\n", ncells);
  // recupere la taille totale de la connectivite
  i += 1;
  E_Int is = i;
  while (buf[i] != '\0' && buf[i] != ' ') { buf2[i-is] = buf[i]; i++; }
  buf2[i-is] = '\0';
  sscanf(buf2, "%d", &size);
  printf("size: %d\n", size);
  cells = new int [size];
  fread(cells, sizeof(int), size, ptrFile);
  if (changeEndian)
  {
    for (E_Int i = 0; i < size; i++) cells[i] = IBE(cells[i]);
  }
  printf("cells=%d\n", cells[0]);
  fgetc(ptrFile); // avoid \n
}

void readCellTypes(FILE* ptrFile, E_Boolean changeEndian, E_Int ncells, int*& cellTypes)
{
  char buf[256];
  readLine(ptrFile, buf);
  cellTypes = new int [ncells];
  fread(cellTypes, sizeof(int), ncells, ptrFile);
  if (changeEndian)
  {
    for (E_Int i = 0; i < ncells; i++) cellTypes[i] = IBE(cellTypes[i]);
  } 
  printf("cellTypes=%d\n", cellTypes[0]);
  fgetc(ptrFile); // avoid \n
}

void readScalar(FILE* ptrFile, E_Boolean changeEndian, E_Int npts, char* varName, float*& field)
{ 
  char buf[256];
  readLine(ptrFile, buf);
  // Recupere le nom du champ
  E_Int i = 8;
  while (buf[i] != '\0' && buf[i] != ' ') { varName[i-8] = buf[i]; i++; }
  varName[i-8] = '\0';
  printf("varname=%s\n", varName);
  // Passe lookup table
  readLine(ptrFile, buf);
  field = new float [npts];
  fread(field, sizeof(float), npts, ptrFile);
  if (changeEndian)
  {
    for (E_Int i = 0; i < npts; i++) field[i] = FBE(field[i]);
  }
  fgetc(ptrFile); // avoid \n
}

void readVector(FILE* ptrFile, E_Boolean changeEndian, E_Int npts, char* varName, float*& field)
{ 
  char buf[256];
  readLine(ptrFile, buf);
  // Recupere le nom du champ
  E_Int i = 8;
  while (buf[i] != '\0' && buf[i] != ' ') { varName[i-8] = buf[i]; i++; }
  varName[i-8] = '\0';
  printf("varname=%s\n", varName);
  field = new float [3*npts];
  fread(field, sizeof(float), 3*npts, ptrFile);
  if (changeEndian)
  {
    for (E_Int i = 0; i < 3*npts; i++) field[i] = FBE(field[i]);
  }
  fgetc(ptrFile); // avoid \n
}


//=============================================================================
/* binvtkread
*/
//=============================================================================
E_Int K_IO::GenIO::binvtkread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Boolean changeEndian = true;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: vtkread: cannot open file %s.\n", file);
    return 1;
  }

  E_Int i; char c;

  // Lecture de la version (jusqu'a newline)
  char buf[256];
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  buf[22] = '\0';
  printf("version: %s\n", buf);
  if (strcmp(buf, "# vtk DataFile Version") != 0) { fclose(ptrFile); return 1; }
  
  // Lecture de la description (jusqu'a newline)
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  printf("description: %s\n", buf);
  
  // Type de format BINARY or ASCII
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i-1] = '\0';
  printf("DATA: %s\n", buf);
  if (strcmp(buf, "BINARY") != 0) { fclose(ptrFile); return 1; }
  // must be binary here
  
  // Type de Data (STRUCTURED_POINTS, STRUCTURED_GRID, RECTILINEAR_GRID, POLYDATA, UNSTRUCTURED_GRID)
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i-1] = '\0';
  printf("TYPE: %s\n", buf);  
  E_Int l = strlen(buf);
  for (E_Int i = 0; i <= l-8; i++) buf[i] = buf[i+8];
  printf("TYPE: %s\n", buf);

  if (strcmp(buf, "FIELD") != 0)
  {
    // Read POINTS
    E_Int npts; float* pts;
    readPoints(ptrFile, changeEndian, npts, pts);
    
    // READ CELLS
    int* cells; int* cellTypes; E_Int ncells; E_Int size;
    readCells(ptrFile, changeEndian, ncells, size, cells);
    readCellTypes(ptrFile, changeEndian, ncells, cellTypes);
    
    // READ POINT_DATA
    E_Int ret = readLine(ptrFile, buf);
    
    // READ SCALARS if ANY
    vector<float*> fields;
    vector<E_Int> type; // 0: scalar, 1: vector
    char varName[256];
    varString = new char [1200];
    strcpy(varString, "x,y,z");
    
    ret = readLine(ptrFile, buf);
    
    if (ret == 0)
    {
      float* ff;
      readScalar(ptrFile, changeEndian, npts, varName, ff);
      strcat(varString, ",");
      strcat(varString, varName);
      fields.push_back(ff);
      type.push_back(0);
    }
    
    if (ret == 0)
    {
      float* ff;
      readVector(ptrFile, changeEndian, npts, varName, ff);
      strcat(varString, ",");
      strcat(varString, varName);
      strcat(varString, "X,");
      strcat(varString, varName);
      strcat(varString, "Y,");
      strcat(varString, varName);
      strcat(varString, "Z");
      fields.push_back(ff);
      type.push_back(1);
    }
    
    if (ret == 0)
    {
      float* ff;
      readScalar(ptrFile, changeEndian, npts, varName, ff);
      strcat(varString, ",");
      strcat(varString, varName);
      fields.push_back(ff);
      type.push_back(0);
    }
    
    E_Int nvars = 0;
    for (size_t i = 0; i < type.size(); i++)
    {
      if (type[i] == 0) nvars += 1;
      else if (type[i] == 1) nvars += 3;
    }
    printf("nvars=%d\n", nvars);
    // Concatenate arrays
    FldArrayF* f = new FldArrayF(npts, 3+nvars);
    E_Float* fx = f->begin(1);
    E_Float* fy = f->begin(2);
    E_Float* fz = f->begin(3);
    
    for (E_Int i = 0; i < npts; i++)
    {
      fx[i] = pts[3*i];
      fy[i] = pts[3*i+1];
      fz[i] = pts[3*i+2];
    }
    delete [] pts;
    
    E_Int np = 4;
    for (size_t n = 0; n < type.size(); n++)
    {
      if (type[n] == 0)
      {
        E_Float* fp = f->begin(np);
        for (E_Int i = 0; i < npts; i++) fp[i] = fields[n][i];
        np += 1;
      }
      else if (type[n] == 1)
      {
        E_Float* fx = f->begin(np);
        E_Float* fy = f->begin(np+1);
        E_Float* fz = f->begin(np+2);
        
        for (E_Int i = 0; i < npts; i++) fx[i] = fields[n][3*i];
        for (E_Int i = 0; i < npts; i++) fy[i] = fields[n][3*i+1];
        for (E_Int i = 0; i < npts; i++) fz[i] = fields[n][3*i+2];
          
        np += 3; 
      }
    }
    
    for (size_t n = 0; n < fields.size(); n++) delete [] fields[n];
    
    // tri suivant les elements types
    E_Int nNODE = 0; // 1
    E_Int nBAR = 0; // 3
    E_Int nTRI = 0; // 5
    E_Int nQUAD = 0; // 9
    E_Int nHEXA = 0; // 12
    E_Int nTETRA = 0; // 10
    E_Int nPENTA = 0; // 13
    E_Int nPYRA = 0; // 14
    for (E_Int i = 0; i < ncells; i++)
    {
      // type de la cellule
      switch (cellTypes[i])
      {
        case 1:
          nNODE++; break;
        case 3:
          nBAR++; break;
        case 5: 
          nTRI++; break;
        case 9:
          nQUAD++; break;
        case 10:
          nTETRA++; break;
        case 12:
          nHEXA++; break;
        case 13:
          nPENTA++; break;
        case 14:
          nPYRA++; break;
        default: ;
      }
    }
    // Allocate
    printf("Elements: BAR=%d TRI=%d QUAD=%d HEXA=%d TETRA=%d PENTA=%d PYRA=%d\n", nBAR, nTRI, nQUAD, nHEXA, nTETRA, nPENTA, nPYRA);
    
    E_Int* cnBAR=NULL; E_Int* cnTRI=NULL; E_Int* cnQUAD=NULL;
    E_Int* cnTETRA=NULL; E_Int* cnPYRA=NULL; E_Int* cnPENTA=NULL;
    E_Int* cnHEXA=NULL;
    FldArrayI* cn;
    
    if (nBAR > 0) 
    { 
      cn = new FldArrayI(nBAR, 2); cnBAR = cn->begin(); 
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(1);
    }
    if (nTRI > 0) 
    { 
      cn = new FldArrayI(nTRI, 3); cnTRI = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(2); 
    }
    if (nQUAD > 0) 
    { 
      cn = new FldArrayI(nQUAD, 4); cnQUAD = cn->begin(); 
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(3);
    }
    if (nHEXA > 0) 
    { 
      cn = new FldArrayI(nHEXA, 8); cnHEXA = cn->begin(); 
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(7);
    }
    if (nTETRA > 0) 
    { 
      cn = new FldArrayI(nTETRA, 4); cnTETRA = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(4); 
    }
    if (nPENTA > 0) 
    { 
      cn = new FldArrayI(nPENTA, 6); cnPENTA = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(6);
    }
    if (nPYRA > 0) 
    { 
      cn = new FldArrayI(nPYRA, 5); cnPYRA = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(5);
    }

    int* ptr = cells;
    E_Int pBAR = 0; E_Int pPENTA = 0;
    E_Int pHEXA = 0; E_Int pTRI = 0; E_Int pQUAD = 0;
    E_Int pPYRA = 0; E_Int pTETRA = 0;
    for (E_Int i = 0; i < ncells; i++)
    {
      // type de la cellule
      switch (cellTypes[i])
      {
        case 3:
          cnBAR[pBAR] = ptr[1]+1;
          cnBAR[pBAR+nBAR] = ptr[2]+1;
          pBAR++;
          break;
        case 5:
          cnTRI[pTRI] = ptr[1]+1;
          cnTRI[pTRI+nTRI] = ptr[2]+1;
          cnTRI[pTRI+2*nTRI] = ptr[3]+1;
          pTRI++;
          break;
        case 9:
          cnQUAD[pQUAD] = ptr[1]+1;
          cnQUAD[pQUAD+nQUAD] = ptr[2]+1;
          cnQUAD[pQUAD+2*nQUAD] = ptr[3]+1;
          cnQUAD[pQUAD+3*nQUAD] = ptr[4]+1;
          pQUAD++;
          break;
        case 10:
          cnTETRA[pTETRA] = ptr[1]+1;
          cnTETRA[pTETRA+nTETRA] = ptr[2]+1;
          cnTETRA[pTETRA+2*nTETRA] = ptr[3]+1;
          cnTETRA[pTETRA+3*nTETRA] = ptr[4]+1;
          pTETRA++;
          break;
        case 12:
          cnHEXA[pHEXA] = ptr[1]+1;
          cnHEXA[pHEXA+nHEXA] = ptr[2]+1;
          cnHEXA[pHEXA+2*nHEXA] = ptr[3]+1;
          cnHEXA[pHEXA+3*nHEXA] = ptr[4]+1;
          cnHEXA[pHEXA+4*nHEXA] = ptr[5]+1;
          cnHEXA[pHEXA+5*nHEXA] = ptr[6]+1;
          cnHEXA[pHEXA+6*nHEXA] = ptr[7]+1;
          cnHEXA[pHEXA+7*nHEXA] = ptr[8]+1;
          pHEXA++;
          break;
        case 13:
          cnPENTA[pPENTA] = ptr[1]+1;
          cnPENTA[pPENTA+nPENTA] = ptr[2]+1;
          cnPENTA[pPENTA+2*nPENTA] = ptr[3]+1;
          cnPENTA[pPENTA+3*nPENTA] = ptr[4]+1;
          cnPENTA[pPENTA+4*nPENTA] = ptr[5]+1;
          cnPENTA[pPENTA+5*nPENTA] = ptr[6]+1;
          pPENTA++;
          break;
        case 14:
          cnPYRA[pPYRA] = ptr[1]+1;
          cnPYRA[pPYRA+nPYRA] = ptr[2]+1;
          cnPYRA[pPYRA+2*nPYRA] = ptr[3]+1;
          cnPYRA[pPYRA+3*nPYRA] = ptr[4]+1;
          cnPYRA[pPYRA+4*nPYRA] = ptr[5]+1;
          pPYRA++;
          break;
        default: ;
      }
      ptr = ptr+(ptr[0]+1);
    }
    
    delete [] cells;
    delete [] cellTypes;  
  }
  
  fclose(ptrFile);

  // Cree le nom de zone
  for (unsigned int i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%d", i);
    zoneNames.push_back(zoneName);
  }
  
  return 0;
}

//=============================================================================
// Only write NODE, BAR, TRI, QUAD, TETRA, HEXA, PYRA, PENTA.
//=============================================================================
E_Int K_IO::GenIO::binvtkwrite(
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
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // NODE, BAR, TRI, QUADS, TETRA, HEXA , PENTA, PYRA, supported
    if (eltType[zone] < 8) nvalidZones++;
    else
      printf("Warning: binvtkwrite: zone %d not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: binvtkwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  printf("Warning: binvtkwrite: not implemented.");
  return 1;
}
