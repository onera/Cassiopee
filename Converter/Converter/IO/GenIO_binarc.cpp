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

// Binary archive (CEDRE) file support

# include "GenIO.h"
# include <stdio.h>
# include <string.h>
# include "Array/Array.h"

using namespace std;
using namespace K_FLD;

// Read a block name jusqu'au caractere nul
void readBlockName(FILE* ptrFile, char* name)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
}

void readVersion(FILE* ptrFile, unsigned char* version)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { version[i] = c; i++; }
  version[i] = '\0';
}

void readTitle(FILE* ptrFile, char* titre, char* date, char* machine, 
               double& tempsZero, double& epsilon)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { titre[i] = c; i++; }
  titre[i] = '\0';
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { date[i] = c; i++; }
  date[i] = '\0';
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { machine[i] = c; i++; }
  machine[i] = '\0';
  fread(&tempsZero, sizeof(double), 1, ptrFile);
  fread(&epsilon, sizeof(double), 1, ptrFile);
}

void readGlobal(FILE* ptrFile, unsigned char& solverType, unsigned char& dimField)
{
  fread(&solverType, sizeof(unsigned char), 1, ptrFile);
  fread(&dimField, sizeof(unsigned char), 1, ptrFile); 
}

void readEspece(FILE* ptrFile, int& numMel, char* nomMel, int& nespeces, char* nomEspece)
{
  int c;
  // numero du melange
  fread(&numMel, sizeof(int), 1, ptrFile);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomMel[i] = c; i++; }
  nomMel[i] = '\0';
  fread(&nespeces, sizeof(int), 1, ptrFile);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomEspece[i] = c; i++; }
  nomEspece[i] = '\0';
}

void readScalar(FILE* ptrFile, int& numGrp, char* nomGrp, int& nelem, 
                char* nomScalar, unsigned char& typeSca)
{
  int c;
  fread(&numGrp, sizeof(int), 1, ptrFile);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomGrp[i] = c; i++; }
  nomGrp[i] = '\0';
  fread(&nelem, sizeof(int), 1, ptrFile);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomScalar[i] = c; i++; }
  nomScalar[i] = '\0';
  fread(&typeSca, sizeof(unsigned char), 1, ptrFile);
}

void readUnit(FILE* ptrFile, char* name, char* unit)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { unit[i] = c; i++; }
  unit[i] = '\0';
}

void readThermo(FILE* ptrFile, unsigned int& read, char* name,  
  unsigned int& nGe, unsigned int& nGr,
  unsigned int& ne, unsigned int& nr, 
  unsigned int*& dthGe, unsigned int*& dthGr,
  unsigned int*& dthe, unsigned int*& dthr)
{
  int c;
  E_Int i = 0;
  fread(&read, sizeof(unsigned int), 1, ptrFile);
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  if (read == 1) return;
  fread(&nGe, sizeof(unsigned int), 1, ptrFile);
  fread(&nGr, sizeof(unsigned int), 1, ptrFile);
  fread(&ne, sizeof(unsigned int), 1, ptrFile);
  fread(&nr, sizeof(unsigned int), 1, ptrFile);
  fread(&dthGe, sizeof(unsigned int), nGe, ptrFile);
  fread(&dthGr, sizeof(unsigned int), nGr, ptrFile);
  fread(&dthe, sizeof(unsigned int), ne, ptrFile);
  fread(&dthr, sizeof(unsigned int), nr, ptrFile); 
}
 
void readStructure(FILE* ptrFile, unsigned int& numabs, unsigned int& numuti,
  char* name, unsigned int& type, unsigned int& situ, int& numMel, int& numGrp)
{
  int c;
  E_Int i = 0;
  fread(&numabs, sizeof(unsigned int), 1, ptrFile);
  fread(&numuti, sizeof(unsigned int), 1, ptrFile);
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  fread(&type, sizeof(unsigned int), 1, ptrFile);
  fread(&situ, sizeof(unsigned int), 1, ptrFile);
  fread(&numMel, sizeof(int), 1, ptrFile);
  fread(&numGrp, sizeof(int), 1, ptrFile);
}

void readDomutil(FILE* ptrFile, unsigned int& numUti, char* name)
{
  int c;
  E_Int i = 0;
  fread(&numUti, sizeof(unsigned int), 1, ptrFile);
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
}


//=============================================================================
/* 
   arcread
   IN: file: file name,
   OUT: varString: variables string
   OUT: structField: field for each structured zones,
   OUT: ni, nj, nk: number of points of each structured zones,
   OUT: unstructField: field for each unstructured zones, 
   OUT: connectivity: connectivity for each unstructured zones,
   OUT: eltType: eltType for each unstructured zones.

   eltType is:
   1: BAR
   2: TRI
   3: QUAD
   4: TETRA
   7: HEXA
   8: NGON
   return 1 if failure.
   return 0 if ok.

*/
//=============================================================================
E_Int K_IO::GenIO::arcread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connectivity,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: arcread: cannot open file %s.\n", file);
    return 1;
  }

  printf("Error: arcread: not implemented.\n");
  return 0;
}

//=============================================================================
/*
  This routine enables binary tecplot of field.
*/
//=============================================================================
E_Int K_IO::GenIO::arcwrite(char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames)
{
  printf("Error: arcwrite: not implemented.\n");
  return 1;
}

