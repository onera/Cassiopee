/*    
    Copyright 2013-2020 Onera.

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

// for dbx
#include <iostream>

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
  printf("blockname: %s\n", name);
}

void readVersion(FILE* ptrFile, unsigned char* version)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { version[i] = c; i++; }
  version[i] = '\0';
  printf("version: %u\n", *version);
}

void readTitle(FILE* ptrFile, char* titre, char* date, char* machine, 
               double& tempsZero, double& epsilon)
{
  int c; E_Int i;
  if (titre != NULL)
  {
    i = 0;
    while ((c = fgetc(ptrFile)) != '\0')
    { titre[i] = c; i++; }
    titre[i] = '\0';
  }
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { date[i] = c; i++; }
  date[i] = '\0';
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { machine[i] = c; i++; }
  machine[i] = '\0';
  fread(&tempsZero, sizeof(double), 1, ptrFile); tempsZero = DBE(tempsZero);
  fread(&epsilon, sizeof(double), 1, ptrFile); epsilon = DBE(epsilon);
  
  if (titre != NULL) printf("titre: %s\n", titre);
  printf("date: %s\n", date);
  printf("machine: %s\n", machine);
  printf("temps0: %g\n", tempsZero);
  printf("epsilon: %g\n", epsilon);
}

void readGlobal(FILE* ptrFile, unsigned char& solverType, unsigned char& dimField)
{
  fread(&solverType, sizeof(unsigned char), 1, ptrFile);
  fread(&dimField, sizeof(unsigned char), 1, ptrFile);
  printf("solverType: %u\n", solverType);
  printf("dimField: %u\n", dimField);
}

void readEspece(FILE* ptrFile, int& numMel, char* nomMel, int& nespeces, char* nomEspece)
{
  int c;
  // numero du melange
  fread(&numMel, sizeof(int), 1, ptrFile); numMel = IBE(numMel);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomMel[i] = c; i++; }
  nomMel[i] = '\0';
  fread(&nespeces, sizeof(int), 1, ptrFile); nespeces = IBE(nespeces);
  i = 0;
  for (E_Int j = 0; j < nespeces; j++)
  {
    while ((c = fgetc(ptrFile)) != '\0')
    { nomEspece[i] = c; i++; }
    nomEspece[i] = '\0';
  }
  printf("numMel: %d\n", numMel);
  printf("nomMel: %s\n", nomMel);
  printf("nbreEspeces: %d\n", nespeces);
  for (E_Int i = 0; i < nespeces; i++) printf("%d %s\n", i, nomEspece);
}

void readScalar(FILE* ptrFile, int& numGrp, char* nomGrp, int& nelem, 
                char* nomScalar, unsigned char& typeSca)
{
  int c;
  fread(&numGrp, sizeof(int), 1, ptrFile); numGrp = IBE(numGrp);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomGrp[i] = c; i++; }
  nomGrp[i] = '\0';
  fread(&nelem, sizeof(int), 1, ptrFile); nelem = IBE(nelem);
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
  printf("name %s\n", name);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { unit[i] = c; i++; }
  unit[i] = '\0';
  printf("unit %s\n", unit);
}

void readThermo(FILE* ptrFile, unsigned int& read, char* name,  
  unsigned int& nGe, unsigned int& nGr,
  unsigned int& ne, unsigned int& nr, 
  unsigned int*& dthGe, unsigned int*& dthGr,
  unsigned int*& dthe, double*& dthr)
{
  fread(&read, sizeof(unsigned int), 1, ptrFile); read = UIBE(read);
  printf("thermo: read %u\n", read);
  
  int c;
  E_Int i = 0;
  while ( (c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  printf("thermo name : %s\n", name);
  if (read == 1) return;
  fread(&nGe, sizeof(unsigned int), 1, ptrFile); nGe = UIBE(nGe);
  fread(&nGr, sizeof(unsigned int), 1, ptrFile); nGr = UIBE(nGr);
  fread(&ne, sizeof(unsigned int), 1, ptrFile); ne = UIBE(ne);
  fread(&nr, sizeof(unsigned int), 1, ptrFile); nr = UIBE(nr);
  printf("nGe: %d %d %d %d\n", nGe, nGr, ne, nr);
  dthGe = new unsigned int [nGe];
  dthGr = new unsigned int [nGr];
  dthe = new unsigned int [ne];
  dthr = new double [nr];
  fread(&dthGe, sizeof(unsigned int), nGe, ptrFile); 
  fread(&dthGr, sizeof(unsigned int), nGr, ptrFile);
  fread(&dthe, sizeof(unsigned int), ne, ptrFile);
  fread(&dthr, sizeof(double), nr, ptrFile);
}
 
void readStructure(FILE* ptrFile, unsigned int& numabs, unsigned int& numuti,
  char* name, unsigned char& type, unsigned char& situ, int& numMel, int& numGrp)
{
  int c;
  E_Int i = 0;
  fread(&numabs, sizeof(unsigned int), 1, ptrFile); numabs = UIBE(numabs);
  fread(&numuti, sizeof(unsigned int), 1, ptrFile); numuti = UIBE(numuti);
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  fread(&type, sizeof(unsigned char), 1, ptrFile); 
  fread(&situ, sizeof(unsigned char), 1, ptrFile); 
  fread(&numMel, sizeof(int), 1, ptrFile); numMel = IBE(numMel);
  fread(&numGrp, sizeof(int), 1, ptrFile); numGrp = IBE(numGrp);
  printf("Structure: %u %u\n", numabs, numuti);
  printf("nom dom: %s\n", name);
  printf("type : %u, situ = %u\n", type, situ);
}

void readDomutil(FILE* ptrFile, unsigned int& numUti, char* nomUti)
{
  int c;
  E_Int i = 0;
  fread(&numUti, sizeof(unsigned int), 1, ptrFile); numUti = UIBE(numUti);
  while ((c = fgetc(ptrFile)) != '\0')
  { nomUti[i] = c; i++; }
  nomUti[i] = '\0';
  printf("nom du domaine %s\n", nomUti);
}

void readNumerotation(FILE* ptrFile, unsigned int& numabs, double& temps,
  int& nsom, int& nfac, int& nflm, int &nflp, int& nceli, int*& isomu,
  int*& ifacu, int*& icelu)
{
  fread(&numabs, sizeof(unsigned int), 1, ptrFile); numabs = UIBE(numabs);
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  fread(&nsom, sizeof(int), 1, ptrFile); nsom = IBE(nsom);
  fread(&nfac, sizeof(int), 1, ptrFile); nfac = IBE(nfac);
  fread(&nflm, sizeof(int), 1, ptrFile); nflm = IBE(nflm);
  fread(&nflp, sizeof(int), 1, ptrFile); nflp = IBE(nflp);
  fread(&nceli, sizeof(int), 1, ptrFile); nceli = IBE(nceli);
  printf("%d %d %d\n", nsom, nfac, nflm);
  isomu = new int [nsom];
  fread(isomu, sizeof(int), nsom, ptrFile);
  for (E_Int i = 0; i < nsom; i++) isomu[i] = IBE(isomu[i]);

  ifacu = new int [nfac];
  fread(ifacu, sizeof(int), nfac, ptrFile);
  for (E_Int i = 0; i < nfac; i++) ifacu[i] = IBE(ifacu[i]);
  
  E_Int size = nceli+nflm+nflp;
  icelu = new int [size];
  fread(icelu, sizeof(int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) icelu[i] = IBE(icelu[i]);
}

void readVariables(FILE* ptrFile, unsigned int& numUti, unsigned int& nvar, char* nomVar,
  char* nomGrp, char* nomCat, unsigned char& glob, unsigned int& elmin,
  unsigned int& nelem, double*& data)
{
  fread(&numUti, sizeof(unsigned int), 1, ptrFile); numUti = UIBE(numUti);
  fread(&nvar, sizeof(unsigned int), 1, ptrFile); nvar = UIBE(nvar);
  
}


void readConnexion(FILE* ptrFile, unsigned int& nd, double& temps, unsigned int& elmin,
  unsigned int& nelem, unsigned int*& num, unsigned char* nbpoints, unsigned int*& numpoints)
{
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  fread(&elmin, sizeof(unsigned int), 1, ptrFile); elmin = UIBE(elmin);
  fread(&nelem, sizeof(unsigned int), 1, ptrFile); nelem = UIBE(nelem);
  
  E_Int size = 2*nelem;
  num = new unsigned int [size];
  fread(num, sizeof(unsigned int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) num[i] = UIBE(num[i]);
  
  nbpoints = new unsigned char [nelem];
  fread(nbpoints, sizeof(unsigned char), nelem, ptrFile);
  
  E_Int nbPointsTot = 0;
  for (size_t i = 0; i < nelem; i++) nbPointsTot += (E_Int)nbpoints[i];
    
  size = nbPointsTot;
  numpoints = new unsigned int [size];
  fread(numpoints, sizeof(unsigned int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) numpoints[i] = UIBE(numpoints[i]);
  
}

// Lit un maillage structure
void readMaillage0(FILE* ptrFile)
{
  
}

// Lit un maillage non structure
void readMaillage1(FILE* ptrFile, unsigned int& nd)
{
    
}

//=============================================================================
/* 
   arcread - read binary cedre format (archive)
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

  //printf("Error: arcread: not implemented.\n");
  unsigned char version;
  char name[256];
  char fileName[256];
  char unitName[256];
  char nomUti[256];
  char titre[256];
  char date[256];
  char machine[256];
  double tempsZero;
  double epsilon;
  unsigned char solverType;
  unsigned char dimField;
  int numMel;
  char nomMel[256];
  int nespeces;
  char nomEspece[256*10];
  int numGrp;
  char nomGrp[256];
  int nelem;
  char nomScalar[256];
  unsigned char typeSca;
  char unit[256];
  unsigned int read;
  unsigned int numUti;
  unsigned int numabs;
  unsigned int numuti;
  unsigned char type;
  unsigned char situ;
  unsigned int nGe, nGr, ne, nr;
  unsigned int* dthGe=NULL;
  unsigned int *dthGr=NULL;
  unsigned int *dthe=NULL;
  double *dthr=NULL;
  double temps;
  int nsom;
  int nfac;
  int nflm;
  int nflp;
  int nceli;
  int* isomu = NULL;
  int* ifacu = NULL;
  int* icelu = NULL;
  unsigned int nd;
  unsigned int elmin;
  unsigned int nelems;
  unsigned int* num = NULL;
  unsigned char* nbpoints = NULL;
  unsigned int* numpoints = NULL;
  
  readBlockName(ptrFile, name);
  readVersion(ptrFile, &version);
  
  while (! feof(ptrFile))
  {
    readBlockName(ptrFile, name);
    if (strcmp(name, "VERSION") == 0) readVersion(ptrFile, &version);
    else if (strcmp(name, "TITRE") == 0) readTitle(ptrFile, titre, date, machine, tempsZero, epsilon);
    else if (strcmp(name, "SANS TITRE") == 0) readTitle(ptrFile, NULL, date, machine, tempsZero, epsilon);
    else if (strcmp(name, "GLOBAL") == 0) readGlobal(ptrFile, solverType, dimField);
    else if (strcmp(name, "ESPECES") == 0) readEspece(ptrFile, numMel, nomMel, nespeces, nomEspece);
    else if (strcmp(name, "SCALAIRES") == 0) readScalar(ptrFile, numGrp, nomGrp, nelem, nomScalar, typeSca);
    else if (strcmp(name, "UNITE") == 0) readUnit(ptrFile, unitName, unit);
    else if (strcmp(name, "THERMODYNAMIQUE") == 0) readThermo(ptrFile, read, fileName, nGe, nGr, ne, nr, dthGe, dthGr, dthe, dthr);
    else if (strcmp(name, "DOMUTIL") == 0) readDomutil(ptrFile, numUti, nomUti);
    else if (strcmp(name, "STRUCTURE") == 0) readStructure(ptrFile, numabs, numuti, name, type, situ, numMel, numGrp);
    else if (strcmp(name, "NUMEROTATION") == 0) readNumerotation(ptrFile, numabs, temps, nsom, nfac, nflm, nflp, nceli, isomu, ifacu, icelu);
    else if (strcmp(name, "CONNEXION") == 0) readConnexion(ptrFile, nd, temps, elmin, nelems, num, nbpoints, numpoints);

    else { printf("Block pas encore implemente: %s\n", name); break; }
  }
  delete [] dthGe; delete [] dthGr; delete [] dthe; delete [] dthr;
  delete [] isomu; delete [] ifacu; delete [] icelu;
  delete [] num; delete [] nbpoints; delete [] numpoints;
  return 0;
}

//=============================================================================
/*
  This routine enables binary write of archive file.
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

