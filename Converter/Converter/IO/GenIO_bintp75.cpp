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

// Binary tecplot v75 file support

#include "GenIO.h"
#include <stdio.h>
#include <stdlib.h>
# include <string.h>

#include "Def/DefFunction.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Read zone header.
   return 1 if pb else return 0 

   IN: ptrFile: opened stream,
   
   if zone is structured then returns :
   OUT: dim, ni, nj, nk
   OUT: zoneName
   OUT: dataPacking: 0 (block), 1 (point)

   if zone is unstructured, returns:
   OUT: npts: number of points
   OUT: nelts: number of elements
   OUT: eltType: type of elements 
   (2: TRI, 3: QUAD, 4: TETRA, 7: HEXA).
   OUT: zoneName
   OUT: dataPacking: 0 (block), 1 (point)
   OUT: geometries (1D fields)

   v75 compatible.
*/
//=============================================================================
E_Int K_IO::GenIO::readZoneHeader75(FILE* ptrFile, E_Int& dim, 
                                    E_Int& ni, E_Int& nj, E_Int& nk,
                                    E_Int& npts, E_Int& nelts, E_Int& eltType,
                                    char* zoneName, E_Int& dataPacking,
                                    vector<FldArrayF*>& geom)
{
  char dummy[BUFSIZE+1];
  float a;
  int ib;
  E_Int i; E_Int fmt;

  /* Constants */
  E_Int si = sizeof(int);

  /* Zone marker (299) */
  fread(&a, sizeof(float), 1, ptrFile);
  if (K_FUNC::fEqualZero(a - 299.) == false) return 1;
  
  /* Zone name */
  i = 0; ib = 1;
  while (ib != 0)
  {
    fread(&ib, si, 1, ptrFile);
    dummy[i] = (char)ib;
    i++;
  }
  strcpy(zoneName, dummy);
  
  /* Format of zone */
  fread(&ib, si, 1, ptrFile);
  fmt = ib;
  if (ib <= 1) dataPacking = ib;
  else dataPacking = ib-2;

  /* zone color */
  fread(&ib, si, 1, ptrFile);
  
  switch (fmt)
  {
    case 0: // structure, format block
    case 1: // structure, format point
      /* ni, nj, nk */
      fread(&ib, si, 1, ptrFile);
      if (ib > 1) dim = 1;
      ni = ib;
      fread(&ib, si, 1, ptrFile);
      if (ib > 1) dim = 2;
      nj = ib;
      fread(&ib, si, 1, ptrFile);
      if (ib > 1) dim = 3;
      nk = ib;
      break;
   
    case 2: // non structure, format feblock
    case 3: // non structure, format fepoint
      fread(&ib, si, 1, ptrFile);
      npts = ib;
      fread(&ib, si, 1, ptrFile);
      nelts = ib;
      fread(&ib, si, 1, ptrFile);
      if (ib == 4) eltType = 1; // BAR
      else if (ib == 0) eltType = 2; // TRI
      else if (ib == 1) eltType = 3; // QUAD
      else if (ib == 2) eltType = 4; // TETRA
      else if (ib == 3) eltType = 7; // HEXA
      break;

    default:
      printf("#readZoneHeader: unknown type of elements.\n");
      exit(0);
  }
  return 0;
}

//=============================================================================
/* Read zone header
   return 1 if pb else return 0.
   This is the v75 version.
   if fmt is structured then returns dim, ni,nj ,nk
   else return npts: number of points
               nelts: number of elements
               eltType: type of elements 
               (2: TRI, 3: QUAD, 4: TETRA, 7: HEXA)
*/
//=============================================================================
E_Int K_IO::GenIO::readZoneHeader75CE(
  FILE* ptr_file, E_Int& dim, 
  E_Int& ni, E_Int& nj, E_Int& nk,
  E_Int& npts, E_Int& nelts, E_Int& eltType,
  char* zoneName, E_Int& dataPacking, vector<FldArrayF*>& geom)
{
  char dummy[BUFSIZE+1];
  float a;
  int ib;
  E_Int i;
  E_Int fmt;// format : block, point, feblock, fepoint

  /* Constants */
  E_Int si = sizeof(int);
  
  /* Zone marker (299) */
  fread(&a, sizeof(float), 1, ptr_file);
  a = FBE(a);
  if (K_FUNC::fEqualZero(a - 299.) == false) return 1;
  
  /* Zone name */
  i = 0; ib = 1;
  while (ib != 0 && i < BUFSIZE)
  {
    fread(&ib, si, 1, ptr_file);
    dummy[i] = (char)IBE(ib);
    i++;
  }
  if (i == BUFSIZE) dummy[BUFSIZE] = '\0';
  strcpy(zoneName, dummy);
  
  /* Format de la zone */
  fread(&ib, si, 1, ptr_file);
  ib = IBE(ib);
  fmt = ib;
  if ( ib <= 1 ) dataPacking = ib;
  else dataPacking = ib-2;

  /* zone color */
  fread(&ib, si, 1, ptr_file);

  switch (fmt)
  {
    case 0: //block
    case 1: //point
      /* ni, nj, nk */
      fread(&ib, si, 1, ptr_file);
      ib = IBE(ib);
      if (ib > 1) dim = 1;
      ni = ib;
      fread(&ib, si, 1, ptr_file);
      ib = IBE(ib);
      if (ib > 1) dim = 2;
      nj = ib;
      fread(&ib, si, 1, ptr_file);
      ib = IBE(ib);
      if (ib > 1) dim = 3;
      nk = ib;
      break;
   
    case 2://feblock
    case 3://fepoint
      fread(&ib, si, 1, ptr_file);
      npts = IBE(ib);
      fread(&ib, si, 1, ptr_file);
      nelts = IBE(ib);
      fread(&ib, si, 1, ptr_file);
      if (IBE(ib) == 0) eltType = 2; // TRI
      else if (IBE(ib) == 1) eltType = 3; // QUAD
      else if (IBE(ib) == 2) eltType = 4; // TETRA
      else if (IBE(ib) == 3) eltType = 7; // HEXA
      break;
    default:
      printf("#readzoneheader: unknown type of elements.\n");
      exit(0); // to fix
  }  
  return 0;
}

//=============================================================================
/* 
   readData pour les zones structurees
   IN: ptrFile: ptr on file being read
   IN: ni, nj, nk: size of zone.
   IN: dataPacking: 0 (block), 1 (point)
   OUT: f: field read. Must be already dimensioned.
 */
//=============================================================================
E_Int K_IO::GenIO::readData75(FILE* ptrFile, E_Int ni, E_Int nj, E_Int nk, 
                              E_Int dataPacking,
                              FldArrayF& f)
{
  float a;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = f.getNfld();
  E_Int si = sizeof(int);
  E_Int npts = ni*nj*nk;
      
  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  
  fread(&ib, si, 1, ptrFile); /* repeat */
  if (ib != 0)
  {
    printf("readData: repeat is not supported.\n");
  }
  
  // Variable type
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile);
    if (ib == 1) sizer = 4;
  }
  
  /* Read dump */
  if (sizer == 4 && dataPacking == 0) // block
  {
    float* buf = new float [npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(float), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i,n+1) = buf[i];
    }
    delete [] buf;
  }
  else if (sizer == 8 && dataPacking == 0) // block
  {
    for (n = 0; n < nfield; n++)
    {
      fread(f.begin(n+1), sizeof(E_Float), npts, ptrFile);
    } 
  }
  else if (sizer == 4 && dataPacking == 1) // point
  {
    float* buf = new float[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = buf[i];
    }
    delete [] buf;
  }
  else if (sizer == 8 && dataPacking == 1) // point
  {
    double* buf = new double[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(E_Float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = buf[i];
    }
    delete [] buf;
  }

  return 0;
}

//=============================================================================
/* 
   readData pour les zones non structurees.
   IN: ptrFile: ptr on file being read
   IN: dataPacking: 0 (block), 1 (point)
   OUT: f: field read. Must be already dimensioned.
        c: connectivity read. Must be already dimensioned.
 */
//=============================================================================
E_Int K_IO::GenIO::readData75(FILE* ptrFile, 
                              E_Int dataPacking, 
                              FldArrayF& f, FldArrayI& c)
{
  float a;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = f.getNfld();
  E_Int npts = f.getSize();
  E_Int nelts = c.getSize();
  E_Int eltType = c.getNfld();
  E_Int si = sizeof(int);

  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  
  fread(&ib, si, 1, ptrFile); /* repeat */
  if (ib != 0)
  {
    printf("readData: repeat is not supported.\n");
  }
  
  // Variable type
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile);
    if (ib == 1) sizer = 4;
  }

  /* Read dump */
  if (sizer == 4 && dataPacking == 0) // block
  {
    float* buf = new float[npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(float), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i, n+1) = buf[i];
    }
    delete [] buf;
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }
    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = buf2[i];
    }
    delete [] buf2;
  }
  else if (sizer == 8 && dataPacking == 0) // block
  {
    double* buf = new double[npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(double), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i, n+1) = buf[i];
    } 
    delete [] buf;
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }
    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = buf2[i];
    }
    delete [] buf2;
  }
  else if (sizer == 4 && dataPacking == 1) //point
  {
    float* buf = new float[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = buf[i];
    }
    delete [] buf;
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }
    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = buf2[i];
    }
    delete [] buf2;
  }
  else if (sizer == 8 && dataPacking == 1) //point
  {
    double* buf = new double[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(double), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = buf[i];
    }
    delete [] buf;
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }
    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = buf2[i];
    }
    delete [] buf2;
  }

  return 0;
}

//=============================================================================
/* 
   readData pour les zones structurees (endian conversion)
   IN: ptrFile: ptr on file being read
   IN: ni, nj, nk: size of zone.
   IN: dataPacking: 0 (block), 1 (point)
   OUT: f: field read. Must be already dimensioned.
 */
//=============================================================================
E_Int K_IO::GenIO::readData75CE(FILE* ptrFile, E_Int ni, E_Int nj, E_Int nk, 
                                E_Int dataPacking,
                                FldArrayF& f)
{
  float a;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = f.getNfld();
  E_Int si = sizeof(int);

  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  
  fread(&ib, si, 1, ptrFile); /* repeat */
  ib = IBE(ib);
  if (ib != 0)
  {
    printf("readData: repeat is not supported.\n");
  }
  
  // variable type
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib == 1) sizer = 4;
  }
  
  /* Read dump */
  E_Int npts = ni*nj*nk;
  if (sizer == 4 && dataPacking == 0) // block
  {
    float* buf = new float[npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(float), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i, n+1) = FBE(buf[i]);
    }
    delete [] buf;
  }
  else if (sizer == 8 && dataPacking == 0) // block
  {
    double* buf = new double[npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(E_Float), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i, n+1) = DBE(buf[i]);
    }
    delete [] buf;
  }
  else if (sizer == 4 && dataPacking == 1) // point
  {
    float* buf = new float[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = FBE(buf[i]);
    }
    delete [] buf;
  }
  else if (sizer == 8 && dataPacking == 1) // point
  {
    double* buf = new double[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(E_Float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = DBE(buf[i]);
    }
    delete [] buf;
  }

  return 0;
}

//=============================================================================
/* 
   readData pour les zones non-structurees (endian conversion)
   IN: ptrFile: ptr on file being read
   IN: dataPacking: 0 (block), 1 (point)
   OUT: f: field read. Must be already dimensioned.
        c: connectivity. Must be already dimensioned.
 */
//=============================================================================
E_Int K_IO::GenIO::readData75CE(FILE* ptrFile, 
                                E_Int dataPacking,
                                FldArrayF& f, FldArrayI& c)
{
  float a;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = f.getNfld();
  E_Int npts = f.getSize();
  E_Int nelts = c.getSize();
  E_Int eltType = c.getNfld();
                           
  E_Int si = sizeof(int);

  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  
  fread(&ib, si, 1, ptrFile); /* repeat */
  ib = IBE(ib);
  if (ib != 0)
  {
    printf("readData: repeat is not supported.\n");
  }
  
  // variable type
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib == 1) sizer = 4;
  }
  
  /* Read dump */
  if (sizer == 4 && dataPacking == 0) // block
  {
    float* buf = new float[npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(float), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i, n+1) = buf[i];
    }
    delete [] buf;

    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }

    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = IBE(buf2[i]);
    }
    delete [] buf2;
  }
  else if (sizer == 8 && dataPacking == 0) // block
  {
    double* buf = new double[npts];
    for (n = 0; n < nfield; n++)
    {
      fread(buf, sizeof(E_Float), npts, ptrFile);
      for (i = 0; i < npts; i++) f(i, n+1) = DBE(buf[i]);
    }
    delete [] buf;

    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }

    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = IBE(buf2[i]);
    }
    delete [] buf2;
  }
  else if (sizer == 4 && dataPacking == 1)
  {
    float* buf = new float[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = FBE(buf[i]);
    }
    delete [] buf;
    
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }

    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = IBE(buf2[i]);
    }
    delete [] buf2;
  }
  else if (sizer == 8 && dataPacking == 1)
  {
    double* buf = new double[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(E_Float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) f(n, i+1) = DBE(buf[i]);
    }
    delete [] buf;

    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      printf("readData: repeat is not supported.\n");
    }

    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      for (i = 0; i < eltType; i++) c(n, i+1) = IBE(buf2[i]);
    }
    delete [] buf2;
  }
  return 0;
}

//=============================================================================
/*
  This routine enables binary tecplot format of structured and unstructured
  grids.
  IN: file: file name
  IN: varString: strings of vars
  IN: structField: field defined on structured grids
  IN: ni, nj, nk: dimension of structured grids
  IN: unstructField: field defined on unstructured grids.
  IN: connect: connectivity of unstructured grids.
  IN: eltType: element type : 
  1 (BAR), 2 (TRI), 3 (QUAD), 4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA)
  This routine is v75, format point.
  return 1 if failed
*/
//=============================================================================
E_Int K_IO::GenIO::tecwrite75(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector <FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  FILE* ptrFile;
  float a;         // must be basic C types !!
  int   ib;
  char cop[BUFSIZE+1];
  char titre[] = "field";
  FldArrayF buffer;
  E_Int no, nol, p, i, j, n, nt, nv, nf;

  // Constants
  E_Int si = sizeof(int);
  E_Int sf = sizeof(float);

  // Number of variables to write
  E_Int n_vartot = 0;
  if (structField.size() > 0)
    n_vartot = structField[0]->getNfld();
  else if (unstructField.size() > 0)
    n_vartot = unstructField[0]->getNfld();
  else
  {
    printf("Warning: tecwrite: no field given.\n");
    return 1;
  }

  // Check if all elements have the same number of fields
  E_Int structFieldSize = structField.size();
  for (E_Int n = 0; n < structFieldSize; n++)
  {
    E_Int nvar = structField[n]->getNfld();
    if (nvar != n_vartot)
    {
      printf("Warning: tecwrite: number of variables differs for structured field: %d.\n", n+1);
      return 1;
    }
  }
  
  E_Int unstructFieldSize = unstructField.size();
  for (E_Int n = 0; n < unstructFieldSize; n++)
  {
    E_Int nvar = unstructField[n]->getNfld();
    if (nvar != n_vartot)
    {
      printf("Warning: tecwrite: number of variables differs for unstructured field: %d.\n",n+1);
      return 1;
    }
  }
  
  // Start of write
  ptrFile = fopen(file, "wb");
  if (ptrFile == NULL)
  {
    printf("Warning: tecwrite: cannot open file %s.\n", file);
    return 1;
  }

  // Version number
  char version[20];
  strcpy(version,"#!TDV75 ");
  fwrite(version, sizeof(char), 8, ptrFile);
    
  // Magic number
  ib = 1;
  fwrite(&ib, si, 1, ptrFile);
  
  // Title
  ib = titre[0];
  i = 0;
  while (titre[i] != '\0')
  {
    fwrite(&ib, si, 1, ptrFile);
    i++;
    ib = titre[i];
  }
  ib = 0;
  fwrite(&ib, si, 1, ptrFile);
    
  // Number of variables
  ib = n_vartot;
  fwrite(&ib, si, 1, ptrFile);
  
  // Variables
  i = 0;
  j = varString[i];
  while (j != '\0')
  {
    if (j != ',')
    {
      if (j != '"')
      {
        ib = j;
        fwrite(&ib, si, 1, ptrFile);	
      }
    }
    else
    {
      ib = 0;
      fwrite(&ib, si, 1, ptrFile);
    }
    i++;
    j = varString[i];
  }
  ib = 0;
  fwrite(&ib, si, 1, ptrFile);
  
  no = 0;  
  while (no < structFieldSize + unstructFieldSize)
  {
    // Write zone header

    // Zone name
    strcpy(cop, zoneNames[no]);
    
    a = 299.; // zone marker
    fwrite(&a, sf, 1, ptrFile);
      
    ib = cop[0];
    i = 0;
    while (cop[i] != '\0')
    {
      fwrite(&ib, si, 1, ptrFile);
      i++;
      ib = cop[i];
    }
    ib = 0; /* null terminated */
    fwrite(&ib, si, 1, ptrFile);
    if (no < structFieldSize)
    {
      ib = 1; /* point */
      fwrite(&ib, si, 1, ptrFile);	
      ib = -1; /* color auto */
      fwrite(&ib, si, 1, ptrFile);
      ib = ni[no]; /* imax,jmax,kmax */
      fwrite(&ib, si, 1, ptrFile);
      ib = nj[no];
      fwrite(&ib, si, 1, ptrFile);	
      ib = nk[no];
      fwrite(&ib, si, 1, ptrFile);
    }
    else
    {
      nol = no - structFieldSize;
      ib = 3; /* fepoint */
      fwrite(&ib, si, 1, ptrFile);
      ib = -1; /* color auto */
      fwrite(&ib, si, 1, ptrFile); 
      ib = unstructField[nol]->getSize(); /* Number of points */
      fwrite(&ib, si, 1, ptrFile);
      ib = connect[nol]->getSize(); /* Number of elts */
      fwrite(&ib, si, 1, ptrFile);	
      /* Type of elts */
      switch (eltType[nol])
      {
        case 1: // BAR
          ib = 4;
          break;
        case 2:  // TRI
          ib = 0;
          break;
        case 3: // QUAD
          ib = 1;
          break;
        case 4: // TETRA
          ib = 2;
          break;
        case 5: // PYRA - FIX as HEXA
          ib = 3; // Dans ce cas, on trace des hexa degeneres
          break;
        case 6: // PENTA - FIX as HEXA
          ib = 3; // Dans ce cas, on trace des hexa degeneres
          break;
        case 7: // HEXA
          ib = 3;
          break;
        default:
          printf("Warning: tecwrite: wrong element type.\n");
          return 1;
      }
      fwrite(&ib, si, 1, ptrFile);
    }
    
    no++;
  }
  
  /*-------END OF HEADER ZONE---------*/
  a = 357.0; /* secret number */
  fwrite(&a, sf, 1, ptrFile);
  /*-------START OF DATA ZONE---------*/

  /* Structured zones */
  no = 0;
  while (no < structFieldSize)
  {
    FldArrayF& f = *structField[no];

    a = 299.;
    fwrite(&a, sf, 1, ptrFile);
    ib = 0; /* repeat 0 */
    fwrite(&ib, si, 1, ptrFile);
    
    // format of data
    for (i = 0; i < n_vartot; i++)
    {
#ifdef E_DOUBLEREAL
      ib = 2;
#else
      ib = 1;
#endif
      fwrite(&ib, si, 1, ptrFile);
    }

    E_Int ni1 = ni[no];
    E_Int nj1 = nj[no];
    E_Int nk1 = nk[no];
    
    // Buffer preparation
    buffer.resize(ni1*nj1*nk1*n_vartot);
    E_Float* bufferp = buffer.begin();

    for (n = 0; n < ni1*nj1*nk1; n++)
    {
      p = n * n_vartot;
      for (nf = 1; nf <= n_vartot; nf++)
        bufferp[p+nf-1] = f(n, nf);
    }
    
    // buffer write
    fwrite(bufferp, sizeof(E_Float), buffer.getSize(), ptrFile);
    no++;
  }

  /* Unstructured zones */
  no = 0;
  while (no < unstructFieldSize)
  {
    FldArrayF& f = *unstructField[no];
    FldArrayI& c = *connect[no];

    a = 299.;
    fwrite(&a, sf, 1, ptrFile);
    ib = 0; /* repeat 0 */
    fwrite(&ib, si, 1, ptrFile);
    
    // format of data
    for (i = 0; i < n_vartot; i++)
    {
#ifdef E_DOUBLEREAL
      ib = 2;
#else
      ib = 1;
#endif
      fwrite(&ib, si, 1, ptrFile);
    }

    // Buffer preparation
    nt = f.getSize();
    buffer.resize(nt * n_vartot);
    E_Float* bufferp = buffer.begin();

    for (n = 0; n < nt; n++)
    {
      p = n * n_vartot;
      for (nf = 1; nf <= n_vartot; nf++) bufferp[p+nf-1] = f(n, nf);
    }

    // buffer write
    fwrite(bufferp, sizeof(E_Float), buffer.getSize(), ptrFile);
    
    // Connectivity
    nt = c.getSize();
    nv = c.getNfld();
    int* bufferi;
    
    if (eltType[no] == 5) // FIX pour PYRA as HEXA
    {
      bufferi = new int[nt * 8];
      for (n = 0; n < nt; n++)
      {
        p = n * 8;
        bufferi[p  ] = c(n, 1);
        bufferi[p+1] = c(n, 2);
        bufferi[p+2] = c(n, 3);
        bufferi[p+3] = c(n, 4);
        bufferi[p+4] = c(n, 5);
        bufferi[p+5] = c(n, 5);
        bufferi[p+6] = c(n, 5);
        bufferi[p+7] = c(n, 5);
      }
      nv = 8;
    }
    else if (eltType[no] == 6) // FIX pour PENTA as HEXA
    {
      bufferi = new int[nt * 8];
      for (n = 0; n < nt; n++)
      {
        p = n * 8;
        bufferi[p  ] = c(n, 1);
        bufferi[p+1] = c(n, 2);
        bufferi[p+2] = c(n, 2);
        bufferi[p+3] = c(n, 3);
        bufferi[p+4] = c(n, 4);
        bufferi[p+5] = c(n, 5);
        bufferi[p+6] = c(n, 5);
        bufferi[p+7] = c(n, 6);
      }
      nv = 8;
    }
    else // CAS standard
    {
      bufferi = new int[nt * nv];

      for (n = 0; n < nt; n++)
      {
        p = n * nv;
        for (nf = 1; nf <= nv; nf++) bufferi[p+nf-1] = c(n, nf);
      }
    }

    ib = 0; // no repeat
    fwrite(&ib, si, 1, ptrFile);

    fwrite(bufferi, si, nt*nv, ptrFile);
    delete [] bufferi;

    no++;
  }

  fclose(ptrFile);
  return 0;
}
