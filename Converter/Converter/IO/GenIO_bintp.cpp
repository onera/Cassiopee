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

// Binary tecplot file support

# include "GenIO.h"
# include <stdio.h>
# include <string.h>
# include "Array/Array.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/*
  Tecstat.
  Read nfield, varString and set the _convertEndian variable. 
  IN: file: file name
  OUT: varString: filled variables. Must be allocated by the calling 
  routine.
  OUT: nfield: number of variables in file.

  Return 1 if file doesnt exist.
  Return 0 if ok.
  This routine is v75 to v112 compatible.
*/
//=============================================================================
E_Int K_IO::GenIO::tecstat(char* file, char* varString,
                           E_Int& nfield)
{
  FILE* ptrFile;
  char version[9];

  // Check endianess
  E_Int ret = tecCheckEndian(file);
  if (ret == -1)
  {
    printf("Warning: tecstat: can not open file %s.\n", file);
    return 1;
  }
  else if (ret == 1) _convertEndian = true;
  else _convertEndian = false;
  
  /* Opening */
  ptrFile = fopen(file, "rb");
  
  /* Read file header */
  if (_convertEndian == false)
    ret = readHeader(ptrFile, nfield, varString, version);
  else
    ret = readHeaderCE(ptrFile, nfield, varString, version);
  if (ret == 0) return 1;

  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* 
   Read binary tecplot file header.
   IN: ptrFile: opened file ptr
   OUT: nfield: the number of field in the file
   OUT: varString: filled varString. Must be allocated by the calling 
   routine.
   OUT: version: format of file.
   This routine is v75 to v112 compatible.
   Retourne 1 (success), 0 (failed).
*/
//=============================================================================
E_Int K_IO::GenIO::readHeader(FILE *ptrFile, E_Int& nfield, char*& varString, 
                              char* version)
{
  E_Int i, j;
  int ib;

  /* Constants */
  E_Int si = sizeof(int);
  
  /* Version */
  char c[9];
  fread(c, sizeof(char), 8, ptrFile); // version
  c[8] = '\0';
  strcpy(version, c);
  E_Int nversion = numeralVersion(version);
  if (nversion == -1) return 0;

  /* endian check 1 */
  fread(&ib, si, 1, ptrFile);
  
  /* FileType (FULL, ...): ajoute depuis la version 111 */
  if (nversion >= 111) fread(&ib, si, 1, ptrFile);

  /* Title */
  i = 0; ib = 1;
  while (ib != 0)
  {
    fread(&ib, si, 1, ptrFile);
    i++;
  }
  
  /* Number of variables */
  fread(&ib, si, 1, ptrFile);
  nfield = ib;
  varString = new char [nfield*K_ARRAY::VARNAMELENGTH];

  /* Variables name */
  varString[0] = '\0';
  j = 0; i = 0;
  while (j < nfield)
  {
    ib = 1;
    while (ib != 0)
    {
      fread(&ib, si, 1, ptrFile);
      if (ib != 0) varString[i] = char(ib);
      else varString[i] = ','; 
      i++;
    }
    j++;
  }
  varString[i-1] = '\0';
  return 1;
}

//=============================================================================
/* 
   tecread
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

   Notes sur les versions:
   Header: depend de la version du fichier. Est le meme de v75 a v108. 
   Un champ fileType en plus depuis v112. 
   La lecture du header des zones depend de la version du fichier. On utilise
   la version 108 pour les versions futures.
*/
//=============================================================================
E_Int K_IO::GenIO::tecread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connectivity,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  FILE* ptrFile;
  E_Int error, nfield, zone, no, zoneStruct, zoneUnstruct;
  char version[9];
  E_Int ni1, nj1, nk1, dim;
  E_Int et = -1;
  vector<E_Int> etl; vector<E_Int> numFacesl; 
  vector<E_Int> numFaceNodesl; vector<E_Int> neltsl;
  vector<E_Int> numBoundaryFacesl; vector<E_Int> numBoundaryConnectionsl;
  E_Int npts, nelts, numFaces, numFaceNodes;
  E_Int numBoundaryFaces, numBoundaryConnections;

  // Check endianess
  E_Int ret = tecCheckEndian(file);
  if (ret == -1)
  {
    printf("Warning: tecread: can not open file %s.\n", file);
    return 1;
  }
  else if (ret == 1) _convertEndian = true;
  else _convertEndian = false;

  /* Opening */
  ptrFile = fopen(file, "rb");

  /* Read file header */
  if (_convertEndian == false)
    ret = readHeader(ptrFile, nfield, varString, version);
  else
    ret = readHeaderCE(ptrFile, nfield, varString, version);
  if (ret == 0) return 1;

  /* Verification des Versions */
  E_Int vers = numeralVersion(version);
  if (vers == 75 || vers == 108 || vers == 112) { ; /* supported */ }
  else if (vers < 75) 
  { printf("Warning: tecread: the file version %d is not really supported. Trying to read with %d. ", vers, 75); vers = 75; }
  else if (vers < 100)
  { printf("Warning: tecread: the file version %d is not really supported. Trying to read with %d. ", vers, 75); vers = 75; }
  else if (vers < 108)
  { printf("Warning: tecread: the file version %d is not really supported. Trying to read with %d. ", vers, 108); vers = 108; }
  else if (vers < 111)
  { printf("Warning: tecread: the file version %d is not really supported. Trying to read with %d. ", vers, 108); vers = 108; }
  else if (vers > 112)
  { printf("Warning: tecread: the file version %d is not really supported. Trying to read with %d. ", vers, 112); vers = 112; }

  /* Local vector for structured and unstructured zones names */
  vector<char*> structZoneNames, unstructZoneNames;

  /* Look for zones */
  vector<FldArrayF*> geom;
  vector<E_Int> loc;
  E_Int dataPacking=1;
  E_Int strand=0; E_Float time=0.;
  no = 0;
  error = 0;

  FldArrayF* f; FldArrayI* c;

  while (error == 0)
  {
    et = -1;
    char* zoneName = new char[BUFSIZE+1];
    switch (vers)
    {
      case 75:
        if (_convertEndian == false)
          error = readZoneHeader75(ptrFile, dim, ni1, nj1, nk1, 
                                   npts, nelts, et, zoneName, dataPacking,
                                   geom);
        else
          error = readZoneHeader75CE(ptrFile, dim, ni1, nj1, nk1, 
                                     npts, nelts, et, zoneName, dataPacking,
                                     geom);
        if (error == 0) 
        {
          if (et == -1) structZoneNames.push_back(zoneName);
          else unstructZoneNames.push_back(zoneName);
        }
        break;
        
      case 108:
      case 111:
      case 112:
        if (_convertEndian == false)
          error = readZoneHeader108(vers, nfield, ptrFile, dim, 
                                    ni1, nj1, nk1, 
                                    npts, nelts, numFaces, numFaceNodes,
                                    numBoundaryFaces, numBoundaryConnections,
                                    et, zoneName, 
                                    dataPacking, strand, time, 
                                    loc, geom);
        else
          error = readZoneHeader108CE(vers, nfield, ptrFile, dim, 
                                      ni1, nj1, nk1, 
                                      npts, nelts, numFaces, numFaceNodes,
                                      numBoundaryFaces, numBoundaryConnections,
                                      et, zoneName, 
                                      dataPacking, strand, time, 
                                      loc, geom);
        if (error == 0) 
        {
          if (et == -1) structZoneNames.push_back(zoneName);
          else 
          {
            unstructZoneNames.push_back(zoneName);
            numFacesl.push_back(numFaces);
            numFaceNodesl.push_back(numFaceNodes);
            numBoundaryFacesl.push_back(numBoundaryFaces);
            numBoundaryConnectionsl.push_back(numBoundaryConnections);
            neltsl.push_back(nelts);
          }
        }
        break;
    
      default:
        error = 1;
    }
    if (error != 0) break;
    
    // Concatenation of structured and unstructed zones names lists
    zoneNames = structZoneNames;
    zoneNames.insert(zoneNames.end(),unstructZoneNames.begin(),unstructZoneNames.end());

    switch (et)
    {
      case 1: // BAR
        f = new FldArrayF(npts, nfield);
        c = new FldArrayI(nelts,2);
        unstructField.push_back(f);
        connectivity.push_back(c);
        eltType.push_back(1); // "BAR"
        break;
      case 2: // TRI
        f = new FldArrayF(npts, nfield);
        c = new FldArrayI(nelts,3);
        unstructField.push_back(f);
        connectivity.push_back(c);
        eltType.push_back(2); // "TRI"
        break;
      case 3: // QUAD
        f = new FldArrayF(npts, nfield);
        c = new FldArrayI(nelts,4);
        unstructField.push_back(f);
        connectivity.push_back(c);
        eltType.push_back(3); // "QUAD"
        break;
      case 4: // TETRA
        f = new FldArrayF(npts, nfield);
        c = new FldArrayI(nelts,4);
        unstructField.push_back(f);
        connectivity.push_back(c);
        eltType.push_back(4); // "TETRA"
        break;
      case 7: // HEXA 
        f = new FldArrayF(npts, nfield);
        c = new FldArrayI(nelts,8);
        unstructField.push_back(f);
        connectivity.push_back(c);
        eltType.push_back(7); // "HEXA"
        break;
      case 8: // NGON 
        f = new FldArrayF(npts, nfield);
        c = new FldArrayI(1); // ne peut pas etre dimensionne
        unstructField.push_back(f);
        connectivity.push_back(c);
        eltType.push_back(8); // "NGON"
        break;   

      default: // structure
        ni.push_back(ni1);
        nj.push_back(nj1);
        nk.push_back(nk1);
        f = new FldArrayF(ni1*nj1*nk1, nfield);
        structField.push_back(f);
        break;
    }
    etl.push_back(et);
    no++;
  }

  // Data section
  zone = 0;
  zoneStruct = 0;
  zoneUnstruct = 0;

  while (zone < no)
  {
    if (etl[zone] == -1) // structured
    {
      FldArrayF& f = *structField[zoneStruct];
      switch (vers)
      {
        case 75:
          if (_convertEndian == false)
            readData75(ptrFile, 
                       ni[zoneStruct], nj[zoneStruct], nk[zoneStruct],
                       dataPacking, f);
          else
            readData75CE(ptrFile, 
                         ni[zoneStruct], nj[zoneStruct], nk[zoneStruct], 
                         dataPacking, f);
          break;
          
        case 108:
        case 111:
        case 112:
          if (_convertEndian == false)
            readData108(ptrFile, 
                        ni[zoneStruct], nj[zoneStruct], nk[zoneStruct], 
                        dataPacking, f);
          else
            readData108CE(ptrFile, 
                          ni[zoneStruct], nj[zoneStruct], nk[zoneStruct], 
                          dataPacking, f);
          break;
          
        default:;
          break;
      }
      zoneStruct++;
    }
    else
    { 
      // unstructured
      FldArrayF& f1 = *unstructField[zoneUnstruct];
      FldArrayI& c1 = *connectivity[zoneUnstruct];
      switch (vers)
      {
        case 75:
          if (_convertEndian == false)
            readData75(ptrFile, dataPacking, f1, c1);
          else
            readData75CE(ptrFile, dataPacking, f1, c1);
          break;
          
        case 108:
        case 111:
        case 112:
          if (_convertEndian == false)
            readData108(ptrFile, dataPacking, etl[zone], 
                        numFacesl[zoneUnstruct], numFaceNodesl[zoneUnstruct],
                        numBoundaryFacesl[zoneUnstruct],  
                        numBoundaryConnectionsl[zoneUnstruct],
                        neltsl[zoneUnstruct], f1, c1);
          else
            readData108CE(ptrFile, dataPacking, etl[zone],
                          numFacesl[zoneUnstruct], numFaceNodesl[zoneUnstruct],
                          numBoundaryFacesl[zoneUnstruct],  
                          numBoundaryConnectionsl[zoneUnstruct],
                          neltsl[zoneUnstruct], f1, c1);
          break;
          
        default:;
          break;
      }
      zoneUnstruct++;
    }
    zone++;
  }
  fclose(ptrFile);

  // Add geometries if any
  E_Int geomSize = geom.size();
  if (geomSize != 0)
  {
    if (varString[0] == '\0') strcpy(varString, "x,y,z");
    if (nfield > 3)
    {
      for (E_Int i = 0; i < geomSize; i++)
      {
        geom[i]->reAlloc(geom[i]->getSize(), nfield);
        for (E_Int ind = 0; ind < geom[i]->getSize(); ind++)
          for (E_Int nfld = 4; nfld <= nfield; nfld++)
            (*geom[i])(ind, nfld) = 0.;
      }
    }
    if (nfield < 3 && nfield > 0)
    {
      printf("Warning: tecread: the number of variables for zones can not store a geometry. 3 At least 3 variables are needed. Geometry not read.\n");
      for (E_Int i = 0; i < geomSize; i++) delete geom[i];
      geom.clear();
    }
    else
    {
      vector<char*> vars;
      K_ARRAY::extractVars(varString, vars);
      if ( strcmp(vars[0], "x") !=0 && strcmp(vars[0], "CoordinateX") !=0)
        printf(
          "Warning: tecread: first geometry variable is set to: %s\n",
          vars[0]);
      if ( strcmp(vars[1], "y") !=0 && strcmp(vars[1], "CoordinateY") !=0)
        printf(
          "Warning: tecread: second geometry variable is set to: %s\n",
          vars[1]);
      if ( strcmp(vars[2], "z") !=0 && strcmp(vars[2], "CoordinateZ") !=0)
        printf(
          "Warning: tecread: third geometry variable is set to: %s\n",
          vars[2]);
      E_Int varsSize = vars.size();
      for (E_Int i = 0; i < varsSize; i++)
        delete [] vars[i];
    }
  }

  for (E_Int i = 0; i < geomSize; i++)
  {
    ni.push_back(geom[i]->getSize());
    nj.push_back(1);
    nk.push_back(1);
    structField.push_back(geom[i]);
  }

  return 0;
}

//=============================================================================
/*
  This routine enables binary tecplot of field.
*/
//=============================================================================
E_Int K_IO::GenIO::tecwrite(char* file, char* dataFmt, char* varString,
                            E_Int ni, E_Int nj, E_Int nk,
                            const FldArrayF& coord,
                            const FldArrayF& field,
			    vector<char*>& zoneNames)
{
  assert(coord.getSize() == field.getSize());
  vector<E_Int> vni; vector<E_Int> vnj; vector<E_Int> vnk;
  vni.push_back(ni); vnj.push_back(nj); vnk.push_back(nk);
  vector<FldArrayF*> vfield;
  E_Int np = coord.getSize();
  E_Int nfc = coord.getNfld();
  E_Int nff = field.getNfld();
  FldArrayF f(np, nfc+nff);
  for (E_Int j = 1; j <= nfc; j++)
  {
    E_Float* fp = f.begin(j);
    E_Float* coordp = (E_Float*)coord.begin(j);
    for (E_Int i = 0; i < np; i++) fp[i] = coordp[i];
  }
  for (E_Int j = 1; j <= nff; j++)
  {
    E_Float* fp = f.begin(j+nfc);
    E_Float* fieldp = (E_Float*)field.begin(j);
    for (E_Int i = 0; i < np; i++) fp[i] = fieldp[i];
  }
  vfield.push_back(&f);
  vector<FldArrayF*> dummy; vector<FldArrayI*> dummy2; vector<E_Int> dummy3;
  return tecwrite108(file, dataFmt, varString, vni, vnj, vnk, vfield, 
                     dummy, dummy2, dummy3, zoneNames); 
}

//=============================================================================
/*
  This routine enables binary tecplot format of structured grids.
*/
//=============================================================================
E_Int K_IO::GenIO::tecwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector <FldArrayF*>& structField,
  vector<char*>& zoneNames)
{
  vector<FldArrayF*> dummy;
  vector<FldArrayI*> dummy2;
  vector<E_Int> dummy3;
  return tecwrite108(file, dataFmt, varString, ni, nj, nk, structField, 
                     dummy, dummy2, dummy3, zoneNames); 
}

//=============================================================================
/*
  This routine enables binary tecplot format of unstructured grids.
*/
//=============================================================================
E_Int K_IO::GenIO::tecwrite(char* file, char* dataFmt, char* varString,
			    vector<FldArrayF*>& unstructField, 
			    vector<FldArrayI*>& connect,
			    vector<E_Int>& eltType,
                            vector<char*>& zoneNames) 
{
  vector<FldArrayF*> dummy;
  vector<E_Int> d1; vector<E_Int> d2; vector<E_Int> d3;  
  return tecwrite108(file, dataFmt, varString, d1, d2, d3, dummy, 
                     unstructField, connect, eltType, zoneNames);
}

//=============================================================================
/*
  This routine enables binary tecplot format.
*/
//=============================================================================
E_Int K_IO::GenIO::tecwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames)
{
  return tecwrite108(file, dataFmt, varString, ni, nj, nk, structField,  
                     unstructField, connect, eltType, zoneNames);
}
