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

// Binary tecplot v108 file support

//#define SENTINELLE *(double*)&0xffefffffffffffff
#define SENTINELLE -1.79769e+308

# include "GenIO.h"
# include <stdio.h>
# include <string.h>
# include <stdint.h>
# include "Def/DefFunction.h"
# include <math.h>
# include "kcore.h"

using namespace std;
using namespace K_FLD;

// Retourne les begin pour les nodes et les centres
void getBeginFromLoc(vector<E_Int>& loc, vector<E_Int>& beginNodes, 
                     vector<E_Int>& beginCenters)
{
  E_Int nvars = loc.size();
  beginNodes.resize(nvars);
  beginCenters.resize(nvars);
  E_Int nfldNodes, nfldCenters;
  nfldNodes = 1; nfldCenters = 1;
  for (size_t i = 0; i < loc.size(); i++)
  {
    beginNodes[i] = nfldNodes;
    beginCenters[i] = nfldCenters;
    if (loc[i] == 0) { nfldNodes++; }
    else { nfldCenters++; }
  }
}

//=============================================================================
/* Read zone header.
   return 1 if pb else return 0.

   IN: version: tecplot version (108...)
   IN: nfield: number of fields
   IN: ptrFile: ptr on file

   if zone is structured, the returns:
   OUT: dim: 1,2,3
   OUT: ni, nj, nk: size of current zone
   OUT: zoneName: name of zone
   OUT: dataPacking: 0 (block), 1 (point)

   if zone is unstructured, returns:
   OUT: npts: number of points
   OUT: nelts: number of elements
   OUT: numFaces: only for NGONs
   OUT: numFaceNodes: only for NGONS
   OUT: eltType: type of elements
   (1: BAR, 2: TRI, 3: QUAD, 4: TETRA, 7: HEXA).
   OUT: zoneName
   OUT: dataPacking: 0 (block), 1 (point)
   OUT: strand: strand number registered in zone
   OUT: time: time registered in zone
   OUT: rawlocal: 1 si connect est donnee par le parent element
   OUT: loc: for each field, 0 (node), 1(center)
   OUT: geometries (1D fields)
   v108-112 compatible.
*/
//=============================================================================
E_Int K_IO::GenIO::readZoneHeader108(
  E_Int version, E_Int nfield,
  FILE* ptrFile, E_Int& dim,
  E_Int& ni, E_Int& nj, E_Int& nk,
  E_Int& npts, E_Int& nelts,
  E_Int& numFaces, E_Int& numFaceNodes,
  E_Int& numBoundaryFaces, E_Int& numBoundaryConnections,
  E_Int& eltType, E_Int& rawlocal,
  char* zoneName, E_Int& dataPacking,
  E_Int& strand, E_Float& time,
  vector<E_Int>& loc,
  vector<FldArrayF*>& geom)
{
  char dummy[BUFSIZE+1];
  float a;
  double t;
  int ib;
  E_Int i, elt;

  numFaces = -1; numFaceNodes = -1;
  numBoundaryFaces = 0; numBoundaryConnections = 0;

  /* Constants */
  E_Int si = sizeof(int);

  /* Zone marker (299) */
  zonemarker:
  fread(&a, sizeof(float), 1, ptrFile);
  if (K_FUNC::fEqualZero(a - 799.) == true)
  {
    // 799 means undocumented Aux data zone: discarded
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 899.) == true)
  {
    // 899 means Variable Aux data : discarded
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 399.) == true)
  {
    // Geometry
    FldArrayF* field;
    E_Int ret = readGeom108(ptrFile, field);
    if (ret == 1) geom.push_back(field);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 499.) == true)
  {
    // Text : discarded
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 299.) == false) return 1;

  /* Zone name */
  i = 0; ib = 1;
  while (ib != 0 && i < BUFSIZE)
  {
    fread(&ib, si, 1, ptrFile);
    dummy[i] = (char)ib;
    i++;
  }
  if (i == BUFSIZE) dummy[BUFSIZE] = '\0';
  strcpy(zoneName, dummy);

  if (version > 101)
  {
    /* Parent zone: not used. */
    fread(&ib, si, 1, ptrFile);

    /* Strand id. */
    fread(&ib, si, 1, ptrFile);
    strand = ib;

    /* Solution time. */
    fread(&t, sizeof(double), 1, ptrFile);
    time = t;
  }

  /* Zone color: not used. */
  fread(&ib, si, 1, ptrFile);
  
  /* Zone type: checked. */
  fread(&ib, si, 1, ptrFile);
  elt = ib;
  if (elt < 0 || elt > 7)
  {
    printf("Warning: readZoneHeader: those elements are unknown.\n");
    return 1;
  }
  
  /* Data packing: used. Disappear in v112. */
  dataPacking = 0;
  if (version < 112)
  {
    fread(&ib, si, 1, ptrFile);
    dataPacking = ib;
    if (ib != 0 && ib != 1)
    {
      printf("Warning: readZoneHeader: incorrect data packing.\n");
      return 1;
    }
  }

  /* Var location */
  fread(&ib, si, 1, ptrFile);
  loc.resize(nfield);
  for (E_Int i = 0; i < nfield; i++) loc[i] = 0;
  if (ib != 0)
  {
    // var location is specified
    E_Int center;
    for (E_Int i = 0; i < nfield; i++)
    {
      center = 0;
      fread(&ib, si, 1, ptrFile); if (ib == 1) center = 1;
      loc[i] = center;
    }
  }

  if (version > 101)
  {
    /* Raw local 1 to 1 face neighbors supplied: not supported. */
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      printf("Warning: readZoneHeader: raw local faces not supported.\n");
      rawlocal = 1;
    }
  }

  /* Misc user defined faces: not supported. */
  fread(&ib, si, 1, ptrFile);
  if (ib != 0)
  {
    fread(&ib, si, 1, ptrFile); // face mode (local, 1-to-1)
    if (elt != 0) fread(&ib, si, 1, ptrFile);
    printf("Warning: readZoneHeader: user defined faces not supported.\n");
  }

  switch (elt)
  {
    case 0: // Structure
      /* Size of zone: ni, nj, nk */
      fread(&ib, si, 1, ptrFile);
      if (ib > 1) dim = 1;
      ni = ib;
      fread(&ib, si, 1, ptrFile);
      if (ib > 1) dim = 2;
      nj = ib;
      fread(&ib, si, 1, ptrFile);
      if (ib > 1) dim = 3;
      nk = ib;
      eltType = -1;
      break;

    case 1: // BAR (tecplot type)
    case 2: // TRI
    case 3: // QUAD
    case 4: // TETRA
    case 5: // HEXA
      fread(&ib, si, 1, ptrFile);
      npts = ib;
      fread(&ib, si, 1, ptrFile);
      nelts = ib;
      if (elt == 1) eltType = 1; // BAR
      else if (elt == 2) eltType = 2; // TRI
      else if (elt == 3) eltType = 3; // QUAD
      else if (elt == 4) eltType = 4; // TETRA
      else if (elt == 5) eltType = 7; // HEXA
      fread(&ib, si, 1, ptrFile); // cell dim
      fread(&ib, si, 1, ptrFile);
      fread(&ib, si, 1, ptrFile);
      break;

    case 6: // POLYGON
    case 7: // POLYHEDRA
      fread(&ib, si, 1, ptrFile);
      npts = ib;
      fread(&ib, si, 1, ptrFile);
      numFaces = ib;
      fread(&ib, si, 1, ptrFile);
      numFaceNodes = ib;
      fread(&ib, si, 1, ptrFile);
      numBoundaryFaces = ib;
      fread(&ib, si, 1, ptrFile);
      numBoundaryConnections = ib;
      fread(&ib, si, 1, ptrFile);
      nelts = ib;
      eltType = 8; // NGON
      fread(&ib, si, 1, ptrFile); // cell dim
      fread(&ib, si, 1, ptrFile);
      fread(&ib, si, 1, ptrFile);
      break;
  }

  /* Auxiliary data: not supported. */
  fread(&ib, si, 1, ptrFile);
  while (ib == 1)
  {
    printf("Warning: readZoneHeader: auxiliary data is not supported.\n");
    // read name string
    fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);

    // read value format
    fread(&ib, si, 1, ptrFile);

    // read value string
    fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);

    // read next
    fread(&ib, si, 1, ptrFile);
  }

  return 0;
}

//=============================================================================
/* Read zone header.
   return 1 if pb else return 0
   v108-112 compatible.
*/
//=============================================================================
E_Int K_IO::GenIO::readZoneHeader108CE(
  E_Int version, E_Int nfield,
  FILE* ptrFile, E_Int& dim,
  E_Int& ni, E_Int& nj, E_Int& nk,
  E_Int& npts, E_Int& nelts,
  E_Int& numFaces, E_Int& numFaceNodes,
  E_Int& numBoundaryFaces, E_Int& numBoundaryConnections,
  E_Int& eltType, E_Int& rawlocal, 
  char* zoneName, E_Int& dataPacking,
  E_Int& strand, E_Float& time,
  vector<E_Int>& loc,
  vector<FldArrayF*>& geom)
{
  char dummy[BUFSIZE+1];
  float a;
  double t;
  int ib;
  E_Int i, elt;

  numFaces = -1; numFaceNodes = -1;
  numBoundaryFaces = 0; numBoundaryConnections = 0;

  /* Constants */
  E_Int si = sizeof(int);

  /* Zone marker (299) */
  zonemarker:
  fread(&a, sizeof(float), 1, ptrFile); a = FBE(a);

  if (K_FUNC::fEqualZero(a - 799.) == true)
  {
    // 799 means undocumented Aux data zone: discarded
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 899.) == true)
  {
    // 899 means Variable Aux data : discarded
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 399.) == true)
  {
    // Geometry
    FldArrayF* field;
    E_Int ret = readGeom108CE(ptrFile, field);
    if (ret == 1) geom.push_back(field);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 499.) == true)
  {
    // Text : discarded
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    fread(&t, sizeof(double), 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    ib = fread(&ib, si, 1, ptrFile);
    while (ib != 0) fread(&ib, si, 1, ptrFile);
    goto zonemarker;
  }
  else if (K_FUNC::fEqualZero(a - 299.) == false) return 1;

  /* Zone name */
  i = 0; ib = 1;
  while (ib != 0 && i < BUFSIZE)
  {
    fread(&ib, si, 1, ptrFile);
    ib = IBE(ib);
    dummy[i] = (char)ib;
    i++;
  }
  if (i == BUFSIZE) dummy[BUFSIZE] = '\0';
  strcpy(zoneName, dummy);
  
  if (version > 101)
  {
    /* Parent zone: not used */
    fread(&ib, si, 1, ptrFile);
  
    /* Strand id */
    fread(&ib, si, 1, ptrFile);
    strand = IBE(ib);
  
    /* Solution time */
    fread(&t, sizeof(double), 1, ptrFile);
    time = DBE(t);
  }

  /* Zone color: not used */
  fread(&ib, si, 1, ptrFile);
  
  /* Zone type: checked */
  fread(&ib, si, 1, ptrFile);
  elt = IBE(ib);
  if (elt < 0 || elt > 7)
  {
    printf("Warning: readZoneHeader: the element type is unknown.\n");
    return 1;
  }

  /* data packing. Disappear in v112.  */
  dataPacking = 0;
  if (version < 112)
  {
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    dataPacking = ib;
    if (ib != 0 && ib != 1)
    {
      printf("Warning: readZoneHeader: incorrect data packing.\n");
      return 1;
    }
  }

  /* Var location */
  fread(&ib, si, 1, ptrFile); ib = IBE(ib);
  loc.resize(nfield);
  for (E_Int i = 0; i < nfield; i++) loc[i] = 0;
  if (ib != 0)
  {
    // var location is specified
    E_Int center;
    for (E_Int i = 0; i < nfield; i++)
    {
      center = 0;
      fread(&ib, si, 1, ptrFile); ib = IBE(ib); if (ib == 1) center = 1;
      loc[i] = center;
    }
  }

  if (version > 101)
  {
    /* Raw local: not supported */
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      printf("Warning: readZoneHeader: raw local faces not supported.\n");
      rawlocal = 1;
    }
  }

  /* Misc user defined faces: not supported */
  fread(&ib, si, 1, ptrFile); ib = IBE(ib);
  if (ib != 0)
  {
    fread(&ib, si, 1, ptrFile); // face mode (local, 1-to-1)
    if (elt != 0) fread(&ib, si, 1, ptrFile);
    printf("Warning: readZoneHeader: user defined faces not supported.\n");
  }
  
  switch (elt)
  {
    case 0: // Structure
      /* Size of zone: ni, nj, nk */
      fread(&ib, si, 1, ptrFile); ib = IBE(ib);
      if (ib > 1) dim = 1;
      ni = ib;
      fread(&ib, si, 1, ptrFile); ib = IBE(ib);
      if (ib > 1) dim = 2;
      nj = ib;
      fread(&ib, si, 1, ptrFile); ib = IBE(ib);
      if (ib > 1) dim = 3;
      nk = ib;
      break;

    case 1: // BAR (tecplot type)
    case 2: // TRI
    case 3: // QUAD
    case 4: // TETRA
    case 5: // HEXA
      fread(&ib, si, 1, ptrFile);
      npts = IBE(ib);
      fread(&ib, si, 1, ptrFile);
      nelts = IBE(ib);
      if (elt == 1) eltType = 1; // BAR
      else if (elt == 2) eltType = 2; // TRI
      else if (elt == 3) eltType = 3; // QUAD
      else if (elt == 4) eltType = 4; // TETRA
      else if (elt == 5) eltType = 7; // HEXA
      fread(&ib, si, 1, ptrFile); // cell dim
      fread(&ib, si, 1, ptrFile);
      fread(&ib, si, 1, ptrFile);
      break;

    case 6: // POLYGON
    case 7: // POLYHEDRA
      fread(&ib, si, 1, ptrFile);
      npts = IBE(ib);
      fread(&ib, si, 1, ptrFile);
      numFaces = IBE(ib);
      fread(&ib, si, 1, ptrFile);
      numFaceNodes = IBE(ib);
      fread(&ib, si, 1, ptrFile);
      numBoundaryFaces = IBE(ib);
      fread(&ib, si, 1, ptrFile);
      numBoundaryConnections = IBE(ib);
      fread(&ib, si, 1, ptrFile);
      nelts = IBE(ib);
      eltType = 8; // NGON
      fread(&ib, si, 1, ptrFile); // cell dim
      fread(&ib, si, 1, ptrFile);
      fread(&ib, si, 1, ptrFile);
      break;
  }

  /* Auxiliary data: not supported. */
  fread(&ib, si, 1, ptrFile); ib = IBE(ib);
  if (ib == 1)
  {
    printf("Warning: readZoneHeader: auxiliary data is not supported.\n");
    // read name string
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    while (ib != 0) { fread(&ib, si, 1, ptrFile); ib = IBE(ib); }

    // read value format
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);

    // read value string
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    while (ib != 0) { fread(&ib, si, 1, ptrFile); ib = IBE(ib); }

    // read next
    fread(&ib, si, 1, ptrFile);
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
E_Int K_IO::GenIO::readData108(E_Int version, FILE* ptrFile, 
                               E_Int ni, E_Int nj, E_Int nk,
                               E_Int dataPacking, vector<E_Int>& loc,
                               FldArrayF* f, FldArrayF* fc)
{
  float a;
  int ib;
  double t;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = loc.size();
  if (nfield == 0) nfield = f->getNfld();
  E_Int si = sizeof(int);
  E_Int npts = ni*nj*nk;
  E_Int nelts = K_FUNC::E_max((ni-1),1)*K_FUNC::E_max((nj-1),1)*K_FUNC::E_max((nk-1),1);

  vector<E_Int> beginNodes; vector<E_Int> beginCenters;
  getBeginFromLoc(loc, beginNodes, beginCenters);
    
  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  if (K_FUNC::fEqualZero(a - 299.) == false) return 1;
  
  /* Type des variables */
  vector<int> varType(nfield);
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile); // variables type
    if (ib == 1) sizer = 4;
    varType[i] = ib; // 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
  }

  /* Passive variables */
  E_Int* passive = new E_Int[nfield];
  for (E_Int i = 0; i < nfield; i++) passive[i] = 0;

  if (version > 101)
  {
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      for (E_Int i = 0; i < nfield; i++)
      { 
        fread(&ib, si, 1, ptrFile);
        passive[i] = ib;
      }
    }
  }

  /* Sharing variables */
  fread(&ib, si, 1, ptrFile);
  if (ib != 0)
  {
    for (E_Int i = 0; i < nfield; i++) fread(&ib, si, 1, ptrFile);
    printf("Warning: this file has sharing variables. Not supported.\n");
  }

  /* Sharing zone */
  fread(&ib, si, 1, ptrFile);
  if (ib != -1)
  {
    printf("Warning: this file has sharing connectivity. Not supported.\n");
  }

  if (version > 101)
  {
    /* Min-Max since no sharing and no passive. */
    for (n = 0; n < nfield; n++)
    {
      if (passive[n] == 0)
      {
        fread(&t, sizeof(double), 1, ptrFile);
        fread(&t, sizeof(double), 1, ptrFile);
      }
    }
  }

  /* Read dump */
  E_Int size;
  if (dataPacking == 0) // block all types
  {
    for (n = 0; n < nfield; n++)
    {
      float* buf1=NULL; double* buf2=NULL; int64_t* buf3=NULL;
      int32_t* buf4=NULL; int8_t* buf5=NULL; 
      E_Float* fp;
      if (loc[n] == 1) { size = nelts; fp = fc->begin(beginCenters[n]); }
      else { size = npts; fp = f->begin(beginNodes[n]); }

      if (varType[n] == 1) buf1 = new float[size];
      else if (varType[n] == 2) buf2 = new double[size];
      else if (varType[n] == 3) buf3 = new int64_t [size];
      else if (varType[n] == 4) buf4 = new int32_t [size];
      else if (varType[n] == 5) buf5 = new int8_t [size];
      else printf("Warning: unknow type of variable.\n");

      if (passive[n] == 0)
      {
        if (varType[n] == 1) 
        {  
          fread(buf1, sizeof(float), size, ptrFile);
          for (i = 0; i < size; i++) fp[i] = buf1[i];
          delete [] buf1;
        }    
        else if (varType[n] == 2) 
        {
          fread(buf2, sizeof(double), size, ptrFile);
          for (i = 0; i < size; i++) fp[i] = buf2[i];
          delete [] buf2;
        }
        else if (varType[n] == 3) 
        {
          fread(buf3, sizeof(int64_t), size, ptrFile);
          for (i = 0; i < size; i++) fp[i] = buf3[i];
          delete [] buf3;
        }
        else if (varType[n] == 4) 
        {
          fread(buf4, sizeof(int32_t), size, ptrFile);
          for (i = 0; i < size; i++) fp[i] = buf4[i];
          delete [] buf4;
        }
        else if (varType[n] == 5) 
        {
          fread(buf5, sizeof(int8_t), size, ptrFile);
          for (i = 0; i < size; i++) fp[i] = buf5[i];
          delete [] buf5;
        }
      }
      else
      {
        for (E_Int i = 0; i < size; i++) fp[i] = 0.;
      }
    }    
  }
  else if (dataPacking == 1) // point all types (forcement node)
  {
    FldArrayF& fp = *f;
    
    // rip passives
    int nfieldLoc = 0;
    vector<int> index(nfield); // index of non passive field
    int bsize = 0; 
    vector<int> pos(nfield+1); // pos of non passive field in buffer
    for (i = 0; i < nfield; i++)
    {
      if (passive[i] == 0) 
      {
        if (varType[i] == 1) { pos[nfieldLoc] = bsize; bsize += 4; }
        else if (varType[i] == 2) { pos[nfieldLoc] = bsize; bsize += 8;  }
        else if (varType[i] == 3) { pos[nfieldLoc] = bsize; bsize += 8;  }
        else if (varType[i] == 4) { pos[nfieldLoc] = bsize; bsize += 4;  }
        else if (varType[i] == 5) { pos[nfieldLoc] = bsize; bsize += 1;  }
        index[nfieldLoc] = i; nfieldLoc += 1;
      }
    }
    char* buf = new char[bsize];

    for (n = 0; n < npts; n++)
    {
      for (i = 0; i < nfield; i++) fp(n, i+1) = 0.;
      fread(buf, sizeof(char), bsize, ptrFile);
      for (i = 0; i < nfieldLoc; i++) 
      {
        switch (varType[index[i]])
        {
          case 1:
          {
            float* ptr = (float*)&buf[pos[i]]; 
            fp(n, index[i]+1) = ptr[0];
          }
          break;
          case 2:
          {
            double* ptr = (double*)&buf[pos[i]]; 
            fp(n, index[i]+1) = ptr[0];
          }
          break;
          case 3:
          {
            int64_t* ptr = (int64_t*)&buf[pos[i]]; 
            fp(n, index[i]+1) = ptr[0];
          }
          break;
          case 4:
          {
            int32_t* ptr = (int32_t*)&buf[pos[i]]; 
            fp(n, index[i]+1) = ptr[0];
          }
          break;
          case 5:
          {
            int8_t* ptr = (int8_t*)&buf[pos[i]]; 
            fp(n, index[i]+1) = ptr[0];
          }
          break;
        }
      }    
    }
    delete [] buf;
  }
  else if (dataPacking == 1 && sizer == 4) // point old code
  {
    float* buf = new float[nfield];
    FldArrayF& fp = *f;
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) fp(n, i+1) = buf[i];
    }
    delete [] buf;
  }
  else if (dataPacking == 1 && sizer == 8) // point old code
  {
    FldArrayF& fp = *f;
    double* buf = new double[nfield];
    for (n = 0; n < npts; n++)
    {
      fread(buf, sizeof(E_Float), nfield, ptrFile);
      for (i = 0; i < nfield; i++) fp(n, i+1) = buf[i];
    }
    delete [] buf;
  }

  delete [] passive;
  return 0;
}
//=============================================================================
/*
   readData pour les zones non-structurees.
   IN: ptrFile: ptr on file being read
   IN: dataPacking: 0 (block), 1 (point)
   IN: et: type d'element (utilise seult pour les NGON)
   IN: numFaces: nbre de faces (uniquement pour les NGON)
   IN: numFaceNodes: nbre de noeuds par face pour toutes les faces (NGON)
   IN: ne: nbre d'elements (NGON)
   OUT: f: field read. Must be already dimensioned.
        c: connectivity read. Must be already dimensioned (sauf pour les NGONs)
 */
//=============================================================================
E_Int K_IO::GenIO::readData108(
  E_Int version,
  FILE* ptrFile,
  E_Int dataPacking, vector<E_Int>& loc, E_Int et,
  E_Int numFaces, E_Int numFaceNodes,
  E_Int numBoundaryFaces, E_Int numBoundaryConnections, 
  E_Int ne, E_Int rawlocal,
  FldArrayF* f, FldArrayI& c, FldArrayF* fc)
{
  float a;
  double t;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = loc.size();
  if (nfield == 0) nfield = f->getNfld(); 
  E_Int npts;
  if (f != NULL) npts = f->getSize();
  else npts = 0;
  E_Int nelts = c.getSize();
  E_Int eltType = c.getNfld();
  E_Int si = sizeof(int);

  vector<E_Int> beginNodes; vector<E_Int> beginCenters;
  getBeginFromLoc(loc, beginNodes, beginCenters);
    
  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  if (K_FUNC::fEqualZero(a - 299.) == false) return 1;

  /* Type des variables */
  vector<int> varType(nfield);
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile); // variables type
    if (ib == 1) sizer = 4;
    varType[i] = ib; // 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
  }

  /* Passive variables */
  E_Int* passive = new E_Int[nfield];
  for (E_Int i = 0; i < nfield; i++) passive[i] = 0;
  
  if (version > 101)
  {
    fread(&ib, si, 1, ptrFile);
    if (ib != 0)
    {
      for (E_Int i = 0; i < nfield; i++)
      {
        fread(&ib, si, 1, ptrFile);
        passive[i] = ib;
      }
    }
  }

  /* Sharing variables */
  fread(&ib, si, 1, ptrFile);
  if (ib != 0)
  {
    for (E_Int i = 0; i < nfield; i++) fread(&ib, si, 1, ptrFile);
    printf("Warning: this file has sharing variables. Not supported.\n");
  }

  /* Sharing connectivity */
  fread(&ib, si, 1, ptrFile);
  if (ib != -1)
  {
    printf("Warning: this file has sharing connectivity. Not supported.\n");
  }

  if (version > 101)
  {
    /* Min-Max for no sharing and no passive. */
    for (n = 0; n < nfield; n++)
    {
      if (passive[n] == 0)
      {
        fread(&t, sizeof(double), 1, ptrFile);
        fread(&t, sizeof(double), 1, ptrFile);
      }
    }
  }

  /* Read dump */
  E_Int size;
  if (dataPacking == 0) // block all types
  {
    for (n = 0; n < nfield; n++)
    {
      float* buf1=NULL; double* buf2=NULL; int64_t* buf3=NULL;
      int32_t* buf4=NULL; int8_t* buf5=NULL;
      
      E_Float* fp;
      if (loc[n] == 1) { size = nelts; fp = fc->begin(beginCenters[n]); }
      else { size = npts; fp = f->begin(beginNodes[n]); }

      if (varType[n] == 1) buf1 = new float[size];
      else if (varType[n] == 2) buf2 = new double[size];
      else if (varType[n] == 3) buf3 = new int64_t [size];
      else if (varType[n] == 4) buf4 = new int32_t [size];
      else if (varType[n] == 5) buf5 = new int8_t [size];
      else printf("Warning: unknown type of variable.\n");

      if (passive[n] == 0)
      {
        if (varType[n] == 1) 
        {
            fread(buf1, sizeof(float), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = buf1[i];
            delete [] buf1;
        }    
        else if (varType[n] == 2) 
        {
            fread(buf2, sizeof(double), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = buf2[i];
            delete [] buf2;
        }
        else if (varType[n] == 3) 
        {
            fread(buf3, sizeof(int64_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = buf3[i];
            delete [] buf3;
        }
        else if (varType[n] == 4) 
        {
            fread(buf4, sizeof(int32_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = buf4[i];
            delete [] buf4;
        }
        else if (varType[n] == 5) 
        {
            fread(buf5, sizeof(int8_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = buf5[i];
            delete [] buf5;
        }
      }
      else
      {
        for (E_Int i = 0; i < size; i++) fp[i] = 0.;
      }
    } 
  }
  else if (dataPacking == 1) // point all types (forcement node)
  {
    FldArrayF& fp = *f;
    
    // rip passives
    int nfieldLoc = 0;
    vector<int> index(nfield); // index of non passive field
    int bsize = 0; 
    vector<int> pos(nfield+1); // pos of non passive field in buffer
    for (i = 0; i < nfield; i++)
    {
        if (passive[i] == 0) 
        {
            if (varType[i] == 1) { pos[nfieldLoc] = bsize; bsize += 4; }
            else if (varType[i] == 2) { pos[nfieldLoc] = bsize; bsize += 8;  }
            else if (varType[i] == 3) { pos[nfieldLoc] = bsize; bsize += 8;  }
            else if (varType[i] == 4) { pos[nfieldLoc] = bsize; bsize += 4;  }
            else if (varType[i] == 5) { pos[nfieldLoc] = bsize; bsize += 1;  }
            index[nfieldLoc] = i; nfieldLoc += 1;
        }
    }
    char* buf = new char[bsize];

    for (n = 0; n < npts; n++)
    {
       for (i = 0; i < nfield; i++) fp(n, i+1) = 0.;
       fread(buf, sizeof(char), bsize, ptrFile);
       for (i = 0; i < nfieldLoc; i++) 
       {
           switch (varType[index[i]])
           {
               case 1:
               {
                float* ptr = (float*)&buf[pos[i]]; 
                fp(n, index[i]+1) = ptr[0];
               }
               break;
               case 2:
               {
                double* ptr = (double*)&buf[pos[i]]; 
                fp(n, index[i]+1) = ptr[0];
               }
               break;
               case 3:
               {
                int64_t* ptr = (int64_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = ptr[0];
               }
               break;
               case 4:
               {
                int32_t* ptr = (int32_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = ptr[0];
               }
               break;
               case 5:
               {
                int8_t* ptr = (int8_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = ptr[0];
               }
               break;
           }
       }    
    }
    delete [] buf;
  }
  else if (dataPacking == 1 && sizer == 4) // point old code
  {
    float* buf = new float[nfield];
    FldArrayF& fp = *f;
    for (n = 0; n < npts; n++)
    {
        fread(buf, sizeof(float), nfield, ptrFile);
        for (i = 0; i < nfield; i++) fp(n, i+1) = buf[i];
    }
    delete [] buf;
  }
  else if (dataPacking == 1 && sizer == 8) // point old code
  {
    FldArrayF& fp = *f;
    double* buf = new double[nfield];
    for (n = 0; n < npts; n++)
    {
        fread(buf, sizeof(E_Float), nfield, ptrFile);
        for (i = 0; i < nfield; i++) fp(n, i+1) = buf[i];
    }
    delete [] buf;
  }

  delete [] passive;
  
  // Connectivity
  if (et != 8) // elements basiques
  {
    int* buf2 = new int[eltType];
    for (n = 0; n < nelts; n++)
    {
      fread(buf2, si, eltType, ptrFile);
      if (version <= 101)
      { for (i = 0; i < eltType; i++) c(n, i+1) = buf2[i]; }
      else // la numerotation a change en version 108!
      { for (i = 0; i < eltType; i++) c(n, i+1) = buf2[i]+1; }
    }
    delete [] buf2;
  }
  else // NGON
  {
    // read face offset
    int* offset = new int[numFaces+1];
    if (numFaceNodes != 2*numFaces) // volumic
    {
      fread(offset, si, numFaces+1, ptrFile);
    }
    else // surfacic
    {
      for (E_Int i = 0; i < numFaces+1; i++) offset[i] = 2*i;
    }
    // faceNodes
    int* faceNodes = new int[numFaceNodes];
    fread(faceNodes, si, numFaceNodes, ptrFile);

    // left elements
    int* left = new int[numFaces];
    fread(left, si, numFaces, ptrFile);

    // right elements
    int* right = new int[numFaces];
    fread(right, si, numFaces, ptrFile);

    // skip boundary data
    if (numBoundaryFaces != 0)
    {
      for (E_Int i = 0; i < numBoundaryFaces; i++)
        fread(&ib, si, 1, ptrFile);
      for (E_Int i = 0; i < numBoundaryConnections; i++)
        fread(&ib, si, 1, ptrFile);
      for (E_Int i = 0; i < numBoundaryConnections; i++)
        fread(&ib, si, 1, ptrFile);
    }

    // Construit la connectivite array
    int* count = new int[numFaces];
    E_Int size = 0;
    for (E_Int i = 0; i < numFaces; i++)
    {
      count[i] = offset[i+1]-offset[i];
      size += count[i]+1;
    }
    // build FN
    FldArrayI faces(size);
    E_Int* ptr = faces.begin();
    E_Int n;
    E_Int p = 0;
    for (E_Int i = 0; i < numFaces; i++)
    {
      n = count[i];
      ptr[0] = n;
      for (E_Int j = 0; j < n; j++)
      {
        ptr[j+1] = faceNodes[p]+1; p++;
      }
      ptr += n+1;
    }
    delete [] count; delete [] offset; delete [] faceNodes;

    // cFE
    FldArrayI cFE(numFaces, 2);
    E_Int* cFE1 = cFE.begin(1); E_Int* cFE2 = cFE.begin(2);
    for (E_Int i = 0; i < numFaces; i++)
    {
      cFE1[i] = left[i]+1; cFE2[i] = right[i]+1;
    }
    delete [] left; delete [] right;

    // Construit la connectivite elts->faces
    FldArrayI cEF;
    K_CONNECT::connectFE2EF(cFE, ne, cEF);

    // Fusionne les connectivites
    c.malloc(faces.getSize()+cEF.getSize()+4);
    E_Int* cp = c.begin();
    cp[0] = numFaces;
    cp[1] = faces.getSize();
    ptr = faces.begin();
    for (E_Int i = 0; i < faces.getSize(); i++)
    {
      cp[i+2] = ptr[i];
    }
    int pt = 2+faces.getSize();
    cp[pt] = ne;
    cp[pt+1] = cEF.getSize();
    ptr = cEF.begin();
    for (E_Int i = 0; i < cEF.getSize(); i++)
    {
      cp[pt+2+i] = ptr[i];
    }
  }
  return 0;
}

//=============================================================================
/*
   readData pour les zones structurees (endian conversion)
   IN: ptrFile: ptr on file being read
   IN: ni, nj, nk: size of zone.
   IN: dataPacking: 0 (block), 1 (point)
   IN: et: type d'element (utilise seult pour les NGON)
   IN: numFaces: nbre de faces (uniquement pour les NGON)
   IN: numFaceNodes: nbre de noeuds par face pour toutes les faces (NGON)
   OUT: f: field read. Must be already dimensioned.
 */
//=============================================================================
E_Int K_IO::GenIO::readData108CE(E_Int version, FILE* ptrFile, 
                                 E_Int ni, E_Int nj, E_Int nk,
                                 E_Int dataPacking, vector<E_Int>& loc,
                                 FldArrayF* f, FldArrayF* fc)
{
  float a;
  double t;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = loc.size();
  if (nfield == 0) nfield = f->getNfld();
  E_Int si = sizeof(int);
  E_Int npts = ni*nj*nk;
  E_Int nelts = K_FUNC::E_max((ni-1),1)*K_FUNC::E_max((nj-1),1)*K_FUNC::E_max((nk-1),1);

  vector<E_Int> beginNodes; vector<E_Int> beginCenters;
  getBeginFromLoc(loc, beginNodes, beginCenters);  
  
  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  a = FBE(a);
  if (K_FUNC::fEqualZero(a - 299.) == false) return 1;

  /* Type des variables */
  vector<int> varType(nfield);
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile); // variables type
    ib = IBE(ib);
    if (ib == 1) sizer = 4;
    varType[i] = ib;
  }

  /* Passive variables */
  E_Int* passive = new E_Int[nfield];
  for (E_Int i = 0; i < nfield; i++) passive[i] = 0;
  if (version > 101)
  {
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      for (E_Int i = 0; i < nfield; i++) 
      {
        fread(&ib, si, 1, ptrFile);
        passive[i] = IBE(ib);
      }
    }
  }

  /* Sharing variables */
  fread(&ib, si, 1, ptrFile); ib = IBE(ib);
  if (ib != 0)
  {
    for (E_Int i = 0; i < nfield; i++) fread(&ib, si, 1, ptrFile);
    printf("Warning: this file has sharing variables. Not supported.\n");
  }

  /* Sharing connectivity  */
  fread(&ib, si, 1, ptrFile); ib = IBE(ib);
  if (ib != -1)
  {
    printf("Warning: this file has sharing connectivity. Not supported.\n");
  }

  if (version > 101)
  {
    /* Min-Max since no sharing and no passive. */
    for (n = 0; n < nfield; n++)
    {
      if (passive[n] == 0)
      {
        fread(&t, sizeof(double), 1, ptrFile);
        fread(&t, sizeof(double), 1, ptrFile);
      }
    }
  }

  /* Read dump */
  E_Int size;
  if (dataPacking == 0) // block all types
  {
    for (n = 0; n < nfield; n++)
    {
      float* buf1=NULL; double* buf2=NULL; int64_t* buf3=NULL;
      int32_t* buf4=NULL; int8_t* buf5=NULL; 
      E_Float* fp;
      if (loc[n] == 1) { size = nelts; fp = fc->begin(beginCenters[n]); }
      else { size = npts; fp = f->begin(beginNodes[n]); }
      if (varType[n] == 1) buf1 = new float[size];
      else if (varType[n] == 2) buf2 = new double[size];
      else if (varType[n] == 3) buf3 = new int64_t [size];
      else if (varType[n] == 4) buf4 = new int32_t [size];
      else if (varType[n] == 5) buf5 = new int8_t [size];
      else printf("Warning: unknown type of variable.\n");

      if (passive[n] == 0)
      {
        if (varType[n] == 1) 
        {
            fread(buf1, sizeof(float), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = FBE(buf1[i]);
            delete [] buf1;
        }    
        else if (varType[n] == 2) 
        {
            fread(buf2, sizeof(double), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = DBE(buf2[i]);
            delete [] buf2;
        }
        else if (varType[n] == 3)
        {
            fread(buf3, sizeof(int64_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = LBE(buf3[i]);
            delete [] buf3;
        }
        else if (varType[n] == 4) 
        {
            fread(buf4, sizeof(int32_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = IBE(buf4[i]);
            delete [] buf4;
        }
        else if (varType[n] == 5) 
        {
            fread(buf5, sizeof(int8_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = SBE(buf5[i]);
            delete [] buf5;
        }
      }
      else
      {
        for (E_Int i = 0; i < size; i++) fp[i] = 0.;
      }
    }    
  }
  else if (dataPacking == 1) // point all types (forcement node)
  {
    FldArrayF& fp = *f;
    
    // rip passives
    int nfieldLoc = 0;
    vector<int> index(nfield); // index of non passive field
    int bsize = 0; 
    vector<int> pos(nfield+1); // pos of non passive field in buffer
    for (i = 0; i < nfield; i++)
    {
        if (passive[i] == 0) 
        {
            if (varType[i] == 1) { pos[nfieldLoc] = bsize; bsize += 4; }
            else if (varType[i] == 2) { pos[nfieldLoc] = bsize; bsize += 8;  }
            else if (varType[i] == 3) { pos[nfieldLoc] = bsize; bsize += 8;  }
            else if (varType[i] == 4) { pos[nfieldLoc] = bsize; bsize += 4;  }
            else if (varType[i] == 5) { pos[nfieldLoc] = bsize; bsize += 1;  }
            index[nfieldLoc] = i; nfieldLoc += 1;
        }
    }
    char* buf = new char[bsize];

    for (n = 0; n < npts; n++)
    {
       for (i = 0; i < nfield; i++) fp(n, i+1) = 0.;
       fread(buf, sizeof(char), bsize, ptrFile);
       for (i = 0; i < nfieldLoc; i++) 
       {
           switch (varType[index[i]])
           {
               case 1:
               {
                float* ptr = (float*)&buf[pos[i]]; 
                fp(n, index[i]+1) = FBE(ptr[0]);
               }
               break;
               case 2:
               {
                double* ptr = (double*)&buf[pos[i]]; 
                fp(n, index[i]+1) = DBE(ptr[0]);
               }
               break;
               case 3:
               {
                int64_t* ptr = (int64_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = LBE(ptr[0]);
               }
               break;
               case 4:
               {
                int32_t* ptr = (int32_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = IBE(ptr[0]);
               }
               break;
               case 5:
               {
                int8_t* ptr = (int8_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = SBE(ptr[0]);
               }
               break;
           }
       }    
    }
    delete [] buf;
  }
  else if (sizer == 4 && dataPacking == 1) // point old code
  {
    vector<float> buf(nfield);
    FldArrayF& fp = *f;
    for (n = 0; n < npts; n++)
    {
        fread(&buf[0], sizeof(float), nfield, ptrFile);
        for (i = 0; i < nfield; i++) fp(n, i+1) = FBE(buf[i]);
    }
  }
  else if (sizer == 8 && dataPacking == 1) // point old code
  {
    vector<double> buf(nfield);
    FldArrayF& fp = *f;
    for (n = 0; n < npts; n++)
    {
        fread(&buf[0], sizeof(E_Float), nfield, ptrFile);
        for (i = 0; i < nfield; i++) fp(n, i+1) = DBE(buf[i]);
    }
  }

  delete [] passive;
  
  return 0;
}
//=============================================================================
/*
   readData pour les zones non-structurees (endian conversion)
   IN: ptrFile : ptr on file being read
   IN: dataPacking : 0 (block), 1 (point)
   IN: et: element type
   OUT: f: field read. Must be already dimensioned.
        c: connectivity read. Must be already dimensioned.
 */
//=============================================================================
E_Int K_IO::GenIO::readData108CE(
  E_Int version, FILE* ptrFile,
  E_Int dataPacking, vector<E_Int>& loc, E_Int et,
  E_Int numFaces, E_Int numFaceNodes,
  E_Int numBoundaryFaces, E_Int numBoundaryConnections,
  E_Int ne, E_Int rawlocal,
  FldArrayF* f, FldArrayI& c, FldArrayF* fc)
{
  float a;
  double t;
  int ib;
  E_Int i, n;
  E_Int sizer = 8;
  E_Int nfield = loc.size();
  if (nfield == 0) nfield = f->getNfld();
  E_Int npts;
  if (f != NULL) npts = f->getSize();
  else npts = 0; 
  E_Int nelts = c.getSize();
  E_Int eltType = c.getNfld();
  E_Int si = sizeof(int);
  
  vector<E_Int> beginNodes; vector<E_Int> beginCenters;
  getBeginFromLoc(loc, beginNodes, beginCenters);
  
  // Read zone separator
  fread(&a, sizeof(float), 1, ptrFile); // 299.
  a = FBE(a);
  if (K_FUNC::fEqualZero(a - 299.) == false) return 1;

  /* Type des variables */
  vector<int> varType(nfield);
  for (i = 0; i < nfield; i++)
  {
    fread(&ib, si, 1, ptrFile); ib = IBE(ib); // variable type
    if (ib == 1) sizer = 4;
    varType[i] = ib;
  }
  
  E_Int* passive = new E_Int[nfield];
  for (E_Int i = 0; i < nfield; i++) passive[i] = 0;
  if (version > 101)
  {
    /* Passive variables */
    fread(&ib, si, 1, ptrFile); ib = IBE(ib);
    if (ib != 0)
    {
      for (E_Int i = 0; i < nfield; i++) 
      {
        fread(&ib, si, 1, ptrFile);
        passive[i] = IBE(ib);
      }
    }
  }

  /* Sharing variables */
  fread(&ib, si, 1, ptrFile); ib = IBE(ib);
  if (ib != 0)
  {
    for (E_Int i = 0; i < nfield; i++) fread(&ib, si, 1, ptrFile);
    printf("Warning: this file has sharing variables. Not supported.\n");
  }

  /* Sharing zone connectivity */
  fread(&ib, si, 1, ptrFile);
  if (ib != -1)
  {
    printf("This file has sharing connectivity. Not supported.\n");
  }

  if (version > 101)
  {
    /* Min-Max since no sharing and no passive. */
    for (n = 0; n < nfield; n++)
    {
      if (passive[n] == 0)
      {
        fread(&t, sizeof(double), 1, ptrFile);
        fread(&t, sizeof(double), 1, ptrFile);
      }
    }
  }

  /* Read dump */
  E_Int size;
  if (dataPacking == 0) // block all types
  {
    for (n = 0; n < nfield; n++)
    {
      float* buf1=NULL; double* buf2=NULL; int64_t* buf3=NULL;
      int32_t* buf4=NULL; int8_t* buf5=NULL;
      E_Float* fp;
      if (loc[n] == 1) { size = nelts; fp = fc->begin(beginCenters[n]); }
      else { size = npts; fp = f->begin(beginNodes[n]); }
      if (varType[n] == 1) buf1 = new float[size];
      else if (varType[n] == 2) buf2 = new double[size];
      else if (varType[n] == 3) buf3 = new int64_t [size];
      else if (varType[n] == 4) buf4 = new int32_t [size];
      else if (varType[n] == 5) buf5 = new int8_t [size];
      else printf("Warning: unknown type of variable.\n");

      if (passive[n] == 0)
      {
        if (varType[n] == 1) 
        {
            fread(buf1, sizeof(float), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = FBE(buf1[i]);
            delete [] buf1;
        }    
        else if (varType[n] == 2) 
        {
            fread(buf2, sizeof(double), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = DBE(buf2[i]);
            delete [] buf2;
        }
        else if (varType[n] == 3)
        {
            fread(buf3, sizeof(int64_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = LBE(buf3[i]);
            delete [] buf3;
        }
        else if (varType[n] == 4) 
        {
            fread(buf4, sizeof(int32_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = IBE(buf4[i]);
            delete [] buf4;
        }
        else if (varType[n] == 5) 
        {
            fread(buf5, sizeof(int8_t), size, ptrFile);
            for (i = 0; i < size; i++) fp[i] = SBE(buf5[i]);
            delete [] buf5;
        }
      }
      else
      {
        for (E_Int i = 0; i < size; i++) fp[i] = 0.;
      }
    }    
  }
  else if (dataPacking == 1) // point all types (forcement node)
  {
    FldArrayF& fp = *f;
    
    // rip passives
    int nfieldLoc = 0;
    vector<int> index(nfield); // index of non passive field
    int bsize = 0; 
    vector<int> pos(nfield+1); // pos of non passive field in buffer
    for (i = 0; i < nfield; i++)
    {
        if (passive[i] == 0) 
        {
            if (varType[i] == 1) { pos[nfieldLoc] = bsize; bsize += 4; }
            else if (varType[i] == 2) { pos[nfieldLoc] = bsize; bsize += 8;  }
            else if (varType[i] == 3) { pos[nfieldLoc] = bsize; bsize += 8;  }
            else if (varType[i] == 4) { pos[nfieldLoc] = bsize; bsize += 4;  }
            else if (varType[i] == 5) { pos[nfieldLoc] = bsize; bsize += 1;  }
            index[nfieldLoc] = i; nfieldLoc += 1;
        }
    }
    char* buf = new char[bsize];

    for (n = 0; n < npts; n++)
    {
       for (i = 0; i < nfield; i++) fp(n, i+1) = 0.;
       fread(buf, sizeof(char), bsize, ptrFile);
       for (i = 0; i < nfieldLoc; i++) 
       {
           switch (varType[index[i]])
           {
               case 1:
               {
                float* ptr = (float*)&buf[pos[i]]; 
                fp(n, index[i]+1) = FBE(ptr[0]);
               }
               break;
               case 2:
               {
                double* ptr = (double*)&buf[pos[i]]; 
                fp(n, index[i]+1) = DBE(ptr[0]);
               }
               break;
               case 3:
               {
                int64_t* ptr = (int64_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = LBE(ptr[0]);
               }
               break;
               case 4:
               {
                int32_t* ptr = (int32_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = IBE(ptr[0]);
               }
               break;
               case 5:
               {
                int8_t* ptr = (int8_t*)&buf[pos[i]]; 
                fp(n, index[i]+1) = SBE(ptr[0]);
               }
               break;
           }
       }    
    }
    delete [] buf;
  }

  else if (sizer == 4 && dataPacking == 1) // old point code
  {
    vector<float> buf(nfield);
    FldArrayF& fp = *f;
    for (n = 0; n < npts; n++)
    {
        fread(&buf[0], sizeof(float), nfield, ptrFile);
        for (i = 0; i < nfield; i++) fp(n, i+1) = FBE(buf[i]);
    }
  }
  else if (sizer == 8 && dataPacking == 1) // point
  {
    vector<double> buf(nfield);
    FldArrayF& fp = *f;
    for (n = 0; n < npts; n++)
    {
        fread(&buf[0], sizeof(E_Float), nfield, ptrFile);
        for (i = 0; i < nfield; i++) fp(n, i+1) = DBE(buf[i]);
    }
  }

  delete [] passive;
  
  // Connectivity
  if (et != 8) // elements basiques
  {
    vector<int> buf2(eltType);
    for (n = 0; n < nelts; n++)
    {
      fread(&buf2[0], si, eltType, ptrFile);
      if (version <= 101)
      { for (i = 0; i < eltType; i++) c(n, i+1) = IBE(buf2[i]); }
      else  // la numerotation a change en version 108!
      { for (i = 0; i < eltType; i++) c(n, i+1) = IBE(buf2[i])+1; }
    }
    if (rawlocal == 1)
    {
      int b;
      for (n = 0; n < nelts; n++)
      {
        fread(&b, si, 1, ptrFile);
      }
    }
  }
  else // NGON
  {
    // read face offset
    vector<int> offset(numFaces+1);
    if (numFaceNodes != 2*numFaces) // volumic
    {
      fread(&offset[0], si, numFaces+1, ptrFile);
    }
    else // surfacic
    {
      for (E_Int i = 0; i < numFaces+1; i++) offset[i] = 2*i;
    }
    // faceNodes
    vector<int> faceNodes(numFaceNodes);
    fread(&faceNodes[0], si, numFaceNodes, ptrFile);

    // left elements
    vector<int> left(numFaces);
    //int* left = new int[numFaces];
    fread(&left[0], si, numFaces, ptrFile);

    // right elements
    vector<int> right(numFaces);
    fread(&right[0], si, numFaces, ptrFile);

    // skip boundary data
    if (numBoundaryFaces != 0)
    {
      for (E_Int i = 0; i < numBoundaryFaces; i++)
        fread(&ib, si, 1, ptrFile);
      for (E_Int i = 0; i < numBoundaryConnections; i++)
        fread(&ib, si, 1, ptrFile);
      for (E_Int i = 0; i < numBoundaryConnections; i++)
        fread(&ib, si, 1, ptrFile);
    }

    // Construit la connectivite array
    vector<int> count(numFaces);
    //int* count = new int[numFaces];
    E_Int size = 0;
    for (E_Int i = 0; i < numFaces; i++)
    {
      count[i] = IBE(offset[i+1])-IBE(offset[i]);
      size += count[i]+1;
    }
    // build FN
    FldArrayI faces(size);
    E_Int* ptr = faces.begin();
    E_Int n;
    E_Int p = 0;
    for (E_Int i = 0; i < numFaces; i++)
    {
      n = count[i];
      ptr[0] = n;
      for (E_Int j = 0; j < n; j++)
      {
        ptr[j+1] = IBE(faceNodes[p])+1; p++;
      }
      ptr += n+1;
    }
    vector<int>().swap(count); vector<int>().swap(offset);
    //delete [] count; delete [] offset;

    // cFE
    FldArrayI cFE(numFaces, 2);
    E_Int* cFE1 = cFE.begin(1); E_Int* cFE2 = cFE.begin(2);
    for (E_Int i = 0; i < numFaces; i++)
    {
      cFE1[i] = IBE(left[i])+1; cFE2[i] = IBE(right[i])+1;
    }
    vector<int>().swap(left); vector<int>().swap(right);
    //delete [] left; delete [] right;

    // Construit la connectivite elts->faces
    FldArrayI cEF;
    K_CONNECT::connectFE2EF(cFE, ne, cEF);

    // Fusionne les connectivites
    c.malloc(faces.getSize()+cEF.getSize()+4);
    E_Int* cp = c.begin();
    cp[0] = numFaces;
    cp[1] = faces.getSize();
    ptr = faces.begin();
    for (E_Int i = 0; i < faces.getSize(); i++)
    {
      cp[i+2] = ptr[i];
    }
    int pt = 2+faces.getSize();
    cp[pt] = ne;
    cp[pt+1] = cEF.getSize();
    ptr = cEF.begin();
    for (E_Int i = 0; i < cEF.getSize(); i++)
    {
      cp[pt+2+i] = ptr[i];
    }
  }

  //if (sizer == 8 && dataPacking == 0) convertEndianField(*f);
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::readGeom108(FILE* ptrFile,
                               FldArrayF*& f)
{
  int ival;
  double dval;
  float fval;
  E_Int position;
  double coord0[3];
  E_Float x0, y0, z0;
  E_Int fieldType = 1;
  E_Int geomType = 0;
  fread(&ival, sizeof(int), 1, ptrFile); // position coordsys
  position = ival;
  fread(&ival, sizeof(int), 1, ptrFile); // Scope
  fread(&ival, sizeof(int), 1, ptrFile); // Draw order
  fread(coord0, sizeof(double), 3, ptrFile); // X0, Y0, Z0
  x0 = coord0[0]; y0 = coord0[1]; z0 = coord0[2];
  fread(&ival, sizeof(int), 1, ptrFile); // Zone
  fread(&ival, sizeof(int), 1, ptrFile); // Color
  fread(&ival, sizeof(int), 1, ptrFile); // Fill color
  fread(&ival, sizeof(int), 1, ptrFile); // is filled?
  fread(&ival, sizeof(int), 1, ptrFile); // GeomType
  geomType = ival;
  fread(&ival, sizeof(int), 1, ptrFile); // Line pattern
  fread(&dval, sizeof(double), 1, ptrFile); // Pattern length
  fread(&dval, sizeof(double), 1, ptrFile); // Line thickness
  fread(&ival, sizeof(int), 1, ptrFile); // Num ellipse pts
  fread(&ival, sizeof(int), 1, ptrFile); // Arrowhead style
  fread(&ival, sizeof(int), 1, ptrFile); // Arrowhead attachment
  fread(&dval, sizeof(double), 1, ptrFile); // Arrow head Size
  fread(&dval, sizeof(double), 1, ptrFile); // Arrow angle
  fread(&ival, sizeof(int), 1, ptrFile); // Macro
  fread(&ival, sizeof(int), 1, ptrFile); // Data type
  fieldType = ival;
  fread(&ival, sizeof(int), 1, ptrFile); // Clipping

  switch (geomType)
  {
    case 0: // polyline
    {
      fread(&ival, sizeof(int), 1, ptrFile); //npl = ival;
      E_Int npts = 0;
      fread(&ival, sizeof(int), 1, ptrFile); npts = ival;
      if (fieldType == 1)
      {
        f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&fval, sizeof(float), 1, ptrFile); (*f)(i,1) = x0 + fval;
        }
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&fval, sizeof(float), 1, ptrFile); (*f)(i,2) = y0 + fval;
        }
        if (position == 4)
        {
          for (E_Int i = 0; i < npts; i++)
          {
            fread(&fval, sizeof(float), 1, ptrFile); (*f)(i,3) = z0 + fval;
          }
        }
      }
      else
      {
        f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&dval, sizeof(double), 1, ptrFile); (*f)(i,1) = x0 + dval;
        }
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&dval, sizeof(double), 1, ptrFile); (*f)(i,2) = y0 + dval;
        }
        if (position == 4)
        {
          for (E_Int i = 0; i < npts; i++)
          {
            fread(&dval, sizeof(double), 1, ptrFile); (*f)(i,3) = z0 + dval;
          }
        }
      }
    }
    return 1;

    case 1: // rectangle
    {
      f = new FldArrayF(5, 3); f->setAllValuesAtNull();
      E_Float Rx, Ry;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        Rx = fval;
        fread(&fval, sizeof(float), 1, ptrFile);
        Ry = fval;
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        Rx = dval;
        fread(&dval, sizeof(double), 1, ptrFile);
        Ry = dval;
      }
      (*f)(0,1) = x0; (*f)(0,2) = y0; (*f)(0,3) = z0;
      (*f)(1,1) = x0+Rx; (*f)(1,2) = y0; (*f)(1,3) = z0;
      (*f)(2,1) = x0+Rx; (*f)(2,2) = y0+Ry; (*f)(2,3) = z0;
      (*f)(3,1) = x0; (*f)(3,2) = y0+Ry; (*f)(3,3) = z0;
      (*f)(4,1) = x0; (*f)(4,2) = y0; (*f)(4,3) = z0;
    }
    return 1;

    case 2: // square
    {
      f = new FldArrayF(5, 3); f->setAllValuesAtNull();
      E_Float R;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        R = fval;
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        R = dval;
      }
      (*f)(0,1) = x0; (*f)(0,2) = y0; (*f)(0,3) = z0;
      (*f)(1,1) = x0+R; (*f)(1,2) = y0; (*f)(1,3) = z0;
      (*f)(2,1) = x0+R; (*f)(2,2) = y0+R; (*f)(2,3) = z0;
      (*f)(3,1) = x0; (*f)(3,2) = y0+R; (*f)(3,3) = z0;
      (*f)(4,1) = x0; (*f)(4,2) = y0; (*f)(4,3) = z0;
    }
    return 1;

    case 3: // circle
    {
      E_Float teta;
      E_Float pi = 4*atan(1.);
      E_Int npts = 100;
      f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
      E_Float R;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        R = fval;
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        R = dval;
      }

      for (E_Int i = 0; i < npts; i++)
      {
        teta = i*2*pi/(npts - 1.);
        (*f)(i,1) = x0 + R*cos(teta);
        (*f)(i,2) = y0 + R*sin(teta);
      }
    }
    return 1;

    case 4: // ellipse
    {
      E_Float teta;
      E_Float pi = 4*atan(1.);
      E_Int npts = 100;
      f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
      E_Float Rx, Ry;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        Rx = fval;
        fread(&fval, sizeof(float), 1, ptrFile);
        Ry = fval;
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        Rx = dval;
        fread(&dval, sizeof(double), 1, ptrFile);
        Ry = dval;
      }

      for (E_Int i = 0; i < npts; i++)
      {
        teta = i*2*pi/(npts - 1.);
        (*f)(i,1) = x0 + Rx*cos(teta);
        (*f)(i,2) = y0 + Ry*sin(teta);
      }
    }
    return 1;

    default:
      printf("Warning: readGeom: this kind of geometry is unknown: " SF_D_ ".\n",
             geomType);
  }
  return 0; // nothing created
}
//=============================================================================
E_Int K_IO::GenIO::readGeom108CE(FILE* ptrFile,
                                 FldArrayF*& f)
{
  int ival;
  double dval;
  float fval;
  E_Int position;
  double coord0[3];
  E_Float x0, y0, z0;
  E_Int fieldType = 1;
  E_Int geomType = 0;
  fread(&ival, sizeof(int), 1, ptrFile); // position coordsys
  position = IBE(ival);
  fread(&ival, sizeof(int), 1, ptrFile); // Scope
  fread(&ival, sizeof(int), 1, ptrFile); // Draw order
  fread(coord0, sizeof(double), 3, ptrFile); // X0, Y0, Z0
  x0 = DBE(coord0[0]); y0 = DBE(coord0[1]); z0 = DBE(coord0[2]);
  fread(&ival, sizeof(int), 1, ptrFile); // Zone
  fread(&ival, sizeof(int), 1, ptrFile); // Color
  fread(&ival, sizeof(int), 1, ptrFile); // Fill color
  fread(&ival, sizeof(int), 1, ptrFile); // is filled?
  fread(&ival, sizeof(int), 1, ptrFile); // GeomType
  geomType = IBE(ival);
  fread(&ival, sizeof(int), 1, ptrFile); // Line pattern
  fread(&dval, sizeof(double), 1, ptrFile); // Pattern length
  fread(&dval, sizeof(double), 1, ptrFile); // Line thickness
  fread(&ival, sizeof(int), 1, ptrFile); // Num ellipse pts
  fread(&ival, sizeof(int), 1, ptrFile); // Arrowhead style
  fread(&ival, sizeof(int), 1, ptrFile); // Arrowhead attachment
  fread(&dval, sizeof(double), 1, ptrFile); // Arrow head Size
  fread(&dval, sizeof(double), 1, ptrFile); // Arrow angle
  fread(&ival, sizeof(int), 1, ptrFile); // Macro
  fread(&ival, sizeof(int), 1, ptrFile); // Data type
  fieldType = IBE(ival);
  fread(&ival, sizeof(int), 1, ptrFile); // Clipping

  switch (geomType)
  {
    case 0: // polyline
    {
      fread(&ival, sizeof(int), 1, ptrFile); //npl = IBE(ival);
      E_Int npts = 0;
      fread(&ival, sizeof(int), 1, ptrFile); npts = IBE(ival);
      if (fieldType == 1)
      {
        f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&fval, sizeof(float), 1, ptrFile);
          (*f)(i,1) = x0 + FBE(fval);
        }
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&fval, sizeof(float), 1, ptrFile);
          (*f)(i,2) = y0 + FBE(fval);
        }
        if (position == 4)
        {
          for (E_Int i = 0; i < npts; i++)
          {
            fread(&fval, sizeof(float), 1, ptrFile);
            (*f)(i,3) = z0 + FBE(fval);
          }
        }
      }
      else
      {
        f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&dval, sizeof(double), 1, ptrFile);
          (*f)(i,1) = x0 + DBE(dval);
        }
        for (E_Int i = 0; i < npts; i++)
        {
          fread(&dval, sizeof(double), 1, ptrFile);
          (*f)(i,2) = y0 + DBE(dval);
        }
        if (position == 4)
        {
          for (E_Int i = 0; i < npts; i++)
          {
            fread(&dval, sizeof(double), 1, ptrFile);
            (*f)(i,3) = z0 + DBE(dval);
          }
        }
      }
    }
    return 1;

    case 1: // rectangle
    {
      f = new FldArrayF(5, 3); f->setAllValuesAtNull();
      E_Float Rx, Ry;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        Rx = FBE(fval);
        fread(&fval, sizeof(float), 1, ptrFile);
        Ry = FBE(fval);
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        Rx = DBE(dval);
        fread(&dval, sizeof(double), 1, ptrFile);
        Ry = DBE(dval);
      }
      (*f)(0,1) = x0; (*f)(0,2) = y0; (*f)(0,3) = z0;
      (*f)(1,1) = x0+Rx; (*f)(1,2) = y0; (*f)(1,3) = z0;
      (*f)(2,1) = x0+Rx; (*f)(2,2) = y0+Ry; (*f)(2,3) = z0;
      (*f)(3,1) = x0; (*f)(3,2) = y0+Ry; (*f)(3,3) = z0;
      (*f)(4,1) = x0; (*f)(4,2) = y0; (*f)(4,3) = z0;
    }
    return 1;

    case 2: // square
    {
      f = new FldArrayF(5, 3); f->setAllValuesAtNull();
      E_Float R;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        R = FBE(fval);
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        R = DBE(dval);
      }
      (*f)(0,1) = x0; (*f)(0,2) = y0; (*f)(0,3) = z0;
      (*f)(1,1) = x0+R; (*f)(1,2) = y0; (*f)(1,3) = z0;
      (*f)(2,1) = x0+R; (*f)(2,2) = y0+R; (*f)(2,3) = z0;
      (*f)(3,1) = x0; (*f)(3,2) = y0+R; (*f)(3,3) = z0;
      (*f)(4,1) = x0; (*f)(4,2) = y0; (*f)(4,3) = z0;
    }
    return 1;

    case 3: // circle
    {
      E_Float teta;
      E_Float pi = 4*atan(1.);
      E_Int npts = 100;
      f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
      E_Float R;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        R = FBE(fval);
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        R = DBE(dval);
      }

      for (E_Int i = 0; i < npts; i++)
      {
        teta = i*2*pi/(npts - 1.);
        (*f)(i,1) = x0 + R*cos(teta);
        (*f)(i,2) = y0 + R*sin(teta);
      }
    }
    return 1;

    case 4: // ellipse
    {
      E_Float teta;
      E_Float pi = 4*atan(1.);
      E_Int npts = 100;
      f = new FldArrayF(npts, 3); f->setAllValuesAtNull();
      E_Float Rx, Ry;
      if (fieldType == 1)
      {
        fread(&fval, sizeof(float), 1, ptrFile);
        Rx = FBE(fval);
        fread(&fval, sizeof(float), 1, ptrFile);
        Ry = FBE(fval);
      }
      else
      {
        fread(&dval, sizeof(double), 1, ptrFile);
        Rx = DBE(dval);
        fread(&dval, sizeof(double), 1, ptrFile);
        Ry = DBE(dval);
      }

      for (E_Int i = 0; i < npts; i++)
      {
        teta = i*2*pi/(npts - 1.);
        (*f)(i,1) = x0 + Rx*cos(teta);
        (*f)(i,2) = y0 + Ry*sin(teta);
      }
    }
    return 1;

    default:
      printf("Warning: readGeom: this kind of geometry is unknown: " SF_D_ ".\n",
             geomType);
  }
  return 0; // nothing created
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
  IN: eltType: element type:
  1 (BAR), 2 (TRI), 3 (QUAD), 4 (TETRA), 5 (PYRA), 6 (PENTA),
  7 (HEXA), 8 (NGON)
  This routine is v112 (sinon pas de NGON)
  return 1 if failed
*/
//=============================================================================
E_Int K_IO::GenIO::tecwrite108(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector <FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{ 
  FILE* ptrFile;
  float a;         // must be basic C types !!
  int   ib;
  char cop[BUFSIZE+1];
  char titre[] = "field";
  FldArrayF buffer;
  double t;
  E_Int no, nol, p, i, j, n, nt, nv, nf;

  // Constants
  E_Int si = sizeof(int);
  E_Int sf = sizeof(float);

  // Number of variables to write
  E_Int nfield = 0;
  if (structField.size() > 0) nfield = structField[0]->getNfld();
  else if (unstructField.size() > 0) nfield = unstructField[0]->getNfld();
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
    if (nvar != nfield)
    {
      printf("Warning: tecwrite: number of variables differs for structured field: " SF_D_ ".\n", n+1);
      return 1;
    }
  }

  E_Int unstructFieldSize = unstructField.size();
  for (E_Int n = 0; n < unstructFieldSize; n++)
  {
    E_Int nvar = unstructField[n]->getNfld();
    if (nvar != nfield)
    {
      printf("Warning: tecwrite: number of variables differs for unstructured field: " SF_D_ ".\n",n+1);
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
  strcpy(version,"#!TDV112");
  fwrite(version, sizeof(char), 8, ptrFile);

  // Magic number for endian identification
  ib = 1;
  fwrite(&ib, si, 1, ptrFile);

  // file type (ajoute en 111)
  ib = 0; // FULL
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
  ib = nfield;
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

  // To keep track of type change due to tecplot limitations
  E_Int changeME2HEXA = 0;
  E_Int changePENTA2HEXA = 0;

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

    // Parent zone
    ib = -1; /* no */
    fwrite(&ib, si, 1, ptrFile);

    // Strand id
    ib = -1; /* -2: assigned by tecplot, -1: static, other:valid */
    fwrite(&ib, si, 1, ptrFile);

    // time
    t = 0.;
    fwrite(&t, sizeof(double), 1, ptrFile);

    // unused
    ib = -1;
    fwrite(&ib, si, 1, ptrFile);

    // type de zones
    if (no < structFieldSize)
    {
      ib = 0; // structured
    }
    else
    {
      nol = no - structFieldSize;
      FldArrayI* c = connect[nol];
      E_Int nc = c->getNConnect();
      /* Type of elts */
      switch (eltTypes[nol][0])
      {
        case 1: // BAR
          ib = 1;
          break;
        case 2:  // TRI
          ib = 2;
          if (nc > 1) { ib = 3; } // FIX as QUAD for ME
          break;
        case 3: // QUAD
          ib = 3;
          break;
        case 4: // TETRA
          ib = 4;
          if (nc > 1) { ib = 5; } // FIX as HEXA for ME
          break;
        case 5: // PYRA - FIX as HEXA
          ib = 5; // Dans ce cas, on trace des hexa degeneres
          break;
        case 6: // PENTA - FIX as HEXA
          ib = 5; // Dans ce cas, on trace des hexa degeneres
          break;
        case 7: // HEXA
          ib = 5;
          break;
        case 8: // NGON
        {
          E_Int* ngon = c->getNGon();
          E_Int* indPG = c->getIndPG();
          E_Int nf; c->getFace(0, nf, ngon, indPG);
          if (nf > 2) ib = 7; // polyhedron
          else ib = 6; // polygon
        }
        break;
        default:
          printf("Warning: tecwrite: wrong element type.\n");
          return 1;
      }
    }
    fwrite(&ib, si, 1, ptrFile); // write type

    // data packing (supp in 112)
    //ib = 0; // block
    //fwrite(&ib, si, 1, ptrFile);

    // var location
    ib = 0; /* dont specify */
    fwrite(&ib, si, 1, ptrFile);

    // raw local
    ib = 0; /* dont specify */
    fwrite(&ib, si, 1, ptrFile);

    // user defined faces
    ib = 0; /* dont specify */
    fwrite(&ib, si, 1, ptrFile);

    if (no < structFieldSize)
    {
      /* ni, nj, nk */
      ib = ni[no]; fwrite(&ib, si, 1, ptrFile);
      ib = nj[no]; fwrite(&ib, si, 1, ptrFile);
      ib = nk[no]; fwrite(&ib, si, 1, ptrFile);
    }
    else
    {
      nol = no - structFieldSize;
      FldArrayI* c = connect[nol];
      E_Int nc = c->getNConnect();
      
      if (eltTypes[nol][0] != 8) // elements basiques
      {
        // numPts
        ib = unstructField[nol]->getSize();
        fwrite(&ib, si, 1, ptrFile);

        // num elts
        E_Int nelts = 0;
        for (E_Int n = 0; n < nc; n++) nelts += c->getConnect(n)->getSize(); 
        ib = nelts;
        //ib = c->getSize();
        fwrite(&ib, si, 1, ptrFile);

        // cellDim
        ib = 0; fwrite(&ib, si, 1, ptrFile);
        ib = 0; fwrite(&ib, si, 1, ptrFile);
        ib = 0; fwrite(&ib, si, 1, ptrFile);
      }
      else // NGON
      {
        // numPts
        ib = unstructField[nol]->getSize();
        fwrite(&ib, si, 1, ptrFile);
        // num faces
        ib = c->getNFaces(); fwrite(&ib, si, 1, ptrFile);
        // numFacesNodes
        E_Int isNGon = c->getNGonType();
        E_Int size = c->getSizeNGon();
        if (isNGon != 3) size -= c->getNFaces();
        ib = size; fwrite(&ib, si, 1, ptrFile);
        // Boundary faces
        ib = 0; fwrite(&ib, si, 1, ptrFile);
        // Boundary connections
        ib = 0; fwrite(&ib, si, 1, ptrFile);
        // num elts
        ib = c->getNElts(); fwrite(&ib, si, 1, ptrFile);
        // cellDim
        ib = 0; fwrite(&ib, si, 1, ptrFile);
        ib = 0; fwrite(&ib, si, 1, ptrFile);
        ib = 0; fwrite(&ib, si, 1, ptrFile);
      }
    }

    // auxiliary data
    ib = 0;
    fwrite(&ib, si, 1, ptrFile);

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

    // format of data
    for (i = 0; i < nfield; i++)
    {
#ifdef E_DOUBLEREAL
      ib = 2;
#else
      ib = 1;
#endif
      fwrite(&ib, si, 1, ptrFile);
    }

    // passive
    ib = 0;
    fwrite(&ib, si, 1, ptrFile);

    // sharing
    ib = 0;
    fwrite(&ib, si, 1, ptrFile);

    // share connect
    ib = -1;
    fwrite(&ib, si, 1, ptrFile);

    /* Min-Max since no sharing and no passive. */
    for (n = 0; n < nfield; n++)
    {
      E_Float* fp = f.begin(n+1);
      E_Int nt = f.getSize();
      E_Float fmin = K_CONST::E_MAX_FLOAT;
      E_Float fmax = -K_CONST::E_MAX_FLOAT;
      for (E_Int i = 0; i < nt; i++)
      { fmin = K_FUNC::E_min(fmin, fp[i]);
        fmax = K_FUNC::E_max(fmax, fp[i]); }
      t = fmin;
      fwrite(&t, sizeof(double), 1, ptrFile);
      t = fmax;
      fwrite(&t, sizeof(double), 1, ptrFile);
    }

    // field
    for (n = 1; n <= nfield; n++)
      fwrite(f.begin(n), sizeof(E_Float), f.getSize(), ptrFile);
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

    // format of data
    for (i = 0; i < nfield; i++)
    {
#ifdef E_DOUBLEREAL
      ib = 2;
#else
      ib = 1;
#endif
      fwrite(&ib, si, 1, ptrFile);
    }

    // passive
    ib = 0;
    fwrite(&ib, si, 1, ptrFile);

    // sharing
    ib = 0;
    fwrite(&ib, si, 1, ptrFile);

    // share connect
    ib = -1;
    fwrite(&ib, si, 1, ptrFile);

    // min max
    for (n = 0; n < nfield; n++)
    {
      E_Float* fp = f.begin(n+1);
      E_Int nt = f.getSize();
      E_Float fmin = K_CONST::E_MAX_FLOAT;
      E_Float fmax = -K_CONST::E_MAX_FLOAT;
      for (E_Int i = 0; i < nt; i++)
      { fmin = K_FUNC::E_min(fmin, fp[i]);
        fmax = K_FUNC::E_max(fmax, fp[i]); }
      t = fmin;
      fwrite(&t, sizeof(double), 1, ptrFile);
      t = fmax;
      fwrite(&t, sizeof(double), 1, ptrFile);
    }

    // field
    for (n = 1; n <= nfield; n++)
      fwrite(f.begin(n), sizeof(E_Float), f.getSize(), ptrFile);

    // Connectivity
    E_Int nc = c.getNConnect();
    int* bufferi;
    E_Int sizet = 0;
    if (eltTypes[no][0] == 8) // NGON
    {
      E_Int* ngon = c.getNGon();
      E_Int* indPG = c.getIndPG();
      E_Int numFaces = c.getNFaces();
      E_Int isNGon = c.getNGonType();
      E_Int size = c.getSizeNGon();
      if (isNGon == 3) size += numFaces;
      E_Int nf; c.getFace(0, nf, ngon, indPG);
      if (nf > 2) sizet = numFaces+1 + (size-numFaces) + 2*numFaces;
      else sizet = (size-numFaces) + 2*numFaces;
    }
    else if (nc > 1) // ME
    {
      for (E_Int n = 0; n < nc; n++)
      {
        FldArrayI* cl = c.getConnect(n);
        switch (eltTypes[no][n])
        {
          case 2: sizet += 4*cl->getSize(); break; // ALL FIX
          case 3: sizet += 4*cl->getSize(); break;
          case 4: sizet += 8*cl->getSize(); break;
          case 5: sizet += 8*cl->getSize(); break;
          case 6: sizet += 8*cl->getSize(); break;
          case 7: sizet += 8*cl->getSize(); break;
          default: break;
        }
      }
    }
    else // BE
    {
      switch (eltTypes[no][0])
      {
        case 1: sizet += 2*c.getSize(); break; // BAR
        case 2: sizet += 3*c.getSize(); break; // TRI
        case 3: sizet += 4*c.getSize(); break; // QUAD
        case 4: sizet += 4*c.getSize(); break; // TETRA
        case 5: sizet += 8*c.getSize(); break; // PYRA->HEXA
        case 6: sizet += 8*c.getSize(); break; // PENTA->HEXA
        case 7: sizet += 8*c.getSize(); break; // HEXA
        default: break;
      }
    }
    
    bufferi = new int [sizet];

    if (eltTypes[no][0] == 8) // NGON
    {
      E_Int* ngon = c.getNGon();
      E_Int* indPG = c.getIndPG();
      E_Int numFaces = c.getNFaces();
      E_Int isNGon = c.getNGonType();
      E_Int size = c.getSizeNGon();
      if (isNGon == 3) size += numFaces;
      E_Int nf; c.getFace(0, nf, ngon, indPG);
      if (nf > 2) nt = numFaces+1 + (size-numFaces) + 2*numFaces;
      else nt = (size-numFaces) + 2*numFaces;
      
      // face offset
      int* ptri = bufferi;
      if (nf > 2) // only for volumic
      {
        ptri[0] = 0;
        for (E_Int i = 0; i < numFaces; i++)
        {
          c.getFace(i, n, ngon, indPG);
          ptri[i+1] = ptri[i] + n;
        }
        ptri += numFaces+1;
      }

      // face nodes
      for (E_Int i = 0; i < numFaces; i++)
      {
        E_Int* face = c.getFace(i, n, ngon, indPG);
        for (E_Int j = 0; j < n; j++) ptri[j] = face[j]-1;
        ptri += n;
      }

      FldArrayI cFE;
      K_CONNECT::connectNG2FE(c, cFE);
      E_Int* cFE1 = cFE.begin(1);
      E_Int* cFE2 = cFE.begin(2);

      // left
      for (E_Int i = 0; i < numFaces; i++) ptri[i] = cFE1[i]-1;
      ptri += numFaces;
      // right
      for (E_Int i = 0; i < numFaces; i++) ptri[i] = cFE2[i]-1;
    }
    else // ME and BE
    {
      int* ptri = bufferi;
      for (E_Int n = 0; n < nc; n++)
      {
        FldArrayI& cl = *(c.getConnect(n));
        nt = cl.getSize(); nv = cl.getNfld();
        if (eltTypes[no][n] == 5) // FIX pour PYRA as HEXA
        {
          changePENTA2HEXA += 1;
          for (i = 0; i < nt; i++)
          {
            p = i * 8;
            ptri[p  ] = cl(i, 1)-1;
            ptri[p+1] = cl(i, 2)-1;
            ptri[p+2] = cl(i, 3)-1;
            ptri[p+3] = cl(i, 4)-1;
            ptri[p+4] = cl(i, 5)-1;
            ptri[p+5] = cl(i, 5)-1;
            ptri[p+6] = cl(i, 5)-1;
            ptri[p+7] = cl(i, 5)-1;
          }
          ptri += nt*8;
        }
        else if (eltTypes[no][n] == 6) // FIX pour PENTA as HEXA
        {
          changePENTA2HEXA += 1;
          for (i = 0; i < nt; i++)
          {
            p = i * 8;
            ptri[p  ] = cl(i, 1)-1;
            ptri[p+1] = cl(i, 2)-1;
            ptri[p+2] = cl(i, 2)-1;
            ptri[p+3] = cl(i, 3)-1;
            ptri[p+4] = cl(i, 4)-1;
            ptri[p+5] = cl(i, 5)-1;
            ptri[p+6] = cl(i, 5)-1;
            ptri[p+7] = cl(i, 6)-1;
          }
          ptri += nt*8;
        }
        else if (eltTypes[no][n] == 2 && nc > 1) // FIX pour TRI as QUAD for ME
        {
          changeME2HEXA += 1;
          for (i = 0; i < nt; i++)
          {
            p = i * 4;
            ptri[p  ] = cl(i, 1)-1;
            ptri[p+1] = cl(i, 2)-1;
            ptri[p+2] = cl(i, 3)-1;
            ptri[p+3] = cl(i, 3)-1;
          }
          ptri += nt*4;
        }
        else if (eltTypes[no][n] == 4 && nc > 1) // FIX pour TETRA as HEXA for ME
        {
          changeME2HEXA += 1;
          for (i = 0; i < nt; i++)
          {
            p = i * 8;
            ptri[p  ] = cl(i, 1)-1;
            ptri[p+1] = cl(i, 2)-1;
            ptri[p+2] = cl(i, 3)-1;
            ptri[p+3] = cl(i, 3)-1;
            ptri[p+4] = cl(i, 4)-1;
            ptri[p+5] = cl(i, 4)-1;
            ptri[p+6] = cl(i, 4)-1;
            ptri[p+7] = cl(i, 4)-1;
          }
          ptri += nt*8;
        }
        else // CAS standard
        {
          for (i = 0; i < nt; i++)
          {
            p = i * nv;
            for (nf = 1; nf <= nv; nf++) ptri[p+nf-1] = cl(i, nf)-1;
          }
          ptri += nt*nv;
        }
      }
    }

    fwrite(bufferi, si, sizet, ptrFile);
    delete [] bufferi;

    no++;
  }

  if (changeME2HEXA > 0) printf("Warning: tecwrite: I changed some multi-elements to HEXA or QUAD.\n");
  if (changePENTA2HEXA > 0) printf("Warning: tecwrite: I changed PENTA or PYRA to HEXA.\n");
  
  fclose(ptrFile);
  return 0;
}
