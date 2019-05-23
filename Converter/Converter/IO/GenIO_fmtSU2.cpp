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

// SU2 (Stanford) file support

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
/* su2read
   Lit une zone comme plusieurs zones en elements basiques.
*/
//=============================================================================
E_Int K_IO::GenIO::su2read(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames,
  vector<FldArrayI*>& BCFaces, vector<char*>& BCNames)
{
  E_Int res; E_Int ti;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: su2read: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de la dimension
  res = readGivenKeyword(ptrFile, "NDIME=");
  res = readInt(ptrFile, ti, -1);
  E_Int dim = ti;
  if (dim != 2 && dim != 1 && dim != 3)
    printf("Warning: su2read: dimension is strange (%d).\n", ti);

  // Lecture du nombre d'elements
  res = readGivenKeyword(ptrFile, "NELEM=");
  res = readInt(ptrFile, ti, -1);
  E_Int ne = ti;
  //printf("ne=%d, res=%d\n", ne, res);

  //===========================================================================
  // Lecture polyedrique
  //============================================================================
  // On compte les faces et les elements
  /*
  E_Int nfaces = 0;
  E_Int sizeFN = 0;
  E_Int nelts = ne;
  E_Int sizeEF = 0;
  
  E_LONG pos = KFTELL(ptrFile);

  E_Int c = 0;
  while (c < ne)
  { 
    res = readInt(ptrFile, ti, -1); // type d'element
    switch (ti)
    { 
      case 3: //BAR
        nfaces += 2;
        sizeFN += 4;
        sizeEF += 3;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 5: // TRI
        nfaces += 3;
        sizeFN += 9;
        sizeEF += 4;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 9: // QUAD
        nfaces += 4;
        sizeFN += 12;
        sizeEF += 5;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 10: // TETRA
        nfaces += 4;
        sizeFN += 16;
        sizeEF += 5;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 12: // HEXA
        nfaces += 6;
        sizeFN += 30;
        sizeEF += 7;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 13: // PENTA
        nfaces += 5;
        sizeFN += 20;
        sizeEF += 6;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 14: // PYRA
        nfaces += 5;
        sizeFN += 17;
        sizeEF += 6;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;
    }
    res = readInt(ptrFile, ti, -1); // index
    c++;
  }

  co = new FldArrayI(2+sizeFN+2+sizeEF);
  cf = co->begin();
  cf[0] = nfaces;
  cf[1] = sizeFN;
  ce = co+2+sizeFN;
  ce[0] = nelts;
  ce[1] = sizeEF;

  c = 0; pf = 0; pe = 0; face = 1;
  while (c < ne)
  { 
    res = readInt(ptrFile, ti, -1); // type d'element
    switch (ti)
    { 
      case 3: //BAR
        cf[pf] = 1;
        res = readInt(ptrFile, ti, -1); cf[pf+1] = ti+1; pf += 2;
        cf[pf] = 1;
        res = readInt(ptrFile, ti, -1); cf[pf+1] = ti+1; pf += 2;
        ce[pe] = 2;
        ce[pe+1] = face+1;
        ce[pe+2] = face+2; face += 2;
        break;

      case 5: // TRI
        res = readInt(ptrFile, ti1, -1);
        res = readInt(ptrFile, ti2, -1);
        res = readInt(ptrFile, ti3, -1);
        cf[pf] = 2;
        cf[pf+1] = ti1+1; cf[pf+2] = ti2+1; pf += 3;
        cf[pf] = 2;
        cf[pf+1] = ti2+1; cf[pf+2] = ti3+1; pf += 3;
        cf[pf] = 2;
        cf[pf+1] = ti3+1; cf[pf+2] = ti1+1; pf += 3;
        ce[pe] = 3;
        ce[pe+1] = face+1;
        ce[pe+2] = face+2; 
        ce[pe+3] = face+3; face += 3;
        break;

      case 9: // QUAD
        res = readInt(ptrFile, ti1, -1);
        res = readInt(ptrFile, ti2, -1);
        res = readInt(ptrFile, ti3, -1);
        res = readInt(ptrFile, ti4, -1);
        cf[pf] = 2;
        cf[pf+1] = ti1+1; cf[pf+2] = ti2+1; pf += 3;
        cf[pf] = 2;
        cf[pf+1] = ti2+1; cf[pf+2] = ti3+1; pf += 3;
        cf[pf] = 2;
        cf[pf+1] = ti3+1; cf[pf+2] = ti4+1; pf += 3;
        cf[pf] = 2;
        cf[pf+1] = ti4+1; cf[pf+2] = ti1+1; pf += 3;
        ce[pe] = 4;
        ce[pe+1] = face+1;
        ce[pe+2] = face+2; 
        ce[pe+3] = face+3; 
        ce[pe+4] = face+4; face += 4;
        break;

      case 10: // TETRA
        res = readInt(ptrFile, ti1, -1);
        res = readInt(ptrFile, ti2, -1);
        res = readInt(ptrFile, ti3, -1);
        res = readInt(ptrFile, ti4, -1);
        cf[pf] = 3;
        cf[pf+1] = ti1+1; cf[pf+2] = ti2+1; cf[pf+3] = ti3+1; pf += 4;
        cf[pf] = 3;
        cf[pf+1] = ti1+1; cf[pf+2] = ti2+1; cf[pf+3] = ti4+1; pf += 4;
        cf[pf] = 3;
        cf[pf+1] = ti2+1; cf[pf+2] = ti3+1; cf[pf+3] = ti4+1; pf += 4;
        cf[pf] = 3;
        cf[pf+1] = ti3+1; cf[pf+2] = ti1+1; cf[pf+3] = ti4+1; pf += 4;
        ce[pe] = 4;
        ce[pe+1] = face+1;
        ce[pe+2] = face+2; 
        ce[pe+3] = face+3; 
        ce[pe+4] = face+4; face += 4;
        break;

      case 12: // HEXA
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 13: // PENTA
        nfaces += 5;
        sizeFN += 20;
        sizeEF += 6;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 14: // PYRA
        nfaces += 5;
        sizeFN += 17;
        sizeEF += 6;
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;
    }
    res = readInt(ptrFile, ti, -1); // index
    c++;
  }
  */

  //============================================================================
  // Lecture par elements basiques
  //============================================================================
  E_Int BARS = 0;
  E_Int TRIS = 0;
  E_Int QUADS = 0;
  E_Int TETRAS = 0;
  E_Int HEXAS = 0;
  E_Int PENTAS = 0;
  E_Int PYRAS = 0;
  
  //if (res == 1) skipLine(ptrFile);
  E_LONG pos = KFTELL(ptrFile);

  E_Int c = 0;
  while (c < ne)
  { 
    res = readInt(ptrFile, ti, -1); // type d'element
    //printf("%d\n", ti);
    switch (ti)
    { 
      case 3: //BAR
        BARS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 5: // TRI
        TRIS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 9: // QUAD
        QUADS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 10: // TETRA
        TETRAS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 12: // HEXA
        HEXAS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 13: // PENTA
        PENTAS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;

      case 14: // PYRA
        PYRAS++; 
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        res = readInt(ptrFile, ti, -1);
        break;
    }
    res = readInt(ptrFile, ti, -1); // index
    c++;
  }

  KFSEEK(ptrFile, pos, SEEK_SET);
  FldArrayI *cBAR, *cTRI, *cQUAD, *cTETRA, *cPENTA, *cPYRA, *cHEXA;
  E_Int *cBAR1, *cBAR2, *cTRI1, *cTRI2, *cTRI3, *cQUAD1, *cQUAD2;
  E_Int *cQUAD3, *cQUAD4, *cTETRA1, *cTETRA2, *cTETRA3, *cTETRA4;
  E_Int *cHEXA1, *cHEXA2, *cHEXA3, *cHEXA4, *cHEXA5, *cHEXA6, *cHEXA7, *cHEXA8;
  E_Int *cPENTA1, *cPENTA2, *cPENTA3, *cPENTA4, *cPENTA5, *cPENTA6;
  E_Int *cPYRA1, *cPYRA2, *cPYRA3, *cPYRA4, *cPYRA5;
  cBAR1 = NULL; cBAR2 = NULL; cTRI1 = NULL; cTRI2 = NULL; cTRI3 = NULL;
  cQUAD1 = NULL; cQUAD2 = NULL; cQUAD3 = NULL; cQUAD4 = NULL;
  cTETRA1 = NULL; cTETRA2 = NULL; cTETRA3 = NULL; cTETRA4 = NULL;
  cHEXA1 = NULL; cHEXA2 = NULL; cHEXA3 = NULL; cHEXA4 = NULL;
  cHEXA5 = NULL; cHEXA6 = NULL; cHEXA7 = NULL; cHEXA8 = NULL;
  cPENTA1 = NULL; cPENTA2 = NULL; cPENTA3 = NULL; cPENTA4 = NULL;
  cPENTA5 = NULL; cPENTA6 = NULL;
  cPYRA1 = NULL; cPYRA2 = NULL; cPYRA3 = NULL; cPYRA4 = NULL; cPYRA5 = NULL;

  if (BARS > 0) 
  {
    cBAR = new FldArrayI(BARS, 2);
    cBAR1 = cBAR->begin(1); cBAR2 = cBAR->begin(2);
  }
  if (TRIS > 0) 
  {
    cTRI = new FldArrayI(TRIS, 3);
    cTRI1 = cTRI->begin(1); cTRI2 = cTRI->begin(2); cTRI3 = cTRI->begin(3);
  }
  if (QUADS > 0) 
  {
    cQUAD = new FldArrayI(QUADS, 4);
    cQUAD1 = cQUAD->begin(1); cQUAD2 = cQUAD->begin(2); 
    cQUAD3 = cQUAD->begin(3); cQUAD4 = cQUAD->begin(4);
  }
  if (TETRAS > 0) 
  {
    cTETRA = new FldArrayI(TETRAS, 4);
    cTETRA1 = cTETRA->begin(1); cTETRA2 = cTETRA->begin(2); 
    cTETRA3 = cTETRA->begin(3); cTETRA4 = cTETRA->begin(4);
  }
  if (HEXAS > 0) 
  {
    cHEXA = new FldArrayI(HEXAS, 8);
    cHEXA1 = cHEXA->begin(1); cHEXA2 = cHEXA->begin(2);
    cHEXA3 = cHEXA->begin(3); cHEXA4 = cHEXA->begin(4);
    cHEXA5 = cHEXA->begin(5); cHEXA6 = cHEXA->begin(6);
    cHEXA7 = cHEXA->begin(7); cHEXA8 = cHEXA->begin(8);
  }
  if (PENTAS > 0) 
  {
    cPENTA = new FldArrayI(PENTAS, 6);
    cPENTA1 = cPENTA->begin(1); cPENTA2 = cPENTA->begin(2);
    cPENTA3 = cPENTA->begin(3); cPENTA4 = cPENTA->begin(4);
    cPENTA5 = cPENTA->begin(5); cPENTA6 = cPENTA->begin(6);
  }
  if (PYRAS > 0) 
  {
    cPYRA = new FldArrayI(PYRAS, 5);
    cPYRA1 = cPYRA->begin(1); cPYRA2 = cPYRA->begin(2);
    cPYRA3 = cPYRA->begin(3); cPYRA4 = cPYRA->begin(4);
    cPYRA5 = cPYRA->begin(5);
  }
  c = 0;
  BARS = 0; TRIS = 0; QUADS = 0; TETRAS = 0; HEXAS = 0; PENTAS = 0; PYRAS = 0;
  while (c < ne)
  { 
    res = readInt(ptrFile, ti, -1); // type d'element
    switch (ti)
    { 
      case 3: // BAR
        res = readInt(ptrFile, ti, -1); cBAR1[BARS] = ti+1;
        res = readInt(ptrFile, ti, -1); cBAR2[BARS] = ti+1;
        BARS++;
        break;

      case 5: // TRI
        res = readInt(ptrFile, ti, -1); cTRI1[TRIS] = ti+1;
        res = readInt(ptrFile, ti, -1); cTRI2[TRIS] = ti+1;
        res = readInt(ptrFile, ti, -1); cTRI3[TRIS] = ti+1;
        TRIS++;
        break;

      case 9: // QUAD
        res = readInt(ptrFile, ti, -1); cQUAD1[QUADS] = ti+1;
        res = readInt(ptrFile, ti, -1); cQUAD2[QUADS] = ti+1;
        res = readInt(ptrFile, ti, -1); cQUAD3[QUADS] = ti+1;
        res = readInt(ptrFile, ti, -1); cQUAD4[QUADS] = ti+1;
        QUADS++;
        break;

      case 10: // TETRA
        res = readInt(ptrFile, ti, -1); cTETRA1[TETRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cTETRA2[TETRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cTETRA3[TETRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cTETRA4[TETRAS] = ti+1;
        TETRAS++; 
        break;

      case 12: // HEXA
        res = readInt(ptrFile, ti, -1); cHEXA1[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA2[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA3[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA4[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA5[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA6[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA7[HEXAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cHEXA8[HEXAS] = ti+1;
        HEXAS++;
        break;

      case 13: // PENTA
        res = readInt(ptrFile, ti, -1); cPENTA1[PENTAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPENTA2[PENTAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPENTA3[PENTAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPENTA4[PENTAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPENTA5[PENTAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPENTA6[PENTAS] = ti+1;
        PENTAS++; 
        break;

      case 14: // PYRA
        res = readInt(ptrFile, ti, -1); cPYRA1[PYRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPYRA2[PYRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPYRA3[PYRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPYRA4[PYRAS] = ti+1;
        res = readInt(ptrFile, ti, -1); cPYRA5[PYRAS] = ti+1;
        PYRAS++;
        break;
    }
    res = readInt(ptrFile, ti, -1); // index of cell (unused here)
    c++;
  }
  // Lecture Vertices (Global)
  res = readGivenKeyword(ptrFile, "NPOIN=");
  res = readInt(ptrFile, ti, -1);
  E_Int np = ti;
  //printf("np %d - res=%d\n", np, res);
  // skip - parfois il semble y avoir un autre entier
  if (res == 1) skipLine(ptrFile);
  FldArrayF f(np, 3); E_Float fx, fy, fz;
  for (E_Int i = 0; i < np; i++)
  {
    fy = 0.; fz = 0.;
    if (dim == 1)
    {
      res = readDouble(ptrFile, fx, -1);
      res = readInt(ptrFile, ti, -1);
    }
    else if (dim == 2)
    { 
      res = readDouble(ptrFile, fx, -1);
      res = readDouble(ptrFile, fy, -1);
      res = readInt(ptrFile, ti, -1);
    }
    else 
    {
      res = readDouble(ptrFile, fx, -1);
      res = readDouble(ptrFile, fy, -1);
      res = readDouble(ptrFile, fz, -1);
      res = readInt(ptrFile, ti, -1);
    }
    //printf("%d\n", res); 
    if (res == 1) skipLine(ptrFile);
    f(ti,1)= fx; f(ti,2) = fy; f(ti,3) = fz;
    //printf("%f %f %f\n", f(ti,1), f(ti,2), f(ti,3));
  }

  // Formation des vertices de chacun
  if (BARS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cBAR);
    eltType.push_back(1);
  }
  if (TRIS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cTRI);
    eltType.push_back(2);
  }
  if (QUADS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cQUAD);
    eltType.push_back(3);
  }
  if (TETRAS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cTETRA);
    eltType.push_back(4);
  }
  if (HEXAS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cHEXA);
    eltType.push_back(7);
  }
  if (PENTAS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cPENTA);
    eltType.push_back(6);
  }
   if (PYRAS > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    connect.push_back(cPYRA);
    eltType.push_back(5);
  }

  // Cree le nom de zone
  for (unsigned int i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%d", i);
    zoneNames.push_back(zoneName);
  }
  //printf("%d %d\n", unstructField.size(), connect.size());

  varString = new char [8];
  strcpy(varString, "x,y,z");

  //====================================
  // Lecture des conditions aux limites
  //====================================
  E_Int nbnds, nfaces;
  res = readGivenKeyword(ptrFile, "NMARK=");
  if (res == 0) { fclose(ptrFile); return 0; }
  res = readInt(ptrFile, nbnds, -1);
  if (nbnds == 0) { fclose(ptrFile); return 0; }

  vector< vector<E_Int> > cVF(np);
  if (TRIS > 0) K_CONNECT::connectEV2VF(*cTRI, "TRI", cVF);
  if (QUADS > 0) K_CONNECT::connectEV2VF(*cQUAD, "QUAD", cVF);
  if (TETRAS > 0) K_CONNECT::connectEV2VF(*cTETRA, "TETRA", cVF);
  if (HEXAS > 0) K_CONNECT::connectEV2VF(*cHEXA, "HEXA", cVF);
  if (PENTAS > 0) K_CONNECT::connectEV2VF(*cPENTA, "PENTA", cVF);
  if (PYRAS > 0) K_CONNECT::connectEV2VF(*cPYRA, "PYRA", cVF);

  pos = KFTELL(ptrFile);
  res = readGivenKeyword(ptrFile, "MARKER_TAG=");
  E_Int inds[4]; E_Int indf;
  char buf[BUFSIZE];

  E_Int nfacetot = 0; res = 1;
  while (res == 1) // compte les faces
  {
    res = readWord(ptrFile, buf); //printf("%s\n", buf);
    res = readGivenKeyword(ptrFile, "MARKER_ELEMS=");
    res = readInt(ptrFile, nfaces, -1);
    nfacetot += nfaces;
    
    for (E_Int i = 0; i < nfaces; i++)
    { 
      res = readInt(ptrFile, ti, -1); // type d'element
      switch (ti)
      { 
        case 3: // BAR
          res = readInt(ptrFile, inds[0], -1);
          res = readInt(ptrFile, inds[1], -1);
          break;

        case 5: // TRI
          res = readInt(ptrFile, inds[0], -1);
          res = readInt(ptrFile, inds[1], -1);
          res = readInt(ptrFile, inds[2], -1);
          break;

        case 9: // QUAD
          res = readInt(ptrFile, inds[0], -1);
          res = readInt(ptrFile, inds[1], -1);
          res = readInt(ptrFile, inds[2], -1);
          res = readInt(ptrFile, inds[3], -1);
          break;
      }
      if (res == 1) skipLine(ptrFile);
    }
    res = readGivenKeyword(ptrFile, "MARKER_TAG=");
  }

  KFSEEK(ptrFile, pos, SEEK_SET);
  res = readGivenKeyword(ptrFile, "MARKER_TAG=");
  //printf("nafectot=%d\n", nfacetot);
  char* names = new char [BUFSIZE*nfacetot];
  BCNames.push_back(names);
  FldArrayI* faceList = new FldArrayI (nfacetot);
  BCFaces.push_back(faceList);
  E_Int* facep = faceList->begin();
  E_Int lenbuf;

  c = 0;
  E_Int p = 0; res = 1;
  while (res == 1) // trouve des BCS
  {
    res = readWord(ptrFile, buf);
    lenbuf = strlen(buf); //printf("%s %d\n", buf, lenbuf);

    res = readGivenKeyword(ptrFile, "MARKER_ELEMS=");
    res = readInt(ptrFile, nfaces, -1);
    //printf("%s %d\n", buf, nfaces);

    for (E_Int i = 0; i < nfaces; i++)
    { 
      res = readInt(ptrFile, ti, -1); // type d'element
      for (E_Int k = 0; k < lenbuf; k++) names[k+p] = buf[k];
      p += lenbuf; names[p] = '\0'; p++;

      switch (ti)
      { 
        case 3: // BAR
          res = readInt(ptrFile, inds[0], -1);
          res = readInt(ptrFile, inds[1], -1);
          indf = K_CONNECT::identifyFace(inds, 2, cVF);
          facep[c] = std::max(indf,1);
          break;

        case 5: // TRI
          res = readInt(ptrFile, inds[0], -1);
          res = readInt(ptrFile, inds[1], -1);
          res = readInt(ptrFile, inds[2], -1);
          indf = K_CONNECT::identifyFace(inds, 3, cVF);
          facep[c] = std::max(indf,1);
          break;

        case 9: // QUAD
          res = readInt(ptrFile, inds[0], -1);
          res = readInt(ptrFile, inds[1], -1);
          res = readInt(ptrFile, inds[2], -1);
          res = readInt(ptrFile, inds[3], -1);
          indf = K_CONNECT::identifyFace(inds, 4, cVF);
          facep[c] = std::max(indf,1);
          break;
      }
      if (res == 1) skipLine(ptrFile);
      c++;
    }
    res = readGivenKeyword(ptrFile, "MARKER_TAG=");
  }
  //printf("final: %d %d\n", c, nfacetot); 

  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Ecrit des zones si elles sont en elements basiques
//=============================================================================
E_Int K_IO::GenIO::su2write(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames, 
  PyObject* BCFaces)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // triangles, quads, tetra, hexa supported
    if (eltType[zone] == 1 || eltType[zone] == 2 || eltType[zone] == 3 ||
        eltType[zone] == 4 || eltType[zone] == 5 || eltType[zone] == 6 ||
        eltType[zone] == 7)
      nvalidZones++;
    else
      printf("Warning: su2write: zone %d not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: su2write: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  char format1[43]; char format2[85]; char format3[127]; 
  char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';
  
  // Build format for data
  sprintf(format1,"\t%s ", dataFmtl);
  sprintf(format2,"\t%s \t%s ", dataFmt, dataFmtl);
  sprintf(format3,"\t%s \t%s \t%s ", dataFmt, dataFmt, dataFmtl);
  strcat(format1,"\t%d \n");
  strcat(format2,"\t%d \n");
  strcat(format3,"\t%d \n");
  //printf("format=%s\n", format3);

  // Concatenate all vertices in one field
  E_Int size = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    size += unstructField[i]->getSize();
  }
  FldArrayF* vertices = new FldArrayF(size, 3);
  FldArrayF& v = *vertices;
  E_Float* v1 = v.begin(1);
  E_Float* v2 = v.begin(2);
  E_Float* v3 = v.begin(3);
  E_Int c = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayF& field = *unstructField[i];
    for (E_Int n = 0; n < field.getSize(); n++)
    {
      v1[c] = field(n, posx);
      v2[c] = field(n, posy);
      v3[c] = field(n, posz); c++;
    }
  }

  // Connectivite par elts
  E_Int shift = 0;
  vector<FldArrayI*> connectBar;
  vector<FldArrayI*> connectTri;
  vector<FldArrayI*> connectQuad;
  vector<FldArrayI*> connectTetra;
  vector<FldArrayI*> connectPyra;
  vector<FldArrayI*> connectPenta;
  vector<FldArrayI*> connectHexa;

  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayI& cn = *connect[i];
    E_Int elt = eltType[i];
    if (elt == 1) // Bar
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
      }
      connectBar.push_back(cpp);
    }
    else if (elt == 2) // Tri
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
        cp(n,3) = cn(n,3) + shift;
      }
      connectTri.push_back(cpp);
    }
    else if (elt == 3) // quads
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
        cp(n,3) = cn(n,3) + shift;
        cp(n,4) = cn(n,4) + shift;
      }
      connectQuad.push_back(cpp);
    }
    else if (elt == 4) // tetra
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
        cp(n,3) = cn(n,3) + shift;
        cp(n,4) = cn(n,4) + shift;
      }
      connectTetra.push_back(cpp);
    }
    else if (elt == 5) // pyra
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
        cp(n,3) = cn(n,3) + shift;
        cp(n,4) = cn(n,4) + shift;
        cp(n,5) = cn(n,5) + shift;
      }
      connectPyra.push_back(cpp);
    }
    else if (elt == 6) // penta
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
        cp(n,3) = cn(n,3) + shift;
        cp(n,4) = cn(n,4) + shift;
        cp(n,5) = cn(n,5) + shift;
        cp(n,6) = cn(n,6) + shift;
      }
      connectPenta.push_back(cpp);
    }
    else if (elt == 7) // hexa
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
        cp(n,3) = cn(n,3) + shift;
        cp(n,4) = cn(n,4) + shift;
        cp(n,5) = cn(n,1) + shift;
        cp(n,6) = cn(n,2) + shift;
        cp(n,7) = cn(n,3) + shift;
        cp(n,8) = cn(n,4) + shift;
      }
      connectHexa.push_back(cpp);
    }
    shift += unstructField[i]->getSize();
  }

  E_Int connectBarSize = connectBar.size();
  E_Int connectTriSize = connectTri.size();
  E_Int connectQuadSize = connectQuad.size();
  E_Int connectTetraSize = connectTetra.size();
  E_Int connectPyraSize = connectPyra.size();
  E_Int connectPentaSize = connectPenta.size();
  E_Int connectHexaSize = connectHexa.size();

  // Find dim
  E_Int dim = 1;
  if (connectTriSize > 0 || connectQuadSize > 0) dim = 2;
  if (connectTetraSize > 0 || connectPyraSize > 0 || 
      connectPentaSize > 0 || connectHexaSize > 0) dim = 3;

  // Find number of elements
  E_Int ne = 0;
  for (E_Int i = 0; i < connectBarSize; i++) ne += connectBar[i]->getSize();
  for (E_Int i = 0; i < connectTriSize; i++) ne += connectTri[i]->getSize();
  for (E_Int i = 0; i < connectQuadSize; i++) ne += connectQuad[i]->getSize();
  for (E_Int i = 0; i < connectTetraSize; i++) ne += connectTetra[i]->getSize();
  for (E_Int i = 0; i < connectPyraSize; i++) ne += connectPyra[i]->getSize();
  for (E_Int i = 0; i < connectPentaSize; i++) ne += connectPenta[i]->getSize();
  for (E_Int i = 0; i < connectHexaSize; i++) ne += connectHexa[i]->getSize();

  // Ecriture
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL)
  {
    printf("Warning: su2write: I can't open file %s.\n", file);
    return 1;  
  }

  fprintf(ptrFile, "NDIME= %d\n", dim);
  fprintf(ptrFile, "NELEM= %d\n", ne);

  c = 0;
  for (E_Int i = 0; i < connectBarSize; i++) 
  {
    FldArrayI& cp = *connectBar[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "3 \t%d \t%d \t%d \n", cp(i,1)-1, cp(i,2)-1, c); c++;
    }
  }

  for (E_Int i = 0; i < connectTriSize; i++) 
  {
    FldArrayI& cp = *connectTri[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "5 \t%d \t%d \t%d \t%d \n", 
              cp(i,1)-1, cp(i,2)-1, cp(i,3)-1, c); c++;
    }
  }
 
  for (E_Int i = 0; i < connectQuadSize; i++) 
  {
    FldArrayI& cp = *connectQuad[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "9 \t%d \t%d \t%d \t%d \t%d \n", 
              cp(i,1)-1, cp(i,2)-1, cp(i,3)-1, cp(i,4)-1, c); c++;
    }
  }
  
  for (E_Int i = 0; i < connectTetraSize; i++) 
  {
    FldArrayI& cp = *connectTetra[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "10 \t%d \t%d \t%d \t%d \t%d \n", 
              cp(i,1)-1, cp(i,2)-1, cp(i,3)-1, cp(i,4)-1, c); c++;
    }
  }

  for (E_Int i = 0; i < connectPyraSize; i++) 
  {
    FldArrayI& cp = *connectPyra[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "14 \t%d \t%d \t%d \t%d \t%d \t%d \n", 
              cp(i,1)-1, cp(i,2)-1, cp(i,3)-1, cp(i,4)-1, cp(i,5)-1, c); c++;
    }
  }
  
  for (E_Int i = 0; i < connectPentaSize; i++) 
  {
    FldArrayI& cp = *connectPenta[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "13 \t%d \t%d \t%d \t%d \t%d \t%d \t%d \n", 
              cp(i,1)-1, cp(i,2)-1, cp(i,3)-1, cp(i,4)-1, 
              cp(i,5)-1, cp(i,6)-1, c); c++;
    }
  }
  
  for (E_Int i = 0; i < connectHexaSize; i++) 
  {
    FldArrayI& cp = *connectHexa[i];
    for (E_Int i = 0; i < cp.getSize(); i++)
    {
      fprintf(ptrFile, "12 \t%d \t%d \t%d \t%d \t%d \t%d \t%d \n", 
              cp(i,1)-1, cp(i,2)-1, cp(i,3)-1, cp(i,4)-1, 
              cp(i,5)-1, cp(i,6)-1, c); c++;
    }
  }
  
  // Vertices
  c = 0;
  fprintf(ptrFile, "NPOIN= %d\n", v.getSize());
  if (dim == 1)
  {
    for (E_Int i = 0; i < v.getSize(); i++)
    { fprintf(ptrFile, format1, v(i,1), c); c++; }
  }
  else if (dim == 2)
  {
    for (E_Int i = 0; i < v.getSize(); i++)
    { fprintf(ptrFile, format2, v(i,1), v(i,2), c); c++; }
  }
  else if (dim == 3)
  {
    for (E_Int i = 0; i < v.getSize(); i++)
    { fprintf(ptrFile, format3, v(i,1), v(i,2), v(i,3), c); c++; }
  }

  delete vertices;
  
  for (E_Int i = 0; i < connectBarSize; i++) delete connectBar[i];
  for (E_Int i = 0; i < connectTriSize; i++) delete connectTri[i];
  for (E_Int i = 0; i < connectQuadSize; i++) delete connectQuad[i];
  for (E_Int i = 0; i < connectTetraSize; i++) delete connectTetra[i];
  for (E_Int i = 0; i < connectPyraSize; i++) delete connectPyra[i];
  for (E_Int i = 0; i < connectPentaSize; i++) delete connectPenta[i];
  for (E_Int i = 0; i < connectHexaSize; i++) delete connectHexa[i];

  // BC (if any)
  E_Int BCFacesSize = 0;
  if (PyList_Check(BCFaces) == true) BCFacesSize = PyList_Size(BCFaces);
  if (BCFacesSize > 0) // il y a des BCs
  {
    E_Int indFace, inde, nof;
    IMPORTNUMPY;
    shift = 0;
    E_Int face[6][4];
    for (E_Int i = 0; i < nzone; i++)
    {
      FldArrayI& cn = *connect[i];
      E_Int elt = eltType[i];
      //E_Int ne = cn.getSize();
      // nf: nbre de faces, nn : nbre de noeud par face
      E_Int eltBnd = 1; E_Int nf = 1; E_Int nn = 1;
      switch (elt)
      {
        case 2: // TRI
          eltBnd = 3; nf = 3; nn = 2;
          face[0][0] = 1; face[0][1] = 2;
          face[1][0] = 2; face[1][1] = 3;
          face[2][0] = 3; face[2][1] = 1;
          break;
        case 3: // QUAD
          eltBnd = 3; nf = 4; nn = 2;
          face[0][0] = 1; face[0][1] = 2;
          face[1][0] = 2; face[1][1] = 3;
          face[2][0] = 3; face[2][1] = 4;
          face[3][0] = 4; face[3][1] = 1;
          break;
        case 4: // TETRA
          eltBnd = 5; nf = 4; nn = 3;
          face[0][0] = 1; face[0][1] = 3; face[0][2] = 2;
          face[1][0] = 1; face[1][1] = 2; face[1][2] = 4;
          face[2][0] = 2; face[2][1] = 3; face[2][2] = 4;
          face[3][0] = 3; face[3][1] = 1; face[3][2] = 4;
          break;
        case 5: // PYRA
          eltBnd = 5; nf = 5; nn = 3; // vary
          face[0][0] = 1; face[0][1] = 4; face[0][2] = 3;
          face[1][0] = 3; face[1][1] = 2; face[1][2] = 1;
          face[2][0] = 1; face[2][1] = 2; face[2][2] = 5; 
          face[3][0] = 2; face[3][1] = 3; face[3][2] = 5;
          face[4][0] = 3; face[4][1] = 4; face[4][2] = 5; face[4][3] = 2;
          break;
        case 6: // PENTA
          eltBnd = 5; nf = 5; nn = 4; // vary
          face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4;
          face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5;
          face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6;
          face[3][0] = 1; face[3][1] = 3; face[3][2] = 2;
          face[4][0] = 4; face[4][1] = 5; face[4][2] = 6;
          break; // vary
          
        case 7: // HEXA
          eltBnd = 9; nf = 6; nn = 4;
          face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; face[0][3] = 2;
          face[1][0] = 1; face[1][1] = 2; face[1][2] = 6; face[1][3] = 5;
          face[2][0] = 2; face[2][1] = 3; face[2][2] = 7; face[2][3] = 6;
          face[3][0] = 3; face[3][1] = 4; face[3][2] = 8; face[3][3] = 7;
          face[4][0] = 1; face[4][1] = 5; face[4][2] = 8; face[4][3] = 4;
          face[5][0] = 5; face[5][1] = 6; face[5][2] = 7; face[5][3] = 8;
          break; 
      }
      PyObject* BCs = PyList_GetItem(BCFaces, i);
      E_Int size = PyList_Size(BCs);
      fprintf(ptrFile, "NMARK= %d\n", size/2);
      for (E_Int j = 0; j < size/2; j++) // marker differents
      {
        char* name = NULL; 
        PyObject* o = PyList_GetItem(BCs, 2*j);
        if (PyString_Check(o)) name = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(o)) name = PyBytes_AsString(PyUnicode_AsUTF8String(o));
#endif
        fprintf(ptrFile, "MARKER_TAG= %s\n", name);
        PyArrayObject* array = (PyArrayObject*)PyList_GetItem(BCs, 2*j+1);
        int* ptr = (int*)PyArray_DATA(array);
        E_Int np = PyArray_SIZE(array);
        fprintf(ptrFile, "MARKER_ELEMS= %d\n", np);
        for (E_Int j = 0; j < np; j++)
        {
          indFace = ptr[j];
          inde = (indFace-1)/nf;
          nof = (indFace-1)-inde*nf;
          fprintf(ptrFile, "%d ", eltBnd);
          if (elt == 5) // PYRA
          {
            if (nof == 4) // base
              for (E_Int i = 0; i < 4; i++) fprintf(ptrFile, "%d ", cn(inde,face[nof][i]));
            else
              for (E_Int i = 0; i < nn; i++) fprintf(ptrFile, "%d ", cn(inde,face[nof][i]));
          }
          else if (elt == 6) // PENTA
          {
            if (nof == 3 || nof == 4) 
              for (E_Int i = 0; i < 3; i++) fprintf(ptrFile, "%d ", cn(inde,face[nof][i]));
            else
              for (E_Int i = 0; i < nn; i++) fprintf(ptrFile, "%d ", cn(inde,face[nof][i]));
          }
          else
          {
            for (E_Int i = 0; i < nn; i++) fprintf(ptrFile, "%d ", cn(inde,face[nof][i]));
          }
          fprintf(ptrFile, "\n");
        }
      }
      shift += unstructField[i]->getSize();
    } // pour chaque zone
  }

  fclose(ptrFile);
  return 0;
}
