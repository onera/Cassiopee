/*    
    Copyright 2013-2024 Onera.

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

// Formated Mesh (INRIA) file support

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
/* meshread
   Manque: prise en compte des commentaires
   Nettoyage de f a la fin (vertex reellement utilises)
*/
//=============================================================================
E_Int K_IO::GenIO::meshread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t; E_Int ti;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: meshread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de l'entete
  char buf[BUFSIZE];
  res = readWord(ptrFile, buf);
  if (res == -1) { fclose(ptrFile); return 1; }
  if (strcmp(buf, "MeshVersionFormatted") != 0) { fclose(ptrFile); return 1; }

  res = readDouble(ptrFile, t, -1);
  if (t != 1)
    printf("Warning: meshread: version number is not really supported.\n");
  res = readGivenKeyword(ptrFile, "DIMENSION");
  res = readDouble(ptrFile, t, -1);
  if (t != 3)
  {
    printf("Warning: meshread: only 3D files are supported.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Lecture Vertices (Global)
  res = readGivenKeyword(ptrFile, "VERTICES");
  res = readDouble(ptrFile, t, -1);
  E_Int nb = E_Int(t);
  FldArrayF f(nb, 3);
  for (E_Int i = 0; i < nb; i++)
  {
    res = readDouble(ptrFile, t, -1); f(i,1) = t;
    res = readDouble(ptrFile, t, -1); f(i,2) = t;
    res = readDouble(ptrFile, t, -1); f(i,3) = t;
    res = readDouble(ptrFile, t, -1);  // ref discarded
    //printf("%f %f %f\n", f(i,1), f(i,2), f(i,3));
  }

  E_Int st;
  res = readWord(ptrFile, buf);

  E_Boolean foundEdge = false;
  E_Boolean foundTri = false;
  E_Boolean foundQuad = false;
  E_Boolean foundTetra = false;
  E_Boolean foundHexa = false;
  E_Boolean foundPenta = false;
  E_Boolean foundPyra = false;

  while (strcmp(buf, "End") != 0 && res >= 1)
  { 
    if (strcmp(buf, "Edges") == 0)
    {
      // Connectivity Edges (Global) -> 1 zone
      res = readInt(ptrFile, st, -1);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 2);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      connect.push_back(cn);
      eltType.push_back(1);
      foundEdge = true;
    }
    else if (strcmp(buf, "Triangles") == 0)
    {
      // Connectivity triangles (Global) -> 1 zone
      res = readInt(ptrFile, st, -1);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 3);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      connect.push_back(cn);
      eltType.push_back(2);
      foundTri = true;
    }
    else if (strcmp(buf, "Quadrilaterals") == 0)
    {
      // Connectivity Quadrangles (Global) -> 1 zone
      res = readInt(ptrFile, st, -1);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 4);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
        //printf("%d %d %d %d\n", c(i,1), c(i,2), c(i,3), c(i,4));
      }
      connect.push_back(cn);
      eltType.push_back(3);
      foundQuad = true;
    }
    else if (strcmp(buf, "Tetrahedra") == 0)
    {
      // Connectivity Tetrahedra (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 4);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      connect.push_back(cn);
      eltType.push_back(4);
      foundTetra = true;
    }
    else if (strcmp(buf, "Hexahedra") == 0)
    {
      // Connectivity Hexahedra (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 8);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); c(i,5) = ti;
        res = readInt(ptrFile, ti, -1); c(i,6) = ti;
        res = readInt(ptrFile, ti, -1); c(i,7) = ti;
        res = readInt(ptrFile, ti, -1); c(i,8) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      connect.push_back(cn);
      eltType.push_back(7);
      foundHexa = true;
    }
    else if (strcmp(buf, "Prisms") == 0)
    {
      // Connectivity Penta (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 6);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); c(i,5) = ti;
        res = readInt(ptrFile, ti, -1); c(i,6) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      connect.push_back(cn);
      eltType.push_back(6);
      foundPenta = true;
    }
    else if (strcmp(buf, "Pyramids") == 0)
    {
      // Connectivity Penta (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 5);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); c(i,5) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      connect.push_back(cn);
      eltType.push_back(5);
      foundPyra = true;
    }
    next: ;
    res = readWord(ptrFile, buf);
    //strcpy(buf, "End");
  }
  
  // Formation des vertices de chacun
  if (foundEdge)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (foundTri)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (foundQuad)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (foundTetra)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (foundHexa)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (foundPenta)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (foundPyra)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  // Cree le nom de zone
  for (size_t i = 0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%ld", i);
    zoneNames.push_back(zoneName);
  }
  //printf("%ld %ld\n", unstructField.size(), connect.size());

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);

  return 0;
}

E_Int K_IO::GenIO::meshread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<vector<E_Int> >& eltType, vector<char*>& zoneNames, E_Int api)
{
  E_Int res;
  E_Float t; E_Int ti;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: meshread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de l'entete
  char buf[BUFSIZE];
  res = readWord(ptrFile, buf);
  if (res == -1) { fclose(ptrFile); return 1; }
  if (strcmp(buf, "MeshVersionFormatted") != 0) { fclose(ptrFile); return 1; }

  res = readDouble(ptrFile, t, -1);
  if (t != 1)
    printf("Warning: meshread: version number is not really supported.\n");
  res = readGivenKeyword(ptrFile, "DIMENSION");
  res = readDouble(ptrFile, t, -1);
  if (t != 3)
  {
    printf("Warning: meshread: only 3D files are supported.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Lecture Vertices (Global)
  res = readGivenKeyword(ptrFile, "VERTICES");
  res = readDouble(ptrFile, t, -1);
  E_Int nb = E_Int(t);
  FldArrayF f(nb, 3);
  for (E_Int i = 0; i < nb; i++)
  {
    res = readDouble(ptrFile, t, -1); f(i,1) = t;
    res = readDouble(ptrFile, t, -1); f(i,2) = t;
    res = readDouble(ptrFile, t, -1); f(i,3) = t;
    res = readDouble(ptrFile, t, -1);  // ref discarded
  }

  E_Int st;
  res = readWord(ptrFile, buf);

  vector<FldArrayI*> tmpConnect;
  E_Boolean foundTopo[7] = {0};
  if (api == 3) eltType.resize(1);
  else eltType.clear();

  while (strcmp(buf, "End") != 0 && res >= 1)
  { 
    if (strcmp(buf, "Edges") == 0)
    {
      // Connectivity Edges (Global) -> 1 zone
      res = readInt(ptrFile, st, -1);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 2);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      
      if (api == 3)
      {
        eltType[0].push_back(1);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({1});
        connect.push_back(cn);
      }
      foundTopo[0] = true;
    }
    else if (strcmp(buf, "Triangles") == 0)
    {
      // Connectivity triangles (Global) -> 1 zone
      res = readInt(ptrFile, st, -1);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 3);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      if (api == 3)
      {
        eltType[0].push_back(2);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({2});
        connect.push_back(cn);
      }
      foundTopo[1] = true;
    }
    else if (strcmp(buf, "Quadrilaterals") == 0)
    {
      // Connectivity Quadrangles (Global) -> 1 zone
      res = readInt(ptrFile, st, -1);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 4);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      if (api == 3)
      {
        eltType[0].push_back(3);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({3});
        connect.push_back(cn);
      }
      foundTopo[2] = true;
    }
    else if (strcmp(buf, "Tetrahedra") == 0)
    {
      // Connectivity Tetrahedra (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 4);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      if (api == 3)
      {
        eltType[0].push_back(4);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({4});
        connect.push_back(cn);
      }
      foundTopo[3] = true;
    }
    else if (strcmp(buf, "Hexahedra") == 0)
    {
      // Connectivity Hexahedra (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 8);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); c(i,5) = ti;
        res = readInt(ptrFile, ti, -1); c(i,6) = ti;
        res = readInt(ptrFile, ti, -1); c(i,7) = ti;
        res = readInt(ptrFile, ti, -1); c(i,8) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      if (api == 3)
      {
        eltType[0].push_back(7);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({7});
        connect.push_back(cn);
      }
      foundTopo[4] = true;
    }
    else if (strcmp(buf, "Prisms") == 0)
    {
      // Connectivity Penta (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 6);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); c(i,5) = ti;
        res = readInt(ptrFile, ti, -1); c(i,6) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      if (api == 3)
      {
        eltType[0].push_back(6);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({6});
        connect.push_back(cn);
      }
      foundTopo[5] = true;
    }
    else if (strcmp(buf, "Pyramids") == 0)
    {
      // Connectivity Penta (Global) -> 1 zone
      res = readDouble(ptrFile, t, -1);
      st = E_Int(t);
      if (st == 0) goto next;
      FldArrayI* cn = new FldArrayI(st, 5);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        res = readInt(ptrFile, ti, -1); c(i,1) = ti;
        res = readInt(ptrFile, ti, -1); c(i,2) = ti;
        res = readInt(ptrFile, ti, -1); c(i,3) = ti;
        res = readInt(ptrFile, ti, -1); c(i,4) = ti;
        res = readInt(ptrFile, ti, -1); c(i,5) = ti;
        res = readInt(ptrFile, ti, -1); // discarded
      }
      if (api == 3)
      {
        eltType[0].push_back(5);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({5});
        connect.push_back(cn);
      }
      foundTopo[6] = true;
    }
    next: ;
    res = readWord(ptrFile, buf);
  }

  // Cree la varString
  varString = new char [8];
  strcpy(varString, "x,y,z");
  
  // Formation des vertices de chacun
  if (api == 3) // Create ME
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
    
    E_Int nfld = f.getNfld(); E_Int npts = f.getSize();
    E_Int nc = tmpConnect.size(); vector<E_Int> nepc(nc);
    for (E_Int ic = 0; ic < nc; ic++) nepc[ic] = tmpConnect[ic]->getSize();
    char eltString[256]; vector<E_Int> dummy(1);
    K_ARRAY::typeId2eltString(eltType[0], 0, eltString, dummy);

    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nepc,
                                         eltString, 0, api);
    FldArrayI* cn2; FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    #pragma omp parallel
    {
      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *tmpConnect[ic]->getConnect(0);
        FldArrayI& cm2 = *(cn2->getConnect(ic));
        #pragma omp for
        for (E_Int i = 0; i < cm2.getSize(); i++)
          for (E_Int j = 1; j <= cm2.getNfld(); j++)
            cm2(i,j) = cm(i,j);
      }
    }
    connect.push_back(cn2); 
    
    for (E_Int ic = 0; ic < nc; ic++) delete tmpConnect[ic];
    tmpConnect.clear();
    delete f2;
  }
  else // Create BEs
  {
    for (E_Int i = 0; i < 7; i++)
    { 
      if (foundTopo[i])
      {
        FldArrayF* an = new FldArrayF(f);
        unstructField.push_back(an);
      }
    }
  }
  
  // Cree le nom de zone
  for (size_t i = 0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%ld", i);
    zoneNames.push_back(zoneName);
  }

  fclose(ptrFile);

  return 0;
}

//=============================================================================
// Only write triangles, quads, tetra, hexa meshes. Others are discarded.
//=============================================================================
E_Int K_IO::GenIO::meshwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames,
  vector<vector<E_Int> >* colors)
{
  E_Int nzone = unstructField.size();

  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // triangles, quads, tetra, hexa, edges, supported
    if (eltType[zone] == 1 || eltType[zone] == 2 || eltType[zone] == 3 || 
        eltType[zone] == 4 || eltType[zone] == 5 || eltType[zone] == 6 || eltType[zone] == 7) 
      nvalidZones++;
    else
#ifdef E_DOUBLEINT
      printf("Warning: meshwrite: zone %ld not written (not a valid element type: %ld).", zone, eltType[zone]);
#else
      printf("Warning: meshwrite: zone %d not written (not a valid element type: %d).", zone, eltType[zone]);
#endif
  }

  if (nvalidZones == 0) return 1;

  // Check if this is the same zone with mutiple connectivity (based on coords)
  E_Int uniqueZone = 1;
  E_Int npts = -1;
  for (E_Int i = 0; i < nzone; i++)
  {
    E_Int n = unstructField[i]->getSize();
    if (npts == -1) npts = n;
    if (npts != n) { uniqueZone = 0; break; }
  }
  if (uniqueZone == 1) printf("[unique zone]...");

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: meshwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

 char format1[40]; char dataFmtl[40];
 strcpy(dataFmtl, dataFmt);
 int l = strlen(dataFmt); 
 if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format for data
  sprintf(format1,"%s%s%s", dataFmt, dataFmt, dataFmtl);
#ifdef E_DOUBLEINT
  strcat(format1," %ld\n");
#else
  strcat(format1," %d\n");
#endif

  // Concatenate all vertices in one field
  FldArrayF* vertices;
  if (uniqueZone == 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < nzone; i++)
    {
      size += unstructField[i]->getSize();
    }
    vertices = new FldArrayF(size, 3);
    FldArrayF& v = *vertices;
    E_Int c = 0;
    for (E_Int i = 0; i < nzone; i++)
    {
      FldArrayF& field = *unstructField[i];
      for (E_Int n = 0; n < field.getSize(); n++)
      {
        v(c, 1) = field(n, posx);
        v(c, 2) = field(n, posy);
        v(c, 3) = (posz>0) ? field(n, posz) : 0.; c++;
      }
    }
  }
  else
  {
    // shared coords
    E_Int size = unstructField[0]->getSize();
    vertices = new FldArrayF(size, 3);
    FldArrayF& v = *vertices;
    FldArrayF& field = *unstructField[0];
    for (E_Int n = 0; n < field.getSize(); n++)
    {
      v(n, 1) = field(n, posx);
      v(n, 2) = field(n, posy);
      v(n, 3) = (posz>0) ? field(n, posz) : 0.;
    }
  }
  FldArrayF& v = *vertices;
  
  // Connectivite par elts
  E_Int shift = 0;
  vector<FldArrayI*> connectEdge;
  vector<FldArrayI*> connectTri;
  vector<FldArrayI*> connectQuad;
  vector<FldArrayI*> connectTetra;
  vector<FldArrayI*> connectHexa;
  vector<FldArrayI*> connectPenta;
  vector<FldArrayI*> connectPyra;

  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayI& cn = *connect[i];
    E_Int elt = eltType[i];
    if (elt == 1) // Edge
    {
      FldArrayI* cpp = new FldArrayI(cn);
      FldArrayI& cp = *cpp;
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        cp(n,1) = cn(n,1) + shift;
        cp(n,2) = cn(n,2) + shift;
      }
      connectEdge.push_back(cpp);
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
        cp(n,5) = cn(n,5) + shift;
        cp(n,6) = cn(n,6) + shift;
        cp(n,7) = cn(n,7) + shift;
        cp(n,8) = cn(n,8) + shift;
      }
      connectHexa.push_back(cpp);
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
    if (uniqueZone == 0) shift += unstructField[i]->getSize();
  }

  // Ecriture
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: meshwrite: I can't open file %s.\n", file);
    return 1;  
  }

  fprintf(ptrFile, "MeshVersionFormatted\n1\n");
  fprintf(ptrFile, "Dimension\n3\n");
#ifdef E_DOUBLEINT
  fprintf(ptrFile, "Vertices\n%ld\n", v.getSize());
#else
  fprintf(ptrFile, "Vertices\n%d\n", v.getSize());
#endif
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, v(i,1), v(i,2), v(i,3), 0);
  
  E_Int connectEdgeSize = connectEdge.size();
  E_Int c;
  if (connectEdgeSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectEdgeSize; i++)
      size = size + connectEdge[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Edges\n%ld\n", size);
#else
    fprintf(ptrFile, "Edges\n%d\n", size);
#endif
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 1)
      {
        FldArrayI& cp = *connectEdge[c];
        if (!colors)
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld\n", cp(i,1), cp(i,2), c);
#else
            fprintf(ptrFile, "%d %d %d\n", cp(i,1), cp(i,2), c);
#endif
        else
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld\n", cp(i,1), cp(i,2), (*colors)[c][i]);
#else
            fprintf(ptrFile, "%d %d %d\n", cp(i,1), cp(i,2), (*colors)[c][i]);
#endif
        c++;
      }
    }
  }

  E_Int connectTriSize = connectTri.size();
  if (connectTriSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectTriSize; i++)
      size = size + connectTri[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Triangles\n%ld\n", size);
#else
    fprintf(ptrFile, "Triangles\n%d\n", size);
#endif
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 2)
      {
        FldArrayI& cp = *connectTri[c];
        if (!colors)
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld\n", cp(i,1), cp(i,2), cp(i,3), c);
#else
            fprintf(ptrFile, "%d %d %d %d\n", cp(i,1), cp(i,2), cp(i,3), c);
#endif
        else
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld\n", cp(i,1), cp(i,2), cp(i,3), (*colors)[c][i]);
#else
            fprintf(ptrFile, "%d %d %d %d\n", cp(i,1), cp(i,2), cp(i,3), (*colors)[c][i]);
#endif
        c++;
      }
    }
  }

  E_Int connectQuadSize = connectQuad.size();
  if (connectQuadSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectQuadSize; i++)
      size = size + connectQuad[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Quadrilaterals\n%ld\n", size);
#else
    fprintf(ptrFile, "Quadrilaterals\n%d\n", size);
#endif
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 3)
      {
        FldArrayI& cp = *connectQuad[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
          fprintf(ptrFile, "%ld %ld %ld %ld %ld\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#else
        fprintf(ptrFile, "%d %d %d %d %d\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#endif
        c++;
      }
    }
  }

  E_Int connectTetraSize = connectTetra.size();
  if (connectTetraSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectTetraSize; i++)
      size = size + connectTetra[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Tetrahedra\n%ld\n", size);
#else
    fprintf(ptrFile, "Tetrahedra\n%d\n", size);
#endif

    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 4)
      {
        FldArrayI& cp = *connectTetra[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
          fprintf(ptrFile, "%ld %ld %ld %ld %ld\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#else
          fprintf(ptrFile, "%d %d %d %d %d\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#endif
        c++;
      }
    }
  }

  E_Int connectHexaSize = connectHexa.size();
  if (connectHexaSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectHexaSize; i++)
      size = size + connectHexa[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Hexahedra\n%ld\n", size);
#else
    fprintf(ptrFile, "Hexahedra\n%d\n", size);
#endif
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 7)
      {
        FldArrayI& cp = *connectHexa[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
          fprintf(ptrFile, "%ld %ld %ld %ld %ld %ld %ld %ld %ld\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), cp(i,6), cp(i,7), cp(i,8), c);
#else
          fprintf(ptrFile, "%d %d %d %d %d %d %d %d %d\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), cp(i,6), cp(i,7), cp(i,8), c);
#endif
        c++;
      }
    }
  }
  E_Int connectPentaSize = connectPenta.size();
  if (connectPentaSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectPentaSize; i++)
      size = size + connectPenta[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Prisms\n%ld\n", size);
#else
    fprintf(ptrFile, "Prisms\n%d\n", size);
#endif
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 6)
      {
        FldArrayI& cp = *connectPenta[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
          fprintf(ptrFile, "%ld %ld %ld %ld %ld %ld %ld\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), cp(i,6), c);
#else
          fprintf(ptrFile, "%d %d %d %d %d %d %d\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), cp(i,6), c);
#endif
        c++;
      }
    }
  }
  E_Int connectPyraSize = connectPyra.size();
  if (connectPyraSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectPyraSize; i++)
      size = size + connectPyra[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Pyramids\n%ld\n", size);
#else
    fprintf(ptrFile, "Pyramids\n%d\n", size);
#endif
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 5)
      {
        FldArrayI& cp = *connectPyra[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
          fprintf(ptrFile, "%ld %ld %ld %ld %ld %ld\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), c);
#else
          fprintf(ptrFile, "%d %d %d %d %d %d\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), c);
#endif
        c++;
      }
    }
  }
  delete vertices;
  for (E_Int i = 0; i < connectTriSize; i++)
    delete connectTri[i];
  for (E_Int i = 0; i < connectQuadSize; i++)
    delete connectQuad[i];
  for (E_Int i = 0; i < connectTetraSize; i++)
    delete connectTetra[i];
  for (E_Int i = 0; i < connectHexaSize; i++)
    delete connectHexa[i];
  for (E_Int i = 0; i < connectPentaSize; i++)
    delete connectPenta[i];
  for (E_Int i = 0; i < connectPyraSize; i++)
    delete connectPyra[i];

  fclose(ptrFile);
  return 0;
}

E_Int K_IO::GenIO::meshwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<vector<E_Int> >& eltType,
  vector<char*>& zoneNames,
  vector<vector<E_Int> >* colors)
{
  E_Int nzones = unstructField.size();
  E_Int nvalidZones = 0;
  vector<E_Boolean> isZoneValid(nzones);

  for (E_Int zn = 0; zn < nzones; zn++)
  {
    isZoneValid [zn] = false;
    vector<E_Int>& eltTypeZn = eltType[zn];
    E_Int nvalidEltTypes = 0;
    for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
    {
      // triangles, quads, tetra, hexa, edges, supported
      if (eltTypeZn[ic] == 1 || eltTypeZn[ic] == 2 || eltTypeZn[ic] == 3 || 
          eltTypeZn[ic] == 4 || eltTypeZn[ic] == 5 || eltTypeZn[ic] == 6 ||
          eltTypeZn[ic] == 7) 
        nvalidEltTypes++;
      else
#ifdef E_DOUBLEINT
        printf("Warning: meshwrite: zone %ld not written (not a valid element "
               "type: %ld).", zn, eltTypeZn[ic]);
#else
        printf("Warning: meshwrite: zone %d not written (not a valid element "
               "type: %d).", zn, eltTypeZn[ic]);
#endif
    }
    if (nvalidEltTypes == (int)eltTypeZn.size())
    {
      isZoneValid[zn] = true; nvalidZones++;
    }
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: meshwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

 char format1[40]; char dataFmtl[40];
 strcpy(dataFmtl, dataFmt);
 int l = strlen(dataFmt); 
 if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format for data
  sprintf(format1,"%s%s%s", dataFmt, dataFmt, dataFmtl);
#ifdef E_DOUBLEINT
  strcat(format1," %ld\n");
#else
  strcat(format1," %d\n");
#endif

  // Connectivite par elts
  E_Int c = 0;
  vector<FldArrayI*> connectBar;
  vector<FldArrayI*> connectTri;
  vector<FldArrayI*> connectQuad;
  vector<FldArrayI*> connectTetra;
  vector<FldArrayI*> connectHexa;
  vector<FldArrayI*> connectPenta;
  vector<FldArrayI*> connectPyra;
  
  // Compute offsets and total size before the parallel block
  E_Int cmpt = 0, size = 0;
  vector<E_Int> shift(nzones);
  E_Int cmptBE[8];
  for (E_Int i = 0; i < 8; i++) cmptBE[i] = 0;
  vector<vector<E_Int> > connIdSrc(nzones);
  vector<vector<E_Int> > connIdTgt(nzones);
  for (E_Int zn = 0; zn < nzones; zn++)
  {
    shift[zn] = size;
    vector<E_Int>& eltTypeZn = eltType[zn];
    for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
    {
      E_Int elt = eltTypeZn[ic];
      connIdSrc[zn].push_back(cmpt); cmpt++;
      if (isZoneValid[zn])
      {
        if (elt == 1) connectBar.push_back(NULL);
        else if (elt == 2) connectTri.push_back(NULL);
        else if (elt == 3) connectQuad.push_back(NULL);
        else if (elt == 4) connectTetra.push_back(NULL);
        else if (elt == 7) connectHexa.push_back(NULL);
        else if (elt == 6) connectPenta.push_back(NULL);
        else if (elt == 5) connectPyra.push_back(NULL);
        connIdTgt[zn].push_back(cmptBE[elt]);
        cmptBE[elt]++;
      }
    }
    if (isZoneValid[zn]) size += unstructField[zn]->getSize();
  }

  E_Int connectBarSize = connectBar.size();
  E_Int connectTriSize = connectTri.size();
  E_Int connectQuadSize = connectQuad.size();
  E_Int connectTetraSize = connectTetra.size();
  E_Int connectPyraSize = connectPyra.size();
  E_Int connectPentaSize = connectPenta.size();
  E_Int connectHexaSize = connectHexa.size();

  // Concatenate all vertices in one field
  FldArrayF* vertices;
  vertices = new FldArrayF(size,3);
  FldArrayF& v = *vertices;

  #pragma omp parallel
  {
    for (E_Int zn = 0; zn < nzones; zn++)
    {
      if (not isZoneValid[zn]) continue;
      
      // Field
      E_Int offset = shift[zn];
      FldArrayF& field = *unstructField[zn];
      if (posx > 0 && posy > 0 && posz > 0)
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          v(offset+n,1) = field(n,posx);
          v(offset+n,2) = field(n,posy);
          v(offset+n,3) = field(n,posz);
        }
      }
      else if (posx > 0 && posy > 0)
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          v(offset+n,1) = field(n,posx);
          v(offset+n,2) = field(n,posy);
          v(offset+n,3) = 0.;
        }
      }
      else if (posx > 0 && posz > 0)
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          v(offset+n,1) = field(n,posx);
          v(offset+n,2) = 0.;
          v(offset+n,3) = field(n,posz);
        }
      }
      else if (posy > 0 && posz > 0)
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          v(offset+n,1) = 0.;
          v(offset+n,2) = field(n,posy);
          v(offset+n,3) = field(n,posz);
        }
      }
      else
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          v(offset+n,1) = field(n,posx);
          v(offset+n,2) = 0.;
          v(offset+n,3) = 0.;
        }
      }

      // Connectivities
      const vector<E_Int>& eltTypeZn = eltType[zn];
      const vector<E_Int>& connIdSrcz = connIdSrc[zn];
      const vector<E_Int>& connIdTgtz = connIdTgt[zn];

      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        E_Int elt = eltTypeZn[ic];
        FldArrayI& cn = *connect[connIdSrcz[ic]];
        
        if (elt == 1) // Edge
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
          }
          connectBar[connIdTgtz[ic]] = cpp;
        }
        else if (elt == 2) // Tri
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
            cp(n,3) = cn(n,3) + shift[zn];
          }
          connectTri[connIdTgtz[ic]] = cpp;
        }
        else if (elt == 3) // quads
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
            cp(n,3) = cn(n,3) + shift[zn];
            cp(n,4) = cn(n,4) + shift[zn];
          }
          connectQuad[connIdTgtz[ic]] = cpp;
        }
        else if (elt == 4) // tetra
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
            cp(n,3) = cn(n,3) + shift[zn];
            cp(n,4) = cn(n,4) + shift[zn];
          }
          connectTetra[connIdTgtz[ic]] = cpp;
        }
        else if (elt == 7) // hexa
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
            cp(n,3) = cn(n,3) + shift[zn];
            cp(n,4) = cn(n,4) + shift[zn];
            cp(n,5) = cn(n,5) + shift[zn];
            cp(n,6) = cn(n,6) + shift[zn];
            cp(n,7) = cn(n,7) + shift[zn];
            cp(n,8) = cn(n,8) + shift[zn];
          }
          connectHexa[connIdTgtz[ic]] = cpp;
        }
        else if (elt == 6) // penta
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
            cp(n,3) = cn(n,3) + shift[zn];
            cp(n,4) = cn(n,4) + shift[zn];
            cp(n,5) = cn(n,5) + shift[zn];
            cp(n,6) = cn(n,6) + shift[zn];
          }
          connectPenta[connIdTgtz[ic]] = cpp;
        }
        else if (elt == 5) // pyra
        {
          FldArrayI* cpp = new FldArrayI(cn);
          FldArrayI& cp = *cpp;
          #pragma omp for
          for (E_Int n = 0; n < cn.getSize(); n++)
          {
            cp(n,1) = cn(n,1) + shift[zn];
            cp(n,2) = cn(n,2) + shift[zn];
            cp(n,3) = cn(n,3) + shift[zn];
            cp(n,4) = cn(n,4) + shift[zn];
            cp(n,5) = cn(n,5) + shift[zn];
          }
          connectPyra[connIdTgtz[ic]] = cpp;
        }
      }
    }
  }

  // Ecriture
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: meshwrite: I can't open file %s.\n", file);
    return 1;  
  }

  fprintf(ptrFile, "MeshVersionFormatted\n1\n");
  fprintf(ptrFile, "Dimension\n3\n");
#ifdef E_DOUBLEINT
  fprintf(ptrFile, "Vertices\n%ld\n", v.getSize());
#else
  fprintf(ptrFile, "Vertices\n%d\n", v.getSize());
#endif
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, v(i,1), v(i,2), v(i,3), 0);
  
  if (connectBarSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectBarSize; i++)
      size = size + connectBar[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Edges\n%ld\n", size);
#else
    fprintf(ptrFile, "Edges\n%d\n", size);
#endif
    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 1)
        {
          FldArrayI& cp = *connectBar[c];
          if (!colors)
            for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
              fprintf(ptrFile, "%ld %ld %ld\n", cp(i,1), cp(i,2), c);
#else
              fprintf(ptrFile, "%d %d %d\n", cp(i,1), cp(i,2), c);
#endif
          else
            for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
              fprintf(ptrFile, "%ld %ld %ld\n", cp(i,1), cp(i,2), (*colors)[c][i]);
#else
              fprintf(ptrFile, "%d %d %d\n", cp(i,1), cp(i,2), (*colors)[c][i]);
#endif
          c++;
        }
      }
    }
  }

  if (connectTriSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectTriSize; i++)
      size = size + connectTri[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Triangles\n%ld\n", size);
#else
    fprintf(ptrFile, "Triangles\n%d\n", size);
#endif
    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 2)
        {
          FldArrayI& cp = *connectTri[c];
          if (!colors)
            for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
              fprintf(ptrFile, "%ld %ld %ld %ld\n", cp(i,1), cp(i,2), cp(i,3), c);
#else
              fprintf(ptrFile, "%d %d %d %d\n", cp(i,1), cp(i,2), cp(i,3), c);
#endif
          else
            for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
              fprintf(ptrFile, "%ld %ld %ld %ld\n", cp(i,1), cp(i,2), cp(i,3), (*colors)[c][i]);
#else
              fprintf(ptrFile, "%d %d %d %d\n", cp(i,1), cp(i,2), cp(i,3), (*colors)[c][i]);
#endif
          c++;
        }
      }
    }
  }

  if (connectQuadSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectQuadSize; i++)
      size = size + connectQuad[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Quadrilaterals\n%ld\n", size);
#else
    fprintf(ptrFile, "Quadrilaterals\n%d\n", size);
#endif
    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 3)
        {
          FldArrayI& cp = *connectQuad[c];
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld %ld\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#else
          fprintf(ptrFile, "%d %d %d %d %d\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#endif
          c++;
        }
      }
    }
  }

  if (connectTetraSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectTetraSize; i++)
      size = size + connectTetra[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Tetrahedra\n%ld\n", size);
#else
    fprintf(ptrFile, "Tetrahedra\n%d\n", size);
#endif

    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 4)
        {
          FldArrayI& cp = *connectTetra[c];
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld %ld\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#else
            fprintf(ptrFile, "%d %d %d %d %d\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
#endif
          c++;
        }
      }
    }
  }

  if (connectHexaSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectHexaSize; i++)
      size = size + connectHexa[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Hexahedra\n%ld\n", size);
#else
    fprintf(ptrFile, "Hexahedra\n%d\n", size);
#endif
    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 7)
        {
          FldArrayI& cp = *connectHexa[c];
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld %ld %ld %ld %ld %ld\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                    cp(i,5), cp(i,6), cp(i,7), cp(i,8), c);
#else
            fprintf(ptrFile, "%d %d %d %d %d %d %d %d %d\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                    cp(i,5), cp(i,6), cp(i,7), cp(i,8), c);
#endif
          c++;
        }
      }
    }
  }
  
  if (connectPentaSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectPentaSize; i++)
      size = size + connectPenta[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Prisms\n%ld\n", size);
#else
    fprintf(ptrFile, "Prisms\n%d\n", size);
#endif
    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 6)
        {
          FldArrayI& cp = *connectPenta[c];
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld %ld %ld %ld\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                    cp(i,5), cp(i,6), c);
#else
            fprintf(ptrFile, "%d %d %d %d %d %d %d\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                    cp(i,5), cp(i,6), c);
#endif
          c++;
        }
      }
    }
  }
  
  if (connectPyraSize != 0)
  {
    size = 0;
    for (E_Int i = 0; i < connectPyraSize; i++)
      size = size + connectPyra[i]->getSize();
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "Pyramids\n%ld\n", size);
#else
    fprintf(ptrFile, "Pyramids\n%d\n", size);
#endif
    c = 0;
    for (E_Int zn = 0; zn < nzones; zn++) 
    {
      if (not isZoneValid[zn]) continue;
      vector<E_Int>& eltTypeZn = eltType[zn];
      for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
      {
        if (eltTypeZn[ic] == 5)
        {
          FldArrayI& cp = *connectPyra[c];
          for (E_Int i = 0; i < cp.getSize(); i++)
#ifdef E_DOUBLEINT
            fprintf(ptrFile, "%ld %ld %ld %ld %ld %ld\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                    cp(i,5), c);
#else
            fprintf(ptrFile, "%d %d %d %d %d %d\n", 
                    cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                    cp(i,5), c);
#endif
          c++;
        }
      }
    }
  }
  
  delete vertices;
  for (E_Int i = 0; i < connectTriSize; i++)
    delete connectTri[i];
  for (E_Int i = 0; i < connectQuadSize; i++)
    delete connectQuad[i];
  for (E_Int i = 0; i < connectTetraSize; i++)
    delete connectTetra[i];
  for (E_Int i = 0; i < connectHexaSize; i++)
    delete connectHexa[i];
  for (E_Int i = 0; i < connectPentaSize; i++)
    delete connectPenta[i];
  for (E_Int i = 0; i < connectPyraSize; i++)
    delete connectPyra[i];

  fclose(ptrFile);
  return 0;
}
