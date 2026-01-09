/*    
    Copyright 2013-2026 ONERA.

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
# include "String/kstring.h"
# include <vector>
#include <unordered_map>
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

  E_Bool foundEdge = false;
  E_Bool foundTri = false;
  E_Bool foundQuad = false;
  E_Bool foundTetra = false;
  E_Bool foundHexa = false;
  E_Bool foundPenta = false;
  E_Bool foundPyra = false;

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
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }
  //printf(SF_D2_ "\n", unstructField.size(), connect.size());

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
  E_Int npts = E_Int(t);
  FldArrayF f(npts, 3);
  for (E_Int i = 0; i < npts; i++)
  {
    for (E_Int j = 1; j <= 3; j++)
    {
      res = readDouble(ptrFile, t, -1); f(i,j) = t;
    }
    res = readDouble(ptrFile, t, -1); // ref discarded
  }

  E_Int st;
  res = readWord(ptrFile, buf);

  const E_Int ncmax = 8;
  vector<FldArrayI*> tmpConnect;
  vector<E_Bool> topoFound(ncmax, false);
  
  eltType.clear();
  if (api == 3) eltType.resize(1);

  // Create a map between basic elements and element numbers
  E_Int elt, nvpeElt;
  std::unordered_map<std::string, E_Int> beMap;
  beMap["Edges"] = 1;
  beMap["Triangles"] = 2;
  beMap["Quadrilaterals"] = 3;
  beMap["Tetrahedra"] = 4;
  beMap["Pyramids"] = 5;
  beMap["Prisms"] = 6;
  beMap["Hexahedra"] = 7;

  vector<E_Int> nvpe(ncmax);
  nvpe[1] = 2; nvpe[2] = 3; nvpe[3] = 4;
  nvpe[4] = 4; nvpe[5] = 5; nvpe[6] = 6; nvpe[7] = 8;

  while (strcmp(buf, "End") != 0 && res >= 1)
  { 
    elt = beMap[buf];
    nvpeElt = nvpe[elt];
    // Connectivity `beMap[buf]` (Global) -> 1 zone
    res = readInt(ptrFile, st, -1);
    if (st != 0)
    {
      FldArrayI* cn = new FldArrayI(st, nvpeElt);
      FldArrayI& c = *cn;
      for (E_Int i = 0; i < st; i++)
      {
        for (E_Int j = 1; j <= nvpeElt; j++)
        {
          res = readInt(ptrFile, ti, -1); c(i,j) = ti;
        }
        res = readInt(ptrFile, ti, -1); // discarded
      }
      
      if (api == 3)
      {
        eltType[0].push_back(elt);
        tmpConnect.push_back(cn);
      }
      else
      {
        eltType.push_back({elt});
        connect.push_back(cn);
      }
      topoFound[elt] = true;
    }
    res = readWord(ptrFile, buf);
  }

  // Cree la varString
  varString = new char [8];
  strcpy(varString, "x,y,z");
  
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
    for (size_t i = 1; i < topoFound.size(); i++)
    { 
      if (topoFound[i])
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
    sprintf(zoneName, "Zone%zu", i);
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
      printf("Warning: meshwrite: zone " SF_D_ " not written (not a valid "
             "element type: " SF_D_ ").", zone, eltType[zone]);
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
  strcat(format1," " SF_D_ "\n");

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
  fprintf(ptrFile, "Vertices\n" SF_D_ "\n", v.getSize());
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, v(i,1), v(i,2), v(i,3), 0);
  
  E_Int connectEdgeSize = connectEdge.size();
  E_Int c;
  if (connectEdgeSize != 0)
  {
    E_Int size = 0;
    for (E_Int i = 0; i < connectEdgeSize; i++)
      size = size + connectEdge[i]->getSize();
    fprintf(ptrFile, "Edges\n" SF_D_ "\n", size);
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 1)
      {
        FldArrayI& cp = *connectEdge[c];
        if (!colors)
          for (E_Int i = 0; i < cp.getSize(); i++)
            fprintf(ptrFile, SF_D3_ "\n", cp(i,1), cp(i,2), c);
        else
          for (E_Int i = 0; i < cp.getSize(); i++)
            fprintf(ptrFile, SF_D3_ "\n", cp(i,1), cp(i,2), (*colors)[c][i]);
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
    fprintf(ptrFile, "Triangles\n" SF_D_ "\n", size);
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 2)
      {
        FldArrayI& cp = *connectTri[c];
        if (!colors)
          for (E_Int i = 0; i < cp.getSize(); i++)
            fprintf(ptrFile, SF_D4_ "\n", cp(i,1), cp(i,2), cp(i,3), c);
        else
          for (E_Int i = 0; i < cp.getSize(); i++)
            fprintf(ptrFile, SF_D4_ "\n", cp(i,1), cp(i,2), cp(i,3), (*colors)[c][i]);
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
    fprintf(ptrFile, "Quadrilaterals\n" SF_D_ "\n", size);
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 3)
      {
        FldArrayI& cp = *connectQuad[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
        fprintf(ptrFile, SF_D5_ "\n", cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
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
    fprintf(ptrFile, "Tetrahedra\n" SF_D_ "\n", size);

    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 4)
      {
        FldArrayI& cp = *connectTetra[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
          fprintf(ptrFile, SF_D5_ "\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4), c);
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
    fprintf(ptrFile, "Hexahedra\n" SF_D_ "\n", size);
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 7)
      {
        FldArrayI& cp = *connectHexa[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
          fprintf(ptrFile, SF_D9_ "\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), cp(i,6), cp(i,7), cp(i,8), c);
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
    fprintf(ptrFile, "Prisms\n" SF_D_ "\n", size);
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 6)
      {
        FldArrayI& cp = *connectPenta[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
          fprintf(ptrFile, SF_D7_ "\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), cp(i,6), c);
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
    fprintf(ptrFile, "Pyramids\n" SF_D_ "\n", size);
    c = 0;
    for (E_Int i = 0; i < nzone; i++) 
    {
      if (eltType[i] == 5)
      {
        FldArrayI& cp = *connectPyra[c];
        for (E_Int i = 0; i < cp.getSize(); i++)
          fprintf(ptrFile, SF_D6_ "\n", 
                  cp(i,1), cp(i,2), cp(i,3), cp(i,4),
                  cp(i,5), c);
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
  // Get number of valid zones, ie, zones containing element types that are all
  // valid. NB: this format supports one zone only
  E_Int nzones = unstructField.size();
  E_Int nvalidZones = 0;
  E_Int zoneId = -1;

  for (E_Int zn = 0; zn < nzones; zn++)
  {
    vector<E_Int>& eltTypeZn = eltType[zn];
    E_Int nvalidEltTypes = 0;
    for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
    {
      // All 1D, 2D, and 3D basic elements are supported
      if (eltTypeZn[ic] >= 1 && eltTypeZn[ic] <= 7) nvalidEltTypes++;
      else
        printf("Warning: meshwrite: zone " SF_D_ " not written (not a valid element "
               "type: " SF_D_ ").", zn, eltTypeZn[ic]);
    }
    if (nvalidEltTypes == (E_Int)eltTypeZn.size())
    {
      nvalidZones++;
      if (zoneId == -1) zoneId = zn;
    }
  }

  if (nvalidZones == 0) return 1;
  else if (nvalidZones > 1)
    printf("Warning: meshwrite: monozone format, only the first valid zone "
           "will be written: zone #" SF_D_ ".\n", zoneId+1);

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: meshwrite: zone does not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Build format for data
  char format1[40]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';
  
  sprintf(format1, "%s%s%s", dataFmt, dataFmt, dataFmtl);
  strcat(format1, " " SF_D_ "\n");

  // Connectivite par elts
  E_Int npts = unstructField[zoneId]->getSize();
  vector<FldArrayI*> connectBE(8, NULL);
  vector<E_Int> nvpe(8);
  nvpe[1] = 2; nvpe[2] = 3; nvpe[3] = 4;
  nvpe[4] = 4; nvpe[5] = 5; nvpe[6] = 6; nvpe[7] = 8;

  // Concatenate all vertices in one field
  FldArrayF* vertices;
  vertices = new FldArrayF(npts,3);
  FldArrayF& v = *vertices;
  vector<E_Int> posCoords; posCoords.reserve(6);
  if (posx > 0) {posCoords.push_back(1); posCoords.push_back(posx);}
  if (posy > 0) {posCoords.push_back(2); posCoords.push_back(posy);}
  if (posz > 0) {posCoords.push_back(3); posCoords.push_back(posz);}

  #pragma omp parallel
  {
    E_Int ind1, ind2;
    
    // Field
    FldArrayF& field = *unstructField[zoneId];
    if (posx > 0 && posy > 0 && posz > 0)
    {
      #pragma omp for
      for (E_Int n = 0; n < field.getSize(); n++)
      {
        v(n,1) = field(n,posx);
        v(n,2) = field(n,posy);
        v(n,3) = field(n,posz);
      }
    }
    else
    {
      #pragma omp for
      for (E_Int n = 0; n < field.getSize(); n++)
        for (E_Int j = 1; j <= 3; j++)
          v(n,j) = 0.;

      #pragma omp for
      for (E_Int n = 0; n < field.getSize(); n++)
        for (size_t j = 0; j < posCoords.size(); j+=2)
        {
          ind1 = posCoords[j];
          ind2 = posCoords[j+1];
          v(n,ind1) = field(n,ind2);
        }
    }
  }

  // Connectivities
  const vector<E_Int>& eltTypeZn = eltType[zoneId];
  for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
  {
    E_Int elt = eltTypeZn[ic];
    FldArrayI& cn = *connect[zoneId]->getConnect(ic);
    FldArrayI* cpp = new FldArrayI(cn);
    FldArrayI& cp = *cpp;
    for (E_Int n = 0; n < cn.getSize(); n++)
      for (E_Int j = 1; j <= nvpe[elt]; j++)
        cp(n,j) = cn(n,j);
    connectBE[elt] = cpp;
  }

  // Ecriture
  E_Int c = 0;
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: meshwrite: can't open file %s.\n", file);
    return 1;  
  }

  fprintf(ptrFile, "MeshVersionFormatted\n1\n");
  fprintf(ptrFile, "Dimension\n3\n");
  fprintf(ptrFile, "Vertices\n" SF_D_ "\n", v.getSize());
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, v(i,1), v(i,2), v(i,3), 0);

  const char* eltTypeMesh[8];
  eltTypeMesh[1] = "Edges";
  eltTypeMesh[2] = "Triangles";
  eltTypeMesh[3] = "Quadrilaterals";
  eltTypeMesh[4] = "Tetrahedra";
  eltTypeMesh[5] = "Pyramids";
  eltTypeMesh[6] = "Prisms";
  eltTypeMesh[7] = "Hexahedra";
  
  for (size_t elt = 1; elt < connectBE.size(); elt++)
  {
    if (connectBE[elt] == NULL) continue;

    E_Int size = connectBE[elt]->getSize();
    fprintf(ptrFile, "%s\n" SF_D_ "\n", eltTypeMesh[elt], size);

    FldArrayI& cp = *connectBE[elt];
    if (!colors)
    {
      for (E_Int i = 0; i < cp.getSize(); i++)
      {
        for (E_Int j = 1; j <= nvpe[elt]; j++)
          fprintf(ptrFile, SF_D_ " ", cp(i,j));
        fprintf(ptrFile, SF_D_ "\n", c);
      }
    }
    else
    {
      for (E_Int i = 0; i < cp.getSize(); i++)
      {
        for (E_Int j = 1; j <= nvpe[elt]; j++)
          fprintf(ptrFile, SF_D_ " ", cp(i,j));
        fprintf(ptrFile, SF_D_ "\n", (*colors)[c][i]);
      }
    }
    c++;
  }
  
  delete vertices;
  for (size_t i = 0; i < connectBE.size(); i++)
    if (connectBE[i] != NULL) delete connectBE[i];
  connectBE.clear();

  fclose(ptrFile);
  return 0;
}
