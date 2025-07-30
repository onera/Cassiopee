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

// Formated gmsh (2.2) file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include <vector>
# include <unordered_set>
# include <unordered_map>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* gmshread
*/
//=============================================================================
E_Int K_IO::GenIO::gmshread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t; E_Int ti; E_Int type;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: gmshread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de l'entete
  char buf[BUFSIZE];
  res = readWord(ptrFile, buf);
  if (res == -1) { fclose(ptrFile); return 1; }
  if (strcmp(buf, "$MeshFormat") != 0) { fclose(ptrFile); return 1; }
  res = readDouble(ptrFile, t, -1);
  //printf("version " SF_F_ "\n", t);
  res = readInt(ptrFile, ti, -1); type = ti;
  if (type == 1) return 1; // c'est un binary file
  //printf("file type " SF_D_ "\n", ti);
  res = readInt(ptrFile, ti, -1);
  //printf("data size " SF_D_ "\n", ti);
  res = readGivenKeyword(ptrFile, "$ENDMESHFORMAT");
  
  // Lecture Vertices (Global)
  res = readGivenKeyword(ptrFile, "$NODES");
  res = readInt(ptrFile, ti, -1);
  E_Int nb = E_Int(ti);
  //printf("Number of nodes " SF_D_ "\n", nb);
  FldArrayF f(nb, 3);
  E_Float* f1 = f.begin(1); E_Float* f2 = f.begin(2); E_Float* f3 = f.begin(3);
  FldArrayI indirNodes(nb);
  for (E_Int i = 0; i < nb; i++)
  {
    res = readInt(ptrFile, ti, -1); indirNodes[i] = ti;
    res = readDouble(ptrFile, t, -1); f1[i] = t;
    res = readDouble(ptrFile, t, -1); f2[i] = t;
    res = readDouble(ptrFile, t, -1); f3[i] = t;
    //printf(SF_F3_ "\n", f(i,1), f(i,2), f(i,3));
  }
  //res = readGivenKeyword(ptrFile, "$ENDNODES"); // pas obligatoire?

  /* Elements by zone type */
  res = readGivenKeyword(ptrFile, "$ELEMENTS");
  res = readInt(ptrFile, ti, -1);
  E_Int ne = E_Int(ti); // Global
  printf("Number of elements " SF_D_ "\n", ne);
  FldArrayI indirElements(ne);
  // declarations
#include "GenIO_gmsh3.h"
  E_Int nDiscard = 0;

  E_Int tagl, ind;
  /* Compte les elements par type */
  E_LONG pos = KFTELL(ptrFile);
#define READI readInt(ptrFile, ti, -1)
#include "GenIO_gmsh1.h"
printf("Elements BAR=" SF_D_ " TRI=" SF_D_ " QUAD=" SF_D_ " TETRA=" SF_D_ " HEXA=" SF_D_ " PENTA=" SF_D_ " PYRA=" SF_D_ " NODES=" SF_D_ "\n", 
       nBAR, nTRI, nQUAD, nTETRA, nHEXA, nPENTA, nPYRA, nNODE);
printf("Elements BAR_3=" SF_D_ " TRI_6=" SF_D_ " QUAD_9=" SF_D_ " TETRA_10=" SF_D_ " HEXA_27=" SF_D_ " PENTA_18=" SF_D_ " PYRA_14=" SF_D_ "\n",
       nBAR_3, nTRI_6, nQUAD_9, nTETRA_10, nHEXA_27, nPENTA_18, nPYRA_14);

  /* Allocations */
  E_Bool fo = true;
#include "GenIO_gmsh4.h"
  
  /* Lecture reelle des elements par type */
  KFSEEK(ptrFile, pos, SEEK_SET);
#include "GenIO_gmsh2.h"

  /* Nodes duplications */
#include "GenIO_gmsh5.h"

  // Cree le nom des zones
  //printf("Number of zones " SF_D_ "\n", unstructField.size());
  for (size_t i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }
  //printf("sizes: " SF_D2_ "\n", unstructField.size(), connect.size());

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);

  return 0;
}

E_Int K_IO::GenIO::gmshread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<vector<E_Int> >& eltType, vector<char*>& zoneNames, E_Int api)
{
  E_Int res;
  E_Float t; E_Int ti; E_Int type;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: gmshread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de l'entete
  char buf[BUFSIZE];
  res = readWord(ptrFile, buf);
  if (res == -1) { fclose(ptrFile); return 1; }
  if (strcmp(buf, "$MeshFormat") != 0) { fclose(ptrFile); return 1; }
  res = readDouble(ptrFile, t, -1);
  res = readInt(ptrFile, ti, -1); type = ti;
  if (type == 1) return 1; // c'est un binary file
  res = readInt(ptrFile, ti, -1);
  res = readGivenKeyword(ptrFile, "$ENDMESHFORMAT");
  
  // Lecture Vertices (Global)
  res = readGivenKeyword(ptrFile, "$NODES");
  res = readInt(ptrFile, ti, -1);
  E_Int npts = E_Int(ti);
  FldArrayF f(npts, 3);
  E_Float* f1 = f.begin(1); E_Float* f2 = f.begin(2); E_Float* f3 = f.begin(3);
  FldArrayI indirNodes(npts);
  for (E_Int i = 0; i < npts; i++)
  {
    res = readInt(ptrFile, ti, -1); indirNodes[i] = ti;
    res = readDouble(ptrFile, t, -1); f1[i] = t;
    res = readDouble(ptrFile, t, -1); f2[i] = t;
    res = readDouble(ptrFile, t, -1); f3[i] = t;
  }
  res = readGivenKeyword(ptrFile, "$ENDNODES");

  /* Elements by zone type */
  res = readGivenKeyword(ptrFile, "$ELEMENTS");
  res = readInt(ptrFile, ti, -1);
  E_Int ne = E_Int(ti); // Global
  //printf("Total number of elements (all zones and connectivities): " SF_D_ "\n", ne);
  FldArrayI indirElements(ne);

  // Declarations
  // Correspondance between Gmsh element type and (CGNS elttype, number of vertex per element)
  std::unordered_map<E_Int, std::pair<E_Int, E_Int> > beGmsh;
  beGmsh[1] = std::make_pair(1, 2); // BAR
  beGmsh[2] = std::make_pair(2, 3); // TRI
  beGmsh[3] = std::make_pair(3, 4); // QUAD
  beGmsh[4] = std::make_pair(4, 4); // TETRA
  beGmsh[5] = std::make_pair(7, 8); // HEXA
  beGmsh[6] = std::make_pair(6, 6); // PENTA
  beGmsh[7] = std::make_pair(5, 5); // PYRA
  beGmsh[8] = std::make_pair(-1, 3); // 3-node BAR (second order)
  beGmsh[9] = std::make_pair(-1, 6); // 6-node TRI (second order)
  beGmsh[10] = std::make_pair(-1, 9); // 9-node QUAD (second order)
  beGmsh[11] = std::make_pair(-1, 10); // 10-node TETRA (second order)
  beGmsh[12] = std::make_pair(-1, 27); // 27-node HEXA (second order)
  beGmsh[13] = std::make_pair(-1, 18); // 18-node PENTA (second order)
  beGmsh[14] = std::make_pair(-1, 14); // 14-node PYRA (second order)
  beGmsh[15] = std::make_pair(0, 1); // NODE
  beGmsh[16] = std::make_pair(-1, 8); // 8-node QUAD (second order)
  beGmsh[17] = std::make_pair(-1, 20); // 20-node HEXA (second order)
  beGmsh[18] = std::make_pair(-1, 15); // 15-node PENTA (second order)
  beGmsh[19] = std::make_pair(-1, 13); // 13-node PYRA (second order)
  beGmsh[20] = std::make_pair(-1, 9); // 9-node TRI (third order)
  beGmsh[21] = std::make_pair(-1, 10); // 10-node TRI (third order)
  beGmsh[22] = std::make_pair(-1, 12); // 12-node TRI (fourth order)
  beGmsh[23] = std::make_pair(-1, 15); // 15-node TRI (fourth order)
  beGmsh[24] = std::make_pair(-1, 15); // 15-node incomplete TRI (fifth order) -> pas en CGNS
  beGmsh[25] = std::make_pair(-1, 21); // 21-node TRI (fifth order) -> pas en CGNS
  beGmsh[26] = std::make_pair(-1, 4); // 4-node BAR (third order)
  beGmsh[27] = std::make_pair(-1, 5); // 5-node BAR (fourth order)
  beGmsh[28] = std::make_pair(-1, 6); // 6-node BAR (fifth order)-> pas en CGNS
  beGmsh[29] = std::make_pair(-1, 20); // 20-node TETRA (third order)
  beGmsh[30] = std::make_pair(-1, 35); // 35-node TETRA (fourth order)
  beGmsh[31] = std::make_pair(-1, 56); // 56-node TETRA (fifth order)-> pas en CGNS
  beGmsh[92] = std::make_pair(-1, 64); // 64-node HEXA (third order)
  beGmsh[93] = std::make_pair(-1, 125); // 125-node HEXA (fourth order)

  std::unordered_set<E_Int> uniqueZoneSet; // unique zone indices 
  std::unordered_map<std::pair<E_Int, E_Int>, E_Int, pairHash> neltsBEMap; // number of elements of that BE in zone

  E_Int etGmsh, tagl, tagzn, ind;
  /* Compte les elements par type */
  E_LONG pos = KFTELL(ptrFile);
#define READI readInt(ptrFile, ti, -1)

  for (E_Int i = 0; i < ne; i++)
  {
#ifdef BINARY
    READI; etGmsh = ti; // type d'element
    READI;  // follow
    READI; tagl = ti; // ntag
    READI; // Indirection
    for (E_Int j = 0; j < tagl-1; j++) READI; // ntags
    READI; tagzn = ti; // zone index
#else
    READI; // Indirection
    READI; etGmsh = ti; // type d'element
    READI; tagl = ti; // tags
    for (E_Int j = 0; j < tagl-1; j++) READI; // ntags
    READI; tagzn = ti; // zone index
#endif
    uniqueZoneSet.insert(tagzn);
    neltsBEMap[std::make_pair(tagzn, etGmsh)]++;
    for (E_Int j = 0; j < beGmsh[etGmsh].second; j++) READI;
  }
  
  const E_Int nzones = uniqueZoneSet.size();
  vector<vector<E_Int> > nepc(nzones); // number of elements of each type per zone
  eltType.clear(); eltType.resize(nzones);

  E_Int cmpt = 0;
  vector<E_Int> icCmpt(nzones, 0);
  std::unordered_map<E_Int, E_Int> zoneIdMap;
  vector<std::unordered_map<E_Int, E_Int> > connIdMap(nzones);

  for (const auto& entry : neltsBEMap)
  {
    E_Int zn = entry.first.first;
    auto res = zoneIdMap.insert(std::make_pair(zn, -1));
    if (res.first->second == -1)
    {
      res.first->second = cmpt; cmpt++;
    }
    zn = res.first->second;
    etGmsh = entry.first.second;
    eltType[zn].push_back(beGmsh[etGmsh].first);
    nepc[zn].push_back(entry.second);

    connIdMap[zn].insert(std::make_pair(etGmsh, icCmpt[zn]));
    icCmpt[zn]++;
  }

  /* Allocations */
  varString = new char [8];
  strcpy(varString, "x,y,z");

  for (E_Int zn = 0; zn < nzones; zn++)
  {
    if (eltType[zn][0] == 0) // NODE
    {
      connect.push_back(new FldArrayI());
    }
    else
    {
      char eltString[256]; vector<E_Int> dummy(1);
      K_ARRAY::typeId2eltString(eltType[zn], 0, eltString, dummy);
      
      PyObject* tpl = K_ARRAY::buildArray3(3, varString, npts, nepc[zn],
                                          eltString, 0, api);
      FldArrayI* cn2; FldArrayF* f2;
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      connect.push_back(cn2);
      delete f2;
    }
  }

  /* Lecture reelle des elements par type */
  vector<vector<E_Int> > etCmpt(nzones);
  for (E_Int zn = 0; zn < nzones; zn++) etCmpt[zn].resize(icCmpt[zn], 0);

  KFSEEK(ptrFile, pos, SEEK_SET);
  for (E_Int i = 0; i < ne; i++)
  {
#ifdef BINARY
    READI; etGmsh = ti; // type d'element
    READI;  // follow
    READI; tagl = ti; // ntag
    READI; // Indirection
    for (E_Int j = 0; j < tagl-1; j++) READI; // ntags
    READI; tagzn = ti; // zone index
#else
    READI; // Indirection
    READI; etGmsh = ti; // type d'element
    READI; tagl = ti; // tags
    for (E_Int j = 0; j < tagl-1; j++) READI; // ntags
    READI; tagzn = ti; // zone index
#endif
    E_Int zn = zoneIdMap[tagzn];
    E_Int ic = connIdMap[zn][etGmsh];
    E_Int nvpe = beGmsh[etGmsh].second;
    ind = etCmpt[zn][ic];
    for (E_Int j = 1; j <= nvpe; j++)
      { READI; (*connect[zn]->getConnect(ic))(ind,j) = indirNodes[ti-1]; }
    etCmpt[zn][ic]++;
  }

  for (E_Int zn = 0; zn < nzones; zn++)
  {
    // Cree le nom de la zone
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, zn);
    zoneNames.push_back(zoneName);
    
    if (eltType[zn][0] == 0) // NODE
    {
      FldArrayF* an = new FldArrayF(nepc[zn][0],3);
      for (E_Int i = 0; i < nepc[zn][0]; i++) 
      {
        ind = (*connect[zn]->getConnect(0))(i,1)-1;
        ind = indirNodes[ind]-1;
        (*an)(i,1) = f1[ind]; (*an)(i,2) = f2[ind]; (*an)(i,3) = f3[ind];
      }
      unstructField.push_back(an);
    }
    else
    {
      FldArrayF* an = new FldArrayF(f);
      unstructField.push_back(an);
    }
  }
  
  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Only write NODE, BAR, TRI, QUAD, TETRA, HEXA, PYRA, PENTA.
//=============================================================================
E_Int K_IO::GenIO::gmshwrite(
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
      printf("Warning: gmshwrite: zone " SF_D_ " not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: gmshwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;
    
  char format1[40]; char format2[40]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format for data
  strcpy(format1, SF_D_ " ");
  sprintf(format2,"%s%s%s", dataFmt, dataFmt, dataFmtl);
  strcat(format1, format2);
  strcat(format1, "\n");

  // Concatenate all vertices in one field
  E_Int size = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    size += unstructField[i]->getSize();
  }
  FldArrayF* vertices = new FldArrayF(size, 3);
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

  // Ecriture
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL)
  {
    printf("Warning: gmshwrite: I can't open file %s.\n", file);
    return 1;
  }

  fprintf(ptrFile, "$MeshFormat\n");
  fprintf(ptrFile, "2.2 0 8\n");
  fprintf(ptrFile, "$EndMeshFormat\n");

  fprintf(ptrFile, "$Nodes\n");
  fprintf(ptrFile, SF_D_ "\n", v.getSize());
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, i+1, v(i,1), v(i,2), v(i,3));
  fprintf(ptrFile, "$EndNodes\n");

  fprintf(ptrFile, "$Elements\n");
  size = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    size += connect[i]->getSize();
  }
  fprintf(ptrFile, SF_D_ "\n", size);

  // Connectivite par elts par rapport a la definition globale des vertices
  E_Int shift = 0;
  E_Int ind = 1;
  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayI& cn = *connect[i];
    E_Int elt = eltType[i];
    switch(elt)
    {
      case 0: // NODE
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 15 2 0 14 " SF_D_ "\n", ind, ind);
          ind++;
        }
      }
      break;
      case 1: // BAR
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 1 2 0 14 " SF_D2_ "\n", ind, cn(n,1)+shift, cn(n,2)+shift);
          ind++;
        }
      }
      break;
      case 2: // TRI
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 2 2 0 14 " SF_D3_ "\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift);
          ind++;
        }
      }
      break;
      case 3: // QUAD
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 3 2 0 14 " SF_D4_ "\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,cn(n,4)+shift);
          ind++;
        }
      }
      break;
      case 4: // TETRA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 4 2 0 14 " SF_D4_ "\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,cn(n,4)+shift);
          ind++;
        }
      }
      break;
      case 5: // PYRA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 7 2 0 14 " SF_D5_ "\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,
                  cn(n,4)+shift,cn(n,5)+shift);
          ind++;
        }
      }
      break;
      case 6: // PENTA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 6 2 0 14 " SF_D6_ "\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,
                  cn(n,4)+shift,cn(n,5)+shift,cn(n,6)+shift);
          ind++;
        }
      }
      break;
      case 7: // HEXA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, SF_D_ " 5 2 0 14 " SF_D8_ "\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,
                  cn(n,4)+shift,cn(n,5)+shift,cn(n,6)+shift,
                  cn(n,7)+shift,cn(n,8)+shift);
          ind++;
        }
      }
      break;
    }
    shift += unstructField[i]->getSize();
  }
  fprintf(ptrFile, "$EndElements\n");

  fclose(ptrFile);
  return 0;
}

E_Int K_IO::GenIO::gmshwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<vector<E_Int> >& eltType,
  vector<char*>& zoneNames)
{
  E_Int nzones = unstructField.size();
  E_Int nvalidZones = 0;
  vector<E_Bool> isZoneValid(nzones, false);
  for (E_Int zn = 0; zn < nzones; zn++)
  {
    vector<E_Int>& eltTypeZn = eltType[zn];
    E_Int nvalidEltTypes = 0;
    // All 0D-3D linear basic elements are supported
    for (size_t ic = 0; ic < eltTypeZn.size(); ic++)
    {
      if (eltTypeZn[ic] < 8) nvalidEltTypes++;
      else
        printf("Warning: gmshwrite: zone " SF_D_ " not written (not a valid element "
               "type: " SF_D_ ").", zn, eltTypeZn[ic]);
    }
    if (nvalidEltTypes == (E_Int)eltTypeZn.size())
      { nvalidZones++; isZoneValid[zn] = true; }
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: gmshwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;
    
  char format1[40]; char format2[40]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format for data
  strcpy(format1, SF_D_ " ");
  sprintf(format2,"%s%s%s", dataFmt, dataFmt, dataFmtl);
  strcat(format1, format2);
  strcat(format1, "\n");

  // Connectivite par elts
  E_Int npts = 0;
  vector<E_Int> nvpe(8);
  nvpe[0] = 1; nvpe[1] = 2; nvpe[2] = 3; nvpe[3] = 4;
  nvpe[4] = 4; nvpe[5] = 5; nvpe[6] = 6; nvpe[7] = 8;

  // Concatenate all vertices in one field
  vector<E_Int> voffsets(nzones+1, 0);
  for (E_Int zn = 0; zn < nzones; zn++)
  {
    if (not isZoneValid[zn]) continue;
    voffsets[zn+1] = voffsets[zn] + unstructField[zn]->getSize();
  }
  npts = voffsets[nzones];

  FldArrayF* vertices;
  vertices = new FldArrayF(npts,3);
  FldArrayF& v = *vertices;
  vector<E_Int> posCoords; posCoords.reserve(6);
  if (posx > 0) {posCoords.push_back(1); posCoords.push_back(posx);}
  if (posy > 0) {posCoords.push_back(2); posCoords.push_back(posy);}
  if (posz > 0) {posCoords.push_back(3); posCoords.push_back(posz);}

  #pragma omp parallel
  {
    E_Int ind1, ind2, ind3;
    
    // Field
    for (E_Int zn = 0; zn < nzones; zn++)
    {
      if (not isZoneValid[zn]) continue;
      FldArrayF& field = *unstructField[zn];
      if (posx > 0 && posy > 0 && posz > 0)
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          ind3 = voffsets[zn] + n;
          v(ind3,1) = field(n,posx);
          v(ind3,2) = field(n,posy);
          v(ind3,3) = field(n,posz);
        }
      }
      else
      {
        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          ind3 = voffsets[zn] + n;
          for (E_Int j = 1; j <= 3; j++)
            v(ind3,j) = 0.;
        }

        #pragma omp for
        for (E_Int n = 0; n < field.getSize(); n++)
        {
          ind3 = voffsets[zn] + n;
          for (size_t j = 0; j < posCoords.size(); j+=2)
          {
            ind1 = posCoords[j];
            ind2 = posCoords[j+1];
            v(ind3,ind1) = field(n,ind2);
          }
        }
      }
    }
  }

  // Ecriture
  vector<E_Int> eltNoGmsh(8);
  eltNoGmsh[0] = 15; eltNoGmsh[1] = 1; eltNoGmsh[2] = 2; eltNoGmsh[3] = 3;
  eltNoGmsh[4] = 4; eltNoGmsh[5] = 7; eltNoGmsh[6] = 6; eltNoGmsh[7] = 5;

  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL)
  {
    printf("Warning: gmshwrite: can't open file %s.\n", file);
    return 1;
  }

  fprintf(ptrFile, "$MeshFormat\n");
  fprintf(ptrFile, "2.2 0 8\n");
  fprintf(ptrFile, "$EndMeshFormat\n");

  fprintf(ptrFile, "$Nodes\n");
  fprintf(ptrFile, SF_D_ "\n", v.getSize());
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, i+1, v(i,1), v(i,2), v(i,3));
  fprintf(ptrFile, "$EndNodes\n");

  fprintf(ptrFile, "$Elements\n");
  E_Int nelts = 0;
  for (E_Int zn = 0; zn < nzones; zn++)
  {
    if (not isZoneValid[zn]) continue;
    E_Int nc = connect[zn]->getNConnect();
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *connect[zn]->getConnect(ic);
      nelts += cm.getSize();
    }
  }
  fprintf(ptrFile, SF_D_ "\n", nelts);

  // Connectivite par elts par rapport a la definition globale des vertices
  E_Int shift = 0;
  E_Int ind = 1;
  E_Int elt, nvpeElt, eltNo;
  for (E_Int zn = 0; zn < nzones; zn++)
  {
    if (not isZoneValid[zn]) continue;
    vector<E_Int>& eltTypeZn = eltType[zn];
    E_Int nc = connect[zn]->getNConnect();
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *connect[zn]->getConnect(ic);
      elt = eltTypeZn[ic];
      nvpeElt = nvpe[elt];
      eltNo = eltNoGmsh[elt];
      nelts = cm.getSize();

      if (elt == 0)
      {
        for (E_Int n = 0; n < nelts; n++)
        {
          fprintf(ptrFile, SF_D2_ " 2 0 " SF_D2_ "\n", ind, eltNo, zn, ind);
          ind++;
        }
      }
      else
      {
        for (E_Int n = 0; n < nelts; n++)
        {
          fprintf(ptrFile, SF_D2_ " 2 0 " SF_D_, ind, eltNo, zn);
          for (E_Int j = 1; j <= nvpeElt; j++)
            fprintf(ptrFile, " " SF_D_, cm(n,j) + shift);
          fprintf(ptrFile, "\n");
          ind++;
        }
      }
    }
    shift += unstructField[zn]->getSize();
  }
  fprintf(ptrFile, "$EndElements\n");

  fclose(ptrFile);
  return 0;
}
