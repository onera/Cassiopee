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

// Binary gmsh (2.2) file support

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
/* bingmshread
*/
//=============================================================================
E_Int K_IO::GenIO::bingmshread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  int ti; double v[3];
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: gmshread: cannot open file %s.\n", file);
    return 1;
  }

  E_Int sizeChar = sizeof(char);
  E_Int sizeInt = sizeof(int);
  E_Int sizeFloat = sizeof(double);
  int type;

  // Lecture de l'entete (12 chars)
  char buf[256];
  fread(buf, sizeChar, 12, ptrFile);
  buf[11] = '\0';
  //printf("%s\n", buf);
  if (strcmp(buf, "$MeshFormat") != 0) { fclose(ptrFile); return 1; }

  /* version , ... en chars */
  E_Int i; char c;
  i = 0; c = '1';
  while (c != ' ' && c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  //printf("version %s\n", buf);
  i = 0; c = '1';
  while (c != ' ' && c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  type = atoi(buf);
  if (type == 0) return 1; // formatted
  //printf("file type %s\n", buf);
  i = 0; c = '1';
  while (c != ' ' && c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  //printf("data size %s\n", buf);
  /* one for endianess */
  fread(&ti, sizeInt, 1, ptrFile);
  E_Boolean changeEndian = true;
  if (ti == 1) changeEndian = false; // TO DO
  fread(buf, sizeChar, 16, ptrFile); /* \n + $endMeshFormat */
  buf[16] = '\0';

  // Lecture Vertices (Global)
  fread(buf, sizeChar, 7, ptrFile); /* $Nodes */
  buf[7] = '\0';

  /* Nombre de noeuds , ... en chars */
  i = 0; c = '1';
  while (c != ' ' && c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  E_Int nb = atoi(buf);
  //printf("Number of nodes %d\n", nb);
  FldArrayF f(nb, 3);
  E_Float* f1 = f.begin(1); E_Float* f2 = f.begin(2); E_Float* f3 = f.begin(3);
  FldArrayI indirNodes(nb);
  if (changeEndian == false)
  {
    for (E_Int i = 0; i < nb; i++)
    {
      fread(&ti, sizeInt, 1, ptrFile); indirNodes[i] = ti;
      fread(v, sizeFloat, 3, ptrFile); 
      f1[i] = v[0]; f2[i] = v[1]; f3[i] = v[2];
      //printf("%f %f %f\n", f(i,1), f(i,2), f(i,3));
    }
  }
  else
  {
    for (E_Int i = 0; i < nb; i++)
    {
      fread(&ti, sizeInt, 1, ptrFile); indirNodes[i] = IBE(ti);
      fread(v, sizeFloat, 3, ptrFile); 
      f1[i] = DBE(v[0]); f2[i] = DBE(v[1]); f3[i] = DBE(v[2]);
      //printf("%f %f %f\n", f(i,1), f(i,2), f(i,3));
    }
  }
  fread(buf, sizeChar, 11, ptrFile); /* \n + $endNodes */
  buf[11] = '\0';

  /* Elements by zone type */
  fread(buf, sizeChar, 10, ptrFile); /* $elements */
  buf[10] = '\0';
  i = 0; c = '1';
  while (c != ' ' && c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  //printf("%s\n", buf);
  E_Int ne = atoi(buf); // Global
  //printf("Number of elements %d\n", ne);
  FldArrayI indirElements(ne);
  // declarations
#include "GenIO_gmsh3.h"
  
  E_Int nDiscard = 0;

  E_Int tagl, ind;
  /* Compte les elements par type */
  E_LONG pos = KFTELL(ptrFile);
#define BINARY
  if (changeEndian == false)
  {
#define READI fread(&ti, sizeInt, 1, ptrFile)
#include "GenIO_gmsh1.h"
  }
  else
  {
#undef READI
#define READI { fread(&ti, sizeInt, 1, ptrFile); ti = IBE(ti); }
#include "GenIO_gmsh1.h"
  }

  //printf("Elements BAR=%d TRI=%d QUAD=%d TETRA=%d HEXA=%d NODES=%d\n", 
  //       nBAR, nTRI, nQUAD, nTETRA, nHEXA, nNODE);

  /* Allocations */
  E_Boolean fo = true;
#include "GenIO_gmsh4.h"
  
  /* Lecture reelle des elements par type */
  KFSEEK(ptrFile, pos, SEEK_SET);
  if (changeEndian == false)
  {
#undef READI
#define READI fread(&ti, sizeInt, 1, ptrFile)
#include "GenIO_gmsh2.h"
  }
  else
  {
#undef READI
#define READI { fread(&ti, sizeInt, 1, ptrFile); ti = IBE(ti); }
#include "GenIO_gmsh2.h"
  }
  /* Nodes duplications */
#include "GenIO_gmsh5.h"

  // Cree le nom des zones
  //printf("Number of zones %d\n", unstructField.size());
  for (size_t i = 0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }
  //printf("sizes: %d %d\n", unstructField.size(), connect.size());

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);

  return 0;
}

//=============================================================================
// Only write NODE, BAR, TRI, QUAD, TETRA, HEXA, PYRA, PENTA.
//=============================================================================
E_Int K_IO::GenIO::bingmshwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // NODE, BAR, TRI, QUADS, TETRA, HEXA , PENTA, PYRA, supported
    if (eltTypes[zone][0] < 8) nvalidZones++;
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
  FILE* ptrFile = fopen(file, "wb");
  if (ptrFile == NULL)
  {
    printf("Warning: gmshwrite: I can't open file %s.\n", file);
    return 1;
  }
  
  char buf[256];
  double val; int type; int follow; int tagl; int vali;
  E_Int sizeChar = sizeof(char);
  E_Int sizeInt = sizeof(int);
  E_Int sizeFloat = sizeof(double);

  fwrite("$MeshFormat\n", sizeChar, 12, ptrFile);
  fwrite("2.2 1 8\n", sizeChar, 8, ptrFile);
  int one = 1;
  fwrite(&one, sizeInt, 1, ptrFile); 
  fwrite("\n", sizeChar, 1, ptrFile);
  fwrite("$EndMeshFormat\n", sizeChar, 15, ptrFile);

  fwrite("$Nodes\n", sizeChar, 7, ptrFile);
  sprintf(buf, SF_D_ "\n", v.getSize());
  fwrite(buf, sizeChar, strlen(buf), ptrFile);

  /* Nodes */
  for (E_Int i = 0; i < v.getSize(); i++)
  {
    vali = i+1; fwrite(&vali, sizeInt, 1, ptrFile);
    val = v(i,1); fwrite(&val, sizeFloat, 1, ptrFile);
    val = v(i,2); fwrite(&val, sizeFloat, 1, ptrFile);
    val = v(i,3); fwrite(&val, sizeFloat, 1, ptrFile);
  }
  fwrite("\n$EndNodes\n", sizeChar, 11, ptrFile);

  fwrite("$Elements\n", sizeChar, 10, ptrFile);
  size = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    size += connect[i]->getSize();
  }
  sprintf(buf, SF_D_ "\n", size);
  fwrite(buf, sizeChar, strlen(buf), ptrFile);

  // Connectivite par elts par rapport a la definition globale des vertices
  E_Int shift = 0;
  E_Int ind = 1;
  int bnd[2]; bnd[0] = 0; bnd[1] = 14;

  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayI& cn = *connect[i];
    E_Int elt = eltTypes[i][0];
    switch(elt)
    {
      case 0: // NODE
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 15; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 1: // BAR
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 1; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 2: // TRI
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 2; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,3)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 3: // QUAD
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 3; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,3)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,4)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 4: // TETRA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 4; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,3)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,4)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 5: // PYRA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 7; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,3)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,4)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 6: // PENTA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 6; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,3)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,4)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,5)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,6)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
      case 7: // HEXA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          type = 5; fwrite(&type, sizeInt, 1, ptrFile);
          follow = 1; fwrite(&follow, sizeInt, 1, ptrFile);
          tagl = 2; fwrite(&tagl, sizeInt, 1, ptrFile);
          fwrite(&ind, sizeInt, 1, ptrFile);
          fwrite(bnd, sizeInt, 2, ptrFile);
          vali = cn(n,1)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,2)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,3)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,4)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,5)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,6)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,7)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          vali = cn(n,8)+shift; fwrite(&vali, sizeInt, 1, ptrFile);
          ind++;
        }
      }
      break;
    }
    shift += unstructField[i]->getSize();
  }
  fprintf(ptrFile, "\n$EndElements\n");

  fclose(ptrFile);
  return 0;
}
