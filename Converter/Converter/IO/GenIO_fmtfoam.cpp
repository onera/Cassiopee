/*    
    Copyright 2013-2023 Onera.

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

// Formated OpenFoam file support

# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

#if defined(_WIN32)
# include <direct.h>
#endif

using namespace K_FLD;
using namespace std;

// Return true if path already exists
bool dirExist(char* path)
{
#if defined(_WIN32)
  struct _stat info;
  if (_stat(path, &info) != 0) return false;
  return (info.st_mode & _S_IFDIR) != 0;
#else 
  struct stat info;
  if (stat(path, &info) != 0) return false;
  return (info.st_mode & S_IFDIR) != 0;
#endif
}

// Create a directory at given path if it doesnt exist already
// return 0: OK, 1: FAILED
E_Int createDir(char* path)
{
  if (dirExist(path) == true) { printf("exists : %s\n", path); return 0; }
#if defined(_WIN32)
  int ret = _mkdir(path);
#else
  mode_t mode = 0755;
  int ret = mkdir(path, mode);
#endif
  if (ret != 0) return 1;
  return 0;
}

// create a directory structure for foam
// return 0: OK, 1:FAILED
E_Int createSimpleFoamStructure(char* path)
{
  E_Int ret;
  char fullPath[1024];
  // create directory dirName
  ret = createDir(path);
  if (ret == 1) return 1;

  // constant
  strcpy(fullPath, path);
  strcat(fullPath, "/constant");
  ret = createDir(fullPath);
  if (ret == 1) return 1;

  // constant/polymesh
  strcpy(fullPath, path);
  strcat(fullPath, "/constant/polyMesh");
  ret = createDir(fullPath);
  if (ret == 1) return 1;

  // 0
  // 0/epsilon 0/k 0/mut ...
  
  // system
  // system/fvScheme
  // system/fvSolution...
  return 0;
}

//=============================================================================
/* foamread */
// Manque: commentaires
// Nettoyage de f a la fin (vertex reellement utilises)
//=============================================================================
E_Int K_IO::GenIO::foamread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames,
  vector<FldArrayI*>& BCFaces, vector<char*>& BCNames)
{
  printf("\nfoamread\n");

  // Read points
  FldArrayF* f = new FldArrayF();
  foamReadPoints(file, *f);
  //printf("points\n");
  //for (E_Int i = 0; i < f.getSize(); i++) printf("%g %g %g\n", f(i,1), f(i,2), f(i,3));

  // Read NGON
  FldArrayI cNGon; E_Int nfaces;
  foamReadFaces(file, nfaces, cNGon);
  //printf("NGON\n");
  //fflush(stdout);
  //for (E_Int i = 0; i < cNGON.getSize(); i++) printf("%d ", cNGON[i]);
  //printf("\n"); fflush(stdout);

  // Allocate PE
  FldArrayI PE(nfaces, 2);
  foamReadOwner(file, PE);
  //printf("PE 1 (owner)\n");
  //fflush(stdout);
  //for (E_Int i = 0; i < PE.getSize(); i++) printf("%d ", PE(i,1));
  //printf("\n"); fflush(stdout);

  foamReadNeighbour(file, PE);
  //printf("PE 2 (neighbour)\n");
  //fflush(stdout);
  //for (E_Int i = 0; i < PE.getSize(); i++) printf("%d ", PE(i,2));
  //printf("\n"); fflush(stdout);

  // compute NFace
  FldArrayI cNFace; E_Int nelts;
  K_CONNECT::connectFE2NFace(PE, cNFace, nelts);

  // Merge
  E_Int sizeNGon = cNGon.getSize();
  E_Int sizeNFace = cNFace.getSize();
  FldArrayI* cn = new FldArrayI(4+sizeNGon+sizeNFace);
  E_Int* cnp = cn->begin();
  
  cnp[0] = nfaces;
  cnp[1] = cNGon.getSize();
  cnp += 2;
  for (E_Int i = 0; i < sizeNGon; i++) cnp[i] = cNGon[i];
  cnp += sizeNGon;
  cnp[0] = nelts;
  cnp[1] = sizeNFace;
  cnp += 2;
  for (E_Int i = 0; i < sizeNFace; i++) cnp[i] = cNFace[i];
  
  // push in output
  unstructField.push_back(f);
  connect.push_back(cn);
  eltType.push_back(8); // NGon
  char* zoneName = new char [128];
  strcpy(zoneName, "toto");
  zoneNames.push_back(zoneName);
  varString = new char [16];
  strcpy(varString, "x,y,z");

  return 0;
}

//=============================================================================
// Write Mesh point coordinates
//=============================================================================
E_Int K_IO::GenIO::foamWritePoints(char* file, FldArrayF& f)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/points");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\"\n");
  fprintf(ptrFile, "    class       vectorField;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      points;\n");
  fprintf(ptrFile, "}\n");

  E_Int npts = f.getSize();
  E_Float* x = f.begin(1);
  E_Float* y = f.begin(2);
  E_Float* z = f.begin(3);
  
#ifdef E_DOUBLEINT
  fprintf(ptrFile, "%lld\n", npts);
#else
  fprintf(ptrFile, "%d\n", npts);
#endif
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < npts; i++)
  {
    fprintf(ptrFile, "(%.18g %.18g %.18g)\n", x[i], y[i], z[i]);
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

E_Int K_IO::GenIO::foamReadPoints(char* file, FldArrayF& f)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/points");
  FILE* ptrFile = fopen(fullPath, "r");

  E_Int ret;
  readGivenKeyword(ptrFile, "FOAMFILE");
  for (E_Int i = 0; i < 9; i++) skipLine(ptrFile);

  // Passe comments
  char buf[1024]; E_Int l;
  E_Boolean cont = true;
  while (cont)
  {
    readline(ptrFile, buf, 1024); printf("buf=%s\n", buf);
    l = strlen(buf);
    if (l >= 2 && buf[0] == '/' && buf[1] == '/') continue;
    if (l >= 2 && buf[0] == '/' && buf[1] == '*') continue;
    if (l < 2) continue;
    cont = false;
  }

  // Readint in buf
  E_Int npts; E_Int pos=0;
  readInt(buf, 1024, pos, npts);
  printf("npts=%d\n", npts);

  f.malloc(npts, 3);
  E_Float* x = f.begin(1);
  E_Float* y = f.begin(2);
  E_Float* z = f.begin(3);

  skipLine(ptrFile); // (

  for (E_Int i = 0; i < npts; i++)
  {
    //ret = fgetc(ptrFile); // (
    //ret = readDouble(ptrFile, x[i]);
    //ret = readDouble(ptrFile, y[i]);
    //ret = readDouble(ptrFile, z[i]);
    //ret = fgetc(ptrFile); // )

    readline(ptrFile, buf, 1024); pos = 1;
    ret = readDouble(buf, 1024, pos, x[i]);
    ret = readDouble(buf, 1024, pos, y[i]);
    ret = readDouble(buf, 1024, pos, z[i]);
    
  }
  fclose(ptrFile);

  return 0;
}

//=============================================================================
// face indices (NGON)
//=============================================================================
E_Int K_IO::GenIO::foamWriteFaces(char* file, FldArrayI& cn)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/faces");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\"\n");
  fprintf(ptrFile, "    class       faceList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      faces;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = cn.getNFaces();
  E_Int* cnp = cn.getNGon();
  fprintf(ptrFile, "%d\n", nfaces);
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < nfaces; i++)
  {
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "%ld(", cnp[0]);
    for (E_Int i = 1; i <= cnp[0]; i++) fprintf(ptrFile, "%ld ", cnp[i]-1);
    fprintf(ptrFile, ")\n");
#else
    fprintf(ptrFile, "%d(", cnp[0]);
    for (E_Int i = 1; i <= cnp[0]; i++) fprintf(ptrFile, "%d ", cnp[i]-1);
#endif
    fprintf(ptrFile, ")\n");
    cnp += cnp[0]+1;
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

//===========================================================
E_Int K_IO::GenIO::foamReadFaces(char* file, E_Int& nfaces, FldArrayI& cn)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/faces");
  FILE* ptrFile = fopen(fullPath, "r");

  E_Int ret;
  readGivenKeyword(ptrFile, "FOAMFILE");
  for (E_Int i = 0; i < 9; i++) skipLine(ptrFile);

  // Passe comments
  char buf[1024]; E_Int l;
  E_Boolean cont = true;
  while (cont)
  {
    readline(ptrFile, buf, 1024); printf("buf=%s\n", buf);
    l = strlen(buf);
    if (l >= 2 && buf[0] == '/' && buf[1] == '/') continue;
    if (l >= 2 && buf[0] == '/' && buf[1] == '*') continue;
    if (l < 2) continue;
    cont = false;
  }

  // Readint in buf
  E_Int pos=0; 
  readInt(buf, 1024, pos, nfaces);

  skipLine(ptrFile);

  // Find sizeNGon
  E_LONG fpos = KFTELL(ptrFile);
  E_Int sizeNGon = 0; E_Int nf;
  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readline(ptrFile, buf, 1024); pos = 0;
    readInt(buf, 1024, pos, nf);
    sizeNGon += nf+1;
  }
  KFSEEK(ptrFile, fpos, SEEK_SET);
  //printf("sizeNGon=%d\n", sizeNGon); 
  
  cn.malloc(sizeNGon);
  E_Int* cnp = cn.begin();
  
  E_Int val; sizeNGon = 0;
  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readline(ptrFile, buf, 1024); pos = 0;
    readInt(buf, 1024, pos, nf);
    cnp[0] = nf; 
    for (E_Int e = 0; e < nf; e++)
    {
      ret = readInt(buf, 1024, pos, val);
      cnp[1+e] = val+1;
    }
    cnp += nf+1;
    sizeNGon += nf+1;
  }

  printf("sizeNGon=%d, nfaces=%d\n", sizeNGon, nfaces);

  return 0;
}

//=============================================================================
// All faces (left=owner)
//=============================================================================
E_Int K_IO::GenIO::foamWriteOwner(char* file, FldArrayI& PE)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/owner");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\"\n");
  fprintf(ptrFile, "    class       labelList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      owner;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = PE.getSize();
#ifdef E_DOUBLEINT
  fprintf(ptrFile, "%ld\n", nfaces);
#else
  fprintf(ptrFile, "%d\n", nfaces);
#endif
  fprintf(ptrFile, "(\n");

  E_Int c = 0; E_Int val;
  for (E_Int i = 0; i < nfaces; i++)
  {
    val = PE(i, 1);
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "%ld\n", val-1);
#else
    fprintf(ptrFile, "%d\n", val-1);
#endif
    //if (c > 10) { fprintf(ptrFile, "\n"); c = 0; }
    c += 1;
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::foamReadOwner(char* file, FldArrayI& PE)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/owner");
  FILE* ptrFile = fopen(fullPath, "r");

  E_Int ret;
  readGivenKeyword(ptrFile, "FOAMFILE");
  for (E_Int i = 0; i < 9; i++) skipLine(ptrFile);

  // Passe comments
  char buf[1024]; E_Int l;
  /*
  E_Boolean cont = true;
  while (cont)
  {
    readline(ptrFile, buf, 1024);
    l = strlen(buf);
    if (l >= 2 && buf[0] == '/' && buf[1] == '/') continue;
    if (l >= 2 && buf[0] == '/' && buf[1] == '*') continue;
    if (l < 2) continue;
    cont = false;
  }
  */

  // Readint in buf
  E_Int nfaces; E_Int val;
  E_Int pos=0;
  readline(ptrFile, buf, 1024);
  //printf("BUFF=%s\n", buf);
  readInt(buf, 1024, pos, nfaces);
  //printf("NFACES=%d\n", nfaces);

  skipLine(ptrFile);
  
  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, val); PE(i, 1) = val+1;
    //printf("%d\n", val);
  }
  return 0;
}
//=============================================================================
// Internal faces only (right=neighbour)
//=============================================================================
E_Int K_IO::GenIO::foamWriteNeighbour(char* file, FldArrayI& PE)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/neighbour");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\"\n");
  fprintf(ptrFile, "    class       labelList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      neighbour;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = PE.getSize();
  // Count internal faces
  E_Int count = 0; E_Int val1, val2;
  for (E_Int i = 0; i < nfaces; i++)
  {
    val1 = PE(i,1); val2 = PE(i,2);
    if (val1 > 0 && val2 > 0) count += 1; 
  }

#ifdef E_DOUBLEINT
  fprintf(ptrFile, "%ld\n", count);
#else
  fprintf(ptrFile, "%d\n", count);
#endif
  fprintf(ptrFile, "(\n");

  E_Int c = 0;
  for (E_Int i = 0; i < nfaces; i++)
  {
    val1 = PE(i,1); val2 = PE(i,2);
    if (val1 > 0 && val2 > 0)
    {
#ifdef E_DOUBLEINT
      fprintf(ptrFile, "%ld\n", val2-1);
#else
      fprintf(ptrFile, "%d\n", val2-1);
#endif
      //if (c > 10) { fprintf(ptrFile, "\n"); c = 0; }
      c += 1;
    }
  }
  assert(c == count);
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::foamReadNeighbour(char* file, FldArrayI& PE)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/neighbour");
  FILE* ptrFile = fopen(fullPath, "r");

    E_Int ret;
  readGivenKeyword(ptrFile, "FOAMFILE");
  for (E_Int i = 0; i < 9; i++) skipLine(ptrFile);

  // Passe comments
  char buf[1024]; E_Int l;
  /*
  E_Boolean cont = true;
  while (cont)
  {
    readline(ptrFile, buf, 1024);
    l = strlen(buf);
    if (l >= 2 && buf[0] == '/' && buf[1] == '/') continue;
    if (l >= 2 && buf[0] == '/' && buf[1] == '*') continue;
    if (l < 2) continue;
    cont = false;
  }
  */

  // Readint in buf
  E_Int nfaces; E_Int val;
  E_Int pos=0;
  //printf("buf=%s\n", buf);
  readline(ptrFile, buf, 1024);
  readInt(buf, 1024, pos, nfaces);
  //printf("NNEI=%d\n", nfaces);

  skipLine(ptrFile);

  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, val); PE(i, 2) = val+1;
  }

  // tag exterior faces
  for (E_Int i = nfaces; i < PE.getSize(); i++)
  {
    PE(i,2) = 0; // exterior
  }
  return 0;
}

//=============================================================================
// Write to open foam format
//=============================================================================
E_Int K_IO::GenIO::foamwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames,
  PyObject* BCFaces)
{
  createSimpleFoamStructure(file);

  E_Int nzone = unstructField.size();

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: foamwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;


  // constant/polyMesh/points
  // constant/polyMesh/faces
  // constant/polyMesh/neighbour
  // constant/polyMesh/owner
  // constant/polyMesh/boundary
  // constant/polyMesh/faceZones
  // constant/polyMesh/cellZones
    
  // limited to one NGON zone
  if (nzone == 0) 
  {
    printf("Warning: foamwrite: no unstructured zone in input. Nothing written.\n");
    return 1;
  }

  // find first NGON zone
  E_Int no = -1;
  for (E_Int i = 0; i < nzone; i++)
  {
    E_Int et = eltType[i];
    if (et == 8) { no = i; break; }
  }
  if (no == -1)
  {
    printf("Warning: foamwrite: no NGON zone found in input. Nothing written.\n");
    return 1;
  }

  FldArrayF& field = *unstructField[no];
  FldArrayI& cn = *connect[no];

  // Compute PE
  E_Int nfaces = cn.getNFaces();
  E_Int nelts = cn.getNElts();
  FldArrayI cFE(nfaces,2);
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);
  E_Int* ptrNF = cn.getNFace();
  E_Int face, nf;
  for (E_Int i = 0; i < nfaces*2; i++) cFE[i] = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    nf = ptrNF[0]; // nb de face pour l'elt
    for (E_Int j = 1; j <= nf; j++)
    {
      face = ptrNF[j]-1;
	  if (facesp1[face] == 0) facesp1[face] = i+1;
      else facesp2[face] = i+1;
    }
    ptrNF += nf+1;
  }

  foamWritePoints(file, field);
  foamWriteFaces(file, cn);
  foamWriteOwner(file, cFE);
  foamWriteNeighbour(file, cFE);

  return 0;
}
