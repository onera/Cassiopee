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
// return 0: ok, 1: FAILED
E_Int createDir(char* path)
{
  if (dirExist(path) == true) return 0;
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
  printf("foamread\n");

  // Read points
  FldArrayF f;
  foamReadPoints(file, f);
  for (E_Int i = 0; i < f.getSize(); i++) printf("%g\n", f(i,1));

  // Read NGON
  FldArrayI cNGON; E_Int nfaces;
  foamReadFaces(file, nfaces, cNGON);

  // Allocate PE
  FldArrayI PE(nfaces, 2);
  foamReadOwner(file, PE);
  foamReadNeighbour(file, PE);

  // Merge in a single connect
  FldArrayI cNFace; E_Int nelts;
  K_CONNECT::connectFE2NFace(PE, cNFace, nelts);

  return 0;
}

//=============================================================================
// Mesh point coordinates
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
    fprintf(ptrFile, "(17.18%g 17.18%g 17.18%g)\n", x[i], y[i], z[i]);
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

  E_Int npts;
  ret = -1;
  while (ret == -1) ret = readInt(ptrFile, npts);

  char buf[1024];
  readline(ptrFile, buf, 1024); printf("buf=%s\n", buf);
  printf("npts=%d\n", npts);

  f.malloc(npts, 3);
  E_Float* x = f.begin(1);
  E_Float* y = f.begin(2);
  E_Float* z = f.begin(3);

  skipLine(ptrFile);

  for (E_Int i = 0; i < npts; i++)
  {
    ret = fgetc(ptrFile); // (
    ret = readDouble(ptrFile, x[i]);
    ret = readDouble(ptrFile, y[i]);
    ret = readDouble(ptrFile, z[i]);
    ret = fgetc(ptrFile); // )
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
  fprintf(ptrFile, "    class       vectorField;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      faces;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = cn.getNFaces();
  E_Int* cnp = cn.getNGon();
  fprintf(ptrFile, "%d\n", nfaces);
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < nfaces; i++)
  {
    fprintf(ptrFile, "%d(", cnp[0]);
    for (E_Int i = 1; i <= cnp[0]; i++) fprintf(ptrFile, "%d ", cnp[i]);
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
  FILE* ptrFile = fopen(fullPath, "w");

  E_Int ret; char buf[1024];
  for (E_Int i = 0; i < 8; i++) ret = readline(ptrFile, buf, 1024);

  ret = readInt(ptrFile, nfaces);
  ret = readline(ptrFile, buf, 1024);

  // Compte la taille du tableau
  E_LONG pos = KFTELL(ptrFile);
  E_Int sizeNGon = 0; E_Int val;
  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, val);
    sizeNGon += val+1;
    skipLine(ptrFile);
  }

  cn.malloc(sizeNGon);
  E_Int* cnp = cn.begin();
  KFSEEK(ptrFile, pos, SEEK_SET);

  E_Int nf;
  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, nf);
    cnp[0] = val;
    fgetc(ptrFile);
    for (E_Int e = 0; e < nf; e++)
    {
      ret = readInt(ptrFile, val);
      cnp[1+e] = val;
    }
    cnp += nf+1;
  }

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
  fprintf(ptrFile, "    class       vectorField;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      owner;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = PE.getSize();
  fprintf(ptrFile, "%d\n", nfaces);
  fprintf(ptrFile, "(\n");

  E_Int c = 0; E_Int val;
  for (E_Int i = 0; i < nfaces; i++)
  {
    val = PE(i,1);
    if (val == 0) val = -1;
    fprintf(ptrFile, "%d ", val);
    if (c > 10)
    {
      fprintf(ptrFile, "\n"); c = 0;
    }
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
  FILE* ptrFile = fopen(fullPath, "w");

  E_Int ret; char buf[1024];
  for (E_Int i = 0; i < 8; i++) ret = readline(ptrFile, buf, 1024);

  E_Int nfaces; E_Int val;
  ret = readInt(ptrFile, nfaces);
  ret = readline(ptrFile, buf, 1024);

  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, val);
    PE(i, 1) = val;
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
  fprintf(ptrFile, "    class       vectorField;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      neighbour;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = PE.getSize();
  fprintf(ptrFile, "%d\n", nfaces);
  fprintf(ptrFile, "(\n");

  // Count internal faces
  E_Int count = 0; E_Int val1, val2;
  for (E_Int i = 0; i < nfaces; i++)
  {
    val1 = PE(i,1); val2 = PE(i,2);
    if (val1 > 0 && val2 > 0) count += 1; 
  }

  E_Int c = 0;
  for (E_Int i = 0; i < count; i++)
  {
    val1 = PE(i,1); val2 = PE(i,2);
    if (val1 > 0 && val2 > 0)
    {
      fprintf(ptrFile, "%d ", val2);
      if (c > 10)
      {
        fprintf(ptrFile, "\n"); c = 0;
      }
      c += 1;
    }
  }
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
  FILE* ptrFile = fopen(fullPath, "w");

  E_Int ret; char buf[1024];
  for (E_Int i = 0; i < 8; i++) ret = readline(ptrFile, buf, 1024);

  E_Int nfaces; E_Int val;
  ret = readInt(ptrFile, nfaces);
  ret = readline(ptrFile, buf, 1024);

  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, val);
    PE(i, 2) = val;
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
