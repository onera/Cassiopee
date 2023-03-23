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
# include <queue>

# include <dirent.h>

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


// Return true if path is a file
E_Int fileExist(char* path)
{
#if defined(_WIN32)
  struct _stat info;
  if (_stat(path, &info) != 0) return false;
  return (info.st_mode & _S_IFREG) != 0;
#else 
  struct stat info;
  if (stat(path, &info) != 0) return false;
  return (info.st_mode & S_IFREG) != 0;
#endif
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

E_Int K_IO::GenIO::readScalarField(char *file, FldArrayF& f, E_Int idx)
{
  FILE* ptrFile = fopen(file, "r");
  assert(ptrFile);

  readGivenKeyword(ptrFile, "<SCALAR>");

  skipLine(ptrFile);

  char buf[1024];
  readline(ptrFile, buf, 1024);
  
  E_Int ncells; E_Int pos=0;
  readInt(buf, 1024, pos, ncells);

  E_Float* fld = f.begin(idx);

  skipLine(ptrFile);

  for (E_Int i = 0; i < ncells; i++)
  {
    readline(ptrFile, buf, 1024); pos = 0;
    readDouble(buf, 1024, pos, fld[i]);
  }
  fclose(ptrFile);

  return ncells;
}

E_Int K_IO::GenIO::readVectorField(char *file, FldArrayF& f, E_Int idx)
{
  FILE* ptrFile = fopen(file, "r");
  assert(ptrFile);

  readGivenKeyword(ptrFile, "<VECTOR>");
  
  skipLine(ptrFile);

  char buf[1024];
  readline(ptrFile, buf, 1024);
  
  E_Int ncells; E_Int pos=0;
  readInt(buf, 1024, pos, ncells);

  E_Float* fldx = f.begin(idx);
  E_Float* fldy = f.begin(idx+1);
  E_Float* fldz = f.begin(idx+2);

  skipLine(ptrFile);

  for (E_Int i = 0; i < ncells; i++)
  {
    readline(ptrFile, buf, 1024); pos = 1;
    readDouble(buf, 1024, pos, fldx[i]);
    readDouble(buf, 1024, pos, fldy[i]);
    readDouble(buf, 1024, pos, fldz[i]);
  }
  fclose(ptrFile);

  return ncells;
}

E_Int K_IO::GenIO::readTensorField(char *file, FldArrayF& f, E_Int idx)
{
  FILE* ptrFile = fopen(file, "r");
  assert(ptrFile);

  readGivenKeyword(ptrFile, "<TENSOR>");
  
  skipLine(ptrFile);

  char buf[1024];
  readline(ptrFile, buf, 1024);
  
  E_Int ncells; E_Int pos=0;
  readInt(buf, 1024, pos, ncells);

  E_Float* fldxx = f.begin(idx);
  E_Float* fldxy = f.begin(idx+1);
  E_Float* fldxz = f.begin(idx+2);
  E_Float* fldyx = f.begin(idx+3);
  E_Float* fldyy = f.begin(idx+4);
  E_Float* fldyz = f.begin(idx+5);
  E_Float* fldzx = f.begin(idx+6);
  E_Float* fldzy = f.begin(idx+7);
  E_Float* fldzz = f.begin(idx+8);

  skipLine(ptrFile);

  for (E_Int i = 0; i < ncells; i++)
  {
    readline(ptrFile, buf, 1024); pos = 1;
    readDouble(buf, 1024, pos, fldxx[i]);
    readDouble(buf, 1024, pos, fldxy[i]);
    readDouble(buf, 1024, pos, fldxz[i]);
    readDouble(buf, 1024, pos, fldyx[i]);
    readDouble(buf, 1024, pos, fldyy[i]);
    readDouble(buf, 1024, pos, fldyz[i]);
    readDouble(buf, 1024, pos, fldzx[i]);
    readDouble(buf, 1024, pos, fldzy[i]);
    readDouble(buf, 1024, pos, fldzz[i]);
  }
  fclose(ptrFile);

  return ncells;
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
  vector<FldArrayI*>& BCFaces, vector<char*>& BCNames,
  char*& varStringc,
  vector<FldArrayF*>& centerStructField,
  vector<FldArrayF*>& centerUnstructField)
{
  puts("\n");

  // Read points
  FldArrayF* f = new FldArrayF();
  foamReadPoints(file, *f);

  // Read NGON
  FldArrayI cNGon; E_Int nfaces;
  foamReadFaces(file, nfaces, cNGon);

  // Allocate PE
  FldArrayI PE(nfaces, 2);
  foamReadOwner(file, PE);
  foamReadNeighbour(file, PE);

  // compute NFace
  FldArrayI cNFace; E_Int nelts;
  K_CONNECT::connectFE2NFace(PE, cNFace, nelts);
#ifdef E_DOUBLEINT
  printf("cells: %ld\n", nelts);
#else
  printf("cells: %d\n", nelts);
#endif

  // Read fields
  foamReadFields(file, centerUnstructField, nelts, varStringc);

  // Read boundary
  foamReadBoundary(file, BCFaces, BCNames);

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
  strcpy(zoneName, "FOAM_ZONE");
  zoneNames.push_back(zoneName);
      
  varString = new char [16];
  // add field names
  strcpy(varString, "x,y,z");

  puts("");

  return 0;
}



#define MAX_FIELDS 20

E_Int K_IO::GenIO::foamReadFields(char *file, std::vector<FldArrayF*>& centerUnstructField,
  E_Int ncells, char*& varStringc)
{
  // identify latest time output folder
  DIR *d;
  struct dirent *dir;
  d = opendir(file);

  char fullPath[1024];
  char path[1024];
  strcpy(fullPath, file);
  E_Float ret;
  E_Float max_time = 0;
  while ((dir = readdir(d))) {
    strcpy(path, fullPath);
    strcat(path, "/");
    strcat(path, dir->d_name);
    if (dirExist(path)) {
      ret = atof(dir->d_name);
      max_time = std::max(max_time, ret);
    }
  }
  closedir(d);

  // loop again and get fullPath
  bool found = false;
  char dir_name[260] = {0};
  dir_name[0] = '\0';
  d = opendir(file);
  while ((dir = readdir(d))) {
    strcpy(path, fullPath);
    strcat(path, "/");
    strcat(path, dir->d_name);
    if (dirExist(path)) {
      ret = atof(dir->d_name);
      if (ret == max_time) {
        sprintf(dir_name, "%s", dir->d_name);
        found = true;
      }
      if (found) break;
    }
    if (found) break;
  }
  closedir(d);


  fullPath[0] = '\0';
  strcat(fullPath, file);
  strcat(fullPath, "/");
  strcat(fullPath, dir_name);
  
  d = opendir(fullPath);
  assert(d);
  varStringc = new char[1024];
  varStringc[0] = '\0';

  E_Int nflds = 0;
  E_Int size = 0;
  E_Int field_type[MAX_FIELDS] = {1}; // 1: scalar, 2: vector, 3: tensor
  char field_name[MAX_FIELDS][128];

  while ((dir = readdir(d))) {
    strcpy(path, fullPath);
    strcat(path, "/");
    strcat(path, dir->d_name);
    if (fileExist(path)) {
      // skip .swp files
      if (dir->d_name[0] == '.') continue;
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, dir->d_name);
      
      // read only volScalarFields and volVectorFields
      FILE *fh = fopen(path, "r");
      assert(fh);
      E_Int ret = readGivenKeyword(fh, "VOLSCALARFIELD", "VOLVECTORFIELD", "VOLTENSORFIELD");
      if (ret == 1) { // volScalarField
        strcpy(field_name[size], dir->d_name);

        field_type[size] = 1;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, ",");

        nflds++;
        size++;
      } else if (ret == 2) { // volVectorField
        strcpy(field_name[size], dir->d_name);

        field_type[size] = 2;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "x,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "y,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "z,");
        nflds++;

        size++;
      } else if (ret == 3) { // volTensorField
        strcpy(field_name[size], dir->d_name);

        field_type[size] = 3;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "xx,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "xy,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "xz,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "yx,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "yy,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "yz,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "zx,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "zy,");
        nflds++;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, "zz,");
        nflds++;

        size++;
      }
      if (size == MAX_FIELDS) {
        fprintf(stderr, "Warning: foamread: Trying to read more that maximum number of fields (%d). Aborting.\n", MAX_FIELDS);
        exit(1);
      }
      fclose(fh);
    }
  }
  closedir(d);

  // delete the last comma
  varStringc[strlen(varStringc)-1] = '\0';

  FldArrayF *F = new FldArrayF();
  F->malloc(ncells, nflds);

  d = opendir(fullPath);
  assert(d);

  E_Int idx = 1;
  for (E_Int fld = 0; fld < size; fld++) {
    if (field_type[fld] == 1) {
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, field_name[fld]);
      E_Int ret = readScalarField(path, *F, idx); 
      assert(ret == ncells);
      idx++;
      printf("Info: foamread: reading scalar field %s\n", field_name[fld]);
    } else if (field_type[fld] == 2) {
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, field_name[fld]);
      E_Int ret = readVectorField(path, *F, idx); 
      assert(ret == ncells);
      idx += 3;
      printf("Info: foamread: reading vector field %s\n", field_name[fld]);
    } else if (field_type[fld] == 3) {
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, field_name[fld]);
      E_Int ret = readTensorField(path, *F, idx); 
      assert(ret == ncells);
      idx += 9; 
      printf("Info: foamread: reading tensor field %s\n", field_name[fld]);
    } else {
      assert(false);
    }
  }

  printf("Info: foamread: done reading fields.\n");

  centerUnstructField.push_back(F);

  return 0;
}

E_Int K_IO::GenIO::foamReadPoints(char* file, FldArrayF& f)
{
  char fullPath[1024];
  fullPath[0] = '\0';
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/points");
  FILE* ptrFile = fopen(fullPath, "r");

  E_Int ret;
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;
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

  // Readint in buf
  E_Int npts; E_Int pos=0;
  readInt(buf, 1024, pos, npts);

  f.malloc(npts, 3);

  E_Float* x = f.begin(1);
  E_Float* y = f.begin(2);
  E_Float* z = f.begin(3);

  skipLine(ptrFile); // (

  for (E_Int i = 0; i < npts; i++)
  {
    readline(ptrFile, buf, 1024); pos = 1;
    readDouble(buf, 1024, pos, x[i]);
    readDouble(buf, 1024, pos, y[i]);
    readDouble(buf, 1024, pos, z[i]);
  }
  fclose(ptrFile);

#ifdef E_DOUBLEINT
  printf("points: %ld\n", f.getSize());
#else
  printf("points: %d\n", f.getSize());
#endif

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
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;
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
  fclose(ptrFile);

#ifdef E_DOUBLEINT
  printf("faces: %ld\n", nfaces);
#else
  printf("faces: %d\n", nfaces);
#endif

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
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;

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

  // Readint in buf
  E_Int nfaces; E_Int val;
  E_Int pos=0;
  readInt(buf, 1024, pos, nfaces);

  skipLine(ptrFile);
  
  for (E_Int i = 0; i < nfaces; i++)
  {
    ret = readInt(ptrFile, val); PE(i, 1) = val+1;
  }

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
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;

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

  // Readint in buf
  E_Int nfaces; E_Int val;
  E_Int pos=0;
  readInt(buf, 1024, pos, nfaces);

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

#ifdef E_DOUBLEINT
  printf("internal faces: %ld\n", nfaces);
#else
  printf("internal faces: %d\n", nfaces);
#endif

  fclose(ptrFile);

  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::foamReadBoundary(char* file, std:: vector<FldArrayI*>& BCFaces,
  std::vector<char*>& BCNames)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/boundary");
  FILE* ptrFile = fopen(fullPath, "r");

  E_Int ret;
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;

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

  // Readint in buf
  E_Int nBC;
  E_Int pos=0;
  readInt(buf, 1024, pos, nBC);

  skipLine(ptrFile);

  // extract total number and names of boundary faces
#define BCSTRINGMAXSIZE 50

  char **bcnames = (char **) malloc(nBC * sizeof(char *));
  char **type = (char **) malloc(nBC * sizeof(char *));
  for (E_Int i = 0; i < nBC; i++) {
    bcnames[i] = (char *) malloc(BCSTRINGMAXSIZE);
    type[i] = (char *) malloc(BCSTRINGMAXSIZE);
  }

  E_Int *nFaces = (E_Int *) malloc(nBC * sizeof(E_Int));
  E_Int *startFace = (E_Int *) malloc(nBC * sizeof(E_Int));

  for (E_Int i = 0; i < nBC; i++)
  {
    // name
    readWord(ptrFile, bcnames[i]);

    // type
    readGivenKeyword(ptrFile, "TYPE");
    readWord(ptrFile, type[i]);
    type[i][strlen(type[i])-1] = '\0';
    strcat(bcnames[i], "@");
    strcat(bcnames[i], type[i]);

    // nFaces
    readGivenKeyword(ptrFile, "NFACES");
    readWord(ptrFile, buf);
    buf[strlen(buf)-1] = '\0';
    nFaces[i] = convertString2Int(buf);

    // startFace
    readGivenKeyword(ptrFile, "STARTFACE");
    readWord(ptrFile, buf);
    buf[strlen(buf)-1] = '\0';
    startFace[i] = convertString2Int(buf);

    skipLine(ptrFile);
  }

  fclose(ptrFile);

  E_Int nboundaryfaces = 0;
  for (E_Int i = 0; i < nBC; i++) nboundaryfaces += nFaces[i];
  FldArrayI *faces = new FldArrayI(nboundaryfaces);
  char *names = new char [nboundaryfaces * BCSTRINGMAXSIZE];
  E_Int *facesp = faces->begin();
  E_Int le;
  E_Int c = 0;

  E_Int k = 0;
  for (E_Int i = 0; i < nBC; i++) {
    printf("Reading %s\n", bcnames[i]);
    le = strlen(bcnames[i]);
    le = K_FUNC::E_min(le, BCSTRINGMAXSIZE-1);
    for (E_Int j = 1; j <= nFaces[i]; j++) {
      facesp[k++] = startFace[i] + j;
      for (E_Int l = 0; l < le; l++) { names[c+l] = bcnames[i][l]; }
      c += le;
      names[c] = '\0'; c++;
    }
  }
  assert(k == nboundaryfaces);

  BCFaces.push_back(faces);
  BCNames.push_back(names);

  for (E_Int i = 0; i < nBC; i++) {
    free(bcnames[i]);
    free(type[i]);
  }
  free(bcnames);
  free(type);
  free(startFace);
  free(nFaces);

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
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
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

//=============================================================================
// face indices (NGON)
//=============================================================================
E_Int K_IO::GenIO::foamWriteFaces(char* file, const ngon_t<K_FLD::IntArray>& NG,
  const std::vector<E_Int>& faces)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/faces");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       faceList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      faces;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = NG.PGs.size();

#ifdef E_DOUBLEINT
  fprintf(ptrFile, "%ld\n", nfaces);
#else
  fprintf(ptrFile, "%d\n", nfaces);
#endif
  fprintf(ptrFile, "(\n");

  E_Int idx, stride;
  for (E_Int i = 0; i < nfaces; i++)
  {
    idx = faces[i];
    const E_Int *pN = NG.PGs.get_facets_ptr(idx);
    stride = NG.PGs.stride(idx);
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "%ld(", stride);
    for (E_Int k = 0; k < stride; k++) fprintf(ptrFile, "%ld ", pN[k]-1);
#else
    fprintf(ptrFile, "%d(", stride);
    for (E_Int k = 0; k < stride; k++) fprintf(ptrFile, "%d ", pN[k]-1);
#endif
    fprintf(ptrFile, ")\n");
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);

  return 0;
}



//=============================================================================
// All faces (left=owner)
//=============================================================================
E_Int K_IO::GenIO::foamWriteOwner(char* file, const K_FLD::IntArray& F2E, const std::vector<E_Int>& faces)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/owner");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       labelList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      owner;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = F2E.cols();
#ifdef E_DOUBLEINT
  fprintf(ptrFile, "%ld\n", nfaces);
#else
  fprintf(ptrFile, "%d\n", nfaces);
#endif
  fprintf(ptrFile, "(\n");

  E_Int own, idx;
  for (E_Int i = 0; i < nfaces; i++)
  {
    idx = faces[i];
    own = F2E(0, idx);
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "%ld\n", own);
#else
    fprintf(ptrFile, "%d\n", own);
#endif
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}


//=============================================================================
// Internal faces only (right=neighbour)
//=============================================================================
E_Int K_IO::GenIO::foamWriteNeighbour(char* file, const K_FLD::IntArray& F2E, const std::vector<E_Int>& faces,
  const E_Int ninternal_faces)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/neighbour");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       labelList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      neighbour;\n");
  fprintf(ptrFile, "}\n");

  //E_Int nfaces = F2E.cols();

#ifdef E_DOUBLEINT
  fprintf(ptrFile, "%ld\n", ninternal_faces);
#else
  fprintf(ptrFile, "%d\n", ninternal_faces);
#endif
  fprintf(ptrFile, "(\n");

  E_Int idx, nei;
  for (E_Int i = 0; i < ninternal_faces; i++)
  {
    idx = faces[i];
    nei = F2E(1, idx);
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "%ld\n", nei);
#else
    fprintf(ptrFile, "%d\n", nei);
#endif
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::foamWriteBoundary(char* file, const std::vector<char*>& bc_names,
  const std::vector<E_Int>& bc_nfaces, const std::vector<E_Int>& bc_startfaces)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/boundary");
  FILE* ptrFile = fopen(fullPath, "w");

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       polyBoundaryMesh;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      boundary;\n");
  fprintf(ptrFile, "}\n");

  fprintf(ptrFile, "%ld\n", bc_names.size());
  fprintf(ptrFile, "(\n");

  for (size_t i = 0; i < bc_names.size(); i++)
  {
    char *token = strtok(bc_names[i], "@");
    char name[strlen(token)+1];
    strcpy(name, token);

    fprintf(ptrFile, "    %s\n", token);
    fprintf(ptrFile, "    {\n");

    token = strtok(NULL, "@");
    if (token == NULL) {
      fprintf(stderr, "No type for BC %s, defaulting to wall.\n", bc_names[i]);
      fprintf(ptrFile, "        type            %s;\n", "wall");
      fprintf(ptrFile, "        physicalType    %s;\n", "wall");
    } else {
      fprintf(ptrFile, "        type            %s;\n", token);
      fprintf(ptrFile, "        physicalType    %s;\n", token);
    }
    
#ifdef E_DOUBLEINT
    fprintf(ptrFile, "        nFaces          %ld;\n", bc_nfaces[i]);
    fprintf(ptrFile, "        startFace       %ld;\n", bc_startfaces[i]);
#else
    fprintf(ptrFile, "        nFaces          %d;\n", bc_nfaces[i]);
    fprintf(ptrFile, "        startFace       %d;\n", bc_startfaces[i]);
#endif
    fprintf(ptrFile, "    }\n");
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
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

  std::cout << "nzone = " << nzone << std::endl;

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

  foamWritePoints(file, field);

  K_FLD::IntArray CN(cn);
  K_FLD::FloatArray CRD(field); 

  ngon_t<K_FLD::IntArray> NG(CN);
  NG.flag_externals(1);

  DELAUNAY::Triangulator dt;
  bool has_been_reversed;
  ngon_t<K_FLD::IntArray>::reorient_skins(dt, CRD, NG, has_been_reversed);

  // F2E
  ngon_unit neighbors;
  K_FLD::IntArray F2E;
  NG.build_ph_neighborhood(neighbors);
  NG.build_F2E(neighbors, F2E);

  std::vector<E_Int> faces;
  std::vector<uint8_t> marked(NG.PGs.size(), 0);
  //E_Int nfaces = NG.PGs.size();

  E_Int PGi, stride;
  std::vector<E_Int> neis;
  std::vector<E_Int> pgs;

  for (E_Int PHi = 0; PHi < NG.PHs.size(); PHi++) {
    const E_Int *pF = NG.PHs.get_facets_ptr(PHi);
    stride = NG.PHs.stride(PHi);

    neis.clear();
    pgs.clear();
    
    for (E_Int j = 0; j < stride; j++) {
      PGi = pF[j] - 1;
      
      if (F2E(0, PGi) == IDX_NONE || F2E(1, PGi) == IDX_NONE) continue;

      if (!marked[PGi]) {
        neis.push_back(NEIGHBOR(PHi, F2E, PGi));
        pgs.push_back(PGi);
        marked[PGi] = 1;
        if (F2E(0,PGi) > F2E(1,PGi)) {
          std::swap(F2E(0,PGi), F2E(1,PGi));
          E_Int *pN = NG.PGs.get_facets_ptr(PGi);
          std::reverse(pN,pN+4);
        }
      }
    }

    // Note (Imad) : internal faces are sorted in increasing order of corresponding neighbours (upper triangular ordering)
    std::vector<E_Int> order(neis.size());
    std::iota(order.begin(), order.end(), 0); // init
    std::sort(order.begin(), order.end(), [&](E_Int i, E_Int j){return neis[i] < neis[j];});

    for (size_t i = 0; i < order.size(); i++) {
      E_Int index = order[i];
      faces.push_back(pgs[index]);
    }
  }

  E_Int ninternal_faces = faces.size();

  std::cout << "internal faces: " << ninternal_faces << std::endl;

  // BC
  E_Int BCFacesSize = 0;
  if (PyList_Check(BCFaces)) BCFacesSize = PyList_Size(BCFaces);

  std::vector<E_Int> start_face_per_bc;
  std::vector<E_Int> nfaces_per_bc;
  std::vector<char*> name_per_bc;

  if (BCFacesSize > 0) {
    E_Int indFace;
    IMPORTNUMPY;

    PyObject* BCs = PyList_GetItem(BCFaces, 0);
    E_Int size = PyList_Size(BCs);

    if (size == 0) {
      printf("Warning: foamwrite: requires boundary patches.\n");
      //exit(1);
    }

    E_Int np;
    E_Int *ptr;
    char *name;
    for (E_Int j = 0; j < size/2; j++) {
      name = NULL;
      PyObject *o = PyList_GetItem(BCs, 2*j);
      
      if (PyString_Check(o)) name = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(o)) name = (char *)PyUnicode_AsUTF8(o);
#endif
      name_per_bc.push_back(0);
      name_per_bc[j] = new char[strlen(name) + 1];
      //strncpy(name_per_bc[j], name, strlen(name));
      strcpy(name_per_bc[j], name);
      PyArrayObject *array = (PyArrayObject *) PyList_GetItem(BCs, 2*j+1);
      ptr = (E_Int *) PyArray_DATA(array);
      np = PyArray_SIZE(array); // number of faces in current boundary
      nfaces_per_bc.push_back(np);
      start_face_per_bc.push_back(faces.size());
      for (E_Int k = 0; k < np; k++) {
        indFace = ptr[k]-1;
        faces.push_back(indFace);
      }
    }
  }

  foamWriteFaces(file, NG, faces);
  foamWriteOwner(file, F2E, faces);
  foamWriteNeighbour(file, F2E, faces, ninternal_faces);
  foamWriteBoundary(file, name_per_bc, nfaces_per_bc, start_face_per_bc);

  for (size_t i = 0; i < name_per_bc.size(); i++)
    delete [] name_per_bc[i];

  return 0;
}
