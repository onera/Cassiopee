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

// Formated OpenFoam file support

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "GenIO.h"
#include "Array/Array.h"
# include "String/kstring.h"
#include <vector>
#include "Def/DefFunction.h"
#include "Connect/connect.h"
#include <queue>
#include "Metric/metric.h"
#include "Math/math.h"

#include <dirent.h>

#if defined(_WIN32)
#include <direct.h>
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
  if (dirExist(path) == true) { printf("exists: %s\n", path); return 0; }
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

E_Int K_IO::GenIO::readScalarField(char *file, FldArrayF& f, E_Int idx,
  E_Int *owner, const std::vector<FldArrayI *> &BCFaces,
  std::vector<FldArrayF *> &BCFields, const std::vector<E_Float> &delta,
  E_Int nifaces, const std::vector<E_Int> &indir)
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

  // Boundary fields
  E_Int ret = readGivenKeyword(ptrFile, "BOUNDARYFIELD");
  if (ret != 1) printf("fmt_foam: cant find BOUNDARYFIELD.\n");
  assert(ret == 1);
  skipLine(ptrFile);
  
  E_Int nbnd = indir.size()-1;

  for (E_Int i = 0; i < nbnd; i++) 
  {
    //const char *bcname = BCNames[i];
    E_Int bcsize = indir[i+1] - indir[i];
    E_Int *bcfaces = BCFaces[0]->begin() + indir[i];
    E_Float *bcf = BCFields[i]->begin(idx);

    // name
    char bcname[256];
    readWord(ptrFile, bcname);

    // type
    char type[256];
    ret = readGivenKeyword(ptrFile, "TYPE");
    if (ret != 1) printf("INFO: foamread: cant find TYPE.\n");
    assert(ret == 1);
    readWord(ptrFile, type);

    if (strcmp(type, "fixedGradient") == 0) {
      puts(type);
      E_Int ret = readGivenKeyword(ptrFile, "gradient");
      ret &= readGivenKeyword(ptrFile, "uniform");
      //assert(ret == 1);

      char buf[256];
      readWord(ptrFile, buf);
      char *endptr;

      E_Float gradVal = strtod(buf, &endptr);

      // Make sure a conversion took place
      assert(gradVal != 0.0 && buf != endptr);

      // valb = valb + gradient * delta
      for (E_Int j = 0; j < bcsize; j++) {
        E_Int face = bcfaces[j]-1;
        E_Int own = owner[face]-1;
        bcf[j] = fld[own] + delta[face-nifaces] * gradVal;
      }
    }

    else if (strcmp(type, "zeroGradient") == 0) {
      puts(type);
      // Copy cell values
      for (E_Int j = 0; j < bcsize; j++) {
        E_Int face = bcfaces[j]-1;
        E_Int own = owner[face]-1;
        bcf[j] = fld[own];
      }
    }

    else if (strcmp(type, "empty") == 0) {
      puts(type);
      // Copy cell values
      for (E_Int j = 0; j < bcsize; j++) {
        E_Int face = bcfaces[j]-1;
        E_Int own = owner[face]-1;
        bcf[j] = fld[own];
      }
    }

    else if (strcmp(type, "inletOutlet")       == 0 ||
             strcmp(type, "calculated")        == 0 ||
             strcmp(type, "totalPressure")     == 0 ||
             strcmp(type, "fixedFluxPressure") == 0 ||
             strcmp(type, "waveTransmissive")  == 0 ||
             strcmp(type, "freestreamPressure") == 0) {

      puts(type);

      // Value can be uniform <val> or nonuniform List<scalar>
      readGivenKeyword(ptrFile, "VALUE");
      char buf[128];
      readWord(ptrFile, buf);

      if (strcmp(buf, "uniform") == 0) {
        readWord(ptrFile, buf);
        E_Float val = strtod(buf, NULL);
        printf("Info: foamread: uniform " SF_F_ "\n", val);
        for (E_Int j = 0; j < bcsize; j++) {
          bcf[i] = val;
        }
      } else {
        // Read list after nonuniform<scalar>
        readGivenKeyword(ptrFile, "<SCALAR>");
        skipLine(ptrFile);

        // Make sure you read a number equal to bcsize
        char buf[256];
        readline(ptrFile, buf, 256);

        E_Int nfaces; E_Int pos=0;
        readInt(buf, 1024, pos, nfaces);
        assert(nfaces == bcsize);

        printf("Info: foamread: non uniform list of size " SF_D_ ".\n", nfaces);

        skipLine(ptrFile); // (

        for (E_Int j = 0; j < nfaces; j++)
        {
          readline(ptrFile, buf, 1024); pos = 0;
          readDouble(buf, 1024, pos, bcf[j]);
        }
      }
    }

    else if (strcmp(type, "symmetryPlane") == 0) {
      puts(type);
      // Copy cell values
      for (E_Int j = 0; j < bcsize; j++) {
        E_Int face = bcfaces[j]-1;
        E_Int own = owner[face]-1;
        bcf[j] = fld[own];
      }
    }

    else if (strcmp(type, "fixedValue")        == 0 ||
             strcmp(type, "uniformFixedValue") == 0) {
      puts(type);
      readGivenKeyword(ptrFile, "UNIFORM");
      char buf[128];
      readWord(ptrFile, buf);
      E_Float val = strtod(buf, NULL);
      printf("Info: foamread: uniform " SF_F_ ".\n", val);
      for (E_Int j = 0; j < bcsize; j++) {
        bcf[j] = val;
      }
    }

    else if (strcmp(type, "inletOutlet") == 0 ||
             strcmp(type, "free")) {
      puts(type);

      // Value can be uniform <val> or nonuniform List<scalar>
      readGivenKeyword(ptrFile, "VALUE");
      char buf[128];
      readWord(ptrFile, buf);

      if (strcmp(buf, "uniform") == 0) {
        readWord(ptrFile, buf);
        E_Float val = strtod(buf, NULL);
        printf("Info: foamread: uniform " SF_F_ ".\n", val);
        for (E_Int j = 0; j < bcsize; j++) {
          bcf[i] = val;
        }
      } else {
        // Read list after nonuniform<scalar>
        readGivenKeyword(ptrFile, "<SCALAR>");
        skipLine(ptrFile);

        // Make sure you read a number equal to bcsize
        char buf[256];
        readline(ptrFile, buf, 256);

        E_Int nfaces; E_Int pos=0;
        readInt(buf, 1024, pos, nfaces);
        assert(nfaces == bcsize);

        printf("Info: foamread: non uniform list of size " SF_D_ ".\n", nfaces);

        skipLine(ptrFile); // (

        for (E_Int j = 0; j < nfaces; j++)
        {
          readline(ptrFile, buf, 1024); pos = 0;
          readDouble(buf, 1024, pos, bcf[j]);
        }
      }
    }

    else if (strcmp(type, "externalWallHeatFluxTemperature") == 0 ||
             strcmp(type, "humidityTemperatureCoupledMixed") == 0 ||
             strcmp(type, "kqRWallFunction")                 == 0 ||
             strcmp(type, "omegaWallFunction")               == 0) {
      // skip valueFraction, go to value
      readGivenKeyword(ptrFile, "<VALUE>");
      readGivenKeyword(ptrFile, "<SCALAR>");
      skipLine(ptrFile);
      
      // Make sure you read a number equal to bcsize
      char buf[256];
      readline(ptrFile, buf, 256);
      E_Int nfaces; E_Int pos=0;
      readInt(buf, 1024, pos, nfaces);
      assert(nfaces == bcsize);
      
      printf("Info: foamread: non uniform list of size " SF_D_ ".\n", nfaces);
      skipLine(ptrFile); // (
      
      for (E_Int j = 0; j < nfaces; j++)
      {
        readline(ptrFile, buf, 1024); pos = 0;
        readDouble(buf, 1024, pos, bcf[j]);
      }
    }

    else if (strcmp(type, "compressible::alphatWallFunction")) {
      puts(type);
      readGivenKeyword(ptrFile, "UNIFORM");
      char buf[128];
      readWord(ptrFile, buf);
      E_Float val = strtod(buf, NULL);
      printf("Info: foamread: uniform " SF_F_ ".\n", val);
      for (E_Int j = 0; j < bcsize; j++) {
        bcf[j] = val;
      }
    }

    else {
      printf("Warning: foamread: unsupported boundary condition %s.\n", type);
      return 0;
    }

    readGivenKeyword(ptrFile, "}");
    skipLine(ptrFile);
    puts("");
  }

  fclose(ptrFile);

  return ncells;
}

E_Int K_IO::GenIO::readVectorField(char *file, FldArrayF& f, E_Int idx,
  E_Int *owner, const std::vector<FldArrayI *> &BCFaces,
  std::vector<FldArrayF *> &BCFields, const std::vector<E_Float> &delta,
  E_Int nifaces, const std::vector<E_Int> &indir)
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

  // Boundary fields
  E_Int ret = readGivenKeyword(ptrFile, "BOUNDARYFIELD");
  if (ret != 1) printf("INFO: foamread: cant find BOUNDARYFIELD.\n");
  assert(ret == 1);
  skipLine(ptrFile);
  
  E_Int nbnd = indir.size()-1;

  for (E_Int i = 0; i < nbnd; i++) 
  {
    //const char *bcname = BCNames[i];
    E_Int bcsize = indir[i+1] - indir[i];
    E_Int *bcfaces = BCFaces[0]->begin() + indir[i];
    E_Float *bcfx = BCFields[i]->begin(idx);
    E_Float *bcfy = BCFields[i]->begin(idx+1);
    E_Float *bcfz = BCFields[i]->begin(idx+2);

    // name
    char bcname[256];
    readWord(ptrFile, bcname);

    // type
    char type[256];
    ret = readGivenKeyword(ptrFile, "TYPE");
    assert(ret == 1);
    readWord(ptrFile, type);

    if (strcmp(type, "zeroGradient")  == 0 ||
        strcmp(type, "empty")         == 0 ||
        strcmp(type, "symmetryPlane") == 0) {
      puts(type);
      // Copy cell values
      for (E_Int j = 0; j < bcsize; j++) {
        E_Int face = bcfaces[j]-1;
        E_Int own = owner[face]-1;
        bcfx[j] = fldx[own];
        bcfy[j] = fldy[own];
        bcfz[j] = fldz[own];
      }
    }
    
    else if (strcmp(type, "fixedValue")              == 0 ||
             strcmp(type, "uniformNormalFixedValue") == 0) {
      puts(type);
      E_Int ret = readGivenKeyword(ptrFile, "UNIFORM");
      assert(ret == 1);
      ret &= readGivenKeyword(ptrFile, "(");
      assert(ret == 1);

      char buf[256];
      
      readWord(ptrFile, buf);
      E_Float valx = strtod(buf, NULL);

      readWord(ptrFile, buf);
      E_Float valy = strtod(buf, NULL);

      readWord(ptrFile, buf);
      E_Float valz = strtod(buf, NULL);
      printf("Info: foamread: uniform (" SF_F3_ ").\n", valx, valy, valz);

      for (E_Int j = 0; j < bcsize; j++) {
        bcfx[j] = valx;
        bcfy[j] = valy;
        bcfz[j] = valz;
      }
    }

    else if (strcmp(type, "noSlip") == 0) {
      puts(type);
      // Assumes non-moving wall
      for (E_Int j = 0; j < bcsize; j++) {
        bcfx[j] = 0.0;
        bcfy[j] = 0.0;
        bcfz[j] = 0.0;
      }
    }

    else if (strcmp(type, "calculated") == 0 ||
             strcmp(type, "freestreamVelocity") == 0) {
      readGivenKeyword(ptrFile, "<VECTOR>");
  
      skipLine(ptrFile);

      char buf[1024];
      readline(ptrFile, buf, 1024);
  
      E_Int nfaces; E_Int pos=0;
      readInt(buf, 1024, pos, nfaces);
      assert(nfaces == bcsize);

      printf("Info: foamread: " SF_D_ " faces.", nfaces);

      skipLine(ptrFile);

      for (E_Int j = 0; j < nfaces; j++)
      {
        readline(ptrFile, buf, 1024); pos = 1;
        readDouble(buf, 1024, pos, bcfx[j]);
        readDouble(buf, 1024, pos, bcfy[j]);
        readDouble(buf, 1024, pos, bcfz[j]);
      }
    }
    
    else {
      printf("Warning: foamread: unsupported boundary condition %s.\n", type);
      return 0;
    }

    readGivenKeyword(ptrFile, "}");
    skipLine(ptrFile);
    puts("");
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

static
void computeDeltaCoeffs(FldArrayI &cn, E_Int *owner, E_Int nifaces,
  E_Int nfaces, E_Float *fcenters, E_Float *fareas, E_Float *cx, E_Float *cy,
  E_Float *cz, std::vector<E_Float> &delta)
{
  E_Int nbfaces = nfaces - nifaces;
  for (E_Int i = 0; i < nbfaces; i++) {
    E_Int face = i + nifaces;
    E_Float *FC = &fcenters[3*face];
    E_Float *NORMAL = &fareas[3*face];

    // Vector pointing from cell center to face center
    E_Int own = owner[face]-1;
    E_Float l[3] = {FC[0]-cx[own], FC[1]-cy[own], FC[2]-cz[own]}; 

    // Face area
    E_Float area = K_MATH::norm(NORMAL, 3);

    delta[i] = K_MATH::dot(NORMAL, l, 3) / area;
  }
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
  vector<FldArrayF*>& BCFields,
  char*& varStringc,
  vector<FldArrayF*>& centerStructField,
  vector<FldArrayF*>& centerUnstructField)
{
  puts("\n");

  // Read points
  FldArrayF* f = new FldArrayF();
  foamReadPoints(file, *f);
  E_Float *px = f->begin(1);
  E_Float *py = f->begin(2);
  E_Float *pz = f->begin(3);

  // Read NGON
  FldArrayI cNGon; E_Int nfaces;
  foamReadFaces(file, nfaces, cNGon);

  // Allocate PE
  FldArrayI PE(nfaces, 2);
  foamReadOwner(file, PE);
  E_Int nifaces = foamReadNeighbour(file, PE);
  E_Int nbfaces = nfaces - nifaces;

  // compute NFace
  FldArrayI cNFace; E_Int nelts;
  K_CONNECT::connectFE2NFace(PE, cNFace, nelts);
  printf("Info: foamread: ncells: " SF_D_ ".\n", nelts);

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

  E_Int *indPH = cn->getIndPH();
  E_Int *nface = cn->getNFace();

  // Compute boundary face centers
  std::vector<E_Float> fcenters(3*nfaces);
  std::vector<E_Float> fareas(3*nfaces);
  K_METRIC::compute_face_centers_and_areas(*cn, px, py, pz, &fcenters[0],
    &fareas[0]);

  // Compute boundary cell centers
  std::vector<E_Float> cx(nelts);
  std::vector<E_Float> cy(nelts);
  std::vector<E_Float> cz(nelts);
  E_Int *owner = PE.begin(1);
  for (E_Int i = nifaces; i < nfaces; i++) {
    E_Int own = owner[i]-1;
    E_Int nf = -1;
    E_Int *pf = cn->getElt(own, nf, nface, indPH);
    E_Float vol = 0.0;
    K_METRIC::compute_cell_center_and_volume(own, nf, pf, px, py, pz,
      &fcenters[0], &fareas[0], owner, cx[own], cy[own], cz[own], vol);
  }

  // Compute deltaCoeffs only on boundary faces
  std::vector<E_Float> delta(nbfaces, 0.0);
  computeDeltaCoeffs(*cn, owner, nifaces, nfaces, &fcenters[0], &fareas[0],
    &cx[0], &cy[0], &cz[0], delta);

  // Read boundary
  std::vector<E_Int> indir;
  foamReadBoundary(file, BCFaces, BCNames, indir);

  // Read fields
  foamReadFields(file, centerUnstructField, nelts, varStringc, BCNames, BCFaces,
    BCFields, owner, delta, nifaces, indir);

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
#define SCALARFIELD 1
#define VECTORFIELD 2
#define TENSORFIELD 3

E_Int K_IO::GenIO::foamReadFields(char *file,
  std::vector<FldArrayF*>& centerUnstructField,
  E_Int ncells, char*& varStringc,
  const std::vector<char *> &BCNames,
  const std::vector<FldArrayI*> &BCFaces,
  std::vector<FldArrayF*> &BCFields, E_Int *owner,
  const std::vector<E_Float> &delta, E_Int nifaces,
  const std::vector<E_Int> &indir)
{
  // identify latest time output folder
  DIR *d;
  struct dirent *dir;
  d = opendir(file);

  char fullPath[1024];
  char path[1024];
  strcpy(fullPath, file);
  E_Float ret;
  E_Float max_time = -1.0;
  while ((dir = readdir(d))) {
    strcpy(path, fullPath);
    strcat(path, "/");
    strcat(path, dir->d_name);
    puts(path);
    if (dirExist(path)) {
      char *badptr = NULL;
      ret = strtod(dir->d_name, &badptr);
      if (*badptr == '\0') {
        max_time = std::max(max_time, ret);
      }
    }
  }
  closedir(d);

  // loop again and get fullPath
  bool found = false;
  char dir_name[260] = {0};
  dir_name[0] = '\0';
  d = opendir(file);
  while ((dir = readdir(d)) && !found) {
    strcpy(path, fullPath);
    strcat(path, "/");
    strcat(path, dir->d_name);
    if (dirExist(path)) {
      char *badptr = NULL;
      ret = strtod(dir->d_name, &badptr);
      if (*badptr == '\0' && ret == max_time) {
        sprintf(dir_name, "%s", dir->d_name);
        found = true;
        break;
      }
    }
  }
  closedir(d);

  if (max_time == -1.0) {
    printf("Warning: foamread: no time directory found. Skipping field read.\n");
    return 0;
  } else if (max_time == 0) {
    printf("Warning: foamread: skipping time directory 0.\n");
    return 0;
  } else {
    printf("INFO: foamread: reading fields from time directory %s\n", dir_name);
  }

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
      if (K_STRING::cmp(dir->d_name, 9, "processor") == 0) continue;
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, dir->d_name);
      
      // read only volScalarFields and volVectorFields
      FILE *fh = fopen(path, "r");
      assert(fh);
      E_Int ret = readGivenKeyword(fh, "VOLSCALARFIELD", "VOLVECTORFIELD", "VOLTENSORFIELD");
      if (ret == SCALARFIELD) { // volScalarField
        strcpy(field_name[size], dir->d_name);

        field_type[size] = SCALARFIELD;

        strcat(varStringc, field_name[size]);
        strcat(varStringc, ",");

        nflds++;
        size++;
      } else if (ret == VECTORFIELD) { // volVectorField
        strcpy(field_name[size], dir->d_name);

        field_type[size] = VECTORFIELD;

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
      } else if (ret == TENSORFIELD) { // volTensorField
        strcpy(field_name[size], dir->d_name);

        field_type[size] = TENSORFIELD;

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
      if (size == MAX_FIELDS) 
      {
        printf("Warning: foamread: trying to read more that maximum number of fields (%d). Aborting.\n", MAX_FIELDS);
        exit(1); // must return 
      }
      fclose(fh);
    }
  }
  closedir(d);

  printf("Info: foamread: nflds: " SF_D_ ".\n", nflds);

  // delete the last comma
  varStringc[strlen(varStringc)-1] = '\0';

  FldArrayF *F = new FldArrayF();
  F->malloc(ncells, nflds);

  BCFields.resize(indir.size()-1);
  for (size_t i = 0; i < indir.size()-1; i++)
    BCFields[i] = new FldArrayF(indir[i+1]-indir[i], nflds);

  d = opendir(fullPath);
  assert(d);

  E_Int idx = 1;
  for (E_Int fld = 0; fld < size; fld++) {
    if (field_type[fld] == SCALARFIELD) {
      printf("Info: foamread: reading scalar field %s\n", field_name[fld]);
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, field_name[fld]);
      E_Int ret = readScalarField(path, *F, idx, owner, BCFaces, BCFields,
        delta, nifaces, indir);
      if (ret != ncells) printf("INFO: foamread: wrong size for scalar field.\n");
      assert(ret == ncells);
      idx++;
    } else if (field_type[fld] == VECTORFIELD) {
      printf("INFO: foamread: reading vector field %s\n", field_name[fld]);
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, field_name[fld]);
      E_Int ret = readVectorField(path, *F, idx, owner, BCFaces, BCFields,
        delta, nifaces, indir);
      if (ret != ncells) printf("INFO: foamread: wrong size for vector field.\n");
      assert(ret == ncells);
      idx += 3;
    } else if (field_type[fld] == TENSORFIELD) {
      printf("Info: foamread: reading tensor field %s\n", field_name[fld]);
      char path[1028];
      strcpy(path, fullPath);
      strcat(path, "/");
      strcat(path, field_name[fld]);
      E_Int ret = readTensorField(path, *F, idx);
      if (ret != ncells) printf("INFO: foamread: wrong size for tensor field.\n");
      assert(ret == ncells);
      idx += 9;
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
  assert(ptrFile);

  E_Int ret;
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  if (ret != 1) printf("INFO: foamread: missing FOAMFILE.\n");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;
  E_Bool cont = true;
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

  printf("Info: foamread: points: " SF_D_ "\n", f.getSize());

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
  if (ret != 1) printf("INFO: foamread: missing FOAMFILE.\n");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;
  E_Bool cont = true;
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

  printf("Info: foamread: faces: " SF_D_ "\n", nfaces);

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
  if (ret != 1) printf("INFO: foamread: missing FOAMFILE.\n");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;

  E_Bool cont = true;
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
  if (ret != 1) printf("INFO: foamread: missing FOAMFILE.\n");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;

  E_Bool cont = true;
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
  E_Int nifaces; E_Int val;
  E_Int pos=0;
  readInt(buf, 1024, pos, nifaces);

  skipLine(ptrFile);

  for (E_Int i = 0; i < nifaces; i++)
  {
    ret = readInt(ptrFile, val); PE(i, 2) = val+1;
  }

  // tag exterior faces
  for (E_Int i = nifaces; i < PE.getSize(); i++)
  {
    PE(i,2) = 0; // exterior
  }

  printf("Info: foamread: internal faces: " SF_D_ "\n", nifaces);

  fclose(ptrFile);

  return nifaces;
}

#define BCSTRINGMAXSIZE 50
//=============================================================================
E_Int K_IO::GenIO::foamReadBoundary(char* file,
  std::vector<FldArrayI*>& BCFaces, std::vector<char*>& BCNames,
  std::vector<E_Int>& indir)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/boundary");
  FILE* ptrFile = fopen(fullPath, "r");
  if (ptrFile == NULL)
  {
    printf("Info: foamread: can not open constant/polyMesh/boundary.\n");
    return 1;
  }
  E_Int ret;
  ret = readGivenKeyword(ptrFile, "FOAMFILE");
  if (ret != 1) printf("INFO: foamread: missing FOAMFILE.\n");
  assert(ret == 1);
  ret = readGivenKeyword(ptrFile, "}");
  assert(ret == 1);

  // Passe comments
  char buf[1024]; E_Int l;

  E_Bool cont = true;
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
    strcat(bcnames[i], "@");
    strcat(bcnames[i], type[i]);

    // nFaces
    readGivenKeyword(ptrFile, "NFACES");
    readWord(ptrFile, buf);
    nFaces[i] = convertString2Int(buf);

    // startFace
    readGivenKeyword(ptrFile, "STARTFACE");
    readWord(ptrFile, buf);
    startFace[i] = convertString2Int(buf);

    skipLine(ptrFile);
    skipLine(ptrFile);
  }

  fclose(ptrFile);
  
  indir.resize(nBC+1);
  indir[0] = 0;

  /*
  for (E_Int i = 0; i < nBC; i++) {
    // Copy name
    char *NAME = (char *)malloc(BCSTRINGMAXSIZE);
    strcpy(NAME, bcnames[i]);

    // Alloc
    FldArrayI *faces = new FldArrayI(nFaces[i]);

    // Fill
    E_Int *ptr = faces->begin();
    for (E_Int j = 1; j <= nFaces[i]; j++)
      ptr[j] = startFace[i] + j;
    
    // Push
    BCNames.push_back(NAME);
    BCFaces.push_back(faces);
  }
  */


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
    indir[i+1] = indir[i] + nFaces[i];
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
  
  fprintf(ptrFile, SF_D_ "\n", npts);
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < npts; i++)
  {
    fprintf(ptrFile, "(%.15g %.15g %.15g)\n", x[i], y[i], z[i]);
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
// face indices (NGON)
//=============================================================================
E_Int K_IO::GenIO::foamWriteFaces(char* file, K_FLD::FldArrayI &cn,
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

  E_Int nfaces = faces.size();

  fprintf(ptrFile, SF_D_ "\n", nfaces);
  fprintf(ptrFile, "(\n");

  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();

  for (E_Int i = 0; i < nfaces; i++)
  {
    E_Int face = faces[i];
    E_Int stride = -1;
    E_Int *pN = cn.getFace(face, stride, ngon, indPG);

    fprintf(ptrFile, SF_D_ "(", stride);
    for (E_Int k = 0; k < stride; k++) fprintf(ptrFile, SF_D_ " ", pN[k]-1);
    fprintf(ptrFile, ")\n");
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);

  return 0;
}

//=============================================================================
// All faces (left=owner)
//=============================================================================
E_Int K_IO::GenIO::foamWriteOwner(char* file, const std::vector<E_Int> &owner,
  const std::vector<E_Int> &faces)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/owner");
  FILE* ptrFile = fopen(fullPath, "w");
  if (ptrFile == NULL) 
  {
    printf("foamWrite: can not open constant/polyMesh/owner\n");
    return 1;
  }
  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       labelList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      owner;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = owner.size();
  fprintf(ptrFile, SF_D_ "\n", nfaces);
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int own = owner[faces[i]];
    fprintf(ptrFile, SF_D_ "\n", own);
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}


//=============================================================================
// right=neighbour
//=============================================================================
E_Int K_IO::GenIO::foamWriteNeighbour(char* file,
  const std::vector<E_Int>& neigh, const std::vector<E_Int>& faces)
{
  char fullPath[1024];
  strcpy(fullPath, file);
  strcat(fullPath, "/constant/polyMesh/neighbour");
  FILE* ptrFile = fopen(fullPath, "w");
  if (ptrFile == NULL) 
  {
    printf("foamWrite: can not open constant/polyMesh/neighbour\n");
    return 1;
  }

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       labelList;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      neighbour;\n");
  fprintf(ptrFile, "}\n");

  E_Int nfaces = neigh.size();
  fprintf(ptrFile, SF_D_ "\n", nfaces);
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < nfaces; i++) 
  {
    E_Int nei = neigh[faces[i]];
    fprintf(ptrFile, SF_D_ "\n", nei);
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
  if (ptrFile == NULL) 
  {
    printf("foamWrite: can not open constant/polyMesh/boundary\n");
    return 1;
  }

  fprintf(ptrFile, "FoamFile\n");
  fprintf(ptrFile, "{\n");
  fprintf(ptrFile, "    version     2.0;\n");
  fprintf(ptrFile, "    format      ascii;\n");
  fprintf(ptrFile, "    arch        \"LSB;label=32;scalar=64\";\n");
  fprintf(ptrFile, "    class       polyBoundaryMesh;\n");
  fprintf(ptrFile, "    location    \"constant/polyMesh\";\n");
  fprintf(ptrFile, "    object      boundary;\n");
  fprintf(ptrFile, "}\n");
  
  E_Int nbcs = bc_names.size();
  fprintf(ptrFile, SF_D_ "\n", nbcs);
  fprintf(ptrFile, "(\n");

  for (E_Int i = 0; i < nbcs; i++)
  {
    char *token = strtok(bc_names[i], "@"); // token is first part of bc_names
    
    char name[256];
    strcpy(name, token);

    fprintf(ptrFile, "    %s\n", token);
    fprintf(ptrFile, "    {\n");

    token = strtok(NULL, "@"); // token is now second part of bc_names (CGNS type if any)

    if (token == NULL) {
      fprintf(stderr, "No type for BC %s, defaulting to wall.\n", bc_names[i]);
      fprintf(ptrFile, "        type            %s;\n", "wall");
      fprintf(ptrFile, "        physicalType    %s;\n", "wall");
    } else {
      fprintf(ptrFile, "        type            %s;\n", token);
      fprintf(ptrFile, "        physicalType    %s;\n", token);
    }
    
    fprintf(ptrFile, "        nFaces          " SF_D_ ";\n", bc_nfaces[i]);
    fprintf(ptrFile, "        startFace       " SF_D_ ";\n", bc_startfaces[i]);
    fprintf(ptrFile, "    }\n");
  }
  fprintf(ptrFile, ")\n");
  fclose(ptrFile);
  return 0;
}

static
E_Int get_neighbor(E_Int cell, E_Int face, E_Int *owner, E_Int *neigh)
{
  assert(cell == owner[face] || cell == neigh[face]);
  return cell == owner[face] ? neigh[face] : owner[face];
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
  vector< vector<E_Int> >& eltTypes,
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
    E_Int et = eltTypes[i][0];
    if (et == 8) { no = i; break; }
  }
  if (no == -1)
  {
    printf("Warning: foamwrite: no NGON zone found in input. Nothing written.\n");
    return 1;
  }

  FldArrayF& field = *unstructField[no];
  FldArrayI& cn = *connect[no];

  // Write points
  foamWritePoints(file, field);

  // Build owners and neighbors
  E_Float *px = field.begin(posx);
  E_Float *py = field.begin(posy);
  E_Float *pz = field.begin(posz);
  K_CONNECT::orient_boundary_ngon(px, py, pz, cn);

  E_Int nfaces = cn.getNFaces();
  std::vector<E_Int> owner(nfaces), neigh(nfaces);

  K_CONNECT::build_parent_elements_ngon(cn, &owner[0], &neigh[0]);

  // Renumber faces
  std::vector<E_Int> faces;
  std::vector<E_Int> marked(nfaces, 0);

  E_Int PGi;
  std::vector<E_Int> neis;
  std::vector<E_Int> pgs;

  E_Int ncells = cn.getNElts();
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();

  for (E_Int PHi = 0; PHi < ncells; PHi++) 
  {
    E_Int stride = -1;
    const E_Int *pF = cn.getElt(PHi, stride, nface, indPH);

    neis.clear();
    pgs.clear();
    
    for (E_Int j = 0; j < stride; j++) 
    {
      PGi = pF[j] - 1;
      
      assert(owner[PGi] != -1);

      if (neigh[PGi] == -1) continue;

      if (!marked[PGi]) {
        neis.push_back(get_neighbor(PHi, PGi, &owner[0], &neigh[0]));
        pgs.push_back(PGi);
        marked[PGi] = 1;
        if (owner[PGi] > neigh[PGi]) {
          E_Int stride = -1;
          E_Int *pn = cn.getFace(PGi, stride, ngon, indPG);
          std::swap(owner[PGi], neigh[PGi]);
          std::reverse(pn, pn+stride);
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

  //E_Int ninternal_faces = faces.size();
  //printf("Internal faces = " SF_D_ "\n", ninternal_faces);

  // BC
  E_Int BCFacesSize = 0;
  if (PyList_Check(BCFaces)) BCFacesSize = PyList_Size(BCFaces);
  //printf("BCFacesSize = " SF_D_ "\n", BCFacesSize);

  std::vector<E_Int> start_face_per_bc;
  std::vector<E_Int> nfaces_per_bc;
  std::vector<char*> name_per_bc;

  if (BCFacesSize > 0) {
    E_Int indFace;
    IMPORTNUMPY;

    PyObject* BCs = PyList_GetItem(BCFaces, 0);
    E_Int size = PyList_Size(BCs);

    if (size == 0) 
    {
      printf("Warning: foamwrite: requires boundary patches.\n");
      return 1;
    }

    E_Int np;
    E_Int *ptr;
    char *name;
    for (E_Int j = 0; j < size/2; j++) {
      name = NULL;
      PyObject *o = PyList_GetItem(BCs, 2*j);
      
#if PY_VERSION_HEX >= 0x03000000
      if (PyUnicode_Check(o)) name = (char *)PyUnicode_AsUTF8(o);
#else
      if (PyString_Check(o)) name = PyString_AsString(o);
#endif
      else 
      {
        printf("Error: foamwrite: bad bcname no " SF_D_ ".\n", j);
        return 1;
      }

      name_per_bc.push_back(0);
      name_per_bc[j] = new char[strlen(name) + 1];
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

  //printf("total faces: %zu\n", faces.size());
  //printf("cn faces: " SF_D_ "\n", nfaces);

  foamWriteOwner(file, owner, faces);
  foamWriteNeighbour(file, neigh, faces);
  foamWriteFaces(file, cn, faces);
  foamWriteBoundary(file, name_per_bc, nfaces_per_bc, start_face_per_bc);

  for (size_t i = 0; i < name_per_bc.size(); i++)
    delete [] name_per_bc[i];

  return 0;
}
