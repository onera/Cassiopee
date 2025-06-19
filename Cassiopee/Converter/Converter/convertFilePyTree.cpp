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
// conversion file CGNS / pyTrees CGNS

#ifdef _MPI
#if defined(_WIN64)
#define __int64 long long
#endif
#endif
#include "converter.h"
#include "kcore.h"
#include "IO/GenIO.h"
#ifdef _HDF5
#include "hdf5.h" // necessary for definition of H5_HAVE_PARALLEL
#endif

// ============================================================================
/* Convert file to pyTree */
// ============================================================================
PyObject* K_CONVERTER::convertFile2PyTree(PyObject* self, PyObject* args)
{
  char* fileName; char* format; 
  PyObject* skeletonData; PyObject* dataShape; 
  PyObject* links; PyObject* skipTypes;
  E_Int readIntMode;
  if (!PYPARSETUPLE_(args, SS_ OOOO_ I_, &fileName, &format, &skeletonData, 
                     &dataShape, &links, &skipTypes, &readIntMode))
    return NULL;
  
  if (dataShape == Py_None) { dataShape = NULL; }
  if (links == Py_None) { links = NULL; }
  if (skipTypes == Py_None) { skipTypes = NULL; }

  E_Int l = strlen(format);
  char* myFormat = new char [l+1]; strcpy(myFormat, format);
  if (strcmp(myFormat, "bin_cgns") == 0) strcpy(myFormat, "bin_hdf");

  // Get skeleton data if any
  int skeleton=0; int maxFloatSize=5; int maxDepth=-1;
  if (skeletonData != Py_None)
  {
    skeleton = 1;
    if (PyList_Check(skeletonData) == true)
    {
      if (PyList_Size(skeletonData) == 2)
      {
        maxFloatSize = (int)PyInt_AsLong(PyList_GetItem(skeletonData, 0));
        maxDepth = (int)PyInt_AsLong(PyList_GetItem(skeletonData, 1));
      }
    }
  }

  if (skeletonData == Py_None) printf("Reading %s (%s)...", fileName, myFormat);
  else printf("Reading %s (%s, skeleton)...", fileName, myFormat);
  fflush(stdout);

  PyObject* tree; E_Int ret;
  if (strcmp(myFormat, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsread(fileName, tree, skeleton, maxFloatSize, maxDepth);
  else if (strcmp(myFormat, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsread(fileName, tree, dataShape, links, skeleton, maxFloatSize, 
                                                  maxDepth, readIntMode, skipTypes);
  else
    ret = K_IO::GenIO::getInstance()->adfcgnsread(fileName, tree, skeleton, maxFloatSize, maxDepth);
  printf("done.\n");
  delete [] myFormat;

  if (ret == 1)
  {
    PyErr_SetString(PyExc_IOError,
                    "convertFile2PyTree: fail to read.");
    return NULL;
  }

  return tree;
}

// ============================================================================
/* Convert file to pyTree */
// ============================================================================
PyObject* K_CONVERTER::convertFile2PyTreeFromPath(PyObject* self, PyObject* args)
{
  char* fileName;
  char* format; PyObject* paths; E_Int readIntMode;
  if (!PYPARSETUPLE_(args, SS_ O_ I_, &fileName, &format, &paths, &readIntMode))
    return NULL;

  E_Int l = strlen(format);
  char* myFormat = new char [l+1]; strcpy(myFormat, format);
  if (strcmp(myFormat, "bin_cgns") == 0) strcpy(myFormat, "bin_adf");
  
  PyObject* ret = NULL;
  ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPaths(fileName, paths, 1.e6, -1, readIntMode);
  printf("done.\n");
 
  delete [] myFormat;

  return ret;
}

// ============================================================================
/* Convert pyTree to file */
// ============================================================================
PyObject* K_CONVERTER::convertPyTree2File(PyObject* self, PyObject* args)
{
  char* fileName; char* format;
  PyObject* t; PyObject* links;
  E_Int isize; E_Int rsize;
  if (!PYPARSETUPLE_(args, O_ SS_ O_ II_, &t, &fileName, &format, &links, 
                    &isize, &rsize)) return NULL;

  printf("Writing %s (%s)...", fileName, format);
  fflush(stdout);

  int writeIntMode = 0;
  if (isize == 4) writeIntMode = 1; // force best i4/i8 write
  int writeRealMode = 0;
  if (rsize == 4) writeRealMode = 1; // force r4 write

  if (strcmp(format, "bin_cgns") == 0)
    K_IO::GenIO::getInstance()->hdfcgnswrite(fileName, t, links, 
      writeIntMode, writeRealMode);
  else if (strcmp(format, "bin_hdf") == 0)
    K_IO::GenIO::getInstance()->hdfcgnswrite(fileName, t, links,
      writeIntMode, writeRealMode);  
  else if (strcmp(format, "bin_adf") == 0)
    K_IO::GenIO::getInstance()->adfcgnswrite(fileName, t);
  else
    K_IO::GenIO::getInstance()->hdfcgnswrite(fileName, t, links,
      writeIntMode, writeRealMode);
  printf("done.\n");

  Py_INCREF(Py_None);
  return Py_None;
}

// ============================================================================
/* Lit des noeuds definis dans Filter (partiellement) - hdf only */
// ============================================================================
PyObject* K_CONVERTER::convertFile2PartialPyTree(PyObject* self, PyObject* args)
{
  char* fileName;
  char* format;
  PyObject* skeletonData;
  PyObject* mpi4pyCom;
  PyObject* filter; // dictionnaire des slices
  E_Int readIntMode;
  if (!PYPARSETUPLE_(args, SS_ OOO_ I_, &fileName, &format, &skeletonData,
                        &mpi4pyCom, &filter, &readIntMode))
    return NULL;
  
  E_Int l = strlen(format);
  char* myFormat = new char [l+1]; strcpy(myFormat, format);
  if (strcmp(myFormat, "bin_cgns") == 0) strcpy(myFormat, "bin_hdf");
  if (strcmp(myFormat, "bin_hdf") != 0)
  { printf("convertFile2PartialPyTree: only for HDF.\n"); return NULL; }

  //if (skeletonData == Py_None) printf("Reading %s (%s)...", fileName, myFormat);
  //else printf("Reading %s (%s, skeleton)...", fileName, myFormat);
  //fflush(stdout);
  printf("Reading %s (%s, partial)...", fileName, myFormat);

  PyObject* ret;
  ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPathsPartial(fileName, readIntMode, filter, mpi4pyCom);
  printf("done.\n");
  delete [] myFormat;
  return ret;
}

// ============================================================================
/* Convert file to partial pyTree - hdf only */
// ============================================================================
PyObject* K_CONVERTER::convertPyTree2FilePartial(PyObject* self, PyObject* args)
{
  char* fileName;
  char* format;
  int skeleton;
  PyObject* mpi4pyCom;
  PyObject* Filter;
  PyObject* t;
  PyObject* skeletonData;
  if (!PYPARSETUPLE_(args, O_ SS_ OOO_, &t, &fileName, &format, &skeletonData, &mpi4pyCom, &Filter)) return NULL;

  //printf("Writing (partial) %s (%s)...", fileName, format);
  //fflush(stdout); 

  if (skeletonData != Py_None) {skeleton = 0;}
  else {skeleton = 1;}
  
#ifdef _MPI
  // CB: Je conserve l'option skeleton si compile avec MPI mais sans HDF parallele
#if defined(H5_HAVE_PARALLEL)
  /** Dans le cas MPI+HDF parallele, on cree les dataSpaces en parallele - Pas besoin de Skeleton **/
  skeleton = 0;
#endif
  K_IO::GenIO::getInstance()->hdfcgnsWritePathsPartial(fileName, t, Filter, skeleton, mpi4pyCom);

#else
  /* En sequentiel */
  // skeleton = 1;
  K_IO::GenIO::getInstance()->hdfcgnsWritePathsPartial(fileName, t, Filter, skeleton, mpi4pyCom);
#endif

  printf("done.\n");

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Lit les paths specifies dans le fichier file (seult HDF).
// Retourne une liste d'objets pythons contenant les noeuds pointes par les
// chemins 
//=============================================================================
PyObject* K_CONVERTER::readPyTreeFromPaths(PyObject* self, PyObject* args)
{
  char* fileName; char* format; 
  E_Int maxFloatSize; E_Int maxDepth; E_Int readIntMode; 
  PyObject* paths; PyObject* skipTypes; PyObject* dataShape; PyObject* mpi4pyCom;
  if (!PYPARSETUPLE_(args, S_ O_ S_ III_ OOO_,
                     &fileName, &paths, &format, 
                     &maxFloatSize, &maxDepth, &readIntMode, 
                     &dataShape, &skipTypes, &mpi4pyCom)) return NULL;
  
  if (skipTypes == Py_None) skipTypes = NULL;
  if (dataShape == Py_None) dataShape = NULL;
  PyObject* ret = NULL;
  if (K_STRING::cmp(format, "bin_cgns") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPaths(fileName, paths, maxFloatSize, maxDepth, readIntMode, dataShape, skipTypes, mpi4pyCom);
  else if (K_STRING::cmp(format, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPaths(fileName, paths, maxFloatSize, maxDepth, readIntMode, dataShape, skipTypes, mpi4pyCom);
  else if (K_STRING::cmp(format, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsReadFromPaths(fileName, paths, maxFloatSize, maxDepth);
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "readPyTreeFromPaths: unknown file format.");
    return NULL;
  }
  return ret;
} 

//=============================================================================
// Ecrit les paths specifies dans le fichier file (ADF/HDF).
//=============================================================================
PyObject* K_CONVERTER::writePyTreePaths(PyObject* self, PyObject* args)
{
  char* fileName; char* format; E_Int maxDepth; E_Int mode;
  PyObject* paths; PyObject* nodeList; PyObject* links;
  E_Int isize, rsize;
  if (!PYPARSETUPLE_(args, S_ OO_ S_ II_ O_ II_,
                     &fileName, &nodeList, &paths, &format, &maxDepth, &mode, &links, &isize, &rsize))
    return NULL;
  
  if (links == Py_None) { links = NULL; }
  
  E_Int ret = 1;

  if (K_STRING::cmp(format, "bin_cgns") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsWritePaths(fileName, nodeList, paths, links, maxDepth, mode, isize, rsize);
  else if (K_STRING::cmp(format, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsWritePaths(fileName, nodeList, paths, links, maxDepth, mode, isize, rsize);
  else if (K_STRING::cmp(format, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsWritePaths(fileName, nodeList, paths, maxDepth, mode);
  if (ret == 1) return NULL; // exceptions deja levees

  Py_INCREF(Py_None);
  return Py_None;
} 

//=============================================================================
// Delete des paths du fichier (ADF/HDF).
//=============================================================================
PyObject* K_CONVERTER::deletePyTreePaths(PyObject* self, PyObject* args)
{
  char* fileName; char* format;
  PyObject* paths;
  if (!PYPARSETUPLE_(args, S_ O_ S_, &fileName, &paths, &format))
    return NULL;
  
  E_Int ret = 1;
  if (K_STRING::cmp(format, "bin_cgns") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsDeletePaths(fileName, paths);
  else if (K_STRING::cmp(format, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsDeletePaths(fileName, paths);
  else if (K_STRING::cmp(format, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsDeletePaths(fileName, paths);
  if (ret == 1) return NULL; // exceptions deja levees

  Py_INCREF(Py_None);
  return Py_None;
}
