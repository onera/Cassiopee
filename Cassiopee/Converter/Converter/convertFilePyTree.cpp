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
  E_Int readMode;
  if (!PYPARSETUPLE_(args, SS_ OOOO_ I_, &fileName, &format, &skeletonData, 
                     &dataShape, &links, &skipTypes, &readMode))
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
                                                  maxDepth, readMode, skipTypes);
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
  char* format; PyObject* paths; E_Int readMode;
  if (!PYPARSETUPLE_(args, SS_ O_ I_, &fileName, &format, &paths, &readMode))
    return NULL;

  E_Int l = strlen(format);
  char* myFormat = new char [l+1]; strcpy(myFormat, format);
  if (strcmp(myFormat, "bin_cgns") == 0) strcpy(myFormat, "bin_adf");
  
  PyObject* ret = NULL;
  ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPaths(fileName, paths, 1.e6, -1, readMode);
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
  if (!PYPARSETUPLE_(args, O_ SS_ O_, &t, &fileName, &format, &links)) return NULL;

  printf("Writing %s (%s)...", fileName, format);
  fflush(stdout);

  if (strcmp(format, "bin_cgns") == 0)
    K_IO::GenIO::getInstance()->hdfcgnswrite(fileName, t, links);
  else if (strcmp(format, "bin_hdf") == 0)
    K_IO::GenIO::getInstance()->hdfcgnswrite(fileName, t, links);  
  else if (strcmp(format, "bin_adf") == 0)
    K_IO::GenIO::getInstance()->adfcgnswrite(fileName, t);
  else
    K_IO::GenIO::getInstance()->hdfcgnswrite(fileName, t, links);
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
  E_Int readMode;
  if (!PYPARSETUPLE_(args, SS_ OOO_ I_, &fileName, &format, &skeletonData,
                        &mpi4pyCom, &filter, &readMode))
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
  ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPathsPartial(fileName, readMode, filter, mpi4pyCom);
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

