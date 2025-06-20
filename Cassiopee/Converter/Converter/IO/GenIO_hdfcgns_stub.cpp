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

// Binary HDF CGNS file support
# include "GenIO.h"

//=============================================================================
/* 
   hdfcgnsread
*/
//=============================================================================
E_Int K_IO::GenIO::hdfcgnsread(char* file, PyObject*& tree, PyObject* dataShape, 
                               PyObject* links, int skeleton, int maxFloatSize, 
                               int maxDepth, int readIntMode, PyObject* skipTypes)
{
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return 0;
}

//=============================================================================
PyObject* K_IO::GenIO::hdfcgnsReadFromPaths(char* file, PyObject* paths,
                                            E_Int maxFloatSize, E_Int maxDepth, E_Int readIntMode,
                                            PyObject* dataShape,
                                            PyObject* skipTypes, PyObject* mpi4pyCom)
{ 
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return NULL;
}
//==============================================================================
PyObject* K_IO::GenIO::hdfcgnsReadFromPathsPartial(char * file, E_Int readIntMode,
                                                   PyObject* Filter,
                                                   PyObject* mpi4pyCom)
{
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return NULL;
}

//=============================================================================
/* 
   hdfcgnswrite
*/
//=============================================================================
E_Int K_IO::GenIO::hdfcgnswrite(char* file, PyObject* tree, PyObject* links,
                                int writeIntMode, int writeRealMode)
{
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::hdfcgnsWritePaths(char* file, PyObject* treeList, 
                                     PyObject* paths,  PyObject* links, 
                                     E_Int maxDepth, E_Int mode,
                                     E_Int isize, E_Int rsize)
{
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return 0;
}
//==============================================================================
E_Int K_IO::GenIO::hdfcgnsWritePathsPartial(char* file, PyObject* tree,
                                            PyObject* Filter, int skeleton,
                                            PyObject* mpi4pyCom)
{
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return 0;
}
//===============================================================================
E_Int K_IO::GenIO::hdfcgnsDeletePaths(char* file,
                                      PyObject* paths)
{
  printf("Error: Converter has been installed without CGNS/HDF support.\n");
  printf("Error: please install libhdf5 first for CGNS/HDF support.\n");
  return 0;
} 
