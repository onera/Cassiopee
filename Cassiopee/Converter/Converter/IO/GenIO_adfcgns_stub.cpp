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

// Binary ADF CGNS file support
# include "GenIO.h"

//=============================================================================
/* 
   adfcgnsread
*/
//=============================================================================
E_Int K_IO::GenIO::adfcgnsread(char* file, PyObject*& tree, int skeleton,
                               int maxFloatSize, int maxDepth)
{
  printf("Error: Converter has been installed without CGNS/ADF support.\n");
  printf("Error: please install libcgns first for CGNS/ADF support.\n");
  return 0;
}

//=============================================================================
/* 
   adfcgnswrite
*/
//=============================================================================
E_Int K_IO::GenIO::adfcgnswrite(char* file, PyObject* tree)
{
  printf("Error: Converter has been installed without CGNS/ADF support.\n");
  printf("Error: please install libcgns first for CGNS/ADF support.\n");
  return 0;
}

//=============================================================================
PyObject* K_IO::GenIO::adfcgnsReadFromPaths(char* file, PyObject* paths,
                                            E_Int maxFloatSize, E_Int maxDepth)
{ 
  printf("Error: Converter has been installed without CGNS/ADF support.\n");
  printf("Error: please install libcgns first for CGNS/ADF support.\n");
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::adfcgnsWritePaths(char* file, PyObject* treeList, 
                                     PyObject* paths, E_Int recursive, E_Int mode)
{
  printf("Error: Converter has been installed without CGNS/ADF support.\n");
  printf("Error: please install libcgns first for CGNS/ADF support.\n");
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::adfcgnsDeletePaths(char* file, PyObject* paths)
{
  printf("Error: Converter has been installed without CGNS/ADF support.\n");
  printf("Error: please install libcgns first for CGNS/ADF support.\n");
  return 0;
}
