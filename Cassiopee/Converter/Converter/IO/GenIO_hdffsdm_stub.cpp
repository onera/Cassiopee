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

// Binary HDF FSDM (.h5) file support
# include "GenIO.h"

//=============================================================================
/* 
   hdffsdmread
*/
//=============================================================================
E_Int K_IO::GenIO::hdffsdmread(char* file, PyObject*& tree)
{
  printf("Error: Converter has been installed without HDF/FSDM support.\n");
  printf("Error: please install libhdf5 first for HDF/FSDM support.\n");
  return 0;
}

//=============================================================================
/* 
   hdffsdmwrite
*/
//=============================================================================
E_Int K_IO::GenIO::hdffsdmwrite(char* file, PyObject* tree, PyObject* links)
{
  printf("Error: Converter has been installed without HDF/FSDM support.\n");
  printf("Error: please install libhdf5 first for HDF/FSDM support.\n");
  return 0;
}

