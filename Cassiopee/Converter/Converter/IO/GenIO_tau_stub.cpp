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

// Binary NETCDF3 TAU (.grid) file support
# include "GenIO.h"

//=============================================================================
/* 
   tauread
*/
//=============================================================================
E_Int K_IO::GenIO::tauread(char* file, PyObject*& tree)
{
  printf("Error: Converter has been installed without NETCDF/TAU support.\n");
  printf("Error: please install netcdf first for NETCDF/TAU support.\n");
  return 0;
}

//=============================================================================
/* 
   tauwrite
*/
//=============================================================================
E_Int K_IO::GenIO::tauwrite(char* file, PyObject* tree)
{
  printf("Error: Converter has been installed without NETCDF/TAU support.\n");
  printf("Error: please install netcdf first for NETCDF/TAU support.\n");
  return 0;
}

