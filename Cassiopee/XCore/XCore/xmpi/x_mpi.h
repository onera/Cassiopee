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
#ifndef _XCORE_XMPI_X_MPI_H
#define _XCORE_XMPI_X_MPI_H

#  if defined( _MPI )
// Pour etre compatible avec Microsoft MPI sous MSys 2
#    if defined(_WIN64)
#    include <cstdint>
typedef __int64 int64_t;
#    endif
//
#  include <mpi.h>
#  endif

#endif
