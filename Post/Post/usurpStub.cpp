/*    
    Copyright 2013-2019 Onera.

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
# include "post.h"
# include "stdio.h"

using namespace std;

// ============================================================================
/* USURP */
// ============================================================================
PyObject* K_POST::usurpF(PyObject* self, PyObject* args)
{ 
  PyObject* blkArrays; PyObject* ibArrays;
  
  if (!PyArg_ParseTuple(args, "OO", &blkArrays, &ibArrays)) return NULL;
  
  printf("Warning: usurp: this function is not available in this distribution (no f90 on this platform).\n"); 
  
  Py_INCREF(Py_None);
  return Py_None;
}
