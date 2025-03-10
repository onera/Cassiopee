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

# ifndef _COMPRESSOR_COMPRESSOR_H_
# define _COMPRESSOR_COMPRESSOR_H_

# include "kcore.h"

namespace K_COMPRESSOR
{ 
  PyObject* deltaIndex(PyObject* self, PyObject* args);
  PyObject* writeUnsteadyCoefs(PyObject* self, PyObject* args);
  PyObject* py_cellN_compress(PyObject *self, PyObject *args);
  PyObject* py_cellN_uncompress(PyObject *self, PyObject *args);
  PyObject* py_fpc_compress(PyObject *self, PyObject *args);
  PyObject* py_fpc_uncompress(PyObject *self, PyObject *args);
  PyObject* py_indices_compress(PyObject* self, PyObject* args);
  PyObject* py_indices_uncompress(PyObject* self, PyObject* args);
  PyObject* py_ngon_indices_compress(PyObject* self, PyObject* args);
  PyObject* py_ngon_indices_uncompress(PyObject* self, PyObject* args);
}
#endif
