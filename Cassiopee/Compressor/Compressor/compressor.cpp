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
#define K_ARRAY_UNIQUE_SYMBOL
#include "compressor.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pycompressor [] =
{
  {"deltaIndex", K_COMPRESSOR::deltaIndex, METH_VARARGS},
  {"writeUnsteadyCoefs", K_COMPRESSOR::writeUnsteadyCoefs, METH_VARARGS},
  {"compressCellN", K_COMPRESSOR::py_cellN_compress, METH_VARARGS},
  {"uncompressCellN", K_COMPRESSOR::py_cellN_uncompress, METH_VARARGS},
  {"compressFpc", K_COMPRESSOR::py_fpc_compress, METH_VARARGS},
  {"uncompressFpc", K_COMPRESSOR::py_fpc_uncompress, METH_VARARGS},
  {"compressIndices", K_COMPRESSOR::py_indices_compress, METH_VARARGS},
  {"uncompressIndices", K_COMPRESSOR::py_indices_uncompress, METH_VARARGS},
  {"compressNGonIndices", K_COMPRESSOR::py_ngon_indices_compress, METH_VARARGS},
  {"uncompressNGonIndices", K_COMPRESSOR::py_ngon_indices_uncompress, METH_VARARGS},
  {"compressFpc", K_COMPRESSOR::py_fpc_compress, METH_VARARGS},
  {"uncompressFpc", K_COMPRESSOR::py_fpc_uncompress, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "compressor",
        NULL,
        -1,
        Pycompressor
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_compressor();
  PyMODINIT_FUNC PyInit_compressor()
#else
  PyMODINIT_FUNC initcompressor();
  PyMODINIT_FUNC initcompressor()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("compressor", Pycompressor);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}
