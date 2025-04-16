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

// Convert file tau / pyTree CGNS

#include "converter.h"
#include "kcore.h"
#include "IO/GenIO.h"

// ============================================================================
/* Convert file to pyTree */
// ============================================================================
PyObject* K_CONVERTER::convertFile2PyTreeTau(PyObject* self, PyObject* args)
{
  char* fileName; char* myFormat;
  if (!PYPARSETUPLE_(args, SS_, &fileName, &myFormat)) return NULL;

  PyObject* tree = NULL;
  printf("Reading %s (%s)...", fileName, myFormat);
  E_Int ret = K_IO::GenIO::getInstance()->tauread(fileName, tree);
  printf("done.\n");

  if (ret == 1)
  {
    PyErr_SetString(PyExc_IOError, "convertFile2PyTree: fail to read.");
    return NULL;
  }

  return tree;
}
