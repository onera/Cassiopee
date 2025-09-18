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

#include "cplot.h"
#include "Data.h"
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

//=============================================================================
/* set file name in state */
//=============================================================================
PyObject* K_CPLOT::setFileName(PyObject* self, PyObject* args)
{
  char* fileName;
  if (!PYPARSETUPLE_(args, S_, &fileName)) return NULL;

  // Recuperation du container de donnees
  Data* d = Data::getInstance();

  // Switch - Dangerous zone protegee par _state.lock
  d->ptrState->syncDisplay();
  strcpy(d->ptrState->file, fileName);
  return Py_BuildValue("l", 0);
}
