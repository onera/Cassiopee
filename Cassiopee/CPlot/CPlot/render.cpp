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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

//=============================================================================
/* 
   Force render
 */
//=============================================================================
PyObject* K_CPLOT::render(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  E_Int count = 5000;
  while (d->ptrState->render == 1 && count > 0)
  {
    static struct timeval tv;
    tv.tv_sec = 0;
    tv.tv_usec = 5;
    select(0, NULL, NULL, NULL, &tv);
    count--;
  }

  d->ptrState->render = 1;
  if (d->ptrState->offscreen > 0) d->ptrState->shootScreen = 1;
  return Py_BuildValue("i", KSUCCESS);
}
