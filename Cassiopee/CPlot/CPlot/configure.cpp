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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"

#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

using namespace K_FLD;

//========================================================================
/*
  configure : Permet de choisir entre un rendu en display list ou bien en
              VBO. Par defaut, le rendu est en VBO
*/
//========================================================================
PyObject* K_CPLOT::configure(PyObject* self, PyObject* args)
{
  int useRender;
  if (!PyArg_ParseTuple(args, "i", &useRender)) return NULL;
  switch(useRender) 
  {
  case 0:
    Data::_renderID = Data::Direct;
    break;
  case 1:
    Data::_renderID = Data::VBO;
    break;
  case 2:
    Data::_renderID = Data::DL;
    break;
  default:
    std::cerr << "Render id pas encore disponible.\n";
    PyErr_SetString(PyExc_TypeError, "Render ID not yet implemanted");
    return NULL;
  }
  Py_RETURN_NONE;
}
