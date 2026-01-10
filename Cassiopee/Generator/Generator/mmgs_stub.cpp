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

// Remaillage surfacique avec mmgs

#include "generator.h"
#include "MMGS/mmgs.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC; 

// ============================================================================
/* MMGS
   IN: maillage TRI
   IN: eventuellement metric ou solution
   IN: 
   IN: 
   OUT: maillage TRI remaille. */
// ============================================================================
PyObject* K_GENERATOR::mmgs(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_TypeError,
		  "mmgs: Generator was not installed with mmgs.");
  return NULL;
}
