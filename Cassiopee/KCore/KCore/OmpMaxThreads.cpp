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

//=============================================================================
/* Retourne le nbre max de threads utilise par openMP (OMP_NUM_THREADS) */
//=============================================================================
PyObject* K_KCORE::getOmpMaxThreads(PyObject* self, PyObject* args)
{
  int i = 1;
#ifdef _OPENMP
  i = omp_get_max_threads();
#endif
  return Py_BuildValue("i", i);
}

//=============================================================================
/* Change le nbre max de threads utilise par openMP (OMP_NUM_THREADS) */
//=============================================================================
PyObject* K_KCORE::setOmpMaxThreads(PyObject* self, PyObject* args)
{
  int i;
  if (!PyArg_ParseTuple(args, "i", &i)) return NULL;
#ifdef _OPENMP
  omp_set_num_threads(i);
#endif
  return Py_None;
}
