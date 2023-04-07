/*    
    Copyright 2013-2023 Onera.

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
// tester for openacc

#include "kcore.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>


//==============================================================================
PyObject* K_KCORE::testerAcc(PyObject* self, PyObject* args)
{
    E_Int n, i;

    n = 100;
    double* c = new double [n];

    #pragma acc data copy(c)
    {
        #pragma acc parallel loop
        for (i = 0; i < n; i++) c[i] = 2;
    }

    printf("c=%f\n", c[0]);

    delete [] c; 
    Py_INCREF(Py_None);
    return Py_None;
}
