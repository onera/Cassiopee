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

// Les differents tests
# define EXTARITH       0

//==============================================================================
PyObject* K_KCORE::testerAcc(PyObject* self, PyObject* args)
{
    E_Int n, n2, i, j, k;

    n= 100; n2= n*n;

    double* a = new double [n2];
    double* b = new double [n2];
    double* c = new double [n2];

    for (E_Int i = 0; i < n2 ; ++i)
    { a[i]=1; b[i]=100; c[i]=0; }

    #pragma acc data copyin(a,b) copy(c)
    {
        #pragma acc parallel loop collapse(3)
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                for (k = 0; k < n; ++k)
                {
                    c[i + n*j] = a[i + k*n] * b[ k +j*n] + c[i+j*n];
                }
            }
       }
    }

    printf("c=%f\n", c[0]);

    delete [] a; delete [] b; delete [] c; 
    Py_INCREF(Py_None);
    return Py_None;
}
