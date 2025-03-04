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
// tester for accelerators (gpu with openacc or openmp)

#include "kcore.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

#ifdef _OPENACC
#include "openacc.h"
#endif

//==============================================================================
PyObject* K_KCORE::testerAcc(PyObject* self, PyObject* args)
{
    E_Int n;

    n = 100;
    double* c = new double [n];
    for (E_Int i = 0; i < n; i++) c[i] = 0.;

    // tester for openACC
#ifdef _OPENACC
    // get connected device type (must be 4=gpu)
    acc_device_t devtype = acc_get_device_type();
    printf("devtype=%d - %d (must be 4)\n", devtype, acc_device_nvidia); fflush(stdout);
    // force init
    acc_init(devtype);
    // Get the number of connected device
    // generally 2 : cpu + gpu
    int ndev = acc_get_num_devices(devtype);
    printf("ndev=%d\n", ndev); fflush(stdout);
#endif
    
#ifdef _OPENACC
    #pragma acc data copy(c)
    {
        #pragma acc parallel
        {
            #pragma acc loop
            for (int i = 0; i < n; i++) c[i] = 2.;
        }
    }
#endif

    // tester for openmp5

    printf("c=%f (must be 2)\n", c[0]);

    delete [] c; 
    Py_INCREF(Py_None);
    return Py_None;
}
