/*    
    Copyright 2013-2019 Onera.

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

# include "connector.h"
#include "stub.h"

//=============================================================================
/* Calcule et stocke les coefficients d'interpolation par intersection
   OUT: [donorBlks,donorInd1D, donorType, coefs, extrapInd1D, orphanInd1D] 
        donorBlks: no du blk donneur, démarre à 0
        donorInd1D: indice global (structure), de l elt (NS) du donneur
        donorType: type d interpolation effectué localement
        coefs: coefficients d interpolation, stockés selon le type
        extrapInd1D: indices des pts extrapolés
        orphanInd1D: indices des pts orphelins */
//=============================================================================
PyObject* K_CONNECTOR::setInterpDataCons(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}
