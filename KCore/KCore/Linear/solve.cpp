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

#include "linear.h"
#include "String/kstring.h"
#include <string.h>
#include <stdlib.h>

using namespace K_FLD;

extern "C"
{
  void testdbg_();
  void kgaussj_(E_Float* A, const E_Int& n, 
                E_Float* B, const E_Int& m, 
                E_Int& ok, E_Int* indxc, E_Int* indxr, E_Int* ipiv);
}

//=============================================================================
// IN: A(n,n)
// IN: B(n,m) m seconds membres
// IN: method : "gauss" only for now
// OUT: x(n,m) m solutions
// Return 0 (OK), 1 (Failed)
//=============================================================================
E_Int K_LINEAR::solve(E_Int n, E_Int m, 
                      E_Float* A, E_Float* B, E_Float* x, 
                      char* method)
{
  if (K_STRING::cmp(method, "gauss") == 0)
  {
    E_Int ok;
    // Copie A et B
    E_Float* Ac = new E_Float[n * n];
    memcpy(Ac, A, n*n*sizeof(E_Float));
    memcpy(x, B, n*m*sizeof(E_Float));
    E_Int* indxc = new E_Int[n];
    E_Int* indxr = new E_Int[n];
    E_Int* ipiv = new E_Int[n];
    kgaussj_(Ac, n, x, m, ok, indxc, indxr, ipiv);
    delete [] indxc; delete [] indxr; delete [] ipiv;
    if (ok == 1) return 0;
    else return 1;
  }
  else return 1;
  return 0;
}
