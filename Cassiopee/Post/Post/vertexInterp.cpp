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

# include "post.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Interpolation function for isoLine, isoSurf, isoSurfMC */
//==============================================================================
void K_POST::vertexInterp(E_Int nfld, E_Float value,
                          FldArrayF& f, E_Int poscellN,
                          E_Float f0, E_Float f1,
                          E_Int ind0, E_Int ind1,
                          FldArrayF& fiso, E_Int& npts)
{
  E_Float alpha, alpha1, val;
  E_Float df = f1-f0;
  if (K_FUNC::fEqualZero(df) == true) alpha = 1.; 
  else alpha = (value-f0)/df;
  alpha1 = 1.-alpha;
  for (E_Int j = 1; j <= nfld; j++)
  {
    val = alpha1*f(ind0, j)+alpha*f(ind1, j);
    fiso(npts, j) = val;
  }
  if (poscellN != 0)
  {
    if (f(ind0, poscellN) == 0. || f(ind1, poscellN) == 0.)
      fiso(npts, poscellN) = 0.;
  }
  npts++;
}
