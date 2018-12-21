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

// Cree des elts discrets a partir de definition CAO

# include "GenIO.h"
# include <vector>

using namespace K_FLD;
using namespace std;

//=============================================================================
void K_IO::GenIO::createPoint(vector<FldArrayF*>& unstructField,
                              vector<FldArrayI*>& connect,
                              vector<E_Int>& eltType,
                              E_Float density,
                              E_Float x, E_Float y, E_Float z)
{
  FldArrayF* f = new FldArrayF(1, 3);
  FldArrayI* cn = new FldArrayI();

  FldArrayF& fp = *f;
  fp(0,1) = x; fp(0,2) = y; fp(0,3) = z;

  unstructField.push_back(f);
  connect.push_back(cn);
  eltType.push_back(0);
}

//=============================================================================
void K_IO::GenIO::createLine(
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk, 
  E_Float density,
  E_Float x1, E_Float y1, E_Float z1,
  E_Float x2, E_Float y2, E_Float z2)
{
  E_Float len = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
  E_Int n = E_Int(len*density)+2;
  
  FldArrayF* f = new FldArrayF(n, 3);

  E_Float* xp = f->begin(1);
  E_Float* yp = f->begin(2);
  E_Float* zp = f->begin(3);

  E_Float delta = 1./n;
  for (E_Int i = 0; i < n; i++)
  {
    xp[i] = x1 + i*delta;
    yp[i] = y1 + i*delta;
    zp[i] = z1 + i*delta;
  }
  structField.push_back(f);
  ni.push_back(n); nj.push_back(1); nk.push_back(1);
}
