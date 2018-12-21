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

#include "Zone1D.h"

//=============================================================================
Zone1D::Zone1D(char* varString, K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn)
{
  _np = f.getSize();
  _nv = f.getNfld();
  _ne = cn.getSize();
  if (cn.getNfld() != 2) printf("Warning: invalid array type for 1D display.\n");
  // copie
  E_Float* fp = f.begin();
  _f = new double [_np*_nv];
  for (E_Int i = 0; i < _np*_nv; i++) _f[i] = fp[i];
  E_Int* cnp = cn.begin();
  _cn = new int [_ne*2];
  for (E_Int i = 0; i < _ne*2; i++) _cn[i] = cnp[i];

  std::vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  int varSize = vars.size();
  if (varSize != _nv) printf("Warning: incoherent array varstring.\n");

  // Allocation of var fields
  _varNames = new char* [_nv];
  for (E_Int n = 0; n < _nv; n++)
  {
    _varNames[n] = new char [strlen(vars[n])+1];
    strcpy(_varNames[n], vars[n]);
    delete [] vars[n];
  }
}

//=============================================================================
Zone1D::~Zone1D()
{
  delete [] _f;
  delete [] _cn;
  for (E_Int i = 0; i < _nv; i++) delete [] _varNames[i];
  delete [] _varNames;
}
