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
#include "../Data.h"

//============================================================================
/*
  addPressureVariable.
  This function is a add-a-variable plugin function.
  Add pressure variable to every blocks if possible.
*/
//============================================================================
void addPressureVariable()
{
  int i, nf, ni, nj, nk, found, ind;
  double ro, rou, rov, row, roE, roe;
  int pro, prou, prov, prow, proE;
  char varName[] = "pressure";
  double gamma = 1.4;

  for (i = 0; i < _numberOfZones; i++)
  {
    ni = _szones[i]->ni;
    nj = _szones[i]->nj;
    nk = _szones[i]->nk;

    // Does this variable already exists?
    found = 0;
    for (nf = 0; nf < _szones[i]->nfield; nf++)
    {
      if (strcmp( varName, _szones[i]->varnames[nf] ) == 0 )
      {
        found = 1;
        break;
      }
    }

    if (found == 0) // if not, create it
    {
      nf = _szones[i]->nfield;
      _szones[i]->f[nf-3] = (double*)malloc( sizeof(double)* ni * nj * nk );
      strcpy(_szones[i]->varnames[nf], varName);
      _szones[i]->nfield++;
    }
    
    // Check if depending variables exists
    pro = checkVariable(i, "ro");
    prou = checkVariable(i, "rou");
    prov = checkVariable(i, "rov");
    prow = checkVariable(i, "row");
    proE = checkVariable(i, "roE");

    if (pro == -1 || prou == -1 || prov == -1 || prow == -1 || proE == -1)
    {
      setVariableToZero(ni, nj, nk, _szones[i]->f[nf-3]);
    }
    else
    {
      // Compute variable stored at nf
      for (ind = 0; ind < ni*nj*nk; ind++)
      {
        ro = _szones[i]->f[pro-3][ind];
        rou = _szones[i]->f[prou-3][ind];
        rov = _szones[i]->f[prov-3][ind];
        row = _szones[i]->f[prow-3][ind];
        roE = _szones[i]->f[proE-3][ind];
        
        roe = roE - 0.5*(rou*rou +  rov*rov + row*row)/ro;
        _szones[i]->f[nf-3][ind] = (gamma-1)*roe;
      }
    }
    
    findFMinMax(_szones[i]);
    
  }
  globFMinMax(_szones, _numberOfStructZones);
  printTmpMessage("Variable created.");
}
