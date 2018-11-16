/*    
    Copyright 2013-2018 Onera.

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

//=============================================================================
/*
  Change a HO connectivity to a LO connectivity.
  Two modes: coarse or fine.
*/
//=============================================================================
E_Int K_CONNECT::connectHO2LO(const char* eltTypeHO,
                              FldArrayI& cEVHO,
                              char* eltTypeLO,
                              FldArrayI& cEVLO)
{
  
}

E_Int K_CONNECT::connectHO2LOCoarse(const char* eltTypeHO,
                                    FldArrayI& cEVHO,
                                    char* eltTypeLO,
                                    FldArrayI& cEVLO)
{
  E_Int nskip = 0;
  E_Int nkeep = 0;
  if (eltTypeHO[0] == 'B' && eltTypeHO[1] == 'A' && eltTypeHO[2] == 'R')
  {
    // BAR
    strcpy(eltType, "BAR");
    if (eltTypeHO[4] == '3')
    {
      nkeep = 2; nskip = 1;
    }
    else if (eltTypeHO[4] == '5')
    {
      nkeep = 2; nskip = 3;
    }
  }
  else if (eltTypeHO[0] == 'T' && eltTypeHO[1] == 'R' && eltTypeHO[2] == 'I')
  {
    // TRI
    strcpy(eltType, "TRI");
    if (eltTypeHO[4] == '3')
    {
      nkeep = 2; nskip = 1;
    }
    else if (eltTypeHO[4] == '5')
    {
      nkeep = 2; nskip = 3;
    }
  }

  E_Int* ptHO = cEVHO.begin();
  E_Int size = cEVHO.size();
  E_Int nelts = int(size/(nkeep+nskip));
  E_Int nsize = nelts*nkeep;
  cEVLO.malloc(nelts, nkeep);
  E_Int* ptLO = cEVLO.begin();
  for (E_Int i = 0; i < nelts; i++)
  {
    for (E_Int j = 0; j < nkeep; j++)
    { ptLO[i+j*nkeep] = ptHO[i+j*nkeep]; }
  }
}

E_Int K_CONNECT::connectHO2LOFine(const char* eltTypeHO,
                                  FldArrayI& cEVHO,
                                  char* eltTypeLO,
                                  FldArrayI& cEVLO)
{
  
}
