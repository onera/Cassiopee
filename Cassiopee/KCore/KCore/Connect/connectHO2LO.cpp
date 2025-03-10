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
#include "Connect/connect.h"
#include "String/kstring.h"
using namespace K_FLD;

E_Int connectHO2LOCoarse(const char* eltTypeHO,
                         FldArrayI& cEVHO,
                         FldArrayI& cEVLO);
E_Int connectHO2LOFine(const char* eltTypeHO,
                       FldArrayI& cEVHO,
                       FldArrayI& cEVLO);

//=============================================================================
/*
  Change a HO connectivity to a LO connectivity.
  Two modes: coarse or fine.
  mode=0: return only the main points of element: TRI_6 becomes TRI
  mode=1: descretise high order elements: return TRI or TETRA. 
*/
//=============================================================================
E_Int K_CONNECT::connectHO2LO(const char* eltTypeHO,
                              FldArrayI& cEVHO,
                              FldArrayI& cEVLO, 
                              E_Int mode)
{
  if (mode == 0) return connectHO2LOCoarse(eltTypeHO, cEVHO, cEVLO);
  else if (mode == 1) return connectHO2LOFine(eltTypeHO, cEVHO, cEVLO);
  return 0;
}

// Converti une connectivite HO par elements en connectivite LO coarse
// cEVLO doit etre deja alloue
E_Int connectHO2LOCoarse(const char* eltTypeHO,
                         FldArrayI& cEVHO,
                         FldArrayI& cEVLO)
{
  // Calcul des strides
  
  E_Int nelts = cEVHO.getSize();
  E_Int nvpe = cEVLO.getNfld();
  // keep first
#pragma omp parallel for
  for (E_Int i = 0; i < nelts; i++) 
    for (E_Int n = 1; n <= nvpe; n++) cEVLO(i,n) = cEVHO(i,n);

  return 1;
}

// Retourne toujours une connectivite TRI ou TETRA
E_Int connectHO2LOFine(const char* eltTypeHO,
                       FldArrayI& cEVHO,
                       FldArrayI& cEVLO)
{
  E_Int nelts = cEVHO.getSize();
  
  // BAR
  if (K_STRING::cmp((char*)eltTypeHO, 3, "BAR_3") == 0)
  {
    // nptsF = npts + nelts
    // neltsF = 2*nelts 
#pragma omp parallel for
    for (E_Int i = 0; i < nelts; i++) 
    {
      cEVLO(2*i,1) = cEVHO(i,1);
      cEVLO(2*i,2) = cEVHO(i,3);
      cEVLO(2*i+1,1) = cEVHO(i,3);
      cEVLO(2*i+1,2) = cEVHO(i,2);
    }
  }
  return 1;
}
