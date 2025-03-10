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
#include <string.h>
#include "Array/Array.h"
#include "String/kstring.h"

//=============================================================================
// Analyse a list of typeIds and loc and return eltString and a list of nvpe
//=============================================================================
E_Int K_ARRAY::typeId2eltString(const std::vector<E_Int>& typeId, E_Int loc,
  char* eltString, std::vector<E_Int>& nvpe)
{
  E_Int ret = 1;
  E_Int size = typeId.size();
  nvpe.resize(size);
  char* eltStringi = new char[28];
  
  strcpy(eltString, "");
  for (E_Int i = 0; i < size; i++)
  {
    ret = typeId2eltString(typeId[i], loc, eltStringi, nvpe[i]);
    if (ret == 0) break;
    strcat(eltString, eltStringi);
    if (i < size-1) strcat(eltString, ",");
  }
  
  delete[] eltStringi;
  return ret;
}

//=============================================================================
// Analyse a typeId and loc and return eltString and nvpe
//=============================================================================
E_Int K_ARRAY::typeId2eltString(E_Int typeId, E_Int loc, char* eltString, E_Int& nvpe)
{  
  E_Int ret = 1;
  switch (typeId)
  {
    case 0:
      strcpy(eltString, "NODE"); nvpe = 1; break;
    
    case 1:
      strcpy(eltString, "BAR"); nvpe = 2; break;
      
    case 2:
      strcpy(eltString, "TRI"); nvpe = 3; break;
    
    case 3:
      strcpy(eltString, "QUAD"); nvpe = 4; break;
    
    case 4:
      strcpy(eltString, "TETRA"); nvpe = 4; break;
    
    case 5:
      strcpy(eltString, "PYRA"); nvpe = 5; break;
    
    case 6:
      strcpy(eltString, "PENTA"); nvpe = 6; break;
  
    case 7:
      strcpy(eltString, "HEXA"); nvpe = 8; break;
      
    case 8:
      strcpy(eltString, "NGON"); nvpe = 1; break;
      
    case 10:
      strcpy(eltString, "BAR_3"); nvpe = 3; break;
      
    case 11: 
      strcpy(eltString, "TRI_6"); nvpe = 6; break;
      
    case 12:
      strcpy(eltString, "QUAD_8"); nvpe = 8; break;
      
    case 13:
      strcpy(eltString, "QUAD_9"); nvpe = 9; break;
    
    case 14:
      strcpy(eltString, "TETRA_10"); nvpe = 10; break;
    
    case 15:
      strcpy(eltString, "PYRA_14"); nvpe = 14; break;
    
    case 16:
      strcpy(eltString, "PENTA_15"); nvpe = 15; break;
    
    case 17:
      strcpy(eltString, "PENTA_18"); nvpe = 18; break;
    
    case 18:
      strcpy(eltString, "HEXA_20"); nvpe = 20; break;
    
    case 19:
      strcpy(eltString, "HEXA_27"); nvpe = 27; break;
    
    case 20:
      strcpy(eltString, "PYRA_13"); nvpe = 13; break;
    
    case 30:
      strcpy(eltString, "BAR_4"); nvpe = 4; break;
    
    case 31:
      strcpy(eltString, "TRI_9"); nvpe = 9; break;
    
    case 32:
      strcpy(eltString, "TRI_10"); nvpe = 10; break;
    
    case 33:
      strcpy(eltString, "QUAD_12"); nvpe = 12; break;
    
    case 34:
      strcpy(eltString, "QUAD_16"); nvpe = 16; break;
    
    case 35:
      strcpy(eltString, "TETRA_16"); nvpe = 16; break;
    
    case 36:
      strcpy(eltString, "TETRA_20"); nvpe = 20; break;
    
    case 37:
      strcpy(eltString, "PYRA_21"); nvpe = 21; break;
    
    case 38:
      strcpy(eltString, "PYRA_29"); nvpe = 29; break;
    
    case 39:
      strcpy(eltString, "PYRA_30"); nvpe = 30; break;
    
    case 40:
      strcpy(eltString, "HEXA_32"); nvpe = 32; break;
    
    case 41:
      strcpy(eltString, "HEXA_56"); nvpe = 56; break;
    
    case 42:
      strcpy(eltString, "HEXA_64"); nvpe = 64; break;
      
    default:
      ret = 0; break;
  }
          
  if (loc == 1) strcat(eltString, "*");
  
  return ret;
}
