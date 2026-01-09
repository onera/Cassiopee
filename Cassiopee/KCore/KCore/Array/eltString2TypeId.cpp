/*    
    Copyright 2013-2026 ONERA.

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

//=====================================================================================
// Analyse an element string of the type "TRI_6" and return "TRI", nvpe, loc and typeId
//=====================================================================================
E_Int K_ARRAY::eltString2TypeId(char* eltString, char* eltType, E_Int& nvpe, 
  E_Int& loc, E_Int& typeId)
{
  E_Int ret = 1;
  char n[128];
  E_Int l = strlen(eltString);
  if (eltString[l-1] == '*') { loc = 1; l = l-1; }
  else loc = 0;
  typeId = 0;
  
  if (K_STRING::cmp(eltString, 5, "TETRA") == 0)
  {
    strncpy(eltType, eltString, 5); eltType[5] = '\0';
    if (l == 5) nvpe = 4;
    else
    {
      strncpy(n, eltString+6, l-6); n[l-6] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 4:
        typeId = 4; break;
      case 10:
        typeId = 14; break;
      case 16:
        typeId = 35; break;
      case 20:
        typeId = 36; break;
      default:
        ret = 0; break;
    }
  }
  else if (K_STRING::cmp(eltString, 4, "NODE") == 0)
  {
    strncpy(eltType, eltString, 4); eltType[4] = '\0';
    nvpe = 1;
    typeId = 0;
  }
  else if (K_STRING::cmp(eltString, 4, "QUAD") == 0)
  {
    if (l == 4) nvpe = 4;
    else
    {
      strncpy(eltType, eltString, 4); eltType[4] = '\0';
      strncpy(n, eltString+5, l-5); n[l-5] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 4:
        typeId = 3; break;
      case 8:
        typeId = 12; break;
      case 9:
        typeId = 13; break;
      case 12:
        typeId = 33; break;
      case 16:
        typeId = 34; break;
      case 25:
        typeId = 54; break;
      default:
        ret = 0; break;
    }
  }
  else if (K_STRING::cmp(eltString, 4, "HEXA") == 0)
  {
    if (l == 4) nvpe = 8;
    else
    {
      strncpy(eltType, eltString, 4); eltType[4] = '\0';
      strncpy(n, eltString+5, l-5); n[l-5] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 8:
        typeId = 7; break;
      case 20:
        typeId = 18; break;
      case 27:
        typeId = 19; break;
      case 32:
        typeId = 40; break;
      case 56:
        typeId = 41; break;
      case 64:
        typeId = 42; break;
      default:
        ret = 0; break;  
    }
  }
  else if (K_STRING::cmp(eltString, 4, "PYRA") == 0)
  {
    if (l == 4) nvpe = 5;
    else
    {
      strncpy(eltType, eltString, 4); eltType[4] = '\0';
      strncpy(n, eltString+5, l-5); n[l-5] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 5:
        typeId = 5; break;
      case 14:
        typeId = 15; break;
      case 13:
        typeId = 20; break;
      case 21:
        typeId = 37; break;
      case 29:
        typeId = 38; break;
      case 30:
        typeId = 39; break;
      default:
        ret = 0; break;   
    }
  }
  else if (K_STRING::cmp(eltString, 5, "PENTA") == 0)
  {
    if (l == 5) nvpe = 6;
    else
    {
      strncpy(eltType, eltString, 5); eltType[5] = '\0';
      strncpy(n, eltString+6, l-6); n[l-6] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 6:
        typeId = 6; break;
      case 15:
        typeId = 16; break;
      case 18:
        typeId = 17; break;
      default:
        ret = 0; break;
    }
  }
  else if (K_STRING::cmp(eltString, 4, "NGON") == 0)
  {
    strncpy(eltType, eltString, 4); eltType[4] = '\0';
    nvpe = 0;
    typeId = 8;
  }
  else if (K_STRING::cmp(eltString, 3, "TRI") == 0)
  {
    if (l == 3) nvpe = 3;
    else
    {
      strncpy(eltType, eltString, 3); eltType[3] = '\0';
      strncpy(n, eltString+4, l-4); n[l-4] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 3:
        typeId = 2; break;
      case 6:
        typeId = 11; break;
      case 9:
        typeId = 31; break;
      case 10:
        typeId = 32; break;
      case 12:
        typeId = 52; break;
      case 15:
        typeId = 53; break;
      default:
        ret = 0; break;  
    }
  }
  else if (K_STRING::cmp(eltString, 3, "BAR") == 0)
  {
    if (l == 3) nvpe = 2;
    else
    {
      strncpy(eltType, eltString, 3); eltType[3] = '\0';
      strncpy(n, eltString+4, l-4); n[l-4] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 2:
        typeId = 1; break;
      case 3:
        typeId = 10; break;
      case 4:
        typeId = 30; break;
      case 5:
        typeId = 51; break;
      default:
        ret = 0; break;    
    }
  }
  return ret;
}
