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

//=============================================================================
/* 
   Find zones that are blanked by the current 
   plugin blanking function.
*/
//=============================================================================
void Data::findBlankedZones()
{
  int i, j;
  char* var;
  Zone* z;

  for (i = 0; i < _numberOfZones; i++)
  {
    z = _zones[i];
    z->blank = 0;
    if (_pref.blanking != NULL)
    {
      for (j = 0; j < _nfield; j++)
      {
        if (z->nfield > j)
        {
          var = z->varnames[j];
          if (strcmp(var, _pref.blanking->varName) == 0) z->blank = j+1;
        }
      }
    }
  }
}
