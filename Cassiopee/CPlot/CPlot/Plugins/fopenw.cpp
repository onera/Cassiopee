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

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include "Data.h"

// fopen for write with check of dirname existence (one level only)
FILE* Data::fopenw(const char* path, const char* mode)
{
  E_Int end = strlen(path);
  E_Int start;
  FILE* r = NULL;
  for (start = end-1; start >= 0; start--)
  {
    if (path[start] == '/' || path[start] == '\\') break;
  }
  if (start > 0) // dirname exists
  {
    char* dirname = new char [end+1];
    for (E_Int i = 0; i < start; i++) dirname[i] = path[i];
    dirname[start] = '\0';
    
    r = fopen(dirname, "rb");
    if (r == NULL) 
    { 
      printf("Info: CPlot: directory %s not found. creating.\n", dirname); 
#ifdef _WIN32
      mkdir(dirname);
#else
      mkdir(dirname, 0755);
#endif
    }
    else fclose(r);
    delete [] dirname;
  }

  r = fopen(path, mode);
  return r;
}
