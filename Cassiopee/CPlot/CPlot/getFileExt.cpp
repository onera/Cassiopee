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

#include "cplot.h"

//=============================================================================
// Retourne l'extension du fichier
// Il faut deja avoir allouer ext (de taille 5)
// si il n'y a pas d'extension dans file ou si l'extension est trop 
// superieure a 5 characteres, retourne ext = "None"
// sinon retourne l'extension
//=============================================================================
void getFileExt(char* file, char* ext)
{
  E_Int l = strlen(file);
  char* buf = new char[l];
  E_Int c = 0;
  for (E_Int i = l-1; i >=0; i--)
  {
    if (file[i] == '.') break;
    else
    {
      buf[c] = file[i]; c++;
    }
  }
  if (c < l && c < 5)
  {
    for (E_Int i = 0; i < c; i++)
      ext[i] = buf[c-1-i];
    ext[c]='\0';
  }
  else
    strcpy(ext, "None");
  delete [] buf;
  //printf("ext = %s\n", ext);
}
