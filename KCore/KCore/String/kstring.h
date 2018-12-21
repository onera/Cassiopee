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

#ifndef _KCORE_STRING_H_
#define _KCORE_STRING_H_
#include "Def/DefTypes.h"

namespace K_STRING
{
  /* strcmp avec une chaine non terminee par \0
     IN: s1: chaine non terminee par \0
     IN: s: taille de s1
     IN: s2: chaine a comparer (terminee par \0)
     OUT: 0: si identique, 1: sinon */
  E_Int cmp(char* s1, E_Int s, char* s2);
  E_Int cmp(char* s1, E_Int s, const char* s2);
  /* strcmp(plus rapide)
     IN: s1: chaine terminee par \0
     IN: s2: chaine a comparer (terminee par \0)
     OUT: 0: si identique, 1: sinon */
  E_Int cmp(char* s1, char* s2);
  E_Int cmp(char* s1, const char* s2);
  E_Int cmp(const char* s1, char* s2);
  E_Int cmp(const char* s1, const char* s2);
  /* cpy pour s1 non null terminated
     IN: s1: chaine (non null terminated)
     IN: s: taille de s1
     OUT: s2: chaine null terminated (deja allouee a une taille >= s+1) */
  void cpy(char* s2, char* s1, E_Int s);
}

#endif
