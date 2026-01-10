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
# include "String/kstring.h"

//==============================================================================
// Sorte de strcmp pour s1 non null terminated
// IN: s1: chaine (non null terminated)
// IN: s: taille de s1
// IN: s2: chaine a comparer (null terminated)
// OUT: 0: si identique, 1: sinon
//==============================================================================
#define CMPI E_Int c = 0;						\
  while (c < s && s2[c] != '\0') { if (s1[c] != s2[c]) return 1; c++; }	\
  if (s2[c] == '\0') return 0;                                          \
  return 1;

E_Int K_STRING::cmp(char* s1, E_Int s, char* s2)
{
  CMPI;
}

E_Int K_STRING::cmp(char* s1, E_Int s, const char* s2)
{
  CMPI;
}

E_Int K_STRING::cmp(const char* s1, E_Int s, const char* s2)
{
  CMPI;
}

E_Int K_STRING::cmp(const char* s1, E_Int s, char* s2)
{
  CMPI;
}

//==============================================================================
// strcmp mais plus rapide
// IN: s1: chaine (null terminated)
// IN: s2: chaine a comparer (null terminated)
// OUT: 0: si identique, 1: sinon
//==============================================================================
#define CMP char* pt1 = (char*)s1; char* pt2 = (char*)s2;       \
  while (*pt1 != '\0' && *pt2 != '\0') { \
  if (*pt1 != *pt2) return 1; \
  pt1++; pt2++; } \
  if (*pt1 == *pt2) return 0; \
  else return 1;

E_Int K_STRING::cmp(char* s1, char* s2)
{
  CMP;
}
E_Int K_STRING::cmp(char* s1, const char* s2)
{
  CMP;
}
E_Int K_STRING::cmp(const char* s1, char* s2)
{
  CMP;
}
E_Int K_STRING::cmp(const char* s1, const char* s2)
{
  CMP;
}

//==============================================================================
// Sorte de strcpy pour s1 non null terminated
// IN: s1: chaine (non null terminated)
// IN: s: taille de s1
// si end=true, add \0 to s2 (doit etre deja allouee a une taille >= s+1)
// OUT: s2: chaine null terminated ou non suivant end 
//==============================================================================
void K_STRING::cpy(char* s2, const char* s1, E_Int s, bool end)
{
  for (E_Int i = 0; i < s; i++) s2[i] = s1[i];
  if (end == true) s2[s] = '\0';
}
