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

#ifndef _KCORE_STRING_H_
#define _KCORE_STRING_H_
#include "Def/DefTypes.h"

// C format specifiers:
//   - int and long int: d, ld (or lld on Windows)
//   - simple and double precision float: f, lf
#if defined E_DOUBLEINT
#if defined _WIN32
#define SF_D_ "%lld"
#define SF_D2_ "%lld %lld"
#define SF_D3_ "%lld %lld %lld"
#define SF_D4_ "%lld %lld %lld %lld"
#define SF_D5_ "%lld %lld %lld %lld %lld"
#define SF_D6_ "%lld %lld %lld %lld %lld %lld"
#define SF_D7_ "%lld %lld %lld %lld %lld %lld %lld"
#define SF_D8_ "%lld %lld %lld %lld %lld %lld %lld %lld"
#define SF_D9_ "%lld %lld %lld %lld %lld %lld %lld %lld %lld"
#define SF_W5D_ "%5lld"
#define SF_W6D_ "%6lld"
#else
#define SF_D_ "%ld"
#define SF_D2_ "%ld %ld"
#define SF_D3_ "%ld %ld %ld"
#define SF_D4_ "%ld %ld %ld %ld"
#define SF_D5_ "%ld %ld %ld %ld %ld"
#define SF_D6_ "%ld %ld %ld %ld %ld %ld"
#define SF_D7_ "%ld %ld %ld %ld %ld %ld %ld"
#define SF_D8_ "%ld %ld %ld %ld %ld %ld %ld %ld"
#define SF_D9_ "%ld %ld %ld %ld %ld %ld %ld %ld %ld"
#define SF_W5D_ "%5ld"
#define SF_W6D_ "%6ld"
#endif
#else
#define SF_D_ "%d"
#define SF_D2_ "%d %d"
#define SF_D3_ "%d %d %d"
#define SF_D4_ "%d %d %d %d"
#define SF_D5_ "%d %d %d %d %d"
#define SF_D6_ "%d %d %d %d %d %d"
#define SF_D7_ "%d %d %d %d %d %d %d"
#define SF_D8_ "%d %d %d %d %d %d %d %d"
#define SF_D9_ "%d %d %d %d %d %d %d %d %d"
#define SF_W5D_ "%5d"
#define SF_W6D_ "%6d"
#endif

#if defined E_DOUBLEREAL
#define SF_F_ "%lf"
#define SF_F2_ "%lf %lf"
#define SF_F3_ "%lf %lf %lf"
#define SF_F4_ "%lf %lf %lf %lf"
#define SF_F5_ "%lf %lf %lf %lf %lf"
#define SF_F6_ "%lf %lf %lf %lf %lf %lf"
#define SF_F7_ "%lf %lf %lf %lf %lf %lf %lf"
#define SF_F8_ "%lf %lf %lf %lf %lf %lf %lf %lf"
#else
#define SF_F_ "%f"
#define SF_F2_ "%f %f"
#define SF_F3_ "%f %f %f"
#define SF_F4_ "%f %f %f %f"
#define SF_F5_ "%f %f %f %f %f"
#define SF_F6_ "%f %f %f %f %f %f"
#define SF_F7_ "%f %f %f %f %f %f %f"
#define SF_F8_ "%f %f %f %f %f %f %f %f"
#endif

namespace K_STRING
{
  /* strcmp avec une chaine non terminee par \0
     IN: s1: chaine non terminee par \0
     IN: s: taille de s1
     IN: s2: chaine a comparer (terminee par \0)
     OUT: 0: si identique, 1: sinon */
  E_Int cmp(char* s1, E_Int s, char* s2);
  E_Int cmp(char* s1, E_Int s, const char* s2);
  /* strcmp (plus rapide)
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
     IN: end: si end=true, ajoute null terminated a s2 (doit etre deja alloue a une taille >= s+1)
     OUT: s2: chaine copiee de s1 */
  void cpy(char* s2, const char* s1, E_Int s, bool end=true);
}

#endif
