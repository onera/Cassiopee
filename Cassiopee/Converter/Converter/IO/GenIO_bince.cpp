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

// Binary tecplot routines with endian conversion

# include <string.h>
# include "GenIO.h"
# include "Array/Array.h"

using namespace K_FLD;

//=============================================================================
/* 
   Read file header (endian conversion).
   This routine is compatible from v75 to v112.
   Retourne 1 (success), 0 (failed).
*/
//=============================================================================
E_Int K_IO::GenIO::readHeaderCE(FILE *ptrFile, E_Int& nfield, 
                               char*& varString, char* version)
{
  E_Int i, j;
  int ib;
  
  /* Constants */
  E_Int si = sizeof(int);

  /* Version */
  char c[9];
  fread(c, sizeof(char), 8, ptrFile); // version
  c[8] = '\0';
  strcpy(version, c);
  E_Int nversion = numeralVersion(version);
  if (nversion == -1) return 0;

  /* Endian check 1 */
  fread(&ib, si, 1, ptrFile);

  /* FileType: ajoute depuis la version 112 */
  if (nversion >= 112) fread(&ib, si, 1, ptrFile);

  /* Title */
  i = 0; ib = 1;
  while (ib != 0)
  {
    fread(&ib, si, 1, ptrFile);
    i++;
  }
  
  /* Number of variables */
  fread(&ib, si, 1, ptrFile);
  nfield = IBE(ib);
  varString = new char [nfield*K_ARRAY::VARNAMELENGTH];

  /* Variables name */
  varString[0] = '\0';
  j = 0; i = 0;
  while (j < nfield)
  {
    ib = 1;
    while (ib != 0)
    {
      fread(&ib, si, 1, ptrFile);
      if (ib != 0)
        varString[i] = char(IBE(ib));
      else
        varString[i] = ','; 
      i++;
    }
    j++;
  }
  varString[i-1] = '\0';
  return 1;
}
