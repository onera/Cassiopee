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

// Endian conversion / detection

#include "GenIO.h"

using namespace K_FLD;

//=============================================================================
/* Return endianess of the current machine:
   1: big endian, 0: little endian. */
//=============================================================================
E_Int K_IO::GenIO::machineEndianess()
{
  int t = 1;
  int r = t<<1;
  if (r == 2) return 0;
  else return 1;
}
//=============================================================================
/* Change endianess of a 2 bytes value */
//=============================================================================
char* K_IO::GenIO::conv2(char* x)
{ 
  static int l = 0;
  char* cl = (char*)&l;
  cl[0] = x[1];
  cl[1] = x[0];
  return (char*)&l;
}

//=============================================================================
/* Change endianess of a 4 bytes value */
//=============================================================================
char* K_IO::GenIO::conv4(char* x)
{ 
  static int l = 0;
  char* cl = (char*)&l;
  cl[0] = x[3];
  cl[1] = x[2];
  cl[2] = x[1];
  cl[3] = x[0];
  return (char*)&l;
}

//=============================================================================
/* Change endianess of a 8 bytes value */
//=============================================================================
char* K_IO::GenIO::conv8(char* x)
{ 
  static E_LONG l = 0;
  char* cl = (char*)&l;
  cl[0] = x[7];
  cl[1] = x[6];
  cl[2] = x[5];
  cl[3] = x[4];
  cl[4] = x[3];
  cl[5] = x[2];
  cl[6] = x[1];
  cl[7] = x[0];
  return (char*)&l;
}

//=============================================================================
/* Convert endians in a float field */
//=============================================================================
void K_IO::GenIO::convertEndianField(FldArrayF& field)
{
  E_Int np = field.getSize();
  E_Int nf = field.getNfld();

  for (E_Int n = 0; n < nf; n++)
  {
    E_Float* f = field.begin(n+1);
    for (E_Int i = 0; i < np; i++)
    {
#ifdef E_DOUBLEREAL
      f[i] = DBE(f[i]);
#else
      f[i] = FBE(f[i]);
#endif
    }
  }
}
//=============================================================================
/* Convert endians in an integer field */
//=============================================================================
void K_IO::GenIO::convertEndianField(FldArrayI& field)
{
  E_Int np = field.getSize();
  E_Int nf = field.getNfld();

  for (E_Int n = 0; n < nf; n++)
  {
    E_Int* f = field.begin(n+1);
    for (E_Int i = 0; i < np; i++)
    {
#ifdef E_DOUBLEINT
      f[i] = LBE(f[i]);
#else
      f[i] = IBE(f[i]);
#endif
    }
  }
}
