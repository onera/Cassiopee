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

#include "GenIO.h"
#include <stdlib.h>
#include <string.h>

#include "Def/DefFunction.h"

using namespace K_FUNC;
using namespace K_FLD;

K_IO::GenIO* K_IO::GenIO::_instance = NULL;

//=============================================================================
K_IO::GenIO* K_IO::GenIO::getInstance()
{
  if (_instance == NULL) _instance = new K_IO::GenIO;
  return _instance;
}

//=============================================================================
K_IO::GenIO::GenIO()
{
  _convertEndian = false;
  _intLength = 8;
  _realLength = 8;
}

//=============================================================================
K_IO::GenIO::~GenIO()
{
}

//=============================================================================
/*
  Check if the endian in file is the same as the endian of the current
  machine for binary tecplot files or file too short to read 9 bytes.
  Return -1 if it fails.
  Return 0 if the endian are the same.
  Return 1 otherwise.
*/
//=============================================================================
E_Int K_IO::GenIO::tecCheckEndian(char* file)
{
  int ib, ret;
  FILE* ptrFile;
  
  /* Open file */
  ptrFile = fopen(file, "rb");
  if (ptrFile == NULL) return -1;
  
  /* Version */
  char c[9];
  ret = fread(c, sizeof(char), 8, ptrFile);
  if (ret < 8) return -1;

  /* Cette valeur doit etre 1 */
  ret = fread(&ib, sizeof(int), 1, ptrFile);
  if (ret < 1) return -1;

  fclose(ptrFile);
  if (ib == 1) return 0;
  else return 1;
}

//=============================================================================
/*
  Check if the endian in file is the same as the endian of the current
  machine for v3d binary files.
  Return -1 if file doesnt exist.
  Return -2 if format is not recognised.
  Return 0 if the endian are the same and integers are i4.
  Return 1 if the endian are the same and integers are i8.
  Return 2 if the endian are different and integers are i4.
  Return 3 if the endian are different and integers are i8.
*/
//=============================================================================
E_Int K_IO::GenIO::v3dCheckEndian(char* file)
{
  FILE* ptrFile;
  int ib, rb, ret;
  E_Int si = sizeof(int);

  /* Open file */
  ptrFile = fopen(file, "rb");
  if (ptrFile == NULL) return -1;

  // In ifort binary file, first record is the size of next
  // record coded on 4 bytes.
  ret = fread(&ib, si, 1, ptrFile);
  if (ret < 1) return -1;

  fclose(ptrFile);

  if (ib < 100)
  {
    _convertEndian = false;
    _intLength = ib;
    if (ib == 4) return 0;
    else if (ib == 8) return 1;
    else return -2;
  }
  else
  {
    rb = IBE(ib);
    _convertEndian = true;
    _intLength = rb;
    if (rb == 4) return 2;
    else if (rb == 8) return 3;
    else return -2;
  }
}

//=============================================================================
E_Bool K_IO::GenIO::getConvertEndian()
{
  return _convertEndian;
}

//=============================================================================
/* Retourne le no de version a partir d'une chaine style : "#!TDV108"
   Retourne -1 si la chaine ne peut pas etre interpretee correctement.
*/
//=============================================================================
E_Int K_IO::GenIO::numeralVersion(char* version)
{
  char number[4];
  if (strlen(version) < 8) return -1;
  if (version[0] != '#' || version[1] != '!' || version[2] != 'T' ||
      version[3] != 'D' || version[4] != 'V') return -1;
  number[0] = version[5]; number[1] = version[6]; number[2] = version[7];
  number[3] = '\0';
  E_Int ret = atoi(number);
  return ret;
}

