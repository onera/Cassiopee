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
// Formated Nasa plot3d file support

# include <stdio.h>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include "../converter.h"

using namespace K_ARRAY;
using namespace std;
using namespace K_FLD;

//=============================================================================
/* fp3dread
 Read formatted plot3d file. 
 IN: file: file name
 OUT: varString: var names in string
 OUT: ni,nj,nk: number of points
 OUT: field: field.
 This format supports only structured blocks.
 Only for grid coordinates.
*/
//=============================================================================
E_Int K_IO::GenIO::fp3dread(
  char* file, char*& varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field, vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: fp3dread: cannot open file %s.\n", file);
    return 1;
  }

  E_Int numzone = 0;
  E_Int dim = 3;
  E_Float value; E_Int valueInt; E_Int ret;

  // Lit le nombre de blocks
  ret = readInt(ptrFile, valueInt);
  if (ret == -1 || ret == 2) { fclose(ptrFile); return 1; }
  if (ret == 1) skipLine(ptrFile);
  E_Int nb = valueInt;
  if (nb <= 0) { fclose(ptrFile); return 1; }

  // Lit les dimensions des blocs
  E_Int* dims = new E_Int [3*nb];
  
  dim = 1;

  /*
  char buf[1024]; E_Int n = 0; 
  for (E_Int i = 0; i < nb; i++)
  {
    E_Int ret = readline(ptrFile, buf, 1024);
    E_Int pos = 0;
    ret = readInt(buf, 1024, pos, valueInt);
    if (ret > 0 && valueInt > 0) 
    {
      dim = 1;
      dims[n] = valueInt; n++;
      ret = readInt(buf, 1024, pos, valueInt);
      if (ret > 0 && valueInt > 0) 
      {
        dim = 2;
        dims[n] = valueInt; n++;
        ret = readInt(buf, 1024, pos, valueInt);
        if (ret > 0 && valueInt > 0)
        {dims[n] = valueInt; n++; dim=3;}
      }
    }
  }
  */

  E_Int nr = 0; E_LONG fpos; 
  fpos = KFTELL(ptrFile);
  ret = readInt(ptrFile, valueInt);
  while (ret >= 0 && valueInt > 0 && nr < 3*nb)
  { dims[nr] = valueInt; nr += 1; fpos = KFTELL(ptrFile); ret = readInt(ptrFile, valueInt);}
  KFSEEK(ptrFile, fpos, SEEK_SET);
  if (nr == 3*nb) dim = 3;
  else if (nr == 2*nb) dim = 2;
  else dim = 1;
  
  //printf("dim=" SF_D_ "\n", dim);
  //for (E_Int i = 0; i < nb*dim; i++) printf(SF_D_ "\n", dims[i]);
  
  if (dim == 1)
  {
    for (E_Int i = 0; i < nb; i++)
    {
      ni.push_back(dims[i]); nj.push_back(1); nk.push_back(1);
    }
  }
  else if (dim == 2)
  {
    for (E_Int i = 0; i < nb; i++)
    {
      ni.push_back(dims[2*i]); nj.push_back(dims[2*i+1]); nk.push_back(1);
    }
  }
  else // dim = 3
  {
    for (E_Int i = 0; i < nb; i++)
    {
      ni.push_back(dims[3*i]); 
      nj.push_back(dims[3*i+1]); 
      nk.push_back(dims[3*i+2]);
    }
  }
  delete [] dims;

  // Read zones
  while (numzone < nb)
  {
    E_Int im = ni[numzone];
    E_Int jm = nj[numzone];
    E_Int km = nk[numzone];
    //printf("ni=" SF_D3_ "\n", im,jm,km);
    FldArrayF* fp = new FldArrayF(im*jm*km, 3);
    fp->setAllValuesAtNull();
    field.push_back(fp);
    FldArrayF& f = *fp;
    E_Float* xp = f.begin(1);
    E_Float* yp = f.begin(2);
    E_Float* zp = f.begin(3);

    // Cree les noms des zones
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, numzone);
    zoneNames.push_back(zoneName);

    // Read block
    for (E_Int i = 0; i < im*jm*km; i++)
    {
      readDouble(ptrFile, value);
      xp[i] = value;
    }
    if (dim > 1)
    {
      for (E_Int i = 0; i < im*jm*km; i++)
      {
        readDouble(ptrFile, value);
        yp[i] = value;
      }
    }
    if (dim > 2)
    {
      for (E_Int i = 0; i < im*jm*km; i++)
      {
        readDouble(ptrFile, value);
        zp[i] = value;
      }
    }
    numzone++;
  }

  

  varString = new char [K_ARRAY::VARSTRINGLENGTH];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* fp3dwrite 
   IN: file: file name
   IN: varString: var names in string format
   IN: ni,nj,nk: number of points for each block
   IN: field.
*/
//=============================================================================
E_Int K_IO::GenIO::fp3dwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field,
  vector<char*>& zoneNames)
{
  E_Int nzone = field.size();
  if (nzone <= 0) return 0;

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "w");

  if (ptrFile == NULL)
  {
    printf("Warning: fp3dwrite: cannot open file %s.\n", file);
    return 1;
  }

  // Write header (number of blocks)
  E_Int nb = field.size();
  fprintf(ptrFile, SF_W5D_ "\n", nb);
  
  // Write im,jm,km for each block
  for (E_Int i = 0; i < nb; i++)
  {
    fprintf(ptrFile, SF_W5D_ " " SF_W5D_ " " SF_W5D_, ni[i], nj[i], nk[i]);
  }
  fprintf(ptrFile, "\n");

  // Build writing data format
  char format1[41], format2[82], format3[122];
  
  char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  // length of dataFmt
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format{i}, i=1,6
  sprintf(format1, "%s\n", dataFmtl);
  sprintf(format2, "%s%s\n", dataFmt, dataFmtl);
  sprintf(format3, "%s%s%s\n", dataFmt, dataFmt, dataFmtl);

  // Write zone
  for (E_Int cnt = 0; cnt < nb; cnt++)
  {
    FldArrayF& f = *field[cnt];
    E_Int nijk = ni[cnt]*nj[cnt]*nk[cnt];
    for (E_Int n = 1; n <= 3; n++)
    {
      E_Float* fn = f.begin(n);
      E_Int i = 0;
      while (i < nijk)
      {
        if (i+2 < nijk)
        {
          fprintf(ptrFile, format3, 
                  fn[i], fn[i+1], fn[i+2]);
          i += 3;
        }
        else break;
      }
      if (nijk-i == 2)
        fprintf(ptrFile, format2,
                fn[i], fn[i+1]);
      else if (nijk-i == 1)
        fprintf(ptrFile, format1, fn[i]);
    }
  }

  fclose(ptrFile);
  return 0;
}
