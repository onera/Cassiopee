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
// Formatted Onera v3d file support

# include <stdio.h>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "../converter.h"

using namespace K_ARRAY;
using namespace std;
using namespace K_FLD;

//=============================================================================
/* fv3dread
 Read formatted v3d file. 
 IN: file: file name
 OUT: varString: var names in string
 OUT: ni,nj,nk: number of points
 OUT: field: field.
 Only supports structured blocks.
*/
//=============================================================================
E_Int K_IO::GenIO::fv3dread(
  char* file, char*& varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field, vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: fv3dread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [K_ARRAY::VARSTRINGLENGTH];

  char t[256]; E_Int ret; E_Int block;
  char var[256]; // variable courante

  E_Int numzone = 0;

  // Read nfld
  E_Float value;
  E_Int valueInt;
  ret = readInt(ptrFile, valueInt);
  if (ret == -1 || ret == 2) { fclose(ptrFile); return 1; }
  E_Int nfld = valueInt;
  if (nfld <= 0) { fclose(ptrFile); return 1; }

  // Read varstring
  block = 0;
  fscanf(ptrFile, "%s", t);

  block:
  if (strcmp(t, "va") == 0) fscanf(ptrFile, "%s", t);
  if (block == 0) strcpy(varString, t);
  strcpy(var, t);

  // Read format
  fscanf(ptrFile, "%s", t);
  E_Int formatLength = getFormatLength(t);
  
  // Read im, jm, km
  ret = readInt(ptrFile, valueInt); // no var
  if (strcmp(var, "x") != 0 && strcmp(var, "CoordinateX") != 0 &&
      strcmp(var, "y") != 0 && strcmp(var, "CoordinateY") != 0 &&
      strcmp(var, "z") != 0 && strcmp(var, "CoordinateZ") != 0)
  {
    ret = readInt(ptrFile, valueInt);
  }
  ret = readInt(ptrFile, valueInt);
  E_Int im = valueInt;
  ret = readInt(ptrFile, valueInt);
  E_Int jm = valueInt;
  ret = readInt(ptrFile, valueInt);
  E_Int km = valueInt;
  ni.push_back(im); nj.push_back(jm); nk.push_back(km);
  if (ret == 1) // garbage
  { int c = fgetc(ptrFile);
    while (c != '\n' && c != '\r') c = fgetc(ptrFile);
  }
  
  FldArrayF* fp = new FldArrayF(im*jm*km, nfld);
  field.push_back(fp);
  FldArrayF& f = *fp;
  E_Float* xp = f.begin(1);

  // Cree les noms des zones
  char* zoneName = new char [128];
  sprintf(zoneName, "Zone" SF_D_, numzone); numzone++;
  zoneNames.push_back(zoneName);

  // Read block
  for (E_Int i = 0; i < im*jm*km; i++)
  {
    readDouble(ptrFile, value, formatLength);
    xp[i] = value;
  }

  // Read other vars
  for (E_Int nvar = 2; nvar <= nfld; nvar++)
  {
    fscanf(ptrFile, "%s", t);
    if (strcmp(t, "va") == 0) fscanf(ptrFile, "%s", t);

    if (block == 0)
    {
      strcat(varString, ",");
      strcat(varString, t);
    }
    fscanf(ptrFile, "%s", t);
    ret = readDouble(ptrFile, value);
    ret = readDouble(ptrFile, value);
    ret = readDouble(ptrFile, value);
    ret = readDouble(ptrFile, value);
    if (ret == 1) // garbage
    { int c = fgetc(ptrFile);
      while (c != '\n' && c != '\r') c = fgetc(ptrFile);
    }

    // Read block
    xp = f.begin(nvar);
    for (E_Int i = 0; i < im*jm*km; i++)
    {
      readDouble(ptrFile, value, formatLength);
      xp[i] = value;
    }
  }

  ret = fscanf(ptrFile, "%s", t); 
  if (ret != EOF)
  {
    block++;
    goto block;
  }

  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* fv3dwrite 
   IN: file: file name
   IN: varString: var names in string format
   IN: ni,nj,nk: number of points for each block
   IN: field.
*/
//=============================================================================
E_Int K_IO::GenIO::fv3dwrite(
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
    printf("Warning: fv3dwrite: cannot open file %s.\n", file);
    return 1;
  }

  // Write header (number of variables)
  fprintf(ptrFile, SF_W5D_ "\n", field[0]->getNfld());
  vector<char*> vars;
  extractVars(varString, vars);
  
  // Replace CoordinateX, CoordinateY, CoordinateZ with x,y,z
  for (size_t i = 0; i < vars.size(); i++)
  {
    if (strcmp(vars[i], "CoordinateX") == 0) strcpy(vars[i], "x");
    if (strcmp(vars[i], "CoordinateY") == 0) strcpy(vars[i], "y");
    if (strcmp(vars[i], "CoordinateZ") == 0) strcpy(vars[i], "z");
  }

  // Build writing data format
  char format1[30], format2[60], format3[90], format4[120], 
    format5[150], format6[180];

  char dataFmtl[29];
  strcpy(dataFmtl, dataFmt);
  // length of dataFmt
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format{i}, i=1,6
  sprintf(format1, "%s\n", dataFmtl);
  sprintf(format2, "%s%s\n", dataFmt, dataFmtl);
  sprintf(format3, "%s%s%s\n", dataFmt, dataFmt, dataFmtl);
  sprintf(format4, "%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmtl);
  sprintf(format5, "%s%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmt, 
          dataFmtl);
  sprintf(format6, "%s%s%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmt, 
          dataFmt, dataFmtl);

  // Checking format and getting informations for header
  // dataFmtl = %[width].[precision][specifier]
  // specifier is the only non-optional parameter
  //char flags[20];
  char* width = new char[20];
  char* precision = new char[20];
  char* specifier = new char[20];
  // length of dataFmtl
  int len = strlen(dataFmtl); 
  // specifier: last character
  specifier[0] = dataFmtl[len-1];
  specifier[1] ='\0';
  int ipch = 0;
  while (dataFmtl[ipch] != '.' && ipch < len) ipch++;

  if (ipch != len)
  {
    for (E_Int i = 1; i < ipch; i++)
    {
      width[i-1] = dataFmtl[i];
    }
    width[ipch-1] = '\0';
    if (dataFmtl[len-2] == '.')
    {
      // default value
      precision[0] ='9'; precision[1] = '\0';
    }
    else
    {
      for (E_Int i = 0; i < len-ipch-2; i++)
      {
        precision[i] = dataFmtl[ipch+1+i];
      }
      precision[len-ipch-2] = '\0';
    }
  }
  else
  {
    precision[0] = '9';
    precision[1] = '\0';
    // if last character of dataFmtl is '%', then width has default value
    // otherwise, width has value between dataFmtl[1] and dataFmtl[l-1]
    if (dataFmtl[len-2] == '%')
    {
      // default value
      //width[0] = '1';
      //width[1] = '6';
      width[0] = '\0';
    }
    else
    {
      for (E_Int i = 1; i < len-1; i++)
      {
        width[i-1] = dataFmtl[i];
      }
      width[len-2] = '\0';
    }
  }

  // Write zone
  for (E_Int cnt = 0; cnt < nzone; cnt++)
  {
    FldArrayF& f = *field[cnt];
    E_Int nijk = ni[cnt]*nj[cnt]*nk[cnt];
  
    for (E_Int n = 1; n <= f.getNfld(); n++)
    {
      char* varsn1 = vars[n-1];
      E_Float* fn = f.begin(n);
      if (strcmp(varsn1, "x") == 0 || 
          strcmp(varsn1, "CoordinateX") == 0 ||
          strcmp(varsn1, "y") == 0 || 
          strcmp(varsn1, "CoordinateY") == 0 ||
          strcmp(varsn1, "z") == 0 || 
          strcmp(varsn1, "CoordinateZ") == 0)
      {
        fprintf(ptrFile, "%s                   6%s%s.%s\n", 
                varsn1, specifier, width, precision);
        fprintf(ptrFile, SF_W6D_ SF_W6D_ SF_W6D_ SF_W6D_ "\n",
                E_Int(1), ni[cnt], nj[cnt], nk[cnt]);
      }
      else 
      {
        fprintf(ptrFile, "%s                   6%s%s.%s\n", 
                varsn1, specifier, width, precision);
        fprintf(ptrFile, SF_W6D_ SF_W6D_ SF_W6D_ SF_W6D_ SF_W6D_ "\n",
                n, E_Int(1), ni[cnt], nj[cnt], nk[cnt]);
      }
    
      E_Int i = 0;
      while (i < nijk)
      {
        if (i+5 < nijk)
        {
          fprintf(ptrFile, format6, 
                  fn[i], fn[i+1], fn[i+2],
                  fn[i+3], fn[i+4], fn[i+5]);
          i += 6;
        }
        else break;
      }
      if (nijk-i == 5)
        fprintf(ptrFile, format5,
                fn[i], fn[i+1], fn[i+2],
                fn[i+3], fn[i+4]);
      else if (nijk-i == 4)
        fprintf(ptrFile, format4, 
                fn[i], fn[i+1], fn[i+2], fn[i+3]);
      else if (nijk-i == 3)
      {
        fprintf(ptrFile, format3, 
                fn[i], fn[i+1], fn[i+2]);
      }
      else if (nijk-i == 2)
        fprintf(ptrFile, format2, fn[i], fn[i+1]);
      else if (nijk-i == 1) fprintf(ptrFile, format1, fn[i]);
    }
  }

  // delete arrays
  delete [] width;
  delete [] precision;
  delete [] specifier;
  for (size_t i = 0; i < vars.size(); i++) delete [] vars[i];

  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Return the length of each character if formate string is:
   formatString = 6eX.Y.
   Return X
*/
//=============================================================================
E_Int K_IO::GenIO::getFormatLength(char* formatString)
{
  char number[256];
  E_Int is = -1;
  E_Int ie = -1;
  E_Int c = 0;
  E_Int lenString = strlen(formatString);
  for (c = 0; c < lenString; c++)
  {
    if (formatString[c] == 'e') is = c+1;
    if (formatString[c] == '.' && is != -1) ie = c;
  }
  if (is != -1 && ie != -1 && ie-is > 0)
  {
    for (c = is; c < ie; c++) number[c-is] = formatString[c];
    number[ie-is] = '\0';
    E_Int value = (E_Int)strtod(number, NULL);
    return value;
  }
  else return -1;
}
