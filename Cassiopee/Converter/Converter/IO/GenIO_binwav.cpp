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

// Binary wav (RIFF) file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   Write array in the wav (audio) form.
   - Only first structured or unstructured array is taken into account.
   - Vars are taken to be:
   Time: time, t. If not found, first variable.
   Field: Pressure, Density, second variable if any.
*/
//=============================================================================
E_Int K_IO::GenIO::wavwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = structField.size();
  
  // Look for field
  FldArrayF* a = NULL;
  if (nzone >= 1)
  {
    a = structField[0];
  }
  else 
  {
    nzone = unstructField.size();
    if (nzone > 1) a = unstructField[0];
    else return 1; // no valid field
  }

  // Look for variables
  E_Int time = K_ARRAY::isTimePresent(varString);
  if (time == -1) time = 0;
  E_Int field = K_ARRAY::isPressurePresent(varString);
  if (field == -1) field = K_ARRAY::isDensityPresent(varString);
  if (field == -1) field = 1;
  if (field+1 > a->getNfld()) return 1; // not enough field
  time = time+1;
  field = field+1;

  //printf("time detected as variable %d\n", time);
  //printf("field detected as variable %d\n", field);

  // Calcul de min Deltat
  E_Int N = a->getSize();
  FldArrayF& f = *a;
  E_Float Deltat = K_CONST::E_MAX_FLOAT;
  for (E_Int i = 0; i < N-1; i++)
  {
    Deltat = K_FUNC::E_min( Deltat, f(i+1, time)-f(i, time) );
  }
  Deltat = K_FUNC::E_max(Deltat, 1.e-12);
  //printf("Time step is %f\n", Deltat);

  // Calcul de fmin, fmax
  E_Float dmin = K_CONST::E_MAX_FLOAT;
  E_Float dmax = -K_CONST::E_MAX_FLOAT;
  for (E_Int i = 0; i < N-1; i++)
  {
    dmin = K_FUNC::E_min( dmin, f(i, field) );
    dmax = K_FUNC::E_max( dmax, f(i, field) );
  }

  // Echantillonnage
  E_Int s = E_Int( (f(N-1, time) - f(0, time))/Deltat );
  unsigned char* sig = new unsigned char[s];
  E_Int l = 1;
  E_Float d, t, alpha;
  
  for (E_Int i = 0; i < s; i++)
  {
    t = f(0, time) + Deltat*i;
    while (f(l, time) < t) l++;
    alpha = (t - f(l-1, time))/(f(l, time)-f(l-1, time));
    d = (1. - alpha) * f(l-1, field) + alpha * f(l, field);

    sig[i] = char((dmin + (d - dmin)/(dmax - dmin))*255);
  }

  // Check machine endianess
  E_Int endianess = machineEndianess();

  // Ouverture fichier
  FILE* ptrFile = fopen(file, "wb");
  if (ptrFile == NULL) 
  {
    printf("Warning: wavwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Ecriture du fichier
  char fmt[5];
  strcpy(fmt, "RIFF");
  fwrite(fmt, sizeof(char), 4, ptrFile);
  int size = 4 + (8 + 16) + (8 + s);
  if (endianess == 1) size = IBE(size);
  fwrite(&size, sizeof(int), 1, ptrFile);
  strcpy(fmt, "WAVE");
  fwrite(fmt, sizeof(char), 4, ptrFile);
  
  // Fmt chunk
  strcpy(fmt, "fmt ");
  fwrite(fmt, sizeof(char), 4, ptrFile);
  size = 16;
  if (endianess == 1) size = IBE(size);
  fwrite(&size, sizeof(int), 1, ptrFile);
  short sc = 1; // format tag (no compression = PCM)
  if (endianess == 1) sc = SBE(sc);
  fwrite(&sc, sizeof(short), 1, ptrFile);
  sc = 1; // number of channels
  if (endianess == 1) sc = SBE(sc);
  fwrite(&sc, sizeof(short), 1, ptrFile);
  int smpl = int(1./Deltat); // sample rate
  if (endianess == 1) smpl = IBE(smpl);
  fwrite(&smpl, sizeof(int), 1, ptrFile);
  smpl = int(1./Deltat); // bytes per second
  //printf("bytes per second %d %d\n", smpl, endianess);
  if (endianess == 1) smpl = IBE(smpl);
  fwrite(&smpl, sizeof(int), 1, ptrFile);
  sc = 1; // block align
  if (endianess == 1) sc = SBE(sc);
  fwrite(&sc, sizeof(short), 1, ptrFile);
  sc = 8; // bits per sample
  if (endianess == 1) sc = SBE(sc);
  fwrite(&sc, sizeof(short), 1, ptrFile);

  // Data chunk
  strcpy(fmt, "data");
  fwrite(fmt, sizeof(char), 4, ptrFile);
  size = s;
  if (endianess == 1) size = IBE(size);
  fwrite(&size, sizeof(int), 1, ptrFile);
  if (endianess == 1)
  {
    for (E_Int i = 0; i < s; i++) sig[i] = SBE(sig[i]);
  }
  fwrite(sig, sizeof(char), s, ptrFile);

  delete [] sig;
  fclose(ptrFile);
  return 0;
}
