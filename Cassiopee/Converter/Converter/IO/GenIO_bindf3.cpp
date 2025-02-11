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

// Binary df3 (Povray density file) file support

#include "GenIO.h"
#include "Array/Array.h"
#include "CompGeom/compGeom.h"
#include <vector>
#include <stdio.h>
#include "Interp/Interp.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Write a density file from arrays
//=============================================================================
E_Int K_IO::GenIO::df3write( 
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  // Find posx, posy...
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: df3write: no coordinates in array.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Find density
  E_Int posDens = K_ARRAY::isDensityPresent(varString);
  if (posDens == -1)
  {
    printf("Warning: df3write: no density in array.\n");
    return 1;
  }
  posDens++;

  //Build position vectors for getInterpolationCell
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  vector<E_Int> posct;
  E_Int s = structField.size();
  for (E_Int n = 0; n < s; n++)
  {
    posxt.push_back(posx); posyt.push_back(posy); poszt.push_back(posz);
    posct.push_back(-1);
  }
  // Bounding box globale
  E_Float xmax, ymax, zmax, xmin, ymin, zmin;
  K_COMPGEOM::globalBoundingBox(posxt, posyt, poszt, structField, 
                                xmin, ymin, zmin, xmax, ymax, zmax);

  // Build density grid
  E_Int nig, njg, nkg;
  nig = 100; njg = 100; nkg = 100;
  FldArrayF* d = new FldArrayF(nig*njg*nkg);
  FldArrayF& density = *d;

  // Build Adt for all blocks
  vector<K_INTERP::InterpData*> adts;
  vector<void*> a1; vector<void*> a2; vector<void*> a3; vector<void*> a4;
  E_Int isBuilt;
  for (E_Int n = 0; n < s; n++)
  {
    E_Int npts = structField[n]->getSize();
    K_INTERP::InterpAdt* adt = 
      new K_INTERP::InterpAdt(npts, 
                              structField[n]->begin(posxt[n]),
                              structField[n]->begin(posyt[n]),
                              structField[n]->begin(poszt[n]),
                              &ni[n], &nj[n], &nk[n], isBuilt);
    adts.push_back(adt);
    a1.push_back(&ni[n]); a2.push_back(&nj[n]); a3.push_back(&nk[n]);
    a4.push_back(NULL);
  }
  
  // Find density value for all points of density grid
  E_Float x, y, z;
  E_Float hi = (xmax - xmin)/(nig-1);
  E_Float hj = (ymax - ymin)/(njg-1);
  E_Float hk = (zmax - zmin)/(nkg-1);
  FldArrayI indi(1);
  FldArrayF cf(8);
  FldArrayI tmpIndi(1); FldArrayF tmpCf(8);
  E_Int ind;
  E_Float voli = 0.;
  E_Int noblk = 0;
  E_Int type = 0;
  for (E_Int k = 0; k < nkg; k++)
    for (E_Int j = 0; j < njg; j++)
      for (E_Int i = 0; i < nig; i++)
      {
        //ind = i + j*nig + k*nig*njg;
        ind = k + j*nkg + i*nkg*njg;
        x = xmin + i * hi;
        y = ymin + j * hj;
        z = zmin + k * hk;
        short found = K_INTERP::getInterpolationCell(
          x, y, z, adts, structField, 
          a1, a2, a3, a4,
          posxt, posyt, poszt, posct,
          voli, indi, cf, tmpIndi, tmpCf, type, noblk);
        if (found > 0)
        {
          FldArrayF& field0 = *structField[noblk-1];
          K_INTERP::compOneInterpolatedValue(
            indi.begin(), cf, field0.begin(posDens), 
            a1[noblk-1], a2[noblk-1], a3[noblk-1],
            type, density[ind]);
        }
        else density[ind] = -K_CONST::E_MAX_FLOAT;
      }

  // Delete Adt
  for (E_Int n = 0; n < s; n++)
  {
    delete adts[n];
    if (a3[n] == NULL)
    { delete (FldArrayI*)a1[n]; }
    else {delete (E_Int*)a1[n]; delete (E_Int*)a2[n]; delete (E_Int*)a3[n]; }
  }
  a4.clear();
  // Traitement sur density
  E_Float dmax = -K_CONST::E_MAX_FLOAT;
  E_Float dmin = K_CONST::E_MAX_FLOAT;
  for (E_Int i = 0; i < density.getSize(); i++)
  {
    //density[i] = K_FUNC::E_abs(density[i]);
    dmax = K_FUNC::E_max(density[i], dmax);
    if (density[i] > -K_CONST::E_MAX_FLOAT+1.) 
      dmin = K_FUNC::E_min(density[i], dmin);
  }
  if (dmin > K_CONST::E_MAX_FLOAT-1.) dmin = 0.;
  if (K_FUNC::fEqualZero(dmin-dmax) == true) dmax = dmin + 1.;
  for (E_Int i = 0; i < density.getSize(); i++)
  {
    if (density[i] <= -K_CONST::E_MAX_FLOAT+1.) density[i] = dmin;
  }
  unsigned int* buf = new unsigned int[nig*njg*nkg];
  E_Float val;
  for (E_Int i = 0; i < density.getSize(); i++)
  {
    val = (density[i]-dmin) / (dmax - dmin) * 4294967295.;
    buf[i] = (unsigned int)val;
    //buf[i] = density[i] / dmax * 255;
  }

  // Write density grid (df3)
  FILE* ptrFile = fopen(file, "wb");
  if (ptrFile == NULL) 
  {
    printf("Warning: df3write: I can't open file %s.\n", file);
    return 1;
  }

  E_Int e = GenIO::getInstance()->machineEndianess();

  // Size of grid density
  short n[3];
  if (e == 0)
  {
    n[0] = SBE(nig); n[1] = SBE(njg); n[2] = SBE(nkg);
    for (E_Int i = 0; i < density.getSize(); i++) buf[i] = IBE(buf[i]);
  }
  else
  {
     n[0] = nig; n[1] = njg; n[2] = nkg;
  }
  fwrite(n, sizeof(short), 3, ptrFile);
  fwrite(buf, sizeof(unsigned int), nig*njg*nkg, ptrFile);

  delete [] buf;
  fclose(ptrFile);
  return 0;
}
