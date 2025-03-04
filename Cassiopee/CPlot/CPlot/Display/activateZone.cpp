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
#include "../Data.h"

//=============================================================================
// These functions are called when activating and deactivating a zone
// This triggers a graphical effect defined in those functions.
//=============================================================================
void Data::activateZone()
{
  E_Int i, ind;
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);
  Zone* z = _zones[ptrState->selectedZone-1];
  E_Int npts = z->npts;
  double dx, dy, dz; // deltas
  double xc, yc, zc;

  // Center of zone
  xc = z->xc; yc = z->yc; zc = z->zc;

  // Temporary storage of point coordinates
  E_Int size = sizeof(double) * npts;
  double *vx = (double*)malloc(size);
  double *vy = (double*)malloc(size);
  double *vz = (double*)malloc(size);
  memcpy(vx, z->x, size);
  memcpy(vy, z->y, size);
  memcpy(vz, z->z, size);

  for (i = 0; i < npts; i++)
  {
    z->x[i] = xc; z->y[i] = yc; z->z[i] = zc;
  }
  
  for (i = 0; i < nroll; i++)
  { 
    for (ind = 0; ind < npts; ind++)
    {
      dx = vx[ind] - xc;
      dy = vy[ind] - yc;
      dz = vz[ind] - zc;
      
      z->x[ind] += dx*froll;
      z->y[ind] += dy*froll;
      z->z[ind] += dz*froll;
    }

    displayBB();
    display();
  }

  // Recuperation des coordonnees reelles
  memcpy(z->x, vx, size);
  memcpy(z->y, vy, size);
  memcpy(z->z, vz, size);

  free(vx); free(vy); free(vz);
}

//=============================================================================
void Data::deactivateZone()
{
  E_Int i, ind;
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);
  Zone* z = _zones[ptrState->selectedZone-1];
  E_Int npts = z->npts;
  double dx, dy, dz; // deltas
  double xc, yc, zc;

  // Center of zone
  xc = z->xc;
  yc = z->yc;
  zc = z->zc;

  // Temporary storage of point coordinates
  E_Int size = sizeof(double) * npts;
  double *vx = (double*)malloc(size);
  double *vy = (double*)malloc(size);
  double *vz = (double*)malloc(size);
  memcpy(vx, z->x, size);
  memcpy(vy, z->y, size);
  memcpy(vz, z->z, size);

  for (i = 0; i < nroll; i++)
  { 
    for (ind = 0; ind < npts; ind++)
    {
      dx = vx[ind] - xc;
      dy = vy[ind] - yc;
      dz = vz[ind] - zc;
      
      z->x[ind] = z->x[ind] - dx*froll;
      z->y[ind] = z->y[ind] - dy*froll;
      z->z[ind] = z->z[ind] - dz*froll;
    }

    displayBB();
    display();
  }

  // Recuperation des coordonnees reelles
  memcpy(z->x, vx, size);
  memcpy(z->y, vy, size);
  memcpy(z->z, vz, size);

  free(vx); free(vy); free(vz);
}
