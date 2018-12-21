/*
    Copyright 2013-2019 Onera.

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

#include "Zone.h"
#include "Functions.h"
#include <math.h>
#include <stdio.h>

//=============================================================================
// Find min and max of zone coordinates.
// Store it in zone->xmin...
//=============================================================================
void findMinMax(Zone* zone)
{
  int npts = zone->npts;
  double xminl = 1.e6;
  double xmaxl = -1.e6;
  double yminl = 1.e6;
  double ymaxl = -1.e6;
  double zminl = 1.e6;
  double zmaxl = -1.e6;

  double* x = zone->x;
  double* y = zone->y;
  double* z = zone->z;

//#pragma omp parallel default(shared)
  {
    int i;
    double xminloc = 1.e6;
    double xmaxloc = -1.e6;
    double yminloc = 1.e6;
    double ymaxloc = -1.e6;
    double zminloc = 1.e6;
    double zmaxloc = -1.e6;
    double xx, yy, zz;

//#pragma omp for nowait
    for (i = 0; i < npts; i++)
    {
      xx = x[i];
      xminloc = MIN(xminloc, xx);
      xmaxloc = MAX(xmaxloc, xx);
      yy = y[i];
      yminloc = MIN(yminloc, yy);
      ymaxloc = MAX(ymaxloc, yy);
      zz = z[i];
      zminloc = MIN(zminloc, zz);
      zmaxloc = MAX(zmaxloc, zz);
    }
//#pragma omp critical
    {
      xminl = MIN(xminl, xminloc);
      xmaxl = MAX(xmaxl, xmaxloc);
      yminl = MIN(yminl, yminloc);
      ymaxl = MAX(ymaxl, ymaxloc);
      zminl = MIN(zminl, zminloc);
      zmaxl = MAX(zmaxl, zmaxloc);
    }
  }

  zone->xmin = xminl;
  zone->xmax = xmaxl;
  zone->ymin = yminl;
  zone->ymax = ymaxl;
  zone->zmin = zminl;
  zone->zmax = zmaxl;

  // Adjust dim
  double dz = zone->zmax - zone->zmin;
  if (dz < 1.e-12 && zone->dim == 3) zone->dim = 2;

  // To avoid pbs
  double dx = zone->xmax - zone->xmin;
  if (dx < 1.e-12) { zone->xmin += -1.e-4; zone->xmax += 1.e-4; }
  double dy = zone->ymax - zone->ymin;
  if (dy < 1.e-12) { zone->ymin += -1.e-4; zone->ymax += 1.e-4; }

  // Compute Di, Dj, Dk
  double np = pow(npts, 0.333);
  double l = (dx + dy + dz)*0.333;
  if (dz < 1.e-6) { l = (dx+dy+1.e-6)*0.5; np = pow(npts, 0.5); }
  if (dy < 1.e-6) { l = (dx+dz+1.e-6)*0.5; np = pow(npts, 0.5); }
  if (dx < 1.e-6) { l = (dy+dz+1.e-6)*0.5; np = pow(npts, 0.5); }

  zone->Di = np/l; // Devrait etre la longueur de l'edge max en i
  zone->Dj = np/l;
  zone->Dk = np/l;
}

//============================================================================
/*
  Find min and max of zone fields
  Store it in zone->minf...
*/
//============================================================================
void findFMinMax(Zone* zone)
{
  int npts = zone->npts;
  int nf = zone->nfield;

/* avec openmp
  double* fmin = new double[nf];
  double* fmax = new double[nf];
  for (int j = 0; j < nf; j++) { fmin[j] = 1.e6; fmax[j] = -1.e6; }

#pragma omp parallel default(shared)
  {
    int i, j;
    double fi; double* f;
    double* fminloc = new double[nf];
    double* fmaxloc = new double[nf];

    for (j = 0; j < nf; j++)
    {
      fminloc[j] = 1.e6; fmaxloc[j] = -1.e6;
      f = zone->f[j];

#pragma omp for nowait
      for (i = 0; i < npts; i++)
      {
        fi = f[i];
        fminloc[j] = MIN(fminloc[j], fi);
        fmaxloc[j] = MAX(fmaxloc[j], fi);
      }
    }
#pragma omp critical
    {
      for (j = 0; j < nf; j++)
      {
        fmin[j] = MIN(fmin[j], fminloc[j]);
        fmax[j] = MAX(fmax[j], fmaxloc[j]);
      }
    }
    delete [] fminloc; delete [] fmaxloc;
  }

  for (int j = 0; j < nf; j++)
  {
    zone->minf[j] = fmin[j]; zone->maxf[j] = fmax[j];
    //printf("zone %d: %f %f\n", j, fmin[j], fmax[j]);
  }
  delete [] fmin; delete [] fmax;
*/

  // sans openmp
  int i, j;
  double* f; double fi;
  for (j = 0; j < nf; j++)
  {
    f = zone->f[j];
    zone->minf[j] = 1.e6; zone->maxf[j] = -1.e6;
    for (i = 0; i < npts; i++)
    {
      fi = f[i];
      zone->minf[j] = MIN(zone->minf[j], fi);
      zone->maxf[j] = MAX(zone->maxf[j], fi);
    }
  }

}

//=============================================================================
// Compute the global min max coord of all zones
// The min max of each zone must have been computed before.
// IN: zones: All Zones
// IN: nz: number of zones
//=============================================================================
void globMinMax(Zone** zones, int nz,
                double& xmin, double& xmax,
                double& ymin, double& ymax,
                double& zmin, double& zmax,
                double& epsup, double& epsstrafe, double& dmoy)
{
  int i;
  Zone* z;
  xmin = 1.e6; ymin = 1.e6; zmin = 1.e6;
  xmax = -1.e6; ymax = -1.e6; zmax = -1.e6;

  for (i = 0; i < nz; i++)
  {
    z = zones[i];
    xmin = MIN(xmin, z->xmin);
    ymin = MIN(ymin, z->ymin);
    zmin = MIN(zmin, z->zmin);
    xmax = MAX(xmax, z->xmax);
    ymax = MAX(ymax, z->ymax);
    zmax = MAX(zmax, z->zmax);
  }

  dmoy = (xmax-xmin)+(ymax-ymin)+(zmax-zmin);
  dmoy = dmoy * 0.3;
  epsup = dmoy * 0.002;
  epsstrafe = dmoy*0.001;

  //printf("coord max: %f %f %f\n", xmax, ymax, zmax);
  //printf("coord min: %f %f %f\n", xmin, ymin, zmin);
}

//============================================================================
/*
  Compute the global min max of fields for all zones
  The min max of each zone must have been computed before
*/
//============================================================================
void globFMinMax(Zone** zones, int nz,
                 double* minf, double* maxf)
{
  int i, j, nf;

  if (nz > 0) nf = zones[0]->nfield;
  else nf = 0;

  for (j = 0; j < nf; j++)
  {
    maxf[j] = -1.e6; minf[j] = 1.e6;
    for (i = 0; i < nz; i++)
    {
      Zone* z = zones[i];
      if (j < z->nfield)
      {
        maxf[j] = MAX(maxf[j], z->maxf[j]);
        minf[j] = MIN(minf[j], z->minf[j]);
      }
    }
    //printf("glob:%d: %f %f\n", j, minf[j], maxf[j]);
  }
}
