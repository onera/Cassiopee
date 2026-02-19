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
#include "../Data.h"
#include <math.h>

//=============================================================================
// Look-for plugins
//=============================================================================

//=============================================================================
/*
  Look for active zone.
  Recherche la (les) zone(s) selectionnees.
  Fit le eye.
*/
//=============================================================================
void lookForActiveZone(Data* data)
{
  double xcam, ycam, zcam, xeye, yeye, zeye;
  E_Int nz;
  double xmax, ymax, zmax, xmin, ymin, zmin;
  double k, d, dx, dy, dz, dd, dbb;

  ViewInfo& view = data->_view;

  // Current camera and eye position
  xcam = view.xcam; ycam = view.ycam; zcam = view.zcam;
  xeye = view.xeye; yeye = view.yeye; zeye = view.zeye;

  // Bounding box of active zone
  nz = data->ptrState->selectedZone;
  if (nz == 0)
  {
    // global fit 
    data->initCam(); 
    data->ptrState->farClip = 1; 
    return;
  }
  
  // Calcul la BB de l'ensemble de la selection
  xmin = +1.e+20; ymin = +1.e+20; zmin = +1.e+20;
  xmax = -1.e+20; ymax = -1.e+20; zmax = -1.e+20;
  for (E_Int i = 0; i < data->_numberOfZones; i++)
  {
    Zone* z = data->_zones[i];
    if (z->selected == 1)
    {
      xmax = MAX(xmax, z->xmax);
      ymax = MAX(ymax, z->ymax);
      zmax = MAX(zmax, z->zmax);
      xmin = MIN(xmin, z->xmin);
      ymin = MIN(ymin, z->ymin);
      zmin = MIN(zmin, z->zmin);
    }
  }

  dx = xmax-xmin; dy = ymax-ymin; dz = zmax-zmin;
  d = MAX(dx, dy); d = MAX(d, dz);

  dx = xcam-xeye; dy = ycam-yeye; dz = zcam-zeye;
  dd = dx*dx+dy*dy+dz*dz;
  dd = sqrt(dd);
  dd = MAX(dd, 1.e-6);

  //dbb = data->dist2BB(
  //  data->ptrState->activePointX,  data->ptrState->activePointY,  
  //  data->ptrState->activePointZ, xmin, ymin, zmin, xmax, ymax, zmax);
  double dxl = data->ptrState->activePointX - 0.5*(xmax + xmin);
  double dyl = data->ptrState->activePointY - 0.5*(ymax + ymin);
  double dzl = data->ptrState->activePointZ - 0.5*(zmax + zmin);
  dbb = sqrt(dxl*dxl + dyl*dyl + dzl*dzl);

  // suivie des axes principaux? (si on est en mode xy, xz, yz..., c'est
  // bien de le rester)
  int suiviDir = -1;
  if (fabs(dx) < 1.e-6 && fabs(dy) < 1.e-6) suiviDir = 1; // xy
  else if (fabs(dx) < 1.e-6 && fabs(dz) < 1.e-6) suiviDir = 2; // xz
  else if (fabs(dy) < 1.e-6 && fabs(dz) < 1.e-6) suiviDir = 3; // yz
  
  // New eye except if activePoint is near block center
  // then look to active point
  if (dbb > 0.001*d)
  {
    view.xeye = 0.5*(xmax + xmin);
    view.yeye = 0.5*(ymax + ymin);
    view.zeye = 0.5*(zmax + zmin);
    data->ptrState->farClip = 1;
  }
  else
  {
    view.xeye = data->ptrState->activePointX;
    view.yeye = data->ptrState->activePointY;
    view.zeye = data->ptrState->activePointZ;
  }

  // Compute the k factor
  k = (1.7*d)/dd;
  if (k < 1.)
  {
    // Bouge la camera que si on est tres loin
    view.xcam = view.xeye + k*( xcam-xeye );
    view.ycam = view.yeye + k*( ycam-yeye );
    view.zcam = view.zeye + k*( zcam-zeye );
    data->ptrState->farClip = 1;
  }

  if (suiviDir > 0) // cam doit suivre
  {
    switch (suiviDir)
    {
      case 1:
        view.xcam = view.xeye;
        view.ycam = view.yeye;
        break;
      case 2:
        view.xcam = view.xeye;
        view.zcam = view.zeye;
        break;
      case 3:
        view.ycam = view.yeye;
        view.zcam = view.zeye;
        break;
    }
  }
  // Ajouter le choix du cote par rapport a la BB globale?
}

//=============================================================================
void lookForMaxValue(Data* data)
{
  double xcam, ycam, zcam, xeye, yeye, zeye;
  E_Int nz, nfld;
  Zone* z;
  double* f;
  E_Int npts, v, ind, n;
  E_Int imax, jmax, kmax, indmax;
  double xmax, ymax, zmax, xmin, ymin, zmin;
  double kfact, d, dx, dy, dz, dd;
  double val;

  CPlotState& state = *data->ptrState;

  // Active field
  if (state.mode >= 3) nfld = state.mode-3;
  else return;

  // Active zone
  nz = state.selectedZone;
  if (nz == 0)
  {
    val = data->_zones[0]->maxf[nfld];
    // Find min max par domaine
    for (n = 1; n < data->_numberOfZones; n++)
    {
      if (data->_zones[n]->maxf[nfld] > val)
      {
        val = data->_zones[n]->maxf[nfld];
        nz = n;
      }
    }
    state.selectedZone = nz+1;
  }
  else nz = nz-1;

  // Find Max in active zone
  z = data->_zones[nz];
  npts = z->npts;
  f = z->f[nfld];

  val = f[0];
  indmax = 0;
  imax = 0; jmax = 0; kmax = 0;

  for (ind = 0; ind < npts; ind++)
  {
    if (f[ind] > val)
    {
      val = f[ind];
      indmax = ind;
      imax = 0; jmax = 0; kmax = 0;
    }
  }

  // Set the active point
  state.activePointX = z->x[indmax];
  state.activePointY = z->y[indmax];
  state.activePointZ = z->z[indmax];
  state.activePointI = imax+1;
  state.activePointJ = jmax+1;
  state.activePointK = kmax+1;
  for (v = 0; v < z->nfield; v++)
    state.activePointF[v] = z->f[v][indmax];

  // Eventually change displayed planes
  
  // Current camera and eye position
  ViewInfo& view = data->_view;
  xcam = view.xcam;
  ycam = view.ycam;
  zcam = view.zcam;

  xeye = view.xeye;
  yeye = view.yeye;
  zeye = view.zeye;

  xmax = z->xmax;
  ymax = z->ymax;
  zmax = z->zmax;
  
  xmin = z->xmin;
  ymin = z->ymin;
  zmin = z->zmin;

  dx = xmax-xmin;
  dy = ymax-ymin;
  dz = zmax-zmin;
  
  d = MAX(dx, dy );
  d = MAX(d, dz);

  dx = xcam-xeye;
  dy = ycam-yeye;
  dz = zcam-zeye;
  dd = dx*dx+dy*dy+dz*dz;
  dd = sqrt(dd);
  dd = MAX(dd, 1.e-6);

  // New eye
  view.xeye = state.activePointX;
  view.yeye = state.activePointY;
  view.zeye = state.activePointZ;
  
  // Compute the k factor
  kfact = (1.*d)/dd;
  
  // New cam
  view.xcam = view.xeye + kfact*( xcam-xeye );
  view.ycam = view.yeye + kfact*( ycam-yeye );
  view.zcam = view.zeye + kfact*( zcam-zeye );

  data->ptrState->farClip = 1;
}

//=============================================================================
void lookForMinValue(Data* data)
{
  double xcam, ycam, zcam, xeye, yeye, zeye;
  E_Int nz, nfld;
  Zone* z;
  double* f;
  E_Int npts, v, ind, n;
  E_Int imin, jmin, kmin, indmin;
  double xmax, ymax, zmax, xmin, ymin, zmin;
  double kfact, d, dx, dy, dz, dd;
  double val;

  CPlotState& state = *data->ptrState;

  // Active field
  if (state.mode >= 3) nfld = state.mode - 3;
  else return;

  // Active zone
  nz = state.selectedZone;
  if (nz == 0)
  {
    val = data->_zones[0]->maxf[nfld];
    // Find min max par domaine
    for (n = 1; n < data->_numberOfZones; n++)
    {
      if (data->_zones[n]->maxf[nfld] > val)
      {
        val = data->_zones[n]->maxf[nfld];
        nz = n;
      }
    }
    state.selectedZone = nz+1;
  }
  else nz = nz-1;

  // Find Min in active zone
  z = data->_zones[nz];
  npts = z->npts;
  f = z->f[nfld];

  val = f[0];
  indmin = 0;
  imin = 0; jmin = 0; kmin = 0;

  for (ind = 0; ind < npts; ind++)
  {
    if (f[ind] < val)
    {
      val = f[ind];
      indmin = ind;
      imin = 0; jmin = 0; kmin = 0;
    }
  }

  // Set the active point
  state.activePointX = z->x[indmin];
  state.activePointY = z->y[indmin];
  state.activePointZ = z->z[indmin];
  state.activePointI = imin+1;
  state.activePointJ = jmin+1;
  state.activePointK = kmin+1;
  for (v = 0; v < z->nfield; v++)
    state.activePointF[v] = z->f[v][indmin];

  // Eventually change displayed planes
  
  // Current camera and eye position
  ViewInfo& view = data->_view;
  xcam = view.xcam;
  ycam = view.ycam;
  zcam = view.zcam;

  xeye = view.xeye;
  yeye = view.yeye;
  zeye = view.zeye;

  xmax = z->xmax;
  ymax = z->ymax;
  zmax = z->zmax;
  
  xmin = z->xmin;
  ymin = z->ymin;
  zmin = z->zmin;

  dx = xmax-xmin;
  dy = ymax-ymin;
  dz = zmax-zmin;
  
  d = MAX(dx, dy );
  d = MAX(d, dz);

  dx = xcam-xeye;
  dy = ycam-yeye;
  dz = zcam-zeye;
  dd = dx*dx+dy*dy+dz*dz;
  dd = sqrt(dd);
  dd = MAX(dd, 1.e-6);

  // New eye
  view.xeye = state.activePointX;
  view.yeye = state.activePointY;
  view.zeye = state.activePointZ;
  
  // Compute the k factor
  kfact = (1.*d)/dd;
  
  // New cam
  view.xcam = view.xeye + kfact*( xcam - xeye );
  view.ycam = view.yeye + kfact*( ycam - yeye );
  view.zcam = view.zeye + kfact*( zcam - zeye );

  data->ptrState->farClip = 1;
}
