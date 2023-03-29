/*    
    Copyright 2013-2023 Onera.

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
#include "../ViewInfo.h"
#include "../Zone.h"
#include <math.h>
#include <stdio.h>

#define TOL 1.e-6

//=============================================================================
// Compute frustum planes
//=============================================================================
void computeFrustumPlanes(ViewInfo& view)
{
  double dirx, diry, dirz, normi, tx, ty, tz;
  double Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz;
  double ncx, ncy, ncz;
  double xcam, ycam, zcam, Nx, Ny, Nz, Px, Py, Pz;

  xcam = view.xcam;
  ycam = view.ycam;
  zcam = view.zcam;

  // vecteur Z
  Zx = xcam - view.xeye;
  Zy = ycam - view.yeye;
  Zz = zcam - view.zeye;
  normi = 1./sqrt(Zx*Zx + Zy*Zy + Zz*Zz);
  Zx = Zx * normi;
  Zy = Zy * normi;
  Zz = Zz * normi;

  // vecteur X
  dirx = view.dirx;
  diry = view.diry;
  dirz = view.dirz;
  Xx = diry * Zz - dirz * Zy;
  Xy = dirz * Zx - dirx * Zz;
  Xz = dirx * Zy - diry * Zx;
  // dir et Z sont deja normes
  
  // Vecteur Y
  Yx = Zy * Xz - Zz * Xy;
  Yy = Zz * Xx - Zx * Xz;
  Yz = Zx * Xy - Zy * Xx;
  
  // Centres des plans
  ncx = xcam - view.nearD * Zx;
  ncy = ycam - view.nearD * Zy;
  ncz = zcam - view.nearD * Zz;
  
  //fcx = xcam - view.farD * Zx;
  //fcy = ycam - view.farD * Zy;
  //fcz = zcam - view.farD * Zz;
  
  // Near plane
  //Px = ncx; Py = ncy; Pz = ncz;
  //view.NearD = Zx*Px + Zy*Py + Zz*Pz;
  //view.NearNx = -Zx; view.NearNy = -Zy; view.NearNz = -Zz;
  
  // Far plane
  //Px = fcx; Py = fcy; Pz = fcz;
  //view.FarD = -Zx*Px - Zy*Py - Zz*Pz;
  //view.FarNx = Zx; view.FarNy = Zy; view.FarNz = Zz;

  // Top plane
  tx = ncx + Yx*view.nh - xcam;
  ty = ncy + Yy*view.nh - ycam;
  tz = ncz + Yz*view.nh - zcam;
  normi = 1./sqrt(tx*tx + ty*ty + tz*tz);
  tx = tx * normi;
  ty = ty * normi;
  tz = tz * normi;
  Nx = ty * Xz - tz * Xy;
  Ny = tz * Xx - tx * Xz;
  Nz = tx * Xy - ty * Xx;
  Px = ncx + Yx * view.nh;
  Py = ncy + Yy * view.nh;
  Pz = ncz + Yz * view.nh;
  view.TopD = -Nx*Px - Ny*Py - Nz*Pz;
  view.TopNx = Nx;
  view.TopNy = Ny;
  view.TopNz = Nz;
  
  // Bottom plane
  tx = ncx - Yx*view.nh - xcam;
  ty = ncy - Yy*view.nh - ycam;
  tz = ncz - Yz*view.nh - zcam;
  normi = 1./sqrt(tx*tx + ty*ty + tz*tz);
  tx = tx * normi;
  ty = ty * normi;
  tz = tz * normi;
  Nx = ty * Xz - tz * Xy;
  Ny = tz * Xx - tx * Xz;
  Nz = tx * Xy - ty * Xx;
  Px = ncx - Yx * view.nh;
  Py = ncy - Yy * view.nh;
  Pz = ncz - Yz * view.nh;
  view.BottomD = Nx*Px + Ny*Py + Nz*Pz;
  view.BottomNx = -Nx;
  view.BottomNy = -Ny;
  view.BottomNz = -Nz;

  // Left plane
  tx = ncx - Xx*view.nw - xcam;
  ty = ncy - Xy*view.nw - ycam;
  tz = ncz - Xz*view.nw - zcam;
  normi = 1./sqrt(tx*tx + ty*ty + tz*tz);
  tx = tx * normi;
  ty = ty * normi;
  tz = tz * normi;
  Nx = ty * Yz - tz * Yy;
  Ny = tz * Yx - tx * Yz;
  Nz = tx * Yy - ty * Yx;
  Px = ncx - Xx * view.nw;
  Py = ncy - Xy * view.nw;
  Pz = ncz - Xz * view.nw;
  view.LeftD = -Nx*Px - Ny*Py - Nz*Pz;
  view.LeftNx = Nx;
  view.LeftNy = Ny;
  view.LeftNz = Nz;
  
  // Right plane
  tx = ncx + Xx*view.nw - xcam;
  ty = ncy + Xy*view.nw - ycam;
  tz = ncz + Xz*view.nw - zcam;
  normi = 1./sqrt(tx*tx + ty*ty + tz*tz);
  tx = tx * normi;
  ty = ty * normi;
  tz = tz * normi;
  Nx = ty * Yz - tz * Yy;
  Ny = tz * Yx - tx * Yz;
  Nz = tx * Yy - ty * Yx;
  Px = ncx + Xx * view.nw;
  Py = ncy + Xy * view.nw;
  Pz = ncz + Xz * view.nw;
  view.RightD = Nx*Px + Ny*Py + Nz*Pz;
  view.RightNx = -Nx;
  view.RightNy = -Ny;
  view.RightNz = -Nz;
}

//=============================================================================
// Return 1 if zone is in frustum
// 0 otherwise
//=============================================================================
E_Int isInFrustum(Zone* z, ViewInfo& view)
{
#ifdef __MESA__
  // to avoid abusive clipping in osmesa
  return 1;
#endif

  E_Int out;
  double dist1, dist2, dist3, dist4;

  // Bounding box of zone
  double xmin = z->xmin;
  double xmax = z->xmax;
  double ymin = z->ymin;
  double ymax = z->ymax;
  double zmin = z->zmin;
  double zmax = z->zmax;
  double bbx[8]; double bby[8]; double bbz[8];
  bbx[0] = xmin; bby[0] = ymin; bbz[0] = zmin;
  bbx[1] = xmax; bby[1] = ymin; bbz[1] = zmin;
  bbx[2] = xmax; bby[2] = ymax; bbz[2] = zmin;
  bbx[3] = xmin; bby[3] = ymax; bbz[3] = zmin;
  bbx[4] = xmin; bby[4] = ymin; bbz[4] = zmax;
  bbx[5] = xmax; bby[5] = ymin; bbz[5] = zmax;
  bbx[6] = xmax; bby[6] = ymax; bbz[6] = zmax;
  bbx[7] = xmin; bby[7] = ymax; bbz[7] = zmax;

  double leftD = view.LeftD;
  double leftNx = view.LeftNx;
  double leftNy = view.LeftNy;
  double leftNz = view.LeftNz;
  double rightD = view.RightD;
  double rightNx = view.RightNx;
  double rightNy = view.RightNy;
  double rightNz = view.RightNz;
  double topD = view.TopD;
  double topNx = view.TopNx;
  double topNy = view.TopNy;
  double topNz = view.TopNz;
  double bottomD = view.BottomD;
  double bottomNx = view.BottomNx;
  double bottomNy = view.BottomNy;
  double bottomNz = view.BottomNz;

  // all points in
  for (E_Int k = 0; k < 8; k++)
  {
    dist1 = leftNx*bbx[k] + leftNy*bby[k] + leftNz*bbz[k] + leftD;
    dist2 = rightNx*bbx[k] + rightNy*bby[k] + rightNz*bbz[k] + rightD;
    dist3 = topNx*bbx[k] + topNy*bby[k] + topNz*bbz[k] + topD;
    dist4 = bottomNx*bbx[k] + bottomNy*bby[k] + bottomNz*bbz[k] + bottomD;
    if (dist1 >= -TOL && dist2 >= -TOL && dist3 >= -TOL && dist4 >= -TOL) return 1;
  }

  // all points left of left plane
  out = 0;
  for (E_Int k = 0; k < 8; k++)
  {
    dist1 = leftNx*bbx[k] + leftNy*bby[k] + leftNz*bbz[k] + leftD;
    if (dist1 < -TOL) out++;
  }
  if (out == 8) return 0;
  
  // all points roght of righ plane
  out = 0;
  for (E_Int k = 0; k < 8; k++)
  {
    dist2 = rightNx*bbx[k] + rightNy*bby[k] + rightNz*bbz[k] + rightD;
    if (dist2 < -TOL) out++;
  }
  if (out == 8) return 0;
  
  // all points top of top plane
  out = 0;
  for (E_Int k = 0; k < 8; k++)
  {
    dist3 = topNx*bbx[k] + topNy*bby[k] + topNz*bbz[k] + topD;
    if (dist3 < -TOL) out++;
  }
  if (out == 8) return 0;
  
  // all points bottom of bottom plane
  out = 0;
  for (E_Int k = 0; k < 8; k++)
  {
    dist4 = bottomNx*bbx[k] + bottomNy*bby[k] + bottomNz*bbz[k] + bottomD;
    if (dist4 < -TOL) out++;
  }
  if (out == 8) return 0;

  return 1;
}
