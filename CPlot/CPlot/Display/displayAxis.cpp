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

#include "../Data.h"

//=============================================================================
// Display the axis bottom/right
// Je ne touche pas a la matrice de projection
// Je cree juste un petit view port, j'affiche le repere a la position eye
//=============================================================================
void Data::displayAxis()
{
  // Matrice de vue
  double viewMatrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, viewMatrix);

  //glMatrixMode(GL_MODELVIEW);
  //glPushMatrix();

  glViewport(_view.w-200,0,200,200);
  glDisable(GL_DEPTH_TEST);

  // Position
  double x1 = _view.xeye;
  double y1 = _view.yeye;
  double z1 = _view.zeye;
  // vecteur view
  double vx = _view.xcam-x1;
  double vy = _view.ycam-y1;
  double vz = _view.zcam-z1;
  double dist = sqrt(vx*vx + vy*vy + vz*vz);
  dist = MAX(dist, 1.e-10);
  vx = vx/dist; vy = vy/dist; vz = vz/dist;

  // vecteur normal
  double nx, ny, nz;
  if (vx*vx > 1.e-10)
  { ny = 1.; nz = 0.; nx = vy / vx; }
  else if (vy*vy > 1.e-10)
  { nx = 1.; nz = 0.; ny = vx / vy; }
  else if (vz*vz > 1.e-10)
  { nx = 1.; ny = 0.; nz = vx / vz; }
  else { nx = 1.; ny = 0.; nz = 0.; }
  double n = 1./sqrt(nx*nx+ny*ny+nz*nz);
  nx = n*nx; ny = n*ny; nz = n*nz;

  double tx, ty, tz;
  tx = vy*nz-vz*ny; 
  ty = vz*nx-vx*nz;
  tz = vx*ny-vy*nx;
  //printf("t: %f %f %f\n", tx,ty,tz);
  //printf("n: %f %f %f\n", nx,ny,nz);

  double l = 0.3*dist;

  glColor3f(1., 0, 0);
  glLineWidth(2.);

  glBegin(GL_LINES);
  glColor3f(1., 0, 0); // X
  glVertex3d(x1, y1, z1);
  glVertex3d(x1+l, y1, z1);
  glColor3f(0., 1, 0); // Y
  glVertex3d(x1, y1, z1);
  glVertex3d(x1, y1+l, z1);
  glColor3f(0.1, 0.1, 1); // Z
  glVertex3d(x1, y1, z1);
  glVertex3d(x1, y1, z1+l);
  glEnd();

  
  E_Float offset = l*0.03; // decalage par rapport a l'axe
  E_Float shOffset = offset*0.3; // decalage pour l'ombre
  char sX[5]; strcpy(sX, "X");
  char sY[5]; strcpy(sY, "Y");
  char sZ[5]; strcpy(sZ, "Z");
  double xp, yp, zp;
  double offnx, offny, offnz, offtx, offty, offtz, offntx, offnty, offntz;
  offnx = shOffset*(-vx+nx);
  offny = shOffset*(-vy+ny);
  offnz = shOffset*(-vz+nz);
  offtx = shOffset*(-vx+tx);
  offty = shOffset*(-vy+ty);
  offtz = shOffset*(-vz+nz);
  offntx = shOffset*(-vx+tx+nx);
  offnty = shOffset*(-vy+ty+ny);
  offntz = shOffset*(-vz+tz+nz);

  // X
  vx = 0; vy = 0; vz = 0;
  xp = x1+l+4*offset; yp = y1; zp = z1;
  glColor3f(1., 0.5, 0.5);
  renderBitmapString(xp+offnx, yp+offny, zp+offnz, FONT3, sX); // +n
  renderBitmapString(xp+offtx, yp+offty, zp+offtz, FONT3, sX); // +t
  renderBitmapString(xp+offntx, yp+offnty, zp+offntz, FONT3, sX); // +t+n
  glColor3f(1., 0., 0.); // X
  renderBitmapString(xp, yp, zp, FONT3, sX);

  // Y
  xp = x1; yp = y1+l+offset; zp = z1;
  glColor3f(0.5, 1, 0.5);
  renderBitmapString(xp+offnx, yp+offny, zp+offnz, FONT3, sY); // +n
  renderBitmapString(xp+offtx, yp+offty, zp+offtz, FONT3, sY); // +t
  renderBitmapString(xp+offntx, yp+offnty, zp+offntz, FONT3, sY); // +t+n
  glColor3f(0., 1, 0); // Y
  renderBitmapString(xp, yp, zp, FONT3, sY);

  // Z
  xp = x1; yp = y1; zp = z1+l+offset;
  glColor3f(0.5, 0.5, 1);
  renderBitmapString(xp+offnx, yp+offny, zp+offnz, FONT3, sZ); // +n
  renderBitmapString(xp+offtx, yp+offty, zp+offtz, FONT3, sZ); // +t
  renderBitmapString(xp+offntx, yp+offnty, zp+offntz, FONT3, sZ); // +t+n
  glColor3f(0., 0., 1.); // Z
  renderBitmapString(xp, yp, zp, FONT3, sZ);
  
  glEnable(GL_DEPTH_TEST);
  
  glLineWidth(1.0);
  glColor3f(1., 1., 1.);
  
  //glPopMatrix();
  glViewport(0,0,_view.w,_view.h);  
}
