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
#include <math.h>

//=============================================================================
// Retourne les deux caracteres affiches en 1D a partir d'un varName
//=============================================================================
void Data::getCharsFromVarName(char* varName, char& c1, char& c2)
{
  E_Int l = strlen(varName);
  if (l == 0) { c1 = '\0'; c2 = '\0' ; return; }
  if (l == 1) { c1 = varName[0]; c2 = '\0'; return; }
  char lastChar = varName[l-1];

  if (strcmp(varName, "CoordinateX") == 0)
  {
    c1 = 'X'; c2= '\0';
  }
  else if (strcmp(varName, "CoordinateY") == 0)
  {
    c1 = 'Y'; c2= '\0';
  }
  else if (strcmp(varName, "CoordinateZ") == 0)
  {
    c1 = 'Z'; c2= '\0';
  }
  else if (lastChar == 'X' || lastChar == 'x')
  {
    c1 = varName[0]; c2 = lastChar;
  }
  else if (lastChar == 'Y' || lastChar == 'y')
  {
    c1 = varName[0]; c2 = lastChar;
  }
  else if (lastChar == 'Z' || lastChar == 'z')
  {
    c1 = varName[0]; c2 = lastChar;
  }
  else
  {
    c1 = varName[0];
    if (l > 1) c2 = varName[1];
    else c2 = '\0';
  }
}

//=============================================================================
// display plots
//=============================================================================
void Data::displayPlots()
{
  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix(); glLoadIdentity();
  
  unsigned int ns = _slots1D.size();
  for (unsigned int i = 0; i < ns; i++) displayPlot(_slots1D[i]);

  // Put back the previous projection
  glPopMatrix();
  resetPerspectiveProjection(); 
}

//=============================================================================
// display plot
// - display background
// - display axis
// - display title
// - display plot (for each 1D zone)
//=============================================================================
void Data::displayPlot(Slot1D* s)
{
  // Coord ecran du background
  E_Float delta = 10.;
  E_Float dx = ((_view.w-delta)*1./ptrState->gridSizeI);
  E_Float dy = ((_view.h-delta)*1./ptrState->gridSizeJ);
  E_Float posx = s->_gridPosI*dx+delta*0.5;
  E_Float posy = s->_gridPosJ*dy+delta*0.5;

  // display 1D background (carre blanc)
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_DEPTH_TEST);

  glColor4f(1., 1., 1., s->_bgBlend);
  glBegin(GL_QUADS);
  glVertex3d(posx, posy, 0);
  glVertex3d(posx+dx, posy, 0);
  glVertex3d(posx+dx, posy+dy, 0);
  glVertex3d(posx, posy+dy, 0);
  glEnd();

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  // display 1D axis
  plot1DAxis(s, posx, posy, dx, dy, s->_bgBlend);

  // display Title

  // display plot
  for (unsigned int i = 0; i < s->_zones.size(); i++)
    plotZone(s, s->_zones[i], posx, posy, dx, dy, s->_var1[i], s->_var2[i]);

}

//=============================================================================
void Data::plotZone(Slot1D* s, Zone1D* z, E_Float posx, E_Float posy,
                    E_Float dx, E_Float dy, int var1, int var2)
{
  E_Int ne = z->_ne;
  int c1, c2;
  int* c = z->_cn;
  double* f1 = &(z->_f[var1*(z->_np)]);
  double* f2 = &(z->_f[var2*(z->_np)]);
  // scale
  double borderx = 0.1*dx;
  double bordery = 0.1*dy;
  double r1min = s->_r1min;
  double r1max = s->_r1max;
  double r2min = s->_r2min;
  double r2max = s->_r2max;
  if (fabs(r1min-r1max) < 1.e-10) { r1min = r1min-1.; r1max = r1max+1.; }
  if (fabs(r2min-r2max) < 1.e-10) { r2min = r2min-1.; r2max = r2max+1.; }
  double a1 = (dx-2*borderx)/(r1max-r1min);
  double b1 = (posx+borderx)-a1*r1min;
  double a2 = -(dy-2*bordery)/(r2max-r2min);
  double b2 = (posy+bordery)-a2*r2max;
 
  glDisable(GL_DEPTH_TEST);
  glColor4f(1., 0.11, 0.11, 1.); // red foreground color
  //glColor4f(249./255., 243./255., 52./255., 1.); // yellow foreground color
  glLineWidth(3.);
  glBegin(GL_LINES);
  for (E_Int i = 0; i < ne; i++)
  {
    c1 = c[i]-1; c2 = c[ne+i]-1;
    //printf("%f %f\n", a1*f1[c1]+b1, a2*f2[c1]+b2);
    glVertex3d(a1*f1[c1]+b1, a2*f2[c1]+b2, 0.);
    glVertex3d(a1*f1[c2]+b1, a2*f2[c2]+b2, 0.);
  }
  glEnd();

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glColor4f(1., 1., 1., 0.25); // White shadow en dessous
  glLineWidth(3.);
  glBegin(GL_LINES);
  for (E_Int i = 0; i < ne; i++)
  {
    c1 = c[i]-1; c2 = c[ne+i]-1;
    //printf("%f %f\n", a1*f1[c1]+b1, a2*f2[c1]+b2);
    glVertex3d(a1*f1[c1]+b1, a2*f2[c1]+b2+2, 0.);
    glVertex3d(a1*f1[c2]+b1, a2*f2[c2]+b2+2, 0.);
  }
  glEnd();
  
  /*
  glColor4f(0., 0., 0., 0.25); // Black shadow au dessus
  glBegin(GL_LINES);
  glLineWidth(2.);
  for (E_Int i = 0; i < ne; i++)
  {
    c1 = c[i]-1; c2 = c[ne+i]-1;
    //printf("%f %f\n", a1*f1[c1]+b1, a2*f2[c1]+b2);
    glVertex3d(a1*f1[c1]+b1, a2*f2[c1]+b2-2, 0.);
    glVertex3d(a1*f1[c2]+b1, a2*f2[c2]+b2-2, 0.);
  }
  glEnd();
  */
  glDisable(GL_BLEND);

  glLineWidth(1.);

  glEnable(GL_DEPTH_TEST);
  glColor3f(1., 1., 1.);

  // display active point on 1D plot
  display1DActivePoint(s, z, posx, posy, dx, dy, var1, var2);
}

//=============================================================================
// plot axis
// Essaie d'afficher un nombre raisonable de ticks, correspondant
// a des valeurs remarquables
//=============================================================================
void Data::plot1DAxis(Slot1D* s, E_Float posx, E_Float posy,
                      E_Float dx, E_Float dy, E_Float blend)
{
  E_Float x, y;
  
  float fc, bc; // foreground color, background color
  if (blend < 0.3 && ptrState->bgColor == 0)
  { fc = 0.3; bc = 1.; }
  else {fc = 1.; bc = 0.3; }

  char figure[20];
  double r1min = s->_r1min;
  double r1max = s->_r1max;
  double r2min = s->_r2min;
  double r2max = s->_r2max;
  if (fabs(r1min-r1max) < 1.e-10) { r1min = r1min-1.; r1max = r1max+1.; }
  if (fabs(r2min-r2max) < 1.e-10) { r2min = r2min-1.; r2max = r2max+1.; }

  char c11 = s->_var1NameC1[0];
  char c12 = s->_var1NameC2[0];
  char c21 = s->_var2NameC1[0];
  char c22 = s->_var2NameC2[0];

  double borderx = 0.1*dx;
  double bordery = 0.1*dy;
  //double deltax = dx-2*borderx;
  //double deltay = dy-2*bordery;

  E_Float delta1 = r1max-r1min;
  E_Float delta2 = r2max-r2min;
  double a1 = (dx-2*borderx)/delta1;
  double b1 = (posx+borderx)-a1*r1min;
  double a2 = -(dy-2*bordery)/delta2;
  double b2 = (posy+bordery)-a2*r2max;

  //E_Float fontShift = 0.9;

  // Guess rangeTick
  E_Float tx = getTick(r1min, r1max);
  E_Float ty = getTick(r2min, r2max);
  //printf("tx %f %f\n", tx, ty);

  glDisable(GL_DEPTH_TEST);
  glLineWidth(3.);
  glColor4f(fc, fc, fc, 1.); // foreground color
  glBegin(GL_LINES);
  // x ticks
  glVertex3d(posx+borderx, posy+dy-bordery,0.);
  glVertex3d(posx+dx-borderx, posy+dy-bordery,0.);
  E_Int imin = E_Int(r1min/tx);
  E_Int imax = E_Int(r1max/tx);
  for (E_Int i = imin; i < imax; i++)
  {
    x = a1*i*tx+b1;
    //printf("x=%f\n", i*tx);
    glVertex3d(x, posy+dy-bordery, 0.);
    glVertex3d(x, posy+dy-bordery-5, 0.);
  }
  //x = a1*imin*tx+b1;
  //sprintf(figure, "%f", imin*tx);
  //printf("figure=%s\n", figure);
  //renderBitmapString(x, posy+dy-bordery, 0., FONT1, figure);

  // y ticks
  glVertex3d(posx+borderx, posy+dy-bordery,0.);
  glVertex3d(posx+borderx, posy+bordery,0.);
  imin = E_Int(r2min/ty);
  imax = E_Int(r2max/ty);
  for (E_Int i = imin; i < imax; i++)
  {
    y = a2*i*ty+b2;
    glVertex3d(posx+borderx, y, 0.);
    glVertex3d(posx+borderx+5, y, 0.);
  }
  glEnd();

  // shadows
  glColor4f(bc, bc, bc, 1.); // shadow color
  glBegin(GL_LINES);
  // x ticks
  glVertex3d(posx+borderx, posy+dy-bordery+2,0.);
  glVertex3d(posx+dx-borderx, posy+dy-bordery+2,0.);
  imin = E_Int(r1min/tx);
  imax = E_Int(r1max/tx);
  for (E_Int i = imin; i < imax; i++)
  {
    x = a1*i*tx+b1;
    //printf("x=%f\n", i*tx);
    glVertex3d(x, posy+dy-bordery+2., 0.);
    glVertex3d(x, posy+dy-bordery-5+2., 0.);
  }
  //x = a1*imin*tx+b1;
  //sprintf(figure, "%f", imin*tx);
  //printf("figure=%s\n", figure);
  //renderBitmapString(x, posy+dy-bordery, 0., FONT1, figure);

  // y ticks
  glVertex3d(posx+borderx+2., posy+dy-bordery,0.);
  glVertex3d(posx+borderx+2., posy+bordery,0.);
  imin = E_Int(r2min/ty);
  imax = E_Int(r2max/ty);
  for (E_Int i = imin; i < imax; i++)
  {
    y = a2*i*ty+b2;
    glVertex3d(posx+borderx+2., y, 0.);
    glVertex3d(posx+borderx+5+2., y, 0.);
  }
  glEnd();

  glLineWidth(1.);

  // x figures
  E_Int ly = FONTSIZE1;
  E_Int lx;
  imin = E_Int(r1min/tx);
  imax = E_Int(r1max/tx);
  for (E_Int i = imin; i < imax; i+=5)
  {
    x = a1*i*tx+b1;
    sprintf(figure, "%g", i*tx);
    lx = textWidth(FONT1, figure);
    renderStringWithShadow(x-lx*0.5, posy+dy-bordery+2+ly, 0., FONT1, figure,
                           fc, fc, fc, 1.,
                           bc, bc, bc, 1.);
  }

  // x var name
  figure[0] = c11; figure[1] = c12; figure[2] = '\0';
  lx = textWidth(FONT1, figure);
  renderStringWithShadow(posx+dx-borderx+lx, posy+dy-bordery+2+ly, 0., 
                         FONT1, figure, fc, fc, fc, 1.,
                         bc, bc, bc, 1.);

  // y figures
  imin = E_Int(r2min/ty);
  imax = E_Int(r2max/ty);
  for (E_Int i = imin; i < imax; i+=5)
  {
    y = a2*i*ty+b2;
    sprintf(figure, "%g", i*ty);
    lx = textWidth(FONT1, figure);
    renderStringWithShadow(posx+borderx-lx-3, y+ly*0.5, 0., FONT1, figure,
                           fc, fc, fc, 1.,
                           bc, bc, bc, 1.);
  }
  
  // y var name
  figure[0] = c21; figure[1] = c22; figure[2] = '\0';
  lx = textWidth(FONT1, figure);
  renderStringWithShadow(posx+borderx-3*lx, posy+bordery-ly, 0., FONT1, figure,
                         fc, fc, fc, 1.,
                         bc, bc, bc, 1.);
  glEnable(GL_DEPTH_TEST);
  glColor4f(1., 1., 1., 1.);
}

//=============================================================================
E_Float Data::getTick(E_Float rmin, E_Float rmax)
{
  E_Float delta = (rmax-rmin)/10.;
  E_Float where = log10(delta);
  E_Int wi = (E_Int)floor(where);
  //printf("delta=%f where=%f; wi=%d\n", delta, where, wi);
  E_Float tick = pow(10., wi);
  return tick;
} 

//=============================================================================
// Determine les indices du point actif dans la courbe
// Retourne 1 (trouve) et 0 (FAIL)
//=============================================================================
int Data::getActivePointIndex(Zone1D* z, int var1, int var2,
                              int& e1, int& e2, double& alpha)
{
  // active point
  int nz = ptrState->selectedZone;
  if (nz == 0) return 0; // FAIL
 
  double xP = ptrState->activePointX;
  double yP = ptrState->activePointY;
  double zP = ptrState->activePointZ;
  double* fv1 = &(z->_f[var1*(z->_np)]);
  double* fv2 = &(z->_f[var2*(z->_np)]);

  // Zone 1D
  
  // Check coordinates
  // si var1 ou var2 est une coordonnee, interpolation lineaire directe
  // suivant la coordonnee
  char* v1 = z->_varNames[var1];
  char* v2 = z->_varNames[var2];
  double f = 0.; double* fv = NULL;
  if (strcmp(v1, "CoordinateX") == 0 || strcmp(v1, "x") == 0) 
  { f = xP; fv = fv1; }
  else if (strcmp(v1, "CoordinateY") == 0 || strcmp(v1, "y") == 0) 
  { f = yP; fv = fv1; }
  else if (strcmp(v1, "CoordinateZ") == 0 || strcmp(v1, "z") == 0) 
  { f = zP; fv = fv1; }
  else if (strcmp(v2, "CoordinateX") == 0 || strcmp(v2, "x") == 0) 
  { f = xP; fv = fv2; }
  else if (strcmp(v2, "CoordinateY") == 0 || strcmp(v2, "y") == 0) 
  { f = yP; fv = fv2; }
  else if (strcmp(v2, "CoordinateZ") == 0 || strcmp(v2, "z") == 0) 
  { f = zP; fv = fv2; }
  if (fv != NULL) // une des variables est une coord
  {
    // Check index of f in fv
    E_Int ne = z->_ne;
    int* c = z->_cn;
    int c1, c2;
    for (E_Int i = 0; i < ne; i++)
    {
      c1 = c[i]-1; c2 = c[ne+i]-1;
      //printf("%f %f %f\n", fv[c1], fv[c2], f);
      if (fv[c1] < f+1.e-10 && fv[c2] > f-1.e-10) { e1 = c1; e2 = c2; alpha = (f-fv[c1])/(fv[c2]-fv[c1]); return 1; }
      if (fv[c1] > f-1.e-10 && fv[c2] < f+1.e-10) { e1 = c1; e2 = c2; alpha = (f-fv[c1])/(fv[c2]-fv[c1]); return 1; }
    }
    e1 = 0; e2 = 0; alpha = 0.; // not found
    return 0;
  }

  // autres cas
  int iP = ptrState->activePointI;
  int jP = ptrState->activePointJ;
  int kP = ptrState->activePointK;

  int ind;
  if (nz < _numberOfStructZones) 
  {
    StructZone* zs = (StructZone*)z;
    ind = (iP-1)+(jP-1)*zs->ni+(kP-1)*(zs->ni*zs->nj);
  }
  else ind = iP;

  double f1 = fv1[ind];
  double f2 = fv2[ind];
  
  E_Int ne = z->_ne;
  int* c = z->_cn;
  int c1, c2;
  for (E_Int i = 0; i < ne; i++)
  {
    c1 = c[i]-1; c2 = c[ne+i]-1;
    //printf("%f %f %f\n", fv[c1], fv[c2], f);
    if (fv1[c1] < f1+1.e-10 && fv1[c2] > f1-1.e-10 && fv2[c1] < f2+1.e-10 && fv2[c2] > f2-1.e-10) 
    { e1 = c1; e2 = c2; alpha = (f1-fv1[c1])/(fv1[c2]-fv1[c1]); return 1; }
  }
  e1 = 0; e2 = 0; alpha = 0.; // not found
  return 0;

  e1 = 0; e2 = 0; alpha = 0.; // not found
  return 0; // FAIL
}

//=============================================================================
int Data::display1DActivePoint(Slot1D* s, Zone1D* z, 
                               E_Float posx, E_Float posy,
                               E_Float dx, E_Float dy,
                               int var1, int var2)
{
  int e1, e2;
  double alpha;
  int ret = getActivePointIndex(z, var1, var2, e1, e2, alpha);
  if (ret == 0) return 0; // FAIL

  double* f1 = &(z->_f[var1*(z->_np)]);
  double* f2 = &(z->_f[var2*(z->_np)]);

  // scale
  double borderx = 0.1*dx;
  double bordery = 0.1*dy;
  double r1min = s->_r1min;
  double r1max = s->_r1max;
  double r2min = s->_r2min;
  double r2max = s->_r2max;
  double a1 = (dx-2*borderx)/(r1max-r1min);
  double b1 = (posx+borderx)-a1*r1min;
  double a2 = -(dy-2*bordery)/(r2max-r2min);
  double b2 = (posy+bordery)-a2*r2max;
  
  glDisable(GL_DEPTH_TEST);
  glColor4f(1., 0.11, 0.11, 1.); // foreground
  double x1 = a1*f1[e1]+b1;
  double y1 = a2*f2[e1]+b2;
  double x2 = a1*f1[e2]+b1;
  double y2 = a2*f2[e2]+b2;
  double xp = (1.-alpha)*x1+alpha*x2;
  double yp = (1.-alpha)*y1+alpha*y2;
  //double xp = x1; double yp = y1;

  glBegin(GL_QUADS);
  glVertex3d(xp-5, yp-5, 0.);
  glVertex3d(xp-5, yp+5, 0.);
  glVertex3d(xp+5, yp+5, 0.);
  glVertex3d(xp+5, yp-5, 0.);
  glEnd();
  glEnable(GL_DEPTH_TEST);
  glColor3f(1., 1., 1.);
  return 1;
}

//=============================================================================
int Data::link2View(Zone1D* z, int var1, int var2, 
                    E_Float& r1min, E_Float& r1max, 
                    E_Float& r2min, E_Float& r2max)
{
  // Check coordinates
  // var1 ou var2 doit etre une coordonnee
  char* v1 = z->_varNames[var1];
  char* v2 = z->_varNames[var2];
  int f = 0;
  if (strcmp(v1, "CoordinateX") == 0 || strcmp(v1, "x") == 0) 
  { f = 1; }
  else if (strcmp(v1, "CoordinateY") == 0 || strcmp(v1, "y") == 0) 
  { f = 2; }
  else if (strcmp(v1, "CoordinateZ") == 0 || strcmp(v1, "z") == 0) 
  { f = 3; }
  else if (strcmp(v2, "CoordinateX") == 0 || strcmp(v2, "x") == 0) 
  { f = 4; }
  else if (strcmp(v2, "CoordinateY") == 0 || strcmp(v2, "y") == 0) 
  { f = 5;  }
  else if (strcmp(v2, "CoordinateZ") == 0 || strcmp(v2, "z") == 0) 
  { f = 6; }
  if (f == 0) return 0; // FAIL
  
  // View
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLdouble winX, winY, winZ;
  GLdouble posX, posY, posZ;
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);
  double rmin, rmax;
  rmin = 1.e6; rmax = -1.e6;

  winX = 0.; winY = 0.; winZ = 0.001; // winZ=0, near clipping plane?
  gluUnProject(winX, winY, winZ, modelview, projection, viewport, 
               &posX, &posY, &posZ);
  //printf("%f %f %f\n", posX, posY, posZ);
  switch (f)
  {
    case 1:
    case 4:
      rmin = MIN(rmin, posX); rmax = MAX(rmax, posX);
      break;
    case 2:
    case 5:
      rmin = MIN(rmin, posY); rmax = MAX(rmax, posY);
      break;
    case 3:
    case 6:
      rmin = MIN(rmin, posZ); rmax = MAX(rmax, posZ);
      break;
  }
  winX = _view.w; winY = 0.;
  gluUnProject(winX, winY, winZ, modelview, projection, viewport, 
               &posX, &posY, &posZ);
  //printf("%f %f %f\n", posX, posY, posZ);
  switch (f)
  {
    case 1: 
    case 4:
      rmin = MIN(rmin, posX); rmax = MAX(rmax, posX);
      break;
    case 2:
    case 5:
      rmin = MIN(rmin, posY); rmax = MAX(rmax, posY);
      break;
    case 3:
    case 6:
      rmin = MIN(rmin, posZ); rmax = MAX(rmax, posZ);
      break;
  }
  winX = 0.; winY = _view.h;
  gluUnProject(winX, winY, winZ, modelview, projection, viewport, 
               &posX, &posY, &posZ);
  //printf("%f %f %f\n", posX, posY, posZ);
  switch (f)
  {
    case 1: 
    case 4:
      rmin = MIN(rmin, posX); rmax = MAX(rmax, posX);
      break;
    case 2:
    case 5:
      rmin = MIN(rmin, posY); rmax = MAX(rmax, posY);
      break;
    case 3:
    case 6:
      rmin = MIN(rmin, posZ); rmax = MAX(rmax, posZ);
      break;
  }
  winX = _view.w; winY = _view.h;
  gluUnProject(winX, winY, winZ, modelview, projection, viewport, 
               &posX, &posY, &posZ);
  //printf("%f %f %f\n", posX, posY, posZ);
  switch (f)
  {
    case 1: 
    case 4:
      rmin = MIN(rmin, posX); rmax = MAX(rmax, posX);
      break;
    case 2:
    case 5:
      rmin = MIN(rmin, posY); rmax = MAX(rmax, posY);
      break;
    case 3:
    case 6:
      rmin = MIN(rmin, posZ); rmax = MAX(rmax, posZ);
      break;
  }
  // final
  switch (f)
  {
    case 1:
    case 2:
    case 3:
      r1min = rmin; r1max = rmax; break;
    case 4: 
    case 5:
    case 6:
      r2min = rmin; r2max = rmax; break;
  }
  return 1;
}
