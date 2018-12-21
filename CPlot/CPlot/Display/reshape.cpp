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
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

//=============================================================================
// Called when window is reshaped
//=============================================================================
void reshape(int w, int h)
{
  Data* d = Data::getInstance();
  d->ptrState->render = 0;
  d->ptrState->syncDisplay( 0 );
  if (d->ptrState->offscreen == 0)
  { d->_view.w = w; d->_view.h = h; }

  // Set view angle here
  d->_view.tang = tan(ANG2RAD * d->_view.angle * 0.5);
  d->farClipping();

  /* pour Ivan
  printf("CPlot reshape: coucou \n");
  double alpha = 0.08;
  double dx = (d->_view.xeye - d->_view.xcam)*alpha;
  double dy = (d->_view.yeye - d->_view.ycam)*alpha;
  double dz = (d->_view.zeye - d->_view.zcam)*alpha;
  double di = sqrt(dx*dx+dy*dy+dz*dz);
  d->adaptiveClipping(di);
  */
}

//=============================================================================
void Data::farClipping()
{
  _view.tang = tan(ANG2RAD * _view.angle*0.5);
  _view.clipping = 0;
  glViewport(0, 0, (GLsizei) _view.w, (GLsizei) _view.h);
  double farD = 40000000.;
  double nearD = 0.5*epsup;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  _view.ratio = (double)_view.w/(double)_view.h;
  if (_view.angle > 0)
    gluPerspective(_view.angle, _view.ratio, nearD, farD);
  else // neg angle means ortho
  {
    // a creuser
    //double f = tan(PIS2 + ANG2RAD*_view.angle*0.5);
    //double t = 1./f; 
    //double r = _view.ratio/f;
    //glOrtho(-r, r, -t, t, nearD, farD);
    //farD = 10000; nearD = -0.1;
    double f = tan(PIS2 + ANG2RAD*_view.angle*0.5);
    double h = f*(farD+nearD)*0.25;
    glOrtho(-h*_view.ratio, h*_view.ratio, -h, h, -nearD, -farD);
  }
  glMatrixMode(GL_MODELVIEW);
  _view.nearD = nearD;
  _view.farD = farD;
  _view.nh = nearD * _view.tang;
  _view.nw = _view.nh * _view.ratio;
  _view.fh = farD * _view.tang;
  _view.fw = _view.fh * _view.ratio;
  ptrState->render = 1;
}

//=============================================================================
void Data::closeClipping()
{
   _view.clipping = 1;
   //double farD = 1000;
   double farD = std::max(10000., 10000.*epsup);
   double nearD = epsup * 5.e-3;
   //printf("clipping close %f %f.\n", nearD, farD);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (_view.angle > 0)
     gluPerspective(_view.angle, _view.ratio, nearD, farD);
   //else gluOrtho();
   glMatrixMode(GL_MODELVIEW);
   _view.nearD = nearD;
   _view.farD = farD;
   _view.nh = nearD * _view.tang;
   _view.nw = _view.nh * _view.ratio;
   _view.fh = farD * _view.tang;
   _view.fw = _view.fh * _view.ratio;
   ptrState->render = 1;
}

//=============================================================================
void Data::veryCloseClipping()
{
   _view.clipping = 2;
   //double farD = 1.;
   double farD = std::max(1000., 1000.*epsup);
   double nearD = epsup * 5.e-4;
   //printf("clipping very close %f %f.\n", nearD, farD);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (_view.angle > 0)
     gluPerspective(_view.angle, _view.ratio, nearD, farD);
   //else gluOtho();
   glMatrixMode(GL_MODELVIEW);
   _view.nearD = nearD;
   _view.farD = farD;
   _view.nh = nearD * _view.tang;
   _view.nw = _view.nh * _view.ratio;
   _view.fh = farD * _view.tang;
   _view.fw = _view.fh * _view.ratio;
   ptrState->render = 1;
}

//=============================================================================
void Data::veryVeryCloseClipping()
{
   _view.clipping = 3;
   //double farD = 0.01;
   double farD = std::max(100., 100.*epsup);
   double nearD = epsup * 5.e-5;
   //printf("clipping very very close.\n");
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (_view.angle > 0)
     gluPerspective(_view.angle, _view.ratio, nearD, farD);
   //else gluOrtho();
   glMatrixMode(GL_MODELVIEW);
   _view.nearD = nearD;
   _view.farD = farD;
   _view.nh = nearD * _view.tang;
   _view.nw = _view.nh * _view.ratio;
   _view.fh = farD * _view.tang;
   _view.fw = _view.fh * _view.ratio;
   ptrState->render = 1;
}

//=============================================================================
void Data::adaptiveClipping(double d)
{
  _view.clipping = -1;
  double nearD, farD;
  nearD = d*0.1; farD = nearD + 40000000.;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (_view.angle > 0)
    gluPerspective(_view.angle, _view.ratio, nearD, farD);
  //else gluOrtho();
  glMatrixMode(GL_MODELVIEW);
  _view.nearD = nearD;
  _view.farD = farD;
  _view.nh = nearD * _view.tang;
  _view.nw = _view.nh * _view.ratio;
  _view.fh = farD * _view.tang;
  _view.fw = _view.fh * _view.ratio;
  ptrState->render = 1;
}
