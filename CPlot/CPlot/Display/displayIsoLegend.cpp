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
// Display legend
// IN: dir: 0=vertical, 1=horizontal
//=============================================================================
void Data::displayIsoLegend(int dir)
{
  if (ptrState->mode == MESH) return;
  if (ptrState->mode == SOLID) return;
  if (_numberOfZones < 1) return;
  if (_zones[0]->nfield < 1) return;

  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();

  // Colormap: colors
  float r, g, b;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _pref.colorMap->f;
  double fmax, fmin;
  int niso = ptrState->niso;
  int mod = ptrState->scalarField;
  if (_niso[mod] == -1)
  {
    fmax = maxf[ptrState->scalarField];
    fmin = minf[ptrState->scalarField];
  }
  else
  {
    fmax = _isoMax[ptrState->scalarField];
    fmin = _isoMin[ptrState->scalarField];
  }
  double deltai = 1./(niso*1.);

  int x1 = _view.w-30;
  int y1 = _view.h-2;
  int x2 = _view.w-2;
  int y2 = 50;
  double delta = (y2-y1)/(niso*1.);

  // Rectangle de fond bleu
  /*
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glDisable(GL_DEPTH_TEST);
  glColor4f(0.0, 0.0, 1., 0.3);
  glBegin(GL_QUADS);
  glVertex3d(x1-64, y1+2, -0.1);
  glVertex3d(x2, y1+2, -0.1);
  glVertex3d(x2, y2-20, -0.1);
  glVertex3d(x1-64, y2-20, -0.1);
  glEnd();
  glDisable(GL_BLEND); */

  // Rectangle des couleurs
  glBegin(GL_QUADS);
  for (int i = 0; i < niso; i++)
  {
    getrgb(this, i*deltai, &r, &g, &b);
    glColor4f(r, g, b, 1.);
    glVertex3d(x1, y1+i*delta, 0);
    glVertex3d(x2, y1+i*delta, 0);

    getrgb(this, (i+1)*deltai, &r, &g, &b);  
    glColor4f(r, g, b, 1.);
    glVertex3d(x2, y1+(i+1)*delta, 0);
    glVertex3d(x1, y1+(i+1)*delta, 0);
  }
  glEnd();

  // varname
  char text[128];
  void* font = GLUT_BITMAP_HELVETICA_18;
  //void* font = GLUT_BITMAP_TIMES_ROMAN_24;

  sprintf(text, "%s", _zones[0]->varnames[ptrState->scalarField]);
  int y = y2-5;
  int x = x2 - textWidth(font, text) - 2;
  renderStringWithShadow(x, y, 0, font, text,
                         1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 0.0, 1.0);

  // scale
#define FORMAT "%4.3g"
  sprintf(text, FORMAT, fmin);
  y = y1-5;
  x = x1 - textWidth(font, text) - 2;
  renderStringWithShadow(x, y, 0, font, text,
                         1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 0.0, 1.0);
  sprintf(text, FORMAT, fmax);
  y = y2+17;
  x = x1 - textWidth(font, text) - 2;
  renderStringWithShadow(x, y, 0, font, text,
                         1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 0.0, 1.0);
  sprintf(text, FORMAT, 0.5*(fmax+fmin));
  y = (y1+y2)/2;
  x = x1 - textWidth(font, text) - 2;
  renderStringWithShadow(x, y, 0, font, text,
                         1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 0.0, 1.0);
  sprintf(text, FORMAT, 0.25*fmax+0.75*fmin);
  y = int(0.75*y1+0.25*y2);
  x = x1 - textWidth(font, text) - 2;
  renderStringWithShadow(x, y, 0, font, text,
                         1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 0.0, 1.0);
  sprintf(text, FORMAT, 0.75*fmax+0.25*fmin);
  y = int(0.25*y1+0.75*y2);
  x = x1 - textWidth(font, text) - 2;
  renderStringWithShadow(x, y, 0, font, text,
                         1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 0.0, 1.0);
  glPopMatrix();
  // Put back the previous projection
  resetPerspectiveProjection();
  
  glColor4f(1., 1., 1., 1.);
}
